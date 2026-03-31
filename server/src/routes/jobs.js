const express = require('express');
const multer = require('multer');
const path = require('path');
const fs = require('fs').promises;
const { v4: uuidv4 } = require('uuid');
const { getQueue } = require('../services/queue');

const router = express.Router();

// Ensure temp directories exist
const UPLOAD_DIR = path.join(__dirname, '../../temp/uploads');
const OUTPUT_DIR = path.join(__dirname, '../../temp/outputs');

(async () => {
  await fs.mkdir(UPLOAD_DIR, { recursive: true });
  await fs.mkdir(OUTPUT_DIR, { recursive: true });
})();

// Configure multer for file uploads
const storage = multer.diskStorage({
  destination: UPLOAD_DIR,
  filename: (req, file, cb) => {
    const uniqueName = `${uuidv4()}-${file.originalname}`;
    cb(null, uniqueName);
  },
});

const upload = multer({
  storage,
  limits: {
    fileSize: 500 * 1024 * 1024, // 500MB limit
  },
  fileFilter: (req, file, cb) => {
    // Accept .xyz and .lammpstrj files
    const ext = path.extname(file.originalname).toLowerCase();
    if (ext === '.xyz' || ext === '.lammpstrj' || ext === '.trj') {
      cb(null, true);
    } else {
      cb(new Error('Only .xyz and .lammpstrj files are allowed'));
    }
  },
});

// Submit a new job
router.post('/', upload.single('inputFile'), async (req, res) => {
  try {
    if (!req.file) {
      return res.status(400).json({ error: 'No file uploaded' });
    }

    let { args = '' } = req.body;

    // Reject manual -nt specification
    if (/(?:^|\s)-nt\s*\d+(?=\s|$)/i.test(args)) {
      // Clean up uploaded file before returning error
      try { await fs.unlink(req.file.path); } catch {}
      return res.status(400).json({ error: '网页版不允许指定核数，请不要包含 -nt <number> 参数' });
    }

    // Force -nt 1
    args = args ? `${args.trim()} -nt 1` : '-nt 1';

    const jobId = uuidv4();
    const outputDir = path.join(OUTPUT_DIR, jobId);
    await fs.mkdir(outputDir, { recursive: true });

    const queue = getQueue();
    const job = await queue.add(
      'analysis',
      {
        inputFile: req.file.path,
        originalName: req.file.originalname,
        args,
        outputDir,
        uploadedAt: new Date().toISOString(),
      },
      {
        jobId,
      }
    );

    // Get queue status for better UX
    const [waitingCount, activeCount] = await Promise.all([
      queue.getWaitingCount(),
      queue.getActiveCount(),
    ]);

    res.status(201).json({
      jobId: job.id,
      status: 'queued',
      position: waitingCount,
      runningCount: activeCount,
      message: 'Job submitted successfully',
    });
  } catch (error) {
    console.error('Error submitting job:', error);
    res.status(500).json({ error: error.message });
  }
});

// Get job status
router.get('/:id', async (req, res) => {
  try {
    const { id } = req.params;
    const queue = getQueue();
    const job = await queue.getJob(id);

    if (!job) {
      return res.status(404).json({ error: 'Job not found' });
    }

    const state = await job.getState();
    const progress = job.progress || 0;

    res.json({
      jobId: job.id,
      status: state,
      progress,
      log: job.data.log || '',
      originalName: job.data.originalName,
      args: job.data.args,
      createdAt: job.timestamp,
      processedAt: job.processedOn,
      completedAt: job.finishedOn,
      result: job.returnvalue,
      failedReason: job.failedReason,
    });
  } catch (error) {
    console.error('Error getting job status:', error);
    res.status(500).json({ error: error.message });
  }
});

// SSE stream for real-time updates
router.get('/:id/stream', async (req, res) => {
  try {
    const { id } = req.params;
    const queue = getQueue();
    const job = await queue.getJob(id);

    if (!job) {
      return res.status(404).json({ error: 'Job not found' });
    }

    // Set up SSE headers
    res.setHeader('Content-Type', 'text/event-stream');
    res.setHeader('Cache-Control', 'no-cache');
    res.setHeader('Connection', 'keep-alive');
    res.setHeader('X-Accel-Buffering', 'no'); // Disable nginx buffering

    // Send initial state
    const state = await job.getState();
    res.write(`event: status\ndata: ${JSON.stringify({ status: state, progress: job.progress || 0 })}\n\n`);

    // Set up interval to check progress
    let lastProgress = job.progress || 0;
    let lastLog = '';
    const logFilePath = path.join(job.data.outputDir, 'job.log');

    const interval = setInterval(async () => {
      try {
        const updatedJob = await queue.getJob(id);
        if (!updatedJob) {
          clearInterval(interval);
          res.write(`event: error\ndata: ${JSON.stringify({ message: 'Job disappeared' })}\n\n`);
          res.end();
          return;
        }

        const currentState = await updatedJob.getState();
        const currentProgress = updatedJob.progress || 0;

        // Read log from disk if available, else fallback to job data
        let currentLog = updatedJob.data.log || '';
        try {
          currentLog = await fs.readFile(logFilePath, 'utf-8');
        } catch {
          // file not ready yet, use job data fallback
        }

        // Send progress update if changed
        if (currentProgress !== lastProgress) {
          res.write(`event: progress\ndata: ${JSON.stringify({ progress: currentProgress })}\n\n`);
          lastProgress = currentProgress;
        }

        // Send log update if changed
        if (currentLog !== lastLog) {
          let newLog;
          if (currentLog.length < lastLog.length) {
            newLog = currentLog;
          } else {
            newLog = currentLog.slice(lastLog.length);
          }
          if (newLog) {
            res.write(`event: log\ndata: ${JSON.stringify({ log: newLog })}\n\n`);
          }
          lastLog = currentLog;
        }

        // Check if job is completed or failed
        if (currentState === 'completed') {
          res.write(`event: complete\ndata: ${JSON.stringify({ result: updatedJob.returnvalue })}\n\n`);
          clearInterval(interval);
          res.end();
        } else if (currentState === 'failed') {
          res.write(`event: error\ndata: ${JSON.stringify({ error: updatedJob.failedReason })}\n\n`);
          clearInterval(interval);
          res.end();
        }
      } catch (err) {
        console.error('SSE error:', err);
        clearInterval(interval);
        res.end();
      }
    }, 500); // Check every 500ms

    // Clean up on client disconnect
    req.on('close', () => {
      clearInterval(interval);
      res.end();
    });
  } catch (error) {
    console.error('Error setting up SSE:', error);
    res.status(500).json({ error: error.message });
  }
});

// Download results
router.get('/:id/download', async (req, res) => {
  try {
    const { id } = req.params;
    const queue = getQueue();
    const job = await queue.getJob(id);

    if (!job) {
      return res.status(404).json({ error: 'Job not found' });
    }

    const state = await job.getState();
    if (state !== 'completed') {
      return res.status(400).json({ error: 'Job not completed yet' });
    }

    const zipPath = path.join(job.data.outputDir, `results-${id}.zip`);

    // Check if file exists
    try {
      await fs.access(zipPath);
    } catch {
      return res.status(404).json({ error: 'Result file not found' });
    }

    res.setHeader('Content-Type', 'application/zip');
    res.setHeader('Content-Disposition', `attachment; filename="reax-tools-results-${id}.zip"`);
    res.sendFile(zipPath);
  } catch (error) {
    console.error('Error downloading results:', error);
    res.status(500).json({ error: error.message });
  }
});

// Get specific output file content
router.get('/:id/files/:filename', async (req, res) => {
  try {
    const { id, filename } = req.params;
    const queue = getQueue();
    const job = await queue.getJob(id);

    if (!job) {
      return res.status(404).json({ error: 'Job not found' });
    }

    const filePath = path.join(job.data.outputDir, filename);

    // Security check: ensure file is within output directory
    const resolvedPath = path.resolve(filePath);
    const resolvedOutputDir = path.resolve(job.data.outputDir);
    if (!resolvedPath.startsWith(resolvedOutputDir)) {
      return res.status(403).json({ error: 'Access denied' });
    }

    try {
      const content = await fs.readFile(filePath, 'utf-8');
      res.json({ filename, content });
    } catch {
      res.status(404).json({ error: 'File not found' });
    }
  } catch (error) {
    console.error('Error getting file:', error);
    res.status(500).json({ error: error.message });
  }
});

// List all jobs (with pagination)
router.get('/', async (req, res) => {
  try {
    const queue = getQueue();
    const { page = 1, limit = 20 } = req.query;

    // Get jobs from different states
    const [completed, failed, waiting, active] = await Promise.all([
      queue.getCompleted(0, 100),
      queue.getFailed(0, 100),
      queue.getWaiting(0, 100),
      queue.getActive(0, 100),
    ]);

    const allJobs = [...active, ...waiting, ...completed, ...failed];

    // Sort by timestamp desc
    allJobs.sort((a, b) => b.timestamp - a.timestamp);

    // Paginate
    const start = (page - 1) * limit;
    const paginatedJobs = allJobs.slice(start, start + parseInt(limit));

    const jobs = await Promise.all(
      paginatedJobs.map(async (job) => ({
        jobId: job.id,
        status: await job.getState(),
        progress: job.progress || 0,
        originalName: job.data.originalName,
        args: job.data.args,
        createdAt: job.timestamp,
        completedAt: job.finishedOn,
      }))
    );

    res.json({
      jobs,
      pagination: {
        page: parseInt(page),
        limit: parseInt(limit),
        total: allJobs.length,
      },
    });
  } catch (error) {
    console.error('Error listing jobs:', error);
    res.status(500).json({ error: error.message });
  }
});

// Cancel/delete job
router.delete('/:id', async (req, res) => {
  try {
    const { id } = req.params;
    const queue = getQueue();
    const job = await queue.getJob(id);

    if (!job) {
      return res.status(404).json({ error: 'Job not found' });
    }

    // Remove from queue
    await job.remove();

    // Clean up files
    try {
      await fs.rm(job.data.outputDir, { recursive: true, force: true });
      await fs.unlink(job.data.inputFile).catch(() => {});
    } catch (err) {
      console.error('Error cleaning up files:', err);
    }

    res.json({ message: 'Job removed' });
  } catch (error) {
    console.error('Error removing job:', error);
    res.status(500).json({ error: error.message });
  }
});

module.exports = router;
