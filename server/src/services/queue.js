const { Queue, Worker, Job } = require('bullmq');
const IORedis = require('ioredis');
const path = require('path');
const { spawn, exec } = require('child_process');
const util = require('util');
const fs = require('fs').promises;
const fsSync = require('fs');
const JSZip = require('jszip');

const execPromise = util.promisify(exec);

// Redis connection
const redisConnection = new IORedis({
  host: process.env.REDIS_HOST || 'localhost',
  port: process.env.REDIS_PORT || 6379,
  maxRetriesPerRequest: null,
  enableReadyCheck: false,
});

// Queue instance
let analysisQueue;

// Initialize queue
function initializeQueue() {
  analysisQueue = new Queue('analysis', {
    connection: redisConnection,
    defaultJobOptions: {
      attempts: 2,
      backoff: {
        type: 'exponential',
        delay: 5000,
      },
      removeOnComplete: {
        age: 24 * 3600, // Keep completed jobs for 24 hours
        count: 100,
      },
      removeOnFail: {
        age: 7 * 24 * 3600, // Keep failed jobs for 7 days
      },
    },
  });

  // Start worker
  startWorker();

  return analysisQueue;
}

// Get queue instance
function getQueue() {
  if (!analysisQueue) {
    throw new Error('Queue not initialized');
  }
  return analysisQueue;
}

// Estimate total frames for progress tracking
async function getTotalFrames(filePath) {
  try {
    const ext = path.extname(filePath).toLowerCase();
    if (ext === '.xyz') {
      const fileHandle = await fs.open(filePath, 'r');
      const buf = Buffer.alloc(64);
      const { bytesRead } = await fileHandle.read(buf, 0, 64, 0);
      await fileHandle.close();
      const firstLine = buf.toString('utf8', 0, bytesRead).split('\n')[0].trim();
      const atomCount = parseInt(firstLine);
      if (isNaN(atomCount) || atomCount <= 0) return 0;
      const { stdout } = await execPromise(`wc -l < "${filePath}"`);
      const totalLines = parseInt(stdout.trim());
      return Math.max(1, Math.floor(totalLines / (atomCount + 2)));
    } else if (ext === '.lammpstrj' || ext === '.trj') {
      const { stdout } = await execPromise(`grep -c 'ITEM: TIMESTEP' "${filePath}"`);
      return parseInt(stdout.trim()) || 0;
    }
  } catch (err) {
    console.warn('[Frame Count] Failed:', err.message);
  }
  return 0;
}

// Start worker process
function startWorker() {
  const concurrency = 1; // Web版只允许单任务计算
  const timeout = parseInt(process.env.JOB_TIMEOUT) || 5 * 60 * 1000; // 5 minutes

  const worker = new Worker(
    'analysis',
    async (job) => {
      const { inputFile, args, outputDir } = job.data;
      const binaryPath = process.env.REAX_TOOLS_BINARY || path.join(__dirname, '../../../bin/reax_tools');

      // Pre-compute total frames for progress tracking
      const totalFrames = await getTotalFrames(inputFile);
      const logFilePath = path.join(outputDir, 'job.log');
      try { await fs.writeFile(logFilePath, ''); } catch {}

      return new Promise((resolve, reject) => {
        const startTime = Date.now();
        let stdout = '';
        let stderr = '';
        let latestFrame = 0;

        // Parse args
        const argArray = args ? args.split(/\s+/).filter(Boolean) : [];
        const processArgs = ['-o', outputDir, '-f', inputFile, ...argArray];

        console.log(`[Job ${job.id}] Starting: ${binaryPath} ${processArgs.join(' ')}`);

        // Use stdbuf to force line-buffered stdout for real-time streaming
        let spawnCmd = binaryPath;
        let spawnArgs = processArgs;
        if (process.platform !== 'win32') {
          try {
            fsSync.accessSync('/usr/bin/stdbuf', fsSync.constants.X_OK);
            spawnCmd = 'stdbuf';
            spawnArgs = ['-oL', binaryPath, ...processArgs];
          } catch {
            // stdbuf not available, fall back to direct spawn
          }
        }

        const child = spawn(spawnCmd, spawnArgs, {
          timeout: timeout,
          killSignal: 'SIGTERM',
        });

        child.stdout.on('data', (data) => {
          const text = data.toString();
          stdout += text;

          // Write to disk log for SSE streaming
          try { fsSync.appendFileSync(logFilePath, text); } catch {}

          // Parse frame progress from reax_tools table output
          // Matches lines like "1        2592     2146     446 ..."
          const frameRegex = /^\s*(\d+)\s+\d+/gm;
          let match;
          while ((match = frameRegex.exec(text)) !== null) {
            const frameNum = parseInt(match[1]);
            if (frameNum > latestFrame) latestFrame = frameNum;
          }

          if (latestFrame > 0 && totalFrames > 0) {
            const progress = Math.min(Math.round((latestFrame / totalFrames) * 100), 99);
            job.updateProgress(progress).catch(() => {});
          }

          // Keep reduced log in job data for completed summary
          job.updateData({
            ...job.data,
            log: stdout.slice(-4000),
          }).catch(() => {});
        });

        child.stderr.on('data', (data) => {
          const text = data.toString();
          stderr += text;
          try { fsSync.appendFileSync(logFilePath, `[stderr] ${text}`); } catch {}
        });

        child.on('close', async (code) => {
          const duration = Date.now() - startTime;

          if (code === 0) {
            try {
              await job.updateProgress(100).catch(() => {});
              // Generate static images using reax_plot.py and Graphviz
              await generateResultImages(outputDir, job.id);

              // Create ZIP of results
              const zipPath = await createResultZip(outputDir, job.id);

              // Collect generated image names
              const allFiles = await fs.readdir(outputDir);
              const images = allFiles.filter(f => f.endsWith('.png'));

              resolve({
                success: true,
                duration,
                outputDir,
                zipPath,
                stdout,
                images,
              });
            } catch (zipError) {
              reject(new Error(`Analysis completed but failed to create ZIP: ${zipError.message}`));
            }
          } else {
            reject(new Error(`Process exited with code ${code}: ${stderr || stdout}`));
          }
        });

        child.on('error', (error) => {
          reject(new Error(`Failed to start process: ${error.message}`));
        });
      });
    },
    {
      connection: redisConnection,
      concurrency: concurrency,
    }
  );

  worker.on('completed', (job, result) => {
    console.log(`[Job ${job.id}] Completed in ${result.duration}ms`);
  });

  worker.on('failed', (job, err) => {
    console.error(`[Job ${job.id}] Failed:`, err.message);
  });

  console.log(`Worker started with concurrency=${concurrency}, timeout=${timeout}ms`);
}

// Create ZIP file of results
async function createResultZip(outputDir, jobId) {
  const zip = new JSZip();

  // Read all files in output directory
  const files = await fs.readdir(outputDir);

  for (const file of files) {
    const filePath = path.join(outputDir, file);
    const stat = await fs.stat(filePath);

    if (stat.isFile()) {
      const content = await fs.readFile(filePath);
      zip.file(file, content);
    }
  }

  const zipPath = path.join(outputDir, `results-${jobId}.zip`);
  const zipContent = await zip.generateAsync({ type: 'nodebuffer' });
  await fs.writeFile(zipPath, zipContent);

  return zipPath;
}

// Generate result images using reax_plot.py and Graphviz dot
async function generateResultImages(outputDir, jobId) {
  const plotScript = process.env.REAX_PLOT_SCRIPT || path.join(__dirname, '../../../bin/reax_plot.py');

  // 1. Run reax_plot.py for CSV charts
  try {
    console.log(`[Job ${jobId}] Generating plots with reax_plot.py...`);
    await execPromise(`python3 "${plotScript}" -d "${outputDir}"`, { timeout: 60000 });
  } catch (err) {
    console.warn(`[Job ${jobId}] reax_plot.py failed:`, err.message);
  }

  // 2. Render DOT files with Graphviz
  const dotFiles = ['reaction_flow_full.dot', 'reaction_flow_simplified.dot'];
  for (const dotFile of dotFiles) {
    const dotPath = path.join(outputDir, dotFile);
    try {
      await fs.access(dotPath);
      const baseName = path.basename(dotFile, '.dot');
      const outPath = path.join(outputDir, `${baseName}.png`);
      console.log(`[Job ${jobId}] Rendering ${dotFile} with Graphviz...`);
      await execPromise(`dot -Tpng "${dotPath}" -o "${outPath}"`, { timeout: 60000 });
    } catch {
      // File doesn't exist, skip
    }
  }
}

module.exports = {
  initializeQueue,
  getQueue,
  redisConnection,
};
