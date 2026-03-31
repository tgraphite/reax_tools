const express = require('express');
const cors = require('cors');
const compression = require('compression');
const path = require('path');
require('dotenv').config();

const jobRoutes = require('./routes/jobs');
const { initializeQueue } = require('./services/queue');
const { startCleanupScheduler } = require('./utils/fileCleanup');

const app = express();
const PORT = process.env.PORT || 3000;

// Middleware
app.use(cors());
app.use(compression());
app.use(express.json());

// Static files for uploads and outputs (temporary access)
app.use('/temp', express.static(path.join(__dirname, '../temp')));
app.use('/outputs', express.static(path.join(__dirname, '../temp/outputs')));

// Routes
app.use('/api/jobs', jobRoutes);

// Health check
app.get('/api/health', (req, res) => {
  res.json({ status: 'ok', timestamp: new Date().toISOString() });
});

// Error handler
app.use((err, req, res, next) => {
  console.error('Error:', err);
  res.status(500).json({ error: err.message || 'Internal server error' });
});

// Initialize services
async function startServer() {
  try {
    // Initialize BullMQ queue
    await initializeQueue();
    console.log('Queue initialized');

    // Start cleanup scheduler
    startCleanupScheduler();
    console.log('Cleanup scheduler started');

    // Start server
    app.listen(PORT, () => {
      console.log(`ReaxTools Server running on port ${PORT}`);
      console.log(`API available at http://localhost:${PORT}/api`);
    });
  } catch (error) {
    console.error('Failed to start server:', error);
    process.exit(1);
  }
}

startServer();
