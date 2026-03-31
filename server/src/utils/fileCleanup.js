const fs = require('fs').promises;
const path = require('path');

const UPLOAD_DIR = path.join(__dirname, '../../temp/uploads');
const OUTPUT_DIR = path.join(__dirname, '../../temp/outputs');

// Cleanup files older than 24 hours
async function cleanupOldFiles() {
  const maxAge = 24 * 60 * 60 * 1000; // 24 hours
  const now = Date.now();

  async function cleanDirectory(dir) {
    try {
      const entries = await fs.readdir(dir, { withFileTypes: true });

      for (const entry of entries) {
        const fullPath = path.join(dir, entry.name);

        try {
          const stats = await fs.stat(fullPath);
          const age = now - stats.mtime.getTime();

          if (age > maxAge) {
            if (entry.isDirectory()) {
              await fs.rm(fullPath, { recursive: true, force: true });
              console.log(`Cleaned up old directory: ${fullPath}`);
            } else {
              await fs.unlink(fullPath);
              console.log(`Cleaned up old file: ${fullPath}`);
            }
          }
        } catch (err) {
          console.error(`Error checking ${fullPath}:`, err.message);
        }
      }
    } catch (err) {
      console.error(`Error reading directory ${dir}:`, err.message);
    }
  }

  await cleanDirectory(UPLOAD_DIR);
  await cleanDirectory(OUTPUT_DIR);
}

// Start periodic cleanup
function startCleanupScheduler() {
  // Run immediately on start
  cleanupOldFiles().catch(console.error);

  // Then every hour
  setInterval(() => {
    cleanupOldFiles().catch(console.error);
  }, 60 * 60 * 1000);

  console.log('File cleanup scheduler started (runs every hour)');
}

module.exports = {
  cleanupOldFiles,
  startCleanupScheduler,
};
