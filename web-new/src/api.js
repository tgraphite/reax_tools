/**
 * API Client for ReaxTools Web
 */

const API_BASE_URL = 'http://localhost:3000/api';

/**
 * Submit a new analysis job
 */
async function submitJob(file, args = '') {
  const formData = new FormData();
  formData.append('inputFile', file);
  formData.append('args', args);

  const response = await fetch(`${API_BASE_URL}/jobs`, {
    method: 'POST',
    body: formData,
  });

  if (!response.ok) {
    const error = await response.json();
    throw new Error(error.error || 'Failed to submit job');
  }

  return response.json();
}

/**
 * Get job status
 */
async function getJobStatus(jobId) {
  const response = await fetch(`${API_BASE_URL}/jobs/${jobId}`);

  if (!response.ok) {
    const error = await response.json();
    throw new Error(error.error || 'Failed to get job status');
  }

  return response.json();
}

/**
 * Get list of all jobs
 */
async function getJobs(page = 1, limit = 20) {
  const response = await fetch(`${API_BASE_URL}/jobs?page=${page}&limit=${limit}`);

  if (!response.ok) {
    const error = await response.json();
    throw new Error(error.error || 'Failed to get jobs');
  }

  return response.json();
}

/**
 * Delete a job
 */
async function deleteJob(jobId) {
  const response = await fetch(`${API_BASE_URL}/jobs/${jobId}`, {
    method: 'DELETE',
  });

  if (!response.ok) {
    const error = await response.json();
    throw new Error(error.error || 'Failed to delete job');
  }

  return response.json();
}

/**
 * Get specific file content from job results
 */
async function getJobFile(jobId, filename) {
  const response = await fetch(`${API_BASE_URL}/jobs/${jobId}/files/${filename}`);

  if (!response.ok) {
    const error = await response.json();
    throw new Error(error.error || 'Failed to get file');
  }

  return response.json();
}

/**
 * Download job results as ZIP
 */
function downloadResults(jobId) {
  window.open(`${API_BASE_URL}/jobs/${jobId}/download`, '_blank');
}

/**
 * Create SSE connection for real-time job updates
 */
function createJobStream(jobId, callbacks) {
  const { onProgress, onLog, onComplete, onError, onStatus } = callbacks;
  const eventSource = new EventSource(`${API_BASE_URL}/jobs/${jobId}/stream`);

  eventSource.addEventListener('status', (e) => {
    const data = JSON.parse(e.data);
    if (onStatus) onStatus(data);
  });

  eventSource.addEventListener('progress', (e) => {
    const data = JSON.parse(e.data);
    if (onProgress) onProgress(data.progress);
  });

  eventSource.addEventListener('log', (e) => {
    const data = JSON.parse(e.data);
    if (onLog) onLog(data.log);
  });

  eventSource.addEventListener('complete', (e) => {
    const data = JSON.parse(e.data);
    if (onComplete) onComplete(data.result);
    eventSource.close();
  });

  eventSource.addEventListener('error', (e) => {
    const data = JSON.parse(e.data);
    if (onError) onError(data.error || 'Stream error');
    eventSource.close();
  });

  // Handle connection errors
  eventSource.onerror = (e) => {
    if (eventSource.readyState === EventSource.CLOSED) {
      if (onError) onError('Connection closed');
    }
  };

  return {
    close: () => eventSource.close(),
  };
}

/**
 * Load sample file
 */
async function loadSampleFile() {
  const response = await fetch('../web/reaxff.xyz');
  if (!response.ok) {
    throw new Error('Failed to load sample file');
  }
  const blob = await response.blob();
  return new File([blob], 'reaxff.xyz', { type: 'text/plain' });
}

// Export API functions
window.ReaxAPI = {
  submitJob,
  getJobStatus,
  getJobs,
  deleteJob,
  getJobFile,
  downloadResults,
  createJobStream,
  loadSampleFile,
};
