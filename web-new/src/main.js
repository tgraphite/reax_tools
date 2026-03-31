/**
 * Main application logic for ReaxTools Web Dashboard
 */

// Global state
const state = {
  currentJobId: null,
  jobStream: null,
  uploadedFile: null,
  jobs: [],
};

// DOM Elements
const elements = {};

/**
 * Initialize application
 */
document.addEventListener('DOMContentLoaded', () => {
  // Cache DOM elements
  cacheElements();

  // Initialize theme
  initTheme();

  // Setup event listeners
  setupEventListeners();

  // Load job history
  loadJobHistory();

  // Start queue status polling
  startQueuePolling();

  console.log('ReaxTools Web Dashboard initialized');
});

/**
 * Cache DOM elements
 */
function cacheElements() {
  elements.fileInput = document.getElementById('fileInput');
  elements.dropzone = document.getElementById('dropzone');
  elements.dropzoneContent = document.getElementById('dropzoneContent');
  elements.fileSelected = document.getElementById('fileSelected');
  elements.selectedFilename = document.getElementById('selectedFilename');
  elements.changeFile = document.getElementById('changeFile');
  elements.argsInput = document.getElementById('argsInput');
  elements.submitForm = document.getElementById('submitForm');
  elements.submitBtn = document.getElementById('submitBtn');
  elements.logViewer = document.getElementById('logViewer');
  elements.progressContainer = document.getElementById('progressContainer');
  elements.progressBar = document.getElementById('progressBar');
  elements.progressText = document.getElementById('progressText');
  elements.jobStatus = document.getElementById('jobStatus');
  elements.statusBadge = document.getElementById('statusBadge');
  elements.resultsActions = document.getElementById('resultsActions');
  elements.downloadBtn = document.getElementById('downloadBtn');
  elements.viewResultsBtn = document.getElementById('viewResultsBtn');
  elements.chartsSection = document.getElementById('chartsSection');
  elements.imagesGrid = document.getElementById('imagesGrid');
  elements.jobHistoryTable = document.getElementById('jobHistoryTable');
  elements.queueCount = document.getElementById('queueCount');
  elements.alertContainer = document.getElementById('alertContainer');
  elements.darkModeToggle = document.getElementById('darkModeToggle');
}

/**
 * Initialize dark mode
 */
function initTheme() {
  const savedTheme = localStorage.getItem('darkMode');
  if (savedTheme === 'true' || (!savedTheme && window.matchMedia('(prefers-color-scheme: dark)').matches)) {
    document.documentElement.classList.add('dark');
  }
}

/**
 * Toggle dark mode
 */
function toggleDarkMode() {
  const isDark = document.documentElement.classList.toggle('dark');
  localStorage.setItem('darkMode', isDark);
}

/**
 * Setup event listeners
 */
function setupEventListeners() {
  // Dark mode toggle
  elements.darkModeToggle?.addEventListener('click', toggleDarkMode);

  // File dropzone
  elements.dropzone?.addEventListener('click', () => elements.fileInput?.click());
  elements.dropzone?.addEventListener('dragover', handleDragOver);
  elements.dropzone?.addEventListener('dragleave', handleDragLeave);
  elements.dropzone?.addEventListener('drop', handleDrop);
  elements.fileInput?.addEventListener('change', handleFileSelect);
  elements.changeFile?.addEventListener('click', (e) => {
    e.stopPropagation();
    resetFileSelection();
  });

  // Form submission
  elements.submitForm?.addEventListener('submit', handleSubmit);

  // Results actions
  elements.viewResultsBtn?.addEventListener('click', scrollToCharts);
  elements.downloadBtn?.addEventListener('click', handleDownload);

}

/**
 * Handle drag over
 */
function handleDragOver(e) {
  e.preventDefault();
  elements.dropzone?.classList.add('active');
}

/**
 * Handle drag leave
 */
function handleDragLeave(e) {
  e.preventDefault();
  elements.dropzone?.classList.remove('active');
}

/**
 * Handle file drop
 */
function handleDrop(e) {
  e.preventDefault();
  elements.dropzone?.classList.remove('active');

  const files = e.dataTransfer?.files;
  if (files?.length > 0) {
    selectFile(files[0]);
  }
}

/**
 * Handle file select
 */
function handleFileSelect(e) {
  const file = e.target.files?.[0];
  if (file) {
    selectFile(file);
  }
}

/**
 * Select a file
 */
function selectFile(file) {
  const validExtensions = ['.xyz', '.lammpstrj', '.trj'];
  const ext = '.' + file.name.split('.').pop().toLowerCase();

  if (!validExtensions.includes(ext)) {
    showAlert('Invalid file type. Please upload .xyz or .lammpstrj files.', 'error');
    return;
  }

  state.uploadedFile = file;
  elements.selectedFilename.textContent = `${file.name} (${formatFileSize(file.size)})`;
  elements.dropzoneContent.classList.add('hidden');
  elements.fileSelected.classList.remove('hidden');
  elements.submitBtn.disabled = false;
}

/**
 * Reset file selection
 */
function resetFileSelection() {
  state.uploadedFile = null;
  elements.fileInput.value = '';
  elements.dropzoneContent.classList.remove('hidden');
  elements.fileSelected.classList.add('hidden');
  elements.submitBtn.disabled = true;
}

/**
 * Format file size
 */
function formatFileSize(bytes) {
  if (bytes === 0) return '0 Bytes';
  const k = 1024;
  const sizes = ['Bytes', 'KB', 'MB', 'GB'];
  const i = Math.floor(Math.log(bytes) / Math.log(k));
  return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
}

/**
 * Handle form submission
 */
async function handleSubmit(e) {
  e.preventDefault();

  if (!state.uploadedFile) {
    showAlert('Please select a file first', 'error');
    return;
  }

  try {
    // Reset UI
    resetJobUI();
    elements.submitBtn.disabled = true;
    elements.submitBtn.innerHTML = `
      <svg class="animate-spin -ml-1 mr-2 h-5 w-5 text-white" fill="none" viewBox="0 0 24 24">
        <circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle>
        <path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
      </svg>
      提交中...
    `;

    // Submit job
    const args = elements.argsInput.value.trim();
    const result = await ReaxAPI.submitJob(state.uploadedFile, args);

    state.currentJobId = result.jobId;
    let queueMsg = `任务已提交!`;
    if (result.runningCount > 0) {
      queueMsg += ` 当前有 ${result.runningCount} 个任务正在计算，前方排队: ${result.position}`;
    } else if (result.position > 0) {
      queueMsg += ` 前方排队: ${result.position}`;
    } else {
      queueMsg += ` 即将开始计算`;
    }
    showAlert(queueMsg, 'success');

    // Start monitoring
    startJobMonitoring(result.jobId);

    // Refresh job history
    loadJobHistory();

  } catch (error) {
    showAlert(`提交失败: ${error.message}`, 'error');
    elements.submitBtn.disabled = false;
    elements.submitBtn.innerHTML = `
      <svg class="w-5 h-5 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M14.752 11.168l-3.197-2.132A1 1 0 0010 9.87v4.263a1 1 0 001.555.832l3.197-2.132a1 1 0 000-1.664z"></path>
        <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M21 12a9 9 0 11-18 0 9 9 0 0118 0z"></path>
      </svg>
      开始分析
    `;
  }
}

/**
 * Reset job UI
 */
function resetJobUI() {
  elements.logViewer.innerHTML = '';
  elements.progressContainer.classList.add('hidden');
  elements.progressBar.style.width = '0%';
  elements.progressText.textContent = '0%';
  elements.jobStatus.classList.add('hidden');
  elements.resultsActions.classList.add('hidden');
  elements.chartsSection.classList.add('hidden');
  ReaxCharts.clearCharts();

  if (state.jobStream) {
    state.jobStream.close();
    state.jobStream = null;
  }
}

/**
 * Start monitoring a job
 */
function startJobMonitoring(jobId) {
  elements.jobStatus.classList.remove('hidden');
  elements.progressContainer.classList.remove('hidden');

  state.jobStream = ReaxAPI.createJobStream(jobId, {
    onStatus: (data) => {
      updateStatusBadge(data.status);
    },
    onProgress: (progress) => {
      elements.progressBar.style.width = `${progress}%`;
      elements.progressText.textContent = `${progress}%`;
    },
    onLog: (log) => {
      appendLog(log);
    },
    onComplete: (result) => {
      handleJobComplete(jobId, result);
    },
    onError: (error) => {
      handleJobError(error);
    },
  });
}

/**
 * Update status badge
 */
function updateStatusBadge(status) {
  const badgeClasses = {
    waiting: 'badge-queued',
    queued: 'badge-queued',
    active: 'badge-running',
    running: 'badge-running',
    completed: 'badge-completed',
    failed: 'badge-failed',
  };

  const badgeTexts = {
    waiting: '等待中',
    queued: '等待中',
    active: '运行中',
    running: '运行中',
    completed: '已完成',
    failed: '失败',
  };

  elements.statusBadge.className = `badge ${badgeClasses[status] || 'badge-queued'}`;
  elements.statusBadge.textContent = badgeTexts[status] || status;
}

/**
 * Append log message
 */
function appendLog(message) {
  // Convert newlines to <br> tags for proper rendering
  const lines = message.split('\n');
  lines.forEach(line => {
    if (line.trim()) {
      const div = document.createElement('div');
      div.className = 'text-gray-300 mb-0.5 font-mono text-sm';
      div.textContent = line;
      elements.logViewer.appendChild(div);
    }
  });
  elements.logViewer.scrollTop = elements.logViewer.scrollHeight;
}

/**
 * Handle job completion
 */
async function handleJobComplete(jobId, result) {
  elements.submitBtn.disabled = false;
  elements.submitBtn.innerHTML = `
    <svg class="w-5 h-5 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M14.752 11.168l-3.197-2.132A1 1 0 0010 9.87v4.263a1 1 0 001.555.832l3.197-2.132a1 1 0 000-1.664z"></path>
      <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M21 12a9 9 0 11-18 0 9 9 0 0118 0z"></path>
    </svg>
    开始分析
  `;

  elements.progressBar.style.width = '100%';
  elements.progressText.textContent = '100%';
  updateStatusBadge('completed');

  elements.resultsActions.classList.remove('hidden');
  elements.downloadBtn.href = `http://localhost:3000/api/jobs/${jobId}/download`;

  showAlert('分析完成!', 'success');
  loadJobHistory();

  // Auto-load results
  await loadResults(jobId);
}

/**
 * Handle job error
 */
function handleJobError(error) {
  elements.submitBtn.disabled = false;
  elements.submitBtn.innerHTML = `
    <svg class="w-5 h-5 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M14.752 11.168l-3.197-2.132A1 1 0 0010 9.87v4.263a1 1 0 001.555.832l3.197-2.132a1 1 0 000-1.664z"></path>
      <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M21 12a9 9 0 11-18 0 9 9 0 0118 0z"></path>
    </svg>
    开始分析
  `;

  updateStatusBadge('failed');
  showAlert(`分析失败: ${error}`, 'error');
  loadJobHistory();
}

/**
 * Load and display results (static images)
 */
async function loadResults(jobId) {
  try {
    ReaxCharts.clearCharts();

    // Get job info to know generated images
    const jobInfo = await ReaxAPI.getJobStatus(jobId);
    const images = jobInfo.result?.images || [];
    const baseUrl = `http://localhost:3000/outputs/${jobId}`;

    // Sort images for consistent order
    images.sort();

    if (images.length === 0) {
      elements.imagesGrid.innerHTML = `
        <div class="col-span-full text-center py-8 text-gray-500">
          未生成可视化图片
        </div>
      `;
    } else {
      elements.imagesGrid.innerHTML = images.map(img => {
        const title = img.replace(/\.png$/i, '');
        return `
          <div class="dashboard-card p-4">
            <h3 class="font-medium text-gray-900 dark:text-white mb-3 text-sm">${title}</h3>
            <img src="${baseUrl}/${img}" alt="${title}" class="w-full rounded border border-gray-200 dark:border-gray-700">
          </div>
        `;
      }).join('');
    }

    // Load key molecules report
    try {
      const mdResponse = await ReaxAPI.getJobFile(jobId, 'key_molecules_flow.md');
      ReaxCharts.renderMarkdown('keyMoleculesReport', mdResponse.content);
    } catch (e) {
      console.warn('Failed to load key molecules report:', e);
      document.getElementById('keyMoleculesReport').innerHTML = '';
    }

    elements.chartsSection.classList.remove('hidden');

  } catch (error) {
    console.error('Error loading results:', error);
  }
}

/**
 * Handle download button
 */
function handleDownload(e) {
  if (!state.currentJobId) {
    e.preventDefault();
    showAlert('No results to download', 'error');
  }
}

/**
 * Scroll to charts section
 */
function scrollToCharts() {
  elements.chartsSection?.scrollIntoView({ behavior: 'smooth' });
}

/**
 * Load job history
 */
async function loadJobHistory() {
  try {
    const response = await ReaxAPI.getJobs(1, 20);
    state.jobs = response.jobs;
    renderJobHistory();
  } catch (error) {
    console.error('Failed to load job history:', error);
  }
}

/**
 * Render job history table
 */
function renderJobHistory() {
  if (!state.jobs.length) {
    elements.jobHistoryTable.innerHTML = `
      <tr>
        <td colspan="6" class="text-center py-8 text-gray-500">
          暂无任务记录
        </td>
      </tr>
    `;
    return;
  }

  elements.jobHistoryTable.innerHTML = state.jobs.map(job => {
    const statusClasses = {
      waiting: 'badge-queued',
      queued: 'badge-queued',
      active: 'badge-running',
      running: 'badge-running',
      completed: 'badge-completed',
      failed: 'badge-failed',
    };

    const statusTexts = {
      waiting: '等待中',
      queued: '等待中',
      active: '运行中',
      running: '运行中',
      completed: '已完成',
      failed: '失败',
    };

    const date = new Date(job.createdAt).toLocaleString('zh-CN');
    const canView = job.status === 'completed' || job.status === 'active' || job.status === 'waiting' || job.status === 'running' || job.status === 'queued';

    return `
      <tr>
        <td class="font-mono text-xs">${job.jobId.slice(0, 8)}...</td>
        <td class="max-w-xs truncate" title="${job.originalName}">${job.originalName}</td>
        <td><span class="badge ${statusClasses[job.status] || 'badge-queued'}">${statusTexts[job.status] || job.status}</span></td>
        <td>
          <div class="flex items-center gap-2">
            <div class="w-16 h-2 bg-gray-200 rounded-full overflow-hidden">
              <div class="h-full bg-primary-500" style="width: ${job.progress}%"></div>
            </div>
            <span class="text-xs">${job.progress}%</span>
          </div>
        </td>
        <td class="text-xs text-gray-500">${date}</td>
        <td>
          <div class="flex gap-1">
            ${canView ? `<button onclick="viewJobResults('${job.jobId}')" class="text-xs btn-secondary py-1 px-2">查看</button>` : ''}
            ${job.status === 'completed' ? `<button onclick="downloadJobResults('${job.jobId}')" class="text-xs btn-primary py-1 px-2">下载</button>` : ''}
            <button onclick="deleteJob('${job.jobId}')" class="text-xs btn-danger py-1 px-2">删除</button>
          </div>
        </td>
      </tr>
    `;
  }).join('');
}

/**
 * View job results
 */
async function viewJobResults(jobId) {
  state.currentJobId = jobId;
  const jobInfo = await ReaxAPI.getJobStatus(jobId);

  if (jobInfo.status === 'completed') {
    await loadResults(jobId);
    scrollToCharts();
  } else {
    // For active/waiting jobs, start monitoring and scroll to the real-time panel
    resetJobUI();
    startJobMonitoring(jobId);
    document.querySelector('.dashboard-card.lg\\:col-span-2')?.scrollIntoView({ behavior: 'smooth' });
  }
}

/**
 * Download job results
 */
function downloadJobResults(jobId) {
  ReaxAPI.downloadResults(jobId);
}

/**
 * Delete a job
 */
async function deleteJob(jobId) {
  if (!confirm('确定要删除这个任务吗?')) return;

  try {
    await ReaxAPI.deleteJob(jobId);
    showAlert('任务已删除', 'success');
    loadJobHistory();
  } catch (error) {
    showAlert(`删除失败: ${error.message}`, 'error');
  }
}

/**
 * Start queue status polling
 */
function startQueuePolling() {
  const updateQueueStatus = async () => {
    try {
      const response = await ReaxAPI.getJobs(1, 100);
      const waitingCount = response.jobs.filter(j => j.status === 'queued' || j.status === 'waiting').length;
      const runningCount = response.jobs.filter(j => j.status === 'running' || j.status === 'active').length;
      elements.queueCount.textContent = waitingCount;

      // Update queue dot color based on load
      const queueDot = document.getElementById('queueDot');
      const queuePulse = document.getElementById('queuePulse');
      if (queueDot && queuePulse) {
        if (runningCount > 0) {
          queueDot.className = 'relative inline-flex rounded-full h-3 w-3 bg-yellow-500';
          queuePulse.className = 'animate-ping absolute inline-flex h-full w-full rounded-full bg-yellow-400 opacity-75';
        } else {
          queueDot.className = 'relative inline-flex rounded-full h-3 w-3 bg-green-500';
          queuePulse.className = 'animate-ping absolute inline-flex h-full w-full rounded-full bg-green-400 opacity-75';
        }
      }
    } catch (error) {
      console.error('Failed to update queue status:', error);
    }
  };

  updateQueueStatus();
  setInterval(updateQueueStatus, 5000);
}

/**
 * Show alert message
 */
function showAlert(message, type = 'info') {
  const alertColors = {
    success: 'bg-green-100 border-green-400 text-green-700 dark:bg-green-900 dark:border-green-700 dark:text-green-100',
    error: 'bg-red-100 border-red-400 text-red-700 dark:bg-red-900 dark:border-red-700 dark:text-red-100',
    warning: 'bg-yellow-100 border-yellow-400 text-yellow-700 dark:bg-yellow-900 dark:border-yellow-700 dark:text-yellow-100',
    info: 'bg-blue-100 border-blue-400 text-blue-700 dark:bg-blue-900 dark:border-blue-700 dark:text-blue-100',
  };

  const alert = document.createElement('div');
  alert.className = `px-4 py-3 rounded-lg border ${alertColors[type]} animate-fade-in`;
  alert.innerHTML = `<p class="text-sm">${message}</p>`;

  elements.alertContainer.appendChild(alert);

  setTimeout(() => {
    alert.remove();
  }, 5000);
}

/**
 * Load sample file
 */
async function loadSampleFile() {
  try {
    const file = await ReaxAPI.loadSampleFile();
    selectFile(file);
    showAlert('示例文件已加载', 'success');
  } catch (error) {
    showAlert(`加载示例失败: ${error.message}`, 'error');
  }
}

/**
 * Show help modal
 */
function showHelp() {
  document.getElementById('helpModal').classList.remove('hidden');
}

/**
 * Close help modal
 */
function closeHelp() {
  document.getElementById('helpModal').classList.add('hidden');
}

// Global function exports for onclick handlers
window.viewJobResults = viewJobResults;
window.downloadJobResults = downloadJobResults;
window.deleteJob = deleteJob;
window.loadSampleFile = loadSampleFile;
window.showHelp = showHelp;
window.closeHelp = closeHelp;
