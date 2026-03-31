/**
 * Simplified chart/image module for ReaxTools Web
 * Now only handles markdown rendering and clearing the image grid.
 */

/**
 * Render markdown content
 */
function renderMarkdown(elementId, markdownText) {
  const element = document.getElementById(elementId);
  if (!element) return;

  // Simple markdown to HTML conversion
  let html = markdownText
    // Headers
    .replace(/^### (.*$)/gim, '<h3 class="text-lg font-semibold mt-4 mb-2">$1</h3>')
    .replace(/^## (.*$)/gim, '<h2 class="text-xl font-semibold mt-6 mb-3">$1</h2>')
    .replace(/^# (.*$)/gim, '<h1 class="text-2xl font-bold mt-6 mb-4">$1</h1>')
    // Bold and italic
    .replace(/\*\*\*(.*?)\*\*\*/g, '<strong><em>$1</em></strong>')
    .replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>')
    .replace(/\*(.*?)\*/g, '<em>$1</em>')
    // Code
    .replace(/`([^`]+)`/g, '<code class="bg-gray-100 dark:bg-gray-800 px-1 rounded text-sm">$1</code>')
    // Lists
    .replace(/^\s*[-*+]\s+(.*$)/gim, '<li class="ml-4">$1</li>')
    // Links
    .replace(/\[([^\]]+)\]\(([^)]+)\)/g, '<a href="$2" class="text-primary-500 hover:underline" target="_blank">$1</a>')
    // Paragraphs
    .replace(/\n\n/g, '</p><p class="mb-2">')
    // Line breaks
    .replace(/\n/g, '<br>');

  element.innerHTML = `<div class="prose dark:prose-invert max-w-none">${html}</div>`;
}

/**
 * Clear all generated content
 */
function clearCharts() {
  const imagesGrid = document.getElementById('imagesGrid');
  if (imagesGrid) {
    imagesGrid.innerHTML = '';
  }

  const report = document.getElementById('keyMoleculesReport');
  if (report) {
    report.innerHTML = '';
  }
}

// Export chart functions
window.ReaxCharts = {
  renderMarkdown,
  clearCharts,
};
