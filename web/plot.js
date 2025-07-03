// Plot.js - Extensible plotting module for ReaxTools Web
// Uses Plotly.js for interactive charts

class ReaxToolsPlotter {
    constructor() {
        this.charts = new Map(); // Store chart instances
        this.chartContainer = null;
    }

    // Initialize the plotter with a container
    init(containerId) {
        this.chartContainer = document.getElementById(containerId);
        if (!this.chartContainer) {
            console.error(`Container with id '${containerId}' not found`);
            return false;
        }
        
        return true;
    }

    // Parse CSV data from string
    parseCSV(csvText) {
        const lines = csvText.trim().split('\n');
        if (lines.length < 2) {
            throw new Error('CSV data must have at least header and one data row');
        }

        const headers = lines[0].split(',');
        const data = [];
        
        for (let i = 1; i < lines.length; i++) {
            const values = lines[i].split(',');
            if (values.length === headers.length) {
                data.push(values.map(v => parseFloat(v) || 0));
            }
        }

        return { headers, data };
    }

    // Create a line chart from CSV data
    createLineChart(csvText, chartId, title = 'Data Visualization') {
        try {
            const { headers, data } = this.parseCSV(csvText);
            
            // Prepare traces for each column
            const traces = headers.map((header, index) => ({
                x: Array.from({ length: data.length }, (_, i) => i + 1), // Time steps
                y: data.map(row => row[index]),
                type: 'scatter',
                mode: 'lines+markers',
                name: header,
                line: { width: 2 },
                marker: { size: 4 }
            }));

            const layout = {
                title: title,
                xaxis: { title: 'Time Step' },
                yaxis: { title: 'Count' },
                hovermode: 'closest',
                showlegend: true,
                legend: { x: 0, y: 1 },
                margin: { l: 50, r: 50, t: 50, b: 50 }
            };

            const config = {
                responsive: true,
                displayModeBar: true,
                modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d'],
                displaylogo: false
            };

            // Create or update chart
            if (this.charts.has(chartId)) {
                Plotly.react(chartId, traces, layout, config);
            } else {
                Plotly.newPlot(chartId, traces, layout, config);
                this.charts.set(chartId, true);
            }

            return true;
        } catch (error) {
            console.error('Error creating line chart:', error);
            return false;
        }
    }

    // Create a bar chart from CSV data
    createBarChart(csvText, chartId, title = 'Data Visualization') {
        try {
            const { headers, data } = this.parseCSV(csvText);
            
            // Use the last row of data for bar chart
            const lastRow = data[data.length - 1] || [];
            
            const trace = {
                x: headers,
                y: lastRow,
                type: 'bar',
                marker: {
                    color: 'rgba(58,200,225,0.6)',
                    line: {
                        color: 'rgba(58,200,225,1.0)',
                        width: 1
                    }
                }
            };

            const layout = {
                title: title,
                xaxis: { title: 'Categories' },
                yaxis: { title: 'Count' },
                margin: { l: 50, r: 50, t: 50, b: 50 }
            };

            const config = {
                responsive: true,
                displayModeBar: true,
                modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d'],
                displaylogo: false
            };

            if (this.charts.has(chartId)) {
                Plotly.react(chartId, [trace], layout, config);
            } else {
                Plotly.newPlot(chartId, [trace], layout, config);
                this.charts.set(chartId, true);
            }

            return true;
        } catch (error) {
            console.error('Error creating bar chart:', error);
            return false;
        }
    }

    // Auto-detect chart type and create appropriate visualization
    autoCreateChart(csvText, chartId, title = 'Data Visualization') {
        const { headers, data } = this.parseCSV(csvText);
        
        // Special handling for key_molecules_reactions.csv
        if (title === 'key_molecules_reactions.csv') {
            return this.createKeyMoleculesVisualization(csvText, chartId, title);
        }
        
        // If we have time series data (multiple rows), create line chart
        if (data.length > 1) {
            return this.createLineChart(csvText, chartId, title);
        } else {
            // If single row, create bar chart
            return this.createBarChart(csvText, chartId, title);
        }
    }

    // Create specialized visualization for key_molecules_reactions.csv
    createKeyMoleculesVisualization(csvText, chartId, title = 'Key Molecules Reactions') {
        try {
            const container = document.getElementById(chartId);
            if (!container) {
                console.error(`Container with id '${chartId}' not found`);
                return false;
            }

            // Clear container
            container.innerHTML = '';

            // Parse CSV data
            const lines = csvText.trim().split('\n');
            if (lines.length < 2) {
                throw new Error('CSV data must have at least header and one data row');
            }

            const headers = lines[0].split(',');
            const data = [];
            
            for (let i = 1; i < lines.length; i++) {
                const values = lines[i].split(',');
                if (values.length === headers.length) {
                    const row = {};
                    headers.forEach((header, index) => {
                        row[header.trim()] = values[index] ? values[index].trim() : '';
                    });
                    data.push(row);
                }
            }

            // Create container for charts and table
            const chartsContainer = document.createElement('div');
            chartsContainer.style.cssText = `
                display: flex;
                flex-direction: column;
                gap: 20px;
                width: 100%;
            `;

            // Create bar chart for atom transfers
            const chartDiv = document.createElement('div');
            chartDiv.id = `${chartId}-bar`;
            chartDiv.style.width = '100%';
            chartDiv.style.height = '400px';
            chartsContainer.appendChild(chartDiv);

            // Add charts container to main container first
            container.appendChild(chartsContainer);

            // Prepare data for bar chart
            const molecules = data.map(row => row.molecule || 'Unknown');
            const inTransfers = data.map(row => parseFloat(row['in atom transfer']) || 0);
            const outTransfers = data.map(row => parseFloat(row['out atom transfer']) || 0);

            const traces = [
                {
                    x: molecules,
                    y: inTransfers,
                    type: 'bar',
                    name: 'In Atom Transfer',
                    marker: {
                        color: 'rgba(58,200,225,0.6)',
                        line: {
                            color: 'rgba(58,200,225,1.0)',
                            width: 1
                        }
                    }
                },
                {
                    x: molecules,
                    y: outTransfers,
                    type: 'bar',
                    name: 'Out Atom Transfer',
                    marker: {
                        color: 'rgba(255,99,132,0.6)',
                        line: {
                            color: 'rgba(255,99,132,1.0)',
                            width: 1
                        }
                    }
                }
            ];

            const layout = {
                title: 'Key Molecules - Atom Transfer Analysis',
                xaxis: { 
                    title: 'Molecules',
                    tickangle: -45
                },
                yaxis: { title: 'Atom Transfer Count' },
                barmode: 'group',
                margin: { l: 50, r: 50, t: 80, b: 100 },
                showlegend: true,
                legend: { x: 0, y: 1 }
            };

            const config = {
                responsive: true,
                displayModeBar: true,
                modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d'],
                displaylogo: false
            };

            // Now create the plot after DOM element is added
            Plotly.newPlot(`${chartId}-bar`, traces, layout, config);
            this.charts.set(`${chartId}-bar`, true);

            // Create HTML table for reaction sources and destinations
            const tableContainer = document.createElement('div');
            tableContainer.style.cssText = `
                margin-top: 20px;
                padding: 15px;
                background-color: #f8f9fa;
                border-radius: 5px;
                border: 1px solid #dee2e6;
            `;

            const tableTitle = document.createElement('h4');
            tableTitle.textContent = 'Key Molecules - Reaction Sources and Destinations';
            tableTitle.style.cssText = `
                margin: 0 0 15px 0;
                color: #34495e;
                font-size: 1.1rem;
                font-weight: bold;
            `;
            tableContainer.appendChild(tableTitle);

            const table = document.createElement('table');
            table.className = 'key-mol-table';

            // Create table header
            const thead = document.createElement('thead');
            const headerRow = document.createElement('tr');
            const headers2 = ['Molecule', 'From 1', 'From 2', 'From 3', 'From 4', 'From 5', 'To 1', 'To 2', 'To 3', 'To 4', 'To 5'];
            headers2.forEach(header => {
                const th = document.createElement('th');
                th.textContent = header;
                headerRow.appendChild(th);
            });
            thead.appendChild(headerRow);
            table.appendChild(thead);

            // Create table body
            const tbody = document.createElement('tbody');
            data.forEach((row, index) => {
                const tr = document.createElement('tr');
                // Add molecule name
                const tdMolecule = document.createElement('td');
                tdMolecule.textContent = row.molecule || '';
                tdMolecule.className = 'molecule';
                tr.appendChild(tdMolecule);
                // Add from columns
                for (let i = 1; i <= 5; i++) {
                    const td = document.createElement('td');
                    const value = row[`from ${i}`] || '';
                    td.textContent = value;
                    td.className = 'from';
                    tr.appendChild(td);
                }
                // Add to columns
                for (let i = 1; i <= 5; i++) {
                    const td = document.createElement('td');
                    const value = row[`to ${i}`] || '';
                    td.textContent = value;
                    td.className = 'to';
                    tr.appendChild(td);
                }
                tbody.appendChild(tr);
            });
            table.appendChild(tbody);
            tableContainer.appendChild(table);

            // Add scrollable container for table
            const scrollContainer = document.createElement('div');
            scrollContainer.style.cssText = `
                max-height: 400px;
                overflow-y: auto;
                border: 1px solid #dee2e6;
                border-radius: 5px;
            `;
            scrollContainer.appendChild(table);
            tableContainer.appendChild(scrollContainer);

            chartsContainer.appendChild(tableContainer);

            return true;
        } catch (error) {
            console.error('Error creating key molecules visualization:', error);
            return false;
        }
    }

    // Create multiple charts from different CSV files and DOT files
    createMultipleCharts(files, containerId) {
        if (!this.chartContainer) {
            this.init(containerId);
        }

        const chartContainer = this.chartContainer;
        chartContainer.innerHTML = ''; // Clear existing charts

        files.forEach((file, index) => {
            // 为每个文件创建一个主 chartDiv
            const chartDiv = document.createElement('div');
            chartDiv.id = `chart-${index}`;
            chartDiv.style.width = '100%';
            chartDiv.style.minHeight = '400px';
            chartDiv.style.marginBottom = '20px';
            chartContainer.appendChild(chartDiv);

            // Network (DOT) 文件
            if (file.name.toLowerCase().endsWith('.dot')) {
                // 独立 network 容器
                const networkDiv = document.createElement('div');
                networkDiv.id = `chart-${index}-network`;
                networkDiv.style.width = '100%';
                networkDiv.style.minHeight = '400px';
                chartDiv.appendChild(networkDiv);
                if (window.reaxToolsPlotNetwork && window.reaxToolsPlotNetwork.createSigmaGraph) {
                    if (window.graphlibDot) {
                        try {
                            const graph = window.graphlibDot.read(file.content);
                            const nodes = graph.nodes().map(id => {
                                const attr = graph.node(id) || {};
                                return Object.assign({ id }, attr);
                            });
                            const edges = graph.edges().map(e => {
                                const attr = graph.edge(e) || {};
                                return Object.assign({ id: e.name || `${e.v}->${e.w}`, source: e.v, target: e.w }, attr);
                            });
                            const graphData = { nodes, edges };
                            window.reaxToolsPlotNetwork.createSigmaGraph(graphData, networkDiv.id, file.name, this.charts);
                        } catch (err) {
                            console.error('DOT parsing failed:', err);
                            networkDiv.innerHTML = '<p style="color:red">DOT文件解析失败，无法显示网络图。</p>';
                        }
                    } else {
                        networkDiv.innerHTML = '<p style="color:red">graphlib-dot 未加载，无法解析DOT文件。</p>';
                    }
                } else {
                    console.error('plot_network.js not loaded or createSigmaGraph not found');
                }
            } else if (file.name === 'key_molecules_reactions.csv') {
                // 独立 key molecules 容器
                const keyMolDiv = document.createElement('div');
                keyMolDiv.id = `chart-${index}-keymol`;
                keyMolDiv.style.width = '100%';
                keyMolDiv.style.minHeight = '400px';
                chartDiv.appendChild(keyMolDiv);
                this.createKeyMoleculesVisualization(file.content, keyMolDiv.id, file.name);
            } else {
                // 普通 CSV 图表独立容器
                const csvDiv = document.createElement('div');
                csvDiv.id = `chart-${index}-csv`;
                csvDiv.style.width = '100%';
                csvDiv.style.minHeight = '400px';
                chartDiv.appendChild(csvDiv);
                this.autoCreateChart(file.content, csvDiv.id, file.name);
            }
        });
    }

    // Clear all charts
    clearCharts() {
        this.charts.forEach((_, chartId) => {
            Plotly.purge(chartId);
        });
        this.charts.clear();
        if (this.chartContainer) {
            this.chartContainer.innerHTML = '';
        }
    }

    // Export chart as PNG
    exportChart(chartId, filename = 'chart.png') {
        if (this.charts.has(chartId)) {
            Plotly.downloadImage(chartId, {
                format: 'png',
                filename: filename,
                height: 600,
                width: 800
            });
        }
    }
}

// Global plotter instance
window.reaxToolsPlotter = new ReaxToolsPlotter(); 