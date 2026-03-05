// Plot.js - Extensible plotting module for ReaxTools Web
// Uses Plotly.js for interactive charts

const filename_title = {
	"species_count.csv": "Species count versus time",
	"bond_count.csv": "Bond type count versus time",
	"atom_bonded_num_count.csv": "Atom type count versus time",
	"ring_count.csv": "Ring size count versus time",
	"reactions.dot": "Reaction Network (Main component)",
	"key_molecules_reactions.csv": "Key molecules analysis",
};

class ReaxToolsPlotter {
	constructor() {
		this.charts = new Map(); // Store chart instances
		this.chartContainer = null;
		this.chartData = new Map(); // Store chart data for persistence
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

	// Save chart data to localStorage
	saveChartData() {
		try {
			const chartDataArray = Array.from(this.chartData.entries());
			localStorage.setItem(
				"reaxTools_chart_data",
				JSON.stringify(chartDataArray)
			);
			console.log("Chart data saved to localStorage");
		} catch (error) {
			console.error("Failed to save chart data:", error);
		}
	}

	// Load chart data from localStorage
	loadChartData() {
		try {
			const saved = localStorage.getItem("reaxTools_chart_data");
			if (saved) {
				const chartDataArray = JSON.parse(saved);
				this.chartData = new Map(chartDataArray);
				console.log("Chart data loaded from localStorage");
				return true;
			}
		} catch (error) {
			console.error("Failed to load chart data:", error);
		}
		return false;
	}

	// Clear chart data from localStorage
	clearChartData() {
		try {
			localStorage.removeItem("reaxTools_chart_data");
			this.chartData.clear();
			console.log("Chart data cleared from localStorage");
		} catch (error) {
			console.error("Failed to clear chart data:", error);
		}
	}

	// Parse CSV data from string
	parseCSV(csvText) {
		const lines = csvText.trim().split("\n");
		if (lines.length < 2) {
			throw new Error("CSV data must have at least header and one data row");
		}

		const headers = lines[0].split(",");
		const data = [];

		for (let i = 1; i < lines.length; i++) {
			const values = lines[i].split(",");
			if (values.length === headers.length) {
				data.push(values.map((v) => parseFloat(v) || 0));
			}
		}

		return { headers, data };
	}

	// Create a line chart from CSV data
	createLineChart(csvText, chartId, title = "Data Visualization") {
		try {
			const { headers, data } = this.parseCSV(csvText);
			const defaultTitleFont = { family: "Consolas", size: 20 };
			const xlabel = "Frame";
			const ylabel = "Value";

			// Prepare traces for each column
			const traces = headers.map((header, index) => ({
				x: Array.from({ length: data.length }, (_, i) => i + 1), // Time steps
				y: data.map((row) => row[index]),
				type: "scatter",
				mode: "lines+markers",
				name: header,
				line: { width: 4 },
				marker: { size: 8 },
			}));

			const layout = {
				// title: {text: title, font: defaultTitleFont},
				xaxis: {
					title: { text: xlabel, font: defaultTitleFont },
					tickfont: defaultTitleFont,
				},
				yaxis: {
					title: { text: ylabel, font: defaultTitleFont },
					tickfont: defaultTitleFont,
				},
				hovermode: "closest",
				showlegend: true,
				legend: { x: -0.15, y: 1, font: { size: 18 } },
				margin: { l: 50, r: 50, t: 50, b: 50 },
				hoverlabel: { font: { family: "Consolas", size: 24 } },
			};

			const config = {
				responsive: true,
				displayModeBar: true,
				modeBarButtonsToRemove: ["pan2d", "lasso2d", "select2d"],
				displaylogo: false,
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
			console.error("Error creating line chart:", error);
			return false;
		}
	}

	// Create a bar chart from CSV data
	createBarChart(csvText, chartId, title = "Data Visualization") {
		try {
			const { headers, data } = this.parseCSV(csvText);
			const defaultTitleFont = { family: "Consolas", size: 20 };
			const xlabel = "Frame";
			const ylabel = "Value";

			// Use the last row of data for bar chart
			const lastRow = data[data.length - 1] || [];

			const trace = {
				x: headers,
				y: lastRow,
				type: "bar",
				marker: {
					color: "rgba(58,200,225,0.6)",
					line: {
						color: "rgba(58,200,225,1.0)",
						width: 1,
					},
				},
			};

			const layout = {
				title: title,
				xaxis: { title: "Categories" },
				yaxis: { title: "Count" },
				margin: { l: 50, r: 50, t: 50, b: 50 },
			};

			const config = {
				responsive: true,
				displayModeBar: true,
				modeBarButtonsToRemove: ["pan2d", "lasso2d", "select2d"],
				displaylogo: false,
			};

			if (this.charts.has(chartId)) {
				Plotly.react(chartId, [trace], layout, config);
			} else {
				Plotly.newPlot(chartId, [trace], layout, config);
				this.charts.set(chartId, true);
			}

			return true;
		} catch (error) {
			console.error("Error creating bar chart:", error);
			return false;
		}
	}

	// Auto-detect chart type and create appropriate visualization
	autoCreateChart(csvText, chartId, title = "Data Visualization") {
		const { headers, data } = this.parseCSV(csvText);

		// Special handling for key_molecules_reactions.csv
		if (title === "key_molecules_reactions.csv") {
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
	createKeyMoleculesVisualization(
		csvText,
		chartId,
		title = "Key Molecules Reactions"
	) {
		try {
			const container = document.getElementById(chartId);
			if (!container) {
				console.error(`Container with id '${chartId}' not found`);
				return false;
			}

			// Clear container
			container.innerHTML = "";

			// Parse CSV data
			const lines = csvText.trim().split("\n");
			if (lines.length < 2) {
				throw new Error("CSV data must have at least header and one data row");
			}

			const headers = lines[0].split(",");
			const data = [];

			for (let i = 1; i < lines.length; i++) {
				const values = lines[i].split(",");
				if (values.length === headers.length) {
					const row = {};
					headers.forEach((header, index) => {
						row[header.trim()] = values[index] ? values[index].trim() : "";
					});
					data.push(row);
				}
			}

			// Create container for charts and table
			const chartsContainer = document.createElement("div");
			chartsContainer.className = "key-mol-charts-container";

			// Create bar chart for atom transfers
			const chartDiv = document.createElement("div");
			chartDiv.id = `${chartId}-bar`;
			chartDiv.className = "key-mol-chart-div";
			chartsContainer.appendChild(chartDiv);

			// Add charts container to main container first
			container.appendChild(chartsContainer);

			// Prepare data for bar chart
			const molecules = data.map((row) => row.molecule || "Unknown");
			const inTransfers = data.map(
				(row) => parseFloat(row["in atom transfer"]) || 0
			);
			const outTransfers = data.map(
				(row) => parseFloat(row["out atom transfer"]) || 0
			);

			const defaultTitleFont = { family: "Consolas", size: 20 };
			const in_color = "#2980b9";
			const out_color = "#e74c3c";

			const traces = [
				{
					x: molecules,
					y: inTransfers,
					type: "bar",
					name: "In Atom Transfer",
					marker: { color: in_color },
				},
				{
					x: molecules,
					y: outTransfers,
					type: "bar",
					name: "Out Atom Transfer",
					marker: { color: out_color },
				},
			];

			const layout = {
				xaxis: {
					title: { text: "Molecules", font: defaultTitleFont },
					tickfont: defaultTitleFont,
					tickangle: -45,
				},
				yaxis: {
					title: { text: "Flux", font: defaultTitleFont },
					tickfont: defaultTitleFont
				},
				barmode: "group",
				margin: { l: 100, r: 50, t: 80, b: 100 },
				showlegend: true,
				legend: { x: 0.85, y: 1, font: { size: 18 } },
				hoverlabel: { font: { family: "Consolas", size: 24 } },
			};

			const config = {
				responsive: true,
				displayModeBar: true,
				modeBarButtonsToRemove: ["pan2d", "lasso2d", "select2d"],
				displaylogo: false,
			};

			// Now create the plot after DOM element is added
			Plotly.newPlot(`${chartId}-bar`, traces, layout, config);
			this.charts.set(`${chartId}-bar`, true);

			// Create HTML table for reaction sources and destinations
			const tableContainer = document.createElement("div");
			tableContainer.className = "key-mol-table-container";

			const tableTitle = document.createElement("h4");
			tableTitle.textContent =
				"Key Molecules: Reaction Sources and Destinations";
			tableTitle.className = "key-mol-table-title";
			tableContainer.appendChild(tableTitle);

			const table = document.createElement("table");
			table.className = "key-mol-table";

			// Create table header
			const thead = document.createElement("thead");
			const headerRow = document.createElement("tr");
			const headers2 = [
				"Molecule",
				"From 1",
				"From 2",
				"From 3",
				"From 4",
				"From 5",
				"To 1",
				"To 2",
				"To 3",
				"To 4",
				"To 5",
			];
			headers2.forEach((header) => {
				const th = document.createElement("th");
				th.textContent = header;
				headerRow.appendChild(th);
			});
			thead.appendChild(headerRow);
			table.appendChild(thead);

			// Create table body
			const tbody = document.createElement("tbody");
			data.forEach((row, index) => {
				const tr = document.createElement("tr");
				// Add molecule name
				const tdMolecule = document.createElement("td");
				tdMolecule.textContent = row.molecule || "";
				tdMolecule.className = "molecule";
				tr.appendChild(tdMolecule);
				// Add from columns
				for (let i = 1; i <= 5; i++) {
					const td = document.createElement("td");
					const value = row[`from ${i}`] || "";
					td.textContent = value;
					td.className = "from";
					tr.appendChild(td);
				}
				// Add to columns
				for (let i = 1; i <= 5; i++) {
					const td = document.createElement("td");
					const value = row[`to ${i}`] || "";
					td.textContent = value;
					td.className = "to";
					tr.appendChild(td);
				}
				tbody.appendChild(tr);
			});
			table.appendChild(tbody);
			tableContainer.appendChild(table);

			// Add scrollable container for table
			const scrollContainer = document.createElement("div");
			scrollContainer.className = "key-mol-scroll-container";
			scrollContainer.appendChild(table);
			tableContainer.appendChild(scrollContainer);

			chartsContainer.appendChild(tableContainer);

			return true;
		} catch (error) {
			console.error("Error creating key molecules visualization:", error);
			return false;
		}
	}

	// Create multiple charts from different CSV files and DOT files
	createMultipleCharts(files, containerId) {
		if (!this.chartContainer) {
			this.init(containerId);
		}

		const chartContainer = this.chartContainer;
		chartContainer.innerHTML = ""; // Clear existing charts

		// 在 createMultipleCharts(files, containerId) 方法开头添加如下排序逻辑
		const orderedFiles = [];
		const used = new Set();
		Object.keys(filename_title).forEach((name) => {
			const file = files.find(
				(f) => f.name.toLowerCase() === name.toLowerCase()
			);
			if (file) {
				orderedFiles.push(file);
				used.add(file);
			}
		});
		files.forEach((f) => {
			if (!used.has(f)) orderedFiles.push(f);
		});
		// 后续遍历 orderedFiles 替代 files
		orderedFiles.forEach((file, index) => {
			// 为每个文件创建一个主 chartDiv
			const chartDiv = document.createElement("div");
			chartDiv.id = `chart-${index}`;
			chartDiv.className = "chart-main-div";
			chartContainer.appendChild(chartDiv);

			// Title
			const titleDiv = document.createElement("div");
			titleDiv.className = "chart-suptitle";
			titleDiv.textContent = filename_title[file.name.toLowerCase()];
			chartDiv.appendChild(titleDiv);

			// Network (DOT) 文件
			if (file.name.toLowerCase().endsWith("reactions.dot")) {
				// 独立 network 容器
				const networkDiv = document.createElement("div");
				networkDiv.id = `chart-${index}-network`;
				networkDiv.className = "chart-network-div";
				chartDiv.appendChild(networkDiv);
				if (
					window.reaxToolsPlotNetwork &&
					window.reaxToolsPlotNetwork.createSigmaGraph
				) {
					if (window.graphlibDot) {
						try {
							const graph = window.graphlibDot.read(file.content);
							const nodes = graph.nodes().map((id) => {
								const attr = graph.node(id) || {};
								return Object.assign({ id }, attr);
							});
							const edges = graph.edges().map((e) => {
								const attr = graph.edge(e) || {};
								return Object.assign(
									{ id: e.name || `${e.v}->${e.w}`, source: e.v, target: e.w },
									attr
								);
							});
							const graphData = { nodes, edges };
							window.reaxToolsPlotNetwork.createSigmaGraph(
								graphData,
								networkDiv.id,
								file.name,
								this.charts
							);

							// Save chart data
							this.chartData.set(`chart-${index}-network`, {
								type: "network",
								name: file.name,
								data: graphData,
							});
						} catch (err) {
							console.error("DOT parsing failed:", err);
							networkDiv.innerHTML =
								'<p class="error-message">Error: DOT parse failed, can not create network chart.</p>';
						}
					} else {
						networkDiv.innerHTML =
							'<p class="error-message">Error: graphlib-dot load failed，can not parse DOT file.</p>';
					}
				} else {
					console.error(
						"plot_network.js not loaded or createSigmaGraph not found"
					);
				}
			} else if (file.name === "key_molecules_reactions.csv") {
				// 独立 key molecules 容器
				const keyMolDiv = document.createElement("div");
				keyMolDiv.id = `chart-${index}-keymol`;
				keyMolDiv.className = "chart-keymol-div";
				chartDiv.appendChild(keyMolDiv);
				this.createKeyMoleculesVisualization(
					file.content,
					keyMolDiv.id,
					file.name
				);

				// Save chart data
				this.chartData.set(`chart-${index}-keymol`, {
					type: "key_molecules",
					name: file.name,
					content: file.content,
				});
			} else {
				// 普通 CSV 图表独立容器
				const csvDiv = document.createElement("div");
				csvDiv.id = `chart-${index}-csv`;
				csvDiv.className = "chart-csv-div";
				chartDiv.appendChild(csvDiv);
				this.autoCreateChart(file.content, csvDiv.id, file.name);

				// Save chart data
				this.chartData.set(`chart-${index}-csv`, {
					type: "csv",
					name: file.name,
					content: file.content,
				});
			}
		});

		// Save all chart data
		this.saveChartData();
	}

	// Restore charts from saved data
	restoreCharts(containerId) {
		if (!this.loadChartData()) {
			return false;
		}

		if (!this.chartContainer) {
			this.init(containerId);
		}

		const chartContainer = this.chartContainer;
		chartContainer.innerHTML = ""; // Clear existing charts

		let index = 0;
		this.chartData.forEach((data, chartId) => {
			// Create main chart div
			const chartDiv = document.createElement("div");
			chartDiv.id = `chart-${index}`;
			chartDiv.className = "chart-main-div";
			chartContainer.appendChild(chartDiv);

			const titleDiv = document.createElement("div");
			titleDiv.className = "chart-suptitle";
			titleDiv.textContent = filename_title[data.name.toLowerCase()];
			chartDiv.appendChild(titleDiv);

			if (data.type === "network") {
				// Restore network chart
				const networkDiv = document.createElement("div");
				networkDiv.id = chartId;
				networkDiv.className = "chart-network-div";
				chartDiv.appendChild(networkDiv);

				if (
					window.reaxToolsPlotNetwork &&
					window.reaxToolsPlotNetwork.createSigmaGraph
				) {
					window.reaxToolsPlotNetwork.createSigmaGraph(
						data.data,
						networkDiv.id,
						data.name,
						this.charts
					);
				}
			} else if (data.type === "key_molecules") {
				// Restore key molecules chart
				const keyMolDiv = document.createElement("div");
				keyMolDiv.id = chartId;
				keyMolDiv.className = "chart-keymol-div";
				chartDiv.appendChild(keyMolDiv);

				this.createKeyMoleculesVisualization(
					data.content,
					keyMolDiv.id,
					data.name
				);
			} else if (data.type === "csv") {
				// Restore CSV chart
				const csvDiv = document.createElement("div");
				csvDiv.id = chartId;
				csvDiv.className = "chart-csv-div";
				chartDiv.appendChild(csvDiv);
				this.autoCreateChart(data.content, csvDiv.id, data.name);
			}

			index++;
		});

		return true;
	}

	// Clear all charts
	clearCharts() {
		this.charts.forEach((_, chartId) => {
			Plotly.purge(chartId);
		});
		this.charts.clear();
		this.chartData.clear();
		this.clearChartData(); // Clear from localStorage
		if (this.chartContainer) {
			this.chartContainer.innerHTML = "";
		}
	}

	// Export chart as PNG
	exportChart(chartId, filename = "chart.png") {
		if (this.charts.has(chartId)) {
			Plotly.downloadImage(chartId, {
				format: "png",
				filename: filename,
				height: 600,
				width: 800,
			});
		}
	}
}

// Global plotter instance
window.reaxToolsPlotter = new ReaxToolsPlotter();
