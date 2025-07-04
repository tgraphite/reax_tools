// ===== Utility Functions =====
/**
 * Get the minimum and maximum value from a numeric array.
 * @param {number[]} seriesData - Array of numbers.
 * @returns {{minValue: number, maxValue: number}}
 */
function getMinMaxValues(seriesData) {
    let minValue = 0;
    let maxValue = 0;
    let value = 0;
    for (let i = 0; i < seriesData.length; i++) {
        value = seriesData[i];
        if (value < minValue) minValue = value;
        if (value > maxValue) maxValue = value;
    }
    if (maxValue == minValue)
        throw new Error("function getMinMaxValues: result min == max.");
    return { minValue, maxValue };
}

function linearScale(value, minValue, maxValue) {
    if (maxValue <= minValue)
        throw new Error("function linearScale: max <= min.");

    return (value - minValue) / (maxValue - minValue);
}

function log2Scale(value, minValue, maxValue) {
    if (minValue <= 0 || maxValue <= 0) {
        console.log("min %d max %d", minValue, maxValue);
        throw new Error("function log2Scale: min or max <= 0.");
    }

    let scaled =
        (Math.log2(value) - Math.log2(minValue)) /
        (Math.log2(maxValue) - Math.log2(minValue));
    return Math.max(0, Math.min(1, scaled));
}

/**
 * Convert hex color to RGB array.
 * @param {string} hex - Hex color string.
 * @returns {number[]}
 */
function hexToRgb(hex) {
    hex = hex.replace("#", "");
    if (hex.length === 3)
        hex = hex
            .split("")
            .map((x) => x + x)
            .join("");
    const num = parseInt(hex, 16);
    return [(num >> 16) & 255, (num >> 8) & 255, num & 255];
}

/**
 * Convert RGB array to hex color string.
 * @param {number[]} rgb - Array of [r, g, b].
 * @returns {string}
 */
function rgbToHex([r, g, b]) {
    return "#" + [r, g, b].map((x) => x.toString(16).padStart(2, "0")).join("");
}

/**
 * Interpolate color for node balance.
 * @param {number} balance -
 * @param {number} minbalance -
 * @param {number} maxbalance -
 * @param {string} lowColor - Hex color for low.
 * @param {string} midColor - Hex color for mid.
 * @param {string} highColor - Hex color for high.
 * @returns {string}
 */
function getBalanceColor(
    balance,
    minBalance,
    maxBalance,
    lowColor,
    midColor,
    highColor
) {
    const scale = (maxBalance - minBalance) / 8 || 10;
    let scaled = 0.5 + Math.atan(balance / scale) / Math.PI;

    let colorA, colorB, t;
    if (scaled <= 0.5) {
        colorA = hexToRgb(lowColor);
        colorB = hexToRgb(midColor);
        t = scaled / 0.5;
    } else {
        colorA = hexToRgb(midColor);
        colorB = hexToRgb(highColor);
        t = (scaled - 0.5) / 0.5;
    }
    const rgb = colorA.map((a, i) => Math.round(a * (1 - t) + colorB[i] * t));
    return rgbToHex(rgb);
}

/**
 * Interpolate color for edge weight.
 * @param {number} weight - Edge weight.
 * @param {number} minEdgeWeight - Minimum edge weight.
 * @param {number} maxEdgeWeight - Maximum edge weight.
 * @param {string} lowColor - Hex color for low.
 * @param {string} highColor - Hex color for high.
 * @returns {string}
 */
function getEdgeColor(
    weight,
    minEdgeWeight,
    maxEdgeWeight,
    lowColor,
    highColor
) {
    let scaled = log2Scale(weight, minEdgeWeight, maxEdgeWeight);
    const rgbA = hexToRgb(lowColor);
    const rgbB = hexToRgb(highColor);
    const rgb = rgbA.map((a, i) =>
        Math.round(a * (1 - scaled) + rgbB[i] * scaled)
    );
    return rgbToHex(rgb);
}

/**
 * Interpolate width for edge weight.
 * @param {number} weight - Edge weight.
 * @param {number} minEdgeWeight - Minimum edge weight.
 * @param {number} maxEdgeWeight - Maximum edge weight.
 * @param {number} lowWidth - Width for low.
 * @param {number} highWidth - Width for high.
 * @returns {number}
 */
function getEdgeWidth(
    weight,
    minEdgeWeight,
    maxEdgeWeight,
    lowWidth,
    highWidth
) {
    let scaled = log2Scale(weight, minEdgeWeight, maxEdgeWeight);
    const width = (highWidth - lowWidth) * scaled + lowWidth;
    return width;
}

/**
 * Interpolate size for node degree.
 * @param {number} weight - Node degree.
 * @param {number} minNodeDegree - Minimum node degree.
 * @param {number} maxNodeDegree - Maximum node degree.
 * @param {number} lowSize - Size for low.
 * @param {number} highSize - Size for high.
 * @returns {number}
 */
function getNodeSize(degree, minNodeDegree, maxNodeDegree, lowSize, highSize) {
    let scaled = log2Scale(degree, minNodeDegree, maxNodeDegree);
    const size = (highSize - lowSize) * scaled + lowSize;
    return size;
}

// ===== Color configuration for easy theme adjustment =====
const COLORS = {
    title: "#34495e",
    border: "#3498db",
    background: "#fff",
    nodeDefault: "#eee",
    sigmaFail: "#666",
    balanceLow: "#FF4F0F", // reactant
    balanceMid: "#06923E", // intermediate
    balanceHigh: "#3674B5", // product
    edgeLow: "#EFEFE0", // min edge weight color
    edgeHigh: "#819A91", // max edge weight color
    nodeSelected: "#ff9800", // selected node color
    edgeSelected: "#e53935", // selected edge color
    edgeDefault: "#E6E6E6",
};

// ===== Size configuration for easy adjustment =====
const SIZES = {
    containerHeight: 1000, // px
    borderRadius: 5, // px
    borderWidth: 2, // px
    edgeWidthScale: 1.5, // edge width scale factor
    edgeMinWidth: 4, // px, min edge width
    edgeMaxWidth: 12, // px, max edge width
    edgeLabelSize: 18, // px
    nodeMinSize: 12, // px, min node size
    nodeMaxSize: 30, // px, max node size
    nodeLabelFontSize: 20, // px, node label font size
};

// 兼容 sigma.js UMD 导出
if (
    typeof window.sigma === "undefined" &&
    typeof window.Sigma !== "undefined"
) {
    window.sigma = window.Sigma;
}

function createSigmaGraph(
    graphData,
    chartId,
    title = "Reaction Network",
    chartsMap
) {
    // graphData should be { nodes: [...], edges: [...] }
    const container = document.getElementById(chartId);
    if (!container) {
        console.error(`Container with id '${chartId}' not found`);
        return false;
    }
    container.innerHTML = "";
    // Sigma container
    const sigmaContainer = document.createElement("div");
    sigmaContainer.className = "sigma-container";
    sigmaContainer.style.height = `${SIZES.containerHeight}px`;
    sigmaContainer.id = chartId + "_sigma";
    container.appendChild(sigmaContainer);
    // Info box below sigmaContainer
    let infoBox = container.querySelector(`#${chartId}-info`);
    if (!infoBox) {
        infoBox = document.createElement("div");
        infoBox.id = `${chartId}-info`;
        infoBox.className = "network-info-box";
        // Two columns for info
        const leftCol = document.createElement("div");
        leftCol.className = "network-info-col left";
        infoBox.appendChild(leftCol);
        const rightCol = document.createElement("div");
        rightCol.className = "network-info-col right";
        infoBox.appendChild(rightCol);
        container.appendChild(infoBox);
    }
    // Info note below infoBox
    let infoNote = container.querySelector(`#${chartId}-info-note`);
    if (!infoNote) {
        infoNote = document.createElement("p");
        infoNote.id = `${chartId}-info-note`;
        infoNote.className = "network-info-note";
        infoNote.textContent =
            "Node color: relative inflow/outflow ratio. More outflow = redder, more inflow = bluer. Node size: total flow. Edge color: reaction frequency (darker = more frequent). Click a node to highlight related reactions and show details.";
        container.appendChild(infoNote);
    }
    // Sigma.js rendering
    if (window.sigma && window.graphology) {
        try {
            const Graph = window.graphology;
            const Sigma = window.sigma;
            const g = new Graph();

            let minEdgeWeight = Infinity,
                maxEdgeWeight = -Infinity;
            let minDegree = Infinity,
                maxDegree = -Infinity;
            let minBalance = Infinity,
                maxBalance = -Infinity;
            const inWeights = {},
                outWeights = {};

            graphData.nodes.forEach((node) => {
                inWeights[node.id] = 0;
                outWeights[node.id] = 0;
            });

            graphData.edges.forEach((edge) => {
                let weight = 1;
                if (edge.label) {
                    const match = String(edge.label).match(/R=(\d+(\.\d+)?)/);
                    if (match) weight = parseFloat(match[1]);
                }
                weight = Number(weight);
                outWeights[edge.source] += weight;
                inWeights[edge.target] += weight;
                if (weight < minEdgeWeight) minEdgeWeight = weight;
                if (weight > maxEdgeWeight) maxEdgeWeight = weight;
            });

            const nodeBalances = {},
                nodeDegrees = {};
            graphData.nodes.forEach((node) => {
                const inW = inWeights[node.id];
                const outW = outWeights[node.id];
                const total = inW + outW;
                let balance = 0.5;
                if (total > 0) balance = inW - outW;
                nodeBalances[node.id] = balance;
                if (balance < minBalance) minBalance = balance;
                if (balance > maxBalance) maxBalance = balance;
                const degree = inW + outW;
                nodeDegrees[node.id] = degree;
                if (degree < minDegree) minDegree = degree;
                if (degree > maxDegree) maxDegree = degree;
            });

            // Sort nodes by balance
            const sortedNodes = [...graphData.nodes].sort(
                (a, b) => nodeBalances[a.id] - nodeBalances[b.id]
            );
            const N = sortedNodes.length;
            const radius = 100; // Node circle radius
            // Generate N evenly distributed positions on a circle
            const positions = [];
            for (let i = 0; i < N; i++) {
                const angle = -Math.PI / 2 + (2 * Math.PI * i) / N;
                positions.push({
                    angle: angle,
                    x: radius * Math.cos(angle),
                    y: radius * Math.sin(angle),
                });
            }
            // Sort positions by angle distance from 12 o'clock
            positions.sort((a, b) => {
                const base = -Math.PI / 2;
                const diffA = Math.PI - Math.abs(Math.abs(a.angle - base) - Math.PI);
                const diffB = Math.PI - Math.abs(Math.abs(b.angle - base) - Math.PI);
                return diffB - diffA;
            });
            // Create node objects
            const customNodes = sortedNodes.map((node, i) => {
                const degree = nodeDegrees[node.id];
                let size = SIZES.nodeMinSize;
                if (maxDegree > minDegree) {
                    size =
                        SIZES.nodeMinSize +
                        ((degree - minDegree) * (SIZES.nodeMaxSize - SIZES.nodeMinSize)) /
                        (maxDegree - minDegree);
                }
                const inW = inWeights[node.id];
                const outW = outWeights[node.id];
                const total = inW + outW;
                balance = inW - outW;
                // Use plain object for node
                return {
                    id: node.id,
                    label: node.label,
                    x: positions[i].x,
                    y: positions[i].y,
                    degree,
                    balance,
                };
            });
            // Create edge objects
            const customEdges = graphData.edges.map((edge) => {
                let weight = 1;
                if (edge.label) {
                    const match = String(edge.label).match(/R=(\d+(\.\d+)?)/);
                    if (match) weight = parseFloat(match[1]);
                }
                return {
                    source: edge.source,
                    target: edge.target,
                    label: edge.label ? String(edge.label) : undefined,
                    weight,
                    type: "arrow",
                    arrowSize: SIZES.edgeMinWidth,
                };
            });
            // Add nodes and edges to graph
            customNodes.forEach((node) => {
                g.addNode(node.id, { ...node });
            });
            customEdges.forEach((edge) => {
                g.addEdge(edge.source, edge.target, { ...edge });
            });
            const renderer = new Sigma(g, sigmaContainer);
            // Highlight logic
            let highlightedNodes = new Set();
            let highlightedEdges = new Set();
            let selectedNodeLabel = null;
            let selectedNodeId = null;
            function clearHighlight() {
                highlightedNodes.clear();
                highlightedEdges.clear();
                renderer.refresh();
            }
            renderer.on("enterNode", ({ node }) => {
                const nodeAttr = g.getNodeAttributes(node);
                selectedNodeLabel = nodeAttr.label || node;
                selectedNodeId = node;
                highlightedNodes = new Set([node, ...g.neighbors(node)]);
                highlightedEdges = new Set(g.edges(node));
                renderer.refresh();
                // Show related edges in sidebar
                const infoList = g.edges(node).map((eid) => {
                    const data = g.getEdgeAttributes(eid);
                    const sourceLabel =
                        g.getNodeAttributes(data.source).label || data.source;
                    const targetLabel =
                        g.getNodeAttributes(data.target).label || data.target;
                    let at = "";
                    if (data.label) {
                        const m = String(data.label).match(/AT=(\d+(?:\.\d+)?)/);
                        if (m) at = m[1];
                    }
                    let r = "";
                    if (data.label) {
                        const m = String(data.label).match(/R=(\d+(?:\.\d+)?)/);
                        if (m) r = m[1];
                    }
                    return {
                        from: sourceLabel,
                        fromId: data.source,
                        to: targetLabel,
                        toId: data.target,
                        reactions: r,
                        atomTransfer: at,
                    };
                });
                updateSidebar(infoList);
            });
            renderer.on("clickStage", () => {
                selectedNodeLabel = null;
                selectedNodeId = null;
                clearHighlight();
                updateSidebar([]);
            });
            // Node/edge reducer for highlight effect
            renderer.setSetting("nodeReducer", (node, data) => {
                let dynamicColor = getBalanceColor(
                    data.balance,
                    minBalance,
                    maxBalance,
                    COLORS.balanceLow,
                    COLORS.balanceMid,
                    COLORS.balanceHigh
                );
                let dynamicSize = getNodeSize(
                    data.degree,
                    minDegree,
                    maxDegree,
                    SIZES.nodeMinSize,
                    SIZES.nodeMaxSize
                );

                if (!highlightedNodes.size) {
                    return {
                        ...data,
                        color: dynamicColor,
                        size: dynamicSize,
                    };
                }
                // Something seleted
                return {
                    ...data,
                    color: highlightedNodes.has(node) ? dynamicColor : COLORS.nodeDefault,
                    size: dynamicSize,
                    zIndex: highlightedNodes.has(node) ? 50 : 0,
                };
            });
            renderer.setSetting("edgeReducer", (edge, data) => {
                let weight = data.weight;
                let dynamicWidth = getEdgeWidth(
                    weight,
                    minEdgeWeight,
                    maxEdgeWeight,
                    SIZES.edgeMinWidth,
                    SIZES.edgeMaxWidth
                );
                let dynamicColor = getEdgeColor(
                    weight,
                    minEdgeWeight,
                    maxEdgeWeight,
                    COLORS.edgeLow,
                    COLORS.edgeHigh
                );

                if (!highlightedEdges.size) {
                    return {
                        ...data,
                        color: dynamicColor,
                        size: dynamicWidth,
                        weight,
                        label: "",
                        zIndex: weight,
                    };
                }
                // Something selected
                return {
                    ...data,
                    color: highlightedEdges.has(edge)
                        ? COLORS.edgeSelected
                        : COLORS.edgeDefault,
                    size: dynamicWidth,
                    zIndex: highlightedEdges.has(edge) ? 100 : 0,
                    weight,
                    label: highlightedEdges.has(edge) ? data.label : "",
                };
            });
            // Only show edge labels when highlighted
            renderer.setSetting("renderEdgeLabels", true);
            renderer.setSetting("zIndex", true);
            renderer.setSetting("edgeLabelSize", SIZES.edgeLabelSize);
            renderer.setSetting("edgeLabelColor", "inherit");
            renderer.setSetting("labelSize", SIZES.nodeLabelFontSize);
            // Edge color when not selected
            renderer.setSetting("edgeColor", (edge, data) => {
                return getEdgeColor(
                    data.weight,
                    minEdgeWeight,
                    maxEdgeWeight,
                    COLORS.edgeLow,
                    COLORS.edgeHigh
                );
            });
            renderer.setSetting("edgeLabelReducer", () => null);
            if (chartsMap) chartsMap.set(chartId, true);
        } catch (error) {
            console.error("Sigma.js rendering failed:", error);
            sigmaContainer.innerHTML =
                '<p class="sigma-error-message">Sigma.js rendering failed</p>';
        }
    } else {
        sigmaContainer.innerHTML =
            '<p class="sigma-error-message">Sigma.js not loaded</p>';
    }
    return true;

    /**
     * Update the sidebar with edge/node info.
     * @param {Array} infoList - List of info objects.
     */
    function updateSidebar(infoList) {
        const infoBox = container.querySelector(`#${chartId}-info`);
        if (!infoBox) return;
        // Render a table for info data
        function renderTable(data) {
            const table = document.createElement("table");
            table.className = "key-mol-table";
            const thead = document.createElement("thead");
            const headerRow = document.createElement("tr");
            ["From", "To", "Reactions", "Atom transfer"].forEach((h) => {
                const th = document.createElement("th");
                th.textContent = h;
                headerRow.appendChild(th);
            });
            thead.appendChild(headerRow);
            table.appendChild(thead);
            const tbody = document.createElement("tbody");
            data.forEach((row) => {
                let from = "",
                    to = "",
                    reactions = "",
                    atomTransfer = "";
                if (typeof row === "string") {
                    const match = row.match(
                        /^(.*?) -> (.*?)\s+R=(\d+(?:\.\d+)?)\s*AT=(\d+(?:\.\d+)?)/
                    );
                    if (match) {
                        from = match[1];
                        to = match[2];
                        reactions = match[3];
                        atomTransfer = match[4];
                    }
                } else if (typeof row === "object") {
                    from = row.from || "";
                    to = row.to || "";
                    reactions = row.reactions || "";
                    atomTransfer = row.atomTransfer || "";
                }
                const tr = document.createElement("tr");
                [from, to, reactions, atomTransfer].forEach((val) => {
                    const td = document.createElement("td");
                    td.textContent = val;
                    tr.appendChild(td);
                });
                tbody.appendChild(tr);
            });
            table.appendChild(tbody);
            return table;
        }
        // Render infoList as two-column tables
        const leftCol = infoBox.querySelector(".network-info-col.left");
        const rightCol = infoBox.querySelector(".network-info-col.right");
        leftCol.innerHTML = "";
        rightCol.innerHTML = "";
        if (Array.isArray(infoList) && infoList.length > 0) {
            const mid = Math.ceil(infoList.length / 2);
            const leftData = infoList.slice(0, mid);
            const rightData = infoList.slice(mid);
            if (leftData.length > 0) {
                leftCol.appendChild(renderTable(leftData));
                leftCol.classList.remove("hidden");
                leftCol.classList.add("visible");
            } else {
                leftCol.classList.remove("visible");
                leftCol.classList.add("hidden");
            }
            if (rightData.length > 0) {
                rightCol.appendChild(renderTable(rightData));
                rightCol.classList.remove("hidden");
                rightCol.classList.add("visible");
            } else {
                rightCol.classList.remove("visible");
                rightCol.classList.add("hidden");
            }
        } else {
            leftCol.classList.remove("visible");
            leftCol.classList.add("hidden");
            rightCol.classList.remove("visible");
            rightCol.classList.add("hidden");
        }
    }
}

// UMD-style export for browser global usage
if (typeof window !== "undefined") {
    window.reaxToolsPlotNetwork = {
        createSigmaGraph,
    };
}
