// Network graph (DOT) rendering utilities for ReaxTools Web
// Only code related to DOT/Graphviz network visualization is here

// ===== Color configuration for easy theme adjustment =====
const COLORS = {
    title: '#34495e',
    border: '#3498db',
    background: '#fff',
    nodeDefault: '#eee',
    sigmaFail: '#666',
    balanceLow: '#FF4F0F',    // reactant
    balanceMid: '#06923E',    // intermediate
    balanceHigh: '#3674B5',   // product
    edgeLow: '#EEEEEE',        // min edge weight color
    edgeHigh: '#333333',       // max edge weight color
    nodeSelected: '#ff9800',   // selected node color
    edgeSelected: '#e53935',   // selected edge color
    edgeDefault: '#E6E6E6'
};

// ===== Size configuration for easy adjustment =====
const SIZES = {
    containerHeight: 1000,      // px
    borderRadius: 5,           // px
    borderWidth: 2,            // px
    edgeWidthScale: 1.5,       // edge width scale factor
    edgeMinWidth: 6,           // px, min edge width
    edgeLabelSize: 18,         // px
    nodeMinSize: 12,           // px, min node size
    nodeMaxSize: 36,           // px, max node size
    nodeLabelFontSize: 24,     // px, node label font size
};

// 兼容 sigma.js UMD 导出
if (typeof window.sigma === 'undefined' && typeof window.Sigma !== 'undefined') {
    window.sigma = window.Sigma;
}

// ===== Custom Node and Edge Classes =====
class CustomNode {
    constructor(obj) {
        Object.assign(this, obj);
    }
}
class CustomEdge {
    constructor(obj) {
        Object.assign(this, obj);
    }
}

function createSigmaGraph(graphData, chartId, title = 'Reaction Network', chartsMap) {
    // graphData should be { nodes: [...], edges: [...] }
    const container = document.getElementById(chartId);
    if (!container) {
        console.error(`Container with id '${chartId}' not found`);
        return false;
    }
    container.innerHTML = '';
    // Title
    const titleDiv = document.createElement('p');
    titleDiv.textContent = title;
    titleDiv.style.textAlign = 'center';
    titleDiv.style.marginBottom = '15px';
    titleDiv.style.marginTop = '0';
    titleDiv.style.color = COLORS.title;
    titleDiv.style.fontSize = '1.1rem';
    titleDiv.style.fontWeight = 'bold';
    container.appendChild(titleDiv);
    // Sigma container
    const sigmaContainer = document.createElement('div');
    sigmaContainer.style.cssText = `width: 100%;height: ${SIZES.containerHeight}px;overflow: hidden;border: ${SIZES.borderWidth}px solid ${COLORS.border};border-radius: ${SIZES.borderRadius}px;background-color: ${COLORS.background};padding: 0;margin-bottom: 10px;position: relative;`;
    sigmaContainer.id = chartId + '_sigma';
    container.appendChild(sigmaContainer);
    // 在 sigmaContainer 之后插入 infoBox
    let infoBox = container.querySelector(`#${chartId}-info`);
    if (!infoBox) {
        infoBox = document.createElement('div');
        infoBox.id = `${chartId}-info`;
        infoBox.style.cssText = `
            display: flex;
            flex-direction: row;
            gap: 20px;
            width: 100%;
            margin-top: 10px;
        `;
        // 左右两栏
        const leftCol = document.createElement('div');
        leftCol.className = 'network-info-col left';
        leftCol.style.cssText = 'flex:1; min-width:0; overflow-x:auto;';
        infoBox.appendChild(leftCol);
        const rightCol = document.createElement('div');
        rightCol.className = 'network-info-col right';
        rightCol.style.cssText = 'flex:1; min-width:0; overflow-x:auto;';
        infoBox.appendChild(rightCol);
        container.appendChild(infoBox);
    }
    // 新增：在 infoBox 下方插入说明文字
    let infoNote = container.querySelector(`#${chartId}-info-note`);
    if (!infoNote) {
        infoNote = document.createElement('div');
        infoNote.id = `${chartId}-info-note`;
        infoNote.style.cssText = 'margin-top: 8px; color: #888; font-size: 0.95em;';
        infoNote.textContent = '节点颜色：分子相对流入/流出比，流出占比越高越红，流入占比越高越蓝。节点大小：总流入流出数。\n边的颜色：反应的发生次数。次数越多颜色越深。\n点击特定节点高亮相关反应和上下游分子，并展示表格。';
        container.appendChild(infoNote);
    }
    // Sigma.js rendering
    if (window.sigma && window.graphology) {
        try {
            const Graph = window.graphology;
            const Sigma = window.sigma;
            const g = new Graph();
            // ===== 1. 计算所有动态缩放属性的 min/max =====
            let minEdgeWeight = Infinity, maxEdgeWeight = -Infinity;
            let minDegree = Infinity, maxDegree = -Infinity;
            let minBalance = Infinity, maxBalance = -Infinity;
            const inWeights = {}, outWeights = {};
            graphData.nodes.forEach(node => {
                inWeights[node.id] = 0;
                outWeights[node.id] = 0;
            });
            graphData.edges.forEach(edge => {
                let weight = 1;
                if (edge.label) {
                    const match = String(edge.label).match(/R=(\d+(\.\d+)?)/);
                    if (match) weight = parseFloat(match[1]);
                }
                weight = Number(weight);
                // 调试输出每条边的label和解析出的weight
                outWeights[edge.source] += weight;
                inWeights[edge.target] += weight;
                if (weight < minEdgeWeight) minEdgeWeight = weight;
                if (weight > maxEdgeWeight) maxEdgeWeight = weight;
            });
            const nodeBalances = {}, nodeDegrees = {};
            graphData.nodes.forEach(node => {
                const inW = inWeights[node.id];
                const outW = outWeights[node.id];
                const total = inW + outW;
                let balance = 0.5;
                if (total > 0) balance = inW / total;
                nodeBalances[node.id] = balance;
                if (balance < minBalance) minBalance = balance;
                if (balance > maxBalance) maxBalance = balance;
                const degree = inW + outW;
                nodeDegrees[node.id] = degree;
                if (degree < minDegree) minDegree = degree;
                if (degree > maxDegree) maxDegree = degree;
            });

            // 2. 计算平衡度并设置颜色
            function getBalanceColor(balance, lowColor, midColor, highColor) {
                // Helper: hex to rgb
                function hexToRgb(hex) {
                    hex = hex.replace('#', '');
                    if (hex.length === 3) hex = hex.split('').map(x => x + x).join('');
                    const num = parseInt(hex, 16);
                    return [num >> 16 & 255, num >> 8 & 255, num & 255];
                }
                // Helper: rgb to hex
                function rgbToHex([r, g, b]) {
                    return '#' + [r, g, b].map(x => x.toString(16).padStart(2, '0')).join('');
                }
                let colorA, colorB, t;
                if (balance <= 0.5) {
                    colorA = hexToRgb(lowColor);
                    colorB = hexToRgb(midColor);
                    t = balance / 0.5;
                } else {
                    colorA = hexToRgb(midColor);
                    colorB = hexToRgb(highColor);
                    t = (balance - 0.5) / 0.5;
                }
                const rgb = colorA.map((a, i) => Math.round(a * (1 - t) + colorB[i] * t));
                return rgbToHex(rgb);
            }

            // 边颜色插值函数
            function getEdgeColor(weight, lowColor, highColor) {
                // Helper: hex to rgb
                function hexToRgb(hex) {
                    hex = hex.replace('#', '');
                    if (hex.length === 3) hex = hex.split('').map(x => x + x).join('');
                    const num = parseInt(hex, 16);
                    return [num >> 16 & 255, num >> 8 & 255, num & 255];
                }
                // Helper: rgb to hex
                function rgbToHex([r, g, b]) {
                    return '#' + [r, g, b].map(x => x.toString(16).padStart(2, '0')).join('');
                }
                weight = Number(weight);
                let t = 0;
                if (maxEdgeWeight > minEdgeWeight && weight > 0 && minEdgeWeight > 0) {
                    t = (Math.log2(weight) - Math.log2(minEdgeWeight)) / (Math.log2(maxEdgeWeight) - Math.log2(minEdgeWeight));
                    t = Math.max(0, Math.min(1, t));
                }
                // 调试输出每次调用的weight和t
                const rgbA = hexToRgb(lowColor);
                const rgbB = hexToRgb(highColor);
                const rgb = rgbA.map((a, i) => Math.round(a * (1 - t) + rgbB[i] * t));
                return rgbToHex(rgb);
            }

            // 按平衡度排序节点
            const sortedNodes = [...graphData.nodes].sort((a, b) => nodeBalances[a.id] - nodeBalances[b.id]);
            const N = sortedNodes.length;
            const radius = 100; // You can also move this to SIZES if you want to adjust the node circle radius

            // 生成N个均分圆上的角度，0号为12点（-90°），顺时针
            const positions = [];
            for (let i = 0; i < N; i++) {
                const angle = -Math.PI / 2 + (2 * Math.PI * i) / N;
                positions.push({
                    angle: angle,
                    x: radius * Math.cos(angle),
                    y: radius * Math.sin(angle)
                });
            }
            // 按与12点方向的角度距离排序（12点为-90°，6点为+90°）
            positions.sort((a, b) => {
                const base = -Math.PI / 2; // 12点方向
                const diffA = Math.PI - Math.abs(Math.abs(a.angle - base) - Math.PI);
                const diffB = Math.PI - Math.abs(Math.abs(b.angle - base) - Math.PI);
                return diffB - diffA;
            });

            // 生成CustomNode实例
            const customNodes = sortedNodes.map((node, i) => {
                const degree = nodeDegrees[node.id];
                let size = SIZES.nodeMinSize;
                if (maxDegree > minDegree) {
                    size = SIZES.nodeMinSize + (degree - minDegree) * (SIZES.nodeMaxSize - SIZES.nodeMinSize) / (maxDegree - minDegree);
                }
                const inW = inWeights[node.id];
                const outW = outWeights[node.id];
                const total = inW + outW;
                let balance = 0.5;
                if (total > 0) balance = inW / total;
                let t = 0.5;
                if (maxBalance > minBalance) {
                    t = (balance - minBalance) / (maxBalance - minBalance);
                }
                const color = getBalanceColor(t, COLORS.balanceLow, COLORS.balanceMid, COLORS.balanceHigh);
                // Remove labelAlignment and angle
                return new CustomNode({
                    id: node.id,
                    label: node.label,
                    x: positions[i].x,
                    y: positions[i].y,
                    degree,
                    balance,
                    color,
                    size
                });
            });

            // 生成CustomEdge实例
            const customEdges = graphData.edges.map(edge => {
                let weight = 1;
                if (edge.label) {
                    const match = String(edge.label).match(/R=(\d+(\.\d+)?)/);
                    if (match) weight = parseFloat(match[1]);
                }
                return new CustomEdge({
                    source: edge.source,
                    target: edge.target,
                    label: edge.label ? String(edge.label) : undefined,
                    weight,
                    type: 'arrow',
                    arrowSize: SIZES.edgeMinWidth
                });
            });

            // 批量添加节点和边
            customNodes.forEach(node => {
                g.addNode(node.id, { ...node });
            });
            customEdges.forEach(edge => {
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
            renderer.on('clickNode', ({ node }) => {
                const nodeAttr = g.getNodeAttributes(node);
                selectedNodeLabel = nodeAttr.label || node;
                selectedNodeId = node;
                highlightedNodes = new Set([node, ...g.neighbors(node)]);
                highlightedEdges = new Set(g.edges(node));
                renderer.refresh();
                // 侧边栏显示与该节点相关的边，节点id换成label
                const infoList = g.edges(node).map(eid => {
                    const data = g.getEdgeAttributes(eid);
                    const sourceLabel = g.getNodeAttributes(data.source).label || data.source;
                    const targetLabel = g.getNodeAttributes(data.target).label || data.target;
                    let at = '';
                    if (data.label) {
                        const m = String(data.label).match(/AT=(\d+(?:\.\d+)?)/);
                        if (m) at = m[1];
                    }
                    let r = '';
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
                        atomTransfer: at
                    };
                });
                updateSidebar(infoList);
            });
            renderer.on('clickEdge', ({ edge }) => {
                selectedNodeLabel = null;
                selectedNodeId = null;
                const source = g.source(edge);
                const target = g.target(edge);
                highlightedNodes = new Set([source, target]);
                highlightedEdges = new Set([edge]);
                renderer.refresh();
                const data = g.getEdgeAttributes(edge);
                const sourceLabel = g.getNodeAttributes(data.source).label || data.source;
                const targetLabel = g.getNodeAttributes(data.target).label || data.target;
                let at = '';
                if (data.label) {
                    const m = String(data.label).match(/AT=(\d+(?:\.\d+)?)/);
                    if (m) at = m[1];
                }
                let r = '';
                if (data.label) {
                    const m = String(data.label).match(/R=(\d+(?:\.\d+)?)/);
                    if (m) r = m[1];
                }
                updateSidebar([{ from: sourceLabel, fromId: data.source, to: targetLabel, toId: data.target, reactions: r, atomTransfer: at }]);
            });
            renderer.on('clickStage', () => {
                selectedNodeLabel = null;
                selectedNodeId = null;
                clearHighlight();
                updateSidebar([]);
            });
            // Node/edge reducer for highlight effect
            renderer.setSetting('nodeReducer', (node, data) => {
                if (!highlightedNodes.size) return data;
                return Object.assign({}, data, {
                    color: highlightedNodes.has(node) ? data.color : COLORS.nodeDefault,
                    zIndex: highlightedNodes.has(node) ? 100 : 0
                });
            });
            renderer.setSetting('edgeReducer', (edge, data) => {
                // 直接用 data.weight
                let weight = data.weight;
                if (!highlightedEdges.size) {
                    // Always return dynamic color when not highlighted
                    return {
                        ...data,
                        color: getEdgeColor(weight, COLORS.edgeLow, COLORS.edgeHigh),
                        size: SIZES.edgeMinWidth,
                        weight // 保证 weight 字段传递下去
                    };
                }
                // Only return highlight color when highlighting
                return {
                    ...data,
                    color: highlightedEdges.has(edge) ? COLORS.edgeSelected : COLORS.edgeDefault,
                    size: SIZES.edgeMinWidth,
                    weight // 保证 weight 字段传递下去
                };
            });
            // 只在高亮时显示边label
            renderer.setSetting('renderEdgeLabels', false);
            renderer.setSetting('edgeLabelSize', SIZES.edgeLabelSize);
            renderer.setSetting('edgeLabelColor', 'inherit');
            renderer.setSetting('nodeLabelSize', SIZES.nodeLabelFontSize);

            // 未选中时的边颜色
            renderer.setSetting('edgeColor', (edge, data) => {
                // 用 data.weight 作为权重
                return getEdgeColor(data.weight, COLORS.edgeLow, COLORS.edgeHigh);
            });

            renderer.setSetting('edgeLabelReducer', () => null);

            if (chartsMap) chartsMap.set(chartId, true);
        } catch (error) {
            console.error('Sigma.js rendering failed:', error);
            sigmaContainer.innerHTML = '<p style="text-align: center; padding: 20px; color: ' + COLORS.sigmaFail + '; font-style: italic;">Sigma.js 渲染失败</p>';
        }
    } else {
        sigmaContainer.innerHTML = '<p style="text-align: center; padding: 20px; color: ' + COLORS.sigmaFail + '; font-style: italic;">Sigma.js 未加载</p>';
    }
    return true;

    function updateSidebar(infoList) {
        const infoBox = container.querySelector(`#${chartId}-info`);
        if (!infoBox) return;
        // 渲染表格函数
        function renderTable(data) {
            const table = document.createElement('table');
            table.className = 'key-mol-table';
            const thead = document.createElement('thead');
            const headerRow = document.createElement('tr');
            ['From', 'To', 'Reactions', 'Atom transfer'].forEach(h => {
                const th = document.createElement('th');
                th.textContent = h;
                headerRow.appendChild(th);
            });
            thead.appendChild(headerRow);
            table.appendChild(thead);
            const tbody = document.createElement('tbody');
            data.forEach(row => {
                // row 可能是字符串或对象
                let from = '', to = '', reactions = '', atomTransfer = '';
                if (typeof row === 'string') {
                    const match = row.match(/^(.*?) -> (.*?)\s+R=(\d+(?:\.\d+)?)\s*AT=(\d+(?:\.\d+)?)/);
                    if (match) {
                        from = match[1];
                        to = match[2];
                        reactions = match[3];
                        atomTransfer = match[4];
                    }
                } else if (typeof row === 'object') {
                    from = row.from || '';
                    to = row.to || '';
                    reactions = row.reactions || '';
                    atomTransfer = row.atomTransfer || '';
                }
                const tr = document.createElement('tr');
                [from, to, reactions, atomTransfer].forEach(val => {
                    const td = document.createElement('td');
                    td.textContent = val;
                    tr.appendChild(td);
                });
                tbody.appendChild(tr);
            });
            table.appendChild(tbody);
            return table;
        }
        // 渲染所有 infoList 为两栏表格
        const leftCol = infoBox.querySelector('.network-info-col.left');
        const rightCol = infoBox.querySelector('.network-info-col.right');
        leftCol.innerHTML = '';
        rightCol.innerHTML = '';
        if (Array.isArray(infoList) && infoList.length > 0) {
            const mid = Math.ceil(infoList.length / 2);
            const leftData = infoList.slice(0, mid);
            const rightData = infoList.slice(mid);
            if (leftData.length > 0) {
                leftCol.appendChild(renderTable(leftData));
                leftCol.style.display = '';
            } else {
                leftCol.style.display = 'none';
            }
            if (rightData.length > 0) {
                rightCol.appendChild(renderTable(rightData));
                rightCol.style.display = '';
            } else {
                rightCol.style.display = 'none';
            }
        } else {
            leftCol.style.display = 'none';
            rightCol.style.display = 'none';
        }
    }
}



// UMD-style export for browser global usage
if (typeof window !== 'undefined') {
    window.reaxToolsPlotNetwork = {
        createSigmaGraph
    };
}


