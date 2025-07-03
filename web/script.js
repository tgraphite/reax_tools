// 全局变量：WASM 模块实例和上传的文件对象
let wasmModule = null; // 存储加载后的 WASM 模块
let uploadedFile = null; // 存储用户上传的文件

// 监听文件选择事件
// 当用户选择文件后，保存文件对象，并在页面上显示文件名和大小
// fileInput 是 <input type="file" id="fileInput"> 元素
// output 是用于显示信息的 <div id="output"> 元素

document.getElementById('fileInput').addEventListener('change', (e) => {
    uploadedFile = e.target.files[0]; // 获取用户选择的文件
    if (!uploadedFile) return;
    document.getElementById('output').innerHTML =
        `已选择文件: ${uploadedFile.name} (${(uploadedFile.size / 1024).toFixed(2)} KB)`;
});

// "选择示例文件" 链接的点击事件
document.getElementById('loadSampleBtn').addEventListener('click', async (e) => {
    e.preventDefault(); // 阻止链接的默认行为
    try {
        const response = await fetch('./reaxff.xyz');
        if (!response.ok) {
            throw new Error(`无法加载示例文件: ${response.statusText}`);
        }
        const fileContent = await response.blob();
        uploadedFile = new File([fileContent], "reaxff.xyz", { type: "text/plain" });

        if (!uploadedFile) return;
        document.getElementById('output').innerHTML =
            `已选择文件: ${uploadedFile.name} (${(uploadedFile.size / 1024).toFixed(2)} KB)`;

        // 清空文件输入框，避免混淆
        document.getElementById('fileInput').value = '';

    } catch (error) {
        alert(`加载示例文件失败: ${error.message}`);
        console.error(error);
    }
});

// 加载 WASM 模块的函数
// 这里没有用 MODULARIZE，直接等待 onRuntimeInitialized 回调
function loadWasmModule() {
    return new Promise((resolve, reject) => {
        // 如果已经加载过 Module，直接返回
        if (wasmModule) {
            resolve(wasmModule);
            return;
        }

        const logOutput = document.getElementById('log-output');
        // 配置 Module 对象，设置回调和文件定位
        window.Module = {
            print: text => {
                console.log(`[WASM]: ${text}`);
                if (logOutput) {
                    // 处理 \r 覆盖
                    if (text.includes('\r')) {
                        // 拆分成多段（有可能一行有多个\r）
                        const parts = text.split('\r');
                        // 取最后一段作为最新内容
                        const last = parts[parts.length - 1];
                        // 取当前内容，去掉最后一行
                        let lines = logOutput.textContent.split('\n');
                        // 如果最后一行为空，去掉
                        if (lines.length > 0 && lines[lines.length - 1] === '') lines.pop();
                        // 替换最后一行
                        lines[lines.length - 1] = last;
                        logOutput.textContent = lines.join('\n') + '\n';
                    } else {
                        logOutput.textContent += text + '\n';
                    }
                    logOutput.scrollTop = logOutput.scrollHeight; // 自动滚到最下方
                }
            },
            printErr: text => {
                console.error(`[WASM Error]: ${text}`);
                if (logOutput) {
                    logOutput.textContent += `[ERROR] ${text}\n`;
                    logOutput.scrollTop = logOutput.scrollHeight; // 自动滚到最下方
                }
            },
            onRuntimeInitialized: function () {
                wasmModule = this; // 保存模块实例
                // 确保 /output 目录存在（WASM 虚拟文件系统）
                if (!wasmModule.FS.analyzePath('/output').exists) {
                    wasmModule.FS.mkdir('/output');
                }
                resolve(wasmModule); // 加载完成，返回 Module
            },
            // 定位 .wasm 文件的路径
            locateFile: (path) => {
                if (path.endsWith('.wasm')) return './wasm_main.wasm';
                return path;
            }
        };
        // 动态加载 wasm_main.js（Emscripten 生成的 glue 代码）
        const script = document.createElement('script');
        script.src = './wasm_main.js';
        script.onerror = reject;
        document.head.appendChild(script);
    });
}

// 读取CSV文件并创建图表
function createChartsFromFiles(files) {
    if (!window.reaxToolsPlotter) {
        console.error('Plotter not initialized');
        return;
    }

    console.log('Creating charts from files:', files.map(f => f.name));

    // Initialize plotter with charts area
    window.reaxToolsPlotter.init('charts-area');

    // Create charts for each file
    const chartFiles = files.map(file => ({
        name: file.name,
        content: file.content
    }));

    window.reaxToolsPlotter.createMultipleCharts(chartFiles, 'charts-area');
}

// 运行按钮点击事件
// runBtn 是 <button id="runBtn"> 元素
// output 是用于显示结果的 <div id="output"> 元素

document.getElementById('runBtn').onclick = async function () {
    if (!uploadedFile) {
        alert('请先选择一个文件');
        return;
    }
    const outputDiv = document.getElementById('output');
    const logOutput = document.getElementById('log-output');
    outputDiv.innerHTML = "处理中...";
    logOutput.textContent = ''; // Clear previous output

    // Clear previous charts
    if (window.reaxToolsPlotter) {
        window.reaxToolsPlotter.clearCharts();
    }

    try {
        // 加载 WASM 模块（只加载一次）
        wasmModule = await loadWasmModule();

        // 清理旧的输出文件
        if (wasmModule.FS.analyzePath('/output').exists) {
            const oldFiles = wasmModule.FS.readdir('/output');
            oldFiles.forEach(file => {
                if (file !== '.' && file !== '..') {
                    wasmModule.FS.unlink(`/output/${file}`);
                }
            });
        }
        
        // 将上传的文件写入 WASM 虚拟文件系统（FS）
        const fileContent = new Uint8Array(await uploadedFile.arrayBuffer());
        const inputFilePath = uploadedFile.name; // 保持原名
        wasmModule.FS.writeFile(inputFilePath, fileContent);

        // 构造命令行参数
        const userArgsText = document.getElementById('argsInput').value.trim();
        const userArgs = userArgsText ? userArgsText.split(/\s+/).filter(Boolean) : [];
        const args = ['program', '--traj', inputFilePath, ...userArgs]; // -f 或 --traj 都可以
        
        // 为每个参数分配内存，并写入字符串
        const argvPtrs = args.map(s => {
            const len = wasmModule.lengthBytesUTF8(s) + 1; // 字符串长度+1
            const ptr = wasmModule._malloc(len); // 分配内存
            wasmModule.stringToUTF8(s, ptr, len); // 写入字符串
            return ptr;
        });
        // 分配 argv 数组的内存（每个指针4字节）
        const argvPtr = wasmModule._malloc(argvPtrs.length * 4);
        // 写入每个参数的指针到 argvPtr
        argvPtrs.forEach((ptr, i) => {
            wasmModule.setValue(argvPtr + i * 4, ptr, 'i32');
        });

        console.log('javascript send args:', args);

        // 调用 C++ 的 cpp_main 函数
        // 返回值为处理结果（0=成功，非0=失败）
        const ret = wasmModule.ccall(
            'cpp_main', // C++ 函数名
            'number',   // 返回类型
            ['number', 'number'], // 参数类型
            [args.length, argvPtr] // 参数值
        );
        console.log('cpp_main 返回值:', ret);
        

        // 释放分配的内存
        argvPtrs.forEach(ptr => wasmModule._free(ptr));
        wasmModule._free(argvPtr);

        // 读取处理结果文件
        const outputFiles = wasmModule.FS.readdir('/output').filter(p => p !== '.' && p !== '..');

        // 获取 log-output 内容
        const logOutputText = logOutput.textContent || '';

        // 获取东八区时间字符串
        function getBeijingTimeISO() {
            const now = new Date();
            const beijingOffset = 8 * 60; // 东八区
            const localOffset = now.getTimezoneOffset();
            const beijingTime = new Date(now.getTime() + (beijingOffset + localOffset) * 60000);
            return beijingTime.toISOString().replace('T', ' ').replace('Z', '');
        }
        const beijingTime = getBeijingTimeISO();

        // 发送日志（只在 WASM 处理完成后）
        await fetch('/api/log', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ time: beijingTime, logOutput: logOutputText })
        });

        if (outputFiles.length === 0) {
            outputDiv.innerHTML = `<p>处理完成，但没有生成输出文件。</p>`;
        } else {
            // 分离用于绘图的CSV文件和DOT文件，以及其他所有文件
            const plotFiles = [];
            const allFiles = [];

            // 只绘制species_count.csv bond_count.csv文件和reactions.dot文件等，其他不绘制。
            outputFiles.forEach(fileName => {
                console.log(`Processing output file: ${fileName}`);
                if (fileName.toLowerCase().endsWith('count.csv') || fileName.toLowerCase().endsWith('reactions.dot') || fileName === 'key_molecules_reactions.csv') {
                    try {
                        const fileData = wasmModule.FS.readFile(`/output/${fileName}`);
                        const content = new TextDecoder().decode(fileData);
                        plotFiles.push({ name: fileName, content: content });
                        console.log(`Added ${fileName} for visualization`);
                    } catch (error) {
                        console.error(`Error reading file ${fileName}:`, error);
                    }
                }
                // 所有文件都添加到下载列表中
                allFiles.push(fileName);
            });

            console.log(`Found ${plotFiles.length} files for visualization:`, plotFiles.map(f => f.name));

            // 创建下载链接 - 所有文件都打包
            let downloadLinks = '';
            if (allFiles.length === 1) {
                // 单个文件，直接提供下载
                const fileName = allFiles[0];
                const resultData = wasmModule.FS.readFile(`/output/${fileName}`);
                const blob = new Blob([resultData], { type: 'application/octet-stream' });
                const url = URL.createObjectURL(blob);
                downloadLinks = `<a href="${url}" download="${fileName}">另存为结果文件 (${fileName})</a>`;
            } else {
                // 多个文件，打包成 zip（包含所有文件）
                const zip = new JSZip();
                allFiles.forEach(fileName => {
                    const fileData = wasmModule.FS.readFile(`/output/${fileName}`);
                    zip.file(fileName, fileData);
                });

                const zipBlob = await zip.generateAsync({ type: 'blob' });
                const url = URL.createObjectURL(zipBlob);
                const dateStr = new Date().toISOString().replace(/:/g, '-').split('.')[0];
                downloadLinks = `<a href="${url}" download="output-${dateStr}.zip">另存所有结果文件 (.zip)</a>`;
            }

            // 显示处理结果
            let resultHtml = `<p>处理完成！`;
            resultHtml += `一共生成并打包了 ${allFiles.length} 个文件。`;
            if (plotFiles.length > 0) {
                resultHtml += `自动生成了 ${plotFiles.length} 个页面图表。`;
            }
            resultHtml += `</p>`;
            
            if (downloadLinks) {
                resultHtml += downloadLinks;
            }

            outputDiv.innerHTML = resultHtml;

            // 如果有用于绘图的文件，创建图表
            if (plotFiles.length > 0) {
                createChartsFromFiles(plotFiles);
            }
        }

    } catch (error) {
        outputDiv.innerHTML = `<p style="color:red">出错: ${error.toString()}</p>`;
        console.error(error);
    }
};

// 添加图表控制按钮事件监听器
document.addEventListener('DOMContentLoaded', function() {
    // 清除所有图表按钮
    document.getElementById('clearChartsBtn').addEventListener('click', function() {
        if (window.reaxToolsPlotter) {
            window.reaxToolsPlotter.clearCharts();
        }
    });

    // 导出所有图表按钮
    document.getElementById('exportAllChartsBtn').addEventListener('click', function() {
        if (window.reaxToolsPlotter && window.reaxToolsPlotter.charts) {
            window.reaxToolsPlotter.charts.forEach((_, chartId) => {
                window.reaxToolsPlotter.exportChart(chartId, `${chartId}.png`);
            });
        }
    });
});
