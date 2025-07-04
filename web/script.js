//
let wasmModule = null; //
let uploadedFile = null; //

// Storage keys for localStorage
const STORAGE_KEYS = {
    CALCULATION_RESULTS: 'reaxTools_calculation_results',
    UPLOADED_FILE: 'reaxTools_uploaded_file',
    USER_ARGS: 'reaxTools_user_args',
    LOG_OUTPUT: 'reaxTools_log_output'
};

// Save calculation results to localStorage
function saveCalculationResults(outputFiles, plotFiles, logOutputText, userArgs) {
    try {
        const results = {
            timestamp: Date.now(),
            outputFiles: outputFiles,
            plotFiles: plotFiles,
            logOutput: logOutputText,
            userArgs: userArgs
        };
        localStorage.setItem(STORAGE_KEYS.CALCULATION_RESULTS, JSON.stringify(results));
        console.log('Calculation results saved to localStorage');
    } catch (error) {
        console.error('Failed to save calculation results:', error);
    }
}

// Load calculation results from localStorage
function loadCalculationResults() {
    try {
        const saved = localStorage.getItem(STORAGE_KEYS.CALCULATION_RESULTS);
        if (saved) {
            return JSON.parse(saved);
        }
    } catch (error) {
        console.error('Failed to load calculation results:', error);
    }
    return null;
}

// Save uploaded file info to localStorage
function saveUploadedFileInfo(file) {
    try {
        const fileInfo = {
            name: file.name,
            size: file.size,
            type: file.type,
            lastModified: file.lastModified
        };
        localStorage.setItem(STORAGE_KEYS.UPLOADED_FILE, JSON.stringify(fileInfo));
    } catch (error) {
        console.error('Failed to save uploaded file info:', error);
    }
}

// Load uploaded file info from localStorage
function loadUploadedFileInfo() {
    try {
        const saved = localStorage.getItem(STORAGE_KEYS.UPLOADED_FILE);
        if (saved) {
            return JSON.parse(saved);
        }
    } catch (error) {
        console.error('Failed to load uploaded file info:', error);
    }
    return null;
}

// Save user arguments to localStorage
function saveUserArgs(args) {
    try {
        localStorage.setItem(STORAGE_KEYS.USER_ARGS, args);
    } catch (error) {
        console.error('Failed to save user args:', error);
    }
}

// Load user arguments from localStorage
function loadUserArgs() {
    try {
        return localStorage.getItem(STORAGE_KEYS.USER_ARGS) || '';
    } catch (error) {
        console.error('Failed to load user args:', error);
        return '';
    }
}

// Clear all stored data
function clearStoredData() {
    try {
        // Clear localStorage
        Object.values(STORAGE_KEYS).forEach(key => {
            localStorage.removeItem(key);
        });

        // Clear virtual file system if WASM module is available
        if (wasmModule && wasmModule.FS) {
            // Clear output directory
            if (wasmModule.FS.analyzePath('/output').exists) {
                const outputFiles = wasmModule.FS.readdir('/output');
                outputFiles.forEach(file => {
                    if (file !== '.' && file !== '..') {
                        try {
                            wasmModule.FS.unlink(`/output/${file}`);
                        } catch (error) {
                            console.warn(`Failed to delete output file ${file}:`, error);
                        }
                    }
                });
            }

            // Clear input files (files uploaded by user)
            if (uploadedFile) {
                try {
                    wasmModule.FS.unlink(uploadedFile.name);
                } catch (error) {
                    console.warn(`Failed to delete input file ${uploadedFile.name}:`, error);
                }
                uploadedFile = null; // Reset uploaded file reference
            }

            // Clear any other files that might exist in root directory
            const rootFiles = wasmModule.FS.readdir('/');
            rootFiles.forEach(file => {
                if (file !== '.' && file !== '..' && file !== 'output' && file !== 'dev' && file !== 'proc' && file !== 'tmp') {
                    try {
                        const stat = wasmModule.FS.stat(`/${file}`);
                        if (stat.isFile()) {
                            wasmModule.FS.unlink(`/${file}`);
                        }
                    } catch (error) {
                        console.warn(`Failed to delete file ${file}:`, error);
                    }
                }
            });
        }

        // Clear chart data and charts
        if (window.reaxToolsPlotter) {
            window.reaxToolsPlotter.clearChartData();
            window.reaxToolsPlotter.clearCharts();
        }

        console.log('All stored data cleared including virtual file system and charts');
    } catch (error) {
        console.error('Failed to clear stored data:', error);
    }
}

// Restore calculation results and display them
async function restoreCalculationResults() {
    const results = loadCalculationResults();
    if (!results) {
        return false;
    }

    const outputDiv = document.getElementById('output');
    const logOutput = document.getElementById('log-output');
    const argsInput = document.getElementById('argsInput');

    // Restore log output
    if (results.logOutput) {
        logOutput.textContent = results.logOutput;
    }

    // Restore user arguments
    if (results.userArgs) {
        argsInput.value = results.userArgs;
    }

    // Display results
    if (results.outputFiles && results.outputFiles.length > 0) {
        let resultHtml = `<p>Computation Finished. `;
        resultHtml += `Created ${results.outputFiles.length} files. `;
        if (results.plotFiles && results.plotFiles.length > 0) {
            resultHtml += `Created ${results.plotFiles.length} charts on web page.`;
        }
        resultHtml += `</p>`;

        // Add restore info
        const restoreTime = new Date(results.timestamp).toLocaleString('zh-CN');
        resultHtml += `<p class="restore-info">Previous results restored: (${restoreTime})</p>`;

        outputDiv.innerHTML = resultHtml;

        // Initialize plotter and restore charts
        if (window.reaxToolsPlotter) {
            window.reaxToolsPlotter.init('charts-area');

            // Try to restore charts from saved data first
            if (!window.reaxToolsPlotter.restoreCharts('charts-area')) {
                // If no saved chart data, recreate charts from plot files
                if (results.plotFiles && results.plotFiles.length > 0) {
                    createChartsFromFiles(results.plotFiles);
                }
            }
        }

        return true;
    }

    return false;
}

//
//
//
//

document.getElementById('fileInput').addEventListener('change', (e) => {
    uploadedFile = e.target.files[0]; //
    if (!uploadedFile) return;

    // Save file info to localStorage
    saveUploadedFileInfo(uploadedFile);

    document.getElementById('output').innerHTML =
        `Selected file: ${uploadedFile.name} (${(uploadedFile.size / 1024).toFixed(2)} KB)`;
});

//
document.getElementById('loadSampleBtn').addEventListener('click', async (e) => {
    e.preventDefault(); //
    try {
        const response = await fetch('./reaxff.xyz');
        if (!response.ok) {
            throw new Error(`Can not load example file: ${response.statusText}`);
        }
        const fileContent = await response.blob();
        uploadedFile = new File([fileContent], "reaxff.xyz", { type: "text/plain" });

        if (!uploadedFile) return;

        // Save file info to localStorage
        saveUploadedFileInfo(uploadedFile);

        document.getElementById('output').innerHTML =
            `Seleted file: ${uploadedFile.name} (${(uploadedFile.size / 1024).toFixed(2)} KB)`;

        //
        document.getElementById('fileInput').value = '';

    } catch (error) {
        alert(`Failed to load file: ${error.message}`);
        console.error(error);
    }
});

// Page load event handler to restore data
document.addEventListener('DOMContentLoaded', async function () {
    // Try to restore calculation results first
    const restored = await restoreCalculationResults();

    if (restored) {
        console.log('Calculation results restored from localStorage');
    } else {
        // If no calculation results, try to restore file info and args
        const fileInfo = loadUploadedFileInfo();
        const userArgs = loadUserArgs();

        if (fileInfo) {
            document.getElementById('output').innerHTML =
                `Seleted file: ${fileInfo.name} (${(fileInfo.size / 1024).toFixed(2)} KB)`;
        }

        if (userArgs) {
            document.getElementById('argsInput').value = userArgs;
        }
    }

    // Add chart control button event listeners
    //
    document.getElementById('clearBtn').addEventListener('click', function () {
        // Clear charts first
        if (window.reaxToolsPlotter) {
            window.reaxToolsPlotter.clearCharts();
        }
        // Clear all stored data including virtual file system
        clearStoredData();
        // Clear output display and UI elements
        document.getElementById('output').innerHTML = '';
        document.getElementById('log-output').textContent = '';
        document.getElementById('argsInput').value = '';
        document.getElementById('fileInput').value = '';
        uploadedFile = null; // Reset uploaded file reference

        // Refresh the browser to ensure complete reset
        window.location.reload();
    });

    //
    document.getElementById('exportAllChartsBtn').addEventListener('click', function () {
        if (window.reaxToolsPlotter && window.reaxToolsPlotter.charts) {
            window.reaxToolsPlotter.charts.forEach((_, chartId) => {
                window.reaxToolsPlotter.exportChart(chartId, `${chartId}.png`);
            });
        }
    });

});

//
//
function loadWasmModule() {
    return new Promise((resolve, reject) => {
        //
        if (wasmModule) {
            resolve(wasmModule);
            return;
        }

        const logOutput = document.getElementById('log-output');
        //
        window.Module = {
            print: text => {
                console.log(`[WASM]: ${text}`);
                if (logOutput) {
                    //
                    if (text.includes('\r')) {
                        //
                        const parts = text.split('\r');
                        //
                        const last = parts[parts.length - 1];
                        //
                        let lines = logOutput.textContent.split('\n');
                        //
                        if (lines.length > 0 && lines[lines.length - 1] === '') lines.pop();
                        //
                        lines[lines.length - 1] = last;
                        logOutput.textContent = lines.join('\n') + '\n';
                    } else {
                        logOutput.textContent += text + '\n';
                    }
                    logOutput.scrollTop = logOutput.scrollHeight; //
                }
            },
            printErr: text => {
                console.error(`[WASM Error]: ${text}`);
                if (logOutput) {
                    logOutput.textContent += `[ERROR] ${text}\n`;
                    logOutput.scrollTop = logOutput.scrollHeight; //
                }
            },
            onRuntimeInitialized: function () {
                wasmModule = this; //
                //
                if (!wasmModule.FS.analyzePath('/output').exists) {
                    wasmModule.FS.mkdir('/output');
                }
                resolve(wasmModule); //
            },
            //
            locateFile: (path) => {
                if (path.endsWith('.wasm')) return './wasm_main.wasm';
                return path;
            }
        };
        //
        const script = document.createElement('script');
        script.src = './wasm_main.js';
        script.onerror = reject;
        document.head.appendChild(script);
    });
}

//
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

//
//
//

document.getElementById('runBtn').onclick = async function () {
    if (!uploadedFile) {
        alert('Please choose input file');
        return;
    }
    const outputDiv = document.getElementById('output');
    const logOutput = document.getElementById('log-output');
    const userArgsText = document.getElementById('argsInput').value.trim();

    outputDiv.innerHTML = "Computing...";
    logOutput.textContent = ''; // Clear previous output

    // Save user arguments
    saveUserArgs(userArgsText);

    // Clear previous charts
    if (window.reaxToolsPlotter) {
        window.reaxToolsPlotter.clearCharts();
    }

    try {
        //
        wasmModule = await loadWasmModule();

        //
        if (wasmModule.FS.analyzePath('/output').exists) {
            const oldFiles = wasmModule.FS.readdir('/output');
            oldFiles.forEach(file => {
                if (file !== '.' && file !== '..') {
                    wasmModule.FS.unlink(`/output/${file}`);
                }
            });
        }

        //
        const fileContent = new Uint8Array(await uploadedFile.arrayBuffer());
        const inputFilePath = uploadedFile.name; //
        wasmModule.FS.writeFile(inputFilePath, fileContent);

        //
        const userArgs = userArgsText ? userArgsText.split(/\s+/).filter(Boolean) : [];
        const args = ['program', '--traj', inputFilePath, ...userArgs]; //

        //
        const argvPtrs = args.map(s => {
            const len = wasmModule.lengthBytesUTF8(s) + 1; //
            const ptr = wasmModule._malloc(len); //
            wasmModule.stringToUTF8(s, ptr, len); //
            return ptr;
        });
        //
        const argvPtr = wasmModule._malloc(argvPtrs.length * 4);
        //
        argvPtrs.forEach((ptr, i) => {
            wasmModule.setValue(argvPtr + i * 4, ptr, 'i32');
        });

        //
        //
        const ret = wasmModule.ccall(
            'cpp_main', //
            'number',   //
            ['number', 'number'], //
            [args.length, argvPtr] //
        );

        //
        argvPtrs.forEach(ptr => wasmModule._free(ptr));
        wasmModule._free(argvPtr);

        //
        const outputFiles = wasmModule.FS.readdir('/output').filter(p => p !== '.' && p !== '..');

        //
        const logOutputText = logOutput.textContent || '';

        //
        function getBeijingTimeISO() {
            const now = new Date();
            const beijingOffset = 8 * 60; //
            const localOffset = now.getTimezoneOffset();
            const beijingTime = new Date(now.getTime() + (beijingOffset + localOffset) * 60000);
            return beijingTime.toISOString().replace('T', ' ').replace('Z', '');
        }
        const beijingTime = getBeijingTimeISO();

        //
        await fetch('/api/log', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ time: beijingTime, logOutput: logOutputText })
        });

        if (outputFiles.length === 0) {
            outputDiv.innerHTML = `<p>Computation terminated, but not data output.</p>`;
        } else {
            //
            const plotFiles = [];
            const allFiles = [];

            //
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
                //
                allFiles.push(fileName);
            });

            console.log(`Found ${plotFiles.length} files for visualization:`, plotFiles.map(f => f.name));

            //
            let downloadLinks = '';
            if (allFiles.length === 1) {
                //
                const fileName = allFiles[0];
                const resultData = wasmModule.FS.readFile(`/output/${fileName}`);
                const blob = new Blob([resultData], { type: 'application/octet-stream' });
                const url = URL.createObjectURL(blob);
                downloadLinks = `<a href="${url}" download="${fileName}">Save result file (${fileName})</a>`;
            } else {
                //
                const zip = new JSZip();
                allFiles.forEach(fileName => {
                    const fileData = wasmModule.FS.readFile(`/output/${fileName}`);
                    zip.file(fileName, fileData);
                });

                const zipBlob = await zip.generateAsync({ type: 'blob' });
                const url = URL.createObjectURL(zipBlob);
                const dateStr = new Date().toISOString().replace(/:/g, '-').split('.')[0];
                downloadLinks = `<a href="${url}" download="output-${dateStr}.zip">Save all result files (.zip)</a>`;
            }

            //
            let resultHtml = `<p>Computation finished. `;
            resultHtml += `Created ${allFiles.length} files. `;
            if (plotFiles.length > 0) {
                resultHtml += `Created ${plotFiles.length} charts on web page.`;
            }
            resultHtml += `</p>`;

            if (downloadLinks) {
                resultHtml += downloadLinks;
            }

            outputDiv.innerHTML = resultHtml;

            // Save calculation results to localStorage
            saveCalculationResults(allFiles, plotFiles, logOutputText, userArgsText);

            //
            if (plotFiles.length > 0) {
                createChartsFromFiles(plotFiles);
            }
        }

    } catch (error) {
        outputDiv.innerHTML = `<p class="error-message">Error: ${error.toString()}</p>`;
        console.error(error);
    }
};
