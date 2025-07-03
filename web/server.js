// Simple Node.js server for logging time and IP
const express = require('express');
const fs = require('fs');
const path = require('path');
const axios = require('axios'); // For IP geolocation
const app = express();
const PORT = 3000;

app.set('trust proxy', true); // 信任代理，才能正确获取 X-Forwarded-For

app.use(express.json());

app.post('/api/log', async (req, res) => {
    // Get Beijing time (UTC+8)
    function getBeijingTime() {
        const now = new Date();
        return now.toLocaleString('zh-CN', { hour12: false, timeZone: 'Asia/Shanghai' });
    }
    const time = getBeijingTime();
    const logOutput = req.body.logOutput || '';
    // Get client IP
    let ip = req.ip;
    let location = '';
    try {
        // Use ip-api.com for geolocation
        const geoRes = await axios.get(`http://ip-api.com/json/${ip}?lang=en`);
        if (geoRes.data && geoRes.data.status === 'success') {
            location = `${geoRes.data.country || ''} ${geoRes.data.regionName || ''} ${geoRes.data.city || ''}`.trim();
        }
    } catch (e) {
        // Ignore geolocation errors
    }
    // Special separator for log analysis
    const separator = '###LOG###';
    const logLine = `${separator} ${time} - ${ip}${location ? ' (' + location + ')' : ''}\n${logOutput}\n---\n`;
    const logPath = path.join(__dirname, 'access.log');
    fs.appendFile(logPath, logLine, err => {
        if (err) {
            console.error('Failed to write log:', err);
            return res.status(500).send('Log error');
        }
        res.send('Logged');
    });
});

app.listen(PORT, () => {
    console.log(`Log server running at http://localhost:${PORT}`);
});
