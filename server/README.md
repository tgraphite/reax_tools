# ReaxTools Web 2.0

反应分子动力学分析平台 - 基于服务器的队列化任务处理系统

## 架构

```
┌─────────────┐      ┌─────────────┐      ┌─────────────┐
│   前端      │──────▶│   API 服务  │──────▶│  任务队列   │
│ (Tailwind)  │◀──────│  (Express)  │◀──────│  (BullMQ)   │
└─────────────┘      └─────────────┘      └──────┬──────┘
                                                  │
                                                  ▼
                                           ┌─────────────┐
                                           │  Worker进程 │
                                           │ (C++二进制) │
                                           └─────────────┘
```

## 快速开始

### 方式一：Docker Compose (推荐)

```bash
docker-compose up -d
```

访问 http://localhost

### 方式二：本地开发

**前提条件:**
- Node.js 18+
- Redis 7+
- 编译好的 reax_tools 二进制文件

**步骤:**

1. 启动 Redis:
```bash
docker run -d -p 6379:6379 redis:7-alpine
```

2. 启动开发服务器:
```bash
bash start-dev.sh
```

或者直接:
```bash
cd server && npm install && npm start
cd web-new && npm install && npm run build
# 然后打开 web-new/index.html
```

## API 接口

| 方法 | 路径 | 说明 |
|------|------|------|
| POST | /api/jobs | 提交新任务 |
| GET | /api/jobs/:id | 查询任务状态 |
| GET | /api/jobs/:id/stream | SSE 实时日志流 |
| GET | /api/jobs/:id/download | 下载结果 ZIP |
| GET | /api/jobs | 任务列表 |
| DELETE | /api/jobs/:id | 删除任务 |

## 特性

- **队列系统**: 支持多任务排队，并发控制
- **实时日志**: Server-Sent Events 实时推送分析进度
- **数据大屏**: 多面板布局，物种/键/环统计图表
- **反应网络**: Sigma.js 可视化物质转化关系
- **暗色模式**: 支持深色主题
- **响应式设计**: 适配桌面和移动设备

## 配置

编辑 `server/.env`:

```
PORT=3000
REDIS_HOST=localhost
REDIS_PORT=6379
MAX_CONCURRENT_JOBS=2
JOB_TIMEOUT=300000
```

## 文件结构

```
reax_tools/
├── server/          # Node.js 后端
│   ├── src/
│   │   ├── index.js         # 服务器入口
│   │   ├── routes/jobs.js   # API 路由
│   │   └── services/        # 队列服务
│   └── .env
├── web-new/         # 新前端 (Tailwind)
│   ├── index.html
│   ├── src/
│   │   ├── main.js          # 主逻辑
│   │   ├── api.js           # API 客户端
│   │   └── charts.js        # 图表组件
│   └── dist/styles.css      # 编译后的 CSS
└── docker-compose.yml
```

## 从旧版本迁移

旧版本使用 WASM 在浏览器中运行，新版本改为服务器端执行：

| 特性 | 旧版 (WASM) | 新版 (Server) |
|------|-------------|---------------|
| 计算位置 | 浏览器本地 | 服务器 |
| 支持文件大小 | <100MB | <500MB |
| 并发处理 | 单任务 | 多任务队列 |
| 离线使用 | 支持 | 不支持 |
| 启动速度 | 慢 (需加载WASM) | 快 |

## 许可证

同主项目
