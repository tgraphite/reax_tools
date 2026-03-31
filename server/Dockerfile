# Production Dockerfile for ReaxTools Server
FROM node:20-alpine

WORKDIR /app

# Install Python3, matplotlib, pandas, numpy, and Graphviz for image generation
RUN apk add --no-cache python3 py3-matplotlib py3-pandas py3-numpy graphviz coreutils

# Install dependencies
COPY package*.json ./
RUN npm ci --only=production

# Copy source code
COPY . .

# Create temp directories
RUN mkdir -p temp/uploads temp/outputs

# Expose port
EXPOSE 3000

# Health check
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
  CMD node -e "require('http').get('http://localhost:3000/api/health', (r) => r.statusCode === 200 ? process.exit(0) : process.exit(1))"

# Start server
CMD ["node", "src/index.js"]
