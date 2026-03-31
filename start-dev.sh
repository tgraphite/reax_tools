#!/bin/bash

# ReaxTools Web Development Startup Script

echo "==============================================="
echo "  ReaxTools Web 2.0 - Development Server"
echo "==============================================="
echo ""

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if Redis is running
if ! redis-cli ping > /dev/null 2>&1; then
    echo -e "${YELLOW}⚠ Redis is not running. Starting Redis...${NC}"

    # Try to start Redis with Docker
    if command -v docker &> /dev/null; then
        docker run -d --name reax-tools-redis -p 6379:6379 redis:7-alpine > /dev/null 2>&1 || true
        sleep 2
    fi

    # Check again
    if ! redis-cli ping > /dev/null 2>&1; then
        echo -e "${RED}✗ Failed to start Redis. Please start Redis manually.${NC}"
        echo "  Docker: docker run -d -p 6379:6379 redis:7-alpine"
        echo "  Or install Redis locally"
        exit 1
    fi
fi

echo -e "${GREEN}✓ Redis is running${NC}"

# Check if binary exists
if [ ! -f "bin/reax_tools" ]; then
    echo -e "${RED}✗ reax_tools binary not found at bin/reax_tools${NC}"
    echo "  Please build the binary first: bash build_local.sh"
    exit 1
fi

echo -e "${GREEN}✓ reax_tools binary found${NC}"

# Install server dependencies if needed
if [ ! -d "server/node_modules" ]; then
    echo -e "${YELLOW}⚠ Server dependencies not found. Installing...${NC}"
    cd server && npm install && cd ..
fi

# Install frontend dependencies if needed
if [ ! -d "web-new/node_modules" ]; then
    echo -e "${YELLOW}⚠ Frontend dependencies not found. Installing...${NC}"
    cd web-new && npm install && cd ..
fi

# Build frontend CSS if needed
if [ ! -f "web-new/dist/styles.css" ]; then
    echo -e "${YELLOW}⚠ Frontend CSS not built. Building...${NC}"
    cd web-new && npm run build && cd ..
fi

echo ""
echo "==============================================="
echo "  Starting services..."
echo "==============================================="
echo ""

# Create temp directories
mkdir -p server/temp/uploads server/temp/outputs

# Start server
echo -e "${GREEN}Starting API server on http://localhost:3000${NC}"
cd server && npm start &
SERVER_PID=$!
cd ..

# Wait a moment for server to start
sleep 2

# Check if server started
if ! kill -0 $SERVER_PID 2>/dev/null; then
    echo -e "${RED}✗ Server failed to start${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Server started (PID: $SERVER_PID)${NC}"
echo ""

echo "==============================================="
echo "  ReaxTools Web is ready!"
echo "==============================================="
echo ""
echo "  Frontend:  file:///$(pwd)/web-new/index.html"
echo "  API:       http://localhost:3000/api"
echo "  Health:    http://localhost:3000/api/health"
echo ""
echo "Press Ctrl+C to stop all services"
echo ""

# Handle shutdown
cleanup() {
    echo ""
    echo "Shutting down..."
    kill $SERVER_PID 2>/dev/null
    exit 0
}

trap cleanup INT TERM

# Wait
wait $SERVER_PID
