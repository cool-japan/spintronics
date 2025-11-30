#!/bin/bash
# Simple HTTP server for WASM demo

PORT=${1:-8080}

echo "🌐 Starting HTTP server on port $PORT..."
echo "📱 Open http://localhost:$PORT in your browser"
echo ""
echo "Press Ctrl+C to stop"
echo ""

if command -v python3 &> /dev/null; then
    python3 -m http.server $PORT
elif command -v python &> /dev/null; then
    python -m SimpleHTTPServer $PORT
else
    echo "❌ Python not found!"
    echo "Install Python or use another HTTP server"
    exit 1
fi
