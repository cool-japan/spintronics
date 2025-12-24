#!/bin/bash
# Run the Spintronics Interactive Demo Server

set -e

echo "🚀 Starting Spintronics Interactive Demos..."
echo "📍 Server will be available at http://localhost:3000"
echo ""

# Build and run in release mode for better performance
cargo run --release

# Alternative development mode with auto-reload:
# cargo watch -x 'run -p spintronics-demo'
