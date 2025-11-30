#!/bin/bash
# Build script for spintronics WASM demo

set -e

echo "🦀 Building Spintronics WASM package..."
echo ""

# Check if wasm-pack is installed
if ! command -v wasm-pack &> /dev/null; then
    echo "❌ wasm-pack not found!"
    echo "Install it with: cargo install wasm-pack"
    exit 1
fi

# Build the WASM package
echo "📦 Compiling Rust to WebAssembly..."
wasm-pack build --features wasm --no-default-features

# Copy to demo directory
echo "📋 Copying to wasm-demo/..."
rm -rf wasm-demo/pkg
cp -r pkg wasm-demo/

# Get package size
WASM_SIZE=$(du -h wasm-demo/pkg/spintronics_bg.wasm | cut -f1)
echo ""
echo "✅ Build complete!"
echo "   WASM size: $WASM_SIZE"
echo ""
echo "🚀 To run the demo:"
echo "   cd wasm-demo"
echo "   python3 -m http.server 8080"
echo "   "
echo "   Then open: http://localhost:8080"
echo ""
