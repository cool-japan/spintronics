#!/bin/bash
# Build script for Spintronics WASM
# Publishes as @cooljapan/spintronics on npm

set -e

SCOPE="cooljapan"
NPM_NAME="@${SCOPE}/spintronics"

echo "Building Spintronics WASM (${NPM_NAME})..."
echo ""

# Check if wasm-pack is installed
if ! command -v wasm-pack &> /dev/null; then
    echo "Error: wasm-pack is not installed"
    echo "Install it with: cargo install wasm-pack"
    exit 1
fi

# Build for web target
echo "Building for web target..."
wasm-pack build --scope "${SCOPE}" --target web --out-dir pkg-web --release -- --features wasm --no-default-features

# Build for Node.js target
echo "Building for Node.js target..."
wasm-pack build --scope "${SCOPE}" --target nodejs --out-dir pkg-node --release -- --features wasm --no-default-features

# Build for bundler target (webpack, etc.)
echo "Building for bundler target..."
wasm-pack build --scope "${SCOPE}" --target bundler --out-dir pkg-bundler --release -- --features wasm --no-default-features

# Rename package to @cooljapan/spintronics in all generated package.json files
for dir in pkg-web pkg-node pkg-bundler; do
    if [ -f "${dir}/package.json" ]; then
        sed -i '' "s|\"@${SCOPE}/spintronics\"|\"${NPM_NAME}\"|g" "${dir}/package.json"
        echo "  Updated ${dir}/package.json -> ${NPM_NAME}"
    fi
done

# Copy web build to demo directory for local testing
echo ""
echo "Copying web build to wasm-demo/..."
rm -rf wasm-demo/pkg
cp -r pkg-web wasm-demo/pkg

# Get package size
WASM_SIZE=$(du -h pkg-bundler/spintronics_bg.wasm | cut -f1)
echo ""
echo "Build complete!"
echo "  WASM size: ${WASM_SIZE}"
echo ""
echo "Output directories:"
echo "  - pkg-web/      : For use in web browsers"
echo "  - pkg-node/     : For use in Node.js"
echo "  - pkg-bundler/  : For use with webpack/rollup/etc."
echo ""
echo "To publish (bundler target recommended for npm):"
echo "  cd pkg-bundler && npm publish --access public"
echo ""
echo "To run the demo:"
echo "  cd wasm-demo"
echo "  python3 -m http.server 8080"
echo "  Then open: http://localhost:8080"
