# Spintronics WASM Demo

Interactive browser-based spintronics simulations powered by Rust and WebAssembly.

**npm package**: [`@cooljapan/spintronics`](https://www.npmjs.com/package/@cooljapan/spintronics)

```bash
npm install @cooljapan/spintronics
```

## Features

### 🎯 Single Spin Dynamics
- Real-time LLG equation solver
- Interactive magnetic field control
- Adjustable damping parameter
- 3D visualization with x-y projection and z-component indicator

### 🌊 Magnon Propagation
- 1D spin chain simulation
- Wave excitation and propagation
- Configurable chain length and coupling
- Color-coded visualization

### ⚡ Spin Hall Effect Calculator
- Charge-to-spin current conversion
- Adjustable spin Hall angle
- Instant calculation for device design

## Building the Demo

### Prerequisites

1. Install Rust (if not already installed):
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2. Install wasm-pack:
```bash
cargo install wasm-pack
```

### Build Steps

1. Navigate to the project root:
```bash
cd /path/to/spintronics
```

2. Build the WASM package:
```bash
wasm-pack build --features wasm --target web
```

This will create a `pkg/` directory containing:
- `spintronics_bg.wasm` - The compiled WebAssembly binary
- `spintronics.js` - JavaScript bindings
- `spintronics.d.ts` - TypeScript definitions

3. Copy the pkg directory to wasm-demo:
```bash
cp -r pkg wasm-demo/
```

4. Serve the demo with a local HTTP server:
```bash
cd wasm-demo
python3 -m http.server 8080
```

Or using Node.js:
```bash
npx http-server -p 8080
```

5. Open your browser and navigate to:
```
http://localhost:8080
```

## Quick Build Script

We provide a convenience script to automate the build process:

```bash
./build-wasm.sh
```

## Performance Notes

- The simulation runs entirely client-side (no server required)
- Physics calculations are compiled to WebAssembly for near-native speed
- Typical performance: 1000+ simulation steps per second
- Best experienced on modern browsers (Chrome 90+, Firefox 88+, Safari 15+, Edge 90+)

## Physics Background

### LLG Equation
The Landau-Lifshitz-Gilbert equation governs magnetization dynamics:

```
dm/dt = -γ/(1+α²) [m × H_eff + α m × (m × H_eff)]
```

where:
- `m` - normalized magnetization direction
- `γ` - gyromagnetic ratio (2.21×10⁵ m/(A·s) for electrons)
- `α` - Gilbert damping constant
- `H_eff` - effective magnetic field

### Spin Hall Effect
Converts charge current into transverse spin current:

```
j_s = θ_SH × j_c
```

where:
- `j_s` - spin current density
- `j_c` - charge current density
- `θ_SH` - spin Hall angle (material-dependent, ~0.07 for Pt)

### Magnon Propagation
Spin waves in magnetic materials described by coupled LLG equations
with exchange interaction between neighboring spins.

## Troubleshooting

### WASM module fails to load
- Ensure you're serving the page via HTTP(S), not file://
- Check browser console for detailed error messages
- Verify pkg/ directory exists and contains .wasm file

### Simulation is slow
- Close other browser tabs
- Disable browser extensions
- Try a different browser (Chrome generally has best WASM performance)

### Canvas not rendering
- Check if JavaScript is enabled
- Verify browser supports Canvas API
- Check browser console for errors

## Technical Details

### WASM Module Size
- Uncompressed: ~2-3 MB
- Gzip compressed: ~500-700 KB
- Typical load time: <1 second on broadband

### Browser Compatibility
- Chrome/Edge 90+ ✅
- Firefox 88+ ✅
- Safari 15+ ✅
- Opera 76+ ✅

### Limitations
- No multithreading (WASM threads are experimental)
- No file I/O (browser security sandbox)
- Limited to client-side resources

## License

Apache-2.0

## References

Based on research by Prof. Eiji Saitoh's group:
- E. Saitoh et al., Appl. Phys. Lett. 88, 182509 (2006)
- K. Uchida et al., Nature 455, 778-781 (2008)
