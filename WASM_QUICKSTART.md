# WebAssembly Quick Start Guide

Get your spintronics simulations running in the browser in under 5 minutes!

## Prerequisites

```bash
# Install Rust (if not already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install wasm-pack
cargo install wasm-pack
```

## Build & Run

### Option 1: Use the build script (recommended)

```bash
# Build WASM package
./build-wasm.sh

# Start the demo server
cd wasm-demo
./serve.sh

# Open browser to http://localhost:8080
```

### Option 2: Manual build

```bash
# Build WASM
wasm-pack build --features wasm --no-default-features

# Copy to demo directory
cp -r pkg wasm-demo/

# Serve
cd wasm-demo
python3 -m http.server 8080
```

## What You'll See

The demo includes three interactive simulations:

### 1. 🎯 Single Spin Precession
Watch a magnetic moment precess around an applied field in real-time.
- Adjust field components (Hx, Hy, Hz)
- Control damping parameter (α)
- See live magnetization components

### 2. 🌊 Magnon Propagation
Observe spin waves traveling through a 1D spin chain.
- Excite waves at any position
- Adjust amplitude and chain length
- Watch magnons propagate

### 3. ⚡ Spin Hall Effect
Calculate spin current generation from charge current.
- Set charge current density
- Adjust spin Hall angle
- Instant results

## Performance

- **WASM size**: 62 KB (optimized)
- **Load time**: <1 second
- **Frame rate**: 60 FPS
- **Simulation rate**: 1000+ steps/second

## Browser Requirements

- Chrome 90+
- Firefox 88+
- Safari 15+
- Edge 90+

## Using in Your Own Project

### JavaScript Example

```javascript
import init, { SpinSimulator } from './pkg/spintronics.js';

async function runSimulation() {
    // Initialize WASM module
    await init();

    // Create simulator
    const sim = new SpinSimulator();

    // Set magnetic field (A/m)
    sim.set_field(1000, 0, 10000);

    // Run simulation loop
    for (let i = 0; i < 1000; i++) {
        sim.step(0.01);  // 0.01 ns time step

        // Get magnetization
        const mx = sim.get_mx();
        const my = sim.get_my();
        const mz = sim.get_mz();

        console.log(`Step ${i}: m = (${mx}, ${my}, ${mz})`);
    }
}

runSimulation();
```

### TypeScript Example

```typescript
import init, { SpinSimulator, SpinChain } from './pkg/spintronics';

async function magnonDemo() {
    await init();

    const chain = new SpinChain(50, 1e-21);
    chain.excite_wave(25, 0.5, 3.0);

    setInterval(() => {
        chain.step(0.01);

        // Export to JSON for visualization
        const data = chain.export_json();
        updateVisualization(JSON.parse(data));
    }, 16); // ~60 FPS
}
```

## API Reference

### SpinSimulator

```javascript
// Constructor
const sim = new SpinSimulator();
const sim = SpinSimulator.with_material(ms, alpha, gamma);

// Methods
sim.set_field(hx, hy, hz);           // Set external field (A/m)
sim.set_magnetization(mx, my, mz);   // Set initial magnetization
sim.step(dt_ns);                     // Advance by dt nanoseconds
sim.run_steps(n, dt_ns);             // Run multiple steps
sim.reset();                          // Reset to initial state

// Getters
sim.get_mx();                         // x-component (normalized)
sim.get_my();                         // y-component (normalized)
sim.get_mz();                         // z-component (normalized)
sim.get_time();                       // Current time (seconds)
sim.get_magnetization();              // Vec3 object
sim.get_effective_field();            // Effective field (A/m)
```

### SpinChain

```javascript
// Constructor
const chain = new SpinChain(n_spins, exchange_j);

// Methods
chain.set_spin(i, mx, my, mz);       // Set spin at index
chain.excite_wave(pos, amp, width);  // Excite spin wave
chain.step(dt_ns);                   // Time evolution
chain.export_json();                  // Export all spins

// Getters
chain.get_mx(i);                      // x-component at index i
chain.get_my(i);                      // y-component at index i
chain.get_mz(i);                      // z-component at index i
chain.length();                       // Number of spins
chain.get_time();                     // Current time
```

### SpinHallCalculator

```javascript
// Constructor
const calc = new SpinHallCalculator();          // Platinum
const calc = SpinHallCalculator.with_angle(θ);  // Custom

// Methods
calc.spin_current(jc, theta, phi);    // Calculate js (A/m²)
calc.get_angle();                      // Get θ_SH
```

## Troubleshooting

### WASM fails to load
- Serve via HTTP(S), not `file://`
- Check browser console for errors
- Ensure pkg/ directory exists

### Simulation is slow
- Close other tabs
- Use Chrome for best performance
- Reduce time step or chain length

### Build fails
- Ensure wasm32 target: `rustup target add wasm32-unknown-unknown`
- Update wasm-pack: `cargo install --force wasm-pack`
- Check Rust version: `rustc --version` (need 1.70+)

## Next Steps

- Customize the HTML/CSS in `wasm-demo/index.html`
- Add more visualization (WebGL, Canvas2D)
- Integrate with React/Vue/Svelte
- Build a full simulation dashboard

## Resources

- [WASM Bindgen Book](https://rustwasm.github.io/wasm-bindgen/)
- [WebAssembly.org](https://webassembly.org/)
- [Spintronics Library Docs](https://docs.rs/spintronics)

Happy simulating! 🧲✨
