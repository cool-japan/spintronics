# Spintronics Interactive Web Demonstrations

**HTMX + Rust** powered real-time physics simulations

## Overview

This subcrate provides interactive web demonstrations of spintronics physics phenomena using a modern, server-side rendering approach:

- **Backend**: Axum (Rust web framework)
- **Frontend**: HTMX (hypermedia-driven interactions)
- **Templates**: Askama (type-safe Rust templates)
- **Physics**: Full `spintronics` library (server-side calculations)

## Features

### 🌀 LLG Magnetization Dynamics
- Real-time Landau-Lifshitz-Gilbert equation solver
- Interactive parameter controls (damping, applied field, initial state)
- Trajectory visualization on Bloch sphere
- Single-step and continuous animation modes

### ⚡ Spin Pumping Calculator
- Reproduces Saitoh 2006 APL experiment
- YIG/Pt, Py/Pt, CoFeB/Pt material systems
- Calculates spin current and ISHE voltage
- Automatic parameter updates via HTMX

### 🔬 Materials Explorer
- Comprehensive ferromagnet database
- Compare M_s, α, A_ex, K across materials
- Detailed explanations for YIG, Permalloy, CoFeB, Fe, Co, Ni

### 🎯 Magnetic Skyrmion Visualizer
- Real-time magnetization field rendering
- Tunable radius, helicity (Néel/Bloch), chirality (CW/CCW)
- Topological charge calculation
- Educational content on DMI and skyrmion physics

## Quick Start

### Prerequisites
- Rust 1.70.0 or later
- Tokio async runtime

### Run the Server

```bash
cd demo
cargo run --release
```

Then open http://localhost:3000 in your browser.

### Development Mode

```bash
# With hot-reload (requires cargo-watch)
cargo watch -x 'run'
```

## Architecture

### HTMX-Powered Interactivity

The demo uses HTMX for seamless server-client communication without JavaScript frameworks:

```html
<!-- Example: LLG simulation trigger -->
<button hx-get="/api/llg/simulate"
        hx-include="#llg-form"
        hx-target="#llg-result"
        hx-trigger="click">
    Run Simulation
</button>
```

When clicked:
1. Form parameters sent to `/api/llg/simulate`
2. Server runs physics calculation (Rust)
3. Server returns HTML fragment
4. HTMX swaps content into `#llg-result`
5. **No page reload, no JSON parsing, no JavaScript**

### Server-Side Rendering Benefits

- **Full Library Access**: Use all `spintronics` features (FEM, HDF5, advanced solvers)
- **No WASM Limitations**: No need to compile to WebAssembly
- **Type Safety**: Templates checked at compile time (Askama)
- **Simple Deploy**: Single binary, no build steps
- **Progressive Enhancement**: Works without JavaScript enabled

### Project Structure

```
demo/
├── Cargo.toml              # Dependencies (axum, htmx, askama)
├── src/
│   ├── main.rs            # Server setup, routes
│   ├── demos.rs           # Physics calculations & HTML generation
│   └── templates.rs       # Askama template definitions
├── templates/             # HTML templates
│   ├── base.html          # Base layout with nav/footer
│   ├── index.html         # Landing page
│   ├── llg.html           # LLG dynamics demo
│   ├── spin_pumping.html  # Spin pumping calculator
│   ├── materials.html     # Materials explorer
│   └── skyrmion.html      # Skyrmion visualizer
└── static/
    └── css/
        └── style.css      # Modern, physics-inspired styling
```

## API Endpoints

### GET `/`
Home page with demo cards

### GET `/llg`
LLG dynamics interactive demo

### POST `/api/llg/simulate`
Parameters: `mx, my, mz, hx, hy, hz, alpha, dt_ps, steps`
Returns: HTML with trajectory SVG

### POST `/api/llg/step`
Parameters: Same as simulate
Returns: HTML with current magnetization state

### GET `/api/spin-pumping/calculate`
Parameters: `material, frequency_ghz, h_rf`
Returns: HTML with spin current and voltage results

### GET `/api/materials/compare`
Returns: HTML table comparing all materials

### POST `/api/skyrmion/visualize`
Parameters: `radius_nm, helicity, chirality`
Returns: HTML with SVG magnetization field

## Customization

### Adding New Demos

1. **Create template** in `templates/your_demo.html`:
```html
{% extends "base.html" %}
{% block content %}
<!-- Your interactive content -->
{% endblock %}
```

2. **Add route handler** in `main.rs`:
```rust
.route("/your-demo", get(your_demo_page))
.route("/api/your-demo/calculate", post(your_demo_calculate))
```

3. **Implement physics** in `demos.rs`:
```rust
pub async fn your_demo_calculate(Query(params): Query<YourRequest>) -> impl IntoResponse {
    // Use spintronics library
    let result = calculate_physics(params);
    // Return HTML fragment
    Html(generate_html(&result))
}
```

### Styling

Edit `static/css/style.css` to customize:
- Color scheme (CSS variables in `:root`)
- Layout (grid, flexbox)
- Animations and transitions
- Responsive breakpoints

## Performance

- **Server-side calculations**: Leverage full Rust optimization
- **Minimal client load**: Only HTML transferred (no heavy JS bundles)
- **Efficient updates**: HTMX sends only changed fragments
- **Concurrent requests**: Tokio async runtime handles parallel demos

## Deployment

### Local Binary
```bash
cargo build --release
./target/release/spintronics-demo
```

### Docker
```dockerfile
FROM rust:1.70 as builder
WORKDIR /app
COPY . .
RUN cargo build --release -p spintronics-demo

FROM debian:bookworm-slim
COPY --from=builder /app/target/release/spintronics-demo /usr/local/bin/
COPY --from=builder /app/demo/templates /templates
COPY --from=builder /app/demo/static /static
CMD ["spintronics-demo"]
```

### Cloud Deploy
- **Fly.io**: `fly launch` (zero-config Rust support)
- **Railway**: Connect GitHub, auto-deploy
- **Shuttle.rs**: Native Rust PaaS

## Educational Use

These demos are designed for:
- **Physics courses**: Visualize abstract concepts (precession, damping, topology)
- **Research onboarding**: New group members learn spintronics interactively
- **Conference demos**: Live physics at poster sessions
- **Self-study**: Progressive difficulty (basic → advanced)

## References

The physics implementations are validated against:

1. **E. Saitoh et al., APL 88, 182509 (2006)** - Spin pumping & ISHE
2. **K. Uchida et al., Nature 455, 778 (2008)** - Spin Seebeck effect
3. **I. M. Miron et al., Nature 476, 189 (2011)** - Spin-orbit torque
4. **S. Woo et al., Nat. Mater. 15, 501 (2016)** - Room-temperature skyrmions

## License

MIT OR Apache-2.0

## Contributing

Contributions welcome! Ideas for new demos:
- Antiferromagnetic THz dynamics
- Magnonic crystal band structures
- Topological Hall effect mapping
- Reservoir computing with spin waves

## Support

For issues or questions:
- GitHub: https://github.com/cool-japan/spintronics/issues
- Documentation: https://docs.rs/spintronics

---

**Built with ❤️ by COOLJAPAN OÜ (Team KitaSan)**
