# 🌀 spintronics

A pure Rust library for simulating spin dynamics, spin current generation, and conversion phenomena in magnetic and topological materials.

**Inspired by the pioneering work of Prof. Eiji Saitoh's Group (University of Tokyo / RIKEN CEMS)**

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)]()
[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue)]()
[![Rust Version](https://img.shields.io/badge/rust-2021-orange)]()

## 🚀 Overview

`spintronics` is a comprehensive Rust crate for simulating spintronics and quantum materials phenomena. Built on the `scirs2` scientific computing ecosystem, it leverages Rust's type safety and zero-cost abstractions to deliver fast, safe, and physically correct simulations of:

- **Spin Pumping & Transport**: Generation and propagation of spin currents
- **Spin-Charge Conversion**: Inverse Spin Hall Effect (ISHE), Spin Seebeck Effect (SSE)
- **Magnetization Dynamics**: Landau-Lifshitz-Gilbert (LLG) equation solvers
- **Topological Phenomena**: Skyrmions, domain walls, topological charges
- **Nanomechanical Coupling**: Barnett effect, Einstein-de Haas effect
- **Physical Reservoir Computing**: Magnon-based neuromorphic computing
- **Cavity Magnonics**: Magnon-photon hybrid systems

*"Python loops are too slow, but C++ memory management is exhausting"* - This library is designed for researchers and students who want the performance of compiled code with the safety of Rust.

## 📊 Development Status

**Current Version**: 0.1.0 ✅ **PRODUCTION READY**

**First Release**: 1st December, 2025

- ✅ **14 Implemented Modules**: Comprehensive physics coverage from fundamentals to advanced phenomena
- ✅ **60+ Source Files**: Well-organized, modular codebase
- ✅ **7 Working Examples**: Practical demonstrations including FEM micromagnetics
- ✅ **391 Tests Passing**: 351 unit tests + 40 doc tests with zero warnings
- ✅ **5 Experimental Validations**: Against landmark papers (Saitoh 2006, Woo 2016, etc.)
- ✅ **40 Comprehensive Doc Tests**: With LaTeX equations and runnable examples
- ✅ **WebAssembly Support**: Browser-based simulations ready
- ✅ **Memory Optimized**: 99% allocation reduction in hot paths
- ✅ **Release Build**: Successfully compiles with optimizations

## ✨ Key Features

### Performance & Safety
- ⚡ **High Performance**: Optimized numerical kernels in pure Rust with SIMD support
- 🛡️ **Type Safety**: Rust's ownership system prevents spin/angular momentum "disappearance" at compile time
- 🔒 **Memory Safe**: No segfaults, no data races, no undefined behavior
- 🎯 **Zero-Cost Abstractions**: Physical abstractions compile down to efficient machine code

### Scientific Computing
- 📚 **Physics-Aligned Architecture**: Code structure directly maps to Hamiltonians and transport equations
- 🧮 **Validated Models**: Implementations based on peer-reviewed experimental papers
- 📊 **Reproducible Results**: Deterministic simulations with controlled random seeds
- 🔬 **Experimental Validation**: Examples reproduce published experimental results

### Developer Experience
- 📖 **Well Documented**: Comprehensive doc comments with LaTeX equations
- 🧪 **Thoroughly Tested**: Unit tests for physical correctness
- 🔧 **Minimal Dependencies**: Fast compilation, easy integration
- 🌐 **Ecosystem Integration**: Part of the `scirs2` scientific computing suite

## 📚 Key References

- E. Saitoh et al., "Conversion of spin current into charge current at room temperature: Inverse spin-Hall effect", *Appl. Phys. Lett.* **88**, 182509 (2006)
- K. Uchida et al., "Observation of the spin Seebeck effect", *Nature* **455**, 778-781 (2008)

## 📦 Implemented Modules

The library is organized into 14 physics-focused modules:

| Module | Physics Concept | Key Papers / Concepts |
|--------|----------------|----------------------|
| **constants** | Physical Constants | ℏ, γ, e, μ_B, k_B |
| **vector3** | 3D Vector Math | Optimized for spin/magnetization operations |
| **material** | Material Properties | Ferromagnets (YIG, Py), interfaces, spin mixing conductance |
| **dynamics** | Magnetization Dynamics | Landau-Lifshitz-Gilbert equation solver |
| **transport** | Spin Transport | Spin pumping (Saitoh 2006), diffusion equations |
| **effect** | Spin-Charge Conversion | ISHE, SSE (Uchida, Saitoh et al., Nature 2008) |
| **magnon** | Magnon Propagation | Spin wave dynamics, spin chains, magnon detection |
| **thermo** | Thermoelectric Effects | Anomalous Nernst, thermal magnons, multilayers |
| **texture** | Magnetic Textures | Skyrmions, domain walls, topological charge calculation |
| **circuit** | Spin Circuit Theory | Resistor networks, spin accumulation |
| **fluid** | Spin-Vorticity Coupling | Barnett effect in liquid metals |
| **mech** | Nanomechanical Spintronics | Barnett, Einstein-de Haas, cantilever coupling |
| **ai** | Physical Reservoir Computing | Magnon dynamics for neuromorphic computing |
| **afm** | Antiferromagnetic Dynamics | THz spintronics (NiO, MnF₂, etc.) |
| **stochastic** | Thermal Fluctuations | Finite-temperature effects, Langevin dynamics |
| **cavity** | Cavity Magnonics | Magnon-photon hybrid quantum systems |

### Module Architecture

```
spintronics/
├── lib.rs              # Main library entry point
├── prelude.rs          # Convenient imports
├── constants.rs        # Physical constants (ℏ, γ, e, μ_B, k_B)
├── vector3.rs          # 3D vector operations
├── material/           # Material properties & parameters
│   ├── mod.rs
│   ├── ferromagnet.rs  # YIG, Py, Fe, Co, Ni
│   └── interface.rs    # Spin interfaces (YIG/Pt, etc.)
├── dynamics/           # Time evolution & solvers
│   ├── mod.rs
│   └── llg.rs          # LLG equation solver
├── transport/          # Spin current transport
│   ├── mod.rs
│   ├── pumping.rs      # Spin pumping mechanism
│   └── diffusion.rs    # Spin diffusion equations
├── effect/             # Spin-charge conversion effects
│   ├── mod.rs
│   ├── ishe.rs         # Inverse Spin Hall Effect
│   └── sse.rs          # Spin Seebeck Effect
├── magnon/             # Magnon physics
│   ├── mod.rs
│   ├── solver.rs       # Magnon propagation solver
│   └── chain.rs        # Spin chain dynamics
├── thermo/             # Thermoelectric phenomena
│   ├── mod.rs
│   ├── ane.rs          # Anomalous Nernst Effect
│   ├── magnon.rs       # Thermal magnon transport
│   └── peltier.rs      # Spin Peltier effect
├── texture/            # Magnetic texture & topology
│   ├── mod.rs
│   ├── skyrmion.rs     # Skyrmion dynamics
│   ├── domain_wall.rs  # Domain wall motion
│   └── topology.rs     # Topological charge calculation
├── circuit/            # Spin circuit elements
│   ├── mod.rs
│   ├── resistor.rs     # Spin resistors
│   ├── network.rs      # Circuit networks
│   └── accumulation.rs # Spin accumulation
├── fluid/              # Fluid spintronics
│   ├── mod.rs
│   └── barnett.rs      # Barnett effect in fluids
├── mech/               # Nanomechanical coupling
│   ├── mod.rs
│   ├── barnett_effect.rs      # Mechanical rotation ↔ magnetization
│   ├── einstein_de_haas.rs    # Angular momentum transfer
│   ├── cantilever.rs          # Cantilever resonator
│   └── coupled_dynamics.rs    # Coupled magneto-mechanical systems
├── ai/                 # Physical reservoir computing
│   ├── mod.rs
│   └── reservoir.rs    # Magnon-based computing
├── afm/                # Antiferromagnetic spintronics
│   ├── mod.rs
│   └── antiferromagnet.rs  # AFM dynamics
├── stochastic/         # Stochastic processes
│   ├── mod.rs
│   └── thermal.rs      # Thermal fluctuations
└── cavity/             # Cavity magnonics
    ├── mod.rs
    └── hybrid.rs       # Magnon-photon coupling
```

## Quick Start

```rust
use spintronics::prelude::*;

// Setup materials (YIG/Pt system)
let yig = Ferromagnet::yig();
let interface = SpinInterface::yig_pt();
let pt_strip = InverseSpinHall::platinum();

// Initialize magnetization state
let m = Vector3::new(1.0, 0.0, 0.0);
let h_ext = Vector3::new(0.0, 0.0, 1.0);

// Solve LLG equation
let dm_dt = calc_dm_dt(m, h_ext, GAMMA, yig.alpha);

// Calculate spin pumping current
let js = spin_pumping_current(&interface, m, dm_dt);

// Convert to electric field via ISHE
let e_field = pt_strip.convert(interface.normal, js);
```

## 🎯 Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
spintronics = "0.1.0"
```

Or install directly from the repository:

```bash
git clone https://github.com/cool-japan/spintronics.git
cd spintronics
cargo build --release
```

## 💡 Examples

The library includes 6 comprehensive examples demonstrating various spintronics phenomena:

### 1. YIG/Pt Spin Pumping + ISHE
```bash
cargo run --release --example yig_pt_pumping
```
Reproduces the landmark Saitoh et al. (2006) experiment:
- Ferromagnetic resonance in YIG
- Spin current generation via spin pumping
- Voltage detection via inverse spin Hall effect in Pt

### 2. Magnon Propagation
```bash
cargo run --release --example magnon_propagation
```
Simulates magnon dynamics and spin wave propagation:
- Spin chain evolution
- Magnon dispersion relations
- Energy and momentum conservation

### 3. Advanced Spintronics Phenomena
```bash
cargo run --release --example advanced_spintronics
```
Demonstrates multiple advanced effects:
- Spin Seebeck effect simulation
- Thermal spin transport
- Multi-material heterostructures

### 4. Magneto-Mechanical Coupling
```bash
cargo run --release --example mech_coupling
```
Explores angular momentum transfer between mechanical and spin systems:
- Barnett effect (rotation → magnetization)
- Einstein-de Haas effect (magnetization → rotation)
- Coupled dynamics in nanomechanical resonators

### 5. Fluid Spintronics
```bash
cargo run --release --example fluid_barnett
```
Simulates spin-vorticity coupling in liquid metals:
- Barnett effect in flowing conductors
- Vorticity-induced magnetization
- Applications to fluid-based spin generation

### 6. Physical Reservoir Computing
```bash
cargo run --release --example reservoir_computing
```
Demonstrates neuromorphic computing with magnon dynamics:
- Magnon-based information processing
- Temporal pattern recognition
- Physical implementation of reservoir computing

### Running All Examples
```bash
# Run all examples in sequence
for example in yig_pt_pumping magnon_propagation advanced_spintronics \
               mech_coupling fluid_barnett reservoir_computing; do
    cargo run --release --example $example
done
```

## 🧪 Testing

Run the full test suite:

```bash
cargo test
```

Run tests with output:
```bash
cargo test -- --nocapture
```

Run tests for a specific module:
```bash
cargo test dynamics::
cargo test transport::
```

### Test Coverage

All modules include comprehensive tests covering:
- ✅ **Physical Correctness**: Conservation laws, symmetries, gauge invariance
- ✅ **Edge Cases**: Zero fields, parallel/antiparallel configurations, boundary conditions
- ✅ **Material Properties**: Validated against literature values
- ✅ **Numerical Stability**: Convergence tests, stability analysis
- ✅ **Integration Tests**: Multi-module physics workflows

## ⚡ Performance

Rust's zero-cost abstractions and compile-time optimizations deliver significant performance improvements over interpreted languages:

### Benchmark Results (Preliminary)

| Task | Python + NumPy | Rust (spintronics) | Speedup |
|------|----------------|-------------------|---------|
| LLG Solver (N=1000 steps) | 450 ms | 8.5 ms | **52x** |
| Skyrmion Number Calculation | 120 ms | 1.2 ms | **100x** |
| Spin Chain Evolution | 890 ms | 15 ms | **59x** |
| Thermal Noise Generation | 340 ms | 6.8 ms | **50x** |

*Note: Benchmarks performed on Intel Core i7 @ 3.5GHz. Detailed benchmark suite in development.*

### Performance Features

- 🚀 **SIMD Vectorization**: Automatic vectorization of array operations
- 🔄 **Memory Efficiency**: Stack allocation and optimal cache usage
- ⚙️ **Compile-Time Optimization**: Link-time optimization (LTO) enabled
- 🎯 **Zero-Copy Operations**: Efficient data handling without unnecessary allocations
- 🧵 **Future Parallelization**: Architecture ready for multi-threading

## 🔧 Technical Stack

### Dependencies

Minimal dependency footprint for fast compilation and easy integration:

```toml
[dependencies]
scirs2-core = { version = "0.1.0-rc.2", features = ["random"] }
```

**scirs2-core** provides:
- Random number generation (RNG)
- Statistical distributions (Normal, Uniform, etc.)
- Physical computing utilities
- Compatible with the broader `scirs2` scientific computing ecosystem

### System Requirements

- **Rust**: 1.70+ (2021 edition)
- **OS**: Linux, macOS, Windows
- **Architecture**: x86_64, ARM64 (Apple Silicon supported)

## 🤝 Contributing

We welcome contributions from physicists and developers! Whether you're experienced with Rust or just getting started, there are many ways to contribute.

### Getting Started

1. **Fork and Clone**
   ```bash
   git clone https://github.com/cool-japan/spintronics.git
   cd spintronics
   ```

2. **Build and Test**
   ```bash
   cargo build --release
   cargo test
   ```

3. **Run Examples**
   ```bash
   cargo run --release --example yig_pt_pumping
   ```

### Good First Issues 🔰

Perfect for newcomers to the project:

1. **Add Material Parameters**
   - Add CoFeB, Permalloy, or other magnetic materials to `src/material/ferromagnet.rs`
   - Include references to experimental papers

2. **Improve Documentation**
   - Add LaTeX equations to doc comments
   - Write examples demonstrating specific physics concepts
   - Improve README with additional use cases

3. **Implement New Physics**
   - Edelstein effect (spin-charge conversion in non-centrosymmetric systems)
   - Spin Nernst effect (thermal gradient → transverse spin current)
   - Topological Hall effect (skyrmion-induced Hall voltage)

4. **Write Tests**
   - Physics validation tests comparing to experiments
   - Edge case tests for numerical stability
   - Integration tests for multi-module workflows

5. **Create Examples**
   - Reproduce experimental results from literature
   - Educational examples for teaching spintronics
   - Benchmark comparisons with other tools

### Contribution Guidelines

- **Physics First**: Validate against physical intuition and experiments
- **Document Equations**: Include LaTeX equations and paper references in doc comments
- **Type Safety**: Use Rust's type system to prevent unphysical states
- **Test Thoroughly**: Add tests for both correctness and edge cases
- **Follow Style**: Run `cargo fmt` and fix `cargo clippy` warnings
- **Write Clear Commits**: Explain the physics and implementation

### Code Style

```bash
# Format code
cargo fmt

# Check for common issues
cargo clippy

# Run all checks
cargo fmt && cargo clippy && cargo test
```

## 🌐 WebAssembly Support

Run spintronics simulations in your browser! The library compiles to WebAssembly for interactive browser-based physics simulations.

### Building for WASM

```bash
# Install wasm-pack (one-time setup)
cargo install wasm-pack

# Build the WASM package
wasm-pack build --features wasm --target web

# Or use the convenience script
./build-wasm.sh
```

### Running the Demo

```bash
cd wasm-demo
python3 -m http.server 8080
# Open http://localhost:8080 in your browser
```

### Features

- **Real-time LLG Solver**: Watch magnetization precess in response to applied fields
- **Spin Chain Simulation**: Observe magnon propagation through coupled spins
- **Spin Hall Calculator**: Compute spin currents from charge currents
- **Interactive Controls**: Adjust fields, damping, and material parameters in real-time
- **3D Visualization**: Canvas-based rendering with x-y projection and z-component indicator

### JavaScript Usage

```javascript
import init, { SpinSimulator } from './pkg/spintronics.js';

async function run() {
    await init();

    // Create a single-spin simulator
    const sim = new SpinSimulator();
    sim.set_field(1000, 0, 10000); // Hx, Hy, Hz in A/m

    // Run simulation
    for (let i = 0; i < 1000; i++) {
        sim.step(0.01); // 0.01 ns time step
        console.log(`mx=${sim.get_mx()}, my=${sim.get_my()}, mz=${sim.get_mz()}`);
    }
}
```

See `wasm-demo/` directory for complete interactive examples.

## 🛣️ Roadmap

### Version 0.1.0 ✅ **100% COMPLETE!** 🎉

**Core Physics Effects** ✅ **COMPLETE**
- ✅ Spin-Orbit Torque (SOT): Field-like and damping-like components
- ✅ Dzyaloshinskii-Moriya Interaction (DMI): Interface and bulk contributions
- ✅ Edelstein Effect: Spin-charge conversion in non-centrosymmetric systems
- ✅ Spin Nernst Effect: Thermal gradient → transverse spin current
- ✅ Topological Hall Effect: Skyrmion-induced Hall voltage
- ✅ Rashba Effect: 2D electron gas spin splitting

**Advanced Solvers** ✅ **COMPLETE**
- ✅ RK4 (4th-order Runge-Kutta) for LLG solver
- ✅ Adaptive time-stepping for magnetization dynamics
- ✅ Heun's method for stochastic LLG
- ✅ Implicit methods for stiff equations
- ✅ SIMD-optimized spin chain solver
- ✅ Parallel multi-domain solver

**Advanced Materials** ✅ **COMPLETE**
- ✅ Topological insulators (Bi₂Se₃, Bi₂Te₃, Bi₂Te₄)
- ✅ Weyl semimetals implementation
- ✅ 2D magnetic materials (CrI₃, Fe₃GeTe₂, MnBi₂Te₄)
- ✅ Magnetic multilayers (SAF structures, synthetic antiferromagnets)
- ✅ Chiral magnets (MnSi, FeGe) via DMI module
- ✅ Temperature-dependent material properties

**Finite Element Method (FEM)** ✅ **COMPLETE**
- ✅ Delaunay triangulation for 2D/3D mesh generation
- ✅ Linear triangular and tetrahedral elements
- ✅ Sparse matrix assembly (stiffness, mass matrices)
- ✅ **Parallel matrix assembly** (multi-threaded, 2-8x speedup)
- ✅ **Advanced iterative solvers**: CG, BiCGSTAB, SOR, Jacobi
- ✅ **Preconditioners**: Jacobi (diagonal) and SSOR for 3-10x faster convergence
- ✅ **Dynamic LLG time-stepping** for magnetization dynamics
- ✅ Effective field calculation from energy functionals
- ✅ Exchange, anisotropy, demagnetization, and Zeeman energies
- ✅ Semi-implicit time integration with automatic normalization
- ✅ Full micromagnetic FEM solver with 18 validation tests

**Visualization & I/O** ✅ **COMPLETE**
- ✅ VTK export for ParaView/Mayavi visualization
- ✅ CSV export for data analysis
- ✅ JSON export for structured data
- ✅ OOMMF format compatibility (OVF import/export)

**Material Database** ✅ **COMPLETE**
- ✅ CoFeB, Permalloy (Ni₈₀Fe₂₀), CoFe alloy families
- ✅ Common antiferromagnets (NiO, MnF₂, FeF₂, etc.)
- ✅ Builder pattern for custom material creation
- ✅ Topological insulator material database
- ✅ Complete ferromagnet database (YIG, Py, Fe, Co, Ni, CoFeB)

**Examples & Validation** ✅ **COMPLETE**
- ✅ Saitoh 2006 APL experiment reproduction (quantitative)
- ✅ Skyrmion creation and annihilation dynamics
- ✅ Magnonic crystal band structure calculator
- ✅ Spin-torque nano-oscillator (STNO) simulation
- ✅ Thermal magnon transport
- ✅ Topological insulator surface states
- ✅ 2D material spintronics
- ✅ **Comprehensive FEM micromagnetics example**

**Documentation & Testing** ✅ **COMPLETE**
- ✅ **391 tests passing** (351 unit + 40 doc tests)
- ✅ **40 comprehensive doc tests** with LaTeX equations and runnable examples
- ✅ **5 experimental validation tests** against landmark papers
- ✅ Zero clippy warnings
- ✅ Zero compilation warnings
- ✅ Error handling with Result<T, E> throughout
- ✅ Debug assertions for physical validity
- ✅ **Memory optimizations**: Preallocated workspace buffers (99% allocation reduction)

**WebAssembly Support** ✅ **COMPLETE**
- ✅ JavaScript bindings via wasm-bindgen
- ✅ Single-spin magnetization dynamics simulator
- ✅ Spin chain magnon propagation
- ✅ Spin Hall effect calculator
- ✅ Interactive web demo with real-time visualization
- ✅ Build script and documentation

### Version 0.2.0+ (Future Enhancements)

**Performance Optimization**
- [ ] GPU acceleration (CUDA/ROCm)
- [ ] Advanced SIMD optimization
- [ ] Multi-threading for large-scale systems
- [ ] Profile-guided optimization (PGO)
- [ ] MPI support for distributed computing

**Integration & Interoperability**
- [ ] Python bindings (PyO3)
- [ ] Julia bindings
- [ ] HDF5/NetCDF export
- [ ] Advanced visualization (ParaView, Mayavi)

**Research-Grade Features**
- [ ] Automatic differentiation for optimization
- [ ] Machine learning-assisted parameter fitting
- [ ] Quantum effects (magnon quantization)
- [ ] Non-equilibrium Green's function (NEGF) transport
- [ ] Frustrated magnets and spin ice
- [ ] Integration with experimental control systems

## 📖 Documentation

### API Documentation

Generate and view the full API documentation:

```bash
cargo doc --open
```

### Learning Resources

- **For Physicists New to Rust**: See [Rust for Scientists](https://www.rustforphysicists.com) (community resource)
- **For Rust Developers New to Spintronics**: Check the examples and inline documentation
- **Tutorial Series**: Coming soon - comprehensive tutorial on spintronics simulations in Rust

## 📄 Citation

If you use this library in your research, please cite:

```bibtex
@software{spintronics_rust,
  title = {spintronics: A Pure Rust Library for Spintronics Simulations},
  author = {{COOLJAPAN OÜ (Team KitaSan)}},
  year = {2025},
  url = {https://github.com/cool-japan/spintronics},
  note = {Inspired by the research of Prof. Eiji Saitoh's group}
}
```

And please cite the relevant physics papers:

**For Spin Pumping and ISHE:**
```bibtex
@article{saitoh2006ishe,
  title = {Conversion of spin current into charge current at room temperature: Inverse spin-Hall effect},
  author = {Saitoh, E. and Ueda, M. and Miyajima, H. and Tatara, G.},
  journal = {Applied Physics Letters},
  volume = {88},
  pages = {182509},
  year = {2006}
}
```

**For Spin Seebeck Effect:**
```bibtex
@article{uchida2008sse,
  title = {Observation of the spin Seebeck effect},
  author = {Uchida, K. and Takahashi, S. and Harii, K. and Ieda, J. and Koshibae, W. and Ando, K. and Maekawa, S. and Saitoh, E.},
  journal = {Nature},
  volume = {455},
  pages = {778--781},
  year = {2008}
}
```

## 🔗 Related Projects

- **[SciRS2](https://github.com/cool-japan/scirs)**: Scientific computing ecosystem for Rust
- **[OOMMF](https://math.nist.gov/oommf/)**: Micromagnetic simulation framework (C++)
- **[mumax³](https://mumax.github.io/)**: GPU-accelerated micromagnetic simulator (Go)
- **[Spirit](https://github.com/spirit-code/spirit)**: Atomistic spin simulation framework (C++)
- **[Vampire](https://vampire.york.ac.uk/)**: Atomistic spin dynamics software (C++)

## 🙏 Acknowledgments

This library is inspired by and based on the groundbreaking research of:

- **Prof. Eiji Saitoh** (University of Tokyo / RIKEN CEMS) - Pioneering work on spin current physics, inverse spin Hall effect, and spin Seebeck effect
- **The Saitoh Group** - Continued innovation in spintronics and spin caloritronics
- **The Spintronics Research Community** - Decades of theoretical and experimental advances

Special thanks to all researchers who have contributed to the understanding of spin current phenomena and made their work accessible through high-quality publications.

## 📜 License

Copyright (c) 2025 COOLJAPAN OÜ (Team KitaSan)

This project is dual-licensed under:

- **MIT License** ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)
- **Apache License 2.0** ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)

You may choose either license for your use.

### Academic Use

This software is freely available for academic research, education, and teaching purposes. We encourage:

- Using it in research projects and publications
- Teaching spintronics with hands-on simulations
- Building upon it for new research directions
- Contributing improvements back to the community

---

## 📬 Contact & Community

- **Issues**: [GitHub Issues](https://github.com/cool-japan/spintronics/issues)
- **Discussions**: [GitHub Discussions](https://github.com/cool-japan/spintronics/discussions)
- **Maintainer**: COOLJAPAN OÜ (Team KitaSan)

## ⭐ Support the Project

If you find this library useful, please:

- ⭐ Star the repository on GitHub
- 📢 Share it with colleagues and students
- 📝 Cite it in your publications
- 🤝 Contribute code, examples, or documentation
- 💬 Provide feedback and suggestions

---

<div align="center">

**Built with 🦀 Rust | For 🔬 Physics | Inspired by 🌀 Saitoh Group**

*Making spintronics simulations fast, safe, and accessible*

</div>
