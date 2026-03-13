# Examples

This directory contains **25 example programs** demonstrating various features of the spintronics library, organized by difficulty level.

## 📚 Basic Examples (Beginner-friendly)

Start here if you're new to the library or spintronics simulations.

### Core Physics
- **`yig_pt_pumping.rs`** - Basic spin pumping in YIG/Pt bilayers
  - Demonstrates: Ferromagnetic resonance, spin pumping, ISHE detection
  - Physics: Saitoh et al. (2006) experiment reproduction
  - Difficulty: ⭐ Beginner
  - Run: `cargo run --example yig_pt_pumping`

- **`saitoh_2006_ishe.rs`** - Inverse Spin Hall Effect validation
  - Demonstrates: Spin current to charge current conversion
  - Physics: Quantitative validation against experimental data
  - Difficulty: ⭐ Beginner
  - Run: `cargo run --example saitoh_2006_ishe`

### Dynamics
- **`llg_dynamics_visualization.rs`** - LLG equation solver basics
  - Demonstrates: Magnetization precession, damping, RK4 solver
  - Physics: Landau-Lifshitz-Gilbert dynamics
  - Difficulty: ⭐ Beginner
  - Run: `cargo run --example llg_dynamics_visualization`

---

## 🔬 Intermediate Examples

These examples introduce more complex physics and multi-component systems.

### Magnon Physics
- **`magnon_propagation.rs`** - Spin wave propagation in chains
  - Demonstrates: Magnon dispersion, group velocity, damping
  - Physics: Spin wave dynamics in 1D systems
  - Difficulty: ⭐⭐ Intermediate
  - Run: `cargo run --example magnon_propagation`

### Topological Textures
- **`skyrmion_visualization.rs`** - Skyrmion creation and visualization
  - Demonstrates: Magnetic skyrmions, topological charge, DMI
  - Physics: Néel and Bloch skyrmions
  - Difficulty: ⭐⭐ Intermediate
  - Run: `cargo run --example skyrmion_visualization`

- **`skyrmion_dynamics.rs`** - Skyrmion motion and switching
  - Demonstrates: Current-driven skyrmion motion, skyrmion Hall effect
  - Physics: Topological spin textures under current
  - Difficulty: ⭐⭐ Intermediate
  - Run: `cargo run --example skyrmion_dynamics`

### Spin-Orbit Effects
- **`spin_torque_oscillator.rs`** - STNO auto-oscillations
  - Demonstrates: Spin-transfer torque, auto-oscillation, phase locking
  - Physics: Non-linear magnetization dynamics
  - Difficulty: ⭐⭐ Intermediate
  - Run: `cargo run --example spin_torque_oscillator`

### 2D Materials
- **`2d_material_spintronics.rs`** - Van der Waals heterostructures
  - Demonstrates: CrI₃, Fe₃GeTe₂, proximity effects
  - Physics: 2D magnetism and spintronics
  - Difficulty: ⭐⭐ Intermediate
  - Run: `cargo run --example 2d_material_spintronics`

### New in v0.3.0
- **`simulation_builder.rs`** - SimulationBuilder fluent API
  - Demonstrates: Builder pattern for setting up complex simulations
  - Physics: Comprehensive spin dynamics configuration
  - Difficulty: ⭐⭐ Intermediate
  - Run: `cargo run --example simulation_builder`

---

## 🚀 Advanced Examples

Complex simulations requiring understanding of multiple physics domains.

### Computational Methods
- **`fem_micromagnetics.rs`** - Finite element micromagnetics
  - Demonstrates: FEM mesh, iterative solvers, energy minimization
  - Physics: Micromagnetic simulations with realistic geometries
  - Difficulty: ⭐⭐⭐ Advanced
  - Requirements: `fem` feature (`cargo run --features fem --example fem_micromagnetics`)

- **`parallel_magnon_dynamics.rs`** - Multi-threaded magnon solver
  - Demonstrates: Parallel computing, multi-domain systems
  - Physics: Large-scale spin wave simulations
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --release --example parallel_magnon_dynamics`

### Collective Phenomena
- **`magnonic_crystal.rs`** - Band structure calculations
  - Demonstrates: Magnon dispersion, band gaps, Brillouin zone
  - Physics: Periodic magnetic structures
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example magnonic_crystal`

### Thermal Effects
- **`thermal_magnon_transport.rs`** - Thermoelectric effects
  - Demonstrates: Spin Seebeck effect, thermal magnon transport
  - Physics: Temperature gradients and spin currents
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example thermal_magnon_transport`

### Topological Materials
- **`topological_insulator.rs`** - Surface states and transport
  - Demonstrates: Topological surface states, Edelstein effect
  - Physics: Bi₂Se₃, spin-momentum locking
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example topological_insulator`

### Emerging Technologies
- **`reservoir_computing.rs`** - Neuromorphic computing with magnons
  - Demonstrates: Physical reservoir computing, pattern recognition
  - Physics: Non-linear magnon dynamics for ML
  - Difficulty: ⭐⭐⭐ Advanced
  - Requirements: `scirs2` feature
  - Run: `cargo run --example reservoir_computing`

- **`advanced_spintronics.rs`** - Multi-effect integration
  - Demonstrates: SOT, DMI, thermal effects combined
  - Physics: Comprehensive spintronics device simulation
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example advanced_spintronics`

### Coupled Systems
- **`mech_coupling.rs`** - Magneto-mechanical coupling
  - Demonstrates: Barnett effect, Einstein-de Haas effect
  - Physics: Coupling between magnetization and rotation
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example mech_coupling`

- **`fluid_barnett.rs`** - Spin-vorticity coupling in fluids
  - Demonstrates: Barnett effect in rotating liquid metals
  - Physics: Fluid dynamics and magnetization
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example fluid_barnett`

### New in v0.3.0
- **`altermagnet_ruo2.rs`** - Altermagnetism in RuO₂
  - Demonstrates: Time-reversal symmetry breaking without net moment, anomalous Hall response
  - Physics: Altermagnetism, RuO₂ material model
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example altermagnet_ruo2`

- **`orbital_hall_effect.rs`** - Orbital Hall effect and orbital currents
  - Demonstrates: Orbital current generation, orbital Hall conductivity
  - Physics: Orbitronics, spin-orbit coupling
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example orbital_hall_effect`

- **`spin_ice_monopoles.rs`** - Frustrated spin ice and magnetic monopoles
  - Demonstrates: Kagome and pyrochlore spin ice, monopole excitations
  - Physics: Frustrated magnets, emergent magnetic monopoles
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example spin_ice_monopoles`

- **`spin_wave_dispersion.rs`** - Spin wave dispersion and magnon dynamics
  - Demonstrates: Dispersion relations, magnon group velocity, band structure
  - Physics: Linear spin wave theory, exchange and dipolar contributions
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example spin_wave_dispersion`

- **`hopfion_3d.rs`** - Hopfion topology in 3D
  - Demonstrates: Hopfion creation, Hopf invariant calculation, 3D topology
  - Physics: Three-dimensional topological solitons
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example hopfion_3d`

- **`magnon_bec.rs`** - Magnon Bose-Einstein condensation
  - Demonstrates: Magnon condensation, coherent magnon states, critical temperature
  - Physics: Quantum magnon statistics, BEC in driven magnetic systems
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example magnon_bec`

- **`magnetoelastic_switching.rs`** - Magnetoelastic coupling and magnetostriction
  - Demonstrates: Strain-induced switching, magnetoelastic energy, magnetostriction
  - Physics: Magnetoelastic coupling, coupled spin-lattice dynamics
  - Difficulty: ⭐⭐⭐ Advanced
  - Run: `cargo run --example magnetoelastic_switching`

---

## 🎯 Learning Path

### Path 1: Basic Spintronics
1. `yig_pt_pumping.rs` - Understand spin pumping basics
2. `saitoh_2006_ishe.rs` - Learn spin-to-charge conversion
3. `llg_dynamics_visualization.rs` - Master magnetization dynamics
4. `magnon_propagation.rs` - Study spin wave physics

### Path 2: Topological Spintronics
1. `llg_dynamics_visualization.rs` - Learn dynamics fundamentals
2. `skyrmion_visualization.rs` - Understand skyrmion structures
3. `skyrmion_dynamics.rs` - Explore current-driven motion
4. `topological_insulator.rs` - Advanced topological effects

### Path 3: Device Engineering
1. `yig_pt_pumping.rs` - Basic device structure
2. `spin_torque_oscillator.rs` - Active spintronic devices
3. `2d_material_spintronics.rs` - Modern material platforms
4. `advanced_spintronics.rs` - Complete device simulation

### Path 4: Computational Methods
1. `llg_dynamics_visualization.rs` - Basic numerical methods
2. `magnon_propagation.rs` - 1D simulations
3. `fem_micromagnetics.rs` - Advanced numerical techniques
4. `parallel_magnon_dynamics.rs` - High-performance computing

### Path 5: New in v0.3.0 - Emerging Quantum Materials
1. `altermagnet_ruo2.rs` - Altermagnetism without net moment
2. `orbital_hall_effect.rs` - Orbital currents and orbitronics
3. `magnon_bec.rs` - Quantum magnon condensation
4. `hopfion_3d.rs` - Three-dimensional topological solitons

### Path 6: Frustrated Magnetism and Complex Textures
1. `spin_ice_monopoles.rs` - Frustrated lattices and emergent monopoles
2. `spin_wave_dispersion.rs` - Magnon dispersion in frustrated systems
3. `magnetoelastic_switching.rs` - Coupled spin-lattice dynamics
4. `hopfion_3d.rs` - 3D topological spin structures

---

## 📋 Feature Requirements

Some examples require optional features to be enabled:

- **FEM examples**: `cargo run --features fem --example fem_micromagnetics`
- **Reservoir computing**: `cargo run --features scirs2 --example reservoir_computing`
- **Python interop**: Enable `python` feature for PyO3 bindings

## 🔧 Building and Running

```bash
# Run a specific example
cargo run --example yig_pt_pumping

# Run with release optimizations (recommended for heavy simulations)
cargo run --release --example parallel_magnon_dynamics

# Run with specific features
cargo run --features fem --example fem_micromagnetics

# Build all examples
cargo build --examples --release
```

## 📖 Additional Resources

- **Main documentation**: `cargo doc --open`
- **API reference**: https://docs.rs/spintronics
- **Tutorial series**: See `docs/tutorials/` (coming soon)
- **Physics background**: See `README.md` for key references

## 🤝 Contributing Examples

We welcome contributions of new examples! When adding an example:

1. Choose appropriate difficulty level
2. Add clear documentation header
3. Include physics references
4. Add entry to this README
5. Ensure it runs without warnings

See `CONTRIBUTING.md` for detailed guidelines.
