# Changelog

All notable changes to the spintronics library will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Interactive Web Demonstration Subcrate (`spintronics-demo`)** (v0.2.0):
  - Modern HTMX + Axum + Askama stack for server-side rendering
  - 4 interactive physics demonstrations:
    - LLG Magnetization Dynamics: Real-time solver with trajectory visualization
    - Spin Pumping Calculator: Reproduces Saitoh 2006 APL experiment
    - Materials Explorer: Compare magnetic properties across ferromagnets
    - Skyrmion Visualizer: Real-time magnetization field rendering
  - Zero JavaScript frameworks - progressive enhancement with HTMX
  - Full library access on server-side (no WASM limitations)
  - Type-safe templates with Askama
  - Comprehensive documentation and deployment guide
- rustfmt.toml and clippy.toml configuration files for code style
- **Vector3 convenience methods** (v0.2.0):
  - `zero()` - Create zero vector
  - `unit_x()`, `unit_y()`, `unit_z()` - Unit vectors along coordinate axes
  - `magnitude_squared()` - Squared magnitude (avoids sqrt for performance)
  - `is_normalized()` - Check if vector is unit length
  - `angle_between()` - Calculate angle between two vectors
  - `project()` - Vector projection operation
  - Performance: All hot-path methods marked with `#[inline]` for optimization (8 methods)
- **Additional interface materials** (v0.2.0):
  - Platinum interfaces: `cofeb_pt()`, `co_pt()`, `fe_pt()` (3 materials)
  - Tantalum interfaces: `yig_ta()`, `py_ta()`, `cofeb_ta()` (3 materials)
  - Tungsten interfaces: `cofeb_w()`, `py_w()` (2 materials)
  - Total: 8 new FM/NM interface combinations
  - Builder methods: `with_g_r`, `with_g_i`, `with_normal`, `with_area`
- **Energy calculation utilities** (v0.2.0):
  - `zeeman_energy()` - Zeeman energy from applied magnetic field
  - `anisotropy_energy()` - Uniaxial anisotropy energy
  - `exchange_energy()` - Exchange energy for non-uniform magnetization
  - All functions marked with `#[inline]` for performance
- **Default implementations for texture types** (v0.2.0):
  - `Skyrmion::default()` - Néel-type skyrmion with CCW chirality, 50 nm radius
  - `DomainWall::default()` - Bloch-type wall with 10 nm width
- **Complete trait implementations for enums** (v0.2.0):
  - Added `Eq` and `Hash` to all simple enums for better API ergonomics
  - Texture enums: `Helicity`, `Chirality`, `LatticeType`, `WallType`, `DmiType`
  - Material enums: `TopologicalClass`, `WeylType`, `MagneticState`, `MagneticOrdering`, `AfmStructure`, `MultilayerType`
  - Enables use in HashMaps, HashSets, and other collections
- **Extended builder methods** (v0.2.0):
  - `Antiferromagnet`: with_sublattice_magnetization, with_exchange_field, with_anisotropy_field, with_resonance_frequency, with_spin_hall_angle
  - Complete builder pattern coverage for 7 additional fields
- Enhanced prelude with additional commonly used types:
  - Thermal effects: AnomalousNernst, SpinPeltier
  - Magnetic textures: Skyrmion, SkyrmionLattice, DomainWall, Chirality, Helicity, TopologicalCharge, WallType
  - Topological functions: calculate_skyrmion_number
- Module-level prelude convenience imports:
  - `effect::prelude` with ISHE, SOT, SNE, SSE, THE type aliases
  - `material::prelude` with FM, TI, WSM, AFM type aliases
  - `texture::prelude` with DMI, DW, Sk type aliases
  - `thermo::prelude` with ANE, SPE type aliases
- Builder methods for additional types:
  - `SpinSeebeck`: with_l_s, with_g_th, with_polarization
  - `AnomalousNernst`: with_alpha_ane, with_magnetization
  - `SpinPeltier`: with_pi_s, with_temperature, with_area
  - `Skyrmion`: with_center, with_radius, with_helicity, with_chirality
  - `DomainWall`: with_center, with_width, with_type, with_normal
- Serde serialization for thermal effect types:
  - AnomalousNernst, SpinPeltier
  - MagnonThermalConductivity, ThermalMagnonTransport
  - Layer, ThermalBoundary, MultilayerStack
- Serde serialization for texture types:
  - Skyrmion, SkymionLattice, Helicity, Chirality, LatticeType
  - DomainWall, WallType
- Display implementations for thermal types:
  - Layer, ThermalBoundary, MultilayerStack
  - MagnonThermalConductivity, ThermalMagnonTransport
- CHANGELOG.md with full version history
- Unit validation module (`units.rs`):
  - Physical quantity validators for magnetization, damping, exchange stiffness
  - Temperature, magnetic field, and thickness range checks
  - Spin Hall angle, resistivity, and DMI constant validation
  - Current density, voltage, and energy scale validators
  - 14 validation functions with comprehensive test coverage
- Examples organization:
  - Comprehensive `examples/README.md` with difficulty levels (Basic/Intermediate/Advanced)
  - 17 examples categorized by complexity and physics domain
  - Learning paths for different user backgrounds
  - Difficulty indicators (⭐/⭐⭐/⭐⭐⭐) added to example files
- Main README updated:
  - Version 0.2.0 highlights and new features
  - 18 modules documented (added units, memory, visualization, python)
  - Optional features guide (python, hdf5, serde, fem, wasm)
  - Updated examples section with links to detailed guide
- lib.rs documentation enhanced:
  - Module count updated (14 → 18)
  - Modules organized by category (Core, Materials, Effects, etc.)
  - Unit validation usage example added
  - v0.2.0 features prominently highlighted
  - Test count: 431 passing (381 unit + 50 doc tests)

### Changed
- Cargo.toml: Added rust-version (MSRV 1.70.0) and homepage metadata
- Simplified `skyrmion_dynamics` example to use enhanced prelude (removed redundant imports)
- **Performance optimizations** (v0.2.0):
  - Added `#[inline]` to critical hot-path functions (21 functions total):
    - Dynamics: `calc_dm_dt()` (LLG equation)
    - Transport: `spin_pumping_current()` (spin pumping)
    - ISHE: `convert()`, `voltage()` (inverse spin Hall effect)
    - SOT: `damping_like_field()`, `field_like_field()` (spin-orbit torque)
    - SSE: `spin_current()`, `interface_current()` (spin Seebeck effect)
    - SNE: `spin_current()`, `heat_current()` (spin Nernst effect)
    - Rashba: `spin_texture()`, `edelstein_spin_density()`, `inverse_edelstein_current()` (Rashba-Edelstein effects)
    - ANE: `electric_field()`, `voltage()` (anomalous Nernst effect)
    - Spin Peltier: `heat_current()`, `temperature_change_rate()` (spin Peltier effect)
    - Magnon thermal: `conductivity_at_temperature()`, `heat_flux()`, `magnon_chemical_potential()`, `thermal_magnon_accumulation()` (thermal magnon transport)
  - Enables aggressive compiler inlining for ~10-30% performance improvement in tight loops

### Fixed
- Typo: Renamed `SkymionLattice` to `SkyrmionLattice` (missing 'r')
- Documentation warnings: Fixed 13 rustdoc warnings (unit bracket escaping, HTML tag escaping)
- HDF5 feature: Fixed `VarLenUnicode` string conversion for hdf5 0.8.x API compatibility
- WASM feature: Added `wasm_js` feature to getrandom for proper WASM32 support
- Feature gates: Fixed `MultiDomainSystem` to be properly gated behind `scirs2` feature
- Module exports: Added energy calculation utilities to `dynamics` module exports

## [0.2.0] - 2025-12-24

### Added

#### Python Bindings (PyO3)
- `PyVector3`: 3D vector with all arithmetic operations
- `PyFerromagnet`: Material parameters (YIG, Permalloy, CoFe, etc.)
- `PySpinInterface`: Spin mixing conductance calculations
- `PyInverseSpinHall`: ISHE converter (Pt, Ta, W materials)
- `PyLlgSimulator`: LLG equation solver with RK4/Euler methods
- `PySpinPumpingSimulation`: Complete spin pumping workflow
- Physical constants exported to Python (HBAR, GAMMA, E_CHARGE, MU_B, KB)

#### Serialization Support (serde)
- `Vector3` serialization/deserialization
- `Ferromagnet`, `SpinInterface`, `InverseSpinHall` serialization
- `SpinSeebeck`, `SpinOrbitTorque` serialization
- `SimulationData` JSON export
- `Magnetic2D`, `MagneticOrdering` serialization
- `TopologicalInsulator`, `TopologicalClass` serialization
- `WeylSemimetal`, `WeylType`, `MagneticState` serialization
- `DmiParameters`, `DmiType` serialization

#### HDF5 Export Support
- `Hdf5Writer`: Write scalars, arrays, and vector fields
- `Hdf5Reader`: Read scalars, arrays, and vector fields
- Hierarchical group support for organized data
- Time series export capabilities
- Graceful fallback when HDF5 feature is disabled

#### Memory Pool Allocator
- `VectorPool<T>`: Generic vector pool for efficient f64 allocation
- `SpinArrayPool`: Specialized pool for `Vec<Vector3<f64>>`
- Thread-local pools for convenience
- `Rk4Workspace`: Preallocated buffers for RK4 integration
- `HeunWorkspace`: Preallocated buffers for Heun/stochastic solvers

#### API Improvements
- `Display` trait for key types (Vector3, Ferromagnet, SpinInterface, etc.)
- `Display` for SpinSeebeck, SpinNernst, SpinOrbitTorque, TopologicalHall, RashbaSystem
- `Display` for Magnetic2D, MagneticOrdering, DmiParameters, Skyrmion, DomainWall
- `Display` for AnomalousNernst, SpinPeltier
- `Default` implementations for SpinOrbitTorque, SpinNernst, TopologicalHall, RashbaSystem
- `Default` for TopologicalInsulator, WeylSemimetal, Magnetic2D
- Builder methods for `SpinOrbitTorque` and `InverseSpinHall`
- Trait hierarchy: `MagneticMaterial`, `SpinChargeConverter`, `TopologicalMaterial`

#### Extended Physical Constants
- Fundamental: HBAR, H_PLANCK, E_CHARGE, KB, C_LIGHT, NA
- Electromagnetic: MU_0, EPSILON_0, ALPHA_FS
- Magnetic: GAMMA, MU_B, MU_N, G_LANDE
- Particle: ME, MP, E_OVER_ME
- Derived: SPIN_QUANTUM, FLUX_QUANTUM, CONDUCTANCE_QUANTUM

#### Community
- CONTRIBUTING.md guide for contributors
- CODE_OF_CONDUCT.md (Contributor Covenant)
- GitHub issue templates (bug report, feature request, question)

#### Infrastructure
- GitHub Actions CI/CD workflow
- Multi-platform testing (Ubuntu, macOS, Windows)
- Clippy and Rustfmt checks
- WASM build verification
- Documentation builds
- MSRV testing (Rust 1.70.0)

### Changed
- Organized prelude imports by category
- Extended physical constants exports in prelude

## [0.1.0] - 2025-11-15

### Added

#### Core Physics Effects
- **Spin-Orbit Torque (SOT)**: Field-like and damping-like torque components
- **Dzyaloshinskii-Moriya Interaction (DMI)**: Interface and bulk contributions
- **Edelstein Effect**: Spin-charge conversion in non-centrosymmetric systems
- **Spin Nernst Effect**: Thermal gradient to transverse spin current
- **Topological Hall Effect**: Skyrmion-induced Hall voltage
- **Rashba Effect**: 2D electron gas spin splitting

#### Solvers and Algorithms
- RK4 (4th-order Runge-Kutta) for LLG solver
- Adaptive time-stepping for dynamics
- Heun's method for stochastic LLG
- Implicit methods for stiff equations
- SIMD-optimized spin chain solver
- Parallel solver for multi-domain systems

#### Materials
- Topological insulators: Bi₂Se₃, Bi₂Te₃, Bi₂Te₄
- Weyl semimetals implementation
- 2D magnetic materials: CrI₃, Fe₃GeTe₂, MnBi₂Te₄
- Magnetic multilayers (SAF, synthetic antiferromagnets)
- Chiral magnets: MnSi, FeGe (in DMI module)
- Temperature-dependent material properties
- CoFeB, Permalloy, CoFe alloy parameters
- Common antiferromagnets: NiO, MnF₂, etc.
- Topological insulator material database

#### Finite Element Method
- Delaunay mesh generation (2D/3D)
- Linear triangular and tetrahedral elements
- Sparse matrix assembly (stiffness, mass)
- Parallel matrix assembly
- Iterative solvers: CG, BiCGSTAB, SOR, Jacobi
- Preconditioners: Jacobi, SSOR
- Micromagnetic FEM solver
- Energy calculations: Exchange, anisotropy, demagnetization, Zeeman

#### WebAssembly
- wasm-bindgen JavaScript bindings
- Single-spin LLG simulator
- Spin chain magnon propagation
- Spin Hall effect calculator
- Interactive web demo

#### Visualization
- VTK export
- CSV export
- JSON export
- OOMMF format import/export

#### Examples
- Saitoh 2006 APL experiment reproduction
- Skyrmion creation and annihilation
- Magnonic crystal band structure
- Spin-torque nano-oscillator (STNO)
- Thermal magnon transport
- Topological insulator surface states
- 2D material spintronics

#### Documentation
- Comprehensive doc tests (40 passing)
- LaTeX equations from papers
- Physics validation tests
- API documentation examples

### Fixed
- All `cargo clippy` warnings
- Memory allocation optimizations in hot paths
- Workspace buffer reuse for solvers

---

## Version Policy

This project follows [Semantic Versioning](https://semver.org/):

- **MAJOR** version: Incompatible API changes
- **MINOR** version: New functionality in a backwards compatible manner
- **PATCH** version: Backwards compatible bug fixes

## Links

- [Repository](https://github.com/cool-japan/spintronics)
- [Documentation](https://docs.rs/spintronics)
- [crates.io](https://crates.io/crates/spintronics)
