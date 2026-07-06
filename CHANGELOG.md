# Changelog

All notable changes to the spintronics library will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.2] - 2026-07-06

### Added

#### Altermagnets (`src/altermagnet/`)
- `AltermagnetBandModel`: momentum-resolved, two-sublattice k·p Bloch-Hamiltonian model
  for collinear altermagnets, with `mnte()`, `crsb()`, and `ruo2()` presets, closed-form
  Berry curvature, and Fermi-sea `crystal_hall_conductivity()` / `spin_hall_conductivity()`
  integrals that are correctly zero without spin-orbit coupling and nonzero (with sign
  flipping under Néel reversal) once an optional `lambda_soc` term is enabled
- `AltermagnetSpinValve`: giant magnetoresistance without ferromagnetism, driven purely by
  the relative crystal-axis angle between two zero-net-moment altermagnetic layers

#### Orbitronics (`src/orbitronics/`)
- `crystal_field` module: real-cubic-harmonic d-orbital angular-momentum operators and
  octahedral/tetrahedral/tetragonal crystal-field + spin-orbit Hamiltonians
- `d_orbital_moment` module: `CrystalFieldModel` with Hund's-rule free-ion terms,
  high-spin/low-spin ground states, orbital quenching and spin-orbit-driven unquenching,
  effective magnetic moments, and 12 preset 3d transition-metal ions spanning Ti³⁺ to
  Cu²⁺, including both high-spin and low-spin configurations (e.g. Fe²⁺, Co³⁺)

#### Frustrated Magnetism (`src/frustrated/`)
- `rvb` module: resonating-valence-bond quantum spin-liquid solver — dimer-covering
  enumeration, Sutherland loop-counting overlaps, a generalized-eigenproblem variational
  ground state (Löwdin canonical orthogonalization), exact diagonalization (dense, plus a
  matrix-free Lanczos solver in the S_z=0 sector for up to 16 sites), and spinon-pair /
  deconfinement diagnostics (53 tests)
- `transport` module (`FrustratedTransport`): per-plaquette scalar spin chirality →
  emergent magnetic field → topological Hall response pipeline for triangular, kagome,
  and pyrochlore lattices
- `FrustratedLattice::set_umbrella_order()`: noncoplanar, canted 120-degree "umbrella"
  spin-texture initializer
- `TopologicalHall::emergent_field_from_solid_angle()` /
  `::hall_resistivity_from_charge_density()`: generalized entry points for arbitrary
  noncoplanar spin textures, not just discrete skyrmions

#### Hopfions (`src/texture/`)
- `hopfion_stability_modes`: collective-coordinate linear-stability (eigenmode) analysis
  of relaxed hopfions, diagonalizing the Hessian of the total energy functional
- `Hopfion::with_profile_and_offset()`: rigidly-translated hopfion ansatz, used to probe
  the translation zero-modes

#### Magnonics (`src/magnon/`)
- `SdwRelaxationDynamics`: time-dependent Ginzburg-Landau (Model A) relaxational dynamics
  for the spin-density-wave order-parameter amplitude, relaxing toward `SdwGapSolver`'s
  self-consistent equilibrium gap
- `ParametricAmplification::with_supercriticality(xi)`: rescale the parametric pump field
  to a chosen operating point relative to threshold

#### Straintronics (`src/mech/`)
- `strain_driven_dynamics`: time-domain LLG driver coupling oscillating strain — surface
  acoustic wave (`SawStrainDrive`) or AC piezoelectric (`PiezoAcStrainDrive`) — into a
  magnetoelastic effective field, integrated with the existing Dormand-Prince adaptive
  solver

#### Materials (`src/material/`)
- `GradedInterface` / `GradingLaw`: `Linear`, `Exponential`, and `ErrorFunction` grading
  laws for spatially-varying Ms/A/K profiles across a compositionally graded or
  interdiffused interface, with optional roughness-like jitter

#### Simulation Infrastructure & I/O
- `Simulation::run_streaming()` (`src/builder/mod.rs`): per-step callback execution for
  large or long-running simulations, without materializing the full trajectory in memory
- Binary OVF format read/write (`src/io/ovf.rs`): `Binary4_1_0` (OVF 1.0, big-endian) and
  `Binary4_2_0` / `Binary8_2_0` (OVF 2.0, little-endian), with control-value validation to
  catch byte-order mismatches or truncated files
- `benches/material_benchmark.rs`, `benches/skyrmion_benchmark.rs`: criterion benchmarks
  for material-preset creation and skyrmion-dynamics simulation cost

#### Python Bindings
- `LlbMaterial`, `LlbSolver`, `OnsagerMatrix`, `SpinCaloritronicsMaterial`,
  `batch_rk4_step`, `batch_rk4_multistep` now fully typed and documented in
  `python/spintronics/__init__.pyi`
- New pytest suite under `py/tests/` covering vectors, materials, LLG/LLB solvers, spin
  pumping, spin Hall, caloritronics, and batch SIMD evolution (168 tests)

### Changed

- `ParametricAmplification::degenerate_from_yig` (`src/magnon/nonlinear.rs`): the pump
  field is now the physics-derived degenerate parametric threshold field
  `h_th = α·ω_p/(γ·Ms)` instead of a hard-coded 0.1 mT placeholder; the preset now sits
  exactly at marginal stability (use the new `with_supercriticality` to move off it)
- Added a workspace-level `[workspace.dependencies]` table (COOLJAPAN workspace policy):
  `scirs2-core`, `scirs2-spatial`, `rayon`, `serde`, `serde_json`, `pyo3`, `numpy`,
  `hdf5`, `base64`, the `wasm-bindgen` family, `criterion`, and `proptest` are now
  declared once and consumed via `.workspace = true` from the root crate, `demo/`, and
  `py/`
- Dependency updates: `scirs2-core`/`scirs2-spatial` 0.4.4 → 0.6.0, `pyo3`/`numpy` 0.25 →
  0.29 (Python bindings migrated to the `Bound<'py, T>` API), `rayon` → 1.12, `criterion`
  → 0.8 (benchmarks migrated from `criterion::black_box` to `std::hint::black_box`),
  `proptest` → 1.11, `getrandom` → 0.4
- `spintronics-demo`: `axum` 0.7 → 0.8, `tower-http` 0.6 → 0.7, `askama` 0.12 → 0.16 —
  `askama_axum` is now a deprecated tombstone, so `demo/src/main.rs` bridges templates to
  axum responses manually via a small `render_template` helper
- `TODO.md`: reconciled roughly 22 items that were already shipped in prior releases
  (deferred integrators, spin-wave-theory modes, disorder/defect models, the type-state
  builder, GPU CPU fallback, autodiff RL/diffusion/quantum-classical modules,
  property-based tests, experimental validations) from unchecked to checked; no
  functional changes

### Fixed

- `CMatrix::hermitian_eigendecomposition` (`src/math/matrix.rs`): the Householder-
  tridiagonalization + implicit-QL implementation returned correct eigenvalues but
  subtly incorrect eigenvectors for any non-diagonal Hermitian matrix of size 3 or
  larger (a phase/sign bookkeeping gap between the accumulated unitary and the
  tridiagonal form). Replaced with a cyclic Jacobi algorithm, verified to machine
  precision up to n=16. This transparently corrects every consumer that relies on
  eigenvectors rather than just eigenvalues, including `topomagnon` edge-mode
  localization, 3D band eigenvectors and Wilson-loop centers, `spinwave::magnonic_crystal`,
  `cavity::polariton`, `magnon::spectral`, and the new orbitronics/RVB/hopfion-stability
  modules
- `SpinPumpingSimulation.run()` (Python binding, `src/python/simulation.rs`) did not
  update the tracked magnetization state after evolving; `get_magnetization()` after
  `.run()` now correctly reflects the evolved state, consistent with
  `LlgSimulator.evolve()`
- Cargo lib-name collision between the root `spintronics` crate and the `spintronics-py`
  cdylib: the latter's `[lib] name` is now `spintronics_native` (maturin's
  `module-name = "spintronics"` keeps the Python-facing import name unchanged)
- 21 hardcoded `/tmp/...` test paths replaced with `std::env::temp_dir()` across
  `src/visualization/{csv,vtk,json,hdf5}.rs` and `src/io/ovf.rs`
- `.gitignore`: anchored the `test_*` / `performance_*` rules to the repository root
  (`/test_*`, `/performance_*`) — the previous unanchored globs were silently excluding
  `py/tests/test_*.py` source files from version control
- `py/pyproject.toml` and `py/python/spintronics/__init__.py`: stale Python package
  version `0.5.0` corrected to `0.3.2`, matching the workspace version
- Build hygiene: 7 `autodiff`-feature examples now declare
  `required-features = ["autodiff"]` so `cargo build --examples` succeeds under default
  features; the `data_export_formats` example's `temp_path` helper gained a precise
  `cfg_attr` to avoid a dead-code warning when none of the `vti`/`netcdf`/`zarr` features
  are enabled
- Removed 3 redundant `#[allow(dead_code)]` / `#[allow(unused_imports)]` attributes that
  were not suppressing any active warning (`src/benchmark.rs`, `src/stochastic/thermal.rs`,
  `src/visualization/xdmf.rs`)
- `KaneMeleModel::hamiltonian_at` Rashba block (`src/topomagnon/qsh.rs`): fixed a
  time-reversal-symmetry violation present whenever `lambda_r != 0` (up to ~0.3 residual
  in `Θ·H(k)·Θ⁻¹ = H(−k)`), caused by assembling the spin-flip block from
  `H[B↑,A↓](k) = -conj(H[A↑,B↓](k))` instead of the correct `-H[A↑,B↓](-k)`; this also
  made the Rashba coupling unphysically vanish at the Dirac point K for every `λ_R`. Gap
  closing at the literature Kane-Mele-Rashba critical ratio `λ_R = 2√3·λ_SO` is now
  reproduced at K
- `KaneMeleModel::z2_invariant` (`src/topomagnon/qsh.rs`) for `λ_R != 0`: replaced the
  TRIM-point Pfaffian method — which used the wrong TRIM points for the honeycomb lattice,
  a non-antiunitary time-reversal operator (missing complex conjugation), and an inherent
  gauge dependence that persisted even after both of those were fixed — with a
  gauge-invariant Wilson-loop / hybrid-Wannier-charge-center calculation
  (`src/topomagnon/wilson.rs`, `src/topomagnon/qsh.rs`), validated against the textbook
  Kane-Mele phase diagram and confirmed gauge-invariant under adversarial eigenbasis
  remixing; the Pfaffian method is retained privately only as a fixed-gauge cross-check

## [0.3.1] - 2026-06-10

### Changed

- `DemagField::compute` (`src/micromagnetics/demag.rs`) — inner O(N²) convolution loop
  refactored to use direct flat kernel-array indexing (`kidx = row_base + u_base − jx`),
  eliminating the former Option-returning `get_n_xx/yy/zz` helper dispatch. Results are
  bit-for-bit identical to the previous implementation (verified by a new regression test).
  When the `parallel` feature is enabled, all target cells are computed concurrently via
  `rayon::into_par_iter`, providing multi-core speedup to the demag convolution.
- `BbhModel::hamiltonian_at(kx, ky)` (`src/topomagnon/hoti.rs`) — return type changed from
  `CMatrix` to `Result<CMatrix>`, propagating matrix-dimension errors through `?` instead of
  panicking on shape mismatch. **Breaking change**: callers that previously used the return
  value directly must now handle the `Result` (add `.unwrap()` or `?`-propagate).
- `StandardProblem3` test (`src/validation/standard_problems/sp3.rs`) —
  `test_large_cube_vortex_stable` reduced from 100 to 25 relaxation steps to stay within the
  CI time budget on the O(N²) Newell-tensor demag path; test comments clarified.
- `scirs2-core` updated from `0.3.1` to `0.4.4`; `default-features = false` added with explicit
  `["std", "array", "random", "parallel"]` feature list to avoid pulling in unused optional
  sub-crate dependencies.
- `scirs2-spatial` updated from `0.3.1` to `0.4.4` with `default-features = false`.
- Workspace version tracking unified: top-level `Cargo.toml` and `demo/Cargo.toml` now use
  `version.workspace = true` instead of an explicit version string.

## [0.3.0] - 2026-03-13

### Added

#### Higher-Order Integrators (`src/dynamics/integrators/`)
- Refactored `integrators.rs` into a multi-file module (`mod.rs`, `rhs_fn.rs`, `dormand_prince.rs`, `symplectic.rs`, `semi_implicit.rs`, `adaptive.rs`, `tests.rs`)
- `DormandPrince45`: Embedded 5(4) adaptive Runge-Kutta integrator
- `DormandPrince87`: Embedded 8(7) high-accuracy Runge-Kutta integrator
- `Yoshida4`: Fourth-order symplectic Yoshida integrator for energy-conserving problems
- `ForestRuth`: Forest-Ruth symplectic integrator
- `SemiImplicit`: Semi-implicit integrator for stiff spin systems
- `AdaptiveIntegrator`: Automatic step-size control wrapper for any error-estimating integrator

#### SimulationBuilder Enhancements (`src/builder/mod.rs`)
- `SolverKind` expanded from 3 to 8 variants: `Rk4`, `Euler`, `Heun`, `Dp45`, `Dp87`, `Yoshida4`, `ForestRuth`, `SemiImplicit`
- New builder methods: `solver_dp45()`, `solver_dp87()`, `solver_yoshida4()`, `solver_forest_ruth()`, `solver_semi_implicit()`

#### LLB Equation (`src/dynamics/llb.rs`)
- `LlbMaterial` with `iron()`, `nickel()`, `cofeb()` presets (Curie temperature, exchange stiffness, damping)
- `LlbSolver` with RK4 integration step and full trajectory `run()`
- `LlbResult` trajectory container with time/magnetization history
- Brillouin function and `equilibrium_magnetization(T)` for finite-temperature physics
- Temperature-dependent damping via longitudinal and transverse coefficients

#### Hopfion Dynamics (`src/texture/hopfion_dynamics.rs`)
- `HopfionDynamicsConfig` with physical parameter validation
- `HopfionDynamicsSolver` with 6-point Laplacian exchange, bulk DMI curl, periodic boundary conditions
- Per-site LLG RK4 integration on 3D spin grid
- Hopf invariant calculation via Berry-connection (Whitehead) integral method
- `HopfionDynamicsResult` with trajectory of Hopf invariant, total energy, and magnetization

#### Caloritronics Module (`src/caloritronics/`)
- `OnsagerMatrix` with `yig_pt()`, `fe_pt()`, `cofeb_pt()` material presets
- `HeatCurrentCalculator` computing Fourier, Peltier, and spin-Peltier contributions
- `SpinCaloritronicsMaterial` with unified `compute_all()` entry point
- `CaloritronicsResult` combining all current contributions
- `AllCurrents` struct for structured output

#### SIMD Batch LLG (`src/simd.rs`)
- `batch_add_scaled()`: SIMD-friendly vector accumulation
- `batch_calc_dm_dt()`: Vectorised LLG right-hand side for N spins
- `batch_evolve_rk4()`: Single RK4 step over a batch of N spins
- `batch_evolve_multi_step()`: Multi-step batch evolution with optional normalisation
- Benchmark: `scalar_rk4_N1024` vs `simd_rk4_N1024` in `benches/llg_benchmark.rs`

### Changed
- `src/dynamics/mod.rs` updated to re-export new integrator types and LLB solver
- `src/lib.rs` doc comments updated; test count corrected to 718
- `src/prelude.rs` extended with `LlbMaterial`, `LlbSolver`, `LlbResult`, `HopfionDynamicsConfig`, `HopfionDynamicsSolver`, `HopfionDynamicsResult`, `OnsagerMatrix`, `AllCurrents`, `HeatCurrentCalculator`, `SpinCaloritronicsMaterial`, `CaloritronicsResult`

### Fixed
- Clippy: removed unnecessary `as f64` casts in integrator tests
- Clippy: replaced `for i in 0..len` with `enumerate()` in SIMD tests
- Clippy: replaced `vec![false; N]` with `[false; N]` in parallel sweep tests
- Clippy: added `#[allow(clippy::needless_range_loop)]` where 3-D indices are genuinely required

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

[0.3.2]: https://github.com/cool-japan/spintronics/releases/tag/v0.3.2
[0.3.1]: https://github.com/cool-japan/spintronics/releases/tag/v0.3.1
[0.3.0]: https://github.com/cool-japan/spintronics/releases/tag/v0.3.0
[0.2.0]: https://github.com/cool-japan/spintronics/releases/tag/v0.2.0
[0.1.0]: https://github.com/cool-japan/spintronics/releases/tag/v0.1.0
