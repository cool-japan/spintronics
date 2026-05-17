# TODO List for Spintronics Library

**Version**: 0.6.0 COMPLETE
**Last Updated**: 2026-05-17 - v0.6.0 released
**Status**: 1200 tests passing, ~65K lines (Rust code: ~53K+)

---

## v0.2.0 - COMPLETE (December 2025)

**Release Date**: December 2025
**Status**: Production-ready, 448 tests passing, 0 warnings
**Highlights**: Python bindings, HDF5 export, memory optimization, enhanced API

### Summary of v0.2.0 Deliverables
- PyO3 bindings for core types (PyVector3, PyFerromagnet, PyLlgSimulator, etc.)
- Serialization support (serde) for all public types
- HDF5 export/import (Hdf5Writer, Hdf5Reader)
- Memory pool allocator (VectorPool, SpinArrayPool, Rk4Workspace, HeunWorkspace)
- Enhanced prelude system with module-level preludes
- Display trait, trait hierarchy, Default implementations, builder methods
- Vector3 enhancements (convenience constructors, magnitude_squared, angle_between, project)
- 8 new FM/NM interface materials (Pt, Ta, W combinations)
- Eq + Hash on all enums, inline attributes on hot-path functions
- Unit validation system (14 validators)
- Interactive web demo subcrate (spintronics-demo)
- GitHub Actions CI/CD, CONTRIBUTING.md, CODE_OF_CONDUCT.md

---

## v0.1.0 - COMPLETE (2025)

**Status**: All planned features implemented and tested.

### Summary of v0.1.0 Deliverables
- Core physics effects: SOT, DMI, Edelstein, Spin Nernst, Topological Hall, Rashba
- Solvers: RK4, adaptive time-stepping, Heun stochastic, implicit, SIMD spin chain, parallel multi-domain
- Materials: topological insulators, Weyl semimetals, 2D magnets, multilayers, chiral magnets
- Material database: ferromagnets, antiferromagnets, interfaces
- FEM with Delaunay mesh, WASM, visualization (VTK/CSV/JSON), OOMMF import/export
- 7 examples reproducing experimental/theoretical results

---

## v0.3.0 - COMPLETE (2026-03-13)

**Released**: 2026-03-13
**Theme**: Advanced Physics Modules, Performance, Simulation Infrastructure
**Tests**: 718 passing, 0 failures, 0 warnings
**Code Size**: ~40K total lines (~30K Rust code)

### Priority 1: Advanced Integrators
- [x] Dormand-Prince RK5 (embedded error estimation) — `DormandPrince45`
- [x] Dormand-Prince RK8 (high-order accuracy) — `DormandPrince87`
- [x] Symplectic integrators for energy conservation
  - [x] Velocity Verlet variant for spin dynamics — `Yoshida4`, `ForestRuth`
  - [x] Partitioned Runge-Kutta methods
- [x] Semi-implicit methods for stiff problems — `SemiImplicit`
  - [ ] Implicit midpoint with Newton iteration (deferred)
  - [ ] Crank-Nicolson for diffusion-dominated systems (deferred)
- [ ] Spectral methods for periodic systems (deferred to v0.4.0)

### Priority 2: Spin Wave Theory
- [x] Analytical dispersion relations for thin films — `spinwave` module
- [x] Magnon dispersion and band structure — `spin_wave_dispersion` example
- [ ] Damon-Eshbach modes for in-plane magnetized films (deferred)
- [ ] Backward volume modes for perpendicular magnetization (deferred)
- [ ] Surface spin waves in semi-infinite media (deferred)
- [ ] Spin wave quantization in nanostructures (deferred)
- [ ] Mode decomposition and spectral analysis (deferred)

### Priority 3: Altermagnets
- [x] RuO2 material model and parameters — `altermagnet` module, `altermagnet_ruo2` example
- [x] Spin splitter effect — implemented in `altermagnet`
- [x] Anomalous Hall effect in altermagnets — implemented
- [ ] CrSb material model (deferred)
- [ ] MnTe material model (deferred)
- [ ] Giant magnetoresistance without ferromagnetism (deferred)

### Priority 4: Orbitronics
- [x] Orbital Hall effect — `orbitronics` module
- [x] Orbital torques — implemented
- [x] Orbital-to-spin conversion — implemented
- [ ] d-orbital magnetism (deferred)

### Priority 5: Frustrated Magnets
- [x] Spin ice models — `frustrated` module, `spin_ice_monopoles` example
- [x] Kagome lattice magnets — implemented in `frustrated`
- [ ] Spin liquids (resonating valence bond states) (deferred)
- [ ] Geometric frustration effects on transport (deferred)

### Priority 6: Hopfions
- [x] 3D topological soliton structure — `texture/hopfion_dynamics.rs`
- [x] Hopf index computation — Berry-connection (Whitehead) method
- [x] Current-driven dynamics — per-site LLG RK4 on 3D grid
- [ ] Stability analysis (deferred)

### Priority 7: Magnon BEC & Spin Density Waves
- [x] Magnon Bose-Einstein condensation — `magnon/bec`, `magnon_bec` example
- [ ] Spin density wave formation and dynamics (deferred)
- [ ] Helical magnets and spirals (deferred)
- [ ] Magnon-magnon interactions (deferred)

### Priority 8: Magnetoelastic Coupling (Straintronics)
- [x] Magnetoelastic coupling tensor — `mech/magnetoelastic`
- [x] Strain-induced anisotropy — implemented
- [ ] Piezoelectric control of magnetism (deferred)
- [ ] Surface acoustic wave (SAW) driven spin dynamics (deferred)

### Priority 9: Advanced Disorder & Defects
- [ ] Random anisotropy model (deferred to v0.4.0)
- [ ] Grain boundary effects in polycrystalline films (deferred)
- [ ] Point defects and pinning sites (deferred)
- [ ] Surface roughness modeling (deferred)
- [ ] Inhomogeneous material parameters (graded interfaces) (deferred)

### Priority 10: SIMD Optimization
- [x] Auto-vectorization hints for vector operations — `simd` module
- [x] `batch_add_scaled`, `batch_calc_dm_dt`, `batch_evolve_rk4`, `batch_evolve_multi_step`
- [x] Benchmark SIMD vs scalar — `benches/llg_benchmark.rs`
- [ ] Explicit SIMD with portable_simd (deferred)
- [ ] Target AVX2/AVX-512 for x86_64 (deferred)
- [ ] NEON optimization for ARM (deferred)

### Priority 11: Enhanced Parallel Computing
- [x] Multi-threading with rayon — `parallel` module
- [x] Parallel domain decomposition
- [x] Parallel parameter sweeps
- [ ] Lock-free data structures (deferred)
- [ ] Work-stealing scheduler (deferred)

### Priority 12: Criterion Benchmarks
- [x] Benchmark suite in `benches/` — `llg_benchmark.rs`
- [x] LLG solver performance (scalar vs SIMD batch)
- [ ] Material creation benchmarks (deferred)
- [ ] Skyrmion dynamics benchmarks (deferred)
- [ ] Automated regression alerts (deferred)

### Priority 13: SimulationBuilder
- [x] SimulationBuilder with method chaining — `builder` module
- [x] 8 SolverKind variants: Rk4, Euler, Heun, Dp45, Dp87, Yoshida4, ForestRuth, SemiImplicit
- [x] Preset configurations — `simulation_builder` example
- [ ] Type-state pattern compile-time validation (deferred)
- [ ] Streaming API for large datasets (deferred)

### Spin Caloritronics (New in v0.3.0)
- [x] `OnsagerMatrix` with yig_pt, fe_pt, cofeb_pt presets
- [x] `HeatCurrentCalculator` (Fourier, Peltier, spin-Peltier contributions)
- [x] `SpinCaloritronicsMaterial` with `compute_all()`
- [x] `CaloritronicsResult` and `AllCurrents`

### LLB Equation (New in v0.3.0)
- [x] `LlbMaterial` with iron(), nickel(), cofeb() presets
- [x] `LlbSolver` with RK4 integration and run()
- [x] Brillouin function and equilibrium_magnetization(T)
- [x] Temperature-dependent longitudinal and transverse damping

### Quantum Effects (v0.4.0 — COMPLETE)
- [x] Magnon quantization in confined geometries — `quantum::ZeroPointFluctuations`
- [x] Zero-point fluctuations at T=0 — `ZeroPointFluctuations::zero_point_amplitude`
- [x] Quantum spin Hall effect in 2D TIs (Kane-Mele model — completed in v0.5.0)
- [x] Magnon-photon coupling strength (cavity QED regime) — `TavisCummings`, `MagnonPolariton`

### Non-Equilibrium Transport (v0.4.0 — COMPLETE)
- [x] Non-equilibrium Green's function (NEGF) formalism — `negf::GreenFunction`
- [x] Keldysh formalism for time-dependent transport — `negf::KeldyshSolver`
- [x] Shot noise in spin transport — `negf::ShotNoise`, Fano factor
- [x] Spin accumulation dynamics with diffusion-drift equations — `negf::SpinAccumulation1D`

---

## v0.4.0 - COMPLETE (2026-05-17)

**Released**: 2026-05-17
**Theme**: Quantum Magnonics, NEGF Transport, Topological Bands, Cavity Extensions
**Tests**: 917 passing, 0 failures, 0 warnings
**Code Size**: ~46K total lines (~36K Rust code)

### Math Primitives
- [x] `Complex` — canonical complex number type promoted from `magnon::bec` — `math::complex`
- [x] `CMatrix` — N×N dense complex matrix, Gauss-Jordan inverse, TQLI eigendecomposition — `math::matrix`

### Quantum Magnonics
- [x] `HolsteinPrimakoff` — linear/quadratic HP transform, YIG/AFM presets — `quantum::holstein_primakoff`
- [x] `BogoliubovTransform` — analytical Bogoliubov diagonalization, vacuum occupation — `quantum::bogoliubov`
- [x] `ZeroPointFluctuations` — zero-point amplitude, Casimir free energy — `quantum::zero_point`

### NEGF Non-Equilibrium Transport
- [x] `Hamiltonian1D`, `LeadSelfEnergy`, `SanchoRubio` — tight-binding + leads — `negf::green_function`
- [x] `GreenFunction`, `TransportCalculator` — Landauer transmission, I-V, DOS — `negf::green_function`
- [x] `KeldyshSolver` — Keldysh lesser/greater GFs, non-equilibrium density — `negf::keldysh`
- [x] `ShotNoise` — zero-freq noise, Fano factor, Johnson-Nyquist — `negf::shot_noise`
- [x] `SpinAccumulation1D` — FTCS + implicit (Thomas) spin diffusion — `negf::accumulation`

### Topological Magnon Bands
- [x] `MagnonBandModel` — Haldane honeycomb, Kagome, square-DMI — `topomagnon::band_model`
- [x] `BerryCurvature` — sum-over-states + finite-diff curvature, BZ integration — `topomagnon::berry_curvature`
- [x] `ChernNumber` — Fukui-Hatsugai-Suzuki discrete method, Wilson loop — `topomagnon::chern_number`
- [x] `EdgeModes` — strip diagonalization, IPR localization, chiral velocity — `topomagnon::edge_modes`
- [x] `MagnonHallConductivity` — Matsumoto-Murakami thermal Hall σ_xy — `topomagnon::magnon_hall`

### Cavity Extensions
- [x] `TavisCummings` — g√N collective coupling, Dicke superradiance — `cavity::tavis_cummings`
- [x] `MagnonPolariton` + `MultiModePolariton` — Hopfield diagonalization — `cavity::polariton`
- [x] `BrillouinScattering`, `MicrowaveToOptical`, `MagnonicFrequencyComb` — `cavity::optomagnonic`

### Random Anisotropy
- [x] `RandomAnisotropy` — Imry-Ma, Harris, Marsaglia axes, LLG field — `material::random_anisotropy`

### Examples (6 new)
- [x] `magnon_zero_point`, `negf_transport`, `tavis_cummings_dicke`
- [x] `magnon_polariton`, `topological_magnon_haldane`, `random_anisotropy_disorder`

## v0.5.0 - COMPLETE (2026-05-17)

**Released**: 2026-05-17
**Theme**: Non-Collinear Magnetism, Multiferroics, Topological QSH, Nonlinear Magnons
**Tests**: 996 passing, 0 failures, 0 warnings
**Code Size**: ~52K total lines (~42K Rust code)

### Non-Collinear Magnetism (`src/noncollinear/`)
- [x] Spin spirals: cycloidal, helical, conical, fan structure — `SpinSpiral`
- [x] Luttinger-Tisza ground-state search — `LuttingerTisza`
- [x] Exchange Fourier transform J(q), frustration ratio, ordering temperature
- [x] TbMnO₃ preset, J1-J2 chain, ferromagnet, antiferromagnet presets

### Multiferroics / Magnetoelectric Coupling (`src/multiferroic/`)
- [x] Linear ME tensor α_ij with Dzyaloshinskii bound — `MagnetoelectricTensor`
- [x] Presets: BiFeO₃, TbMnO₃, Cr₂O₃
- [x] KNB mechanism P ∝ e_ij × (S_i × S_j) — `KnbMechanism`
- [x] Inverse ME: E-field control of magnetisation — `InverseMagnetoelectric`
- [x] Free functions: DM polarization, exchange striction, toroidal moment

### Quantum Spin Hall / Kane-Mele (`src/topomagnon/qsh.rs`)
- [x] Full 4-band Bloch Hamiltonian — `KaneMeleModel`
- [x] Z2 invariant via Fukui-Hatsugai on honeycomb BZ parallelogram
- [x] Rashba and staggered potential phase boundaries
- [x] Helical edge states in strip geometry

### Nonlinear Magnon Physics (`src/magnon/nonlinear.rs`)
- [x] Four-magnon scattering vertex T_kk — `FourMagnonScattering`
- [x] Suhl instability threshold and parametric growth rate
- [x] Parametric amplification gain — `ParametricAmplification`
- [x] Nonlinear FMR linewidth and bistability — `NonlinearFmrLinewidth`

### Examples (4 new)
- [x] `spin_spiral_tbmno3.rs` — LT ground-state + KNB polarization
- [x] `bife_o3_multiferroic.rs` — ME coupling + DM + toroidal moments
- [x] `kane_mele_qsh.rs` — Z2 phase diagram + edge states
- [x] `nonlinear_magnon_suhl.rs` — Suhl instability + parametric amp

---

## v0.6.0 - COMPLETE (2026-05-17)

**Released**: 2026-05-17
**Theme**: Spin Wave Extensions, HOTI, Axion Electrodynamics, Data Formats, ML Autodiff
**Tests**: 1200 passing, 0 failures, 0 warnings
**Code Size**: ~65K total lines (~53K Rust code)

### Spin Wave Theory Extensions
- [x] `DamonEshbachDetailed` — full DE dispersion, non-reciprocity, surface localization (src/spinwave/damon_eshbach.rs)
- [x] `BackwardVolumeMSW` — BVMSW dispersion, negative group velocity, crossover wavevector (src/spinwave/bvmsw.rs)
- [x] `SurfaceSpinWave` — semi-infinite medium, Rado-Weertman boundary condition (src/spinwave/surface.rs)
- [x] `SpectralMagnonSolver` — FFT/CMatrix eigenmodes, DOS, spectral weight, mode decomposition (src/magnon/spectral.rs)

### Advanced Topological Phenomena
- [x] `WilsonLoop` — multi-band Wilson loop, link matrices, Wannier centers, nested polarization
- [x] `BbhModel` — BBH 4-band HOTI, quadrupole moment, topological corner states
- [x] `BreathingKagomeModel` — 3-band kagome, corner Z₃ polarization
- [x] `CornerStateSolver` — OBC×OBC finite-cluster corner state spectrum + IPR localization
- [x] `MagnonBandModel3D` — 3D cubic Haldane + pyrochlore presets
- [x] `AxionElectrodynamics` — 3D Berry-curvature θ-term (Chern-Simons form), α_TME, axion response
- [x] `AxionMagnonPhoton` — axion-mediated magnon-photon coupling, Faraday rotation, cooperativity

### Data Export Formats
- [x] `VtiWriter` — VTK ImageData (XML + base64 binary), feature `vti` (src/visualization/vti.rs)
- [x] `XdmfWriter` — XDMF v3.0 (XML + raw f64 binary), feature `xdmf` (src/visualization/xdmf.rs)
- [x] `NetCdfWriter`/`NetCdfReader` — pure-Rust NetCDF3 Classic (XDR binary), feature `netcdf` (src/visualization/netcdf.rs)
- [x] `ZarrStore`/`ZarrArray` — pure-Rust Zarr v2 (JSON+binary chunks), feature `zarr` (src/visualization/zarr.rs)

### ML Autodiff (`#[cfg(feature = "autodiff")]`)
- [x] `Tape` / `Var<'t>` — reverse-mode AD tape; arithmetic ops + sin/cos/exp/ln/sqrt/tanh/powi/powf
- [x] `Sgd` (with momentum), `Adam` (Kingma & Ba 2014), `LBfgs` (two-loop BFGS)
- [x] `ParameterFitter` — closure-based gradient fitting; `FitResult`
- [x] Differentiable physics: `kittel_frequency_diff`, `zeeman_energy_diff`, `exchange_energy_diff`, `dmi_energy_diff`, `anisotropy_energy_diff`

### Examples (6 new, total 41)
- [x] `damon_eshbach_nonreciprocity.rs`, `backward_volume_magnons.rs`, `hoti_corner_states.rs`
- [x] `axion_magnon_photon.rs`, `data_export_formats.rs`, `autodiff_parameter_fitting.rs`

## v0.7.0 - ROADMAP

**Target**: Q4 2026
**Theme**: Language Bindings, Neural Network Potentials, Extended ML, Experimental Validation

### Language Bindings
- [ ] Julia bindings via julia-rs
- [ ] C/C++ bindings via cbindgen
- [ ] R bindings for statistical analysis

### ML Enhancements (Phase 2)
- [ ] Neural network potentials for exchange interactions (trainable Heisenberg/DMI)
- [ ] ML-based surrogate models for fast parameter space exploration
- [ ] Gradient-based structure optimization (LLG + autodiff)
- [ ] Physics-informed neural networks for spin dynamics

### Advanced Spin Wave Theory
- [ ] Damon-Eshbach in semi-infinite ferromagnets (single-surface non-reciprocity)
- [ ] Experimental validation against Demidov et al. 2006 BLS data
- [ ] Spin wave quantization in nanodisks (combined Bessel + Kalinikos-Slavin)
- [ ] Magnonic crystal band gaps

### Workspace Restructuring (candidate for v1.0.0)
- [ ] Split into workspace with multiple crates:
  - `spintronics-core` — Core physics and materials
  - `spintronics-solver` — Numerical solvers
  - `spintronics-io` — I/O and visualization
  - `spintronics-python` — Python bindings
  - `spintronics-cli` — Command-line tools

---

## v1.0.0 - STABLE RELEASE

**Target**: 2027
**Theme**: API Stabilization, Production-Grade, Full Ecosystem

### API Stabilization
- [ ] Comprehensive API review and stabilization
- [ ] Migration guides for all breaking changes
- [ ] Semantic versioning policy document
- [ ] Backward compatibility testing with semver-check

### GPU Acceleration (Feature-Gated, Not Default)
- [ ] CUDA backend for LLG solver (feature = "cuda")
  - [ ] GPU-accelerated RK4 integration
  - [ ] Parallel magnetization dynamics for large systems (>1M spins)
  - [ ] Benchmark against CPU (target: 10-100x speedup)
- [ ] ROCm support for AMD GPUs (feature = "rocm")
- [ ] Fallback to CPU for systems without GPU
- [ ] Unified API: transparent GPU/CPU selection

### MPI Distributed Computing (Feature-Gated, Not Default)
- [ ] MPI support for distributed computing (feature = "mpi")
  - [ ] Domain decomposition across nodes
  - [ ] Halo exchange for boundary communication
  - [ ] Parallel I/O with parallel HDF5
- [ ] Hybrid MPI + GPU for supercomputer deployment

### Production Features
- [ ] Comprehensive documentation website (mdBook)
- [ ] Multi-language bindings (Python, Julia, C/C++, R, MATLAB)
- [ ] Extensive experimental validation (NIST Standard Problems 1-5)
- [ ] Performance competitive with OOMMF, mumax3

---

## Known Issues & Technical Debt

### Current Known Issues
- [ ] Thermal noise implementation needs validation against experiments
  - Compare with Einstein relation for FMR linewidth
- [ ] Stochastic solver convergence for very small damping (alpha < 0.001)
- [ ] Edge cases in skyrmion number calculation near boundaries
- [ ] Numerical stability in strong exchange limit
- [x] Unit consistency checks across modules (units.rs with 14 validators)

### Technical Debt

#### High Priority
- [ ] Refactor Vector3 to use scirs2-linalg (Breaking change candidate)
- [ ] Profile-guided optimization infrastructure

#### Medium Priority
- [ ] Consolidate similar code patterns across effect modules
- [ ] Clean up redundant type conversions
- [ ] Improve naming consistency across public APIs

#### Low Priority
- [ ] Consider const generics for compile-time dimensions
- [ ] Explore zero-copy serialization with rkyv
- [ ] Split large modules into submodules if they grow

---

## Testing & Quality Infrastructure

### Testing
- [ ] Code coverage with tarpaulin or cargo-llvm-cov (target: >80%)
- [ ] Property-based testing with proptest (conservation laws, symmetries)
- [ ] Fuzzing with cargo-fuzz (numerical solvers, serialization)
- [ ] Mutation testing with cargo-mutants

### Documentation
- [ ] API documentation improvements (cross-linking, more inline examples)
- [ ] Documentation testing in CI (compile examples, check links)
- [ ] Changelog automation with git-cliff

### Release Management
- [ ] Release checklist automation
- [ ] Automated releases to crates.io (tag-based)
- [ ] Backward compatibility testing (semver-check)

---

## Research & Validation

### Experimental Validation
- [ ] Spin Pumping: Validate against Saitoh et al. 2006 APL
- [ ] Spin Seebeck Effect: Compare with Uchida et al. 2008 Nature
- [ ] Skyrmion Size: Validate with Woo 2016 Nat. Mater.
- [ ] Magnon Dispersion: Check against neutron scattering data
- [ ] Thermal Transport: Validate magnon thermal conductivity in YIG
- [ ] Micromagnetic Benchmarks: Cross-check with OOMMF, mumax3

### Literature Review
- [ ] Spin-orbit torques survey (2020-2026)
- [ ] Cavity magnonics latest experiments
- [ ] Antiferromagnetic spintronics developments
- [ ] 2D materials spin-orbitronics
- [ ] Topological magnonics
- [ ] Magnon-based computing (neuromorphic, reservoir)

---

## Documentation & Community

### Tutorials
- [ ] Tutorial 1: Introduction to Spintronics Simulations (YIG/Pt spin pumping)
- [ ] Tutorial 2: Skyrmion Physics (creation, manipulation, topological Hall)
- [ ] Tutorial 3: Thermal Spintronics (Seebeck, Nernst, magnon transport)
- [ ] Tutorial 4: Advanced Topics (AFM THz, cavity magnonics, reservoir computing)

### Community
- [ ] Set up GitHub Discussions forum
- [ ] Write JOSS paper
- [ ] Present at conferences (MMM, Intermag, APS March Meeting)
- [ ] Submit to This Week in Rust

---

## Release Schedule

| Version | Target | Theme | Status |
|---------|--------|-------|--------|
| v0.1.0 | 2025 | Core Physics & Materials | COMPLETE |
| v0.2.0 | Dec 2025 | Python Bindings, HDF5, Memory Optimization | COMPLETE |
| v0.3.0 | 2026-03-13 | Advanced Physics, Performance, Simulation Infrastructure | COMPLETE |
| v0.4.0 | Q4 2026 | Research Features, ML, Ecosystem Expansion | Planned |
| v1.0.0 | 2027 | API Stabilization, Production-Grade | Planned |

---

## Success Metrics

### v0.3.0 Goals (ACHIEVED)
- [x] **Tests**: 718 tests passing (target was 570+)
- [x] **Code Size**: ~40K total lines / ~30K Rust code (target was ~29,500 lines)
- [x] **Performance**: SIMD batch LLG + parallel lattice evolution implemented
- [x] **Examples**: 25 examples (target achieved)
- [x] **Quality**: 0 warnings, 0 unwrap() in production code

### Long-Term Vision (v1.0.0)
- [ ] De facto standard for spintronics simulation in Rust
- [ ] Performance competitive with OOMMF, mumax3
- [ ] Used by multiple research groups
- [ ] Multiple papers citing the library
- [ ] Used in university courses

### Community Targets
- GitHub Stars: 100+ by Q2 2026
- Contributors: 5+ active
- Papers citing: 3+ by end 2026
- Downloads: 2000+ from crates.io

---

## Notes for Contributors

### Development Philosophy
- **Physics First**: Validate against experiments and physical intuition
- **Type Safety**: Use Rust's type system to prevent unphysical states
- **Pure Rust**: Default features must be 100% Pure Rust (C/Fortran deps feature-gated only)
- **Performance**: Profile before optimizing; correctness > speed
- **Documentation**: Every public function with doc comments and physics context
- **Testing**: Tests that verify physical behavior, not just code coverage
- **References**: Cite papers in code comments for implemented equations

### Code Standards
- Formatting: rustfmt.toml (enforced in CI)
- Linting: clippy with warnings as errors
- Testing: >80% code coverage target
- Documentation: All public APIs documented
- No unwrap() usage
- Snake_case naming convention

### Getting Started
1. Read CONTRIBUTING.md
2. Check Issues labeled "good first issue"
3. Join community discussions
4. Start with documentation or test additions
5. Gradually move to feature implementation

---

**Maintained by**: COOLJAPAN OU (Team KitaSan)
**License**: MIT OR Apache-2.0
**Repository**: https://github.com/cool-japan/spintronics
**Contact**: See CONTRIBUTING.md for communication channels
