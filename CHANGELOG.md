# Changelog

All notable changes to the spintronics library will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.9.0] - 2026-05-17

### Added

#### Advanced ML Phase 4 (`#[cfg(feature = "autodiff")]`)

- `LatticeGraph`, `GraphMessagePassingLayer`, `GraphMlp`, `NodeFeatures` — graph neural network potentials for spin Hamiltonians on arbitrary lattice topology, building on v0.8.0 equivariant layers. Each message-passing layer uses `EquivariantLinear` for both the message function `f_msg(s_i, s_j, v_i, v_j, |r_ij|, r̂_ij)` and the update function `f_update(s_i, m_i_s, v_i, m_i_v)`; aggregation by summation preserves SO(3) equivariance; rotation invariance verified to 2.8×10⁻¹⁴ (machine precision). Convenience graph builders: `chain_1d`, `ring_1d`, `square_lattice_2d` (`src/autodiff/graph_nn.rs`)
- `GaussianProcess`, `GpConfig`, `AcquisitionStrategy`, `BayesianOptConfig`, `BayesianOptResult`, `BayesianOptimizer` — Bayesian optimization with RBF-kernel Gaussian Process surrogate (Cholesky + jitter retry) and Expected Improvement / UCB / Posterior Variance acquisition strategies. Self-contained `erf` via Abramowitz–Stegun 7.1.26 for Φ(z); converges on `(x-0.7)²` to within 0.05 (`src/autodiff/bayesian_opt.rs`)

#### More Experimental Validations

- `Garello2013Validation` — Garello et al. Nat. Nanotechnol. 8, 587 (2013): angular-harmonic decomposition of SOT in Pt/Co/AlOx, H_DL ≈ 0.4 mT per 10¹¹ A/m², H_FL/H_DL ≈ 0.3, critical J_c ≈ 5×10¹⁰ A/m² (`src/validation/experimental/garello_2013.rs`)
- `Boona2014Validation` — Boona, Heremans, Goennenwein et al. MRS Bulletin 39, 426 (2014): LSSE in granular YIG/Pt with Seebeck coefficient S_S ≈ 0.5–1.0 μV/K, linear V_LSSE(ΔT) up to 50 K, grain-boundary attenuation factor (`src/validation/experimental/boona_2014.rs`)
- **Total: 7 landmark-paper validations** spanning v0.7.0–v0.9.0 (Demidov 2006, Saitoh 2006, Uchida 2008, Mosendz 2010, Liu 2012, Garello 2013, Boona 2014)

#### GPU Acceleration Skeleton (feature `cuda`)

- `Device` trait — object-safe abstraction for portable LLG simulation: `name`, `is_available`, `max_spins`, `step_llg_rk4`, `step_llg_rk4_multi`, `zeeman_energy` (`src/gpu/mod.rs`)
- `CpuDevice` — always available; wraps the scalar CPU path via `crate::llg::calc_dm_dt` with frozen-field RK4 and per-spin renormalization. Demonstrated at ~92 ns / (spin · step) for 100–1000-spin batches (`src/gpu/cpu.rs`)
- `CudaDevice` (`#[cfg(feature = "cuda")]`) — v0.9.0 skeleton. Constructs successfully with `available=false`; all operations return `numerical_error("CUDA backend: not implemented in v0.9.0 ...")`. Establishes the type signatures and contracts that v1.0.0 CUDA kernels must honor (`src/gpu/cuda.rs`)
- `available_devices()`, `select_best_device()` — device enumeration + auto-selection helpers
- Feature `cuda = []` (no external deps yet — pure skeleton plumbing)

#### Examples (4 new, total 56)

- `graph_nn_lattice.rs` — 9-spin lattice graph with EquivariantLinear message passing; rotation invariance verified to 2.8×10⁻¹⁴
- `bayesian_opt_materials.rs` — BO finds (J*, K*) optimum on 625-point grid within 0.026 in 20 evaluations (vs random baseline)
- `validation_full_suite.rs` — runs all 7 landmark validations; 18/23 quantitative checks pass within 30% relative-error window
- `gpu_device_demo.rs` — enumerates available devices, evolves 100 spins × 200 steps in 1.8 ms (|m| drift 2.2×10⁻¹⁶), scaling sweep over 10–1000 spins

### Changed

- `prelude` re-exports v0.9.0 types: `LatticeGraph`, `GraphMessagePassingLayer`, `GraphMlp`, `NodeFeatures`, `GaussianProcess`, `GpConfig`, `AcquisitionStrategy`, `BayesianOptConfig`, `BayesianOptResult`, `BayesianOptimizer`, `Device`, `CpuDevice`, `CudaDevice` (feature-gated), `available_devices`, `select_best_device`, `Garello2013Validation`, `Boona2014Validation`

### Test Counts

- **Library tests**: 1463 / 1463 passing (was 1403 in v0.8.0; +60 new across graph_nn, bayesian_opt, garello_2013, boona_2014, gpu)
- **Integration property tests**: 27 / 27 passing (unchanged)
- **Doctests**: 103 / 103 passing (5 ignored)
- **Zero clippy warnings** with `-D warnings` on `--all-targets` and full feature set (incl. `cuda`)
- **4 new examples** all run with physically sensible output

## [0.8.0] - 2026-05-17

### Added

#### Stochastic Methods Improvements (`#[cfg(feature = "scirs2")]`)

- `HeunAdaptive` — Heun-Euler embedded pair (orders 1 and 2) with PI-style step controller `dt_new = dt · clamp(safety · (tol/err)^(1/2), min_factor, max_factor)`. Frozen-noise convention preserves Brownian statistics across step rejections via `H_th · √dt_old/√dt_new` rescaling (`src/stochastic/heun_adaptive.rs`)
- `ImplicitMilstein` — implicit-midpoint drift + explicit Wiener-coupled stochastic term with Newton-Raphson 3×3 finite-difference Jacobian and self-contained Gauss elimination. Optional Milstein correction `Δm = ½(−γ/(1+α²))² · (m × σΔW) × σΔW` captures leading off-diagonal precessional contribution; disable to recover implicit Euler-Maruyama (`src/stochastic/implicit_milstein.rs`)
- `PimcSimulation`, `PimcConfig`, `PimcLattice` (Chain1D / Ring1D), `PimcResult` — worldline-based path-integral Monte Carlo for finite-T Heisenberg chains with Trotter-discretized imaginary time, periodic-in-τ rigidity coupling K_rig = 1/Δτ, single-spin Metropolis updates with Marsaglia random unit-vector proposals; estimators: `⟨M_z⟩`, per-site energy, β² · var(E) · N specific heat, β · var(M_z) · N susceptibility (`src/stochastic/pimc.rs`)

#### Advanced ML Phase 3 (`#[cfg(feature = "autodiff")]`)

- `EquivariantLinear`, `EquivariantMlp`, `EquivariantConfig` — Cartesian-tensor O(3)-equivariant NN layers: scalar channels (rotation-invariant) + vector channels (rotation-equivariant) + scalar-vector mixing via `s_j · v_j` (equivariant) and `|v|` (invariant); gated tanh activation preserves equivariance; rotation invariance verified to machine precision (`src/autodiff/equivariant.rs`)
- `random_so3(seed)` — Marsaglia-normal axis + Rodrigues rotation matrix; orthonormality verified to 1e-12
- `rotate_vector(R, v)` — 3×3 matrix-vector product
- `ActiveLearner`, `ActiveLearningConfig`, `QueryStrategy` (UncertaintySampling / QueryByCommittee / RandomBaseline), `ActiveLearnResult` — active learning loop with ensemble bootstrap, predict_with_uncertainty, select_next_query (argmax sum-of-stds), full `fit(oracle, pool)` (`src/autodiff/active_learning.rs`)

#### More Experimental Validations

- `Mosendz2010Validation` — Mosendz et al. PRL 104, 046601 (2010): V_ISHE vs Pt thickness scaling, linewidth enhancement Δα_eff = (γℏ g↑↓_eff)/(4π·M_s·t_F), spin mixing conductance g↑↓ ≈ 2.1×10¹⁹ m⁻², refined θ_SH(Pt) ≈ 0.013 (`src/validation/experimental/mosendz_2010.rs`)
- `Liu2012Validation` — Liu et al. Science 336, 555 (2012): β-Ta spin Hall angle θ_SH ≈ -0.12 (opposite sign to Pt), critical SOT switching current J_c, thickness scaling, polarity convention (`src/validation/experimental/liu_2012.rs`)

#### Property-Based Testing (first `tests/` integration suite)

- `tests/property_conservation.rs` — 12 property tests with 32 cases each: |m|=1 conservation under RK4/Heun/Euler, Zeeman energy conservation at α=0, damping aligns m with H, Larmor reverses under field flip, dm/dt perpendicular to m and h at α=0, cross-product anti-commutativity, triple product cyclic, Lagrange identity, normalize idempotent
- `tests/property_symmetries.rs` — 15 property tests with 32 cases each: rotation matrix orthogonality + determinant, SO(3) invariants (dot product, magnitude), cross-product equivariance, exchange / Zeeman / uniaxial-anisotropy invariance under co-rotation, calc_dm_dt equivariance, time-reversal at α=0, linearity in H, parity, rotation composition + inverse
- Total: **27 property tests × 32 random cases = 864 randomized property trials per `cargo test` run**
- `proptest = "1.6"` added to `[dev-dependencies]`

#### Examples (5 new, total 52)

- `heun_adaptive_thermal_llg.rs` — Permalloy thin film at T = 300 K with `HeunAdaptive`; dt auto-adjusts (96–276 fs) as magnetization warms up
- `pimc_heisenberg_chain.rs` — 1D Heisenberg ring at varying β: ⟨M_z⟩ transitions from paramagnetic (β=0.05) to aligned (β=5.0)
- `equivariant_nn_demo.rs` — 4-spin EquivariantMlp; rotation invariance preserved to 4.4×10⁻¹⁶ (machine precision) over 20 random SO(3) rotations
- `active_learning_demo.rs` — `f(x) = sin(5x)·exp(-x²)` with 20-query budget: RandomBaseline MSE = 0.184; QueryByCommittee MSE = 0.0075 (24× better); UncertaintySampling MSE = 0.0161 (11× better)
- `validation_landmark_suite.rs` — runs all 5 landmark validations (3 from v0.7.0 + 2 new); reports 12/15 quantitative checks pass within 30%-relative-error window

### Changed

- `prelude` now re-exports v0.8.0 types: `HeunAdaptive`, `ImplicitMilstein`, `PimcSimulation`, `PimcConfig`, `PimcLattice`, `PimcResult`, `EquivariantLinear`, `EquivariantMlp`, `EquivariantConfig`, `ActiveLearner`, `ActiveLearningConfig`, `QueryStrategy`, `ActiveLearnResult`, `random_so3`, `rotate_vector`, `Mosendz2010Validation`, `Liu2012Validation`
- Refactored two `assertions_on_constants` clippy warnings (introduced during Phase 3) at the source: const-only assertions moved to module-level `const _: ()` compile-time checks; runtime assertions retained for non-const-evaluable comparisons

### Test Counts

- **Library tests**: 1403 / 1403 passing (was 1329 in v0.7.0; +74 new across stochastic, autodiff, and validation modules)
- **Integration property tests**: 27 / 27 passing (12 conservation + 15 symmetries, 32 cases each, **first `tests/` directory in repo**)
- **Doctests**: 103 / 103 passing
- **Zero clippy warnings** with `-D warnings` on `--all-targets` and full feature set
- **5 new examples** all run with physically sensible output

## [0.7.0] - 2026-05-17

### Added

#### ML Enhancements (Phase 2) — autodiff feature

- `Mlp`, `Layer`, `Activation` (`Relu`, `Tanh`, `Sigmoid`, `Gelu`, `Linear`) — feed-forward neural network with Xavier/He init, full `Var<'t>` tape integration for backprop; LCG RNG for reproducibility (`src/autodiff/neural.rs`)
- `NeuralExchange` — trainable Heisenberg J(r) potential, rescales r to [-1,1] for stable training; `coupling()` for fast inference, `coupling_diff()` for differentiable forward
- `NeuralAnisotropy` — trainable single-ion anisotropy K(m_x, m_y, m_z)
- `LlgPinn`, `PinnTrainer` — physics-informed neural network for LLG dynamics: time derivative via central finite difference on the tape; explicit LL form `dm/dt = -γ/(1+α²)[m×H + α m×(m×H)]`; configurable residual / IC / norm penalty weights (`src/autodiff/pinn.rs`)
- `SpinConfig` (spherical coords, unit-norm built-in), `EnergyFunctional` (Heisenberg + uniaxial anisotropy + Zeeman with user/chain/ring bond presets), `MagneticStructureOptimizer`, `StructureOptResult`, `find_fm_ground_state`, `find_afm_ground_state` (`src/autodiff/structure_opt.rs`)

#### Advanced Spin Wave Theory

- `NanodiskSpinWaves` — confined modes in magnetic disks: radial (Bessel zero) × azimuthal angular quantization, `mode_frequency(n_rad, m_az)`, `mode_profile`, `mode_spectrum(n_max, m_max)`, `group_velocity`, `propagation_length`; presets `yig_nanodisk(R)`, `permalloy_nanodisk(R)`, `cofeb_nanodisk(R)` (`src/spinwave/nanodisk.rs`)
- `MagnonicCrystal1D` / `MagnonicCrystal2D` — periodic magnetic media via plane-wave expansion: `hamiltonian_at(k)` builds Hermitian Hamiltonian in reciprocal-lattice basis, `band_structure(n_kpoints)`, `band_gap(band_idx)`, `group_velocity(band_idx, k)`; constraint n_pw ≤ 31 (1D) / 3 (2D) from `CMatrix::MAX_DIM=64`; presets `yig_pt_alternating(period)`, `nife_cofe_alternating(period)`, `MagnonicCrystal2D::checkerboard()` (`src/spinwave/magnonic_crystal.rs`)
- `SemiInfiniteDamonEshbach` — single-surface DE in semi-infinite media (distinct from thin-film DE and `SurfaceSpinWave`): ω² = (ω_H + Dk²)(ω_H + ω_M + Dk²) + (ω_M/2)², canonical √[ω_H(ω_H+ω_M)] + ω_M/2 in k→0 limit; surface localization, non-reciprocity, group velocity, propagation length; presets `yig_bulk()`, `iron_bulk()`, `permalloy_bulk()` (`src/spinwave/semi_infinite_de.rs`)

#### Stiff/Diffusion Integrators

- `ImplicitMidpointNewton` — A-stable 2nd-order implicit midpoint with Newton-Raphson + finite-difference Jacobian; explicit Euler initial guess; self-contained Gauss elimination with partial pivoting; `with_max_iter`, `with_tol`, `with_fd_step`; implements existing `Integrator` trait for drop-in replacement (`src/dynamics/integrators/implicit_midpoint.rs`)
- `CrankNicolsonDiffusion` — unconditionally stable Crank-Nicolson for 1D diffusion equations: Dirichlet / Neumann (zero-flux) / Periodic boundaries; Thomas algorithm for tridiagonal (Sherman-Morrison for periodic); `step`, `evolve`, `steady_state` methods (`src/dynamics/integrators/crank_nicolson.rs`)
- `SpinDiffusionCrankNicolson` — specialized for ∂μ_s/∂t = D∇²μ_s − μ_s/τ_sf with implicit-trapezoidal sink; `spin_diffusion_length() = √(D τ_sf)`

#### Experimental Validation Refactor

- Refactored `src/validation.rs` → `src/validation/` directory hierarchy (backward-compatible via `pub use parameter_checks::*` in `mod.rs`)
- `src/validation/parameter_checks.rs` — 11 existing `check_*` / `is_valid_*` functions preserved verbatim
- `src/validation/experimental/mod.rs` — `ValidationResult { name, max_relative_error, mean_relative_error, n_points, tolerance, passed }`
- `Demidov2006Validation` — Damon-Eshbach dispersion and non-reciprocity vs BLS data from Demidov et al. PRL 96, 097202 (2006); embedded `K_VALUES`, `OMEGA_GHZ`, `NONRECIP_PCT` const arrays (`src/validation/experimental/demidov_2006.rs`)
- `Saitoh2006Validation` — Pt spin Hall angle, ISHE polarity, linear J_s scaling vs Saitoh et al. APL 88, 182509 (2006); literature range 0.0037–0.013 (`src/validation/experimental/saitoh_2006.rs`)
- `Uchida2008Validation` — Longitudinal Spin Seebeck Effect: linear V_LSSE(ΔT), polarity, order-of-magnitude check vs Uchida et al. Nature 455, 778 (2008) (`src/validation/experimental/uchida_2008.rs`)

#### Examples (6 new, total 47)

- `neural_exchange_training.rs` — Train `NeuralExchange` on synthetic RKKY-like target with Adam + finite-difference gradient (~3× MSE reduction)
- `nanodisk_modes.rs` — YIG 100 nm disk: spectrum table, lowest mode at ~2.2 GHz, radial profile (vanishes at boundary), material comparison
- `implicit_midpoint_stiff_demo.rs` — Stiff LLG (H=10 MA/m, α=0.1): implicit midpoint at dt=T_L/30 remains bounded; explicit Euler diverges; 2nd-order convergence verified to machine precision
- `pinn_llg_solver.rs` — Train `LlgPinn` on Larmor precession, compare against RK4 reference; demonstrates norm constraint conservation
- `magnonic_crystal_bandgap.rs` — 1D NiFe/CoFeB alternating crystal (200 nm period): band structure, first band gap (~57% relative width), group velocity vanishes at BZ edges
- `experimental_validation_demo.rs` — Run all 3 landmark validations: 6/8 quantitative checks pass with 30%-relative-error tolerance

### Changed

- `prelude` now re-exports v0.7.0 types: `Mlp`, `NeuralExchange`, `NeuralAnisotropy`, `LlgPinn`, `PinnTrainer`, `MagneticStructureOptimizer`, `SpinConfig`, `EnergyFunctional`, `NanodiskSpinWaves`, `MagnonicCrystal1D`, `MagnonicCrystal2D`, `SemiInfiniteDamonEshbach`, `ImplicitMidpointNewton`, `CrankNicolsonDiffusion`, `DiffusionBoundary` (alias to avoid clash with `negf::BoundaryCondition`), `SpinDiffusionCrankNicolson`, `Demidov2006Validation`, `Saitoh2006Validation`, `Uchida2008Validation`, `ExperimentalValidationResult`
- `validation` is now a directory module with `parameter_checks` and `experimental` submodules; all `crate::validation::check_*` calls remain backward-compatible
- Fixed pre-existing doctest issues in `src/multiferroic/mod.rs` (bond vector parallel to spin current → zero polarization) and `src/noncollinear/mod.rs` (tolerance too strict for J1-J2 chain spiral detection)

### Test Counts

- **Library tests**: 1329 / 1329 passing (was 1200 in v0.6.0; +129 new)
- **Doctests**: 103 passing (2 pre-existing failures fixed)
- **Zero clippy warnings** with `-D warnings` on `--all-targets` and full feature set
- **6 new examples** all run end-to-end with physically sensible output

## [0.6.0] - 2026-05-17

### Added

#### Spin Wave Theory Extensions

- `DamonEshbachDetailed` — full Kalinikos-Slavin in-plane thin-film dispersion with chiral surface localization, non-reciprocity Δω=ω(+k)−ω(−k), surface localization length, propagation length; presets `yig_film_micron()`, `permalloy_thin()`, `cofeb_strip()` (`src/spinwave/damon_eshbach.rs`)
- `BackwardVolumeMSW` — perpendicular-magnetization BVMSW dispersion, **negative group velocity** below crossover k*, bisection-based `crossover_wavevector()`, cosine thickness mode profiles; presets `yig_yagi(h_perp)`, `permalloy_perp(h_perp)` (`src/spinwave/bvmsw.rs`)
- `SurfaceSpinWave` — semi-infinite medium spin wave: `dispersion_omega(k)`, exchange-corrected penetration depth, Rado-Weertman surface anisotropy shift, exponential field amplitude; presets `bulk_yig()`, `bulk_iron()` (`src/spinwave/surface.rs`)
- `SpectralMagnonSolver` — FFT-based magnon eigenmodes via `CMatrix::hermitian_eigendecomposition`; Bloch-resolved bands for square FM and honeycomb AFM; Lorentzian DOS, spectral weight A(ω,k), mode decomposition (`src/magnon/spectral.rs`, `#[cfg(feature = "scirs2")]`)
- Extended `SpinWaveMode` enum with `SurfaceLocalized` variant and corresponding `mode_profile` dispatch

#### Advanced Topological Physics

- `WilsonLoop` — multi-band Wilson loop operator W_y(kx): `link_matrix` (overlap matrix M_mn=⟨u_m^a|u_n^b⟩), `wilson_unitary_y`, `wannier_centers` (eigenphases of W_y), `nested_polarization`, `polarization_p_x/y`; foundation for HOTI Z₂ invariant (`src/topomagnon/wilson.rs`)
- `BbhModel` — 4-band Benalcazar-Bernevig-Hughes model: `hamiltonian_at(kx,ky)` with anticommuting Γ matrices, `energy_bands`, `band_gap`, `is_higher_order_topological`, `quadrupole_moment` via nested Wilson loop; presets `topological_phase()`, `trivial_phase()` (`src/topomagnon/hoti.rs`)
- `BreathingKagomeModel` — 3-band kagome with topological corner Z₃: `hamiltonian_at`, `corner_polarization` via Wannier centers; presets `topological_kagome()`, `trivial_kagome()`
- `CornerStateSolver` — OBC×OBC finite-cluster diagonalization: `build_cluster_hamiltonian`, `solve_finite_cluster`, `corner_localization` (IPR per corner), `count_corner_states`; MAX cluster lx×ly×4≤64
- `MagnonBandModel3D` — 3D band model with kz extension: `cubic_haldane(j_nn,j_nnn,dmi,h_ext)` and `pyrochlore_topological(j,d)` presets (`src/topomagnon/band_model_3d.rs`)
- `AxionElectrodynamics` — Chern-Simons 3-form integral θ = (1/4π)∫_BZ tr(F−(2/3)A³); `axion_angle()` (0 or π), `topological_magnetoelectric_polarizability()`, `axion_response(E,B)` → emergent (P, M) (`src/topomagnon/axion.rs`)
- `AxionMagnonPhoton` — axion-mediated magnon-photon coupling: `magnon_photon_conversion_efficiency()`, `parity_violation_angle()` (Faraday/Kerr), `cooperativity()`

#### Data Export Formats (v0.6.0)

- `VtiWriter` — VTK ImageData (XML + base64-encoded f32 LE binary appended data); `write_vector_field`, `write_scalar_field`, `write_multi_field`; ParaView-compatible (`src/visualization/vti.rs`, feature `vti`)
- `XdmfWriter` — XDMF v3.0 temporal collection: `add_time_step`, `write(xmf_path, bin_path)`; raw f64 LE binary with seek offsets; ParaView-compatible (`src/visualization/xdmf.rs`, feature `xdmf`)
- `NetCdfWriter` / `NetCdfReader` — pure-Rust NetCDF3 Classic (CDF-1) binary: complete XDR big-endian encoding, dim_list, att_list, var_list, data section; `write_vector_field_cf` for CF-1.10 convention; round-trip reader (`src/visualization/netcdf.rs`, feature `netcdf`)
- `ZarrStore` / `ZarrArray` / `ZarrDtype` — pure-Rust Zarr v2 on-disk store: hand-written `.zarray` JSON metadata, C-order LE binary chunks, N-dimensional chunked I/O, round-trip reader (`src/visualization/zarr.rs`, feature `zarr`)

#### ML Autodiff (`src/autodiff/`, feature `autodiff`)

- `Tape` — reverse-mode AD tape with `push_leaf`, `push_op`, `backward(loss)`, `reset`; fixed-size arrays (no per-op allocation)
- `Var<'t>` — tape-bound variable: `leaf(tape, value)`, `value()`, `grad()`; arithmetic via `std::ops` traits (Add, Sub, Mul, Div, Neg); transcendentals: `sin`, `cos`, `exp`, `ln`, `sqrt`, `tanh`, `powi`, `powf`, `abs`, `recip`; also Var op f64 for constant arithmetic
- `Sgd` — SGD with momentum; `Adam` — Kingma & Ba 2014; `LBfgs` — limited-memory BFGS two-loop recursion (Nocedal 1980)
- `ParameterFitter` — closure-based fitting loop (`fit`); returns `FitResult { final_params, final_loss, n_iterations, converged, loss_history }`
- `kittel_frequency_diff`, `zeeman_energy_diff`, `exchange_energy_diff`, `dmi_energy_diff`, `anisotropy_energy_diff`, `llg_torque_norm_diff` — differentiable physics functions

#### Examples (6 new)

- `damon_eshbach_nonreciprocity.rs` — DE non-reciprocity sweep in YIG/Py/CoFeB thin films
- `backward_volume_magnons.rs` — BVMSW dispersion with negative v_g and thickness mode quantization
- `hoti_corner_states.rs` — BBH phase diagram, quadrupole moment, corner state spectrum, breathing kagome
- `axion_magnon_photon.rs` — 3D axion angle, topological ME polarizability, Faraday rotation
- `data_export_formats.rs` — 3D skyrmion texture export to VTK/VTI/NetCDF/Zarr with round-trip verification
- `autodiff_parameter_fitting.rs` — reverse-mode AD demonstration, YIG α-fitting with Adam, SGD/Adam/L-BFGS comparison

### Changed

- Version bumped `0.5.0` → `0.6.0`
- Test count: 996 → 1200 (+204 new tests)
- Examples count: 35 → 41 (+6 new)
- `src/spinwave/mod.rs` — added `damon_eshbach`, `bvmsw`, `surface` submodules
- `src/magnon/mod.rs` — added `spectral` submodule (behind `scirs2` feature)
- `src/topomagnon/mod.rs` — added `wilson`, `hoti`, `axion`, `band_model_3d` submodules
- `src/visualization/mod.rs` — added feature-gated `vti`, `xdmf`, `netcdf`, `zarr` submodules
- `src/lib.rs` — added `#[cfg(feature = "autodiff")] pub mod autodiff`
- `Cargo.toml` — added `base64 = "0.22"` (optional); added features `vti`, `xdmf`, `netcdf`, `zarr`, `autodiff`

## [0.5.0] - 2026-05-17

### Added

#### Non-Collinear Magnetism (`src/noncollinear/`)
- `SpinSpiral` — cycloidal, helical, conical, and fan-structure spiral types; magnetization profile `M(r)`; spin structure factor; Landau-Lifshitz exchange energy
- `SpiralType` / `SpiralChirality` — enums for classifying spiral order and handedness
- `LuttingerTisza` — exchange Fourier transform J(q); ground-state wavevector search; spiral vs. ferromagnetic classification; frustration ratio; ordering temperature estimate
- `ExchangeInteraction` — bond vector + exchange coupling building block for LT models
- Presets: `SpinSpiral::terbium_manganese_oxide()`, `LuttingerTisza::j1j2_chain()`, `LuttingerTisza::antiferromagnet()`

#### Multiferroic / Magnetoelectric Coupling (`src/multiferroic/`)
- `MagnetoelectricTensor` — linear ME tensor α_ij; presets for BiFeO₃, TbMnO₃, Cr₂O₃; `electric_polarization_from_field`, `magnetization_from_efield`, Dzyaloshinskii bound
- `MultiferroicType` — Type-I, Type-II, Type-III classification
- `KnbMechanism` — KNB formula P ∝ e_ij × (S_i × S_j) with `tbmno3()` preset; chain polarization sum
- `InverseMagnetoelectric` — E-field-induced magnetisation change; energy density
- Free functions: `dzyaloshinskii_moriya_polarization`, `exchange_striction_polarization`, `toroidal_moment`
- Free functions: `spin_current_from_spiral`, `magnon_drag_contribution`

#### Quantum Spin Hall / Kane-Mele (`src/topomagnon/qsh.rs`)
- `KaneMeleModel` — 4-band Bloch Hamiltonian on honeycomb lattice in (A↑,B↑,A↓,B↓) basis; intrinsic SOC λ_SO, Rashba λ_R, staggered potential λ_V
- `z2_invariant()` — Fukui-Hatsugai discrete Chern method on honeycomb BZ parallelogram; TRIM Pfaffian fallback for λ_R≠0
- `band_gap()`, `is_topological()`, `energy_bands(kx,ky)`, `edge_spectrum(n_cells, kx_min, kx_max, n_kx)`
- Presets: `graphene_with_soc`, `topological_phase`, `trivial_phase`

#### Nonlinear Magnon Physics (`src/magnon/nonlinear.rs`)
- `FourMagnonScattering` — Kalinikos-Slavin dispersion; four-magnon coupling T_kk; Suhl first-order instability threshold; parametric growth rate σ = √[(|T|h)² − Δω²]; `yig()` preset
- `ParametricAmplification` — degenerate signal-idler amplifier; threshold field, gain coefficient, signal gain (dB), pump depletion; `degenerate_from_yig()` preset
- `NonlinearFmrLinewidth` — nonlinear linewidth broadening, foldover field, bistability threshold power, power saturation factor; `from_yig()` preset
- Free functions: `magnon_magnon_interaction_energy`, `four_magnon_relaxation_rate`, `suhl_spin_wave_instability_power`

#### Examples (4 new)
- `spin_spiral_tbmno3.rs` (⭐⭐⭐) — TbMnO₃ ground-state search via Luttinger-Tisza; KNB electric polarization; cycloidal vs. helical comparison; spin structure factor
- `bife_o3_multiferroic.rs` (⭐⭐⭐) — BiFeO₃/TbMnO₃/Cr₂O₃ ME database; linear ME effect; DM polarization; exchange striction; toroidal moment; switching energy
- `kane_mele_qsh.rs` (⭐⭐⭐⭐) — Z2 phase diagram; Rashba-driven topological transition; band structure; helical edge states in strip geometry
- `nonlinear_magnon_suhl.rs` (⭐⭐⭐⭐) — YIG Suhl instability threshold; parametric growth rate; nonlinear FMR linewidth; bistability; parametric amplifier gain

### Changed
- Version bumped 0.4.0 → 0.5.0
- `src/prelude.rs`: added exports for all new v0.5.0 types and free functions
- `src/magnon/mod.rs`: extended nonlinear re-exports to include free functions
- Test count 917 → 996 (79 new tests)
- Example count 31 → 35
- lib.rs banner and architecture section updated

## [0.4.0] - 2026-05-17

### Added

#### Math Primitives (`src/math/`)
- `Complex` — lightweight complex number with constants `ZERO`, `ONE`, `I`; methods `new`, `from_real`, `from_polar`, `add`, `sub`, `mul`, `div`, `scale`, `conj`, `neg`, `exp`, `pow_n`, `mul_i`
- `CMatrix` — dense N×N complex matrix (N ≤ 64); Gauss-Jordan inverse with partial pivoting; Householder tridiagonalization + TQLI QL eigendecomposition for Hermitian matrices

#### Quantum Magnonics (`src/quantum/`)
- `HolsteinPrimakoff` — Holstein-Primakoff boson mapping for spin-S systems; linear and quadratic expansion orders; YIG/AFM presets
- `BogoliubovTransform` — analytical diagonalization via `tanh(2θ)=-B/A`; vacuum occupation |v_k|²; ground-state energy; AFM square-lattice preset
- `ZeroPointFluctuations` — zero-point amplitude `√(ℏ/2mω)`, ground-state energy ½Σℏω, Casimir free energy; bridge from `spinwave::QuantizedModes`

#### NEGF Non-Equilibrium Transport (`src/negf/`)
- `Hamiltonian1D` — 1D tight-binding chain with uniform or disordered onsite energies
- `LeadSelfEnergy` / `SanchoRubio` — wide-band and iterative surface Green's function leads
- `GreenFunction` — retarded/advanced/lesser Green's functions, spectral function, DOS, LDOS
- `TransportCalculator` — Landauer transmission T(E), I-V curve, differential conductance
- `KeldyshSolver` — Keldysh Σ<, Σ>, G<, G>, non-equilibrium density
- `ShotNoise` — zero-frequency noise, Fano factor, Johnson-Nyquist thermal noise
- `SpinAccumulation1D` — FTCS and implicit backward-Euler spin diffusion, steady-state, Thomas algorithm

#### Topological Magnon Bands (`src/topomagnon/`)
- `MagnonBandModel` — Haldane honeycomb (2-band), Kagome (3-band), square-DMI; `hamiltonian_at(k)`, `bands()`, `band_gap()`
- `BerryCurvature` — sum-over-states and link-variable curvature; 2D BZ grid integration
- `ChernNumber` — Fukui-Hatsugai-Suzuki discrete method; Wilson loop; total Chern sum validation
- `EdgeModes` — block-tridiagonal strip diagonalization; IPR-based edge localization; chiral velocity
- `MagnonHallConductivity` — Matsumoto-Murakami thermal Hall conductivity; dilogarithm; temperature sweep

#### Cavity Extensions (`src/cavity/`)
- `TavisCummings` — collective coupling g√N, polariton frequencies, mean-field Rabi dynamics, superradiant threshold, cooperativity; `yig_ensemble` preset
- `MagnonPolariton` — Hopfield diagonalization, photon/magnon fractions, vacuum Rabi splitting, strong-coupling detection; `from_yig_cavity` preset
- `MultiModePolariton` — multi-mode Jaynes-Cummings, real-symmetric Jacobi eigendecomposition
- `BrillouinScattering` — scattering rate, Stokes/anti-Stokes frequency shift, cross section
- `OptomagnonicCoupling` — optomagnonic coupling, Kerr shift
- `MicrowaveToOptical` — conversion efficiency 4C_me·C_mo/(1+C_me+C_mo)², bandwidth, impedance matching
- `MagnonicFrequencyComb` — comb spectrum, phase noise, coherence time

#### Random Anisotropy Disorder (`src/material/random_anisotropy.rs`)
- `RandomAnisotropy` — Gaussian/Uniform/FixedAxes distributions; Marsaglia sphere sampling; Imry-Ma correlation length; Harris criterion; LLG coupling via `effective_field()`; `nanocrystalline` preset

#### New Examples (6)
- `magnon_zero_point` — YIG quantized modes → zero-point amplitudes, Casimir free energy, Bogoliubov squeeze
- `negf_transport` — tight-binding chain, Landauer I-V, shot noise/Fano factor, Anderson disorder
- `tavis_cummings_dicke` — √N collective coupling, anticrossing, superradiant threshold, mean-field dynamics
- `magnon_polariton` — anticrossing sweep, Hopfield fractions, coupling regimes, multi-mode polariton
- `topological_magnon_haldane` — Berry curvature, Chern number, DMI phase transition, chiral edge modes
- `random_anisotropy_disorder` — Imry-Ma length, Harris criterion, LLG with disorder dephasing

### Changed
- `cargo check` no longer fails due to `oxiarc-lz4`/`oxiarc-zstd`: `scirs2-core` and `scirs2-spatial` now use `default-features = false`
- `math::Complex` is now the canonical complex number type; `magnon::bec::Complex` is a deprecated re-export
- Test count increased from 718 (v0.3.0) to 917 (v0.4.0)
- 25 examples → 31 examples

### Deprecated
- `spintronics::magnon::bec::Complex` — use `spintronics::math::Complex` instead (deprecated since 0.4.0)

### Fixed
- `CMatrix::hermitian_eigendecomposition` — replaced naive Jacobi with Householder tridiagonalization + implicit QL (Wilkinson shift) to fix convergence for matrices with large diagonal/off-diagonal ratios (e.g. 10 GHz / 100 MHz = 100:1)
- Zero clippy warnings maintained across all modules and examples

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
