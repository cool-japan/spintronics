# TODO List for Spintronics Library

**Version**: 0.3.2
**Last Updated**: 2026-07-06 - v0.3.2 in development
**Status**: 1977 lib + 65 proptest passing = 2042 nextest (main crate, `cargo nextest run --all-features --workspace`) + 118 doctests passing, 5 ignored (`cargo test --all-features --doc`), 0 warnings, ~101K lines (Rust code: ~79.5K+)

---

## Stubs to implement (added 2026-06-12 by /cooljapan-stub-check)

- [x] `spintronics`: `src/magnon/nonlinear.rs` — replaced placeholder pump_h (0.1 mT) with the
  physics-derived degenerate parametric **threshold field** `h_th = √(γ_s·γ_i)/κ = α·ω_p/(γ_gyro·Ms)`
  in `ParametricAmplification::degenerate_from_yig`. Added the `with_supercriticality(ξ)` builder so
  the preset (which now sits at marginal stability) can be moved to any operating point ξ = h_p/h_th.
  3 new tests (threshold identity, α/ω_p/Ms scaling, supercriticality operating point); example
  `nonlinear_magnon_suhl` updated to twice-critical operation. (2026-06-22, P2, done)

### Build hygiene fixes (2026-06-22)
- [x] Added `required-features = ["autodiff"]` to 7 previously-ungated autodiff examples
  (`autodiff_parameter_fitting`, `neural_exchange_training`, `pinn_llg_solver`, `equivariant_nn_demo`,
  `active_learning_demo`, `graph_nn_lattice`, `bayesian_opt_materials`) so `cargo build --examples`
  no longer fails under default features.
- [x] Fixed dead-code warning for `temp_path` in `data_export_formats` example via a precise
  `cfg_attr(not(any(feature = "vti"/"netcdf"/"zarr")), allow(dead_code))`.
- Verified: `cargo clippy --all-targets -- -D warnings` clean for both default and `autodiff` features.

### Binary OVF I/O implemented (2026-06-22 by /stub-check)
- [x] `spintronics`: `src/io/ovf.rs` — implemented binary OVF read/write, replacing the
  `Err("Binary OVF format not yet implemented")` stub. Supports `Binary4_1_0` (OVF 1.0, **big-endian**),
  `Binary4_2_0` and `Binary8_2_0` (OVF 2.0, **little-endian**) per the OOMMF OVF spec; control/check
  values `1234567.0_f32` / `123456789012345.0_f64` validated on read to catch byte-order mismatch;
  byte-oriented reader with truncation guards (text path unchanged). Header writing refactored into
  shared `write_header_1_0/2_0` helpers. +5 round-trip/error tests (io::ovf now 9 tests). File 1162 lines.
- Remaining intentional/blocked stubs (NOT implemented, correctly so): CUDA backend (`src/gpu/cuda.rs`,
  needs `cudarc` + GPU hardware → v1.0.0 roadmap); HDF5 feature-off fallback (`src/visualization/hdf5.rs`,
  correct feature-gating, full impl active under `--features hdf5`).

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
  - [x] Implicit midpoint with Newton iteration (deferred) — shipped: src/dynamics/integrators/implicit_midpoint.rs (test: test_newton_converges_quickly)
  - [x] Crank-Nicolson for diffusion-dominated systems (deferred) — shipped: src/dynamics/integrators/crank_nicolson.rs (test: test_unconditional_stability_large_dt)
- [ ] Spectral methods for periodic systems (deferred to v0.4.0)

### Priority 2: Spin Wave Theory
- [x] Analytical dispersion relations for thin films — `spinwave` module
- [x] Magnon dispersion and band structure — `spin_wave_dispersion` example
- [x] Damon-Eshbach modes for in-plane magnetized films (deferred) — shipped: src/spinwave/damon_eshbach.rs (test: test_nonreciprocity_positive)
- [x] Backward volume modes for perpendicular magnetization (deferred) — shipped: src/spinwave/bvmsw.rs (test: test_backward_group_velocity_small_k)
- [x] Surface spin waves in semi-infinite media (deferred) — shipped: src/spinwave/surface.rs, semi_infinite_de.rs (test: test_penetration_depth_decreases_with_k)
- [x] Spin wave quantization in nanostructures (deferred) — shipped: src/spinwave/quantization.rs, nanodisk.rs (analytic stripe/disk/rectangle + nanodisk radial-azimuthal quantization)
- [x] Mode decomposition and spectral analysis (deferred) — shipped: src/magnon/spectral.rs (test: test_mode_decomposition_uniform_field)

### Priority 3: Altermagnets
- [x] RuO2 material model and parameters — `altermagnet` module, `altermagnet_ruo2` example
- [x] Spin splitter effect — implemented in `altermagnet`
- [x] Anomalous Hall effect in altermagnets — implemented
- [x] CrSb material model (deferred) (done 2026-07-05)
  - **Goal:** Momentum-resolved electronic altermagnet band model with per-spin block-diagonal Hamiltonians, closed-form Berry curvature, Fermi-sea crystal/spin Hall (correctly zero without spin-orbit coupling, nonzero and Neel-sign-flipping with optional SOC), and an altermagnet spin-valve GMR-without-ferromagnetism device. Full design in Batch C1 of the woolly-prancing-squid planning notes (local Claude Code plan file, outside this repo).
  - **Design:** Collinear Neel order along z, no SOC by default so spin is block-diagonal; per-spin H_sigma(k) = eps0(k)*tau0 + d_sigma(k).tau with d_sigma = (Re gamma, -Im gamma, eps_a(k) - sigma*M), altermagnetic form factor eps_a(k) = t_AM*(ak)^2*cos(n*(phi-phi_Neel)) with n = AltermagneticSymmetry::harmonic_order() (d/g/i = 2/4/6). SOC (when enabled) enters via a mixed-k-parity structure factor: b_x/b_z odd in k, b_y even in k (required by the Theta=i*tau_y*K operator identity Theta H(k;M) Theta^-1 = H(-k;-M); a fully-odd b_y would instead make the Neel-reversal Hall response even, incorrectly vanishing). New BlochHamiltonian trait (src/altermagnet/bloch_hamiltonian.rs: `hamiltonian_at`/`n_bands`/`diagonalize_at`) implemented for a per-spin sector via `AltermagnetSpinHamiltonian`, both paths (closed-form eigenvalues and the explicit 2x2 matrix) built from the same `spin_d_vector`/`eps0` so they cannot diverge; plus a standalone KuboBerry finite-difference cross-check calculator (src/altermagnet/kubo_berry.rs) built purely against the new trait. src/topomagnon/* shared code (BerryCurvature, MagnonBandModel, Chern/Wilson/axion/edge-modes) is NOT modified, to avoid Phase-3 test blast radius. Correct symmetry invariant is spin-group E_up(k) = E_down(R k), where R is a proper rotation by pi/n (n = harmonic_order()) -- generalizes to every harmonic order, NOT just the d-wave-specific (kx,ky)->(ky,kx) swap and NOT E_up(k) = E_down(-k).
  - **Files:** src/altermagnet/band_model.rs (new), src/altermagnet/bloch_hamiltonian.rs (new: BlochHamiltonian trait + AltermagnetSpinHamiltonian), src/altermagnet/kubo_berry.rs (new: KuboBerry cross-check calculator), src/altermagnet/spin_valve.rs (new), src/altermagnet/mod.rs (+exports), src/prelude.rs (+exports for AltermagnetBandModel/Band/Spin/SpinBands only -- BlochHamiltonian/AltermagnetSpinHamiltonian/KuboBerry are reached via `spintronics::altermagnet::*`, not re-exported from the prelude), Cargo.toml (+examples), examples/altermagnet_band_structure.rs, examples/altermagnet_spin_valve.rs, tests/property_altermagnet.rs (new)
  - **Prerequisites:** none - scalar Altermagnet/AltermagneticSymmetry presets (materials.rs), MagnonBandModel/BerryCurvature patterns (reference only, not modified), and MagneticMultilayer GMR template (multilayer.rs) all already shipped
  - **Tests:** Hermiticity of H_sigma(k) (unit + property test, via the BlochHamiltonian trait); spin-degeneracy when M=0 or when t_AM=0 (latter = ordinary AFM limit); nodal lines at node_angles(); each spin even in k; spin-group relation E_up(k)=E_down(Rk) via the pi/n rotation, verified for d-wave (unit) and generalized to g-wave/i-wave (unit + property, any harmonic order); Berry curvature band-sum-to-zero (Omega_Lower(k)+Omega_Upper(k)=0, any spin) and oddness under the combined spin+momentum+Neel-order inversion Omega_sigma(k;M) = -Omega_{-sigma}(-k;-M) -- the naive same-spin Omega_sigma(k)=-Omega_sigma(-k) does NOT hold in general once delta_hyb and SOC are both nonzero (verified false by direct computation, residual comparable to the curvature itself), so the combined-flip form -- a direct consequence of the Theta=i*tau_y*K operator identity -- is what is tested; charge Hall approx 0 at lambda_soc=0, nonzero and sign-flipping under Neel reversal at lambda_soc!=0; zero net spin polarization (exact for d-wave at any parameters/Fermi energy/Neel angle; MnTe and CrSb presets individually asserted at their default Fermi energy); closed-form vs KuboBerry finite-difference cross-check agreement (relative agreement with SOC on, both vanish together with SOC off); GMR-without-FM: zero net magnetization in every configuration plus angle-dependent resistance
  - **Risk:** highest of the three hard-physics modules (new trait, new physics domain, Hall-transport correctness subtlety) - mitigated by the BlochHamiltonian/KuboBerry cross-check and the invariant test battery above (property_altermagnet.rs plus band_model.rs/bloch_hamiltonian.rs/kubo_berry.rs unit tests); must NOT modify (and does not modify) src/topomagnon/* shared code
- [x] MnTe material model (deferred) (done 2026-07-05)
  - **Goal:** Momentum-resolved electronic altermagnet band model with per-spin block-diagonal Hamiltonians, closed-form Berry curvature, Fermi-sea crystal/spin Hall (correctly zero without spin-orbit coupling, nonzero and Neel-sign-flipping with optional SOC), and an altermagnet spin-valve GMR-without-ferromagnetism device. Full design in Batch C1 of the woolly-prancing-squid planning notes (local Claude Code plan file, outside this repo).
  - **Design:** Collinear Neel order along z, no SOC by default so spin is block-diagonal; per-spin H_sigma(k) = eps0(k)*tau0 + d_sigma(k).tau with d_sigma = (Re gamma, -Im gamma, eps_a(k) - sigma*M), altermagnetic form factor eps_a(k) = t_AM*(ak)^2*cos(n*(phi-phi_Neel)) with n = AltermagneticSymmetry::harmonic_order() (d/g/i = 2/4/6). SOC (when enabled) enters via a mixed-k-parity structure factor: b_x/b_z odd in k, b_y even in k (required by the Theta=i*tau_y*K operator identity Theta H(k;M) Theta^-1 = H(-k;-M); a fully-odd b_y would instead make the Neel-reversal Hall response even, incorrectly vanishing). New BlochHamiltonian trait (src/altermagnet/bloch_hamiltonian.rs: `hamiltonian_at`/`n_bands`/`diagonalize_at`) implemented for a per-spin sector via `AltermagnetSpinHamiltonian`, both paths (closed-form eigenvalues and the explicit 2x2 matrix) built from the same `spin_d_vector`/`eps0` so they cannot diverge; plus a standalone KuboBerry finite-difference cross-check calculator (src/altermagnet/kubo_berry.rs) built purely against the new trait. src/topomagnon/* shared code (BerryCurvature, MagnonBandModel, Chern/Wilson/axion/edge-modes) is NOT modified, to avoid Phase-3 test blast radius. Correct symmetry invariant is spin-group E_up(k) = E_down(R k), where R is a proper rotation by pi/n (n = harmonic_order()) -- generalizes to every harmonic order, NOT just the d-wave-specific (kx,ky)->(ky,kx) swap and NOT E_up(k) = E_down(-k).
  - **Files:** src/altermagnet/band_model.rs (new), src/altermagnet/bloch_hamiltonian.rs (new: BlochHamiltonian trait + AltermagnetSpinHamiltonian), src/altermagnet/kubo_berry.rs (new: KuboBerry cross-check calculator), src/altermagnet/spin_valve.rs (new), src/altermagnet/mod.rs (+exports), src/prelude.rs (+exports for AltermagnetBandModel/Band/Spin/SpinBands only -- BlochHamiltonian/AltermagnetSpinHamiltonian/KuboBerry are reached via `spintronics::altermagnet::*`, not re-exported from the prelude), Cargo.toml (+examples), examples/altermagnet_band_structure.rs, examples/altermagnet_spin_valve.rs, tests/property_altermagnet.rs (new)
  - **Prerequisites:** none - scalar Altermagnet/AltermagneticSymmetry presets (materials.rs), MagnonBandModel/BerryCurvature patterns (reference only, not modified), and MagneticMultilayer GMR template (multilayer.rs) all already shipped
  - **Tests:** Hermiticity of H_sigma(k) (unit + property test, via the BlochHamiltonian trait); spin-degeneracy when M=0 or when t_AM=0 (latter = ordinary AFM limit); nodal lines at node_angles(); each spin even in k; spin-group relation E_up(k)=E_down(Rk) via the pi/n rotation, verified for d-wave (unit) and generalized to g-wave/i-wave (unit + property, any harmonic order); Berry curvature band-sum-to-zero (Omega_Lower(k)+Omega_Upper(k)=0, any spin) and oddness under the combined spin+momentum+Neel-order inversion Omega_sigma(k;M) = -Omega_{-sigma}(-k;-M) -- the naive same-spin Omega_sigma(k)=-Omega_sigma(-k) does NOT hold in general once delta_hyb and SOC are both nonzero (verified false by direct computation, residual comparable to the curvature itself), so the combined-flip form -- a direct consequence of the Theta=i*tau_y*K operator identity -- is what is tested; charge Hall approx 0 at lambda_soc=0, nonzero and sign-flipping under Neel reversal at lambda_soc!=0; zero net spin polarization (exact for d-wave at any parameters/Fermi energy/Neel angle; MnTe and CrSb presets individually asserted at their default Fermi energy); closed-form vs KuboBerry finite-difference cross-check agreement (relative agreement with SOC on, both vanish together with SOC off); GMR-without-FM: zero net magnetization in every configuration plus angle-dependent resistance
  - **Risk:** highest of the three hard-physics modules (new trait, new physics domain, Hall-transport correctness subtlety) - mitigated by the BlochHamiltonian/KuboBerry cross-check and the invariant test battery above (property_altermagnet.rs plus band_model.rs/bloch_hamiltonian.rs/kubo_berry.rs unit tests); must NOT modify (and does not modify) src/topomagnon/* shared code
- [x] Giant magnetoresistance without ferromagnetism (deferred) (done 2026-07-05)
  - **Goal:** Momentum-resolved electronic altermagnet band model with per-spin block-diagonal Hamiltonians, closed-form Berry curvature, Fermi-sea crystal/spin Hall (correctly zero without spin-orbit coupling, nonzero and Neel-sign-flipping with optional SOC), and an altermagnet spin-valve GMR-without-ferromagnetism device. Full design in Batch C1 of the woolly-prancing-squid planning notes (local Claude Code plan file, outside this repo).
  - **Design:** Collinear Neel order along z, no SOC by default so spin is block-diagonal; per-spin H_sigma(k) = eps0(k)*tau0 + d_sigma(k).tau with d_sigma = (Re gamma, -Im gamma, eps_a(k) - sigma*M), altermagnetic form factor eps_a(k) = t_AM*(ak)^2*cos(n*(phi-phi_Neel)) with n = AltermagneticSymmetry::harmonic_order() (d/g/i = 2/4/6). SOC (when enabled) enters via a mixed-k-parity structure factor: b_x/b_z odd in k, b_y even in k (required by the Theta=i*tau_y*K operator identity Theta H(k;M) Theta^-1 = H(-k;-M); a fully-odd b_y would instead make the Neel-reversal Hall response even, incorrectly vanishing). New BlochHamiltonian trait (src/altermagnet/bloch_hamiltonian.rs: `hamiltonian_at`/`n_bands`/`diagonalize_at`) implemented for a per-spin sector via `AltermagnetSpinHamiltonian`, both paths (closed-form eigenvalues and the explicit 2x2 matrix) built from the same `spin_d_vector`/`eps0` so they cannot diverge; plus a standalone KuboBerry finite-difference cross-check calculator (src/altermagnet/kubo_berry.rs) built purely against the new trait. src/topomagnon/* shared code (BerryCurvature, MagnonBandModel, Chern/Wilson/axion/edge-modes) is NOT modified, to avoid Phase-3 test blast radius. Correct symmetry invariant is spin-group E_up(k) = E_down(R k), where R is a proper rotation by pi/n (n = harmonic_order()) -- generalizes to every harmonic order, NOT just the d-wave-specific (kx,ky)->(ky,kx) swap and NOT E_up(k) = E_down(-k).
  - **Files:** src/altermagnet/band_model.rs (new), src/altermagnet/bloch_hamiltonian.rs (new: BlochHamiltonian trait + AltermagnetSpinHamiltonian), src/altermagnet/kubo_berry.rs (new: KuboBerry cross-check calculator), src/altermagnet/spin_valve.rs (new), src/altermagnet/mod.rs (+exports), src/prelude.rs (+exports for AltermagnetBandModel/Band/Spin/SpinBands only -- BlochHamiltonian/AltermagnetSpinHamiltonian/KuboBerry are reached via `spintronics::altermagnet::*`, not re-exported from the prelude), Cargo.toml (+examples), examples/altermagnet_band_structure.rs, examples/altermagnet_spin_valve.rs, tests/property_altermagnet.rs (new)
  - **Prerequisites:** none - scalar Altermagnet/AltermagneticSymmetry presets (materials.rs), MagnonBandModel/BerryCurvature patterns (reference only, not modified), and MagneticMultilayer GMR template (multilayer.rs) all already shipped
  - **Tests:** Hermiticity of H_sigma(k) (unit + property test, via the BlochHamiltonian trait); spin-degeneracy when M=0 or when t_AM=0 (latter = ordinary AFM limit); nodal lines at node_angles(); each spin even in k; spin-group relation E_up(k)=E_down(Rk) via the pi/n rotation, verified for d-wave (unit) and generalized to g-wave/i-wave (unit + property, any harmonic order); Berry curvature band-sum-to-zero (Omega_Lower(k)+Omega_Upper(k)=0, any spin) and oddness under the combined spin+momentum+Neel-order inversion Omega_sigma(k;M) = -Omega_{-sigma}(-k;-M) -- the naive same-spin Omega_sigma(k)=-Omega_sigma(-k) does NOT hold in general once delta_hyb and SOC are both nonzero (verified false by direct computation, residual comparable to the curvature itself), so the combined-flip form -- a direct consequence of the Theta=i*tau_y*K operator identity -- is what is tested; charge Hall approx 0 at lambda_soc=0, nonzero and sign-flipping under Neel reversal at lambda_soc!=0; zero net spin polarization (exact for d-wave at any parameters/Fermi energy/Neel angle; MnTe and CrSb presets individually asserted at their default Fermi energy); closed-form vs KuboBerry finite-difference cross-check agreement (relative agreement with SOC on, both vanish together with SOC off); GMR-without-FM: zero net magnetization in every configuration plus angle-dependent resistance
  - **Risk:** highest of the three hard-physics modules (new trait, new physics domain, Hall-transport correctness subtlety) - mitigated by the BlochHamiltonian/KuboBerry cross-check and the invariant test battery above (property_altermagnet.rs plus band_model.rs/bloch_hamiltonian.rs/kubo_berry.rs unit tests); must NOT modify (and does not modify) src/topomagnon/* shared code

### Priority 4: Orbitronics
- [x] Orbital Hall effect — `orbitronics` module
- [x] Orbital torques — implemented
- [x] Orbital-to-spin conversion — implemented
- [x] d-orbital magnetism (deferred) (done 2026-07-05)
  - **Goal:** Local d-orbital moment model: crystal-field splitting (octahedral/tetrahedral/tetragonal), Hund's-rule free-ion terms, high-spin/low-spin, orbital quenching plus SOC-driven unquenching, effective moments, ~12 preset ions. Full design in Batch C2 of the woolly-prancing-squid planning notes (local Claude Code plan file, outside this repo).
  - **Design:** Build L_x,L_y,L_z operators in the complex |l=2,m> basis (exact ladder operators) and rotate to the real cubic-harmonic basis via a fixed unitary U - do not hand-type the real matrices in code (use them only as test fixtures). Diagonal crystal-field Hamiltonians for O_h (t2g -0.4*10Dq, eg +0.6*10Dq), T_d (inverted), tetragonal (Ballhausen Ds/Dt). H_SOC = lambda*L.S on the 10x10 orbital-tensor-spin space. Documented fidelity boundary: single-configuration ligand-field + atomic Hund's rules + single-particle SOC, NOT full many-electron multiplet/Tanabe-Sugano CI.
  - **Files:** src/orbitronics/crystal_field.rs (new - operators), src/orbitronics/d_orbital_moment.rs (new - CrystalFieldModel, Hund's rules, moments, presets), src/orbitronics/mod.rs (+exports), examples/d_orbital_local_moments.rs, Cargo.toml (+example)
  - **Prerequisites:** none - OrbitalHallMaterial exists as a loose one-accessor bridge only, no type coupling
  - **Tests:** full Hund's-rule d1-d9 (S,L,J) table including spot checks (d5 to S=5/2,L=0; d2 to S=1,L=3; d7 to S=3/2,L=3); electron-hole symmetry; Lande g-factor; spin-only mu_eff textbook values (Fe3+ 5.92, Ni2+ 2.83, Cr3+ 3.87); operator Hermiticity, [L_x,L_y]=iL_z commutator, L^2=6*I, L_z spectrum {-2..2}; orbital quenching for A/E ground terms (zero diagonal L_i); SOC-driven unquenching for T ground terms; 10Dq to infinity forces low-spin; preset ion table matches expected S/L/mu values
  - **Risk:** moderate - self-contained new module, no shared-code touch, but real risk is physics-correctness (operator algebra bugs) rather than integration blast radius

### Priority 5: Frustrated Magnets
- [x] Spin ice models — `frustrated` module, `spin_ice_monopoles` example
- [x] Kagome lattice magnets — implemented in `frustrated`
- [x] Spin liquids (resonating valence bond states) (deferred) (done 2026-07-05)
  - **Goal:** Resonating-valence-bond quantum module: dimer-covering enumeration, Sutherland loop-counting overlaps, VB-basis Heisenberg matrix elements, generalized-eigenproblem ground state, exact-diagonalization benchmark (dense + Lanczos), spinon/deconfinement diagnostics. Full design in Batch C3 of the woolly-prancing-squid planning notes (local Claude Code plan file, outside this repo).
  - **Design:** Abstract from_bonds(num_sites, bonds, J) API (open validation clusters, sign-free) plus from_lattice(&FrustratedLattice) (classical lattices use periodic BCs that would distort tiny-cluster matchings). Bipartite (positive-sign) path first, then non-bipartite signs (no Monte-Carlo sign problem since this is exact small-D linear algebra). CMatrix::MAX_DIM=64 drives the design: VB basis stays on CMatrix (cap D<=64, explicit error carrying the measured covering count if exceeded - never a silent cap); a matrix-free Lanczos S_z=0 exact-diagonalization backend (N<=16) is required for the benchmark since dense CMatrix ED already exceeds MAX_DIM at N=8 (C(8,4)=70>64). All matrices are real symmetric.
  - **Files:** src/frustrated/rvb/mod.rs, dimer.rs, valence_bond.rs, solver.rs, exact.rs, spinon.rs, diagnostics.rs (all new), src/frustrated/mod.rs (+exports), examples/rvb_triangular_spin_liquid.rs, Cargo.toml (+example)
  - **Prerequisites:** none - FrustratedLattice (neighbors/spins/positions) and CMatrix::hermitian_eigendecomposition already shipped as the foundation
  - **Tests:** 2-site exact singlet E=-3J/4; 4-ring exact E=-2J matching a worked S/H example; covering counts vs known combinatorics (4-ring=2, 2x4 ladder=5, 4x4 square=36); loop-rule overlaps vs direct spin-basis inner products including a non-bipartite cluster; total spin S_tot^2 approx 0; variational ordering E_equal_amplitude >= E_ground_state >= E_exact_diagonalization; dense-vs-Lanczos ED agreement; explicit error (not silent truncation) when covering count exceeds MAX_VB_BASIS=64
  - **Risk:** moderate-high - most numerically subtle (overcomplete non-orthogonal basis, generalized eigenproblem, non-bipartite sign bookkeeping) of the three hard modules; self-contained new subtree, no shared-code touch
  - **Shipped:** src/frustrated/rvb/{mod,dimer,valence_bond,solver,exact,spinon,diagnostics}.rs (53 tests), examples/rvb_triangular_spin_liquid.rs. Bonus fix: while cross-checking exact.rs's total-spin-squared invariant, found and fixed a pre-existing sign/phase bug in `CMatrix::hermitian_eigendecomposition` (src/math/matrix.rs) — eigenvalues were always correct, but eigenvectors were wrong for any non-diagonal matrix of size >=3 (never caught since existing tests only checked eigenvalue correctness or eigenvector orthonormality, not the eigenvalue equation itself, beyond n<=2). Replaced the Householder+QL implementation with a cyclic Jacobi algorithm (phase pre-rotation + real 2x2 rotation per sweep step), verified to machine precision on real and complex Hermitian cases up to n=16; full workspace suite (1989 tests) passes after the fix, including several other consumer modules (topomagnon band models, orbitronics, hopfion stability modes) that were silently relying on eigenvectors this whole time.
- [x] Geometric frustration effects on transport (deferred) (done 2026-07-05) -- shipped: src/frustrated/transport.rs (FrustratedTransport, plaquette chirality/solid-angle/emergent-field), src/frustrated/lattice.rs (set_umbrella_order), src/effect/topological_hall.rs (+emergent_field_from_solid_angle, +hall_resistivity_from_charge_density) (test: test_kagome_umbrella_gives_nonzero_uniform_sign_chirality_and_hall_response)
  - **Goal:** Wire the scalar spin-chirality of a FrustratedLattice configuration into the existing topological-Hall transport machinery, giving a working chirality-to-emergent-field-to-Hall-response pipeline for frustrated (triangular/kagome/pyrochlore) spin textures.
  - **Design:** New function/struct (e.g. FrustratedTransport or a free fn frustration_hall_response) computing per-plaquette scalar chirality chi = S_i . (S_j x S_k) from FrustratedLattice's spins/neighbors/positions, converting to an emergent field, and feeding src/effect/topological_hall.rs's TopologicalHallEffect.
  - **Files:** src/frustrated/lattice.rs or new src/frustrated/transport.rs, src/effect/topological_hall.rs (extend, additive only), src/frustrated/mod.rs, examples/frustrated_transport.rs
  - **Prerequisites:** none - both sides (FrustratedLattice, TopologicalHallEffect) already shipped
  - **Tests:** zero Hall response for a coplanar spin configuration; nonzero for a noncoplanar 120-degree/umbrella configuration; sign reversal under chirality reversal; consistency with existing topological_hall.rs conventions
  - **Risk:** low-moderate - touches effect/topological_hall.rs, but additively (new entry point, no signature changes to existing public API)

### Priority 6: Hopfions
- [x] 3D topological soliton structure — `texture/hopfion_dynamics.rs`
- [x] Hopf index computation — Berry-connection (Whitehead) method
- [x] Current-driven dynamics — per-site LLG RK4 on 3D grid
- [x] Stability analysis (deferred) (done 2026-07-05)
  - **Goal:** Dynamical/eigenmode linear-stability analysis for a relaxed hopfion (the normal-mode spectrum / growth rates around equilibrium), complementing the already-shipped energy-landscape (collapse/expansion-radius) stability in src/texture/hopfion.rs.
  - **Design:** A linear-stability analyzer computing the Hessian of the total-energy functional (HopfionEnergy in hopfion.rs) about a relaxed configuration via finite differences, then diagonalizing it (crate::math::CMatrix::hermitian_eigendecomposition) to get the normal-mode eigenvalue spectrum.
  - **Files:** src/texture/hopfion.rs (extend) or new src/texture/hopfion_stability_modes.rs if hopfion.rs would cross 2000 lines (check with rslines/wc -l first; split via splitrs if needed)
  - **Prerequisites:** none - HopfionEnergy and HopfionStability (energy-landscape variant) already shipped as the foundation
  - **Tests:** positive-definite eigenvalue spectrum at a genuinely stable radius; at least one negative eigenvalue past the known collapse-radius boundary (from existing find_stability_boundaries); zero-modes corresponding to translation invariance (within numerical tolerance)
  - **Risk:** low-moderate - check hopfion.rs current line count before extending; split via splitrs if it would cross 2000 lines

### Priority 7: Magnon BEC & Spin Density Waves
- [x] Magnon Bose-Einstein condensation — `magnon/bec`, `magnon_bec` example
- [x] Spin density wave formation and dynamics (deferred) (done 2026-07-05) — shipped: `src/magnon/spin_density_wave.rs` (`SdwRelaxationDynamics`, test: `test_sdw_relaxation_converges_to_self_consistent_gap`)
  - **Goal:** Time-domain evolution for the spin-density-wave order parameter (amplitude + phase), relaxing toward the existing self-consistent BCS-like gap, complementing the already-shipped static/thermodynamic SDW physics.
  - **Design:** New time-stepping method(s) on or alongside SpinDensityWave/SdwGapSolver (src/magnon/spin_density_wave.rs) implementing TDGL-style relaxational dynamics: d(amplitude)/dt proportional to -dF/d(amplitude) toward the self-consistent gap at temperature T, plus phase dynamics if driven.
  - **Files:** src/magnon/spin_density_wave.rs (extend)
  - **Prerequisites:** none - static SDW machinery (order parameter, gap solver, energies) already shipped and is the foundation
  - **Tests:** relaxation converges to the known self-consistent equilibrium gap value; free energy monotonically decreases during relaxation; amplitude vanishes continuously as T approaches T_N
  - **Risk:** low - additive extension of an existing, well-tested module
- [x] Helical magnets and spirals (deferred) — shipped: src/noncollinear/spiral.rs (test: test_helical_gives_no_polarization)
- [ ] Magnon-magnon interactions (deferred)

### Priority 8: Magnetoelastic Coupling (Straintronics)
- [x] Magnetoelastic coupling tensor — `mech/magnetoelastic`
- [x] Strain-induced anisotropy — implemented
- [x] Piezoelectric control of magnetism (deferred) (done 2026-07-05) — shipped: src/mech/strain_driven_dynamics.rs (StrainDrivenLlgDriver, PiezoAcStrainDrive, test: test_piezo_driven_dynamics_evolves_and_conserves_norm)
  - **Goal:** A time-domain LLG driver coupling the existing periodic strain fields (SAW: SawSource::strain_at_point(x,z,t); piezo: PiezoelectricSubstrate/StraintronicDevice strain) into a magnetoelastic effective field that drives real magnetization dynamics, complementing the shipped steady-state/analytical SAW and piezo-strain models.
  - **Design:** New driver struct/function coupling strain(t) to magnetoelastic anisotropy field (reusing MagnetoelasticMaterial's villari_effective_field_vector / stress_induced_anisotropy machinery in src/mech/magnetoelastic.rs) to time-stepped LLG integration (reusing existing dynamics::llg / integrators infrastructure) over a driven simulation window.
  - **Files:** src/mech/saw.rs (extend) and/or new src/mech/strain_driven_dynamics.rs, integrating with src/dynamics/llg.rs and existing integrators
  - **Prerequisites:** none - all referenced components (SawSource, MagnetoelasticMaterial, LLG integrators) already shipped
  - **Tests:** resonant magnetization response at the acoustic-FMR condition (acoustic_fmr_condition_h_ext); suppressed response off-resonance; net energy pumping into the magnetic subsystem over a drive cycle consistent with resonant_absorption()
  - **Risk:** moderate - integrates across two subsystems (mech/ and dynamics/); keep the coupling as a new entry point, not a modification of either subsystem's existing public API
- [x] Surface acoustic wave (SAW) driven spin dynamics (deferred) (done 2026-07-05) — shipped: src/mech/strain_driven_dynamics.rs (StrainDrivenLlgDriver::from_saw_magnetoacoustics, SawStrainDrive, test: test_resonant_response_exceeds_off_resonant)
  - **Goal:** A time-domain LLG driver coupling the existing periodic strain fields (SAW: SawSource::strain_at_point(x,z,t); piezo: PiezoelectricSubstrate/StraintronicDevice strain) into a magnetoelastic effective field that drives real magnetization dynamics, complementing the shipped steady-state/analytical SAW and piezo-strain models.
  - **Design:** New driver struct/function coupling strain(t) to magnetoelastic anisotropy field (reusing MagnetoelasticMaterial's villari_effective_field_vector / stress_induced_anisotropy machinery in src/mech/magnetoelastic.rs) to time-stepped LLG integration (reusing existing dynamics::llg / integrators infrastructure) over a driven simulation window.
  - **Files:** src/mech/saw.rs (extend) and/or new src/mech/strain_driven_dynamics.rs, integrating with src/dynamics/llg.rs and existing integrators
  - **Prerequisites:** none - all referenced components (SawSource, MagnetoelasticMaterial, LLG integrators) already shipped
  - **Tests:** resonant magnetization response at the acoustic-FMR condition (acoustic_fmr_condition_h_ext); suppressed response off-resonance; net energy pumping into the magnetic subsystem over a drive cycle consistent with resonant_absorption()
  - **Risk:** moderate - integrates across two subsystems (mech/ and dynamics/); keep the coupling as a new entry point, not a modification of either subsystem's existing public API

### Priority 9: Advanced Disorder & Defects
- [x] Random anisotropy model (deferred to v0.4.0) — shipped: src/material/random_anisotropy.rs (test: test_imry_ma_length_positive)
- [x] Grain boundary effects in polycrystalline films (deferred) — shipped: src/material/defects.rs (GrainBoundary), src/material/disorder.rs (GrainStructure)
- [x] Point defects and pinning sites (deferred) — shipped: src/material/defects.rs (DefectSite/DepinningParams, test: coercivity_estimate)
- [x] Surface roughness modeling (deferred) — shipped: src/material/disorder.rs (SurfaceRoughness::generate, test: actual_rms)
- [x] Inhomogeneous material parameters (graded interfaces) (deferred) (done 2026-07-05) — shipped: src/material/disorder.rs (GradedInterface::generate, GradingLaw, test: test_graded_interface_barycenter_symmetric_profile)
  - **Goal:** A GradedInterface model producing spatially-varying Ms/A/K material-parameter profiles across an interface region (linear, exponential, and error-function grading laws), sampled on a simulation grid.
  - **Design:** New struct in src/material/disorder.rs following the sibling generate(...) -> Result<Self> + seeded Xorshift64 pattern used by SurfaceRoughness/GrainStructure. Grading law enum (Linear, Exponential, ErrorFunction) parameterized by transition width and endpoint values; a sampling method returning the parameter value at a given position/depth.
  - **Files:** src/material/disorder.rs (extend), src/material/mod.rs (re-export), examples/graded_interface.rs
  - **Prerequisites:** none - sibling disorder infra already shipped
  - **Tests:** monotonic profile across the transition; correct endpoint values at the domain boundary; correct mean/barycenter for symmetric profiles; invalid-parameter rejection (negative width, non-finite bounds)
  - **Risk:** low - purely additive, follows an established sibling pattern (SurfaceRoughness/GrainStructure) closely; no shared-code touch

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
- [x] Material creation benchmarks (deferred) (done 2026-07-05)
  - **Goal:** Two new criterion benches: material-preset creation cost (Ferromagnet::yig/permalloy/cofeb/iron/cobalt/nickel/cofe) and skyrmion-dynamics simulation cost (SimulationBuilder::skyrmion_dynamics() preset run).
  - **Design:** New benches/material_benchmark.rs and benches/skyrmion_benchmark.rs, harness=false matching existing benches/{llg,spinchain,vector3}_benchmark.rs style; register [[bench]] entries in Cargo.toml.
  - **Files:** benches/material_benchmark.rs (new), benches/skyrmion_benchmark.rs (new), Cargo.toml (+2 [[bench]] entries)
  - **Prerequisites:** none
  - **Tests:** benches must compile cleanly under `cargo bench --no-run` and execute at least one iteration
  - **Risk:** low - purely additive dev-only files, zero production code touched
- [x] Skyrmion dynamics benchmarks (deferred) (done 2026-07-05)
  - **Goal:** Two new criterion benches: material-preset creation cost (Ferromagnet::yig/permalloy/cofeb/iron/cobalt/nickel/cofe) and skyrmion-dynamics simulation cost (SimulationBuilder::skyrmion_dynamics() preset run).
  - **Design:** New benches/material_benchmark.rs and benches/skyrmion_benchmark.rs, harness=false matching existing benches/{llg,spinchain,vector3}_benchmark.rs style; register [[bench]] entries in Cargo.toml.
  - **Files:** benches/material_benchmark.rs (new), benches/skyrmion_benchmark.rs (new), Cargo.toml (+2 [[bench]] entries)
  - **Prerequisites:** none
  - **Tests:** benches must compile cleanly under `cargo bench --no-run` and execute at least one iteration
  - **Risk:** low - purely additive dev-only files, zero production code touched
- [ ] Automated regression alerts (deferred)

### Priority 13: SimulationBuilder
- [x] SimulationBuilder with method chaining — `builder` module
- [x] 8 SolverKind variants: Rk4, Euler, Heun, Dp45, Dp87, Yoshida4, ForestRuth, SemiImplicit
- [x] Preset configurations — `simulation_builder` example
- [x] Type-state pattern compile-time validation (deferred) — shipped: src/builder/mod.rs, states.rs (type-state SimulationBuilder<M,F,S>, compile_fail doctest proof)
- [x] Streaming API for large datasets (deferred) (done 2026-07-05)
  - **Goal:** A streaming/iterator (or callback) simulation-output path that does not retain the full trajectory in memory, for large-N or long-run simulations.
  - **Design:** Extend Simulation/SimulationBuilder (src/builder/mod.rs) with a run_streaming<F: FnMut(usize, &Vector3<f64>, f64)>(&self, callback: F) -> Result<()> method (or a StepIterator implementing Iterator<Item=Result<(f64,Vector3<f64>,f64)>>) that performs the same integration loop as run() but yields/calls back per-step instead of collecting into trajectory/energies Vecs.
  - **Files:** src/builder/mod.rs (extend, additive method only - no breaking change to existing run()/SimulationResult), examples/streaming_simulation.rs
  - **Prerequisites:** none
  - **Tests:** streamed per-step values equal the corresponding entries from the collected run() trajectory on a small cross-check case; memory-behavior test showing no growing Vec accumulation during streaming
  - **Risk:** low-moderate - touches builder/mod.rs (a core, heavily-used file); must be strictly additive (new method(s), zero changes to existing run()/build() signatures or behavior)

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

## v0.7.0 - COMPLETE (2026-05-17)

**Released**: 2026-05-17
**Theme**: ML Phase 2, Advanced Spin Waves, Stiff/Diffusion Integrators, Experimental Validation
**Tests**: 1329 passing, 0 failures, 0 warnings
**Code Size**: ~71K total lines (~58K Rust code)

### ML Enhancements / Phase 2 (`#[cfg(feature = "autodiff")]`)
- [x] `Mlp`, `Layer`, `Activation` enum (`Relu`, `Tanh`, `Sigmoid`, `Gelu`, `Linear`) — feed-forward NN with Xavier/He init (`src/autodiff/neural.rs`)
- [x] `NeuralExchange` — trainable J(r), rescales r to [-1,1] for stable training
- [x] `NeuralAnisotropy` — trainable K(m_x, m_y, m_z) surrogate
- [x] `LlgPinn` / `PinnTrainer` — physics-informed NN for LLG with FD time derivative on tape (`src/autodiff/pinn.rs`)
- [x] `SpinConfig` (spherical coords), `EnergyFunctional`, `MagneticStructureOptimizer`, `find_fm_ground_state`, `find_afm_ground_state` (`src/autodiff/structure_opt.rs`)

### Advanced Spin Wave Theory
- [x] `NanodiskSpinWaves` — radial Bessel × azimuthal modes, mode spectrum, group velocity, propagation length (`src/spinwave/nanodisk.rs`)
- [x] `MagnonicCrystal1D` / `MagnonicCrystal2D` — plane-wave band structure, band gap, group velocity (`src/spinwave/magnonic_crystal.rs`)
- [x] `SemiInfiniteDamonEshbach` — single-surface DE in semi-infinite media (`src/spinwave/semi_infinite_de.rs`)

### Stiff/Diffusion Integrators
- [x] `ImplicitMidpointNewton` — A-stable 2nd-order with finite-diff Jacobian + Gauss elimination (`src/dynamics/integrators/implicit_midpoint.rs`)
- [x] `CrankNicolsonDiffusion` — unconditionally stable, Dirichlet/Neumann/Periodic BC via Thomas algorithm (`src/dynamics/integrators/crank_nicolson.rs`)
- [x] `SpinDiffusionCrankNicolson` — specialized for ∂μ_s/∂t = D∇²μ_s − μ_s/τ_sf
- [x] Resolves v0.3.0-deferred items: implicit midpoint with Newton iteration; Crank-Nicolson for diffusion-dominated systems

### Experimental Validation Refactor
- [x] `src/validation.rs` → `src/validation/` directory (backward-compatible via `pub use parameter_checks::*`)
- [x] `Demidov2006Validation` — DE dispersion + non-reciprocity vs PRL 96, 097202 (2006)
- [x] `Saitoh2006Validation` — Pt spin Hall angle + ISHE polarity + linear scaling vs APL 88, 182509 (2006)
- [x] `Uchida2008Validation` — LSSE linear thermal response + polarity vs Nature 455, 778 (2008)
- [x] `ValidationResult` type with max/mean relative error, n_points, tolerance, pass flag

### Examples (6 new, total 47)
- [x] `neural_exchange_training.rs` — NeuralExchange + Adam (FD gradient demo)
- [x] `nanodisk_modes.rs` — YIG 100 nm disk fundamental at ~2.2 GHz
- [x] `implicit_midpoint_stiff_demo.rs` — Stiff LLG comparison vs explicit Euler
- [x] `pinn_llg_solver.rs` — LlgPinn on Larmor precession
- [x] `magnonic_crystal_bandgap.rs` — 1D NiFe/CoFeB crystal (~57% relative band gap)
- [x] `experimental_validation_demo.rs` — Run all 3 landmark paper validations

---

## v0.8.0 - COMPLETE (2026-05-17)

**Released**: 2026-05-17
**Theme**: Stochastic Methods, Advanced ML Phase 3, More Validations, Property-Based Testing
**Tests**: 1403 lib + 27 proptest + 103 doctests passing, 0 failures, 0 warnings
**Code Size**: ~76K total lines (~63K Rust code)

### Stochastic Methods Improvements (`#[cfg(feature = "scirs2")]`)
- [x] `HeunAdaptive` — Heun-Euler embedded pair with PI controller, frozen-noise rejection retry (`src/stochastic/heun_adaptive.rs`)
- [x] `ImplicitMilstein` — Newton + 3×3 FD Jacobian + optional Milstein correction for multiplicative noise (`src/stochastic/implicit_milstein.rs`)
- [x] `PimcSimulation` — worldline path-integral MC for finite-T Heisenberg chains, Chain1D / Ring1D lattices, Trotter rigidity coupling, Marsaglia spin proposals (`src/stochastic/pimc.rs`)

### Advanced ML Phase 3 (`#[cfg(feature = "autodiff")]`)
- [x] `EquivariantLinear`, `EquivariantMlp`, `EquivariantConfig` — Cartesian-tensor O(3)-equivariant NN layers; rotation invariance preserved to machine precision (`src/autodiff/equivariant.rs`)
- [x] `random_so3`, `rotate_vector` — Marsaglia + Rodrigues helpers for testing equivariance
- [x] `ActiveLearner`, `ActiveLearningConfig`, `QueryStrategy` (UncertaintySampling / QueryByCommittee / RandomBaseline), `ActiveLearnResult` — active learning with ensemble bootstrap and `fit(oracle, pool)` loop (`src/autodiff/active_learning.rs`)

### More Experimental Validations
- [x] `Mosendz2010Validation` — Mosendz et al. PRL 104, 046601 (2010): V_ISHE vs Pt thickness, Δα_eff linewidth enhancement, g↑↓ ≈ 2.1×10¹⁹ m⁻², refined θ_SH ≈ 0.013 (`src/validation/experimental/mosendz_2010.rs`)
- [x] `Liu2012Validation` — Liu et al. Science 336, 555 (2012): β-Ta θ_SH ≈ -0.12 (negative!), critical J_c, polarity, thickness scaling (`src/validation/experimental/liu_2012.rs`)

### Property-Based Testing (first `tests/` integration suite, `proptest = "1.6"`)
- [x] `tests/property_conservation.rs` — 12 properties × 32 cases: |m|=1 under RK4/Heun/Euler, Zeeman energy at α=0, damping alignment, Larmor sign reversal, dm/dt orthogonality, cross-product anti-commutativity, triple product cyclic, Lagrange identity, normalize idempotent
- [x] `tests/property_symmetries.rs` — 15 properties × 32 cases: SO(3) matrix orthogonality + determinant + trace bounded, SO(3) invariants (dot, magnitude), cross-product equivariance, exchange/Zeeman/anisotropy invariance, calc_dm_dt equivariance, time-reversal at α=0, linearity in H, parity, rotation composition + inverse
- [x] Total: 27 property tests × 32 cases = 864 randomized trials per `cargo test` run

### Examples (5 new, total 52)
- [x] `heun_adaptive_thermal_llg.rs` — Permalloy at T=300 K, adaptive dt 96–276 fs
- [x] `pimc_heisenberg_chain.rs` — 1D Heisenberg ring vs β: paramagnetic → aligned
- [x] `equivariant_nn_demo.rs` — 4-spin EquivariantMlp; rotation drift 4.4×10⁻¹⁶
- [x] `active_learning_demo.rs` — QueryByCommittee 24× better than RandomBaseline on sin(5x)·exp(-x²)
- [x] `validation_landmark_suite.rs` — 5 papers / 12 of 15 quantitative checks pass

---

## v0.9.0 - COMPLETE (2026-05-17)

**Released**: 2026-05-17
**Theme**: ML Phase 4 (Graph NN + Bayesian Opt), Garello 2013 + Boona 2014 validations, GPU device abstraction skeleton
**Tests**: 1463 lib + 27 proptest + 103 doctests passing, 0 failures, 0 warnings
**Code Size**: ~80K total lines (~67K Rust code)

### Advanced ML Phase 4 (`#[cfg(feature = "autodiff")]`)
- [x] `LatticeGraph`, `GraphMessagePassingLayer`, `GraphMlp`, `NodeFeatures` — equivariant graph NN message passing on arbitrary lattice topology; chain_1d, ring_1d, square_lattice_2d builders; rotation invariance to 2.8e-14 (`src/autodiff/graph_nn.rs`)
- [x] `GaussianProcess`, `GpConfig`, `BayesianOptimizer`, `BayesianOptConfig`, `AcquisitionStrategy` (EI / UCB / PosteriorVariance), `BayesianOptResult` — BO with RBF-kernel GP, Cholesky + jitter, custom erf via Abramowitz-Stegun (`src/autodiff/bayesian_opt.rs`)

### More Experimental Validations (total: 7 papers)
- [x] `Garello2013Validation` — Nat. Nanotechnol. 8, 587 (2013): angular-harmonic SOT decomposition Pt/Co/AlOx
- [x] `Boona2014Validation` — MRS Bulletin 39, 426 (2014): LSSE in granular YIG/Pt
- Previously: Demidov 2006, Saitoh 2006, Uchida 2008 (v0.7.0); Mosendz 2010, Liu 2012 (v0.8.0)

### GPU Acceleration Skeleton (feature `cuda`)
- [x] `Device` trait — object-safe `Box<dyn Device>` abstraction (`src/gpu/mod.rs`)
- [x] `CpuDevice` — always available; wraps existing CPU LLG with frozen-field RK4 (~92 ns/spin/step) (`src/gpu/cpu.rs`)
- [x] `CudaDevice` — `#[cfg(feature = "cuda")]` skeleton; constructs OK with available=false; all ops return numerical_error pending v1.0.0 CUDA kernels (`src/gpu/cuda.rs`)
- [x] `available_devices()`, `select_best_device()` — auto-enumeration
- [x] Feature `cuda = []` — no external deps yet (pure plumbing)

### Examples (4 new, total 56)
- [x] `graph_nn_lattice.rs` — rotation invariance to 2.8e-14
- [x] `bayesian_opt_materials.rs` — BO converges within 0.026 of true optimum in 20 evals
- [x] `validation_full_suite.rs` — 18/23 checks across 7 landmark papers
- [x] `gpu_device_demo.rs` — 100 spins × 200 steps in 1.8 ms; scaling sweep

---

## v0.5.0 - COMPLETE (2026-05-31)

**Released**: 2026-05-31
**Theme**: SMR/STNO/AOS/SAW Physics, ML Phase 5, RL, Micromagnetics, New Validations
**Tests**: 1689 lib + 42 proptest + 109 doctests passing, 0 failures, 0 warnings
**Code Size**: ~95K total lines (~80K Rust code)

### New Physics Effect Modules
- [x] `SpinHallMagnetoresistance`, `UnidirectionalSmr` — SMR/USMR with Chen (PRB 2013) formula, angular scans, Pt/YIG & W/YIG & Pt/Co presets (`src/effect/smr.rs`)
- [x] `SpinTorqueOscillator`, `SpinTorqueOscillatorConfig` — Slonczewski STT auto-oscillation, threshold current, Slavin–Tiberkevich linewidth, Adler locking, Permalloy preset (`src/effect/stno.rs`)
- [x] `CircularHelicity`, `LaserPulseParams`, `OpticalMagneticMaterial`, `OpticalSwitching`, `OpticalSwitchResult` — Inverse Faraday Effect, HDS, ultrafast demag model; GdFeCo, Co, Py presets (`src/effect/optical_switching.rs`)
- [x] `AcSpinPumping`, `SpinBattery` — backflow-corrected G_r_eff; DC/2ω spin current at FMR; Δα enhancement; ISHE voltage; YIG/Pt preset (`src/transport/ac_pumping.rs`)
- [x] `PiezoSubstrate`, `SawMagnetoelastic`, `SawSource`, `SawMagnetoacoustics`, `SawSpinWaveExcitation` — SAW magnetoacoustics; resonant precession (Lorentzian); acoustic spin pumping; LiNbO₃/GaAs/ZnO presets (`src/mech/saw.rs`)

### ML Phase 5 (`#[cfg(feature = "autodiff")]`)
- [x] `DiffusionModel`, `NoiseSchedule`, `SpinTexture` — DDPM with 2-layer MLP denoiser + manual backprop; Adam training; reverse diffusion; topological charge validation (`src/autodiff/diffusion_model.rs`)
- [x] `QuantumClassicalOptimizer`, `MagnonNeuralNetwork`, `MagnonHamiltonianParams`, `QuantumClassicalResult` — MLP → Bogoliubov (A_k, B_k) → ε_k = √(A_k²−B_k²); central FD gradients; Adam training (`src/autodiff/quantum_classical.rs`)

### RL for SOT Switching
- [x] `SotSwitchingEnv`, `CemPolicy`, `SotRlOptimizer`, `SotRlResult`, `SotSwitchingConfig` — PMA macrospin LLG+SOT environment; CEM for pulse optimization; CoFeB/Pt preset (`src/ai/rl.rs`)

### Micromagnetics Infrastructure
- [x] `NewellTensor`, `DemagField` — analytic Newell (1993) demag; 8-corner f-function; direct O(N²) convolution (`src/micromagnetics/demag.rs`)
- [x] `MicromagneticGrid`, `GridConfig`, `LlgResult` — FD exchange + demag + Zeeman + anisotropy; per-cell LLG RK4 (`src/micromagnetics/grid.rs`)
- [x] `StandardProblem3`, `Sp3Config`, `Sp3Result`, `StableState` — muMAG SP#3 flower↔vortex energy comparison (`src/validation/standard_problems/sp3.rs`)

### Experimental Validations (total: 11 papers)
- [x] `Nakayama2013Validation` — Pt/YIG SMR angular scan + ratio vs PRL 110, 206601 (2013) (`src/validation/experimental/nakayama_2013.rs`)
- [x] `Avci2015Validation` — Pt/Co USMR current linearity + coefficient vs Nat. Phys. 11, 570 (2015) (`src/validation/experimental/avci_2015.rs`)
- [x] `Woo2016Validation` — skyrmion diameter in Pt/CoFeB/MgO vs Nat. Mater. 15, 501 (2016) (`src/validation/experimental/woo_2016.rs`)
- [x] `Cornelissen2015Validation` — nonlocal magnon transport, λ_m=9.4 μm vs Nat. Phys. 11, 1022 (2015) (`src/validation/experimental/cornelissen_2015.rs`)

### Property-Based Testing Phase 2
- [x] `tests/property_physics.rs` — 15 property tests × 32 cases = 480 trials: SMR symmetries, USMR antisymmetry, STNO norm conservation + torque perpendicularity, AC pumping inequalities

### Examples (6 new, total 62)
- [x] `smr_angular_scan.rs` — SMR angular scan + USMR + Nakayama/Avci validations
- [x] `stno_auto_oscillation.rs` — STNO auto-oscillation, linewidth, injection locking
- [x] `ac_spin_pumping_battery.rs` — YIG/Pt spin pumping DC/2ω + ISHE voltage
- [x] `rl_sot_switching.rs` — CEM RL agent for SOT pulse protocol
- [x] `diffusion_skyrmion_gen.rs` — DDPM training + topological charge (autodiff)
- [x] `variational_magnon_nn.rs` — quantum-classical hybrid NN (autodiff)

---

## v1.0.0 - ROADMAP (Stable Release)

**Target**: Q3 2027
**Theme**: Stable API, Workspace Split, Real CUDA Kernels, Language Bindings

### Language Bindings (FFI foundation)
- [ ] C/C++ bindings via cbindgen (foundation for FFI to other languages)
- [ ] Build script `build.rs` for automatic `spintronics.h` generation
- [ ] Example C program calling LLG solver + Vector3
- [ ] Julia bindings via julia-rs or jlrs (builds on C ABI)
- [ ] R bindings via Rcpp (builds on C ABI)

### Workspace Restructuring
- [ ] Split monolithic crate into workspace with multiple crates:
  - `spintronics-core` — Core physics and materials
  - `spintronics-solver` — Numerical solvers (LLG, integrators)
  - `spintronics-spinwave` — Spin wave theory + magnonics
  - `spintronics-topomagnon` — Topological magnon bands + HOTI + axion
  - `spintronics-autodiff` — ML autodiff + neural potentials + PINN + equivariant + graph NN
  - `spintronics-io` — I/O and visualization (VTK, HDF5, NetCDF, Zarr)
  - `spintronics-gpu` — Device abstraction + CUDA / ROCm kernels
  - `spintronics-python` — Python bindings
  - `spintronics-c` — C bindings
  - `spintronics-cli` — Command-line tools

### Real GPU Kernels (drop the v0.9.0 skeleton)
- [ ] CUDA kernels via `cudarc` crate
- [ ] LLG RK4 kernel (one block per spin)
- [ ] Zeeman / exchange / DMI effective field kernels
- [ ] Async device-host transfers with stream synchronisation
- [ ] ROCm equivalent via `rocm-rs`
- [ ] Benchmark vs OOMMF, mumax3 (target 10–100× speedup on >1M spins)

### API Stabilisation
- [ ] Comprehensive API review for v1.0
- [ ] Migration guide for v0.x → v1.0 breaking changes
- [ ] semver-check in CI
- [ ] Long-term support announcements

### Advanced ML Phase 5
- [x] Hybrid quantum-classical NN for magnetic Hamiltonians — shipped: src/autodiff/quantum_classical.rs (test: test_training_reduces_loss)
- [x] Diffusion models for skyrmion lattice generation — shipped: src/autodiff/diffusion_model.rs (test: test_topological_charge_skyrmion)
- [x] Reinforcement learning for SOT switching protocol design — shipped: src/ai/rl.rs (test: test_cem_policy_converges)

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
- [x] Fallback to CPU for systems without GPU — shipped: src/gpu/cpu.rs (CpuDevice, test: test_multi_step_runs_and_norm_preserved)
- [x] Unified API: transparent GPU/CPU selection — shipped: src/gpu/mod.rs (select_best_device())

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
- [x] Property-based testing with proptest (conservation laws, symmetries) — shipped: tests/property_conservation.rs, property_symmetries.rs, property_physics.rs, property_optics_saw.rs (58 proptest properties total)
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
- [x] Spin Pumping: Validate against Saitoh et al. 2006 APL — shipped: src/validation/experimental/saitoh_2006.rs (test: test_ishe_scaling_linear)
- [x] Spin Seebeck Effect: Compare with Uchida et al. 2008 Nature — shipped: src/validation/experimental/uchida_2008.rs (test: test_seebeck_coefficient_in_order_window)
- [x] Skyrmion Size: Validate with Woo 2016 Nat. Mater. — shipped: src/validation/experimental/woo_2016.rs (test: test_validation_passes_30pct)
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
| v0.3.1 | 2026-06-10 | DemagField optimization, hamiltonian_at Result, scirs2 0.5.0 | COMPLETE |
| v0.3.2 | TBD | In development | IN PROGRESS |
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
**License**: Apache-2.0
**Repository**: https://github.com/cool-japan/spintronics
**Contact**: See CONTRIBUTING.md for communication channels
