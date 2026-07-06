//! Resonating valence bond (RVB) quantum spin-liquid module.
//!
//! # Physics background
//!
//! On a geometrically frustrated (or simply strongly quantum) antiferromagnet, the
//! classical Néel-ordered state competes with a purely quantum-mechanical alternative:
//! a superposition of many ways of pairing spins into singlet "valence bonds", first
//! proposed by Anderson as the ground state of the triangular-lattice Heisenberg
//! antiferromagnet (P.W. Anderson, Mater. Res. Bull. 8, 153 (1973)) and developed for
//! anisotropic triangular lattices by Fazekas and Anderson (P. Fazekas, P.W. Anderson,
//! Philos. Mag. 30, 423 (1974)). Such a **resonating valence bond (RVB)** state has no
//! long-range magnetic order at any temperature (a "quantum spin liquid") and supports
//! fractionalized, deconfined spin-1/2 **spinon** excitations rather than ordinary
//! spin-1 magnons.
//!
//! This module implements the standard machinery for putting RVB variational states on
//! a rigorous, checkable footing:
//!
//! - **[`dimer`]** — combinatorial enumeration of dimer (perfect-matching) and
//!   monomer-dimer coverings of a bond graph, restricted to nearest-neighbor valence
//!   bonds (the standard tractable truncation of Liang, Doucot & Anderson,
//!   Phys. Rev. Lett. 61, 365 (1988)).
//! - **[`valence_bond`]** — the non-orthogonal overlap between two valence-bond states
//!   via Sutherland's loop-counting formula (B. Sutherland, Phys. Rev. B 37, 3786
//!   (1988)), and the Heisenberg Hamiltonian's matrix elements in this basis.
//! - **[`solver`]** — assembly of the overlap (`S`) and Hamiltonian (`H`) matrices and
//!   the generalized eigenproblem `Hc = ESc` (Löwdin canonical orthogonalization),
//!   giving the variational RVB ground state.
//! - **[`exact`]** — independent exact-diagonalization benchmarks (dense, and
//!   matrix-free Lanczos in the `S_z=0` sector) that the variational solver is checked
//!   against.
//! - **[`spinon`]** — energetics of a separated spinon pair and a deconfinement
//!   diagnostic.
//! - **[`diagnostics`]** — real-space spin correlations, the static structure factor,
//!   and a long-range-order screening report.
//!
//! # The `CMatrix::MAX_DIM = 64` constraint
//!
//! [`crate::math::CMatrix`] is the only Hermitian eigensolver available in this crate,
//! and it is capped at `MAX_DIM = 64`. This drives the whole module's design:
//!
//! - The valence-bond basis (dimension = number of enumerated dimer coverings) is
//!   capped at [`MAX_VB_BASIS`] `= 64`. Enumeration never silently truncates: exceeding
//!   the cap returns [`crate::error::Error::InvalidParameter`] naming the *actual
//!   measured* covering count.
//! - A dense exact-diagonalization benchmark (the full `2^N × 2^N` Hamiltonian as a
//!   `CMatrix`) is only honest up to `N ≤ 6` (`2^6 = 64`; already `2^8 = 256 > 64`), so
//!   [`exact::ExactDiagonalization`] additionally provides a **matrix-free Lanczos**
//!   solver in the `S_z = 0` sector (bit-encoded basis, sparse Hamiltonian-vector
//!   product, full reorthogonalization) valid up to [`MAX_ED_SITES`] `= 16` sites — the
//!   real benchmark for anything interesting-sized.
//! - All matrices here are real symmetric (the Heisenberg model has no spin-orbit or
//!   Dzyaloshinskii-Moriya terms); they are embedded into `CMatrix` with zero imaginary
//!   parts purely to reuse `hermitian_eigendecomposition`.
//!
//! # Quick start
//!
//! ```rust
//! use spintronics::frustrated::rvb::{RvbSolver, ExactDiagonalization};
//!
//! // 4-site Heisenberg ring: exactly solvable, ground energy -2J.
//! let solver = RvbSolver::from_bonds(4, vec![(0, 1), (1, 2), (2, 3), (3, 0)], 1.0)
//!     .expect("valid bond list");
//! let ground = solver.ground_state().expect("generalized eigenproblem should solve");
//! assert!((ground.energy - (-2.0)).abs() < 1e-9);
//!
//! // Cross-check against a fully independent exact-diagonalization benchmark.
//! let ed = ExactDiagonalization::new(4, vec![(0, 1), (1, 2), (2, 3), (3, 0)], 1.0)
//!     .expect("valid ED instance");
//! let e_exact = ed.dense_ground_energy().expect("dense diagonalization");
//! assert!((ground.energy - e_exact).abs() < 1e-8);
//! ```

pub mod diagnostics;
pub mod dimer;
pub mod exact;
pub mod solver;
pub mod spinon;
pub mod valence_bond;

/// Maximum valence-bond basis dimension, bounded by `CMatrix::MAX_DIM = 64`.
///
/// Dimer-covering enumeration never silently truncates: exceeding this cap returns an
/// [`crate::error::Error::InvalidParameter`] naming the actual measured covering count.
pub const MAX_VB_BASIS: usize = 64;

/// Maximum number of sites for the matrix-free Lanczos exact-diagonalization backend
/// (`S_z = 0` sector). Dense `CMatrix`-based exact diagonalization is honest only up to
/// 6 sites (`2^6 = 64 = CMatrix::MAX_DIM`); Lanczos extends the benchmark to sizes where
/// the valence-bond basis is actually interesting.
pub const MAX_ED_SITES: usize = 16;

/// Relative singular-value cutoff for Löwdin canonical orthogonalization: overlap-matrix
/// eigenvalues below `S_SINGULAR_TOL * s_max` are dropped before projecting the
/// Hamiltonian, since the (generally overcomplete, non-orthogonal) valence-bond basis
/// can have near-null overlap eigenvalues.
pub const S_SINGULAR_TOL: f64 = 1e-10;

pub use diagnostics::{is_spin_liquid, spin_correlation, structure_factor, SpinLiquidReport};
pub use dimer::DimerCovering;
pub use exact::{ExactDiagonalization, SpinBasisState};
pub use solver::{GroundState, RvbSolver, ShortRangeRvb};
pub use spinon::{deconfinement_diagnostic, spinon_pair_energy, SpinonSeparationPoint};
pub use valence_bond::ValenceBondState;
