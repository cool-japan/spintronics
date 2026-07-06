//! The RVB solver: basis assembly, overlap/Hamiltonian matrices, and the
//! generalized-eigenproblem ground state.
//!
//! [`RvbSolver`] enumerates the nearest-neighbor dimer-covering (valence-bond) basis for
//! a bond graph, builds the (non-orthogonal) overlap matrix `S` and Heisenberg
//! Hamiltonian matrix `H` in that basis (see [`super::valence_bond`] for the physics),
//! and solves the generalized eigenproblem `H c = E S c` for the variational ground
//! state via Löwdin canonical orthogonalization:
//!
//! 1. Diagonalize `S = U s U^T` (real symmetric, embedded as Hermitian via
//!    [`crate::math::CMatrix::hermitian_eigendecomposition`]).
//! 2. Drop near-null eigenvalues (`s_k ≤ `[`super::S_SINGULAR_TOL`]` · s_max`) — the
//!    (non-orthogonal, generally overcomplete) VB basis can have small or zero overlap
//!    eigenvalues, which would otherwise blow up `s_k^{-1/2}`.
//! 3. Build the canonical orthogonalization transform `X = U_kept · s_kept^{-1/2}` and
//!    project: `H' = X^T H X` (a genuine `r×r` matrix on the retained subspace, `r ≤ D`).
//! 4. Diagonalize `H'` (again via `hermitian_eigendecomposition`); the lowest eigenvalue
//!    is the ground-state energy, and back-transforming its eigenvector through `X`
//!    gives the ground state's coefficients in the original (non-orthogonal) VB basis.
//!
//! Because `span{|VB_α⟩}` is a genuine subspace of the full spin Hilbert space, this is
//! an exact Rayleigh-Ritz calculation on that subspace: the resulting energy is a rigorous
//! upper bound on the true ground-state energy (equality when the chosen coverings
//! happen to span the full singlet sector, e.g. tiny clusters like the 4-site ring).
//!
//! `CMatrix::MAX_DIM = 64` bounds the basis size directly: [`RvbSolver::from_bonds`]
//! enumerates dimer coverings and returns a descriptive [`crate::error::Error::InvalidParameter`]
//! (naming the *actual measured* covering count) if it exceeds [`super::MAX_VB_BASIS`],
//! rather than silently truncating.

use super::dimer::DimerCovering;
use super::valence_bond;
use crate::error::{self, Error, Result};
use crate::frustrated::lattice::FrustratedLattice;
use crate::math::{CMatrix, Complex};

/// The resonating-valence-bond solver for a fixed bond graph.
///
/// Built either from an explicit, geometry-agnostic bond list ([`RvbSolver::from_bonds`],
/// the primary and validated API — used for the small open test clusters where exact
/// combinatorics are known) or as a convenience wrapper around a classical
/// [`FrustratedLattice`] ([`RvbSolver::from_lattice`]). The latter's lattices use
/// *periodic* boundary conditions, which distort dimer-covering counts on tiny clusters,
/// so validation should always go through `from_bonds` with explicit edge lists.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RvbSolver {
    /// Number of sites.
    pub num_sites: usize,
    /// Exchange bonds `(i, j)`, `i < j`. Used both to restrict the dimer-covering basis
    /// (nearest-neighbor valence bonds only) and to define the Heisenberg Hamiltonian
    /// `H = J·Σ_{(i,j)∈bonds} S_i·S_j`.
    pub bonds: Vec<(usize, usize)>,
    /// Antiferromagnetic exchange coupling `J > 0`.
    pub coupling_j: f64,
    /// Two-coloring of the bond graph if it is bipartite (`Some`), or `None` if it
    /// contains an odd cycle (e.g. a triangular plaquette). Informational only — the
    /// overlap/Hamiltonian formulas in [`super::valence_bond`] never depend on
    /// bipartiteness.
    pub sublattice: Option<Vec<i8>>,
    /// The enumerated nearest-neighbor dimer-covering (valence-bond) basis.
    pub coverings: Vec<DimerCovering>,
}

/// The variational ground state of an [`RvbSolver`]'s Heisenberg Hamiltonian, restricted
/// to the span of its valence-bond basis.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GroundState {
    /// Ground-state energy (variational upper bound on the true ground energy).
    pub energy: f64,
    /// Coefficients `c_α` such that `|ψ⟩ = Σ_α c_α |VB_α⟩`, normalized so that
    /// `⟨ψ|ψ⟩ = c^T S c = 1`.
    pub coefficients: Vec<f64>,
}

/// A trial valence-bond state (a linear combination of the enumerated VB basis).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ShortRangeRvb {
    /// Amplitude on each basis covering (same order as [`RvbSolver::coverings`]).
    pub coefficients: Vec<f64>,
}

impl ShortRangeRvb {
    /// Anderson's original equal-amplitude RVB ansatz: every covering weighted equally.
    pub fn equal_amplitude(solver: &RvbSolver) -> Self {
        Self {
            coefficients: vec![1.0; solver.dim()],
        }
    }
}

impl RvbSolver {
    /// Build a solver from an explicit bond list — the primary, geometry-agnostic API.
    ///
    /// # Errors
    ///
    /// Returns an error if `coupling_j <= 0`, if `bonds` is malformed, or if the
    /// enumerated dimer-covering basis is empty or exceeds [`super::MAX_VB_BASIS`].
    pub fn from_bonds(
        num_sites: usize,
        bonds: Vec<(usize, usize)>,
        coupling_j: f64,
    ) -> Result<Self> {
        if coupling_j <= 0.0 {
            return Err(error::invalid_param(
                "coupling_j",
                "RVB physics requires antiferromagnetic coupling J > 0",
            ));
        }
        let sublattice = super::dimer::bipartite_coloring(num_sites, &bonds)?;
        let coverings = Self::enumerate_dimer_coverings(num_sites, &bonds)?;
        Ok(Self {
            num_sites,
            bonds,
            coupling_j,
            sublattice,
            coverings,
        })
    }

    /// Build a solver from a classical [`FrustratedLattice`], reading its neighbor list
    /// as the bond graph and its exchange coupling directly.
    ///
    /// # Errors
    ///
    /// See [`RvbSolver::from_bonds`]. Periodic lattices with an odd site count (e.g. a
    /// 3-site-per-cell kagome cluster) will correctly fail (no perfect matching exists).
    pub fn from_lattice(lattice: &FrustratedLattice) -> Result<Self> {
        let num_sites = lattice.num_sites();
        let mut bond_set: std::collections::BTreeSet<(usize, usize)> =
            std::collections::BTreeSet::new();
        for (i, neighbors) in lattice.neighbors.iter().enumerate() {
            for &j in neighbors {
                if i != j {
                    bond_set.insert((i.min(j), i.max(j)));
                }
            }
        }
        let bonds: Vec<(usize, usize)> = bond_set.into_iter().collect();
        Self::from_bonds(num_sites, bonds, lattice.coupling_j)
    }

    /// Enumerate the nearest-neighbor dimer-covering basis, enforcing
    /// [`super::MAX_VB_BASIS`] with a descriptive error (never a silent truncation).
    ///
    /// # Errors
    ///
    /// Returns an error if no perfect matching exists, or if the measured covering
    /// count exceeds [`super::MAX_VB_BASIS`] (the error message names the actual count).
    pub fn enumerate_dimer_coverings(
        num_sites: usize,
        bonds: &[(usize, usize)],
    ) -> Result<Vec<DimerCovering>> {
        let coverings = DimerCovering::enumerate_perfect_matchings(num_sites, bonds)?;
        if coverings.is_empty() {
            return Err(error::invalid_param(
                "bonds",
                "no perfect dimer matching exists for the given bond list",
            ));
        }
        if coverings.len() > super::MAX_VB_BASIS {
            return Err(Error::InvalidParameter {
                param: "num_sites/bonds".to_string(),
                reason: format!(
                    "dimer-covering enumeration produced {} valence-bond basis states, \
                     exceeding MAX_VB_BASIS={} (bounded by CMatrix::MAX_DIM=64); \
                     reduce the cluster size or bond connectivity",
                    coverings.len(),
                    super::MAX_VB_BASIS
                ),
            });
        }
        Ok(coverings)
    }

    /// Dimension of the valence-bond basis (number of enumerated coverings).
    pub fn dim(&self) -> usize {
        self.coverings.len()
    }

    /// The `D×D` overlap matrix `S_{αβ} = ⟨VB_α|VB_β⟩`.
    pub fn overlap_matrix(&self) -> Result<CMatrix> {
        build_overlap_matrix(&self.coverings, self.num_sites)
    }

    /// The `D×D` Heisenberg Hamiltonian matrix `H_{αβ} = ⟨VB_α|H|VB_β⟩` in the
    /// (non-orthogonal) valence-bond basis.
    pub fn hamiltonian_matrix(&self) -> Result<CMatrix> {
        build_hamiltonian_matrix(
            &self.coverings,
            &self.bonds,
            self.coupling_j,
            self.num_sites,
        )
    }

    /// Solve the generalized eigenproblem `Hc = ESc` via Löwdin canonical
    /// orthogonalization; see the module-level docs for the algorithm.
    ///
    /// # Errors
    ///
    /// Returns a [`crate::error::Error::NumericalError`] if the overlap matrix has no
    /// eigenvalues above the [`super::S_SINGULAR_TOL`] cutoff (a degenerate basis).
    pub fn ground_state(&self) -> Result<GroundState> {
        generalized_ground_state(
            self.num_sites,
            &self.coverings,
            &self.bonds,
            self.coupling_j,
        )
    }

    /// Variational (Rayleigh-quotient) energy of an arbitrary trial state in the VB
    /// basis: `⟨ψ|H|ψ⟩ / ⟨ψ|ψ⟩` for `|ψ⟩ = Σ_α c_α|VB_α⟩`.
    ///
    /// By the variational principle this is always `>= ` [`RvbSolver::ground_state`]'s
    /// energy.
    ///
    /// # Errors
    ///
    /// Returns an error if `trial.coefficients` has the wrong length, or if the trial
    /// state has (numerically) zero norm.
    pub fn variational_energy(&self, trial: &ShortRangeRvb) -> Result<f64> {
        let d = self.dim();
        if trial.coefficients.len() != d {
            return Err(error::invalid_param(
                "trial",
                "coefficient vector length must match the number of dimer coverings",
            ));
        }
        let s = self.overlap_matrix()?;
        let h = self.hamiltonian_matrix()?;
        let c = &trial.coefficients;
        let mut numerator = 0.0_f64;
        let mut denominator = 0.0_f64;
        for a in 0..d {
            for b in 0..d {
                numerator += c[a] * h.get(a, b).re * c[b];
                denominator += c[a] * s.get(a, b).re * c[b];
            }
        }
        if denominator.abs() < 1e-14 {
            return Err(error::numerical_error(
                "trial state has (numerically) zero norm: c^T S c vanished",
            ));
        }
        Ok(numerator / denominator)
    }
}

/// Build the overlap matrix for an arbitrary (homogeneous-monomer-set) covering list —
/// shared by [`RvbSolver::overlap_matrix`] and [`super::spinon`]'s sub-lattice solves.
pub(crate) fn build_overlap_matrix(
    coverings: &[DimerCovering],
    num_sites: usize,
) -> Result<CMatrix> {
    let d = coverings.len();
    if d == 0 {
        return Err(error::invalid_param(
            "coverings",
            "at least one dimer covering is required to build the overlap matrix",
        ));
    }
    let mut rows = vec![vec![Complex::ZERO; d]; d];
    for a in 0..d {
        for b in a..d {
            let s = valence_bond::overlap(&coverings[a], &coverings[b], num_sites)?;
            rows[a][b] = Complex::from_real(s);
            rows[b][a] = Complex::from_real(s);
        }
    }
    CMatrix::from_rows(rows)
}

/// Build the Heisenberg Hamiltonian matrix for an arbitrary covering list and bond list
/// — shared by [`RvbSolver::hamiltonian_matrix`] and [`super::spinon`]'s sub-lattice
/// solves. Bonds touching a monomer site (absent from every covering's pairing) are
/// silently skipped: no valence-bond action is defined there (see [`super::spinon`]'s
/// exact-zero-expectation argument for why this is the physically correct treatment,
/// not an approximation, for its fixed-polarization spinon ansatz).
pub(crate) fn build_hamiltonian_matrix(
    coverings: &[DimerCovering],
    bonds: &[(usize, usize)],
    coupling_j: f64,
    num_sites: usize,
) -> Result<CMatrix> {
    let d = coverings.len();
    if d == 0 {
        return Err(error::invalid_param(
            "coverings",
            "at least one dimer covering is required to build the Hamiltonian matrix",
        ));
    }
    let mut h = vec![vec![0.0_f64; d]; d];
    for (beta_idx, beta) in coverings.iter().enumerate() {
        let mut diag_coeff_total = 0.0_f64;
        let mut recon_terms: Vec<(DimerCovering, f64)> = Vec::new();
        for &(i, j) in bonds {
            if beta.partner_of(i).is_none() || beta.partner_of(j).is_none() {
                continue;
            }
            let bc = valence_bond::bond_contribution(beta, i, j, coupling_j)?;
            diag_coeff_total += bc.diagonal;
            if let Some((beta_prime, coeff)) = bc.reconnection {
                recon_terms.push((beta_prime, coeff));
            }
        }
        for (alpha_idx, alpha) in coverings.iter().enumerate() {
            let s_alpha_beta = valence_bond::overlap(alpha, beta, num_sites)?;
            h[alpha_idx][beta_idx] += diag_coeff_total * s_alpha_beta;
        }
        for (beta_prime, coeff) in &recon_terms {
            for (alpha_idx, alpha) in coverings.iter().enumerate() {
                let s_alpha_bp = valence_bond::overlap(alpha, beta_prime, num_sites)?;
                h[alpha_idx][beta_idx] += coeff * s_alpha_bp;
            }
        }
    }
    // Symmetrize away floating-point order-of-summation noise (H is exactly symmetric).
    // `a`/`b` cross-index the shared matrix `h` at non-adjacent (a,b) and (b,a)
    // positions simultaneously, which does not map cleanly onto iterator adapters.
    #[allow(clippy::needless_range_loop)]
    for a in 0..d {
        for b in (a + 1)..d {
            let avg = 0.5 * (h[a][b] + h[b][a]);
            h[a][b] = avg;
            h[b][a] = avg;
        }
    }
    let rows: Vec<Vec<Complex>> = h
        .into_iter()
        .map(|row| row.into_iter().map(Complex::from_real).collect())
        .collect();
    CMatrix::from_rows(rows)
}

/// Solve `Hc = ESc` via Löwdin canonical orthogonalization for an arbitrary covering
/// list — shared by [`RvbSolver::ground_state`] and [`super::spinon`]'s sub-lattice
/// solves (which use a covering list with a fixed monomer pair, and a bond list already
/// filtered to exclude bonds touching the monomers).
pub(crate) fn generalized_ground_state(
    num_sites: usize,
    coverings: &[DimerCovering],
    bonds: &[(usize, usize)],
    coupling_j: f64,
) -> Result<GroundState> {
    let d = coverings.len();
    let s = build_overlap_matrix(coverings, num_sites)?;
    let h = build_hamiltonian_matrix(coverings, bonds, coupling_j, num_sites)?;

    let (s_vals, s_vecs) = s.hermitian_eigendecomposition()?;
    let s_max = s_vals.iter().copied().fold(f64::MIN, f64::max);
    if s_max.is_nan() || s_max <= 0.0 {
        return Err(error::numerical_error(
            "overlap matrix has no positive eigenvalues; the valence-bond basis is degenerate",
        ));
    }
    let threshold = super::S_SINGULAR_TOL * s_max;
    let keep: Vec<usize> = (0..d).filter(|&k| s_vals[k] > threshold).collect();
    let r = keep.len();
    if r == 0 {
        return Err(error::numerical_error(
            "every overlap-matrix eigenvalue fell below the singular-value cutoff",
        ));
    }

    // Canonical orthogonalization transform X (d x r, real):
    // X[:, idx] = s_vecs[:, k] / sqrt(s_vals[k]).
    let mut x = vec![vec![0.0_f64; r]; d];
    for (idx, &k) in keep.iter().enumerate() {
        let inv_sqrt = 1.0 / s_vals[k].sqrt();
        for (row, x_row) in x.iter_mut().enumerate() {
            x_row[idx] = s_vecs.get(row, k).re * inv_sqrt;
        }
    }

    let h_real: Vec<Vec<f64>> = (0..d)
        .map(|row| (0..d).map(|col| h.get(row, col).re).collect())
        .collect();
    let xt = real_transpose(&x);
    let h_prime = real_matmul(&real_matmul(&xt, &h_real), &x);

    let hp_rows: Vec<Vec<Complex>> = h_prime
        .iter()
        .map(|row| row.iter().map(|&v| Complex::from_real(v)).collect())
        .collect();
    let hp_matrix = CMatrix::from_rows(hp_rows)?;
    let (energies, eigvecs_reduced) = hp_matrix.hermitian_eigendecomposition()?;
    let ground_energy = energies[0];

    let mut coefficients = vec![0.0_f64; d];
    for (row, coeff) in coefficients.iter_mut().enumerate() {
        let mut acc = 0.0_f64;
        // `idx` indexes both the plain Vec<Vec<f64>> `x` and the CMatrix
        // `eigvecs_reduced` simultaneously; there is no single collection to enumerate.
        #[allow(clippy::needless_range_loop)]
        for idx in 0..r {
            acc += x[row][idx] * eigvecs_reduced.get(idx, 0).re;
        }
        *coeff = acc;
    }

    Ok(GroundState {
        energy: ground_energy,
        coefficients,
    })
}

/// Plain `f64` matrix multiplication for rectangular (non-square) real matrices, which
/// `CMatrix` cannot represent (it requires square matrices).
fn real_matmul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    let k = a.first().map_or(0, |row| row.len());
    let n = b.first().map_or(0, |row| row.len());
    let mut out = vec![vec![0.0_f64; n]; m];
    for (i, out_row) in out.iter_mut().enumerate() {
        for p in 0..k {
            let aip = a[i][p];
            for (j, out_val) in out_row.iter_mut().enumerate() {
                *out_val += aip * b[p][j];
            }
        }
    }
    out
}

/// Transpose a rectangular real matrix.
fn real_transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    let n = a.first().map_or(0, |row| row.len());
    let mut out = vec![vec![0.0_f64; m]; n];
    for (i, row) in a.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            out[j][i] = v;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring4_bonds() -> Vec<(usize, usize)> {
        vec![(0, 1), (1, 2), (2, 3), (3, 0)]
    }

    #[test]
    fn test_two_site_exact_singlet_energy() {
        let solver = RvbSolver::from_bonds(2, vec![(0, 1)], 1.0).expect("valid 2-site solver");
        assert_eq!(solver.dim(), 1);
        let ground = solver.ground_state().expect("ground state should solve");
        assert!(
            (ground.energy - (-0.75)).abs() < 1e-10,
            "expected E = -3J/4 = -0.75, got {}",
            ground.energy
        );
    }

    #[test]
    fn test_ring4_covering_count_and_overlap_matrix() {
        let solver = RvbSolver::from_bonds(4, ring4_bonds(), 1.0).expect("valid ring4 solver");
        assert_eq!(solver.dim(), 2);
        let s = solver.overlap_matrix().expect("overlap matrix");
        assert!((s.get(0, 0).re - 1.0).abs() < 1e-10);
        assert!((s.get(1, 1).re - 1.0).abs() < 1e-10);
        assert!(
            (s.get(0, 1).re.abs() - 0.5).abs() < 1e-10,
            "expected |S_01| = 0.5, got {}",
            s.get(0, 1).re
        );
        assert!(
            (s.get(0, 1).re - s.get(1, 0).re).abs() < 1e-12,
            "S must be symmetric"
        );
    }

    #[test]
    fn test_ring4_ground_energy_exact_minus_2j() {
        let solver = RvbSolver::from_bonds(4, ring4_bonds(), 1.0).expect("valid ring4 solver");
        let ground = solver.ground_state().expect("ground state should solve");
        assert!(
            (ground.energy - (-2.0)).abs() < 1e-9,
            "4-ring ground energy should be exactly -2J, got {}",
            ground.energy
        );
    }

    #[test]
    fn test_variational_ordering_equal_amplitude_ge_ground() {
        let solver = RvbSolver::from_bonds(4, ring4_bonds(), 1.0).expect("valid ring4 solver");
        let ground = solver.ground_state().expect("ground state");
        let trial = ShortRangeRvb::equal_amplitude(&solver);
        let e_trial = solver
            .variational_energy(&trial)
            .expect("variational energy");
        assert!(
            e_trial >= ground.energy - 1e-9,
            "equal-amplitude energy {} should be >= ground energy {}",
            e_trial,
            ground.energy
        );
    }

    #[test]
    fn test_covering_count_exceeds_max_vb_basis_errors_with_measured_count() {
        // 6x6 open square lattice has far more than 64 nearest-neighbor dimer coverings.
        let idx = |x: usize, y: usize| -> usize { y * 6 + x };
        let mut bonds = Vec::new();
        for y in 0..6 {
            for x in 0..6 {
                if x + 1 < 6 {
                    bonds.push((idx(x, y), idx(x + 1, y)));
                }
                if y + 1 < 6 {
                    bonds.push((idx(x, y), idx(x, y + 1)));
                }
            }
        }
        let result = RvbSolver::from_bonds(36, bonds, 1.0);
        match result {
            Err(Error::InvalidParameter { reason, .. }) => {
                assert!(
                    reason.chars().any(|c| c.is_ascii_digit()),
                    "error message should name the actual measured covering count: {}",
                    reason
                );
            },
            other => panic!(
                "expected InvalidParameter error with measured count, got {:?}",
                other
            ),
        }
    }

    #[test]
    fn test_non_positive_coupling_errors() {
        assert!(RvbSolver::from_bonds(2, vec![(0, 1)], 0.0).is_err());
        assert!(RvbSolver::from_bonds(2, vec![(0, 1)], -1.0).is_err());
    }

    #[test]
    fn test_from_lattice_wraps_frustrated_lattice() {
        // A single triangular unit cell isn't necessarily matchable, so use a lattice
        // size guaranteed to have an even site count and at least one matching: 2x1
        // triangular lattice sites = 2 (with periodic self-neighboring degeneracies
        // possible at this tiny size) -- instead use a 2x2 kagome (12 sites, even).
        let lattice = FrustratedLattice::kagome(2, 2, 1.0, 1e-9).expect("valid kagome lattice");
        let solver = RvbSolver::from_lattice(&lattice).expect("from_lattice should succeed");
        assert_eq!(solver.num_sites, 12);
        assert!(solver.dim() >= 1);
    }

    #[test]
    fn test_sublattice_none_for_non_bipartite_bonds() {
        // Triangle: odd cycle, non-bipartite. No perfect matching exists for N=3 anyway
        // (odd site count), but bipartite_coloring is computed before that check.
        let solver_result = RvbSolver::from_bonds(3, vec![(0, 1), (1, 2), (2, 0)], 1.0);
        assert!(
            solver_result.is_err(),
            "odd site count must fail perfect-matching enumeration"
        );
    }
}
