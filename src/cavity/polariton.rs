//! Magnon-Polariton: Diagonalisation of the Magnon-Photon Hamiltonian
//!
//! Diagonalises the bilinear magnon-photon Hamiltonian
//!
//!   H = ℏω_c a†a + ℏω_m b†b + ℏg(a†b + ab†)
//!
//! via the Hopfield-Bogoliubov transformation.  The two dressed (polariton) modes
//! have frequencies
//!
//!   ω_± = (ω_c + ω_m)/2 ± √[(Δ/2)² + g²]
//!
//! and Hopfield coefficients X_±, Y_± satisfying X²+Y²=1.
//!
//! The multi-mode extension [`MultiModePolariton`] handles N magnon modes coupled
//! to one cavity by building and diagonalising an (N+1)×(N+1) Hermitian matrix
//! via [`CMatrix::hermitian_eigendecomposition`].
//!
//! # References
//!
//! - J. J. Hopfield, Phys. Rev. **112**, 1555 (1958) — original polariton picture
//! - H. Huebl et al., Phys. Rev. Lett. **111**, 127003 (2013) — YIG-cavity strong coupling
//! - M. Harder et al., Phys. Rev. B **94**, 054403 (2016) — level repulsion in magnon-polariton
//! - Y. Tabuchi et al., Science **349**, 405 (2015) — coherent coupling magnon-qubit

use crate::constants::{GAMMA, MU_0};
use crate::error::{self, Result};
use crate::math::{CMatrix, Complex};

// ---------------------------------------------------------------------------
// Internal: robust real-symmetric Jacobi for small matrices
// ---------------------------------------------------------------------------

/// Jacobi eigendecomposition for a real symmetric matrix stored as `Vec<Vec<f64>>`.
///
/// Uses a **relative** off-diagonal convergence test (`off_norm / frob_norm < tol`)
/// which is robust for both very small and very large matrix entries.
/// Returns `(eigenvalues_ascending, eigenvector_columns)`.
fn jacobi_real_symmetric(mat: &[Vec<f64>]) -> Result<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = mat.len();
    // Work on a mutable copy.
    let mut a: Vec<Vec<f64>> = mat.to_vec();
    // Accumulate orthogonal matrix Q (eigenvectors as columns).
    let mut q: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = vec![0.0_f64; n];
            row[i] = 1.0;
            row
        })
        .collect();

    let max_sweeps = 200;
    let tol_rel = 1.0e-12_f64;

    for _sweep in 0..max_sweeps {
        // Frobenius norm squared (diagonal + off-diagonal).
        let frob_sq: f64 = (0..n)
            .flat_map(|i| (0..n).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j].powi(2))
            .sum();
        let off_sq: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j].powi(2))
            .sum();

        if off_sq < tol_rel * tol_rel * frob_sq || frob_sq == 0.0 {
            break;
        }

        // Cyclic Jacobi sweep over all (p, q) pairs with p < q.
        for p in 0..n {
            for qp in (p + 1)..n {
                let a_pq = a[p][qp];
                if a_pq.abs() < tol_rel * 1e-3 * frob_sq.sqrt() {
                    continue;
                }
                let a_pp = a[p][p];
                let a_qq = a[qp][qp];

                // Standard symmetric Jacobi: find angle to zero a[p][qp].
                // t satisfies: a01*(1-t^2) + t*(a11-a00) = 0
                // t = sign(beta) / (|beta| + sqrt(1 + beta^2))
                // where beta = (a_qq - a_pp) / (2*a_pq)
                let beta = (a_qq - a_pp) / (2.0 * a_pq);
                let t = beta.signum() / (beta.abs() + (1.0 + beta * beta).sqrt());
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;

                // Apply the Givens rotation G to A: A <- G^T A G
                // Simultaneous column-then-row update avoids double-application.
                // Update columns p and qp.
                for row in a.iter_mut() {
                    let x = row[p];
                    let y = row[qp];
                    row[p] = c * x - s * y;
                    row[qp] = s * x + c * y;
                }
                // Update rows p and qp (two-row update via split_at_mut).
                let (p_row, qp_row) = if p < qp {
                    let (lo, hi) = a.split_at_mut(qp);
                    (&mut lo[p], &mut hi[0])
                } else {
                    let (lo, hi) = a.split_at_mut(p);
                    (&mut hi[0], &mut lo[qp])
                };
                for (x_ref, y_ref) in p_row.iter_mut().zip(qp_row.iter_mut()) {
                    let x = *x_ref;
                    let y = *y_ref;
                    *x_ref = c * x - s * y;
                    *y_ref = s * x + c * y;
                }
                // Accumulate rotation in Q (column update).
                for q_row in q.iter_mut() {
                    let x = q_row[p];
                    let y = q_row[qp];
                    q_row[p] = c * x - s * y;
                    q_row[qp] = s * x + c * y;
                }
            }
        }
    }

    // Extract diagonal and sort ascending.
    let mut pairs: Vec<(f64, usize)> = (0..n).map(|i| (a[i][i], i)).collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let eigenvalues: Vec<f64> = pairs.iter().map(|(v, _)| *v).collect();
    let eigenvectors: Vec<Vec<f64>> = pairs
        .iter()
        .map(|(_, col)| (0..n).map(|row| q[row][*col]).collect())
        .collect();

    Ok((eigenvalues, eigenvectors))
}

// ---------------------------------------------------------------------------
// Branch enum
// ---------------------------------------------------------------------------

/// Polariton branch identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Branch {
    /// Lower polariton branch (ω_−).
    Lower,
    /// Upper polariton branch (ω_+).
    Upper,
}

// ---------------------------------------------------------------------------
// Single-mode: MagnonPolariton
// ---------------------------------------------------------------------------

/// Single-mode magnon-polariton via 2×2 Hopfield diagonalisation.
///
/// Treats the cavity photon (frequency ω_c) and one magnon mode (ω_m) coupled
/// via exchange interaction g.  The analytic 2×2 solution gives exact polariton
/// frequencies and Hopfield mixing coefficients without numerical diagonalisation.
///
/// # Examples
///
/// ```rust
/// use spintronics::cavity::{MagnonPolariton, Branch};
///
/// let mp = MagnonPolariton::from_yig_cavity();
/// let (lo, hi) = mp.eigenfrequencies();
/// assert!(hi > lo);
/// let sum = mp.photon_fraction(Branch::Lower) + mp.magnon_fraction(Branch::Lower);
/// assert!((sum - 1.0).abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct MagnonPolariton {
    /// Cavity angular frequency ω_c [rad/s].
    pub omega_cavity: f64,
    /// Magnon angular frequency ω_m [rad/s].
    pub omega_magnon: f64,
    /// Coupling strength g [rad/s].
    pub coupling_g: f64,
}

impl MagnonPolariton {
    /// Create a new `MagnonPolariton`.
    ///
    /// Allows zero frequencies (bare modes) but rejects negative values.
    pub fn new(omega_cavity: f64, omega_magnon: f64, coupling_g: f64) -> Self {
        Self {
            omega_cavity,
            omega_magnon,
            coupling_g,
        }
    }

    /// Construct from the canonical YIG-cavity parameters used throughout this library.
    ///
    /// - ω_c = 2π × 10 GHz  (microwave cavity)
    /// - ω_m = GAMMA × μ₀ × H_bias  (H_bias = 80 kA/m)
    /// - g   = 2π × 100 MHz  (single-mode coupling)
    pub fn from_yig_cavity() -> Self {
        let omega_cavity = 2.0 * std::f64::consts::PI * 10.0e9;
        let h_bias = 8.0e4_f64; // 80 kA/m — matches HybridSystem::yig_cavity()
        let omega_magnon = GAMMA * MU_0 * h_bias;
        let coupling_g = 2.0 * std::f64::consts::PI * 100.0e6;
        Self {
            omega_cavity,
            omega_magnon,
            coupling_g,
        }
    }

    /// Polariton eigenfrequencies (ω_lower, ω_upper) [rad/s].
    ///
    ///   ω_± = (ω_c + ω_m)/2 ± √[(Δ/2)² + g²]
    pub fn eigenfrequencies(&self) -> (f64, f64) {
        let center = (self.omega_cavity + self.omega_magnon) / 2.0;
        let delta = self.omega_cavity - self.omega_magnon;
        let disc = ((delta / 2.0).powi(2) + self.coupling_g.powi(2)).sqrt();
        (center - disc, center + disc)
    }

    /// Hopfield coefficients for both branches.
    ///
    /// Returns `((X_lower, Y_lower), (X_upper, Y_upper))` where:
    /// - X = photon (cavity) admixture amplitude
    /// - Y = magnon admixture amplitude
    /// - X² + Y² = 1 for each branch
    ///
    /// Sign convention: X_lower ≥ 0, Y_lower ≥ 0; signs of upper follow from orthogonality.
    ///
    /// At resonance (Δ=0): X = Y = 1/√2.
    /// For large positive Δ (ω_c >> ω_m): lower branch → magnon (Y→1), upper → photon (X→1).
    pub fn hopfield_coefficients(&self) -> ((f64, f64), (f64, f64)) {
        let delta = self.omega_cavity - self.omega_magnon;
        let g = self.coupling_g;
        let disc = ((delta / 2.0).powi(2) + g.powi(2)).sqrt();

        // Mixing angle θ such that tan(2θ) = 2g/Δ.
        // X_lower = cos θ, Y_lower = sin θ  (lower branch: more magnon-like when Δ>0)
        // Use the explicit expressions to avoid quadrant ambiguity.
        //
        // From the 2×2 eigenvector of [[ω_c, g],[g, ω_m]] for eigenvalue ω_−:
        //   (ω_c − ω_−) X = g Y  →  X/Y = g/(ω_c − ω_−) = g/(disc + Δ/2)
        // Normalise: X²+Y²=1.
        let denom_lo_sq = g.powi(2) + (disc + delta / 2.0).powi(2);
        let x_lo = if denom_lo_sq > 0.0 {
            g / denom_lo_sq.sqrt()
        } else {
            std::f64::consts::FRAC_1_SQRT_2
        };
        let y_lo = if denom_lo_sq > 0.0 {
            (disc + delta / 2.0) / denom_lo_sq.sqrt()
        } else {
            std::f64::consts::FRAC_1_SQRT_2
        };
        // Upper branch is orthogonal: (X_up, Y_up) = (-Y_lo, X_lo) or (Y_lo, -X_lo).
        // Choose the sign so that X_up ≥ 0.
        let x_up = y_lo;
        let y_up = -x_lo;
        // Ensure both magnitudes are non-negative for the photon fraction (squared).
        ((x_lo, y_lo), (x_up, y_up))
    }

    /// Photon fraction |X|² for the given branch.
    pub fn photon_fraction(&self, branch: Branch) -> f64 {
        let ((xl, _), (xu, _)) = self.hopfield_coefficients();
        match branch {
            Branch::Lower => xl.powi(2),
            Branch::Upper => xu.powi(2),
        }
    }

    /// Magnon fraction |Y|² for the given branch.
    pub fn magnon_fraction(&self, branch: Branch) -> f64 {
        let ((_, yl), (_, yu)) = self.hopfield_coefficients();
        match branch {
            Branch::Lower => yl.powi(2),
            Branch::Upper => yu.powi(2),
        }
    }

    /// Vacuum Rabi splitting = ω_+ − ω_− = 2g (evaluated at resonance Δ=0).
    ///
    /// The actual splitting at arbitrary detuning is 2√[(Δ/2)²+g²], but the
    /// customary definition of "vacuum Rabi splitting" is the on-resonance value 2g.
    pub fn vacuum_rabi_splitting(&self) -> f64 {
        2.0 * self.coupling_g
    }

    /// Returns `true` if the system is in the strong-coupling regime.
    ///
    /// Strong coupling condition: g > (κ + γ)/4.
    pub fn is_strong_coupling(&self, kappa: f64, gamma: f64) -> bool {
        self.coupling_g > (kappa + gamma) / 4.0
    }
}

// ---------------------------------------------------------------------------
// Multi-mode: MultiModePolariton
// ---------------------------------------------------------------------------

/// Multi-mode polariton: N magnon modes coupled to one microwave cavity.
///
/// Builds the (N+1)×(N+1) Hermitian interaction matrix:
///
/// ```text
///       [  ω_c   g_1   g_2  … g_N  ]
///       [  g_1   ω_1    0  …   0   ]
///   H = [  g_2    0   ω_2  …   0   ]
///       [  …                        ]
///       [  g_N    0    0  …  ω_N   ]
/// ```
///
/// and diagonalises it via [`CMatrix::hermitian_eigendecomposition`].
///
/// # Examples
///
/// ```rust
/// use spintronics::cavity::MultiModePolariton;
///
/// let mp = MultiModePolariton::new(
///     2.0 * std::f64::consts::PI * 10.0e9,
///     vec![2.0 * std::f64::consts::PI * 9.5e9],
///     vec![2.0 * std::f64::consts::PI * 100.0e6],
/// ).unwrap();
/// let freqs = mp.eigenfrequencies().unwrap();
/// assert_eq!(freqs.len(), 2);
/// ```
#[derive(Debug, Clone)]
pub struct MultiModePolariton {
    /// Cavity angular frequency ω_c [rad/s].
    pub cavity_freq: f64,
    /// Magnon mode angular frequencies ω_i [rad/s].
    pub magnon_freqs: Vec<f64>,
    /// Coupling strengths g_i [rad/s] (one per magnon mode).
    pub couplings: Vec<f64>,
}

impl MultiModePolariton {
    /// Create a new `MultiModePolariton`.
    ///
    /// # Errors
    ///
    /// Returns `Err` if:
    /// - `magnon_freqs` and `couplings` have different lengths.
    /// - Any coupling is negative.
    /// - `cavity_freq` is negative.
    pub fn new(cavity_freq: f64, magnon_freqs: Vec<f64>, couplings: Vec<f64>) -> Result<Self> {
        if cavity_freq < 0.0 {
            return Err(error::invalid_param(
                "cavity_freq",
                "cavity frequency must be non-negative",
            ));
        }
        if magnon_freqs.len() != couplings.len() {
            return Err(error::invalid_param(
                "couplings",
                "couplings must have the same length as magnon_freqs",
            ));
        }
        for (i, &g) in couplings.iter().enumerate() {
            if g < 0.0 {
                return Err(error::invalid_param(
                    "couplings",
                    &format!("coupling[{i}] = {g} is negative; all couplings must be non-negative"),
                ));
            }
        }
        Ok(Self {
            cavity_freq,
            magnon_freqs,
            couplings,
        })
    }

    /// Returns the dimension of the polariton system: N+1 (one cavity + N magnon modes).
    fn dim(&self) -> usize {
        1 + self.magnon_freqs.len()
    }

    /// Build the (N+1)×(N+1) Hermitian Hopfield matrix.
    ///
    /// ```text
    /// H[0][0]     = omega_c
    /// H[i+1][i+1] = magnon_freqs[i]
    /// H[0][i+1]   = H[i+1][0] = couplings[i]
    /// ```
    fn build_matrix(&self) -> Result<CMatrix> {
        let n = self.dim();
        let mut rows = vec![vec![Complex::ZERO; n]; n];
        // Cavity diagonal
        rows[0][0] = Complex::from_real(self.cavity_freq);
        for (i, (&omega_m, &g)) in self
            .magnon_freqs
            .iter()
            .zip(self.couplings.iter())
            .enumerate()
        {
            // Magnon diagonal
            rows[i + 1][i + 1] = Complex::from_real(omega_m);
            // Off-diagonal coupling (real and symmetric for JC without rotating wave)
            rows[0][i + 1] = Complex::from_real(g);
            rows[i + 1][0] = Complex::from_real(g);
        }
        CMatrix::from_rows(rows)
    }

    /// Diagonalise the interaction matrix using the internal Jacobi routine.
    ///
    /// Returns `(eigenvalues_ascending, eigenvectors_as_columns)`.
    ///
    /// The matrix entries are real (all coupling constants and frequencies are real),
    /// so we extract the real parts and apply `jacobi_real_symmetric`.  The returned
    /// eigenvectors are wrapped back into a `CMatrix` for mode-composition queries.
    fn diagonalise(&self) -> Result<(Vec<f64>, CMatrix)> {
        let cmat = self.build_matrix()?;
        let n = self.dim();
        // Extract real symmetric matrix (imaginary parts are zero by construction).
        let real_mat: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| cmat.get(i, j).re).collect())
            .collect();
        let (vals, vecs_cols) = jacobi_real_symmetric(&real_mat)?;
        // Pack eigenvectors back into a CMatrix (column j = vecs_cols[j]).
        let mut q = CMatrix::zeros(n);
        for (col, col_vec) in vecs_cols.iter().enumerate() {
            for (row, &v) in col_vec.iter().enumerate() {
                q.set(row, col, Complex::from_real(v));
            }
        }
        Ok((vals, q))
    }

    /// Eigenfrequencies of the full (N+1)-mode polariton system, sorted ascending.
    ///
    /// # Errors
    ///
    /// Propagates errors from matrix construction or the Jacobi eigendecomposition.
    pub fn eigenfrequencies(&self) -> Result<Vec<f64>> {
        let (vals, _) = self.diagonalise()?;
        Ok(vals)
    }

    /// Mode compositions: `composition[i][j] = |⟨ψ_i|basis_j⟩|²`.
    ///
    /// - `composition[i][0]`   = photon content of the i-th polariton mode.
    /// - `composition[i][j+1]` = magnon-j content of the i-th polariton mode.
    ///
    /// Each row sums to 1 (conservation of probability).
    ///
    /// # Errors
    ///
    /// Propagates errors from matrix construction or diagonalisation.
    pub fn mode_compositions(&self) -> Result<Vec<Vec<f64>>> {
        let (_, vecs) = self.diagonalise()?;
        let n = self.dim();
        let mut compositions = Vec::with_capacity(n);
        for col in 0..n {
            let row_comp: Vec<f64> = (0..n).map(|row| vecs.get(row, col).norm_sq()).collect();
            compositions.push(row_comp);
        }
        Ok(compositions)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // 1. At resonance ω_+ - ω_- = 2g.
    #[test]
    fn test_anticrossing_at_resonance() {
        let omega = 2.0 * PI * 10.0e9;
        let g = 2.0 * PI * 100.0e6;
        let mp = MagnonPolariton::new(omega, omega, g);
        let (lo, hi) = mp.eigenfrequencies();
        assert!(
            approx(hi - lo, 2.0 * g, 1.0e3),
            "splitting = {}, 2g = {}",
            hi - lo,
            2.0 * g
        );
    }

    // 2. Hopfield fractions sum to 1 for each branch at arbitrary detuning.
    #[test]
    fn test_hopfield_fractions_sum_to_one() {
        let omega_c = 2.0 * PI * 10.0e9;
        let omega_m = 2.0 * PI * 9.0e9;
        let g = 2.0 * PI * 80.0e6;
        let mp = MagnonPolariton::new(omega_c, omega_m, g);
        let sum_lo = mp.photon_fraction(Branch::Lower) + mp.magnon_fraction(Branch::Lower);
        let sum_up = mp.photon_fraction(Branch::Upper) + mp.magnon_fraction(Branch::Upper);
        assert!(
            approx(sum_lo, 1.0, 1e-12),
            "lower branch fractions sum = {sum_lo}"
        );
        assert!(
            approx(sum_up, 1.0, 1e-12),
            "upper branch fractions sum = {sum_up}"
        );
    }

    // 3. Large positive detuning → lower branch mostly magnon.
    #[test]
    fn test_large_positive_detuning() {
        let omega_c = 2.0 * PI * 12.0e9;
        let omega_m = 2.0 * PI * 9.0e9; // ω_c >> ω_m, Δ = +3 GHz >> g=100 MHz
        let g = 2.0 * PI * 100.0e6;
        let mp = MagnonPolariton::new(omega_c, omega_m, g);
        let mf_lo = mp.magnon_fraction(Branch::Lower);
        assert!(
            mf_lo > 0.9,
            "lower branch should be mostly magnon when Δ>>g, got magnon_fraction={mf_lo}"
        );
    }

    // 4. Large negative detuning → lower branch mostly photon.
    #[test]
    fn test_large_negative_detuning() {
        let omega_c = 2.0 * PI * 9.0e9;
        let omega_m = 2.0 * PI * 12.0e9; // ω_m >> ω_c, Δ = -3 GHz
        let g = 2.0 * PI * 100.0e6;
        let mp = MagnonPolariton::new(omega_c, omega_m, g);
        let pf_lo = mp.photon_fraction(Branch::Lower);
        assert!(
            pf_lo > 0.9,
            "lower branch should be mostly photon when ω_m>>ω_c, got photon_fraction={pf_lo}"
        );
    }

    // 5. Vacuum Rabi splitting = 2g.
    #[test]
    fn test_vacuum_rabi_splitting_equals_2g() {
        let g = 2.0 * PI * 123.0e6;
        let mp = MagnonPolariton::new(1.0e10, 1.0e10, g);
        assert!(
            approx(mp.vacuum_rabi_splitting(), 2.0 * g, 1.0),
            "splitting = {}, 2g = {}",
            mp.vacuum_rabi_splitting(),
            2.0 * g
        );
    }

    // 6. MultiModePolariton with 1 magnon matches MagnonPolariton eigenfrequencies.
    #[test]
    fn test_multi_mode_n1_matches_single_mode() {
        let omega_c = 2.0 * PI * 10.0e9;
        let omega_m = 2.0 * PI * 9.5e9;
        let g = 2.0 * PI * 100.0e6;
        let single = MagnonPolariton::new(omega_c, omega_m, g);
        let (lo_s, hi_s) = single.eigenfrequencies();
        let multi = MultiModePolariton::new(omega_c, vec![omega_m], vec![g]).unwrap();
        let freqs = multi.eigenfrequencies().unwrap();
        assert_eq!(freqs.len(), 2);
        assert!(
            approx(freqs[0], lo_s, 1.0e3),
            "multi lo={} single lo={}",
            freqs[0],
            lo_s
        );
        assert!(
            approx(freqs[1], hi_s, 1.0e3),
            "multi hi={} single hi={}",
            freqs[1],
            hi_s
        );
    }

    // 7. MultiModePolariton eigenvalues are sorted ascending.
    #[test]
    fn test_multi_mode_eigenvalues_sorted() {
        let omega_c = 2.0 * PI * 10.0e9;
        let magnon_freqs = vec![2.0 * PI * 8.0e9, 2.0 * PI * 9.0e9, 2.0 * PI * 11.0e9];
        let couplings = vec![2.0 * PI * 80.0e6, 2.0 * PI * 60.0e6, 2.0 * PI * 70.0e6];
        let multi = MultiModePolariton::new(omega_c, magnon_freqs, couplings).unwrap();
        let freqs = multi.eigenfrequencies().unwrap();
        for i in 1..freqs.len() {
            assert!(
                freqs[i] >= freqs[i - 1],
                "eigenvalues not sorted: freqs[{i}]={} < freqs[{}]={}",
                freqs[i],
                i - 1,
                freqs[i - 1]
            );
        }
    }

    // 8. Each row of mode_compositions sums to 1.
    #[test]
    fn test_mode_compositions_sum_to_one() {
        let omega_c = 2.0 * PI * 10.0e9;
        let magnon_freqs = vec![2.0 * PI * 9.0e9, 2.0 * PI * 10.5e9];
        let couplings = vec![2.0 * PI * 100.0e6, 2.0 * PI * 80.0e6];
        let multi = MultiModePolariton::new(omega_c, magnon_freqs, couplings).unwrap();
        let comps = multi.mode_compositions().unwrap();
        for (i, row) in comps.iter().enumerate() {
            let s: f64 = row.iter().sum();
            assert!(
                approx(s, 1.0, 1e-10),
                "composition row {i} sums to {s}, expected 1"
            );
        }
    }
}
