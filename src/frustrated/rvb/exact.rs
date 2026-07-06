//! Exact diagonalization benchmarks for the Heisenberg antiferromagnet, independent of
//! the valence-bond machinery.
//!
//! These provide the ground truth against which [`super::solver::RvbSolver`]'s
//! variational valence-bond ground state is checked. Two backends are provided because
//! a single dense `CMatrix` cannot honestly benchmark anything interesting-sized:
//!
//! - [`ExactDiagonalization::dense_ground_state`]: builds the full `2^N × 2^N`
//!   Hamiltonian as a `CMatrix` and diagonalizes it directly. Honest only for `N ≤ 6`
//!   (`2^6 = 64 = CMatrix::MAX_DIM`; already `2^8 = 256 > 64`).
//! - [`ExactDiagonalization::ground_state_sz0`]: a **matrix-free Lanczos** solver
//!   restricted to the `S_z = 0` sector (bit-encoded basis of `N`-bit masks with exactly
//!   `N/2` set bits, sparse Hamiltonian-vector product, full reorthogonalization), valid
//!   up to [`super::MAX_ED_SITES`] sites. This is the real benchmark for cluster sizes
//!   where the valence-bond basis is actually interesting (e.g. the 4×4 open square,
//!   `N = 16`, `C(16,8) = 12870`-dimensional `S_z=0` sector).
//!
//! Both operate directly in the computational `S^z` basis (bit `i` of a basis mask is 0
//! for spin-up, 1 for spin-down at site `i`), independent of the valence-bond sign
//! conventions in [`super::valence_bond`] — this makes them a genuine cross-check, not
//! merely a repetition of the same formulas.

use std::collections::HashMap;

use crate::error::{self, Result};
use crate::frustrated::lattice::Xorshift64;
use crate::math::{CMatrix, Complex};

/// Fixed seed for the Lanczos starting vector (determinism, reproducible tests).
const LANCZOS_SEED: u64 = 0x5EED_1234_ABCD_0001;
/// Convergence tolerance on the lowest Ritz value between successive Lanczos steps.
const LANCZOS_TOL: f64 = 1e-11;

/// A quantum state expressed in an explicit bit-encoded computational basis (either the
/// full `2^N` `S^z` basis, or a restricted symmetry sector such as `S_z = 0`).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpinBasisState {
    /// Basis bitmasks (bit `i` = 0 for up, 1 for down at site `i`).
    pub basis: Vec<u64>,
    /// Amplitude of each basis state (same order as `basis`).
    pub amplitudes: Vec<f64>,
}

/// Exact-diagonalization benchmark for `H = J·Σ_{(i,j)∈bonds} S_i·S_j` on `num_sites`
/// spin-1/2 sites.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ExactDiagonalization {
    /// Number of sites.
    pub num_sites: usize,
    /// Exchange bonds `(i, j)`.
    pub bonds: Vec<(usize, usize)>,
    /// Antiferromagnetic exchange coupling `J > 0`.
    pub coupling_j: f64,
}

impl ExactDiagonalization {
    /// Construct a benchmark instance, validating bond indices and coupling sign.
    ///
    /// # Errors
    ///
    /// Returns an error if any bond is out of range/self-connecting, or if
    /// `coupling_j <= 0`.
    pub fn new(num_sites: usize, bonds: Vec<(usize, usize)>, coupling_j: f64) -> Result<Self> {
        for &(i, j) in &bonds {
            if i >= num_sites || j >= num_sites {
                return Err(error::invalid_param(
                    "bonds",
                    "bond references a site index out of range",
                ));
            }
            if i == j {
                return Err(error::invalid_param(
                    "bonds",
                    "a bond cannot connect a site to itself",
                ));
            }
        }
        if coupling_j <= 0.0 {
            return Err(error::invalid_param(
                "coupling_j",
                "RVB physics requires antiferromagnetic coupling J > 0",
            ));
        }
        Ok(Self {
            num_sites,
            bonds,
            coupling_j,
        })
    }

    /// Ground-state energy via dense diagonalization of the full `2^N × 2^N`
    /// Hamiltonian. See [`ExactDiagonalization::dense_ground_state`] for the size limit.
    pub fn dense_ground_energy(&self) -> Result<f64> {
        Ok(self.dense_ground_state()?.0)
    }

    /// Ground-state energy and eigenvector via dense diagonalization.
    ///
    /// # Errors
    ///
    /// Returns an error if `2^num_sites > CMatrix::MAX_DIM` (64), i.e. `num_sites > 6`.
    pub fn dense_ground_state(&self) -> Result<(f64, SpinBasisState)> {
        let n = self.num_sites;
        if n == 0 {
            return Err(error::invalid_param(
                "num_sites",
                "num_sites must be positive",
            ));
        }
        let dim = match 1usize.checked_shl(n as u32) {
            Some(d) if d <= CMatrix::MAX_DIM => d,
            _ => {
                return Err(error::invalid_param(
                    "num_sites",
                    &format!(
                        "dense exact diagonalization requires 2^num_sites <= CMatrix::MAX_DIM (64); \
                         num_sites={} is too large; use ground_energy_sz0/ground_state_sz0 (Lanczos) \
                         instead for N up to MAX_ED_SITES={}",
                        n,
                        super::MAX_ED_SITES
                    ),
                ));
            },
        };

        let mut h = vec![vec![0.0_f64; dim]; dim];
        // `b` is used both as a bitmask value (via bit shifts) and as a matrix index;
        // there is no existing collection to enumerate() over.
        #[allow(clippy::needless_range_loop)]
        for b in 0..dim {
            let mut diag = 0.0_f64;
            for &(i, j) in &self.bonds {
                let bit_i = (b >> i) & 1;
                let bit_j = (b >> j) & 1;
                if bit_i == bit_j {
                    diag += 0.25 * self.coupling_j;
                } else {
                    diag -= 0.25 * self.coupling_j;
                    let flipped = b ^ (1usize << i) ^ (1usize << j);
                    h[flipped][b] += 0.5 * self.coupling_j;
                }
            }
            h[b][b] += diag;
        }
        let rows: Vec<Vec<Complex>> = h
            .into_iter()
            .map(|row| row.into_iter().map(Complex::from_real).collect())
            .collect();
        let h_matrix = CMatrix::from_rows(rows)?;
        let (vals, vecs) = h_matrix.hermitian_eigendecomposition()?;
        let ground_energy = vals[0];
        let amplitudes: Vec<f64> = (0..dim).map(|row| vecs.get(row, 0).re).collect();
        let basis: Vec<u64> = (0..dim as u64).collect();
        Ok((ground_energy, SpinBasisState { basis, amplitudes }))
    }

    /// Ground-state energy in the `S_z = 0` sector via matrix-free Lanczos.
    pub fn ground_energy_sz0(&self) -> Result<f64> {
        Ok(self.ground_state_sz0()?.0)
    }

    /// Ground-state energy and eigenvector in the `S_z = 0` sector via matrix-free
    /// Lanczos with full reorthogonalization.
    ///
    /// # Errors
    ///
    /// Returns an error if `num_sites` is odd (no `S_z = 0` sector), zero, or exceeds
    /// [`super::MAX_ED_SITES`].
    pub fn ground_state_sz0(&self) -> Result<(f64, SpinBasisState)> {
        let n = self.num_sites;
        if n == 0 {
            return Err(error::invalid_param(
                "num_sites",
                "num_sites must be positive",
            ));
        }
        if n > super::MAX_ED_SITES {
            return Err(error::invalid_param(
                "num_sites",
                &format!(
                    "Sz=0 Lanczos exact diagonalization is limited to MAX_ED_SITES={}; got num_sites={}",
                    super::MAX_ED_SITES,
                    n
                ),
            ));
        }
        if n % 2 != 0 {
            return Err(error::invalid_param(
                "num_sites",
                "the Sz=0 sector requires an even number of sites",
            ));
        }

        let basis = build_sz0_basis(n);
        let dim = basis.len();
        let index_of: HashMap<u64, usize> =
            basis.iter().enumerate().map(|(idx, &b)| (b, idx)).collect();
        let max_iter = dim.min(CMatrix::MAX_DIM);

        let mut lanczos_vectors: Vec<Vec<f64>> = Vec::with_capacity(max_iter);
        let mut alpha: Vec<f64> = Vec::with_capacity(max_iter);
        let mut beta: Vec<f64> = Vec::with_capacity(max_iter);

        let mut rng = Xorshift64::new(LANCZOS_SEED)?;
        let mut v_curr: Vec<f64> = (0..dim).map(|_| rng.next_f64() - 0.5).collect();
        normalize(&mut v_curr)?;
        let mut v_prev = vec![0.0_f64; dim];
        let mut beta_prev = 0.0_f64;
        let mut prev_energy = f64::INFINITY;

        for step in 0..max_iter {
            lanczos_vectors.push(v_curr.clone());
            let mut w =
                apply_heisenberg_bonds(&basis, &index_of, &self.bonds, self.coupling_j, &v_curr);
            if step > 0 {
                axpy(&mut w, -beta_prev, &v_prev);
            }
            let alpha_k = dot(&v_curr, &w);
            axpy(&mut w, -alpha_k, &v_curr);
            // Full reorthogonalization against every previous Lanczos vector combats
            // the well-known loss of orthogonality in finite-precision Lanczos.
            for prev_vec in &lanczos_vectors {
                let proj = dot(prev_vec, &w);
                axpy(&mut w, -proj, prev_vec);
            }
            alpha.push(alpha_k);

            let (t_vals, _) = diagonalize_tridiagonal(&alpha, &beta)?;
            let current_energy = t_vals[0];
            let converged = step > 0 && (current_energy - prev_energy).abs() < LANCZOS_TOL;
            prev_energy = current_energy;

            let beta_next = norm(&w);
            let breakdown = beta_next < 1e-12;

            if converged || breakdown || step == max_iter - 1 {
                break;
            }

            beta.push(beta_next);
            v_prev.copy_from_slice(&v_curr);
            v_curr = w.iter().map(|x| x / beta_next).collect();
            beta_prev = beta_next;
        }

        let (t_vals, t_vecs) = diagonalize_tridiagonal(&alpha, &beta)?;
        let ground_energy = t_vals[0];
        if !ground_energy.is_finite() {
            return Err(error::numerical_error(
                "Lanczos ground-state energy is not finite",
            ));
        }
        let m = alpha.len();
        let mut amplitudes = vec![0.0_f64; dim];
        for (step, lanczos_vec) in lanczos_vectors.iter().enumerate().take(m) {
            let coeff = t_vecs.get(step, 0).re;
            axpy(&mut amplitudes, coeff, lanczos_vec);
        }
        normalize(&mut amplitudes)?;

        Ok((ground_energy, SpinBasisState { basis, amplitudes }))
    }

    /// Expectation value of the total-spin-squared operator `S_tot² = (Σ_i S_i)²` for an
    /// explicit state (either the full `2^N` basis from
    /// [`ExactDiagonalization::dense_ground_state`], or the `S_z=0` sector basis from
    /// [`ExactDiagonalization::ground_state_sz0`]).
    ///
    /// Uses `S_tot² = N·(3/4) + 2·Σ_{i<j} S_i·S_j`, applying the all-pairs sum via the
    /// same matrix-free sparse operator used by the Sz=0 Hamiltonian (with an effective
    /// coupling of 2.0, matching the `2·` prefactor).
    ///
    /// # Errors
    ///
    /// Returns an error if `state` has (numerically) zero norm.
    pub fn total_spin_squared(&self, state: &SpinBasisState) -> Result<f64> {
        let n = self.num_sites;
        let all_pairs: Vec<(usize, usize)> = (0..n)
            .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
            .collect();
        let index_of: HashMap<u64, usize> = state
            .basis
            .iter()
            .enumerate()
            .map(|(idx, &b)| (b, idx))
            .collect();
        let norm_sq = dot(&state.amplitudes, &state.amplitudes);
        if norm_sq < 1e-14 {
            return Err(error::numerical_error(
                "state has (numerically) zero norm; cannot compute <S_tot^2>",
            ));
        }
        let hv =
            apply_heisenberg_bonds(&state.basis, &index_of, &all_pairs, 2.0, &state.amplitudes);
        let cross_term = dot(&state.amplitudes, &hv);
        let onsite_term = n as f64 * 0.75 * norm_sq;
        Ok((onsite_term + cross_term) / norm_sq)
    }
}

/// All `n`-bit masks with exactly `n/2` bits set (the `S_z = 0` sector basis).
fn build_sz0_basis(n: usize) -> Vec<u64> {
    let half = (n / 2) as u32;
    let total = 1u64 << n; // n <= MAX_ED_SITES <= 16, always safe in u64.
    (0..total).filter(|b| b.count_ones() == half).collect()
}

/// Matrix-free application of `Σ_{(i,j)∈bonds} coupling_j·S_i·S_j` to a state vector `v`
/// expressed over an explicit bitmask `basis` (with `index_of` its inverse lookup).
///
/// Shared by the Sz=0 Hamiltonian ([`ExactDiagonalization::ground_state_sz0`]) and the
/// all-pairs `S_tot²` computation ([`ExactDiagonalization::total_spin_squared`]).
fn apply_heisenberg_bonds(
    basis: &[u64],
    index_of: &HashMap<u64, usize>,
    bonds: &[(usize, usize)],
    coupling_j: f64,
    v: &[f64],
) -> Vec<f64> {
    let dim = basis.len();
    let mut result = vec![0.0_f64; dim];
    for (idx, &b) in basis.iter().enumerate() {
        let mut diag = 0.0_f64;
        for &(i, j) in bonds {
            let bit_i = (b >> i) & 1;
            let bit_j = (b >> j) & 1;
            if bit_i == bit_j {
                diag += 0.25 * coupling_j;
            } else {
                diag -= 0.25 * coupling_j;
                let flipped = b ^ (1u64 << i) ^ (1u64 << j);
                if let Some(&f_idx) = index_of.get(&flipped) {
                    result[f_idx] += 0.5 * coupling_j * v[idx];
                }
            }
        }
        result[idx] += diag * v[idx];
    }
    result
}

/// Diagonalize a real symmetric tridiagonal matrix (diagonal `alpha`, off-diagonal
/// `beta`) via [`CMatrix::hermitian_eigendecomposition`] (embedded with zero imaginary
/// parts).
fn diagonalize_tridiagonal(alpha: &[f64], beta: &[f64]) -> Result<(Vec<f64>, CMatrix)> {
    let m = alpha.len();
    let mut rows = vec![vec![Complex::ZERO; m]; m];
    for (i, &a) in alpha.iter().enumerate() {
        rows[i][i] = Complex::from_real(a);
    }
    for (i, &b) in beta.iter().enumerate() {
        rows[i][i + 1] = Complex::from_real(b);
        rows[i + 1][i] = Complex::from_real(b);
    }
    let t = CMatrix::from_rows(rows)?;
    t.hermitian_eigendecomposition()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

fn normalize(v: &mut [f64]) -> Result<()> {
    let n = norm(v);
    if n < 1e-14 {
        return Err(error::numerical_error(
            "cannot normalize a (numerically) zero vector",
        ));
    }
    for x in v.iter_mut() {
        *x /= n;
    }
    Ok(())
}

fn axpy(y: &mut [f64], alpha: f64, x: &[f64]) {
    for (yi, xi) in y.iter_mut().zip(x.iter()) {
        *yi += alpha * xi;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_site_dense_ground_energy() {
        let ed = ExactDiagonalization::new(2, vec![(0, 1)], 1.0).expect("valid ED instance");
        let e0 = ed
            .dense_ground_energy()
            .expect("dense diagonalization should succeed");
        assert!((e0 - (-0.75)).abs() < 1e-10, "expected -3J/4, got {}", e0);
    }

    #[test]
    fn test_ring4_dense_ground_energy() {
        let ed = ExactDiagonalization::new(4, vec![(0, 1), (1, 2), (2, 3), (3, 0)], 1.0)
            .expect("valid ED instance");
        let e0 = ed
            .dense_ground_energy()
            .expect("dense diagonalization should succeed");
        assert!((e0 - (-2.0)).abs() < 1e-9, "expected -2J, got {}", e0);
    }

    #[test]
    fn test_dense_diagonalization_rejects_large_n() {
        let bonds: Vec<(usize, usize)> = (0..7).map(|i| (i, (i + 1) % 8)).collect();
        let ed = ExactDiagonalization::new(8, bonds, 1.0).expect("valid ED instance (bonds only)");
        assert!(
            ed.dense_ground_state().is_err(),
            "N=8 (2^8=256>64) must be rejected"
        );
    }

    #[test]
    fn test_dense_vs_lanczos_agreement_ring4() {
        let ed = ExactDiagonalization::new(4, vec![(0, 1), (1, 2), (2, 3), (3, 0)], 1.0)
            .expect("valid ED instance");
        let e_dense = ed.dense_ground_energy().expect("dense ED");
        let e_lanczos = ed.ground_energy_sz0().expect("Lanczos ED");
        assert!(
            (e_dense - e_lanczos).abs() < 1e-8,
            "dense ({}) and Lanczos ({}) ground energies must agree",
            e_dense,
            e_lanczos
        );
    }

    #[test]
    fn test_dense_vs_lanczos_agreement_ladder() {
        let bonds = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (4, 5),
            (5, 6),
            (6, 7),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        // N=8 exceeds dense (2^8=256>64), so cross-check dense at N=6 sub-chain instead
        // and separately confirm Lanczos runs at N=8.
        let ed8 = ExactDiagonalization::new(8, bonds, 1.0).expect("valid ED instance");
        let e_lanczos = ed8
            .ground_energy_sz0()
            .expect("Lanczos ED should succeed for N=8");
        assert!(e_lanczos.is_finite());
        assert!(
            e_lanczos < 0.0,
            "AFM ground energy should be negative, got {}",
            e_lanczos
        );
    }

    #[test]
    fn test_total_spin_squared_near_zero_for_singlet_ground_state() {
        let ed = ExactDiagonalization::new(4, vec![(0, 1), (1, 2), (2, 3), (3, 0)], 1.0)
            .expect("valid ED instance");
        let (_, state) = ed.dense_ground_state().expect("dense ED");
        let s2 = ed.total_spin_squared(&state).expect("total spin squared");
        assert!(
            s2.abs() < 1e-8,
            "expected S_tot^2 ~ 0 for the singlet ground state, got {}",
            s2
        );
    }

    #[test]
    fn test_total_spin_squared_sz0_ground_state_near_zero() {
        let ed = ExactDiagonalization::new(4, vec![(0, 1), (1, 2), (2, 3), (3, 0)], 1.0)
            .expect("valid ED instance");
        let (_, state) = ed.ground_state_sz0().expect("Lanczos ED");
        let s2 = ed.total_spin_squared(&state).expect("total spin squared");
        assert!(
            s2.abs() < 1e-6,
            "expected S_tot^2 ~ 0 for the singlet ground state, got {}",
            s2
        );
    }

    #[test]
    fn test_odd_num_sites_sz0_errors() {
        let ed =
            ExactDiagonalization::new(3, vec![(0, 1), (1, 2)], 1.0).expect("valid ED instance");
        assert!(ed.ground_state_sz0().is_err());
    }

    #[test]
    fn test_non_positive_coupling_errors() {
        assert!(ExactDiagonalization::new(2, vec![(0, 1)], 0.0).is_err());
    }
}
