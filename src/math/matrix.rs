//! Dense complex matrix (`CMatrix`) for small-dimensional quantum mechanics.
//!
//! Provides an N×N complex matrix type with:
//! - Gauss-Jordan inversion with partial pivoting
//! - Jacobi rotation diagonalization for Hermitian matrices
//! - Standard arithmetic operations: multiply, add, subtract, scale, trace, conjugate-transpose
//!
//! All operations are O(N³) and intended for matrices with N ≤ 64. This is sufficient
//! for magnon band models (2×2 honeycomb, 3×3 kagome), NEGF Hamiltonians (≤ 32 sites),
//! and strip Hamiltonians for edge-mode calculations (≤ 60 sites).
//!
//! # Example
//!
//! ```rust
//! use spintronics::math::{CMatrix, Complex};
//!
//! let h = CMatrix::eye(2);
//! let (vals, vecs) = h.hermitian_eigendecomposition().unwrap();
//! assert!((vals[0] - 1.0).abs() < 1e-12);
//! assert!((vals[1] - 1.0).abs() < 1e-12);
//! ```

use crate::error::{self, Result};
use crate::math::Complex;

/// Dense N×N complex matrix stored in row-major order.
#[derive(Debug, Clone)]
pub struct CMatrix {
    /// Row-major flat storage: data[i*n + j] = M\[i\]\[j\]
    data: Vec<Complex>,
    n: usize,
}

impl CMatrix {
    /// Maximum allowed dimension (runtime enforcement).
    pub const MAX_DIM: usize = 64;

    /// Create an N×N zero matrix.
    pub fn zeros(n: usize) -> Self {
        Self {
            data: vec![Complex::ZERO; n * n],
            n,
        }
    }

    /// Create an N×N identity matrix.
    pub fn eye(n: usize) -> Self {
        let mut m = Self::zeros(n);
        for i in 0..n {
            m.set(i, i, Complex::ONE);
        }
        m
    }

    /// Build from a `Vec<Vec<Complex>>` (rows outer, columns inner).
    ///
    /// Returns an error if the rows are inconsistent or empty.
    pub fn from_rows(rows: Vec<Vec<Complex>>) -> Result<Self> {
        let n = rows.len();
        if n == 0 {
            return Err(error::invalid_param("rows", "matrix must be non-empty"));
        }
        if n > Self::MAX_DIM {
            return Err(error::invalid_param(
                "n",
                "matrix dimension exceeds CMatrix::MAX_DIM (64)",
            ));
        }
        for row in rows.iter() {
            if row.len() != n {
                return Err(error::invalid_param(
                    "rows",
                    "all rows must have length equal to the number of rows (square matrix)",
                ));
            }
        }
        let mut data = Vec::with_capacity(n * n);
        for row in &rows {
            data.extend_from_slice(row);
        }
        Ok(Self { data, n })
    }

    /// Dimension N of this N×N matrix.
    #[inline]
    pub fn n(&self) -> usize {
        self.n
    }

    /// Get element M\[i\]\[j\].
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Complex {
        self.data[i * self.n + j]
    }

    /// Set element M\[i\]\[j\].
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: Complex) {
        self.data[i * self.n + j] = v;
    }

    /// Add `v` to element M\[i\]\[j\] in place.
    #[inline]
    fn add_to(&mut self, i: usize, j: usize, v: Complex) {
        let cur = self.get(i, j);
        self.set(i, j, cur.add(&v));
    }

    /// Trace = Σ M\[i\]\[i\].
    pub fn trace(&self) -> Complex {
        let mut t = Complex::ZERO;
        for i in 0..self.n {
            t = t.add(&self.get(i, i));
        }
        t
    }

    /// Conjugate-transpose (Hermitian adjoint) M†.
    pub fn conj_transpose(&self) -> Self {
        let mut out = Self::zeros(self.n);
        for i in 0..self.n {
            for j in 0..self.n {
                out.set(j, i, self.get(i, j).conj());
            }
        }
        out
    }

    /// Matrix addition A + B.
    pub fn add(&self, other: &Self) -> Result<Self> {
        if self.n != other.n {
            return Err(error::invalid_param(
                "other",
                "matrix dimensions must match for addition",
            ));
        }
        let mut out = Self::zeros(self.n);
        for k in 0..self.data.len() {
            out.data[k] = self.data[k].add(&other.data[k]);
        }
        Ok(out)
    }

    /// Matrix subtraction A − B.
    pub fn sub(&self, other: &Self) -> Result<Self> {
        if self.n != other.n {
            return Err(error::invalid_param(
                "other",
                "matrix dimensions must match for subtraction",
            ));
        }
        let mut out = Self::zeros(self.n);
        for k in 0..self.data.len() {
            out.data[k] = self.data[k].sub(&other.data[k]);
        }
        Ok(out)
    }

    /// Scalar multiplication s·A.
    pub fn scale(&self, s: Complex) -> Self {
        let mut out = self.clone();
        for k in 0..out.data.len() {
            out.data[k] = out.data[k].mul(&s);
        }
        out
    }

    /// Scale by a real scalar.
    pub fn scale_real(&self, s: f64) -> Self {
        self.scale(Complex::from_real(s))
    }

    /// Matrix multiplication A·B (O(N³)).
    pub fn matmul(&self, other: &Self) -> Result<Self> {
        if self.n != other.n {
            return Err(error::invalid_param(
                "other",
                "matrix dimensions must match for multiplication",
            ));
        }
        let n = self.n;
        let mut out = Self::zeros(n);
        for i in 0..n {
            for k in 0..n {
                let aik = self.get(i, k);
                if aik.re == 0.0 && aik.im == 0.0 {
                    continue;
                }
                for j in 0..n {
                    let v = aik.mul(&other.get(k, j));
                    out.add_to(i, j, v);
                }
            }
        }
        Ok(out)
    }

    /// Frobenius norm √(Σ|M\[i\]\[j\]|²).
    pub fn frobenius_norm(&self) -> f64 {
        self.data.iter().map(|c| c.norm_sq()).sum::<f64>().sqrt()
    }

    /// Inverse via Gauss-Jordan elimination with partial pivoting.
    ///
    /// Returns an error if the matrix is singular (pivot < `1e-12 * max_element`).
    pub fn inverse(&self) -> Result<Self> {
        let n = self.n;
        // Augmented matrix [A | I]
        let mut aug: Vec<Vec<Complex>> = (0..n)
            .map(|i| {
                let mut row: Vec<Complex> = (0..n).map(|j| self.get(i, j)).collect();
                for j in 0..n {
                    row.push(if i == j { Complex::ONE } else { Complex::ZERO });
                }
                row
            })
            .collect();

        // Track scale for pivot threshold
        let max_elem = aug
            .iter()
            .flat_map(|row| row.iter().take(n))
            .map(|c| c.norm())
            .fold(0.0_f64, f64::max);
        let pivot_thresh = 1e-14 * max_elem.max(1.0);

        for col in 0..n {
            // Find pivot
            let mut max_row = col;
            let mut max_val = aug[col][col].norm();
            for (row, aug_row) in aug.iter().enumerate().skip(col + 1) {
                let v = aug_row[col].norm();
                if v > max_val {
                    max_val = v;
                    max_row = row;
                }
            }
            if max_val < pivot_thresh {
                return Err(error::numerical_error(
                    "matrix is singular or nearly singular (Gauss-Jordan inverse)",
                ));
            }
            aug.swap(col, max_row);

            // Scale pivot row
            let pivot = aug[col][col];
            let pivot_inv = Complex::ONE.div(&pivot);
            for v in &mut aug[col] {
                *v = v.mul(&pivot_inv);
            }

            // Eliminate column
            for row in 0..n {
                if row == col {
                    continue;
                }
                let factor = aug[row][col];
                if factor.re == 0.0 && factor.im == 0.0 {
                    continue;
                }
                let col_vals: Vec<Complex> = aug[col].clone();
                for (dst, src) in aug[row].iter_mut().zip(col_vals.iter()) {
                    *dst = dst.sub(&factor.mul(src));
                }
            }
        }

        // Extract right half as result
        let mut result = Self::zeros(n);
        for (i, aug_row) in aug.iter().enumerate() {
            for (j, val) in aug_row[n..].iter().enumerate() {
                result.set(i, j, *val);
            }
        }
        Ok(result)
    }

    /// Hermitian eigendecomposition via the cyclic Jacobi eigenvalue algorithm.
    ///
    /// Assumes `self` is Hermitian (A = A†). Repeatedly applies complex Jacobi
    /// rotations to annihilate the largest off-diagonal pair, converging to a diagonal
    /// form:
    ///
    /// 1. For each off-diagonal `A[p][q]` (`p < q`), first apply a diagonal phase
    ///    correction (`column q *= e^{-iφ}`, `row q *= e^{+iφ}`, `φ = arg(A[p][q])`) so
    ///    the pivot becomes real; this reduces the complex Hermitian problem to a real
    ///    symmetric 2×2 rotation at every step.
    /// 2. Apply the standard real Jacobi rotation (stable `tan(θ)` formula, Numerical
    ///    Recipes §11.1) to zero the now-real pivot, updating the two affected
    ///    rows/columns and accumulating the rotation into the eigenvector matrix.
    /// 3. Repeat in cyclic sweeps over all `(p, q)` pairs until the off-diagonal
    ///    Frobenius norm falls below a relative tolerance.
    ///
    /// Returns `(eigenvalues, eigenvectors)` where eigenvalues are sorted ascending and
    /// the k-th column of `eigenvectors` is the k-th eigenvector.
    ///
    /// # Errors
    ///
    /// Returns `NumericalError` if the sweep loop does not converge within a generous
    /// fixed number of sweeps.
    pub fn hermitian_eigendecomposition(&self) -> Result<(Vec<f64>, Self)> {
        hermitian_eig_impl(self)
    }

    /// Build a diagonal matrix from real values.
    pub fn from_diagonal(diag: &[f64]) -> Self {
        let n = diag.len();
        let mut m = Self::zeros(n);
        for (i, &v) in diag.iter().enumerate() {
            m.set(i, i, Complex::from_real(v));
        }
        m
    }

    /// Returns the column `j` as a `Vec<Complex>`.
    pub fn column(&self, j: usize) -> Vec<Complex> {
        (0..self.n).map(|i| self.get(i, j)).collect()
    }

    /// Returns the row `i` as a `Vec<Complex>`.
    pub fn row(&self, i: usize) -> Vec<Complex> {
        (0..self.n).map(|j| self.get(i, j)).collect()
    }
}

// ---------------------------------------------------------------------------
// Hermitian eigendecomposition: cyclic Jacobi eigenvalue algorithm
// ---------------------------------------------------------------------------
//
// NOTE (historical): an earlier version of this file implemented Householder
// tridiagonalization + implicit-QL (Golub & Van Loan §8.3-8.4). That implementation
// computed correct eigenVALUES but subtly incorrect eigenVECTORS for every n >= 3 input
// with a nonzero off-diagonal: each Householder step re-labelled its produced
// sub-diagonal magnitude as +sigma (always non-negative, for the real symmetric
// tridiagonal form expected by QL) while the accumulated unitary Q was built from the
// reflector's *actual* (possibly differently-signed/phased) output, with no compensating
// phase correction applied to Q at intermediate steps (only the very last sub-diagonal
// had such a correction). The result: Q did not actually satisfy `A = Q T Q†`, so
// `hermitian_eigendecomposition` returned eigenvectors failing `H v_k = λ_k v_k` for any
// non-trivial (non-diagonal) matrix of size >= 3, even though `hermitian_eig_impl`'s
// *eigenvalues* were exactly correct (the tridiagonal QL recursion for eigenvalues alone
// does not depend on Q at all). This was never caught by this file's own tests (which
// only checked eigenvalue correctness for n<=3, and eigenvector orthonormality — not the
// eigenvalue equation itself — at n=2), nor exercised by any consumer module that only
// used eigenvalues. It was found while cross-checking `frustrated::rvb`'s exact-
// diagonalization ground state against the total-spin-squared invariant `S_tot^2 = 0`,
// which requires a genuinely correct eigenVECTOR, not just the right eigenvalue.
//
// The cyclic Jacobi algorithm below has no such tridiagonalization phase-bookkeeping to
// get wrong: every rotation is applied *symmetrically and consistently* to the working
// matrix and the accumulated eigenvector matrix in the same step, so `V† A V → diagonal`
// holds by construction at every intermediate sweep, not just at convergence.

/// Cyclic Jacobi eigenvalue algorithm for a complex Hermitian matrix.
///
/// For each off-diagonal pair `(p, q)`, `p < q`:
/// 1. **Phase pre-rotation**: with `A[p][q] = r·e^{iφ}`, scale column `q` by `e^{-iφ}`
///    and row `q` by `e^{+iφ}` (and accumulate the same into the eigenvector matrix `V`
///    via its column `q`). This makes the pivot real (`= r >= 0`) without touching any
///    other pivot, reducing the complex case to a real symmetric 2×2 rotation.
/// 2. **Real Jacobi rotation**: the standard stable `tan(θ)` formula (Numerical Recipes
///    §11.1, `JACOBI`) zeroes the now-real pivot `A[p][q]`, updating `A[p][p]`, `A[q][q]`,
///    every other row/column `k ∉ {p,q}` (`A[k][p], A[k][q]` and their Hermitian mirrors),
///    and accumulating the same rotation into `V`.
///
/// Repeated in cyclic sweeps (all `p<q` pairs) until the off-diagonal Frobenius norm
/// falls below a relative tolerance. Returns `(eigenvalues, eigenvectors)` with
/// eigenvalues sorted ascending and the k-th column of `eigenvectors` the k-th
/// eigenvector — by construction, `V` accumulates exactly the product of all applied
/// unitary rotations, so `A_original = V T V†` (T diagonal) holds at every step, not just
/// asymptotically.
///
/// # Errors
///
/// Returns `NumericalError` if the sweep loop does not converge within `100` sweeps.
///
/// `p`, `q`, `k`, `row`/`col` index simultaneously into both the working matrix `a` and
/// the accumulated eigenvector matrix `v` throughout — clippy's `needless_range_loop` is
/// a false positive for this genuinely index-driven (not merely positional) algorithm.
#[allow(clippy::needless_range_loop)]
fn hermitian_eig_impl(h: &CMatrix) -> Result<(Vec<f64>, CMatrix)> {
    let n = h.n;
    if n == 0 {
        return Ok((vec![], CMatrix::zeros(0)));
    }
    if n == 1 {
        let e = h.get(0, 0).re;
        let mut v = CMatrix::zeros(1);
        v.set(0, 0, Complex::ONE);
        return Ok((vec![e], v));
    }

    // Work copy of A (complex, n×n, dense) and the accumulated eigenvector matrix V
    // (starts as identity).
    let mut a: Vec<Vec<Complex>> = (0..n)
        .map(|i| (0..n).map(|j| h.get(i, j)).collect())
        .collect();
    let mut v: Vec<Vec<Complex>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| if i == j { Complex::ONE } else { Complex::ZERO })
                .collect()
        })
        .collect();

    let frobenius_scale = a
        .iter()
        .flatten()
        .map(|c| c.norm_sq())
        .sum::<f64>()
        .sqrt()
        .max(1.0);
    let tol = 1e-14 * frobenius_scale;
    const MAX_SWEEPS: usize = 100;

    let mut converged = false;
    for _sweep in 0..MAX_SWEEPS {
        let mut off_diag_sq = 0.0_f64;
        for p in 0..n {
            for q in (p + 1)..n {
                off_diag_sq += a[p][q].norm_sq();
            }
        }
        if off_diag_sq.sqrt() < tol {
            converged = true;
            break;
        }

        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p][q];
                if apq.norm() < 1e-300 {
                    continue;
                }

                // Phase pre-rotation: make the (p,q) pivot real and non-negative.
                let phi = apq.phase();
                let phase_q = Complex::from_polar(1.0, -phi);
                for row in a.iter_mut() {
                    row[q] = row[q].mul(&phase_q);
                }
                let phase_q_conj = phase_q.conj();
                for col in 0..n {
                    a[q][col] = a[q][col].mul(&phase_q_conj);
                }
                for row in v.iter_mut() {
                    row[q] = row[q].mul(&phase_q);
                }

                let app = a[p][p].re;
                let aqq = a[q][q].re;
                let r = a[p][q].re;
                if r.abs() < 1e-300 {
                    continue;
                }

                // Standard real-symmetric Jacobi rotation: t = tan(theta) solves
                // r*t^2 - (aqq-app)*t - r = 0 (the root with |t| <= 1, for R = [[c,-s],[s,c]]
                // and new_a[p][p] = app + t*r, new_a[q][q] = aqq - t*r — both derived and
                // verified directly by hand against explicit 2x2 and 3x3 cases). Computed
                // via the numerically stable "large root, then reciprocal" form to avoid
                // cancellation (the magnitude matches Numerical Recipes §11.1; the sign
                // here is fixed to match *this* file's R/diagonal-update convention, which
                // is the opposite of NR's own d[p]-=h/d[q]+=h convention).
                let tau = (aqq - app) / (2.0 * r);
                let t = if tau >= 0.0 {
                    -1.0 / (tau + (1.0 + tau * tau).sqrt())
                } else {
                    1.0 / (-tau + (1.0 + tau * tau).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;

                a[p][p] = Complex::from_real(app + t * r);
                a[q][q] = Complex::from_real(aqq - t * r);
                a[p][q] = Complex::ZERO;
                a[q][p] = Complex::ZERO;

                // Rotate every other row/column k (both Hermitian mirror entries).
                for k in 0..n {
                    if k == p || k == q {
                        continue;
                    }
                    let akp = a[k][p];
                    let akq = a[k][q];
                    let new_akp = akp.scale(c).add(&akq.scale(s));
                    let new_akq = akq.scale(c).sub(&akp.scale(s));
                    a[k][p] = new_akp;
                    a[p][k] = new_akp.conj();
                    a[k][q] = new_akq;
                    a[q][k] = new_akq.conj();
                }

                // Accumulate the same rotation into the eigenvector matrix V.
                for row in v.iter_mut() {
                    let vip = row[p];
                    let viq = row[q];
                    row[p] = vip.scale(c).add(&viq.scale(s));
                    row[q] = viq.scale(c).sub(&vip.scale(s));
                }
            }
        }
    }

    if !converged {
        return Err(crate::error::numerical_error(
            "Jacobi eigenvalue sweep did not converge within 100 sweeps",
        ));
    }

    // Sort eigenvalues ascending and reorder eigenvectors accordingly.
    let diag: Vec<f64> = (0..n).map(|i| a[i][i].re).collect();
    let mut pairs: Vec<(f64, usize)> = diag
        .iter()
        .copied()
        .enumerate()
        .map(|(i, e)| (e, i))
        .collect();
    pairs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));

    let eigenvalues: Vec<f64> = pairs.iter().map(|(e, _)| *e).collect();
    let mut eigenvectors = CMatrix::zeros(n);
    for (col_out, &(_, col_in)) in pairs.iter().enumerate() {
        for row in 0..n {
            eigenvectors.set(row, col_out, v[row][col_in]);
        }
    }

    Ok((eigenvalues, eigenvectors))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn cx(re: f64, im: f64) -> Complex {
        Complex::new(re, im)
    }

    #[test]
    fn test_eye_trace() {
        let m = CMatrix::eye(4);
        assert!((m.trace().re - 4.0).abs() < 1e-14);
        assert!((m.trace().im).abs() < 1e-14);
    }

    #[test]
    fn test_matmul_identity() {
        let a = CMatrix::from_rows(vec![
            vec![cx(1.0, 0.0), cx(2.0, 1.0)],
            vec![cx(3.0, -1.0), cx(4.0, 0.0)],
        ])
        .unwrap();
        let i = CMatrix::eye(2);
        let r = a.matmul(&i).unwrap();
        for row in 0..2 {
            for col in 0..2 {
                let diff = r.get(row, col).sub(&a.get(row, col));
                assert!(diff.norm() < 1e-13);
            }
        }
    }

    #[test]
    fn test_matmul_2x2() {
        let a = CMatrix::from_rows(vec![
            vec![cx(1.0, 0.0), cx(0.0, 1.0)],
            vec![cx(0.0, -1.0), cx(1.0, 0.0)],
        ])
        .unwrap();
        let r = a.matmul(&a).unwrap();
        // a^2 where a = σ_y (Pauli Y)
        // σ_y = [[0,−i],[i,0]] scaled; here [[1,i],[−i,1]] — just verify trace
        assert!(r.trace().re.is_finite());
    }

    #[test]
    fn test_conj_transpose() {
        let a = CMatrix::from_rows(vec![
            vec![cx(1.0, 2.0), cx(3.0, 4.0)],
            vec![cx(5.0, 6.0), cx(7.0, 8.0)],
        ])
        .unwrap();
        let ah = a.conj_transpose();
        assert!((ah.get(0, 1).re - 5.0).abs() < 1e-14);
        assert!((ah.get(0, 1).im + 6.0).abs() < 1e-14);
    }

    #[test]
    fn test_inverse_2x2() {
        let a = CMatrix::from_rows(vec![
            vec![cx(2.0, 0.0), cx(1.0, 0.0)],
            vec![cx(1.0, 0.0), cx(1.0, 0.0)],
        ])
        .unwrap();
        let inv = a.inverse().unwrap();
        let prod = a.matmul(&inv).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(approx_eq(prod.get(i, j).re, expected, 1e-12));
                assert!(approx_eq(prod.get(i, j).im, 0.0, 1e-12));
            }
        }
    }

    #[test]
    fn test_inverse_complex() {
        let a = CMatrix::from_rows(vec![
            vec![cx(1.0, 1.0), cx(0.0, 1.0)],
            vec![cx(0.0, -1.0), cx(1.0, -1.0)],
        ])
        .unwrap();
        let inv = a.inverse().unwrap();
        let prod = a.matmul(&inv).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(approx_eq(prod.get(i, j).re, expected, 1e-12));
            }
        }
    }

    #[test]
    fn test_inverse_singular_errors() {
        let a = CMatrix::from_rows(vec![
            vec![cx(1.0, 0.0), cx(2.0, 0.0)],
            vec![cx(2.0, 0.0), cx(4.0, 0.0)],
        ])
        .unwrap();
        assert!(a.inverse().is_err());
    }

    #[test]
    fn test_hermitian_eigendecomposition_identity() {
        let h = CMatrix::eye(3);
        let (vals, _vecs) = h.hermitian_eigendecomposition().unwrap();
        for v in &vals {
            assert!(approx_eq(*v, 1.0, 1e-10));
        }
    }

    #[test]
    fn test_hermitian_eigendecomposition_2x2_real() {
        // H = [[2, 1],[1, 2]] → eigenvalues 1 and 3
        let h = CMatrix::from_rows(vec![
            vec![cx(2.0, 0.0), cx(1.0, 0.0)],
            vec![cx(1.0, 0.0), cx(2.0, 0.0)],
        ])
        .unwrap();
        let (vals, vecs) = h.hermitian_eigendecomposition().unwrap();
        assert!(approx_eq(vals[0], 1.0, 1e-10));
        assert!(approx_eq(vals[1], 3.0, 1e-10));
        // Verify H v = λ v: compute H * v0 via matmul with vecs
        // (eigenvectors are columns of vecs)
        let hv = h.matmul(&vecs).unwrap();
        // Column 0 of H*vecs should equal vals[0] * column 0 of vecs
        for i in 0..2 {
            let hv_i = hv.get(i, 0);
            let lv_i = vecs.get(i, 0).scale(vals[0]);
            assert!(approx_eq(hv_i.re, lv_i.re, 1e-9));
            assert!(approx_eq(hv_i.im, lv_i.im, 1e-9));
        }
    }

    #[test]
    fn test_hermitian_eigendecomposition_complex_offdiag() {
        // H = [[1, i],[-i, 1]] → eigenvalues 0 and 2
        let h = CMatrix::from_rows(vec![
            vec![cx(1.0, 0.0), cx(0.0, 1.0)],
            vec![cx(0.0, -1.0), cx(1.0, 0.0)],
        ])
        .unwrap();
        let (vals, _) = h.hermitian_eigendecomposition().unwrap();
        assert!(approx_eq(vals[0], 0.0, 1e-9));
        assert!(approx_eq(vals[1], 2.0, 1e-9));
    }

    #[test]
    fn test_eigendecomposition_3x3_diagonal() {
        let h = CMatrix::from_diagonal(&[3.0, 1.0, 2.0]);
        let (vals, _) = h.hermitian_eigendecomposition().unwrap();
        assert!(approx_eq(vals[0], 1.0, 1e-10));
        assert!(approx_eq(vals[1], 2.0, 1e-10));
        assert!(approx_eq(vals[2], 3.0, 1e-10));
    }

    #[test]
    fn test_hermitian_eigendecomposition_satisfies_eigenvalue_equation_3x3() {
        // Regression test: an earlier Householder+QL implementation of
        // `hermitian_eigendecomposition` computed correct eigenVALUES but subtly
        // incorrect eigenVECTORS for every non-diagonal n>=3 Hermitian matrix (a
        // sign/phase bookkeeping bug never caught by this file's own tests, which only
        // checked eigenvalues for n<=3 or eigenvector orthonormality — not the
        // eigenvalue equation itself — at n=2). This 3x3 matrix has an asymmetric
        // diagonal (2, 1.5, 3) specifically because a symmetric-diagonal 2x2/3x3 case
        // cannot distinguish a particular rotation-sign convention error from a correct
        // one (both give the same result when the two diagonal entries being rotated
        // are equal).
        let h = CMatrix::from_rows(vec![
            vec![cx(2.0, 0.0), cx(0.5, 0.0), cx(0.3, 0.0)],
            vec![cx(0.5, 0.0), cx(1.5, 0.0), cx(0.2, 0.0)],
            vec![cx(0.3, 0.0), cx(0.2, 0.0), cx(3.0, 0.0)],
        ])
        .unwrap();
        let (vals, vecs) = h.hermitian_eigendecomposition().unwrap();
        // Eigenvalues ascending.
        assert!(vals[0] < vals[1] && vals[1] < vals[2]);
        // H * v_k = lambda_k * v_k for every column k (the actual eigenvalue equation,
        // not merely orthonormality or a trace/eigenvalue check).
        let hv = h.matmul(&vecs).unwrap();
        // `col` indexes `vals`, `vecs` (via `.get`), and `hv` (via `.get`) simultaneously.
        #[allow(clippy::needless_range_loop)]
        for col in 0..3 {
            for row in 0..3 {
                let lhs = hv.get(row, col);
                let rhs = vecs.get(row, col).scale(vals[col]);
                assert!(
                    approx_eq(lhs.re, rhs.re, 1e-9) && approx_eq(lhs.im, rhs.im, 1e-9),
                    "H v_{} != lambda_{} v_{} at row {}: {:?} vs {:?}",
                    col,
                    col,
                    col,
                    row,
                    lhs,
                    rhs
                );
            }
        }
        // Eigenvectors orthonormal.
        for a in 0..3 {
            for b in 0..3 {
                let dot: Complex = (0..3)
                    .map(|r| vecs.get(r, a).conj().mul(&vecs.get(r, b)))
                    .fold(Complex::ZERO, |acc, c| acc.add(&c));
                let expected = if a == b { 1.0 } else { 0.0 };
                assert!(approx_eq(dot.re, expected, 1e-9) && approx_eq(dot.im, 0.0, 1e-9));
            }
        }
    }

    #[test]
    fn test_eigendecomposition_eigenvectors_orthonormal() {
        let h = CMatrix::from_rows(vec![
            vec![cx(2.0, 0.0), cx(1.0, 0.0)],
            vec![cx(1.0, 0.0), cx(2.0, 0.0)],
        ])
        .unwrap();
        let (_, vecs) = h.hermitian_eigendecomposition().unwrap();
        // Check v0 · v1 = 0, |v0| = 1, |v1| = 1
        let v0: Vec<Complex> = (0..2).map(|i| vecs.get(i, 0)).collect();
        let v1: Vec<Complex> = (0..2).map(|i| vecs.get(i, 1)).collect();
        let dot: Complex = v0
            .iter()
            .zip(v1.iter())
            .map(|(a, b)| a.conj().mul(b))
            .fold(Complex::ZERO, |acc, x| acc.add(&x));
        assert!(dot.norm() < 1e-10);
        let norm0_sq: f64 = v0.iter().map(|c| c.norm_sq()).sum();
        assert!(approx_eq(norm0_sq, 1.0, 1e-10));
    }

    #[test]
    fn test_from_rows_inconsistent_size() {
        let result = CMatrix::from_rows(vec![vec![cx(1.0, 0.0), cx(0.0, 0.0)], vec![cx(0.0, 0.0)]]);
        assert!(result.is_err());
    }

    #[test]
    fn test_add_and_sub() {
        let a = CMatrix::eye(2);
        let b = CMatrix::eye(2);
        let s = a.add(&b).unwrap();
        assert!((s.get(0, 0).re - 2.0).abs() < 1e-14);
        let d = a.sub(&b).unwrap();
        assert!((d.get(0, 0).re).abs() < 1e-14);
    }

    #[test]
    fn test_scale() {
        let m = CMatrix::eye(3);
        let s = m.scale(Complex::new(2.0, 0.0));
        assert!((s.get(0, 0).re - 2.0).abs() < 1e-14);
        assert!((s.get(0, 1).re).abs() < 1e-14);
    }

    #[test]
    fn test_column_and_row() {
        let a = CMatrix::from_rows(vec![
            vec![cx(1.0, 0.0), cx(2.0, 0.0)],
            vec![cx(3.0, 0.0), cx(4.0, 0.0)],
        ])
        .unwrap();
        let col = a.column(1);
        assert!((col[0].re - 2.0).abs() < 1e-14);
        assert!((col[1].re - 4.0).abs() < 1e-14);
        let row = a.row(0);
        assert!((row[0].re - 1.0).abs() < 1e-14);
        assert!((row[1].re - 2.0).abs() < 1e-14);
    }
}
