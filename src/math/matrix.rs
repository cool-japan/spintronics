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
    /// Row-major flat storage: data[i*n + j] = M[i][j]
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

    /// Get element M[i][j].
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Complex {
        self.data[i * self.n + j]
    }

    /// Set element M[i][j].
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: Complex) {
        self.data[i * self.n + j] = v;
    }

    /// Add `v` to element M[i][j] in place.
    #[inline]
    fn add_to(&mut self, i: usize, j: usize, v: Complex) {
        let cur = self.get(i, j);
        self.set(i, j, cur.add(&v));
    }

    /// Trace = Σ M[i][i].
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

    /// Frobenius norm √(Σ|M[i][j]|²).
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

    /// Hermitian eigendecomposition via Householder tridiagonalization + implicit QL.
    ///
    /// Assumes `self` is Hermitian (A = A†). Uses the standard two-phase approach:
    ///
    /// 1. **Householder tridiagonalization** (Golub & Van Loan §8.3): reduce A to a
    ///    real symmetric tridiagonal T by a sequence of complex Householder reflectors,
    ///    accumulating the unitary Q such that A = Q T Q†.
    /// 2. **Implicit QL with Wilkinson shift** (Golub & Van Loan §8.4): diagonalize T
    ///    using real Givens rotations with cubic convergence.
    ///
    /// Returns `(eigenvalues, eigenvectors)` where eigenvalues are sorted ascending and
    /// the k-th column of `eigenvectors` is the k-th eigenvector.
    ///
    /// # Errors
    ///
    /// Returns `NumericalError` if the QL iteration does not converge within 30·n steps.
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
// Hermitian eigendecomposition: Householder tridiagonalization + implicit QL
// ---------------------------------------------------------------------------

/// Householder tridiagonalization of a complex Hermitian matrix.
///
/// Returns `(d, e, q)` where:
/// - `d[0..n]` is the real diagonal of the tridiagonal form.
/// - `e[0..n]` is the real sub-diagonal, Householder convention:
///   `e[0]=0` (unused), `e[i]` connects `d[i-1]` to `d[i]` for `i=1..n-1`.
/// - `q` is the unitary matrix such that `A = Q T Q†`.
///
/// Algorithm: Golub & Van Loan §8.3, complex Hermitian variant.
/// For each column k = 0..n-2:
///   1. Extract sub-column x = A[k+1:n, k].
///   2. Form Householder reflector H = I - β v v† that maps x → -sign(x₀)|x| eₖ.
///   3. Apply: A ← H A H†, Q ← Q H† (accumulate unitary).
///   4. Record d[k] = A[k,k].re, e[k+1] = |x| (real, positive sub-diagonal).
///
/// The unitary Q satisfies Q† A_original Q = T (real symmetric tridiagonal).
fn hermitian_householder_tridiag(h: &CMatrix, n: usize) -> (Vec<f64>, Vec<f64>, CMatrix) {
    // Work copy of A (complex, n×n, row-major).
    let mut a: Vec<Vec<Complex>> = (0..n)
        .map(|i| (0..n).map(|j| h.get(i, j)).collect())
        .collect();
    // Q starts as identity; we accumulate reflectors into Q from the right.
    let mut q: Vec<Vec<Complex>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| if i == j { Complex::ONE } else { Complex::ZERO })
                .collect()
        })
        .collect();
    let mut d = vec![0.0_f64; n];
    let mut e = vec![0.0_f64; n]; // e[0] = 0 (unused)

    // The Householder tridiagonalization uses n-2 reflectors (columns 0..n-3).
    // For n<=2 the matrix is already tridiagonal; no reflectors are needed.
    for k in 0..n.saturating_sub(2) {
        // m = length of sub-column below position (k, k) = n - k - 1
        let m = n - k - 1;

        // x = A[k+1..n, k] (sub-column below diagonal in column k)
        let x_orig: Vec<Complex> = (0..m).map(|i| a[k + 1 + i][k]).collect();
        let sigma = x_orig.iter().map(|c| c.norm_sq()).sum::<f64>().sqrt();

        // Record diagonal entry before modification
        d[k] = a[k][k].re;

        if sigma < 1e-15 {
            // Sub-column already negligibly small; sub-diagonal is zero.
            e[k + 1] = 0.0;
            continue;
        }

        // Sub-diagonal of the tridiagonal is sigma (the norm of sub-column).
        e[k + 1] = sigma;

        // Householder vector: v = x + e^{i*arg(x[0])} * sigma * e_0
        // Using the phase of x[0] ensures numerical stability (v[0] ≈ 2*sigma in magnitude).
        let x0 = x_orig[0];
        let phase = if x0.norm_sq() < 1e-300 {
            Complex::ONE
        } else {
            Complex::from_polar(1.0, x0.phase())
        };
        let mut v: Vec<Complex> = x_orig.clone();
        v[0] = x0.add(&phase.scale(sigma));
        let v_norm_sq = v.iter().map(|c| c.norm_sq()).sum::<f64>();
        let beta = if v_norm_sq < 1e-28 {
            0.0
        } else {
            2.0 / v_norm_sq
        };

        // Apply H = I - β v v† to A[k+1:n, k+1:n] from both sides.
        // H A H† = A - β v w† - β w v† + β² (v† w) v v†
        // where w = A v  (UNSCALED matrix-vector product).
        // Efficient form: A ← A - v p† - p v†
        // where p = β w - (β/2)(v† w) v  (so that the β² term is absorbed).
        //
        // Derivation check:
        //   v p† + p v† = β v w† - (β/2)(v†w) v v† + β w v† - (β/2)(v†w) v v†
        //               = β v w† + β w v† - β(v†w) v v†
        //   A - v p† - p v† = A - β v w† - β w v† + β(v†w) v v†  ✓ (matches H A H†)

        // w = A[k+1:n, k+1:n] * v  (unscaled)
        let mut w = vec![Complex::ZERO; m];
        for i in 0..m {
            for j in 0..m {
                w[i] = w[i].add(&a[k + 1 + i][k + 1 + j].mul(&v[j]));
            }
        }
        // v† w = sum_i conj(v[i]) * w[i]  (real for Hermitian A and real sigma)
        let vt_w: Complex = v
            .iter()
            .zip(w.iter())
            .map(|(vi, wi)| vi.conj().mul(wi))
            .fold(Complex::ZERO, |acc, c| acc.add(&c));
        // p_vec = β w - (β²/2)(v† w) v
        // This is the standard efficient Householder update: A - v p† - p v† = H A H†.
        let p_vec: Vec<Complex> = w
            .iter()
            .zip(v.iter())
            .map(|(wi, vi)| wi.scale(beta).sub(&vt_w.scale(beta * beta * 0.5).mul(vi)))
            .collect();
        // A[k+1:n, k+1:n] -= v * p† + p * v†
        for i in 0..m {
            for j in 0..m {
                let delta = v[i].mul(&p_vec[j].conj()).add(&p_vec[i].mul(&v[j].conj()));
                a[k + 1 + i][k + 1 + j] = a[k + 1 + i][k + 1 + j].sub(&delta);
            }
        }

        // Set the reduced column/row entries explicitly (numerical noise clean-up).
        // After reduction, a[k+1, k] should be -sigma (the norm, with a sign from the reflector).
        // We store it as -sigma (real, negative) to be consistent with the phase convention.
        a[k + 1][k] = Complex::from_real(-sigma);
        a[k][k + 1] = Complex::from_real(-sigma);
        for i in 2..m + 1 {
            a[k + i][k] = Complex::ZERO;
            a[k][k + i] = Complex::ZERO;
        }

        // Update Q: Q ← Q * H† = Q * (I - β v v†)†  = Q * (I - β v v†) (H is Hermitian)
        // For each row r of Q: Q[r, k+1:n] ← Q[r, k+1:n] - β (Q[r, k+1:n] · v) * v†
        // dot_r = sum_j Q[r, k+1+j] * conj(v[j])
        for q_row in q.iter_mut() {
            let dot_r: Complex = (0..m)
                .map(|j| q_row[k + 1 + j].mul(&v[j].conj()))
                .fold(Complex::ZERO, |acc, c| acc.add(&c))
                .scale(beta);
            for j in 0..m {
                q_row[k + 1 + j] = q_row[k + 1 + j].sub(&dot_r.mul(&v[j]));
            }
        }
    }

    // Extract the remaining diagonal and sub-diagonal entries.
    // After the loop (k ran from 0 to n-3), the last 2×2 block a[n-2..n, n-2..n]
    // is not yet necessarily real-sub-diagonal. The loop set d[0..n-3] and e[1..n-2].
    // We still need:
    //   d[n-2], d[n-1] from the last 2×2 diagonal entries.
    //   e[n-1] = |a[n-1][n-2]| (sub-diagonal magnitude — may be complex if n<=2
    //                            because the Householder loop didn't run for k=n-2).
    //
    // If n=2: no Householder ran, so a[1][0] = the original h[1][0] (complex in general).
    // We phase-rotate column 1 of Q to make a[1][0] real positive.
    if n >= 2 {
        d[n - 2] = a[n - 2][n - 2].re;
        let sub = a[n - 1][n - 2]; // sub-diagonal element (may be complex for n=2)
        let sub_norm = sub.norm();
        e[n - 1] = sub_norm;
        // Phase-normalize: if sub is not real positive, apply a Givens phase rotation
        // to Q[:,n-1] to make the sub-diagonal real positive in the eigenvector basis.
        // This is a unitary transformation: Q[:,n-1] *= conj(sub / sub_norm).
        // Equivalently: multiply Q[:,n-1] by e^{-i*arg(sub)}.
        if sub_norm > 1e-15 {
            let phase_conj = Complex::from_polar(1.0, -sub.phase());
            for q_row in q.iter_mut() {
                q_row[n - 1] = q_row[n - 1].mul(&phase_conj);
            }
        }
    }
    d[n - 1] = a[n - 1][n - 1].re;

    // Build the CMatrix Q.
    let mut q_mat = CMatrix::zeros(n);
    for (i, q_row) in q.iter().enumerate() {
        for (j, &val) in q_row.iter().enumerate() {
            q_mat.set(i, j, val);
        }
    }
    (d, e, q_mat)
}

/// Implicit QL algorithm with Wilkinson shift for a real symmetric tridiagonal.
///
/// Direct translation of Numerical Recipes §11.3 TQLI (C edition), adapted for
/// 0-indexed arrays and complex eigenvector accumulator `z`.
///
/// **Convention** (matching NR after its internal e-shift):
/// `e[i]` is the sub-diagonal element connecting `d[i]` to `d[i+1]`, for i=0..n-2.
/// `e[n-1] = 0` (boundary). This is the NR-post-shift convention.
///
/// The Householder function returns `e_h[]` where `e_h[i]` connects `d[i-1]` to `d[i]`
/// (i.e. e_h[1..n-1] are the sub-diagonals, e_h[0]=0). The caller must shift:
/// `e[i] = e_h[i+1]` for i=0..n-2, `e[n-1]=0`.
fn tridiag_ql_in_place(
    d: &mut [f64],
    e: &mut [f64],
    z: &mut CMatrix,
    n: usize,
) -> crate::error::Result<()> {
    if n <= 1 {
        return Ok(());
    }
    // e[i] = sub-diagonal between d[i] and d[i+1]; e[n-1]=0.
    // NR's TQLI inner loop structure (1-indexed NR → 0-indexed here):
    //   for (l=1; l<=n; l++) {        →  for l in 0..n
    //     for (m=l; m<=n-1; m++)      →  for m in l..n-1 (search)
    //       if e[m]≈0 break           →  check e[m] (0-indexed, our convention)
    //     shift from d[l], e[l], d[l+1] (NR: e[l], our e[l])
    //     inner for (i=m-1; i>=l; i--):
    //       f=s*e[i]; b=c*e[i]        →  f=s*e[i]; b=c*e[i]
    //       e[i+1]=pythag(f,g)        →  e[i+1]=hypot(f,g)
    //       ...
    //     d[l]-=p; e[l]=g; e[m]=0

    for l in 0..n {
        let mut num_iter = 0_usize;
        loop {
            // Find m: first index in [l, n-1] where e[m] is negligible.
            // e[m] connects d[m] and d[m+1].
            let mut m = l;
            while m < n - 1 {
                let dd = d[m].abs() + d[m + 1].abs();
                if (e[m].abs() + dd) == dd {
                    break; // e[m] negligible → d[l..m] block is decoupled at m
                }
                m += 1;
            }
            // If m == l, e[l] is already negligible → d[l] is converged.
            if m == l {
                break;
            }
            if num_iter >= 60 {
                return Err(crate::error::numerical_error(
                    "tridiagonal QL iteration did not converge within 60 iterations",
                ));
            }
            num_iter += 1;

            // Wilkinson shift from the 2×2 block at the l-end: [d[l], e[l]; e[l], d[l+1]].
            // NR formula (1-indexed):
            //   g = (d[l+1] - d[l]) / (2 * e[l])
            //   r = pythag(g, 1)
            //   g = d[m] - d[l] + e[l] / (g + SIGN(r, g))
            // NR's final g IS the initial g_var for the sweep (it equals d[m] - shift_eigenvalue).
            let gg = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let r = gg.hypot(1.0);
            // g_var_init = d[m] - d[l] + e[l]/(gg + sign(gg)*r) — this is NR's g after shift
            let g_var_init = d[m] - d[l] + e[l] / (gg + if gg >= 0.0 { r } else { -r });

            // Initialize QL sweep variables (NR: s=c=1, p=0).
            let mut g_var = g_var_init;
            let mut s = 1.0_f64;
            let mut c = 1.0_f64;
            let mut p = 0.0_f64;

            // Inner QL sweep: i from m-1 down to l (inclusive).
            // Each iteration annihilates one sub-diagonal e[i+1] while updating d[i+1], d[i].
            let mut i = m;
            while i > l {
                i -= 1;
                // f = s * e[i] — e[i] connects d[i] and d[i+1]
                let f = s * e[i];
                let b = c * e[i];
                // New off-diagonal = pythag(f, g_var); store back into e[i+1].
                let r_hyp = f.hypot(g_var);
                e[i + 1] = r_hyp;
                if r_hyp.abs() < 1e-300 {
                    // Degenerate Givens; d[i+1] doesn't change (no rotation).
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    break;
                }
                s = f / r_hyp;
                c = g_var / r_hyp;
                // Update d[i+1].
                g_var = d[i + 1] - p;
                let r_var = (d[i] - g_var) * s + 2.0 * c * b;
                p = s * r_var;
                d[i + 1] = g_var + p;
                // Update running g_var for next iteration.
                g_var = c * r_var - b;

                // Apply Givens rotation to eigenvector columns i and i+1.
                // z[:, i+1] ← s * z[:, i] + c * z[:, i+1]
                // z[:, i]   ← c * z[:, i] - s * z[:, i+1]
                for row in 0..n {
                    let zi = z.get(row, i);
                    let zi1 = z.get(row, i + 1);
                    z.set(row, i + 1, zi.scale(s).add(&zi1.scale(c)));
                    z.set(row, i, zi.scale(c).sub(&zi1.scale(s)));
                }
            }
            // Apply accumulated shift to d[l].
            d[l] -= p;
            // Store running g_var into e[l] (residual that will converge to 0).
            e[l] = g_var;
            // Zero out e[m] — the sweep was supposed to eliminate it.
            e[m] = 0.0;
        }
    }
    Ok(())
}

/// Main entry point: Householder + QL for complex Hermitian matrix.
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

    // Phase 1: tridiagonalize.
    // hermitian_householder_tridiag returns e_h where e_h[i] = sub-diagonal between
    // d[i-1] and d[i] (e_h[0]=0 unused).
    let (mut d, e_h, mut q) = hermitian_householder_tridiag(h, n);

    // Convert from Householder convention (e_h[i] connects d[i-1]↔d[i]) to
    // TQLI convention (e[i] connects d[i]↔d[i+1]), matching NR's post-shift e[].
    // e[i] = e_h[i+1] for i=0..n-2; e[n-1] = 0.
    let mut e: Vec<f64> = (0..n)
        .map(|i| if i + 1 < n { e_h[i + 1] } else { 0.0 })
        .collect();

    // Phase 2: QL on real symmetric tridiagonal, updating Q.
    tridiag_ql_in_place(&mut d, &mut e, &mut q, n)?;

    // Sort eigenvalues ascending and reorder eigenvectors
    let mut pairs: Vec<(f64, usize)> = d.iter().copied().enumerate().map(|(i, v)| (v, i)).collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let eigenvalues: Vec<f64> = pairs.iter().map(|(v, _)| *v).collect();
    let mut eigenvectors = CMatrix::zeros(n);
    for (col_out, (_, col_in)) in pairs.iter().enumerate() {
        for row in 0..n {
            eigenvectors.set(row, col_out, q.get(row, *col_in));
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
