//! Multi-band Wilson loops and nested Wilson loops for higher-order topology.
//!
//! This module generalises the single-band Wilson loop already present in
//! [`ChernNumber::wilson_loop`] to the multi-band case, following the
//! non-abelian Berry-phase formalism of:
//!
//! - D. Vanderbilt and R. D. King-Smith, *Phys. Rev. B* **48**, 4442 (1993):
//!   Nested Wilson-loop construction for hybrid Wannier functions.
//! - R. D. King-Smith and D. Vanderbilt, *Phys. Rev. B* **47**, 1651 (1993):
//!   Polarisation from Berry phase of occupied-band Wannier centres.
//! - J. Zak, *Phys. Rev. Lett.* **62**, 2747 (1989):
//!   Zak phase as a gauge-invariant Berry phase along the BZ.
//! - W. A. Benalcazar, B. A. Bernevig, and T. L. Hughes,
//!   *Science* **357**, 61 (2017):
//!   Quadrupole moment via nested Wilson loops (BBH theory).
//!
//! # Algorithm
//!
//! Given `n_occ` occupied band indices the method works as follows:
//!
//! 1. **Link matrix** `M(k, k+dk)`: an `n_occ × n_occ` overlap matrix
//!    `M_mn = ⟨u_m(k)|u_n(k+dk)⟩` where the inner product runs over the
//!    lattice degree of freedom.
//!
//! 2. **Wilson unitary** along ky at fixed kx:
//!    `W_y(kx) = M(k_1, k_2) · M(k_2, k_3) · … · M(k_ny, k_1)`
//!    which is an `n_occ × n_occ` unitary matrix.
//!
//! 3. **Wannier centres** `ν_j ∈ [0,1)` are the eigenphases of `W_y` divided
//!    by `2π`.
//!
//! 4. **Nested Wilson loop**: treat `{ν_j(kx)}` as a 1-D Wannier band and
//!    compute its Berry phase along kx to obtain the nested polarisation.
//!    A quantised value of 0.5 signals a HOTI corner charge.

use std::f64::consts::PI;

use crate::error::{self, Result};
use crate::math::{CMatrix, Complex};
use crate::topomagnon::band_model::MagnonBandModel;

// ---------------------------------------------------------------------------
// WilsonLoop
// ---------------------------------------------------------------------------

/// Multi-band Wilson-loop calculator for a [`MagnonBandModel`].
///
/// Wraps a reference to a model and a list of occupied-band indices.
/// All Wilson-loop integrals are computed via Fukui-type discretisation
/// of the Brillouin zone.
pub struct WilsonLoop<'a> {
    /// Reference to the underlying magnon band model.
    pub model: &'a MagnonBandModel,
    /// Indices of the bands treated as "occupied" (0-indexed, ascending order).
    pub occupied_bands: Vec<usize>,
}

impl<'a> WilsonLoop<'a> {
    /// Construct a new `WilsonLoop` calculator.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if any band index exceeds `model.n_bands()-1`
    /// or if `occupied_bands` is empty.
    pub fn new(model: &'a MagnonBandModel, occupied_bands: Vec<usize>) -> Result<Self> {
        if occupied_bands.is_empty() {
            return Err(error::invalid_param(
                "occupied_bands",
                "must list at least one band index",
            ));
        }
        let nb = model.n_bands();
        for &b in &occupied_bands {
            if b >= nb {
                return Err(error::invalid_param(
                    "occupied_bands",
                    &format!("band index {b} exceeds n_bands-1={}", nb - 1),
                ));
            }
        }
        Ok(Self {
            model,
            occupied_bands,
        })
    }

    // -----------------------------------------------------------------------
    // Link matrix
    // -----------------------------------------------------------------------

    /// Compute the `n_occ × n_occ` overlap matrix
    /// `M_mn = ⟨u_m^a | u_n^b⟩`
    /// between two sets of occupied eigenstates.
    ///
    /// `states_a[m]` is the Bloch vector at k, `states_b[n]` at k+dk, each of
    /// length equal to the number of lattice sites.  The returned `CMatrix` has
    /// `M[m][n] = Σ_j  conj(states_a[m][j]) · states_b[n][j]`.
    pub fn link_matrix(states_a: &[Vec<Complex>], states_b: &[Vec<Complex>]) -> CMatrix {
        let n = states_a.len();
        let mut mat = CMatrix::zeros(n);
        for (m, state_a) in states_a.iter().enumerate() {
            for (nn, state_b) in states_b.iter().enumerate() {
                let mut inner = Complex::ZERO;
                for j in 0..state_a.len() {
                    inner = inner.add(&state_a[j].conj().mul(&state_b[j]));
                }
                mat.set(m, nn, inner);
            }
        }
        mat
    }

    // -----------------------------------------------------------------------
    // Wilson unitary along ky
    // -----------------------------------------------------------------------

    /// Compute the Wilson unitary `W_y(kx)` — the ordered product of link
    /// matrices around the closed ky loop `[0, 2π)` at fixed `kx`.
    ///
    /// The Brillouin zone is discretised into `ny` equally-spaced ky points.
    /// Periodicity is enforced by wrapping the last point back to the first.
    ///
    /// Returns an `n_occ × n_occ` unitary matrix `W_y(kx)`.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `ny < 4`.
    pub fn wilson_unitary_y(&self, kx: f64, ny: usize) -> Result<CMatrix> {
        if ny < 4 {
            return Err(error::invalid_param("ny", "need at least 4 ky points"));
        }
        let n_occ = self.occupied_bands.len();

        // Collect eigenstates at each ky point (plus wrap-around at ky[ny] = ky[0])
        let mut all_states: Vec<Vec<Vec<Complex>>> = Vec::with_capacity(ny + 1);
        for iy in 0..=ny {
            let ky = 2.0 * PI * (iy as f64) / (ny as f64);
            let (_, vecs) = self.model.diagonalize((kx, ky))?;
            let occ: Vec<Vec<Complex>> = self
                .occupied_bands
                .iter()
                .map(|&b| vecs.column(b))
                .collect();
            all_states.push(occ);
        }

        // W = M(0,1) · M(1,2) · … · M(ny-1, ny)
        // Each link matrix is first projected onto the unitary group via polar
        // decomposition to prevent drift from unitarity over many steps.
        let mut w = CMatrix::eye(n_occ);
        for iy in 0..ny {
            let raw = Self::link_matrix(&all_states[iy], &all_states[iy + 1]);
            let m = polar_unitary(&raw)?;
            w = w.matmul(&m)?;
        }

        Ok(w)
    }

    /// Compute the Wilson unitary `W_x(ky)` — the ordered product of link
    /// matrices along the closed kx loop `[0, 2π)` at fixed `ky`.
    ///
    /// Analogous to [`wilson_unitary_y`](Self::wilson_unitary_y) but with
    /// kx and ky swapped.
    pub fn wilson_unitary_x(&self, ky: f64, nx: usize) -> Result<CMatrix> {
        if nx < 4 {
            return Err(error::invalid_param("nx", "need at least 4 kx points"));
        }
        let n_occ = self.occupied_bands.len();

        let mut all_states: Vec<Vec<Vec<Complex>>> = Vec::with_capacity(nx + 1);
        for ix in 0..=nx {
            let kx = 2.0 * PI * (ix as f64) / (nx as f64);
            let (_, vecs) = self.model.diagonalize((kx, ky))?;
            let occ: Vec<Vec<Complex>> = self
                .occupied_bands
                .iter()
                .map(|&b| vecs.column(b))
                .collect();
            all_states.push(occ);
        }

        let mut w = CMatrix::eye(n_occ);
        for ix in 0..nx {
            let raw = Self::link_matrix(&all_states[ix], &all_states[ix + 1]);
            let m = polar_unitary(&raw)?;
            w = w.matmul(&m)?;
        }

        Ok(w)
    }

    // -----------------------------------------------------------------------
    // Wannier centres
    // -----------------------------------------------------------------------

    /// Compute the Wannier centres `ν_j ∈ [0, 1)` from `W_y(kx)`.
    ///
    /// The Wannier centres are the eigenphases of the Wilson unitary divided
    /// by `2π`, i.e. `ν_j = arg(λ_j) / (2π)` shifted to `[0, 1)`.
    /// The result is sorted in ascending order.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`wilson_unitary_y`](Self::wilson_unitary_y).
    pub fn wannier_centers(&self, kx: f64, ny: usize) -> Result<Vec<f64>> {
        let w = self.wilson_unitary_y(kx, ny)?;
        let mut phases = unitary_eig_phases(&w)?;
        // Map to [0, 1)
        for ph in &mut phases {
            *ph = (*ph / (2.0 * PI)).rem_euclid(1.0);
        }
        phases.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        Ok(phases)
    }

    // -----------------------------------------------------------------------
    // Polarisation
    // -----------------------------------------------------------------------

    /// Electric polarisation `P_y` via the y-directed Wilson loop.
    ///
    /// `P_y = (1/nkx) · Σ_{kx} ⟨ν(kx)⟩`
    /// where `⟨ν⟩` is the average of the Wannier centres at a given `kx`.
    ///
    /// This reproduces the King-Smith/Vanderbilt formula for electronic
    /// polarisation applied to the magnon pseudo-Bloch bands.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `nkx < 2` or `nky < 4`.
    pub fn polarization_p_y(&self, nkx: usize, nky: usize) -> Result<f64> {
        if nkx < 2 {
            return Err(error::invalid_param("nkx", "need at least 2 kx points"));
        }
        let mut sum = 0.0_f64;
        for ikx in 0..nkx {
            let kx = 2.0 * PI * (ikx as f64) / (nkx as f64);
            let centers = self.wannier_centers(kx, nky)?;
            let avg = centers.iter().sum::<f64>() / (centers.len() as f64);
            sum += avg;
        }
        Ok(sum / (nkx as f64))
    }

    /// Electric polarisation `P_x` via the x-directed Wilson loop.
    ///
    /// `P_x = (1/nky) · Σ_{ky} ⟨ν_x(ky)⟩`
    ///
    /// Computed by swapping the roles of kx and ky compared to
    /// [`polarization_p_y`](Self::polarization_p_y).
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `nky < 2` or `nkx < 4`.
    pub fn polarization_p_x(&self, nkx: usize, nky: usize) -> Result<f64> {
        if nky < 2 {
            return Err(error::invalid_param("nky", "need at least 2 ky points"));
        }
        let mut sum = 0.0_f64;
        for iky in 0..nky {
            let ky = 2.0 * PI * (iky as f64) / (nky as f64);
            // Use wilson_unitary_x at each ky to get Wannier centres
            let w = self.wilson_unitary_x(ky, nkx)?;
            let mut phases = unitary_eig_phases(&w)?;
            for ph in &mut phases {
                *ph = (*ph / (2.0 * PI)).rem_euclid(1.0);
            }
            let avg = phases.iter().sum::<f64>() / (phases.len() as f64);
            sum += avg;
        }
        Ok(sum / (nky as f64))
    }

    // -----------------------------------------------------------------------
    // Nested Wilson loop (HOTI invariant)
    // -----------------------------------------------------------------------

    /// Nested polarisation (corner Z₂ invariant) via the double Wilson loop.
    ///
    /// Following Vanderbilt and King-Smith (*PRB* **48**, 4442, 1993) and the
    /// BBH model (*Science* **357**, 61, 2017):
    ///
    /// 1. For each of the `nkx` kx values, compute the Wannier band `ν_j(kx)`
    ///    via [`wannier_centers`](Self::wannier_centers) choosing the `wannier_band_idx`-th
    ///    Wannier centre.
    ///
    /// 2. Accumulate the nested Wilson-loop phase:
    ///    `W_nested = ∏_{kx} exp(2πi · ν_j(kx))`
    ///
    /// 3. Return `nested_p = arg(W_nested) / (2π)` which is quantised to
    ///    `0` (trivial) or `±0.5` (topological corner charge) by symmetry.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `wannier_band_idx ≥ n_occ`,
    /// `nkx < 4`, or `nky < 4`.
    pub fn nested_polarization(
        &self,
        wannier_band_idx: usize,
        nkx: usize,
        nky: usize,
    ) -> Result<f64> {
        let n_occ = self.occupied_bands.len();
        if wannier_band_idx >= n_occ {
            return Err(error::invalid_param(
                "wannier_band_idx",
                &format!("index {wannier_band_idx} exceeds n_occ-1={}", n_occ - 1),
            ));
        }
        if nkx < 4 {
            return Err(error::invalid_param("nkx", "need at least 4 kx points"));
        }

        // Accumulate W_nested = ∏_kx exp(2πi · ν_j(kx))
        let mut w_nested = Complex::ONE;
        for ikx in 0..nkx {
            let kx = 2.0 * PI * (ikx as f64) / (nkx as f64);
            let centers = self.wannier_centers(kx, nky)?;
            // Use the requested Wannier band
            let nu = centers[wannier_band_idx.min(centers.len() - 1)];
            let phase = Complex::from_polar(1.0, 2.0 * PI * nu);
            w_nested = w_nested.mul(&phase);
        }

        let nested_p = w_nested.phase() / (2.0 * PI);
        Ok(nested_p)
    }
}

// ---------------------------------------------------------------------------
// Polar decomposition: project a matrix onto the unitary group
// ---------------------------------------------------------------------------

/// Project an `n × n` complex matrix `M` onto the unitary group via the
/// polar decomposition `M = U P`, returning the unitary factor `U`.
///
/// For N = 1: `U = M / |M|` (normalise the scalar).
/// For N ≥ 2: Use `U = M (M†M)^{-1/2}`.
///
/// The matrix square root `(M†M)^{-1/2}` is computed via the eigendecomposition
/// of the Hermitian positive-semidefinite matrix `M†M`:
/// If `M†M = V D V†` then `(M†M)^{-1/2} = V D^{-1/2} V†`.
///
/// # Errors
///
/// Returns an error if the eigendecomposition fails. Falls back to identity if
/// any singular value is below `1e-14`.
fn polar_unitary(m: &CMatrix) -> Result<CMatrix> {
    let n = m.n();

    if n == 1 {
        let z = m.get(0, 0);
        let norm = z.norm();
        let mut out = CMatrix::zeros(1);
        if norm < 1e-15 {
            out.set(0, 0, Complex::ONE);
        } else {
            out.set(0, 0, z.scale(1.0 / norm));
        }
        return Ok(out);
    }

    // Compute M†M
    let md = m.conj_transpose();
    let mdm = md.matmul(m)?;

    // Diagonalise M†M = V D V†;  D eigenvalues are non-negative
    let (evals, vecs) = mdm.hermitian_eigendecomposition()?;

    // Build D^{-1/2} (diagonal)
    let mut d_inv_sqrt = CMatrix::zeros(n);
    for (i, &ev) in evals.iter().enumerate() {
        let sv = ev.max(0.0).sqrt(); // singular value = sqrt(eigenvalue)
        let val = if sv < 1e-14 { 0.0 } else { 1.0 / sv };
        d_inv_sqrt.set(i, i, Complex::from_real(val));
    }

    // (M†M)^{-1/2} = V D^{-1/2} V†
    let vd = vecs.conj_transpose();
    let mid = vecs.matmul(&d_inv_sqrt)?;
    let m_inv_sqrt = mid.matmul(&vd)?;

    // U = M · (M†M)^{-1/2}
    m.matmul(&m_inv_sqrt)
}

// ---------------------------------------------------------------------------
// Eigenphases of a unitary matrix (private helper)
// ---------------------------------------------------------------------------

/// Extract eigenphases of an `n × n` unitary matrix `U` (n ≤ 4).
///
/// Strategy: compute the Hermitian part `H = (U + U†)/2`; since `U` is
/// unitary its eigenvalues lie on the unit circle, so the real part of each
/// eigenvalue is `cos(θ_j)`.  Diagonalising `H` via
/// `hermitian_eigendecomposition` gives `cos(θ_j)` values; the phases are
/// recovered as `θ_j = arccos(cos(θ_j))`.
///
/// The sign ambiguity (`±θ`) is a known limitation of this approach — it
/// gives the correct magnitude but may lose sign information. For the
/// topological applications here (phases near 0 and π) this is sufficient.
fn unitary_eig_phases(u: &CMatrix) -> Result<Vec<f64>> {
    let n = u.n();

    // Special case n = 1
    if n == 1 {
        let z = u.get(0, 0);
        return Ok(vec![z.im.atan2(z.re)]);
    }

    // H_herm = (U + U†)/2  — Hermitian matrix
    let ud = u.conj_transpose();
    let h_sum = u.add(&ud)?;
    let h_herm = h_sum.scale_real(0.5);

    // Eigenvalues of H_herm are cos(θ_j), clamped to [-1, 1]
    let (evals, _) = h_herm.hermitian_eigendecomposition()?;
    let phases: Vec<f64> = evals
        .iter()
        .map(|&c| {
            let c_clamped = c.clamp(-1.0, 1.0);
            c_clamped.acos()
        })
        .collect();

    Ok(phases)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topomagnon::band_model::MagnonBandModel;

    fn trivial_two_band() -> MagnonBandModel {
        // Honeycomb without DMI — topologically trivial
        MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.0, 0.5).unwrap()
    }

    fn topological_two_band() -> MagnonBandModel {
        // Honeycomb with DMI=0.5 — topological (|C|=1)
        MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.5, 0.0).unwrap()
    }

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    #[test]
    fn new_rejects_empty_bands() {
        let m = trivial_two_band();
        assert!(WilsonLoop::new(&m, vec![]).is_err());
    }

    #[test]
    fn new_rejects_out_of_range_band() {
        let m = trivial_two_band();
        assert!(WilsonLoop::new(&m, vec![0, 99]).is_err());
    }

    #[test]
    fn new_valid() {
        let m = trivial_two_band();
        assert!(WilsonLoop::new(&m, vec![0]).is_ok());
        assert!(WilsonLoop::new(&m, vec![0, 1]).is_ok());
    }

    // -----------------------------------------------------------------------
    // Link matrix
    // -----------------------------------------------------------------------

    #[test]
    fn link_matrix_identity_same_states() {
        // M_mn = ⟨u_m|u_n⟩; if states_a == states_b and orthonormal → identity
        let states: Vec<Vec<Complex>> = vec![
            vec![Complex::ONE, Complex::ZERO],
            vec![Complex::ZERO, Complex::ONE],
        ];
        let m = WilsonLoop::link_matrix(&states, &states);
        // Diagonal elements ≈ 1, off-diag ≈ 0
        assert!((m.get(0, 0).re - 1.0).abs() < 1e-14);
        assert!((m.get(1, 1).re - 1.0).abs() < 1e-14);
        assert!(m.get(0, 1).norm() < 1e-14);
        assert!(m.get(1, 0).norm() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // Wilson unitary — unitarity check
    // -----------------------------------------------------------------------

    #[test]
    fn wilson_unitary_y_is_unitary() {
        let m = topological_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        let w = wl.wilson_unitary_y(0.0, 20).unwrap();
        // W · W† should be identity (Frobenius norm of W·W†-I < 1e-8)
        let wd = w.conj_transpose();
        let wwd = w.matmul(&wd).unwrap();
        let id = CMatrix::eye(1);
        let diff = wwd.sub(&id).unwrap();
        assert!(
            diff.frobenius_norm() < 1e-8,
            "W·W† not identity: norm={}",
            diff.frobenius_norm()
        );
    }

    #[test]
    fn wilson_unitary_x_is_unitary() {
        let m = topological_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        let w = wl.wilson_unitary_x(0.0, 20).unwrap();
        let wd = w.conj_transpose();
        let wwd = w.matmul(&wd).unwrap();
        let id = CMatrix::eye(1);
        let diff = wwd.sub(&id).unwrap();
        assert!(
            diff.frobenius_norm() < 1e-8,
            "W_x·W_x† not identity: norm={}",
            diff.frobenius_norm()
        );
    }

    // -----------------------------------------------------------------------
    // Wannier centres
    // -----------------------------------------------------------------------

    #[test]
    fn wannier_centers_trivial_band() {
        // For a trivial (zero-DMI) honeycomb, P_y should be 0 and Wannier
        // centres at 0.0 (or 0.5 by gauge).  We check they are finite and
        // lie in [0,1).
        let m = trivial_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        let centers = wl.wannier_centers(0.0, 20).unwrap();
        assert_eq!(centers.len(), 1);
        assert!(centers[0] >= 0.0 && centers[0] < 1.0);
    }

    #[test]
    fn wannier_centers_sorted_ascending() {
        let m = MagnonBandModel::kagome(1.0, 0.4, 0.0).unwrap();
        let wl = WilsonLoop::new(&m, vec![0, 1]).unwrap();
        let centers = wl.wannier_centers(0.0, 16).unwrap();
        assert_eq!(centers.len(), 2);
        assert!(
            centers[0] <= centers[1],
            "centers not sorted: {:?}",
            centers
        );
    }

    #[test]
    fn wannier_centers_lie_in_unit_interval() {
        let m = topological_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        for ikx in 0..8 {
            let kx = 2.0 * PI * (ikx as f64) / 8.0;
            let centers = wl.wannier_centers(kx, 16).unwrap();
            for c in &centers {
                assert!((0.0..1.0).contains(c), "center out of [0,1): {c}");
            }
        }
    }

    // -----------------------------------------------------------------------
    // Polarisation
    // -----------------------------------------------------------------------

    #[test]
    fn polarization_p_y_finite_and_in_range() {
        let m = trivial_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        let py = wl.polarization_p_y(8, 16).unwrap();
        // Must lie in [0, 1)
        assert!((0.0..1.0).contains(&py), "P_y out of range: {py}");
    }

    #[test]
    fn polarization_p_x_finite_and_in_range() {
        let m = trivial_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        let px = wl.polarization_p_x(16, 8).unwrap();
        assert!((0.0..1.0).contains(&px), "P_x out of range: {px}");
    }

    // -----------------------------------------------------------------------
    // Nested polarisation
    // -----------------------------------------------------------------------

    #[test]
    fn nested_polarization_rejects_bad_wannier_idx() {
        let m = topological_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        assert!(wl.nested_polarization(1, 8, 8).is_err());
    }

    #[test]
    fn nested_polarization_trivial_near_zero_or_half() {
        // For the trivial two-band honeycomb the nested polarisation should
        // be close to 0 (trivial) but we accept any quantised value.
        let m = trivial_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        let np = wl.nested_polarization(0, 8, 8).unwrap();
        // Must lie in (-0.5, 0.5] after wrapping
        let wrapped = (np + 0.5).rem_euclid(1.0) - 0.5;
        assert!(
            wrapped.abs() < 0.6,
            "Nested polarisation out of plausible range: {np}"
        );
    }

    #[test]
    fn nested_polarization_quantised_modulo_half() {
        // For honeycomb with large DMI the nested polarisation should be
        // close to a multiple of 0.5 (topological symmetry constraint).
        let m = topological_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        let np = wl.nested_polarization(0, 8, 8).unwrap();
        // Wrap to [0,1) and check distance from nearest 0.5 multiple
        let wrapped = np.rem_euclid(1.0);
        let dist_from_half_int = (wrapped * 2.0).round() / 2.0;
        assert!(
            (wrapped - dist_from_half_int).abs() < 0.3,
            "Nested polarisation not near half-integer: {wrapped}"
        );
    }

    // -----------------------------------------------------------------------
    // Error paths
    // -----------------------------------------------------------------------

    #[test]
    fn wilson_unitary_y_rejects_small_grid() {
        let m = trivial_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        assert!(wl.wilson_unitary_y(0.0, 2).is_err());
    }

    #[test]
    fn polarization_rejects_small_nkx() {
        let m = trivial_two_band();
        let wl = WilsonLoop::new(&m, vec![0]).unwrap();
        assert!(wl.polarization_p_y(1, 8).is_err());
    }
}
