//! Berry curvature Ω_n(k) for magnon bands.
//!
//! Computes the Berry curvature of the n-th magnon band via the sum-over-states
//! (Kubo) formula, using finite differences to evaluate ∂_k H(k):
//!
//! ```text
//! Ω_n(k) = −2 Im Σ_{m≠n}  ⟨u_n(k)|∂_kx H|u_m(k)⟩ ⟨u_m(k)|∂_ky H|u_n(k)⟩
//!                           ─────────────────────────────────────────────────
//!                                        (ε_n(k) − ε_m(k))²
//! ```
//!
//! Degenerate states (|ε_n − ε_m| < 1 × 10⁻¹⁰) contribute zero to avoid
//! numerical divergences; this is the standard treatment for gapless points.
//!
//! # Integration
//!
//! The Brillouin-zone integral `(1/(2π)) ∫ Ω_n(k) d²k` approximates the
//! Chern number and should agree with the integer returned by [`ChernNumber`]
//! to within ±0.5 for grids coarser than 20×20.
//!
//! [`ChernNumber`]: super::chern_number::ChernNumber

use crate::error::{self, Result};
use crate::math::{CMatrix, Complex};
use std::f64::consts::PI;

use super::band_model::MagnonBandModel;

// ---------------------------------------------------------------------------
// BerryCurvature
// ---------------------------------------------------------------------------

/// Computes the Berry curvature Ω_n(k) of the n-th magnon band.
pub struct BerryCurvature<'a> {
    /// Reference to the magnon band model.
    pub model: &'a MagnonBandModel,
    /// Finite-difference step for ∂_k H.
    pub dk: f64,
}

impl<'a> BerryCurvature<'a> {
    /// Create with default finite-difference step dk = 1 × 10⁻⁴.
    pub fn new(model: &'a MagnonBandModel) -> Self {
        Self { model, dk: 1e-4 }
    }

    /// Override the finite-difference step (builder pattern).
    pub fn with_step(mut self, dk: f64) -> Self {
        self.dk = dk;
        self
    }

    // -----------------------------------------------------------------------
    // Core computation
    // -----------------------------------------------------------------------

    /// Compute Berry curvature Ω_n(k) at a single k-point via the Kubo formula.
    ///
    /// Uses central finite differences:
    /// `∂_kx H ≈ (H(k + dk·x̂) − H(k − dk·x̂)) / (2·dk)`
    ///
    /// Returns `Err` if `band_idx` is out of range.
    pub fn curvature_at(&self, k: (f64, f64), band_idx: usize) -> Result<f64> {
        let nb = self.model.n_bands();
        if band_idx >= nb {
            return Err(error::invalid_param(
                "band_idx",
                "index exceeds number of bands",
            ));
        }

        let (kx, ky) = k;
        let dk = self.dk;

        // Evaluate H at five points: centre, ±dk in x, ±dk in y
        let (evals, vecs) = self.model.diagonalize(k)?;

        let h_px = self.model.hamiltonian_at((kx + dk, ky))?;
        let h_mx = self.model.hamiltonian_at((kx - dk, ky))?;
        let h_py = self.model.hamiltonian_at((kx, ky + dk))?;
        let h_my = self.model.hamiltonian_at((kx, ky - dk))?;

        // ∂_kx H and ∂_ky H (central differences)
        let dh_x = h_px.sub(&h_mx)?.scale_real(1.0 / (2.0 * dk));
        let dh_y = h_py.sub(&h_my)?.scale_real(1.0 / (2.0 * dk));

        // Kubo sum over states m ≠ n
        let eps_n = evals[band_idx];
        let mut omega = 0.0_f64;

        for (m, &eps_m) in evals.iter().enumerate() {
            if m == band_idx {
                continue;
            }
            let denom = eps_n - eps_m;
            if denom.abs() < 1e-10 {
                continue; // degenerate — skip
            }

            // Matrix elements <u_n|∂_kx H|u_m> and <u_m|∂_ky H|u_n>
            let v_n = vecs.column(band_idx);
            let v_m = vecs.column(m);

            let mx = matrix_element(&v_n, &dh_x, &v_m, nb);
            let my = matrix_element(&v_m, &dh_y, &v_n, nb);

            // Contribution: -2 Im( mx * my ) / denom²
            let prod = mx.mul(&my);
            omega += -2.0 * prod.im / (denom * denom);
        }

        Ok(omega)
    }

    // -----------------------------------------------------------------------
    // Grid evaluation
    // -----------------------------------------------------------------------

    /// Compute Ω_n on a uniform `kx_pts × ky_pts` grid over [−π, π]².
    ///
    /// Returns `grid[ikx][iky]` with the Berry curvature at each k-point.
    pub fn compute_grid(
        &self,
        kx_pts: usize,
        ky_pts: usize,
        band_idx: usize,
    ) -> Result<Vec<Vec<f64>>> {
        if band_idx >= self.model.n_bands() {
            return Err(error::invalid_param(
                "band_idx",
                "index exceeds number of bands",
            ));
        }
        if kx_pts < 2 || ky_pts < 2 {
            return Err(error::invalid_param("kx_pts/ky_pts", "must be at least 2"));
        }

        let mut grid = Vec::with_capacity(kx_pts);
        for ix in 0..kx_pts {
            let kx = -PI + 2.0 * PI * (ix as f64) / (kx_pts as f64);
            let mut row = Vec::with_capacity(ky_pts);
            for iy in 0..ky_pts {
                let ky = -PI + 2.0 * PI * (iy as f64) / (ky_pts as f64);
                let omega = self.curvature_at((kx, ky), band_idx)?;
                row.push(omega);
            }
            grid.push(row);
        }
        Ok(grid)
    }

    // -----------------------------------------------------------------------
    // BZ integration
    // -----------------------------------------------------------------------

    /// Integrate the Berry curvature over the Brillouin zone.
    ///
    /// Returns `(1/(2π)) ∫_{BZ} Ω_n(k) d²k` using the trapezoidal rule on the
    /// uniform grid. For an isolated band this should be close to the Chern integer.
    pub fn integrate_brillouin(
        &self,
        kx_pts: usize,
        ky_pts: usize,
        band_idx: usize,
    ) -> Result<f64> {
        let grid = self.compute_grid(kx_pts, ky_pts, band_idx)?;

        let dkx = 2.0 * PI / (kx_pts as f64);
        let dky = 2.0 * PI / (ky_pts as f64);
        let cell_area = dkx * dky;

        let mut sum = 0.0_f64;
        for row in &grid {
            for &omega in row {
                sum += omega;
            }
        }

        // (1/(2π)) * sum * dkx * dky
        Ok(sum * cell_area / (2.0 * PI))
    }

    /// Backward-compatible alias: compute Berry curvature at (kx, ky) for `band`.
    pub fn at(&self, band: usize, kx: f64, ky: f64) -> Result<f64> {
        self.curvature_at((kx, ky), band)
    }

    /// Backward-compatible alias: integrate over BZ on an n×n grid.
    pub fn chern_number_approx(&self, band: usize, n_grid: usize) -> Result<f64> {
        self.integrate_brillouin(n_grid, n_grid, band)
    }
}

// ---------------------------------------------------------------------------
// Helper: matrix element ⟨v_n|M|v_m⟩
// ---------------------------------------------------------------------------

/// Compute ⟨v_n | M | v_m⟩ where v_n, v_m are column vectors of length `n`.
fn matrix_element(v_bra: &[Complex], m: &CMatrix, v_ket: &[Complex], n: usize) -> Complex {
    // Mv_ket = M * v_ket  (n × 1 result)
    let mut mv = vec![Complex::ZERO; n];
    for (i, mv_i) in mv.iter_mut().enumerate() {
        *mv_i = v_ket
            .iter()
            .enumerate()
            .fold(Complex::ZERO, |acc, (j, &vkj)| {
                acc.add(&m.get(i, j).mul(&vkj))
            });
    }
    // ⟨v_bra|mv⟩ = Σ_i conj(v_bra[i]) * mv[i]
    v_bra
        .iter()
        .zip(mv.iter())
        .fold(Complex::ZERO, |acc, (&bra_i, &mv_i)| {
            acc.add(&bra_i.conj().mul(&mv_i))
        })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topomagnon::band_model::MagnonBandModel;

    #[test]
    fn curvature_at_returns_finite() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.2, 0.1).unwrap();
        let bc = BerryCurvature::new(&m);
        let omega = bc.curvature_at((0.3, 0.4), 0).unwrap();
        assert!(omega.is_finite(), "Berry curvature must be finite");
    }

    #[test]
    fn curvature_antisymmetry() {
        // For time-reversal invariant model (h_ext=0, dmi=0): Ω_n(k) = -Ω_n(-k)
        // We test this approximately at a generic k-point.
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.0, 0.0).unwrap();
        let bc = BerryCurvature::new(&m);
        let k = (0.5, 0.3);
        let omega_k = bc.curvature_at(k, 0).unwrap();
        let omega_mk = bc.curvature_at((-k.0, -k.1), 0).unwrap();
        // Should satisfy Ω(k) ≈ -Ω(-k) for TR-symmetric model
        assert!(
            (omega_k + omega_mk).abs() < 0.1,
            "TR antisymmetry violated: Ω(k)={}, Ω(-k)={}",
            omega_k,
            omega_mk
        );
    }

    #[test]
    fn curvature_zero_for_trivial_band() {
        // Square lattice with no DMI: single-band, curvature is identically zero
        let m = MagnonBandModel::square_dmi(1.0, 0.0, 0.0).unwrap();
        let bc = BerryCurvature::new(&m);
        let omega = bc.curvature_at((0.5, 0.7), 0).unwrap();
        // 1×1 matrix: only one band, no off-diagonal terms → omega must be exactly 0
        assert!(
            omega.abs() < 1e-10,
            "1-band curvature must be zero, got {}",
            omega
        );
    }

    #[test]
    fn compute_grid_size_correct() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.1, 0.0).unwrap();
        let bc = BerryCurvature::new(&m);
        let grid = bc.compute_grid(10, 12, 0).unwrap();
        assert_eq!(grid.len(), 10);
        assert_eq!(grid[0].len(), 12);
    }

    #[test]
    fn integrate_matches_chern_within_half() {
        // For honeycomb with large DMI the integral should be close to ±1
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.5, 0.0).unwrap();
        let bc = BerryCurvature::new(&m);
        let integral = bc.integrate_brillouin(25, 25, 0).unwrap();
        // Should be within 0.4 of an integer (Berry curvature integral → Chern number)
        let nearest_int = integral.round();
        assert!(
            (integral - nearest_int).abs() < 0.4,
            "BZ integral {} not close to integer",
            integral
        );
    }

    #[test]
    fn curvature_opposite_sign_bands() {
        // For 2-band model: Ω_0 + Ω_1 should integrate to ~0 (sum rule)
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.3, 0.0).unwrap();
        let bc = BerryCurvature::new(&m);
        let i0 = bc.integrate_brillouin(15, 15, 0).unwrap();
        let i1 = bc.integrate_brillouin(15, 15, 1).unwrap();
        assert!(
            (i0 + i1).abs() < 0.3,
            "Band sum rule violated: Ω_0+Ω_1 integral = {}",
            i0 + i1
        );
    }

    #[test]
    fn invalid_band_idx_errors() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.1, 0.0).unwrap();
        let bc = BerryCurvature::new(&m);
        assert!(bc.curvature_at((0.0, 0.0), 5).is_err());
        assert!(bc.compute_grid(5, 5, 5).is_err());
        assert!(bc.integrate_brillouin(5, 5, 5).is_err());
    }

    #[test]
    fn with_step_sets_dk() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.1, 0.0).unwrap();
        let bc = BerryCurvature::new(&m).with_step(1e-5);
        assert!((bc.dk - 1e-5).abs() < 1e-20);
    }
}
