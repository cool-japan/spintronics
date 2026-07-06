//! Finite-difference Kubo-formula Berry curvature: a numerical cross-check.
//!
//! [`KuboBerry`] independently evaluates the standard interband Kubo sum
//!
//! ```text
//! Omega_n(k) = -2 Im sum_{m != n} <n|dH/dkx|m><m|dH/dky|n> / (E_n - E_m)^2
//! ```
//!
//! against any [`BlochHamiltonian`], finite-differencing `dH/dk` directly from
//! [`BlochHamiltonian::hamiltonian_at`] (no analytic derivative required). This
//! mirrors `crate::topomagnon::berry_curvature::BerryCurvature` (same formula,
//! same finite-difference construction) but is built purely against the new
//! [`BlochHamiltonian`] trait, so it works for the altermagnet band model
//! without touching -- or depending on -- any `src/topomagnon/*` code.
//!
//! `KuboBerry` exists solely as an independent numerical cross-check of
//! [`AltermagnetBandModel::berry_curvature`](crate::altermagnet::band_model::AltermagnetBandModel::berry_curvature)'s
//! closed-form result; it is not a replacement for it, and is deliberately kept
//! lean (finite differences only, no caching or grid machinery).

use crate::error::{self, Result};
use crate::math::{CMatrix, Complex};

use super::bloch_hamiltonian::BlochHamiltonian;

/// Numerical Kubo-formula Berry-curvature calculator for any [`BlochHamiltonian`].
pub struct KuboBerry<'a> {
    /// The Bloch Hamiltonian being probed.
    pub hamiltonian: &'a dyn BlochHamiltonian,
    /// Finite-difference step for `dH/dk`.
    pub dk: f64,
}

impl<'a> KuboBerry<'a> {
    /// Create with the default finite-difference step `dk = 1e-4`.
    pub fn new(hamiltonian: &'a dyn BlochHamiltonian) -> Self {
        Self {
            hamiltonian,
            dk: 1e-4,
        }
    }

    /// Override the finite-difference step (builder pattern).
    pub fn with_step(mut self, dk: f64) -> Self {
        self.dk = dk;
        self
    }

    /// Berry curvature `Omega_n(kx, ky)` of band `band_idx` via the Kubo formula.
    ///
    /// Evaluates `hamiltonian_at` on a five-point finite-difference stencil
    /// (centre plus `+/- dk` along `kx` and `ky`), diagonalizes the centre point
    /// via [`BlochHamiltonian::diagonalize_at`], and sums the interband Kubo
    /// formula over all bands `m != band_idx`. Near-degenerate pairs
    /// (`|E_n - E_m| < 1e-10`) are skipped to avoid numerical divergence, the
    /// same convergence guard used by
    /// `crate::topomagnon::berry_curvature::BerryCurvature`.
    ///
    /// # Errors
    ///
    /// Returns an error if `band_idx >= n_bands()`, or if the underlying
    /// [`BlochHamiltonian`] calls fail.
    pub fn curvature_at(&self, kx: f64, ky: f64, band_idx: usize) -> Result<f64> {
        let nb = self.hamiltonian.n_bands();
        if band_idx >= nb {
            return Err(error::invalid_param(
                "band_idx",
                "index exceeds number of bands",
            ));
        }
        let dk = self.dk;

        let (evals, vecs) = self.hamiltonian.diagonalize_at(kx, ky)?;

        let h_px = self.hamiltonian.hamiltonian_at(kx + dk, ky)?;
        let h_mx = self.hamiltonian.hamiltonian_at(kx - dk, ky)?;
        let h_py = self.hamiltonian.hamiltonian_at(kx, ky + dk)?;
        let h_my = self.hamiltonian.hamiltonian_at(kx, ky - dk)?;

        // d/dkx H and d/dky H via central finite differences.
        let dh_x = h_px.sub(&h_mx)?.scale_real(1.0 / (2.0 * dk));
        let dh_y = h_py.sub(&h_my)?.scale_real(1.0 / (2.0 * dk));

        let eps_n = evals[band_idx];
        let mut omega = 0.0_f64;

        for (m, &eps_m) in evals.iter().enumerate() {
            if m == band_idx {
                continue;
            }
            let denom = eps_n - eps_m;
            if denom.abs() < 1e-10 {
                continue; // near-degenerate -- skip to avoid a numerical blow-up
            }

            let v_n = vecs.column(band_idx);
            let v_m = vecs.column(m);

            let mx = matrix_element(&v_n, &dh_x, &v_m, nb);
            let my = matrix_element(&v_m, &dh_y, &v_n, nb);

            // Contribution: -2 Im(mx * my) / denom^2.
            let prod = mx.mul(&my);
            omega += -2.0 * prod.im / (denom * denom);
        }

        Ok(omega)
    }
}

/// Compute `<v_bra|M|v_ket>` for column vectors of dimension `n`.
fn matrix_element(v_bra: &[Complex], m: &CMatrix, v_ket: &[Complex], n: usize) -> Complex {
    let mut mv = vec![Complex::ZERO; n];
    for (i, mv_i) in mv.iter_mut().enumerate() {
        *mv_i = v_ket
            .iter()
            .enumerate()
            .fold(Complex::ZERO, |acc, (j, &vkj)| {
                acc.add(&m.get(i, j).mul(&vkj))
            });
    }
    v_bra
        .iter()
        .zip(mv.iter())
        .fold(Complex::ZERO, |acc, (&bra_i, &mv_i)| {
            acc.add(&bra_i.conj().mul(&mv_i))
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::altermagnet::band_model::{AltermagnetBandModel, Band, Spin};
    use crate::altermagnet::bloch_hamiltonian::AltermagnetSpinHamiltonian;
    use crate::altermagnet::materials::AltermagneticSymmetry;

    fn soc_model() -> AltermagnetBandModel {
        AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.4,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid")
    }

    fn no_soc_model() -> AltermagnetBandModel {
        AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.0,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid")
    }

    #[test]
    fn curvature_at_returns_finite() {
        let model = soc_model();
        let h = AltermagnetSpinHamiltonian::new(&model, Spin::Up);
        let kubo = KuboBerry::new(&h);
        let omega = kubo.curvature_at(0.4, 0.15, 0).expect("valid");
        assert!(omega.is_finite());
    }

    #[test]
    fn invalid_band_idx_errors() {
        let model = soc_model();
        let h = AltermagnetSpinHamiltonian::new(&model, Spin::Up);
        let kubo = KuboBerry::new(&h);
        assert!(kubo.curvature_at(0.0, 0.0, 2).is_err());
        assert!(kubo.curvature_at(0.0, 0.0, 5).is_err());
    }

    #[test]
    fn with_step_sets_dk() {
        let model = soc_model();
        let h = AltermagnetSpinHamiltonian::new(&model, Spin::Up);
        let kubo = KuboBerry::new(&h).with_step(1e-5);
        assert!((kubo.dk - 1e-5).abs() < 1e-20);
    }

    /// Cross-check: `KuboBerry` (numerical) agrees with the closed-form
    /// `berry_curvature` at several k-points, for SOC on (nonzero curvature
    /// expected).
    #[test]
    fn kubo_matches_closed_form_with_soc() {
        let model = soc_model();
        for spin in [Spin::Up, Spin::Down] {
            let h = AltermagnetSpinHamiltonian::new(&model, spin);
            let kubo = KuboBerry::new(&h);
            for &(kx, ky) in &[(0.4, 0.15), (0.9, -0.6), (-0.3, 0.5)] {
                // band index 0 = Lower (ascending eigenvalues), 1 = Upper.
                let closed_lower = model.berry_curvature(kx, ky, spin, Band::Lower);
                let closed_upper = model.berry_curvature(kx, ky, spin, Band::Upper);
                let numeric_lower = kubo.curvature_at(kx, ky, 0).expect("valid");
                let numeric_upper = kubo.curvature_at(kx, ky, 1).expect("valid");
                let scale = closed_lower.abs().max(closed_upper.abs()).max(1e-6);
                assert!(
                    (closed_lower - numeric_lower).abs() / scale < 1e-3,
                    "lower band mismatch at k=({kx},{ky}), spin={spin:?}: closed={closed_lower}, kubo={numeric_lower}"
                );
                assert!(
                    (closed_upper - numeric_upper).abs() / scale < 1e-3,
                    "upper band mismatch at k=({kx},{ky}), spin={spin:?}: closed={closed_upper}, kubo={numeric_upper}"
                );
            }
        }
    }

    /// Cross-check: both the closed form and `KuboBerry` agree that curvature
    /// vanishes identically once SOC is switched off.
    #[test]
    fn kubo_and_closed_form_both_vanish_without_soc() {
        let model = no_soc_model();
        let h = AltermagnetSpinHamiltonian::new(&model, Spin::Up);
        let kubo = KuboBerry::new(&h);
        for &(kx, ky) in &[(0.4, 0.15), (0.9, -0.6)] {
            let closed = model.berry_curvature(kx, ky, Spin::Up, Band::Lower);
            let numeric = kubo.curvature_at(kx, ky, 0).expect("valid");
            assert!(closed.abs() < 1e-9, "closed form should vanish: {closed}");
            assert!(numeric.abs() < 1e-6, "Kubo cross-check should vanish: {numeric}");
        }
    }
}
