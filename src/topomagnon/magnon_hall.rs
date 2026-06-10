//! Matsumoto-Murakami magnon thermal Hall conductivity.
//!
//! Implements the Matsumoto-Murakami formula for the off-diagonal thermal
//! conductivity κ_xy of magnons in a topological magnon insulator:
//!
//! ```text
//! κ_xy = −(k_B / ℏ) · (1/(2π)²) · Σ_{n} ∫_{BZ} c₂(n_B(ε_n(k))) · Ω_n(k) d²k
//! ```
//!
//! where:
//! - `n_B(ε) = 1/(exp(ε/k_B T) − 1)` is the Bose-Einstein distribution,
//! - `Ω_n(k)` is the Berry curvature of band n,
//! - `c₂(x) = (1+x)(ln((1+x)/x))² − (ln x)² − 2·Li₂(−x)` is the Matsumoto-Murakami
//!   transport Berry factor (function of Bose occupation number x = n_B).
//!
//! # References
//!
//! - R. Matsumoto, S. Murakami, *Phys. Rev. Lett.* **106**, 197202 (2011).
//! - R. Matsumoto, S. Murakami, *Phys. Rev. B* **84**, 184406 (2011).
//!
//! # Units
//!
//! Energy input to the model is assumed in meV. The conductivity is returned
//! in SI units W/(K·m) after applying the physical prefactor with ℏ and k_B.
//! For pure unit-system-agnostic usage, treat the result as proportional to
//! the Hall angle.

use std::f64::consts::PI;

use super::band_model::MagnonBandModel;
use super::berry_curvature::BerryCurvature;
use crate::constants::{HBAR, KB};
use crate::error::{self, Result};

// Conversion factor: 1 meV = 1.602176634e-22 J
const MEV_TO_J: f64 = 1.602_176_634e-22;

// ---------------------------------------------------------------------------
// MagnonHallConductivity
// ---------------------------------------------------------------------------

/// Computes the magnon thermal Hall conductivity κ_xy via the Matsumoto-Murakami formula.
pub struct MagnonHallConductivity<'a> {
    /// Reference to the magnon band model.
    pub model: &'a MagnonBandModel,
    /// Temperature \[K\].
    pub temperature: f64,
}

impl<'a> MagnonHallConductivity<'a> {
    /// Create a magnon Hall conductivity calculator.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `temperature ≤ 0`.
    pub fn new(model: &'a MagnonBandModel, temperature: f64) -> Result<Self> {
        if temperature <= 0.0 {
            return Err(error::invalid_param("temperature", "must be positive"));
        }
        Ok(Self { model, temperature })
    }

    // -----------------------------------------------------------------------
    // Matsumoto-Murakami c₂ function
    // -----------------------------------------------------------------------

    /// Compute the transport Berry factor c₂(x) where x = n_B is the Bose occupation.
    ///
    /// ```text
    /// c₂(x) = (1+x)·(ln((1+x)/x))² − (ln x)² − 2·Li₂(−x)
    /// ```
    ///
    /// Special cases:
    /// - High temperature (ε/k_BT < 1×10⁻⁴): c₂ → π²/3 (classical limit).
    /// - Low temperature (ε/k_BT > 50): c₂ → 0 (exponentially suppressed magnons).
    #[inline]
    pub fn c2_function(energy: f64, temperature: f64) -> f64 {
        if energy <= 0.0 || temperature <= 0.0 {
            return 0.0;
        }
        let x = energy * MEV_TO_J / (KB * temperature); // dimensionless ε/(k_BT)
        if x < 1e-4 {
            return PI * PI / 3.0; // classical limit
        }
        if x > 50.0 {
            return 0.0; // quantum/low-T limit
        }
        // n_B = 1/(e^x - 1)
        let n_b = 1.0 / (x.exp() - 1.0);
        // c₂(n_B) = (1+n_B)·(ln((1+n_B)/n_B))² − (ln n_B)² − 2·Li₂(−n_B)
        let ln_nb = n_b.ln();
        let ln_1p_nb = (1.0 + n_b).ln();
        let term1 = (1.0 + n_b) * (ln_1p_nb - ln_nb).powi(2);
        let term2 = -(ln_nb * ln_nb);
        let li2 = Self::dilogarithm(-n_b);
        term1 + term2 - 2.0 * li2
    }

    // -----------------------------------------------------------------------
    // Dilogarithm Li₂(z)
    // -----------------------------------------------------------------------

    /// Compute the dilogarithm Li₂(z) = −∫₀ᶻ ln(1−t)/t dt.
    ///
    /// For z ∈ (−∞, 0) — the range arising from z = −n_B < 0:
    /// - z ∈ (−1, 0): power series Li₂(z) = Σ z^n/n² (up to 100 terms or tol 1e-12).
    /// - z < −1: functional equation Li₂(z) = −π²/6 − ½(ln(−z))² − Li₂(1/z).
    pub fn dilogarithm(z: f64) -> f64 {
        if z == 0.0 {
            return 0.0;
        }
        if z > 0.0 && (z - 1.0).abs() < 1e-14 {
            return PI * PI / 6.0; // Li₂(1) = π²/6
        }
        if z < -1.0 {
            // Functional equation: Li₂(z) = -π²/6 - ½(ln(-z))² - Li₂(1/z)
            let lnmz = (-z).ln();
            return -(PI * PI / 6.0) - 0.5 * lnmz * lnmz - Self::dilogarithm(1.0 / z);
        }
        if (z - 1.0).abs() < 1e-10 {
            // Near z=1: Li₂(1) = π²/6
            return PI * PI / 6.0;
        }
        // Power series for |z| ≤ 1: Li₂(z) = Σ_{n=1}^∞ z^n / n²
        let mut sum = 0.0_f64;
        let mut z_pow = z;
        for n in 1_usize..=100 {
            let term = z_pow / ((n * n) as f64);
            sum += term;
            z_pow *= z;
            if term.abs() < 1e-12 * (1.0 + sum.abs()) {
                break;
            }
        }
        sum
    }

    // -----------------------------------------------------------------------
    // Main computation
    // -----------------------------------------------------------------------

    /// Compute the magnon thermal Hall conductivity κ_xy.
    ///
    /// Evaluates the Matsumoto-Murakami sum over all bands and k-points on an
    /// `nx × ny` uniform grid. Returns κ_xy in W/(K·m).
    ///
    /// The sign of κ_xy is determined by the Berry curvature distribution and
    /// hence by the sign of the DMI / external field.
    pub fn compute(&self, nx: usize, ny: usize) -> Result<f64> {
        if nx < 2 || ny < 2 {
            return Err(error::invalid_param("nx/ny", "need at least 2 grid points"));
        }

        let nb = self.model.n_bands();
        let bc = BerryCurvature::new(self.model);
        let t = self.temperature;

        let dkx = 2.0 * PI / (nx as f64);
        let dky = 2.0 * PI / (ny as f64);
        let cell_area = dkx * dky;

        let mut sum = 0.0_f64;

        for band in 0..nb {
            for ix in 0..nx {
                let kx = -PI + (ix as f64) * dkx;
                for iy in 0..ny {
                    let ky = -PI + (iy as f64) * dky;
                    let (evals, _) = self.model.diagonalize((kx, ky))?;
                    let eps = evals[band]; // energy in meV (or whatever units j_nn is in)
                    let c2 = Self::c2_function(eps, t);
                    let omega = bc.curvature_at((kx, ky), band)?;
                    sum += c2 * omega;
                }
            }
        }

        // κ_xy = −(k_B / ℏ) · (1/(2π)²) · Σ_{n,k} c₂(n_B(ε_n(k))) · Ω_n(k) · Δkx·Δky
        // The factor k_B/ℏ has units W/(K·m) when Ω is in units of a² (lattice area).
        // We include a·² through model.a_lattice²
        let a = self.model.a_lattice;
        let prefactor = -(KB / HBAR) / (2.0 * PI * 2.0 * PI) / (a * a);
        Ok(prefactor * sum * cell_area)
    }

    // -----------------------------------------------------------------------
    // Temperature sweep
    // -----------------------------------------------------------------------

    /// Compute κ_xy at multiple temperatures.
    ///
    /// Returns a vector of `(T, κ_xy)` pairs.
    pub fn temperature_sweep(
        &self,
        temperatures: &[f64],
        nx: usize,
        ny: usize,
    ) -> Result<Vec<(f64, f64)>> {
        let mut results = Vec::with_capacity(temperatures.len());
        for &t in temperatures {
            if t <= 0.0 {
                return Err(error::invalid_param(
                    "temperatures",
                    "all temperatures must be positive",
                ));
            }
            let mhc = Self::new(self.model, t)?;
            let kxy = mhc.compute(nx, ny)?;
            results.push((t, kxy));
        }
        Ok(results)
    }

    // -----------------------------------------------------------------------
    // Thermal Hall resistivity
    // -----------------------------------------------------------------------

    /// Compute the magnon thermal Hall resistivity ρ_xy = 1 / κ_xy.
    ///
    /// Returns `NumericalError` if κ_xy = 0 (trivial model with no Berry curvature).
    pub fn thermal_hall_resistivity(&self, nx: usize, ny: usize) -> Result<f64> {
        let kxy = self.compute(nx, ny)?;
        if kxy.abs() < 1e-50 {
            return Err(error::numerical_error(
                "thermal Hall conductivity is zero — model has no Berry curvature",
            ));
        }
        Ok(1.0 / kxy)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::topomagnon::band_model::MagnonBandModel;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn c2_classical_limit() {
        // High temperature (x = ε/(k_BT) → 0): c₂ → π²/3
        let c2 = MagnonHallConductivity::c2_function(1e-10, 1e6); // extremely high T
        assert!(
            approx(c2, PI * PI / 3.0, 1e-6),
            "Classical c₂ limit: expected {}, got {}",
            PI * PI / 3.0,
            c2
        );
    }

    #[test]
    fn c2_quantum_limit() {
        // Low temperature (x → ∞): c₂ → 0
        let c2 = MagnonHallConductivity::c2_function(1000.0, 0.001); // very low T
        assert!(c2.abs() < 1e-6, "Quantum c₂ limit: expected ~0, got {}", c2);
    }

    #[test]
    fn dilogarithm_at_zero() {
        let li2 = MagnonHallConductivity::dilogarithm(0.0);
        assert!(approx(li2, 0.0, 1e-15), "Li₂(0) should be 0, got {}", li2);
    }

    #[test]
    fn dilogarithm_half() {
        // Li₂(1/2) ≈ 0.5822... (known value: π²/12 − (ln 2)²/2)
        let li2 = MagnonHallConductivity::dilogarithm(0.5);
        let expected = PI * PI / 12.0 - 0.5 * (2.0_f64.ln()).powi(2);
        assert!(
            approx(li2, expected, 1e-10),
            "Li₂(1/2): expected {:.6}, got {:.6}",
            expected,
            li2
        );
    }

    #[test]
    fn compute_zero_with_zero_dmi() {
        // Square lattice: no Berry curvature → κ_xy = 0
        let m = MagnonBandModel::square_dmi(1.0, 0.0, 0.0).unwrap();
        let mhc = MagnonHallConductivity::new(&m, 100.0).unwrap();
        let kxy = mhc.compute(8, 8).unwrap();
        assert!(
            kxy.abs() < 1e-30,
            "κ_xy should be zero for trivial band, got {}",
            kxy
        );
    }

    #[test]
    fn compute_nonzero_with_dmi() {
        // Honeycomb with DMI: Berry curvature → non-zero κ_xy
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.4, 0.0).unwrap();
        let mhc = MagnonHallConductivity::new(&m, 50.0).unwrap();
        let kxy = mhc.compute(10, 10).unwrap();
        assert!(
            kxy.abs() > 1e-20,
            "κ_xy should be non-zero with DMI, got {}",
            kxy
        );
    }

    #[test]
    fn compute_sign_flips_with_dmi_sign() {
        let m_pos = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.4, 0.0).unwrap();
        let m_neg = MagnonBandModel::honeycomb_haldane(1.0, 0.0, -0.4, 0.0).unwrap();
        let kxy_pos = MagnonHallConductivity::new(&m_pos, 50.0)
            .unwrap()
            .compute(8, 8)
            .unwrap();
        let kxy_neg = MagnonHallConductivity::new(&m_neg, 50.0)
            .unwrap()
            .compute(8, 8)
            .unwrap();
        assert!(
            kxy_pos * kxy_neg < 0.0,
            "κ_xy sign should flip with DMI sign: {} vs {}",
            kxy_pos,
            kxy_neg
        );
    }

    #[test]
    fn temperature_sweep_correct_length() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.3, 0.0).unwrap();
        let temps = vec![10.0, 50.0, 100.0, 200.0];
        let mhc = MagnonHallConductivity::new(&m, 100.0).unwrap();
        let sweep = mhc.temperature_sweep(&temps, 8, 8).unwrap();
        assert_eq!(sweep.len(), 4);
        for (i, &t) in temps.iter().enumerate() {
            assert!((sweep[i].0 - t).abs() < 1e-10);
        }
    }

    #[test]
    fn invalid_temperature_rejected() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.1, 0.0).unwrap();
        assert!(MagnonHallConductivity::new(&m, -1.0).is_err());
        assert!(MagnonHallConductivity::new(&m, 0.0).is_err());
    }
}
