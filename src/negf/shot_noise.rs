//! Shot noise and current fluctuations in mesoscopic conductors.
//!
//! Implements the Landauer–Büttiker shot noise theory for a coherent 1D conductor.
//!
//! # Shot Noise Formula (Blanter–Büttiker)
//!
//! The zero-frequency noise power spectral density is:
//!
//! ```text
//! S(0) = (2e²/h) ∫ { T(E)[f_L(1−f_L) + f_R(1−f_R)]
//!                   + T(E)(1−T(E))(f_L − f_R)² } dE
//! ```
//!
//! The first term is thermal (equilibrium) noise; the second is the shot
//! (partition) noise arising from quantum transmission fluctuations.
//!
//! The Fano factor F = S_shot / (2eI) interpolates between F=1 (full Poissonian
//! shot noise for a tunnel junction, T→0) and F=0 (noiseless ballistic limit, T=1).
//!
//! # References
//!
//! - Ya. M. Blanter, M. Büttiker, Phys. Rep. **336**, 1 (2000)
//! - M. Büttiker, Phys. Rev. Lett. **65**, 2901 (1990)

use super::green_function::TransportCalculator;
use crate::constants::{E_CHARGE, H_PLANCK, KB};
use crate::error::Result;

// ============================================================================
// ShotNoise
// ============================================================================

/// Shot noise calculator for a Landauer–Büttiker coherent conductor.
///
/// Wraps a [`TransportCalculator`] and computes zero-frequency noise power,
/// thermal noise, partition noise, and the Fano factor.
#[derive(Debug, Clone)]
pub struct ShotNoise {
    /// The underlying transport calculator with energy grid.
    pub transport: TransportCalculator,
}

impl ShotNoise {
    /// Create a new `ShotNoise` calculator.
    pub fn new(transport: TransportCalculator) -> Self {
        Self { transport }
    }

    /// Helper: SI prefactor 2e²/h \[S\] for noise calculations.
    #[inline]
    fn noise_prefactor() -> f64 {
        2.0 * E_CHARGE * E_CHARGE / H_PLANCK
    }

    /// Total zero-frequency noise power S(0) \[A²/Hz\].
    ///
    /// Combines thermal fluctuations and shot (partition) noise:
    ///
    /// ```text
    /// S = (2e²/h) ∫ { T(E)[f_L(1−f_L) + f_R(1−f_R)] + T(E)(1−T(E))(f_L−f_R)² } dE
    /// ```
    ///
    /// # Errors
    ///
    /// Propagates errors from the transmission evaluation.
    pub fn noise_zero_freq(&self, v_bias: f64, temperature: f64) -> Result<f64> {
        let mu_l = v_bias / 2.0;
        let mu_r = -v_bias / 2.0;
        let tc = &self.transport;
        let de = (tc.e_max - tc.e_min) / (tc.n_energy - 1) as f64;
        let prefactor = Self::noise_prefactor();

        let mut values = Vec::with_capacity(tc.n_energy);
        for k in 0..tc.n_energy {
            let e = tc.e_min + k as f64 * de;
            let t = tc.transmission(e)?;
            let fl = TransportCalculator::fermi_dirac(e, mu_l, temperature);
            let fr = TransportCalculator::fermi_dirac(e, mu_r, temperature);
            // Thermal term
            let thermal_term = t * (fl * (1.0 - fl) + fr * (1.0 - fr));
            // Partition (shot) term
            let df = fl - fr;
            let shot_term = t * (1.0 - t) * df * df;
            values.push(thermal_term + shot_term);
        }

        let integral = trapezoid_integrate(&values, de);
        Ok(prefactor * integral)
    }

    /// Shot (partition) noise only S_shot \[A²/Hz\].
    ///
    /// ```text
    /// S_shot = (2e²/h) ∫ T(E)(1−T(E))(f_L − f_R)² dE
    /// ```
    ///
    /// # Errors
    ///
    /// Propagates errors from the transmission evaluation.
    pub fn shot_noise_only(&self, v_bias: f64, temperature: f64) -> Result<f64> {
        let mu_l = v_bias / 2.0;
        let mu_r = -v_bias / 2.0;
        let tc = &self.transport;
        let de = (tc.e_max - tc.e_min) / (tc.n_energy - 1) as f64;
        let prefactor = Self::noise_prefactor();

        let mut values = Vec::with_capacity(tc.n_energy);
        for k in 0..tc.n_energy {
            let e = tc.e_min + k as f64 * de;
            let t = tc.transmission(e)?;
            let fl = TransportCalculator::fermi_dirac(e, mu_l, temperature);
            let fr = TransportCalculator::fermi_dirac(e, mu_r, temperature);
            let df = fl - fr;
            values.push(t * (1.0 - t) * df * df);
        }

        let integral = trapezoid_integrate(&values, de);
        Ok(prefactor * integral)
    }

    /// Johnson–Nyquist (thermal) noise at zero bias S_th = 4 k_B T G(V=0) \[A²/Hz\].
    ///
    /// ```text
    /// S_th = 4 k_B T · (e²/h) · T(E_F)
    /// ```
    ///
    /// # Errors
    ///
    /// Propagates errors from the transmission evaluation at the Fermi level.
    pub fn thermal_noise(&self, temperature: f64) -> Result<f64> {
        let e_f = self.transport.gf.sigma_l.e_fermi;
        let t = self.transport.transmission(e_f)?;
        // G = (e²/h) * T, so S_th = 4 kB T G
        let conductance = (E_CHARGE * E_CHARGE / H_PLANCK) * t;
        // KB [J/K] * T [K] = thermal energy in Joules
        Ok(4.0 * KB * temperature * conductance)
    }

    /// Fano factor F = S_shot / (2eI) (dimensionless).
    ///
    /// Interpolates between:
    /// - F = 1 for a tunnel junction (T → 0)
    /// - F = 0 for a perfect ballistic channel (T = 1)
    ///
    /// Returns `1.0` if current |I| < 1e-30 A (vanishing current limit).
    ///
    /// # Errors
    ///
    /// Propagates errors from [`shot_noise_only`](Self::shot_noise_only) or current calculation.
    pub fn fano_factor(&self, v_bias: f64, temperature: f64) -> Result<f64> {
        let s_shot = self.shot_noise_only(v_bias, temperature)?;
        let current = self.transport.current(v_bias, temperature)?;
        let abs_current = current.abs();
        if abs_current < 1e-30 {
            // Vanishing current: return Poissonian limit
            return Ok(1.0);
        }
        Ok(s_shot / (2.0 * E_CHARGE * abs_current))
    }

    /// Partition noise at T=0 \[A²/Hz\].
    ///
    /// At zero temperature, only states between μ_R and μ_L contribute:
    ///
    /// ```text
    /// S_part = (2e²/h) ∫_{μ_R}^{μ_L} T(E)(1−T(E)) dE
    /// ```
    ///
    /// For negative bias (v_bias < 0), the limits are swapped and the result
    /// is still non-negative.
    ///
    /// # Errors
    ///
    /// Propagates errors from the transmission evaluation.
    pub fn partition_noise(&self, v_bias: f64) -> Result<f64> {
        let mu_l = v_bias / 2.0; // chem. pot. of left lead
        let mu_r = -v_bias / 2.0; // chem. pot. of right lead
                                  // Integration window: [min(mu_l, mu_r), max(mu_l, mu_r)]
        let e_lo = mu_l.min(mu_r);
        let e_hi = mu_l.max(mu_r);

        if (e_hi - e_lo).abs() < 1e-15 {
            return Ok(0.0);
        }

        let tc = &self.transport;
        let n_pts = tc.n_energy;
        let de = (e_hi - e_lo) / (n_pts - 1) as f64;
        let prefactor = Self::noise_prefactor();

        let mut values = Vec::with_capacity(n_pts);
        for k in 0..n_pts {
            let e = e_lo + k as f64 * de;
            let t = tc.transmission(e)?;
            values.push(t * (1.0 - t));
        }

        let integral = trapezoid_integrate(&values, de);
        Ok(prefactor * integral)
    }
}

// ============================================================================
// Helpers
// ============================================================================

/// Trapezoid rule integration.
fn trapezoid_integrate(values: &[f64], de: f64) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let n = values.len();
    let mut sum = 0.0;
    for k in 0..(n - 1) {
        sum += 0.5 * (values[k] + values[k + 1]) * de;
    }
    sum
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::negf::green_function::{
        GreenFunction, Hamiltonian1D, LeadSelfEnergy, TransportCalculator,
    };

    /// Helper: build a ShotNoise calculator for an n-site chain with given coupling Γ.
    fn make_noise(n: usize, gamma: f64) -> ShotNoise {
        let h = Hamiltonian1D::from_uniform(n, 0.0, 1.0).expect("valid");
        let sl = LeadSelfEnergy::new(gamma, 0.0).expect("valid");
        let sr = LeadSelfEnergy::new(gamma, 0.0).expect("valid");
        let gf = GreenFunction::new(h, sl, sr, 1e-3).expect("valid");
        let tc = TransportCalculator::new(gf, -3.0, 3.0, 200).expect("valid");
        ShotNoise::new(tc)
    }

    // ------------------------------------------------------------------ fano limits

    #[test]
    fn test_fano_open_channel() {
        // Large coupling → T ≈ 1 at E=0 → Fano ≈ 0
        let sn = make_noise(2, 5.0);
        let fano = sn.fano_factor(0.5, 300.0).expect("ok");
        // Fano should be small for open channel
        assert!(
            fano < 0.5,
            "Fano = {} (should be small for open channel)",
            fano
        );
    }

    #[test]
    fn test_fano_tunneling() {
        // Very small coupling → T ≈ 0 → Fano ≈ 1
        let sn = make_noise(5, 0.01);
        let fano = sn.fano_factor(0.5, 300.0).expect("ok");
        // Fano should approach 1 for tunnel junction
        assert!(
            fano > 0.5,
            "Fano = {} (should approach 1 for tunnel junction)",
            fano
        );
    }

    // ------------------------------------------------------------------ positivity

    #[test]
    fn test_thermal_noise_positive() {
        let sn = make_noise(3, 0.3);
        let s_th = sn.thermal_noise(300.0).expect("ok");
        assert!(s_th >= 0.0, "thermal noise must be non-negative: {}", s_th);
    }

    #[test]
    fn test_noise_positive() {
        let sn = make_noise(3, 0.3);
        let s = sn.noise_zero_freq(0.5, 300.0).expect("ok");
        assert!(s >= 0.0, "total noise must be non-negative: {}", s);
    }

    #[test]
    fn test_total_noise_ge_shot_only() {
        // Total noise ≥ shot noise (thermal term is non-negative)
        let sn = make_noise(4, 0.4);
        let s_total = sn.noise_zero_freq(0.3, 300.0).expect("ok");
        let s_shot = sn.shot_noise_only(0.3, 300.0).expect("ok");
        assert!(
            s_total >= s_shot - 1e-25,
            "total noise ({}) must be ≥ shot noise ({})",
            s_total,
            s_shot
        );
    }

    // ------------------------------------------------------------------ scaling

    #[test]
    fn test_shot_noise_scales_with_v() {
        // Shot noise should increase with voltage
        let sn = make_noise(4, 0.3);
        let s_low = sn.shot_noise_only(0.1, 10.0).expect("ok");
        let s_high = sn.shot_noise_only(0.5, 10.0).expect("ok");
        assert!(
            s_high >= s_low,
            "shot noise should increase with voltage: s_low={}, s_high={}",
            s_low,
            s_high
        );
    }

    #[test]
    fn test_partition_noise_zero_at_v0() {
        // At V=0, integration window [0, 0] → zero
        let sn = make_noise(4, 0.3);
        let s = sn.partition_noise(0.0).expect("ok");
        assert!(
            s.abs() < 1e-20,
            "partition noise at V=0 should be zero: {}",
            s
        );
    }

    #[test]
    fn test_fano_in_zero_one() {
        // Fano factor should be in [0, 1] (or close, allowing small numerical errors)
        let sn = make_noise(4, 0.3);
        let fano = sn.fano_factor(0.4, 300.0).expect("ok");
        assert!(fano >= -0.01, "Fano < 0: {}", fano);
        assert!(fano <= 1.5, "Fano > 1: {}", fano);
    }

    #[test]
    fn test_noise_higher_at_higher_t() {
        // At zero bias, higher T → more thermal noise
        let sn = make_noise(3, 0.3);
        let s_cold = sn.noise_zero_freq(0.0, 50.0).expect("ok");
        let s_hot = sn.noise_zero_freq(0.0, 600.0).expect("ok");
        assert!(
            s_hot >= s_cold,
            "noise should increase with temperature: cold={}, hot={}",
            s_cold,
            s_hot
        );
    }
}
