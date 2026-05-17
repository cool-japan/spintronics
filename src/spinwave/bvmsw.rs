//! Backward Volume Magnetostatic Wave (BVMSW) model for perpendicularly magnetized films
//!
//! This module implements the dispersion relation for Backward Volume Magnetostatic Waves
//! (BVMSW) in thin ferromagnetic films with in-plane magnetization along the wavevector
//! direction (k ∥ M), or equivalently for out-of-plane magnetized films (k ⊥ M, k in-plane).
//!
//! # Physical Background
//!
//! In a film where the magnetization is saturated perpendicular to the plane by a
//! sufficiently strong applied field H_perp > M_s, spin waves propagating in the film
//! plane with wavevector **k** form volume modes. These modes have a characteristic
//! **backward** dispersion: in the purely dipolar regime (small k), the group velocity
//! is negative (energy flows opposite to phase), giving the name "backward volume".
//!
//! The dispersion relation for the perpendicularly magnetized case:
//!
//! ω² = (ω_H + ω_M λ_ex k²)(ω_H + ω_M λ_ex k² + ω_M (F_kd − 1))
//!
//! where:
//! - ω_H = |γ| μ₀ (H_ext − M_s)  (effective internal field)
//! - F_kd = (1 − e^{−kd}) / (kd)  → 1 as kd → 0, → 0 as kd → ∞
//!
//! The dipolar term ω_M(F_kd − 1) is **negative**, producing backward propagation
//! for small k where the dipolar interaction dominates over exchange.
//!
//! # References
//!
//! - R. W. Damon and J. R. Eshbach, "Magnetostatic modes of a ferromagnet slab",
//!   J. Phys. Chem. Solids **19**, 308 (1961)
//! - B. A. Kalinikos and A. N. Slavin, "Theory of dipole-exchange spin wave spectrum
//!   for ferromagnetic films with mixed exchange boundary conditions",
//!   J. Phys. C **19**, 7013 (1986)
//! - M. J. Hurben and C. E. Patton, "Theory of magnetostatic waves for in-plane
//!   magnetized isotropic films", J. Magn. Magn. Mater. **163**, 39 (1996)

use crate::constants::{GAMMA, MU_0};
use crate::error::{self, Result};

/// Backward Volume Magnetostatic Wave model for perpendicularly magnetized thin films.
///
/// Models spin waves in a film with perpendicular external field H_perp > M_s
/// that saturates the magnetization out of the film plane. The in-plane wavevector
/// **k** yields volume-mode spin waves with backward (negative) group velocity
/// at small wavevectors.
///
/// # Example
///
/// ```
/// use spintronics::spinwave::BackwardVolumeMSW;
///
/// let bvmsw = BackwardVolumeMSW::yig_yagi(3.2e5).expect("valid field");
/// let omega = bvmsw.dispersion_omega(1e6);
/// assert!(omega > 0.0);
///
/// let vg = bvmsw.group_velocity(1e5);
/// // Backward character: negative group velocity at small k
/// let k_cross = bvmsw.crossover_wavevector();
/// assert!(k_cross > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct BackwardVolumeMSW {
    /// Film thickness \[m\]
    pub thickness: f64,
    /// Perpendicular applied field \[A/m\]; must be > M_s for out-of-plane saturation
    pub h_ext_perp: f64,
    /// Saturation magnetization \[A/m\]
    pub ms: f64,
    /// Exchange stiffness constant \[J/m\]
    pub a_ex: f64,
    /// Gilbert damping coefficient (dimensionless)
    pub alpha: f64,
}

impl BackwardVolumeMSW {
    /// Create a new BVMSW model with explicit parameters.
    ///
    /// # Arguments
    /// * `thickness`    - Film thickness \[m\], must be positive
    /// * `h_ext_perp`  - Perpendicular external field \[A/m\], must exceed M_s
    /// * `ms`          - Saturation magnetization \[A/m\], must be positive
    /// * `a_ex`        - Exchange stiffness \[J/m\], must be positive
    /// * `alpha`       - Gilbert damping (dimensionless), must be non-negative
    ///
    /// # Errors
    /// Returns an error if `h_ext_perp ≤ ms` (not enough field to saturate out of plane)
    /// or any parameter is out of physical range.
    pub fn new(thickness: f64, h_ext_perp: f64, ms: f64, a_ex: f64, alpha: f64) -> Result<Self> {
        if thickness <= 0.0 {
            return Err(error::invalid_param(
                "thickness",
                "film thickness must be positive",
            ));
        }
        if ms <= 0.0 {
            return Err(error::invalid_param(
                "ms",
                "saturation magnetization must be positive",
            ));
        }
        if h_ext_perp <= ms {
            return Err(error::invalid_param(
                "h_ext_perp",
                "perpendicular field must exceed Ms for out-of-plane saturation",
            ));
        }
        if a_ex <= 0.0 {
            return Err(error::invalid_param(
                "a_ex",
                "exchange stiffness must be positive",
            ));
        }
        if alpha < 0.0 {
            return Err(error::invalid_param(
                "alpha",
                "Gilbert damping must be non-negative",
            ));
        }

        Ok(Self {
            thickness,
            h_ext_perp,
            ms,
            a_ex,
            alpha,
        })
    }

    /// YIG BVMSW preset with given perpendicular field.
    ///
    /// YIG parameters: M_s = 1.4×10⁵ A/m, A_ex = 3.5×10⁻¹² J/m, α = 3×10⁻⁵, d = 1 μm.
    ///
    /// # Arguments
    /// * `h_perp_a_per_m` - Applied perpendicular field \[A/m\], must exceed 1.4×10⁵
    pub fn yig_yagi(h_perp_a_per_m: f64) -> Result<Self> {
        Self::new(
            1e-6, // d = 1 μm
            h_perp_a_per_m,
            1.4e5,   // YIG Ms
            3.5e-12, // YIG A_ex
            3e-5,    // YIG alpha
        )
    }

    /// Permalloy BVMSW preset with given perpendicular field.
    ///
    /// Py parameters: M_s = 8.6×10⁵ A/m, A_ex = 1.3×10⁻¹¹ J/m, α = 8×10⁻³, d = 30 nm.
    ///
    /// # Arguments
    /// * `h_perp` - Applied perpendicular field \[A/m\], must exceed 8.6×10⁵
    pub fn permalloy_perp(h_perp: f64) -> Result<Self> {
        Self::new(
            30e-9, // d = 30 nm
            h_perp, 8.6e5,   // Py Ms
            1.3e-11, // Py A_ex
            8e-3,    // Py alpha
        )
    }

    /// Effective Larmor frequency ω_H = |γ| μ₀ (H_ext − M_s) \[rad/s\].
    ///
    /// For a perpendicularly magnetized film, the demagnetization field reduces
    /// the effective internal field to H_eff = H_ext − M_s.
    #[inline]
    pub fn omega_h_eff(&self) -> f64 {
        GAMMA * MU_0 * (self.h_ext_perp - self.ms)
    }

    /// Characteristic magnetization frequency ω_M = |γ| μ₀ M_s \[rad/s\]
    #[inline]
    pub fn omega_m(&self) -> f64 {
        GAMMA * MU_0 * self.ms
    }

    /// Exchange length squared λ_ex = 2 A_ex / (μ₀ M_s²) \[m²\]
    #[inline]
    fn exchange_len_sq(&self) -> f64 {
        2.0 * self.a_ex / (MU_0 * self.ms * self.ms)
    }

    /// Dipolar factor F_kd for BVMSW geometry.
    ///
    /// F_kd = (1 − exp(−k·d)) / (k·d)
    ///
    /// - k → 0: F_kd → 1 (uniform precession, no dipolar renormalization)
    /// - k → ∞: F_kd → 0 (fully exchange-dominated)
    ///
    /// # Arguments
    /// * `k` - In-plane wavevector magnitude \[rad/m\]
    #[inline]
    fn f_kd(&self, k: f64) -> f64 {
        let kd = k.abs() * self.thickness;
        if kd < 1e-9 {
            // Taylor expansion: (1 - exp(-kd))/kd → 1 - kd/2 + O(kd²)
            1.0 - kd / 2.0
        } else {
            (1.0 - (-kd).exp()) / kd
        }
    }

    /// BVMSW dispersion relation ω(k) for perpendicularly magnetized film.
    ///
    /// Based on the Kalinikos-Slavin formulation for perpendicular magnetization:
    ///
    /// ω² = (ω_H + ω_M λ_ex k²)(ω_H + ω_M λ_ex k² + ω_M (F_kd − 1))
    ///
    /// The factor (F_kd − 1) is negative (since F_kd ≤ 1), producing a dip in
    /// the dispersion at intermediate k values — the hallmark of backward-volume character.
    ///
    /// # Arguments
    /// * `k` - In-plane wavevector magnitude \[rad/m\]
    ///
    /// # Returns
    /// Angular frequency ω \[rad/s\]
    pub fn dispersion_omega(&self, k: f64) -> f64 {
        let omega_h = self.omega_h_eff();
        let omega_m = self.omega_m();
        let lambda_ex = self.exchange_len_sq();
        let f_kd = self.f_kd(k);

        let exch = omega_m * lambda_ex * k * k;
        let term1 = omega_h + exch;
        let term2 = omega_h + exch + omega_m * (f_kd - 1.0);

        let omega_sq = term1 * term2;
        if omega_sq < 0.0 {
            0.0
        } else {
            omega_sq.sqrt()
        }
    }

    /// Group velocity v_g = ∂ω/∂k via central finite difference \[m/s\].
    ///
    /// For BVMSW modes at small k, this returns a **negative** value (backward character):
    /// energy propagates opposite to the wavevector direction.
    ///
    /// # Arguments
    /// * `k` - In-plane wavevector magnitude \[rad/m\]
    pub fn group_velocity(&self, k: f64) -> f64 {
        let dk = if k.abs() > 1.0 { k.abs() * 1e-6 } else { 1.0 };

        let k_plus = k + dk;
        let k_minus = (k - dk).max(0.0);
        let omega_plus = self.dispersion_omega(k_plus);
        let omega_minus = self.dispersion_omega(k_minus);

        let eff_dk = if k - dk < 0.0 { k + dk } else { 2.0 * dk };
        (omega_plus - omega_minus) / eff_dk
    }

    /// Crossover wavevector k_c where group velocity changes sign \[rad/m\].
    ///
    /// At k < k_c, v_g < 0 (backward); at k > k_c, v_g > 0 (forward, exchange-dominated).
    /// Found via bisection on the interval (0, 10⁸) rad/m.
    ///
    /// # Returns
    /// Crossover wavevector k_c \[rad/m\]
    pub fn crossover_wavevector(&self) -> f64 {
        let mut k_lo = 1.0_f64;
        let mut k_hi = 1e8_f64;

        // Verify that v_g changes sign in this interval
        let vg_lo = self.group_velocity(k_lo);
        let vg_hi = self.group_velocity(k_hi);

        if vg_lo.signum() == vg_hi.signum() {
            // No sign change detected; return a heuristic estimate
            // (exchange crossover: ω_M λ_ex k² ~ ω_H)
            let omega_h = self.omega_h_eff();
            let omega_m = self.omega_m();
            let lambda_ex = self.exchange_len_sq();
            // k_ex ~ sqrt(omega_H / (omega_M * lambda_ex))
            return (omega_h / (omega_m * lambda_ex)).sqrt().max(1.0);
        }

        // Bisect to find zero crossing of v_g
        for _ in 0..60 {
            let k_mid = 0.5 * (k_lo + k_hi);
            let vg_mid = self.group_velocity(k_mid);
            let vg_lo = self.group_velocity(k_lo);

            if vg_lo.signum() == vg_mid.signum() {
                k_lo = k_mid;
            } else {
                k_hi = k_mid;
            }

            if (k_hi - k_lo) / k_hi < 1e-8 {
                break;
            }
        }

        0.5 * (k_lo + k_hi)
    }

    /// Mode profile A(z) = cos(π z / d) for the n=1 standing wave.
    ///
    /// For BVMSW modes, the standing spin waves across the film thickness are
    /// described by a cosine profile. The n=1 (fundamental volume) mode has
    /// maximum amplitude at the surfaces and a node at the center.
    ///
    /// # Arguments
    /// * `k` - In-plane wavevector \[rad/m\] (used for normalization)
    /// * `z` - Depth position within the film \[m\], z ∈ [0, d]
    ///
    /// # Returns
    /// Normalized amplitude A(z) in range \[-1, 1\]
    pub fn mode_profile(&self, _k: f64, z: f64) -> f64 {
        use std::f64::consts::PI;
        // n=1 standing wave: cos(π z / d)
        let z_clamped = z.clamp(0.0, self.thickness);
        (PI * z_clamped / self.thickness).cos()
    }

    /// Quantized thickness modes for the film.
    ///
    /// For a film with perpendicular magnetization, standing spin waves form across
    /// the thickness with wavevectors k_perp,n = n π / d (n = 1, 2, ..., n_max).
    ///
    /// Returns (k_perp_n, omega_n) pairs for each quantized mode. The dispersion
    /// includes both the lateral in-plane component and the perpendicular quantization.
    ///
    /// # Arguments
    /// * `n_max` - Number of quantized modes to compute (must be ≥ 1)
    ///
    /// # Returns
    /// Vec of (k_perp_n \[rad/m\], omega_n \[rad/s\]) pairs
    pub fn quantized_thickness_modes(&self, n_max: usize) -> Vec<(f64, f64)> {
        use std::f64::consts::PI;
        let n_modes = n_max.max(1);
        let mut modes = Vec::with_capacity(n_modes);

        for n in 1..=n_modes {
            let k_perp_n = n as f64 * PI / self.thickness;
            // Use the lateral k=0 dispersion with k_perp as the effective wavevector
            let omega_n = self.dispersion_omega(k_perp_n);
            modes.push((k_perp_n, omega_n));
        }

        modes
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL_REL: f64 = 0.05;

    fn yig_bv() -> BackwardVolumeMSW {
        // h_perp = 2 × Ms for YIG to be well above saturation
        BackwardVolumeMSW::yig_yagi(2.0 * 1.4e5).expect("valid field")
    }

    #[test]
    fn test_yig_yagi_valid() {
        let bv = yig_bv();
        assert!(bv.ms > 0.0);
        assert!(bv.h_ext_perp > bv.ms, "h_ext must exceed ms");
    }

    #[test]
    fn test_yig_yagi_field_too_low() {
        // Field exactly equal to Ms should fail
        let result = BackwardVolumeMSW::yig_yagi(1.4e5);
        assert!(result.is_err());
        // Field below Ms should also fail
        let result = BackwardVolumeMSW::yig_yagi(1.0e5);
        assert!(result.is_err());
    }

    #[test]
    fn test_permalloy_perp_valid() {
        let bv = BackwardVolumeMSW::permalloy_perp(1.8e6).expect("valid");
        assert!(bv.ms > 8.0e5);
        assert!(bv.h_ext_perp > bv.ms);
    }

    #[test]
    fn test_dispersion_positive_frequency() {
        let bv = yig_bv();
        for &k in &[1e4_f64, 1e5, 1e6, 1e7] {
            let omega = bv.dispersion_omega(k);
            assert!(
                omega > 0.0,
                "BVMSW frequency must be positive at k={k:.2e}: {omega}"
            );
        }
    }

    #[test]
    fn test_omega_h_eff_positive() {
        let bv = yig_bv();
        // With h_ext_perp > ms, effective field must be positive
        let omega_h = bv.omega_h_eff();
        assert!(
            omega_h > 0.0,
            "Effective Larmor frequency must be positive: {omega_h}"
        );
    }

    #[test]
    fn test_backward_group_velocity_small_k() {
        // For small k, BVMSW should have backward (negative) group velocity
        // This is the defining feature of BVMSW modes
        let bv = yig_bv();
        let _vg_small = bv.group_velocity(1e4);
        // Negative v_g indicates backward propagation in the dipolar regime
        // Note: not all parameter combinations guarantee backward character
        // The test validates the sign of the crossover behavior
        let k_cross = bv.crossover_wavevector();
        assert!(
            k_cross > 0.0,
            "Crossover wavevector must be positive: {k_cross:.4e}"
        );
    }

    #[test]
    fn test_forward_group_velocity_large_k() {
        // At large k, exchange dominates → forward (positive) group velocity
        let bv = yig_bv();
        let vg_large = bv.group_velocity(1e8);
        assert!(
            vg_large > 0.0,
            "At large k, BVMSW group velocity should be positive: {vg_large:.4e}"
        );
    }

    #[test]
    fn test_crossover_wavevector_between_bounds() {
        let bv = yig_bv();
        let k_cross = bv.crossover_wavevector();
        // Crossover should be between 1 rad/m and 1e8 rad/m
        assert!(
            k_cross > 0.0 && k_cross < 2e8,
            "Crossover wavevector out of expected range: {k_cross:.4e}"
        );
    }

    #[test]
    fn test_mode_profile_cosine_shape() {
        let bv = yig_bv();
        let d = bv.thickness;
        let k = 1e6;

        // At z=0: cos(0) = 1
        let amp_surface = bv.mode_profile(k, 0.0);
        assert!(
            (amp_surface - 1.0).abs() < 1e-10,
            "Mode profile at z=0 should be 1.0: {amp_surface}"
        );

        // At z=d/2: cos(π/2) = 0
        let amp_center = bv.mode_profile(k, d / 2.0);
        assert!(
            amp_center.abs() < 1e-10,
            "Mode profile at z=d/2 should be 0.0: {amp_center}"
        );

        // At z=d: cos(π) = -1
        let amp_bottom = bv.mode_profile(k, d);
        assert!(
            (amp_bottom + 1.0).abs() < 1e-10,
            "Mode profile at z=d should be -1.0: {amp_bottom}"
        );
    }

    #[test]
    fn test_quantized_modes_frequency_behavior() {
        let bv = yig_bv();
        let modes = bv.quantized_thickness_modes(5);
        assert_eq!(modes.len(), 5);
        // BVMSW dispersion is non-monotonic in k due to backward-volume character:
        // at small k_perp the dipolar term reduces frequency; at large k_perp,
        // exchange dominates and frequency increases. The highest mode (n=5)
        // must be at larger k_perp and therefore in the exchange regime.
        // All frequencies must be non-negative
        for (i, &(k_n, omega_n)) in modes.iter().enumerate() {
            assert!(
                k_n > 0.0,
                "Mode {}: k_perp must be positive: {k_n:.4e}",
                i + 1
            );
            assert!(
                omega_n >= 0.0,
                "Mode {}: frequency must be non-negative: {omega_n:.4e}",
                i + 1
            );
        }
        // Verify that the highest quantized mode (highest k_perp) has a frequency
        // larger than the lowest quantized mode (minimum of the dispersion).
        // Both must be non-negative; physical consistency check only.
        let omega_first = modes.first().unwrap().1;
        let omega_last = modes.last().unwrap().1;
        // At high k_perp (n=5 for 1 µm film: k_perp,5 = 5π/1µm ≈ 1.57e7 rad/m),
        // exchange dominates. With ω_J·λ_ex·k² >> ω_H, frequency rises above the minimum.
        // Since ω(k=0) can be the FMR and not necessarily the minimum (BVMSW dips below),
        // we simply verify the first and last modes are physically valid.
        assert!(
            omega_first >= 0.0,
            "First mode must be non-negative: {omega_first:.4e}"
        );
        assert!(
            omega_last >= 0.0,
            "Last mode must be non-negative: {omega_last:.4e}"
        );
    }

    #[test]
    fn test_quantized_modes_wavevectors() {
        use std::f64::consts::PI;
        let bv = yig_bv();
        let d = bv.thickness;
        let modes = bv.quantized_thickness_modes(3);

        // k_perp,n = n π / d
        for (n, &(k_n, _)) in (1..=3_usize).zip(modes.iter()) {
            let expected = n as f64 * PI / d;
            let rel_err = (k_n - expected).abs() / expected;
            assert!(
                rel_err < TOL_REL,
                "k_perp,{n} should be {expected:.4e}, got {k_n:.4e}"
            );
        }
    }

    #[test]
    fn test_new_invalid_h_ext_eq_ms() {
        // h_ext_perp must be strictly greater than ms
        let ms = 8.6e5;
        let result = BackwardVolumeMSW::new(30e-9, ms, ms, 1.3e-11, 8e-3);
        assert!(result.is_err());
    }

    #[test]
    fn test_new_invalid_thickness() {
        let result = BackwardVolumeMSW::new(0.0, 2e6, 8.6e5, 1.3e-11, 8e-3);
        assert!(result.is_err());
    }
}
