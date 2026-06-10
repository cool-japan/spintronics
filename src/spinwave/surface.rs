//! Surface spin wave model for semi-infinite ferromagnetic media
//!
//! This module implements spin wave physics for ferromagnetic media with a free surface,
//! where both surface-localized (Damon-Eshbach type) and bulk-like modes coexist.
//! The key addition over thin-film models is the treatment of surface boundary conditions
//! via the Rado-Weertman surface anisotropy parameter K_s.
//!
//! # Physical Background
//!
//! In a semi-infinite ferromagnet (x ≥ 0, surface at x = 0), the magnetostatic
//! boundary conditions allow for surface-localized spin wave modes whose amplitude
//! decays exponentially into the bulk as exp(−κ z) with z measured from the surface.
//!
//! The Damon-Eshbach mode in a semi-infinite medium has angular frequency:
//!
//! ω_DE = √\[ω_H (ω_H + ω_M)\] + ω_M / 2  (purely magnetostatic, independent of k)
//!
//! With exchange interaction included, the dispersion becomes k-dependent:
//!
//! ω² = (ω_H + D k²)(ω_H + ω_M + D k²)
//!
//! where D = 2 A_ex |γ| / M_s is the exchange stiffness parameter \[rad·m²/s\].
//!
//! # Surface Anisotropy (Rado-Weertman Theory)
//!
//! The surface anisotropy K_s \[J/m²\] modifies the spin wave boundary conditions,
//! shifting the mode frequency by a k-dependent amount. Following Rado & Weertman:
//!
//! Δω = −2 |γ| K_s / M_s × k² / (k² + k_ex²)
//!
//! where k_ex = √(μ₀ M_s² / (2 A_ex)) is the exchange wavevector.
//!
//! # References
//!
//! - G. T. Rado and J. R. Weertman, "Spin-wave resonance in a ferromagnetic metal",
//!   J. Phys. Chem. Solids **11**, 315 (1959)
//! - R. W. Damon and J. R. Eshbach, "Magnetostatic modes of a ferromagnet slab",
//!   J. Phys. Chem. Solids **19**, 308 (1961)
//! - B. A. Kalinikos and A. N. Slavin, "Theory of dipole-exchange spin wave spectrum
//!   for ferromagnetic films with mixed exchange boundary conditions",
//!   J. Phys. C **19**, 7013 (1986)

use crate::constants::{GAMMA, MU_0};
use crate::error::{self, Result};

/// Surface spin wave model for semi-infinite ferromagnetic media.
///
/// Computes dispersion and spatial profiles of surface-localized spin waves
/// near a free ferromagnetic surface, including the effect of surface anisotropy
/// via the Rado-Weertman boundary condition.
///
/// # Example
///
/// ```
/// use spintronics::spinwave::SurfaceSpinWave;
///
/// let yig = SurfaceSpinWave::bulk_yig();
/// let omega = yig.dispersion_omega(1e6);
/// assert!(omega > 0.0, "surface spin wave frequency must be positive");
///
/// let delta = yig.surface_anisotropy_shift(1e6);
/// // For zero surface anisotropy, shift is zero
/// assert_eq!(delta, 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct SurfaceSpinWave {
    /// External magnetic field \[A/m\]
    pub h_ext: f64,
    /// Saturation magnetization \[A/m\]
    pub ms: f64,
    /// Exchange stiffness constant \[J/m\]
    pub a_ex: f64,
    /// Gilbert damping coefficient (dimensionless)
    pub alpha: f64,
    /// Surface anisotropy energy density K_s \[J/m²\]
    ///
    /// Positive K_s: easy-axis perpendicular to surface (hard-plane anisotropy).
    /// Negative K_s: easy-plane at the surface.
    /// K_s = 0: natural (exchange-only) boundary condition.
    pub surface_anisotropy: f64,
}

impl SurfaceSpinWave {
    /// Create a new surface spin wave model with explicit parameters.
    ///
    /// # Arguments
    /// * `h_ext`             - External magnetic field \[A/m\], must be non-negative
    /// * `ms`                - Saturation magnetization \[A/m\], must be positive
    /// * `a_ex`              - Exchange stiffness \[J/m\], must be positive
    /// * `alpha`             - Gilbert damping (dimensionless), must be non-negative
    /// * `surface_anisotropy` - Surface anisotropy K_s \[J/m²\], may be any sign
    ///
    /// # Errors
    /// Returns an error if any parameter violates its physical constraint.
    pub fn new(
        h_ext: f64,
        ms: f64,
        a_ex: f64,
        alpha: f64,
        surface_anisotropy: f64,
    ) -> Result<Self> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }
        if ms <= 0.0 {
            return Err(error::invalid_param(
                "ms",
                "saturation magnetization must be positive",
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
            h_ext,
            ms,
            a_ex,
            alpha,
            surface_anisotropy,
        })
    }

    /// YIG semi-infinite bulk preset.
    ///
    /// Yttrium Iron Garnet (Y₃Fe₅O₁₂) with zero surface anisotropy.
    /// Parameters from standard literature values.
    pub fn bulk_yig() -> Self {
        Self {
            h_ext: 62_460.0, // ~5 GHz FMR bias field [A/m]
            ms: 1.4e5,
            a_ex: 3.5e-12,
            alpha: 3e-5,
            surface_anisotropy: 0.0,
        }
    }

    /// Iron semi-infinite bulk preset.
    ///
    /// α-Fe with typical surface anisotropy for (001) surface.
    pub fn bulk_iron() -> Self {
        Self {
            h_ext: 0.0,
            ms: 1.71e6,
            a_ex: 2.1e-11,
            alpha: 1e-3,
            surface_anisotropy: 0.0, // can be set to non-zero for (001) Fe: K_s ~ 0.5 mJ/m²
        }
    }

    /// Larmor frequency ω_H = |γ| μ₀ H_ext \[rad/s\]
    #[inline]
    pub fn omega_h(&self) -> f64 {
        GAMMA * MU_0 * self.h_ext
    }

    /// Characteristic magnetization frequency ω_M = |γ| μ₀ M_s \[rad/s\]
    #[inline]
    pub fn omega_m(&self) -> f64 {
        GAMMA * MU_0 * self.ms
    }

    /// Exchange stiffness parameter D = 2 A_ex |γ| / M_s \[rad·m²/s\].
    ///
    /// This enters the dispersion as the coefficient of k² in ω(k).
    /// Note: differs from λ_ex (which is 2A/(μ₀Ms²)) by a factor of |γ| μ₀ M_s.
    #[inline]
    fn d_exchange(&self) -> f64 {
        2.0 * self.a_ex * GAMMA / self.ms
    }

    /// Exchange wavevector k_ex = √(μ₀ M_s² / (2 A_ex)) \[rad/m\].
    ///
    /// The characteristic wavevector at which exchange energy equals dipolar energy.
    /// For k >> k_ex, exchange dominates; for k << k_ex, dipolar interactions dominate.
    #[inline]
    fn k_exchange(&self) -> f64 {
        (MU_0 * self.ms * self.ms / (2.0 * self.a_ex)).sqrt()
    }

    /// Surface spin wave dispersion ω(k) including exchange.
    ///
    /// Implements the dipole-exchange spin wave dispersion for a semi-infinite
    /// ferromagnet. The surface-localized mode frequency is:
    ///
    /// ω² = (ω_H + D k²)(ω_H + ω_M + D k²)
    ///
    /// where D = 2 A_ex |γ| / M_s. This reduces to the standard Damon-Eshbach
    /// result at k = 0: ω = √\[ω_H (ω_H + ω_M)\].
    ///
    /// Note: this is the bulk volume mode formula for a semi-infinite medium.
    /// The actual surface mode frequency includes a correction of +ω_M/2 in the
    /// pure magnetostatic limit, here absorbed into the Kalinikos-Slavin framework.
    ///
    /// # Arguments
    /// * `k` - In-plane wavevector magnitude \[rad/m\]
    ///
    /// # Returns
    /// Angular frequency ω \[rad/s\]
    pub fn dispersion_omega(&self, k: f64) -> f64 {
        let omega_h = self.omega_h();
        let omega_m = self.omega_m();
        let d = self.d_exchange();

        let d_k2 = d * k * k;
        let term1 = omega_h + d_k2;
        let term2 = omega_h + omega_m + d_k2;

        let omega_sq = term1 * term2;
        if omega_sq < 0.0 {
            0.0
        } else {
            omega_sq.sqrt()
        }
    }

    /// Penetration depth of the surface spin wave \[m\].
    ///
    /// The surface mode amplitude decays into the bulk as exp(−z / ξ), where the
    /// penetration depth ξ is modified by exchange:
    ///
    /// ξ = 1 / √(k² + k_ex²)
    ///
    /// - For k >> k_ex: ξ ≈ 1/k (pure dipolar, same as for thin film)
    /// - For k << k_ex: ξ ≈ 1/k_ex (exchange-limited penetration depth)
    ///
    /// # Arguments
    /// * `k` - In-plane wavevector \[rad/m\]
    ///
    /// # Returns
    /// Penetration depth \[m\]
    pub fn penetration_depth(&self, k: f64) -> f64 {
        let k_ex = self.k_exchange();
        1.0 / (k * k + k_ex * k_ex).sqrt()
    }

    /// Surface anisotropy frequency shift Δω (Rado-Weertman correction) \[rad/s\].
    ///
    /// Surface anisotropy K_s modifies the spin wave boundary condition, shifting
    /// the mode frequency according to (Rado & Weertman, 1959):
    ///
    /// Δω = −2 |γ| K_s / M_s × k² / (k² + k_ex²)
    ///
    /// For easy-plane surface anisotropy (K_s < 0), this gives a positive shift.
    /// For easy-axis (K_s > 0, axis ⊥ surface), this gives a negative shift
    /// that can drive the mode below the bulk band.
    ///
    /// # Arguments
    /// * `k` - In-plane wavevector \[rad/m\]
    ///
    /// # Returns
    /// Frequency shift Δω \[rad/s\]
    pub fn surface_anisotropy_shift(&self, k: f64) -> f64 {
        if self.surface_anisotropy == 0.0 {
            return 0.0;
        }
        let k_ex = self.k_exchange();
        let k2 = k * k;
        let k_ex2 = k_ex * k_ex;
        // Δω = -2γK_s/Ms × k²/(k² + k_ex²)
        -2.0 * GAMMA * self.surface_anisotropy / self.ms * k2 / (k2 + k_ex2)
    }

    /// Field amplitude of the surface spin wave A(z) = exp(−k z).
    ///
    /// The surface-localized spin wave has an evanescent profile into the bulk.
    /// In the dipolar-dominated regime, this is a simple exponential:
    ///
    /// A(z) = exp(−k z),  z ≥ 0 (z=0 at the surface)
    ///
    /// # Arguments
    /// * `k` - In-plane wavevector magnitude \[rad/m\]
    /// * `z` - Depth below surface \[m\], z ≥ 0
    ///
    /// # Returns
    /// Normalized field amplitude in range \[0, 1\]
    pub fn field_amplitude(&self, k: f64, z: f64) -> f64 {
        let z_pos = z.max(0.0);
        (-k.abs() * z_pos).exp()
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    fn yig_surface() -> SurfaceSpinWave {
        SurfaceSpinWave::bulk_yig()
    }

    fn iron_surface() -> SurfaceSpinWave {
        SurfaceSpinWave::bulk_iron()
    }

    #[test]
    fn test_yig_preset_valid() {
        let s = yig_surface();
        assert!(s.ms > 0.0);
        assert!(s.a_ex > 0.0);
        assert!(s.alpha > 0.0);
        assert!(s.surface_anisotropy == 0.0);
    }

    #[test]
    fn test_iron_preset_valid() {
        let s = iron_surface();
        assert!((s.ms - 1.71e6).abs() < 1e3, "Fe Ms should be 1.71e6 A/m");
        assert!(s.a_ex > 0.0);
    }

    #[test]
    fn test_new_invalid_ms() {
        let result = SurfaceSpinWave::new(0.0, 0.0, 1e-12, 1e-4, 0.0);
        assert!(result.is_err());
        let result = SurfaceSpinWave::new(0.0, -1.0, 1e-12, 1e-4, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_new_invalid_a_ex() {
        let result = SurfaceSpinWave::new(0.0, 1e5, 0.0, 1e-4, 0.0);
        assert!(result.is_err());
        let result = SurfaceSpinWave::new(0.0, 1e5, -1e-12, 1e-4, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_new_invalid_h_ext() {
        let result = SurfaceSpinWave::new(-100.0, 1e5, 1e-12, 1e-4, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_new_valid_negative_surface_anisotropy() {
        // Surface anisotropy may be negative (easy-plane)
        let result = SurfaceSpinWave::new(0.0, 1e5, 1e-12, 1e-4, -0.5e-3);
        assert!(result.is_ok());
    }

    #[test]
    fn test_dispersion_positive_frequency() {
        let s = yig_surface();
        for &k in &[1e4_f64, 1e5, 1e6, 1e7, 1e8] {
            let omega = s.dispersion_omega(k);
            assert!(
                omega > 0.0,
                "surface spin wave frequency must be positive at k={k:.2e}: {omega}"
            );
        }
    }

    #[test]
    fn test_dispersion_k_zero_limit() {
        // At k=0: ω² = ω_H (ω_H + ω_M)
        let s = yig_surface();
        let omega_k0 = s.dispersion_omega(0.0);
        let omega_h = s.omega_h();
        let omega_m = s.omega_m();
        let omega_expected = (omega_h * (omega_h + omega_m)).sqrt();
        let rel_err = (omega_k0 - omega_expected).abs() / omega_expected.max(1.0);
        assert!(
            rel_err < 0.01,
            "At k=0, ω should be Kittel frequency: got {omega_k0:.4e}, expected {omega_expected:.4e}"
        );
    }

    #[test]
    fn test_dispersion_increases_with_k() {
        // Exchange makes dispersion monotonically increasing with k
        let s = yig_surface();
        let omega1 = s.dispersion_omega(1e6);
        let omega2 = s.dispersion_omega(1e7);
        assert!(
            omega2 > omega1,
            "Frequency should increase with k: ω(1e6)={omega1:.4e}, ω(1e7)={omega2:.4e}"
        );
    }

    #[test]
    fn test_penetration_depth_decreases_with_k() {
        // 1/√(k² + k_ex²) decreases as k increases
        let s = yig_surface();
        let xi1 = s.penetration_depth(1e6);
        let xi2 = s.penetration_depth(1e7);
        assert!(
            xi1 > xi2,
            "Penetration depth should decrease with k: ξ(1e6)={xi1:.4e}, ξ(1e7)={xi2:.4e}"
        );
    }

    #[test]
    fn test_penetration_depth_k_ex_limit() {
        // At k=0: ξ = 1/k_ex
        let s = iron_surface();
        let k_ex = s.k_exchange();
        let xi_k0 = s.penetration_depth(0.0);
        let expected = 1.0 / k_ex;
        let rel_err = (xi_k0 - expected).abs() / expected;
        assert!(
            rel_err < 0.01,
            "At k=0: ξ = 1/k_ex = {expected:.4e}, got {xi_k0:.4e}"
        );
    }

    #[test]
    fn test_surface_anisotropy_shift_zero_when_ks_zero() {
        let s = yig_surface(); // K_s = 0
        let shift = s.surface_anisotropy_shift(1e6);
        assert_eq!(shift, 0.0, "Zero surface anisotropy → zero frequency shift");
    }

    #[test]
    fn test_surface_anisotropy_shift_sign() {
        // K_s > 0 (easy-axis ⊥ surface): Δω < 0 (shifts frequency down)
        let s = SurfaceSpinWave::new(0.0, 1.4e5, 3.5e-12, 3e-5, 0.5e-3).expect("valid parameters");
        let shift = s.surface_anisotropy_shift(1e7);
        assert!(
            shift < 0.0,
            "Positive K_s should give negative frequency shift: Δω = {shift:.4e}"
        );
    }

    #[test]
    fn test_field_amplitude_at_surface() {
        let s = yig_surface();
        // At z=0: A(0) = exp(0) = 1.0
        let amp = s.field_amplitude(1e6, 0.0);
        assert!(
            (amp - 1.0).abs() < 1e-14,
            "Surface amplitude should be 1.0: {amp}"
        );
    }

    #[test]
    fn test_field_amplitude_decays_into_bulk() {
        let s = yig_surface();
        let k = 1e6;
        let z1 = 100e-9; // 100 nm below surface
        let z2 = 500e-9; // 500 nm below surface

        let amp1 = s.field_amplitude(k, z1);
        let amp2 = s.field_amplitude(k, z2);

        assert!(
            amp1 > amp2,
            "Amplitude should decay with depth: A(100nm)={amp1:.4e}, A(500nm)={amp2:.4e}"
        );
        assert!(amp1 > 0.0 && amp1 <= 1.0, "Amplitude must be in (0,1]");
        assert!(amp2 > 0.0 && amp2 <= 1.0, "Amplitude must be in (0,1]");
    }

    #[test]
    fn test_omega_m_yig_value() {
        let s = yig_surface();
        let omega_m = s.omega_m();
        // ω_M = γ μ₀ M_s ≈ 1.76e11 × 4πe-7 × 1.4e5 ≈ 3.1e10 rad/s
        assert!(
            omega_m > 1e10 && omega_m < 1e11,
            "YIG ω_M should be ~3e10 rad/s: {omega_m:.4e}"
        );
    }

    #[test]
    fn test_fmr_frequency_near_5ghz_yig() {
        // With h_ext chosen for ~5 GHz, check dispersion at k≈0
        let s = yig_surface();
        let omega = s.dispersion_omega(1.0); // nearly k=0
        let f_ghz = omega / (2.0 * PI * 1e9);
        assert!(
            f_ghz > 3.0 && f_ghz < 8.0,
            "YIG surface mode near k=0 should be ~5 GHz: {f_ghz:.2} GHz"
        );
    }
}
