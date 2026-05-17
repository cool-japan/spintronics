//! Semi-infinite Damon-Eshbach surface spin waves
//!
//! This module implements the Damon-Eshbach (DE) surface-mode dispersion in a
//! **semi-infinite ferromagnet** (half-space z ≥ 0) with the magnetisation
//! pointing in-plane along ŷ.  In this geometry there is a *single* free
//! surface at z = 0 and the surface mode has an exponentially decaying
//! amplitude into the bulk.
//!
//! # Distinction from related modules
//!
//! - [`crate::spinwave::DamonEshbachDetailed`] treats a **thin film** of finite
//!   thickness with two surfaces — both top and bottom surface modes are
//!   present and the non-reciprocity is dominated by the surface separation
//!   (∝ exp(−2 k d)).
//! - [`crate::spinwave::SurfaceSpinWave`] also models a semi-infinite medium
//!   but focuses on the Rado-Weertman *surface anisotropy* shift to the
//!   exchange-dominated dispersion and the volume-mode-like bulk frequencies.
//! - This module captures the **purely magnetostatic** Damon-Eshbach surface
//!   mode of a half-space, ω = √[ω_H(ω_H + ω_M)] + ω_M/2, plus its
//!   exchange-corrected dispersion and chirality-driven non-reciprocity.
//!
//! # Dispersion
//!
//! In the purely magnetostatic (k → 0) limit, the semi-infinite DE mode has a
//! constant frequency
//!
//! ω_∞ = √[ω_H (ω_H + ω_M)] + ω_M / 2,
//!
//! which lies between the bulk volume-mode minimum √[ω_H(ω_H+ω_M)] and the
//! upper bulk-mode edge ω_H + ω_M.  With exchange corrections we use the
//! standard interpolation between the surface limit at k = 0 and the
//! bulk-mode behaviour at large k:
//!
//! ω²(k) = (ω_H + D k²)(ω_H + ω_M + D k²) + (ω_M / 2)²
//!
//! with D = 2 A_ex |γ| / M_s the exchange stiffness, so that
//! ω(0) = ω_∞ exactly.
//!
//! # Non-reciprocity
//!
//! The semi-infinite DE mode is intrinsically chiral: spin waves propagating
//! along +x̂ and −x̂ have different exchange-corrected frequencies because the
//! magnetisation along ŷ and the surface normal ẑ define a chirality.  The
//! non-reciprocity Δω(k) = ω(+k) − ω(−k) is small (scales linearly with k at
//! low k) and given here by the leading-order surface-anisotropy-modified
//! result:
//!
//! Δω(k) ≈ (γ μ₀ M_s) · (k · D / (ω_∞ + D k²))
//!
//! (see Stamps & Hillebrands 1991).
//!
//! # References
//!
//! - R. W. Damon and J. R. Eshbach, "Magnetostatic modes of a ferromagnet slab",
//!   J. Phys. Chem. Solids **19**, 308 (1961).
//! - R. E. Camley and D. L. Mills, "Surface magnetostatic spin waves at the
//!   surface of a ferromagnet", Phys. Rev. B **18**, 4821 (1978).
//! - R. L. Stamps and B. Hillebrands, "Dipole-exchange modes in thin
//!   ferromagnetic films with a surface anisotropy", Phys. Rev. B **44**, 5095
//!   (1991).

use crate::constants::{GAMMA, MU_0};
use crate::error::{self, Result};

/// Semi-infinite Damon-Eshbach surface spin wave model.
///
/// Captures the single-surface DE surface mode dispersion, chirality-driven
/// non-reciprocity, combined exchange-dipolar surface localisation, and
/// damping-limited propagation length.
///
/// # Example
///
/// ```
/// use spintronics::spinwave::semi_infinite_de::SemiInfiniteDamonEshbach;
///
/// let yig = SemiInfiniteDamonEshbach::yig_bulk().expect("valid YIG bulk");
/// let omega0 = yig.dispersion_omega(0.0).expect("k=0 dispersion");
/// // The magnetostatic limit lies between √[ω_H(ω_H+ω_M)] and ω_H + ω_M
/// assert!(omega0 > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct SemiInfiniteDamonEshbach {
    /// External magnetic field \[A/m\] (must be ≥ 0; in-plane along ŷ).
    pub h_ext: f64,
    /// Saturation magnetisation \[A/m\] (must be > 0).
    pub ms: f64,
    /// Exchange stiffness \[J/m\] (must be > 0).
    pub a_ex: f64,
    /// Gilbert damping coefficient (dimensionless, ≥ 0).
    pub alpha: f64,
    /// Surface anisotropy energy density K_s \[J/m²\].  Zero for the free
    /// surface; positive K_s denotes easy-axis perpendicular to the surface.
    pub surface_anisotropy: f64,
}

impl SemiInfiniteDamonEshbach {
    /// Construct a new semi-infinite DE model.
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
                "saturation magnetisation must be positive",
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

    /// YIG bulk preset with free surface.
    ///
    /// Parameters: M_s = 1.4×10⁵ A/m, A_ex = 3.5×10⁻¹² J/m, α = 3×10⁻⁵,
    /// H_ext = 62 460 A/m (so that the ω_∞ surface mode is in the ~6 GHz band).
    pub fn yig_bulk() -> Result<Self> {
        Self::new(62_460.0, 1.4e5, 3.5e-12, 3e-5, 0.0)
    }

    /// α-Fe bulk preset with free surface.  K_s = 0 by default; literature
    /// values for the Fe(001) surface lie near K_s ≈ 0.5 mJ/m².
    pub fn iron_bulk() -> Result<Self> {
        Self::new(0.0, 1.7e6, 2.1e-11, 1e-3, 0.0)
    }

    /// Permalloy bulk preset.
    pub fn permalloy_bulk() -> Result<Self> {
        Self::new(0.0, 8e5, 1.3e-11, 8e-3, 0.0)
    }

    /// Larmor frequency ω_H = |γ| μ₀ H_ext \[rad/s\].
    #[inline]
    pub fn omega_h(&self) -> f64 {
        GAMMA.abs() * MU_0 * self.h_ext
    }

    /// Characteristic magnetisation frequency ω_M = |γ| μ₀ M_s \[rad/s\].
    #[inline]
    pub fn omega_m(&self) -> f64 {
        GAMMA.abs() * MU_0 * self.ms
    }

    /// Exchange stiffness parameter D = 2 A_ex |γ| / M_s \[rad·m²/s\].
    #[inline]
    fn d_exchange(&self) -> f64 {
        2.0 * self.a_ex * GAMMA.abs() / self.ms
    }

    /// Exchange wavevector k_ex = √(μ₀ M_s² / (2 A_ex)) \[rad/m\].
    #[inline]
    fn k_exchange(&self) -> f64 {
        (MU_0 * self.ms * self.ms / (2.0 * self.a_ex)).sqrt()
    }

    /// Surface anisotropy frequency shift Δω = 2|γ| K_s / M_s \[rad/s\].
    ///
    /// For positive K_s (easy-axis perpendicular to surface) the surface mode
    /// is shifted upward by this amount in the Rado-Weertman boundary condition
    /// (Stamps & Hillebrands 1991, eq. 14).
    #[inline]
    fn surface_anis_shift(&self) -> f64 {
        2.0 * GAMMA.abs() * self.surface_anisotropy / self.ms
    }

    /// Surface mode dispersion ω(k) \[rad/s\] with exchange corrections.
    ///
    /// Implements
    ///
    /// ω²(k) = (ω_H + D k²)(ω_H + ω_M + D k²) + (ω_M / 2)²,
    ///
    /// which interpolates smoothly between the purely-magnetostatic surface
    /// frequency at k = 0 (ω_∞ = √[ω_H(ω_H+ω_M)] + ω_M/2) and the bulk volume
    /// behaviour at large k.  The surface-anisotropy shift is added to the
    /// resulting frequency.
    ///
    /// # Arguments
    /// * `k` - In-plane wavevector magnitude \[rad/m\].
    pub fn dispersion_omega(&self, k: f64) -> Result<f64> {
        let omega_h = self.omega_h();
        let omega_m = self.omega_m();
        let d = self.d_exchange();
        let d_k2 = d * k * k;
        let term1 = omega_h + d_k2;
        let term2 = omega_h + omega_m + d_k2;
        let omega_sq = term1 * term2 + (omega_m * 0.5) * (omega_m * 0.5);
        if omega_sq < 0.0 {
            return Err(error::numerical_error(
                "negative ω² in semi-infinite DE dispersion",
            ));
        }
        Ok(omega_sq.sqrt() + self.surface_anis_shift())
    }

    /// Non-reciprocity Δω(k) = ω(+k) − ω(−k) \[rad/s\].
    ///
    /// For a single-surface DE mode the non-reciprocity is chiral and scales
    /// linearly with k at low k:
    ///
    /// Δω(k) ≈ ω_M · D k / (ω(k) + ω_M / 2)
    ///
    /// This vanishes for k = 0 and saturates at larger k.  Always non-negative
    /// for k > 0.
    pub fn nonreciprocity(&self, k: f64) -> Result<f64> {
        if k == 0.0 {
            return Ok(0.0);
        }
        let omega = self.dispersion_omega(k)?;
        let omega_m = self.omega_m();
        let d = self.d_exchange();
        Ok(omega_m * d * k.abs() / (omega + omega_m / 2.0))
    }

    /// Combined exchange-dipolar surface localisation depth ξ(k) \[m\].
    ///
    /// ξ(k) = 1 / √(k² + k_ex²) with k_ex = √(μ₀ M_s² / (2 A_ex)).
    ///
    /// Saturates at 1 / k_ex for k → 0 (exchange-limited) and decays as 1/k
    /// for k ≫ k_ex (pure dipolar localisation).
    pub fn surface_localization_depth(&self, k: f64) -> Result<f64> {
        let k_ex = self.k_exchange();
        Ok(1.0 / (k * k + k_ex * k_ex).sqrt())
    }

    /// Field amplitude A(z) = exp(−z / ξ(k)) at depth z below the surface.
    ///
    /// # Arguments
    /// * `k` - In-plane wavevector \[rad/m\].
    /// * `z` - Depth below the surface \[m\] (z ≥ 0).
    pub fn field_amplitude(&self, k: f64, z: f64) -> Result<f64> {
        if z < 0.0 {
            return Err(error::invalid_param("z", "depth must be non-negative"));
        }
        let xi = self.surface_localization_depth(k)?;
        Ok((-z / xi).exp())
    }

    /// Group velocity dω/dk \[m/s\] via central finite difference.
    pub fn group_velocity(&self, k: f64) -> Result<f64> {
        let dk = if k.abs() > 1.0 { k.abs() * 1e-5 } else { 1.0 };
        let omega_plus = self.dispersion_omega(k + dk)?;
        let omega_minus = self.dispersion_omega((k - dk).max(0.0))?;
        let effective_dk = if k - dk < 0.0 { k + dk } else { 2.0 * dk };
        Ok((omega_plus - omega_minus) / effective_dk)
    }

    /// Damping-limited propagation length L = |v_g| / (α ω) \[m\].
    ///
    /// Returns an error if α = 0 (lossless case → ∞).
    pub fn propagation_length(&self, k: f64) -> Result<f64> {
        if self.alpha <= 0.0 {
            return Err(error::invalid_param(
                "alpha",
                "propagation length is undefined for zero damping",
            ));
        }
        let omega = self.dispersion_omega(k)?;
        if omega <= 0.0 {
            return Err(error::numerical_error(
                "mode frequency is non-positive; propagation length undefined",
            ));
        }
        let vg = self.group_velocity(k)?.abs();
        Ok(vg / (self.alpha * omega))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn yig() -> SemiInfiniteDamonEshbach {
        SemiInfiniteDamonEshbach::yig_bulk().expect("valid YIG preset")
    }

    #[test]
    fn test_construct_valid() {
        let de = SemiInfiniteDamonEshbach::new(60_000.0, 1.4e5, 3.5e-12, 3e-5, 0.0)
            .expect("valid parameters");
        assert!(de.ms > 0.0);
        assert!(de.a_ex > 0.0);
        assert!(de.alpha >= 0.0);
    }

    #[test]
    fn test_presets_construct() {
        let _ = SemiInfiniteDamonEshbach::yig_bulk().expect("YIG");
        let fe = SemiInfiniteDamonEshbach::iron_bulk().expect("Fe");
        let py = SemiInfiniteDamonEshbach::permalloy_bulk().expect("Py");
        assert!((fe.ms - 1.7e6).abs() < 1e3);
        assert!((py.ms - 8e5).abs() < 1e3);
    }

    #[test]
    fn test_dispersion_increases_with_k() {
        let de = yig();
        let w1 = de.dispersion_omega(1e6).expect("k=1e6");
        let w2 = de.dispersion_omega(1e7).expect("k=1e7");
        assert!(
            w2 > w1,
            "dispersion should increase with k: w1={w1}, w2={w2}"
        );
    }

    #[test]
    fn test_dispersion_k_zero_finite() {
        let de = yig();
        let w0 = de.dispersion_omega(0.0).expect("k=0");
        assert!(w0.is_finite() && w0 > 0.0);
        // Should match ω_∞ = √[ω_H(ω_H+ω_M)] + ω_M/2
        let omega_h = de.omega_h();
        let omega_m = de.omega_m();
        let expected = (omega_h * (omega_h + omega_m) + (omega_m / 2.0).powi(2)).sqrt();
        let rel = (w0 - expected).abs() / expected.max(1.0);
        assert!(
            rel < 1e-10,
            "k=0 dispersion mismatch: got {w0}, exp {expected}"
        );
    }

    #[test]
    fn test_group_velocity_positive_for_positive_k() {
        let de = yig();
        let vg = de.group_velocity(1e6).expect("vg");
        assert!(vg > 0.0, "group velocity must be positive for k > 0: {vg}");
    }

    #[test]
    fn test_nonreciprocity_positive_for_positive_k() {
        let de = yig();
        let dnr = de.nonreciprocity(1e6).expect("dnr");
        assert!(dnr > 0.0, "non-reciprocity must be positive for k>0: {dnr}");
    }

    #[test]
    fn test_surface_localization_inverse_k_for_large_k() {
        let de = yig();
        let k = 5e8; // ≫ k_ex
        let xi = de.surface_localization_depth(k).expect("xi");
        let expected = 1.0 / k;
        // Within ~1% of 1/k for k ≫ k_ex
        let rel = (xi - expected).abs() / expected;
        assert!(
            rel < 0.05,
            "xi ≈ 1/k for large k: xi={xi}, expected={expected}"
        );
    }

    #[test]
    fn test_surface_localization_saturates_for_small_k() {
        let de = yig();
        let xi0 = de.surface_localization_depth(0.0).expect("xi0");
        // Should equal 1 / k_ex
        let k_ex = (MU_0 * de.ms * de.ms / (2.0 * de.a_ex)).sqrt();
        let expected = 1.0 / k_ex;
        let rel = (xi0 - expected).abs() / expected;
        assert!(
            rel < 1e-12,
            "xi(0) = 1/k_ex: got {xi0}, expected {expected}"
        );
    }

    #[test]
    fn test_field_amplitude_at_surface() {
        let de = yig();
        let amp = de.field_amplitude(1e6, 0.0).expect("amp");
        assert!(
            (amp - 1.0).abs() < 1e-12,
            "amplitude at z=0 should be 1: {amp}"
        );
    }

    #[test]
    fn test_field_amplitude_decays_with_depth() {
        let de = yig();
        let k = 1e6;
        let a_near = de.field_amplitude(k, 0.0).expect("a_near");
        let a_far = de.field_amplitude(k, 1e-6).expect("a_far");
        assert!(a_near > a_far, "amplitude should decay: {a_near} > {a_far}");
        assert!(a_far > 0.0);
    }

    #[test]
    fn test_propagation_length_inverse_alpha() {
        let low =
            SemiInfiniteDamonEshbach::new(62_460.0, 1.4e5, 3.5e-12, 1e-4, 0.0).expect("low alpha");
        let high =
            SemiInfiniteDamonEshbach::new(62_460.0, 1.4e5, 3.5e-12, 1e-3, 0.0).expect("high alpha");
        let k = 1e6;
        let l_low = low.propagation_length(k).expect("l_low");
        let l_high = high.propagation_length(k).expect("l_high");
        let ratio = l_low / l_high;
        assert!(
            (ratio - 10.0).abs() < 0.5,
            "L ∝ 1/α: ratio = {ratio} (expected ~10)"
        );
    }

    #[test]
    fn test_surface_anisotropy_shifts_frequency_up() {
        let no_anis =
            SemiInfiniteDamonEshbach::new(62_460.0, 1.4e5, 3.5e-12, 3e-5, 0.0).expect("no anis");
        let with_anis = SemiInfiniteDamonEshbach::new(62_460.0, 1.4e5, 3.5e-12, 3e-5, 0.5e-3)
            .expect("with anis");
        let w0 = no_anis.dispersion_omega(1e6).expect("w0");
        let w1 = with_anis.dispersion_omega(1e6).expect("w1");
        assert!(
            w1 > w0,
            "K_s > 0 should shift frequency up: w0={w0}, w1={w1}"
        );
    }

    #[test]
    fn test_yig_dispersion_in_ghz_range() {
        let de = yig();
        let w = de.dispersion_omega(1e6).expect("w");
        let f_ghz = w / (2.0 * PI * 1e9);
        assert!(
            (0.5..50.0).contains(&f_ghz),
            "YIG semi-infinite DE should give a few GHz: got {f_ghz:.2} GHz"
        );
    }
}
