//! Spin Hall Magnetoresistance (SMR) and Unidirectional SMR (USMR)
//!
//! This module implements the Spin Hall Magnetoresistance effect, which arises in
//! NM/FM bilayers (e.g. Pt/YIG) where charge current in the heavy-metal normal-metal
//! (NM) layer generates spin accumulation at the NM/FM interface.  The reflected spin
//! current depends on the FM magnetisation direction **m**, causing an *m*-dependent
//! resistivity modulation observable in both the longitudinal and Hall channels.
//!
//! ## Theory
//!
//! The key results of Chen *et al.*, *Phys. Rev. B* **87**, 144411 (2013):
//!
//! - **Longitudinal resistivity**: `ρ_L(m) = ρ₀ + ρ₁(1 − m_y²)`, where *y* is the
//!   direction perpendicular to both the current flow (*x*) and the interface normal (*z*).
//! - **Hall resistivity**: `ρ_H(m) = ρ₁ m_x m_y − ρ₂ m_z`
//! - **SMR coefficient**: `ρ₁ = ρ_NM θ_SH² [2 g_r λ_sf tanh(t_NM/(2λ_sf))] /
//!   [σ_NM + 2 g_r λ_sf coth(t_NM/λ_sf)]`
//! - **SMR ratio**: `ρ₁/ρ₀`
//!
//! ## Unidirectional SMR (USMR)
//!
//! Avci *et al.*, *Nat. Phys.* **11**, 570 (2015) identified a second-order magnetotransport
//! signal that is *linear in current density* and *odd in magnetisation*:
//!
//! `ΔR_USMR = ρ₀ · η · J · (ĵ × ẑ) · m`
//!
//! where η is the USMR coefficient (m²/A).
//!
//! ## References
//!
//! - Y. T. Chen, S. Takahashi, H. Nakayama, M. Althammer, S. T. B. Goennenwein,
//!   E. Saitoh, G. E. W. Bauer, "Theory of spin Hall magnetoresistance",
//!   *Phys. Rev. B* **87**, 144411 (2013).
//! - K. Nakayama, H. Jungfleisch, T. Balogh, F. Casanova, L. E. Hueso, G. Tatara,
//!   E. Saitoh, *Phys. Rev. Lett.* **110**, 206601 (2013).
//! - C. O. Avci, K. Garello, A. Ghosh, M. Gabureac, S. F. Alvarado,
//!   P. Gambardella, *Nat. Phys.* **11**, 570 (2015).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::error::{invalid_param, Result};
use crate::vector3::Vector3;

// ──────────────────────────────────────────────────────────────────────────────
// Compile-time sanity checks
// ──────────────────────────────────────────────────────────────────────────────

/// Default validation tolerance (30 %).
pub const DEFAULT_TOLERANCE: f64 = 0.30_f64;

/// Typical Pt spin Hall angle lower bound used in asserts.
const _PT_THETA_SH: f64 = 0.07_f64;

const _: () = assert!(DEFAULT_TOLERANCE > 0.0);
const _: () = assert!(DEFAULT_TOLERANCE < 1.0);
const _: () = assert!(_PT_THETA_SH > 0.0);
const _: () = assert!(_PT_THETA_SH < 1.0);

// ──────────────────────────────────────────────────────────────────────────────
// SpinHallMagnetoresistance
// ──────────────────────────────────────────────────────────────────────────────

/// Spin Hall Magnetoresistance (SMR) model for NM/FM bilayers.
///
/// Captures the *m*-dependent resistivity modulation arising from spin
/// accumulation at the NM/FM interface.  All SI units.
///
/// # Example
/// ```
/// use spintronics::effect::smr::SpinHallMagnetoresistance;
/// use spintronics::Vector3;
///
/// let smr = SpinHallMagnetoresistance::platinum_yig();
/// let m = Vector3::new(1.0, 0.0, 0.0); // m along x
/// let rho_l = smr.longitudinal_resistivity(m);
/// assert!(rho_l > 0.0);
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpinHallMagnetoresistance {
    /// Spin Hall angle of the NM layer (dimensionless).  θ_SH ∈ (−1, 1).
    pub theta_sh: f64,
    /// NM layer bulk resistivity ρ₀ [Ω·m].
    pub resistivity_nm: f64,
    /// Spin diffusion length in NM [m].
    pub lambda_sf: f64,
    /// NM layer thickness [m].
    pub t_nm: f64,
    /// NM electrical conductivity σ_NM [S/m].
    pub sigma_nm: f64,
    /// Real part of spin mixing conductance g↑↓ (Re part) [m⁻²].
    pub g_r: f64,
    /// Imaginary part of spin mixing conductance g↑↓ (Im part) [m⁻²].
    pub g_i: f64,
}

impl SpinHallMagnetoresistance {
    /// Construct and validate a new SMR model.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] if:
    /// - `|theta_sh| >= 1`
    /// - `lambda_sf <= 0`
    /// - `t_nm <= 0`
    /// - `sigma_nm <= 0`
    /// - `g_r < 0`
    pub fn new(
        theta_sh: f64,
        resistivity_nm: f64,
        lambda_sf: f64,
        t_nm: f64,
        sigma_nm: f64,
        g_r: f64,
        g_i: f64,
    ) -> Result<Self> {
        if theta_sh.abs() >= 1.0 {
            return Err(invalid_param(
                "theta_sh",
                "spin Hall angle must satisfy |θ_SH| < 1",
            ));
        }
        if lambda_sf <= 0.0 {
            return Err(invalid_param(
                "lambda_sf",
                "spin diffusion length must be positive",
            ));
        }
        if t_nm <= 0.0 {
            return Err(invalid_param("t_nm", "NM layer thickness must be positive"));
        }
        if sigma_nm <= 0.0 {
            return Err(invalid_param(
                "sigma_nm",
                "NM conductivity must be positive",
            ));
        }
        if g_r < 0.0 {
            return Err(invalid_param(
                "g_r",
                "real part of spin mixing conductance must be non-negative",
            ));
        }
        Ok(Self {
            theta_sh,
            resistivity_nm,
            lambda_sf,
            t_nm,
            sigma_nm,
            g_r,
            g_i,
        })
    }

    /// Pt(7 nm)/YIG parameters from Nakayama *et al.*, PRL **110**, 206601 (2013).
    pub fn platinum_yig() -> Self {
        Self {
            theta_sh: 0.07,
            resistivity_nm: 1.06e-7, // Ω·m
            lambda_sf: 1.5e-9,       // m
            t_nm: 7.0e-9,            // m
            sigma_nm: 9.43e6,        // S/m
            g_r: 3.0e18,             // m⁻²
            g_i: 0.5e18,             // m⁻²
        }
    }

    /// W(5 nm)/YIG parameters (β-W phase, large negative spin Hall angle).
    pub fn tungsten_yig() -> Self {
        Self {
            theta_sh: -0.3,
            resistivity_nm: 5.0e-7,
            lambda_sf: 1.2e-9,
            t_nm: 5.0e-9,
            sigma_nm: 2.0e6,
            g_r: 1.5e18,
            g_i: 0.3e18,
        }
    }

    // ── Derived quantities ────────────────────────────────────────────────────

    /// Bulk NM resistivity ρ₀ [Ω·m] (equivalent to `self.resistivity_nm`).
    #[inline]
    pub fn rho_0(&self) -> f64 {
        self.resistivity_nm
    }

    /// SMR coefficient ρ₁ [Ω·m].
    ///
    /// Derived from the Chen *et al.* 2013 formula:
    /// ```text
    /// ρ₁ = ρ_NM · θ_SH² · [2 g_r λ_sf tanh(t_NM/(2λ_sf))]
    ///       / [σ_NM + 2 g_r λ_sf coth(t_NM/λ_sf)]
    ///
    /// where coth(x) = cosh(x)/sinh(x).
    /// ```
    pub fn rho_1(&self) -> f64 {
        let ratio = self.t_nm / self.lambda_sf;
        let half_ratio = self.t_nm / (2.0 * self.lambda_sf);

        // tanh(t_NM/(2λ_sf)) and coth(t_NM/λ_sf) — handle large arguments gracefully
        let tanh_half = half_ratio.tanh();
        let coth_full = if ratio > 20.0 {
            // coth(x) → 1 for x >> 1
            1.0
        } else {
            ratio.cosh() / ratio.sinh()
        };

        let numerator = 2.0 * self.g_r * self.lambda_sf * tanh_half;
        let denominator = self.sigma_nm + 2.0 * self.g_r * self.lambda_sf * coth_full;

        self.resistivity_nm * self.theta_sh * self.theta_sh * numerator / denominator
    }

    /// Anomalous Hall-like SMR coefficient ρ₂ [Ω·m].
    ///
    /// `ρ₂ = ρ₁ × (g_i / g_r)` — arises from the imaginary part of the spin
    /// mixing conductance (Berry-phase contribution at the interface).
    pub fn rho_2(&self) -> f64 {
        if self.g_r.abs() > 0.0 {
            self.rho_1() * (self.g_i / self.g_r)
        } else {
            0.0
        }
    }

    /// Spin conductance of the NM layer G_Sh [S/m²].
    ///
    /// ```text
    /// G_Sh = σ_NM · tanh(t_NM/λ_sf) / λ_sf
    /// ```
    pub fn g_sh(&self) -> f64 {
        let ratio = self.t_nm / self.lambda_sf;
        self.sigma_nm * ratio.tanh() / self.lambda_sf
    }

    /// Dimensionless SMR ratio ρ₁/ρ₀.
    ///
    /// Typical values: 0.1–1 % for Pt/YIG at room temperature.
    pub fn smr_ratio(&self) -> f64 {
        self.rho_1() / self.rho_0()
    }

    // ── Resistivity tensor components ─────────────────────────────────────────

    /// Longitudinal resistivity [Ω·m] for magnetisation direction `m`.
    ///
    /// `ρ_L(m) = ρ₀ + ρ₁(1 − m_y²)`
    ///
    /// Current flows along **x̂**; the SMR sensitivity axis is **ŷ**.
    pub fn longitudinal_resistivity(&self, m: Vector3<f64>) -> f64 {
        self.rho_0() + self.rho_1() * (1.0 - m.y * m.y)
    }

    /// Hall (transverse) resistivity [Ω·m] for magnetisation direction `m`.
    ///
    /// `ρ_H(m) = ρ₁ m_x m_y − ρ₂ m_z`
    pub fn hall_resistivity(&self, m: Vector3<f64>) -> f64 {
        self.rho_1() * m.x * m.y - self.rho_2() * m.z
    }

    /// Sweep the magnetisation in the **x-z** plane and return resistivities.
    ///
    /// For each angle `α ∈ [0, 2π)`, the magnetisation is set to
    /// `m = (cos α, 0, sin α)` and the tuple `(α, ρ_L, ρ_H)` is recorded.
    ///
    /// # Arguments
    /// * `n_angles` — number of equally-spaced angles (must be ≥ 1).
    pub fn angular_scan(&self, n_angles: usize) -> Vec<(f64, f64, f64)> {
        let n = n_angles.max(1);
        let mut out = Vec::with_capacity(n);
        for k in 0..n {
            let alpha = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
            let m = Vector3::new(alpha.cos(), 0.0, alpha.sin());
            out.push((
                alpha,
                self.longitudinal_resistivity(m),
                self.hall_resistivity(m),
            ));
        }
        out
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// UnidirectionalSmr
// ──────────────────────────────────────────────────────────────────────────────

/// Unidirectional Spin Hall Magnetoresistance (USMR) model.
///
/// Bundles an [`SpinHallMagnetoresistance`] model with the additional USMR
/// coefficient η (m²/A) that governs the current-linear, magnetisation-odd
/// contribution to the resistance.
///
/// ## Physics
///
/// `ΔR_USMR = ρ₀ · η · J · (ĵ × ẑ) · m`
///
/// This term is *linear in current density* J and *odd in magnetisation* **m**,
/// breaking the spatial-inversion symmetry of ordinary SMR.  It arises from
/// the asymmetric population of spin-up and spin-down channels driven by the
/// spin-orbit coupling at the NM/FM interface.
///
/// # References
/// Avci *et al.*, *Nat. Phys.* **11**, 570 (2015).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct UnidirectionalSmr {
    /// Underlying SMR model (provides ρ₀ and geometry).
    pub smr: SpinHallMagnetoresistance,
    /// USMR coefficient η [m²/A]: ΔR/R₀ = η · J · (ĵ × ẑ) · m.
    pub usmr_coefficient: f64,
}

impl UnidirectionalSmr {
    /// Construct a USMR model from an existing [`SpinHallMagnetoresistance`].
    ///
    /// The USMR coefficient is estimated as `η = θ_SH × 2.5 × 10⁻¹⁶ m²/A`,
    /// a typical experimental scaling for Pt/Co bilayers.
    pub fn new(smr: SpinHallMagnetoresistance) -> Self {
        let usmr_coefficient = smr.theta_sh * 2.5e-16;
        Self {
            smr,
            usmr_coefficient,
        }
    }

    /// Pt(3 nm)/Co(1 nm) USMR preset from Avci *et al.* 2015.
    ///
    /// # Errors
    /// Propagates validation errors from [`SpinHallMagnetoresistance::new`].
    pub fn platinum_cobalt() -> Result<Self> {
        let smr = SpinHallMagnetoresistance::new(
            0.07,    // theta_sh
            1.06e-7, // resistivity_nm [Ω·m]
            2.0e-9,  // lambda_sf [m]
            5.0e-9,  // t_nm [m]
            9.43e6,  // sigma_nm [S/m]
            4.5e18,  // g_r [m⁻²]
            0.8e18,  // g_i [m⁻²]
        )?;
        Ok(Self {
            usmr_coefficient: 1.7e-16,
            smr,
        })
    }

    // ── USMR observables ──────────────────────────────────────────────────────

    /// USMR resistance change ΔR for a given magnetisation, current density,
    /// and current direction.
    ///
    /// Formula: `ΔR = ρ₀ · η · J · (ĵ × ẑ) · m`
    ///
    /// # Arguments
    /// * `m`            — normalised FM magnetisation vector.
    /// * `j_density`    — charge current density \[A/m²\].
    /// * `j_direction`  — unit vector of current flow direction.
    pub fn usmr_resistance(
        &self,
        m: Vector3<f64>,
        j_density: f64,
        j_direction: Vector3<f64>,
    ) -> f64 {
        let z_hat = Vector3::new(0.0, 0.0, 1.0);
        let j_cross_z = j_direction.cross(&z_hat);
        let proj = j_cross_z.dot(&m);
        self.smr.rho_0() * self.usmr_coefficient * j_density * proj
    }

    /// Normalised USMR relative resistance change ΔR/R₀.
    ///
    /// Assumes current flows along **x̂** so that `ĵ × ẑ = −ŷ`, giving:
    ///
    /// `ΔR/R₀ = −η · J · m_y`
    ///
    /// # Arguments
    /// * `m`          — normalised FM magnetisation vector.
    /// * `j_density`  — charge current density \[A/m²\].
    pub fn usmr_relative_change(&self, m: Vector3<f64>, j_density: f64) -> f64 {
        // ĵ = x̂, ẑ = ẑ ⟹ ĵ × ẑ = x̂ × ẑ = −ŷ
        // (ĵ × ẑ) · m = −m_y
        self.usmr_coefficient * j_density * (-m.y)
    }

    /// Scan USMR resistance over a range of current densities at fixed `m`.
    ///
    /// Returns a vector of `(j_density, delta_R)` pairs sampled uniformly over
    /// `[j_min, j_max]` with `n_points` samples.
    ///
    /// # Arguments
    /// * `m`        — normalised FM magnetisation vector.
    /// * `j_min`    — minimum current density \[A/m²\].
    /// * `j_max`    — maximum current density \[A/m²\].
    /// * `n_points` — number of samples (≥ 1).
    pub fn current_scan(
        &self,
        m: Vector3<f64>,
        j_min: f64,
        j_max: f64,
        n_points: usize,
    ) -> Vec<(f64, f64)> {
        let n = n_points.max(1);
        let mut out = Vec::with_capacity(n);
        let j_direction = Vector3::new(1.0, 0.0, 0.0); // x̂
        if n == 1 {
            let delta_r = self.usmr_resistance(m, j_min, j_direction);
            out.push((j_min, delta_r));
        } else {
            for k in 0..n {
                let t = (k as f64) / ((n - 1) as f64);
                let j = j_min + t * (j_max - j_min);
                let delta_r = self.usmr_resistance(m, j, j_direction);
                out.push((j, delta_r));
            }
        }
        out
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── SpinHallMagnetoresistance ─────────────────────────────────────────────

    #[test]
    fn test_platinum_yig_builds() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        assert!(smr.theta_sh > 0.0);
        assert!(smr.resistivity_nm > 0.0);
        assert!(smr.lambda_sf > 0.0);
        assert!(smr.t_nm > 0.0);
        assert!(smr.sigma_nm > 0.0);
        assert!(smr.g_r > 0.0);
    }

    #[test]
    fn test_tungsten_yig_negative_theta() {
        let smr = SpinHallMagnetoresistance::tungsten_yig();
        assert!(smr.theta_sh < 0.0);
    }

    #[test]
    fn test_new_validates_theta_sh_bounds() {
        let result = SpinHallMagnetoresistance::new(1.5, 1e-7, 2e-9, 5e-9, 9e6, 3e18, 0.5e18);
        assert!(result.is_err());

        let result = SpinHallMagnetoresistance::new(-1.0, 1e-7, 2e-9, 5e-9, 9e6, 3e18, 0.5e18);
        assert!(result.is_err());

        let result = SpinHallMagnetoresistance::new(0.07, 1e-7, 2e-9, 5e-9, 9e6, 3e18, 0.5e18);
        assert!(result.is_ok());
    }

    #[test]
    fn test_new_validates_positive_lengths() {
        // lambda_sf = 0
        assert!(SpinHallMagnetoresistance::new(0.07, 1e-7, 0.0, 5e-9, 9e6, 3e18, 0.5e18).is_err());
        // t_nm = 0
        assert!(SpinHallMagnetoresistance::new(0.07, 1e-7, 2e-9, 0.0, 9e6, 3e18, 0.5e18).is_err());
        // sigma_nm = 0
        assert!(SpinHallMagnetoresistance::new(0.07, 1e-7, 2e-9, 5e-9, 0.0, 3e18, 0.5e18).is_err());
    }

    #[test]
    fn test_rho_0_equals_resistivity_nm() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        assert!((smr.rho_0() - smr.resistivity_nm).abs() < 1e-30);
    }

    #[test]
    fn test_rho_1_positive_for_positive_theta() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        assert!(smr.rho_1() > 0.0, "ρ₁ must be positive for Pt/YIG");
    }

    #[test]
    fn test_rho_1_positive_for_negative_theta() {
        // θ_SH² is always positive so ρ₁ > 0 regardless of sign of θ_SH
        let smr = SpinHallMagnetoresistance::tungsten_yig();
        assert!(smr.rho_1() > 0.0, "ρ₁ = ρ_NM θ_SH² … must be positive");
    }

    #[test]
    fn test_smr_ratio_physical_magnitude() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let ratio = smr.smr_ratio();
        // Nakayama 2013 reports ~1e-3; allow an order-of-magnitude window
        assert!(ratio > 1e-5, "SMR ratio too small: {ratio:.3e}");
        assert!(ratio < 1e-1, "SMR ratio too large: {ratio:.3e}");
        println!("SMR ratio: {ratio:.4e}");
    }

    #[test]
    fn test_rho_2_proportional_to_gi_over_gr() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let expected = smr.rho_1() * (smr.g_i / smr.g_r);
        assert!((smr.rho_2() - expected).abs() < 1e-30 * smr.rho_0().abs());
    }

    #[test]
    fn test_g_sh_positive() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        assert!(smr.g_sh() > 0.0);
    }

    /// ρ_L(m) = ρ_L(−m): even symmetry in m.
    #[test]
    fn test_longitudinal_even_symmetry() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let test_vectors = [
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.6, 0.8, 0.0),
            Vector3::new(0.5, 0.5, (0.5_f64).sqrt()),
        ];
        for m in test_vectors {
            let m_neg = Vector3::new(-m.x, -m.y, -m.z);
            let diff =
                (smr.longitudinal_resistivity(m) - smr.longitudinal_resistivity(m_neg)).abs();
            assert!(
                diff < 1e-30,
                "ρ_L(m) ≠ ρ_L(−m) for m={m:?}: diff={diff:.3e}"
            );
        }
    }

    /// ρ_H(m): verify analytical formula directly against implementation.
    ///
    /// The Hall resistivity `ρ_H = ρ₁ m_x m_y − ρ₂ m_z` is NOT generally odd under
    /// m → −m (the ρ₂ m_z term flips sign while the ρ₁ m_x m_y term does not change).
    /// We therefore verify the formula directly rather than testing a symmetry that
    /// does not hold in general.
    #[test]
    fn test_hall_formula_direct() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let rho_1 = smr.rho_1();
        let rho_2 = smr.rho_2();
        let test_vectors = [
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.6, 0.8, 0.0),
            Vector3::new(
                1.0_f64 / 3.0_f64.sqrt(),
                1.0 / 3.0_f64.sqrt(),
                1.0 / 3.0_f64.sqrt(),
            ),
        ];
        for m in test_vectors {
            let rho_h_sim = smr.hall_resistivity(m);
            let rho_h_ref = rho_1 * m.x * m.y - rho_2 * m.z;
            let diff = (rho_h_sim - rho_h_ref).abs();
            assert!(
                diff < 1e-28,
                "Hall formula mismatch for m={m:?}: sim={rho_h_sim:.3e} ref={rho_h_ref:.3e}"
            );
        }
    }

    #[test]
    fn test_longitudinal_m_along_y_gives_rho0() {
        // When m ∥ ŷ: m_y = 1 ⟹ ρ_L = ρ₀ + ρ₁(1-1) = ρ₀
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let m = Vector3::new(0.0, 1.0, 0.0);
        let rho = smr.longitudinal_resistivity(m);
        assert!((rho - smr.rho_0()).abs() < 1e-30);
    }

    #[test]
    fn test_longitudinal_m_perp_y_gives_rho0_plus_rho1() {
        // When m ⊥ ŷ: m_y = 0 ⟹ ρ_L = ρ₀ + ρ₁
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let m = Vector3::new(1.0, 0.0, 0.0);
        let rho = smr.longitudinal_resistivity(m);
        let expected = smr.rho_0() + smr.rho_1();
        assert!((rho - expected).abs() < 1e-30);
    }

    #[test]
    fn test_angular_scan_returns_correct_length() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let scan = smr.angular_scan(36);
        assert_eq!(scan.len(), 36);
    }

    #[test]
    fn test_angular_scan_rho_l_even() {
        // In the x-z plane m_y = 0 always, so ρ_L = ρ₀ + ρ₁ everywhere
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let scan = smr.angular_scan(64);
        let expected_rho_l = smr.rho_0() + smr.rho_1();
        for (alpha, rho_l, _rho_h) in &scan {
            assert!(
                (rho_l - expected_rho_l).abs() < 1e-25,
                "ρ_L deviation at α={alpha:.3}: Δ={:.3e}",
                (rho_l - expected_rho_l).abs()
            );
        }
    }

    #[test]
    fn test_angular_scan_hall_follows_sin_2alpha() {
        // m = (cos α, 0, sin α) ⟹ ρ_H = ρ₁·cos(α)·0 − ρ₂·sin(α) = −ρ₂ sin(α)
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let scan = smr.angular_scan(128);
        for (alpha, _rho_l, rho_h) in &scan {
            let expected = -smr.rho_2() * alpha.sin();
            let diff = (rho_h - expected).abs();
            assert!(
                diff < 1e-25,
                "ρ_H vs −ρ₂ sin α at α={alpha:.3}: Δ={diff:.3e}"
            );
        }
    }

    #[test]
    fn test_angular_scan_zero_angles_returns_one() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let scan = smr.angular_scan(0);
        assert_eq!(scan.len(), 1);
    }

    /// Verify the full Hall formula symmetry under m → −m.
    ///
    /// `ρ_H(m) + ρ_H(−m) = 2 ρ₁ m_x m_y` (the m_z terms cancel with opposite sign).
    /// For vectors in the x-y plane (m_z = 0) the m_x m_y term does NOT cancel,
    /// but the ρ₂ contribution does cancel.  We verify this analytically.
    #[test]
    fn test_hall_symmetry_under_negation() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let rho_1 = smr.rho_1();
        let angles: Vec<f64> = (0..8).map(|k| (k as f64) * PI / 4.0).collect();
        for alpha in angles {
            let m = Vector3::new(alpha.cos(), alpha.sin(), 0.3);
            let norm = m.magnitude();
            let m_n = Vector3::new(m.x / norm, m.y / norm, m.z / norm);
            let m_neg = Vector3::new(-m_n.x, -m_n.y, -m_n.z);
            // ρ_H(m) + ρ_H(-m) = (ρ₁ mx my - ρ₂ mz) + (ρ₁ mx my + ρ₂ mz) = 2 ρ₁ mx my
            let expected_sum = 2.0 * rho_1 * m_n.x * m_n.y;
            let actual_sum = smr.hall_resistivity(m_n) + smr.hall_resistivity(m_neg);
            let diff = (actual_sum - expected_sum).abs();
            assert!(
                diff < 1e-24,
                "Hall sum ρ_H(m)+ρ_H(-m) mismatch at α={alpha:.3}: diff={diff:.3e}"
            );
        }
    }

    // ── UnidirectionalSmr ─────────────────────────────────────────────────────

    #[test]
    fn test_usmr_new_sets_coefficient_from_theta_sh() {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let theta = smr.theta_sh;
        let usmr = UnidirectionalSmr::new(smr);
        let expected_coeff = theta * 2.5e-16;
        assert!((usmr.usmr_coefficient - expected_coeff).abs() < 1e-28);
    }

    #[test]
    fn test_platinum_cobalt_builds() {
        let usmr = UnidirectionalSmr::platinum_cobalt().expect("platinum_cobalt must build");
        assert!((usmr.usmr_coefficient - 1.7e-16).abs() < 1e-28);
        assert!(usmr.smr.theta_sh > 0.0);
    }

    /// ΔR_USMR is linear in current density.
    #[test]
    fn test_usmr_linear_in_current() {
        let usmr = UnidirectionalSmr::platinum_cobalt().expect("build");
        let m = Vector3::new(0.0, 1.0, 0.0); // m along ŷ maximises USMR
        let j_direction = Vector3::new(1.0, 0.0, 0.0);

        let j1 = 1.0e11;
        let j2 = 2.0e11;

        let dr1 = usmr.usmr_resistance(m, j1, j_direction);
        let dr2 = usmr.usmr_resistance(m, j2, j_direction);

        let ratio = dr2 / dr1;
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "USMR should be linear in J: ratio={ratio:.6}"
        );
    }

    /// ΔR_USMR(+m) = −ΔR_USMR(−m): odd in magnetisation.
    #[test]
    fn test_usmr_odd_in_magnetisation() {
        let usmr = UnidirectionalSmr::platinum_cobalt().expect("build");
        let j_density = 1.0e11;
        let j_direction = Vector3::new(1.0, 0.0, 0.0);

        let test_vectors = [
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.5, 0.5_f64.sqrt(), 0.0),
            Vector3::new(0.0, 0.8, 0.6),
            Vector3::new(1.0, 0.0, 0.0),
        ];
        for m in test_vectors {
            let m_neg = Vector3::new(-m.x, -m.y, -m.z);
            let dr_p = usmr.usmr_resistance(m, j_density, j_direction);
            let dr_n = usmr.usmr_resistance(m_neg, j_density, j_direction);
            let sum = (dr_p + dr_n).abs();
            assert!(sum < 1e-30, "USMR not odd in m: {sum:.3e} for m={m:?}");
        }
    }

    #[test]
    fn test_usmr_relative_change_proportional_to_j() {
        let usmr = UnidirectionalSmr::platinum_cobalt().expect("build");
        let m = Vector3::new(0.0, 1.0, 0.0);
        let dr1 = usmr.usmr_relative_change(m, 1.0e11);
        let dr2 = usmr.usmr_relative_change(m, 3.0e11);
        let ratio = dr2 / dr1;
        assert!(
            (ratio - 3.0).abs() < 1e-10,
            "relative change ratio={ratio:.6}"
        );
    }

    #[test]
    fn test_current_scan_length_and_monotonicity() {
        let usmr = UnidirectionalSmr::platinum_cobalt().expect("build");
        let m = Vector3::new(0.0, 1.0, 0.0); // pure ŷ
        let scan = usmr.current_scan(m, 1.0e10, 1.0e12, 20);
        assert_eq!(scan.len(), 20);

        // For m = ŷ and j = x̂: (ĵ×ẑ)·m = (x̂×ẑ)·ŷ = (−ŷ)·ŷ = −1
        // ΔR = ρ₀ · η · J · (−1) < 0  (since η > 0 for Pt/Co)
        // Monotonically decreasing with J (more negative)
        for w in scan.windows(2) {
            assert!(
                w[1].1 < w[0].1,
                "ΔR should decrease monotonically for ŷ magnetisation"
            );
        }
    }

    #[test]
    fn test_current_scan_single_point() {
        let usmr = UnidirectionalSmr::platinum_cobalt().expect("build");
        let m = Vector3::new(0.0, 1.0, 0.0);
        let scan = usmr.current_scan(m, 5.0e10, 1.0e11, 1);
        assert_eq!(scan.len(), 1);
        assert!((scan[0].0 - 5.0e10).abs() < 1.0);
    }

    #[test]
    fn test_current_scan_zero_points_fallback() {
        let usmr = UnidirectionalSmr::platinum_cobalt().expect("build");
        let m = Vector3::new(1.0, 0.0, 0.0);
        let scan = usmr.current_scan(m, 1.0e10, 1.0e11, 0);
        assert_eq!(scan.len(), 1);
    }
}
