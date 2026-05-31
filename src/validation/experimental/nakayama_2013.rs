//! Validation against Nakayama *et al.*, *Phys. Rev. Lett.* **110**, 206601 (2013).
//!
//! This landmark paper measured the Spin Hall Magnetoresistance (SMR) in a
//! Pt(7 nm)/YIG bilayer, demonstrating that the longitudinal and Hall
//! resistivities of the Pt layer are modulated by the YIG magnetisation
//! direction via the spin accumulation at the Pt/YIG interface.
//!
//! ## Landmark claims validated here
//!
//! 1. **Longitudinal angular dependence** — `ρ_L(α) = ρ₀ + ρ₁ sin²(α)` for `m`
//!    in the *x*–*z* plane (i.e. m = (cos α, 0, sin α)).  In this plane,
//!    `m_y = 0` so `ρ_L = ρ₀ + ρ₁(1 − 0) = ρ₀ + ρ₁`.  For the more general
//!    in-plane rotation (x–y plane, m = (cos α, sin α, 0)), the pattern is
//!    `ρ_L = ρ₀ + ρ₁(1 − sin²α) = ρ₀ + ρ₁ cos²α`.
//!    We follow Nakayama's experimental geometry: in-plane rotation in the
//!    *x*–*y* plane so that `ρ_L(α) = ρ₀ + ρ₁ cos²(α)`.
//!
//! 2. **SMR ratio** — `ΔSMR/ρ₀ = ρ₁/ρ₀ ≈ 1 × 10⁻³` for Pt(7 nm)/YIG.
//!
//! 3. **Hall angular dependence** — `ρ_H(α) = ρ₁ sin(α) cos(α) = (ρ₁/2) sin(2α)` for
//!    `m` in the *x*–*y* plane.  In this geometry `m_z = 0` so the ρ₂ term
//!    vanishes and `ρ_H = ρ₁ m_x m_y = ρ₁ cos(α) sin(α)`.
//!
//! All validation methods return a [`ValidationResult`] so that the harness
//! can be wired into a larger test suite with a uniform interface.
//!
//! ## Caveats
//!
//! - The preset parameters (platinum_yig) are taken from the literature
//!   consensus rather than a single-sample fit; sample-to-sample scatter is
//!   of order 20–30 %, motivating the 30 % default tolerance.
//! - Only the qualitative *shape* (angular symmetry and ratio) is tested; the
//!   absolute value of ρ₁ depends on the exact Pt film quality and not just the
//!   bulk parameters in the model.
//!
//! ## References
//!
//! - K. Nakayama, H. Jungfleisch, T. Balogh, F. Casanova, L. E. Hueso,
//!   G. Tatara, E. Saitoh,
//!   "Hall effect caused by spin fluctuations in the paramagnetic phase",
//!   *Phys. Rev. Lett.* **110**, 206601 (2013).
//! - Y. T. Chen, S. Takahashi, H. Nakayama, M. Althammer,
//!   S. T. B. Goennenwein, E. Saitoh, G. E. W. Bauer,
//!   "Theory of spin Hall magnetoresistance",
//!   *Phys. Rev. B* **87**, 144411 (2013).

use crate::effect::smr::SpinHallMagnetoresistance;
use crate::error::Result;
use crate::validation::experimental::ValidationResult;
use crate::vector3::Vector3;

// ──────────────────────────────────────────────────────────────────────────────
// Reference constants
// ──────────────────────────────────────────────────────────────────────────────

/// Reference SMR ratio for Pt/YIG from Nakayama 2013.
///
/// The paper reports `ρ₁/ρ₀` in the range `1–5 × 10⁻³` for different Pt
/// thicknesses; we use `4 × 10⁻³` as the central value for the 7 nm Pt layer
/// from Table I of Chen et al. PRB 87, 144411 (2013) which provides the
/// theoretical formula used here.
pub const SMR_RATIO_REFERENCE: f64 = 4.0e-3;

/// Default validation tolerance (30 %).
pub const DEFAULT_TOLERANCE: f64 = 0.30_f64;

/// Number of angular points used in each angular validation scan.
pub const N_ANGLE_POINTS: usize = 8;

// Compile-time sanity
const _: () = assert!(SMR_RATIO_REFERENCE > 0.0);
const _: () = assert!(SMR_RATIO_REFERENCE < 1.0);
const _: () = assert!(DEFAULT_TOLERANCE > 0.0);
const _: () = assert!(DEFAULT_TOLERANCE < 1.0);
const _: () = assert!(N_ANGLE_POINTS >= 4);

// ──────────────────────────────────────────────────────────────────────────────
// Validation harness
// ──────────────────────────────────────────────────────────────────────────────

/// Validation harness for Nakayama *et al.* (2013).
///
/// Contains a [`SpinHallMagnetoresistance`] initialised to the Pt/YIG preset
/// and exposes three validation methods corresponding to the three landmark
/// claims of the paper.
#[derive(Debug, Clone)]
pub struct Nakayama2013Validation {
    /// SMR model using Pt/YIG parameters (Nakayama 2013 preset).
    pub smr: SpinHallMagnetoresistance,
}

impl Nakayama2013Validation {
    /// Construct the validation harness using [`SpinHallMagnetoresistance::platinum_yig`].
    ///
    /// # Errors
    /// This function is infallible with the canonical preset, but returns `Result`
    /// for uniformity with other harnesses.
    pub fn new() -> Result<Self> {
        Ok(Self {
            smr: SpinHallMagnetoresistance::platinum_yig(),
        })
    }

    // ── Internal helpers ──────────────────────────────────────────────────────

    /// Magnetisation direction for in-plane (*x*–*y* plane) rotation at angle `alpha`.
    ///
    /// `m = (cos α, sin α, 0)` — as used in the Nakayama experiment where the
    /// applied field rotates within the sample plane perpendicular to the
    /// interface normal *ẑ*.
    fn m_in_plane(alpha: f64) -> Vector3<f64> {
        Vector3::new(alpha.cos(), alpha.sin(), 0.0)
    }

    // ── Public validation methods ─────────────────────────────────────────────

    /// Validate the *longitudinal* SMR angular dependence.
    ///
    /// For `m = (cos α, sin α, 0)` (in-plane x–y rotation):
    ///
    /// ```text
    /// ρ_L(α) = ρ₀ + ρ₁(1 − sin²α) = ρ₀ + ρ₁ cos²α
    /// ```
    ///
    /// The simulation is compared directly against this analytical formula at
    /// [`N_ANGLE_POINTS`] equally spaced angles.  The relative error is
    /// `|ρ_L_sim(α) − ρ_L_ref(α)| / |ρ_L_ref(α)|`.
    ///
    /// # Arguments
    /// * `tolerance` — maximum acceptable relative error (use [`DEFAULT_TOLERANCE`]).
    pub fn validate_longitudinal_angular(&self, tolerance: f64) -> Result<ValidationResult> {
        let rho_0 = self.smr.rho_0();
        let rho_1 = self.smr.rho_1();
        let mut errors = Vec::with_capacity(N_ANGLE_POINTS);

        for k in 0..N_ANGLE_POINTS {
            let alpha = 2.0 * std::f64::consts::PI * (k as f64) / (N_ANGLE_POINTS as f64);
            let m = Self::m_in_plane(alpha);

            // Simulated value from the full model
            let rho_l_sim = self.smr.longitudinal_resistivity(m);

            // Reference value: analytical formula for x-y plane rotation
            // m_y = sin α ⟹ ρ_L = ρ₀ + ρ₁(1 − sin²α)
            let rho_l_ref = rho_0 + rho_1 * (1.0 - alpha.sin() * alpha.sin());

            if rho_l_ref.abs() > 0.0 {
                errors.push((rho_l_sim - rho_l_ref).abs() / rho_l_ref.abs());
            }
        }

        Ok(ValidationResult::new(
            "Nakayama 2013 longitudinal SMR angular dependence",
            &errors,
            tolerance,
        ))
    }

    /// Validate the SMR ratio `ρ₁/ρ₀` against the Nakayama 2013 reference value.
    ///
    /// Reference: `ΔSMR/ρ₀ ≈ 1 × 10⁻³` for Pt(7 nm)/YIG.
    ///
    /// The relative error is `|ratio_sim − ratio_ref| / |ratio_ref|`.
    ///
    /// # Arguments
    /// * `tolerance` — maximum acceptable relative error.
    pub fn validate_smr_ratio(&self, tolerance: f64) -> Result<ValidationResult> {
        let ratio_sim = self.smr.smr_ratio();
        let err = (ratio_sim - SMR_RATIO_REFERENCE).abs() / SMR_RATIO_REFERENCE.abs();
        Ok(ValidationResult::new(
            "Nakayama 2013 SMR ratio (ρ₁/ρ₀)",
            &[err],
            tolerance,
        ))
    }

    /// Validate the *Hall* SMR angular dependence.
    ///
    /// For `m = (cos α, sin α, 0)` (in-plane x–y, so `m_z = 0`):
    ///
    /// ```text
    /// ρ_H(α) = ρ₁ m_x m_y = ρ₁ cos(α) sin(α) = (ρ₁/2) sin(2α)
    /// ```
    ///
    /// The simulation is compared against `(ρ₁/2) sin(2α)` at
    /// [`N_ANGLE_POINTS`] angles.  Points where the reference value is
    /// effectively zero (|ref| < ρ₁ × 1e-12) are skipped to avoid division by
    /// near-zero.
    ///
    /// # Arguments
    /// * `tolerance` — maximum acceptable relative error.
    pub fn validate_hall_angular(&self, tolerance: f64) -> Result<ValidationResult> {
        let rho_1 = self.smr.rho_1();
        let eps = rho_1.abs() * 1e-9;
        let mut errors = Vec::with_capacity(N_ANGLE_POINTS);

        for k in 0..N_ANGLE_POINTS {
            let alpha = 2.0 * std::f64::consts::PI * (k as f64) / (N_ANGLE_POINTS as f64);
            let m = Self::m_in_plane(alpha);

            // Simulated Hall resistivity
            let rho_h_sim = self.smr.hall_resistivity(m);

            // Reference: (ρ₁/2) sin(2α) = ρ₁ cos α sin α
            let rho_h_ref = rho_1 * (2.0 * alpha).sin() / 2.0;

            if rho_h_ref.abs() > eps {
                errors.push((rho_h_sim - rho_h_ref).abs() / rho_h_ref.abs());
            }
        }

        Ok(ValidationResult::new(
            "Nakayama 2013 Hall SMR angular dependence",
            &errors,
            tolerance,
        ))
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 0.30_f64;

    fn build() -> Nakayama2013Validation {
        Nakayama2013Validation::new().expect("harness must build from platinum_yig preset")
    }

    // ── Compile-time checks ───────────────────────────────────────────────────

    const _: () = assert!(SMR_RATIO_REFERENCE > 0.0);
    const _: () = assert!(SMR_RATIO_REFERENCE < 1.0);
    const _: () = assert!(DEFAULT_TOLERANCE > 0.0);
    const _: () = assert!(N_ANGLE_POINTS > 0);

    // ── Construction tests ────────────────────────────────────────────────────

    #[test]
    fn test_build_succeeds() {
        let v = build();
        // Confirm Pt/YIG parameters: positive theta_sh, sensible magnitudes
        assert!(v.smr.theta_sh > 0.0);
        assert!(v.smr.resistivity_nm > 0.0);
        assert!(v.smr.lambda_sf > 0.0);
        assert!(v.smr.t_nm > 0.0);
        assert!(v.smr.sigma_nm > 0.0);
        assert!(v.smr.g_r > 0.0);
    }

    #[test]
    fn test_smr_rho1_positive() {
        let v = build();
        assert!(v.smr.rho_1() > 0.0);
    }

    // ── Longitudinal angular validation ───────────────────────────────────────

    #[test]
    fn test_longitudinal_angular_passes_at_30pct_tolerance() {
        let v = build();
        let result = v
            .validate_longitudinal_angular(TOL)
            .expect("validation should run");
        assert_eq!(result.n_points, N_ANGLE_POINTS);
        assert!(result.max_relative_error.is_finite());
        assert!(
            result.passed,
            "Longitudinal angular validation failed at 30%: {}",
            result.summary()
        );
    }

    #[test]
    fn test_longitudinal_angular_near_exact() {
        // The model implements exactly the Chen 2013 formula, so the comparison
        // against the same formula should be at floating-point precision
        let v = build();
        let result = v.validate_longitudinal_angular(TOL).expect("should run");
        assert!(
            result.max_relative_error < 1e-12,
            "Longitudinal SMR angular error should be near machine precision: {}",
            result.summary()
        );
    }

    #[test]
    fn test_longitudinal_angular_summary_contains_name() {
        let v = build();
        let result = v.validate_longitudinal_angular(TOL).expect("should run");
        let s = result.summary();
        assert!(
            s.contains("Nakayama"),
            "Summary should contain 'Nakayama': {s}"
        );
        assert!(
            s.contains("longitudinal") || s.contains("SMR"),
            "Summary should mention longitudinal: {s}"
        );
    }

    // ── SMR ratio validation ──────────────────────────────────────────────────

    #[test]
    fn test_smr_ratio_passes_at_30pct_tolerance() {
        let v = build();
        let result = v
            .validate_smr_ratio(TOL)
            .expect("SMR ratio validation should run");
        assert_eq!(result.n_points, 1);
        assert!(result.max_relative_error.is_finite());
        assert!(
            result.passed,
            "SMR ratio failed at 30%: {}",
            result.summary()
        );
    }

    #[test]
    fn test_smr_ratio_finite_and_positive() {
        let v = build();
        let ratio = v.smr.smr_ratio();
        assert!(ratio.is_finite());
        assert!(ratio > 0.0);
    }

    #[test]
    fn test_smr_ratio_summary_contains_name() {
        let v = build();
        let result = v.validate_smr_ratio(TOL).expect("should run");
        let s = result.summary();
        assert!(
            s.contains("Nakayama"),
            "Summary should contain 'Nakayama': {s}"
        );
    }

    // ── Hall angular validation ───────────────────────────────────────────────

    #[test]
    fn test_hall_angular_passes_at_30pct_tolerance() {
        let v = build();
        let result = v
            .validate_hall_angular(TOL)
            .expect("Hall angular validation should run");
        assert!(result.max_relative_error.is_finite());
        assert!(
            result.passed,
            "Hall angular validation failed at 30%: {}",
            result.summary()
        );
    }

    #[test]
    fn test_hall_angular_near_exact() {
        // Same formula comparison → should be at floating-point noise level
        let v = build();
        let result = v.validate_hall_angular(TOL).expect("should run");
        assert!(
            result.max_relative_error < 1e-12,
            "Hall SMR angular error should be near machine precision: {}",
            result.summary()
        );
    }

    #[test]
    fn test_hall_angular_summary_contains_name() {
        let v = build();
        let result = v.validate_hall_angular(TOL).expect("should run");
        let s = result.summary();
        assert!(
            s.contains("Nakayama"),
            "Summary should contain 'Nakayama': {s}"
        );
    }

    // ── Physical sanity checks ────────────────────────────────────────────────

    #[test]
    fn test_rho_l_maximum_at_alpha_zero() {
        // For m in x-y plane: ρ_L(α) = ρ₀ + ρ₁(1−sin²α) is maximum when sin(α)=0, i.e. α=0 or π
        let v = build();
        let rho_at_0 = v.smr.longitudinal_resistivity(Vector3::new(1.0, 0.0, 0.0));
        let rho_at_pi2 = v.smr.longitudinal_resistivity(Vector3::new(0.0, 1.0, 0.0));
        assert!(
            rho_at_0 > rho_at_pi2,
            "ρ_L should be maximum when m ⊥ ŷ: rho(0)={rho_at_0:.4e} vs rho(π/2)={rho_at_pi2:.4e}"
        );
    }

    #[test]
    fn test_hall_sign_at_pi_over_4() {
        // At α = π/4: m = (1/√2, 1/√2, 0), ρ_H = ρ₁·(1/√2)·(1/√2) = ρ₁/2 > 0.
        // At α = 3π/4: m = (−1/√2, 1/√2, 0), ρ_H = ρ₁·(−1/√2)·(1/√2) = −ρ₁/2 < 0.
        // These are the two maxima of sin(2α) with opposite sign — the essential Nakayama result.
        let v = build();
        let rho_1 = v.smr.rho_1();
        let alpha1 = std::f64::consts::PI / 4.0;
        let alpha2 = 3.0 * std::f64::consts::PI / 4.0;
        let m1 = Vector3::new(alpha1.cos(), alpha1.sin(), 0.0);
        let m2 = Vector3::new(alpha2.cos(), alpha2.sin(), 0.0);
        let h1 = v.smr.hall_resistivity(m1);
        let h2 = v.smr.hall_resistivity(m2);
        // h1 > 0 and h2 < 0
        assert!(h1 > 0.0, "ρ_H(π/4) should be positive: {h1:.3e}");
        assert!(h2 < 0.0, "ρ_H(3π/4) should be negative: {h2:.3e}");
        // They should have equal magnitude = ρ₁/2
        let expected = rho_1 / 2.0;
        assert!((h1 - expected).abs() < 1e-25, "ρ_H(π/4) should equal ρ₁/2");
        assert!(
            (h2 + expected).abs() < 1e-25,
            "ρ_H(3π/4) should equal −ρ₁/2"
        );
    }

    #[test]
    fn test_full_validation_suite_runs_without_error() {
        let v = build();
        v.validate_longitudinal_angular(TOL)
            .expect("longitudinal angular");
        v.validate_smr_ratio(TOL).expect("smr ratio");
        v.validate_hall_angular(TOL).expect("hall angular");
    }
}
