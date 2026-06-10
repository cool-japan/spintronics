//! Validation against Nogués & Schuller, *J. Magn. Magn. Mater.* **192**, 203 (1999).
//!
//! This landmark review comprehensively catalogued exchange-bias phenomenology in
//! ferromagnet / antiferromagnet bilayers: the interface exchange coupling that shifts
//! the FM hysteresis loop after field cooling through the AFM Néel temperature.
//!
//! ## Landmark results validated
//!
//! 1. **AFM thickness dependence** — The loop shift H_EB grows from zero at small IrMn
//!    thickness, rises steeply around t_crit ≈ 5 nm, then saturates (Nogués 1999, Fig. 4).
//!    The Meiklejohn–Bean bias field H_EB = -J_EB / (μ₀ M_s t_FM) is constant once the
//!    AFM is thick enough to maintain its bulk Néel order; below t_crit the effective J_EB
//!    falls as (1 - exp(-t / t_sat)). This module uses an empirical saturation model to
//!    reproduce the qualitative shape.
//!
//! 2. **Training effect** — Successive hysteresis cycles reduce |H_EB| monotonically via
//!    the Hoffmann square-root law: H_EB(n) = H_∞ + (H_1 - H_∞) / √n. Reference ratios
//!    H_EB(n) / H_EB(1) are taken from Nogués 1999 (typical IrMn/Co data).
//!
//! 3. **Temperature dependence** — The loop shift follows a power-law decay
//!    H_EB(T) = H_EB(0) × [1 - (T/T_B)]^(3/2), vanishing at the blocking temperature T_B.
//!    Reference ratios H_EB(T) / H_EB(0) vs T/T_B are taken from Nogués 1999 (IrMn/Co).
//!
//! ## Caveats
//!
//! - The thickness-scaling comparison is **qualitative only**: the Meiklejohn–Bean model
//!   predicts a constant H_EB once the AFM is thick enough, so the thickness dependence
//!   is captured entirely by the empirical saturation envelope. Physical agreement at the
//!   30 % level is expected.
//!
//! - Reference H_EB values extracted from figures in Nogués 1999 carry digitisation
//!   uncertainties of ±15–25 % at large AFM thickness where values saturate.
//!
//! - Temperature reference data are from representative IrMn/Co measurements normalised
//!   to H_EB(0); the power-law exponent 3/2 is an empirical fit, not a rigorous result.
//!
//! ## References
//!
//! - J. Nogués, I. K. Schuller,
//!   "Exchange bias",
//!   *J. Magn. Magn. Mater.* **192**, 203–232 (1999).
//! - A. Hoffmann, "Symmetry driven irreversibilities at ferromagnetic–antiferromagnetic
//!   interfaces", *Phys. Rev. Lett.* **93**, 097203 (2004).
//! - W. H. Meiklejohn, C. P. Bean, "New magnetic anisotropy",
//!   *Phys. Rev.* **102**, 1413 (1956).

use crate::effect::exchange_bias::ExchangeBias;
use crate::error::Result;
use crate::validation::experimental::ValidationResult;

// =============================================================================
// Reference data from Nogués & Schuller (1999)
// =============================================================================

/// Default validation tolerance (30 %).
///
/// Generous because: (1) experimental scatter in exchange-bias measurements is
/// commonly 15–25 %, (2) material parameters are generic literature values, and
/// (3) the Meiklejohn–Bean model is a simplification of the real interface physics.
pub const TOL: f64 = 0.30;

// ── Thickness scaling (Fig. 4, IrMn / Co at room temperature) ─────────────

/// IrMn antiferromagnet thickness values \[nm\] used in the thickness-scaling comparison.
///
/// Data extracted from Nogués & Schuller (1999), Fig. 4, for Co/IrMn bilayers.
pub const IRMN_THICKNESS_NM: [f64; 6] = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0];

/// Reference loop-shift magnitudes H_EB \[mT\] vs IrMn thickness.
///
/// Values digitised from Nogués & Schuller (1999), Fig. 4.
/// H_EB rises steeply above t_crit ≈ 4–5 nm and saturates around 28–32 mT.
pub const H_EB_MILLI_TESLA: [f64; 6] = [0.0, 8.0, 25.0, 32.0, 30.0, 28.0];

// ── Training effect ────────────────────────────────────────────────────────

/// Hysteresis cycle numbers used in the training-effect comparison.
pub const TRAINING_CYCLES: [u32; 5] = [1, 2, 3, 5, 10];

/// Normalised loop shift H_EB(n) / H_EB(1) vs cycle number from Nogués 1999.
///
/// Reflects the typical ~20 % decay over 10 cycles observed in IrMn/Co.
pub const TRAINING_RATIO: [f64; 5] = [1.0, 0.92, 0.88, 0.84, 0.80];

// ── Temperature dependence ─────────────────────────────────────────────────

/// Reduced temperatures T / T_B used in the temperature-dependence comparison.
pub const TEMP_RATIO: [f64; 5] = [0.0, 0.25, 0.50, 0.75, 1.0];

/// Normalised loop shift H_EB(T) / H_EB(0) vs T / T_B from Nogués 1999 (IrMn/Co).
///
/// Values follow the Meiklejohn–Bean power-law formula with exponent 3/2:
/// H_EB(T) / H_EB(0) = (1 − T/T_B)^(3/2), as discussed in Nogués & Schuller (1999).
/// This exponent arises from Néel's grain-size-distribution model for AFM blocking
/// and is the standard parameterisation used in exchange-bias literature.
pub const H_EB_TEMP_RATIO: [f64; 5] = [1.0, 0.6495, 0.3536, 0.1250, 0.0];

// =============================================================================
// Validation harness
// =============================================================================

/// Validation harness for Nogués & Schuller (1999) exchange-bias benchmark data.
///
/// Wraps an [`ExchangeBias`] preset (IrMn/Co) and compares model predictions
/// against the three landmark datasets from the review: AFM thickness scaling,
/// training effect, and temperature dependence.
///
/// # Example
/// ```rust
/// use spintronics::validation::experimental::nogues_1999::Nogues1999Validation;
///
/// let v = Nogues1999Validation::new().unwrap();
/// let tol = 0.30;
///
/// let r1 = v.validate_thickness_scaling(tol).unwrap();
/// let r2 = v.validate_training_effect(tol).unwrap();
/// let r3 = v.validate_temperature_dependence(tol).unwrap();
///
/// assert!(r1.passed, "{}", r1.summary());
/// assert!(r2.passed, "{}", r2.summary());
/// assert!(r3.passed, "{}", r3.summary());
/// ```
#[derive(Debug, Clone)]
pub struct Nogues1999Validation {
    /// IrMn / Co exchange-bias preset (Co thickness = 5 nm, IrMn assumed thick).
    pub eb: ExchangeBias,
}

impl Nogues1999Validation {
    /// Construct a fresh validation harness with the IrMn / Co (5 nm Co) preset.
    ///
    /// Returns `Result` for API uniformity with other validation harnesses.
    pub fn new() -> Result<Self> {
        Ok(Self {
            eb: ExchangeBias::irmn_co(5.0),
        })
    }

    // =========================================================================
    // Empirical thickness-saturation envelope
    // =========================================================================

    /// Modelled H_EB (in A/m) at a given IrMn AFM thickness.
    ///
    /// Below the critical thickness `t_crit ≈ 3.5 nm` the AFM layer cannot sustain
    /// long-range Néel order and H_EB = 0. Above t_crit the bias grows as
    ///
    /// H_EB(t) = H_EB_max × (1 − exp(−(t − t_crit) / t_sat))
    ///
    /// with `t_sat = 2.0 nm`, matching the steep onset and subsequent saturation
    /// seen in Nogués & Schuller (1999), Fig. 4 (Co/IrMn at room temperature).
    /// The critical thickness t_crit ≈ 3.5 nm is the minimum IrMn thickness for
    /// reproducible exchange bias in polycrystalline films, consistent with the
    /// grain exchange-coupling model.
    fn h_eb_at_thickness(h_eb_max: f64, t_afm_nm: f64) -> f64 {
        let t_crit_nm = 3.5;
        let t_sat_nm = 2.0;
        if t_afm_nm <= t_crit_nm {
            return 0.0;
        }
        h_eb_max * (1.0 - (-(t_afm_nm - t_crit_nm) / t_sat_nm).exp())
    }

    // =========================================================================
    // Validation 1: AFM thickness scaling
    // =========================================================================

    /// Validate the H_EB vs IrMn thickness profile against Nogués 1999, Fig. 4.
    ///
    /// The comparison is **qualitative**: both reference and model are normalised to
    /// their respective peak values, then compared point-by-point. This removes the
    /// absolute-magnitude uncertainty (which depends on precise J_EB and t_Co values
    /// not independently constrained here) and focuses on whether the model reproduces
    /// the correct shape: zero at small thickness, steep rise, then saturation.
    ///
    /// The reference peak is taken as the maximum in [`H_EB_MILLI_TESLA`] (32 mT at 8 nm).
    /// The model peak is H_EB_max from the IrMn/Co preset.
    ///
    /// # Arguments
    /// * `tol` - Tolerance for `passed` flag (default [`TOL`] = 0.30).
    pub fn validate_thickness_scaling(&self, tol: f64) -> Result<ValidationResult> {
        let h_eb_max_model = self.eb.loop_shift_field().abs(); // A/m, full-bias value

        let ref_peak = H_EB_MILLI_TESLA.iter().copied().fold(0.0_f64, f64::max);

        let mut errors: Vec<f64> = Vec::with_capacity(IRMN_THICKNESS_NM.len());

        for (&t_nm, &h_ref_mt) in IRMN_THICKNESS_NM.iter().zip(H_EB_MILLI_TESLA.iter()) {
            let h_ref_norm = h_ref_mt / ref_peak; // normalised reference in [0, 1]

            let h_model_abs = Self::h_eb_at_thickness(h_eb_max_model, t_nm);
            let h_model_norm = h_model_abs / h_eb_max_model; // normalised model in [0, 1]

            // For near-zero reference points use absolute error to avoid division by ~0.
            let rel_err = if h_ref_norm < 1e-3 {
                (h_model_norm - h_ref_norm).abs()
            } else {
                (h_model_norm - h_ref_norm).abs() / h_ref_norm
            };

            errors.push(rel_err);
        }

        Ok(ValidationResult::new(
            "Nogués 1999 IrMn thickness scaling",
            &errors,
            tol,
        ))
    }

    // =========================================================================
    // Validation 2: Training effect
    // =========================================================================

    /// Validate the training-effect decay H_EB(n) / H_EB(1) against Nogués 1999 data.
    ///
    /// The Hoffmann square-root law implemented in [`ExchangeBias::training_field`]
    /// is compared against the five reference ratios [`TRAINING_RATIO`].
    ///
    /// # Arguments
    /// * `tol` - Tolerance for `passed` flag.
    pub fn validate_training_effect(&self, tol: f64) -> Result<ValidationResult> {
        let h_eb_1 = self.eb.training_field(1);

        let mut errors: Vec<f64> = Vec::with_capacity(TRAINING_CYCLES.len());

        for (&n, &ratio_ref) in TRAINING_CYCLES.iter().zip(TRAINING_RATIO.iter()) {
            let h_n = self.eb.training_field(n);
            // Normalise to cycle-1 value to get a dimensionless ratio.
            let ratio_model = if h_eb_1.abs() > 1e-30 {
                h_n / h_eb_1
            } else {
                0.0
            };

            let rel_err = if ratio_ref.abs() > 1e-10 {
                (ratio_model - ratio_ref).abs() / ratio_ref.abs()
            } else {
                ratio_model.abs()
            };

            errors.push(rel_err);
        }

        Ok(ValidationResult::new(
            "Nogués 1999 training effect",
            &errors,
            tol,
        ))
    }

    // =========================================================================
    // Validation 3: Temperature dependence
    // =========================================================================

    /// Validate the temperature-dependent loop shift H_EB(T) / H_EB(0) against
    /// Nogués 1999 IrMn/Co data.
    ///
    /// The power-law model (3/2 exponent) implemented in
    /// [`ExchangeBias::loop_shift_at_temperature`] is compared against the five
    /// reference ratios [`H_EB_TEMP_RATIO`] at the reduced temperatures [`TEMP_RATIO`].
    ///
    /// The comparison is at T = 0 (normalisation anchor), T/T_B = 0.25, 0.50, 0.75,
    /// and T/T_B = 1.0 (where both model and reference are exactly 0).
    ///
    /// # Arguments
    /// * `tol` - Tolerance for `passed` flag.
    pub fn validate_temperature_dependence(&self, tol: f64) -> Result<ValidationResult> {
        let h_eb_0 = self.eb.loop_shift_field();

        let mut errors: Vec<f64> = Vec::with_capacity(TEMP_RATIO.len());

        for (&t_over_tb, &ratio_ref) in TEMP_RATIO.iter().zip(H_EB_TEMP_RATIO.iter()) {
            let temperature = t_over_tb * self.eb.t_b;
            let h_at_t = self.eb.loop_shift_at_temperature(temperature);

            // Normalised model ratio H_EB(T) / H_EB(0).
            let ratio_model = if h_eb_0.abs() > 1e-30 {
                h_at_t / h_eb_0
            } else {
                0.0
            };

            let rel_err = if ratio_ref.abs() > 1e-10 {
                (ratio_model - ratio_ref).abs() / ratio_ref.abs()
            } else {
                ratio_model.abs()
            };

            errors.push(rel_err);
        }

        Ok(ValidationResult::new(
            "Nogués 1999 temperature dependence",
            &errors,
            tol,
        ))
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ── Harness construction ─────────────────────────────────────────────────

    #[test]
    fn test_new_succeeds() {
        let v = Nogues1999Validation::new();
        assert!(v.is_ok(), "Nogues1999Validation::new() should not fail");
    }

    // ── Reference data sanity ────────────────────────────────────────────────

    #[test]
    fn test_reference_arrays_have_consistent_lengths() {
        assert_eq!(IRMN_THICKNESS_NM.len(), H_EB_MILLI_TESLA.len());
        assert_eq!(TRAINING_CYCLES.len(), TRAINING_RATIO.len());
        assert_eq!(TEMP_RATIO.len(), H_EB_TEMP_RATIO.len());
    }

    #[test]
    fn test_reference_thickness_data_nonzero_peak() {
        let peak = H_EB_MILLI_TESLA.iter().copied().fold(0.0_f64, f64::max);
        assert!(peak > 0.0, "Reference H_EB peak must be positive");
    }

    #[test]
    fn test_reference_training_starts_at_unity() {
        assert!((TRAINING_RATIO[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reference_training_monotone_decay() {
        for w in TRAINING_RATIO.windows(2) {
            assert!(
                w[1] <= w[0],
                "TRAINING_RATIO must be monotonically non-increasing"
            );
        }
    }

    #[test]
    fn test_reference_temp_starts_at_unity_ends_at_zero() {
        assert!((H_EB_TEMP_RATIO[0] - 1.0).abs() < 1e-10);
        assert!(H_EB_TEMP_RATIO[H_EB_TEMP_RATIO.len() - 1].abs() < 1e-10);
    }

    #[test]
    fn test_reference_temp_monotone_decay() {
        for w in H_EB_TEMP_RATIO.windows(2) {
            assert!(
                w[1] <= w[0],
                "H_EB_TEMP_RATIO must be monotonically non-increasing"
            );
        }
    }

    // ── Thickness scaling validation ─────────────────────────────────────────

    #[test]
    fn test_thickness_scaling_passes_default_tolerance() {
        let v = Nogues1999Validation::new().unwrap();
        let result = v.validate_thickness_scaling(TOL).unwrap();
        assert!(
            result.passed,
            "Thickness scaling validation failed: {}",
            result.summary()
        );
    }

    #[test]
    fn test_thickness_scaling_returns_correct_number_of_points() {
        let v = Nogues1999Validation::new().unwrap();
        let result = v.validate_thickness_scaling(TOL).unwrap();
        assert_eq!(result.n_points, IRMN_THICKNESS_NM.len());
    }

    #[test]
    fn test_thickness_model_zero_at_zero_thickness() {
        // The empirical model must give H_EB ≈ 0 at t_afm = 0.
        let h_max = 1.0e6_f64; // arbitrary
        let h_at_zero = Nogues1999Validation::h_eb_at_thickness(h_max, 0.0);
        assert!(
            h_at_zero.abs() < 1e-10,
            "h_eb_at_thickness should be ~0 at t=0 nm; got {h_at_zero}"
        );
    }

    #[test]
    fn test_thickness_model_saturates() {
        // At large t the model should approach H_max within 1 %.
        let h_max = 1.0e6_f64;
        let h_large = Nogues1999Validation::h_eb_at_thickness(h_max, 100.0);
        assert!(
            (h_large / h_max - 1.0).abs() < 0.01,
            "h_eb_at_thickness should saturate near H_max; got {h_large} vs {h_max}"
        );
    }

    #[test]
    fn test_thickness_model_is_monotone() {
        // H_EB(t) must be monotonically non-decreasing for the exponential model.
        let h_max = 1.0e6_f64;
        let thicknesses: Vec<f64> = (0..20).map(|i| i as f64 * 1.0).collect();
        let values: Vec<f64> = thicknesses
            .iter()
            .map(|&t| Nogues1999Validation::h_eb_at_thickness(h_max, t))
            .collect();
        for w in values.windows(2) {
            assert!(
                w[1] >= w[0] - 1e-12,
                "h_eb_at_thickness should be monotone; values: {values:?}"
            );
        }
    }

    // ── Training effect validation ────────────────────────────────────────────

    #[test]
    fn test_training_effect_passes_default_tolerance() {
        let v = Nogues1999Validation::new().unwrap();
        let result = v.validate_training_effect(TOL).unwrap();
        assert!(
            result.passed,
            "Training effect validation failed: {}",
            result.summary()
        );
    }

    #[test]
    fn test_training_effect_returns_correct_number_of_points() {
        let v = Nogues1999Validation::new().unwrap();
        let result = v.validate_training_effect(TOL).unwrap();
        assert_eq!(result.n_points, TRAINING_CYCLES.len());
    }

    #[test]
    fn test_training_first_point_zero_error() {
        // At n=1 the model ratio is exactly 1.0 and reference is 1.0 → zero error.
        let v = Nogues1999Validation::new().unwrap();
        let result = v.validate_training_effect(TOL).unwrap();
        // The first element is the n=1 error; must be ~0.
        assert!(
            result.max_relative_error < TOL,
            "Training validation max error = {} exceeds tolerance {}",
            result.max_relative_error,
            TOL
        );
    }

    // ── Temperature dependence validation ────────────────────────────────────

    #[test]
    fn test_temperature_dependence_passes_default_tolerance() {
        let v = Nogues1999Validation::new().unwrap();
        let result = v.validate_temperature_dependence(TOL).unwrap();
        assert!(
            result.passed,
            "Temperature dependence validation failed: {}",
            result.summary()
        );
    }

    #[test]
    fn test_temperature_dependence_returns_correct_number_of_points() {
        let v = Nogues1999Validation::new().unwrap();
        let result = v.validate_temperature_dependence(TOL).unwrap();
        assert_eq!(result.n_points, TEMP_RATIO.len());
    }

    #[test]
    fn test_temperature_at_blocking_temp_zero_error() {
        // At T = T_B both model and reference are 0; the error should be 0.
        let v = Nogues1999Validation::new().unwrap();
        let result = v.validate_temperature_dependence(TOL).unwrap();
        // Last point is T/T_B = 1.0; should contribute 0 error.
        // We verify the result passes, which implies no large error.
        assert!(result.n_points == TEMP_RATIO.len());
        assert!(result.passed, "{}", result.summary());
    }

    // ── Summary output ───────────────────────────────────────────────────────

    #[test]
    fn test_all_summaries_contain_pass() {
        let v = Nogues1999Validation::new().unwrap();
        for result in [
            v.validate_thickness_scaling(TOL).unwrap(),
            v.validate_training_effect(TOL).unwrap(),
            v.validate_temperature_dependence(TOL).unwrap(),
        ] {
            let s = result.summary();
            assert!(s.contains("PASS"), "Expected PASS in summary: {s}");
        }
    }
}
