//! Validation against Avci *et al.*, *Nat. Phys.* **11**, 570 (2015).
//!
//! This landmark paper identified and characterised the **Unidirectional Spin
//! Hall Magnetoresistance (USMR)** in Pt(3 nm)/Co(1 nm) bilayers.  Unlike
//! ordinary SMR, the USMR signal is:
//!
//! 1. **Linear in charge current density** J.
//! 2. **Odd in magnetisation** **m** (i.e. ΔR(+m) = −ΔR(−m)).
//!
//! These two symmetry properties, together with the measured coefficient
//! η ≈ 1.7 × 10⁻¹⁶ m²/A for Pt/Co, are the hallmarks that this harness
//! validates.
//!
//! ## Landmark claims validated
//!
//! 1. **Current linearity** — `ΔR/R₀ ∝ J` (positive slope for Pt/Co).
//! 2. **Magnetisation parity** — `ΔR(+m_y) = −ΔR(−m_y)`.
//! 3. **Coefficient magnitude** — `η ≈ 1.7 × 10⁻¹⁶ m²/A`.
//!
//! ## Reference data
//!
//! The embedded reference data ([`CURRENT_DENSITY`], [`USMR_REL_CHANGE`]) are
//! derived from the linear fit in Avci 2015 (Fig. 2) with η = 1.7 × 10⁻¹⁶ m²/A
//! and m along −ŷ so that the USMR relative change is positive for positive J.
//!
//! ## Caveats
//!
//! - The sample-to-sample scatter in η reported by Avci is ~20 %; the default
//!   30 % tolerance accommodates this.
//! - The model uses generic Pt/Co parameters rather than the exact parameters
//!   of the Avci sample; absolute-value agreement at the 30 % level is
//!   sufficient to confirm the correct physics.
//!
//! ## References
//!
//! - C. O. Avci, K. Garello, A. Ghosh, M. Gabureac, S. F. Alvarado,
//!   P. Gambardella,
//!   "Unidirectional spin Hall magnetoresistance in ferromagnet/normal-metal bilayers",
//!   *Nat. Phys.* **11**, 570 (2015).
//! - C. O. Avci, K. Garello, J. Mendil, A. Ghosh, N. Blasakis, M. Gabureac,
//!   M. Trassin, M. Fiebig, P. Gambardella,
//!   "Magnetoresistance of heavy and light metal/ferromagnetic metal bilayers",
//!   *Appl. Phys. Lett.* **107**, 192405 (2015).

use crate::effect::smr::UnidirectionalSmr;
use crate::error::Result;
use crate::validation::experimental::ValidationResult;
use crate::vector3::Vector3;

// ──────────────────────────────────────────────────────────────────────────────
// Reference constants and data
// ──────────────────────────────────────────────────────────────────────────────

/// Charge current densities (A/m²) at which ΔR/R₀ is tabulated.
///
/// Span the experimentally probed linear regime of Avci 2015 (Fig. 2).
pub const CURRENT_DENSITY: &[f64] = &[5.0e10, 1.0e11, 1.5e11, 2.0e11];

/// Reference USMR relative resistance change ΔR/R₀ (dimensionless).
///
/// Derived from `ΔR/R₀ = η × J` with `η = USMR_COEFFICIENT_AVCI`
/// and `m` pointing along `−ŷ` (maximises |ΔR/R₀|, positive sign for Pt/Co).
/// Values are exactly linear in J by construction; the validation method
/// tests that the model also produces this linear scaling.
pub const USMR_REL_CHANGE: &[f64] = &[8.5e-6, 1.7e-5, 2.55e-5, 3.4e-5];

/// Central USMR coefficient η for Pt(3 nm)/Co(1 nm) from Avci 2015 [m²/A].
pub const USMR_COEFFICIENT_AVCI: f64 = 1.7e-16;

// Compile-time sanity
const _: () = assert!(CURRENT_DENSITY.len() == USMR_REL_CHANGE.len());
const _: () = assert!(CURRENT_DENSITY.len() >= 2);
const _: () = assert!(USMR_COEFFICIENT_AVCI > 0.0);

// Confirm reference data is consistent: ΔR/R₀[i] = η × J[i] × 1 (m_eff = 1)
// We cannot evaluate floating-point arithmetic in const context exactly, so
// just assert that the first entry is plausible.
const _: () = assert!(USMR_REL_CHANGE[0] > 0.0);
const _: () = assert!(USMR_REL_CHANGE[0] < 1.0e-3);

/// Default validation tolerance (30 %).
pub const DEFAULT_TOLERANCE: f64 = 0.30_f64;

const _: () = assert!(DEFAULT_TOLERANCE > 0.0);
const _: () = assert!(DEFAULT_TOLERANCE < 1.0);

// ──────────────────────────────────────────────────────────────────────────────
// Validation harness
// ──────────────────────────────────────────────────────────────────────────────

/// Validation harness for Avci *et al.* (2015).
///
/// Contains a [`UnidirectionalSmr`] initialised to the Pt/Co preset and
/// exposes three validation methods corresponding to the three landmark claims
/// of the paper.
#[derive(Debug, Clone)]
pub struct Avci2015Validation {
    /// USMR model using Pt(3 nm)/Co(1 nm) parameters.
    pub usmr: UnidirectionalSmr,
}

impl Avci2015Validation {
    /// Construct the validation harness using [`UnidirectionalSmr::platinum_cobalt`].
    ///
    /// # Errors
    /// Propagates errors from [`UnidirectionalSmr::platinum_cobalt`].
    pub fn new() -> Result<Self> {
        Ok(Self {
            usmr: UnidirectionalSmr::platinum_cobalt()?,
        })
    }

    // ── Public validation methods ─────────────────────────────────────────────

    /// Validate that the simulated USMR relative change `ΔR/R₀` scales
    /// **linearly** with current density.
    ///
    /// ## Method
    ///
    /// For each `J_k` in [`CURRENT_DENSITY`] we compute
    /// `(ΔR/R₀)_sim = usmr_relative_change(m, J_k)` with `m = −ŷ`
    /// (maximises the USMR signal in the Avci geometry where `ĵ = x̂`).
    ///
    /// The simulated curve is rescaled to match the reference at the highest
    /// current density.  The residual per-point relative error then measures
    /// *deviation from linearity*, not absolute magnitude.  Because both the
    /// model and reference are rigorously linear in J, this test should pass
    /// at floating-point precision.
    ///
    /// The absolute magnitude comparison is handled separately by
    /// [`Self::validate_coefficient_magnitude`].
    ///
    /// # Arguments
    /// * `tolerance` — maximum acceptable per-point relative error.
    pub fn validate_current_linearity(&self, tolerance: f64) -> Result<ValidationResult> {
        let n = CURRENT_DENSITY.len();
        // m = −ŷ so that (ĵ × ẑ)·m = (x̂ × ẑ)·(−ŷ) = (−ŷ)·(−ŷ) = 1 > 0
        // giving positive ΔR/R₀ for positive J (matching the Avci sign convention).
        let m = Vector3::new(0.0, -1.0, 0.0);

        let mut sim_vals = Vec::with_capacity(n);
        for &j in CURRENT_DENSITY {
            sim_vals.push(self.usmr.usmr_relative_change(m, j));
        }

        // Rescale at the highest-current anchor point
        let anchor_sim = sim_vals[n - 1];
        let anchor_ref = USMR_REL_CHANGE[n - 1];
        let scale = if anchor_sim.abs() > 0.0 {
            anchor_ref / anchor_sim
        } else {
            1.0
        };

        let mut errors = Vec::with_capacity(n);
        for (k, &sim_val) in sim_vals.iter().enumerate() {
            let ref_val = USMR_REL_CHANGE[k];
            let rescaled = sim_val * scale;
            if ref_val.abs() > 0.0 {
                errors.push((rescaled - ref_val).abs() / ref_val.abs());
            }
        }

        Ok(ValidationResult::new(
            "Avci 2015 USMR linearity in J",
            &errors,
            tolerance,
        ))
    }

    /// Validate that the USMR is **odd in magnetisation**: `ΔR(+m_y) = −ΔR(−m_y)`.
    ///
    /// ## Method
    ///
    /// For a test current density `j_test`, compute `ΔR` for `m = +ŷ` and
    /// `m = −ŷ`.  The two values should be equal in magnitude and opposite in
    /// sign.  Returns `Ok(true)` if the sign test passes, `Ok(false)` otherwise.
    ///
    /// # Arguments
    /// * `j_test` — test current density \[A/m²\] (e.g. `1e11`).
    pub fn validate_magnetization_odd(&self, j_test: f64) -> Result<bool> {
        let m_pos = Vector3::new(0.0, 1.0, 0.0);
        let m_neg = Vector3::new(0.0, -1.0, 0.0);

        let dr_pos = self.usmr.usmr_relative_change(m_pos, j_test);
        let dr_neg = self.usmr.usmr_relative_change(m_neg, j_test);

        // Odd parity: ΔR(+m) + ΔR(−m) = 0
        let parity_sum = (dr_pos + dr_neg).abs();
        // The sign of dr_pos and dr_neg should differ
        let opposite_sign = dr_pos * dr_neg < 0.0;
        Ok(opposite_sign && parity_sum < 1e-20 * (dr_pos.abs() + dr_neg.abs() + 1e-30))
    }

    /// Validate the USMR coefficient magnitude against the Avci 2015 reference.
    ///
    /// Compares `self.usmr.usmr_coefficient` with [`USMR_COEFFICIENT_AVCI`].
    /// The relative error is `|η_sim − η_ref| / η_ref`.
    ///
    /// # Arguments
    /// * `tolerance` — maximum acceptable relative error.
    pub fn validate_coefficient_magnitude(&self, tolerance: f64) -> Result<ValidationResult> {
        let eta_sim = self.usmr.usmr_coefficient;
        let eta_ref = USMR_COEFFICIENT_AVCI;
        let err = (eta_sim - eta_ref).abs() / eta_ref.abs();
        Ok(ValidationResult::new(
            "Avci 2015 USMR coefficient η (m²/A)",
            &[err],
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

    fn build() -> Avci2015Validation {
        Avci2015Validation::new().expect("Avci harness must build from platinum_cobalt preset")
    }

    // ── Compile-time checks ───────────────────────────────────────────────────

    const _: () = assert!(CURRENT_DENSITY.len() == USMR_REL_CHANGE.len());
    const _: () = assert!(CURRENT_DENSITY.len() >= 2);
    const _: () = assert!(USMR_COEFFICIENT_AVCI > 0.0);
    const _: () = assert!(DEFAULT_TOLERANCE > 0.0);

    // ── Construction ─────────────────────────────────────────────────────────

    #[test]
    fn test_build_succeeds() {
        let v = build();
        assert!(v.usmr.usmr_coefficient > 0.0);
        assert!(v.usmr.smr.theta_sh > 0.0);
        assert!(v.usmr.smr.resistivity_nm > 0.0);
        assert!(v.usmr.smr.lambda_sf > 0.0);
    }

    // ── Reference data consistency ────────────────────────────────────────────

    #[test]
    fn test_reference_data_length() {
        assert_eq!(CURRENT_DENSITY.len(), USMR_REL_CHANGE.len());
        assert_eq!(CURRENT_DENSITY.len(), 4);
    }

    #[test]
    fn test_reference_data_positive_and_monotone() {
        for &j in CURRENT_DENSITY {
            assert!(j > 0.0, "current density must be positive");
        }
        for &dr in USMR_REL_CHANGE {
            assert!(dr > 0.0, "USMR_REL_CHANGE must be positive");
        }
        for w in CURRENT_DENSITY.windows(2) {
            assert!(w[1] > w[0], "current densities must be strictly increasing");
        }
        for w in USMR_REL_CHANGE.windows(2) {
            assert!(w[1] > w[0], "USMR_REL_CHANGE must be strictly increasing");
        }
    }

    #[test]
    fn test_reference_data_linear_ratio() {
        // By construction USMR_REL_CHANGE[i] / J[i] should equal η approximately
        for k in 0..CURRENT_DENSITY.len() {
            let eta_k = USMR_REL_CHANGE[k] / CURRENT_DENSITY[k];
            let rel_err = (eta_k - USMR_COEFFICIENT_AVCI).abs() / USMR_COEFFICIENT_AVCI;
            assert!(
                rel_err < 1e-9,
                "Reference data inconsistent with USMR_COEFFICIENT_AVCI at index {k}: rel_err={rel_err:.3e}"
            );
        }
    }

    // ── Current linearity ─────────────────────────────────────────────────────

    #[test]
    fn test_current_linearity_passes_30pct() {
        let v = build();
        let result = v
            .validate_current_linearity(TOL)
            .expect("linearity validation should run");
        assert_eq!(result.n_points, CURRENT_DENSITY.len());
        assert!(result.max_relative_error.is_finite());
        assert!(
            result.passed,
            "USMR linearity failed at 30%: {}",
            result.summary()
        );
    }

    #[test]
    fn test_current_linearity_near_exact() {
        // Both model and reference are linear in J by construction;
        // after rescaling the per-point error should be at machine precision.
        let v = build();
        let result = v.validate_current_linearity(TOL).expect("should run");
        assert!(
            result.max_relative_error < 1e-12,
            "Linearity should be near machine precision: {}",
            result.summary()
        );
    }

    #[test]
    fn test_current_linearity_n_points_correct() {
        let v = build();
        let result = v.validate_current_linearity(TOL).expect("should run");
        assert_eq!(result.n_points, CURRENT_DENSITY.len());
    }

    #[test]
    fn test_linearity_summary_contains_avci() {
        let v = build();
        let result = v.validate_current_linearity(TOL).expect("should run");
        let s = result.summary();
        assert!(s.contains("Avci"), "Summary must contain 'Avci': {s}");
    }

    // ── Magnetisation odd parity ──────────────────────────────────────────────

    #[test]
    fn test_magnetization_odd_at_1e11() {
        let v = build();
        let ok = v
            .validate_magnetization_odd(1.0e11)
            .expect("odd parity check should run");
        assert!(ok, "USMR must be odd in magnetisation at J=1e11 A/m²");
    }

    #[test]
    fn test_magnetization_odd_at_various_currents() {
        let v = build();
        for &j in CURRENT_DENSITY {
            let ok = v.validate_magnetization_odd(j).expect("should run");
            assert!(ok, "USMR must be odd in m at J={j:.2e} A/m²");
        }
    }

    #[test]
    fn test_magnetization_sign_opposite() {
        // Direct check: ΔR(+ŷ) and ΔR(−ŷ) must have opposite signs
        let v = build();
        let j = 1.0e11;
        let m_pos = Vector3::new(0.0, 1.0, 0.0);
        let m_neg = Vector3::new(0.0, -1.0, 0.0);
        let dr_pos = v.usmr.usmr_relative_change(m_pos, j);
        let dr_neg = v.usmr.usmr_relative_change(m_neg, j);
        assert!(
            dr_pos * dr_neg < 0.0,
            "USMR for ±ŷ must have opposite signs: {dr_pos:.3e}, {dr_neg:.3e}"
        );
    }

    // ── Coefficient magnitude ─────────────────────────────────────────────────

    #[test]
    fn test_coefficient_magnitude_passes_30pct() {
        let v = build();
        let result = v
            .validate_coefficient_magnitude(TOL)
            .expect("coefficient validation should run");
        assert_eq!(result.n_points, 1);
        assert!(result.max_relative_error.is_finite());
        assert!(
            result.passed,
            "USMR coefficient magnitude failed at 30%: {}",
            result.summary()
        );
    }

    #[test]
    fn test_coefficient_magnitude_n_points_is_one() {
        let v = build();
        let result = v.validate_coefficient_magnitude(TOL).expect("should run");
        assert_eq!(result.n_points, 1);
    }

    #[test]
    fn test_coefficient_magnitude_summary_contains_avci() {
        let v = build();
        let result = v.validate_coefficient_magnitude(TOL).expect("should run");
        let s = result.summary();
        assert!(s.contains("Avci"), "Summary must contain 'Avci': {s}");
    }

    #[test]
    fn test_coefficient_preset_matches_reference() {
        // The platinum_cobalt preset sets usmr_coefficient = 1.7e-16 exactly
        let v = build();
        let rel_err =
            (v.usmr.usmr_coefficient - USMR_COEFFICIENT_AVCI).abs() / USMR_COEFFICIENT_AVCI;
        assert!(
            rel_err < 1e-12,
            "platinum_cobalt coefficient must equal USMR_COEFFICIENT_AVCI: rel_err={rel_err:.3e}"
        );
    }

    // ── Full suite run ────────────────────────────────────────────────────────

    #[test]
    fn test_full_validation_suite_runs_without_error() {
        let v = build();
        v.validate_current_linearity(TOL).expect("linearity");
        v.validate_magnetization_odd(1.0e11).expect("odd parity");
        v.validate_coefficient_magnitude(TOL).expect("magnitude");
    }

    // ── Physical sanity checks ────────────────────────────────────────────────

    #[test]
    fn test_usmr_positive_for_m_neg_y() {
        // For ĵ = x̂, m = −ŷ: (ĵ × ẑ) = −ŷ, (−ŷ)·(−ŷ) = 1 > 0
        // usmr_relative_change uses −m.y, so for m.y = −1: −(−1) = +1 > 0
        let v = build();
        let m = Vector3::new(0.0, -1.0, 0.0);
        let dr = v.usmr.usmr_relative_change(m, 1.0e11);
        assert!(dr > 0.0, "ΔR/R₀ must be positive for m=−ŷ, J>0: {dr:.3e}");
    }

    #[test]
    fn test_usmr_zero_for_m_along_x() {
        // m along x̂: m.y = 0 ⟹ ΔR/R₀ = η J (−m.y) = 0
        let v = build();
        let m = Vector3::new(1.0, 0.0, 0.0);
        let dr = v.usmr.usmr_relative_change(m, 1.0e11);
        assert!(dr.abs() < 1e-30, "ΔR/R₀ should be zero for m=x̂: {dr:.3e}");
    }

    #[test]
    fn test_current_scan_length_matches_n_points() {
        let v = build();
        let m = Vector3::new(0.0, -1.0, 0.0);
        let scan = v.usmr.current_scan(
            m,
            CURRENT_DENSITY[0],
            *CURRENT_DENSITY.last().unwrap(),
            CURRENT_DENSITY.len(),
        );
        assert_eq!(scan.len(), CURRENT_DENSITY.len());
    }
}
