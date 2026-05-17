//! Validation against landmark experimental papers.
//!
//! Each submodule corresponds to a single landmark paper, with reference data
//! embedded as `const` arrays. The validation result includes maximum and mean
//! relative errors, plus a pass/fail flag based on a configurable tolerance.
//!
//! A "passing" validation does not certify that the simulation reproduces the
//! experiment to within experimental precision; it certifies that the model
//! behaves on the correct order of magnitude and exhibits the right qualitative
//! trends. The default tolerance is intentionally generous (30 %) because:
//!
//! - Experimental papers carry their own systematic uncertainties (often 10–30 %).
//! - Material parameters in the simulation are taken from generic literature
//!   values, not from the specific sample reported in the paper.
//! - Several models implemented here are simplified (e.g. magnetostatic limit,
//!   uniform mode approximations) compared to the full theory in the paper.
//!
//! To tighten the tolerance for a specific use case, supply a smaller value
//! when calling a validation method.
//!
//! ## Landmark Papers Validated
//!
//! | Submodule | Effect | Paper | Year |
//! |---|---|---|---|
//! | [`demidov_2006`] | Damon-Eshbach BLS in YIG | PRL **96**, 097202 | 2006 |
//! | [`saitoh_2006`] | Inverse Spin Hall Effect in Pt/Py | APL **88**, 182509 | 2006 |
//! | [`uchida_2008`] | Longitudinal Spin Seebeck Effect in Pt/YIG | Nature **455**, 778 | 2008 |
//! | [`mosendz_2010`] | Quantitative spin pumping in Py/Pt | PRL **104**, 046601 | 2010 |
//! | [`liu_2012`] | SOT switching with β-Ta in Ta/CoFeB/MgO | Science **336**, 555 | 2012 |
//! | [`garello_2013`] | Quantitative SOT decomposition in Pt/Co/AlO_x | Nat. Nanotechnol. **8**, 587 | 2013 |
//! | [`boona_2014`] | LSSE in granular YIG/Pt | MRS Bulletin **39**, 426 | 2014 |

pub mod boona_2014;
pub mod demidov_2006;
pub mod garello_2013;
pub mod liu_2012;
pub mod mosendz_2010;
pub mod saitoh_2006;
pub mod uchida_2008;

/// Generic result of an experimental validation: relative-error metrics + pass flag.
///
/// All errors are dimensionless relative errors `|sim − exp| / |exp|`. The
/// summary string is human-readable and suitable for direct printing in
/// `cargo test -- --nocapture` output.
#[derive(Debug, Clone)]
pub struct ValidationResult {
    /// Short human-readable name of the validation (e.g. `"Demidov 2006 dispersion"`).
    pub name: String,
    /// Maximum relative error across all comparison points.
    pub max_relative_error: f64,
    /// Arithmetic mean of relative errors over the comparison points.
    pub mean_relative_error: f64,
    /// Number of comparison points used.
    pub n_points: usize,
    /// Tolerance against which `max_relative_error` is compared.
    pub tolerance: f64,
    /// `true` iff `max_relative_error <= tolerance`.
    pub passed: bool,
}

impl ValidationResult {
    /// Build a [`ValidationResult`] from a slice of per-point relative errors.
    ///
    /// # Arguments
    /// * `name`     - Human-readable label.
    /// * `errors`   - Per-point relative errors `|sim − exp| / |exp|`.
    /// * `tolerance` - Threshold for `passed` (compared against the max).
    ///
    /// An empty `errors` slice yields `max = mean = 0.0`, `n_points = 0`, and
    /// `passed = true` (vacuously). This makes it safe to chain even when a
    /// validation is skipped for missing data.
    pub fn new(name: impl Into<String>, errors: &[f64], tolerance: f64) -> Self {
        let n_points = errors.len();
        let (max_relative_error, mean_relative_error) = if n_points == 0 {
            (0.0, 0.0)
        } else {
            let max = errors.iter().copied().fold(0.0_f64, f64::max);
            let mean = errors.iter().sum::<f64>() / (n_points as f64);
            (max, mean)
        };
        let passed = max_relative_error <= tolerance;
        Self {
            name: name.into(),
            max_relative_error,
            mean_relative_error,
            n_points,
            tolerance,
            passed,
        }
    }

    /// Produce a one-line human-readable summary, e.g.
    /// `"Demidov 2006: max=0.12, mean=0.07 over 8 pts (tol=0.30) -> PASS"`.
    pub fn summary(&self) -> String {
        format!(
            "{}: max={:.4}, mean={:.4} over {} pts (tol={:.2}) -> {}",
            self.name,
            self.max_relative_error,
            self.mean_relative_error,
            self.n_points,
            self.tolerance,
            if self.passed { "PASS" } else { "FAIL" }
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validation_result_pass() {
        let errors = vec![0.05, 0.10, 0.15, 0.08];
        let result = ValidationResult::new("test", &errors, 0.30);
        assert_eq!(result.n_points, 4);
        assert!((result.max_relative_error - 0.15).abs() < 1e-12);
        assert!((result.mean_relative_error - 0.095).abs() < 1e-12);
        assert!(result.passed);
    }

    #[test]
    fn test_validation_result_fail() {
        let errors = vec![0.05, 0.40];
        let result = ValidationResult::new("test", &errors, 0.30);
        assert!(!result.passed);
        assert!((result.max_relative_error - 0.40).abs() < 1e-12);
    }

    #[test]
    fn test_validation_result_empty() {
        let errors: Vec<f64> = vec![];
        let result = ValidationResult::new("empty", &errors, 0.30);
        assert_eq!(result.n_points, 0);
        assert_eq!(result.max_relative_error, 0.0);
        assert_eq!(result.mean_relative_error, 0.0);
        assert!(result.passed); // vacuously true
    }

    #[test]
    fn test_summary_format() {
        let errors = vec![0.10, 0.20];
        let result = ValidationResult::new("foo", &errors, 0.30);
        let s = result.summary();
        assert!(s.contains("foo"));
        assert!(s.contains("PASS"));
        assert!(s.contains("max="));
        assert!(s.contains("mean="));
    }
}
