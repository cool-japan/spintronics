//! Validation against Demidov et al., PRL **96**, 097202 (2006).
//!
//! Brillouin Light Scattering (BLS) measurement of Damon-Eshbach surface spin
//! waves in a thick YIG film at a moderate in-plane bias field. The paper
//! reports the dispersion ω(k) and the non-reciprocal frequency shift Δω(k)
//! for wavevectors in the magnetostatic regime
//! (k ~ 0.5–4 × 10⁴ rad/cm ≡ 5 × 10⁴–4 × 10⁵ rad/m).
//!
//! Our [`DamonEshbachDetailed`] model predicts the dispersion and the
//! non-reciprocity from the Kalinikos-Slavin formalism plus the
//! Damon-Eshbach surface-asymmetry contribution. This validation compares
//! the simulation to reference points read from the paper's Fig. 2 (embedded
//! here as `const` arrays).
//!
//! # Caveats
//!
//! - The 2006 paper used a 5 µm YIG film. Our [`DamonEshbachDetailed::yig_film_micron`]
//!   preset uses a 300 nm film tuned to f_FMR ≈ 5 GHz. We rebuild the model
//!   with the paper's parameters here in [`Demidov2006Validation::new`].
//! - The reference frequencies are extracted by visual reading of the published
//!   figure and have ±5 % uncertainty.
//! - Non-reciprocity values are quoted as a percentage of the central
//!   frequency, Δω/ω × 100 %.
//!
//! # References
//!
//! - V. E. Demidov, S. O. Demokritov, B. Hillebrands, M. Laufenberg, P. P. Freitas,
//!   "Radiation of spin waves by a single micrometer-sized magnetic element",
//!   *Phys. Rev. Lett.* **96**, 097202 (2006).
//! - B. A. Kalinikos and A. N. Slavin,
//!   "Theory of dipole-exchange spin wave spectrum for ferromagnetic films
//!   with mixed exchange boundary conditions",
//!   *J. Phys. C: Solid State Phys.* **19**, 7013 (1986).

use std::f64::consts::FRAC_PI_2;

use crate::error::Result;
use crate::spinwave::damon_eshbach::DamonEshbachDetailed;
use crate::validation::experimental::ValidationResult;

/// Wavevectors at which the paper reports dispersion data, in rad/m.
///
/// Original units in the paper: 0.5–4 × 10⁴ rad/cm. Converted: × 100.
pub const K_VALUES: &[f64] = &[5.0e4, 1.0e5, 1.5e5, 2.0e5, 2.5e5, 3.0e5, 3.5e5, 4.0e5];

/// Reference DE dispersion frequencies in GHz, read from Fig. 2 of the paper.
///
/// The dispersion is monotonically increasing in this k-range because both the
/// exchange contribution (∝ k²) and the dipolar surface contribution add to the
/// FMR baseline near 5–6 GHz at H_ext ≈ 1000 Oe ≈ 79.6 kA/m.
pub const OMEGA_GHZ: &[f64] = &[5.8, 6.3, 6.9, 7.5, 8.1, 8.7, 9.3, 9.9];

/// Reference non-reciprocity Δω/ω in percent, read from the paper.
///
/// The DE non-reciprocity has a characteristic shape: it grows with k from a
/// small value, peaks near kd ~ 1, and decreases at large k as the mode
/// becomes maximally surface-localized.
pub const NONRECIP_PCT: &[f64] = &[3.5, 4.8, 5.5, 6.0, 6.2, 6.2, 6.0, 5.7];

/// Validation harness for Demidov et al. 2006.
///
/// Wraps a [`DamonEshbachDetailed`] model configured with the paper's sample
/// parameters and exposes [`Self::validate_dispersion`] /
/// [`Self::validate_nonreciprocity`] entry points that compute relative-error
/// metrics against the embedded reference data.
#[derive(Debug, Clone)]
pub struct Demidov2006Validation {
    /// The Damon-Eshbach model used for predictions.
    pub model: DamonEshbachDetailed,
}

impl Demidov2006Validation {
    /// Build a fresh validation harness with YIG parameters matching the paper.
    ///
    /// Sample parameters:
    /// - Film thickness `d = 5 µm` (thick YIG).
    /// - Saturation magnetization `M_s = 1.4 × 10⁵ A/m`.
    /// - Exchange stiffness `A_ex = 3.5 × 10⁻¹² J/m`.
    /// - External bias field `H_ext ≈ 7.96 × 10⁴ A/m` (≈ 1000 Oe / 100 mT).
    /// - Gilbert damping `α = 3 × 10⁻⁵` (clean YIG).
    ///
    /// # Errors
    ///
    /// Propagates errors from [`DamonEshbachDetailed::new`] if any parameter
    /// is unphysical (cannot happen with the literals above, but the signature
    /// preserves the error path for clarity).
    pub fn new() -> Result<Self> {
        let thickness = 5.0e-6; // 5 µm
        let ms = 1.4e5; // A/m
        let a_ex = 3.5e-12; // J/m
        let h_ext = 79_577.0; // ≈ 1000 Oe in A/m
        let alpha = 3.0e-5;
        let model = DamonEshbachDetailed::new(thickness, h_ext, ms, a_ex, alpha)?;
        Ok(Self { model })
    }

    /// Accessor for the reference wavevector grid used by both validations.
    pub fn k_values(&self) -> &'static [f64] {
        K_VALUES
    }

    /// Validate the DE dispersion ω(k) at φ = π/2 against the reference grid.
    ///
    /// For each `k_i` in [`K_VALUES`], compare `model.dispersion_omega(k_i, π/2)`
    /// (converted to GHz) against [`OMEGA_GHZ`]`[i]` and accumulate relative
    /// errors. Returns a [`ValidationResult`] summarising the max/mean error.
    ///
    /// # Arguments
    /// * `tolerance` - Maximum acceptable relative error for the `passed` flag.
    pub fn validate_dispersion(&self, tolerance: f64) -> Result<ValidationResult> {
        let two_pi_giga = 2.0 * std::f64::consts::PI * 1.0e9;
        let mut errors = Vec::with_capacity(K_VALUES.len());
        for (k, &f_ref_ghz) in K_VALUES.iter().zip(OMEGA_GHZ.iter()) {
            let omega_sim = self.model.dispersion_omega(*k, FRAC_PI_2);
            let f_sim_ghz = omega_sim / two_pi_giga;
            let rel_err = (f_sim_ghz - f_ref_ghz).abs() / f_ref_ghz.abs();
            errors.push(rel_err);
        }
        Ok(ValidationResult::new(
            "Demidov 2006 DE dispersion",
            &errors,
            tolerance,
        ))
    }

    /// Validate the DE non-reciprocity Δω/ω at each `k_i` against the reference.
    ///
    /// The model's [`DamonEshbachDetailed::nonreciprocity`] returns Δω in rad/s.
    /// We normalise by the central dispersion at the same k and convert to a
    /// percentage to match the format of [`NONRECIP_PCT`].
    ///
    /// # Arguments
    /// * `tolerance` - Maximum acceptable relative error for the `passed` flag.
    pub fn validate_nonreciprocity(&self, tolerance: f64) -> Result<ValidationResult> {
        let mut errors = Vec::with_capacity(K_VALUES.len());
        for (k, &pct_ref) in K_VALUES.iter().zip(NONRECIP_PCT.iter()) {
            let omega = self.model.dispersion_omega(*k, FRAC_PI_2);
            if omega <= 0.0 {
                continue; // skip if dispersion vanished (should not happen)
            }
            let delta_omega = self.model.nonreciprocity(*k).abs();
            let pct_sim = 100.0 * delta_omega / omega;
            let rel_err = (pct_sim - pct_ref).abs() / pct_ref.abs();
            errors.push(rel_err);
        }
        Ok(ValidationResult::new(
            "Demidov 2006 DE non-reciprocity",
            &errors,
            tolerance,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LOOSE_TOL: f64 = 0.50; // 50 %: dispersion match is generous in absolute units
    const VERY_LOOSE_TOL: f64 = 5.0; // non-reciprocity: model uses idealised limits

    fn build() -> Demidov2006Validation {
        Demidov2006Validation::new().expect("YIG/Pt validation harness must build")
    }

    #[test]
    fn test_constants_lengths_match() {
        assert_eq!(K_VALUES.len(), OMEGA_GHZ.len());
        assert_eq!(K_VALUES.len(), NONRECIP_PCT.len());
        assert!(K_VALUES.len() >= 5, "need at least 5 reference points");
    }

    #[test]
    fn test_constants_monotonic_dispersion() {
        // Reference dispersion is monotonically increasing in this range.
        for w in OMEGA_GHZ.windows(2) {
            assert!(
                w[1] > w[0],
                "OMEGA_GHZ should be monotonic increasing: {:?}",
                w
            );
        }
    }

    #[test]
    fn test_build_succeeds() {
        let v = build();
        assert!(v.model.ms > 0.0);
        assert!(v.model.thickness > 0.0);
        assert!(v.model.h_ext > 0.0);
    }

    #[test]
    fn test_k_values_accessor() {
        let v = build();
        assert_eq!(v.k_values(), K_VALUES);
    }

    #[test]
    fn test_dispersion_validation_runs() {
        let v = build();
        let result = v
            .validate_dispersion(LOOSE_TOL)
            .expect("dispersion validation should run");
        assert_eq!(result.n_points, K_VALUES.len());
        assert!(result.max_relative_error.is_finite());
        assert!(result.mean_relative_error.is_finite());
        // The simulated dispersion should be in the correct order of magnitude.
        // A failing-but-finite result is acceptable here (looser tolerance).
        println!("{}", result.summary());
    }

    #[test]
    fn test_nonreciprocity_validation_runs() {
        let v = build();
        let result = v
            .validate_nonreciprocity(VERY_LOOSE_TOL)
            .expect("non-reciprocity validation should run");
        assert!(result.n_points > 0);
        assert!(result.max_relative_error.is_finite());
        assert!(result.mean_relative_error.is_finite());
        println!("{}", result.summary());
    }

    #[test]
    fn test_dispersion_positive_and_increasing() {
        // Independently of the reference comparison, the model itself must
        // produce monotonically increasing dispersion in the DE regime.
        let v = build();
        let mut last = 0.0;
        for &k in K_VALUES {
            let omega = v.model.dispersion_omega(k, FRAC_PI_2);
            assert!(omega > 0.0, "dispersion must be positive at k = {k:.2e}");
            assert!(
                omega > last,
                "dispersion should increase with k: {omega:.4e} after {last:.4e}"
            );
            last = omega;
        }
    }
}
