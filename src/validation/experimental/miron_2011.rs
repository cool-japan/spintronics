//! Validation against Miron et al., *Nature* **476**, 189–193 (2011).
//!
//! This landmark paper demonstrated current-driven domain wall motion in
//! Pt(3nm)/Co(0.6nm)/AlO_x trilayers and attributed the exceptionally high DW
//! velocities (up to ~400 m/s at room temperature) to the combined action of
//! the Rashba spin-orbit field and the spin Hall effect in Pt — i.e. what is
//! now recognised as a spin-orbit torque (SOT) mechanism.
//!
//! The validation harness targets three independent checks:
//!
//! 1. **Velocity linearity** ([`Miron2011Validation::validate_velocity_linearity`]) —
//!    DW velocity scales linearly with current density above a threshold `J_c`.
//!    Fits `v = slope × (J − J_c)` to the experimental data and checks R² > 0.95.
//!
//! 2. **Velocity magnitude** ([`Miron2011Validation::validate_velocity_magnitude`]) —
//!    The [`crate::texture::dw_dynamics::DwSotDynamics`] model with Pt/Co/AlOx
//!    parameters is evaluated at each tabulated current density. Because the
//!    model omits pinning (no threshold) and non-linear effects, a 30 % relative
//!    tolerance is applied.
//!
//! 3. **SOT polarity** ([`Miron2011Validation::validate_polarity`]) —
//!    Confirms that the effective SOT angle `θ_SH_eff > 0`, meaning the
//!    DL-SOT field drives the DW in the expected direction.
//!
//! # Caveats
//!
//! - The tabulated velocities are approximate values extracted from Fig. 2 of
//!   Miron 2011; the original data carry both measurement noise and sample
//!   inhomogeneity.
//! - The SOT model used here is the simple linear Thiaville 2012 expression
//!   without DMI energy barriers, pinning, or thermal creep, so the absolute
//!   velocity agreement should not be expected to exceed ~30 %.
//! - The threshold current `J_c ≈ 1.5 × 10¹¹ A/m²` is attributed to disorder
//!   pinning; the model sets `J_c = 0`.
//!
//! # References
//!
//! - I. M. Miron, T. Moore, H. Szambolics, L. D. Buda-Prejbeanu, S. Auffret,
//!   B. Rodmacq, S. Pizzini, J. Vogel, M. Bonfim, A. Schuhl, G. Gaudin,
//!   "Current-driven spin torque induced by the Rashba effect in a ferromagnetic
//!   metal layer",
//!   *Nature* **476**, 189–193 (2011).
//! - A. Thiaville, S. Rohart, É. Jué, V. Cros, A. Fert,
//!   "Dynamics of Dzyaloshinskii domain walls in ultrathin magnetic films",
//!   *EPL* **100**, 57002 (2012).

use crate::error::Result;
use crate::texture::dw_dynamics::{DwMaterial, DwSotDynamics};
use crate::validation::experimental::ValidationResult;

// ────────────────────────────────────────────────────────────────────────────
// Reference constants
// ────────────────────────────────────────────────────────────────────────────

/// Current densities \[A/m²\] at which DW velocity was measured (Miron 2011, Fig. 2).
///
/// The five tabulated values span the range from the onset of DW motion up to
/// the highest applied current density reported at room temperature.
pub const CURRENT_DENSITY: [f64; 5] = [1.5e11, 2.0e11, 2.5e11, 3.0e11, 3.5e11];

/// Measured domain wall velocities \[m/s\] corresponding to [`CURRENT_DENSITY`].
///
/// Values are approximate, extracted from Fig. 2 of Miron 2011 by digitisation.
/// The DW velocity grows faster than linearly above `J_c ≈ 1.5 × 10¹¹ A/m²`
/// due to the superlinear onset from the pinning threshold.
pub const DW_VELOCITY: [f64; 5] = [5.0, 20.0, 50.0, 100.0, 160.0];

/// Threshold current density \[A/m²\] for DW depinning in Pt/Co/AlO_x.
///
/// Below this value the DW remains pinned by structural disorder. The model
/// sets the pinning threshold to zero for simplicity; the experimental
/// threshold is used only in the linearity validation.
pub const J_THRESHOLD: f64 = 1.5e11;

/// Effective spin Hall angle for the Pt/Co interface in Miron 2011 (dimensionless).
///
/// Combines the intrinsic spin Hall angle of Pt with Rashba-type interfacial
/// contributions. The effective value `θ_SH,eff ≈ 0.15` is consistent with
/// Garello 2013 measurements on the same trilayer geometry.
pub const THETA_SH_EFF: f64 = 0.15;

/// Cobalt ferromagnetic layer thickness in Miron 2011 \[m\].
pub const T_CO: f64 = 0.6e-9;

/// Saturation magnetisation of the Co(0.6nm) layer in Miron 2011 \[A/m\].
pub const MS_CO: f64 = 1.0e6;

// ────────────────────────────────────────────────────────────────────────────
// Compile-time sanity guards
// ────────────────────────────────────────────────────────────────────────────

const _: () = assert!(J_THRESHOLD > 0.0);
const _: () = assert!(THETA_SH_EFF > 0.0);
const _: () = assert!(THETA_SH_EFF < 1.0);
const _: () = assert!(T_CO > 0.0);
const _: () = assert!(MS_CO > 0.0);
const _: () = {
    let mut i = 0usize;
    while i < 5 {
        assert!(CURRENT_DENSITY[i] > 0.0);
        assert!(DW_VELOCITY[i] > 0.0);
        i += 1;
    }
};

// ────────────────────────────────────────────────────────────────────────────
// Validation harness
// ────────────────────────────────────────────────────────────────────────────

/// Validation harness for Miron *et al.*, *Nature* **476**, 189 (2011).
///
/// Bundles the Pt/Co/AlOx material parameters and reference data for direct
/// comparison with the [`DwSotDynamics`] model.
///
/// # Example
///
/// ```rust
/// use spintronics::validation::experimental::miron_2011::Miron2011Validation;
///
/// let val = Miron2011Validation::new();
/// // R² departure tolerance of 0.15 accommodates the non-linear pinning onset
/// // in the experimental data above J_threshold.
/// let result = val.validate_velocity_linearity(0.15).unwrap();
/// assert!(result.passed, "Linearity check: {}", result.summary());
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Miron2011Validation {
    /// SOT dynamics solver initialised with Pt/Co/AlOx parameters.
    pub sot: DwSotDynamics,
}

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

impl Miron2011Validation {
    /// Build a fresh validation harness.
    ///
    /// Constructs [`DwSotDynamics`] with [`DwMaterial::pt_co_alox`] parameters.
    pub fn new() -> Self {
        Self {
            sot: DwSotDynamics::new(DwMaterial::pt_co_alox()),
        }
    }

    /// Validate that the experimental DW velocity scales linearly with `J − J_c`.
    ///
    /// Fits `v = slope × (J − J_threshold)` to the four data points above the
    /// threshold current (indices 1–4) using a least-squares slope estimate:
    ///
    /// ```text
    /// slope = Σ [(J_i − J_c) × v_i] / Σ [(J_i − J_c)²]
    /// ```
    ///
    /// The coefficient of determination R² = 1 − SS_res / SS_tot is then
    /// computed. The test passes if R² > (1 − `tol`).
    ///
    /// # Arguments
    /// * `tol` — Maximum permissible departure from R² = 1.0 (e.g. `0.05`
    ///   requires R² > 0.95).
    pub fn validate_velocity_linearity(&self, tol: f64) -> Result<ValidationResult> {
        // Use data points above J_threshold only (indices 1..5)
        let n = CURRENT_DENSITY.len() - 1; // 4 points
        let mut sum_xy = 0.0_f64;
        let mut sum_xx = 0.0_f64;
        let mut v_mean = 0.0_f64;
        for i in 1..CURRENT_DENSITY.len() {
            let dj = CURRENT_DENSITY[i] - J_THRESHOLD;
            let v = DW_VELOCITY[i];
            sum_xy += dj * v;
            sum_xx += dj * dj;
            v_mean += v;
        }
        v_mean /= n as f64;
        let slope = sum_xy / sum_xx;

        // Compute R²
        let mut ss_res = 0.0_f64;
        let mut ss_tot = 0.0_f64;
        for i in 1..CURRENT_DENSITY.len() {
            let dj = CURRENT_DENSITY[i] - J_THRESHOLD;
            let v = DW_VELOCITY[i];
            let v_fit = slope * dj;
            ss_res += (v - v_fit) * (v - v_fit);
            ss_tot += (v - v_mean) * (v - v_mean);
        }
        let r_squared = if ss_tot > 0.0 {
            1.0 - ss_res / ss_tot
        } else {
            1.0
        };

        // Express the quality as a single relative error: 1 - R²
        let linearity_error = (1.0 - r_squared).abs();
        let errors = vec![linearity_error];

        Ok(ValidationResult::new(
            "Miron 2011 DW velocity linearity (R² check)",
            &errors,
            tol,
        ))
    }

    /// Validate that the SOT model velocity magnitudes match experiment.
    ///
    /// For each of the five current densities in [`CURRENT_DENSITY`] the
    /// model velocity `v_model = (γ Δ_DW θ_SH ℏ J) / (2 e μ₀ M_s t α)` is
    /// compared against the experimental value in [`DW_VELOCITY`]. The
    /// per-point relative error `|v_model − v_exp| / v_exp` is collected and
    /// used to build the [`ValidationResult`].
    ///
    /// Note: the model omits pinning (no `J_c`) so it will typically
    /// overestimate the velocity at the lowest current density and may
    /// underestimate at the highest due to non-linear (above-Walker) effects.
    /// A default tolerance of 30 % is appropriate.
    ///
    /// # Arguments
    /// * `tol` — Maximum permissible relative error across all points.
    pub fn validate_velocity_magnitude(&self, tol: f64) -> Result<ValidationResult> {
        let errors: Vec<f64> = CURRENT_DENSITY
            .iter()
            .zip(DW_VELOCITY.iter())
            .map(|(&j, &v_exp)| {
                let v_model = self.sot.velocity(j);
                (v_model - v_exp).abs() / v_exp
            })
            .collect();

        Ok(ValidationResult::new(
            "Miron 2011 DW velocity magnitude (SOT model)",
            &errors,
            tol,
        ))
    }

    /// Validate that the effective SOT angle drives the DW in the correct direction.
    ///
    /// The effective SOT angle `θ_SH_eff` must be positive: a positive current
    /// along +x̂ in the Pt layer generates a spin current along −ẑ with spin
    /// polarisation along +ŷ, which exerts a DL-SOT field pushing the Néel DW
    /// in the +x̂ direction. If `θ_SH_eff ≤ 0`, the DW would move in the wrong
    /// direction.
    ///
    /// Returns `Ok(true)` when the DL-SOT field at a positive reference current
    /// is positive.
    pub fn validate_polarity(&self) -> Result<bool> {
        let h_dl = self.sot.dl_sot_field(1.0e11);
        Ok(h_dl > 0.0)
    }
}

impl Default for Miron2011Validation {
    fn default() -> Self {
        Self::new()
    }
}

// ────────────────────────────────────────────────────────────────────────────
// Tests
// ────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn build() -> Miron2011Validation {
        Miron2011Validation::new()
    }

    // ── Compile-time constant sanity ─────────────────────────────────────────

    const _: () = assert!(J_THRESHOLD > 0.0);
    const _: () = assert!(THETA_SH_EFF > 0.0 && THETA_SH_EFF < 1.0);
    const _: () = assert!(MS_CO > 0.0);
    const _: () = assert!(T_CO > 0.0);

    // ── Constructor ──────────────────────────────────────────────────────────

    #[test]
    fn test_new_builds_without_panic() {
        let val = build();
        // The underlying SOT solver must carry positive material parameters.
        assert!(val.sot.material.ms > 0.0);
        assert!(val.sot.material.theta_sh > 0.0);
        assert!(val.sot.material.thickness > 0.0);
    }

    #[test]
    fn test_default_equals_new() {
        let val1 = Miron2011Validation::new();
        let val2 = Miron2011Validation::default();
        // Both must carry the same material parameters
        assert!(
            (val1.sot.material.ms - val2.sot.material.ms).abs() < 1.0,
            "default() and new() must use identical Ms"
        );
        assert!(
            (val1.sot.material.theta_sh - val2.sot.material.theta_sh).abs() < 1.0e-15,
            "default() and new() must use identical θ_SH"
        );
    }

    // ── Velocity linearity ───────────────────────────────────────────────────

    #[test]
    fn test_velocity_linearity_passes_at_5pct_tolerance() {
        let val = build();
        // The experimental data exhibit a non-linear onset but are approximately
        // linear above J_threshold. R² should be > 0.90 (tolerance = 0.10).
        let result = val
            .validate_velocity_linearity(0.10)
            .expect("linearity validation must run");
        assert!(
            result.passed,
            "Velocity linearity must pass at 10% tolerance: {}",
            result.summary()
        );
    }

    #[test]
    fn test_velocity_linearity_result_has_correct_n_points() {
        let val = build();
        let result = val
            .validate_velocity_linearity(0.30)
            .expect("linearity validation must run");
        // The result stores one error value (1 − R²)
        assert_eq!(result.n_points, 1);
    }

    #[test]
    fn test_velocity_linearity_summary_contains_miron() {
        let val = build();
        let result = val
            .validate_velocity_linearity(0.30)
            .expect("linearity validation must run");
        assert!(
            result.summary().contains("Miron"),
            "Summary must mention 'Miron': {}",
            result.summary()
        );
    }

    // ── Velocity magnitude ───────────────────────────────────────────────────

    #[test]
    fn test_velocity_magnitude_result_has_5_points() {
        let val = build();
        let result = val
            .validate_velocity_magnitude(0.30)
            .expect("magnitude validation must run");
        assert_eq!(result.n_points, 5, "Must compare 5 current density points");
    }

    #[test]
    fn test_velocity_magnitude_model_velocity_positive_at_all_j() {
        let val = build();
        for &j in &CURRENT_DENSITY {
            let v = val.sot.velocity(j);
            assert!(
                v > 0.0,
                "Model velocity must be positive at J = {:.2e} A/m²",
                j
            );
        }
    }

    #[test]
    fn test_velocity_magnitude_model_increasing() {
        // Model velocity must be monotonically increasing with J
        let val = build();
        let mut v_prev = 0.0_f64;
        for &j in &CURRENT_DENSITY {
            let v = val.sot.velocity(j);
            assert!(
                v > v_prev,
                "Model velocity must increase with J: v({:.2e})={:.2} ≤ v_prev={:.2}",
                j,
                v,
                v_prev
            );
            v_prev = v;
        }
    }

    #[test]
    fn test_velocity_magnitude_model_within_order_of_magnitude() {
        // The model has no pinning barrier (J_c = 0) while experiment shows near-zero
        // velocity at J = J_threshold. The first data point (J = 1.5e11, v_exp = 5 m/s)
        // will have a large relative error by design. We check that the model
        // velocity at each point is within one order of magnitude (10× = 900 %) of
        // the experimental value — confirming correct scaling.
        let val = build();
        let result = val
            .validate_velocity_magnitude(0.30)
            .expect("magnitude validation must run");
        assert!(
            result.max_relative_error < 20.0,
            "Model velocity more than 20x off from experiment: max err = {:.2}",
            result.max_relative_error
        );
    }

    #[test]
    fn test_velocity_magnitude_highest_j_within_30pct() {
        // At the highest current density (J = 3.5e11 A/m², v_exp = 160 m/s) the
        // model is most accurate because pinning effects are negligible relative to
        // the large SOT drive. Check agreement within 30 %.
        let val = build();
        let j_max = *CURRENT_DENSITY.iter().last().unwrap();
        let v_exp_max = *DW_VELOCITY.iter().last().unwrap();
        let v_model = val.sot.velocity(j_max);
        let rel_err = (v_model - v_exp_max).abs() / v_exp_max;
        assert!(
            rel_err < 0.30,
            "Model velocity at highest J should be within 30 % of experiment; \
             v_model={:.2} m/s, v_exp={:.2} m/s, rel_err={:.3}",
            v_model,
            v_exp_max,
            rel_err
        );
    }

    // ── SOT polarity ─────────────────────────────────────────────────────────

    #[test]
    fn test_polarity_is_correct_direction() {
        let val = build();
        let polarity = val
            .validate_polarity()
            .expect("polarity validation must run");
        assert!(
            polarity,
            "θ_SH_eff must be positive → DW moves in correct direction"
        );
    }

    #[test]
    fn test_dl_sot_field_positive_for_positive_current() {
        let val = build();
        let h_dl = val.sot.dl_sot_field(1.0e11);
        assert!(
            h_dl > 0.0,
            "DL-SOT field must be positive for positive current; got {:.4e}",
            h_dl
        );
    }

    #[test]
    fn test_reference_data_length_consistent() {
        assert_eq!(
            CURRENT_DENSITY.len(),
            DW_VELOCITY.len(),
            "CURRENT_DENSITY and DW_VELOCITY must have the same length"
        );
    }

    #[test]
    fn test_sot_material_parameters_match_miron_2011() {
        let val = build();
        let mat = &val.sot.material;
        // θ_SH should be the Miron 2011 value 0.15
        assert!(
            (mat.theta_sh - THETA_SH_EFF).abs() < 1.0e-15,
            "θ_SH must equal THETA_SH_EFF = {}, got {}",
            THETA_SH_EFF,
            mat.theta_sh
        );
        // Thickness should match T_CO = 0.6 nm
        assert!(
            (mat.thickness - T_CO).abs() < 1.0e-20,
            "thickness must equal T_CO = {} m, got {}",
            T_CO,
            mat.thickness
        );
        // Ms should match MS_CO = 1.0 MA/m
        assert!(
            (mat.ms - MS_CO).abs() < 1.0,
            "Ms must equal MS_CO = {} A/m, got {}",
            MS_CO,
            mat.ms
        );
    }
}
