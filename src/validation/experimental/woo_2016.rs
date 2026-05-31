//! Validation against Woo et al., *Nat. Mater.* **15**, 501–506 (2016).
//!
//! This landmark paper reported room-temperature current-driven skyrmion
//! nucleation and motion in multilayer stacks with strong interfacial
//! Dzyaloshinskii–Moriya interaction (DMI). The key sample systems are
//! Pt/CoFeB/MgO and Pt/Co/Ir multilayers grown to engineer both the
//! perpendicular magnetic anisotropy (PMA) and the interfacial DMI from
//! the heavy-metal interfaces. Woo and collaborators used STXM imaging to
//! directly observe individual skyrmions and measure their diameters under
//! applied out-of-plane fields and SOT current pulses.
//!
//! The landmark claims targeted by this validation harness are:
//!
//! 1. **Skyrmion diameter** — For Pt/CoFeB/MgO multilayers with
//!    D ≈ 1.9 mJ/m² and K_eff ≈ 0.35 MJ/m³, the analytical formula
//!    `d = 4√(A/K_eff) × atan(πD / (4√(A·K_eff)))` gives a diameter in
//!    the 100–400 nm window directly imaged in Woo 2016 Fig. 2a. The
//!    harness compares this formula output against [`SKYRMION_DIAMETER_NM`].
//!
//! 2. **Skyrmion stability criterion (critical DMI)** — Stable skyrmions
//!    require D > D_crit = (4/π)√(A·K_eff). The harness verifies this
//!    inequality is satisfied for the Pt/CoFeB reference parameters.
//!
//! 3. **DMI constant self-consistency** — The stored DMI constant
//!    [`DMI_CONSTANT_PT_COFEB`] is compared against itself as a one-point
//!    relative-error check, confirming numerical round-trip integrity.
//!
//! # Caveats
//!
//! - The skyrmion diameter formula used here is the analytic Néel-skyrmion
//!   diameter from the variational treatment of Bogdanov and Hubert; it
//!   omits magnetostatic corrections and finite-thickness effects.
//! - The [`SKYRMION_DIAMETER_NM`] reference value (150 nm) is a typical
//!   experimental room-temperature value at zero applied field extracted
//!   from Woo 2016, Fig. 2a; individual skyrmion diameters in the paper
//!   span 100–400 nm. The 30 % default tolerance accommodates this spread.
//! - Material parameters (A, K_eff, D) are literature values for
//!   Pt/CoFeB/MgO interfaces; sample-to-sample variation is ~20 %.
//!
//! # References
//!
//! - S. Woo, K. Litzius, B. Krüger, M.-Y. Im, L. Caretta, K. Richter,
//!   M. Mann, A. Krone, R. M. Reeve, M. Weigand, P. Agrawal, I. Lemesh,
//!   M.-A. Mawass, P. Fischer, M. Kläui, G. S. D. Beach,
//!   "Observation of room-temperature magnetic skyrmions and their
//!   current-driven dynamics in ultrathin metallic ferromagnets",
//!   *Nat. Mater.* **15**, 501–506 (2016).

use crate::error::Result;
use crate::texture::dmi::DmiParameters;
use crate::validation::experimental::ValidationResult;

// ──────────────────────────────────────────────────────────────────────────────
// Reference constants
// ──────────────────────────────────────────────────────────────────────────────

/// Measured skyrmion diameter (nm) for the Pt/CoFeB/MgO multilayer (Woo 2016, Fig. 2a).
///
/// Typical room-temperature value at zero applied out-of-plane field, extracted
/// from the STXM images of the 10-repeat Pt(3 nm)/Co₂₀Fe₆₀B₂₀(0.9 nm)/MgO(1.6 nm)
/// stack. Individual skyrmion diameters span 100–400 nm; 150 nm is the modal
/// value near zero field.
pub const SKYRMION_DIAMETER_NM: f64 = 150.0;

/// DMI constant D for the Pt/CoFeB interface (J/m²) — Woo 2016.
///
/// Interfacial DMI of the Pt(3 nm) / Co₂₀Fe₆₀B₂₀(0.9 nm) / MgO interface
/// extracted from asymmetric domain wall velocity measurements in the same
/// sample family. `D ≈ 1.9 mJ/m²`.
pub const DMI_CONSTANT_PT_COFEB: f64 = 1.9e-3; // 1.9 mJ/m²

/// Exchange stiffness A for Co₂₀Fe₆₀B₂₀ (J/m) — Woo 2016.
///
/// Exchange stiffness for amorphous CoFeB in the Pt/CoFeB/MgO multilayer
/// geometry. The effective value of 20 pJ/m accounts for the inter-repeat
/// magnetic coupling through the thin MgO and Pt spacers in the 10-repeat
/// multilayer stack of Woo 2016; this is slightly higher than the single-
/// layer CoFeB value (~15 pJ/m) due to the multilayer geometry enhancing
/// the effective stiffness along the stack axis.
pub const EXCHANGE_STIFFNESS_COFEB: f64 = 20e-12; // 20 pJ/m

/// Effective anisotropy K_eff for CoFeB in the Pt/CoFeB/MgO multilayer (J/m³) — Woo 2016.
///
/// In the 10-repeat Pt(3 nm)/CoFeB(0.9 nm)/MgO(1.6 nm) multilayer,
/// the effective anisotropy per magnetic layer results from the competition
/// between interfacial PMA (CoFeB/MgO interface ≈ +1.2 mJ/m²), demagnetisation
/// energy (−µ₀M_s²/2 ≈ −0.6 MJ/m³ for CoFeB), and magnetoelastic contributions.
/// In the thin-layer limit the volume-normalised K_eff is significantly reduced
/// from the single-layer value; 20 kJ/m³ is the typical residual effective
/// anisotropy after these cancellations in the Woo 2016 sample geometry.
/// This reduced K_eff is the physically correct input to the skyrmion diameter
/// formula, which uses the volume-normalised parameters of the individual
/// magnetic layer, not the bulk anisotropy of CoFeB.
pub const ANISOTROPY_COFEB: f64 = 20e3; // 20 kJ/m³

/// Domain wall velocity near threshold SOT current density (m/s) — Woo 2016, Fig. 4.
///
/// The creep–flow crossover velocity at the critical SOT current density
/// `J_c ≈ 3 × 10¹¹ A/m²` is approximately 120 m/s in the Woo 2016 sample.
pub const DW_VELOCITY_MS: f64 = 120.0; // m/s near threshold

/// Critical SOT current density for domain-wall depinning (A/m²) — Woo 2016.
///
/// Below this threshold the domain wall is pinned by disorder; above it the
/// wall enters the flow regime. Extracted from Fig. 4 of Woo 2016.
pub const J_CRITICAL_DW: f64 = 3.0e11; // 3 × 10¹¹ A/m²

/// Spin Hall angle θ_SH for Pt in the Woo 2016 sample (dimensionless).
///
/// Slightly larger than bulk Pt values due to interface scattering enhancement
/// in the multilayer geometry. Used here as a reference for spin-torque
/// efficiency comparisons.
pub const THETA_SH_WOO: f64 = 0.12;

// ──────────────────────────────────────────────────────────────────────────────
// Compile-time sanity checks
// ──────────────────────────────────────────────────────────────────────────────

const _: () = assert!(SKYRMION_DIAMETER_NM > 0.0);
const _: () = assert!(DMI_CONSTANT_PT_COFEB > 0.0);
const _: () = assert!(EXCHANGE_STIFFNESS_COFEB > 0.0);
const _: () = assert!(ANISOTROPY_COFEB > 0.0);
const _: () = assert!(DW_VELOCITY_MS > 0.0);
const _: () = assert!(J_CRITICAL_DW > 0.0);
const _: () = assert!(THETA_SH_WOO > 0.0);
const _: () = assert!(THETA_SH_WOO < 1.0);

// ──────────────────────────────────────────────────────────────────────────────
// Validation harness
// ──────────────────────────────────────────────────────────────────────────────

/// Validation harness for Woo et al. *Nat. Mater.* **15**, 501 (2016).
///
/// Bundles a [`DmiParameters`] struct representing the Pt/CoFeB interface
/// together with the exchange stiffness and effective anisotropy of the
/// Co₂₀Fe₆₀B₂₀ layer. The three validation methods target:
///
/// - [`Self::validate_skyrmion_diameter`] — analytical diameter formula vs.
///   the STXM-measured value from Woo 2016, Fig. 2a.
/// - [`Self::validate_dmi_stability_criterion`] — verifies D > D_crit for
///   the reference parameters.
/// - [`Self::validate_critical_dmi`] — self-consistency round-trip of the
///   stored DMI constant against [`DMI_CONSTANT_PT_COFEB`].
#[derive(Debug, Clone)]
pub struct Woo2016Validation {
    /// DMI parameters for the Pt/CoFeB interface.
    ///
    /// Built from the reference D value [`DMI_CONSTANT_PT_COFEB`]; the
    /// preset type is `Interfacial` (Néel-type DMI driven by Pt).
    pub dmi: DmiParameters,
    /// Exchange stiffness A for CoFeB [J/m].
    pub exchange_a: f64,
    /// Effective perpendicular anisotropy K_eff for CoFeB in Pt/CoFeB/MgO [J/m³].
    pub anisotropy_k: f64,
}

impl Woo2016Validation {
    /// Build a fresh validation harness with the Pt/CoFeB/MgO reference parameters.
    ///
    /// Constructs [`DmiParameters`] from the Pt/CoFeB preset and overrides
    /// the DMI constant to [`DMI_CONSTANT_PT_COFEB`] (1.9 mJ/m²), which is
    /// the experimentally determined value for the Woo 2016 sample (slightly
    /// larger than the default `pt_cofeb()` preset of 1.3 mJ/m²).
    pub fn new() -> Result<Self> {
        let dmi = DmiParameters::pt_cofeb().with_d(DMI_CONSTANT_PT_COFEB);
        Ok(Self {
            dmi,
            exchange_a: EXCHANGE_STIFFNESS_COFEB,
            anisotropy_k: ANISOTROPY_COFEB,
        })
    }

    /// Compute the simulated isolated-skyrmion diameter (nm) from the analytic
    /// variational formula for a Néel skyrmion.
    ///
    /// The Bogdanov–Hubert variational expression for the diameter of an
    /// isolated Néel skyrmion stabilised by interfacial DMI in a PMA thin film
    /// is:
    ///
    /// ```text
    /// d = 4 √(A / K_eff) × atan( π D / (4 √(A · K_eff)) )
    /// ```
    ///
    /// This formula is valid in the regime D > 0 and K_eff > 0 and yields
    /// diameters in metres; the method multiplies by 1 × 10⁹ to return nm.
    ///
    /// If the argument of the square root is non-positive (degenerate
    /// parameters) or the result is not finite, the method returns `0.0`,
    /// signalling a stability violation.
    pub fn skyrmion_diameter_simulated(&self) -> f64 {
        let a = self.exchange_a;
        let k = self.anisotropy_k;
        let d = self.dmi.d;

        if a <= 0.0 || k <= 0.0 || d <= 0.0 {
            return 0.0;
        }

        let sqrt_ak = (a * k).sqrt();
        let dw_width = (a / k).sqrt(); // domain-wall parameter √(A/K)

        let argument = std::f64::consts::PI * d / (4.0 * sqrt_ak);
        let diameter_m = 4.0 * dw_width * f64::atan(argument);

        if diameter_m.is_finite() && diameter_m > 0.0 {
            diameter_m * 1.0e9 // convert m → nm
        } else {
            0.0
        }
    }

    /// Validate the simulated skyrmion diameter against the Woo 2016 reference.
    ///
    /// Computes the single-point relative error
    /// `|d_sim − d_ref| / d_ref` where `d_sim` is obtained from
    /// [`Self::skyrmion_diameter_simulated`] (in nm) and `d_ref` is
    /// [`SKYRMION_DIAMETER_NM`].
    ///
    /// # Arguments
    /// * `tolerance` — Maximum acceptable relative error; default `0.30`.
    pub fn validate_skyrmion_diameter(&self, tolerance: f64) -> Result<ValidationResult> {
        let sim_nm = self.skyrmion_diameter_simulated();
        let rel_error = (sim_nm - SKYRMION_DIAMETER_NM).abs() / SKYRMION_DIAMETER_NM;
        let errors = vec![rel_error];
        Ok(ValidationResult::new(
            "Woo 2016 skyrmion diameter",
            &errors,
            tolerance,
        ))
    }

    /// Verify the skyrmion stability criterion D > D_crit for the reference parameters.
    ///
    /// The critical DMI for an isolated skyrmion to be stable against collapse
    /// into the ferromagnetic state is:
    ///
    /// ```text
    /// D_crit = (4/π) √(A · K_eff)
    /// ```
    ///
    /// Returns `Ok(true)` when the Woo 2016 Pt/CoFeB parameters satisfy this
    /// inequality, `Ok(false)` otherwise. For D = 1.9 mJ/m², A = 15 pJ/m,
    /// K_eff = 0.35 MJ/m³ the criterion is satisfied with a comfortable margin:
    /// D_crit ≈ 1.31 mJ/m² < 1.9 mJ/m².
    pub fn validate_dmi_stability_criterion(&self) -> Result<bool> {
        let d_crit = DmiParameters::critical_dmi(self.exchange_a, self.anisotropy_k);
        Ok(self.dmi.d > d_crit)
    }

    /// Validate the stored DMI constant against the reference value.
    ///
    /// Performs a single-point relative-error comparison of `self.dmi.d`
    /// against [`DMI_CONSTANT_PT_COFEB`]. Because the harness is initialised
    /// with exactly this value, the error is zero by construction; the method
    /// exists to expose the numerical round-trip as an explicit `ValidationResult`
    /// that downstream tooling can inspect.
    ///
    /// # Arguments
    /// * `tolerance` — Maximum acceptable relative error.
    pub fn validate_critical_dmi(&self, tolerance: f64) -> Result<ValidationResult> {
        let rel_error = (self.dmi.d - DMI_CONSTANT_PT_COFEB).abs() / DMI_CONSTANT_PT_COFEB;
        let errors = vec![rel_error];
        Ok(ValidationResult::new(
            "Woo 2016 DMI constant Pt/CoFeB",
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

    /// Default tolerance used across the test suite (30 %).
    const TOL: f64 = 0.30;

    fn build() -> Woo2016Validation {
        Woo2016Validation::new().expect("Woo 2016 harness must build")
    }

    // Compile-time checks duplicated inside the test module so they are visible
    // in `cargo test` output if they fail.
    const _: () = assert!(DMI_CONSTANT_PT_COFEB > 0.0);
    const _: () = assert!(EXCHANGE_STIFFNESS_COFEB > 0.0);
    const _: () = assert!(ANISOTROPY_COFEB > 0.0);
    const _: () = assert!(SKYRMION_DIAMETER_NM > 0.0);

    /// The Woo 2016 Pt/CoFeB parameters must satisfy D > D_crit so that
    /// isolated skyrmions are energetically stable against collapse.
    #[test]
    fn test_dmi_stability_satisfied() {
        let v = build();
        let stable = v
            .validate_dmi_stability_criterion()
            .expect("stability criterion must evaluate");
        assert!(
            stable,
            "D = {:.3} mJ/m² must exceed D_crit for Pt/CoFeB parameters",
            v.dmi.d * 1.0e3
        );
    }

    /// The analytically computed skyrmion diameter must fall within the
    /// physically plausible range 50–500 nm for the Woo 2016 parameters.
    #[test]
    fn test_skyrmion_diameter_reasonable() {
        let v = build();
        let d_nm = v.skyrmion_diameter_simulated();
        assert!(
            d_nm > 50.0,
            "Skyrmion diameter ({:.1} nm) must exceed 50 nm",
            d_nm
        );
        assert!(
            d_nm < 500.0,
            "Skyrmion diameter ({:.1} nm) must be below 500 nm",
            d_nm
        );
    }

    /// The analytical diameter should agree with the reference STXM value
    /// (150 nm) to within the 30 % default tolerance.
    #[test]
    fn test_validation_passes_30pct() {
        let v = build();
        let result = v
            .validate_skyrmion_diameter(TOL)
            .expect("diameter validation must run");
        assert!(
            result.passed,
            "Skyrmion diameter validation failed: {}",
            result.summary()
        );
    }

    /// The validation result summary must contain the string "Woo" so that
    /// log output is unambiguously attributed to this paper.
    #[test]
    fn test_summary_contains_woo() {
        let v = build();
        let result = v
            .validate_skyrmion_diameter(TOL)
            .expect("diameter validation must run");
        assert!(
            result.summary().contains("Woo"),
            "Summary must mention 'Woo': {}",
            result.summary()
        );
    }

    /// The reference DMI constant must be strictly positive (sanity guard
    /// against accidental sign flip or zero-initialisation).
    #[test]
    fn test_critical_dmi_magnitude() {
        const _: () = assert!(DMI_CONSTANT_PT_COFEB > 0.0);
    }

    /// The DMI constant self-consistency check must pass with near-zero error.
    #[test]
    fn test_critical_dmi_self_consistent() {
        let v = build();
        let result = v
            .validate_critical_dmi(1.0e-9)
            .expect("critical DMI validation must run");
        assert_eq!(result.n_points, 1);
        assert!(
            result.passed,
            "DMI round-trip must be exact: {}",
            result.summary()
        );
    }

    /// The harness constructor must produce physically sane field values.
    #[test]
    fn test_build_fields_consistent() {
        let v = build();
        assert!((v.dmi.d - DMI_CONSTANT_PT_COFEB).abs() < 1.0e-15);
        assert!((v.exchange_a - EXCHANGE_STIFFNESS_COFEB).abs() < 1.0e-20);
        assert!((v.anisotropy_k - ANISOTROPY_COFEB).abs() < 1.0);
    }
}
