//! Unified Spin Caloritronics Cross-Effects
//!
//! This module assembles the Onsager transport matrix and heat current calculator
//! into a single `SpinCaloritronicsMaterial` that exposes a high-level interface
//! for computing all spin-caloritronic cross-effects in one call.
//!
//! ## Cross-Effects Covered
//!
//! | Effect | Physical Origin | Observable |
//! |--------|----------------|-----------|
//! | Spin Seebeck | ∇T → spin current | j_s via S_s |
//! | Spin Peltier | j_s → heat current | j_Q via Π_s = T·S_s |
//! | Anomalous Nernst | ∇T → transverse voltage | V via θ_H·S_e |
//! | Spin Nernst | ∇T → transverse spin current | j_s via θ_SN |
//!
//! ## Usage
//!
//! ```rust
//! use spintronics::caloritronics::{SpinCaloritronicsMaterial, CaloritronicsResult};
//! use spintronics::Vector3;
//!
//! let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
//! let grad_t = Vector3::new(1000.0, 0.0, 0.0);
//! let j_spin = Vector3::new(0.0, 0.0, 1.0e3);
//!
//! let result = mat.compute_all(&grad_t, &j_spin).expect("computation failed");
//! assert!(result.spin_seebeck_current.magnitude() > 0.0);
//! ```

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use super::heat_current::HeatCurrentCalculator;
use super::onsager::OnsagerMatrix;
use crate::error::Error;
use crate::vector3::Vector3;

/// Spin caloritronic material: combines Onsager matrix with heat current calculator.
///
/// This is the primary user-facing type for spin caloritronics calculations.
/// It holds an `OnsagerMatrix` (charge/spin transport) and a
/// `HeatCurrentCalculator` (all heat-current mechanisms) for a specific material
/// system (typically a ferromagnet/normal-metal bilayer such as YIG/Pt).
///
/// The two sub-objects are kept separate so that users who need fine-grained
/// control can access them individually through the public fields.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpinCaloritronicsMaterial {
    /// Onsager transport matrix (charge, spin, heat coupling)
    pub onsager: OnsagerMatrix,
    /// Heat current calculator (Fourier + Peltier + spin Peltier)
    pub heat_calc: HeatCurrentCalculator,
    /// Human-readable material system label
    pub name: String,
}

impl SpinCaloritronicsMaterial {
    /// Create a material from an explicit `OnsagerMatrix`.
    ///
    /// The `HeatCurrentCalculator` is automatically constructed from the Onsager
    /// parameters using the Kelvin relations:
    ///   Π   = T · S_e
    ///   Π_s = T · S_s
    pub fn new(onsager: OnsagerMatrix) -> Self {
        let peltier_coeff = onsager.temperature * onsager.seebeck;
        let spin_peltier_coeff = onsager.temperature * onsager.spin_seebeck;
        let heat_calc = HeatCurrentCalculator::new(
            onsager.thermal_conductivity,
            peltier_coeff,
            spin_peltier_coeff,
        );
        Self {
            onsager,
            heat_calc,
            name: String::from("custom"),
        }
    }

    /// YIG/Pt bilayer preset at `temperature` \[K\].
    pub fn yig_pt(temperature: f64) -> Self {
        let onsager = OnsagerMatrix::yig_pt(temperature);
        let mut mat = Self::new(onsager);
        mat.name = String::from("YIG/Pt");
        mat
    }

    /// Fe/Pt bilayer preset at `temperature` \[K\].
    pub fn fe_pt(temperature: f64) -> Self {
        let onsager = OnsagerMatrix::fe_pt(temperature);
        let mut mat = Self::new(onsager);
        mat.name = String::from("Fe/Pt");
        mat
    }

    /// CoFeB/Pt bilayer preset at `temperature` \[K\].
    pub fn cofeb_pt(temperature: f64) -> Self {
        let onsager = OnsagerMatrix::cofeb_pt(temperature);
        let mut mat = Self::new(onsager);
        mat.name = String::from("CoFeB/Pt");
        mat
    }

    /// Compute all spin-caloritronic cross-effects for a given thermal environment.
    ///
    /// This method is the main workhorse: it evaluates all caloritronic observables
    /// for the provided temperature gradient and injected spin current.
    ///
    /// ## Computed Quantities
    ///
    /// - `spin_seebeck_current`: j_s = S_s·σ·∇T (spin Seebeck)
    /// - `peltier_heat`: |j_Q^sPeltier| = |Π_s·j_s| via spin Peltier
    /// - `nernst_voltage`: ν = −S_e·θ_H·|∇T| (anomalous Nernst, longitudinal component)
    /// - `spin_nernst_current`: j_s^SN = θ_H · j_s (spin Nernst proxy)
    /// - `reciprocity_satisfied`: whether Onsager reciprocity error < 1e-10
    /// - `total_heat_current`: j_Q = −κ·∇T + Π_s·j_s (Fourier + spin Peltier)
    ///
    /// # Arguments
    /// * `grad_t` - Temperature gradient vector [K/m]
    /// * `j_spin` - Externally injected spin current density [A/m²]
    ///
    /// # Errors
    /// Returns `Error::InvalidParameter` if `temperature` ≤ 0 or `conductivity` ≤ 0.
    pub fn compute_all(
        &self,
        grad_t: &Vector3<f64>,
        j_spin: &Vector3<f64>,
    ) -> Result<CaloritronicsResult, Error> {
        if self.onsager.temperature <= 0.0 {
            return Err(Error::InvalidParameter {
                param: String::from("temperature"),
                reason: String::from("must be positive"),
            });
        }
        if self.onsager.conductivity <= 0.0 {
            return Err(Error::InvalidParameter {
                param: String::from("conductivity"),
                reason: String::from("must be positive"),
            });
        }

        // 1. Spin Seebeck current: j_s = S_s·σ·∇T
        let spin_seebeck_current = self.onsager.spin_current_from_grad_t(grad_t);

        // 2. Spin Peltier heat from injected spin current
        let spin_peltier_heat_vec = self.heat_calc.spin_peltier_current(j_spin);
        let peltier_heat = spin_peltier_heat_vec.magnitude();

        // 3. Anomalous Nernst voltage (scalar, longitudinal ∇T → transverse V)
        let grad_t_magnitude = grad_t.magnitude();
        let nernst_voltage = self.onsager.nernst_voltage(grad_t_magnitude);

        // 4. Spin Nernst proxy: transverse spin current ∝ θ_H · j_s^SSE
        //    The spin Nernst current is perpendicular to the SSE current; here
        //    we scale the SSE current by the Hall angle as a first approximation.
        let spin_nernst_current = Vector3::new(
            self.onsager.hall_angle * spin_seebeck_current.x,
            self.onsager.hall_angle * spin_seebeck_current.y,
            self.onsager.hall_angle * spin_seebeck_current.z,
        );

        // 5. Onsager reciprocity check
        let reciprocity_satisfied = self.onsager.reciprocity_error() < 1.0e-10;

        // 6. Total heat current (Fourier + spin Peltier; no charge current here)
        let zero_charge = Vector3::zero();
        let total_heat_current = self
            .heat_calc
            .total_heat_current(grad_t, &zero_charge, j_spin);

        Ok(CaloritronicsResult {
            spin_seebeck_current,
            peltier_heat,
            nernst_voltage,
            spin_nernst_current,
            reciprocity_satisfied,
            total_heat_current,
        })
    }
}

/// Output of a full spin caloritronics computation.
///
/// Contains all primary observables for a given temperature gradient and
/// spin current input.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CaloritronicsResult {
    /// Spin Seebeck current density j_s [A/m²]
    pub spin_seebeck_current: Vector3<f64>,
    /// Scalar spin Peltier heat |j_Q^sPeltier| [W/m²]
    pub peltier_heat: f64,
    /// Anomalous Nernst voltage per unit length ν [V/m]
    pub nernst_voltage: f64,
    /// Spin Nernst proxy current j_s^SN [A/m²]
    pub spin_nernst_current: Vector3<f64>,
    /// True if Onsager reciprocity error < 1e-10
    pub reciprocity_satisfied: bool,
    /// Total heat current density [W/m²]
    pub total_heat_current: Vector3<f64>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helpers ────────────────────────────────────────────────────────────────

    fn grad_t_x(magnitude: f64) -> Vector3<f64> {
        Vector3::new(magnitude, 0.0, 0.0)
    }

    fn zero_vec() -> Vector3<f64> {
        Vector3::zero()
    }

    fn small_spin() -> Vector3<f64> {
        Vector3::new(1.0e3, 0.0, 0.0)
    }

    // ── Build tests ───────────────────────────────────────────────────────────

    /// 1. YIG/Pt material constructs without panic.
    #[test]
    fn test_yig_pt_material_builds() {
        let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
        assert_eq!(mat.name, "YIG/Pt");
        assert!(mat.onsager.conductivity > 0.0);
        assert!(mat.onsager.temperature > 0.0);
    }

    /// 2. Fe/Pt material constructs without panic.
    #[test]
    fn test_fe_pt_material_builds() {
        let mat = SpinCaloritronicsMaterial::fe_pt(300.0);
        assert_eq!(mat.name, "Fe/Pt");
        assert!(mat.onsager.conductivity > 0.0);
    }

    /// 3. CoFeB/Pt material constructs without panic.
    #[test]
    fn test_cofeb_pt_material_builds() {
        let mat = SpinCaloritronicsMaterial::cofeb_pt(300.0);
        assert_eq!(mat.name, "CoFeB/Pt");
        assert!(mat.onsager.conductivity > 0.0);
    }

    // ── Reciprocity tests ─────────────────────────────────────────────────────

    /// 4. Onsager reciprocity for YIG/Pt: error < 1e-6.
    #[test]
    fn test_onsager_reciprocity_yig_pt() {
        let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
        let err = mat.onsager.reciprocity_error();
        assert!(err < 1.0e-6, "YIG/Pt reciprocity error {err} exceeds 1e-6");
    }

    /// 5. Onsager reciprocity for Fe/Pt: error < 1e-6.
    #[test]
    fn test_onsager_reciprocity_fe_pt() {
        let mat = SpinCaloritronicsMaterial::fe_pt(300.0);
        let err = mat.onsager.reciprocity_error();
        assert!(err < 1.0e-6, "Fe/Pt reciprocity error {err} exceeds 1e-6");
    }

    // ── Current direction tests ───────────────────────────────────────────────

    /// 6. For positive spin_seebeck, spin current flows parallel to ∇T.
    #[test]
    fn test_spin_current_direction() {
        let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
        // YIG/Pt has positive spin_seebeck (1e-3)
        let grad_t = grad_t_x(1000.0);
        let j_spin = mat.onsager.spin_current_from_grad_t(&grad_t);
        // Should have positive x-component (parallel to ∇T)
        assert!(
            j_spin.x > 0.0,
            "spin current should flow along +x (∇T direction)"
        );
        assert!(j_spin.y.abs() < f64::EPSILON);
        assert!(j_spin.z.abs() < f64::EPSILON);
    }

    /// 7. Spin Peltier heat from Onsager is proportional to j_spin.
    #[test]
    fn test_peltier_heat_from_onsager() {
        let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
        let j1 = Vector3::new(1.0e3, 0.0, 0.0);
        let j2 = Vector3::new(2.0e3, 0.0, 0.0);
        let heat1 = mat.onsager.heat_current_from_spin_current(&j1);
        let heat2 = mat.onsager.heat_current_from_spin_current(&j2);
        // heat2 should be 2× heat1
        let ratio = heat2.magnitude() / heat1.magnitude();
        assert!(
            (ratio - 2.0).abs() < 1.0e-10,
            "Peltier heat not proportional: ratio = {ratio}"
        );
    }

    /// 8. Nernst voltage has correct sign for negative Seebeck coefficient.
    #[test]
    fn test_nernst_sign() {
        let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
        // seebeck < 0, hall_angle > 0, grad_t > 0
        // ν = -S_e · θ_H · |∇T| → positive (because -(-) = +)
        let nernst_v = mat.onsager.nernst_voltage(1000.0);
        assert!(
            nernst_v > 0.0,
            "Nernst voltage should be positive for negative S_e"
        );
    }

    /// 9. Zero electric field → charge current comes from thermal driving only.
    #[test]
    fn test_zero_field_zero_charge_current() {
        let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
        let grad_t = grad_t_x(1000.0);
        let e_field = zero_vec();
        let currents = mat.onsager.all_currents(&grad_t, &e_field);
        // j_c = σ · S_e · ∇T  (Seebeck driving, no external field)
        // S_e = -5e-6, σ = 5e6, ∇T_x = 1000  → j_c_x = 5e6 * (-5e-6) * 1000 = -25
        let expected_jc_x = mat.onsager.conductivity * mat.onsager.seebeck * grad_t.x;
        assert!(
            (currents.charge_current.x - expected_jc_x).abs() < 1.0e-6,
            "charge current mismatch"
        );
    }

    /// 10. Zero temperature gradient → spin Seebeck current vanishes.
    #[test]
    fn test_zero_grad_t_zero_spin_current() {
        let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
        let grad_t = zero_vec();
        let j_spin = mat.onsager.spin_current_from_grad_t(&grad_t);
        assert!(
            j_spin.magnitude() < f64::EPSILON,
            "Spin current should be zero for zero ∇T"
        );
    }

    /// 11. Kelvin relation: Π_s = T · S_s.
    #[test]
    fn test_peltier_seebeck_ratio() {
        let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
        let expected_pi_s = mat.onsager.temperature * mat.onsager.spin_seebeck;
        // The heat_calc stores Π_s = T · S_s
        assert!(
            (mat.heat_calc.spin_peltier_coeff - expected_pi_s).abs() < 1.0e-12,
            "Kelvin relation violated: Π_s = {}, T·S_s = {}",
            mat.heat_calc.spin_peltier_coeff,
            expected_pi_s
        );
    }

    /// 12. compute_all with nonzero ∇T gives nonzero spin_seebeck_current.
    #[test]
    fn test_compute_all_nonzero_for_grad_t() {
        let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
        let grad_t = grad_t_x(500.0);
        let j_spin = zero_vec();
        let result = mat
            .compute_all(&grad_t, &j_spin)
            .expect("compute_all failed");
        assert!(
            result.spin_seebeck_current.magnitude() > 0.0,
            "spin Seebeck current should be nonzero"
        );
    }

    /// 13. compute_all reciprocity_satisfied flag is true for preset materials.
    #[test]
    fn test_compute_all_reciprocity_flag() {
        for (name, mat) in [
            ("YIG/Pt", SpinCaloritronicsMaterial::yig_pt(300.0)),
            ("Fe/Pt", SpinCaloritronicsMaterial::fe_pt(300.0)),
            ("CoFeB/Pt", SpinCaloritronicsMaterial::cofeb_pt(300.0)),
        ] {
            let grad_t = grad_t_x(100.0);
            let j_spin = small_spin();
            let result = mat
                .compute_all(&grad_t, &j_spin)
                .expect("compute_all failed");
            assert!(
                result.reciprocity_satisfied,
                "{name}: reciprocity_satisfied should be true"
            );
        }
    }

    /// 14. HeatCurrentCalculator fourier_current = −κ · ∇T.
    #[test]
    fn test_heat_calculator_fourier() {
        let kappa = 46.0_f64;
        let calc = HeatCurrentCalculator::new(kappa, 0.0, 0.0);
        let grad_t = Vector3::new(1000.0, 0.0, 0.0);
        let j_q = calc.fourier_current(&grad_t);
        // Expected: x = -46.0 * 1000.0 = -46_000
        assert!(
            (j_q.x - (-kappa * 1000.0)).abs() < 1.0e-8,
            "Fourier current x mismatch: got {}, expected {}",
            j_q.x,
            -kappa * 1000.0
        );
        assert!(j_q.y.abs() < f64::EPSILON);
        assert!(j_q.z.abs() < f64::EPSILON);
    }

    /// 15. total_heat_current combines all three mechanisms correctly.
    #[test]
    fn test_total_heat_current_combines_all_mechanisms() {
        let kappa = 46.0_f64;
        let pi_charge = -1.5e-3_f64;
        let pi_spin = 0.3_f64;
        let calc = HeatCurrentCalculator::new(kappa, pi_charge, pi_spin);

        let grad_t = Vector3::new(1000.0, 0.0, 0.0);
        let j_charge = Vector3::new(1.0e6, 0.0, 0.0);
        let j_spin = Vector3::new(1.0e3, 0.0, 0.0);

        let total = calc.total_heat_current(&grad_t, &j_charge, &j_spin);

        // Expected: (-κ·∇T_x) + (Π·j_c_x) + (Π_s·j_s_x)
        let expected_x = (-kappa * 1000.0) + (pi_charge * 1.0e6) + (pi_spin * 1.0e3);

        assert!(
            (total.x - expected_x).abs() < 1.0e-6,
            "total heat current x mismatch: got {}, expected {}",
            total.x,
            expected_x
        );

        // Individual contributions
        let fourier = calc.fourier_current(&grad_t);
        let peltier = calc.peltier_current(&j_charge);
        let spin_peltier = calc.spin_peltier_current(&j_spin);

        // Sum should match total
        let sum_x = fourier.x + peltier.x + spin_peltier.x;
        assert!(
            (total.x - sum_x).abs() < 1.0e-8,
            "total should equal sum of components"
        );
    }
}
