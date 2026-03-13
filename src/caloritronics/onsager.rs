//! Onsager Transport Matrix for Spin Caloritronics
//!
//! This module implements the Onsager reciprocal relations coupling charge,
//! spin, and heat currents in magnetic heterostructures. The Onsager formalism
//! provides a unified thermodynamic framework for all caloritronic cross-effects.
//!
//! ## Physical Background
//!
//! In linear response theory, generalized fluxes (currents) J_i are related
//! to thermodynamic forces (gradients) X_j by the Onsager matrix L:
//!
//!   J_i = Σ_j L_ij * X_j
//!
//! The Onsager reciprocity theorem states: L_ij = T * L_ji (for time-reversal
//! symmetric systems), ensuring thermodynamic consistency.
//!
//! For spin caloritronics, the relevant currents are:
//! - Charge current density j_c [A/m²]
//! - Spin current density j_s [A/m²] (angular momentum current)
//! - Heat current density j_Q [W/m²]
//!
//! And the corresponding forces:
//! - Electric field E [V/m]
//! - Spin chemical potential gradient ∇μ_s [V/m]
//! - Temperature gradient ∇T/T [1/m]
//!
//! ## Key References
//!
//! - L. Onsager, "Reciprocal Relations in Irreversible Processes", Phys. Rev. 37, 405 (1931)
//! - G. E. W. Bauer et al., "Spin caloritronics", Nat. Mater. 11, 391 (2012)
//! - K. Uchida et al., "Spin Seebeck insulator", Nat. Mater. 9, 894 (2010)

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::vector3::Vector3;

/// Onsager transport matrix for spin caloritronics.
///
/// Encodes the linear-response coupling between charge, spin, and heat currents
/// via Onsager reciprocity. This is the central object of spin caloritronics,
/// unifying the Seebeck effect, spin Seebeck effect, Peltier effect, spin Peltier
/// effect, and Nernst/Ettingshausen effects into one framework.
///
/// The 2×2 Onsager matrix for charge-heat coupling (ignoring spin for a moment):
///
///   \[j_c\]   \[σ      σ·S \]   \[E     \]
///   \[j_Q\] = \[Π·σ    κ+Π·S·σ\] \[∇T    \]
///
/// where Π = T·S is the Peltier coefficient (Kelvin relation).
///
/// The spin Seebeck coefficient S_s couples spin current to temperature gradient:
///   j_s = S_s · σ · ∇T
///
/// # Material Presets
///
/// ```rust
/// use spintronics::caloritronics::OnsagerMatrix;
///
/// let yig_pt = OnsagerMatrix::yig_pt(300.0);
/// let fe_pt  = OnsagerMatrix::fe_pt(300.0);
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OnsagerMatrix {
    /// Temperature \[K\]
    pub temperature: f64,
    /// Electrical conductivity σ [S/m]
    pub conductivity: f64,
    /// Seebeck coefficient S_e [V/K]
    ///
    /// Negative for electron-dominated transport (most metals)
    pub seebeck: f64,
    /// Spin Seebeck coefficient S_s [A/(m·K)]
    ///
    /// Governs thermal generation of spin current: j_s = S_s · σ · ∇T
    pub spin_seebeck: f64,
    /// Anomalous Hall angle θ_H \[dimensionless\]
    ///
    /// Relates transverse (Hall) conductivity to longitudinal: σ_xy = θ_H · σ_xx
    pub hall_angle: f64,
    /// Thermal conductivity κ [W/(m·K)]
    pub thermal_conductivity: f64,
}

impl OnsagerMatrix {
    /// Create a new Onsager matrix with explicit parameters.
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]
    /// * `conductivity` - Electrical conductivity σ [S/m]
    /// * `seebeck` - Seebeck coefficient S_e [V/K]
    /// * `spin_seebeck` - Spin Seebeck coefficient S_s [A/(m·K)]
    /// * `hall_angle` - Anomalous Hall/Spin Hall angle θ_H
    /// * `thermal_conductivity` - Thermal conductivity κ [W/(m·K)]
    pub fn new(
        temperature: f64,
        conductivity: f64,
        seebeck: f64,
        spin_seebeck: f64,
        hall_angle: f64,
        thermal_conductivity: f64,
    ) -> Self {
        Self {
            temperature,
            conductivity,
            seebeck,
            spin_seebeck,
            hall_angle,
            thermal_conductivity,
        }
    }

    /// YIG/Pt bilayer at given temperature.
    ///
    /// YIG (Y₃Fe₅O₁₂) is the canonical spin Seebeck insulator: negligible
    /// charge current but large spin Seebeck coefficient due to long magnon
    /// mean free path. Pt has strong spin-orbit coupling (large Hall angle).
    ///
    /// Values are approximate room-temperature parameters for the composite system.
    pub fn yig_pt(temperature: f64) -> Self {
        Self {
            temperature,
            conductivity: 5.0e6,
            seebeck: -5.0e-6,
            spin_seebeck: 1.0e-3,
            hall_angle: 0.1,
            thermal_conductivity: 46.0,
        }
    }

    /// Fe/Pt bilayer at given temperature.
    ///
    /// Iron has large spontaneous magnetization and moderate Seebeck coefficient.
    /// Combined with the Pt overlayer, this gives a strong spin Peltier response.
    pub fn fe_pt(temperature: f64) -> Self {
        Self {
            temperature,
            conductivity: 1.0e7,
            seebeck: -15.0e-6,
            spin_seebeck: 5.0e-3,
            hall_angle: 0.08,
            thermal_conductivity: 80.0,
        }
    }

    /// CoFeB/Pt bilayer at given temperature.
    ///
    /// CoFeB is the preferred magnetic electrode in magnetic tunnel junctions
    /// and spin-orbit torque devices. It has a large anomalous Nernst coefficient.
    pub fn cofeb_pt(temperature: f64) -> Self {
        Self {
            temperature,
            conductivity: 8.0e6,
            seebeck: -10.0e-6,
            spin_seebeck: 3.0e-3,
            hall_angle: 0.06,
            thermal_conductivity: 60.0,
        }
    }

    /// Check Onsager reciprocity: L_ij = T · L_ji.
    ///
    /// For the charge–heat subblock, the off-diagonal elements must satisfy:
    ///
    ///   L_cQ = σ · S_e   (Seebeck)
    ///   L_Qc = Π · σ = T · S_e · σ   (Peltier, Kelvin relation Π = T·S)
    ///
    /// Reciprocity requires L_Qc = T · L_cQ, i.e. T·σ·S_e = T·σ·S_e ✓
    ///
    /// For the spin–heat channel:
    ///   L_sQ = S_s · σ
    ///   L_Qs = Π_s · σ = T · S_s · σ
    ///   Reciprocity: L_Qs = T · L_sQ ✓
    ///
    /// This method computes the relative deviation from exact reciprocity.
    /// For ideal (analytically constructed) matrices this is zero by construction,
    /// but numerical evaluation lets tests confirm the property.
    ///
    /// # Returns
    /// Relative reciprocity error (dimensionless, should be ≈ 0)
    pub fn reciprocity_error(&self) -> f64 {
        // Charge-heat coupling: L_cQ = σ·S_e, L_Qc = T·σ·S_e
        // Reciprocity: L_Qc == T · L_cQ  →  T·σ·S = T·σ·S  (exact by construction)
        //
        // Spin-heat coupling: L_sQ = S_s·σ, L_Qs = T·S_s·σ
        // Reciprocity: L_Qs == T · L_sQ  (exact by construction)
        //
        // We verify this numerically to expose any floating-point discrepancy.
        let l_cq = self.conductivity * self.seebeck;
        let l_qc = self.temperature * self.conductivity * self.seebeck;
        let l_sq = self.spin_seebeck * self.conductivity;
        let l_qs = self.temperature * self.spin_seebeck * self.conductivity;

        // Relative error = |L_Qc - T·L_cQ| / |L_Qc| + |L_Qs - T·L_sQ| / |L_Qs|
        let charge_err = if l_qc.abs() > f64::EPSILON {
            (l_qc - self.temperature * l_cq).abs() / l_qc.abs()
        } else {
            0.0
        };
        let spin_err = if l_qs.abs() > f64::EPSILON {
            (l_qs - self.temperature * l_sq).abs() / l_qs.abs()
        } else {
            0.0
        };
        charge_err + spin_err
    }

    /// Compute spin current density from a temperature gradient (spin Seebeck effect).
    ///
    /// j_s = S_s · σ · ∇T   [A/m²]
    ///
    /// The spin current flows parallel to the temperature gradient, with magnitude
    /// determined by the spin Seebeck coefficient and conductivity.
    pub fn spin_current_from_grad_t(&self, grad_t: &Vector3<f64>) -> Vector3<f64> {
        let prefactor = self.spin_seebeck * self.conductivity;
        Vector3::new(
            prefactor * grad_t.x,
            prefactor * grad_t.y,
            prefactor * grad_t.z,
        )
    }

    /// Compute heat current from a spin current (spin Peltier effect).
    ///
    /// The spin Peltier coefficient is given by the Kelvin relation:
    ///   Π_s = T · S_s
    ///
    /// Heat current: j_Q = Π_s · j_s = T · S_s · j_s   [W/m²]
    pub fn heat_current_from_spin_current(&self, j_spin: &Vector3<f64>) -> Vector3<f64> {
        let pi_s = self.temperature * self.spin_seebeck;
        Vector3::new(pi_s * j_spin.x, pi_s * j_spin.y, pi_s * j_spin.z)
    }

    /// Compute the anomalous (or ordinary) Nernst voltage.
    ///
    /// The Nernst coefficient relates the transverse electric field to the
    /// longitudinal temperature gradient via the Hall conductivity:
    ///
    ///   ν = -S_e · (σ_xy / σ_xx) ≈ -S_e · tan(θ_H) ≈ -S_e · θ_H
    ///
    /// For small Hall angles tan(θ_H) ≈ θ_H.
    ///
    /// # Arguments
    /// * `grad_t` - Longitudinal temperature gradient [K/m]
    ///
    /// # Returns
    /// Transverse Nernst voltage [V/m] (voltage per unit length)
    pub fn nernst_voltage(&self, grad_t: f64) -> f64 {
        // ν = -S_e · θ_H · grad_T
        -self.seebeck * self.hall_angle * grad_t
    }

    /// Compute all three currents simultaneously for given fields.
    ///
    /// This is the full Onsager response calculation including:
    /// - Charge current from electric field (Ohm) and temperature gradient (Seebeck)
    /// - Spin current from temperature gradient (spin Seebeck)
    /// - Heat current from electric field (Peltier) and temperature gradient (Fourier)
    ///
    /// Equations:
    ///   j_c = σ·E + σ·S_e·∇T         [A/m²]
    ///   j_s = S_s·σ·∇T                [A/m²]
    ///   j_Q = Π·σ·E - κ·∇T           [W/m²]
    ///       = T·S_e·σ·E - κ·∇T
    ///
    /// # Arguments
    /// * `grad_t` - Temperature gradient vector [K/m]
    /// * `e_field` - Electric field vector [V/m]
    pub fn all_currents(&self, grad_t: &Vector3<f64>, e_field: &Vector3<f64>) -> AllCurrents {
        // Charge current: j_c = σ·(E + S_e·∇T)
        let seebeck_force_x = self.seebeck * grad_t.x;
        let seebeck_force_y = self.seebeck * grad_t.y;
        let seebeck_force_z = self.seebeck * grad_t.z;
        let charge_current = Vector3::new(
            self.conductivity * (e_field.x + seebeck_force_x),
            self.conductivity * (e_field.y + seebeck_force_y),
            self.conductivity * (e_field.z + seebeck_force_z),
        );

        // Spin current: j_s = S_s·σ·∇T
        let spin_current = self.spin_current_from_grad_t(grad_t);

        // Heat current: j_Q = T·S_e·σ·E - κ·∇T  (Peltier + Fourier)
        let peltier_coeff = self.temperature * self.seebeck;
        let heat_current = Vector3::new(
            peltier_coeff * self.conductivity * e_field.x - self.thermal_conductivity * grad_t.x,
            peltier_coeff * self.conductivity * e_field.y - self.thermal_conductivity * grad_t.y,
            peltier_coeff * self.conductivity * e_field.z - self.thermal_conductivity * grad_t.z,
        );

        AllCurrents {
            charge_current,
            spin_current,
            heat_current,
        }
    }
}

/// All three physical currents computed by the Onsager matrix.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AllCurrents {
    /// Charge current density j_c [A/m²]
    pub charge_current: Vector3<f64>,
    /// Spin current density j_s [A/m²]
    pub spin_current: Vector3<f64>,
    /// Heat current density j_Q [W/m²]
    pub heat_current: Vector3<f64>,
}
