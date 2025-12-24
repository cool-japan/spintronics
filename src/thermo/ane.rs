//! Anomalous Nernst Effect (ANE)
//!
//! The ANE generates a transverse electric field in response to a temperature
//! gradient in magnetic materials. Unlike the normal Nernst effect, ANE does
//! not require an external magnetic field.
//!
//! E_ANE = α_ANE × (∇T × m)
//!
//! where α_ANE is the anomalous Nernst coefficient and m is magnetization direction.

use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::vector3::Vector3;

/// Anomalous Nernst Effect calculator
///
/// The ANE is particularly strong in materials with large spin-orbit coupling
/// and is being explored for thermoelectric energy harvesting.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AnomalousNernst {
    /// Anomalous Nernst coefficient [V/(K·m)]
    ///
    /// Typical values:
    /// - Fe: ~1-3 μV/K
    /// - CoFeB: ~3-5 μV/K
    /// - Heusler alloys: up to 6 μV/K
    pub alpha_ane: f64,

    /// Magnetization direction (normalized)
    pub magnetization: Vector3<f64>,
}

impl Default for AnomalousNernst {
    fn default() -> Self {
        Self {
            alpha_ane: 3.0e-6, // 3 μV/K (typical for CoFeB)
            magnetization: Vector3::new(0.0, 0.0, 1.0),
        }
    }
}

impl AnomalousNernst {
    /// Create ANE properties for iron
    pub fn iron() -> Self {
        Self {
            alpha_ane: 2.0e-6,
            magnetization: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Create ANE properties for CoFeB
    pub fn cofeb() -> Self {
        Self {
            alpha_ane: 4.0e-6,
            magnetization: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Create ANE properties for Heusler alloy
    pub fn heusler() -> Self {
        Self {
            alpha_ane: 6.0e-6,
            magnetization: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Calculate ANE electric field
    ///
    /// # Arguments
    /// * `grad_t` - Temperature gradient vector \[K/m\]
    ///
    /// # Returns
    /// Electric field vector [V/m]
    ///
    /// # Physical Interpretation
    /// The ANE creates a transverse voltage perpendicular to both the
    /// temperature gradient and magnetization, enabling thermal energy
    /// harvesting without external magnetic fields.
    #[inline]
    pub fn electric_field(&self, grad_t: Vector3<f64>) -> Vector3<f64> {
        // E_ANE = α_ANE × (∇T × m)
        let cross = grad_t.cross(&self.magnetization);
        cross * self.alpha_ane
    }

    /// Calculate voltage across a sample
    ///
    /// # Arguments
    /// * `grad_t` - Temperature gradient \[K/m\]
    /// * `width` - Sample width perpendicular to E field \[m\]
    ///
    /// # Returns
    /// Voltage \[V\]
    #[inline]
    pub fn voltage(&self, grad_t: Vector3<f64>, width: f64) -> f64 {
        let e_field = self.electric_field(grad_t);
        e_field.magnitude() * width
    }

    /// Calculate power factor for thermoelectric figure of merit
    ///
    /// Power factor = α²/ρ, where α is Seebeck coefficient
    /// For ANE: PF = α_ANE² / ρ
    pub fn power_factor(&self, resistivity: f64) -> f64 {
        self.alpha_ane.powi(2) / resistivity
    }

    /// Builder method to set ANE coefficient
    pub fn with_alpha_ane(mut self, alpha_ane: f64) -> Self {
        self.alpha_ane = alpha_ane;
        self
    }

    /// Builder method to set magnetization direction
    pub fn with_magnetization(mut self, magnetization: Vector3<f64>) -> Self {
        self.magnetization = magnetization.normalize();
        self
    }
}

impl fmt::Display for AnomalousNernst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "AnomalousNernst: α_ANE={:.1} μV/K, m=[{:.2},{:.2},{:.2}]",
            self.alpha_ane * 1e6,
            self.magnetization.x,
            self.magnetization.y,
            self.magnetization.z
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ane_perpendicularity() {
        let ane = AnomalousNernst::default();
        let grad_t = Vector3::new(1000.0, 0.0, 0.0); // 1000 K/m in x

        let e_field = ane.electric_field(grad_t);

        // E should be perpendicular to both ∇T and m
        assert!(e_field.dot(&grad_t).abs() < 1e-20);
        assert!(e_field.dot(&ane.magnetization).abs() < 1e-20);
    }

    #[test]
    fn test_ane_zero_gradient() {
        let ane = AnomalousNernst::default();
        let grad_t = Vector3::new(0.0, 0.0, 0.0);

        let e_field = ane.electric_field(grad_t);
        assert!(e_field.magnitude() < 1e-30);
    }

    #[test]
    fn test_ane_parallel_gradient_magnetization() {
        let ane = AnomalousNernst::cofeb();
        let grad_t = Vector3::new(0.0, 0.0, 1000.0); // Parallel to m

        let e_field = ane.electric_field(grad_t);
        // Parallel vectors give zero cross product
        assert!(e_field.magnitude() < 1e-20);
    }

    #[test]
    fn test_material_coefficients() {
        let fe = AnomalousNernst::iron();
        let cofeb = AnomalousNernst::cofeb();
        let heusler = AnomalousNernst::heusler();

        assert!(fe.alpha_ane < cofeb.alpha_ane);
        assert!(cofeb.alpha_ane < heusler.alpha_ane);
    }
}
