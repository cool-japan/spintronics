//! Barnett Effect
//!
//! The Barnett effect: mechanical rotation induces magnetization.
//!
//! When a material rotates with angular velocity Ω, it experiences an
//! effective magnetic field (Barnett field):
//!
//! H_B = -Ω / γ
//!
//! where γ is the gyromagnetic ratio.

use crate::constants::GAMMA;
use crate::vector3::Vector3;

/// Barnett field calculator
///
/// Converts mechanical rotation (angular velocity) into an effective
/// magnetic field that acts on electron spins.
#[derive(Debug, Clone)]
pub struct BarnettField {
    /// Gyromagnetic ratio \[rad/(s·T)\]
    pub gamma: f64,
}

impl Default for BarnettField {
    fn default() -> Self {
        Self { gamma: GAMMA }
    }
}

impl BarnettField {
    /// Create a new Barnett field calculator
    pub fn new(gamma: f64) -> Self {
        Self { gamma }
    }

    /// Calculate Barnett field from angular velocity
    ///
    /// # Arguments
    /// * `omega` - Angular velocity vector \[rad/s\]
    ///
    /// # Returns
    /// Effective magnetic field \[A/m\]
    ///
    /// # Physical Interpretation
    /// A rotating frame creates a fictitious magnetic field that
    /// couples to the electron's magnetic moment, causing spin polarization
    /// without any actual magnetic field.
    ///
    /// For liquid metals rotating at ~1000 rpm, this can create
    /// fields equivalent to ~1 mT.
    pub fn field_from_rotation(&self, omega: Vector3<f64>) -> Vector3<f64> {
        // H_B = -ω / γ
        // Convert to A/m: H \[A/m\] = (ω \[rad/s\] / γ \[rad/(s·T)\]) / μ₀
        let mu0 = 1.256637e-6;
        omega * (-1.0 / (self.gamma * mu0))
    }

    /// Calculate Barnett field from vorticity (fluid mechanics)
    ///
    /// In a fluid, local vorticity ω = ∇ × v acts as local rotation
    ///
    /// # Arguments
    /// * `vorticity` - Fluid vorticity [1/s]
    ///
    /// # Returns
    /// Local Barnett field \[A/m\]
    pub fn field_from_vorticity(&self, vorticity: Vector3<f64>) -> Vector3<f64> {
        // Vorticity is 2 × angular velocity in fluid mechanics
        let omega = vorticity * 0.5;
        self.field_from_rotation(omega)
    }

    /// Calculate spin polarization from rotation
    ///
    /// P = -ω / (γ k_B T) × μ_B
    ///
    /// For small rotations, this gives the fractional spin polarization
    pub fn spin_polarization(&self, omega: Vector3<f64>, temperature: f64) -> f64 {
        use crate::constants::{KB, MU_B};

        if temperature < 1.0 {
            return 0.0;
        }

        let omega_mag = omega.magnitude();
        let thermal_energy = KB * temperature;

        // Polarization ∝ magnetic energy / thermal energy
        (MU_B * omega_mag) / (self.gamma * thermal_energy)
    }

    /// Calculate spin current density from vorticity gradient
    ///
    /// J_s ∝ ∇ω (spin current flows from high to low vorticity)
    ///
    /// This is the key to Saitoh's liquid metal spin current generation
    pub fn spin_current_from_vorticity_gradient(
        &self,
        vorticity_gradient: f64,
        density: f64,
    ) -> f64 {
        use crate::constants::HBAR;

        // Simplified model: J_s ~ (ℏ n / 2) × ∇ω
        // where n is carrier density
        (HBAR * density / 2.0) * vorticity_gradient
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_barnett_field_creation() {
        let barnett = BarnettField::default();
        assert_eq!(barnett.gamma, GAMMA);
    }

    #[test]
    fn test_field_from_rotation() {
        let barnett = BarnettField::default();
        let omega = Vector3::new(0.0, 0.0, 100.0); // 100 rad/s ~ 955 rpm

        let h_field = barnett.field_from_rotation(omega);

        // Field should be anti-parallel to rotation
        assert!(h_field.z < 0.0);
        assert!(h_field.magnitude() > 0.0);
    }

    #[test]
    fn test_zero_rotation() {
        let barnett = BarnettField::default();
        let omega = Vector3::new(0.0, 0.0, 0.0);

        let h_field = barnett.field_from_rotation(omega);
        assert!(h_field.magnitude() < 1e-20);
    }

    #[test]
    fn test_vorticity_to_field() {
        let barnett = BarnettField::default();
        let vorticity = Vector3::new(0.0, 0.0, 1000.0); // 1000 1/s

        let h_field = barnett.field_from_vorticity(vorticity);
        assert!(h_field.magnitude() > 0.0);
    }

    #[test]
    fn test_spin_polarization() {
        let barnett = BarnettField::default();
        let omega = Vector3::new(0.0, 0.0, 1000.0);
        let temperature = 300.0; // K

        let polarization = barnett.spin_polarization(omega, temperature);
        assert!(polarization > 0.0);
        assert!(polarization < 1.0); // Should be small at room temperature
    }
}
