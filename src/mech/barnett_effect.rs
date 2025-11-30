//! Barnett Effect: Rotation induces magnetization
//!
//! When a ferromagnetic sample rotates, it becomes magnetized even without
//! an external magnetic field. This is because rotation in the laboratory
//! frame appears as an effective magnetic field in the rotating frame.

use crate::constants::{GAMMA, HBAR};
use crate::vector3::Vector3;

/// Barnett magnetization calculator
#[derive(Debug, Clone)]
pub struct BarnettMagnetization {
    /// Gyromagnetic ratio \[rad/(s·T)\]
    pub gamma: f64,

    /// Number density of magnetic moments [1/m³]
    pub moment_density: f64,

    /// Saturation magnetization \[A/m\]
    pub ms: f64,
}

impl BarnettMagnetization {
    /// Create for a ferromagnetic material
    ///
    /// # Arguments
    /// * `ms` - Saturation magnetization \[A/m\]
    /// * `moment_density` - Magnetic moment density [1/m³]
    pub fn new(ms: f64, moment_density: f64) -> Self {
        Self {
            gamma: GAMMA,
            moment_density,
            ms,
        }
    }

    /// Create for iron
    pub fn iron() -> Self {
        Self::new(
            1.7e6,  // Ms for Fe
            8.5e28, // Moment density
        )
    }

    /// Create for permalloy
    pub fn permalloy() -> Self {
        Self::new(8.0e5, 5.8e28)
    }

    /// Calculate induced magnetization from rotation
    ///
    /// M_B = (ℏ n / 2) × Ω
    ///
    /// where n is the density of magnetic moments
    ///
    /// # Arguments
    /// * `omega` - Angular velocity \[rad/s\]
    ///
    /// # Returns
    /// Magnetization vector \[A/m\]
    pub fn magnetization_from_rotation(&self, omega: Vector3<f64>) -> Vector3<f64> {
        // M = (ℏ n / 2) × Ω
        omega * (HBAR * self.moment_density / 2.0)
    }

    /// Calculate fractional magnetization
    ///
    /// m = M_B / M_s
    pub fn fractional_magnetization(&self, omega: Vector3<f64>) -> f64 {
        let m_b = self.magnetization_from_rotation(omega);
        m_b.magnitude() / self.ms
    }

    /// Calculate rotation rate needed for saturation
    ///
    /// Ω_sat = 2 M_s / (ℏ n)
    ///
    /// This is typically extremely high (unrealistic for macroscopic samples)
    pub fn saturation_rotation_rate(&self) -> f64 {
        (2.0 * self.ms) / (HBAR * self.moment_density)
    }

    /// Calculate Barnett field strength
    ///
    /// H_B = -Ω / γ
    ///
    /// This is the effective field experienced by spins in the rotating frame
    pub fn barnett_field(&self, omega: Vector3<f64>) -> Vector3<f64> {
        omega * (-1.0 / self.gamma)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_barnett_creation() {
        let barnett = BarnettMagnetization::iron();
        assert!(barnett.ms > 0.0);
    }

    #[test]
    fn test_magnetization_from_rotation() {
        let barnett = BarnettMagnetization::permalloy();
        let omega = Vector3::new(0.0, 0.0, 1000.0); // 1000 rad/s

        let m = barnett.magnetization_from_rotation(omega);
        assert!(m.magnitude() > 0.0);
        assert!(m.z > 0.0); // Should be parallel to rotation
    }

    #[test]
    fn test_zero_rotation() {
        let barnett = BarnettMagnetization::iron();
        let omega = Vector3::new(0.0, 0.0, 0.0);

        let m = barnett.magnetization_from_rotation(omega);
        assert!(m.magnitude() < 1e-50);
    }

    #[test]
    fn test_fractional_magnetization_small() {
        let barnett = BarnettMagnetization::iron();
        let omega = Vector3::new(0.0, 0.0, 100.0); // Moderate rotation

        let frac = barnett.fractional_magnetization(omega);
        // Should be very small for realistic rotation rates
        assert!(frac < 1.0);
        assert!(frac > 0.0);
    }
}
