//! Einstein-de Haas Effect
//!
//! Change in magnetization causes mechanical rotation due to conservation
//! of angular momentum. The "twin" of the Barnett effect.

use crate::constants::HBAR;
use crate::vector3::Vector3;

/// Einstein-de Haas calculator
#[derive(Debug, Clone)]
pub struct EinsteinDeHaas {
    /// Moment of inertia [kg·m²]
    pub moment_of_inertia: f64,

    /// Sample volume \[m³\]
    pub volume: f64,

    /// Density [kg/m³]
    pub density: f64,
}

impl EinsteinDeHaas {
    /// Create new calculator
    pub fn new(moment_of_inertia: f64, volume: f64, density: f64) -> Self {
        Self {
            moment_of_inertia,
            volume,
            density,
        }
    }

    /// Calculate angular velocity from magnetization change
    ///
    /// Conservation of angular momentum:
    /// L_spin + L_orbital = const
    /// ΔL_spin = -ΔL_mech
    ///
    /// Ω = -(ℏ / I) × Δ(M × V)
    ///
    /// # Arguments
    /// * `delta_m` - Change in magnetization \[A/m\]
    ///
    /// # Returns
    /// Angular velocity \[rad/s\]
    pub fn angular_velocity_from_magnetization_change(
        &self,
        delta_m: Vector3<f64>,
    ) -> Vector3<f64> {
        // Total spin angular momentum change
        let delta_l_spin = delta_m * (HBAR * self.volume);

        // Mechanical angular momentum change (opposite)
        let delta_l_mech = delta_l_spin * (-1.0);

        // Ω = L / I
        delta_l_mech * (1.0 / self.moment_of_inertia)
    }

    /// Calculate rotation angle from magnetization reversal
    ///
    /// For a complete magnetization flip (M → -M)
    pub fn rotation_angle_magnetization_flip(&self, ms: f64) -> f64 {
        let delta_m_magnitude = 2.0 * ms;
        let delta_l = HBAR * delta_m_magnitude * self.volume;

        // Assuming the rotation decays exponentially, the total angle is:
        // θ = ∫ ω dt = L / (I × damping_rate)
        // For estimation, use: θ ~ L / I
        delta_l / self.moment_of_inertia
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edh_creation() {
        let edh = EinsteinDeHaas::new(1.0e-15, 1.0e-18, 8000.0);
        assert!(edh.moment_of_inertia > 0.0);
    }

    #[test]
    fn test_angular_velocity_calculation() {
        let edh = EinsteinDeHaas::new(1.0e-15, 1.0e-18, 8000.0);
        let delta_m = Vector3::new(0.0, 0.0, 1.0e5);

        let omega = edh.angular_velocity_from_magnetization_change(delta_m);
        assert!(omega.magnitude() > 0.0);
    }
}
