//! Coupled spin-mechanical dynamics
//!
//! Simultaneous solution of LLG + Newton's equations with angular momentum exchange

use crate::constants::GAMMA;
use crate::vector3::Vector3;

/// Spin-mechanical coupling system
#[derive(Debug, Clone)]
pub struct SpinMechanicalCoupling {
    /// Magnetization (normalized)
    pub magnetization: Vector3<f64>,

    /// Angular velocity \[rad/s\]
    pub omega: Vector3<f64>,

    /// Moment of inertia [kg·m²]
    pub moment_of_inertia: f64,

    /// Gilbert damping
    pub alpha: f64,

    /// Saturation magnetization \[A/m\]
    pub ms: f64,

    /// Sample volume \[m³\]
    pub volume: f64,
}

impl SpinMechanicalCoupling {
    /// Create new coupled system
    pub fn new(
        magnetization: Vector3<f64>,
        moment_of_inertia: f64,
        alpha: f64,
        ms: f64,
        volume: f64,
    ) -> Self {
        Self {
            magnetization: magnetization.normalize(),
            omega: Vector3::new(0.0, 0.0, 0.0),
            moment_of_inertia,
            alpha,
            ms,
            volume,
        }
    }

    /// Evolve coupled system by one time step
    ///
    /// Couples:
    /// 1. LLG with Barnett field: H_eff includes -Ω/γ term
    /// 2. Newton's equation with spin torque: τ = -dL_spin/dt
    pub fn evolve(&mut self, h_ext: Vector3<f64>, dt: f64) {
        // Total field = external + Barnett
        let h_barnett = self.omega * (-1.0 / GAMMA);
        let h_total = h_ext + h_barnett;

        // LLG dynamics
        let m_cross_h = self.magnetization.cross(&h_total);
        let damping = self.magnetization.cross(&m_cross_h) * self.alpha;

        let dm_dt = (m_cross_h + damping) * (-GAMMA / (1.0 + self.alpha * self.alpha));

        // Update magnetization
        self.magnetization = (self.magnetization + dm_dt * dt).normalize();

        // Mechanical torque from spin dynamics
        // τ = -V M_s (dm/dt)
        let mechanical_torque = dm_dt * (-self.volume * self.ms);

        // Angular acceleration
        let alpha_mech = mechanical_torque * (1.0 / self.moment_of_inertia);

        // Update angular velocity
        self.omega = self.omega + alpha_mech * dt;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coupled_system_creation() {
        let system =
            SpinMechanicalCoupling::new(Vector3::new(1.0, 0.0, 0.0), 1.0e-15, 0.01, 1.0e6, 1.0e-18);

        assert!((system.magnetization.magnitude() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_angular_momentum_exchange() {
        let mut system =
            SpinMechanicalCoupling::new(Vector3::new(1.0, 0.0, 0.0), 1.0e-15, 0.01, 1.0e6, 1.0e-18);

        let h_ext = Vector3::new(0.0, 0.0, 1.0e5);

        // Initial rotation should be zero
        assert!(system.omega.magnitude() < 1e-20);

        // Evolve - magnetization will precess, inducing mechanical rotation
        for _ in 0..100 {
            system.evolve(h_ext, 1.0e-13);
        }

        // Some rotation should develop
        assert!(system.omega.magnitude() > 0.0);
    }
}
