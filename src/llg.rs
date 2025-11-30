//! Landau-Lifshitz-Gilbert (LLG) equation solver
//!
//! This module provides a solver for the LLG equation which governs
//! the time evolution of magnetization in magnetic materials.
//!
//! # LLG Equation
//!
//! The LLG equation is:
//! ```text
//! dm/dt = -γ/(1+α²) [m × H_eff + α m × (m × H_eff)]
//! ```
//!
//! where:
//! - m is the normalized magnetization direction
//! - γ is the gyromagnetic ratio
//! - α is the Gilbert damping constant
//! - H_eff is the effective magnetic field

use crate::constants::GAMMA;
use crate::material::Ferromagnet;
use crate::vector3::Vector3;

/// LLG equation solver
pub struct LLGSolver {
    /// Material properties
    pub material: Ferromagnet,
    /// External magnetic field (A/m)
    pub h_ext: Vector3<f64>,
}

impl LLGSolver {
    /// Create a new LLG solver
    pub fn new(material: Ferromagnet) -> Self {
        LLGSolver {
            material,
            h_ext: Vector3::new(0.0, 0.0, 0.0),
        }
    }

    /// Calculate effective field for a given magnetization
    ///
    /// Currently includes:
    /// - External field
    /// - Uniaxial anisotropy field
    ///
    /// # Arguments
    /// * `m` - Magnetization vector (A/m)
    ///
    /// # Returns
    /// Effective field H_eff (A/m)
    pub fn effective_field(&self, m: Vector3<f64>) -> Vector3<f64> {
        let mut h_eff = self.h_ext;

        // Anisotropy field: H_anis = (2K_u/μ₀M_s²)(m·e)e
        if self.material.anisotropy_k.abs() > 1e-10 {
            let mu0 = 4.0 * std::f64::consts::PI * 1e-7;
            let m_norm = m.normalize();
            let m_dot_e = m_norm.dot(&self.material.easy_axis);
            let h_anis_coeff =
                2.0 * self.material.anisotropy_k / (mu0 * self.material.ms * self.material.ms);
            h_eff = h_eff + self.material.easy_axis * (h_anis_coeff * m_dot_e);
        }

        h_eff
    }

    /// Calculate dm/dt from the LLG equation
    ///
    /// # Arguments
    /// * `m` - Magnetization vector (A/m)
    ///
    /// # Returns
    /// Time derivative dm/dt (A/m/s)
    pub fn dmdt(&self, m: Vector3<f64>) -> Vector3<f64> {
        let h_eff = self.effective_field(m);
        let m_norm = m.normalize();

        let gamma = GAMMA; // Use global constant
        let alpha = self.material.alpha;

        // LLG equation:
        // dm/dt = -γ/(1+α²) [m × H_eff + α m × (m × H_eff)]
        let m_cross_h = m_norm.cross(&h_eff);
        let m_cross_m_cross_h = m_norm.cross(&m_cross_h);

        let prefactor = -gamma / (1.0 + alpha * alpha);
        (m_cross_h + m_cross_m_cross_h * alpha) * (prefactor * self.material.ms)
    }

    /// Advance magnetization by one time step using Euler method
    ///
    /// # Arguments
    /// * `m` - Current magnetization (A/m)
    /// * `dt` - Time step (s)
    ///
    /// # Returns
    /// Updated magnetization (A/m)
    pub fn step_euler(&self, m: Vector3<f64>, dt: f64) -> Vector3<f64> {
        let dmdt = self.dmdt(m);
        let m_new = m + dmdt * dt;

        // Renormalize to maintain |m| = Ms
        let mag = m_new.magnitude();
        if mag > 1e-10 {
            m_new * (self.material.ms / mag)
        } else {
            m
        }
    }

    /// Advance magnetization by one time step using RK4 method
    ///
    /// 4th-order Runge-Kutta provides better accuracy than Euler
    ///
    /// # Arguments
    /// * `m` - Current magnetization (A/m)
    /// * `dt` - Time step (s)
    ///
    /// # Returns
    /// Updated magnetization (A/m)
    pub fn step_rk4(&self, m: Vector3<f64>, dt: f64) -> Vector3<f64> {
        let k1 = self.dmdt(m);
        let k2 = self.dmdt(m + k1 * (dt / 2.0));
        let k3 = self.dmdt(m + k2 * (dt / 2.0));
        let k4 = self.dmdt(m + k3 * dt);

        let m_new = m + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);

        // Renormalize
        let mag = m_new.magnitude();
        if mag > 1e-10 {
            m_new * (self.material.ms / mag)
        } else {
            m
        }
    }

    /// Solve LLG equation for multiple time steps
    ///
    /// # Arguments
    /// * `m0` - Initial magnetization (A/m)
    /// * `dt` - Time step (s)
    /// * `n_steps` - Number of steps
    ///
    /// # Returns
    /// Vector of magnetization at each time step
    pub fn solve(&self, m0: Vector3<f64>, dt: f64, n_steps: usize) -> Vec<Vector3<f64>> {
        let mut trajectory = Vec::with_capacity(n_steps + 1);
        trajectory.push(m0);

        let mut m = m0;
        for _ in 0..n_steps {
            m = self.step_rk4(m, dt);
            trajectory.push(m);
        }

        trajectory
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_llg_solver_creation() {
        let material = Ferromagnet::permalloy();
        let solver = LLGSolver::new(material);
        assert_eq!(solver.h_ext.magnitude(), 0.0);
    }

    #[test]
    fn test_effective_field() {
        let material = Ferromagnet::permalloy();
        let ms = material.ms;
        let mut solver = LLGSolver::new(material);
        solver.h_ext = Vector3::new(1000.0, 0.0, 0.0);

        let m = Vector3::new(0.0, 0.0, ms);
        let h_eff = solver.effective_field(m);

        // Should be dominated by external field
        assert!((h_eff.x - 1000.0).abs() < 1e-6);
    }

    #[test]
    fn test_dmdt_perpendicular() {
        let material = Ferromagnet::permalloy();
        let ms = material.ms;
        let mut solver = LLGSolver::new(material);

        // Apply field in x, magnetization in z
        solver.h_ext = Vector3::new(1000.0, 0.0, 0.0);
        let m = Vector3::new(0.0, 0.0, ms);

        let dmdt = solver.dmdt(m);

        // dm/dt should be primarily in y-direction (m × H)
        // The precession term dominates, creating motion in y
        assert!(dmdt.y.abs() > 1e6); // Should be large
                                     // x-component may be small but non-zero due to damping
        assert!(dmdt.x.abs() < dmdt.y.abs() * 0.1); // Much smaller than y
    }

    #[test]
    fn test_step_euler_normalization() {
        let material = Ferromagnet::permalloy();
        let ms = material.ms;
        let mut solver = LLGSolver::new(material);
        solver.h_ext = Vector3::new(1000.0, 0.0, 0.0);

        let m0 = Vector3::new(0.0, 0.0, ms);
        let m1 = solver.step_euler(m0, 1e-12);

        // Magnitude should be preserved
        assert!((m1.magnitude() - ms).abs() < 1e-6);
    }

    #[test]
    fn test_step_rk4_normalization() {
        let material = Ferromagnet::permalloy();
        let ms = material.ms;
        let mut solver = LLGSolver::new(material);
        solver.h_ext = Vector3::new(1000.0, 0.0, 0.0);

        let m0 = Vector3::new(0.0, 0.0, ms);
        let m1 = solver.step_rk4(m0, 1e-12);

        // Magnitude should be preserved
        assert!((m1.magnitude() - ms).abs() < 1e-6);
    }

    #[test]
    fn test_precession() {
        let material = Ferromagnet::permalloy();
        let ms = material.ms;
        let mut solver = LLGSolver::new(material);

        // Field in z, magnetization slightly tilted
        solver.h_ext = Vector3::new(0.0, 0.0, 10000.0);
        let m0 = Vector3::new(0.1, 0.0, 1.0).normalize() * ms;

        // Run for a few steps
        let dt = 1e-12; // 1 ps
        let m1 = solver.step_rk4(m0, dt);

        // Magnetization should have precessed (y-component should appear)
        assert!(m1.y.abs() > 1e-10);
    }

    #[test]
    fn test_solve_trajectory() {
        let material = Ferromagnet::permalloy();
        let ms = material.ms;
        let mut solver = LLGSolver::new(material);
        solver.h_ext = Vector3::new(1000.0, 0.0, 0.0);

        let m0 = Vector3::new(0.0, 0.0, ms);
        let trajectory = solver.solve(m0, 1e-12, 10);

        assert_eq!(trajectory.len(), 11); // Initial + 10 steps

        // All magnetizations should have same magnitude
        for m in &trajectory {
            assert!((m.magnitude() - ms).abs() < 1e-6);
        }
    }
}
