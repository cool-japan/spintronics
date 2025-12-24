//! Python bindings for magnetization dynamics

use pyo3::prelude::*;

use super::materials::PyFerromagnet;
use super::vector::PyVector3;
use crate::constants::GAMMA;
use crate::dynamics::calc_dm_dt;
use crate::material::Ferromagnet;
use crate::vector3::Vector3;

/// LLG (Landau-Lifshitz-Gilbert) equation simulator
///
/// Simulates magnetization dynamics under external fields and torques.
/// The LLG equation describes the time evolution of the magnetization:
///
/// dm/dt = -γ(m × H_eff) + α(m × dm/dt)
///
/// ## Example
/// ```python
/// yig = Ferromagnet.yig()
/// sim = LlgSimulator(yig)
/// sim.set_magnetization(1.0, 0.0, 0.0)
/// sim.set_external_field(0.0, 0.0, 0.1)  # 100 mT in z
///
/// # Evolve for 1 ns with 1000 steps
/// trajectory = sim.evolve(1e-9, 1000)
/// ```
#[pyclass(name = "LlgSimulator")]
pub struct PyLlgSimulator {
    material: Ferromagnet,
    magnetization: Vector3<f64>,
    external_field: Vector3<f64>,
    time: f64,
}

#[pymethods]
impl PyLlgSimulator {
    /// Create a new LLG simulator
    ///
    /// Args:
    ///     material: Ferromagnetic material parameters
    #[new]
    pub fn new(material: &PyFerromagnet) -> Self {
        Self {
            material: material.inner().clone(),
            magnetization: Vector3::new(1.0, 0.0, 0.0),
            external_field: Vector3::new(0.0, 0.0, 0.0),
            time: 0.0,
        }
    }

    /// Set magnetization direction (will be normalized)
    ///
    /// Args:
    ///     mx, my, mz: Magnetization components
    pub fn set_magnetization(&mut self, mx: f64, my: f64, mz: f64) {
        self.magnetization = Vector3::new(mx, my, mz).normalize();
    }

    /// Get current magnetization
    pub fn get_magnetization(&self) -> PyVector3 {
        PyVector3::from_inner(self.magnetization)
    }

    /// Set external magnetic field (Tesla)
    ///
    /// Args:
    ///     hx, hy, hz: Field components in Tesla
    pub fn set_external_field(&mut self, hx: f64, hy: f64, hz: f64) {
        self.external_field = Vector3::new(hx, hy, hz);
    }

    /// Get current external field
    pub fn get_external_field(&self) -> PyVector3 {
        PyVector3::from_inner(self.external_field)
    }

    /// Get current simulation time (seconds)
    #[getter]
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Reset simulation time to zero
    pub fn reset_time(&mut self) {
        self.time = 0.0;
    }

    /// Calculate dm/dt at current state
    pub fn dm_dt(&self) -> PyVector3 {
        let dm_dt = calc_dm_dt(
            self.magnetization,
            self.external_field,
            GAMMA,
            self.material.alpha,
        );
        PyVector3::from_inner(dm_dt)
    }

    /// Perform one RK4 time step
    ///
    /// Args:
    ///     dt: Time step (seconds)
    pub fn step_rk4(&mut self, dt: f64) {
        let m = self.magnetization;
        let h = self.external_field;
        let alpha = self.material.alpha;

        // RK4 integration
        let k1 = calc_dm_dt(m, h, GAMMA, alpha);
        let k2 = calc_dm_dt((m + k1 * (dt / 2.0)).normalize(), h, GAMMA, alpha);
        let k3 = calc_dm_dt((m + k2 * (dt / 2.0)).normalize(), h, GAMMA, alpha);
        let k4 = calc_dm_dt((m + k3 * dt).normalize(), h, GAMMA, alpha);

        self.magnetization = (m + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0)).normalize();
        self.time += dt;
    }

    /// Perform one Euler time step
    ///
    /// Args:
    ///     dt: Time step (seconds)
    pub fn step_euler(&mut self, dt: f64) {
        let dm_dt = calc_dm_dt(
            self.magnetization,
            self.external_field,
            GAMMA,
            self.material.alpha,
        );
        self.magnetization = (self.magnetization + dm_dt * dt).normalize();
        self.time += dt;
    }

    /// Evolve magnetization for specified duration
    ///
    /// Args:
    ///     duration: Total simulation time (seconds)
    ///     n_steps: Number of time steps
    ///     method: Integration method ("rk4" or "euler"), default "rk4"
    ///
    /// Returns:
    ///     List of (time, mx, my, mz) tuples
    #[pyo3(signature = (duration, n_steps, method = "rk4"))]
    pub fn evolve(
        &mut self,
        duration: f64,
        n_steps: usize,
        method: &str,
    ) -> Vec<(f64, f64, f64, f64)> {
        let dt = duration / n_steps as f64;
        let mut trajectory = Vec::with_capacity(n_steps + 1);

        // Record initial state
        trajectory.push((
            self.time,
            self.magnetization.x,
            self.magnetization.y,
            self.magnetization.z,
        ));

        // Evolve
        for _ in 0..n_steps {
            match method {
                "euler" => self.step_euler(dt),
                _ => self.step_rk4(dt),
            }
            trajectory.push((
                self.time,
                self.magnetization.x,
                self.magnetization.y,
                self.magnetization.z,
            ));
        }

        trajectory
    }

    /// Calculate precession frequency (rad/s)
    ///
    /// Returns the Larmor frequency for the current field.
    pub fn precession_frequency(&self) -> f64 {
        GAMMA * self.external_field.magnitude()
    }

    /// Calculate precession period (seconds)
    pub fn precession_period(&self) -> f64 {
        let omega = self.precession_frequency();
        if omega > 0.0 {
            2.0 * std::f64::consts::PI / omega
        } else {
            f64::INFINITY
        }
    }

    /// String representation
    pub fn __repr__(&self) -> String {
        format!(
            "LlgSimulator(m=({:.3}, {:.3}, {:.3}), H=({:.3}, {:.3}, {:.3}), t={:.2e}s)",
            self.magnetization.x,
            self.magnetization.y,
            self.magnetization.z,
            self.external_field.x,
            self.external_field.y,
            self.external_field.z,
            self.time
        )
    }
}
