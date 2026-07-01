//! Python bindings for Landau-Lifshitz-Bloch (LLB) dynamics.
//!
//! The LLB equation extends the LLG equation to finite temperatures,
//! including near and above the Curie temperature. Unlike LLG, the
//! magnitude |m| is NOT conserved — longitudinal relaxation is included.
//!
//! ## Python Usage
//!
//! ```python
//! import spintronics
//!
//! mat = spintronics.LlbMaterial.iron()
//! print(f"Curie temperature: {mat.curie_temp} K")
//! print(f"m_e(300 K) = {mat.equilibrium_magnetization(300.0):.4f}")
//!
//! solver = spintronics.LlbSolver(mat, dt=1e-14, temperature=300.0,
//!                                 h_ext=(0.0, 0.0, 1.0))
//! result = solver.run(m0=(0.8, 0.1, 0.0), num_steps=1000, record_every=10)
//! print(result['equilibrium_m'])
//! ```

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::dynamics::llb::{LlbMaterial, LlbResult, LlbSolver};
use crate::vector3::Vector3;

// ─────────────────────────────────────────────────────────────────────────────
// PyLlbMaterial
// ─────────────────────────────────────────────────────────────────────────────

/// Material parameters for the LLB (Landau-Lifshitz-Bloch) equation.
///
/// Encapsulates Curie temperature, Gilbert damping, quantum spin number,
/// and zero-temperature saturation magnetization needed for finite-temperature
/// spin dynamics.
///
/// ## Preset Materials
///
/// ```python
/// iron  = LlbMaterial.iron()    # T_C = 1043 K, α = 0.01
/// ni    = LlbMaterial.nickel()  # T_C = 631  K, α = 0.064
/// cofeb = LlbMaterial.cofeb()   # T_C = 1000 K, α = 0.005
/// ```
#[pyclass(name = "LlbMaterial", skip_from_py_object)]
#[derive(Clone)]
pub struct PyLlbMaterial {
    inner: LlbMaterial,
}

#[pymethods]
impl PyLlbMaterial {
    /// Create a new LLB material with explicit parameters.
    ///
    /// Args:
    ///     curie_temp: Curie temperature T_C \[K\]
    ///     alpha: Base Gilbert damping constant (dimensionless)
    ///     spin_s: Quantum spin number S (e.g. 0.5 for spin-1/2)
    ///     ms_0: Zero-temperature saturation magnetization [A/m]
    #[new]
    pub fn new(curie_temp: f64, alpha: f64, spin_s: f64, ms_0: f64) -> Self {
        Self {
            inner: LlbMaterial {
                curie_temp,
                alpha,
                spin_s,
                ms_0,
            },
        }
    }

    /// Iron (Fe) LLB material preset (T_C = 1043 K, α = 0.01, S = 1.0).
    #[staticmethod]
    pub fn iron() -> Self {
        Self {
            inner: LlbMaterial::iron(),
        }
    }

    /// Nickel (Ni) LLB material preset (T_C = 631 K, α = 0.064, S = 0.3).
    #[staticmethod]
    pub fn nickel() -> Self {
        Self {
            inner: LlbMaterial::nickel(),
        }
    }

    /// CoFeB LLB material preset (T_C = 1000 K, α = 0.005, S = 0.5).
    #[staticmethod]
    pub fn cofeb() -> Self {
        Self {
            inner: LlbMaterial::cofeb(),
        }
    }

    /// Compute equilibrium magnetization m_e(T) via self-consistent Brillouin equation.
    ///
    /// Returns 0.0 for T >= T_C or T <= 0.
    ///
    /// Args:
    ///     temperature: Temperature \[K\]
    ///
    /// Returns:
    ///     Dimensionless equilibrium magnetization in [0, 1]
    pub fn equilibrium_magnetization(&self, temperature: f64) -> f64 {
        self.inner.equilibrium_magnetization(temperature)
    }

    /// Temperature-dependent longitudinal damping α_∥(T).
    ///
    /// - T < T_C:  α_∥ = α · (2/5 + 3T/(5T_C))
    /// - T ≥ T_C:  α_∥ = 2α · T/(5T_C)
    ///
    /// Args:
    ///     temperature: Temperature \[K\]
    pub fn alpha_parallel(&self, temperature: f64) -> f64 {
        self.inner.alpha_parallel(temperature)
    }

    /// Temperature-dependent transverse damping α_⊥(T).
    ///
    /// - T < T_C:  α_⊥ = α · T/T_C  (minimum 1e-10)
    /// - T ≥ T_C:  α_⊥ = 2α · T/(5T_C)
    ///
    /// Args:
    ///     temperature: Temperature \[K\]
    pub fn alpha_perp(&self, temperature: f64) -> f64 {
        self.inner.alpha_perp(temperature)
    }

    /// Curie temperature T_C \[K\]
    #[getter]
    pub fn curie_temp(&self) -> f64 {
        self.inner.curie_temp
    }

    /// Base Gilbert damping constant α (dimensionless)
    #[getter]
    pub fn alpha(&self) -> f64 {
        self.inner.alpha
    }

    /// Quantum spin number S
    #[getter]
    pub fn spin_s(&self) -> f64 {
        self.inner.spin_s
    }

    /// Zero-temperature saturation magnetization M_s [A/m]
    #[getter]
    pub fn ms_0(&self) -> f64 {
        self.inner.ms_0
    }

    /// Python repr string for LlbMaterial.
    pub fn __repr__(&self) -> String {
        format!(
            "LlbMaterial(curie_temp={:.1}, alpha={:.4}, spin_s={:.2}, ms_0={:.3e})",
            self.inner.curie_temp, self.inner.alpha, self.inner.spin_s, self.inner.ms_0
        )
    }
}

impl PyLlbMaterial {
    /// Access the inner Rust struct
    pub fn inner(&self) -> &LlbMaterial {
        &self.inner
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyLlbSolver
// ─────────────────────────────────────────────────────────────────────────────

/// LLB equation solver with 4th-order Runge-Kutta integration.
///
/// Unlike the standard LLG solver, the LLB solver allows |m| to change
/// over time via longitudinal relaxation. This is essential near and above
/// the Curie temperature T_C.
///
/// ## Example
///
/// ```python
/// mat    = LlbMaterial.iron()
/// solver = LlbSolver(mat, dt=1e-14, temperature=300.0, h_ext=(0.0, 0.0, 1.0))
/// result = solver.run(m0=(0.8, 0.0, 0.0), num_steps=5000, record_every=50)
///
/// import numpy as np
/// import matplotlib.pyplot as plt
/// plt.plot(result['time'], result['m_magnitude'])
/// plt.xlabel("Time [s]")
/// plt.ylabel("|m|")
/// plt.show()
/// ```
#[pyclass(name = "LlbSolver")]
pub struct PyLlbSolver {
    inner: LlbSolver,
}

#[pymethods]
impl PyLlbSolver {
    /// Create a new LLB solver.
    ///
    /// Args:
    ///     material: LlbMaterial with Curie temperature, damping, spin number, M_s
    ///     dt: Integration time step \[s\]
    ///     temperature: Simulation temperature \[K\]
    ///     h_ext: External applied field (hx, hy, hz) \[T\]
    #[new]
    pub fn new(material: &PyLlbMaterial, dt: f64, temperature: f64, h_ext: [f64; 3]) -> Self {
        let h = Vector3::new(h_ext[0], h_ext[1], h_ext[2]);
        Self {
            inner: LlbSolver::new(material.inner.clone(), dt, temperature, h),
        }
    }

    /// Advance the magnetization by one RK4 time step.
    ///
    /// Note: Unlike LLG, |m| is NOT renormalized after each step. The LLB
    /// equation explicitly evolves the magnetization magnitude.
    ///
    /// Args:
    ///     m: Current magnetization vector (mx, my, mz)
    ///
    /// Returns:
    ///     Updated magnetization vector (mx, my, mz) after one dt
    pub fn step(&self, m: [f64; 3]) -> PyResult<[f64; 3]> {
        let m_vec = Vector3::new(m[0], m[1], m[2]);
        let m_new = self
            .inner
            .step(&m_vec)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok([m_new.x, m_new.y, m_new.z])
    }

    /// Run the LLB simulation for num_steps time steps.
    ///
    /// Args:
    ///     m0: Initial magnetization vector (mx, my, mz)
    ///     num_steps: Total number of integration steps
    ///     record_every: Snapshot interval (1 = every step, N = every Nth step)
    ///
    /// Returns:
    ///     dict with keys:
    ///       'mx', 'my', 'mz': list of float — trajectory components
    ///       'm_magnitude': list of float — |m(t)| at each snapshot
    ///       'time': list of float — time stamps \[s\]
    ///       'equilibrium_m': float — m_e(T) for the simulation temperature
    pub fn run<'py>(
        &self,
        py: Python<'py>,
        m0: [f64; 3],
        num_steps: usize,
        record_every: usize,
    ) -> PyResult<Bound<'py, PyDict>> {
        let m0_vec = Vector3::new(m0[0], m0[1], m0[2]);
        let result: LlbResult = self
            .inner
            .run(m0_vec, num_steps, record_every)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let mx: Vec<f64> = result.trajectory.iter().map(|v| v.x).collect();
        let my: Vec<f64> = result.trajectory.iter().map(|v| v.y).collect();
        let mz: Vec<f64> = result.trajectory.iter().map(|v| v.z).collect();

        let dict = PyDict::new(py);
        dict.set_item("mx", mx)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item("my", my)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item("mz", mz)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item("m_magnitude", result.m_magnitude)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item("time", result.time)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item("equilibrium_m", result.equilibrium_m)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(dict)
    }

    /// Integration time step dt \[s\]
    #[getter]
    pub fn dt(&self) -> f64 {
        self.inner.dt
    }

    /// Simulation temperature \[K\]
    #[getter]
    pub fn temperature(&self) -> f64 {
        self.inner.temperature
    }

    /// External field h_ext as (hx, hy, hz) tuple \[T\]
    #[getter]
    pub fn h_ext(&self) -> (f64, f64, f64) {
        (self.inner.h_ext.x, self.inner.h_ext.y, self.inner.h_ext.z)
    }

    /// Gyromagnetic ratio γ [rad/(s·T)]
    #[getter]
    pub fn gamma(&self) -> f64 {
        self.inner.gamma
    }

    /// Python repr string for LlbSolver.
    pub fn __repr__(&self) -> String {
        format!(
            "LlbSolver(T={:.1} K, dt={:.2e} s, H=({:.3}, {:.3}, {:.3}) T)",
            self.inner.temperature,
            self.inner.dt,
            self.inner.h_ext.x,
            self.inner.h_ext.y,
            self.inner.h_ext.z
        )
    }
}
