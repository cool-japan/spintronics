//! Python bindings for high-level simulation workflows

use pyo3::prelude::*;

use super::vector::PyVector3;
use crate::constants::GAMMA;
use crate::dynamics::calc_dm_dt;
use crate::effect::InverseSpinHall;
use crate::material::{Ferromagnet, SpinInterface};
use crate::transport::spin_pumping_current;
use crate::vector3::Vector3;

/// Complete spin pumping simulation
///
/// Simulates the full spin pumping workflow:
/// 1. Magnetization precession under FMR
/// 2. Spin pumping current generation
/// 3. ISHE voltage detection
///
/// ## Example
/// ```python
/// sim = SpinPumpingSimulation()
/// sim.set_fmr_conditions(9.65e9, 0.1)  # 9.65 GHz, 100 mT
/// result = sim.run(1e-9, 1000)
/// print(f"Peak ISHE voltage: {result['peak_voltage']} V")
/// ```
#[pyclass(name = "SpinPumpingSimulation")]
pub struct PySpinPumpingSimulation {
    ferromagnet: Ferromagnet,
    interface: SpinInterface,
    detector: InverseSpinHall,
    magnetization: Vector3<f64>,
    external_field: Vector3<f64>,
    sample_length: f64,
}

#[pymethods]
impl PySpinPumpingSimulation {
    /// Create a new spin pumping simulation with default YIG/Pt system
    #[new]
    pub fn new() -> Self {
        Self {
            ferromagnet: Ferromagnet::yig(),
            interface: SpinInterface::yig_pt(),
            detector: InverseSpinHall::platinum(),
            magnetization: Vector3::new(1.0, 0.0, 0.0),
            external_field: Vector3::new(0.0, 0.0, 0.1),
            sample_length: 5e-3, // 5 mm default
        }
    }

    /// Set sample length for voltage measurement
    pub fn set_sample_length(&mut self, length: f64) {
        self.sample_length = length;
    }

    /// Set external magnetic field (Tesla)
    pub fn set_field(&mut self, hx: f64, hy: f64, hz: f64) {
        self.external_field = Vector3::new(hx, hy, hz);
    }

    /// Set initial magnetization direction
    pub fn set_magnetization(&mut self, mx: f64, my: f64, mz: f64) {
        self.magnetization = Vector3::new(mx, my, mz).normalize();
    }

    /// Set FMR conditions
    ///
    /// Args:
    ///     frequency: Microwave frequency (Hz)
    ///     field: DC magnetic field magnitude (T)
    pub fn set_fmr_conditions(&mut self, _frequency: f64, field: f64) {
        // Set field along z-axis for in-plane FMR geometry
        self.external_field = Vector3::new(0.0, 0.0, field);
        // Initial magnetization in x-y plane
        self.magnetization = Vector3::new(1.0, 0.0, 0.0);
    }

    /// Run the simulation
    ///
    /// Args:
    ///     duration: Total simulation time (seconds)
    ///     n_steps: Number of time steps
    ///
    /// Returns:
    ///     Dictionary with simulation results
    pub fn run(&mut self, duration: f64, n_steps: usize) -> PyResult<PyObject> {
        let dt = duration / n_steps as f64;
        let alpha = self.ferromagnet.alpha;

        let mut times = Vec::with_capacity(n_steps + 1);
        let mut mx_vals = Vec::with_capacity(n_steps + 1);
        let mut my_vals = Vec::with_capacity(n_steps + 1);
        let mut mz_vals = Vec::with_capacity(n_steps + 1);
        let mut js_vals = Vec::with_capacity(n_steps + 1);
        let mut voltage_vals = Vec::with_capacity(n_steps + 1);

        let mut time = 0.0;
        let mut m = self.magnetization;
        let h = self.external_field;

        // Spin current flow direction (from FM to NM, along interface normal)
        let js_flow = self.interface.normal;

        for _ in 0..=n_steps {
            // Record state
            times.push(time);
            mx_vals.push(m.x);
            my_vals.push(m.y);
            mz_vals.push(m.z);

            // Calculate dm/dt
            let dm_dt = calc_dm_dt(m, h, GAMMA, alpha);

            // Calculate spin pumping current
            let js = spin_pumping_current(&self.interface, m, dm_dt);
            let js_mag = js.magnitude();
            js_vals.push(js_mag);

            // Calculate ISHE voltage
            // The spin current flows in the interface normal direction
            // and has polarization given by js vector
            let e_field = self.detector.convert(js_flow, js);
            let voltage = e_field.magnitude() * self.sample_length;
            voltage_vals.push(voltage);

            // RK4 step
            let k1 = calc_dm_dt(m, h, GAMMA, alpha);
            let k2 = calc_dm_dt((m + k1 * (dt / 2.0)).normalize(), h, GAMMA, alpha);
            let k3 = calc_dm_dt((m + k2 * (dt / 2.0)).normalize(), h, GAMMA, alpha);
            let k4 = calc_dm_dt((m + k3 * dt).normalize(), h, GAMMA, alpha);

            m = (m + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0)).normalize();
            time += dt;
        }

        // Calculate statistics
        let peak_voltage = voltage_vals.iter().cloned().fold(0.0_f64, f64::max);
        let avg_voltage: f64 = voltage_vals.iter().sum::<f64>() / voltage_vals.len() as f64;
        let peak_js = js_vals.iter().cloned().fold(0.0_f64, f64::max);

        // Build result dictionary
        Python::with_gil(|py| {
            let dict = pyo3::types::PyDict::new(py);
            dict.set_item("times", times)?;
            dict.set_item("mx", mx_vals)?;
            dict.set_item("my", my_vals)?;
            dict.set_item("mz", mz_vals)?;
            dict.set_item("spin_current", js_vals)?;
            dict.set_item("voltage", voltage_vals)?;
            dict.set_item("peak_voltage", peak_voltage)?;
            dict.set_item("avg_voltage", avg_voltage)?;
            dict.set_item("peak_spin_current", peak_js)?;
            Ok(dict.into())
        })
    }

    /// Get current magnetization
    pub fn get_magnetization(&self) -> PyVector3 {
        PyVector3::from_inner(self.magnetization)
    }

    /// String representation
    pub fn __repr__(&self) -> String {
        format!(
            "SpinPumpingSimulation(material=YIG, detector=Pt, field={:.3} T)",
            self.external_field.magnitude()
        )
    }
}

impl Default for PySpinPumpingSimulation {
    fn default() -> Self {
        Self::new()
    }
}
