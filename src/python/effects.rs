//! Python bindings for spin-charge conversion effects

use pyo3::prelude::*;

use super::vector::PyVector3;
use crate::effect::InverseSpinHall;

/// Inverse Spin Hall Effect converter
///
/// Converts spin current to charge current via spin-orbit coupling.
/// The electric field E is generated perpendicular to both the
/// spin current direction and spin polarization:
///
/// E = rho * theta_SH * (j_s x sigma)
///
/// ## Example
/// ```python
/// pt = InverseSpinHall.platinum()
/// js_flow = Vector3(0.0, 0.0, 1.0)  # Spin current flow direction
/// js_pol = Vector3(1e8, 0.0, 0.0)   # Spin polarization (A/m^2)
///
/// e_field = pt.convert(js_flow, js_pol)
/// print(f"Electric field: {e_field} V/m")
/// ```
#[pyclass(name = "InverseSpinHall")]
#[derive(Clone)]
pub struct PyInverseSpinHall {
    inner: InverseSpinHall,
}

#[pymethods]
impl PyInverseSpinHall {
    /// Create a new ISHE converter
    ///
    /// Args:
    ///     theta_sh: Spin Hall angle (dimensionless)
    ///     rho: Electrical resistivity (Ohm*m)
    #[new]
    pub fn new(theta_sh: f64, rho: f64) -> Self {
        Self {
            inner: InverseSpinHall { theta_sh, rho },
        }
    }

    /// Spin Hall angle (dimensionless)
    #[getter]
    pub fn theta_sh(&self) -> f64 {
        self.inner.theta_sh
    }

    /// Electrical resistivity (Ohm*m)
    #[getter]
    pub fn rho(&self) -> f64 {
        self.inner.rho
    }

    /// Convert spin current to electric field
    ///
    /// Args:
    ///     js_flow: Spin current flow direction (Vector3)
    ///     js_polarization: Spin polarization direction with magnitude (Vector3, A/m^2)
    ///
    /// Returns:
    ///     Electric field vector (V/m)
    pub fn convert(&self, js_flow: &PyVector3, js_polarization: &PyVector3) -> PyVector3 {
        let e_field = self.inner.convert(js_flow.inner(), js_polarization.inner());
        PyVector3::from_inner(e_field)
    }

    /// Calculate voltage across a strip of given width
    ///
    /// Args:
    ///     js_flow: Spin current flow direction (Vector3)
    ///     js_polarization: Spin polarization direction with magnitude (Vector3, A/m^2)
    ///     strip_width: Width of the conducting strip (m)
    ///
    /// Returns:
    ///     Voltage (V)
    pub fn voltage(
        &self,
        js_flow: &PyVector3,
        js_polarization: &PyVector3,
        strip_width: f64,
    ) -> f64 {
        self.inner
            .voltage(js_flow.inner(), js_polarization.inner(), strip_width)
    }

    /// Calculate conversion efficiency (V/W ratio)
    pub fn efficiency(&self) -> f64 {
        self.inner.efficiency()
    }

    /// Create Platinum (Pt) ISHE converter
    ///
    /// Platinum has a large positive spin Hall angle (~0.08)
    /// and is the most commonly used material for ISHE detection.
    #[staticmethod]
    pub fn platinum() -> Self {
        Self {
            inner: InverseSpinHall::platinum(),
        }
    }

    /// Create Tantalum (Ta) ISHE converter
    ///
    /// Tantalum has a large positive spin Hall angle (~0.12).
    #[staticmethod]
    pub fn tantalum() -> Self {
        Self {
            inner: InverseSpinHall::tantalum(),
        }
    }

    /// Create Tungsten (W) ISHE converter
    ///
    /// Tungsten has the largest known spin Hall angle (~-0.30, negative).
    #[staticmethod]
    pub fn tungsten() -> Self {
        Self {
            inner: InverseSpinHall::tungsten(),
        }
    }

    /// String representation
    pub fn __repr__(&self) -> String {
        format!(
            "InverseSpinHall(theta_sh={:.3}, rho={:.2e})",
            self.inner.theta_sh, self.inner.rho
        )
    }
}

impl PyInverseSpinHall {
    /// Get the inner InverseSpinHall
    #[allow(dead_code)]
    pub fn inner(&self) -> &InverseSpinHall {
        &self.inner
    }
}
