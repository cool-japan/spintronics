//! Python bindings for spin caloritronics.
//!
//! Exposes the Onsager transport matrix and spin-caloritronic cross-effects
//! to Python, allowing computation of spin Seebeck, spin Peltier, anomalous
//! Nernst, and spin Nernst responses.
//!
//! ## Python Usage
//!
//! ```python
//! import spintronics
//!
//! # Onsager matrix for YIG/Pt at 300 K
//! mat = spintronics.OnsagerMatrix.yig_pt(300.0)
//! grad_t  = (1000.0, 0.0, 0.0)   # 1000 K/m along x
//! e_field = (0.0, 0.0, 0.0)
//! currents = mat.all_currents(grad_t, e_field)
//! print(f"Spin current: {currents['spin_current']}")
//!
//! # Full spin caloritronics computation
//! scm = spintronics.SpinCaloritronicsMaterial.yig_pt(300.0)
//! result = scm.compute_all(grad_t=(1000.0, 0.0, 0.0), j_spin=(0.0, 0.0, 1e3))
//! print(f"Nernst voltage: {result['nernst_voltage']:.4e} V/m")
//! ```

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::caloritronics::cross_effects::SpinCaloritronicsMaterial;
use crate::caloritronics::onsager::OnsagerMatrix;
use crate::vector3::Vector3;

// ─────────────────────────────────────────────────────────────────────────────
// Helper: convert (f64, f64, f64) tuple → Vector3
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn arr_to_vec(a: [f64; 3]) -> Vector3<f64> {
    Vector3::new(a[0], a[1], a[2])
}

#[inline]
fn vec_to_arr(v: &Vector3<f64>) -> [f64; 3] {
    [v.x, v.y, v.z]
}

// ─────────────────────────────────────────────────────────────────────────────
// PyOnsagerMatrix
// ─────────────────────────────────────────────────────────────────────────────

/// Onsager transport matrix for spin caloritronics.
///
/// Encodes the linear-response coupling between charge, spin, and heat
/// currents in a magnetic heterostructure using Onsager reciprocal relations.
///
/// ## Preset Systems
///
/// ```python
/// yig_pt  = OnsagerMatrix.yig_pt(300.0)   # YIG/Pt bilayer at 300 K
/// fe_pt   = OnsagerMatrix.fe_pt(300.0)    # Fe/Pt bilayer at 300 K
/// cofeb_pt = OnsagerMatrix.cofeb_pt(300.0) # CoFeB/Pt bilayer
/// ```
///
/// ## Key Coefficients
///
/// - `conductivity`: electrical conductivity σ [S/m]
/// - `seebeck`: Seebeck coefficient S_e [V/K] (usually negative for metals)
/// - `spin_seebeck`: spin Seebeck coefficient S_s [A/(m·K)]
/// - `hall_angle`: anomalous Hall angle θ_H (dimensionless)
/// - `thermal_conductivity`: κ [W/(m·K)]
#[pyclass(name = "OnsagerMatrix")]
#[derive(Clone)]
pub struct PyOnsagerMatrix {
    inner: OnsagerMatrix,
}

#[pymethods]
impl PyOnsagerMatrix {
    /// Create an Onsager matrix with explicit parameters.
    ///
    /// Args:
    ///     temperature: Temperature \[K\]
    ///     conductivity: Electrical conductivity σ [S/m]
    ///     seebeck: Seebeck coefficient S_e [V/K]
    ///     spin_seebeck: Spin Seebeck coefficient S_s [A/(m·K)]
    ///     hall_angle: Anomalous Hall angle θ_H (dimensionless)
    ///     thermal_conductivity: Thermal conductivity κ [W/(m·K)]
    #[new]
    pub fn new(
        temperature: f64,
        conductivity: f64,
        seebeck: f64,
        spin_seebeck: f64,
        hall_angle: f64,
        thermal_conductivity: f64,
    ) -> Self {
        Self {
            inner: OnsagerMatrix::new(
                temperature,
                conductivity,
                seebeck,
                spin_seebeck,
                hall_angle,
                thermal_conductivity,
            ),
        }
    }

    /// YIG/Pt bilayer preset at `temperature` \[K\].
    #[staticmethod]
    pub fn yig_pt(temperature: f64) -> Self {
        Self {
            inner: OnsagerMatrix::yig_pt(temperature),
        }
    }

    /// Fe/Pt bilayer preset at `temperature` \[K\].
    #[staticmethod]
    pub fn fe_pt(temperature: f64) -> Self {
        Self {
            inner: OnsagerMatrix::fe_pt(temperature),
        }
    }

    /// CoFeB/Pt bilayer preset at `temperature` \[K\].
    #[staticmethod]
    pub fn cofeb_pt(temperature: f64) -> Self {
        Self {
            inner: OnsagerMatrix::cofeb_pt(temperature),
        }
    }

    /// Temperature \[K\]
    #[getter]
    pub fn temperature(&self) -> f64 {
        self.inner.temperature
    }

    /// Electrical conductivity σ [S/m]
    #[getter]
    pub fn conductivity(&self) -> f64 {
        self.inner.conductivity
    }

    /// Seebeck coefficient S_e [V/K]
    #[getter]
    pub fn seebeck(&self) -> f64 {
        self.inner.seebeck
    }

    /// Spin Seebeck coefficient S_s [A/(m·K)]
    #[getter]
    pub fn spin_seebeck(&self) -> f64 {
        self.inner.spin_seebeck
    }

    /// Anomalous Hall angle θ_H (dimensionless)
    #[getter]
    pub fn hall_angle(&self) -> f64 {
        self.inner.hall_angle
    }

    /// Thermal conductivity κ [W/(m·K)]
    #[getter]
    pub fn thermal_conductivity(&self) -> f64 {
        self.inner.thermal_conductivity
    }

    /// Check Onsager reciprocity: returns relative deviation from L_ij = T·L_ji.
    ///
    /// For analytically constructed matrices this is zero by construction.
    /// Values < 1e-10 indicate reciprocity is satisfied.
    pub fn reciprocity_error(&self) -> f64 {
        self.inner.reciprocity_error()
    }

    /// Compute spin current density from a temperature gradient (spin Seebeck effect).
    ///
    /// j_s = S_s · σ · ∇T  [A/m²]
    ///
    /// Args:
    ///     grad_t: Temperature gradient vector (∂T/∂x, ∂T/∂y, ∂T/∂z) [K/m]
    ///
    /// Returns:
    ///     Spin current density (jx, jy, jz) [A/m²]
    pub fn spin_current_from_grad_t(&self, grad_t: [f64; 3]) -> [f64; 3] {
        let result = self.inner.spin_current_from_grad_t(&arr_to_vec(grad_t));
        vec_to_arr(&result)
    }

    /// Compute heat current from a spin current (spin Peltier effect).
    ///
    /// j_Q = Π_s · j_s = T · S_s · j_s  [W/m²]
    ///
    /// Args:
    ///     j_spin: Spin current density (jx, jy, jz) [A/m²]
    ///
    /// Returns:
    ///     Heat current density (jQx, jQy, jQz) [W/m²]
    pub fn heat_current_from_spin_current(&self, j_spin: [f64; 3]) -> [f64; 3] {
        let result = self
            .inner
            .heat_current_from_spin_current(&arr_to_vec(j_spin));
        vec_to_arr(&result)
    }

    /// Compute the anomalous Nernst voltage per unit length.
    ///
    /// ν = -S_e · θ_H · |∇T|  [V/m]
    ///
    /// Args:
    ///     grad_t: Longitudinal temperature gradient magnitude [K/m]
    ///
    /// Returns:
    ///     Transverse Nernst voltage [V/m]
    pub fn nernst_voltage(&self, grad_t: f64) -> f64 {
        self.inner.nernst_voltage(grad_t)
    }

    /// Compute charge, spin, and heat currents simultaneously.
    ///
    /// Full Onsager response:
    ///   j_c = σ·(E + S_e·∇T)      [A/m²]
    ///   j_s = S_s·σ·∇T            [A/m²]
    ///   j_Q = T·S_e·σ·E - κ·∇T   [W/m²]
    ///
    /// Args:
    ///     grad_t: Temperature gradient (∂T/∂x, ∂T/∂y, ∂T/∂z) [K/m]
    ///     e_field: Electric field (Ex, Ey, Ez) [V/m]
    ///
    /// Returns:
    ///     dict with keys 'charge_current', 'spin_current', 'heat_current'
    ///     (each is a list of 3 floats)
    pub fn all_currents(
        &self,
        py: Python<'_>,
        grad_t: [f64; 3],
        e_field: [f64; 3],
    ) -> PyResult<PyObject> {
        let currents = self
            .inner
            .all_currents(&arr_to_vec(grad_t), &arr_to_vec(e_field));

        let dict = PyDict::new(py);
        dict.set_item(
            "charge_current",
            vec_to_arr(&currents.charge_current).to_vec(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item("spin_current", vec_to_arr(&currents.spin_current).to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item("heat_current", vec_to_arr(&currents.heat_current).to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(dict.into())
    }

    /// Python repr string for OnsagerMatrix.
    pub fn __repr__(&self) -> String {
        format!(
            "OnsagerMatrix(T={:.1} K, σ={:.2e} S/m, S_e={:.2e} V/K, S_s={:.2e} A/(m·K))",
            self.inner.temperature,
            self.inner.conductivity,
            self.inner.seebeck,
            self.inner.spin_seebeck
        )
    }
}

impl PyOnsagerMatrix {
    /// Access the inner Rust struct
    pub fn inner(&self) -> &OnsagerMatrix {
        &self.inner
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PySpinCaloritronicsMaterial
// ─────────────────────────────────────────────────────────────────────────────

/// Unified spin caloritronic material: combines Onsager matrix with heat
/// current calculator for a single high-level computation.
///
/// This is the primary user-facing type for spin caloritronics calculations.
/// Use `compute_all()` to obtain all cross-effects simultaneously.
///
/// ## Preset Systems
///
/// ```python
/// yig_pt   = SpinCaloritronicsMaterial.yig_pt(300.0)
/// fe_pt    = SpinCaloritronicsMaterial.fe_pt(300.0)
/// cofeb_pt = SpinCaloritronicsMaterial.cofeb_pt(300.0)
/// ```
///
/// ## Example
///
/// ```python
/// mat = SpinCaloritronicsMaterial.fe_pt(400.0)
/// result = mat.compute_all(grad_t=(500.0, 0.0, 0.0), j_spin=(0.0, 0.0, 1e4))
/// print(f"Nernst voltage: {result['nernst_voltage']:.4e} V/m")
/// print(f"Spin Seebeck current: {result['spin_seebeck_current']}")
/// ```
#[pyclass(name = "SpinCaloritronicsMaterial")]
#[derive(Clone)]
pub struct PySpinCaloritronicsMaterial {
    inner: SpinCaloritronicsMaterial,
}

#[pymethods]
impl PySpinCaloritronicsMaterial {
    /// YIG/Pt bilayer preset at `temperature` \[K\].
    #[staticmethod]
    pub fn yig_pt(temperature: f64) -> Self {
        Self {
            inner: SpinCaloritronicsMaterial::yig_pt(temperature),
        }
    }

    /// Fe/Pt bilayer preset at `temperature` \[K\].
    #[staticmethod]
    pub fn fe_pt(temperature: f64) -> Self {
        Self {
            inner: SpinCaloritronicsMaterial::fe_pt(temperature),
        }
    }

    /// CoFeB/Pt bilayer preset at `temperature` \[K\].
    #[staticmethod]
    pub fn cofeb_pt(temperature: f64) -> Self {
        Self {
            inner: SpinCaloritronicsMaterial::cofeb_pt(temperature),
        }
    }

    /// Create a material from an Onsager matrix.
    ///
    /// The HeatCurrentCalculator is built automatically via Kelvin relations:
    ///   Π   = T · S_e
    ///   Π_s = T · S_s
    ///
    /// Args:
    ///     onsager: OnsagerMatrix for the desired material system
    #[staticmethod]
    pub fn from_onsager(onsager: &PyOnsagerMatrix) -> Self {
        Self {
            inner: SpinCaloritronicsMaterial::new(onsager.inner.clone()),
        }
    }

    /// Human-readable material system label (e.g. "YIG/Pt", "Fe/Pt").
    #[getter]
    pub fn name(&self) -> &str {
        &self.inner.name
    }

    /// The underlying OnsagerMatrix as a Python object.
    #[getter]
    pub fn onsager(&self) -> PyOnsagerMatrix {
        PyOnsagerMatrix {
            inner: self.inner.onsager.clone(),
        }
    }

    /// Compute all spin-caloritronic cross-effects.
    ///
    /// Evaluates:
    ///   - Spin Seebeck current: j_s = S_s·σ·∇T
    ///   - Spin Peltier heat: |Π_s·j_s|
    ///   - Anomalous Nernst voltage: ν = -S_e·θ_H·|∇T|
    ///   - Spin Nernst current (proxy): θ_H · j_s^SSE
    ///   - Total heat current: −κ·∇T + Π_s·j_s
    ///   - Onsager reciprocity check
    ///
    /// Args:
    ///     grad_t: Temperature gradient (∂T/∂x, ∂T/∂y, ∂T/∂z) [K/m]
    ///     j_spin: Injected spin current density (jx, jy, jz) [A/m²]
    ///
    /// Returns:
    ///     dict with keys:
    ///       'spin_seebeck_current': list\[float\] — j_s [A/m²]
    ///       'peltier_heat': float — |j_Q^sPeltier| [W/m²]
    ///       'nernst_voltage': float — anomalous Nernst ν [V/m]
    ///       'spin_nernst_current': list\[float\] — j_s^SN proxy [A/m²]
    ///       'reciprocity_satisfied': bool — Onsager error < 1e-10
    ///       'total_heat_current': list\[float\] — total j_Q [W/m²]
    pub fn compute_all(
        &self,
        py: Python<'_>,
        grad_t: [f64; 3],
        j_spin: [f64; 3],
    ) -> PyResult<PyObject> {
        let result = self
            .inner
            .compute_all(&arr_to_vec(grad_t), &arr_to_vec(j_spin))
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let dict = PyDict::new(py);
        dict.set_item(
            "spin_seebeck_current",
            vec_to_arr(&result.spin_seebeck_current).to_vec(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item("peltier_heat", result.peltier_heat)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item("nernst_voltage", result.nernst_voltage)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item(
            "spin_nernst_current",
            vec_to_arr(&result.spin_nernst_current).to_vec(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item("reciprocity_satisfied", result.reciprocity_satisfied)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        dict.set_item(
            "total_heat_current",
            vec_to_arr(&result.total_heat_current).to_vec(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(dict.into())
    }

    /// Python repr string for SpinCaloritronicsMaterial.
    pub fn __repr__(&self) -> String {
        format!(
            "SpinCaloritronicsMaterial(name='{}', T={:.1} K)",
            self.inner.name, self.inner.onsager.temperature
        )
    }
}
