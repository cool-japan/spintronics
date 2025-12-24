//! Python bindings for material types

use pyo3::prelude::*;

use super::vector::PyVector3;
use crate::material::{Ferromagnet, SpinInterface};
use crate::vector3::Vector3;

/// Ferromagnetic material properties
///
/// Contains magnetic parameters like saturation magnetization,
/// exchange stiffness, and Gilbert damping.
///
/// ## Example
/// ```python
/// # Use predefined materials
/// yig = Ferromagnet.yig()
/// permalloy = Ferromagnet.permalloy()
/// cofe = Ferromagnet.cofe()
///
/// # Access properties
/// print(f"YIG saturation magnetization: {yig.ms} A/m")
/// print(f"YIG damping: {yig.alpha}")
/// ```
#[pyclass(name = "Ferromagnet")]
#[derive(Clone)]
pub struct PyFerromagnet {
    inner: Ferromagnet,
}

#[pymethods]
impl PyFerromagnet {
    /// Create a new ferromagnet with custom parameters
    ///
    /// Args:
    ///     alpha: Gilbert damping parameter (dimensionless)
    ///     ms: Saturation magnetization (A/m)
    ///     anisotropy_k: Uniaxial anisotropy constant (J/m^3)
    ///     easy_axis: Easy axis direction (tuple of 3 floats)
    ///     exchange_a: Exchange stiffness (J/m)
    #[new]
    #[pyo3(signature = (alpha, ms, anisotropy_k=0.0, easy_axis=(0.0, 0.0, 1.0), exchange_a=1e-11))]
    pub fn new(
        alpha: f64,
        ms: f64,
        anisotropy_k: f64,
        easy_axis: (f64, f64, f64),
        exchange_a: f64,
    ) -> Self {
        Self {
            inner: Ferromagnet {
                alpha,
                ms,
                anisotropy_k,
                easy_axis: Vector3::new(easy_axis.0, easy_axis.1, easy_axis.2).normalize(),
                exchange_a,
            },
        }
    }

    /// Gilbert damping parameter (dimensionless)
    #[getter]
    pub fn alpha(&self) -> f64 {
        self.inner.alpha
    }

    /// Saturation magnetization (A/m)
    #[getter]
    pub fn ms(&self) -> f64 {
        self.inner.ms
    }

    /// Uniaxial anisotropy constant (J/m^3)
    #[getter]
    pub fn anisotropy_k(&self) -> f64 {
        self.inner.anisotropy_k
    }

    /// Easy axis direction (normalized)
    #[getter]
    pub fn easy_axis(&self) -> PyVector3 {
        PyVector3::from_inner(self.inner.easy_axis)
    }

    /// Exchange stiffness (J/m)
    #[getter]
    pub fn exchange_a(&self) -> f64 {
        self.inner.exchange_a
    }

    /// Create YIG (Yttrium Iron Garnet) material
    ///
    /// YIG is a ferrimagnetic insulator with very low damping,
    /// ideal for spin pumping and magnon transport experiments.
    #[staticmethod]
    pub fn yig() -> Self {
        Self {
            inner: Ferromagnet::yig(),
        }
    }

    /// Create Permalloy (Ni80Fe20) material
    ///
    /// Permalloy is a soft magnetic alloy with near-zero
    /// magnetostriction and high permeability.
    #[staticmethod]
    pub fn permalloy() -> Self {
        Self {
            inner: Ferromagnet::permalloy(),
        }
    }

    /// Create CoFe (Cobalt-Iron) alloy
    ///
    /// CoFe has high saturation magnetization, useful for
    /// spin-transfer torque applications.
    #[staticmethod]
    pub fn cofe() -> Self {
        Self {
            inner: Ferromagnet::cofe(),
        }
    }

    /// Create CoFeB (Cobalt-Iron-Boron) alloy
    ///
    /// CoFeB is widely used in magnetic tunnel junctions
    /// due to its perpendicular magnetic anisotropy.
    #[staticmethod]
    pub fn cofeb() -> Self {
        Self {
            inner: Ferromagnet::cofeb(),
        }
    }

    /// Create Iron (Fe) material
    #[staticmethod]
    pub fn iron() -> Self {
        Self {
            inner: Ferromagnet::iron(),
        }
    }

    /// Create Cobalt (Co) material
    #[staticmethod]
    pub fn cobalt() -> Self {
        Self {
            inner: Ferromagnet::cobalt(),
        }
    }

    /// Create Nickel (Ni) material
    #[staticmethod]
    pub fn nickel() -> Self {
        Self {
            inner: Ferromagnet::nickel(),
        }
    }

    /// String representation
    pub fn __repr__(&self) -> String {
        format!(
            "Ferromagnet(alpha={:.4}, ms={:.2e}, anisotropy_k={:.2e}, exchange_a={:.2e})",
            self.inner.alpha, self.inner.ms, self.inner.anisotropy_k, self.inner.exchange_a
        )
    }
}

impl PyFerromagnet {
    /// Get the inner Ferromagnet
    pub fn inner(&self) -> &Ferromagnet {
        &self.inner
    }
}

/// Spin interface between ferromagnet and normal metal
///
/// Contains interface parameters for spin pumping and
/// spin-charge conversion.
///
/// ## Example
/// ```python
/// # Create YIG/Pt interface
/// interface = SpinInterface.yig_pt()
/// print(f"Spin mixing conductance: {interface.g_r} 1/(Ohm*m^2)")
/// ```
#[pyclass(name = "SpinInterface")]
#[derive(Clone)]
pub struct PySpinInterface {
    inner: SpinInterface,
}

#[pymethods]
impl PySpinInterface {
    /// Create a new spin interface
    ///
    /// Args:
    ///     g_r: Real part of spin mixing conductance (1/(Ohm*m^2))
    ///     g_i: Imaginary part of spin mixing conductance (1/(Ohm*m^2))
    ///     normal: Interface normal direction (tuple of 3 floats)
    ///     area: Interface area (m^2)
    #[new]
    #[pyo3(signature = (g_r, g_i=0.0, normal=(0.0, 1.0, 0.0), area=1e-12))]
    pub fn new(g_r: f64, g_i: f64, normal: (f64, f64, f64), area: f64) -> Self {
        Self {
            inner: SpinInterface {
                g_r,
                g_i,
                normal: Vector3::new(normal.0, normal.1, normal.2),
                area,
            },
        }
    }

    /// Real part of spin mixing conductance (1/(Ohm*m^2))
    #[getter]
    pub fn g_r(&self) -> f64 {
        self.inner.g_r
    }

    /// Imaginary part of spin mixing conductance (1/(Ohm*m^2))
    #[getter]
    pub fn g_i(&self) -> f64 {
        self.inner.g_i
    }

    /// Interface normal direction
    #[getter]
    pub fn normal(&self) -> PyVector3 {
        PyVector3::from_inner(self.inner.normal)
    }

    /// Interface area (m^2)
    #[getter]
    pub fn area(&self) -> f64 {
        self.inner.area
    }

    /// Create YIG/Pt interface
    ///
    /// The YIG/Pt interface is the canonical system for
    /// spin pumping and ISHE detection experiments.
    #[staticmethod]
    pub fn yig_pt() -> Self {
        Self {
            inner: SpinInterface::yig_pt(),
        }
    }

    /// Create Permalloy/Pt interface
    #[staticmethod]
    pub fn py_pt() -> Self {
        Self {
            inner: SpinInterface::py_pt(),
        }
    }

    /// String representation
    pub fn __repr__(&self) -> String {
        format!(
            "SpinInterface(g_r={:.2e}, g_i={:.2e}, area={:.2e})",
            self.inner.g_r, self.inner.g_i, self.inner.area
        )
    }
}

impl PySpinInterface {
    /// Get the inner SpinInterface
    pub fn inner(&self) -> &SpinInterface {
        &self.inner
    }
}
