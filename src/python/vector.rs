//! Python bindings for Vector3

use pyo3::prelude::*;

use crate::vector3::Vector3;

/// A 3D vector for spintronics calculations
///
/// This class wraps the Rust `Vector3<f64>` type for use in Python.
///
/// ## Example
/// ```python
/// v1 = Vector3(1.0, 0.0, 0.0)
/// v2 = Vector3(0.0, 1.0, 0.0)
/// cross = v1.cross(v2)
/// print(f"Cross product: ({cross.x}, {cross.y}, {cross.z})")
/// ```
#[pyclass(name = "Vector3")]
#[derive(Clone)]
pub struct PyVector3 {
    inner: Vector3<f64>,
}

#[pymethods]
impl PyVector3 {
    /// Create a new 3D vector
    #[new]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            inner: Vector3::new(x, y, z),
        }
    }

    /// X component
    #[getter]
    pub fn x(&self) -> f64 {
        self.inner.x
    }

    /// Y component
    #[getter]
    pub fn y(&self) -> f64 {
        self.inner.y
    }

    /// Z component
    #[getter]
    pub fn z(&self) -> f64 {
        self.inner.z
    }

    /// Set X component
    #[setter]
    pub fn set_x(&mut self, value: f64) {
        self.inner.x = value;
    }

    /// Set Y component
    #[setter]
    pub fn set_y(&mut self, value: f64) {
        self.inner.y = value;
    }

    /// Set Z component
    #[setter]
    pub fn set_z(&mut self, value: f64) {
        self.inner.z = value;
    }

    /// Calculate the dot product with another vector
    pub fn dot(&self, other: &PyVector3) -> f64 {
        self.inner.dot(&other.inner)
    }

    /// Calculate the cross product with another vector
    pub fn cross(&self, other: &PyVector3) -> PyVector3 {
        PyVector3 {
            inner: self.inner.cross(&other.inner),
        }
    }

    /// Calculate the magnitude (length) of the vector
    pub fn magnitude(&self) -> f64 {
        self.inner.magnitude()
    }

    /// Return a normalized (unit) vector in the same direction
    pub fn normalize(&self) -> PyVector3 {
        PyVector3 {
            inner: self.inner.normalize(),
        }
    }

    /// Add two vectors
    pub fn __add__(&self, other: &PyVector3) -> PyVector3 {
        PyVector3 {
            inner: self.inner + other.inner,
        }
    }

    /// Subtract two vectors
    pub fn __sub__(&self, other: &PyVector3) -> PyVector3 {
        PyVector3 {
            inner: self.inner - other.inner,
        }
    }

    /// Multiply vector by scalar
    pub fn __mul__(&self, scalar: f64) -> PyVector3 {
        PyVector3 {
            inner: self.inner * scalar,
        }
    }

    /// Right multiply (scalar * vector)
    pub fn __rmul__(&self, scalar: f64) -> PyVector3 {
        PyVector3 {
            inner: self.inner * scalar,
        }
    }

    /// String representation
    pub fn __repr__(&self) -> String {
        format!(
            "Vector3({:.6}, {:.6}, {:.6})",
            self.inner.x, self.inner.y, self.inner.z
        )
    }

    /// Convert to tuple
    pub fn to_tuple(&self) -> (f64, f64, f64) {
        (self.inner.x, self.inner.y, self.inner.z)
    }

    /// Convert to list
    pub fn to_list(&self) -> Vec<f64> {
        vec![self.inner.x, self.inner.y, self.inner.z]
    }

    /// Create unit vector in x direction
    #[staticmethod]
    pub fn unit_x() -> PyVector3 {
        PyVector3 {
            inner: Vector3::new(1.0, 0.0, 0.0),
        }
    }

    /// Create unit vector in y direction
    #[staticmethod]
    pub fn unit_y() -> PyVector3 {
        PyVector3 {
            inner: Vector3::new(0.0, 1.0, 0.0),
        }
    }

    /// Create unit vector in z direction
    #[staticmethod]
    pub fn unit_z() -> PyVector3 {
        PyVector3 {
            inner: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Create zero vector
    #[staticmethod]
    pub fn zero() -> PyVector3 {
        PyVector3 {
            inner: Vector3::new(0.0, 0.0, 0.0),
        }
    }
}

impl PyVector3 {
    /// Get the inner Vector3
    pub fn inner(&self) -> Vector3<f64> {
        self.inner
    }

    /// Create from inner Vector3
    pub fn from_inner(inner: Vector3<f64>) -> Self {
        Self { inner }
    }
}
