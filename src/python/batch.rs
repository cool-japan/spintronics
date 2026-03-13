//! Python/numpy bindings for SIMD batch LLG evolution.
//!
//! Exposes the Structure-of-Arrays (SoA) SIMD batch processing functions to
//! Python via numpy arrays. The (N,3) numpy array layout (AoS) is converted
//! to/from the internal SoA `SimdBatch` layout transparently.
//!
//! ## Python Usage
//!
//! ```python
//! import numpy as np
//! import spintronics
//!
//! N = 1024
//! # N spins along x, field along z
//! m = np.zeros((N, 3))
//! m[:, 0] = 1.0
//!
//! h = np.zeros((N, 3))
//! h[:, 2] = 0.1   # 100 mT
//!
//! alpha, gamma, dt = 0.01, 1.76085963e11, 1e-13
//!
//! # One RK4 step
//! m1 = spintronics.batch_rk4_step(m, h, alpha=alpha, gamma=gamma, dt=dt)
//!
//! # 100 RK4 steps at once (more efficient than looping)
//! m100 = spintronics.batch_rk4_multistep(m, h, alpha=alpha, gamma=gamma,
//!                                          dt=dt, num_steps=100)
//! ```

use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::simd::{batch_evolve_multi_step, batch_evolve_rk4, SimdBatch};
use crate::vector3::Vector3;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers: (N,3) numpy array ↔ SimdBatch SoA layout
// ─────────────────────────────────────────────────────────────────────────────

/// Convert a readonly (N, 3) numpy f64 array into a `SimdBatch` (SoA layout).
///
/// Validates that the array has exactly 2 dimensions and 3 columns.
fn numpy_to_simd(arr: &PyReadonlyArray2<f64>) -> PyResult<SimdBatch> {
    let shape = arr.shape();
    if shape.len() != 2 {
        return Err(PyValueError::new_err(format!(
            "expected 2D array, got {}D",
            shape.len()
        )));
    }
    if shape[1] != 3 {
        return Err(PyValueError::new_err(format!(
            "expected shape (N, 3), got ({}, {})",
            shape[0], shape[1]
        )));
    }

    let n = shape[0];
    let slice = arr
        .as_slice()
        .map_err(|e| PyValueError::new_err(format!("array is not contiguous: {e}")))?;

    let spins: Vec<Vector3<f64>> = (0..n)
        .map(|i| Vector3::new(slice[i * 3], slice[i * 3 + 1], slice[i * 3 + 2]))
        .collect();

    Ok(SimdBatch::from_vector3_slice(&spins))
}

/// Convert a `SimdBatch` back to a numpy (N, 3) f64 array.
fn simd_to_numpy<'py>(py: Python<'py>, batch: &SimdBatch) -> Bound<'py, PyArray2<f64>> {
    let n = batch.len();
    let mut flat = Vec::with_capacity(n * 3);
    for i in 0..n {
        flat.push(batch.x[i]);
        flat.push(batch.y[i]);
        flat.push(batch.z[i]);
    }
    // ndarray-style: reshape flat Vec into (N,3)
    let arr = numpy::ndarray::Array2::from_shape_vec((n, 3), flat).expect("shape is always (n, 3)");
    arr.into_pyarray(py)
}

// ─────────────────────────────────────────────────────────────────────────────
// Public Python functions
// ─────────────────────────────────────────────────────────────────────────────

/// Evolve N spins for one RK4 time step using SIMD batch processing.
///
/// Each spin evolves independently under its own effective field (no coupling).
/// The result is normalized so that |m_i| = 1 for every spin.
///
/// Args:
///     m: (N, 3) numpy array of magnetization vectors (row-major, each row is m_i)
///     h_eff: (N, 3) numpy array of effective field vectors [T]
///     alpha: Gilbert damping constant (dimensionless)
///     gamma: Gyromagnetic ratio [rad/(s·T)], default ≈ 1.761e11
///     dt: Integration time step [s]
///
/// Returns:
///     (N, 3) numpy array of updated magnetization vectors (unit length)
///
/// Raises:
///     ValueError: If array shapes are incompatible or lengths differ
#[pyfunction]
pub fn batch_rk4_step<'py>(
    py: Python<'py>,
    m: PyReadonlyArray2<'py, f64>,
    h_eff: PyReadonlyArray2<'py, f64>,
    alpha: f64,
    gamma: f64,
    dt: f64,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let m_batch = numpy_to_simd(&m)?;
    let h_batch = numpy_to_simd(&h_eff)?;

    if m_batch.len() != h_batch.len() {
        return Err(PyValueError::new_err(format!(
            "m and h_eff must have the same number of rows: {} vs {}",
            m_batch.len(),
            h_batch.len()
        )));
    }

    let m_new = batch_evolve_rk4(&m_batch, &h_batch, alpha, gamma, dt)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    Ok(simd_to_numpy(py, &m_new))
}

/// Evolve N spins for `num_steps` RK4 time steps using SIMD batch processing.
///
/// `h_eff` is held constant throughout (static field approximation).
/// This is more efficient than calling `batch_rk4_step` in a loop because
/// no Python-to-Rust conversion overhead is incurred at each step.
///
/// The result is normalized so that |m_i| = 1 after all steps.
///
/// Args:
///     m: (N, 3) numpy array of magnetization vectors (initial state)
///     h_eff: (N, 3) numpy array of effective field vectors [T] (constant)
///     alpha: Gilbert damping constant (dimensionless)
///     gamma: Gyromagnetic ratio [rad/(s·T)]
///     dt: Integration time step [s]
///     num_steps: Number of RK4 iterations to perform
///
/// Returns:
///     (N, 3) numpy array of final magnetization vectors (unit length)
///
/// Raises:
///     ValueError: If array shapes are incompatible or lengths differ
#[pyfunction]
pub fn batch_rk4_multistep<'py>(
    py: Python<'py>,
    m: PyReadonlyArray2<'py, f64>,
    h_eff: PyReadonlyArray2<'py, f64>,
    alpha: f64,
    gamma: f64,
    dt: f64,
    num_steps: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let m_batch = numpy_to_simd(&m)?;
    let h_batch = numpy_to_simd(&h_eff)?;

    if m_batch.len() != h_batch.len() {
        return Err(PyValueError::new_err(format!(
            "m and h_eff must have the same number of rows: {} vs {}",
            m_batch.len(),
            h_batch.len()
        )));
    }

    let m_new = batch_evolve_multi_step(m_batch, &h_batch, alpha, gamma, dt, num_steps)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    Ok(simd_to_numpy(py, &m_new))
}
