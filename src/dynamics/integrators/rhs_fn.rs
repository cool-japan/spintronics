//! Core data types, the [`Integrator`] trait, and vector arithmetic helpers.

use crate::error::{Error, Result};
use crate::vector3::Vector3;

// ---------------------------------------------------------------------------
// Core data types
// ---------------------------------------------------------------------------

/// Output produced by a single integrator step.
#[derive(Debug, Clone)]
pub struct IntegratorOutput {
    /// The state vector after the step.
    pub new_state: Vec<Vector3<f64>>,
    /// Estimated local truncation error (if the method supports it).
    pub error_estimate: Option<f64>,
    /// Suggested next time step (if the method supports adaptive control).
    pub suggested_dt: Option<f64>,
}

/// The right-hand-side function type used by integrators.
///
/// Given a state slice and current time, returns the derivative for each
/// component.
pub type RhsFn<'a> = dyn Fn(&[Vector3<f64>], f64) -> Vec<Vector3<f64>> + 'a;

// ---------------------------------------------------------------------------
// Integrator trait
// ---------------------------------------------------------------------------

/// Trait for numerical integrators that advance a system of `Vector3` ODEs.
pub trait Integrator {
    /// Advance the state by one step.
    ///
    /// # Arguments
    /// * `state` - Current state vector
    /// * `t`     - Current time
    /// * `dt`    - Time step size
    /// * `f`     - Right-hand-side function  dy/dt = f(y, t)
    ///
    /// # Errors
    /// Returns an error if the integration step fails (e.g. NaN detected).
    fn step(
        &mut self,
        state: &[Vector3<f64>],
        t: f64,
        dt: f64,
        f: &RhsFn<'_>,
    ) -> Result<IntegratorOutput>;
}

// ---------------------------------------------------------------------------
// Helper: vector arithmetic on slices
// ---------------------------------------------------------------------------

/// y = a + b * scalar
#[inline]
pub(super) fn vec_add_scaled(a: &[Vector3<f64>], b: &[Vector3<f64>], s: f64) -> Vec<Vector3<f64>> {
    a.iter()
        .zip(b.iter())
        .map(|(&ai, &bi)| ai + bi * s)
        .collect()
}

/// sum of multiple k-vectors with given coefficients: base + sum(c_i * k_i)
#[inline]
pub(super) fn vec_combine(
    base: &[Vector3<f64>],
    ks: &[&[Vector3<f64>]],
    coeffs: &[f64],
) -> Vec<Vector3<f64>> {
    let n = base.len();
    let mut out = base.to_vec();
    for i in 0..n {
        for (k, &c) in ks.iter().zip(coeffs.iter()) {
            out[i] = out[i] + k[i] * c;
        }
    }
    out
}

/// Compute the maximum norm of the difference between two state vectors.
#[inline]
pub(super) fn max_error_norm(a: &[Vector3<f64>], b: &[Vector3<f64>]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(&ai, &bi)| (ai - bi).magnitude())
        .fold(0.0_f64, f64::max)
}

/// Check a state for NaN and return an error if found.
pub(super) fn check_nan(state: &[Vector3<f64>]) -> Result<()> {
    for v in state {
        if v.x.is_nan() || v.y.is_nan() || v.z.is_nan() {
            return Err(Error::NumericalError {
                description: "NaN detected in integrator output".to_string(),
            });
        }
    }
    Ok(())
}
