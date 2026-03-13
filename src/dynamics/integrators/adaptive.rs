//! Adaptive step-size controller wrapping any [`Integrator`] with error estimates.

use super::rhs_fn::{Integrator, IntegratorOutput, RhsFn};
use crate::error::{Error, Result};
use crate::vector3::Vector3;

// =========================================================================
// Adaptive step size controller
// =========================================================================

/// Adaptive step size controller wrapping any integrator that provides
/// error estimates.
///
/// Uses a PI controller (proportional-integral) for smooth step size
/// adjustment, avoiding the oscillations of a pure proportional controller.
///
/// # Algorithm
///
/// The step size is updated according to:
///
/// ```text
/// dt_new = dt * S * (tol / err_n)^alpha * (err_{n-1} / tol)^beta
/// ```
///
/// where `S` is a safety factor (0.9), `alpha` and `beta` are the PI gains
/// tuned for the method order, and `err_n`, `err_{n-1}` are the current and
/// previous error estimates.
pub struct AdaptiveIntegrator<I: Integrator> {
    /// The underlying fixed-step integrator.
    pub integrator: I,
    /// Desired local error tolerance.
    pub tolerance: f64,
    /// Safety factor (default 0.9).
    pub safety: f64,
    /// Minimum allowed step size.
    pub dt_min: f64,
    /// Maximum allowed step size.
    pub dt_max: f64,
    /// Maximum step size growth factor per step.
    pub max_factor: f64,
    /// Minimum step size shrink factor per step.
    pub min_factor: f64,
    /// Order of the error estimator (for PI controller gains).
    pub order: f64,
    /// Previous error (for PI controller).
    prev_error: f64,
}

impl<I: Integrator> AdaptiveIntegrator<I> {
    /// Create a new adaptive integrator.
    ///
    /// # Arguments
    /// * `integrator` - The underlying fixed-step integrator.
    /// * `tolerance`  - Desired local truncation error bound.
    /// * `order`      - Order of the lower-order embedded method (e.g., 4 for DP45).
    pub fn new(integrator: I, tolerance: f64, order: f64) -> Self {
        Self {
            integrator,
            tolerance,
            safety: 0.9,
            dt_min: 1e-20,
            dt_max: 1e10,
            max_factor: 5.0,
            min_factor: 0.2,
            order,
            prev_error: tolerance,
        }
    }

    /// Set the minimum allowed step size.
    pub fn with_dt_min(mut self, dt_min: f64) -> Self {
        self.dt_min = dt_min;
        self
    }

    /// Set the maximum allowed step size.
    pub fn with_dt_max(mut self, dt_max: f64) -> Self {
        self.dt_max = dt_max;
        self
    }

    /// Set the safety factor.
    pub fn with_safety(mut self, safety: f64) -> Self {
        self.safety = safety;
        self
    }

    /// Perform an adaptive step, rejecting and retrying if the error is too large.
    ///
    /// Returns the accepted output together with the actual step size used.
    ///
    /// # Arguments
    /// * `state` - Current state vector.
    /// * `t`     - Current time.
    /// * `dt`    - Initial step size guess.
    /// * `f`     - Right-hand-side function.
    ///
    /// # Errors
    /// Returns an error if the step size shrinks below `dt_min` without
    /// achieving the tolerance, or if the integrator itself fails.
    pub fn adaptive_step(
        &mut self,
        state: &[Vector3<f64>],
        t: f64,
        dt: f64,
        f: &RhsFn<'_>,
    ) -> Result<(IntegratorOutput, f64)> {
        let mut current_dt = dt.clamp(self.dt_min, self.dt_max);

        loop {
            let output = self.integrator.step(state, t, current_dt, f)?;

            let error = match output.error_estimate {
                Some(e) if e > 0.0 => e,
                Some(_) => {
                    // Error is zero or negative — accept and grow
                    let new_dt = (current_dt * self.max_factor).clamp(self.dt_min, self.dt_max);
                    self.prev_error = self.tolerance;
                    let mut accepted = output;
                    accepted.suggested_dt = Some(new_dt);
                    return Ok((accepted, current_dt));
                },
                None => {
                    // Integrator provides no error estimate — accept as-is
                    return Ok((output, current_dt));
                },
            };

            // PI controller gains
            let alpha_pi = 0.7 / (self.order + 1.0);
            let beta_pi = 0.4 / (self.order + 1.0);

            let factor = self.safety
                * (self.tolerance / error).powf(alpha_pi)
                * (self.prev_error / self.tolerance).powf(beta_pi);

            let factor = factor.clamp(self.min_factor, self.max_factor);

            if error <= self.tolerance {
                // Accept step
                self.prev_error = error;
                let new_dt = (current_dt * factor).clamp(self.dt_min, self.dt_max);
                let mut accepted = output;
                accepted.suggested_dt = Some(new_dt);
                return Ok((accepted, current_dt));
            }

            // Reject step — shrink dt and retry
            current_dt = (current_dt * factor).clamp(self.dt_min, self.dt_max);

            if current_dt <= self.dt_min {
                return Err(Error::NumericalError {
                    description: format!(
                        "Adaptive integrator: step size {:.2e} reached minimum {:.2e} \
                         without achieving tolerance {:.2e} (error = {:.2e})",
                        current_dt, self.dt_min, self.tolerance, error
                    ),
                });
            }
        }
    }
}
