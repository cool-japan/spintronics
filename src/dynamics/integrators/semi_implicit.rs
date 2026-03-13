//! Semi-implicit integrator for stiff ODE systems.

use super::rhs_fn::{check_nan, max_error_norm, Integrator, IntegratorOutput, RhsFn};
use crate::error::Result;
use crate::vector3::Vector3;

// =========================================================================
// 4. Semi-implicit integrator for stiff problems
// =========================================================================

/// Semi-implicit integrator using the implicit midpoint rule.
///
/// Solves stiff ODE systems (large exchange coupling, small damping) via
/// fixed-point iteration on the implicit midpoint equation:
///
///   y_{n+1} = y_n + dt * f( (y_n + y_{n+1})/2, t + dt/2 )
///
/// The implicit midpoint rule is A-stable (unconditionally stable for linear
/// problems) and preserves quadratic invariants, making it suitable for
/// near-Hamiltonian systems with mild dissipation.
pub struct SemiImplicit {
    /// Maximum number of fixed-point iterations.
    pub max_iterations: usize,
    /// Convergence tolerance for the fixed-point iteration.
    pub tolerance: f64,
}

impl SemiImplicit {
    /// Create a new semi-implicit integrator.
    ///
    /// # Arguments
    /// * `max_iterations` - Maximum fixed-point iterations per step.
    /// * `tolerance`      - Convergence threshold (max-norm of successive iterates).
    pub fn new(max_iterations: usize, tolerance: f64) -> Self {
        Self {
            max_iterations,
            tolerance,
        }
    }
}

impl Default for SemiImplicit {
    fn default() -> Self {
        Self::new(50, 1e-12)
    }
}

impl Integrator for SemiImplicit {
    fn step(
        &mut self,
        state: &[Vector3<f64>],
        t: f64,
        dt: f64,
        f: &RhsFn<'_>,
    ) -> Result<IntegratorOutput> {
        let n = state.len();
        let t_mid = t + dt * 0.5;

        // Initial guess: explicit Euler
        let f0 = f(state, t);
        let mut y_new: Vec<Vector3<f64>> = state
            .iter()
            .zip(f0.iter())
            .map(|(&si, &fi)| si + fi * dt)
            .collect();

        let mut converged = false;
        for _ in 0..self.max_iterations {
            // Midpoint
            let y_mid: Vec<Vector3<f64>> = state
                .iter()
                .zip(y_new.iter())
                .map(|(&si, &yi)| (si + yi) * 0.5)
                .collect();

            let f_mid = f(&y_mid, t_mid);

            let y_next: Vec<Vector3<f64>> = state
                .iter()
                .zip(f_mid.iter())
                .map(|(&si, &fi)| si + fi * dt)
                .collect();

            let diff = max_error_norm(&y_next, &y_new);
            y_new = y_next;

            if diff < self.tolerance {
                converged = true;
                break;
            }
        }

        check_nan(&y_new)?;

        if !converged {
            // Still return the best estimate but report the issue through a
            // large error estimate so the adaptive wrapper can shrink the step.
            let residual = {
                let y_mid: Vec<Vector3<f64>> = state
                    .iter()
                    .zip(y_new.iter())
                    .map(|(&si, &yi)| (si + yi) * 0.5)
                    .collect();
                let f_mid = f(&y_mid, t_mid);
                let y_check: Vec<Vector3<f64>> = state
                    .iter()
                    .zip(f_mid.iter())
                    .map(|(&si, &fi)| si + fi * dt)
                    .collect();
                max_error_norm(&y_check, &y_new)
            };

            return Ok(IntegratorOutput {
                new_state: y_new,
                error_estimate: Some(residual),
                suggested_dt: Some(dt * 0.5),
            });
        }

        // Estimate error via the difference between implicit midpoint and
        // explicit Euler (a cheap 1st-order method).
        let y_euler: Vec<Vector3<f64>> = (0..n).map(|i| state[i] + f0[i] * dt).collect();
        let error = max_error_norm(&y_new, &y_euler);

        Ok(IntegratorOutput {
            new_state: y_new,
            error_estimate: Some(error),
            suggested_dt: None,
        })
    }
}
