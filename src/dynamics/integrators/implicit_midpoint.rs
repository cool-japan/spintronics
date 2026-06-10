//! Implicit midpoint integrator with Newton-Raphson + finite-difference Jacobian.
//!
//! Solves
//!
//!   y_{n+1} = y_n + dt · f((y_n + y_{n+1})/2, t_n + dt/2)
//!
//! by treating the implicit equation
//!
//!   F(y) = y - y_n - dt · f((y_n + y)/2, t_n + dt/2) = 0
//!
//! as a root-finding problem for y = y_{n+1}, and solving it with full
//! Newton-Raphson iteration. The Jacobian J = ∂F/∂y is constructed numerically
//! via central finite differences on the right-hand side, evaluated at the
//! current midpoint guess:
//!
//!   J(y) = I - (dt/2) · ∂f/∂y_mid |_{y_mid = (y_n + y)/2}
//!
//! For stiff systems (e.g. the Landau-Lifshitz-Gilbert equation with strong
//! damping α near 1 or with strong exchange/anisotropy fields), explicit
//! methods such as RK4 require dt of order 1/(γ H_eff), while the implicit
//! midpoint rule is A-stable and remains stable for arbitrarily large dt.
//! The trade-off is the cost of Newton iteration (typically 2–5 iterations
//! per step) and finite-difference Jacobian construction
//! (2 · 3N right-hand-side evaluations for N spin DOFs using central
//! differences).
//!
//! The implicit midpoint rule is second-order accurate, symmetric in time,
//! symplectic for canonical Hamiltonian systems, and exactly preserves
//! quadratic invariants — in particular |m|² for the LLG equation in its
//! Landau-Lifshitz form. This makes it a popular choice for long-time
//! integration of micromagnetic systems.
//!
//! # References
//! - Hairer, E. and Wanner, G., *Solving Ordinary Differential Equations II*,
//!   Springer, Chapter IV (1996).
//! - Mentink, J. H. et al., *Stable and fast semi-implicit integration of the
//!   stochastic Landau-Lifshitz equation*, J. Phys.: Condens. Matter **22**,
//!   176001 (2010).

use super::rhs_fn::{check_nan, Integrator, IntegratorOutput, RhsFn};
use crate::error::{numerical_error, Result};
use crate::vector3::Vector3;

// =========================================================================
// ImplicitMidpointNewton
// =========================================================================

/// Implicit midpoint integrator using Newton-Raphson with a numerical
/// (finite-difference) Jacobian.
///
/// Unlike the [`super::SemiImplicit`] fixed-point variant, this integrator
/// can handle problems where the fixed-point iteration fails to converge
/// (typically when |λ dt/2| ≥ 1 for a stiff eigenmode). It builds and
/// factorises the full Jacobian at each Newton step, giving quadratic
/// convergence near the root.
pub struct ImplicitMidpointNewton {
    /// Maximum number of Newton iterations per step.
    pub max_newton_iter: usize,
    /// Convergence tolerance on the Newton residual (L²-norm of F(y)).
    pub newton_tol: f64,
    /// Finite-difference step h used to estimate ∂f/∂y. Central differences
    /// give an error of O(h²) on the Jacobian, so values around the square
    /// root of machine epsilon (≈ 1.5e-8) are a good default.
    pub fd_step: f64,
    /// Counter for how many `step` calls have been issued (for diagnostics).
    pub max_step_count: usize,
}

impl ImplicitMidpointNewton {
    /// Construct an implicit midpoint integrator with the default settings:
    /// `max_newton_iter = 20`, `newton_tol = 1e-12`, `fd_step = 1e-8`.
    pub fn new() -> Self {
        Self {
            max_newton_iter: 20,
            newton_tol: 1e-12,
            fd_step: 1e-8,
            max_step_count: 0,
        }
    }

    /// Set the maximum number of Newton iterations.
    pub fn with_max_iter(mut self, n: usize) -> Self {
        self.max_newton_iter = n;
        self
    }

    /// Set the Newton convergence tolerance (L²-norm of F).
    pub fn with_tol(mut self, tol: f64) -> Self {
        self.newton_tol = tol;
        self
    }

    /// Set the finite-difference step used for the Jacobian.
    pub fn with_fd_step(mut self, h: f64) -> Self {
        self.fd_step = h;
        self
    }
}

impl Default for ImplicitMidpointNewton {
    fn default() -> Self {
        Self::new()
    }
}

impl Integrator for ImplicitMidpointNewton {
    fn step(
        &mut self,
        state: &[Vector3<f64>],
        t: f64,
        dt: f64,
        f: &RhsFn<'_>,
    ) -> Result<IntegratorOutput> {
        let n = state.len();
        let dof = 3 * n;
        let t_mid = t + 0.5 * dt;
        let half_dt = 0.5 * dt;

        if self.fd_step.abs() < 1e-15 {
            return Err(numerical_error(
                "ImplicitMidpointNewton: fd_step is too small (must be >= 1e-15)",
            ));
        }
        if !dt.is_finite() || dt == 0.0 {
            return Err(numerical_error(
                "ImplicitMidpointNewton: dt must be finite and non-zero",
            ));
        }

        // Initial Newton guess via explicit Euler.
        let f0 = f(state, t);
        let mut y_new: Vec<Vector3<f64>> = state
            .iter()
            .zip(f0.iter())
            .map(|(&si, &fi)| si + fi * dt)
            .collect();

        let state_flat = vector3_to_flat(state);
        let mut residual_norm = f64::INFINITY;
        let mut converged = false;

        for _ in 0..self.max_newton_iter {
            // Midpoint y_mid = (y_n + y_k) / 2
            let y_mid: Vec<Vector3<f64>> = state
                .iter()
                .zip(y_new.iter())
                .map(|(&si, &yi)| (si + yi) * 0.5)
                .collect();

            let f_mid = f(&y_mid, t_mid);

            // F(y_k) = y_k - y_n - dt * f(y_mid, t_mid)
            let mut residual = vec![0.0_f64; dof];
            let y_flat = vector3_to_flat(&y_new);
            let f_flat = vector3_to_flat(&f_mid);
            for i in 0..dof {
                residual[i] = y_flat[i] - state_flat[i] - dt * f_flat[i];
            }

            residual_norm = l2_norm(&residual);
            if residual_norm < self.newton_tol {
                converged = true;
                break;
            }

            // Build Jacobian J = I - (dt/2) * df/dy_mid via central differences.
            let mut jacobian = vec![0.0_f64; dof * dof];
            build_jacobian(f, &y_mid, t_mid, self.fd_step, dof, half_dt, &mut jacobian);

            // Solve J · delta = -F  for the Newton step "delta".
            let mut rhs = vec![0.0_f64; dof];
            for i in 0..dof {
                rhs[i] = -residual[i];
            }

            let delta = match gauss_solve(jacobian, rhs, dof) {
                Some(d) => d,
                None => {
                    return Err(numerical_error(
                        "ImplicitMidpointNewton: singular Jacobian during Newton solve",
                    ));
                },
            };

            // y_{k+1} = y_k + delta
            for i in 0..dof {
                let new_val = y_flat[i] + delta[i];
                if !new_val.is_finite() {
                    return Err(numerical_error(
                        "ImplicitMidpointNewton: non-finite Newton update",
                    ));
                }
            }
            let mut updated_flat = vec![0.0_f64; dof];
            for i in 0..dof {
                updated_flat[i] = y_flat[i] + delta[i];
            }
            y_new = flat_to_vector3(&updated_flat);
        }

        check_nan(&y_new)?;

        self.max_step_count = self.max_step_count.saturating_add(1);

        if !converged {
            return Err(numerical_error(&format!(
                "ImplicitMidpointNewton: Newton failed to converge in {} iterations \
                 (residual = {:.3e}, tol = {:.3e})",
                self.max_newton_iter, residual_norm, self.newton_tol,
            )));
        }

        Ok(IntegratorOutput {
            new_state: y_new,
            error_estimate: Some(residual_norm),
            suggested_dt: None,
        })
    }
}

// =========================================================================
// Helpers
// =========================================================================

/// Flatten a `Vec<Vector3<f64>>` of length N into a `Vec<f64>` of length 3N
/// using the layout `[x0, y0, z0, x1, y1, z1, ...]`.
fn vector3_to_flat(v: &[Vector3<f64>]) -> Vec<f64> {
    let mut out = Vec::with_capacity(3 * v.len());
    for vi in v {
        out.push(vi.x);
        out.push(vi.y);
        out.push(vi.z);
    }
    out
}

/// Inverse of [`vector3_to_flat`].
fn flat_to_vector3(flat: &[f64]) -> Vec<Vector3<f64>> {
    debug_assert!(flat.len() % 3 == 0);
    let n = flat.len() / 3;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let off = 3 * i;
        out.push(Vector3::new(flat[off], flat[off + 1], flat[off + 2]));
    }
    out
}

/// Euclidean (L²) norm of a flat slice.
fn l2_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Build the Newton Jacobian J = I - half_dt · ∂f/∂y_mid at the supplied
/// midpoint state. The result is written into `jacobian` in row-major order
/// (J[row, col] = jacobian[row * dof + col]).
fn build_jacobian(
    f: &RhsFn<'_>,
    y_mid: &[Vector3<f64>],
    t_mid: f64,
    fd_step: f64,
    dof: usize,
    half_dt: f64,
    jacobian: &mut [f64],
) {
    // Adaptive finite-difference step: h_j = fd_step * max(|y_j|, 1).
    let y_mid_flat = vector3_to_flat(y_mid);

    for j in 0..dof {
        let h = fd_step * y_mid_flat[j].abs().max(1.0);
        let two_h = 2.0 * h;

        // y_mid + h e_j
        let mut yp = y_mid_flat.clone();
        yp[j] += h;
        let f_plus = f(&flat_to_vector3(&yp), t_mid);
        let f_plus_flat = vector3_to_flat(&f_plus);

        // y_mid - h e_j
        let mut ym = y_mid_flat.clone();
        ym[j] -= h;
        let f_minus = f(&flat_to_vector3(&ym), t_mid);
        let f_minus_flat = vector3_to_flat(&f_minus);

        // Column j of df/dy: (f_plus - f_minus) / (2h)
        // Then J[:, j] = -half_dt · column   (the identity contribution is
        // added below).
        for row in 0..dof {
            let dfi = (f_plus_flat[row] - f_minus_flat[row]) / two_h;
            jacobian[row * dof + j] = -half_dt * dfi;
        }
    }

    // Add the identity matrix.
    for i in 0..dof {
        jacobian[i * dof + i] += 1.0;
    }
}

/// Solve A x = b via Gauss elimination with partial pivoting.
///
/// `a` is the row-major n×n matrix, `b` is the right-hand side. Returns
/// `Some(x)` on success or `None` if the matrix is numerically singular.
fn gauss_solve(mut a: Vec<f64>, mut b: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    if n == 0 {
        return Some(Vec::new());
    }
    debug_assert_eq!(a.len(), n * n);
    debug_assert_eq!(b.len(), n);

    // Forward elimination with partial pivoting.
    for k in 0..n {
        // Find pivot row.
        let mut pivot_row = k;
        let mut pivot_val = a[k * n + k].abs();
        for r in (k + 1)..n {
            let v = a[r * n + k].abs();
            if v > pivot_val {
                pivot_val = v;
                pivot_row = r;
            }
        }
        if pivot_val < 1e-30 {
            return None;
        }

        if pivot_row != k {
            // Swap rows k and pivot_row in A.
            for c in 0..n {
                a.swap(k * n + c, pivot_row * n + c);
            }
            b.swap(k, pivot_row);
        }

        // Eliminate below the pivot.
        let inv_pivot = 1.0 / a[k * n + k];
        for r in (k + 1)..n {
            let factor = a[r * n + k] * inv_pivot;
            if factor == 0.0 {
                continue;
            }
            for c in k..n {
                a[r * n + c] -= factor * a[k * n + c];
            }
            b[r] -= factor * b[k];
        }
    }

    // Back substitution.
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = b[i];
        for c in (i + 1)..n {
            sum -= a[i * n + c] * x[c];
        }
        let diag = a[i * n + i];
        if diag.abs() < 1e-30 {
            return None;
        }
        x[i] = sum / diag;
    }
    Some(x)
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GAMMA;

    // -- Helper RHS functions --------------------------------------------

    /// dy/dt = -y  (linear decay)
    fn exponential_rhs(state: &[Vector3<f64>], _t: f64) -> Vec<Vector3<f64>> {
        state.iter().map(|&v| v * (-1.0)).collect()
    }

    /// dy/dt = -lambda y  (Dahlquist test)
    fn make_dahlquist_rhs(lambda: f64) -> impl Fn(&[Vector3<f64>], f64) -> Vec<Vector3<f64>> {
        move |state: &[Vector3<f64>], _t: f64| state.iter().map(|&v| v * (-lambda)).collect()
    }

    /// dy/dt = (-y.y, y.x, 0) - alpha y  (damped precession)
    fn damped_precession_rhs(state: &[Vector3<f64>], _t: f64) -> Vec<Vector3<f64>> {
        let alpha = 0.1;
        state
            .iter()
            .map(|v| Vector3::new(-v.y, v.x, 0.0) + *v * (-alpha))
            .collect()
    }

    /// Landau-Lifshitz form for a single spin in a constant field B along z:
    ///   dm/dt = -gamma m × B - alpha gamma m × (m × B)
    fn make_llg_rhs(
        b_field: Vector3<f64>,
        alpha: f64,
    ) -> impl Fn(&[Vector3<f64>], f64) -> Vec<Vector3<f64>> {
        move |state: &[Vector3<f64>], _t: f64| {
            state
                .iter()
                .map(|m| {
                    let prec = m.cross(&b_field) * (-GAMMA);
                    let damp = m.cross(&m.cross(&b_field)) * (-alpha * GAMMA);
                    prec + damp
                })
                .collect()
        }
    }

    // -- Test 1: constructor / builder methods ---------------------------

    #[test]
    fn test_builder_methods() {
        let integ = ImplicitMidpointNewton::new()
            .with_max_iter(50)
            .with_tol(1e-14)
            .with_fd_step(1e-7);
        assert_eq!(integ.max_newton_iter, 50);
        assert!((integ.newton_tol - 1e-14).abs() < 1e-30);
        assert!((integ.fd_step - 1e-7).abs() < 1e-20);
    }

    // -- Test 2: default values ------------------------------------------

    #[test]
    fn test_default_values() {
        let integ = ImplicitMidpointNewton::default();
        assert_eq!(integ.max_newton_iter, 20);
        assert!((integ.newton_tol - 1e-12).abs() < 1e-30);
        assert!((integ.fd_step - 1e-8).abs() < 1e-20);
        assert_eq!(integ.max_step_count, 0);
    }

    // -- Test 3: A-stability on Dahlquist test ---------------------------

    #[test]
    fn test_a_stability_dahlquist() {
        // For lambda = -100, an explicit RK4 needs dt < 2.78/100 ~ 0.028.
        // The implicit midpoint rule is A-stable, so dt = 1.0 must still
        // produce a bounded, monotonically decaying result.
        let lambda = 100.0;
        let rhs = make_dahlquist_rhs(lambda);
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let dt = 1.0;
        let mut integ = ImplicitMidpointNewton::new();

        let mut state = y0;
        for _ in 0..5 {
            let out = integ.step(&state, 0.0, dt, &rhs).expect("step ok");
            state = out.new_state;
            assert!(
                state[0].x.is_finite(),
                "A-stability: result must remain finite"
            );
            assert!(
                state[0].x.abs() <= 1.0,
                "A-stability: |y| must not grow, got {}",
                state[0].x.abs()
            );
        }
    }

    // -- Test 4: 2nd-order convergence -----------------------------------

    #[test]
    fn test_second_order_convergence() {
        // dy/dt = -y, integrate to t = 1.0, halve dt and check the error
        // scales as dt^2.
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let t_end = 1.0_f64;
        let analytical = y0[0] * (-t_end).exp();
        let dts = [0.1_f64, 0.05, 0.025];
        let mut errors = Vec::new();

        for &dt in &dts {
            let mut integ = ImplicitMidpointNewton::new();
            let n_steps = (t_end / dt).round() as usize;
            let mut state = y0.clone();
            let mut t = 0.0;
            for _ in 0..n_steps {
                let out = integ
                    .step(&state, t, dt, &exponential_rhs)
                    .expect("step ok");
                state = out.new_state;
                t += dt;
            }
            let err = (state[0] - analytical).magnitude();
            errors.push(err);
        }

        // Ratio of consecutive errors should be ~4 (since dt halves).
        let order1 = (errors[0] / errors[1]).ln() / (dts[0] / dts[1]).ln();
        let order2 = (errors[1] / errors[2]).ln() / (dts[1] / dts[2]).ln();

        assert!(
            order1 > 1.7,
            "2nd-order convergence (1-2): got order {:.2}",
            order1
        );
        assert!(
            order2 > 1.7,
            "2nd-order convergence (2-3): got order {:.2}",
            order2
        );
    }

    // -- Test 5: |m| conservation over one Larmor period ----------------

    #[test]
    fn test_llg_larmor_norm_conservation() {
        // For a constant external field B along z and zero damping the
        // Landau-Lifshitz equation reduces to undamped precession at the
        // Larmor frequency omega = gamma B. The implicit midpoint rule
        // preserves the quadratic invariant |m|^2 exactly.
        let b_field = Vector3::new(0.0, 0.0, 1.0);
        let rhs = make_llg_rhs(b_field, 0.0);
        let m0 = vec![Vector3::new(1.0, 0.0, 0.0).normalize()];

        let omega = GAMMA * b_field.magnitude();
        let period = 2.0 * std::f64::consts::PI / omega;
        let n_steps = 200;
        let dt = period / n_steps as f64;
        let mut integ = ImplicitMidpointNewton::new();

        let mut state = m0;
        let mut t = 0.0;
        for _ in 0..n_steps {
            let out = integ.step(&state, t, dt, &rhs).expect("LLG step ok");
            state = out.new_state;
            t += dt;
            let mag = state[0].magnitude();
            assert!(
                (mag - 1.0).abs() < 1e-9,
                "|m| should stay on unit sphere, got |m| = {}",
                mag
            );
        }
    }

    // -- Test 6: energy-like behaviour on a symplectic problem ------------

    #[test]
    fn test_quadratic_invariant_preservation() {
        // For the harmonic oscillator dq/dt = p, dp/dt = -q the implicit
        // midpoint rule preserves the energy E = (q^2 + p^2)/2 to machine
        // precision. Here we encode (q, p) as two Vector3 entries.
        fn harmonic_rhs(state: &[Vector3<f64>], _t: f64) -> Vec<Vector3<f64>> {
            let half = state.len() / 2;
            let mut d = vec![Vector3::zero(); state.len()];
            for i in 0..half {
                d[i] = state[half + i];
                d[half + i] = state[i] * (-1.0);
            }
            d
        }

        let state0 = vec![
            Vector3::new(1.0, 0.0, 0.0), // q
            Vector3::new(0.0, 0.0, 0.0), // p
        ];
        let e0 = 0.5 * (state0[0].magnitude_squared() + state0[1].magnitude_squared());

        let mut integ = ImplicitMidpointNewton::new();
        let mut state = state0;
        let dt = 0.05_f64;
        for _ in 0..1000 {
            let out = integ.step(&state, 0.0, dt, &harmonic_rhs).expect("step ok");
            state = out.new_state;
        }
        let e_final = 0.5 * (state[0].magnitude_squared() + state[1].magnitude_squared());
        let drift = ((e_final - e0) / e0).abs();
        assert!(drift < 1e-10, "Energy drift too large: {:.3e}", drift);
    }

    // -- Test 7: stiff Dahlquist (where explicit RK4 would blow up) -----

    #[test]
    fn test_stiff_dahlquist_stable_at_large_dt() {
        // With lambda = 1e6, RK4 needs dt < 2.78e-6. The implicit method
        // should accept dt = 1e-3 and produce a tiny, finite y(t).
        let lambda = 1.0e6;
        let rhs = make_dahlquist_rhs(lambda);
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let dt = 1.0e-3;
        let mut integ = ImplicitMidpointNewton::new();

        let mut state = y0;
        for _ in 0..10 {
            let out = integ.step(&state, 0.0, dt, &rhs).expect("step ok");
            state = out.new_state;
            assert!(state[0].x.is_finite());
            assert!(state[0].x.abs() < 1.0, "should decay");
        }
    }

    // -- Test 8: damped precession matches analytical ---------------------

    #[test]
    fn test_damped_precession_analytical() {
        // Analytical solution for dy/dt = (-y.y, y.x, 0) - 0.1 y starting at
        // (1, 0, 0): y(t) = e^{-0.1 t} (cos t, sin t, 0).
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let t_end = 1.0_f64;
        let dt = 0.01_f64;
        let n_steps = (t_end / dt).round() as usize;
        let mut integ = ImplicitMidpointNewton::new();

        let mut state = y0;
        let mut t = 0.0;
        for _ in 0..n_steps {
            let out = integ
                .step(&state, t, dt, &damped_precession_rhs)
                .expect("damped precession step ok");
            state = out.new_state;
            t += dt;
        }

        let analytical = Vector3::new(
            (-0.1 * t_end).exp() * t_end.cos(),
            (-0.1 * t_end).exp() * t_end.sin(),
            0.0,
        );
        let err = (state[0] - analytical).magnitude();
        assert!(err < 5e-4, "damped precession error: {:.3e}", err);
    }

    // -- Test 9: Newton converges within a few iterations ---------------

    #[test]
    fn test_newton_converges_quickly() {
        // For a mild problem the Newton solve should always converge well
        // within the default 20 iterations; verify the residual reported by
        // the integrator is below the tolerance.
        let y0 = vec![Vector3::new(1.0, 0.5, -0.25)];
        let mut integ = ImplicitMidpointNewton::new();
        let out = integ
            .step(&y0, 0.0, 0.01, &exponential_rhs)
            .expect("step ok");
        let residual = out.error_estimate.expect("integrator must report residual");
        assert!(
            residual < 1e-12,
            "Newton residual {:.3e} should be below tol",
            residual
        );
    }

    // -- Test 10: returns Err when Newton fails to converge --------------

    #[test]
    fn test_newton_failure_returns_err() {
        // Force failure by demanding a tiny tolerance with a single Newton
        // iteration on a non-trivial problem.
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let mut integ = ImplicitMidpointNewton::new()
            .with_max_iter(1)
            .with_tol(1e-30);
        let result = integ.step(&y0, 0.0, 0.01, &damped_precession_rhs);
        assert!(
            result.is_err(),
            "Newton with max_iter=1 and impossible tol should fail"
        );
    }

    // -- Test 11: rejects degenerate fd_step ------------------------------

    #[test]
    fn test_rejects_tiny_fd_step() {
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let mut integ = ImplicitMidpointNewton::new().with_fd_step(1e-20);
        let result = integ.step(&y0, 0.0, 0.01, &exponential_rhs);
        assert!(result.is_err(), "tiny fd_step must be rejected");
    }

    // -- Test 12: Gauss elimination sanity check --------------------------

    #[test]
    fn test_gauss_solve_2x2() {
        // 2x + y = 5
        //  x + 3y = 10  =>  x = 1, y = 3
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let b = vec![5.0, 10.0];
        let x = gauss_solve(a, b, 2).expect("non-singular");
        assert!((x[0] - 1.0).abs() < 1e-12);
        assert!((x[1] - 3.0).abs() < 1e-12);
    }

    // -- Test 13: vector3 <-> flat round-trip ----------------------------

    #[test]
    fn test_vector3_flat_roundtrip() {
        let v = vec![Vector3::new(1.0, 2.0, 3.0), Vector3::new(-4.0, 5.5, -6.25)];
        let flat = vector3_to_flat(&v);
        assert_eq!(flat, vec![1.0, 2.0, 3.0, -4.0, 5.5, -6.25]);
        let back = flat_to_vector3(&flat);
        assert_eq!(back, v);
    }
}
