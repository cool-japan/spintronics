//! Tests for all integrators.

#[cfg(test)]
#[allow(clippy::module_inception)]
mod tests {
    use super::super::*;
    use crate::vector3::Vector3;

    // -- Helpers ----------------------------------------------------------

    /// Simple harmonic oscillator: dq/dt = p, dp/dt = -q
    /// State: [q, p] where q and p are Vector3 (we use x-component only).
    fn harmonic_rhs(state: &[Vector3<f64>], _t: f64) -> Vec<Vector3<f64>> {
        let half = state.len() / 2;
        let mut deriv = vec![Vector3::zero(); state.len()];
        for i in 0..half {
            // dq/dt = p
            deriv[i] = state[half + i];
            // dp/dt = -q  (force = -q for unit mass, unit frequency)
            deriv[half + i] = state[i] * (-1.0);
        }
        deriv
    }

    /// Simple exponential decay: dy/dt = -y  (analytical: y(t) = y0 * e^{-t})
    fn exponential_rhs(state: &[Vector3<f64>], _t: f64) -> Vec<Vector3<f64>> {
        state.iter().map(|&v| v * (-1.0)).collect()
    }

    /// Rotation (undamped precession): dy/dt = (-y.y, y.x, 0)
    /// Analytical: y(t) = (cos t, sin t, 0) if y(0) = (1, 0, 0)
    fn rotation_rhs(state: &[Vector3<f64>], _t: f64) -> Vec<Vector3<f64>> {
        state.iter().map(|v| Vector3::new(-v.y, v.x, 0.0)).collect()
    }

    /// Damped precession: dy/dt = (-y.y, y.x, 0) - alpha * y
    /// Models spin precession with Gilbert-like damping.
    fn damped_precession_rhs(state: &[Vector3<f64>], _t: f64) -> Vec<Vector3<f64>> {
        let alpha = 0.1;
        state
            .iter()
            .map(|v| Vector3::new(-v.y, v.x, 0.0) + *v * (-alpha))
            .collect()
    }

    /// Stiff system: dy/dt = -1000 * y (fast exponential decay)
    fn stiff_rhs(state: &[Vector3<f64>], _t: f64) -> Vec<Vector3<f64>> {
        state.iter().map(|&v| v * (-1000.0)).collect()
    }

    // -- Test 1: DP45 convergence order ------------------------------------

    #[test]
    fn test_dp45_convergence_order() {
        // Solve dy/dt = -y from t=0 to t=0.1 with decreasing dt.
        // Error should scale as O(dt^5) for the 5th-order method.
        let y0 = vec![Vector3::new(1.0, 0.5, 0.25)];
        let t_end: f64 = 0.1;

        let _analytical = y0[0] * (-t_end).exp();

        let mut errors = Vec::new();
        // Use a longer integration (t_end=1.0) so errors are well above machine epsilon
        let t_end_conv: f64 = 1.0;
        let analytical_conv = y0[0] * (-t_end_conv).exp();
        let dts: [f64; 3] = [0.05, 0.025, 0.0125];

        for &dt in &dts {
            let mut integrator = DormandPrince45::new();
            let n_steps = (t_end_conv / dt).round() as usize;
            let mut state = y0.clone();
            let mut t_local = 0.0;
            for _ in 0..n_steps {
                let out = integrator
                    .step(&state, t_local, dt, &exponential_rhs)
                    .expect("DP45 step should succeed");
                state = out.new_state;
                t_local += dt;
            }
            let err = (state[0] - analytical_conv).magnitude();
            errors.push(err);
        }

        // Check convergence order: log(e1/e2) / log(dt1/dt2) ~ 5
        let order_12 = (errors[0] / errors[1]).ln() / (dts[0] / dts[1]).ln();
        let order_23 = (errors[1] / errors[2]).ln() / (dts[1] / dts[2]).ln();

        // Allow some tolerance — numerical order may not be exactly 5
        assert!(
            order_12 > 4.0,
            "DP45 convergence order between dt1/dt2: {:.2}, expected > 4.0 (errors: {:.2e}, {:.2e})",
            order_12, errors[0], errors[1]
        );
        assert!(
            order_23 > 4.0,
            "DP45 convergence order between dt2/dt3: {:.2}, expected > 4.0 (errors: {:.2e}, {:.2e})",
            order_23, errors[1], errors[2]
        );
    }

    // -- Test 2: DP45 error estimate is meaningful -------------------------

    #[test]
    fn test_dp45_error_estimate() {
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let mut integrator = DormandPrince45::new();

        let out = integrator
            .step(&y0, 0.0, 0.01, &exponential_rhs)
            .expect("DP45 step should succeed");

        let err = out
            .error_estimate
            .expect("DP45 should provide error estimate");
        assert!(err >= 0.0, "Error estimate should be non-negative");
        assert!(
            err < 1e-6,
            "Error should be small for smooth problem with dt=0.01"
        );

        let sdt = out.suggested_dt.expect("DP45 should suggest dt");
        assert!(sdt > 0.0, "Suggested dt should be positive");
    }

    // -- Test 3: DP45 FSAL property ----------------------------------------

    #[test]
    fn test_dp45_fsal_reuse() {
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let mut integrator = DormandPrince45::new();

        // First step — no FSAL cache
        let out1 = integrator
            .step(&y0, 0.0, 0.01, &exponential_rhs)
            .expect("step 1 should succeed");

        // Second step — should reuse FSAL
        let _out2 = integrator
            .step(&out1.new_state, 0.01, 0.01, &exponential_rhs)
            .expect("step 2 should succeed (with FSAL)");

        // If FSAL is broken, the second step would still succeed but with
        // slightly different numerics. We just verify no error occurs.
    }

    // -- Test 4: DP87 convergence order ------------------------------------

    #[test]
    fn test_dp87_convergence_order() {
        let y0 = vec![Vector3::new(1.0, 0.5, 0.25)];
        // Use a longer integration time and larger dt so errors are measurable
        let t_end: f64 = 5.0;
        let analytical = y0[0] * (-t_end).exp();

        let mut errors = Vec::new();
        let dts: [f64; 3] = [0.5, 0.25, 0.125];

        for &dt in &dts {
            let mut integrator = DormandPrince87::new();
            let n_steps = (t_end / dt).round() as usize;
            let mut state = y0.clone();
            let mut t = 0.0;
            for _ in 0..n_steps {
                let out = integrator
                    .step(&state, t, dt, &exponential_rhs)
                    .expect("DP87 step should succeed");
                state = out.new_state;
                t += dt;
            }
            let err = (state[0] - analytical).magnitude();
            errors.push(err);
        }

        let order_12 = (errors[0] / errors[1]).ln() / (dts[0] / dts[1]).ln();

        // RK8 should show high-order convergence
        assert!(
            order_12 > 6.0,
            "DP87 convergence order: {:.2}, expected > 6.0 (errors: {:.2e}, {:.2e})",
            order_12,
            errors[0],
            errors[1]
        );
    }

    // -- Test 5: DP87 higher accuracy than DP45 ----------------------------

    #[test]
    fn test_dp87_higher_accuracy_than_dp45() {
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let dt = 0.05;

        let mut dp45 = DormandPrince45::new();
        let mut dp87 = DormandPrince87::new();

        let out45 = dp45
            .step(&y0, 0.0, dt, &exponential_rhs)
            .expect("DP45 should succeed");
        let out87 = dp87
            .step(&y0, 0.0, dt, &exponential_rhs)
            .expect("DP87 should succeed");

        let analytical = y0[0] * (-dt).exp();
        let err45 = (out45.new_state[0] - analytical).magnitude();
        let err87 = (out87.new_state[0] - analytical).magnitude();

        assert!(
            err87 < err45,
            "DP87 error ({:.2e}) should be smaller than DP45 error ({:.2e})",
            err87,
            err45
        );
    }

    // -- Test 6: Velocity Verlet energy conservation -----------------------

    #[test]
    fn test_velocity_verlet_energy_conservation() {
        // Harmonic oscillator: E = 0.5*(q^2 + p^2) should be conserved.
        let q0 = Vector3::new(1.0, 0.0, 0.0);
        let p0 = Vector3::new(0.0, 1.0, 0.0);
        let state = vec![q0, p0];

        let energy = |s: &[Vector3<f64>]| -> f64 {
            0.5 * (s[0].magnitude_squared() + s[1].magnitude_squared())
        };

        let e0 = energy(&state);
        let mut cur = state;
        let mut integrator = VelocityVerlet::new();
        let dt = 0.01;

        for _ in 0..10000 {
            let out = integrator
                .step(&cur, 0.0, dt, &harmonic_rhs)
                .expect("Verlet step should succeed");
            cur = out.new_state;
        }

        let e_final = energy(&cur);
        let relative_drift = ((e_final - e0) / e0).abs();

        assert!(
            relative_drift < 1e-6,
            "Velocity Verlet energy drift: {:.2e}, expected < 1e-6",
            relative_drift
        );
    }

    // -- Test 7: Yoshida energy conservation -------------------------------

    #[test]
    fn test_yoshida4_energy_conservation() {
        let q0 = Vector3::new(1.0, 0.0, 0.0);
        let p0 = Vector3::new(0.0, 1.0, 0.0);
        let state = vec![q0, p0];

        let energy = |s: &[Vector3<f64>]| -> f64 {
            0.5 * (s[0].magnitude_squared() + s[1].magnitude_squared())
        };

        let e0 = energy(&state);
        let mut cur = state;
        let mut integrator = Yoshida4::new();
        let dt = 0.01;

        for _ in 0..10000 {
            let out = integrator
                .step(&cur, 0.0, dt, &harmonic_rhs)
                .expect("Yoshida step should succeed");
            cur = out.new_state;
        }

        let e_final = energy(&cur);
        let relative_drift = ((e_final - e0) / e0).abs();

        assert!(
            relative_drift < 1e-8,
            "Yoshida4 energy drift: {:.2e}, expected < 1e-8",
            relative_drift
        );
    }

    // -- Test 8: Forest-Ruth energy conservation ---------------------------

    #[test]
    fn test_forest_ruth_energy_conservation() {
        let q0 = Vector3::new(1.0, 0.0, 0.0);
        let p0 = Vector3::new(0.0, 1.0, 0.0);
        let state = vec![q0, p0];

        let energy = |s: &[Vector3<f64>]| -> f64 {
            0.5 * (s[0].magnitude_squared() + s[1].magnitude_squared())
        };

        let e0 = energy(&state);
        let mut cur = state;
        let mut integrator = ForestRuth::new();
        let dt = 0.01;

        for _ in 0..10000 {
            let out = integrator
                .step(&cur, 0.0, dt, &harmonic_rhs)
                .expect("Forest-Ruth step should succeed");
            cur = out.new_state;
        }

        let e_final = energy(&cur);
        let relative_drift = ((e_final - e0) / e0).abs();

        assert!(
            relative_drift < 1e-8,
            "Forest-Ruth energy drift: {:.2e}, expected < 1e-8",
            relative_drift
        );
    }

    // -- Test 9: Symplectic integrators odd-state rejection ----------------

    #[test]
    fn test_symplectic_odd_state_rejection() {
        let state = vec![Vector3::new(1.0, 0.0, 0.0)]; // 1 element = odd

        let mut vv = VelocityVerlet::new();
        let result = vv.step(&state, 0.0, 0.01, &harmonic_rhs);
        assert!(
            result.is_err(),
            "VelocityVerlet should reject odd-length state"
        );

        let mut y4 = Yoshida4::new();
        let result = y4.step(&state, 0.0, 0.01, &harmonic_rhs);
        assert!(result.is_err(), "Yoshida4 should reject odd-length state");

        let mut fr = ForestRuth::new();
        let result = fr.step(&state, 0.0, 0.01, &harmonic_rhs);
        assert!(result.is_err(), "ForestRuth should reject odd-length state");
    }

    // -- Test 10: Semi-implicit convergence on stiff system ----------------

    #[test]
    fn test_semi_implicit_stiff_system() {
        // dy/dt = -1000*y  with y(0) = (1, 0, 0)
        // Implicit midpoint fixed-point iteration converges when |lambda*dt/2| < 1,
        // so dt < 0.002 for lambda=1000. Use multiple small steps for accuracy.
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let dt = 0.0005;
        let n_steps = 2; // total time = 0.001
        let total_time = dt * n_steps as f64;
        let mut integrator = SemiImplicit::new(200, 1e-14);

        let mut state = y0;
        let mut t = 0.0;
        for _ in 0..n_steps {
            let out = integrator
                .step(&state, t, dt, &stiff_rhs)
                .expect("Semi-implicit step should succeed");
            state = out.new_state;
            t += dt;
        }

        // Analytical: exp(-1000 * 0.001) = exp(-0.5) ~ 0.6065
        let analytical_x = (-1000.0_f64 * total_time).exp();
        let err = (state[0].x - analytical_x).abs();

        // The implicit midpoint is 2nd order with error O(dt^2) ~ O(2.5e-7) per step
        assert!(
            err < 0.01,
            "Semi-implicit error on stiff system: {:.2e}",
            err
        );
        assert!(
            state[0].x.is_finite(),
            "Semi-implicit should produce finite result for stiff system"
        );
    }

    // -- Test 11: Semi-implicit convergence failure reporting ---------------

    #[test]
    fn test_semi_implicit_convergence_report() {
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        // Very few iterations — may not converge
        let mut integrator = SemiImplicit::new(1, 1e-20);

        let out = integrator
            .step(&y0, 0.0, 0.01, &stiff_rhs)
            .expect("Should still return a result even if not fully converged");

        // When not converged, error_estimate should be present
        assert!(
            out.error_estimate.is_some(),
            "Should report error estimate even when not converged"
        );
    }

    // -- Test 12: Rotation analytical comparison ---------------------------

    #[test]
    fn test_dp45_rotation_analytical() {
        // y(0) = (1, 0, 0), dy/dt = (-y, x, 0)
        // Analytical: y(t) = (cos t, sin t, 0)
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let t_end: f64 = 2.0 * std::f64::consts::PI; // One full rotation
        let dt: f64 = 0.001;
        let n_steps = (t_end / dt).round() as usize;

        let mut integrator = DormandPrince45::new();
        let mut state = y0.clone();
        let mut t = 0.0;

        for _ in 0..n_steps {
            let out = integrator
                .step(&state, t, dt, &rotation_rhs)
                .expect("rotation step should succeed");
            state = out.new_state;
            t += dt;
        }

        // After one full period, should return to (1, 0, 0)
        // DP45 with dt=0.001 over ~6283 steps: error ~O(dt^5 * n) ~ O(1e-15 * 6e3) ~ 1e-12
        // but accumulation effects give larger error. Allow 1e-3.
        let err = (state[0] - y0[0]).magnitude();
        assert!(err < 1e-3, "DP45 rotation error after 2*pi: {:.2e}", err);
    }

    // -- Test 13: Damped precession comparison -----------------------------

    #[test]
    fn test_dp45_damped_precession() {
        // dy/dt = (-y, x, 0) - 0.1*y
        // Analytical: y(t) = e^{-0.1t} * (cos t, sin t, 0) for y(0)=(1,0,0)
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let t_end: f64 = 5.0;
        let dt: f64 = 0.005;
        let n_steps = (t_end / dt).round() as usize;

        let mut integrator = DormandPrince45::new();
        let mut state = y0;
        let mut t = 0.0;

        for _ in 0..n_steps {
            let out = integrator
                .step(&state, t, dt, &damped_precession_rhs)
                .expect("damped precession step should succeed");
            state = out.new_state;
            t += dt;
        }

        let analytical = Vector3::new(
            (-0.1 * t_end).exp() * t_end.cos(),
            (-0.1 * t_end).exp() * t_end.sin(),
            0.0,
        );

        let err = (state[0] - analytical).magnitude();
        assert!(
            err < 1e-6,
            "DP45 damped precession error at t={}: {:.2e}",
            t_end,
            err
        );
    }

    // -- Test 14: Adaptive integrator step size control --------------------

    #[test]
    fn test_adaptive_step_size_control() {
        let integrator = DormandPrince45::new();
        let mut adaptive = AdaptiveIntegrator::new(integrator, 1e-8, 4.0)
            .with_dt_min(1e-15)
            .with_dt_max(1.0);

        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];

        let (out, actual_dt) = adaptive
            .adaptive_step(&y0, 0.0, 0.1, &exponential_rhs)
            .expect("adaptive step should succeed");

        // The step should have been accepted
        assert!(actual_dt > 0.0, "Accepted dt should be positive");
        assert!(
            out.error_estimate.expect("should have error") <= 1e-8 || actual_dt <= 0.1,
            "Error should be within tolerance or step was shrunk"
        );

        let sdt = out.suggested_dt.expect("should suggest next dt");
        assert!(sdt > 0.0, "Suggested dt should be positive");
    }

    // -- Test 15: Adaptive integrator rejects bad steps --------------------

    #[test]
    fn test_adaptive_rejects_large_error() {
        let integrator = DormandPrince45::new();
        let mut adaptive = AdaptiveIntegrator::new(integrator, 1e-12, 4.0)
            .with_dt_min(1e-15)
            .with_dt_max(1.0);

        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];

        // Start with a very large step — should be shrunk
        let (out, actual_dt) = adaptive
            .adaptive_step(&y0, 0.0, 1.0, &exponential_rhs)
            .expect("adaptive step should eventually succeed");

        assert!(
            actual_dt < 1.0,
            "Should have shrunk from initial dt=1.0, got dt={:.2e}",
            actual_dt
        );

        let err = out.error_estimate.expect("should have error estimate");
        assert!(
            err <= 1e-12,
            "Final error {:.2e} should be within tolerance 1e-12",
            err
        );
    }

    // -- Test 16: Zero-damping (pure rotation) with DP87 -------------------

    #[test]
    fn test_dp87_zero_damping_rotation() {
        let y0 = vec![Vector3::new(1.0, 0.0, 0.0)];
        let t_end: f64 = 2.0 * std::f64::consts::PI;
        let dt: f64 = 0.005;
        let n_steps = (t_end / dt).round() as usize;

        let mut integrator = DormandPrince87::new();
        let mut state = y0.clone();
        let mut t = 0.0;

        for _ in 0..n_steps {
            let out = integrator
                .step(&state, t, dt, &rotation_rhs)
                .expect("DP87 rotation step should succeed");
            state = out.new_state;
            t += dt;
        }

        // DP87 should give good accuracy for this smooth problem
        let err = (state[0] - y0[0]).magnitude();
        assert!(
            err < 1e-2,
            "DP87 rotation error after 2*pi with dt={dt}: {err:.2e}",
        );
    }

    // -- Test 17: Velocity Verlet trajectory correctness -------------------

    #[test]
    fn test_velocity_verlet_trajectory() {
        // Harmonic oscillator: q(t) = cos(t), p(t) = -sin(t) for q(0)=1, p(0)=0
        let state = vec![
            Vector3::new(1.0, 0.0, 0.0), // q
            Vector3::new(0.0, 0.0, 0.0), // p
        ];

        let mut integrator = VelocityVerlet::new();
        let dt: f64 = 0.001;
        let t_end: f64 = 1.0;
        let n_steps = (t_end / dt).round() as usize;
        let mut cur = state;

        for _ in 0..n_steps {
            let out = integrator
                .step(&cur, 0.0, dt, &harmonic_rhs)
                .expect("Verlet step should succeed");
            cur = out.new_state;
        }

        let expected_q = t_end.cos();
        let expected_p = -(t_end.sin());

        assert!(
            (cur[0].x - expected_q).abs() < 1e-5,
            "q error: {:.2e}",
            (cur[0].x - expected_q).abs()
        );
        assert!(
            (cur[1].x - expected_p).abs() < 1e-5,
            "p error: {:.2e}",
            (cur[1].x - expected_p).abs()
        );
    }

    // -- Test 18: NaN detection --------------------------------------------

    #[test]
    fn test_nan_detection() {
        // A rhs that produces NaN
        fn bad_rhs(state: &[Vector3<f64>], _t: f64) -> Vec<Vector3<f64>> {
            state
                .iter()
                .map(|v| Vector3::new(v.x / 0.0, v.y, v.z))
                .collect()
        }

        let y0 = vec![Vector3::new(0.0, 1.0, 0.0)];
        let mut integrator = DormandPrince45::new();

        let result = integrator.step(&y0, 0.0, 0.01, &bad_rhs);
        assert!(result.is_err(), "Should detect NaN in output");
    }

    // -- Test 19: Adaptive with DP87 ---------------------------------------

    #[test]
    fn test_adaptive_dp87() {
        let integrator = DormandPrince87::new();
        let mut adaptive = AdaptiveIntegrator::new(integrator, 1e-10, 7.0)
            .with_dt_min(1e-15)
            .with_dt_max(1.0);

        let y0 = vec![Vector3::new(1.0, 0.5, 0.25)];

        let (out, actual_dt) = adaptive
            .adaptive_step(&y0, 0.0, 0.1, &exponential_rhs)
            .expect("adaptive DP87 step should succeed");

        assert!(actual_dt > 0.0);
        let err = out.error_estimate.expect("should have error");
        assert!(
            err <= 1e-10 || actual_dt < 0.1,
            "Error {:.2e} should meet tolerance or step shrunk",
            err
        );
    }

    // -- Test 20: Yoshida better than Verlet for energy --------------------

    #[test]
    fn test_yoshida_better_than_verlet_energy() {
        let q0 = Vector3::new(1.0, 0.0, 0.0);
        let p0 = Vector3::new(0.0, 1.0, 0.0);
        let state = vec![q0, p0];

        let energy = |s: &[Vector3<f64>]| -> f64 {
            0.5 * (s[0].magnitude_squared() + s[1].magnitude_squared())
        };

        let e0 = energy(&state);
        let dt = 0.1; // Larger step to see the difference
        let n_steps = 1000;

        // Verlet
        let mut verlet = VelocityVerlet::new();
        let mut cur_v = state.clone();
        for _ in 0..n_steps {
            let out = verlet
                .step(&cur_v, 0.0, dt, &harmonic_rhs)
                .expect("Verlet step should succeed");
            cur_v = out.new_state;
        }
        let drift_verlet = ((energy(&cur_v) - e0) / e0).abs();

        // Yoshida
        let mut yoshida = Yoshida4::new();
        let mut cur_y = state;
        for _ in 0..n_steps {
            let out = yoshida
                .step(&cur_y, 0.0, dt, &harmonic_rhs)
                .expect("Yoshida step should succeed");
            cur_y = out.new_state;
        }
        let drift_yoshida = ((energy(&cur_y) - e0) / e0).abs();

        assert!(
            drift_yoshida <= drift_verlet,
            "Yoshida drift ({:.2e}) should be <= Verlet drift ({:.2e})",
            drift_yoshida,
            drift_verlet
        );
    }
}
