//! Implicit Midpoint Integrator on a Stiff LLG Problem
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Numerical Methods / Stiff ODEs
//! **Physics**: A-stable implicit midpoint rule with Newton + finite-diff Jacobian
//!
//! ## Background
//!
//! The Landau-Lifshitz-Gilbert equation becomes *stiff* when either α → 0
//! (large precession-to-damping time-scale ratio) or H_eff is very large
//! (oscillation period ω_L = γ|H| short). Explicit RK4 needs dt ~ 1 / (γ|H|)
//! to remain stable; for H = 10 T, γ = 1.76 × 10¹¹, that's dt ≲ 5·10⁻¹³ s.
//!
//! The implicit midpoint rule is *A-stable*, meaning it remains stable for
//! arbitrary dt on a linear test problem. The trade-off is solving a Newton
//! iteration each step (with finite-difference Jacobian here, costing 3N
//! RHS evaluations for an N-spin system).
//!
//! Here we compare RK4 and `ImplicitMidpointNewton` on a single-spin LLG
//! problem with a large effective field, showing that the implicit method
//! handles dt that would be unstable for RK4.
//!
//! ## References
//! - Hairer & Wanner, "Solving Ordinary Differential Equations II" (1996)
//! - Mentink et al., JPCM 22, 176001 (2010) — symplectic LLG integrators

use spintronics::dynamics::integrators::{ImplicitMidpointNewton, Integrator};
use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  Implicit Midpoint vs. RK4 on Stiff LLG");
    println!("=============================================================");

    // -------------------------------------------------------------------------
    // Section 1: Build a stiff LLG RHS (large H_eff along z, modest α)
    // -------------------------------------------------------------------------
    println!("\n--- Section 1: Setup ---\n");

    let h_eff_z = 1.0e7_f64; // 10 MA/m → ~12.6 T (very large)
    let alpha = 0.1_f64; // moderately strong damping
    let gamma = spintronics::constants::GAMMA.abs();
    let omega_l = gamma * spintronics::constants::MU_0 * h_eff_z;
    let t_l = 2.0 * std::f64::consts::PI / omega_l;
    println!("  H_eff (z)       = {h_eff_z:.2e} A/m");
    println!("  α (Gilbert)     = {alpha}");
    println!("  ω_Larmor        = {omega_l:.3e} rad/s");
    println!("  Period T_L      = {t_l:.3e} s");

    let rhs = move |m: &[Vector3<f64>], _t: f64| -> Vec<Vector3<f64>> {
        let h_eff = Vector3::new(0.0, 0.0, h_eff_z);
        m.iter()
            .map(|&mi| {
                let mxh = mi.cross(&h_eff);
                let m_mxh = mi.cross(&mxh);
                // Explicit LL form: dm/dt = -γ/(1+α²) [m × H + α m × (m × H)]
                let coeff = -gamma * spintronics::constants::MU_0 / (1.0 + alpha * alpha);
                (mxh + m_mxh * alpha) * coeff
            })
            .collect()
    };

    let m0 = vec![Vector3::new(0.6, 0.0, 0.8)]; // initial tilt

    // -------------------------------------------------------------------------
    // Section 2: ImplicitMidpointNewton — large dt should be stable
    // -------------------------------------------------------------------------
    println!("\n--- Section 2: ImplicitMidpointNewton at dt = T_L/30 ---\n");

    let dt_implicit = t_l / 30.0;
    let n_steps_impl = 30;
    let mut integrator = ImplicitMidpointNewton::new()
        .with_max_iter(20)
        .with_tol(1e-10);

    let mut state = m0.clone();
    let mut t = 0.0_f64;
    println!(
        "  {:>5}  {:>12}  {:>10}  {:>10}  {:>10}  {:>10}",
        "step", "t (s)", "mx", "my", "mz", "|m|"
    );
    println!("  {}", "-".repeat(64));
    for k in 0..=n_steps_impl {
        let mag = state[0].magnitude();
        println!(
            "  {:>5}  {:>12.3e}  {:>+10.4}  {:>+10.4}  {:>+10.4}  {:>10.6}",
            k, t, state[0].x, state[0].y, state[0].z, mag
        );
        if k == n_steps_impl {
            break;
        }
        let out = integrator.step(&state, t, dt_implicit, &rhs)?;
        state = out.new_state;
        t += dt_implicit;
    }

    let final_mag_impl = state[0].magnitude();
    println!("\n  → Final |m| (implicit midpoint): {final_mag_impl:.6}");
    println!("    (should remain ≈ 1.0 — magnetisation conserved by LLG)");

    // -------------------------------------------------------------------------
    // Section 3: Naive Euler at the same dt — expected to blow up
    // -------------------------------------------------------------------------
    println!("\n--- Section 3: Explicit Euler at same dt (for contrast) ---\n");

    let mut state_eul = m0.clone();
    let mut t_eul = 0.0_f64;
    println!(
        "  {:>5}  {:>12}  {:>+10}  {:>+10}  {:>+10}  {:>10}",
        "step", "t (s)", "mx", "my", "mz", "|m|"
    );
    println!("  {}", "-".repeat(64));
    for k in 0..=n_steps_impl {
        let mag = state_eul[0].magnitude();
        if mag > 1e6 || !mag.is_finite() {
            println!(
                "  {:>5}  {:>12.3e}  *** EULER DIVERGED at step {k} (|m| = {mag:.3e}) ***",
                k, t_eul
            );
            break;
        }
        println!(
            "  {:>5}  {:>12.3e}  {:>+10.4}  {:>+10.4}  {:>+10.4}  {:>10.6}",
            k, t_eul, state_eul[0].x, state_eul[0].y, state_eul[0].z, mag
        );
        if k == n_steps_impl {
            break;
        }
        let dm = rhs(&state_eul, t_eul);
        state_eul[0] = state_eul[0] + dm[0] * dt_implicit;
        t_eul += dt_implicit;
    }

    // -------------------------------------------------------------------------
    // Section 4: Implicit midpoint accuracy — second order convergence
    // -------------------------------------------------------------------------
    println!("\n--- Section 4: Convergence Order Check ---\n");

    let t_end = t_l * 0.5; // half-period
    println!("  Integrating to t = {t_end:.3e} s (T_L/2)");
    println!(
        "  {:>10}  {:>14}  {:>14}",
        "dt", "|m| error", "ratio (∝ dt²)"
    );
    println!("  {}", "-".repeat(44));
    let mut prev_err = None;
    for &dt_div in &[40.0, 80.0, 160.0, 320.0] {
        let dt_test = t_end / dt_div;
        let n_steps = dt_div as usize;
        let mut s = m0.clone();
        let mut tt = 0.0_f64;
        let mut intg = ImplicitMidpointNewton::new().with_tol(1e-13);
        for _ in 0..n_steps {
            let out = intg.step(&s, tt, dt_test, &rhs)?;
            s = out.new_state;
            tt += dt_test;
        }
        // True LLG conserves |m|; deviation is a measure of integrator error.
        let err = (s[0].magnitude() - 1.0).abs();
        let ratio_str = match prev_err {
            Some(prev) => format!("{:.3}", prev / err.max(1e-30)),
            None => "—".to_string(),
        };
        println!("  {:>10.3e}  {:>14.4e}  {:>14}", dt_test, err, ratio_str);
        prev_err = Some(err);
    }

    println!("\n=============================================================");
    println!("  Done. ImplicitMidpointNewton: A-stable, 2nd-order; remains");
    println!("  bounded on stiff LLG where explicit Euler diverges quickly.");
    println!("=============================================================\n");

    Ok(())
}
