//! Automatic Differentiation for Spintronic Parameter Fitting
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Machine Learning / Optimisation
//! **Physics**: Ferromagnetic resonance (Kittel formula), Gilbert damping,
//!               FMR linewidth fitting
//!
//! ## Background
//!
//! Reverse-mode automatic differentiation (AD) tracks every scalar operation
//! on a Wengert-list (tape) and propagates gradients backwards (back-prop)
//! from the loss to each parameter.  This example demonstrates the spintronics
//! AD engine in three scenarios:
//!
//! 1. **Basic tape usage** — illustrate `Var` arithmetic and backward pass.
//! 2. **Damping-parameter fitting** — fit the Gilbert damping constant α from
//!    synthetic FMR linewidth data using the Kittel formula and Adam.
//! 3. **Optimizer comparison** — compare SGD, Adam, and L-BFGS on a simple
//!    2D quadratic loss landscape.
//!
//! ## Gilbert damping and FMR linewidth
//!
//! The half-power FMR linewidth (in field units) is related to the damping α:
//!
//! ```text
//! ΔH = (2α / |γ| μ₀) · ω₀
//! ```
//!
//! where ω₀ is the Kittel resonance frequency.  Given ω₀(H), the linewidth
//! is linear in α, so minimising `Σ (predicted_ΔH - measured_ΔH)²` over α
//! recovers the true damping constant.
//!
//! ## References
//!
//! - Baydin et al., "Automatic Differentiation in Machine Learning: a Survey",
//!   *JMLR* **18**, 1 (2018)
//! - Kingma & Ba, "Adam: A Method for Stochastic Optimization", *ICLR* 2015
//! - Liu & Nocedal, "On the limited memory BFGS method",
//!   *Math. Prog.* **45**, 503 (1989)
//! - Gilbert & Nocedal, "Global convergence properties of conjugate gradient methods",
//!   *SIAM J. Optim.* **2**, 21 (1992)

#![cfg(feature = "autodiff")]

use std::f64::consts::PI;

use spintronics::autodiff::physics_fns::kittel_frequency_diff;
use spintronics::prelude::*;

// Physical constants used throughout this example.
// These are re-exported via the prelude:
//   GAMMA = gyromagnetic ratio |γ| [rad/(s·T)] ≈ 1.7608e11
//   MU_0  = vacuum permeability [H/m] = 4π×10⁻⁷

/// True Gilbert damping constant for YIG (target of the fitting).
const ALPHA_TRUE: f64 = 3.0e-5;

/// Initial guess for the damping — off by a factor of ~3.
const ALPHA_INIT: f64 = 1.0e-4;

/// Saturation magnetisation of YIG [A/m] (fixed, not fitted here).
const MS_YIG: f64 = 1.4e5; // ~140 kA/m

/// Number of synthetic FMR field points.
const N_FIELD: usize = 10;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    // ─────────────────────────────────────────────────────────────────────────
    // Section 1: Demonstrate the tape
    // ─────────────────────────────────────────────────────────────────────────
    println!("======================================================");
    println!(" Section 1: Tape arithmetic and reverse-mode AD");
    println!("======================================================");

    // Compute z = x·y + sin(x) and differentiate w.r.t. x and y.
    //   dz/dx = y + cos(x)
    //   dz/dy = x
    let tape = Tape::new();
    let x = Var::leaf(&tape, 3.0_f64);
    let y = Var::leaf(&tape, 2.0_f64);
    let z = x * y + x.sin();

    tape.backward(z);

    let dz_dx_ad = x.grad();
    let dz_dy_ad = y.grad();

    // Analytic reference.
    let dz_dx_ref = y.value() + x.value().cos(); // y + cos(x)
    let dz_dy_ref = x.value(); // x

    println!("z = x·y + sin(x)  at  x=3, y=2:");
    println!(
        "  forward:  z           = {:.6}  (ref: {:.6})",
        z.value(),
        3.0_f64 * 2.0 + 3.0_f64.sin()
    );
    println!(
        "  backward: dz/dx (AD)  = {:.6}  (analytic: {:.6})",
        dz_dx_ad, dz_dx_ref
    );
    println!(
        "  backward: dz/dy (AD)  = {:.6}  (analytic: {:.6})",
        dz_dy_ad, dz_dy_ref
    );
    println!(
        "  error dz/dx = {:.2e}  |  error dz/dy = {:.2e}",
        (dz_dx_ad - dz_dx_ref).abs(),
        (dz_dy_ad - dz_dy_ref).abs()
    );

    // Record tape statistics.
    println!("Tape: {} values, {} ops", tape.values_len(), tape.ops_len());

    // ─────────────────────────────────────────────────────────────────────────
    // Section 2: Fit the Gilbert damping constant α
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n======================================================");
    println!(" Section 2: Fitting Gilbert damping α from FMR linewidth");
    println!("======================================================");

    // Build the synthetic "experimental" data.
    // H ranges from 50 kA/m to 500 kA/m (10 field values).
    let h_fields: Vec<f64> = (0..N_FIELD)
        .map(|i| 50e3 + (450e3 / (N_FIELD - 1) as f64) * i as f64) // A/m
        .collect();

    // "True" FMR linewidth for each field value:
    //   ω₀(H) = |γ|·μ₀·√(H·(H+M_s))   (Kittel)
    //   ΔH(H) = 2α·ω₀ / (|γ|·μ₀)      [A/m]
    let measured_dh: Vec<f64> = h_fields
        .iter()
        .map(|&h| {
            let omega_0 = GAMMA.abs() * MU_0 * (h * (h + MS_YIG)).sqrt();
            2.0 * ALPHA_TRUE * omega_0 / (GAMMA.abs() * MU_0)
        })
        .collect();

    println!("Synthetic FMR linewidth data ({} field points):", N_FIELD);
    println!(
        "{:>10}  {:>14}  {:>14}",
        "H (kA/m)", "ω₀ (GHz·2π)", "ΔH_meas (A/m)"
    );
    for (h, dh) in h_fields.iter().zip(measured_dh.iter()) {
        let omega0 = GAMMA.abs() * MU_0 * (h * (h + MS_YIG)).sqrt();
        println!(
            "{:>10.1}  {:>14.4}  {:>14.2}",
            h * 1e-3,
            omega0 / (2.0 * PI * 1e9),
            dh
        );
    }

    println!("\nFitting α with Adam (lr=1e-5, 2000 iterations):");
    println!("  true α  = {:.4e}", ALPHA_TRUE);
    println!(
        "  init α  = {:.4e}  (× {:.1} error)",
        ALPHA_INIT,
        ALPHA_INIT / ALPHA_TRUE
    );

    // Build the fitter.  Single parameter: α (scalar).
    let fitter = ParameterFitter::new_adam(
        vec![ALPHA_INIT],
        1e-5,  // lr — small because α lives at ~10⁻⁵
        2000,  // max_iter
        1e-20, // tol (tight — run to convergence)
    );

    // The loss closure builds the Kittel omega for each field and computes
    // ΔH = 2α·ω₀ / (|γ|·μ₀), then minimises Σ (ΔH_pred - ΔH_meas)².
    let h_fields_clone = h_fields.clone();
    let measured_dh_clone = measured_dh.clone();

    let result = fitter.fit(|tape, leaves| {
        let alpha_var = leaves[0]; // the single fitted parameter

        let mut loss = Var::leaf(tape, 0.0_f64);
        for (h_val, &dh_meas) in h_fields_clone.iter().zip(measured_dh_clone.iter()) {
            // Kittel frequency (differentiable w.r.t. ms_var, but ms is fixed here).
            let ms_var = Var::leaf(tape, MS_YIG);
            let h_var = Var::leaf(tape, *h_val);
            let omega0 = kittel_frequency_diff(tape, ms_var, h_var);

            // Predicted linewidth: ΔH = 2α·ω₀ / (|γ|·μ₀)
            let coeff = 2.0 / (GAMMA.abs() * MU_0);
            let dh_pred = alpha_var * omega0 * coeff;

            // Squared residual
            let residual = dh_pred - dh_meas;
            let sq = residual * residual;
            loss = loss + sq;
        }
        loss
    });

    let alpha_fitted = result.final_params[0];
    let rel_error = (alpha_fitted - ALPHA_TRUE).abs() / ALPHA_TRUE;

    println!("\nFit result:");
    println!("  fitted α    = {:.6e}", alpha_fitted);
    println!("  true α      = {:.6e}", ALPHA_TRUE);
    println!(
        "  relative err = {:.2e}  ({:.4}%)",
        rel_error,
        rel_error * 100.0
    );
    println!("  final loss   = {:.4e}", result.final_loss);
    println!("  n_iterations = {}", result.n_iterations);
    println!("  converged    = {}", result.converged);

    // Show loss decay over first/last few iterations.
    if result.loss_history.len() >= 4 {
        println!("\nLoss history:");
        let n = result.loss_history.len();
        for (i, &l) in result.loss_history.iter().enumerate().take(3) {
            println!("  iter {:4}: loss = {:.4e}", i + 1, l);
        }
        println!("  ...");
        for (i, &l) in result.loss_history[n - 3..].iter().enumerate() {
            println!("  iter {:4}: loss = {:.4e}", n - 2 + i, l);
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Section 3: Optimizer comparison on 2D quadratic
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n======================================================");
    println!(" Section 3: Optimizer comparison  f(x,y) = x² + 4y²");
    println!("======================================================");
    println!("  Starting point: (3.0, 2.0)   |  Target: (0.0, 0.0)");
    println!("  (Elongated bowl — Hessian eigenvalues 1 and 4)");
    println!();

    // Build the loss closure separately for each optimizer.
    // The `ParameterFitter::fit` signature requires a higher-ranked lifetime
    // (`for<'a> Fn(&'a Tape, &[Var<'a>]) -> Var<'a>`), so each call gets its
    // own monomorphised closure.

    // ── SGD ──
    let sgd_fitter = ParameterFitter::new_sgd(vec![3.0, 2.0], 0.1, 0.0, 200, 1e-12);
    let sgd_result = sgd_fitter.fit(|_tape, leaves| {
        let x_var = leaves[0];
        let y_var = leaves[1];
        x_var * x_var + y_var * y_var * 4.0
    });
    let sgd_norm = (sgd_result.final_params[0].powi(2) + sgd_result.final_params[1].powi(2)).sqrt();

    // ── Adam ──
    let adam_fitter = ParameterFitter::new_adam(vec![3.0, 2.0], 0.01, 200, 1e-12);
    let adam_result = adam_fitter.fit(|_tape, leaves| {
        let x_var = leaves[0];
        let y_var = leaves[1];
        x_var * x_var + y_var * y_var * 4.0
    });
    let adam_norm =
        (adam_result.final_params[0].powi(2) + adam_result.final_params[1].powi(2)).sqrt();

    // ── L-BFGS ──
    let lbfgs_fitter = ParameterFitter::new_lbfgs(vec![3.0, 2.0], 0.1, 10, 50, 1e-12);
    let lbfgs_result = lbfgs_fitter.fit(|_tape, leaves| {
        let x_var = leaves[0];
        let y_var = leaves[1];
        x_var * x_var + y_var * y_var * 4.0
    });
    let lbfgs_norm =
        (lbfgs_result.final_params[0].powi(2) + lbfgs_result.final_params[1].powi(2)).sqrt();

    // Print comparison table.
    println!(
        "{:<8}  {:>8}  {:>10}  {:>12}  {:>10}",
        "Optimizer", "lr", "n_iters", "||(x,y)||", "converged"
    );
    println!("{}", "-".repeat(58));
    println!(
        "{:<8}  {:>8}  {:>10}  {:>12.4e}  {:>10}",
        "SGD", "0.1", sgd_result.n_iterations, sgd_norm, sgd_result.converged
    );
    println!(
        "{:<8}  {:>8}  {:>10}  {:>12.4e}  {:>10}",
        "Adam", "0.01", adam_result.n_iterations, adam_norm, adam_result.converged
    );
    println!(
        "{:<8}  {:>8}  {:>10}  {:>12.4e}  {:>10}",
        "L-BFGS", "0.1", lbfgs_result.n_iterations, lbfgs_norm, lbfgs_result.converged
    );

    // Final parameter values.
    println!("\nFinal parameters:");
    println!(
        "  SGD:    x={:.4e},  y={:.4e}",
        sgd_result.final_params[0], sgd_result.final_params[1]
    );
    println!(
        "  Adam:   x={:.4e},  y={:.4e}",
        adam_result.final_params[0], adam_result.final_params[1]
    );
    println!(
        "  L-BFGS: x={:.4e},  y={:.4e}",
        lbfgs_result.final_params[0], lbfgs_result.final_params[1]
    );

    // Physical insight.
    println!("\nPhysical insight:");
    println!("  - L-BFGS uses curvature information → converges fastest on smooth problems.");
    println!("  - Adam adapts per-parameter learning rates → robust on ill-conditioned loss.");
    println!("  - SGD with momentum is simplest and most predictable for smooth convex cases.");
    println!("  - For physics fitting (α, D, K_u) with tiny parameters, scale lr accordingly.");

    Ok(())
}
