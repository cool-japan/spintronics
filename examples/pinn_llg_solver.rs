//! Physics-Informed Neural Network for LLG Dynamics
//!
//! **Difficulty**: ⭐⭐⭐⭐
//! **Category**: ML / Physics-informed methods
//! **Physics**: Solve Larmor precession via PINN, compare with RK4 ground truth
//!
//! ## Background
//!
//! Physics-informed neural networks (PINNs) approximate a solution u(t) by a
//! neural network u_θ(t), and *fit* the parameters θ so that the residual of
//! the governing PDE/ODE is small at collocation points. For LLG:
//!     ∂m/∂t = -γ/(1+α²)·[m × H + α·m × (m × H)]
//! the PINN loss is
//!     L = Σ_k |Residual(t_k)|² + λ_IC |m(0) - m_0|²
//!
//! The autodiff tape lets us compute ∂m_θ/∂t analytically (or, here, via
//! central finite differences on the same tape, which avoids nested tapes).
//!
//! This example trains a small MLP on Larmor precession (α=0.05, H || ẑ) and
//! compares predictions against an RK4 reference.
//!
//! ## References
//! - Raissi, Perdikaris, Karniadakis, JCP 378, 686 (2019) — PINN paper
//! - Karniadakis et al., Nat. Rev. Phys. 3, 422 (2021) — review
//! - Baydin et al., JMLR 18, 1 (2018) — autodiff for ML

use spintronics::prelude::*;

fn rk4_step(m: [f64; 3], h_eff: [f64; 3], alpha: f64, gamma: f64, dt: f64) -> [f64; 3] {
    let f = |s: [f64; 3]| {
        let mxh = [
            s[1] * h_eff[2] - s[2] * h_eff[1],
            s[2] * h_eff[0] - s[0] * h_eff[2],
            s[0] * h_eff[1] - s[1] * h_eff[0],
        ];
        let m_mxh = [
            s[1] * mxh[2] - s[2] * mxh[1],
            s[2] * mxh[0] - s[0] * mxh[2],
            s[0] * mxh[1] - s[1] * mxh[0],
        ];
        let coeff = -gamma / (1.0 + alpha * alpha);
        [
            coeff * (mxh[0] + alpha * m_mxh[0]),
            coeff * (mxh[1] + alpha * m_mxh[1]),
            coeff * (mxh[2] + alpha * m_mxh[2]),
        ]
    };
    let k1 = f(m);
    let m2 = [
        m[0] + 0.5 * dt * k1[0],
        m[1] + 0.5 * dt * k1[1],
        m[2] + 0.5 * dt * k1[2],
    ];
    let k2 = f(m2);
    let m3 = [
        m[0] + 0.5 * dt * k2[0],
        m[1] + 0.5 * dt * k2[1],
        m[2] + 0.5 * dt * k2[2],
    ];
    let k3 = f(m3);
    let m4 = [m[0] + dt * k3[0], m[1] + dt * k3[1], m[2] + dt * k3[2]];
    let k4 = f(m4);
    [
        m[0] + (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        m[1] + (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        m[2] + (dt / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
    ]
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  Physics-Informed Neural Network for LLG Larmor Precession");
    println!("=============================================================");

    // -------------------------------------------------------------------------
    // Section 1: Problem setup — Larmor precession in scaled time
    // -------------------------------------------------------------------------
    println!("\n--- Section 1: Setup ---\n");

    let alpha = 0.05_f64;
    // Use SCALED units: γ·μ_0·H = 1 (dimensionless), so one Larmor period = 2π
    let gamma = 1.0_f64;
    let h_eff = [0.0_f64, 0.0, 1.0]; // ẑ
    let initial_m = [1.0_f64, 0.0, 0.0]; // start in x̂
    let t_end = 2.0 * std::f64::consts::PI; // one Larmor period

    println!("  α            = {alpha}");
    println!("  γ            = {gamma} (scaled units)");
    println!("  H_eff        = {h_eff:?}");
    println!("  m₀           = {initial_m:?}");
    println!("  t_end        = {t_end:.4} (one Larmor period)");

    // -------------------------------------------------------------------------
    // Section 2: Build and train PINN
    // -------------------------------------------------------------------------
    println!("\n--- Section 2: Build & Train PINN ---\n");

    // Hidden layers: [16, 16]; mapping t → (mx, my, mz)
    let mut pinn = LlgPinn::new(&[16, 16], h_eff, alpha, gamma, 1234)?;
    println!(
        "  PINN MLP: 1 → 16 → 16 → 3  ({} parameters)",
        pinn.n_params()
    );

    // Collocation points uniformly distributed in [0, t_end]
    let n_coll = 32;
    let t_coll: Vec<f64> = (0..n_coll)
        .map(|i| t_end * (i as f64) / (n_coll as f64 - 1.0))
        .collect();
    println!("  Collocation points: {n_coll}");

    let trainer = PinnTrainer::new(t_coll.clone(), initial_m).with_weights(1.0, 1000.0, 10.0);

    let result = trainer.train(&mut pinn, 2000, 3e-3, OptimizerKind::Adam)?;
    println!("  Trained for {} iterations (Adam)", result.n_iterations);
    println!(
        "  Initial loss: {:.4e}",
        result.loss_history.first().copied().unwrap_or(0.0)
    );
    println!("  Final loss:   {:.4e}", result.final_loss);
    println!("  Converged:    {}", result.converged);

    // -------------------------------------------------------------------------
    // Section 3: Compare PINN prediction with RK4 ground truth
    // -------------------------------------------------------------------------
    println!("\n--- Section 3: PINN vs RK4 Reference ---\n");

    // Generate RK4 reference trajectory
    let n_ref = 16;
    let dt_ref = t_end / (n_ref as f64);
    let mut m_ref = initial_m;
    let mut ref_trajectory = vec![(0.0_f64, m_ref)];
    for _ in 0..n_ref {
        m_ref = rk4_step(m_ref, h_eff, alpha, gamma, dt_ref);
        let t_curr = ref_trajectory.last().unwrap().0 + dt_ref;
        ref_trajectory.push((t_curr, m_ref));
    }

    println!(
        "  {:>8}  {:>9}  {:>9}  {:>9}  {:>9}  {:>9}  {:>9}",
        "t", "mx_RK4", "my_RK4", "mz_RK4", "mx_NN", "my_NN", "mz_NN"
    );
    println!("  {}", "-".repeat(72));
    for &(t_curr, m_true) in ref_trajectory.iter() {
        let m_nn = pinn.predict(t_curr)?;
        println!(
            "  {:>8.3}  {:>+9.4}  {:>+9.4}  {:>+9.4}  {:>+9.4}  {:>+9.4}  {:>+9.4}",
            t_curr, m_true[0], m_true[1], m_true[2], m_nn[0], m_nn[1], m_nn[2]
        );
    }

    // L²-error on a denser grid
    let n_test = 33;
    let mut l2 = 0.0_f64;
    for i in 0..n_test {
        let t_curr = t_end * (i as f64) / (n_test as f64 - 1.0);
        // Get RK4 result by integrating from 0 with fine dt
        let mut m_int = initial_m;
        let n_fine = 200;
        let dt_fine = t_curr / (n_fine as f64);
        if t_curr > 0.0 {
            for _ in 0..n_fine {
                m_int = rk4_step(m_int, h_eff, alpha, gamma, dt_fine);
            }
        }
        let m_nn = pinn.predict(t_curr)?;
        let dx = m_nn[0] - m_int[0];
        let dy = m_nn[1] - m_int[1];
        let dz = m_nn[2] - m_int[2];
        l2 += dx * dx + dy * dy + dz * dz;
    }
    l2 = (l2 / (n_test as f64)).sqrt();
    println!("\n  L² error PINN vs RK4 (averaged over {n_test} pts): {l2:.4e}");

    // -------------------------------------------------------------------------
    // Section 4: Conservation diagnostics
    // -------------------------------------------------------------------------
    println!("\n--- Section 4: |m(t)| Conservation (PINN) ---\n");
    println!("  {:>8}  {:>14}", "t", "|m_NN|");
    println!("  {}", "-".repeat(24));
    for i in 0..=8 {
        let t_curr = t_end * (i as f64) / 8.0;
        let m_nn = pinn.predict(t_curr)?;
        let norm = (m_nn[0].powi(2) + m_nn[1].powi(2) + m_nn[2].powi(2)).sqrt();
        println!("  {:>8.3}  {:>14.6}", t_curr, norm);
    }

    println!("\n=============================================================");
    println!("  Done. PINN trained on Larmor precession. |m| conserved by norm");
    println!("  penalty; full convergence to the exact trajectory remains hard");
    println!("  (L² error ~{l2:.2e}) — typical of unconstrained PINNs with");
    println!("  small networks and few collocation points.");
    println!("=============================================================\n");

    Ok(())
}
