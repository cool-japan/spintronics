//! Reinforcement Learning for SOT Switching Protocol Optimization
//!
//! Uses the Cross-Entropy Method (CEM) to find optimal current pulse sequences
//! for switching perpendicular magnetization (PMA) in a CoFeB/Pt device.
//! The agent maximises a reward that balances switching success with low
//! energy dissipation and short switching time.
//!
//! Run with: cargo run --example rl_sot_switching

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== RL-Optimized SOT Switching Protocol ===\n");

    let env = SotSwitchingEnv::default_cofeb_pt();
    println!("CoFeB/Pt device (PMA):");
    println!("  J_max           = {:.1e} A/m²", env.config.j_max);
    println!("  t_FM            = {:.1} nm", env.config.t_fm * 1e9);
    println!("  H_bias (in-plane) = {:.1} kA/m", env.config.h_bias * 1e-3);
    println!("  n_pulse_steps   = {}", env.config.n_pulse_steps);
    println!("  Max episode steps = {}", env.config.max_steps);
    println!(
        "  Switch threshold (m_z) = {:.2}",
        env.config.switch_threshold
    );

    let mut optimizer = SotRlOptimizer::new(env);
    println!("\nCEM hyperparameters:");
    println!("  Population size = {}", optimizer.population_size);
    println!(
        "  Elite fraction  = {:.0}%",
        optimizer.elite_fraction * 100.0
    );

    // ── Demonstrate single-step SOT dynamics before CEM training ──────────────
    {
        let mut demo_env = SotSwitchingEnv::default_cofeb_pt();
        let m0 = demo_env.reset();
        println!(
            "\nSingle-step dynamics demo from m = ({:.3}, {:.3}, {:.3}):",
            m0.x, m0.y, m0.z
        );
        let j_max = demo_env.config.j_max;
        for step in 0..5 {
            let j = if step % 2 == 0 { j_max } else { -j_max };
            let (m, reward, done) = demo_env.step(j);
            println!(
                "  step {}: j={:+.1e}  m_z={:+.4}  reward={:+.3}{}",
                step + 1,
                j,
                m.z,
                reward,
                if done { " (done)" } else { "" }
            );
            if done {
                break;
            }
        }
    }

    println!("\nTraining CEM agent (20 generations)...");
    let result = optimizer.train(20, 42);

    println!("\nTraining complete:");
    println!("  Best reward:         {:.3}", result.best_reward);
    println!("  Switching achieved:  {}", result.switching_achieved);
    println!(
        "  Initial mean reward (gen 1):  {:.4}",
        result.mean_rewards_per_gen.first().copied().unwrap_or(0.0)
    );
    println!(
        "  Final mean reward  (gen 20):  {:.4}",
        result.mean_rewards_per_gen.last().copied().unwrap_or(0.0)
    );
    // Improvement metric
    let first_r = result.mean_rewards_per_gen.first().copied().unwrap_or(0.0);
    let last_r = result.mean_rewards_per_gen.last().copied().unwrap_or(0.0);
    println!(
        "  Reward improvement: {:.4} → {:.4}  (Δ = {:.4})",
        first_r,
        last_r,
        last_r - first_r
    );

    // Learning curve
    println!("\nLearning curve (mean reward per generation):");
    for (i, &r) in result.mean_rewards_per_gen.iter().enumerate() {
        // Normalise to [0,1] for bar: reward is in [-3, 10]
        let bar_len = (((r + 0.1) / 0.1).max(0.0) * 5.0) as usize;
        let bar: String = "#".repeat(bar_len.min(50));
        println!("  Gen {:2}: {:8.4} |{}", i + 1, r, bar);
    }

    println!(
        "\nBest pulse sequence (J amplitudes, n={}):",
        result.best_policy.len()
    );
    let j_max = optimizer.env.config.j_max;
    for (i, j) in result.best_policy.iter().enumerate() {
        let pct = (j.abs() / j_max) * 100.0;
        let dir = if *j >= 0.0 { "+" } else { "-" };
        println!(
            "  step {:2}: J = {:+.2e} A/m²  ({}{:.0}% of |J_max|)",
            i + 1,
            j,
            dir,
            pct
        );
    }

    // Verification run with best policy
    let mut verify_env = SotSwitchingEnv::default_cofeb_pt();
    let m0 = verify_env.reset();
    println!(
        "\nVerification run from m = ({:.4}, {:.4}, {:.4}):",
        m0.x, m0.y, m0.z
    );
    println!("  (showing all steps up to 10, then final state)");

    let mut m = m0;
    let mut total_r = 0.0_f64;
    let n_show = result.best_policy.len().min(10);
    for (step, &j_amp) in result.best_policy.iter().enumerate() {
        let (new_m, reward, done) = verify_env.step(j_amp);
        m = new_m;
        total_r += reward;
        if step < n_show || done {
            println!(
                "  step {:2}: m_z={:+.4}  reward={:+.4}{}",
                step + 1,
                m.z,
                reward,
                if done { " ← episode done" } else { "" }
            );
        }
        if done {
            break;
        }
    }

    println!("\n=== Summary ===");
    println!("  Final m_z      = {:.4}", m.z);
    println!(
        "  Target m_z     < {:.2}",
        verify_env.config.switch_threshold
    );
    println!(
        "  Switched:        {}",
        m.z < verify_env.config.switch_threshold
    );
    println!("  Total reward   = {:.4}", total_r);
    println!("  Best gen reward = {:.4}", result.best_reward);

    Ok(())
}
