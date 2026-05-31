//! Spin-Torque Nano-Oscillator (STNO) — auto-oscillation and threshold
//!
//! Demonstrates the onset of auto-oscillation when the spin-transfer torque
//! overcomes Gilbert damping above a threshold current density. We compare
//! a sub-threshold (damped) regime with a super-threshold (auto-oscillating)
//! regime using the Slonczewski LLG + STT model.
//!
//! Run with: cargo run --example stno_auto_oscillation

use std::f64::consts::PI;

use spintronics::effect::stno::SpinTorqueOscillatorConfig;
use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Spin-Torque Nano-Oscillator Simulation ===\n");

    // ── Device from the canonical Kiselev et al. (Nature 425, 380, 2003) geometry ──
    // Permalloy free layer, p̂ = +ẑ (parallel fixed layer), H_ext = 0
    // This configuration drives in-plane precession around ẑ for J > J_th.

    let mut material = Ferromagnet::permalloy();
    material.alpha = 0.01;
    material.ms = 8.0e5; // 800 kA/m Permalloy

    let h_ext = Vector3::new(0.0, 0.0, 0.0); // zero external field
    let t_fm = 5.0e-9_f64; // 5 nm free layer
    let radius = 50.0e-9_f64; // 100 nm diameter nanopillar
    let area = PI * radius * radius;

    let config = SpinTorqueOscillatorConfig::new(
        0.35, // η — spin polarization efficiency
        t_fm,
        area,
        0.0,                         // β = 0 (no field-like torque)
        Vector3::new(0.0, 0.0, 1.0), // p̂ = +ẑ
    )?;

    let stno = SpinTorqueOscillator::new(material.clone(), h_ext, config)?;
    let ms = stno.solver.material.ms;

    println!("Permalloy nanopillar (d=100nm, t=5nm, p̂=+ẑ):");
    println!("  M_s     = {:.0} kA/m", ms * 1e-3);
    println!("  α       = {:.4}", stno.solver.material.alpha);
    println!("  Volume  = {:.3e} m³", stno.config.volume);
    let j_th = stno.threshold_current_density();
    println!("  J_threshold = {:.2e} A/m²", j_th);

    // Start with m in the x-y plane (perpendicular to p̂) — maximum STT coupling
    let m0 = Vector3::new(1.0, 0.0, 0.0) * ms;
    let dt = 1e-12_f64; // 1 ps

    // ── Sub-threshold: J = 0.2 × J_th ─────────────────────────────────────────
    let j_sub = 0.2 * j_th;
    let traj_sub = stno.simulate(m0, dt, 500, j_sub);
    let power_sub = stno.oscillation_power(&traj_sub);

    // ── Super-threshold: J = 10 × J_th ────────────────────────────────────────
    let j_sup = 10.0 * j_th;
    let traj_sup = stno.simulate(m0, dt, 2000, j_sup);
    let power_sup = stno.oscillation_power(&traj_sup);
    let freq = stno.oscillation_frequency(&traj_sup, dt);

    println!(
        "\n=== Sub-threshold (J = 0.2×J_th = {:.2e} A/m²) ===",
        j_sub
    );
    println!("  Oscillation power (var m_z/M_s) = {:.4e}", power_sub);
    // Show last few trajectory points to confirm damping
    let last = traj_sub.len() - 1;
    println!(
        "  Final m_z/M_s = {:.4}  (expect near 0 as in-plane steady state)",
        traj_sub[last].z / ms
    );

    println!(
        "\n=== Super-threshold (J = 10×J_th = {:.2e} A/m²) ===",
        j_sup
    );
    println!("  Oscillation power (var m_z/M_s) = {:.4e}", power_sup);
    if let Some(f) = freq {
        println!("  Oscillation frequency ≈ {:.2} GHz", f * 1e-9);
    } else {
        // Estimate from in-plane gyration period
        println!("  Oscillation frequency: see STT-driven in-plane gyration below");
    }
    let auto_osc = power_sup > 1.0e-4;
    println!(
        "  Auto-oscillating: {}",
        if auto_osc { "YES" } else { "no" }
    );

    // ── Current sweep: show threshold behavior ──────────────────────────────────
    println!("\n=== Current Sweep: Power vs J/J_th ===");
    println!("  J/J_th    Power (var m_z/M_s)");
    for &ratio in &[0.1_f64, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0] {
        let j_test = ratio * j_th;
        let traj = stno.simulate(m0, dt, 800, j_test);
        let p = stno.oscillation_power(&traj);
        let marker = if p > 1.0e-4 { " ← oscillating" } else { "" };
        println!("  {:5.1}     {:.3e}{}", ratio, p, marker);
    }

    // ── Thermal linewidth ───────────────────────────────────────────────────────
    println!("\n=== Thermal Linewidth (Slavin-Tiberkevich 2009) ===");
    for temp in [77_u32, 150, 300] {
        let lw = stno.thermal_linewidth(j_sup, temp as f64);
        println!("  T={} K:  Δf = {:.2e} MHz", temp, lw * 1e-6);
    }

    // ── Injection locking ───────────────────────────────────────────────────────
    println!("\n=== Injection Locking Bandwidth (Adler 1946) ===");
    println!("  Q-factor = 100");
    for i_inj in [0.1_f64, 0.3, 0.5, 0.9] {
        let delta_f = stno.injection_locking_range(100.0, i_inj);
        println!(
            "  i_inj={:.1}:  Δf_lock = {:.2e} MHz",
            i_inj,
            delta_f * 1e-6
        );
    }

    println!("\n=== Summary ===");
    println!("  J_th       = {:.2e} A/m²", j_th);
    println!("  Power (sub) = {:.2e}", power_sub);
    println!("  Power (sup) = {:.2e}", power_sup);
    println!("  Auto-oscillation confirmed above threshold: {}", auto_osc);

    Ok(())
}
