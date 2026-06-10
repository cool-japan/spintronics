//! GPU Device Abstraction Demo
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Numerical Methods / Heterogeneous Computing
//! **Physics**: Device-abstracted LLG evolution for portable performance
//!
//! ## Background
//!
//! The v0.9.0 GPU skeleton introduces a `Device` trait that lets simulation
//! code run on CPU or GPU transparently. v0.9.0 ships:
//!   - `CpuDevice` — always available, drop-in baseline (wraps existing CPU path)
//!   - `CudaDevice` (feature = "cuda") — stub for v1.0.0 CUDA kernels
//!
//! This example demonstrates the device-agnostic API on a batch of spins:
//!   1. Enumerate available devices
//!   2. Use `select_best_device()` to pick the best
//!   3. Evolve a batch of N=100 spins under a uniform field for 100 RK4 steps
//!   4. Compare wall time vs. spin count to show scaling characteristics
//!
//! When v1.0.0 lands a real CUDA backend, the same example code will pick
//! `CudaDevice` automatically (if available) and run on GPU.
//!
//! ## References
//! - mumax3 — Vansteenkiste et al., AIP Adv. 4, 107133 (2014)

use std::time::Instant;

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  GPU Device Abstraction Demo");
    println!("=============================================================");

    // -------------------------------------------------------------------------
    // Section 1: Enumerate available devices
    // -------------------------------------------------------------------------
    println!("\n--- Section 1: Available Devices ---\n");

    let devices = available_devices();
    println!(
        "  {:>5}  {:>10}  {:>12}  {:>10}",
        "idx", "name", "available?", "max_spins"
    );
    println!("  {}", "-".repeat(42));
    for (i, dev) in devices.iter().enumerate() {
        let avail = if dev.is_available() { "yes" } else { "no" };
        let max = dev.max_spins();
        let max_str = if max == usize::MAX {
            "∞".to_string()
        } else {
            format!("{max}")
        };
        println!(
            "  {:>5}  {:>10}  {:>12}  {:>10}",
            i,
            dev.name(),
            avail,
            max_str
        );
    }

    // -------------------------------------------------------------------------
    // Section 2: Select best device
    // -------------------------------------------------------------------------
    println!("\n--- Section 2: Selected Device ---\n");

    let device = select_best_device();
    println!("  select_best_device() → '{}'", device.name());
    println!("  is_available()       → {}", device.is_available());

    // -------------------------------------------------------------------------
    // Section 3: Evolve a batch of spins
    // -------------------------------------------------------------------------
    println!("\n--- Section 3: Batch Evolution (100 spins, 200 steps) ---\n");

    let n_spins = 100_usize;
    let n_steps = 200_usize;
    let alpha = 0.05_f64;
    let dt = 1.0e-13_f64;

    // Build batch: spins tilted slightly from +ẑ
    let mut spins: Vec<[f64; 3]> = (0..n_spins)
        .map(|i| {
            let theta = 0.3 + 0.005 * (i as f64);
            let phi = 0.1 * (i as f64);
            [
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            ]
        })
        .collect();
    let h_eff: Vec<[f64; 3]> = vec![[0.0, 0.0, 1.0e6]; n_spins];

    let start = Instant::now();
    device.step_llg_rk4_multi(&mut spins, &h_eff, alpha, dt, n_steps)?;
    let elapsed = start.elapsed();

    let mz_avg: f64 = spins.iter().map(|s| s[2]).sum::<f64>() / n_spins as f64;
    let mag_drift_max: f64 = spins
        .iter()
        .map(|s| ((s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt() - 1.0).abs())
        .fold(0.0_f64, f64::max);

    println!(
        "  Wall time:                {:.3} ms",
        elapsed.as_secs_f64() * 1000.0
    );
    println!("  Final ⟨m_z⟩ (post evolve): {mz_avg:.6}");
    println!("  Max |m| drift from 1.0:   {mag_drift_max:.2e}");
    println!(
        "  Time per spin per step:   {:.3} ns",
        elapsed.as_secs_f64() * 1e9 / (n_spins * n_steps) as f64
    );

    // -------------------------------------------------------------------------
    // Section 4: Zeeman energy via device
    // -------------------------------------------------------------------------
    println!("\n--- Section 4: Zeeman Energy via Device ---\n");

    let ms = 1.4e5_f64;
    let energy = device.zeeman_energy(&spins, &h_eff, ms)?;
    println!("  E_Zeeman = {energy:.4e} J·m³ (Σ_i -μ₀ ms m_i·h_i)");

    // -------------------------------------------------------------------------
    // Section 5: Scaling — compare 10, 100, 1000 spins
    // -------------------------------------------------------------------------
    println!("\n--- Section 5: Scaling Characteristics ---\n");
    println!(
        "  {:>10}  {:>14}  {:>14}",
        "n_spins", "time_total (ms)", "ns / (spin·step)"
    );
    println!("  {}", "-".repeat(42));
    for &n in [10_usize, 100, 500, 1000].iter() {
        let mut batch: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let theta = 0.3 + 0.005 * (i as f64);
                [theta.sin(), 0.0, theta.cos()]
            })
            .collect();
        let h_batch: Vec<[f64; 3]> = vec![[0.0, 0.0, 1.0e6]; n];
        let n_steps_local = 100_usize;
        let t0 = Instant::now();
        device.step_llg_rk4_multi(&mut batch, &h_batch, alpha, dt, n_steps_local)?;
        let dt_total = t0.elapsed().as_secs_f64() * 1000.0;
        let per = dt_total * 1e6 / (n * n_steps_local) as f64;
        println!("  {:>10}  {:>14.3}  {:>14.3}", n, dt_total, per);
    }

    println!("\n=============================================================");
    println!("  Done. The Device trait gives a stable API for v1.0.0 CUDA");
    println!("  kernels; CpuDevice provides the unconditional baseline today.");
    println!("=============================================================\n");

    Ok(())
}
