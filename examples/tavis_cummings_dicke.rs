//! Tavis-Cummings / Dicke Model: Collective Magnon-Photon Coupling
//!
//! **Difficulty**: ⭐⭐⭐ Advanced
//! **Category**: Cavity Magnonics
//! **Physics**: Tavis-Cummings model, Dicke superradiance, collective coupling √N enhancement,
//!              polariton splitting
//!
//! This example demonstrates the Tavis-Cummings (TC) model for N spin-1/2 emitters
//! (magnons) collectively coupled to a single microwave cavity mode. Unlike the
//! single-emitter Jaynes-Cummings model, the collective coupling constant
//!
//!   g_N = g · √N
//!
//! grows with the square root of the ensemble size, enabling macroscopic systems
//! such as YIG spheres (N ~ 10^17 spins) to reach the deep strong-coupling regime.
//!
//! ## What is demonstrated
//!
//! 1. **Collective coupling sweep** – g_N vs N for N ∈ {1, 4, 16, 64, 256, 1024}.
//!    Confirms the √N scaling and reports polariton frequencies, splitting, and cooperativity.
//!
//! 2. **Detuning sweep / anticrossing** – Fixing N = 256, sweeps ω_magnon across the
//!    cavity resonance (0.5 × ω_c → 1.5 × ω_c, 11 points), revealing the characteristic
//!    avoided crossing (level repulsion) of the polariton branches.  Hopfield photon/magnon
//!    fractions are computed analytically from the eigenvector of the 2×2 matrix.
//!
//! 3. **Superradiant phase-transition threshold** – Hepp-Lieb-Wang critical coupling
//!    g_c = √(ω_c · ω_m) / 2.  Checks whether a macroscopic YIG ensemble (N = 10^17)
//!    already exceeds g_c.
//!
//! 4. **Mean-field dynamics at resonance** – Euler integration of the coupled cavity-
//!    magnon equations of motion starting from (|α| = 1, |β| = 0).  The populations
//!    |α|² and |β|² oscillate at the vacuum Rabi frequency 2 g_N, analogous to the
//!    single-spin Rabi problem.
//!
//! ## References
//!
//! - M. Tavis, F. W. Cummings, Phys. Rev. **170**, 379 (1968)
//! - R. H. Dicke, Phys. Rev. **93**, 99 (1954)
//! - K. Hepp, E. H. Lieb, Ann. Phys. **76**, 360 (1973) — superradiant phase transition
//! - Y. Zhang et al., Phys. Rev. Lett. **113**, 156401 (2014) — magnon-polariton strong coupling
//! - H. Huebl et al., Phys. Rev. Lett. **111**, 127003 (2013) — YIG-cavity experiments

use std::f64::consts::PI;

use spintronics::math::Complex;
use spintronics::prelude::*;

// ---------------------------------------------------------------------------
// Helper: Hopfield coefficients for a two-mode coupled oscillator.
//
// For the 2×2 matrix
//
//   H = [[ω_c,  g_N],
//        [g_N,  ω_m]]
//
// the upper polariton eigenvector (cos θ, sin θ) gives photon fraction cos²θ
// and magnon fraction sin²θ, where
//
//   tan 2θ = 2 g_N / (ω_c - ω_m)   (Hopfield mixing angle)
//
// For the lower polariton the fractions swap.
// Returns (photon_fraction_upper, magnon_fraction_lower).
// ---------------------------------------------------------------------------
fn hopfield_fractions(omega_c: f64, omega_m: f64, g_n: f64) -> (f64, f64) {
    let delta = omega_c - omega_m;
    // Mixing angle θ.
    let theta = 0.5_f64 * (2.0 * g_n).atan2(delta);
    let cos2 = theta.cos().powi(2);
    let sin2 = theta.sin().powi(2);
    // Upper polariton: photon fraction = cos²θ (when ω_c > ω_m)
    // Lower polariton: magnon fraction = cos²θ
    // Return (photon fraction of upper branch, magnon fraction of lower branch).
    (cos2, sin2)
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Tavis-Cummings / Dicke Model: Collective Magnon-Photon Coupling ===\n");

    // -----------------------------------------------------------------------
    // 1. Collective coupling sweep: N ∈ {1, 4, 16, 64, 256, 1024}
    // -----------------------------------------------------------------------
    println!("=== 1. Collective Coupling √N Scaling ===");
    println!(
        "  {:>6}  {:>14}  {:>12}  {:>12}  {:>14}  {:>14}",
        "N", "g√N [MHz]", "ω+ [GHz]", "ω- [GHz]", "Splitting[MHz]", "Cooperativity"
    );
    println!("  {}", "-".repeat(82));

    let n_values: &[usize] = &[1, 4, 16, 64, 256, 1024];
    let two_pi = 2.0 * PI;
    let mhz = 1.0e6_f64;
    let ghz = 1.0e9_f64;

    for &n in n_values {
        let tc = TavisCummings::yig_ensemble(n);
        let g_n = tc.collective_coupling();
        let (omega_lo, omega_hi) = tc.polariton_frequencies();
        let splitting = omega_hi - omega_lo;
        let cooperativity = tc.cooperativity();

        println!(
            "  {:>6}  {:>14.3}  {:>12.5}  {:>12.5}  {:>14.3}  {:>14.3e}",
            n,
            g_n / (two_pi * mhz),
            omega_hi / (two_pi * ghz),
            omega_lo / (two_pi * ghz),
            splitting / (two_pi * mhz),
            cooperativity,
        );
    }

    // Verify √N ratio: g_N(4N)/g_N(N) = 2
    let tc_base = TavisCummings::yig_ensemble(64);
    let tc_quad = TavisCummings::yig_ensemble(256);
    let ratio = tc_quad.collective_coupling() / tc_base.collective_coupling();
    println!(
        "\n  √N scaling check: g_N(256)/g_N(64) = {:.6}  (expected 2.000000)",
        ratio
    );

    // -----------------------------------------------------------------------
    // 2. Detuning sweep — anticrossing (N = 256)
    // -----------------------------------------------------------------------
    println!("\n=== 2. Detuning Sweep — Anticrossing (N = 256) ===");
    println!(
        "  {:>8}  {:>10}  {:>10}  {:>14}  {:>18}  {:>18}",
        "Δ/ω_c",
        "ω+ [GHz]",
        "ω- [GHz]",
        "Split [MHz]",
        "Photon frac (upper)",
        "Magnon frac (lower)"
    );
    println!("  {}", "-".repeat(90));

    let n_sweep: usize = 256;
    // Build an ensemble at resonance to read off ω_c and g (per-spin).
    let tc_ref = TavisCummings::yig_ensemble(n_sweep);
    let omega_c_ref = tc_ref.omega_cavity;
    let g_per_spin = tc_ref.coupling_g;
    let kappa_ref = tc_ref.kappa;
    let gamma_ref = tc_ref.gamma_damp;
    let g_n_ref = tc_ref.collective_coupling();

    // Sweep ω_m from 0.5 × ω_c to 1.5 × ω_c in 11 steps.
    let n_det_steps = 11_usize;
    for step in 0..n_det_steps {
        let frac = 0.5 + (step as f64) * 1.0 / ((n_det_steps - 1) as f64);
        let omega_m_sweep = frac * omega_c_ref;
        let tc_sweep = TavisCummings::new(
            n_sweep,
            omega_c_ref,
            omega_m_sweep,
            g_per_spin,
            kappa_ref,
            gamma_ref,
        )?;

        let (omega_lo, omega_hi) = tc_sweep.polariton_frequencies();
        let splitting = omega_hi - omega_lo;
        let detuning_norm = (omega_c_ref - omega_m_sweep) / omega_c_ref;

        let (ph_upper, mag_lower) = hopfield_fractions(omega_c_ref, omega_m_sweep, g_n_ref);

        println!(
            "  {:>8.4}  {:>10.5}  {:>10.5}  {:>14.3}  {:>18.4}  {:>18.4}",
            detuning_norm,
            omega_hi / (two_pi * ghz),
            omega_lo / (two_pi * ghz),
            splitting / (two_pi * mhz),
            ph_upper,
            mag_lower,
        );
    }

    println!(
        "\n  Note: At resonance (Δ/ω_c = 0) splitting = 2g√N = {:.3} MHz",
        2.0 * g_n_ref / (two_pi * mhz)
    );

    // -----------------------------------------------------------------------
    // 3. Superradiant phase-transition threshold
    // -----------------------------------------------------------------------
    println!("\n=== 3. Superradiant Phase-Transition Threshold (Hepp-Lieb-Wang) ===");

    // Reference system at resonance (N=256 YIG parameters)
    let tc_sr = TavisCummings::yig_ensemble(256);
    let g_c = tc_sr.superradiant_threshold();
    let g_n_sr = tc_sr.collective_coupling();

    println!(
        "  ω_c = {:.4} GHz,  ω_m = {:.4} GHz",
        tc_sr.omega_cavity / (two_pi * ghz),
        tc_sr.omega_magnon / (two_pi * ghz),
    );
    println!("  g_c = √(ω_c·ω_m)/2 = {:.6} GHz", g_c / (two_pi * ghz));
    println!("  g_N (N=256) = {:.3} MHz", g_n_sr / (two_pi * mhz));
    println!("  g_N / g_c = {:.4e}   (superradiant if > 1)", g_n_sr / g_c);

    // Check a macroscopic YIG ensemble: N ~ 10^17 (1 mm sphere, ~4×10^17 spins)
    let n_yig_macro: usize = 100_000_000_000_000_000; // 10^17
    let tc_macro = TavisCummings::yig_ensemble(n_yig_macro);
    let g_n_macro = tc_macro.collective_coupling();
    let g_c_macro = tc_macro.superradiant_threshold();
    let above_threshold = g_n_macro > g_c_macro;

    println!("\n  Macroscopic YIG ensemble (N = 10^17):");
    println!(
        "    g_N = {:.4} GHz  (collective coupling)",
        g_n_macro / (two_pi * ghz)
    );
    println!(
        "    g_c = {:.4} GHz  (threshold)",
        g_c_macro / (two_pi * ghz)
    );
    println!("    g_N / g_c = {:.4e}", g_n_macro / g_c_macro);
    println!(
        "    Above Dicke superradiant threshold? {}",
        if above_threshold {
            "YES — deep superradiant regime"
        } else {
            "NO — normal phase"
        }
    );

    // -----------------------------------------------------------------------
    // 4. Mean-field dynamics at resonance (N = 64)
    // -----------------------------------------------------------------------
    println!("\n=== 4. Mean-Field Dynamics at Resonance (N = 64) ===");
    println!("  Equations of motion (lossless, rotating frame):");
    println!("    dα/dt = −iω_c α − i g_N β");
    println!("    dβ/dt = −iω_m β − i g_N α");
    println!("  Initial state: α(0) = 1.0 (all energy in cavity), β(0) = 0");
    println!("  Integration: explicit Euler, dt = 1/(100·ω_c), 200 steps\n");

    let n_dyn: usize = 64;
    let mut tc_dyn = TavisCummings::yig_ensemble(n_dyn);

    // Set ω_m = ω_c for exact resonance (modify via rebuild at resonance).
    tc_dyn = TavisCummings::new(
        n_dyn,
        tc_dyn.omega_cavity,
        tc_dyn.omega_cavity, // resonance: ω_m = ω_c
        tc_dyn.coupling_g,
        tc_dyn.kappa,
        tc_dyn.gamma_damp,
    )?;

    // Set initial conditions: cavity amplitude = 1+0i, magnon amplitude = 0.
    // We use mean_field_evolve with zero drive; seed cavity via first step pre-load.
    // Because the fields are private, we use a trick: seed by one drive step then
    // reset drive to zero. Instead, we expose initial state by relying on the fact
    // that mean_field_evolve accumulates da from the initial zero state with a drive.
    // The cleanest approach: drive with a delta-function kick at t=0.

    let omega_c_dyn = tc_dyn.omega_cavity;
    // Time step: dt = 1 / (100 × ω_c)
    let dt = 1.0 / (100.0 * omega_c_dyn);
    let n_steps: usize = 200;
    // Print interval: every 20 steps → 10 time points (including t=0 via initial kick)
    let print_every: usize = 20;

    // We inject a strong single-step drive at t=0 to load α(0) ≈ 1.
    // drive magnitude chosen so α grows by ≈ 1 in one dt:
    //   α(dt) ≈ α(0) + da(0)·dt ≈ 0 + drive·dt ≈ 1  →  drive = 1/dt
    let drive_kick = Complex::new(1.0 / dt, 0.0);

    // Step 0 (kick): load α
    tc_dyn.mean_field_evolve(dt, drive_kick);

    println!(
        "  {:>8}  {:>8}  {:>14}  {:>14}",
        "step", "t [ps]", "|α|² (photon)", "|β|² (magnon)"
    );
    println!("  {}", "-".repeat(52));

    let t_ps_scale = 1.0e12_f64; // seconds → picoseconds

    // Print state just after kick (step 0 printed as t=dt).
    let alpha0 = tc_dyn.cavity_amplitude();
    let beta0 = tc_dyn.magnon_amplitude();
    println!(
        "  {:>8}  {:>8.3}  {:>14.6}  {:>14.6}",
        0,
        dt * t_ps_scale,
        alpha0.norm_sq(),
        beta0.norm_sq(),
    );

    // Evolve 200 steps with zero drive; print every 20 steps.
    for step in 1..=n_steps {
        tc_dyn.mean_field_evolve(dt, Complex::ZERO);

        if step % print_every == 0 {
            let alpha = tc_dyn.cavity_amplitude();
            let beta = tc_dyn.magnon_amplitude();
            let t_total = (step as f64 + 1.0) * dt; // +1 for the kick step
            println!(
                "  {:>8}  {:>8.3}  {:>14.6}  {:>14.6}",
                step,
                t_total * t_ps_scale,
                alpha.norm_sq(),
                beta.norm_sq(),
            );
        }
    }

    let g_n_dyn = tc_dyn.collective_coupling();
    let rabi_period_ps = PI / g_n_dyn * t_ps_scale; // T_Rabi = π / g_N
    println!(
        "\n  g_N (N=64) = {:.3} MHz   →   Rabi half-period π/g_N = {:.3} ps",
        g_n_dyn / (two_pi * mhz),
        rabi_period_ps,
    );
    println!("  (Energy oscillates between cavity and magnon modes at frequency 2g_N)");

    // -----------------------------------------------------------------------
    // Summary
    // -----------------------------------------------------------------------
    println!("\n=== Summary ===");
    println!("Tavis-Cummings / Dicke model highlights:");
    println!("  - Collective coupling g_N = g·√N: doubling N increases coupling by √2 ≈ 1.414");
    println!("  - At resonance (Δ=0) polariton splitting = 2g_N (vacuum Rabi splitting)");
    println!(
        "  - Hopfield fractions: upper (lower) polariton is photon-like (magnon-like) when ω_c > ω_m"
    );
    println!("  - YIG macroscopic ensemble (N=10^17) exceeds superradiant threshold");
    println!("    (g_N ~ 10^4 THz >> g_c ~ 2.65 GHz for YIG with g=2π×100 MHz per spin)");
    println!(
        "  - Mean-field dynamics show coherent Rabi oscillations between cavity and magnon modes"
    );

    Ok(())
}
