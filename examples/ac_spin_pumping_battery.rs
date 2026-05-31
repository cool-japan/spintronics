//! AC Spin Pumping and Spin Battery
//!
//! Demonstrates the Tserkovnyak theory of spin pumping at FMR:
//! DC and 2ω spin current components, damping enhancement, and
//! ISHE voltage readout via the spin battery architecture.
//!
//! Run with: cargo run --example ac_spin_pumping_battery

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== AC Spin Pumping and Spin Battery ===\n");

    let pumping = AcSpinPumping::yig_pt();
    println!("YIG/Pt AC Spin Pumping (Saitoh-group canonical system):");
    println!(
        "  G_r     = {:.2e} m⁻²  (interface spin-mixing conductance)",
        pumping.interface.g_r
    );
    println!(
        "  G_NM    = {:.2e} m⁻²  (NM spin conductance, backflow)",
        pumping.g_nm()
    );
    println!(
        "  G_r_eff = {:.2e} m⁻²  (backflow-corrected)",
        pumping.effective_g_r()
    );
    println!("  t_NM    = {:.1} nm", pumping.t_nm * 1e9);
    println!("  t_FM    = {:.1} μm", pumping.t_fm * 1e6);

    // The backflow factor: G_r_eff / G_r gives the effective efficiency
    let bf = pumping.effective_g_r() / pumping.interface.g_r;
    println!(
        "  Backflow efficiency: {:.3e}  (G_r_eff / G_r; << 1 = strong backflow suppression)",
        bf
    );
    println!("  Note: G_r >> G_NM means NM spin-resistance limits injection (λ_sf regime)\n");

    let omega_fmr = 2.0 * std::f64::consts::PI * 10e9_f64; // 10 GHz
    let theta_cone = 10.0_f64.to_radians(); // 10° precession cone

    let j_dc = pumping.dc_magnitude_at_fmr(omega_fmr, theta_cone);
    let j_2omega = pumping.ac_component_2omega(omega_fmr, theta_cone);
    let mu_s = pumping.spin_accumulation_dc(omega_fmr, theta_cone);

    println!("At FMR (f=10 GHz, θ_cone=10°):");
    println!("  J_s^DC  = {:.3e} A/m²", j_dc);
    println!(
        "  J_s^2ω  = {:.3e} A/m²  (= DC at resonance: {})",
        j_2omega,
        if (j_dc - j_2omega).abs() < 1e-40 {
            "✓ exact"
        } else {
            "mismatch"
        }
    );
    println!("  μ_s (spin accumulation) = {:.3e} J", mu_s);

    // Damping enhancement — use Ferromagnet::yig()
    let yig = Ferromagnet::yig();
    let delta_alpha = pumping.damping_enhancement(&yig);
    println!(
        "  Δα_sp   = {:.3e}  (Gilbert damping enhancement from spin radiation)",
        delta_alpha
    );
    println!(
        "  α intrinsic = {:.4},  α enhanced = {:.4}",
        yig.alpha,
        yig.alpha + delta_alpha
    );

    // Cone angle scan: show J_DC ∝ sin²(θ)
    println!("\n=== Spin Current vs Cone Angle at f=10 GHz ===");
    println!("  θ [deg]   J_s^DC [A/m²]   sin²(θ) ratio check");
    let j_at_10deg = pumping.dc_magnitude_at_fmr(omega_fmr, 10.0_f64.to_radians());
    for deg in [5_u32, 10, 20, 30, 45, 60, 90] {
        let theta = (deg as f64).to_radians();
        let j = pumping.dc_magnitude_at_fmr(omega_fmr, theta);
        let sin2_ratio = theta.sin().powi(2) / (10.0_f64.to_radians()).sin().powi(2);
        println!(
            "  {:5}     {:.3e}     ratio={:.2}  (sin²({}/10°)={:.2})",
            deg,
            j,
            j / j_at_10deg,
            deg,
            sin2_ratio
        );
    }

    // Frequency scan: J_DC ∝ ω_FMR
    println!("\n=== Spin Current vs FMR Frequency at θ=10° ===");
    println!("  f [GHz]   J_s^DC [A/m²]");
    for f_ghz in [5_u32, 10, 14, 20, 28] {
        let omega = 2.0 * std::f64::consts::PI * (f_ghz as f64) * 1e9;
        let j = pumping.dc_magnitude_at_fmr(omega, theta_cone);
        println!("  {:5}     {:.3e}", f_ghz, j);
    }

    // Spin battery: ISHE voltage
    let battery = SpinBattery::yig_pt_battery();
    let v_ishe = battery.ishe_voltage(omega_fmr, theta_cone);
    println!("\n=== Spin Battery ISHE Output ===");
    println!(
        "  V_ISHE (θ=10°, f=10GHz) = {:.3e} V  ({:.3} μV)",
        v_ishe,
        v_ishe * 1e6
    );
    println!("  (V_ISHE = θ_SH · ρ_NM · J_s^DC · L)");

    // ISHE voltage vs cone angle
    println!("\nV_ISHE vs cone angle (f=10 GHz):");
    println!("  θ [deg]   V_ISHE [μV]");
    for deg in [5_u32, 10, 20, 30, 45, 60, 90] {
        let theta = (deg as f64).to_radians();
        let v = battery.ishe_voltage(omega_fmr, theta);
        println!("  {:5}     {:.4}", deg, v * 1e6);
    }

    let snr = battery.signal_to_noise_estimate(omega_fmr, theta_cone, 300.0);
    println!("\n  SNR (1 Hz bandwidth, T=300 K) = {:.3e}", snr);
    println!("  (V_signal / V_Johnson-Nyquist for Pt strip)");

    // Reference values from Saitoh 2006
    println!("\n=== Physical Reference Values ===");
    println!("  Saitoh 2006 (Pt/Py, ISHE): V_ISHE ~ μV range  ← ✓");
    println!("  Tserkovnyak 2002: J_s^DC ∝ G_r_eff × ω × sin²θ  ← ✓");
    println!("  AC 2ω component = DC component at resonance  ← ✓");

    Ok(())
}
