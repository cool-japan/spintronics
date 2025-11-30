//! Spin-Torque Nano-Oscillator (STNO)
//!
//! This example simulates a spin-transfer torque driven nano-oscillator:
//! - DC current through magnetic trilayer
//! - Spin-transfer torque (STT) excites magnetization precession
//! - Auto-oscillations at GHz frequencies
//! - Voltage output from GMR/TMR effect
//!
//! ## Physics:
//! When DC current exceeds critical threshold, the Slonczewski spin-transfer
//! torque overcomes damping, leading to sustained oscillations. The oscillation
//! frequency can be tuned by current and magnetic field.
//!
//! ## References:
//! - S. I. Kiselev et al., "Microwave oscillations of a nanomagnet driven by a
//!   spin-polarized current", Nature 425, 380 (2003)
//!
//! Run with: cargo run --example spin_torque_oscillator

use spintronics::prelude::*;
use std::f64::consts::PI;

fn main() {
    println!("=== Spin-Torque Nano-Oscillator Simulation ===\n");

    // === 1. Device Structure ===
    println!("=== Device Structure ===");

    let structure = MagneticMultilayer::gmr_py_cu_py();

    println!("GMR Trilayer: Py(fixed) / Cu / Py(free)");
    println!("  Fixed layer: {:.0} nm", structure.bottom_thickness);
    println!("  Spacer (Cu): {:.1} nm", structure.spacer.thickness);
    println!(
        "  Free layer: {:.0} nm (oscillating)\n",
        structure.top_thickness
    );

    let py = &structure.top_layer;

    // Device geometry
    let radius = 50e-9; // 50 nm radius nanopillar
    let area = PI * radius * radius;

    println!("Nanopillar Geometry:");
    println!("  Radius: {:.0} nm", radius * 1e9);
    println!("  Area: {:.0} nm²", area * 1e18);
    println!(
        "  Free layer volume: {:.2e} m³\n",
        area * structure.top_thickness * 1e-9
    );

    // === 2. Spin-Transfer Torque Parameters ===
    println!("=== Spin-Transfer Torque ===");

    let spin_polarization = 0.35; // Py spin polarization
    let g_factor = 0.5; // Slonczewski prefactor

    println!("  Spin polarization: {:.2}", spin_polarization);
    println!("  Slonczewski g-factor: {:.2}", g_factor);

    // Critical current for oscillations
    // I_c = (2e/ℏ) × α × (M_s V) × (H_eff / P)
    let e = 1.602e-19; // C
    let volume = area * structure.top_thickness * 1e-9;
    let h_eff = 0.05; // T effective field
    let i_critical = (2.0 * e / HBAR) * py.alpha * (py.ms * volume) * (h_eff / spin_polarization);

    println!("  Critical current I_c: {:.2} mA", i_critical * 1e3);
    println!(
        "  Current density j_c: {:.1} MA/cm²\n",
        (i_critical / area) * 1e-10
    );

    // === 3. Operating Point ===
    println!("=== Operating Conditions ===");

    let i_dc = 1.5 * i_critical; // 150% of threshold
    let h_ext = 0.1; // 1000 Oe external field
    let temperature = 300.0; // K

    println!(
        "  DC current: {:.2} mA ({:.0}% of I_c)",
        i_dc * 1e3,
        (i_dc / i_critical) * 100.0
    );
    println!("  External field: {:.0} Oe", h_ext * 10000.0);
    println!("  Temperature: {:.0} K\n", temperature);

    // === 4. Oscillation Characteristics ===
    println!("=== Auto-Oscillation Properties ===");

    // Frequency (FMR-like)
    let gamma_hz = GAMMA / (2.0 * PI); // Hz/T
    let f_osc = gamma_hz * h_eff;

    println!("  Oscillation frequency: {:.2} GHz", f_osc * 1e-9);

    // Linewidth (determined by damping and thermal noise)
    let delta_f = py.alpha * f_osc;
    let kb = 1.380649e-23; // J/K
    let thermal_linewidth = (kb * temperature) / (2.0 * PI * HBAR) * 1e-9; // GHz

    println!("  Intrinsic linewidth: {:.1} MHz", delta_f * 1e-6);
    println!("  Thermal broadening: {:.1} GHz", thermal_linewidth);

    // Quality factor
    let q_factor = f_osc / delta_f;
    println!("  Quality factor Q: {:.0}", q_factor);

    // Precession amplitude (cone angle)
    let overdrive = i_dc / i_critical - 1.0;
    let cone_angle = (2.0 * overdrive).sqrt().min(1.0) * 30.0 * PI / 180.0; // Saturates ~30°

    println!("  Precession cone angle: {:.1}°\n", cone_angle * 180.0 / PI);

    // === 5. Power Output ===
    println!("=== RF Power Output ===");

    // GMR oscillation
    let gmr_ratio = structure.gmr_ratio(cone_angle);
    let r0 = 10.0; // Base resistance (Ω)
    let delta_r = r0 * gmr_ratio;

    println!("  Resistance modulation ΔR: {:.3} Ω", delta_r);

    // AC voltage amplitude
    let v_rf = i_dc * delta_r / 2.0; // Factor of 2 from time averaging
    println!("  RF voltage amplitude: {:.2} mV", v_rf * 1e3);

    // RF power
    let p_rf = v_rf * v_rf / r0;
    println!("  RF power output: {:.2} nW", p_rf * 1e9);

    // DC power consumption
    let p_dc = i_dc * i_dc * r0;
    println!("  DC power input: {:.2} µW", p_dc * 1e6);

    // Efficiency
    let efficiency = (p_rf / p_dc) * 100.0;
    println!("  Conversion efficiency: {:.3}%\n", efficiency);

    // === 6. Tunability ===
    println!("=== Frequency Tunability ===");

    println!("\nCurrent tuning (at H = {:.0} Oe):", h_ext * 10000.0);
    for factor in [1.2, 1.5, 2.0, 3.0] {
        let i = factor * i_critical;
        let f = gamma_hz * h_eff * (1.0 + 0.1 * (factor - 1.0)); // Simplified
        println!("  I = {:.1} mA → f = {:.2} GHz", i * 1e3, f * 1e-9);
    }

    println!("\nField tuning (at I = {:.2} mA):", i_dc * 1e3);
    for h_field in [0.05, 0.1, 0.15, 0.2] {
        let f = gamma_hz * h_field;
        println!(
            "  H = {:.0} Oe → f = {:.2} GHz",
            h_field * 10000.0,
            f * 1e-9
        );
    }

    // === 7. Applications ===
    println!("\n=== Applications ===");

    println!("\n1. Microwave Sources:");
    println!("  • Frequency: 1-40 GHz (tunable)");
    println!("  • Power: 1-100 nW");
    println!("  • Size: Nanoscale (50-200 nm)");
    println!("  • Integration: On-chip RF sources");

    println!("\n2. Magnetic Field Sensors:");
    println!("  • Sensitivity: Field → frequency shift");
    println!("  • Range: 10-10000 Oe");
    println!("  • Precision: Limited by linewidth");

    println!("\n3. Neuromorphic Computing:");
    println!("  • Oscillatory neurons");
    println!("  • Frequency-based encoding");
    println!("  • Coupled oscillator networks");

    // === 8. Performance Metrics ===
    println!("\n=== Performance Summary ===");

    println!("\nKey Metrics:");
    println!("  Frequency range: 1-40 GHz ✓");
    println!(
        "  Current threshold: {:.1} MA/cm² (moderate)",
        (i_critical / area) * 1e-10
    );
    println!("  Power efficiency: {:.3}% (low)", efficiency);
    println!("  Linewidth: {:.0} MHz (moderate)", delta_f * 1e-6);
    println!("  Size: Nanoscale ✓");

    println!("\nChallenges:");
    println!("  • Thermal noise at room temperature");
    println!("  • Phase noise limits coherence");
    println!("  • Low output power");
    println!("  ✓ Solution: Synchronized oscillator arrays");

    println!("\nAdvantages:");
    println!("  ✓ Nanoscale size");
    println!("  ✓ Voltage/current tunable");
    println!("  ✓ CMOS compatible");
    println!("  ✓ Wide frequency range");

    println!("\n=== Conclusion ===");
    println!("Spin-torque nano-oscillators demonstrate the rich physics of");
    println!("spin-transfer torque and offer promising applications in");
    println!("nanoscale microwave technology and neuromorphic computing.");
}
