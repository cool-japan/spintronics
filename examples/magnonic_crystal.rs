//! Magnonic Crystal Band Structure Calculator
//!
//! This example demonstrates:
//! - Creating a periodic magnetic structure (magnonic crystal)
//! - Calculating magnon dispersion relation
//! - Finding bandgaps from Bragg scattering
//! - Comparing with photonic/phononic analogs
//!
//! ## Physics:
//! Magnonic crystals are periodic modulations of magnetic properties that
//! create bandgaps for spin waves, analogous to photonic crystals for light.
//! Applications include magnon filters, isolators, and logic devices.
//!
//! ## References:
//! - V. V. Kruglyak et al., "Magnonics", J. Phys. D: Appl. Phys. 43, 264001 (2010)
//! - S. A. Nikitov et al., "Magnonics: a new research area in spintronics and
//!   spin wave electronics", Phys. Usp. 58, 1002 (2015)
//!
//! Run with: cargo run --example magnonic_crystal

use spintronics::prelude::*;
use std::f64::consts::PI;

fn main() {
    println!("=== Magnonic Crystal Band Structure ===\n");

    // === 1. Crystal Structure ===
    println!("=== Magnonic Crystal Design ===");

    // Bilayer system: alternating Py and CoFeB
    let py = Ferromagnet::permalloy();
    let cofeb = Ferromagnet::cofeb();

    println!("Material A: Permalloy");
    println!("  M_s = {:.2} kA/m", py.ms * 1e-3);
    println!("  A_ex = {:.2} pJ/m", py.exchange_a * 1e12);
    println!("  K_u = {:.2} kJ/m³", py.anisotropy_k * 1e-3);

    println!("\nMaterial B: CoFeB");
    println!("  M_s = {:.2} kA/m", cofeb.ms * 1e-3);
    println!("  A_ex = {:.2} pJ/m", cofeb.exchange_a * 1e12);
    println!("  K_u = {:.2} kJ/m³\n", cofeb.anisotropy_k * 1e-3);

    // Periodicity
    let lattice_constant = 100e-9; // 100 nm period
    let d_a = 50e-9; // Py thickness
    let d_b = 50e-9; // CoFeB thickness

    println!("Periodicity:");
    println!("  Lattice constant: {:.0} nm", lattice_constant * 1e9);
    println!("  Py layer: {:.0} nm", d_a * 1e9);
    println!("  CoFeB layer: {:.0} nm\n", d_b * 1e9);

    // === 2. Magnon Dispersion (Uniform Film) ===
    println!("=== Magnon Dispersion Relation ===");

    let h_ext = 0.1; // 1000 Oe external field
    let _gamma = GAMMA; // rad/(s·T)

    println!("External field: {:.0} Oe\n", h_ext * 10000.0);

    // Dispersion for Permalloy
    println!("Permalloy magnons:");
    for k in [0.0, 1e6, 2e6, 5e6, 1e7] {
        let omega = magnon_frequency(&py, h_ext, k);
        let f_ghz = omega / (2.0 * PI * 1e9);
        let v_group = magnon_group_velocity(&py, h_ext, k);
        println!(
            "  k = {:.1e} m⁻¹ → f = {:.2} GHz, v_g = {:.1} km/s",
            k,
            f_ghz,
            v_group * 1e-3
        );
    }

    // Dispersion for CoFeB
    println!("\nCoFeB magnons:");
    for k in [0.0, 1e6, 2e6, 5e6, 1e7] {
        let omega = magnon_frequency(&cofeb, h_ext, k);
        let f_ghz = omega / (2.0 * PI * 1e9);
        let v_group = magnon_group_velocity(&cofeb, h_ext, k);
        println!(
            "  k = {:.1e} m⁻¹ → f = {:.2} GHz, v_g = {:.1} km/s",
            k,
            f_ghz,
            v_group * 1e-3
        );
    }

    // === 3. Brillouin Zone and Bandgaps ===
    println!("\n=== Brillouin Zone Analysis ===");

    let k_bragg = PI / lattice_constant;
    let _k_max = 2.0 * k_bragg;

    println!("First Brillouin zone:");
    println!("  k_Bragg = π/a = {:.2e} m⁻¹", k_bragg);
    println!("  BZ edge: ±{:.2e} m⁻¹\n", k_bragg);

    // Bandgap estimation (simplified)
    let omega_a = magnon_frequency(&py, h_ext, k_bragg);
    let omega_b = magnon_frequency(&cofeb, h_ext, k_bragg);
    let omega_avg = (omega_a + omega_b) / 2.0;
    let delta_omega = (omega_a - omega_b).abs();

    println!("Bandgap at k = k_Bragg:");
    println!("  ω_Py = {:.2} GHz", omega_a / (2.0 * PI * 1e9));
    println!("  ω_CoFeB = {:.2} GHz", omega_b / (2.0 * PI * 1e9));
    println!(
        "  Center frequency: {:.2} GHz",
        omega_avg / (2.0 * PI * 1e9)
    );
    println!(
        "  Bandgap width: {:.3} GHz\n",
        delta_omega / (2.0 * PI * 1e9)
    );

    // === 4. Band Structure Calculation ===
    println!("=== Band Structure (First Two Bands) ===\n");

    println!("Reduced wavevector (k/k_Bragg) | Lower band (GHz) | Upper band (GHz)");
    println!("---------------------------------------------------------------------");

    for i in 0..11 {
        let k_reduced = i as f64 / 10.0;
        let k = k_reduced * k_bragg;

        // Lower band (mostly Py-like)
        let omega_lower = magnon_frequency(&py, h_ext, k);
        // Upper band (mostly CoFeB-like)
        let omega_upper = magnon_frequency(&cofeb, h_ext, k);

        // Apply bandgap near BZ edge
        let gap_coupling = (k / k_bragg).min(1.0);
        let gap_shift = delta_omega * gap_coupling * 0.5;

        let f_lower = (omega_lower - gap_shift) / (2.0 * PI * 1e9);
        let f_upper = (omega_upper + gap_shift) / (2.0 * PI * 1e9);

        println!(
            "  {:.2}                        | {:.3}            | {:.3}",
            k_reduced, f_lower, f_upper
        );
    }

    // === 5. Density of States ===
    println!("\n=== Magnon Density of States ===");

    let dos_py = magnon_dos(&py, h_ext, 1e7);
    let dos_cofeb = magnon_dos(&cofeb, h_ext, 1e7);

    println!("  Py DOS: {:.2e} states/(J·m³)", dos_py);
    println!("  CoFeB DOS: {:.2e} states/(J·m³)", dos_cofeb);
    println!("  DOS vanishes in bandgap\n");

    // === 6. Applications ===
    println!("=== Magnonic Crystal Applications ===");

    println!("\n1. Magnon Filters:");
    println!(
        "  • Passband: {:.2}-{:.2} GHz",
        (omega_avg - delta_omega) / (2.0 * PI * 1e9),
        (omega_avg + delta_omega) / (2.0 * PI * 1e9)
    );
    println!(
        "  • Stopband: {:.2} GHz width",
        delta_omega / (2.0 * PI * 1e9)
    );
    println!("  • Quality factor: {:.1}", omega_avg / delta_omega);

    println!("\n2. Magnonic Logic:");
    println!("  • AND/OR gates using interference");
    println!("  • XOR from phase accumulation");
    println!("  • NOT gate from edge modes");

    println!("\n3. Slow Light Analog:");
    let v_avg = magnon_group_velocity(&py, h_ext, k_bragg * 0.9);
    println!("  • Group velocity near bandgap: {:.1} m/s", v_avg);
    println!("  • Slow magnon regime (v_g << 10 km/s)");
    println!("  • Enhanced nonlinear effects");

    // === 7. Thermal Magnons ===
    println!("\n=== Thermal Magnon Population ===");

    let t = 300.0; // K
    let kb = KB;
    let hbar = HBAR;

    let omega_thermal = kb * t / hbar;
    let f_thermal = omega_thermal / (2.0 * PI * 1e9);

    println!("  Temperature: {:.0} K", t);
    println!("  Thermal energy: {:.2} meV", kb * t / E_CHARGE * 1e3);
    println!("  Thermal frequency: {:.1} GHz", f_thermal);

    let n_magnons_py = thermal_magnon_number(omega_a, t);
    let n_magnons_cofeb = thermal_magnon_number(omega_b, t);

    println!("  Py thermal population: {:.2e}", n_magnons_py);
    println!("  CoFeB thermal population: {:.2e}", n_magnons_cofeb);

    // === 8. Performance Metrics ===
    println!("\n=== Performance Summary ===");

    let q_factor = omega_avg / delta_omega;
    let contrast = delta_omega / omega_avg;
    let transmission_gap: f64 = 0.01; // 1% in stopband
    let transmission_pass: f64 = 0.95; // 95% in passband

    println!("\nKey Metrics:");
    println!(
        "  Operating frequency: {:.1} GHz",
        omega_avg / (2.0 * PI * 1e9)
    );
    println!("  Bandgap/center frequency: {:.1}%", contrast * 100.0);
    println!("  Quality factor: {:.0}", q_factor);
    println!(
        "  Stopband rejection: {:.0} dB",
        -10.0 * transmission_gap.log10()
    );
    println!(
        "  Passband loss: {:.2} dB",
        -10.0 * transmission_pass.log10()
    );

    println!("\nDesign Trade-offs:");
    println!("  Larger Δω → wider bandgap (choose different materials)");
    println!("  Smaller period → higher frequency bandgap");
    println!("  Thicker layers → sharper band edges");

    println!("\n=== Conclusion ===");
    println!("Magnonic crystals enable control of spin wave propagation");
    println!("through engineered band structures, opening pathways for");
    println!("magnon-based signal processing and quantum information.");
}

// Helper functions

fn magnon_frequency(mat: &Ferromagnet, h_ext: f64, k: f64) -> f64 {
    let gamma = GAMMA;
    let mu0 = 4.0 * PI * 1e-7;

    // Dispersion: ω = γ * sqrt(H(H + M_s) + (2A/M_s)k²)
    let h_term = h_ext * (h_ext + mu0 * mat.ms);
    let exchange_term = (2.0 * mat.exchange_a / (mu0 * mat.ms)) * k * k;

    gamma * (h_term + exchange_term).sqrt()
}

fn magnon_group_velocity(mat: &Ferromagnet, h_ext: f64, k: f64) -> f64 {
    // v_g = dω/dk
    let dk = 1e4;
    let omega_plus = magnon_frequency(mat, h_ext, k + dk);
    let omega_minus = magnon_frequency(mat, h_ext, k - dk);
    (omega_plus - omega_minus) / (2.0 * dk)
}

fn magnon_dos(mat: &Ferromagnet, h_ext: f64, k_max: f64) -> f64 {
    // Simplified 3D DOS
    let omega_max = magnon_frequency(mat, h_ext, k_max);
    let volume = 1.0; // Normalize to 1 m³

    // DOS(ω) ~ ω^(1/2) for magnons
    (volume / (2.0 * PI * PI)) * (omega_max.sqrt())
}

fn thermal_magnon_number(omega: f64, temperature: f64) -> f64 {
    // Bose-Einstein distribution
    let kb = KB;
    let hbar = HBAR;
    let x = hbar * omega / (kb * temperature);

    if x > 100.0 {
        0.0
    } else {
        1.0 / (x.exp() - 1.0)
    }
}
