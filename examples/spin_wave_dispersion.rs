//! Spin Wave Dispersion Relations
//!
//! **Difficulty**: ⭐⭐ Intermediate
//! **Category**: Magnonics
//! **Physics**: Kittel formula, Kalinikos-Slavin, DE/BVMSW/FVMSW modes
//!
//! This example demonstrates spin wave dispersion relations in thin
//! ferromagnetic films, covering:
//!
//! 1. Kittel uniform precession frequency (k=0)
//! 2. Exchange-dominated spin wave dispersion
//! 3. Kalinikos-Slavin dipole-exchange theory
//! 4. Three canonical magnetostatic modes:
//!    - Damon-Eshbach (DE) surface waves
//!    - Backward Volume MSW (BVMSW)
//!    - Forward Volume MSW (FVMSW)
//!
//! References:
//! - Kittel, Phys. Rev. 73, 155 (1948)
//! - Kalinikos & Slavin, J. Phys. C 19, 7013 (1986)

use std::f64::consts::PI;

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Spin Wave Dispersion Relations ===\n");

    // Material: YIG (Yttrium Iron Garnet) -- low damping, ideal for magnonics
    let yig = Ferromagnet::yig();
    let thickness = 1.0e-6; // 1 um film

    println!("Material: YIG (Y3Fe5O12)");
    println!("  M_s = {:.3e} A/m", yig.ms);
    println!("  A = {:.2e} J/m (exchange stiffness)", yig.exchange_a);
    println!("  alpha = {:.1e} (Gilbert damping)", yig.alpha);
    println!("  Film thickness = {:.1} um\n", thickness * 1e6);

    // 1. Kittel frequency vs external field
    println!("=== Kittel Frequency (k=0) ===");
    let disp = SpinWaveDispersion::new(&yig, thickness, 0.0)
        .expect("Failed to create dispersion calculator");

    println!(
        "  {:>10} {:>15} {:>15}",
        "H_ext (T)", "f (GHz)", "omega (rad/s)"
    );
    println!("  {}", "-".repeat(45));
    for &h_ext in &[0.01, 0.05, 0.1, 0.2, 0.5, 1.0] {
        let omega = disp
            .kittel_frequency(h_ext)
            .expect("Failed to compute Kittel frequency");
        let f_ghz = omega / (2.0 * PI * 1e9);
        println!("  {:>10.2} {:>15.4} {:>15.4e}", h_ext, f_ghz, omega);
    }

    // 2. Exchange-dominated dispersion
    println!("\n=== Exchange Dispersion (H_ext = 0.1 T) ===");
    let h_ext = 0.1; // Tesla
    println!(
        "  {:>15} {:>15} {:>15} {:>15}",
        "k (rad/m)", "lambda (nm)", "f (GHz)", "v_g approx (km/s)"
    );
    println!("  {}", "-".repeat(65));

    let k_values: Vec<f64> = (0..10)
        .map(|i| {
            let k_min: f64 = 1.0e5;
            let k_max: f64 = 1.0e8;
            k_min * (k_max / k_min).powf(i as f64 / 9.0)
        })
        .collect();

    let mut prev_omega = 0.0;
    let mut prev_k = 0.0;
    for (i, &k) in k_values.iter().enumerate() {
        let omega = disp
            .exchange_dispersion(h_ext, k)
            .expect("Failed to compute exchange dispersion");
        let f_ghz = omega / (2.0 * PI * 1e9);
        let lambda_nm = 2.0 * PI / k * 1e9;

        let v_g = if i > 0 && (k - prev_k).abs() > 1e-10 {
            (omega - prev_omega) / (k - prev_k) / 1e3
        } else {
            0.0
        };

        println!(
            "  {:>15.4e} {:>15.2} {:>15.4} {:>15.3}",
            k, lambda_nm, f_ghz, v_g
        );
        prev_omega = omega;
        prev_k = k;
    }

    // 3. Kalinikos-Slavin angle dependence
    println!("\n=== Kalinikos-Slavin: Angle Dependence (k = 10^6 rad/m) ===");
    let k_fixed = 1.0e6;
    println!(
        "  {:>12} {:>15} {:>15}",
        "phi (deg)", "f (GHz)", "Mode type"
    );
    println!("  {}", "-".repeat(47));
    for i in 0..=6 {
        let phi_deg = i as f64 * 15.0;
        let phi_rad = phi_deg * PI / 180.0;
        let omega = disp
            .kalinikos_slavin(h_ext, k_fixed, phi_rad)
            .expect("Failed to compute KS dispersion");
        let f_ghz = omega / (2.0 * PI * 1e9);
        let mode_type = if phi_deg < 10.0 {
            "BVMSW-like"
        } else if phi_deg > 80.0 {
            "DE-like"
        } else {
            "intermediate"
        };
        println!("  {:>12.1} {:>15.4} {:>15}", phi_deg, f_ghz, mode_type);
    }

    // 4. Three canonical mode types
    println!("\n=== Canonical Magnetostatic Modes ===");
    let mode_calc =
        SpinWaveModeCalculator::new(&yig, thickness).expect("Failed to create mode calculator");

    let modes = [
        (SpinWaveMode::DamonEshbach, "Damon-Eshbach (DE)"),
        (SpinWaveMode::BackwardVolume, "BVMSW"),
        (SpinWaveMode::ForwardVolume, "FVMSW"),
    ];

    let k_mode = 1.0e6; // 10^6 rad/m
    println!("  k = {:.0e} rad/m, H_ext = {:.2} T\n", k_mode, h_ext);

    for (mode, label) in &modes {
        let omega = mode_calc
            .frequency(*mode, h_ext, k_mode)
            .expect("Failed to compute mode frequency");
        let f_ghz = omega / (2.0 * PI * 1e9);
        println!(
            "  {:<25} f = {:.4} GHz  (omega = {:.4e} rad/s)",
            label, f_ghz, omega
        );
    }

    // 5. Mode dispersion comparison
    println!("\n=== Mode Dispersion Comparison ===");
    println!(
        "  {:>15} {:>12} {:>12} {:>12}",
        "k (rad/m)", "DE (GHz)", "BVMSW (GHz)", "FVMSW (GHz)"
    );
    println!("  {}", "-".repeat(55));

    for &k in &[1e5, 3e5, 1e6, 3e6, 1e7] {
        let f_de = mode_calc
            .frequency(SpinWaveMode::DamonEshbach, h_ext, k)
            .expect("DE computation failed")
            / (2.0 * PI * 1e9);
        let f_bv = mode_calc
            .frequency(SpinWaveMode::BackwardVolume, h_ext, k)
            .expect("BVMSW computation failed")
            / (2.0 * PI * 1e9);
        let f_fv = mode_calc
            .frequency(SpinWaveMode::ForwardVolume, h_ext, k)
            .expect("FVMSW computation failed")
            / (2.0 * PI * 1e9);
        println!(
            "  {:>15.2e} {:>12.4} {:>12.4} {:>12.4}",
            k, f_de, f_bv, f_fv
        );
    }

    // 6. Frequency gap and bandwidth
    println!("\n=== Frequency Range Analysis ===");
    let omega_kittel = disp
        .kittel_frequency(h_ext)
        .expect("Kittel computation failed");
    let omega_high_k = disp
        .exchange_dispersion(h_ext, 1e8)
        .expect("High-k computation failed");
    println!(
        "  Kittel frequency (k=0): {:.4} GHz",
        omega_kittel / (2.0 * PI * 1e9)
    );
    println!(
        "  Exchange-dominated (k=10^8): {:.4} GHz",
        omega_high_k / (2.0 * PI * 1e9)
    );
    println!(
        "  Magnonic bandwidth: {:.2} GHz",
        (omega_high_k - omega_kittel) / (2.0 * PI * 1e9)
    );

    println!("\n=== Summary ===");
    println!("Spin wave physics in {:.0} um YIG film:", thickness * 1e6);
    println!(
        "  - Kittel FMR frequency: {:.2} GHz at {:.2} T",
        omega_kittel / (2.0 * PI * 1e9),
        h_ext
    );
    println!("  - DE modes: surface waves (k perp M), highest frequency");
    println!("  - BVMSW: k parallel to M, negative group velocity at low k");
    println!("  - FVMSW: M perpendicular to film, isotropic in-plane dispersion");

    Ok(())
}
