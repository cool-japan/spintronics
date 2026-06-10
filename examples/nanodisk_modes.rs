//! Spin Wave Modes in a Magnetic Nanodisk
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Spin Wave Theory / Confined Geometries
//! **Physics**: Radial Bessel modes × azimuthal whispering-gallery, Kalinikos-Slavin dispersion
//!
//! ## Background
//!
//! In a magnetic nanodisk of radius R with perpendicular magnetization, confined
//! spin wave modes are indexed by `(n_rad, m_az)`:
//!   - `n_rad ≥ 1` : radial quantum number, with k_rad = α_{m_az, n_rad} / R
//!   - `m_az ≥ 0`  : azimuthal index (ψ ∝ exp(i m_az θ))
//!
//! where α_{m, n} is the n-th zero of the Bessel function J_m.
//!
//! Each mode frequency is given by the standard Kalinikos-Slavin dispersion
//! evaluated at the quantized k = k_rad:
//!   ω² = ω_H (ω_H + ω_M · F(k)) + (D · k²)·(ω_H + ω_M·F(k))
//! where F(k) is the dipolar form factor and D = 2|γ| A_ex/M_s the exchange stiffness.
//!
//! For a YIG 100 nm-radius / 20 nm-thick disk the lowest mode (n=1, m=0)
//! sits in the few-GHz range — exactly what spin-wave nano-oscillators target.
//!
//! ## References
//! - Demidov & Demokritov, IEEE Trans. Magn. 51, 0800215 (2015)
//! - Klingler, Chumak, Mewes et al., APL 110, 092409 (2017)

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  Confined Spin Wave Modes in a YIG Nanodisk");
    println!("=============================================================");

    // -------------------------------------------------------------------------
    // Section 1: Construct YIG 100 nm-radius nanodisk
    // -------------------------------------------------------------------------
    println!("\n--- Section 1: YIG Nanodisk (R = 100 nm, t = 20 nm) ---\n");

    let disk = NanodiskSpinWaves::yig_nanodisk(100.0)?;
    println!("  Radius:       {:>10.2e} m", disk.radius);
    println!("  Thickness:    {:>10.2e} m", disk.thickness);
    println!("  H_ext:        {:>10.2e} A/m (perpendicular)", disk.h_ext);
    println!("  M_s:          {:>10.2e} A/m", disk.ms);
    println!("  A_ex:         {:>10.2e} J/m", disk.a_ex);
    println!("  α (Gilbert):  {:>10.2e}", disk.alpha);
    println!("  ω_H:          {:>10.3e} rad/s", disk.omega_h());
    println!("  ω_M:          {:>10.3e} rad/s", disk.omega_m());

    // -------------------------------------------------------------------------
    // Section 2: Sweep (n_rad, m_az) up to (3, 4) and tabulate frequencies
    // -------------------------------------------------------------------------
    println!("\n--- Section 2: Mode Spectrum (n_rad, m_az) ---\n");

    let n_max = 3;
    let m_max = 4;
    let spectrum = disk.mode_spectrum(n_max, m_max)?;

    println!(
        "  {:>6}  {:>6}  {:>12}  {:>12}  {:>12}",
        "n_rad", "m_az", "k (1/m)", "ω (Grad/s)", "f (GHz)"
    );
    println!("  {}", "-".repeat(56));
    for &(n, m, omega) in spectrum.iter().take(12) {
        let k = disk.radial_wavevector(n, m)?;
        let f_ghz = omega / (2.0 * std::f64::consts::PI) * 1e-9;
        println!(
            "  {:>6}  {:>6}  {:>12.3e}  {:>12.3}  {:>12.3}",
            n,
            m,
            k,
            omega * 1e-9,
            f_ghz
        );
    }

    // Lowest mode identification
    let (n_low, m_low, omega_low) = spectrum[0];
    let f_low_ghz = omega_low / (2.0 * std::f64::consts::PI) * 1e-9;
    println!("\n  → Lowest-frequency mode: (n={n_low}, m={m_low}) at f = {f_low_ghz:.3} GHz");

    // -------------------------------------------------------------------------
    // Section 3: Group velocity & propagation length for fundamental
    // -------------------------------------------------------------------------
    println!("\n--- Section 3: Mode Properties (Fundamental) ---\n");

    let vg = disk.group_velocity(1, 0)?;
    let lp = disk.propagation_length(1, 0)?;
    println!("  Fundamental (n=1, m=0):");
    println!("    Group velocity v_g = {vg:.3e} m/s");
    println!(
        "    Propagation length L_p = {lp:.3e} m  (~{:.1} nm)",
        lp * 1e9
    );
    println!(
        "    L_p / R = {:.2}  (>1 → mode survives a full disk radius)",
        lp / disk.radius
    );

    // -------------------------------------------------------------------------
    // Section 4: Mode profile slice for fundamental
    // -------------------------------------------------------------------------
    println!("\n--- Section 4: Radial Profile of Fundamental Mode ---\n");
    println!("  {:>8}  {:>14}", "r/R", "|profile|²");
    println!("  {}", "-".repeat(26));
    for i in 0..=10 {
        let r_norm = i as f64 / 10.0;
        let r = r_norm * disk.radius;
        let psi = disk.mode_profile(1, 0, r, 0.0)?;
        println!("  {:>8.2}  {:>14.4e}", r_norm, psi * psi);
    }

    // -------------------------------------------------------------------------
    // Section 5: Compare across materials (Permalloy, CoFeB) at fixed radius
    // -------------------------------------------------------------------------
    println!("\n--- Section 5: Material Comparison (R = 100 nm, fundamental) ---\n");

    let yig = NanodiskSpinWaves::yig_nanodisk(100.0)?;
    let nife = NanodiskSpinWaves::permalloy_nanodisk(100.0)?;
    let cofeb = NanodiskSpinWaves::cofeb_nanodisk(100.0)?;

    println!(
        "  {:>15}  {:>10}  {:>10}  {:>10}",
        "Material", "f (GHz)", "v_g (m/s)", "L_p (nm)"
    );
    println!("  {}", "-".repeat(50));
    for (label, m) in [("YIG", &yig), ("Permalloy", &nife), ("CoFeB", &cofeb)] {
        let omega = m.mode_frequency(1, 0)?;
        let vg = m.group_velocity(1, 0)?;
        let lp = m.propagation_length(1, 0)?;
        println!(
            "  {:>15}  {:>10.3}  {:>10.2e}  {:>10.1}",
            label,
            omega / (2.0 * std::f64::consts::PI) * 1e-9,
            vg,
            lp * 1e9
        );
    }

    println!("\n=============================================================");
    println!("  Done. YIG nanodisk fundamental mode in low-GHz regime,");
    println!("  consistent with Demidov & Demokritov 2015.");
    println!("=============================================================\n");

    Ok(())
}
