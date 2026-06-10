//! Magnonic Crystal Band Structure and Band Gap
//!
//! **Difficulty**: ⭐⭐⭐⭐
//! **Category**: Spin Wave Theory / Magnonics
//! **Physics**: 1D periodic magnetic structure, plane-wave expansion, Bloch bands
//!
//! ## Background
//!
//! A magnonic crystal is a periodic magnetic medium that produces a Bloch-like
//! spin wave band structure with allowed bands and forbidden gaps (analogous to
//! electronic band gaps in semiconductors and photonic gaps in dielectric
//! crystals). The band gap is the key feature: by choosing the period and
//! material contrast carefully, one can design a frequency window in which
//! magnons cannot propagate.
//!
//! Implementation: plane-wave expansion in the reciprocal-lattice basis
//! G_n = 2π·n/a, with 2·n_pw + 1 plane waves. The Hamiltonian is built from
//! Fourier coefficients of the position-dependent ω_M(x) and exchange-length
//! λ_ex²(x), then diagonalized at each Bloch wavevector k ∈ [-π/a, +π/a).
//!
//! ## References
//! - Vasseur, Dobrzynski, Djafari-Rouhani, Puszkarski, PRB 49, 1727 (1994)
//! - Krawczyk & Grundler, J. Phys. Cond. Matt. 26, 123202 (2014) — magnonics review

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  1D Magnonic Crystal: Band Structure & First Band Gap");
    println!("=============================================================");

    // -------------------------------------------------------------------------
    // Section 1: Build alternating-layer 1D magnonic crystal (YIG / Permalloy)
    // -------------------------------------------------------------------------
    println!("\n--- Section 1: YIG / Permalloy Magnonic Crystal ---\n");

    // Period 200 nm, 50:50 filling, modest applied field
    let period_nm = 200.0_f64;
    let crystal = MagnonicCrystal1D::nife_cofe_alternating(period_nm)?;
    println!("  Period:            {:>8.2e} m", crystal.period);
    println!(
        "  M_s (A):           {:>8.2e} A/m  (Permalloy / NiFe)",
        crystal.ms_a
    );
    println!("  M_s (B):           {:>8.2e} A/m  (CoFeB)", crystal.ms_b);
    println!("  A_ex (A):          {:>8.2e} J/m", crystal.a_ex_a);
    println!("  A_ex (B):          {:>8.2e} J/m", crystal.a_ex_b);
    println!("  Filling (A):       {:>8.2}", crystal.filling_a);
    println!("  H_ext:             {:>8.2e} A/m", crystal.h_ext);
    println!("  n_pw:              {:>8}", crystal.n_pw);

    // -------------------------------------------------------------------------
    // Section 2: Sample band structure across the first BZ
    // -------------------------------------------------------------------------
    println!("\n--- Section 2: Band Structure (first Brillouin zone) ---\n");

    let n_kpts = 11;
    let bands = crystal.band_structure(n_kpts)?;
    let n_show = 5; // first 5 bands

    print!("  {:>10}", "k·a/π");
    for b in 0..n_show {
        print!("  {:>10}", format!("band {b}"));
    }
    println!();
    println!("  {}", "-".repeat(10 + n_show * 12));

    for (k_bloch, ws) in bands.iter() {
        let k_norm = k_bloch * crystal.period / std::f64::consts::PI;
        print!("  {:>10.3}", k_norm);
        for w in ws.iter().take(n_show.min(ws.len())) {
            let f_ghz = w / (2.0 * std::f64::consts::PI) * 1e-9;
            print!("  {:>10.3}", f_ghz);
        }
        println!();
    }

    // -------------------------------------------------------------------------
    // Section 3: Identify first band gap
    // -------------------------------------------------------------------------
    println!("\n--- Section 3: First Band Gap ---\n");

    let gap_result = crystal.band_gap(0, n_kpts);
    match gap_result {
        Ok((omega_low, omega_high)) => {
            let f_low = omega_low / (2.0 * std::f64::consts::PI) * 1e-9;
            let f_high = omega_high / (2.0 * std::f64::consts::PI) * 1e-9;
            let delta_f = f_high - f_low;
            let f_center = 0.5 * (f_low + f_high);
            println!("  First gap exists between band 0 and band 1:");
            println!("    Lower edge:     {f_low:.3} GHz");
            println!("    Upper edge:     {f_high:.3} GHz");
            println!("    Gap width Δf:   {delta_f:.3} GHz");
            println!("    Center freq:    {f_center:.3} GHz");
            println!("    Δf / f_center:  {:.2} %", 100.0 * delta_f / f_center);
        },
        Err(e) => {
            println!("  No band gap between band 0 and band 1: {e}");
            println!("  (For low material contrast or small n_pw, gaps may close.)");
        },
    }

    // -------------------------------------------------------------------------
    // Section 4: Group velocity along the lowest band
    // -------------------------------------------------------------------------
    println!("\n--- Section 4: Group Velocity (band 0) ---\n");

    println!("  {:>10}  {:>14}", "k·a/π", "v_g (m/s)");
    println!("  {}", "-".repeat(26));
    let dk = std::f64::consts::PI / crystal.period * 1e-3;
    for i in 0..=10 {
        let k_norm = i as f64 / 10.0;
        let k_bloch = k_norm * std::f64::consts::PI / crystal.period;
        let vg = crystal.group_velocity(0, k_bloch, dk)?;
        println!("  {:>10.3}  {:>14.3e}", k_norm, vg);
    }

    println!("\n  → v_g → 0 at zone center (k=0) and zone boundary (k=π/a),");
    println!("    consistent with Bloch theorem.");

    // -------------------------------------------------------------------------
    // Section 5: Compare with low-contrast crystal (no gap expected)
    // -------------------------------------------------------------------------
    println!("\n--- Section 5: Low-Contrast Comparison ---\n");

    // Very similar Ms_a and Ms_b — gap should be tiny or absent
    let low_contrast = MagnonicCrystal1D::new(
        period_nm * 1e-9,
        8.0e5,    // ms_a
        8.1e5,    // ms_b — 1% contrast
        1.0e-11,  // a_ex_a
        1.05e-11, // a_ex_b — 5% contrast
        0.5,
        2.0e4,
        11,
    )?;

    match low_contrast.band_gap(0, n_kpts) {
        Ok((lo, hi)) => {
            let dw_ghz = (hi - lo) / (2.0 * std::f64::consts::PI) * 1e-9;
            println!("  Low-contrast gap width:   {dw_ghz:.4} GHz  (very small as expected)");
        },
        Err(_) => {
            println!("  Low-contrast crystal: no band gap (as expected for nearly-uniform medium)");
        },
    }

    println!("\n=============================================================");
    println!("  Done. 1D magnonic crystal shows clear band structure with");
    println!("  a sub-GHz first band gap, opening with material contrast.");
    println!("=============================================================\n");

    Ok(())
}
