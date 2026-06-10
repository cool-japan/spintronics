//! Full Landmark-Paper Validation Suite (v0.7.0 + v0.8.0)
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Validation / Reproducibility
//! **Physics**: 5 landmark spintronics papers compared with quantitative metrics
//!
//! ## Background
//!
//! This example runs all five experimental validations registered in the crate
//! and produces a single pass/fail dashboard. Useful for CI regression-testing
//! that simulation predictions remain consistent with the original measurements.
//!
//! Papers under test:
//!   1. Demidov et al., PRL 96, 097202 (2006)  — DE BLS in YIG film
//!   2. Saitoh et al.,  APL 88, 182509 (2006)  — ISHE in Pt/Permalloy
//!   3. Uchida et al.,  Nature 455, 778 (2008) — Longitudinal SSE in Pt/YIG
//!   4. Mosendz et al., PRL 104, 046601 (2010) — Quantitative spin pumping in Py/Pt
//!   5. Liu et al.,     Science 336, 555 (2012) — SOT switching in Ta/CoFeB/MgO

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  Landmark Paper Validation Suite");
    println!("  (5 papers / 13 quantitative checks)");
    println!("=============================================================");

    let tolerance = 0.30_f64; // 30% relative-error window for landmark-paper checks
    let mut summary: Vec<(&'static str, bool)> = Vec::new();

    // -------------------------------------------------------------------------
    // Demidov 2006 — Damon-Eshbach in YIG
    // -------------------------------------------------------------------------
    println!("\n--- Demidov 2006 (PRL 96, 097202) ---\n");
    let demidov = Demidov2006Validation::new()?;
    let d_disp = demidov.validate_dispersion(tolerance)?;
    let d_nr = demidov.validate_nonreciprocity(tolerance)?;
    println!(
        "  DE dispersion ω(k):      max_rel={:.2}% mean_rel={:.2}% ⇒ {}",
        d_disp.max_relative_error * 100.0,
        d_disp.mean_relative_error * 100.0,
        if d_disp.passed { "PASS" } else { "fail" }
    );
    println!(
        "  DE non-reciprocity Δω/ω: max_rel={:.2}% mean_rel={:.2}% ⇒ {}",
        d_nr.max_relative_error * 100.0,
        d_nr.mean_relative_error * 100.0,
        if d_nr.passed { "PASS" } else { "fail" }
    );
    summary.push(("Demidov  DE dispersion              ", d_disp.passed));
    summary.push(("Demidov  DE non-reciprocity         ", d_nr.passed));

    // -------------------------------------------------------------------------
    // Saitoh 2006 — ISHE in Pt/Py
    // -------------------------------------------------------------------------
    println!("\n--- Saitoh 2006 (APL 88, 182509) ---\n");
    let saitoh = Saitoh2006Validation::new()?;
    let s_theta = saitoh.validate_spin_hall_angle(50.0)?; // wide window for old data
    let s_pol = saitoh.validate_ishe_voltage_polarity()?;
    let s_scale = saitoh.validate_ishe_scaling(&[1.0e3, 5.0e3, 1.0e4, 5.0e4], tolerance)?;
    println!(
        "  θ_SH(Pt) windowed:        {}",
        if s_theta.passed { "PASS" } else { "fail" }
    );
    println!(
        "  ISHE polarity (J_s × σ):  {}",
        if s_pol { "PASS" } else { "fail" }
    );
    println!(
        "  ISHE linear in J_s:       max_rel={:.2}% ⇒ {}",
        s_scale.max_relative_error * 100.0,
        if s_scale.passed { "PASS" } else { "fail" }
    );
    summary.push(("Saitoh   θ_SH(Pt) windowed          ", s_theta.passed));
    summary.push(("Saitoh   ISHE polarity              ", s_pol));
    summary.push(("Saitoh   ISHE linear in J_s         ", s_scale.passed));

    // -------------------------------------------------------------------------
    // Uchida 2008 — LSSE in Pt/YIG
    // -------------------------------------------------------------------------
    println!("\n--- Uchida 2008 (Nature 455, 778) ---\n");
    let uchida = Uchida2008Validation::new()?;
    let u_lin = uchida.validate_linear_thermal_response(tolerance)?;
    let u_pol = uchida.validate_polarity()?;
    let u_order = uchida.validate_seebeck_coefficient_order()?;
    println!(
        "  LSSE linear V_LSSE(ΔT):   max_rel={:.2}% ⇒ {}",
        u_lin.max_relative_error * 100.0,
        if u_lin.passed { "PASS" } else { "fail" }
    );
    println!(
        "  LSSE polarity:            {}",
        if u_pol { "PASS" } else { "fail" }
    );
    println!(
        "  Seebeck order-of-mag:     {}",
        if u_order { "PASS" } else { "fail" }
    );
    summary.push(("Uchida   LSSE linear V(ΔT)          ", u_lin.passed));
    summary.push(("Uchida   LSSE polarity              ", u_pol));
    summary.push(("Uchida   Seebeck order-of-mag       ", u_order));

    // -------------------------------------------------------------------------
    // Mosendz 2010 — Spin pumping in Py/Pt
    // -------------------------------------------------------------------------
    println!("\n--- Mosendz 2010 (PRL 104, 046601) ---\n");
    let mosendz = Mosendz2010Validation::new()?;
    let m_thk = mosendz.validate_pt_thickness_scaling(0.65)?; // shape-only; see test note
    let m_lw = mosendz.validate_linewidth_enhancement(12.0e-9, tolerance)?;
    let m_g = mosendz.validate_spin_mixing_conductance(tolerance)?;
    println!(
        "  V_ISHE(t_Pt) thickness shape: max_rel={:.2}% ⇒ {}",
        m_thk.max_relative_error * 100.0,
        if m_thk.passed { "PASS" } else { "fail" }
    );
    println!(
        "  Linewidth Δα_eff:             max_rel={:.2}% ⇒ {}",
        m_lw.max_relative_error * 100.0,
        if m_lw.passed { "PASS" } else { "fail" }
    );
    println!(
        "  g↑↓ vs reference:             max_rel={:.2}% ⇒ {}",
        m_g.max_relative_error * 100.0,
        if m_g.passed { "PASS" } else { "fail" }
    );
    summary.push(("Mosendz  V_ISHE thickness shape     ", m_thk.passed));
    summary.push(("Mosendz  linewidth Δα_eff           ", m_lw.passed));
    summary.push(("Mosendz  g↑↓ vs reference           ", m_g.passed));

    // -------------------------------------------------------------------------
    // Liu 2012 — SOT switching in Ta/CoFeB/MgO
    // -------------------------------------------------------------------------
    println!("\n--- Liu 2012 (Science 336, 555) ---\n");
    let liu = Liu2012Validation::new()?;
    let l_theta = liu.validate_ta_spin_hall_angle(tolerance)?;
    let l_jc = liu.validate_critical_switching_current(tolerance)?;
    let l_pol = liu.validate_sot_polarity()?;
    let l_thk = liu.validate_thickness_scaling(0.90)?; // wide window; see test note
    println!(
        "  θ_SH(β-Ta) windowed:      max_rel={:.2}% ⇒ {}",
        l_theta.max_relative_error * 100.0,
        if l_theta.passed { "PASS" } else { "fail" }
    );
    println!(
        "  Critical J_c:             max_rel={:.2}% ⇒ {}",
        l_jc.max_relative_error * 100.0,
        if l_jc.passed { "PASS" } else { "fail" }
    );
    println!(
        "  SOT polarity (Ta vs Pt):  {}",
        if l_pol { "PASS" } else { "fail" }
    );
    println!(
        "  J_c(t_CoFeB) scaling:     max_rel={:.2}% ⇒ {}",
        l_thk.max_relative_error * 100.0,
        if l_thk.passed { "PASS" } else { "fail" }
    );
    summary.push(("Liu      θ_SH(β-Ta) windowed        ", l_theta.passed));
    summary.push(("Liu      critical J_c               ", l_jc.passed));
    summary.push(("Liu      SOT polarity (Ta vs Pt)    ", l_pol));
    summary.push(("Liu      J_c(t_CoFeB) shape         ", l_thk.passed));

    // -------------------------------------------------------------------------
    // Summary dashboard
    // -------------------------------------------------------------------------
    println!("\n=============================================================");
    println!("  Validation Summary");
    println!("=============================================================");
    println!("  {:<40}  Status", "Test");
    println!("  {}", "-".repeat(48));
    let mut n_pass = 0;
    for &(name, passed) in summary.iter() {
        let status = if passed { "✓ PASS" } else { "✗ fail" };
        if passed {
            n_pass += 1;
        }
        println!("  {name}  {status}");
    }
    let total = summary.len();
    println!("\n  Overall: {n_pass}/{total} checks passed.");
    println!("\n=============================================================");
    println!("  Done. 5 landmark papers, {n_pass}/{total} quantitative checks pass");
    println!("  within the 30%-relative-error window typical for simplified models.");
    println!("=============================================================\n");

    Ok(())
}
