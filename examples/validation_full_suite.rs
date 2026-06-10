//! Full Landmark Paper Validation Suite (v0.7.0 + v0.8.0 + v0.9.0)
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Validation / Reproducibility
//! **Physics**: 7 landmark spintronics papers, quantitative pass/fail dashboard
//!
//! ## Background
//!
//! Extends `validation_landmark_suite.rs` to include the two v0.9.0 additions
//! (Garello 2013 + Boona 2014), bringing total coverage to 7 landmark papers:
//!
//!   1. Demidov  2006 — DE BLS in YIG film              (v0.7.0)
//!   2. Saitoh   2006 — ISHE in Pt/Permalloy             (v0.7.0)
//!   3. Uchida   2008 — Longitudinal SSE in Pt/YIG       (v0.7.0)
//!   4. Mosendz  2010 — Quantitative spin pumping Py/Pt  (v0.8.0)
//!   5. Liu      2012 — SOT switching Ta/CoFeB/MgO       (v0.8.0)
//!   6. Garello  2013 — SOT angular harmonics Pt/Co/AlOx (v0.9.0)
//!   7. Boona    2014 — LSSE in granular YIG/Pt          (v0.9.0)

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  Full Landmark-Paper Validation Suite (v0.9.0)");
    println!("  7 papers, ~20 quantitative checks");
    println!("=============================================================");

    let tolerance = 0.30_f64;
    let mut summary: Vec<(&'static str, bool)> = Vec::new();

    // -------------------------------------------------------------------------
    // 1. Demidov 2006
    // -------------------------------------------------------------------------
    println!("\n--- Demidov 2006 (PRL 96, 097202) ---\n");
    let demidov = Demidov2006Validation::new()?;
    let d_disp = demidov.validate_dispersion(tolerance)?;
    let d_nr = demidov.validate_nonreciprocity(tolerance)?;
    println!(
        "  DE dispersion ω(k):       max_rel={:.2}% ⇒ {}",
        d_disp.max_relative_error * 100.0,
        if d_disp.passed { "PASS" } else { "fail" }
    );
    println!(
        "  DE non-reciprocity Δω/ω:  max_rel={:.2}% ⇒ {}",
        d_nr.max_relative_error * 100.0,
        if d_nr.passed { "PASS" } else { "fail" }
    );
    summary.push(("Demidov   DE dispersion              ", d_disp.passed));
    summary.push(("Demidov   DE non-reciprocity         ", d_nr.passed));

    // -------------------------------------------------------------------------
    // 2. Saitoh 2006
    // -------------------------------------------------------------------------
    println!("\n--- Saitoh 2006 (APL 88, 182509) ---\n");
    let saitoh = Saitoh2006Validation::new()?;
    let s_theta = saitoh.validate_spin_hall_angle(50.0)?;
    let s_pol = saitoh.validate_ishe_voltage_polarity()?;
    let s_scale = saitoh.validate_ishe_scaling(&[1.0e3, 5.0e3, 1.0e4, 5.0e4], tolerance)?;
    println!(
        "  θ_SH(Pt) windowed:        {}",
        if s_theta.passed { "PASS" } else { "fail" }
    );
    println!(
        "  ISHE polarity:            {}",
        if s_pol { "PASS" } else { "fail" }
    );
    println!(
        "  ISHE linear in J_s:       max_rel={:.2}% ⇒ {}",
        s_scale.max_relative_error * 100.0,
        if s_scale.passed { "PASS" } else { "fail" }
    );
    summary.push(("Saitoh    θ_SH(Pt) windowed          ", s_theta.passed));
    summary.push(("Saitoh    ISHE polarity              ", s_pol));
    summary.push(("Saitoh    ISHE linear in J_s         ", s_scale.passed));

    // -------------------------------------------------------------------------
    // 3. Uchida 2008
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
    summary.push(("Uchida    LSSE linear V(ΔT)          ", u_lin.passed));
    summary.push(("Uchida    LSSE polarity              ", u_pol));
    summary.push(("Uchida    Seebeck order-of-mag       ", u_order));

    // -------------------------------------------------------------------------
    // 4. Mosendz 2010
    // -------------------------------------------------------------------------
    println!("\n--- Mosendz 2010 (PRL 104, 046601) ---\n");
    let mosendz = Mosendz2010Validation::new()?;
    let m_thk = mosendz.validate_pt_thickness_scaling(0.65)?;
    let m_lw = mosendz.validate_linewidth_enhancement(12.0e-9, tolerance)?;
    let m_g = mosendz.validate_spin_mixing_conductance(tolerance)?;
    println!(
        "  V_ISHE thickness shape:   max_rel={:.2}% ⇒ {}",
        m_thk.max_relative_error * 100.0,
        if m_thk.passed { "PASS" } else { "fail" }
    );
    println!(
        "  Linewidth Δα_eff:         max_rel={:.2}% ⇒ {}",
        m_lw.max_relative_error * 100.0,
        if m_lw.passed { "PASS" } else { "fail" }
    );
    println!(
        "  g↑↓ vs reference:         max_rel={:.2}% ⇒ {}",
        m_g.max_relative_error * 100.0,
        if m_g.passed { "PASS" } else { "fail" }
    );
    summary.push(("Mosendz   V_ISHE thickness shape     ", m_thk.passed));
    summary.push(("Mosendz   linewidth Δα_eff           ", m_lw.passed));
    summary.push(("Mosendz   g↑↓ vs reference           ", m_g.passed));

    // -------------------------------------------------------------------------
    // 5. Liu 2012
    // -------------------------------------------------------------------------
    println!("\n--- Liu 2012 (Science 336, 555) ---\n");
    let liu = Liu2012Validation::new()?;
    let l_theta = liu.validate_ta_spin_hall_angle(tolerance)?;
    let l_jc = liu.validate_critical_switching_current(tolerance)?;
    let l_pol = liu.validate_sot_polarity()?;
    let l_thk = liu.validate_thickness_scaling(0.90)?;
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
        "  SOT polarity:             {}",
        if l_pol { "PASS" } else { "fail" }
    );
    println!(
        "  J_c(t_CoFeB) scaling:     max_rel={:.2}% ⇒ {}",
        l_thk.max_relative_error * 100.0,
        if l_thk.passed { "PASS" } else { "fail" }
    );
    summary.push(("Liu       θ_SH(β-Ta) windowed        ", l_theta.passed));
    summary.push(("Liu       critical J_c               ", l_jc.passed));
    summary.push(("Liu       SOT polarity (Ta vs Pt)    ", l_pol));
    summary.push(("Liu       J_c(t_CoFeB) shape         ", l_thk.passed));

    // -------------------------------------------------------------------------
    // 6. Garello 2013 (NEW in v0.9.0)
    // -------------------------------------------------------------------------
    println!("\n--- Garello 2013 (Nat. Nanotechnol. 8, 587) ---\n");
    let garello = Garello2013Validation::new()?;
    let g_lin = garello.validate_damping_like_field_linearity(tolerance)?;
    let g_mag = garello.validate_damping_like_magnitude(0.50)?;
    let g_ratio = garello.validate_fl_dl_ratio(tolerance)?;
    let g_jc = garello.validate_critical_switching_current(0.50)?;
    println!(
        "  H_DL linearity in J:      max_rel={:.2}% ⇒ {}",
        g_lin.max_relative_error * 100.0,
        if g_lin.passed { "PASS" } else { "fail" }
    );
    println!(
        "  H_DL magnitude (per J):   max_rel={:.2}% ⇒ {}",
        g_mag.max_relative_error * 100.0,
        if g_mag.passed { "PASS" } else { "fail" }
    );
    println!(
        "  H_FL : H_DL ratio:        max_rel={:.2}% ⇒ {}",
        g_ratio.max_relative_error * 100.0,
        if g_ratio.passed { "PASS" } else { "fail" }
    );
    println!(
        "  Critical J_c (Pt/Co):     max_rel={:.2}% ⇒ {}",
        g_jc.max_relative_error * 100.0,
        if g_jc.passed { "PASS" } else { "fail" }
    );
    summary.push(("Garello   H_DL linearity             ", g_lin.passed));
    summary.push(("Garello   H_DL magnitude             ", g_mag.passed));
    summary.push(("Garello   H_FL : H_DL ratio          ", g_ratio.passed));
    summary.push(("Garello   critical J_c (Pt/Co)       ", g_jc.passed));

    // -------------------------------------------------------------------------
    // 7. Boona 2014 (NEW in v0.9.0)
    // -------------------------------------------------------------------------
    println!("\n--- Boona 2014 (MRS Bulletin 39, 426) ---\n");
    let boona = Boona2014Validation::new()?;
    let b_lin = boona.validate_linear_thermal_response(tolerance)?;
    let b_win = boona.validate_granular_seebeck_window(tolerance)?;
    let b_pol = boona.validate_polarity()?;
    let b_att = boona.validate_grain_boundary_attenuation(tolerance)?;
    println!(
        "  LSSE linear V(ΔT)≤50K:    max_rel={:.2}% ⇒ {}",
        b_lin.max_relative_error * 100.0,
        if b_lin.passed { "PASS" } else { "fail" }
    );
    println!(
        "  Granular S_S window:      max_rel={:.2}% ⇒ {}",
        b_win.max_relative_error * 100.0,
        if b_win.passed { "PASS" } else { "fail" }
    );
    println!(
        "  LSSE polarity:            {}",
        if b_pol { "PASS" } else { "fail" }
    );
    println!(
        "  Grain-boundary atten.:    max_rel={:.2}% ⇒ {}",
        b_att.max_relative_error * 100.0,
        if b_att.passed { "PASS" } else { "fail" }
    );
    summary.push(("Boona     LSSE linear (≤ 50 K)       ", b_lin.passed));
    summary.push(("Boona     granular S_S window        ", b_win.passed));
    summary.push(("Boona     LSSE polarity              ", b_pol));
    summary.push(("Boona     grain-boundary attenuation ", b_att.passed));

    // -------------------------------------------------------------------------
    // Summary dashboard
    // -------------------------------------------------------------------------
    println!("\n=============================================================");
    println!("  Validation Summary (7 papers, v0.7.0–v0.9.0)");
    println!("=============================================================");
    println!("  {:<40}  Status", "Test");
    println!("  {}", "-".repeat(48));
    let mut n_pass = 0_usize;
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
    println!("  Done. 7 landmark papers, {n_pass}/{total} quantitative checks pass.");
    println!("=============================================================\n");

    Ok(())
}
