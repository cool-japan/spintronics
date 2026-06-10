//! Spin Hall Magnetoresistance — angular scan and USMR
//!
//! Demonstrates the angular dependence of longitudinal and Hall resistivity
//! in Pt/YIG bilayer due to SMR. Also shows USMR linear-in-current behavior.
//!
//! Run with: cargo run --example smr_angular_scan

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Spin Hall Magnetoresistance (SMR) Angular Scan ===\n");

    let smr = SpinHallMagnetoresistance::platinum_yig();
    println!("Pt/YIG SMR parameters:");
    println!("  θ_SH = {:.3}", smr.theta_sh);
    println!("  λ_sf = {:.1} nm", smr.lambda_sf * 1e9);
    println!("  t_NM = {:.1} nm", smr.t_nm * 1e9);
    println!("  ρ₁/ρ₀ (SMR ratio) = {:.2e}", smr.smr_ratio());
    println!("  Expected ≈ 1×10⁻³ for Pt/YIG\n");

    // Angular scan: sweep M in x-z plane (m_y = 0 throughout)
    // In this plane: ρ_L = ρ₀ + ρ₁ (constant), ρ_H = −ρ₂ sin(α)
    println!("=== Longitudinal Resistivity vs Angle (m in x-z plane) ===");
    println!("  α [deg]   ρ_L/ρ₀      ρ_H/ρ₀");
    let n = 8_usize;
    for i in 0..=n {
        let alpha_deg = i as f64 * 360.0 / n as f64;
        let alpha_rad = alpha_deg.to_radians();
        // m in x-z plane: m = (cos α, 0, sin α)
        let m = Vector3::new(alpha_rad.cos(), 0.0, alpha_rad.sin());
        let rho_l = smr.longitudinal_resistivity(m);
        let rho_h = smr.hall_resistivity(m);
        let rho0 = smr.rho_0();
        println!(
            "  {:6.1}     {:.6}    {:.6}",
            alpha_deg,
            rho_l / rho0,
            rho_h / rho0
        );
    }

    // USMR section
    println!("\n=== Unidirectional SMR (USMR) — current-linear term ===");
    let usmr = UnidirectionalSmr::platinum_cobalt()?;
    println!("Pt/Co USMR coefficient: {:.2e} m²/A", usmr.usmr_coefficient);
    let j_vals = [5e10_f64, 1e11, 1.5e11, 2e11];
    // m along +y maximises the USMR signal (ĵ × ẑ) · m = −m_y
    let m_plus = Vector3::new(0.0, 1.0, 0.0);
    let m_minus = Vector3::new(0.0, -1.0, 0.0);
    println!("  J [A/m²]     ΔR/R₀ (+m_y)    ΔR/R₀ (-m_y)   Odd?");
    for &j in &j_vals {
        let dr_plus = usmr.usmr_relative_change(m_plus, j);
        let dr_minus = usmr.usmr_relative_change(m_minus, j);
        let is_odd = (dr_plus + dr_minus).abs() < 1e-15;
        println!(
            "  {:.1e}    {:.2e}        {:.2e}        {}",
            j,
            dr_plus,
            dr_minus,
            if is_odd { "yes" } else { "no" }
        );
    }

    println!("\n=== Nakayama 2013 Validation ===");
    let val = Nakayama2013Validation::new()?;
    let tol = 0.30;
    let r1 = val.validate_longitudinal_angular(tol)?;
    let r2 = val.validate_smr_ratio(tol)?;
    println!("  {}", r1.summary());
    println!("  {}", r2.summary());

    println!("\n=== Avci 2015 USMR Validation ===");
    let avci = Avci2015Validation::new()?;
    let r3 = avci.validate_current_linearity(tol)?;
    let r4 = avci.validate_coefficient_magnitude(tol)?;
    println!("  {}", r3.summary());
    println!("  {}", r4.summary());

    Ok(())
}
