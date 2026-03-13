//! Orbital Hall Effect in Light Metals
//!
//! **Difficulty**: ⭐⭐ Intermediate
//! **Category**: Orbitronics
//! **Physics**: Orbital Hall effect, orbital-to-spin conversion, orbital torques
//!
//! This example demonstrates that the Orbital Hall Effect (OHE) can be giant
//! in light 3d transition metals (Cr, Ti) despite weak spin-orbit coupling,
//! exceeding the Spin Hall Effect (SHE) of heavy metals like Pt.
//!
//! Key insight: OHE originates from orbital texture in k-space, not SOC.
//!
//! References:
//! - Kontani et al., Phys. Rev. Lett. 102, 016601 (2009)
//! - Jo et al., Phys. Rev. B 98, 214405 (2018)

use spintronics::orbitronics::orbital_hall::{
    cr_pt_system, material_ohe_ranking, sigma_oh_to_si, ti_pt_system,
};
use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Orbital Hall Effect: Light Metals vs Heavy Metals ===\n");

    // 1. Material comparison: OHE conductivities
    println!("=== Orbital Hall Conductivity Ranking ===");
    let ranking = material_ohe_ranking();
    println!(
        "{:<6} {:>15} {:>10} {:>12}",
        "Metal", "sigma_OH", "SOC", "OHE/Pt-SHE"
    );
    println!(
        "{:<6} {:>15} {:>10} {:>12}",
        "", "(hbar/e)(Ohm*cm)^-1", "(rel.)", "ratio"
    );
    println!("{}", "-".repeat(50));
    for (name, sigma_oh) in &ranking {
        let ratio = sigma_oh / 2000.0; // Pt SHE reference
        let soc = match *name {
            "Cr" => 0.04,
            "Ti" => 0.02,
            "V" => 0.03,
            "Cu" => 0.01,
            "Al" => 0.005,
            "Pt" => 1.0,
            "W" => 0.8,
            "Mn" => 0.05,
            _ => 0.0,
        };
        println!(
            "{:<6} {:>15.0} {:>10.3} {:>12.2}",
            name, sigma_oh, soc, ratio
        );
    }

    println!("\n  Key: Cr OHE = 5000 >> Pt SHE = 2000, yet Cr SOC << Pt SOC");
    println!("  This proves OHE does not require strong spin-orbit coupling!\n");

    // 2. Compare Cr/Pt and Ti/Pt OHE systems
    println!("=== OHE Systems: Cr/Pt vs Ti/Pt ===");
    let cr_system = cr_pt_system();
    let ti_system = ti_pt_system();

    let e_field = Vector3::new(1.0e5, 0.0, 0.0); // 100 kV/m along x

    let j_orbital_cr = cr_system.orbital_current_density(e_field);
    let j_orbital_ti = ti_system.orbital_current_density(e_field);
    let j_spin_cr = cr_system.effective_spin_current(e_field);
    let j_spin_ti = ti_system.effective_spin_current(e_field);

    println!("  Electric field: E = {:.0e} V/m along x\n", e_field.x);

    println!(
        "  {:<12} {:>18} {:>18} {:>12}",
        "System", "|J_orbital|", "|J_spin_eff|", "Enhancement"
    );
    println!("  {}", "-".repeat(65));
    println!(
        "  {:<12} {:>18.4e} {:>18.4e} {:>12.2}x",
        "Cr/Pt",
        j_orbital_cr.magnitude(),
        j_spin_cr.magnitude(),
        cr_system.enhancement_over_pt_she()
    );
    println!(
        "  {:<12} {:>18.4e} {:>18.4e} {:>12.2}x",
        "Ti/Pt",
        j_orbital_ti.magnitude(),
        j_spin_ti.magnitude(),
        ti_system.enhancement_over_pt_she()
    );

    // 3. Giant OHE despite weak SOC
    println!("\n=== Giant OHE in Light Metals ===");
    let cr = OrbitalHallMaterial::chromium();
    let ti = OrbitalHallMaterial::titanium();
    let pt = OrbitalHallMaterial::platinum();

    println!(
        "  Cr: sigma_OH = {:.0}, SOC = {:.3} (light metal: {})",
        cr.sigma_oh,
        cr.soc_strength,
        cr.is_light_metal()
    );
    println!(
        "  Ti: sigma_OH = {:.0}, SOC = {:.3} (light metal: {})",
        ti.sigma_oh,
        ti.soc_strength,
        ti.is_light_metal()
    );
    println!(
        "  Pt: sigma_OH = {:.0}, SOC = {:.3} (light metal: {})",
        pt.sigma_oh,
        pt.soc_strength,
        pt.is_light_metal()
    );
    println!(
        "\n  Cr OHE / Pt SHE ratio = {:.2}",
        cr.ohe_to_pt_she_ratio()
    );
    println!("  Ti OHE / Pt SHE ratio = {:.2}", ti.ohe_to_pt_she_ratio());

    // 4. Orbital torques comparison
    println!("\n=== Orbital Torques: Cr/Pt/CoFeB ===");
    let ot_cr = OrbitalTorque::cr_pt_cofeb();
    let ot_ti = OrbitalTorque::ti_pt_cofeb();

    let m_perp = Vector3::new(0.0, 0.0, 1.0); // perpendicular magnetization
    let j_charge = 1.0e11; // 10^11 A/m^2
    let current_dir = Vector3::new(1.0, 0.0, 0.0);

    let (tau_dl_cr, tau_fl_cr) = ot_cr.compute_torques(m_perp, j_charge, current_dir);
    let (tau_dl_ti, tau_fl_ti) = ot_ti.compute_torques(m_perp, j_charge, current_dir);

    println!("  J_c = {:.0e} A/m^2, m = z (perpendicular)\n", j_charge);
    println!(
        "  {:<12} {:>15} {:>15} {:>10}",
        "System", "|tau_DL|", "|tau_FL|", "OT/SOT"
    );
    println!("  {}", "-".repeat(57));
    let ratio_cr = ot_cr.ot_to_sot_ratio(j_charge, m_perp, current_dir);
    let ratio_ti = ot_ti.ot_to_sot_ratio(j_charge, m_perp, current_dir);
    println!(
        "  {:<12} {:>15.4e} {:>15.4e} {:>10.3}",
        "Cr/Pt/CoFeB",
        tau_dl_cr.magnitude(),
        tau_fl_cr.magnitude(),
        ratio_cr
    );
    println!(
        "  {:<12} {:>15.4e} {:>15.4e} {:>10.3}",
        "Ti/Pt/CoFeB",
        tau_dl_ti.magnitude(),
        tau_fl_ti.magnitude(),
        ratio_ti
    );

    // 5. Orbital accumulation at interfaces
    println!("\n=== Orbital Accumulation at Interfaces ===");
    let e_mag = 1.0e5; // V/m
    let lambda_l = 10.0e-9; // 10 nm orbital diffusion length

    let acc_cr = cr_system.orbital_accumulation(e_mag, lambda_l);
    let acc_ti = ti_system.orbital_accumulation(e_mag, lambda_l);

    println!(
        "  E = {:.0e} V/m, lambda_L = {:.0} nm",
        e_mag,
        lambda_l * 1e9
    );
    println!("  Cr/Pt accumulation: {:.4e} hbar", acc_cr);
    println!("  Ti/Pt accumulation: {:.4e} hbar", acc_ti);

    // 6. Unit conversion
    println!("\n=== SI Unit Conversion ===");
    let sigma_si = sigma_oh_to_si(5000.0);
    println!(
        "  Cr: sigma_OH = 5000 (hbar/e)(Ohm*cm)^-1 = {:.4e} (hbar/(2e))*S/m",
        sigma_si
    );

    println!("\n=== Summary ===");
    println!("Orbitronics key findings:");
    println!("  - Light metals (Cr, Ti) have giant OHE despite weak SOC");
    println!(
        "  - Cr OHE ({:.0}) exceeds Pt SHE (2000) by {:.1}x",
        cr.sigma_oh,
        cr.ohe_to_pt_she_ratio()
    );
    println!(
        "  - Orbital-to-spin conversion at Cr/Pt interface: eta = {:.2}",
        cr_system.converter.conversion_efficiency
    );
    println!("  - Orbital torques can exceed conventional SOT for optimized stacks");

    Ok(())
}
