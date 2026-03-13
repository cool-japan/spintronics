//! Magnetoelastic Switching and Straintronics
//!
//! **Difficulty**: ⭐⭐ Intermediate
//! **Category**: Straintronics / Multiferroics
//! **Physics**: Magnetostriction, Villari effect, strain-mediated switching, energy comparison
//!
//! This example demonstrates magnetoelastic coupling and strain-mediated
//! magnetization control, including:
//!
//! 1. Magnetostrictive material properties (Terfenol-D, Galfenol, Nickel)
//! 2. Magnetostrictive strain tensor computation
//! 3. Magnetoelastic energy under applied strain
//! 4. Villari effect (stress-induced effective field)
//! 5. Strain vs SOT vs STT switching energy comparison
//! 6. Piezoelectric substrate and critical switching voltage
//!
//! References:
//! - Roy et al., "Hybrid spintronics and straintronics" (2011)
//! - Biswas et al., Phys. Rev. B (2014)

use spintronics::mech::magnetoelastic::{
    isotropic_magnetostrictive_strain, magnetoelastic_energy_density_cubic,
    magnetostrictive_strain_tensor, stress_induced_anisotropy, villari_effective_field_magnitude,
    MagnetoelasticMaterial, PiezoelectricSubstrate, StrainTensor,
};
use spintronics::mech::straintronics::{
    critical_strain_uniaxial, critical_voltage, strain_switching_energy,
    typical_switching_energy_comparison, StraintronicDevice,
};
use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Magnetoelastic Switching and Straintronics ===\n");

    // 1. Material comparison
    println!("=== Magnetostrictive Materials ===");
    let materials = [
        MagnetoelasticMaterial::terfenol_d(),
        MagnetoelasticMaterial::galfenol(),
        MagnetoelasticMaterial::nickel(),
        MagnetoelasticMaterial::cofeb(),
        MagnetoelasticMaterial::iron(),
    ];

    println!(
        "  {:<12} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "Material", "l_100 (ppm)", "l_111 (ppm)", "l_s (ppm)", "Y (GPa)", "M_s (kA/m)"
    );
    println!("  {}", "-".repeat(68));
    for mat in &materials {
        println!(
            "  {:<12} {:>10.0} {:>10.0} {:>10.1} {:>10.0} {:>10.0}",
            mat.name,
            mat.lambda_100 * 1e6,
            mat.lambda_111 * 1e6,
            mat.lambda_s * 1e6,
            mat.youngs_modulus * 1e-9,
            mat.ms * 1e-3
        );
    }

    // 2. Magnetostrictive strain tensor
    println!("\n=== Magnetostrictive Strain (Terfenol-D) ===");
    let terfenol = MagnetoelasticMaterial::terfenol_d();
    let directions = [
        (Vector3::new(1.0, 0.0, 0.0), "[100]"),
        (Vector3::new(0.0, 0.0, 1.0), "[001]"),
        (Vector3::new(1.0, 1.0, 1.0), "[111]"),
    ];

    for (dir, label) in &directions {
        let strain = magnetostrictive_strain_tensor(&terfenol, dir)
            .expect("Failed to compute strain tensor");
        let (exx, eyy, ezz) = strain.diagonal();
        println!(
            "  M along {}: e_xx={:+.1} ppm, e_yy={:+.1} ppm, e_zz={:+.1} ppm",
            label,
            exx * 1e6,
            eyy * 1e6,
            ezz * 1e6
        );
    }

    // 3. Isotropic magnetostrictive strain vs angle
    println!("\n=== Isotropic Strain vs Angle ===");
    println!("  {:>12} {:>15}", "theta (deg)", "strain (ppm)");
    println!("  {}", "-".repeat(30));
    for &theta_deg in &[0.0, 30.0, 45.0, 54.74, 60.0, 90.0] {
        let cos_theta = (theta_deg * std::f64::consts::PI / 180.0).cos();
        let strain = isotropic_magnetostrictive_strain(terfenol.lambda_s, cos_theta);
        println!("  {:>12.2} {:>15.2}", theta_deg, strain * 1e6);
    }
    println!("  (54.74 deg = magic angle where strain vanishes)");

    // 4. Magnetoelastic energy
    println!("\n=== Magnetoelastic Energy Under Applied Strain ===");
    let strain_applied =
        StrainTensor::uniaxial(500.0e-6, 0).expect("Failed to create strain tensor");
    println!(
        "  Applied strain: e_xx = {:.0} ppm (uniaxial along x)\n",
        500.0
    );

    println!(
        "  {:>12} {:>20} {:>15}",
        "M direction", "E_me density (J/m^3)", "Delta E (kJ/m^3)"
    );
    println!("  {}", "-".repeat(50));

    let m_x = Vector3::new(1.0, 0.0, 0.0);
    let m_z = Vector3::new(0.0, 0.0, 1.0);
    let e_along = magnetoelastic_energy_density_cubic(&terfenol, &strain_applied, &m_x)
        .expect("Failed to compute ME energy");
    let e_perp = magnetoelastic_energy_density_cubic(&terfenol, &strain_applied, &m_z)
        .expect("Failed to compute ME energy");

    println!("  {:>12} {:>20.2} {:>15}", "[100] (||)", e_along, "-");
    println!(
        "  {:>12} {:>20.2} {:>15.2}",
        "[001] (perp)",
        e_perp,
        (e_perp - e_along) * 1e-3
    );

    // 5. Villari effect
    println!("\n=== Villari Effect (Stress -> Effective Field) ===");
    println!(
        "  {:>12} {:>15} {:>15}",
        "Stress (MPa)", "H_eff (A/m)", "H_eff (Oe)"
    );
    println!("  {}", "-".repeat(47));
    for &stress_mpa in &[1.0, 10.0, 50.0, 100.0, 500.0] {
        let stress = stress_mpa * 1e6;
        let h_eff = villari_effective_field_magnitude(&terfenol, stress)
            .expect("Failed to compute Villari field");
        let h_oe = h_eff / 79.577; // 1 Oe = 79.577 A/m
        println!("  {:>12.0} {:>15.0} {:>15.0}", stress_mpa, h_eff, h_oe);
    }

    // 6. Stress-induced anisotropy
    println!("\n=== Stress-Induced Anisotropy ===");
    let stress_50mpa = 50.0e6;
    let k_sigma = stress_induced_anisotropy(terfenol.lambda_s, stress_50mpa);
    println!(
        "  Stress = 50 MPa, lambda_s = {:.0} ppm",
        terfenol.lambda_s * 1e6
    );
    println!(
        "  K_sigma = -(3/2) lambda_s * sigma = {:.2} kJ/m^3",
        k_sigma * 1e-3
    );

    // 7. Switching energy comparison
    println!("\n=== Switching Energy Comparison ===");
    let (e_strain, e_sot, e_stt) = typical_switching_energy_comparison();
    println!(
        "  {:<20} {:>15} {:>12}",
        "Mechanism", "Energy (J)", "Relative"
    );
    println!("  {}", "-".repeat(50));
    println!(
        "  {:<20} {:>15.2e} {:>12.1}x",
        "Strain switching", e_strain, 1.0
    );
    println!(
        "  {:<20} {:>15.2e} {:>12.0}x",
        "SOT switching",
        e_sot,
        e_sot / e_strain
    );
    println!(
        "  {:<20} {:>15.2e} {:>12.0}x",
        "STT switching",
        e_stt,
        e_stt / e_strain
    );
    println!(
        "  -> Strain switching is {:.0}x more efficient than SOT!",
        e_sot / e_strain
    );

    // 8. Straintronic device
    println!("\n=== Straintronic Device (Terfenol-D / PMN-PT) ===");
    let substrate = PiezoelectricSubstrate::pmn_pt();
    let ku = 5.0e4; // 50 kJ/m^3 uniaxial anisotropy
    let device = StraintronicDevice::new(terfenol, substrate, 2.0e-9, 100.0e-9, ku)
        .expect("Failed to create straintronic device");

    println!(
        "  Free layer: {} ({:.0} nm x {:.0} nm x {:.1} nm)",
        terfenol.name,
        device.free_layer_width * 1e9,
        device.free_layer_width * 1e9,
        device.free_layer_thickness * 1e9
    );
    println!("  Substrate: {}", substrate.name);
    println!("  K_u = {:.0} kJ/m^3", ku * 1e-3);
    println!("  Volume = {:.2e} m^3", device.free_layer_volume());

    // Critical strain and voltage
    let eps_crit = critical_strain_uniaxial(ku, terfenol.lambda_s, terfenol.youngs_modulus)
        .expect("Failed to compute critical strain");
    println!("\n  Critical strain: {:.2} ppm", eps_crit * 1e6);

    let substrate_thickness = 500.0e-6; // 500 um
    let v_crit = critical_voltage(eps_crit, &substrate, substrate_thickness, 0.9)
        .expect("Failed to compute critical voltage");
    println!(
        "  Critical voltage: {:.2} mV (at eta=0.9, d={:.0} um substrate)",
        v_crit * 1e3,
        substrate_thickness * 1e6
    );

    // 9. Strain switching energy for device
    let m_initial = Vector3::new(1.0, 0.0, 0.0); // along x
    let m_final = Vector3::new(0.0, 1.0, 0.0); // 90 deg switch to y
    let strain_switch = StrainTensor::uniaxial(eps_crit, 0).expect("Failed to create strain");
    let e_switch = strain_switching_energy(&device, &strain_switch, &m_initial, &m_final)
        .expect("Failed to compute switching energy");
    println!(
        "  90-deg switching energy: {:.4e} J ({:.2} aJ)",
        e_switch,
        e_switch * 1e18
    );

    println!("\n=== Summary ===");
    println!("Magnetoelastic switching key results:");
    println!(
        "  - Terfenol-D: giant lambda_s = {:.0} ppm",
        terfenol.lambda_s * 1e6
    );
    println!(
        "  - Strain switching {:.0}x more energy-efficient than SOT",
        e_sot / e_strain
    );
    println!(
        "  - Critical voltage: {:.2} mV for PMN-PT substrate",
        v_crit * 1e3
    );
    println!("  - Ultra-low switching energy: {:.2} aJ", e_switch * 1e18);

    Ok(())
}
