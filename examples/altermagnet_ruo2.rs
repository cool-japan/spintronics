//! RuO2 Altermagnet Properties
//!
//! **Difficulty**: ⭐⭐ Intermediate
//! **Category**: Altermagnetic spintronics
//! **Physics**: Spin splitting, spin-splitter effect, crystal Hall effect
//!
//! This example demonstrates the properties of RuO2, a prototypical d-wave
//! altermagnet. Altermagnets have zero net magnetization like antiferromagnets,
//! but exhibit non-relativistic spin splitting like ferromagnets. We show:
//!
//! 1. d-wave angular symmetry of spin splitting
//! 2. Spin-splitter effect (spin current generation without SOC)
//! 3. Crystal Hall conductivity
//! 4. Temperature dependence of the order parameter
//!
//! Reference: Smejkal et al., Phys. Rev. X 12, 040501 (2022)

use std::f64::consts::PI;

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== RuO2 Altermagnet Properties ===\n");

    // 1. Create RuO2 altermagnet
    let ruo2 = Altermagnet::ruo2();
    println!("Material: {}", ruo2);
    println!("  Crystal structure: {}", ruo2.crystal_structure);
    println!(
        "  Symmetry: {} (harmonic order {})",
        ruo2.symmetry,
        ruo2.symmetry.harmonic_order()
    );
    println!("  Neel temperature T_N = {:.0} K", ruo2.neel_temperature);
    println!("  Max spin splitting = {:.2} eV", ruo2.spin_splitting);
    println!(
        "  Net magnetization = {:.1e} A/m (zero by definition)\n",
        ruo2.net_magnetization()
    );

    // 2. Angular dependence of spin splitting (d-wave symmetry)
    println!("=== d-wave Spin Splitting vs Angle ===");
    println!(
        "{:>10} {:>15} {:>15}",
        "Angle (deg)", "Splitting (eV)", "cos(2*phi)"
    );
    println!("{}", "-".repeat(45));

    let n_angles = 13;
    for i in 0..n_angles {
        let phi_deg = i as f64 * 360.0 / (n_angles - 1) as f64;
        let phi_rad = phi_deg * PI / 180.0;
        let splitting = ruo2.splitting_at_angle(phi_rad);
        let cos2phi = (2.0 * phi_rad).cos();
        println!("{:10.1} {:15.4} {:15.4}", phi_deg, splitting, cos2phi);
    }

    // Show node angles
    let nodes = ruo2.symmetry.node_angles();
    print!("\n  Nodes ({} total): ", nodes.len());
    for node in &nodes {
        print!("{:.1} deg  ", node.to_degrees());
    }
    println!("\n");

    // 3. Transport: spin-splitter effect
    println!("=== Spin-Splitter Transport ===");
    let conductivity = 1.0e6; // S/m (metallic)
    let scattering_time = 1.0e-14; // s
    let transport = AltermagnetTransport::new(&ruo2, conductivity, scattering_time)?;

    let efficiency = transport.spin_splitter_efficiency();
    println!("  Spin-splitter efficiency theta_AM = {:.4}", efficiency);
    println!("  (Compare: Pt spin Hall angle ~ 0.07)\n");

    // Spin current vs charge current
    println!("  Charge current J_c (A/m^2)   Spin current J_s (A/m^2)");
    println!("  {}", "-".repeat(55));
    for &j_c in &[1.0e9, 5.0e9, 1.0e10, 5.0e10, 1.0e11] {
        let j_s = transport.spin_splitter_current(j_c)?;
        println!("  {:>15.2e}              {:>15.4e}", j_c, j_s);
    }

    // 4. Directional spin current (angle-resolved)
    println!("\n=== Angle-Resolved Spin Current (J_c = 10^10 A/m^2) ===");
    let j_c = 1.0e10;
    println!("{:>10} {:>20}", "Angle (deg)", "J_s (A/m^2)");
    println!("{}", "-".repeat(35));
    for i in 0..9 {
        let phi = i as f64 * PI / 8.0;
        let j_s_dir = transport.directional_spin_current(j_c, phi)?;
        println!("{:10.1} {:20.4e}", phi.to_degrees(), j_s_dir);
    }

    // 5. Crystal Hall conductivity
    println!("\n=== Crystal Hall Effect ===");
    let sigma_h = transport.crystal_hall_conductivity(ruo2.lattice_constant)?;
    println!("  Crystal Hall conductivity: {:.4e} S/m", sigma_h);
    println!("  (Arises from Berry curvature, no SOC needed)");

    // 6. Temperature dependence
    println!("\n=== Temperature Dependence ===");
    println!(
        "{:>10} {:>15} {:>15} {:>15}",
        "T (K)", "Order param", "Splitting (eV)", "Efficiency"
    );
    println!("{}", "-".repeat(60));
    for &temp in &[0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 290.0, 300.0, 350.0] {
        let order = ruo2.order_parameter(temp)?;
        let split = ruo2.splitting_at_temperature(temp)?;
        let eff = transport.spin_splitter_efficiency_at_temperature(temp)?;
        println!("{:10.0} {:15.4} {:15.4} {:15.6}", temp, order, split, eff);
    }

    // 7. Compare with other altermagnets
    println!("\n=== Altermagnet Material Comparison ===");
    let materials = [
        Altermagnet::ruo2(),
        Altermagnet::crsb(),
        Altermagnet::mnte(),
        Altermagnet::fe2o3(),
    ];
    println!(
        "{:<8} {:>10} {:>10} {:>12} {:>8}",
        "Name", "Symmetry", "T_N (K)", "Delta (eV)", "Nodes"
    );
    println!("{}", "-".repeat(55));
    for mat in &materials {
        println!(
            "{:<8} {:>10} {:>10.0} {:>12.2} {:>8}",
            mat.name,
            format!("{}", mat.symmetry),
            mat.neel_temperature,
            mat.spin_splitting,
            mat.symmetry.num_nodes(),
        );
    }

    println!("\n=== Summary ===");
    println!("RuO2 demonstrates key altermagnet properties:");
    println!("  - Zero net magnetization with large (1.4 eV) spin splitting");
    println!("  - d-wave angular symmetry: nodes at 45, 135, 225, 315 degrees");
    println!(
        "  - Spin-splitter efficiency ~ {:.1}% (comparable to heavy metals)",
        efficiency * 100.0
    );
    println!("  - No spin-orbit coupling required for spin current generation");

    Ok(())
}
