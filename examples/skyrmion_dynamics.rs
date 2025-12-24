//! Skyrmion Creation and Dynamics
//!
//! **Difficulty**: ⭐⭐ Intermediate
//! **Category**: Topological Textures
//! **Physics**: Magnetic skyrmions, DMI, spin-orbit torque, topological charge
//!
//! This example demonstrates:
//! - Creating a Néel-type skyrmion with DMI
//! - Calculating topological charge
//! - Current-driven skyrmion motion via spin-orbit torque
//! - Skyrmion Hall effect
//!
//! ## Physics:
//! Skyrmions are topologically protected spin textures stabilized by
//! Dzyaloshinskii-Moriya interaction (DMI). They can be moved by
//! electric currents with extremely low threshold densities.
//!
//! Run with: `cargo run --example skyrmion_dynamics`

use std::f64::consts::PI;

use spintronics::prelude::*;

fn main() {
    println!("=== Skyrmion Dynamics Simulation ===\n");

    // Material system: Pt/CoFeB (interfacial DMI)
    let dmi = DmiParameters::pt_cofeb();
    let cofeb = Ferromagnet::cofeb();
    let sot = SpinOrbitTorque::platinum_cofeb();

    println!("Material: Pt/CoFeB");
    println!("  DMI constant: {:.2} mJ/m²", dmi.d * 1000.0);
    println!("  Exchange: {:.2} pJ/m", cofeb.exchange_a * 1e12);
    println!("  Anisotropy: {:.2} kJ/m³", cofeb.anisotropy_k * 1e-3);
    println!(
        "  Spin Hall angle: {:.3}\n",
        sot.effective_spin_hall_efficiency()
    );

    // Calculate skyrmion diameter
    let sk_diameter = dmi.skyrmion_diameter(cofeb.exchange_a);
    println!("Skyrmion Properties:");
    println!("  Diameter: {:.1} nm", sk_diameter);

    // Create skyrmion
    let sk_radius = sk_diameter / 2.0 * 1e-9; // Convert to meters
    let center_nm = (sk_diameter, sk_diameter); // Center in nm
    let center_m = (center_nm.0 * 1e-9, center_nm.1 * 1e-9); // Convert to meters

    println!(
        "  System size: {:.0} × {:.0} nm²",
        sk_diameter * 4.0,
        sk_diameter * 4.0
    );
    println!(
        "  Skyrmion center: ({:.1}, {:.1}) nm\n",
        center_nm.0, center_nm.1
    );

    // === 1. Create Skyrmion ===
    println!("=== Creating Skyrmion ===");

    let skyrmion = Skyrmion::new(
        center_m,
        sk_radius,
        Helicity::Neel,
        Chirality::CounterClockwise,
    );

    println!("  Created Néel-type skyrmion");
    println!("  Center: ({:.1}, {:.1}) nm", center_nm.0, center_nm.1);
    println!("  Radius: {:.1} nm", sk_radius * 1e9);

    // Topological charge
    let q_total = skyrmion.topological_charge;
    println!("  Topological charge Q: {}", q_total);
    println!("  (Q = -1 indicates one skyrmion)\n");

    // Sample magnetization at a few points
    let wall_width = sk_radius / 5.0; // Domain wall width
    let m_center = skyrmion.magnetization_at(center_m.0, center_m.1, wall_width);
    let m_edge = skyrmion.magnetization_at(center_m.0 + sk_radius, center_m.1, wall_width);

    println!(
        "  Magnetization at center: ({:.2}, {:.2}, {:.2})",
        m_center.x, m_center.y, m_center.z
    );
    println!(
        "  Magnetization at edge: ({:.2}, {:.2}, {:.2})",
        m_edge.x, m_edge.y, m_edge.z
    );

    // === 2. Energetics ===
    println!("\n=== Skyrmion Energetics ===");

    // DMI energy favors Néel texture
    let lattice_constant = 0.5; // nm, typical for atomic lattice
    let e_dmi_per_spin = dmi.d * lattice_constant * 1e-9; // J
    println!("  DMI energy scale: {:.2e} J per spin", e_dmi_per_spin);

    // Thermal stability
    let t_room = 300.0; // K
    let kb = 1.380649e-23; // J/K
    let e_thermal = kb * t_room;
    let stability_ratio = (e_dmi_per_spin / e_thermal).abs();
    println!("  E_DMI / k_B T: {:.1}", stability_ratio);

    if stability_ratio > 10.0 {
        println!("  ✓ Thermally stable at room temperature");
    } else {
        println!("  ✗ Thermally unstable");
    }

    // === 3. Current-Driven Motion ===
    println!("\n=== Current-Driven Skyrmion Motion ===");

    let j_current = 1.0e11; // A/m² (10 MA/cm²)
    let current_direction = Vector3::new(1.0, 0.0, 0.0); // Along x

    println!("  Applied current: {:.0} MA/cm²", j_current * 1e-10);
    println!("  Direction: x-axis");

    // Spin-orbit torque at skyrmion edge
    let m_sot = Vector3::new(0.7, 0.0, 0.7).normalize(); // Tilted spin
    let t_dl = sot.damping_like_field(j_current, m_sot, current_direction, cofeb.ms);

    println!("  Damping-like field: {:.2e} T", t_dl.magnitude());

    // Skyrmion velocity (simplified)
    // v = (j_e / e n_s) × (θ_SH / α)
    let e = 1.602e-19; // C
    let n_s = 1.0 / (lattice_constant * 1e-9).powi(2); // spins/m²
    let v_drift = (j_current / (e * n_s)) * (sot.effective_spin_hall_efficiency() / cofeb.alpha);

    println!("  Drift velocity: {:.1} m/s", v_drift);
    println!("  Drift velocity: {:.0} nm/ns\n", v_drift * 1.0);

    // === 4. Skyrmion Hall Effect ===
    println!("=== Skyrmion Hall Effect ===");

    // Skyrmions move at an angle to the current direction
    let _magnus_force_direction = Vector3::new(0.0, 1.0, 0.0); // Perpendicular
    let skyrmion_hall_angle = 5.0 * PI / 180.0; // Typical: 5-20 degrees

    println!("  Magnus force: perpendicular to current");
    println!(
        "  Skyrmion Hall angle: {:.1}°",
        skyrmion_hall_angle * 180.0 / PI
    );
    println!(
        "  Net motion: {:.1}° from current direction\n",
        skyrmion_hall_angle * 180.0 / PI
    );

    // Velocity components
    let v_longitudinal = v_drift * skyrmion_hall_angle.cos();
    let v_transverse = v_drift * skyrmion_hall_angle.sin();

    println!("  Longitudinal velocity: {:.1} m/s", v_longitudinal);
    println!("  Transverse velocity: {:.2} m/s", v_transverse);

    // === 5. Application Potential ===
    println!("\n=== Application Potential ===");

    // Racetrack memory: skyrmion as bit
    let bit_separation = 100.0; // nm
    let write_time = (bit_separation * 1e-9) / v_drift;
    let write_frequency = 1.0 / write_time;

    println!("Racetrack Memory:");
    println!("  Bit separation: {:.0} nm", bit_separation);
    println!("  Write time: {:.2} ns", write_time * 1e9);
    println!("  Write frequency: {:.1} GHz", write_frequency * 1e-9);

    // Energy per bit operation
    let voltage = 1.0; // V
    let resistance = 1000.0; // Ω
    let power = voltage * voltage / resistance;
    let energy_per_bit = power * write_time;

    println!("  Energy per bit: {:.2} fJ", energy_per_bit * 1e15);

    // === 6. Comparison with Domain Walls ===
    println!("\n=== Advantages over Domain Walls ===");
    println!("  ✓ Topological protection (robust)");
    println!("  ✓ Particle-like behavior (well-defined position)");
    println!("  ✓ Low depinning current density");
    println!("  ✗ Skyrmion Hall effect (needs compensation)");

    println!("\n=== Summary ===");
    println!("Skyrmions are promising for:");
    println!("  • Racetrack memory (high density, low power)");
    println!("  • Logic devices (skyrmion-based gates)");
    println!("  • Neuromorphic computing (synaptic weights)");
    println!("\nKey challenge: Suppress skyrmion Hall effect for straight-line motion");
}
