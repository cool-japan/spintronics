//! Barnett Effect in Liquid Metals
//!
//! Demonstrates spin current generation via spin-vorticity coupling in liquid mercury.
//!
//! Based on: Takahashi et al., Nature Physics 12, 52-56 (2016)
//! "Spin hydrodynamic generation"

use spintronics::fluid::{BarnettField, NavierStokes};
use spintronics::vector3::Vector3;

fn main() {
    println!("=== Barnett Effect in Liquid Mercury ===\n");

    // Setup mercury flow with vorticity
    let dt = 1.0e-6; // 1 microsecond time step
    let mercury = NavierStokes::mercury(dt);
    let barnett = BarnettField::default();

    println!("Liquid: Mercury (Hg)");
    println!("  Density: {:.0} kg/m³", mercury.density);
    println!("  Kinematic viscosity: {:.2e} m²/s", mercury.viscosity);
    println!();

    // Simulate vortex flow
    println!("Simulating Taylor-Couette flow (rotating cylinders)...\n");

    let inner_radius = 0.01; // 1 cm
    let outer_radius = 0.02; // 2 cm
    let omega_inner = 100.0; // 100 rad/s

    // Calculate velocity profile at mid-radius
    let r = (inner_radius + outer_radius) / 2.0;
    let v_theta = omega_inner * inner_radius * inner_radius / r;

    println!("Inner cylinder rotation: {:.1} rad/s", omega_inner);
    println!("At r = {:.3} m:", r);
    println!("  Tangential velocity: {:.3} m/s", v_theta);

    // For Taylor-Couette flow, vorticity in z-direction
    let omega_z = 2.0 * omega_inner * inner_radius * inner_radius
        / (outer_radius * outer_radius - inner_radius * inner_radius);

    let vorticity = Vector3::new(0.0, 0.0, omega_z);

    println!("  Vorticity: {:.2} rad/s (z-component)", omega_z);
    println!();

    // Calculate Barnett field
    let h_barnett = barnett.field_from_vorticity(vorticity);

    println!("=== Barnett Field ===");
    println!("H_B = -Ω / (γ μ₀)");
    println!("  Hx: {:.3e} A/m", h_barnett.x);
    println!("  Hy: {:.3e} A/m", h_barnett.y);
    println!("  Hz: {:.3e} A/m", h_barnett.z);
    println!("  |H_B|: {:.3e} A/m", h_barnett.magnitude());
    println!();

    // Convert to magnetic field in Tesla
    let mu0 = 1.256637e-6;
    let b_barnett = h_barnett * mu0;

    println!("Effective magnetic field:");
    println!(
        "  |B_B|: {:.3e} T = {:.3} nT",
        b_barnett.magnitude(),
        b_barnett.magnitude() * 1e9
    );
    println!();

    // Compare with typical lab fields
    println!("Context:");
    println!("  Earth's field: ~50 μT");
    println!("  Typical lab magnet: ~1 T");
    println!(
        "  Barnett field (this flow): {:.3} nT",
        b_barnett.magnitude() * 1e9
    );
    println!();

    // Estimate induced spin polarization
    println!("=== Spin Polarization ===");

    let temperature = 300.0; // Room temperature
    let omega = vorticity * 0.5; // Convert vorticity to angular velocity
    let spin_polarization = barnett.spin_polarization(omega, temperature);
    println!("P = {:.3e}", spin_polarization);
    println!();

    // Higher rotation rate example
    println!("=== Higher Rotation Rate ===");
    let omega_fast = 10000.0; // 10,000 rad/s
    let vorticity_fast = Vector3::new(0.0, 0.0, omega_fast);

    let h_fast = barnett.field_from_vorticity(vorticity_fast);
    let b_fast = h_fast * mu0;

    println!("Ω = {} rad/s", omega_fast);
    println!("  Barnett field: {:.3} μT", b_fast.magnitude() * 1e6);
    let omega_fast_vec = Vector3::new(0.0, 0.0, omega_fast);
    println!(
        "  Spin polarization: {:.3e}",
        barnett.spin_polarization(omega_fast_vec, temperature)
    );
    println!();

    // Reynolds number analysis
    println!("=== Flow Regime ===");
    let re = mercury.reynolds_number(v_theta, outer_radius - inner_radius);
    println!("Reynolds number: {:.1}", re);
    if re < 2300.0 {
        println!("  → Laminar flow (Re < 2300)");
    } else {
        println!("  → Turbulent flow (Re > 2300)");
    }
    println!();

    println!("=== Summary ===");
    println!("Fluid vorticity creates effective magnetic field via Barnett effect.");
    println!("This induces spin polarization in the conduction electrons of liquid metal.");
    println!("Experiment: Measure voltage from ISHE in adjacent Pt layer.");
}
