//! SimulationBuilder: Type-State Pattern for Spintronics Simulations
//!
//! **Difficulty**: ⭐ Beginner
//! **Category**: Simulation Framework
//! **Physics**: LLG dynamics, Zeeman energy, ferromagnetic resonance
//!
//! This example demonstrates the type-state SimulationBuilder pattern
//! for constructing and running micromagnetic simulations. The builder
//! ensures at compile time that all required parameters (material,
//! external field, solver) are set before a simulation can be built.
//!
//! We show:
//! 1. Building a simulation with the type-state pattern
//! 2. Running LLG dynamics with RK4 integration
//! 3. Magnetization precession under external field
//! 4. Energy evolution during relaxation
//! 5. Comparing different materials and field configurations
//!
//! Reference: Gilbert, IEEE Trans. Magn. 40, 3443 (2004)

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== SimulationBuilder: Type-State Pattern ===\n");

    // 1. Basic simulation: YIG in external field
    println!("=== Simulation 1: YIG Magnetization Precession ===");
    let yig = Ferromagnet::yig();
    println!("  Material: YIG");
    println!("  M_s = {:.3e} A/m", yig.ms);
    println!("  alpha = {:.1e} (Gilbert damping)", yig.alpha);

    let h_ext = Vector3::new(0.0, 0.0, 0.1); // 0.1 T along z
    println!(
        "  External field: ({:.2}, {:.2}, {:.2}) T",
        h_ext.x, h_ext.y, h_ext.z
    );

    let mut sim = SimulationBuilder::new()
        .material(yig)
        .external_field(h_ext)
        .solver_rk4()
        .num_steps(1000)
        .build()
        .expect("Failed to build YIG simulation");

    let result = sim.run().expect("Failed to run YIG simulation");

    println!(
        "  Initial m: ({:.4}, {:.4}, {:.4})",
        result.trajectory[0].x, result.trajectory[0].y, result.trajectory[0].z
    );
    println!(
        "  Final m:   ({:.4}, {:.4}, {:.4})",
        result.final_magnetization.x, result.final_magnetization.y, result.final_magnetization.z
    );
    println!("  Steps: {}", result.trajectory.len() - 1);
    println!(
        "  |m_final| = {:.8}",
        result.final_magnetization.magnitude()
    );

    // Energy evolution
    if !result.energies.is_empty() {
        let e_initial = result.energies[0];
        let e_final = result.energies[result.energies.len() - 1];
        println!("  Initial energy density: {:.4e} J/m^3", e_initial);
        println!("  Final energy density:   {:.4e} J/m^3", e_final);
        println!("  Energy change: {:.4e} J/m^3", e_final - e_initial);
    }

    // 2. Permalloy with tilted initial magnetization
    println!("\n=== Simulation 2: Permalloy Relaxation ===");
    let py = Ferromagnet::permalloy();
    println!("  Material: Permalloy (Py)");
    println!("  M_s = {:.3e} A/m", py.ms);
    println!(
        "  alpha = {:.3} (higher damping -> faster relaxation)",
        py.alpha
    );

    let h_z = Vector3::new(0.0, 0.0, 0.5); // 0.5 T along z

    let mut sim_py = SimulationBuilder::new()
        .material(py)
        .external_field(h_z)
        .solver_rk4()
        .num_steps(2000)
        .build()
        .expect("Failed to build Permalloy simulation");

    let result_py = sim_py.run().expect("Failed to run Permalloy simulation");

    println!("  External field: {:.2} T along z", h_z.z);
    println!(
        "  Final m: ({:.6}, {:.6}, {:.6})",
        result_py.final_magnetization.x,
        result_py.final_magnetization.y,
        result_py.final_magnetization.z
    );
    println!("  m_z alignment: {:.6}", result_py.final_magnetization.z);

    // 3. Cobalt: strong magnetization
    println!("\n=== Simulation 3: Cobalt in Strong Field ===");
    let co = Ferromagnet::cobalt();
    println!("  Material: Cobalt (Co)");
    println!("  M_s = {:.3e} A/m", co.ms);

    let h_strong = Vector3::new(0.0, 0.0, 1.0); // 1.0 T along z

    let mut sim_co = SimulationBuilder::new()
        .material(co)
        .external_field(h_strong)
        .solver_rk4()
        .num_steps(500)
        .build()
        .expect("Failed to build Cobalt simulation");

    let result_co = sim_co.run().expect("Failed to run Cobalt simulation");

    println!("  H_ext = 1.0 T along z");
    println!(
        "  Final m: ({:.6}, {:.6}, {:.6})",
        result_co.final_magnetization.x,
        result_co.final_magnetization.y,
        result_co.final_magnetization.z
    );

    // 4. Material comparison table
    println!("\n=== Material Comparison Under Same Field ===");
    let field = Vector3::new(0.0, 0.0, 0.2); // 0.2 T
    let step_count = 1000;

    let ferromagnets = [
        ("YIG", Ferromagnet::yig()),
        ("Permalloy", Ferromagnet::permalloy()),
        ("Cobalt", Ferromagnet::cobalt()),
        ("Iron", Ferromagnet::iron()),
        ("Nickel", Ferromagnet::nickel()),
    ];

    println!("  H_ext = {:.2} T along z, {} steps\n", field.z, step_count);
    println!(
        "  {:<12} {:>10} {:>10} {:>12} {:>12}",
        "Material", "M_s (kA/m)", "alpha", "m_z final", "|dm|"
    );
    println!("  {}", "-".repeat(60));

    for (name, mat) in &ferromagnets {
        let mut sim_mat = SimulationBuilder::new()
            .material(mat.clone())
            .external_field(field)
            .solver_rk4()
            .num_steps(step_count)
            .build()
            .expect("Failed to build simulation");

        let res = sim_mat.run().expect("Failed to run simulation");

        let dm = res.final_magnetization - res.trajectory[0];
        println!(
            "  {:<12} {:>10.0} {:>10.4} {:>12.6} {:>12.6}",
            name,
            mat.ms * 1e-3,
            mat.alpha,
            res.final_magnetization.z,
            dm.magnitude()
        );
    }

    // 5. Field scan: FMR-like frequency determination
    println!("\n=== Field Scan: Precession Dynamics ===");
    println!("  Material: YIG");
    println!(
        "  {:>10} {:>15} {:>15}",
        "H_z (T)", "m_z (final)", "|m_final|"
    );
    println!("  {}", "-".repeat(45));

    let yig2 = Ferromagnet::yig();
    for &h_z_val in &[0.01, 0.05, 0.1, 0.5, 1.0, 2.0] {
        let h = Vector3::new(0.0, 0.0, h_z_val);
        let mut sim_scan = SimulationBuilder::new()
            .material(yig2.clone())
            .external_field(h)
            .solver_rk4()
            .num_steps(500)
            .build()
            .expect("Failed to build scan simulation");

        let res_scan = sim_scan.run().expect("Failed to run scan simulation");

        println!(
            "  {:>10.2} {:>15.6} {:>15.8}",
            h_z_val,
            res_scan.final_magnetization.z,
            res_scan.final_magnetization.magnitude()
        );
    }

    println!("\n=== Summary ===");
    println!("SimulationBuilder type-state pattern:");
    println!("  - Compile-time guarantee: material + field + solver all required");
    println!("  - RK4 integration of LLG equation");
    println!("  - Magnetization conserved: |m| = 1 throughout dynamics");
    println!("  - Higher damping (alpha) -> faster alignment with external field");
    println!("  - Energy decreases monotonically during relaxation");

    Ok(())
}
