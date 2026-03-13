//! 3D Topological Hopfions
//!
//! **Difficulty**: ⭐⭐⭐ Advanced
//! **Category**: Topological Magnetism
//! **Physics**: Hopf fibration, Hopf invariant, toroidal solitons, exchange energy
//!
//! This example demonstrates 3D magnetic hopfions -- topological solitons
//! characterized by the Hopf invariant Q_H, which counts the linking number
//! of preimage curves on S^2. We show:
//!
//! 1. Hopf fibration map S^3 -> S^2
//! 2. Hopfion construction on a 3D grid with toroidal magnetization texture
//! 3. Hopf invariant calculation (topological charge)
//! 4. Exchange energy and total energy of hopfion configurations
//! 5. Comparison of different Hopf charges and profile functions
//!
//! References:
//! - Rybakov et al., APL Mater. 10, 111113 (2022)
//! - Kent et al., Nat. Commun. 12, 1562 (2021)

use spintronics::texture::hopfion::{
    HopfFibration, HopfInvariant, Hopfion, HopfionEnergy, HopfionEnergyParams, ProfileFunction,
};

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== 3D Topological Hopfions ===\n");

    // 1. Hopf fibration: map points from S^3 to S^2
    println!("=== Hopf Fibration Map ===");
    let test_points: Vec<(f64, f64, f64, f64, &str)> = vec![
        (1.0, 0.0, 0.0, 0.0, "(1,0,0,0)"),
        (0.0, 1.0, 0.0, 0.0, "(0,1,0,0)"),
        (0.0, 0.0, 1.0, 0.0, "(0,0,1,0)"),
        (0.0, 0.0, 0.0, 1.0, "(0,0,0,1)"),
    ];
    println!(
        "  {:>15} -> {:>8} {:>8} {:>8}",
        "S^3 point", "n_x", "n_y", "n_z"
    );
    println!("  {}", "-".repeat(50));
    for (a, b, c, d, label) in &test_points {
        let n = HopfFibration::map(*a, *b, *c, *d);
        println!("  {:>15} -> {:8.4} {:8.4} {:8.4}", label, n.x, n.y, n.z);
    }

    // 2. Create hopfion with Q_H = 1
    println!("\n=== Hopfion Construction ===");
    let grid_size = (16, 16, 16);
    let radius = 5.0e-9; // 5 nm
    let hopfion = Hopfion::new(grid_size, radius, 1).expect("Failed to create hopfion");

    println!("  {}", hopfion);
    println!(
        "  Grid: {}x{}x{} = {} cells",
        grid_size.0,
        grid_size.1,
        grid_size.2,
        grid_size.0 * grid_size.1 * grid_size.2
    );
    println!("  Cell size: {:.3} nm", hopfion.cell_size * 1e9);
    println!("  Radius: {:.1} nm", radius * 1e9);
    println!("  Normalized: {}", hopfion.is_normalized(1e-10));

    // Sample magnetization at key points
    println!("\n  Magnetization at selected grid points:");
    let mid = grid_size.0 / 2;
    for (label, ix, iy, iz) in &[
        ("Center", mid, mid, mid),
        ("Corner", 0usize, 0usize, 0usize),
        ("Edge-x", grid_size.0 - 1, mid, mid),
        ("Edge-z", mid, mid, grid_size.2 - 1),
    ] {
        let m = hopfion.magnetization_at(*ix, *iy, *iz);
        println!(
            "    {:<10} ({},{},{}): m = ({:.3}, {:.3}, {:.3}), |m| = {:.6}",
            label,
            ix,
            iy,
            iz,
            m.x,
            m.y,
            m.z,
            m.magnitude()
        );
    }

    // 3. Hopf invariant (topological charge)
    println!("\n=== Hopf Invariant ===");
    let q_h = HopfInvariant::calculate(&hopfion).expect("Failed to calculate Hopf invariant");
    println!("  Q_H = {:.4} (expected: 1)", q_h);
    println!("  Deviation from integer: {:.4}", (q_h - q_h.round()).abs());

    // 4. Energy calculations
    println!("\n=== Energy Analysis ===");
    let exchange_a = 1.0e-11; // 10 pJ/m (typical ferromagnet)
    let e_exchange = HopfionEnergy::exchange_energy(&hopfion, exchange_a)
        .expect("Failed to compute exchange energy");
    println!("  Exchange stiffness A = {:.1e} J/m", exchange_a);
    println!("  Exchange energy: {:.4e} J", e_exchange);

    // Total energy with default parameters
    let params = HopfionEnergyParams::default();
    let e_total =
        HopfionEnergy::total_energy(&hopfion, &params).expect("Failed to compute total energy");
    println!(
        "  Total energy (A={:.0e}, D={:.0e}, B_z={:.1} T): {:.4e} J",
        params.exchange_stiffness, params.dmi_constant, params.external_field.z, e_total
    );

    // 5. Compare different Hopf charges
    println!("\n=== Hopf Charge Comparison ===");
    println!(
        "  {:>5} {:>15} {:>15} {:>15}",
        "Q_H", "E_exchange (J)", "Q_H_calc", "Deviation"
    );
    println!("  {}", "-".repeat(55));

    for charge in [1, 2] {
        let h = Hopfion::new(grid_size, radius, charge).expect("Failed to create hopfion");
        let e_ex = HopfionEnergy::exchange_energy(&h, exchange_a)
            .expect("Failed to compute exchange energy");
        let q = HopfInvariant::calculate(&h).expect("Failed to calculate Hopf invariant");
        println!(
            "  {:>5} {:>15.4e} {:>15.4} {:>15.4}",
            charge,
            e_ex,
            q,
            (q - charge as f64).abs()
        );
    }

    // 6. Compare profile functions
    println!("\n=== Profile Function Comparison ===");
    let profiles = [
        (ProfileFunction::Linear, "Linear"),
        (ProfileFunction::Gaussian { sigma: 3.0e-9 }, "Gaussian"),
        (ProfileFunction::Tanh { wall_width: 3.0e-9 }, "Tanh"),
    ];
    println!("  {:>10} {:>15} {:>12}", "Profile", "E_exchange (J)", "Q_H");
    println!("  {}", "-".repeat(42));
    for (profile, name) in &profiles {
        let h = Hopfion::with_profile(grid_size, radius, 1, *profile)
            .expect("Failed to create hopfion");
        let e_ex = HopfionEnergy::exchange_energy(&h, exchange_a)
            .expect("Failed to compute exchange energy");
        let q = HopfInvariant::calculate(&h).expect("Failed to calculate Hopf invariant");
        println!("  {:>10} {:>15.4e} {:>12.4}", name, e_ex, q);
    }

    // 7. Uniform reference state
    println!("\n=== Uniform vs Hopfion ===");
    let uniform =
        Hopfion::uniform(grid_size, hopfion.cell_size).expect("Failed to create uniform state");
    let e_uniform = HopfionEnergy::exchange_energy(&uniform, exchange_a)
        .expect("Failed to compute uniform energy");
    println!("  Uniform state exchange energy: {:.4e} J", e_uniform);
    println!(
        "  Hopfion excitation energy: {:.4e} J",
        e_exchange - e_uniform
    );

    println!("\n=== Summary ===");
    println!("Hopfion physics:");
    println!(
        "  - 3D topological solitons with Hopf invariant Q_H = {:.2}",
        q_h
    );
    println!(
        "  - Toroidal magnetization texture on {}^3 grid",
        grid_size.0
    );
    println!("  - Exchange energy E_ex = {:.2e} J", e_exchange);
    println!("  - Stabilized by frustrated exchange and DMI");

    Ok(())
}
