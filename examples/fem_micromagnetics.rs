//! Finite Element Micromagnetic Simulation Example
//!
//! **Difficulty**: ⭐⭐⭐ Advanced
//! **Category**: Computational Methods
//! **Physics**: Micromagnetics, FEM, energy minimization
//! **Requirements**: `fem` feature
//!
//! Demonstrates the use of the FEM module for micromagnetic simulations,
//! including mesh generation, energy calculations, and solver usage.
//!
//! Run with:
//! ```bash
//! cargo run --features fem --example fem_micromagnetics
//! ```

use spintronics::fem::{Mesh2D, MicromagneticFEM, SolverParams, SolverType};
use spintronics::material::Ferromagnet;
use spintronics::vector3::Vector3;

fn main() -> Result<(), String> {
    println!("=== Finite Element Micromagnetic Simulation ===\n");

    // 1. Generate a 2D triangular mesh for a thin film
    println!("1. Generating mesh...");
    let width = 200e-9; // 200 nm
    let height = 100e-9; // 100 nm
    let element_size = 10e-9; // 10 nm elements

    let mesh = Mesh2D::rectangle(width, height, element_size)?;
    println!(
        "   - Mesh created: {} nodes, {} elements",
        mesh.n_nodes(),
        mesh.n_elements()
    );

    // 2. Create micromagnetic FEM solver with Permalloy
    println!("\n2. Setting up micromagnetic problem...");
    let material = Ferromagnet::permalloy();
    println!("   - Material: Permalloy");
    println!("   - Saturation magnetization: {:.2e} A/m", material.ms);
    println!("   - Exchange constant: {:.2e} J/m", material.exchange_a);
    println!("   - Damping: {}", material.alpha);

    let ms = material.ms; // Extract for later use
    let mut fem = MicromagneticFEM::new(mesh, material);

    // 3. Calculate initial energies
    println!("\n3. Calculating initial energies...");
    let h_ext = Vector3::new(0.0, 0.0, 10000.0); // 10 kA/m external field in z

    let e_exchange = fem.exchange_energy();
    let e_zeeman = fem.zeeman_energy(h_ext);
    let e_total = fem.total_energy(h_ext);

    println!("   - Exchange energy: {:.3e} J", e_exchange);
    println!("   - Zeeman energy: {:.3e} J", e_zeeman);
    println!("   - Total energy: {:.3e} J", e_total);

    // 4. Include anisotropy effects
    println!("\n4. Including anisotropy...");
    let easy_axis = Vector3::new(0.0, 0.0, 1.0); // Easy axis along z
    let e_anisotropy = fem.uniaxial_anisotropy_energy(easy_axis);
    println!("   - Uniaxial anisotropy energy: {:.3e} J", e_anisotropy);

    let e_total_with_anis = fem.total_energy_with_anisotropy(h_ext, easy_axis);
    println!(
        "   - Total energy (with anisotropy): {:.3e} J",
        e_total_with_anis
    );

    // 5. Include demagnetization effects
    println!("\n5. Including demagnetization...");
    // Thin film demagnetization factors: strong out-of-plane demagnetization
    let demag_factors = Vector3::new(0.0, 0.0, 1.0);
    let e_demag = fem.demagnetization_energy(demag_factors);
    println!("   - Demagnetization energy: {:.3e} J", e_demag);

    let e_complete = fem.total_energy_complete(h_ext, easy_axis, demag_factors);
    println!("   - Complete micromagnetic energy: {:.3e} J", e_complete);

    // 6. Calculate average magnetization
    println!("\n6. Analyzing magnetization state...");
    let m_avg = fem.average_magnetization();
    println!(
        "   - Average magnetization: ({:.2e}, {:.2e}, {:.2e}) A/m",
        m_avg.x, m_avg.y, m_avg.z
    );
    println!("   - Magnitude: {:.2e} A/m", m_avg.magnitude());
    println!(
        "   - Normalized: ({:.3}, {:.3}, {:.3})",
        m_avg.x / ms,
        m_avg.y / ms,
        m_avg.z / ms
    );

    // 7. Demonstrate different solvers
    println!("\n7. Testing different linear solvers...");
    demonstrate_solvers()?;

    // 8. Solve for equilibrium magnetization
    println!("\n8. Solving for equilibrium magnetization...");
    let h_applied = Vector3::new(5000.0, 0.0, 10000.0); // Mixed field
    match fem.solve_equilibrium(h_applied) {
        Ok(m_eq) => {
            let m_eq_avg = average_vector(&m_eq);
            println!(
                "   - Equilibrium magnetization: ({:.2e}, {:.2e}, {:.2e}) A/m",
                m_eq_avg.x, m_eq_avg.y, m_eq_avg.z
            );
            println!("   ✓ Equilibrium solution converged");
        },
        Err(e) => println!("   ✗ Equilibrium solver failed: {}", e),
    }

    println!("\n=== Simulation Complete ===");
    Ok(())
}

/// Demonstrate different linear solver types
fn demonstrate_solvers() -> Result<(), String> {
    use spintronics::fem::assembly::SparseMatrix;

    // Create a simple 3x3 SPD system
    let mut a = SparseMatrix::new(3, 3);
    a.add_entry(0, 0, 4.0);
    a.add_entry(0, 1, 1.0);
    a.add_entry(1, 0, 1.0);
    a.add_entry(1, 1, 4.0);
    a.add_entry(1, 2, 1.0);
    a.add_entry(2, 1, 1.0);
    a.add_entry(2, 2, 4.0);

    let b = vec![5.0, 6.0, 5.0];

    // Test CG solver
    let params_cg = SolverParams {
        max_iter: 100,
        tolerance: 1e-8,
        omega: 1.0,
        preconditioner: spintronics::fem::Preconditioner::None,
    };
    match spintronics::fem::solve_linear_system_with_params(&a, &b, SolverType::CG, params_cg) {
        Ok(x) => println!(
            "   - CG solver: x = [{:.4}, {:.4}, {:.4}]",
            x[0], x[1], x[2]
        ),
        Err(e) => println!("   ✗ CG failed: {}", e),
    }

    // Test BiCGSTAB solver
    match spintronics::fem::solve_linear_system(&a, &b, SolverType::BiCGSTAB) {
        Ok(x) => println!(
            "   - BiCGSTAB solver: x = [{:.4}, {:.4}, {:.4}]",
            x[0], x[1], x[2]
        ),
        Err(e) => println!("   ✗ BiCGSTAB failed: {}", e),
    }

    // Test SOR solver with optimal relaxation
    let params_sor = SolverParams {
        max_iter: 100,
        tolerance: 1e-6,
        omega: 1.5, // Typical optimal value
        preconditioner: spintronics::fem::Preconditioner::None,
    };
    match spintronics::fem::solve_linear_system_with_params(&a, &b, SolverType::SOR, params_sor) {
        Ok(x) => println!(
            "   - SOR solver: x = [{:.4}, {:.4}, {:.4}]",
            x[0], x[1], x[2]
        ),
        Err(e) => println!("   ✗ SOR failed: {}", e),
    }

    Ok(())
}

/// Calculate average of a vector field
fn average_vector(vectors: &[Vector3<f64>]) -> Vector3<f64> {
    let n = vectors.len() as f64;
    if n == 0.0 {
        return Vector3::zero();
    }

    let mut sum = Vector3::zero();
    for v in vectors {
        sum = sum + *v;
    }

    Vector3::new(sum.x / n, sum.y / n, sum.z / n)
}
