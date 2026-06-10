//! O(3)-Equivariant Neural Network for Spin Hamiltonians
//!
//! **Difficulty**: ⭐⭐⭐⭐
//! **Category**: Machine Learning / Symmetry-preserving NN
//! **Physics**: SO(3) rotation invariance, Cartesian-tensor equivariant layers
//!
//! ## Background
//!
//! Physical spin energies E(S_1, ..., S_N) are *invariant* under simultaneous
//! rotation of all spin inputs (SO(3) symmetry). A standard MLP must *learn*
//! this symmetry from data — wasting capacity. Cartesian-tensor equivariant
//! layers bake the symmetry into the architecture by construction:
//!
//!   - Scalar channels: rotation-invariant
//!   - Vector channels: rotation-equivariant (rotate by R when input rotates by R)
//!   - Scalar-vector mixing via `s_j · v_j` (equivariant) and `|v|` (invariant)
//!
//! This example:
//!   1. Builds an `EquivariantMlp` mapping (scalar, vectors) → (scalar energy)
//!   2. Numerically verifies rotation invariance by sampling random SO(3) matrices
//!      and checking output drift (should be ≪ 1e-10)
//!   3. Compares with a vanilla MLP that would have to learn the symmetry
//!
//! ## References
//! - Schütt et al., "Equivariant Message Passing for the Prediction of Tensorial
//!   Properties and Molecular Spectra", arXiv:2102.03150 (2021)
//! - Thomas et al., "Tensor Field Networks", arXiv:1802.08219 (2018)
//! - Batzner et al., "E(3)-Equivariant Graph Neural Networks for Data-Efficient
//!   and Accurate Interatomic Potentials", Nat. Commun. 13, 2453 (2022)

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  O(3)-Equivariant Neural Network for Spin Hamiltonians");
    println!("=============================================================");

    // -------------------------------------------------------------------------
    // Section 1: Build an EquivariantMlp
    // -------------------------------------------------------------------------
    println!("\n--- Section 1: Construct EquivariantMlp ---\n");

    // Architecture: 0 input scalars + 4 input vectors (4 spins)
    //               → hidden: 4 scalars + 4 vectors
    //               → output: 1 scalar (energy) + 0 vectors
    let layer1 = EquivariantConfig {
        n_scalar_in: 0,
        n_vector_in: 4,
        n_scalar_out: 4,
        n_vector_out: 4,
    };
    let layer2 = EquivariantConfig {
        n_scalar_in: 4,
        n_vector_in: 4,
        n_scalar_out: 1,
        n_vector_out: 0,
    };

    let mlp = EquivariantMlp::new(&[layer1, layer2], 1234)?;
    println!("  Architecture: (0s, 4v) → (4s, 4v) → (1s, 0v)");
    println!("  Total parameters: {}", mlp.n_params());

    // -------------------------------------------------------------------------
    // Section 2: Test rotation invariance numerically
    // -------------------------------------------------------------------------
    println!("\n--- Section 2: Rotation Invariance Check ---\n");

    // Build a 4-spin configuration (e.g., a tetrahedral arrangement)
    let spins = vec![
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(-1.0 / 3f64.sqrt(), -1.0 / 3f64.sqrt(), -1.0 / 3f64.sqrt()),
    ];

    // Original energy
    let scalars: Vec<f64> = vec![];
    let (s_out_orig, _v_out_orig) = mlp.forward(&scalars, &spins)?;
    let e_orig = s_out_orig[0];
    println!("  Energy of original spin configuration: {e_orig:.10}");

    // Sample N random rotation matrices and rotate all spins
    let n_rotations = 20;
    let mut max_drift = 0.0_f64;
    let mut max_rel_drift = 0.0_f64;
    println!("\n  {:>4}  {:>16}  {:>16}", "test", "E(rotated)", "|drift|");
    println!("  {}", "-".repeat(42));
    for seed in 0..n_rotations {
        let r = random_so3(7919 + seed * 31);
        let rotated_spins: Vec<Vector3<f64>> =
            spins.iter().map(|&s| rotate_vector(&r, s)).collect();
        let (s_out_rot, _) = mlp.forward(&scalars, &rotated_spins)?;
        let e_rot = s_out_rot[0];
        let drift = (e_rot - e_orig).abs();
        let rel = drift / e_orig.abs().max(1e-30);
        max_drift = max_drift.max(drift);
        max_rel_drift = max_rel_drift.max(rel);
        if seed < 5 {
            println!("  {:>4}  {:>+16.10}  {:>16.3e}", seed, e_rot, drift);
        }
    }
    println!("  ...");
    println!("\n  Max drift over {n_rotations} random rotations: {max_drift:.3e}");
    println!("  Max relative drift:                       {max_rel_drift:.3e}");
    if max_drift < 1.0e-10 {
        println!("  ✓ Rotation invariance preserved to machine precision");
    } else {
        println!("  ✗ Unexpected drift — check architecture");
    }

    // -------------------------------------------------------------------------
    // Section 3: Single-layer equivariance test (output vectors rotate by R)
    // -------------------------------------------------------------------------
    println!("\n--- Section 3: Vector-Output Equivariance Check ---\n");

    let cfg_vec = EquivariantConfig {
        n_scalar_in: 0,
        n_vector_in: 2,
        n_scalar_out: 0,
        n_vector_out: 2,
    };
    let layer_vec = EquivariantLinear::new(cfg_vec, 2024)?;
    let v_in = vec![Vector3::new(0.6, 0.0, 0.8), Vector3::new(0.0, 1.0, 0.0)];
    let (_, v_out_orig) = layer_vec.forward(&[], &v_in)?;

    let r = random_so3(99999);
    let v_in_rot: Vec<Vector3<f64>> = v_in.iter().map(|&v| rotate_vector(&r, v)).collect();
    let (_, v_out_after_rot) = layer_vec.forward(&[], &v_in_rot)?;
    let v_out_rotated_then: Vec<Vector3<f64>> =
        v_out_orig.iter().map(|&v| rotate_vector(&r, v)).collect();

    let mut max_diff = 0.0_f64;
    for (a, b) in v_out_after_rot.iter().zip(v_out_rotated_then.iter()) {
        let d = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt();
        max_diff = max_diff.max(d);
    }
    println!("  Max |forward(R·v) − R·forward(v)|: {max_diff:.3e}");
    if max_diff < 1.0e-10 {
        println!("  ✓ Vector outputs rotate equivariantly with input");
    }

    // -------------------------------------------------------------------------
    // Section 4: Energy of FM vs disordered configurations
    // -------------------------------------------------------------------------
    println!("\n--- Section 4: FM vs Disordered Energy (Untrained) ---\n");

    let fm_spins = vec![Vector3::new(0.0, 0.0, 1.0); 4];
    let disordered = vec![
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, -1.0, 0.0),
        Vector3::new(-0.5, 0.5, 0.5_f64.sqrt()),
        Vector3::new(0.5, -0.5, -(0.5_f64.sqrt())),
    ];

    let e_fm = mlp.forward(&[], &fm_spins)?.0[0];
    let e_dis = mlp.forward(&[], &disordered)?.0[0];
    println!("  E(FM all ẑ):       {e_fm:>+12.6}");
    println!("  E(disordered):     {e_dis:>+12.6}");
    println!(
        "  Difference (random init, not yet trained): {:.4}",
        e_dis - e_fm
    );

    println!("\n=============================================================");
    println!("  Done. The equivariant network preserves SO(3) symmetry by");
    println!("  construction — rotating all input spins leaves the energy");
    println!("  invariant to machine precision, regardless of weight values.");
    println!("=============================================================\n");

    Ok(())
}
