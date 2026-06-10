//! Graph Neural Network Potential on a 1D Spin Chain
//!
//! **Difficulty**: ⭐⭐⭐⭐
//! **Category**: Machine Learning / Equivariant message passing
//! **Physics**: Surrogate Hamiltonian for arbitrary spin lattices
//!
//! ## Background
//!
//! Standard neural network potentials require hand-crafted descriptors and a
//! fixed coordination number. A message-passing **graph neural network** (GNN)
//! generalises to arbitrary lattice topology, while the **Cartesian-tensor
//! equivariant** construction from v0.8.0 guarantees that the predicted energy
//! is exactly invariant under simultaneous rotation of every spin.
//!
//! This demo:
//!   1. Builds a 1D periodic ring of `N = 8` sites.
//!   2. Constructs a 2-layer [`GraphMlp`] with 3 scalar + 3 vector channels.
//!   3. Evaluates the energy of three test configurations
//!      (ferromagnetic, antiferromagnetic, 120° spiral) and prints them.
//!   4. Verifies rotation invariance by applying a random SO(3) rotation to
//!      every spin (and to every edge displacement) and checking that the
//!      energy is preserved to ~10⁻¹⁰.
//!
//! ## References
//! - Gilmer et al., "Neural Message Passing for Quantum Chemistry",
//!   *ICML* (2017), arXiv:1704.01212.
//! - Schütt, Unke & Gastegger, "Equivariant Message Passing for the Prediction
//!   of Tensorial Properties and Molecular Spectra", *ICML* (2021),
//!   arXiv:2102.03150.

use spintronics::autodiff::{random_so3, rotate_vector, GraphMlp, LatticeGraph};
use spintronics::vector3::Vector3;

fn main() {
    println!("Graph Neural Network Potential — 1D Ring Demo\n");

    // 1. Build the lattice topology.
    let n_sites = 8;
    let graph = LatticeGraph::ring_1d(n_sites, 1.0).expect("ring construction");
    println!(
        "Lattice: 1D periodic ring with {} sites, {} directed edges",
        graph.n_nodes,
        graph.n_edges(),
    );

    // 2. Build a 2-layer GraphMLP with 3 scalar + 3 vector channels per node.
    let net = GraphMlp::new(2, 3, 3, 0x1234_5678).expect("GraphMlp construction");
    println!(
        "Graph MLP: {} layers, {} total params\n",
        net.layers.len(),
        net.n_params()
    );

    // 3. Evaluate three spin configurations.
    let fm: Vec<Vector3<f64>> = (0..n_sites).map(|_| Vector3::unit_z()).collect();
    let afm: Vec<Vector3<f64>> = (0..n_sites)
        .map(|i| {
            if i % 2 == 0 {
                Vector3::unit_z()
            } else {
                Vector3::unit_z() * -1.0
            }
        })
        .collect();
    let spiral: Vec<Vector3<f64>> = (0..n_sites)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / 3.0;
            Vector3::new(theta.cos(), theta.sin(), 0.0)
        })
        .collect();

    let e_fm = net.energy(&graph, &fm).expect("FM energy");
    let e_afm = net.energy(&graph, &afm).expect("AFM energy");
    let e_spiral = net.energy(&graph, &spiral).expect("spiral energy");
    println!("Configuration energies (untrained random init):");
    println!("  Ferromagnetic (all ẑ):           {:>12.6e}", e_fm);
    println!("  Antiferromagnetic (alternating): {:>12.6e}", e_afm);
    println!("  120° spiral in xy plane:         {:>12.6e}\n", e_spiral);

    // 4. Rotation invariance check.
    let r_mat = random_so3(0x4321_8765);
    let rotated_spins: Vec<Vector3<f64>> = fm.iter().map(|s| rotate_vector(&r_mat, *s)).collect();
    let mut rotated_graph = LatticeGraph::new(graph.n_nodes).expect("rotated graph");
    for &(i, j, r_ij) in &graph.edges {
        rotated_graph
            .add_edge(i, j, rotate_vector(&r_mat, r_ij))
            .expect("rotated edge");
    }
    let e_fm_rot = net
        .energy(&rotated_graph, &rotated_spins)
        .expect("rotated FM");
    let drift = (e_fm - e_fm_rot).abs();
    println!("Rotation-invariance test (FM configuration):");
    println!("  E(original)        = {:>12.6e}", e_fm);
    println!("  E(after R·spins)   = {:>12.6e}", e_fm_rot);
    println!("  |ΔE|               = {:>12.6e}", drift);
    println!(
        "  {} (tolerance 1e-10)",
        if drift < 1.0e-10 { "PASS" } else { "FAIL" },
    );
}
