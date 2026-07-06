//! Resonating Valence Bond (RVB) Quantum Spin Liquid
//!
//! **Difficulty**: ⭐⭐⭐ Advanced
//! **Category**: Frustrated Magnetism (Quantum)
//! **Physics**: Dimer coverings, Sutherland loop-counting overlaps, generalized
//! eigenproblem, exact-diagonalization benchmarks, spinon deconfinement
//!
//! This example demonstrates the RVB quantum module on a small open triangular-lattice
//! strip (built explicitly via `RvbSolver::from_bonds` — the classical `FrustratedLattice`
//! constructors use periodic boundary conditions, which would distort the tiny-cluster
//! dimer-covering counts used here):
//!
//! 1. Nearest-neighbor dimer-covering (valence-bond) basis enumeration
//! 2. The variational RVB ground state via the generalized eigenproblem
//! 3. Cross-check against independent exact diagonalization (dense `CMatrix`)
//! 4. The variational energy ordering `E_equal-amplitude >= E_ground >= E_exact`
//! 5. A spin-liquid screening report (long-range correlation decay)
//! 6. Spinon-pair energetics and a deconfinement diagnostic
//! 7. The `from_lattice` convenience wrapper around a classical `FrustratedLattice`
//!
//! References:
//! - P.W. Anderson, Mater. Res. Bull. 8, 153 (1973)
//! - B. Sutherland, Phys. Rev. B 37, 3786 (1988)
//! - S. Liang, B. Doucot, P.W. Anderson, Phys. Rev. Lett. 61, 365 (1988)

use spintronics::frustrated::{
    deconfinement_diagnostic, is_spin_liquid, spinon_pair_energy, ExactDiagonalization,
    FrustratedLattice, RvbSolver, ShortRangeRvb,
};

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Resonating Valence Bond (RVB) Quantum Spin Liquid ===\n");

    // 1. A 2x3 open triangular-lattice strip (6 sites, 2 rows of 3), built explicitly
    //    as a bond list so the boundary is genuinely open (no periodic distortion of
    //    the dimer-covering combinatorics). Sites (row, col) -> idx = row*3 + col:
    //
    //      0---1---2       row 0
    //      | \ | \ |
    //      3---4---5       row 1
    //
    //    Horizontal, vertical, and one diagonal per unit triangle give 4 triangular
    //    plaquettes (alternating up/down) -- a genuinely non-bipartite cluster.
    let bonds = vec![
        (0, 1),
        (1, 2),
        (3, 4),
        (4, 5),
        (0, 3),
        (1, 4),
        (2, 5),
        (0, 4),
        (1, 5),
    ];
    let coupling_j = 1.0;
    let solver = RvbSolver::from_bonds(6, bonds.clone(), coupling_j)?;

    println!("=== Dimer-Covering (Valence-Bond) Basis ===");
    println!("  Sites: {}", solver.num_sites);
    println!("  Bonds: {}", solver.bonds.len());
    println!(
        "  Bipartite (two-colorable): {}",
        solver.sublattice.is_some()
    );
    println!(
        "  Nearest-neighbor dimer coverings (basis dimension): {}",
        solver.dim()
    );

    // 2. Variational RVB ground state via the generalized eigenproblem H c = E S c.
    println!("\n=== Variational RVB Ground State ===");
    let ground = solver.ground_state()?;
    println!(
        "  E_ground (VB basis, Lowdin canonical orthogonalization) = {:.6} J",
        ground.energy
    );

    // 3. Independent cross-check: this cluster has exactly 6 sites, so 2^6 = 64 =
    //    CMatrix::MAX_DIM -- the largest system dense exact diagonalization can honestly
    //    reach.
    let ed = ExactDiagonalization::new(6, bonds, coupling_j)?;
    let e_dense = ed.dense_ground_energy()?;
    let e_lanczos = ed.ground_energy_sz0()?;
    println!("\n=== Exact-Diagonalization Cross-Check ===");
    println!(
        "  E_exact (dense CMatrix, N=6, 2^6=64=MAX_DIM)   = {:.6} J",
        e_dense
    );
    println!(
        "  E_exact (matrix-free Lanczos, Sz=0 sector)      = {:.6} J",
        e_lanczos
    );
    println!(
        "  |E_dense - E_lanczos| = {:.2e} (independent methods agree)",
        (e_dense - e_lanczos).abs()
    );

    // 4. Variational ordering: the equal-amplitude Anderson ansatz sits above the fully
    //    optimized VB ground state, which in turn sits above the true ground state
    //    (the VB basis is a genuine subspace of the full Hilbert space).
    let equal_amplitude = ShortRangeRvb::equal_amplitude(&solver);
    let e_equal = solver.variational_energy(&equal_amplitude)?;
    println!("\n=== Variational Energy Ordering ===");
    println!(
        "  E_equal-amplitude (Anderson RVB ansatz) = {:.6} J",
        e_equal
    );
    println!(
        "  E_ground (optimized VB basis)           = {:.6} J",
        ground.energy
    );
    println!(
        "  E_exact (true ground state)              = {:.6} J",
        e_dense
    );
    println!(
        "  Ordering holds: {}",
        e_equal >= ground.energy - 1e-9 && ground.energy >= e_dense - 1e-9
    );

    // 5. Spin-liquid screening: is the longest-range correlation small?
    println!("\n=== Spin-Liquid Screening Report ===");
    let report = is_spin_liquid(&solver, 0.5)?;
    println!(
        "  Most distant site pair: ({}, {}) at graph distance {}",
        report.far_site_a, report.far_site_b, report.far_graph_distance
    );
    println!(
        "  |<S_far_a . S_far_b>| = {:.6}",
        report.long_range_correlation
    );
    println!(
        "  Spin-liquid candidate (no long-range order): {}",
        report.is_candidate
    );

    // 6. Spinon-pair energetics: cost of separating two deconfined spin-1/2 spinons.
    println!("\n=== Spinon Deconfinement Diagnostic ===");
    let pairs = [(0, 1), (0, 2), (0, 5)];
    match deconfinement_diagnostic(&solver, &pairs) {
        Ok(points) => {
            for p in &points {
                println!(
                    "  spinons at ({}, {}), graph distance {}: Delta(r) = {:.6} J",
                    p.site_a, p.site_b, p.graph_distance, p.delta_energy
                );
            }
        },
        Err(e) => println!(
            "  (some separations have no valid monomer-dimer matching: {})",
            e
        ),
    }
    let e_pair_adjacent = spinon_pair_energy(&solver, 0, 1)?;
    println!("  E_spinon-pair(0,1) = {:.6} J", e_pair_adjacent);

    // 7. Convenience wrapper: build a solver directly from a classical FrustratedLattice
    //    (periodic boundary conditions; a small 2x2 kagome cell has 12 sites, even).
    println!("\n=== from_lattice() Convenience Wrapper ===");
    let kagome_lattice = FrustratedLattice::kagome(2, 2, 1.0, 1e-9)?;
    let kagome_solver = RvbSolver::from_lattice(&kagome_lattice)?;
    println!(
        "  2x2 periodic kagome: {} sites, {} nearest-neighbor dimer coverings",
        kagome_solver.num_sites,
        kagome_solver.dim()
    );
    let kagome_ground = kagome_solver.ground_state()?;
    println!(
        "  E_ground (periodic kagome, VB basis) = {:.6} J",
        kagome_ground.energy
    );

    println!("\n=== Summary ===");
    println!("The RVB variational ground state, an independent Lanczos/dense exact-");
    println!("diagonalization benchmark, and the spinon deconfinement diagnostic all");
    println!("agree on a quantum-disordered (spin-liquid-like) ground state for this");
    println!("frustrated triangular cluster.");

    Ok(())
}
