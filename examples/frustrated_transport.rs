//! Geometric Frustration Effects on Transport: Chirality-Driven Topological
//! Hall Effect
//!
//! **Difficulty**: ⭐⭐⭐ Advanced
//! **Category**: Frustrated Magnetism / Transport
//! **Physics**: Scalar spin chirality, Berg-Lüscher solid angle, emergent
//!              electrodynamics, topological Hall effect
//!
//! This example wires the scalar spin chirality of a `FrustratedLattice`
//! configuration into the emergent-field topological Hall machinery of
//! `spintronics::effect::topological_hall`, demonstrating the full
//! chirality -> emergent-field -> Hall-response pipeline:
//!
//! 1. The coplanar 120-degree Néel order carries zero scalar spin chirality
//!    (any three coplanar spins have a vanishing scalar triple product) and
//!    therefore drives no topological Hall response.
//! 2. Canting that order into a noncoplanar "umbrella" state switches on a
//!    nonzero chirality on every elementary triangular plaquette -- but
//!    whether this survives as a *net* (lattice-summed) signal is itself a
//!    lattice-geometry-dependent question: on the triangular lattice the
//!    "up" and "down" elementary triangles carry equal-and-opposite
//!    chirality and cancel exactly in the net (like a Néel order canceling
//!    in net magnetization), while on the kagome lattice they carry the
//!    *same* sign and add constructively to a genuinely nonzero net
//!    topological Hall resistivity.
//! 3. A generic (disordered) noncoplanar texture on the triangular lattice
//!    does not have this cancellation.
//! 4. Reversing the global chirality (mirroring one spin component) exactly
//!    reverses the sign of the Hall response.
//!
//! # References
//! - Y. Taguchi et al., "Spin Chirality, Berry Phase, and Anomalous Hall
//!   Effect in a Frustrated Ferromagnet", Science 291, 2573 (2001)
//! - H. Berg & M. Lüscher, Nucl. Phys. B 190, 412 (1981)
//! - A. Van Oosterom & J. Strackee, IEEE Trans. Biomed. Eng. BME-30, 125 (1983)

use spintronics::frustrated::{frustration_hall_response, FrustratedLattice, FrustratedTransport};
use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Geometric Frustration Effects on Transport ===\n");

    let mnsi = TopologicalHall::mnsi();

    // =========================================================================
    // 1. Coplanar 120-degree Neel order: reference case, zero chirality
    // =========================================================================
    println!("=== 1. Coplanar 120-degree Neel Order (Reference) ===");
    let mut lat = FrustratedLattice::triangular(8, 8, 1.0, 3e-10)?;
    lat.set_120_degree_order();

    let transport = FrustratedTransport::from_lattice(&lat)?;
    println!("  Plaquettes found: {}", transport.num_plaquettes());
    println!(
        "  Total chirality: {:.3e} (coplanar spins -> exactly zero)",
        transport.total_chirality()
    );
    println!(
        "  Hall resistivity: {:.3e} Ohm*cm (no topological Hall signal)\n",
        transport.hall_resistivity(&mnsi)
    );

    // =========================================================================
    // 2. Triangular umbrella order: nonzero per-plaquette chirality, but the
    //    "up" and "down" triangles cancel exactly in the net
    // =========================================================================
    println!("=== 2. Triangular Umbrella Order: Per-Plaquette Nonzero, Net Cancels ===");
    println!(
        "  {:>12} {:>14} {:>16} {:>18}",
        "theta (deg)", "mean |chi|", "total chi (net)", "rho_THE (Ohm*cm)"
    );
    println!("  {}", "-".repeat(66));
    for theta_deg in [10.0_f64, 30.0, 45.0, 54.7, 60.0, 75.0, 90.0] {
        let theta = theta_deg.to_radians();
        let mut lat = FrustratedLattice::triangular(8, 8, 1.0, 3e-10)?;
        lat.set_umbrella_order(theta);

        let transport = FrustratedTransport::from_lattice(&lat)?;
        let mean_abs_chi: f64 = transport
            .plaquettes
            .iter()
            .map(|p| p.chirality.abs())
            .sum::<f64>()
            / transport.num_plaquettes() as f64;
        let rho = transport.hall_resistivity(&mnsi);
        println!(
            "  {:>12.1} {:>14.4} {:>16.4e} {:>18.4e}",
            theta_deg,
            mean_abs_chi,
            transport.total_chirality(),
            rho
        );
    }
    println!(
        "  (every plaquette is noncoplanar -- mean |chi| grows with theta -- but the up/down\n\
         \x20  triangles' equal-and-opposite chirality makes the net total ~0 for every theta)\n"
    );

    // =========================================================================
    // 3. A generic (disordered) noncoplanar texture does not cancel
    // =========================================================================
    println!("=== 3. Generic Disordered Noncoplanar Texture (Triangular Lattice) ===");
    let mut disordered_lat = FrustratedLattice::triangular(8, 8, 1.0, 3e-10)?;
    let mut rng = spintronics::frustrated::Xorshift64::new(2024)?;
    for spin in disordered_lat.spins.iter_mut() {
        *spin = rng.random_unit_vector();
    }
    let disordered_transport = FrustratedTransport::from_lattice(&disordered_lat)?;
    println!(
        "  Total chirality: {:.4} (no symmetry-enforced cancellation)",
        disordered_transport.total_chirality()
    );
    println!(
        "  Hall resistivity: {:.4e} Ohm*cm\n",
        disordered_transport.hall_resistivity(&mnsi)
    );

    // =========================================================================
    // 4. Kagome umbrella order: up/down triangles add constructively
    // =========================================================================
    println!("=== 4. Kagome Umbrella Order: Genuine Net Topological Hall Signal ===");
    println!(
        "  {:>12} {:>14} {:>16} {:>18}",
        "theta (deg)", "mean chi", "total chi (net)", "rho_THE (Ohm*cm)"
    );
    println!("  {}", "-".repeat(66));
    for theta_deg in [10.0_f64, 30.0, 45.0, 54.7, 60.0, 75.0, 90.0] {
        let theta = theta_deg.to_radians();
        let mut kagome_lat = FrustratedLattice::kagome(6, 6, 1.0, 5e-10)?;
        kagome_lat.set_umbrella_order(theta);

        let kagome_transport = FrustratedTransport::from_lattice(&kagome_lat)?;
        let rho = kagome_transport.hall_resistivity(&mnsi);
        println!(
            "  {:>12.1} {:>14.4} {:>16.4e} {:>18.4e}",
            theta_deg,
            kagome_transport.mean_chirality(),
            kagome_transport.total_chirality(),
            rho
        );
    }
    println!();

    // =========================================================================
    // 5. Sign reversal under global chirality reversal (kagome baseline)
    // =========================================================================
    println!("=== 5. Sign Reversal Under Global Chirality Reversal (Kagome) ===");
    let mut lat = FrustratedLattice::kagome(6, 6, 1.0, 5e-10)?;
    lat.set_umbrella_order(std::f64::consts::FRAC_PI_3);

    let rho_before = frustration_hall_response(&lat, &mnsi)?;
    println!("  Before mirroring: rho_THE = {:.4e} Ohm*cm", rho_before);

    for spin in lat.spins.iter_mut() {
        spin.z = -spin.z; // mirror one spin component of every site
    }
    let rho_after = frustration_hall_response(&lat, &mnsi)?;
    println!(
        "  After mirroring Sz -> -Sz: rho_THE = {:.4e} Ohm*cm",
        rho_after
    );

    let sign_reversed = (rho_before + rho_after).abs() < 1e-9 * rho_before.abs().max(1e-9);
    println!("  Exact sign reversal confirmed: {}\n", sign_reversed);

    println!("=== Summary ===");
    println!("  - Coplanar (theta=90 deg) 120-degree order: chi=0, no topological Hall signal");
    println!("  - Triangular umbrella order: every plaquette is noncoplanar, but the lattice's");
    println!("    up/down triangle relationship cancels the net response exactly");
    println!("  - A generic disordered noncoplanar texture on the same lattice does not cancel");
    println!("  - Kagome umbrella order: up/down triangles reinforce, giving a genuine nonzero");
    println!("    net topological Hall resistivity");
    println!("  - Global chirality reversal (mirroring Sz) exactly reverses the Hall response");

    Ok(())
}
