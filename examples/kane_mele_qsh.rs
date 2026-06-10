//! Kane-Mele Quantum Spin Hall Insulator: Z2 Topology and Edge States
//!
//! **Difficulty**: ⭐⭐⭐⭐
//! **Category**: Topological Physics
//! **Physics**: Z2 topological invariant; spin-orbit coupling on honeycomb lattice;
//! helical edge states; topological phase transition; Rashba effect on topology.
//!
//! ## Background
//!
//! The Kane-Mele model (2005) describes graphene with spin-orbit coupling.  The
//! Bloch Hamiltonian in the (A↑, B↑, A↓, B↓) basis is:
//!
//! ```text
//! H(k) = ε(k)·σ_x ⊗ I + Δ(k)·σ_y ⊗ I    (hopping)
//!       + λ_SO · ν_ij · I ⊗ σ_z            (intrinsic SOC — topological gap)
//!       + λ_R · (σ × d̂)_z ⊗ I              (Rashba SOC — breaks QSH for large λ_R)
//!       + λ_V · ξ_i · I ⊗ I                 (staggered potential — trivial gap)
//! ```
//!
//! The Z2 invariant distinguishes the QSH phase (ν=1, helical edges) from the
//! trivial insulator (ν=0).  The topological phase requires λ_SO > 0 and
//! λ_R < 2√3 λ_SO.
//!
//! ## References
//!
//! - C. L. Kane & E. J. Mele, "Quantum Spin Hall Effect in Graphene",
//!   *Phys. Rev. Lett.* **95**, 226801 (2005).
//! - C. L. Kane & E. J. Mele, "Z2 Topological Order and the Quantum Spin Hall Effect",
//!   *Phys. Rev. Lett.* **95**, 146802 (2005).
//! - M. Z. Hasan & C. L. Kane, "Colloquium: Topological insulators",
//!   *Rev. Mod. Phys.* **82**, 3045 (2010).
//! - T. Fukui, Y. Hatsugai & H. Suzuki, "Chern numbers in discretized Brillouin zone",
//!   *J. Phys. Soc. Jpn.* **74**, 1674 (2005).

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Kane-Mele Model: Quantum Spin Hall Insulator & Z2 Topology");
    println!("=================================================================\n");

    // -------------------------------------------------------------------------
    // 1. Phase diagram: sweep λ_SO with λ_R = 0
    // -------------------------------------------------------------------------
    println!("--- Phase Diagram: Z2 invariant vs λ_SO (λ_R = 0, λ_V = 0) ---");
    println!(
        "  {:>10}  {:>10}  {:>8}  {:>14}  {:>12}",
        "λ_SO/t", "λ_R/t", "Z2", "Gap (meV)", "Phase"
    );
    let t = 1.0_f64; // hopping (eV scale)
    let a = 2.46e-10_f64;
    let lambda_so_vals = [0.0, 0.02, 0.05, 0.1, 0.2, 0.3];
    for &lso in &lambda_so_vals {
        let model = KaneMeleModel::new(t, lso * t, 0.0, 0.0)?;
        let z2 = model.z2_invariant().unwrap_or(-1);
        let gap = model.band_gap().unwrap_or(0.0) * 1000.0; // meV
        let phase = if z2 == 1 {
            "QSH (topological)"
        } else {
            "Trivial insulator"
        };
        println!(
            "  {:>10.3}  {:>10.3}  {:>8}  {:>14.2}  {:>12}",
            lso, 0.0, z2, gap, phase
        );
    }
    println!();

    // -------------------------------------------------------------------------
    // 2. Rashba-driven topological phase transition
    // -------------------------------------------------------------------------
    println!("--- Rashba Effect: Topological → Trivial Transition ---");
    println!(
        "  λ_SO = 0.1·t fixed.  Critical λ_R = 2√3·λ_SO ≈ {:.3}·t",
        2.0 * 3.0_f64.sqrt() * 0.1
    );
    println!(
        "  {:>10}  {:>10}  {:>8}  {:>14}",
        "λ_R/t", "λ_R/λ_SO", "Z2", "Gap (meV)"
    );
    let lso_fixed = 0.1 * t;
    let rashba_vals = [
        0.0,
        0.1,
        0.2,
        0.3,
        0.35,
        2.0_f64.sqrt() * 3.0_f64.sqrt() * 0.1,
        0.4,
        0.5,
    ];
    for &lr in &rashba_vals {
        let model = KaneMeleModel::new(t, lso_fixed, lr * t, 0.0)?;
        let z2 = model.z2_invariant().unwrap_or(-1);
        let gap = model.band_gap().unwrap_or(0.0) * 1000.0;
        println!(
            "  {:>10.3}  {:>10.3}  {:>8}  {:>14.2}",
            lr,
            lr / 0.1,
            z2,
            gap
        );
    }
    println!();

    // -------------------------------------------------------------------------
    // 3. Staggered potential λ_V: gap inversion
    // -------------------------------------------------------------------------
    println!("--- Staggered Potential λ_V: Trivial gap vs. Topological gap ---");
    println!(
        "  λ_SO = 0.1·t, λ_R = 0.  Critical λ_V = 3√3·λ_SO ≈ {:.3}·t",
        3.0 * 3.0_f64.sqrt() * 0.1
    );
    println!("  {:>10}  {:>8}  {:>14}", "λ_V/t", "Z2", "Gap (meV)");
    let lv_vals = [0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0];
    for &lv in &lv_vals {
        let model = KaneMeleModel::new(t, lso_fixed, 0.0, lv * t)?;
        let z2 = model.z2_invariant().unwrap_or(-1);
        let gap = model.band_gap().unwrap_or(0.0) * 1000.0;
        println!("  {:>10.3}  {:>8}  {:>14.2}", lv, z2, gap);
    }
    println!();

    // -------------------------------------------------------------------------
    // 4. Band structure at key k-points for topological phase
    // -------------------------------------------------------------------------
    let topo = KaneMeleModel::topological_phase();
    println!("--- Band Structure at High-Symmetry Points (Topological Phase) ---");
    println!("  Model: λ_SO = 0.06t, λ_R = 0, λ_V = 0");
    let k_points: &[(&str, f64, f64)] = &[
        ("Γ", 0.0, 0.0),
        ("K", 4.0 * std::f64::consts::PI / (3.0 * a), 0.0),
        ("K'", -4.0 * std::f64::consts::PI / (3.0 * a), 0.0),
        (
            "M",
            std::f64::consts::PI / a,
            std::f64::consts::PI / (a * 3.0_f64.sqrt()),
        ),
    ];
    println!(
        "  {:>8}  {:>12}  {:>12}  {:>12}  {:>12}",
        "k-point", "E₁ (eV)", "E₂ (eV)", "E₃ (eV)", "E₄ (eV)"
    );
    for (name, kx, ky) in k_points {
        match topo.energy_bands(*kx, *ky) {
            Ok(bands) => println!(
                "  {:>8}  {:>12.6}  {:>12.6}  {:>12.6}  {:>12.6}",
                name, bands[0], bands[1], bands[2], bands[3]
            ),
            Err(e) => println!("  {:>8}  ERROR: {}", name, e),
        }
    }
    let gap_topo = topo.band_gap().unwrap_or(0.0);
    let z2_topo = topo.z2_invariant().unwrap_or(-1);
    println!("  Direct gap = {:.4} eV,  Z2 = {}", gap_topo, z2_topo);
    println!("  is_topological() = {}", topo.is_topological());
    println!();

    // -------------------------------------------------------------------------
    // 5. Trivial phase comparison
    // -------------------------------------------------------------------------
    let trivial = KaneMeleModel::trivial_phase();
    let gap_triv = trivial.band_gap().unwrap_or(0.0);
    let z2_triv = trivial.z2_invariant().unwrap_or(-1);
    println!("--- Trivial Insulator Phase (λ_V dominates) ---");
    println!(
        "  Gap = {:.4} eV,  Z2 = {},  is_topological() = {}",
        gap_triv,
        z2_triv,
        trivial.is_topological()
    );
    println!();

    // -------------------------------------------------------------------------
    // 6. Edge state spectrum for topological phase (strip geometry)
    // -------------------------------------------------------------------------
    println!("--- Helical Edge States (Strip Geometry, N=20 cells) ---");
    println!("  (Solving block-tridiagonal Hamiltonian along kx with open BC in y)");
    let n_strip = 16_usize; // MAX_DIM=64 → n_cells*4=64 (exact limit)
    let n_kx = 15_usize;
    let pi = std::f64::consts::PI;
    let kx_min = -pi / a * 0.8;
    let kx_max = pi / a * 0.8;
    match topo.edge_spectrum(n_strip, kx_min, kx_max, n_kx) {
        Ok(spectrum) => {
            println!(
                "  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}",
                "kx·a/π", "E_mid-1", "E_mid", "E_mid+1", "E_mid+2"
            );
            for (kx_val, bands) in &spectrum {
                let n_b = bands.len();
                if n_b >= 4 {
                    let mid = n_b / 2;
                    println!(
                        "  {:>10.4}  {:>10.4}  {:>10.4}  {:>10.4}  {:>10.4}",
                        kx_val * a / pi,
                        bands[mid - 1],
                        bands[mid],
                        bands[mid + 1],
                        bands[mid + 2],
                    );
                }
            }
            // Count kx-points where bands straddle zero energy (in-gap crossings)
            let crossings = spectrum
                .iter()
                .filter(|(_, bands)| {
                    let n_b = bands.len();
                    if n_b < 2 {
                        return false;
                    }
                    let mid = n_b / 2;
                    bands[mid - 1] < 0.0 && bands[mid] > 0.0
                })
                .count();
            println!(
                "\n  Gapless crossings at E=0: {} / {} kx-points",
                crossings, n_kx
            );
            println!("  (Crossings at Fermi level → helical edge conduction)");
        },
        Err(e) => println!("  Edge spectrum: {}", e),
    }
    println!();

    // -------------------------------------------------------------------------
    // 7. Graphene with realistic SOC (tiny but finite)
    // -------------------------------------------------------------------------
    println!("--- Graphene with Realistic SOC ---");
    let graphene = KaneMeleModel::graphene_with_soc(1.0e-3 * t);
    let gap_gr = graphene.band_gap().unwrap_or(0.0) * 1000.0;
    let z2_gr = graphene.z2_invariant().unwrap_or(-1);
    println!("  λ_SO = 1 meV (intrinsic graphene SOC)");
    println!("  Gap  = {:.4} meV", gap_gr);
    println!(
        "  Z2   = {}  (topological gap, too small to observe at RT)",
        z2_gr
    );
    println!("  (Bilayer graphene with adatoms can enhance this to ~10 meV)");

    println!("\n=================================================================");
    println!("  Summary: Kane-Mele model hosts QSH phase for small SOC;");
    println!("  Rashba coupling λ_R > 2√3·λ_SO drives topological → trivial.");
    println!("  Helical edge states are protected by time-reversal symmetry.");
    println!("=================================================================\n");

    Ok(())
}
