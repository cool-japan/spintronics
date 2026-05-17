//! Higher-Order Topological Insulator: BBH Corner States
//!
//! **Difficulty**: ⭐⭐⭐⭐
//! **Category**: Topological Condensed Matter
//! **Physics**: BBH quadrupole insulator, corner states, breathing kagome HOTI
//!
//! ## Background
//!
//! Higher-Order Topological Insulators (HOTIs) are a class of topological phases
//! that host protected boundary modes in codimension greater than one. In 2D, a
//! second-order TI has corner modes rather than edge modes. The paradigmatic
//! example is the Benalcazar-Bernevig-Hughes (BBH) model on the square lattice,
//! which supports a quantised electric quadrupole moment q_xy = 1/2 in the
//! topological phase. This is detected via a nested Wilson loop calculation:
//! the Wannier centres of the occupied bands wind around the BZ and produce a
//! quantised nested polarisation of ±1/2.
//!
//! The BBH Hamiltonian in the 4-band basis (A, B, C, D):
//!   H(k) = (λ_x + γ_x cos kx)Γ₁ + γ_x sin kx·Γ₂
//!         + (λ_y + γ_y cos ky)Γ₃ + γ_y sin ky·Γ₄
//!
//! Topological phase: |γ| > |λ| in both directions.
//! Trivial phase:     |λ| > |γ| in at least one direction.
//!
//! The breathing kagome lattice provides a second HOTI realisation with C₃
//! rotational symmetry and corner states protected by three-fold rotation.
//!
//! ## References
//! - Benalcazar, Bernevig & Hughes, Science 357, 61 (2017)
//! - Ezawa, PRL 120, 026801 (2018)
//! - Liu & Liu, PRB 99, 075129 (2019)

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  Higher-Order Topological Insulator: BBH Corner States");
    println!("=============================================================");

    // -------------------------------------------------------------------------
    // Section 1: BBH Phase Diagram — sweep λ/γ from 0.2 to 2.0
    // Fix γ_x = γ_y = 1.0, vary λ_x = λ_y = λ from 0.2 to 2.0
    // -------------------------------------------------------------------------
    println!("\n--- Section 1: BBH Phase Diagram (γ = 1.0, λ swept) ---\n");

    let gamma: f64 = 1.0;
    let lambda_values: &[f64] = &[0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0];

    println!(
        "  {:>10}  {:>10}  {:>12}  {:>8}  {:>8}",
        "λ_x/γ_x", "λ_y/γ_y", "gap (eV)", "q_xy", "HOTI?"
    );
    println!("  {}", "-".repeat(56));

    for &lam in lambda_values {
        let model = BbhModel::new(lam, gamma, lam, gamma)?;
        let gap = model.band_gap(12, 12);
        // Use a coarser grid for the quadrupole moment to keep the example fast
        let q_xy = model.quadrupole_moment(8, 8)?;
        let is_hoti = model.is_higher_order_topological();

        // Map q_xy to [-0.5, 0.5] for display
        let q_display = {
            let w = q_xy.rem_euclid(1.0);
            if w > 0.5 {
                w - 1.0
            } else {
                w
            }
        };

        let ratio = lam / gamma;
        println!(
            "  {:>10.2}  {:>10.2}  {:>12.4}  {:>8.3}  {:>8}",
            ratio,
            ratio,
            gap,
            q_display,
            if is_hoti { "YES" } else { "no" }
        );
    }

    println!("\n  → Topological transition occurs at λ = γ (ratio = 1.0)");

    // -------------------------------------------------------------------------
    // Section 2: Topological vs Trivial — compare canonical phases
    // -------------------------------------------------------------------------
    println!("\n--- Section 2: Topological vs Trivial BBH Phase ---\n");

    let topo = BbhModel::topological_phase(); // λ=0.5, γ=1.0
    let trivial = BbhModel::trivial_phase(); // λ=1.0, γ=0.5

    for (label, model) in [
        ("Topological (λ=0.5, γ=1.0)", &topo),
        ("Trivial     (λ=1.0, γ=0.5)", &trivial),
    ] {
        let gap = model.band_gap(16, 16);
        let q_xy = model.quadrupole_moment(8, 8)?;
        let q_display = {
            let w = q_xy.rem_euclid(1.0);
            if w > 0.5 {
                w - 1.0
            } else {
                w
            }
        };
        let is_hoti = model.is_higher_order_topological();

        println!("  {}:", label);
        println!("    band gap    = {:.4}  (units of hopping energy)", gap);
        println!(
            "    q_xy        = {:.4}  (quantised: 0=trivial, ±0.5=HOTI)",
            q_display
        );
        println!("    is_HOTI     = {}", is_hoti);
        println!();
    }

    // -------------------------------------------------------------------------
    // Section 3: Corner state spectrum on a 3×3 cluster
    // 3×3×4 orbitals = 36 states ≤ CMatrix::MAX_DIM = 64
    // -------------------------------------------------------------------------
    println!("--- Section 3: Corner State Spectrum (3×3 cluster, 36 sites) ---\n");

    // Use topological parameters λ=0.5, γ=1.0
    let solver = CornerStateSolver::new(3, 3, 0.5, 1.0, 0.5, 1.0)?;
    let (eigvals, eigvecs) = solver.solve_finite_cluster()?;

    let n = eigvals.len(); // = 36
    println!("  Cluster dimension: {n} states (3×3×4 orbitals)");

    // Print the 8 eigenvalues nearest to zero (middle of the spectrum)
    let mid = n / 2;
    let lo = mid.saturating_sub(4);
    let hi = (lo + 8).min(n);

    println!("\n  Mid-gap eigenvalues (indices {lo}..{}):", hi - 1);
    println!("  {:>8}  {:>14}  {:>10}", "index", "energy (t)", "label");
    println!("  {}", "-".repeat(38));

    for (i, &ev) in eigvals.iter().enumerate().take(hi).skip(lo) {
        let label = if ev.abs() < 0.15 { "<-- gap" } else { "" };
        println!("  {:>8}  {:>14.6}  {}", i, ev, label);
    }

    // Count corner states in the gap [-0.5, 0.5]
    // For a 3×3 cluster the corner states are not exactly at zero
    // but remain well within the bulk gap (|E| < 0.5 t).
    let gap_min = -0.5_f64;
    let gap_max = 0.5_f64;
    let n_corner = solver.count_corner_states(&eigvals, gap_min, gap_max);
    println!(
        "\n  Corner states in gap ({}, {}) : {}",
        gap_min, gap_max, n_corner
    );

    // Corner localization for the state(s) nearest to zero energy
    // Find the index of the eigenvalue closest to 0.0
    let zero_idx = (lo..hi)
        .min_by(|&i, &j| {
            eigvals[i]
                .abs()
                .partial_cmp(&eigvals[j].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap_or(mid);

    let w = solver.corner_localization(zero_idx, &eigvecs);
    println!(
        "\n  Corner localization of state #{zero_idx} (E = {:.6} t):",
        eigvals[zero_idx]
    );
    println!("    corner (0,0)           : {:.4}", w[0]);
    println!("    corner (lx-1, 0)       : {:.4}", w[1]);
    println!("    corner (0, ly-1)       : {:.4}", w[2]);
    println!("    corner (lx-1, ly-1)    : {:.4}", w[3]);
    let corner_total: f64 = w.iter().sum();
    println!(
        "    total corner weight    : {:.4}  (1.0 = fully corner-localised)",
        corner_total
    );

    // -------------------------------------------------------------------------
    // Section 4: Breathing Kagome comparison
    // -------------------------------------------------------------------------
    println!("\n--- Section 4: Breathing Kagome HOTI ---\n");

    let topo_kg = BreathingKagomeModel::topological_kagome(); // t_inter=1.0 > t_intra=0.5
    let trivial_kg = BreathingKagomeModel::trivial_kagome(); // t_intra=1.0 > t_inter=0.5

    println!(
        "  {:>22}  {:>8}  {:>10}  {:>12}  {:>12}",
        "Model", "t_intra", "t_inter", "band_gap", "is_topo"
    );
    println!("  {}", "-".repeat(70));

    for (label, model) in [
        ("Topological kagome", &topo_kg),
        ("Trivial kagome", &trivial_kg),
    ] {
        // Band gap: minimum positive eigenvalue over the BZ
        let mut min_pos = f64::INFINITY;
        let nk = 16_usize;
        for ix in 0..nk {
            let kx = 2.0 * std::f64::consts::PI * (ix as f64 + 0.5) / nk as f64;
            for iy in 0..nk {
                let ky = 2.0 * std::f64::consts::PI * (iy as f64 + 0.5) / nk as f64;
                if let Ok(e) = model.energy_bands(kx, ky) {
                    // lowest positive band (index 1 since there are 3 bands and e[0]<0)
                    for &ev in &e {
                        if ev > 0.0 && ev < min_pos {
                            min_pos = ev;
                        }
                    }
                }
            }
        }
        let gap = if min_pos.is_infinite() { 0.0 } else { min_pos };

        println!(
            "  {:>22}  {:>8.2}  {:>10.2}  {:>12.4}  {:>12}",
            label,
            model.t_intra,
            model.t_inter,
            gap,
            if model.is_topological() { "YES" } else { "no" }
        );
    }

    // Corner polarization for topological kagome
    let cp = topo_kg.corner_polarization(8, 8)?;
    println!("\n  Corner polarization (topological kagome) = {:.4}", cp);
    println!("  → Near 0.5 indicates non-trivial corner charge (HOTI)");

    println!("\n=============================================================");
    println!("  Done. BBH model shows quantised q_xy = ±0.5 in the");
    println!("  topological phase with 4 mid-gap corner states visible");
    println!("  in the finite-cluster spectrum.");
    println!("=============================================================\n");

    Ok(())
}
