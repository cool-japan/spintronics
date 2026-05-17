//! Topological Magnon Bands: Haldane Model and Chern Numbers
//!
//! **Difficulty**: ⭐⭐⭐⭐ Expert
//! **Category**: Topological Magnonics
//! **Physics**: Haldane honeycomb model, Berry curvature, Chern number,
//!              topological phase transition, edge modes
//!
//! This example demonstrates the bosonic analogue of the Haldane model applied
//! to magnon quasiparticles in a honeycomb ferromagnet. The Dzyaloshinskii-Moriya
//! interaction (DMI) plays the role of Haldane's imaginary next-nearest-neighbour
//! hopping, opening a topological gap at the K/K' points and endowing the lower
//! magnon band with a non-trivial Chern number C = ±1.
//!
//! We demonstrate:
//!
//! 1. Honeycomb Haldane model setup and high-symmetry-point diagonalization
//! 2. Band structure along the Γ→K→M→Γ high-symmetry path
//! 3. Berry curvature on a 20×20 Brillouin-zone grid
//! 4. Chern numbers via the Fukui-Hatsugai discrete-gauge method
//! 5. Topological phase transition: gap and Chern number vs DMI strength
//! 6. Chiral edge modes in a 20-cell strip geometry
//!
//! References:
//! - Haldane, PRL 61, 2015 (1988)
//! - Fukui, Hatsugai & Suzuki, J. Phys. Soc. Jpn. 74, 1674 (2005)
//! - Owerre, J. Phys. Cond. Matter 28, 386001 (2016)
//! - Matsumoto & Murakami, PRL 106, 197202 (2011)

use spintronics::topomagnon::{BerryCurvature, ChernNumber, EdgeModes, MagnonBandModel};

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Topological Magnon Bands: Haldane Model and Chern Numbers ===\n");

    // -------------------------------------------------------------------------
    // Physical parameters (all energies in meV)
    // -------------------------------------------------------------------------
    let j_nn: f64 = 1.0; // nearest-neighbour exchange [meV]
    let j_nnn: f64 = 0.1; // next-nearest-neighbour exchange [meV]
    let dmi: f64 = 0.3; // DMI strength — topological phase (|D| > |D_c|) [meV]
    let h_ext: f64 = 1e-3; // external Zeeman field [meV]  (≈ 1e-3 T × g μ_B)

    // -------------------------------------------------------------------------
    // 1. Model setup
    // -------------------------------------------------------------------------
    println!("=== 1. Honeycomb Haldane Model Setup ===");
    let model = MagnonBandModel::honeycomb_haldane(j_nn, j_nnn, dmi, h_ext)
        .expect("Failed to construct honeycomb Haldane model");

    println!("  Lattice type    : {:?}", model.lattice);
    println!("  Number of bands : {}", model.n_bands());
    println!("  Lattice constant: {:.2e} m", model.a_lattice);
    println!(
        "  J_NN            : {:.3} meV  (nearest-neighbour exchange)",
        j_nn
    );
    println!(
        "  J_NNN           : {:.3} meV  (next-nearest-neighbour exchange)",
        j_nnn
    );
    println!(
        "  DMI             : {:.3} meV  (Haldane mass / topological gap driver)",
        dmi
    );
    println!("  H_ext           : {:.4} meV  (Zeeman bias)", h_ext);

    // High-symmetry k-points for honeycomb BZ (in reciprocal-lattice units a=1)
    // Γ = (0, 0), K = (4π/3, 0), M = (π, π/√3)
    let pi = std::f64::consts::PI;
    let sqrt3 = 3.0_f64.sqrt();
    let k_gamma = (0.0_f64, 0.0_f64);
    let k_k = (4.0 * pi / 3.0, 0.0_f64);
    let k_m = (pi, pi / sqrt3);

    println!("\n  High-symmetry k-points:");
    println!("    Γ = ({:.4}, {:.4})", k_gamma.0, k_gamma.1);
    println!("    K = ({:.4}, {:.4})", k_k.0, k_k.1);
    println!("    M = ({:.4}, {:.4})", k_m.0, k_m.1);

    println!("\n  Eigenvalues at high-symmetry points:");
    println!(
        "  {:>6}  {:>12}  {:>12}  {:>12}",
        "Point", "E_1 [meV]", "E_2 [meV]", "Gap [meV]"
    );
    println!("  {}", "-".repeat(50));

    for (label, k) in [("Γ", k_gamma), ("K", k_k), ("M", k_m)] {
        let (evals, _) = model
            .diagonalize(k)
            .expect("Failed to diagonalize at k-point");
        let gap = evals[1] - evals[0];
        println!(
            "  {:>6}  {:>12.6}  {:>12.6}  {:>12.6}",
            label, evals[0], evals[1], gap
        );
    }

    // -------------------------------------------------------------------------
    // 2. Band structure along Γ→K→M→Γ path (24 points)
    // -------------------------------------------------------------------------
    println!("\n=== 2. Band Structure along Γ→K→M→Γ ===");

    let n_seg = 8_usize; // points per segment
    let n_total = 3 * n_seg; // 24 total

    // Build path: Γ→K (0..8), K→M (8..16), M→Γ (16..24)
    let mut k_path: Vec<(f64, f64)> = Vec::with_capacity(n_total);
    let mut k_labels: Vec<&str> = Vec::with_capacity(n_total);

    // Segment Γ→K
    for i in 0..n_seg {
        let t = i as f64 / n_seg as f64;
        let kx = k_gamma.0 + t * (k_k.0 - k_gamma.0);
        let ky = k_gamma.1 + t * (k_k.1 - k_gamma.1);
        k_path.push((kx, ky));
        k_labels.push(if i == 0 { "Γ" } else { "" });
    }

    // Segment K→M
    for i in 0..n_seg {
        let t = i as f64 / n_seg as f64;
        let kx = k_k.0 + t * (k_m.0 - k_k.0);
        let ky = k_k.1 + t * (k_m.1 - k_k.1);
        k_path.push((kx, ky));
        k_labels.push(if i == 0 { "K" } else { "" });
    }

    // Segment M→Γ
    for i in 0..n_seg {
        let t = i as f64 / n_seg as f64;
        let kx = k_m.0 + t * (k_gamma.0 - k_m.0);
        let ky = k_m.1 + t * (k_gamma.1 - k_m.1);
        k_path.push((kx, ky));
        k_labels.push(if i == 0 { "M" } else { "" });
    }

    println!(
        "  {:>8}  {:>6}  {:>12}  {:>12}  {:>12}",
        "k_idx", "Label", "E_1 [meV]", "E_2 [meV]", "Gap [meV]"
    );
    println!("  {}", "-".repeat(58));

    for (idx, (&k, label)) in k_path.iter().zip(k_labels.iter()).enumerate() {
        let (evals, _) = model
            .diagonalize(k)
            .expect("Failed to diagonalize along path");
        let gap = evals[1] - evals[0];
        println!(
            "  {:>8}  {:>6}  {:>12.6}  {:>12.6}  {:>12.6}",
            idx, label, evals[0], evals[1], gap
        );
    }

    // -------------------------------------------------------------------------
    // 3. Berry curvature on 20×20 BZ grid for band 0
    // -------------------------------------------------------------------------
    println!("\n=== 3. Berry Curvature on 20×20 BZ Grid (Band 0) ===");

    let nx_bc = 20_usize;
    let ny_bc = 20_usize;

    let berry = BerryCurvature::new(&model);
    let grid = berry
        .compute_grid(nx_bc, ny_bc, 0)
        .expect("Failed to compute Berry curvature grid");

    let mut max_omega = f64::NEG_INFINITY;
    let mut min_omega = f64::INFINITY;
    let mut sum_sq = 0.0_f64;
    let mut sum_omega = 0.0_f64;
    let n_pts = (nx_bc * ny_bc) as f64;

    for row in &grid {
        for &omega in row {
            if omega > max_omega {
                max_omega = omega;
            }
            if omega < min_omega {
                min_omega = omega;
            }
            sum_sq += omega * omega;
            sum_omega += omega;
        }
    }
    let rms_omega = (sum_sq / n_pts).sqrt();
    let mean_omega = sum_omega / n_pts;

    println!("  Grid size        : {}×{}", nx_bc, ny_bc);
    let abs_max = max_omega.abs().max(min_omega.abs());
    println!("  max |Ω|(k)       : {:.6}  (rad²·a²)", abs_max);
    println!("  max  Ω (k)       : {:.6}", max_omega);
    println!("  min  Ω (k)       : {:.6}", min_omega);
    println!("  Mean Ω           : {:.6}", mean_omega);
    println!("  RMS  Ω           : {:.6}", rms_omega);

    // BZ integral ≈ Chern number
    let bz_integral = berry
        .integrate_brillouin(nx_bc, ny_bc, 0)
        .expect("Failed to integrate Berry curvature");
    println!(
        "  BZ integral (≈C) : {:.4}  (should be ≈ integer Chern number)",
        bz_integral
    );

    // -------------------------------------------------------------------------
    // 4. Chern numbers via Fukui-Hatsugai method
    // -------------------------------------------------------------------------
    println!("\n=== 4. Chern Numbers via Fukui-Hatsugai Method (20×20 grid) ===");

    let nx_cn = 20_usize;
    let ny_cn = 20_usize;

    let chern_calc = ChernNumber::new(&model);

    println!("  {:>6}  {:>12}", "Band", "Chern number");
    println!("  {}", "-".repeat(22));

    let mut chern_sum = 0_i32;
    for band in 0..model.n_bands() {
        let c = chern_calc
            .compute(band, nx_cn, ny_cn)
            .expect("Failed to compute Chern number");
        println!("  {:>6}  {:>12}", band, c);
        chern_sum += c;
    }
    println!(
        "  {:>6}  {:>12}  (Nielsen-Ninomiya sum — must be 0)",
        "Total", chern_sum
    );

    // -------------------------------------------------------------------------
    // 5. Topological phase transition: sweep DMI 0 → 0.6 meV (7 points)
    // -------------------------------------------------------------------------
    println!("\n=== 5. Topological Phase Transition: Gap and Chern Number vs DMI ===");
    println!(
        "  {:>12}  {:>14}  {:>12}  {:>12}",
        "DMI [meV]", "Gap [meV]", "Chern_0", "Phase"
    );
    println!("  {}", "-".repeat(58));

    let dmi_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6_f64];
    for &dmi_sweep in &dmi_values {
        let m_sweep = MagnonBandModel::honeycomb_haldane(j_nn, j_nnn, dmi_sweep, h_ext)
            .expect("Failed to construct sweep model");

        // Minimum direct gap over 20×20 BZ grid
        let gap = m_sweep
            .band_gap(0, 1, 20, 20)
            .expect("Failed to compute band gap");

        // Chern number (or trivial if gap too small)
        let chern_phase = if gap < 1e-3 {
            // Near-gapless: skip Chern (unreliable)
            "     -- (gapless)"
        } else {
            let c_sw = ChernNumber::new(&m_sweep)
                .compute(0, 20, 20)
                .expect("Failed to compute Chern for sweep");
            if c_sw != 0 {
                "TOPOLOGICAL"
            } else {
                "TRIVIAL"
            }
        };

        // Chern value for numeric column
        let c_val: Option<i32> = if gap >= 1e-3 {
            Some(
                ChernNumber::new(&m_sweep)
                    .compute(0, 20, 20)
                    .expect("Chern recompute failed"),
            )
        } else {
            None
        };

        let c_str = match c_val {
            Some(c) => format!("{:>12}", c),
            None => format!("{:>12}", "N/A"),
        };

        println!(
            "  {:>12.3}  {:>14.6}  {}  {:>12}",
            dmi_sweep, gap, c_str, chern_phase
        );
    }

    // -------------------------------------------------------------------------
    // 6. Edge modes on a 20-cell strip
    // -------------------------------------------------------------------------
    println!("\n=== 6. Edge Modes on 20-Cell Strip ===");

    let strip_width = 20_usize;
    let edge = EdgeModes::new(&model, strip_width).expect("Failed to create EdgeModes calculator");

    // Determine the strip gap from the strip's own spectrum at kx = π/2.
    // The strip approximation shifts absolute energies slightly from the bulk.
    // We scan the strip spectrum at mid-BZ to find the largest internal gap.
    let modes_mid = edge
        .solve_strip(pi / 2.0)
        .expect("Failed to solve strip at mid BZ");
    let mut strip_evals: Vec<f64> = modes_mid.iter().map(|m| m.frequency).collect();
    strip_evals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Find the largest consecutive gap in the strip spectrum
    let (gap_min, gap_max) = {
        let mut best_gap = 0.0_f64;
        let mut best_lo = 0.0_f64;
        let mut best_hi = 0.0_f64;
        for w in strip_evals.windows(2) {
            let g = w[1] - w[0];
            if g > best_gap {
                best_gap = g;
                best_lo = w[0];
                best_hi = w[1];
            }
        }
        (best_lo, best_hi)
    };
    let gap_center = (gap_min + gap_max) / 2.0;

    // Bulk K-point gap for reference
    let (evals_k, _) = model.diagonalize(k_k).expect("Diagonalize at K failed");

    println!(
        "  Strip width       : {} unit cells ({} modes per kx)",
        strip_width,
        strip_width * model.n_bands()
    );
    println!(
        "  Strip gap center  : {:.6} meV  (from strip spectrum at kx=π/2)",
        gap_center
    );
    println!("  Strip gap window  : [{:.6}, {:.6}] meV", gap_min, gap_max);
    println!("  Bulk gap at K     : {:.6} meV", evals_k[1] - evals_k[0]);
    println!();

    let n_kx = 11_usize;
    let disp = edge
        .dispersion_curve(0.0, pi, n_kx)
        .expect("Failed to compute dispersion curve");

    println!(
        "  {:>8}  {:>14}  {:>14}  {:>14}  {:>14}",
        "kx", "E_edge_min", "E_edge_max", "n_in_gap", "Edge type"
    );
    println!("  {}", "-".repeat(75));

    for (kx, evals) in &disp {
        // Collect in-gap eigenvalues
        let in_gap: Vec<f64> = evals
            .iter()
            .copied()
            .filter(|&e| e >= gap_min && e <= gap_max)
            .collect();

        let n_in_gap = in_gap.len();

        // Edge type summary for this kx
        let edge_modes_at_kx = edge
            .find_in_gap_modes(*kx, gap_min, gap_max)
            .expect("Failed to find in-gap modes");

        let edge_type_str = if edge_modes_at_kx.is_empty() {
            "Bulk only".to_string()
        } else {
            let has_top = edge_modes_at_kx
                .iter()
                .any(|m| m.edge == spintronics::topomagnon::EdgeSide::Top);
            let has_bot = edge_modes_at_kx
                .iter()
                .any(|m| m.edge == spintronics::topomagnon::EdgeSide::Bottom);
            match (has_top, has_bot) {
                (true, true) => "Top+Bottom".to_string(),
                (true, false) => "Top edge".to_string(),
                (false, true) => "Bottom edge".to_string(),
                (false, false) => "In-gap/Bulk".to_string(),
            }
        };

        let e_min_str = if in_gap.is_empty() {
            "      ---".to_string()
        } else {
            format!(
                "{:.6}",
                in_gap.iter().copied().fold(f64::INFINITY, f64::min)
            )
        };
        let e_max_str = if in_gap.is_empty() {
            "      ---".to_string()
        } else {
            format!(
                "{:.6}",
                in_gap.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            )
        };

        println!(
            "  {:>8.4}  {:>14}  {:>14}  {:>14}  {:>14}",
            kx, e_min_str, e_max_str, n_in_gap, edge_type_str
        );
    }

    // -------------------------------------------------------------------------
    // Summary
    // -------------------------------------------------------------------------
    println!("\n=== Summary ===");
    println!("Topological magnon Haldane model on honeycomb lattice:");
    println!(
        "  - J_NN = {:.2} meV, DMI = {:.2} meV → topological phase with |C| = 1",
        j_nn, dmi
    );
    println!(
        "  - Berry curvature peaks at K/K' points (RMS Ω = {:.4})",
        rms_omega
    );
    println!(
        "  - Chern number C_0 = {} (lower band), C_1 = {} (upper band), sum = {}",
        chern_calc.compute(0, nx_cn, ny_cn).expect("C_0 failed"),
        chern_calc.compute(1, nx_cn, ny_cn).expect("C_1 failed"),
        chern_sum
    );
    println!("  - Phase transition: topological gap opens for DMI > critical value");
    let bulk_gap_k = evals_k[1] - evals_k[0];
    println!(
        "  - Chiral edge modes appear in the bulk gap ({:.4} meV wide at K)",
        bulk_gap_k
    );
    println!("  - Magnon Hall effect and non-zero κ_xy arise from non-trivial Chern number");

    Ok(())
}
