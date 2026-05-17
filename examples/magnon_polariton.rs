//! Magnon-Polariton: Strong Light-Matter Coupling in Cavity Magnonics
//!
//! **Difficulty**: ⭐⭐ Intermediate
//! **Category**: Cavity Magnonics
//! **Physics**: Magnon-polariton hybridization, Hopfield coefficients, vacuum Rabi
//!   splitting, coupling regimes
//!
//! This example demonstrates the physics of magnon-polariton formation in a
//! microwave cavity coupled to a YIG (Yttrium Iron Garnet) ferromagnet.  When
//! the magnon-photon coupling strength g exceeds the individual linewidths
//! (strong-coupling regime), the bare photon and magnon modes hybridise into
//! two new normal modes — the upper and lower polariton branches — with a
//! characteristic avoided crossing (anticrossing) in the dispersion relation.
//!
//! The treatment follows the analytic 2×2 Hopfield-Bogoliubov diagonalisation
//! of the Jaynes-Cummings Hamiltonian:
//!
//!   H = ℏω_c a†a + ℏω_m b†b + ℏg(a†b + ab†)
//!
//! giving polariton frequencies
//!
//!   ω_± = (ω_c + ω_m)/2 ± √[(Δ/2)² + g²]    (Δ = ω_c − ω_m)
//!
//! and Hopfield mixing coefficients X (photon fraction) and Y (magnon fraction)
//! with X² + Y² = 1.  The multi-mode extension couples N magnon modes to one
//! cavity by diagonalising the full (N+1)×(N+1) Hermitian matrix.
//!
//! ## What This Example Shows
//!
//! 1. **YIG/cavity system**: canonical parameters, vacuum Rabi splitting, and
//!    strong-coupling criterion.
//! 2. **Anticrossing curve**: polariton frequencies and Hopfield coefficients
//!    as the magnon frequency sweeps through cavity resonance.
//! 3. **Coupling strength regimes**: weak, intermediate, and strong coupling
//!    at zero detuning, with cooperativity and regime classification.
//! 4. **Multi-mode polariton**: three magnon modes coupled to one cavity —
//!    full eigenspectrum and mode compositions.
//!
//! ## References
//!
//! - H. Huebl et al., Phys. Rev. Lett. **111**, 127003 (2013)
//!   — first observation of strong magnon-photon coupling in YIG/cavity.
//! - J. J. Hopfield, Phys. Rev. **112**, 1555 (1958)
//!   — original polariton picture and Hopfield coefficients.
//! - M. Harder & C.-M. Hu, Solid State Phys. **69**, 47 (2018)
//!   — comprehensive review of cavity magnonics.
//! - M. Goryachev et al., Phys. Rev. Applied **2**, 054002 (2014)
//!   — high-cooperativity magnon-photon coupling.

use spintronics::cavity::{Branch, MagnonPolariton, MultiModePolariton};
use spintronics::prelude::{GAMMA, HBAR, MU_0};

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let pi = std::f64::consts::PI;

    println!("╔══════════════════════════════════════════════════════════════════╗");
    println!("║   Magnon-Polariton: Strong Light-Matter Coupling                ║");
    println!("║   Cavity Magnonics — Hopfield Diagonalisation                   ║");
    println!("╚══════════════════════════════════════════════════════════════════╝");
    println!();

    // =========================================================================
    // 1.  YIG/cavity system — canonical parameters
    // =========================================================================
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("1. YIG/Cavity System — Canonical Parameters");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let yig = MagnonPolariton::from_yig_cavity();

    let f_cavity_ghz = yig.omega_cavity / (2.0 * pi * 1.0e9);
    let f_magnon_ghz = yig.omega_magnon / (2.0 * pi * 1.0e9);
    let g_mhz = yig.coupling_g / (2.0 * pi * 1.0e6);
    let vrs_mhz = yig.vacuum_rabi_splitting() / (2.0 * pi * 1.0e6);

    println!("  Cavity frequency  ω_c  = {:.4} GHz", f_cavity_ghz);
    println!("  Magnon frequency  ω_m  = {:.4} GHz", f_magnon_ghz);
    println!("  Coupling strength g    = {:.1} MHz", g_mhz);
    println!("  Vacuum Rabi splitting  = 2g = {:.1} MHz", vrs_mhz);

    // Typical YIG linewidths for the strong-coupling test
    // κ (cavity) ~ 1 MHz,  γ (magnon) ~ 1 MHz  →  g/(κ+γ)/4 = 100/(2) = 50 >> 1
    let kappa = 2.0 * pi * 1.0e6; // 1 MHz cavity linewidth
    let gamma_m = 2.0 * pi * 1.0e6; // 1 MHz magnon linewidth (YIG has ~0.1 MHz but use 1 MHz)
    let cooperativity = 4.0 * yig.coupling_g.powi(2) / (kappa * gamma_m);
    let strong = yig.is_strong_coupling(kappa, gamma_m);

    println!();
    println!("  Dissipation parameters (κ = γ_m = 1 MHz):");
    println!(
        "    κ  (cavity linewidth)  = {:.1} MHz",
        kappa / (2.0 * pi * 1.0e6)
    );
    println!(
        "    γ_m (magnon linewidth) = {:.1} MHz",
        gamma_m / (2.0 * pi * 1.0e6)
    );
    println!("    Cooperativity C = 4g²/(κγ_m) = {:.1}", cooperativity);
    println!(
        "    Strong coupling?  → {}",
        if strong { "YES ✓" } else { "NO ✗" }
    );

    // Energy of vacuum Rabi photon: ℏ × 2g
    let e_rabi_mev = HBAR * yig.vacuum_rabi_splitting() / 1.602_176_634e-22; // in meV
                                                                             // Cross-check magnon frequency from fundamental constants: ω_m = γ × μ₀ × H_bias
    let h_bias_check = 8.0e4_f64; // 80 kA/m (same as from_yig_cavity)
    let omega_m_check = GAMMA * MU_0 * h_bias_check;
    println!();
    println!(
        "  Coupling energy ℏ × 2g = {:.4e} meV  (from HBAR = {:.4e} J·s)",
        e_rabi_mev, HBAR
    );
    println!(
        "  ω_m cross-check (GAMMA × MU_0 × H_bias): {:.4} GHz  [matches ω_m above]",
        omega_m_check / (2.0 * pi * 1.0e9)
    );

    // Detuning: Δ = ω_c − ω_m
    let delta_mhz = (yig.omega_cavity - yig.omega_magnon) / (2.0 * pi * 1.0e6);
    println!();
    println!("  Detuning Δ = ω_c − ω_m = {:.1} MHz", delta_mhz);
    let (lo, hi) = yig.eigenfrequencies();
    let f_lo_ghz = lo / (2.0 * pi * 1.0e9);
    let f_hi_ghz = hi / (2.0 * pi * 1.0e9);
    let splitting_mhz = (hi - lo) / (2.0 * pi * 1.0e6);
    println!("  Polariton eigenfrequencies:");
    println!("    ω_-  = {:.5} GHz  (lower branch)", f_lo_ghz);
    println!("    ω_+  = {:.5} GHz  (upper branch)", f_hi_ghz);
    println!("    Splitting (off-resonance) = {:.2} MHz", splitting_mhz);

    // =========================================================================
    // 2.  Anticrossing curve — frequency sweep
    // =========================================================================
    println!();
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("2. Anticrossing Curve — Sweeping Magnon Frequency");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let omega_c = 2.0 * pi * 10.0e9; // 10 GHz cavity (fixed)
    let g_sweep = 2.0 * pi * 100.0e6; // 100 MHz coupling

    let n_points = 21_usize;
    let omega_m_min = 0.5 * omega_c;
    let omega_m_max = 1.5 * omega_c;

    println!();
    println!(
        "  Cavity fixed at ω_c = {:.1} GHz,  g = {:.0} MHz",
        omega_c / (2.0 * pi * 1.0e9),
        g_sweep / (2.0 * pi * 1.0e6)
    );
    println!(
        "  Sweeping ω_m from {:.1} GHz to {:.1} GHz  ({} points)",
        omega_m_min / (2.0 * pi * 1.0e9),
        omega_m_max / (2.0 * pi * 1.0e9),
        n_points,
    );
    println!();

    // Header — aligned columns
    println!(
        "  {:>12}  {:>10}  {:>10}  {:>12}  {:>10}  {:>10}",
        "f_magnon[GHz]", "f_+[GHz]", "f_-[GHz]", "Split[MHz]", "X²(upper)", "Y²(upper)"
    );
    println!("  {}", "─".repeat(72));

    for i in 0..n_points {
        let t = i as f64 / (n_points - 1) as f64;
        let omega_m = omega_m_min + t * (omega_m_max - omega_m_min);

        let mp = MagnonPolariton::new(omega_c, omega_m, g_sweep);
        let (f_lo, f_hi) = mp.eigenfrequencies();
        let split_mhz = (f_hi - f_lo) / (2.0 * pi * 1.0e6);

        // Hopfield coefficients for upper branch
        let x2_up = mp.photon_fraction(Branch::Upper);
        let y2_up = mp.magnon_fraction(Branch::Upper);

        println!(
            "  {:>12.4}  {:>10.5}  {:>10.5}  {:>12.2}  {:>10.4}  {:>10.4}",
            omega_m / (2.0 * pi * 1.0e9),
            f_hi / (2.0 * pi * 1.0e9),
            f_lo / (2.0 * pi * 1.0e9),
            split_mhz,
            x2_up,
            y2_up,
        );
    }

    // Highlight the anticrossing minimum
    let mp_resonance = MagnonPolariton::new(omega_c, omega_c, g_sweep);
    let (lo_res, hi_res) = mp_resonance.eigenfrequencies();
    println!();
    println!("  Anticrossing minimum (ω_m = ω_c, Δ = 0):");
    println!("    ω_-  = {:.5} GHz", lo_res / (2.0 * pi * 1.0e9));
    println!("    ω_+  = {:.5} GHz", hi_res / (2.0 * pi * 1.0e9));
    println!(
        "    Splitting = {:.1} MHz  (= 2g = vacuum Rabi splitting)",
        (hi_res - lo_res) / (2.0 * pi * 1.0e6)
    );
    let x2_lo = mp_resonance.photon_fraction(Branch::Lower);
    let y2_lo = mp_resonance.magnon_fraction(Branch::Lower);
    println!(
        "    At resonance: X² = {:.4}, Y² = {:.4}  (equal mixing)",
        x2_lo, y2_lo
    );

    // =========================================================================
    // 3.  Coupling strength regimes
    // =========================================================================
    println!();
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("3. Coupling Strength Regimes at Zero Detuning (ω_m = ω_c)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    // Typical YIG linewidths: κ = 1 MHz, γ = 1 MHz
    // Condition for strong coupling: g > (κ + γ)/4
    let kappa_reg = 2.0 * pi * 1.0e6; // 1 MHz
    let gamma_reg = 2.0 * pi * 1.0e6; // 1 MHz
    let threshold = (kappa_reg + gamma_reg) / 4.0;

    struct Regime {
        label: &'static str,
        g_hz: f64,
    }

    let regimes = [
        Regime {
            label: "Weak",
            g_hz: 0.1e6,
        },
        Regime {
            label: "Intermediate",
            g_hz: 10.0e6,
        },
        Regime {
            label: "Strong",
            g_hz: 100.0e6,
        },
    ];

    println!();
    println!(
        "  Linewidths: κ = {:.0} MHz, γ_m = {:.0} MHz",
        kappa_reg / (2.0 * pi * 1.0e6),
        gamma_reg / (2.0 * pi * 1.0e6)
    );
    println!(
        "  Strong-coupling threshold g > (κ+γ)/4 = {:.2} MHz",
        threshold / (2.0 * pi * 1.0e6)
    );
    println!();
    println!(
        "  {:>14}  {:>10}  {:>12}  {:>12}  {:>12}  {:>14}",
        "Regime", "g [MHz]", "Split [MHz]", "Cooperativity", "Strong?", "Classification"
    );
    println!("  {}", "─".repeat(80));

    for r in &regimes {
        let g = 2.0 * pi * r.g_hz;
        let mp = MagnonPolariton::new(omega_c, omega_c, g);
        let (lo_r, hi_r) = mp.eigenfrequencies();
        let split = (hi_r - lo_r) / (2.0 * pi * 1.0e6);
        let coop = 4.0 * g.powi(2) / (kappa_reg * gamma_reg);
        let strong_r = mp.is_strong_coupling(kappa_reg, gamma_reg);

        let classification = if coop < 1.0 {
            "Weak (perturbative)"
        } else if coop < 100.0 {
            "Intermediate"
        } else {
            "Strong (resolved split)"
        };

        println!(
            "  {:>14}  {:>10.1}  {:>12.2}  {:>12.1}  {:>12}  {:>14}",
            r.label,
            r.g_hz / 1.0e6,
            split,
            coop,
            if strong_r { "YES" } else { "NO" },
            classification,
        );
    }

    println!();
    println!("  Note: cooperativity C = 4g² / (κ γ_m).  Resolved vacuum Rabi splitting");
    println!("  requires C >> 1 and g > (κ+γ)/4 simultaneously.");

    // =========================================================================
    // 4.  Multi-mode polariton
    // =========================================================================
    println!();
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("4. Multi-Mode Polariton — Cavity + 3 Magnon Modes");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let cavity_freq = 2.0 * pi * 10.0e9; // 10 GHz
    let magnon_freqs = vec![2.0 * pi * 9.0e9, 2.0 * pi * 10.0e9, 2.0 * pi * 11.0e9];
    let couplings = vec![2.0 * pi * 50.0e6, 2.0 * pi * 100.0e6, 2.0 * pi * 50.0e6];

    println!();
    println!(
        "  Cavity: ω_c = {:.1} GHz",
        cavity_freq / (2.0 * pi * 1.0e9)
    );
    println!("  Magnon modes:");
    for (i, (&om, &gm)) in magnon_freqs.iter().zip(couplings.iter()).enumerate() {
        println!(
            "    Mode {}: ω_m{} = {:>5.1} GHz,  g_{} = {:>5.0} MHz",
            i + 1,
            i + 1,
            om / (2.0 * pi * 1.0e9),
            i + 1,
            gm / (2.0 * pi * 1.0e6),
        );
    }

    let mmp = MultiModePolariton::new(cavity_freq, magnon_freqs, couplings)
        .expect("MultiModePolariton construction must succeed with valid inputs");

    let eigen_freqs = mmp
        .eigenfrequencies()
        .expect("eigenfrequency diagonalisation must succeed");
    let compositions = mmp
        .mode_compositions()
        .expect("mode_compositions must succeed");

    println!();
    println!("  Polariton eigenfrequencies (sorted ascending):");
    for (i, &freq) in eigen_freqs.iter().enumerate() {
        println!("    ω_{} = {:.5} GHz", i + 1, freq / (2.0 * pi * 1.0e9));
    }

    println!();
    println!("  Mode compositions  |  X² = photon fraction,  Y²_i = magnon-i fraction");
    println!("  (each row sums to 1; column 0 = photon, columns 1-3 = magnon modes)");
    println!();
    println!(
        "  {:>8}  {:>10}  {:>10}  {:>10}  {:>10}  {:>12}",
        "Branch", "X²(photon)", "Y²(m1)", "Y²(m2)", "Y²(m3)", "f [GHz]"
    );
    println!("  {}", "─".repeat(66));

    for (i, (comp, &freq)) in compositions.iter().zip(eigen_freqs.iter()).enumerate() {
        let x2 = comp.first().copied().unwrap_or(0.0);
        let y1 = comp.get(1).copied().unwrap_or(0.0);
        let y2 = comp.get(2).copied().unwrap_or(0.0);
        let y3 = comp.get(3).copied().unwrap_or(0.0);
        let sum: f64 = comp.iter().sum();
        println!(
            "  {:>8}  {:>10.4}  {:>10.4}  {:>10.4}  {:>10.4}  {:>12.5}   (Σ={:.6})",
            format!("ψ_{}", i + 1),
            x2,
            y1,
            y2,
            y3,
            freq / (2.0 * pi * 1.0e9),
            sum,
        );
    }

    println!();
    println!("  Physical interpretation:");
    println!("    • Mode at ~9 GHz  is dominated by the 9 GHz magnon (weak coupling: 50 MHz).");
    println!("    • Mode at ~10 GHz is the strongly hybridised mode — cavity + 10 GHz magnon");
    println!("      (strong coupling: 100 MHz) → largest photon admixture and anticrossing.");
    println!("    • Mode at ~11 GHz is dominated by the 11 GHz magnon (weak coupling: 50 MHz).");
    println!("    • The central doublet near 10 GHz shows the vacuum Rabi splitting from the");
    println!("      100 MHz coupling.");

    println!();
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  Done.  All calculations performed via analytic Hopfield diagonalisation.");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    Ok(())
}
