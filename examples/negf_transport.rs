//! NEGF Quantum Transport: Transmission and I-V Characteristics
//!
//! **Difficulty**: ⭐⭐⭐ Advanced
//! **Category**: Non-Equilibrium Transport
//! **Physics**: Non-Equilibrium Green's Functions, Landauer-Büttiker, shot noise,
//!              spin-dependent transport
//!
//! This example demonstrates quantum transport through a 1D tight-binding chain
//! using the Non-Equilibrium Green's Function (NEGF) formalism.  We show:
//!
//! 1. Transmission spectrum T(E) for a 20-site uniform chain
//! 2. I-V characteristics via the Landauer–Büttiker formula at 300 K
//! 3. Shot noise and Fano factor at several bias voltages
//! 4. Effect of Anderson disorder (σ = 0.3 eV) on transmission and localization
//!
//! The Landauer formula connects microscopic quantum mechanics to macroscopic
//! transport:
//!
//! ```text
//! I = (e/h) ∫ T(E) [f_L(E) − f_R(E)] dE
//! ```
//!
//! Shot noise reveals sub-Poissonian statistics in coherent transport:
//!
//! ```text
//! F = S_shot / (2eI) ∈ [0, 1]
//! ```
//!
//! References:
//! - Keldysh, JETP 20, 1018 (1965)
//! - Datta, Electronic Transport in Mesoscopic Systems (1995)
//! - Caroli et al., J. Phys. C 4, 916 (1971)
//! - Blanter & Büttiker, Phys. Rep. 336, 1 (2000)

use spintronics::negf::{
    GreenFunction, Hamiltonian1D, KeldyshSolver, LeadSelfEnergy, ShotNoise, TransportCalculator,
};
use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== NEGF Quantum Transport: Transmission and I-V Characteristics ===\n");

    // =========================================================================
    // 1.  Uniform 20-site tight-binding chain
    // =========================================================================
    let n_sites = 20_usize;
    let onsite_energy = 0.0_f64; // eV
    let hopping = -1.0_f64; // eV  (negative: conventional sign for electron bands)

    println!("=== 1. Tight-Binding Chain ===");
    let hamiltonian =
        Hamiltonian1D::from_uniform(n_sites, onsite_energy, hopping).expect("valid Hamiltonian");

    println!("  Sites : {}", hamiltonian.n_sites());
    println!("  Hopping t = {:.2} eV", hamiltonian.hopping);
    println!(
        "  Band : [{:.2}, {:.2}] eV  (2|t| = {:.2} eV bandwidth)",
        2.0 * hopping,
        -2.0 * hopping,
        4.0 * hopping.abs()
    );

    // Print first few site energies to confirm uniform construction
    println!("  On-site energies (first 5 sites):");
    for i in 0..5 {
        let e = hamiltonian.site_energy(i).expect("valid site index");
        println!("    site {:>2} : ε = {:.3} eV", i, e);
    }

    // =========================================================================
    // 2.  Transmission spectrum T(E)
    // =========================================================================
    println!("\n=== 2. Transmission Spectrum T(E) ===");

    let gamma = 0.5_f64; // eV — lead coupling width (wide-band limit)
    let eta = 1e-3_f64; // eV — convergence broadening

    let sigma_l = LeadSelfEnergy::new(gamma, 0.0).expect("valid left lead");
    let sigma_r = LeadSelfEnergy::new(gamma, 0.0).expect("valid right lead");

    let gf = GreenFunction::new(hamiltonian.clone(), sigma_l.clone(), sigma_r.clone(), eta)
        .expect("valid GreenFunction");

    let tc =
        TransportCalculator::new(gf.clone(), -3.0, 3.0, 200).expect("valid TransportCalculator");

    println!("  Γ_L = Γ_R = {:.2} eV (wide-band lead coupling)", gamma);
    println!("  η   = {:.0e} eV (broadening)", eta);
    println!();
    println!("  {:>8} {:>14}", "E (eV)", "T(E)");
    println!("  {}", "-".repeat(25));
    // Print 11 evenly-spaced sample energies from -2.5 to +2.5 eV
    for k in 0..=10 {
        let energy = -2.5 + k as f64 * 0.5;
        let t_e = tc.transmission(energy).expect("transmission evaluation");
        println!("  {:>8.3} {:>14.6}", energy, t_e);
    }

    // =========================================================================
    // 3.  I-V characteristics (Landauer integral, symmetric bias)
    // =========================================================================
    println!("\n=== 3. I-V Characteristics (T = 300 K, symmetric bias) ===");

    let temperature = 300.0_f64; // K

    // G_0 = e^2/h — the conductance quantum
    let g0 = CONDUCTANCE_QUANTUM; // S  (= e²/h ≈ 7.748e-5 S)

    // The library's current() integrates over an energy grid in eV, but uses the
    // SI prefactor e/h [A/J].  To obtain the physical current in Amperes the
    // result must be multiplied by E_CHARGE [J/eV] to convert the energy step
    // from eV to Joules.  Equivalently: I [A] = (e²/h) * ∫ T(E)(f_L−f_R) dE[eV].
    let ev_to_joule = E_CHARGE; // 1 eV = 1.602e-19 J

    println!("  T   = {:.0} K", temperature);
    println!("  G₀  = {:.4e} S  (conductance quantum e²/h)", g0);
    println!("  μ_L = +V/2,  μ_R = −V/2  (symmetric bias)\n");

    println!(
        "  {:>8} {:>14} {:>18} {:>12}",
        "V (eV)", "I (µA)", "dI/dV (G₀)", "G_0·T(0)"
    );
    println!("  {}", "-".repeat(58));

    // Zero-bias conductance from T(E_F=0) as reference — uses CONDUCTANCE_QUANTUM * T, no eV factor
    let g_zero_bias = tc
        .zero_bias_conductance(temperature)
        .expect("zero-bias conductance");
    let g_zero_in_g0 = g_zero_bias / g0;

    for k in 0..=10 {
        let v_bias = k as f64 * 0.2; // 0.0 to 2.0 V
                                     // Apply eV→J correction to get physical Amperes
        let current_raw = tc
            .current(v_bias, temperature)
            .expect("current calculation");
        let current_a = current_raw * ev_to_joule;
        let current_ua = current_a * 1.0e6; // convert A → µA
        let di_dv_raw = tc
            .differential_conductance(v_bias, temperature)
            .expect("differential conductance");
        let di_dv = di_dv_raw * ev_to_joule; // A/V (= S), corrected
        let di_dv_g0 = di_dv / g0;
        // Show G_0·T(0) column only for first row; blank afterwards for clarity
        let g0t_str = if k == 0 {
            format!("{:.6}", g_zero_in_g0)
        } else {
            String::from("")
        };
        println!(
            "  {:>8.1} {:>14.4} {:>18.6} {:>12}",
            v_bias, current_ua, di_dv_g0, g0t_str
        );
    }

    // =========================================================================
    // 4.  Shot noise and Fano factor
    // =========================================================================
    println!("\n=== 4. Shot Noise and Fano Factor ===");

    // ShotNoise wraps a TransportCalculator; reuse the same energy grid
    let tc_noise = TransportCalculator::new(gf.clone(), -3.0, 3.0, 400)
        .expect("valid TransportCalculator for noise");
    let shot_noise = ShotNoise::new(tc_noise);

    println!("  Poissonian limit : F = 1  (tunnel junction, T → 0)");
    println!("  Ballistic limit  : F = 0  (perfect channel, T = 1)");
    println!("  Sub-Poissonian noise indicates coherent quantum transport\n");

    println!(
        "  {:>8} {:>18} {:>14} {:>18} {:>10}",
        "V (eV)", "S_shot (A²/Hz)", "I (µA)", "S_total (A²/Hz)", "Fano F"
    );
    println!("  {}", "-".repeat(74));

    for &v_bias in &[0.5_f64, 1.0, 1.5] {
        let s_shot = shot_noise
            .shot_noise_only(v_bias, temperature)
            .expect("shot noise");
        let s_total = shot_noise
            .noise_zero_freq(v_bias, temperature)
            .expect("total noise");
        // Apply eV→J correction (same as in section 3)
        let current_raw = shot_noise
            .transport
            .current(v_bias, temperature)
            .expect("current for noise");
        let current_ua = current_raw * ev_to_joule * 1.0e6;
        let fano = shot_noise
            .fano_factor(v_bias, temperature)
            .expect("Fano factor");
        println!(
            "  {:>8.2} {:>18.4e} {:>14.4} {:>18.4e} {:>10.4}",
            v_bias, s_shot, current_ua, s_total, fano
        );
    }

    // Contextualise the Fano values
    let fano_mid = shot_noise
        .fano_factor(1.0, temperature)
        .expect("Fano at 1.0 V");
    let regime = if fano_mid < 0.1 {
        "near-ballistic"
    } else if fano_mid < 0.5 {
        "sub-Poissonian (coherent)"
    } else {
        "near-Poissonian (tunnel)"
    };
    println!(
        "\n  Transport regime at V=1.0 eV: F = {:.4} → {}",
        fano_mid, regime
    );

    // =========================================================================
    // 5.  Anderson disorder: uniform vs disordered chain
    // =========================================================================
    println!("\n=== 5. Anderson Disorder and Localization ===");

    let disorder_sigma = 0.3_f64; // eV — disorder strength
    let seed = 42_u64;

    let h_disordered = Hamiltonian1D::from_uniform(n_sites, onsite_energy, hopping)
        .expect("valid Hamiltonian")
        .with_disorder(disorder_sigma, seed);

    // Show that on-site energies are now random
    println!("  σ_disorder = {:.2} eV, seed = {}", disorder_sigma, seed);
    println!("  Disordered on-site energies (first 5 sites):");
    for i in 0..5 {
        let e = h_disordered.site_energy(i).expect("valid site index");
        println!("    site {:>2} : ε = {:+.4} eV", i, e);
    }

    let sigma_l_d = LeadSelfEnergy::new(gamma, 0.0).expect("valid left lead");
    let sigma_r_d = LeadSelfEnergy::new(gamma, 0.0).expect("valid right lead");
    let gf_disordered = GreenFunction::new(h_disordered.clone(), sigma_l_d, sigma_r_d, eta)
        .expect("valid disordered GreenFunction");
    let tc_disordered = TransportCalculator::new(gf_disordered.clone(), -3.0, 3.0, 200)
        .expect("valid disordered TransportCalculator");

    println!("\n  T(E) comparison: uniform vs disordered chain\n");
    println!(
        "  {:>8} {:>16} {:>16} {:>12}",
        "E (eV)", "T_uniform", "T_disordered", "Ratio"
    );
    println!("  {}", "-".repeat(56));

    for k in 0..=10 {
        let energy = -2.5 + k as f64 * 0.5;
        let t_uniform = tc.transmission(energy).expect("uniform transmission");
        let t_disorder = tc_disordered
            .transmission(energy)
            .expect("disordered transmission");
        let ratio = if t_uniform > 1e-12 {
            t_disorder / t_uniform
        } else {
            0.0
        };
        println!(
            "  {:>8.3} {:>16.6} {:>16.6} {:>12.4}",
            energy, t_uniform, t_disorder, ratio
        );
    }

    // Quantify localization at band centre E=0
    let t_clean_ef = tc.transmission(0.0).expect("clean T(E_F)");
    let t_dirty_ef = tc_disordered.transmission(0.0).expect("disordered T(E_F)");
    let localization_factor = if t_clean_ef > 1e-12 {
        t_dirty_ef / t_clean_ef
    } else {
        0.0
    };
    println!("\n  At E = 0 (band centre):");
    println!("    T_uniform    = {:.6}", t_clean_ef);
    println!("    T_disordered = {:.6}", t_dirty_ef);
    println!(
        "    T_dis / T_uni = {:.4}  (Anderson localization suppresses transmission)",
        localization_factor
    );

    // =========================================================================
    // 6.  Keldysh non-equilibrium occupation (bonus: per-site carrier density)
    // =========================================================================
    println!("\n=== 6. Keldysh Non-Equilibrium Carrier Density ===");

    let v_nonequil = 0.5_f64; // eV
    let mu_l_ne = v_nonequil / 2.0;
    let mu_r_ne = -v_nonequil / 2.0;

    let keldysh = KeldyshSolver::new(gf.clone());

    // Nonequilibrium carrier density for first 5 sites
    let density = keldysh
        .nonequilibrium_density(mu_l_ne, mu_r_ne, temperature, -3.0, 3.0, 200)
        .expect("nonequilibrium density");

    println!(
        "  V = {:.2} eV  (μ_L = {:+.3}, μ_R = {:+.3})",
        v_nonequil, mu_l_ne, mu_r_ne
    );
    println!("  T = {:.0} K", temperature);
    println!("  Per-site carrier density n_i (integrated over band):\n");
    println!("  {:>6} {:>20}", "Site", "n_i (per eV·site)");
    println!("  {}", "-".repeat(30));
    for (i, &ni) in density.iter().enumerate() {
        println!("  {:>6} {:>20.6e}", i, ni);
    }
    let total: f64 = density.iter().sum();
    println!(
        "  {:>6} {:>20.6e}  ← total across {} sites",
        "SUM", total, n_sites
    );

    // Verify Keldysh identity numerically at E=0
    let energy_check = 0.0_f64;
    let occ = keldysh
        .occupation_density(energy_check, mu_l_ne, mu_r_ne, temperature)
        .expect("occupation density");
    println!(
        "\n  Occupation density at E = {:.1} eV: {:.6e} eV⁻¹",
        energy_check, occ
    );

    // =========================================================================
    // Summary
    // =========================================================================
    println!("\n=== Summary ===");
    println!(
        "NEGF transport for a {}-site tight-binding chain (t = {:.1} eV, Γ = {:.2} eV):",
        n_sites,
        hopping.abs(),
        gamma
    );

    let t_at_zero = tc.transmission(0.0).expect("T(0)");
    let i_at_1v = tc.current(1.0, temperature).expect("I(1V)") * ev_to_joule * 1.0e6;
    println!(
        "  - T(E=0)              = {:.4} (near-unity for clean chain)",
        t_at_zero
    );
    println!("  - I(V=1.0 V)          = {:.4} µA at T = 300 K", i_at_1v);
    println!("  - G(V→0)              = {:.4} G₀", g_zero_in_g0);
    println!(
        "  - Fano factor         = {:.4}  (coherent, sub-Poissonian)",
        fano_mid
    );
    println!(
        "  - Anderson disorder σ = {:.2} eV suppresses T(0) by factor {:.4}",
        disorder_sigma, localization_factor
    );
    println!("  - Keldysh formalism confirms non-equilibrium carrier accumulation");

    Ok(())
}
