//! Random Anisotropy Disorder: Imry-Ma Length and LLG Dynamics
//!
//! **Difficulty**: ⭐⭐⭐ Advanced
//! **Category**: Disordered Magnets
//! **Physics**: Random anisotropy, Imry-Ma correlation length, Harris criterion,
//!              LLG with disorder field
//!
//! This example demonstrates the physics of random anisotropy disorder in magnetic
//! materials.  We explore:
//!
//! 1. Nanocrystalline and custom Gaussian random-anisotropy models
//! 2. Spatial correlation functions and Imry-Ma magnetic correlation length
//! 3. Harris criterion relevance for varying sample sizes
//! 4. Random easy-axis statistics and per-site anisotropy energetics
//! 5. Simplified LLG dynamics driven by the random anisotropy effective field,
//!    showing dephasing/demagnetisation of an initially saturated spin ensemble
//! 6. Magnetisation autocorrelation decay across disorder realisations
//!
//! # References
//! - Y. Imry & S. K. Ma, Phys. Rev. Lett. 35, 1399 (1975)
//! - A. B. Harris, J. Phys. C 7, 1671 (1974)
//! - E. M. Chudnovsky & R. A. Serota, Phys. Rev. B 33, 251 (1986)
//! - G. Herzer, IEEE Trans. Magn. 26, 1397 (1990)

use spintronics::material::random_anisotropy::{RandomAnisotropy, RandomAnisotropyDistribution};
use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Random Anisotropy Disorder: Imry-Ma Length and LLG Dynamics ===\n");

    // =========================================================================
    // 1. Model construction
    // =========================================================================
    println!("=== 1. Random Anisotropy Models ===");

    // Nanocrystalline preset (no seed argument — preset fixes seed = 42)
    let nano = RandomAnisotropy::nanocrystalline();
    println!("\nNanocrystalline preset:");
    println!("  k_mean              = {:.3e} J/m³", nano.k_mean);
    println!("  k_std               = {:.3e} J/m³", nano.k_std);
    println!("  k_rms               = {:.3e} J/m³", nano.k_rms());
    println!(
        "  correlation_length  = {:.1} nm",
        nano.correlation_length * 1e9
    );
    println!("  distribution        = {:?}", nano.distribution);
    println!("  seed                = {}", nano.seed);

    // Custom model with Gaussian distribution
    let k_mean = 5e3_f64; // J/m³
    let k_std = 5e2_f64; // J/m³
    let xi_a = 8e-9_f64; // 8 nm grain size
    let mut model = RandomAnisotropy::new(
        k_mean,
        k_std,
        xi_a,
        RandomAnisotropyDistribution::Gaussian,
        1234,
    )?;

    println!("\nCustom Gaussian model:");
    println!("  k_mean              = {:.3e} J/m³", model.k_mean);
    println!("  k_std               = {:.3e} J/m³", model.k_std);
    println!("  k_rms               = {:.3e} J/m³", model.k_rms());
    println!(
        "  correlation_length  = {:.1} nm",
        model.correlation_length * 1e9
    );
    println!("  distribution        = {:?}", model.distribution);

    // =========================================================================
    // 2. Spatial correlation and Imry-Ma length
    // =========================================================================
    println!("\n=== 2. Spatial Correlation and Imry-Ma Length ===");

    let exchange_a = 1e-11_f64; // J/m (typical soft ferromagnet)

    println!(
        "\nCorrelation function C(r) = exp(−r / ξ_a)  [ξ_a = {:.0} nm]:",
        xi_a * 1e9
    );
    println!("  {:>10}   {:>12}", "r [nm]", "C(r)");
    println!("  {}", "-".repeat(26));
    for &r_nm in &[0.0_f64, 1.0, 2.0, 5.0, 10.0] {
        let r = r_nm * 1e-9;
        println!("  {:>10.1}   {:>12.6}", r_nm, model.correlation_function(r));
    }

    let xi_m = model.imry_ma_correlation_length(exchange_a);
    println!("\nExchange stiffness A   = {:.2e} J/m", exchange_a);
    println!(
        "Imry-Ma length  ξ_M    = {:.3e} m  ({:.1} μm)",
        xi_m,
        xi_m * 1e6
    );
    println!("Grain length    ξ_a    = {:.1} nm", xi_a * 1e9);
    println!(
        "Ratio ξ_M / ξ_a       = {:.1e}  (>> 1 for weak disorder / strong exchange)",
        xi_m / xi_a
    );

    // =========================================================================
    // 3. Harris criterion
    // =========================================================================
    println!("\n=== 3. Harris Criterion Relevance ===");
    println!(
        "  Disorder is relevant when sample size L > ξ_a = {:.0} nm",
        xi_a * 1e9
    );
    println!(
        "\n  {:>18}   {:>15}   {:>10}",
        "Sample size (L)", "L / ξ_a", "Relevant?"
    );
    println!("  {}", "-".repeat(50));
    for &n_units in &[10_usize, 50, 100, 500, 1000] {
        let l = n_units as f64 * xi_a;
        let relevant = model.harris_relevant(l);
        println!(
            "  {:>14} ξ_a   {:>15.1}   {:>10}",
            n_units,
            l / xi_a,
            if relevant { "YES" } else { "no" }
        );
    }

    // =========================================================================
    // 4. Random axes, strengths, and energetics
    // =========================================================================
    println!("\n=== 4. Random Axes, Strengths, and Energetics ===");

    let n_sites = 500_usize;
    let axes = model.generate_axes(n_sites);
    let strengths = model.generate_strengths(n_sites);

    println!("\nFirst 5 random easy-axis directions (should be unit vectors):");
    println!(
        "  {:>4}   {:>10}  {:>10}  {:>10}   {:>10}",
        "site", "n_x", "n_y", "n_z", "|n|"
    );
    println!("  {}", "-".repeat(52));
    for (i, n) in axes.iter().take(5).enumerate() {
        let mag = n.magnitude();
        println!(
            "  {:>4}   {:>10.6}  {:>10.6}  {:>10.6}   {:>10.8}",
            i, n.x, n.y, n.z, mag
        );
    }

    let mean_mag: f64 = axes.iter().map(|v| v.magnitude()).sum::<f64>() / n_sites as f64;
    let mean_k: f64 = strengths.iter().sum::<f64>() / n_sites as f64;
    let var_k: f64 =
        strengths.iter().map(|&k| (k - mean_k).powi(2)).sum::<f64>() / (n_sites - 1) as f64;
    let std_k = var_k.sqrt();
    println!("\nAxis and strength statistics (N = {}):", n_sites);
    println!("  mean |axis|   = {:.8} (should be 1.0)", mean_mag);
    println!(
        "  mean K        = {:.4e} J/m³  (expected {:.3e})",
        mean_k, k_mean
    );
    println!(
        "  std K         = {:.4e} J/m³  (expected {:.3e})",
        std_k, k_std
    );

    // Anisotropy energy for m pointing along x, y, z
    let m_x_all: Vec<Vector3<f64>> = (0..n_sites).map(|_| Vector3::unit_x()).collect();
    let m_y_all: Vec<Vector3<f64>> = (0..n_sites).map(|_| Vector3::unit_y()).collect();
    let m_z_all: Vec<Vector3<f64>> = (0..n_sites).map(|_| Vector3::unit_z()).collect();

    let e_x = model.anisotropy_energy(&m_x_all, &axes, &strengths);
    let e_y = model.anisotropy_energy(&m_y_all, &axes, &strengths);
    let e_z = model.anisotropy_energy(&m_z_all, &axes, &strengths);
    println!("\nAnisotropy energy for uniform magnetisation:");
    println!("  E(m ∥ x) = {:.4e} J/m³", e_x);
    println!("  E(m ∥ y) = {:.4e} J/m³", e_y);
    println!("  E(m ∥ z) = {:.4e} J/m³", e_z);

    // Effective anisotropy field magnitude for m = [1,0,0] at each site
    let ms = 8.0e5_f64; // A/m, typical permalloy
    let h_fields = model.effective_field(&m_x_all, &axes, &strengths, ms)?;
    let h_norms: Vec<f64> = h_fields.iter().map(|h| h.magnitude()).collect();
    let mean_h = h_norms.iter().sum::<f64>() / n_sites as f64;
    let max_h = h_norms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_h = h_norms.iter().cloned().fold(f64::INFINITY, f64::min);

    println!(
        "\nEffective anisotropy field |H_eff| for m = [1,0,0]  (Ms = {:.1e} A/m):",
        ms
    );
    println!("  First 5 site magnitudes:");
    for (i, &h) in h_norms.iter().take(5).enumerate() {
        println!("    site {:>3}: |H_eff| = {:.4e} A/m", i, h);
    }
    println!("  Mean |H_eff|  = {:.4e} A/m", mean_h);
    println!("  Min  |H_eff|  = {:.4e} A/m", min_h);
    println!("  Max  |H_eff|  = {:.4e} A/m", max_h);

    // =========================================================================
    // 5. Simplified LLG dynamics with random anisotropy field
    // =========================================================================
    println!("\n=== 5. LLG Dynamics with Random Anisotropy Field ===");
    println!("N = 100 spin sites, H_ext = 0, α = 0.01, dt = 1e-12 s, 50 steps");

    let n_llg = 100_usize;
    let gamma = GAMMA; // gyromagnetic ratio [rad/(s·T)] — from prelude constants
                       // Convert A/m field to Tesla: μ₀ H (T)
    let mu_0 = MU_0;
    let alpha = 0.01_f64;
    let dt = 1e-12_f64; // s

    // Fresh axes/strengths for LLG ensemble (reset to a reproducible sub-seed)
    let mut llg_model = RandomAnisotropy::new(
        k_mean,
        k_std,
        xi_a,
        RandomAnisotropyDistribution::Gaussian,
        9999,
    )?;
    let llg_axes = llg_model.generate_axes(n_llg);
    let llg_strengths = llg_model.generate_strengths(n_llg);

    // Initial state: all spins along x
    let mut spins: Vec<Vector3<f64>> = (0..n_llg).map(|_| Vector3::unit_x()).collect();

    println!(
        "\n  {:>6}   {:>14}   {:>14}   {:>14}",
        "step", "m̄_x", "m̄_y", "m̄_z"
    );
    println!("  {}", "-".repeat(54));

    let print_steps: std::collections::HashSet<usize> = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 49]
        .iter()
        .cloned()
        .collect();

    for step in 0..50_usize {
        // Print mean magnetisation at selected steps
        if print_steps.contains(&step) {
            let mx: f64 = spins.iter().map(|s| s.x).sum::<f64>() / n_llg as f64;
            let my: f64 = spins.iter().map(|s| s.y).sum::<f64>() / n_llg as f64;
            let mz: f64 = spins.iter().map(|s| s.z).sum::<f64>() / n_llg as f64;
            println!("  {:>6}   {:>14.8}   {:>14.8}   {:>14.8}", step, mx, my, mz);
        }

        // Compute anisotropy fields [A/m] → convert to Tesla for LLG
        let h_eff_am = llg_model.effective_field(&spins, &llg_axes, &llg_strengths, ms)?;

        // Euler step for each spin
        for i in 0..n_llg {
            let m = spins[i];
            // H_eff in Tesla
            let h = Vector3::new(
                mu_0 * h_eff_am[i].x,
                mu_0 * h_eff_am[i].y,
                mu_0 * h_eff_am[i].z,
            );
            // dm/dt = -γ (m × H) - α γ (m × (m × H))
            let m_cross_h = m.cross(&h);
            let m_cross_m_cross_h = m.cross(&m_cross_h);
            let dm = Vector3::new(
                -gamma * m_cross_h.x - alpha * gamma * m_cross_m_cross_h.x,
                -gamma * m_cross_h.y - alpha * gamma * m_cross_m_cross_h.y,
                -gamma * m_cross_h.z - alpha * gamma * m_cross_m_cross_h.z,
            );
            let m_new = Vector3::new(m.x + dt * dm.x, m.y + dt * dm.y, m.z + dt * dm.z);
            // Re-normalise to keep on unit sphere
            spins[i] = m_new.normalize();
        }
    }

    // Final mean magnetisation
    let mx_final: f64 = spins.iter().map(|s| s.x).sum::<f64>() / n_llg as f64;
    let my_final: f64 = spins.iter().map(|s| s.y).sum::<f64>() / n_llg as f64;
    let mz_final: f64 = spins.iter().map(|s| s.z).sum::<f64>() / n_llg as f64;
    let m_total = (mx_final * mx_final + my_final * my_final + mz_final * mz_final).sqrt();
    println!("\n  Final mean magnetisation:");
    println!("    m̄_x = {:.8}", mx_final);
    println!("    m̄_y = {:.8}", my_final);
    println!("    m̄_z = {:.8}", mz_final);
    println!(
        "    |m̄| = {:.8}  (< 1 shows dephasing by disorder)",
        m_total
    );

    // =========================================================================
    // 6. Magnetisation autocorrelation decay across disorder realisations
    // =========================================================================
    println!("\n=== 6. Magnetisation Autocorrelation Across Realisations ===");
    println!("  10 independent disorder seeds, 10 LLG steps each (N = 100 sites)");
    println!(
        "\n  {:>12}   {:>14}   {:>10}",
        "realisation", "m_x after 10 steps", "seed"
    );
    println!("  {}", "-".repeat(44));

    let n_real = 10_usize;
    let n_steps = 10_usize;
    let mut mx_values = Vec::with_capacity(n_real);

    for real_idx in 0..n_real {
        let seed_r = 7777_u64 + real_idx as u64 * 31;
        let mut rm = RandomAnisotropy::new(
            k_mean,
            k_std,
            xi_a,
            RandomAnisotropyDistribution::Gaussian,
            seed_r,
        )?;
        let r_axes = rm.generate_axes(n_llg);
        let r_strengths = rm.generate_strengths(n_llg);

        let mut s: Vec<Vector3<f64>> = (0..n_llg).map(|_| Vector3::unit_x()).collect();
        for _ in 0..n_steps {
            let hf = rm.effective_field(&s, &r_axes, &r_strengths, ms)?;
            for i in 0..n_llg {
                let m = s[i];
                let h = Vector3::new(mu_0 * hf[i].x, mu_0 * hf[i].y, mu_0 * hf[i].z);
                let m_cross_h = m.cross(&h);
                let m_cross_m_cross_h = m.cross(&m_cross_h);
                let dm = Vector3::new(
                    -gamma * m_cross_h.x - alpha * gamma * m_cross_m_cross_h.x,
                    -gamma * m_cross_h.y - alpha * gamma * m_cross_m_cross_h.y,
                    -gamma * m_cross_h.z - alpha * gamma * m_cross_m_cross_h.z,
                );
                let m_new = Vector3::new(m.x + dt * dm.x, m.y + dt * dm.y, m.z + dt * dm.z);
                s[i] = m_new.normalize();
            }
        }
        let mx_r: f64 = s.iter().map(|sv| sv.x).sum::<f64>() / n_llg as f64;
        mx_values.push(mx_r);
        println!("  {:>12}   {:>14.8}   {:>10}", real_idx + 1, mx_r, seed_r);
    }

    let mx_mean_real: f64 = mx_values.iter().sum::<f64>() / n_real as f64;
    let mx_var_real: f64 = mx_values
        .iter()
        .map(|&v| (v - mx_mean_real).powi(2))
        .sum::<f64>()
        / (n_real - 1) as f64;
    println!(
        "\n  Ensemble average ⟨m_x⟩ after 10 steps = {:.6}",
        mx_mean_real
    );
    println!(
        "  Ensemble std dev σ(m_x)               = {:.6}",
        mx_var_real.sqrt()
    );
    println!("  (Deviation from 1 confirms disorder-driven dephasing)");

    // =========================================================================
    // Summary
    // =========================================================================
    println!("\n=== Summary ===");
    println!("Random anisotropy disorder physics demonstrated:");
    println!("  - Grain correlation length ξ_a  = {:.0} nm", xi_a * 1e9);
    println!(
        "  - RMS anisotropy            K_rms = {:.3e} J/m³",
        model.k_rms()
    );
    println!(
        "  - Imry-Ma magnetic length   ξ_M   = {:.2e} m = {:.1} μm",
        xi_m,
        xi_m * 1e6
    );
    println!(
        "  - Harris disorder relevant for L > {:.0} nm (grain size)",
        xi_a * 1e9
    );
    println!(
        "  - LLG dephasing: |m̄| = {:.4} after 50 steps  (started at 1.0)",
        m_total
    );
    println!(
        "  - Ensemble ⟨m_x⟩ = {:.4} ± {:.4} across 10 disorder realisations",
        mx_mean_real,
        mx_var_real.sqrt()
    );

    Ok(())
}
