//! Graded (Inhomogeneous) Interfaces: Linear, Exponential, and Error-Function
//! Material-Parameter Profiles
//!
//! **Difficulty**: ⭐⭐ Intermediate
//! **Category**: Disordered Magnets / Multilayers
//! **Physics**: Compositionally graded interfaces, interdiffusion broadening,
//!              Boltzmann-Matano diffusion profiles, interfacial roughness
//!
//! Real magnetic interfaces (e.g. sputtered or MBE-grown multilayers) rarely
//! switch abruptly between two bulk material parameters. Interdiffusion
//! during growth and post-growth annealing smears the transition over a
//! finite width, and interfacial roughness makes that transition position
//! fluctuate laterally from site to site. This example demonstrates
//! `GradedInterface`, which models a scalar material parameter (Ms, exchange
//! stiffness A, or anisotropy K) varying smoothly across such a transition
//! region.
//!
//! We explore:
//!
//! 1. The three grading laws (Linear, Exponential, ErrorFunction) applied to
//!    a saturation-magnetisation step between two ferromagnets
//! 2. Endpoint saturation and the exact point-antisymmetry ("barycenter")
//!    identity shared by all three laws
//! 3. Exchange-stiffness grading across a diffuse interface
//! 4. Roughness-induced lateral jitter of the transition position
//! 5. Robust rejection of invalid parameters
//!
//! # References
//! - J. Crank, "The Mathematics of Diffusion", 2nd ed., Oxford (1975) --
//!   error-function solutions of the diffusion equation
//! - C. Matano, Jpn. J. Phys. 8, 109 (1933) -- Boltzmann-Matano analysis of
//!   interdiffusion profiles
//! - E. E. Fullerton et al., Phys. Rev. B 45, 9292 (1992) -- interdiffusion
//!   at sputtered magnetic multilayer interfaces

use spintronics::material::disorder::{GradedInterface, GradingLaw};

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Graded (Inhomogeneous) Interfaces ===\n");

    // =========================================================================
    // 1. Three grading laws for an Ms step
    // =========================================================================
    println!("=== 1. Saturation Magnetisation Ms Across a Graded Interface ===");

    let ms_left = 4.0e5_f64; // A/m (e.g. a softer ferromagnet)
    let ms_right = 8.0e5_f64; // A/m (e.g. Permalloy-like)
    let width = 2.0e-9_f64; // 2 nm transition
    let center = 0.0_f64;

    let linear = GradedInterface::generate(
        GradingLaw::Linear,
        center,
        width,
        ms_left,
        ms_right,
        1,
        0.0,
        1,
    )?;
    let exponential = GradedInterface::generate(
        GradingLaw::Exponential,
        center,
        width,
        ms_left,
        ms_right,
        1,
        0.0,
        2,
    )?;
    let error_fn = GradedInterface::generate(
        GradingLaw::ErrorFunction,
        center,
        width,
        ms_left,
        ms_right,
        1,
        0.0,
        3,
    )?;

    println!(
        "\nMs(left) = {:.3e} A/m,  Ms(right) = {:.3e} A/m,  width = {:.1} nm\n",
        ms_left,
        ms_right,
        width * 1e9
    );
    println!(
        "  {:>10}   {:>14}   {:>14}   {:>14}",
        "depth [nm]", "Linear", "Exponential", "ErrorFunction"
    );
    println!("  {}", "-".repeat(58));
    for i in -10..=10 {
        let z = i as f64 * width; // depth in units of one width, -10w..+10w
        println!(
            "  {:>10.2}   {:>14.4e}   {:>14.4e}   {:>14.4e}",
            z * 1e9,
            linear.value_at(z),
            exponential.value_at(z),
            error_fn.value_at(z)
        );
    }

    // =========================================================================
    // 2. Endpoint saturation and the barycenter identity
    // =========================================================================
    println!("\n=== 2. Endpoint Saturation and the Barycenter Identity ===");
    let far = 10.0 * width;
    println!(
        "\nAt |depth| = 10 x width = {:.1} nm, all three laws saturate to the bulk values:",
        far * 1e9
    );
    for (name, profile) in [
        ("Linear", &linear),
        ("Exponential", &exponential),
        ("ErrorFunction", &error_fn),
    ] {
        println!(
            "  {:<14}  value({:+.1} nm) = {:.4e} A/m,  value({:+.1} nm) = {:.4e} A/m",
            name,
            -far * 1e9,
            profile.value_at(-far),
            far * 1e9,
            profile.value_at(far)
        );
    }

    println!(
        "\nBarycenter (mean of the two endpoints) = {:.4e} A/m",
        linear.barycenter()
    );
    println!("Point-antisymmetry check: value(c+d) + value(c-d) == Ms_left + Ms_right exactly:");
    for &d_nm in &[0.5_f64, 1.0, 3.0, 8.0] {
        let d = d_nm * 1e-9;
        let sum = error_fn.value_at(center + d) + error_fn.value_at(center - d);
        println!(
            "  d = {:>4.1} nm:  sum = {:.6e} A/m  (expected {:.6e} A/m)",
            d_nm,
            sum,
            ms_left + ms_right
        );
    }

    // =========================================================================
    // 3. Exchange stiffness A graded across a diffuse interface
    // =========================================================================
    println!("\n=== 3. Exchange Stiffness A Across a Diffuse Interface ===");
    let a_left = 1.0e-11_f64; // J/m
    let a_right = 2.0e-11_f64; // J/m
    let a_width = 1.0e-9_f64; // sharper, 1 nm interdiffusion length
    let a_profile = GradedInterface::generate(
        GradingLaw::ErrorFunction,
        5.0e-9, // interface offset from the origin
        a_width,
        a_left,
        a_right,
        1,
        0.0,
        4,
    )?;
    println!(
        "\nA(left) = {:.3e} J/m,  A(right) = {:.3e} J/m,  interface at z0 = 5.0 nm",
        a_left, a_right
    );
    for z_nm in [0.0_f64, 3.0, 5.0, 7.0, 10.0] {
        let z = z_nm * 1e-9;
        println!(
            "  A(z = {:>5.1} nm) = {:.4e} J/m",
            z_nm,
            a_profile.value_at(z)
        );
    }

    // =========================================================================
    // 4. Roughness-induced lateral jitter of the transition position
    // =========================================================================
    println!("\n=== 4. Roughness-Induced Lateral Jitter ===");
    let n_sites = 2000_usize;
    let jitter_rms = 0.3e-9; // 0.3 nm interfacial roughness
    let rough = GradedInterface::generate(
        GradingLaw::ErrorFunction,
        center,
        width,
        ms_left,
        ms_right,
        n_sites,
        jitter_rms,
        99,
    )?;

    let mean_sq: f64 = rough.center_offsets.iter().map(|o| o * o).sum::<f64>() / n_sites as f64;
    println!(
        "\nRequested jitter RMS = {:.3} nm,  realized RMS = {:.3} nm  (N = {} sites)",
        jitter_rms * 1e9,
        mean_sq.sqrt() * 1e9,
        n_sites
    );

    // Local Ms at fixed depth z=0 varies from site to site because the local
    // interface center is jittered even though the nominal center sits at z=0.
    let z_probe = 0.0;
    let mut local_values: Vec<f64> = (0..10).map(|s| rough.value_at_site(z_probe, s)).collect();
    println!("\nMs at z = 0 for the first 10 lateral sites (local center jitters per site):");
    for (s, v) in local_values.iter().enumerate() {
        println!("  site {:>2}: Ms = {:.4e} A/m", s, v);
    }
    local_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    println!(
        "  range across first 10 sites: [{:.4e}, {:.4e}] A/m",
        local_values[0],
        local_values[local_values.len() - 1]
    );

    // =========================================================================
    // 5. Invalid-parameter rejection
    // =========================================================================
    println!("\n=== 5. Invalid-Parameter Rejection ===");
    let bad_width = GradedInterface::generate(
        GradingLaw::Linear,
        0.0,
        -1.0e-9, // negative width: invalid
        ms_left,
        ms_right,
        1,
        0.0,
        1,
    );
    println!(
        "  negative width           -> {}",
        if bad_width.is_err() {
            "rejected (Err)"
        } else {
            "accepted"
        }
    );
    let bad_endpoint = GradedInterface::generate(
        GradingLaw::Linear,
        0.0,
        width,
        f64::NAN, // non-finite endpoint: invalid
        ms_right,
        1,
        0.0,
        1,
    );
    println!(
        "  non-finite endpoint      -> {}",
        if bad_endpoint.is_err() {
            "rejected (Err)"
        } else {
            "accepted"
        }
    );
    let bad_sites =
        GradedInterface::generate(GradingLaw::Linear, 0.0, width, ms_left, ms_right, 0, 0.0, 1);
    println!(
        "  zero lateral sites       -> {}",
        if bad_sites.is_err() {
            "rejected (Err)"
        } else {
            "accepted"
        }
    );

    // =========================================================================
    // Summary
    // =========================================================================
    println!("\n=== Summary ===");
    println!("Graded interface physics demonstrated:");
    println!(
        "  - Ms step {:.2e} -> {:.2e} A/m graded over {:.1} nm (3 laws)",
        ms_left,
        ms_right,
        width * 1e9
    );
    println!("  - Barycenter identity value(c+d)+value(c-d) = Ms_left+Ms_right holds exactly");
    println!(
        "  - A step {:.2e} -> {:.2e} J/m graded over {:.1} nm interdiffusion length",
        a_left,
        a_right,
        a_width * 1e9
    );
    println!(
        "  - Interfacial roughness jitter: requested {:.2} nm, realized {:.2} nm RMS ({} sites)",
        jitter_rms * 1e9,
        mean_sq.sqrt() * 1e9,
        n_sites
    );
    println!("  - Invalid parameters (negative width, non-finite endpoints, zero sites) rejected");

    Ok(())
}
