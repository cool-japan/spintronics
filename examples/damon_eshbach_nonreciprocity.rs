//! Damon-Eshbach Non-Reciprocity in YIG Thin Films
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Spin Wave Physics
//! **Physics**: Non-reciprocal spin wave propagation in in-plane magnetized thin films
//!
//! ## Background
//!
//! The Damon-Eshbach (DE) geometry describes spin waves propagating perpendicular to
//! the in-plane equilibrium magnetization in a ferromagnetic thin film. These surface-
//! localized modes exhibit a hallmark non-reciprocity: waves propagating in the +k and
//! -k directions have different frequencies. This asymmetry arises because the dynamic
//! dipolar field breaks the mirror symmetry between the two film surfaces. The upper
//! surface carries the +k mode while the lower surface carries the -k mode, and they
//! reside at slightly different frequencies. This non-reciprocity enables directional
//! spin wave devices, magnonic diodes, and non-reciprocal phase shifters.
//!
//! The full Kalinikos-Slavin dispersion is used:
//!   ω² = (ω_H + ω_M λ_ex k²)(ω_H + ω_M λ_ex k² + ω_M F_kd sin²φ + ω_M(1−F_kd))
//!
//! where F_kd = 1 − (1 − exp(−kd))/(kd) is the dipolar propagation factor.
//!
//! ## References
//! - Damon & Eshbach, J. Phys. Chem. Solids 19, 308 (1961)
//! - Demidov et al., PRL 96, 097202 (2006)
//! - Hurben & Patton, J. Magn. Magn. Mater. 163, 39 (1996)

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("============================================================");
    println!("  DamonEshbach Nonreciprocity (in-plane magnetized film)");
    println!("============================================================");

    // -------------------------------------------------------------------------
    // Section 1: YIG 300 nm film at μ₀H = 100 mT
    // H_ext = B / μ₀ = 0.1 / (4π×10⁻⁷) ≈ 79,577 A/m
    // -------------------------------------------------------------------------
    println!("\n--- Section 1: YIG 300 nm film at μ₀H = 100 mT ---\n");

    let h_ext_100mt: f64 = 0.1 / MU_0; // 100 mT → A/m ≈ 79.58 kA/m
    let ms_yig: f64 = 1.4e5; // A/m
    let a_ex_yig: f64 = 3.5e-12; // J/m
    let alpha_yig: f64 = 3e-5;
    let d_yig: f64 = 300e-9; // 300 nm

    let yig = DamonEshbachDetailed::new(d_yig, h_ext_100mt, ms_yig, a_ex_yig, alpha_yig)?;

    let omega_h = yig.omega_h();
    let omega_m = yig.omega_m();

    println!(
        "  ω_H / 2π = {:>8.3} GHz   (Larmor frequency)",
        omega_h / (2.0 * std::f64::consts::PI * 1e9)
    );
    println!(
        "  ω_M / 2π = {:>8.3} GHz   (magnetization frequency)",
        omega_m / (2.0 * std::f64::consts::PI * 1e9)
    );
    println!("  thickness = {:>8.1} nm", d_yig * 1e9);
    println!("  α (Gilbert) = {:>8.2e}", alpha_yig);

    // Dispersion table
    let k_values: &[f64] = &[1e5, 1e6, 5e6, 1e7, 5e7];
    println!();
    println!(
        "  {:>12}  {:>10}  {:>10}  {:>12}",
        "k (rad/m)", "f(+k) GHz", "f(-k) GHz", "Δf/f_avg (%)"
    );
    println!("  {}", "-".repeat(50));

    for &k in k_values {
        let omega_pos = yig.dispersion_omega(k, std::f64::consts::FRAC_PI_2);
        // For the anti-reciprocal direction we use the nonreciprocity shift
        let delta_omega = yig.nonreciprocity(k);
        let omega_neg = (omega_pos - delta_omega).max(0.0);

        let f_pos = omega_pos / (2.0 * std::f64::consts::PI * 1e9);
        let f_neg = omega_neg / (2.0 * std::f64::consts::PI * 1e9);
        let f_avg = 0.5 * (f_pos + f_neg);
        let delta_f_pct = if f_avg > 1e-6 {
            100.0 * (f_pos - f_neg).abs() / f_avg
        } else {
            0.0
        };

        println!(
            "  {:>12.3e}  {:>10.4}  {:>10.4}  {:>12.3}",
            k, f_pos, f_neg, delta_f_pct
        );
    }

    // Group velocity at k = 5e6 rad/m
    let k_ref: f64 = 5e6;
    let vg = yig.group_velocity(k_ref, std::f64::consts::FRAC_PI_2);
    println!(
        "\n  Group velocity at k = {:.0e} rad/m : {:>8.1} m/s",
        k_ref, vg
    );

    // Surface localization length at k = 5e6 rad/m
    let loc_len = yig.surface_localization_length(k_ref);
    println!(
        "  Surface localization length       : {:>8.1} nm",
        loc_len * 1e9
    );

    // Propagation length at k = 5e6 rad/m
    let prop_len = yig.propagation_length(k_ref);
    println!(
        "  Propagation length (1/e amplitude): {:>8.1} μm",
        prop_len * 1e6
    );

    // -------------------------------------------------------------------------
    // Section 2: Material comparison at k = 5e6 rad/m
    // -------------------------------------------------------------------------
    println!("\n--- Section 2: Material comparison at k = 5×10⁶ rad/m ---\n");

    // YIG 300 nm – use the custom 100 mT variant above
    let yig_300nm = &yig;

    // Permalloy 20 nm – zero field (h_ext = 0 by default)
    let permalloy = DamonEshbachDetailed::permalloy_thin();

    // CoFeB strip – zero field (h_ext = 0 by default)
    let cofeb = DamonEshbachDetailed::cofeb_strip();

    println!(
        "  {:>16}  {:>10}  {:>14}  {:>16}",
        "Material", "f (GHz)", "Δf/f_avg (%)", "PropLen (μm)"
    );
    println!("  {}", "-".repeat(62));

    struct MatEntry<'a> {
        name: &'a str,
        model: &'a DamonEshbachDetailed,
    }

    let materials = [
        MatEntry {
            name: "YIG 300 nm",
            model: yig_300nm,
        },
        MatEntry {
            name: "Permalloy 20nm",
            model: &permalloy,
        },
        MatEntry {
            name: "CoFeB strip",
            model: &cofeb,
        },
    ];

    for mat in &materials {
        let k = k_ref;
        let omega = mat.model.dispersion_omega(k, std::f64::consts::FRAC_PI_2);
        let delta_omega = mat.model.nonreciprocity(k);
        let omega_neg = (omega - delta_omega).max(0.0);

        let f_ghz = omega / (2.0 * std::f64::consts::PI * 1e9);
        let f_neg_ghz = omega_neg / (2.0 * std::f64::consts::PI * 1e9);
        let f_avg = 0.5 * (f_ghz + f_neg_ghz);
        let delta_pct = if f_avg > 1e-6 {
            100.0 * (f_ghz - f_neg_ghz).abs() / f_avg
        } else {
            0.0
        };

        let prop_len_um = mat.model.propagation_length(k) * 1e6;

        println!(
            "  {:>16}  {:>10.4}  {:>14.3}  {:>16.2}",
            mat.name, f_ghz, delta_pct, prop_len_um
        );
    }

    println!("\n============================================================");
    println!("  Done. Non-reciprocity grows with k (stronger surface");
    println!("  localization) and is largest for low-damping YIG.");
    println!("============================================================\n");

    Ok(())
}
