//! Nonlinear Magnon Physics: Suhl Instability and Parametric Amplification
//!
//! **Difficulty**: ⭐⭐⭐⭐
//! **Category**: Nonlinear Spin Dynamics / Parametric Magnon Processes
//! **Physics**: Four-magnon scattering; Suhl spin-wave instability threshold;
//! parametric amplification of spin waves; nonlinear FMR linewidth broadening;
//! bistability and foldover; Kalinikos-Slavin dispersion.
//!
//! ## Background
//!
//! When the microwave pump field exceeds a threshold h_th, energy is transferred
//! from the uniform precession (k=0) to spin-wave pairs (k, -k) via four-magnon
//! scattering — the **Suhl first-order spin-wave instability**:
//!
//! ```text
//! h_th = μ₀·Δω_k / |T_kk|
//! ```
//!
//! where Δω_k = ω_k - ω_pump/2 is the detuning and T_kk the four-magnon
//! vertex.  Above threshold, parametric gain σ = √[(|T|·h)² - Δω²] grows
//! exponentially.  This is exploited in:
//! - YIG parametric amplifiers (ultra-low-noise microwave amplification)
//! - Magnon Bose-Einstein condensation (spin-wave driven BEC)
//! - Nonlinear FMR spectroscopy (linewidth broadening, foldover)
//!
//! ## References
//!
//! - H. Suhl, "The Theory of Ferromagnetic Resonance at High Signal Powers",
//!   *J. Phys. Chem. Solids* **1**, 209 (1957).
//! - B. Kalinikos & A. Slavin, "Theory of dipole-exchange spin wave spectrum",
//!   *J. Phys. C* **19**, 7013 (1986).
//! - V. E. Zakharov et al., "Spin-wave turbulence beyond the parametric
//!   excitation threshold", *Sov. Phys. Usp.* **17**, 896 (1975).
//! - A. Chumak et al., "Magnon spintronics", *Nat. Phys.* **11**, 453 (2015).

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Nonlinear Magnon Physics: Suhl Instability & Parametric Amp.");
    println!("=================================================================\n");

    // -------------------------------------------------------------------------
    // 1. YIG material parameters
    //    FourMagnonScattering::yig() expects h_ext in Tesla.
    //    B₀ = μ₀ H = 4π×10⁻⁷ × 159.2 kA/m ≈ 0.200 T → f_FMR ≈ 5.6 GHz
    // -------------------------------------------------------------------------
    let b_ext = 0.200_f64; // [T]  applied flux density
    let yig = FourMagnonScattering::yig(b_ext);

    let f_fmr_ghz = yig.magnon_frequency(0.0, 0.0) / (2.0 * std::f64::consts::PI * 1.0e9);
    println!("--- YIG Material: FMR Operating Point ---");
    println!("  B_ext = {:.3} T  (H ≈ 159 kA/m)", b_ext);
    println!("  f_FMR = {:.4} GHz", f_fmr_ghz);
    println!();

    // -------------------------------------------------------------------------
    // 2. Magnon dispersion (Kalinikos-Slavin)
    // -------------------------------------------------------------------------
    println!("--- Magnon Dispersion at H_ext = 159.2 kA/m ---");
    println!(
        "  {:>14}  {:>12}  {:>12}  {:>12}",
        "k (rad/m)", "θ=0° (GHz)", "θ=45° (GHz)", "θ=90° (GHz)"
    );
    let k_vals: &[f64] = &[0.0, 1.0e5, 5.0e5, 1.0e6, 5.0e6, 1.0e7];
    for &k in k_vals {
        let f0 = yig.magnon_frequency(k, 0.0) / (2.0 * std::f64::consts::PI * 1.0e9);
        let f45 = yig.magnon_frequency(k, std::f64::consts::FRAC_PI_4)
            / (2.0 * std::f64::consts::PI * 1.0e9);
        let f90 = yig.magnon_frequency(k, std::f64::consts::FRAC_PI_2)
            / (2.0 * std::f64::consts::PI * 1.0e9);
        println!("  {:>14.3e}  {:>12.4}  {:>12.4}  {:>12.4}", k, f0, f45, f90);
    }
    println!();

    // -------------------------------------------------------------------------
    // 3. Suhl instability threshold vs wavevector
    // -------------------------------------------------------------------------
    // h_th returned in T; convert to A/m and mT for display
    println!("--- Suhl Threshold Field h_th vs Wavevector ---");
    println!(
        "  {:>14}  {:>12}  {:>14}  {:>14}",
        "k (rad/m)", "θ_k", "h_th (mT)", "h_th (A/m)"
    );
    let k_suhl_vals = [1.0e4, 1.0e5, 5.0e5, 1.0e6, 5.0e6, 1.0e7];
    let theta_suhl = std::f64::consts::FRAC_PI_4; // θ = 45° (MSSW geometry)
    for &k in &k_suhl_vals {
        let h_th_t = yig.suhl_threshold_field(k, theta_suhl); // [T]
        let h_th_mt = h_th_t * 1.0e3; // T → mT
        let h_th_am = h_th_t / (4.0 * std::f64::consts::PI * 1.0e-7); // T → A/m
        println!(
            "  {:>14.3e}  {:>12.1}°  {:>14.4}  {:>14.2}",
            k,
            theta_suhl.to_degrees(),
            h_th_mt,
            h_th_am
        );
    }

    // Optimal wavevector for pump at f_pump = f_FMR / 2 (degenerate case)
    let f_pump_hz = yig.magnon_frequency(0.0, 0.0) * 0.999; // near-degenerate
    let k_opt = yig.optimal_k_for_suhl(f_pump_hz);
    println!(
        "\n  Optimal k* for f_pump = {:.4} GHz → k* = {:.4e} rad/m",
        f_pump_hz / (2.0 * std::f64::consts::PI * 1.0e9),
        k_opt
    );
    println!();

    // -------------------------------------------------------------------------
    // 4. Parametric growth rate above threshold
    // -------------------------------------------------------------------------
    println!("--- Parametric Growth Rate σ vs Pump Field ---");
    let k_param = k_opt.max(1.0e5);
    let theta_param = std::f64::consts::FRAC_PI_4;
    let h_th_t = yig.suhl_threshold_field(k_param, theta_param); // [T]
    let h_th_mt = h_th_t * 1.0e3; // mT
    println!(
        "  Threshold field h_th = {:.4} mT at k={:.2e} rad/m",
        h_th_mt, k_param
    );
    println!(
        "  {:>14}  {:>14}  {:>16}  {:>12}",
        "h_pump (T)", "h/h_th", "σ (rad/s)", "σ (GHz)"
    );
    let h_ratios = [0.5, 0.9, 1.0, 1.2, 1.5, 2.0, 5.0];
    for &ratio in &h_ratios {
        let h = h_th_t * ratio;
        let sigma = yig.parametric_growth_rate(h, k_param, theta_param);
        let sigma_ghz = sigma / (2.0 * std::f64::consts::PI * 1.0e9);
        let in_str = if ratio > 1.0 { "ABOVE" } else { "below" };
        println!(
            "  {:>14.4e}  {:>14.2}  {:>16.4e}  {:>12.6}  {}",
            h, ratio, sigma, sigma_ghz, in_str
        );
    }
    println!();

    // -------------------------------------------------------------------------
    // 5. Instability check at fixed pump field
    // -------------------------------------------------------------------------
    let h_applied = h_th_t * 1.5;
    let is_unstable = yig.instability_is_present(h_applied, f_pump_hz);
    println!("--- Instability Check ---");
    println!("  h_pump = 1.5·h_th = {:.4e} T", h_applied);
    println!("  Suhl instability active? → {}", is_unstable);
    println!();

    // -------------------------------------------------------------------------
    // 6. Parametric amplification (degenerate signal-idler)
    //    degenerate_from_yig(omega_pump, alpha, Ms): creates degenerate amp
    //    where signal = idler = omega_pump/2 (both in rad/s).
    // -------------------------------------------------------------------------
    println!("--- Parametric Amplification (Degenerate YIG Amplifier) ---");
    let omega_pump = f_pump_hz; // omega in rad/s (f_pump_hz is actually omega — see magnon_frequency)
    let alpha_pa = 3.0e-5_f64;
    let ms_pa = 1.4e5_f64;
    // Pass pump angular freq omega_pump so signal/idler = omega_pump/2.
    // The preset's pump field is the physics-derived threshold h_th = α·ω_p/(γ·Ms),
    // so it sits exactly at marginal stability. Operate it at twice-critical
    // (ξ = h_p/h_th = 2) to demonstrate above-threshold parametric gain.
    let pa = ParametricAmplification::degenerate_from_yig(omega_pump, alpha_pa, ms_pa)
        .with_supercriticality(2.0)
        .expect("twice-critical pump is a valid operating point");

    let h_thresh = pa.threshold_pump_field();
    let gain_coeff = pa.gain_coefficient();
    let is_above = pa.is_above_threshold();
    // omega_signal = omega_pump/2, displayed in GHz
    let f_sig_ghz = pa.idler_frequency() / (2.0 * std::f64::consts::PI * 1.0e9);
    println!("  Signal/Idler freq = {:.4} GHz", f_sig_ghz);
    println!(
        "  Threshold h_th    = {:.4e} T  ({:.2e} A/m)",
        h_thresh,
        h_thresh / (4.0 * std::f64::consts::PI * 1.0e-7)
    );
    println!("  Gain coefficient  = {:.4e} rad/s", gain_coeff);
    println!("  Above threshold?  = {}", is_above);
    println!(
        "  Phase mismatch    = {:.4e} rad/s",
        pa.phase_mismatch(omega_pump)
    );

    // Show gain vs round-trips with a pump 1.5× above threshold
    println!("\n  Parametric growth: gain coefficient scales as σ = √((gH)²-Δω²)");
    println!("  (Gain diverges exponentially above threshold — shown as finite σ [rad/s])");
    println!(
        "  {:>10}  {:>14}  {:>14}  {:>14}",
        "h/h_th", "g·H (rad/s)", "Δω (rad/s)", "σ (krad/s)"
    );
    let h_th_pa = h_thresh;
    let gamma_pa = alpha_pa * omega_pump / 2.0; // linewidth at signal freq
    let coupling_pa = pa.gain_coefficient(); // gain at current pump
    let _ = coupling_pa; // already displayed
    for &ratio in &[0.5, 0.9, 1.0, 1.5, 2.0, 5.0] {
        let h_pump = h_th_pa * ratio;
        // coupling = GAMMA*Ms/2, drive = coupling*h_pump
        let coupling = GAMMA.abs() * ms_pa / 2.0;
        let drive = coupling * h_pump;
        let disc = drive * drive - gamma_pa * gamma_pa;
        let sigma = if disc > 0.0 { disc.sqrt() } else { 0.0 };
        println!(
            "  {:>10.2}  {:>14.4e}  {:>14.4e}  {:>14.4}",
            ratio,
            drive,
            gamma_pa,
            sigma * 1.0e-3
        );
    }
    println!();

    // -------------------------------------------------------------------------
    // 7. Nonlinear FMR linewidth and bistability
    // -------------------------------------------------------------------------
    println!("--- Nonlinear FMR Linewidth and Bistability ---");
    let nlw = NonlinearFmrLinewidth::from_yig(b_ext);
    let lw_linear = nlw.linear_linewidth_field();
    let fold_field = nlw.foldover_field();
    let bist_power = nlw.bistability_threshold_power();

    // linear_linewidth_field returns [T]; convert to mT and Oe for display
    let lw_mt = lw_linear * 1.0e3;
    let lw_oe = lw_linear * 1.0e4; // 1 T = 1e4 Oe
    println!(
        "  Linear FMR linewidth ΔH₀      = {:.4} mT  ({:.4} Oe)",
        lw_mt, lw_oe
    );
    let fold_mt = fold_field * 1.0e3;
    println!("  Foldover field H_fold          = {:.4} mT", fold_mt);
    println!("  Bistability threshold P_bist   = {:.4e} W", bist_power);

    // Effective linewidth vs pump amplitude (result in T, display in mT)
    println!("\n  Nonlinear linewidth ΔH(h_pump) at several pump levels:");
    println!(
        "  {:>14}  {:>14}  {:>14}",
        "h_pump/h_lin", "ΔH_NL (mT)", "Broadening"
    );
    let pump_ratios = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0];
    for &pr in &pump_ratios {
        let lw_nl = nlw.nonlinear_linewidth_field(pr) * 1.0e3; // T → mT
        let lw_linear_mt = lw_linear * 1.0e3;
        let broaden = if lw_linear_mt > 0.0 {
            lw_nl / lw_linear_mt
        } else {
            1.0
        };
        println!("  {:>14.1}  {:>14.4}  {:>14.4}×", pr, lw_nl, broaden);
    }

    // Power saturation (Suhl-Weiss compression)
    println!("\n  Power saturation factor vs input power:");
    println!("  {:>14}  {:>14}", "P_in (mW)", "η_sat");
    let powers_mw = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0];
    for &p_mw in &powers_mw {
        let eta = nlw.power_saturation_factor(p_mw * 1.0e-3);
        println!("  {:>14.3}  {:>14.6}", p_mw, eta);
    }

    // -------------------------------------------------------------------------
    // 8. Free-function: global instability power threshold
    // -------------------------------------------------------------------------
    let gamma = GAMMA.abs();
    let p_inst = suhl_spin_wave_instability_power(1.4e5, 3.0e-5, b_ext, gamma);
    println!(
        "\n--- Global Suhl Instability Power ---\n  P_inst (YIG, B={:.3} T) = {:.4e} W",
        b_ext, p_inst
    );

    // -------------------------------------------------------------------------
    // 9. Magnon-magnon interaction energy
    // -------------------------------------------------------------------------
    println!("\n--- Four-Magnon Interaction Energy ---");
    let t_kk = 1.0e-26_f64; // representative coupling J
    println!(
        "  {:>10}  {:>10}  {:>20}",
        "n1 (magnons)", "n2 (magnons)", "E_interaction (J)"
    );
    for n1 in [0.0, 1.0, 5.0, 10.0].iter() {
        for n2 in [0.0, 1.0, 5.0].iter() {
            let e_int = magnon_magnon_interaction_energy(*n1, *n2, t_kk);
            println!("  {:>10.0}  {:>10.0}  {:>20.4e}", n1, n2, e_int);
        }
    }

    println!("\n=================================================================");
    println!("  Summary: YIG nonlinear magnon physics demonstrates the Suhl");
    println!("  threshold h_th, parametric gain above threshold, and linewidth");
    println!("  broadening driven by four-magnon scattering at large amplitudes.");
    println!("=================================================================\n");

    Ok(())
}
