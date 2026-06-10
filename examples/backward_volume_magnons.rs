//! Backward Volume Magnetostatic Spin Waves with Negative Group Velocity
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Spin Wave Physics
//! **Physics**: BVMSW dispersion, backward propagation, and quantized thickness modes
//!
//! ## Background
//!
//! Backward Volume Magnetostatic Spin Waves (BVMSW) arise in ferromagnetic thin
//! films magnetized perpendicular to the film plane by a sufficiently strong field
//! H_perp > M_s. In-plane spin waves in this geometry exhibit a characteristic
//! backward dispersion: at small wavevectors the group velocity is negative, meaning
//! energy propagates anti-parallel to the wavevector. This exotic behaviour stems
//! from the dipolar interaction, which produces a negative restoring contribution
//! to the dispersion. At large k the exchange interaction takes over and the group
//! velocity becomes positive again — the crossover wavevector k* separates the
//! dipolar-dominated (backward) regime from the exchange-dominated (forward) regime.
//!
//! The dispersion relation (Kalinikos-Slavin, perpendicular geometry):
//!   ω² = (ω_H + ω_M λ_ex k²)(ω_H + ω_M λ_ex k² + ω_M (F_kd − 1))
//!
//! where ω_H = |γ| μ₀ (H_ext − M_s) is the effective Larmor frequency and
//! F_kd = (1 − exp(−kd)) / (kd) → 1 at k→0 and → 0 at k→∞.
//!
//! ## References
//! - Kalinikos & Slavin, J. Phys. C 19, 7013 (1986)
//! - Schneider et al., APL 92, 022505 (2008)
//! - Damon & Eshbach, J. Phys. Chem. Solids 19, 308 (1961)

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  Backward Volume Magnetostatic Spin Waves (BVMSW)");
    println!("=============================================================");

    // -------------------------------------------------------------------------
    // Section 1: YIG perpendicular field h_perp = 600 kA/m
    // μ₀h ≈ 754 mT, well above μ₀·Ms ≈ 176 mT for YIG (Ms = 1.4×10⁵ A/m)
    // -------------------------------------------------------------------------
    println!("\n--- Section 1: YIG 1 μm film, h_perp = 600 kA/m (μ₀h ≈ 754 mT) ---\n");

    let h_perp: f64 = 600e3; // A/m
    let bv = BackwardVolumeMSW::yig_yagi(h_perp)?;

    let omega_h_eff = bv.omega_h_eff();
    let omega_m = bv.omega_m();
    let h_eff_a_m = h_perp - bv.ms; // H_eff = H_ext - M_s

    println!(
        "  ω_H_eff / 2π = {:>8.3} GHz   (effective Larmor: H_ext - M_s)",
        omega_h_eff / (2.0 * std::f64::consts::PI * 1e9)
    );
    println!(
        "  ω_M   / 2π   = {:>8.3} GHz   (characteristic magnetization freq.)",
        omega_m / (2.0 * std::f64::consts::PI * 1e9)
    );
    println!(
        "  H_eff = H_ext - M_s = {:.3e} A/m = {:.1} kA/m",
        h_eff_a_m,
        h_eff_a_m * 1e-3
    );
    println!("  M_s   = {:.3e} A/m", bv.ms);
    println!("  d     = {:.1} μm", bv.thickness * 1e6);

    // Dispersion table
    let k_values: &[f64] = &[0.0, 1e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7];
    println!();
    println!(
        "  {:>12}  {:>10}  {:>12}  {:>8}",
        "k (rad/m)", "f (GHz)", "v_g (m/s)", "sign(v_g)"
    );
    println!("  {}", "-".repeat(50));

    for &k in k_values {
        let omega = bv.dispersion_omega(k);
        let f_ghz = omega / (2.0 * std::f64::consts::PI * 1e9);

        let vg = if k < 1.0 {
            // Numerical derivative at k=0 via one-sided forward difference
            let dk = 1.0_f64; // 1 rad/m
            let o1 = bv.dispersion_omega(dk);
            (o1 - omega) / dk
        } else {
            bv.group_velocity(k)
        };

        let sign_label = if vg < -1.0 {
            "BACKWARD"
        } else if vg > 1.0 {
            "forward"
        } else {
            "~zero"
        };

        println!(
            "  {:>12.3e}  {:>10.4}  {:>12.2}  {:>8}",
            k, f_ghz, vg, sign_label
        );
    }

    // -------------------------------------------------------------------------
    // Section 2: Crossover wavevector k* where v_g changes sign
    // -------------------------------------------------------------------------
    println!("\n--- Section 2: Crossover wavevector k* (v_g = 0) ---\n");

    let k_star = bv.crossover_wavevector();
    let omega_star = bv.dispersion_omega(k_star);
    let f_star = omega_star / (2.0 * std::f64::consts::PI * 1e9);

    println!("  k*             = {:.4e} rad/m", k_star);
    println!("  f(k*)          = {:.4} GHz", f_star);
    println!(
        "  v_g at 0.5·k*  = {:.3} m/s  (backward region)",
        bv.group_velocity(0.5 * k_star)
    );
    println!(
        "  v_g at 2·k*    = {:.3} m/s  (forward region)",
        bv.group_velocity(2.0 * k_star)
    );

    // -------------------------------------------------------------------------
    // Section 3: Quantized thickness modes n = 1 … 5
    // k_perp,n = n π / d;  dispersion evaluated at lateral k = k_perp,n
    // -------------------------------------------------------------------------
    println!("\n--- Section 3: Quantized thickness modes (n = 1..5) ---\n");

    let modes = bv.quantized_thickness_modes(5);
    println!(
        "  {:>6}  {:>14}  {:>10}",
        "mode n", "k_perp (rad/m)", "f (GHz)"
    );
    println!("  {}", "-".repeat(36));

    for (idx, &(k_perp, omega_n)) in modes.iter().enumerate() {
        let n = idx + 1;
        let f_n = omega_n / (2.0 * std::f64::consts::PI * 1e9);
        println!("  {:>6}  {:>14.4e}  {:>10.4}", n, k_perp, f_n);
    }

    println!("\n=============================================================");
    println!("  Done. At small k the group velocity is negative (BACKWARD).");
    println!("  The crossover k* marks the dipolar-to-exchange transition.");
    println!("=============================================================\n");

    Ok(())
}
