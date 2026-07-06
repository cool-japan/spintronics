//! Magnon Bose-Einstein Condensation
//!
//! **Difficulty**: ⭐⭐⭐ Advanced
//! **Category**: Magnon Physics
//! **Physics**: Magnon BEC, parametric pumping, Suhl threshold, Gross-Pitaevskii
//!
//! This example demonstrates magnon Bose-Einstein condensation (BEC),
//! where bosonic magnon quasiparticles accumulate in the lowest energy state.
//! We show:
//!
//! 1. BEC transition temperature and critical density
//! 2. Condensate fraction vs temperature
//! 3. Healing length and Bogoliubov sound speed
//! 4. Parametric pumping and Suhl instability threshold
//! 5. Magnon growth rate above threshold
//!
//! References:
//! - Demokritov et al., Nature 443, 430 (2006)
//! - Serga et al., Nat. Commun. 5, 3452 (2014)

// This example exercises `spintronics::magnon::bec` (magnon BEC, parametric
// pumping), which is excluded from wasm32 builds (see
// `#[cfg(not(target_arch = "wasm32"))]` on `pub mod magnon;` in
// `src/lib.rs`), so it is a no-op there.
#[cfg(not(target_arch = "wasm32"))]
use spintronics::magnon::bec::{
    bec_temperature, critical_density, magnon_distribution, MagnonCondensate, ParametricPumping,
};
#[cfg(not(target_arch = "wasm32"))]
use spintronics::prelude::*;

#[cfg(not(target_arch = "wasm32"))]
fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Magnon Bose-Einstein Condensation ===\n");

    // Physical parameters for YIG magnons
    let effective_mass = 1.0e-28; // ~0.1 m_e (magnon effective mass in YIG)
    let magnon_density = 1.0e23; // 10^23 /m^3 (high pumping regime)
    let interaction_strength = 1.0e-50; // magnon-magnon interaction [J*m^3]

    // 1. BEC transition temperature
    println!("=== BEC Transition Temperature ===");
    let t_bec =
        bec_temperature(magnon_density, effective_mass).expect("Failed to compute BEC temperature");
    println!("  Magnon density: {:.2e} /m^3", magnon_density);
    println!(
        "  Effective mass: {:.2e} kg ({:.2} m_e)",
        effective_mass,
        effective_mass / ME
    );
    println!("  T_BEC = {:.2} K", t_bec);

    // 2. Critical density at various temperatures
    println!("\n=== Critical Density vs Temperature ===");
    println!(
        "  {:>8} {:>20} {:>15}",
        "T (K)", "n_critical (/m^3)", "n/n_c"
    );
    println!("  {}", "-".repeat(48));
    for &temp in &[50.0, 100.0, 200.0, 300.0, 500.0] {
        let n_c =
            critical_density(effective_mass, temp).expect("Failed to compute critical density");
        let ratio = magnon_density / n_c;
        println!("  {:>8.0} {:>20.4e} {:>15.4e}", temp, n_c, ratio);
    }

    // 3. Condensate fraction vs temperature
    println!("\n=== Condensate Fraction vs Temperature ===");
    println!(
        "  {:>8} {:>12} {:>15} {:>15} {:>12}",
        "T (K)", "T/T_BEC", "Cond. fraction", "n_0 (/m^3)", "n_th (/m^3)"
    );
    println!("  {}", "-".repeat(70));

    let temp_fracs = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 1.0, 1.5, 2.0];
    for &frac in &temp_fracs {
        let temp = frac * t_bec;
        // Avoid exactly 0 for the distribution but allow it for condensate
        let condensate =
            MagnonCondensate::new(temp, magnon_density, effective_mass, interaction_strength)
                .expect("Failed to create condensate");

        println!(
            "  {:>8.2} {:>12.3} {:>15.6} {:>15.4e} {:>12.4e}",
            temp,
            frac,
            condensate.condensate_fraction,
            condensate.condensate_density(),
            condensate.thermal_density()
        );
    }

    // 4. Condensate quantum properties
    println!("\n=== Condensate Quantum Properties ===");
    let temp_low = 0.3 * t_bec;
    let condensate = MagnonCondensate::new(
        temp_low,
        magnon_density,
        effective_mass,
        interaction_strength,
    )
    .expect("Failed to create low-T condensate");

    let t_bec_check = condensate
        .transition_temperature()
        .expect("Failed to get T_BEC");
    println!(
        "  Temperature: {:.2} K ({:.1}% of T_BEC)",
        temp_low,
        100.0 * temp_low / t_bec_check
    );
    println!(
        "  Condensate fraction: {:.4}",
        condensate.condensate_fraction
    );
    println!(
        "  Condensate density: {:.4e} /m^3",
        condensate.condensate_density()
    );
    println!(
        "  Chemical potential: {:.4e} J",
        condensate.chemical_potential
    );

    let xi = condensate
        .healing_length()
        .expect("Failed to compute healing length");
    let c_s = condensate
        .sound_speed()
        .expect("Failed to compute sound speed");
    println!("  Healing length: {:.2} nm", xi * 1e9);
    println!("  Bogoliubov sound speed: {:.2} m/s", c_s);

    // 5. Bose-Einstein distribution
    println!("\n=== Bose-Einstein Distribution ===");
    let temp_dist = 300.0;
    let mu = -1.0e-22; // chemical potential below band minimum
    println!("  T = {:.0} K, mu = {:.2e} J", temp_dist, mu);
    println!("  {:>15} {:>20}", "omega (rad/s)", "n(omega)");
    println!("  {}", "-".repeat(40));
    for &omega_ghz in &[1.0, 5.0, 10.0, 50.0, 100.0] {
        let omega = omega_ghz * 1e9 * 2.0 * std::f64::consts::PI;
        let n = magnon_distribution(omega, temp_dist, mu).expect("Failed to compute distribution");
        println!("  {:>15.2e} {:>20.6e}", omega, n);
    }

    // 6. Parametric pumping
    println!("\n=== Parametric Pumping ===");
    let pump_freq = 2.0 * std::f64::consts::PI * 14.0e9; // 14 GHz pump
    let damping = 1.0e-4; // YIG low damping
    let coupling = 1.0e-23; // parametric coupling strength

    let pumping = ParametricPumping::new(pump_freq, damping, coupling)
        .expect("Failed to create pumping parameters");

    println!(
        "  Pump frequency: {:.2} GHz",
        pump_freq / (2.0 * std::f64::consts::PI * 1e9)
    );
    println!(
        "  Magnon frequency: {:.2} GHz (omega_p/2)",
        pumping.magnon_frequency() / (2.0 * std::f64::consts::PI * 1e9)
    );
    println!("  Gilbert damping: {:.1e}", damping);
    println!("  Suhl threshold: {:.4e}", pumping.suhl_threshold());

    // 7. Growth rate vs pump field
    println!("\n=== Magnon Growth Rate vs Pump Field ===");
    let h_th = pumping.suhl_threshold();
    println!("  {:>12} {:>15} {:>10}", "h/h_th", "Growth rate", "Status");
    println!("  {}", "-".repeat(42));
    for &ratio in &[0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0] {
        let h = ratio * h_th;
        let gamma_rate = pumping.growth_rate(h);
        let status = if gamma_rate > 0.0 {
            "GROWING"
        } else {
            "DAMPED"
        };
        println!("  {:>12.2} {:>15.4e} {:>10}", ratio, gamma_rate, status);
    }

    println!("\n=== Summary ===");
    println!("Magnon BEC physics:");
    println!(
        "  - BEC transition at T_BEC = {:.2} K for n = {:.0e} /m^3",
        t_bec, magnon_density
    );
    println!(
        "  - Condensate fraction = {:.1}% at T = {:.1} K",
        condensate.condensate_fraction * 100.0,
        temp_low
    );
    println!("  - Healing length xi = {:.1} nm", xi * 1e9);
    println!("  - Bogoliubov sound speed c_s = {:.1} m/s", c_s);
    println!("  - Parametric pumping creates magnon pairs at omega_p/2");

    Ok(())
}

#[cfg(target_arch = "wasm32")]
fn main() {}
