//! Spin Ice and Emergent Magnetic Monopoles
//!
//! **Difficulty**: ⭐⭐ Intermediate
//! **Category**: Frustrated Magnetism
//! **Physics**: Pyrochlore lattice, ice rules, Pauling entropy, magnetic monopoles
//!
//! This example demonstrates spin ice physics in Dy2Ti2O7, including:
//! 1. Pyrochlore lattice with Ising spins along local [111] directions
//! 2. Ice rule (2-in, 2-out) constraints on tetrahedra
//! 3. Pauling residual entropy S_P = (R/2) ln(3/2)
//! 4. Emergent magnetic monopole creation energy and Coulomb interactions
//!
//! References:
//! - Bramwell & Gingras, Science 294, 1495 (2001)
//! - Castelnovo, Moessner, Sondhi, Nature 451, 42 (2008)

use spintronics::frustrated::{
    monopole_coulomb_interaction, monopole_creation_energy, monopole_density_factor,
    pauling_entropy_per_spin_kb,
};
use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Spin Ice and Emergent Magnetic Monopoles ===\n");

    // 1. Material parameters for Dy2Ti2O7
    let params = SpinIceParams::dy2ti2o7();
    println!("Material: {}", params.name);
    println!("  J_nn = {:.2} K (exchange coupling)", params.j_nn);
    println!("  D_nn = {:.2} K (dipolar coupling)", params.d_nn);
    println!("  J_eff = J_nn + D_nn = {:.2} K", params.j_eff);
    println!("  Magnetic moment = {:.0} mu_B per ion", params.moment);
    println!(
        "  NN distance = {:.2} Angstrom\n",
        params.nn_distance * 1e10
    );

    // 2. Pauling residual entropy
    println!("=== Pauling Residual Entropy ===");
    let s_pauling = pauling_entropy();
    let s_per_spin = pauling_entropy_per_spin_kb();
    println!("  S_P = (R/2) ln(3/2) = {:.4} J/(mol*K)", s_pauling);
    println!("  s_P = (1/2) ln(3/2) = {:.4} k_B per spin", s_per_spin);
    println!("  (Same as water ice -- hence 'spin ice')\n");

    // 3. Create pyrochlore spin ice system
    println!("=== Pyrochlore Lattice Construction ===");
    let nx = 3;
    let ny = 3;
    let mut ice = SpinIce::new(nx, ny, params.clone()).expect("Failed to create spin ice");
    println!("  Lattice: {}x{} pyrochlore", nx, ny);
    println!("  Number of sites: {}", ice.lattice.num_sites());
    println!("  Number of tetrahedra: {}", ice.tetrahedra.len());

    // 4. Apply ice rules
    println!("\n=== Applying Ice Rules ===");
    let violations_before = ice.count_ice_rule_violations();
    let fraction_before = ice.ice_rule_fraction();
    println!("  Before initialization:");
    println!("    Ice rule violations: {}", violations_before);
    println!("    Fraction satisfied: {:.1}%", fraction_before * 100.0);

    ice.initialize_ice_state(42)
        .expect("Failed to initialize ice state");
    let violations_after = ice.count_ice_rule_violations();
    let fraction_after = ice.ice_rule_fraction();
    println!("  After initialization:");
    println!("    Ice rule violations: {}", violations_after);
    println!("    Fraction satisfied: {:.1}%", fraction_after * 100.0);

    // Check total magnetization
    let mag = ice.total_magnetization();
    println!(
        "    Total magnetization: ({:.3}, {:.3}, {:.3})",
        mag.x, mag.y, mag.z
    );
    println!(
        "    |M| = {:.4} (should be small for ice state)",
        mag.magnitude()
    );

    // 5. Monopole creation energy
    println!("\n=== Magnetic Monopole Energetics ===");
    let e_monopole = monopole_creation_energy(params.j_eff, params.moment, params.nn_distance)
        .expect("Failed to compute monopole energy");

    let e_monopole_kelvin = e_monopole / KB;
    println!("  Monopole creation energy:");
    println!("    E_create = {:.4e} J", e_monopole);
    println!(
        "    E_create = {:.2} K (in temperature units)",
        e_monopole_kelvin
    );
    println!("    (Literature: ~4.35 K for Dy2Ti2O7)");

    // 6. Monopole Coulomb interaction
    println!("\n=== Monopole Coulomb Interaction ===");
    let diamond_const = params.nn_distance * 8.0_f64.sqrt(); // pyrochlore lattice constant
    println!(
        "  Diamond lattice constant a_d = {:.2} Angstrom",
        diamond_const * 1e10
    );
    println!(
        "\n  {:>15} {:>15} {:>15}",
        "Separation (nm)", "V_Coulomb (J)", "V_Coulomb (K)"
    );
    println!("  {}", "-".repeat(50));

    for &r_nm in &[0.5, 1.0, 2.0, 5.0, 10.0, 20.0] {
        let r = r_nm * 1e-9;
        let v_coulomb = monopole_coulomb_interaction(r, params.moment, diamond_const)
            .expect("Failed to compute Coulomb energy");
        println!(
            "  {:>15.1} {:>15.4e} {:>15.3}",
            r_nm,
            v_coulomb,
            v_coulomb / KB
        );
    }

    // 7. Monopole density vs temperature
    println!("\n=== Monopole Density vs Temperature ===");
    println!(
        "  {:>10} {:>20} {:>15}",
        "T (K)", "Boltzmann factor", "Relative density"
    );
    println!("  {}", "-".repeat(50));

    let density_ref =
        monopole_density_factor(10.0, e_monopole).expect("Failed to compute density factor");

    for &temp in &[0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0] {
        let density =
            monopole_density_factor(temp, e_monopole).expect("Failed to compute density factor");
        let relative = if density_ref > 1e-300 {
            density / density_ref
        } else {
            0.0
        };
        println!("  {:>10.1} {:>20.6e} {:>15.6e}", temp, density, relative);
    }

    // 8. Compare Dy2Ti2O7 and Ho2Ti2O7
    println!("\n=== Comparison: Dy2Ti2O7 vs Ho2Ti2O7 ===");
    let ho_params = SpinIceParams::ho2ti2o7();
    let e_ho = monopole_creation_energy(ho_params.j_eff, ho_params.moment, ho_params.nn_distance)
        .expect("Failed to compute Ho monopole energy");

    println!(
        "  {:<12} {:>10} {:>10} {:>10} {:>12}",
        "Material", "J_nn (K)", "D_nn (K)", "J_eff (K)", "E_mono (K)"
    );
    println!("  {}", "-".repeat(60));
    println!(
        "  {:<12} {:>10.2} {:>10.2} {:>10.2} {:>12.2}",
        params.name, params.j_nn, params.d_nn, params.j_eff, e_monopole_kelvin
    );
    println!(
        "  {:<12} {:>10.2} {:>10.2} {:>10.2} {:>12.2}",
        ho_params.name,
        ho_params.j_nn,
        ho_params.d_nn,
        ho_params.j_eff,
        e_ho / KB
    );

    // 9. Frustration parameter
    println!("\n=== Frustration Analysis ===");
    let theta_cw = -1.0; // Curie-Weiss temperature approximation
    let t_order = 0.5; // ordering temperature
    let f_param =
        frustration_parameter(theta_cw, t_order).expect("Failed to compute frustration parameter");
    println!(
        "  Frustration parameter f = |theta_CW| / T_order = {:.1}",
        f_param
    );
    println!("  f > 10: strongly frustrated");
    println!("  Spin ice: f ~ infinity (no long-range ordering)\n");

    println!("=== Summary ===");
    println!("Spin ice physics:");
    println!("  - Pyrochlore lattice with Ising spins along local [111]");
    println!("  - Ice rule (2-in, 2-out) creates macroscopic ground state degeneracy");
    println!(
        "  - Pauling entropy S = {:.4} J/(mol*K) identical to water ice",
        s_pauling
    );
    println!(
        "  - Violations create emergent magnetic monopoles (E ~ {:.1} K)",
        e_monopole_kelvin
    );
    println!("  - Monopoles interact via magnetic Coulomb law V ~ 1/r");

    Ok(())
}
