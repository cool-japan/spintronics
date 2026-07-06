//! Momentum-resolved altermagnet band structure and crystal Hall response.
//!
//! Demonstrates [`AltermagnetBandModel`] for the MnTe (d-wave) and CrSb (g-wave)
//! presets: the non-relativistic spin splitting away from the nodal directions,
//! the compensated (zero) net spin polarization, and the crystal Hall effect —
//! which is zero without spin-orbit coupling (SOC) and becomes nonzero and
//! Néel-sign-reversing once SOC is switched on.
//!
//! Run with: `cargo run --example altermagnet_band_structure`

use spintronics::altermagnet::{AltermagnetBandModel, AltermagneticSymmetry, Spin};
use spintronics::error::Result;
use std::f64::consts::PI;

fn main() -> Result<()> {
    println!("=== Altermagnet momentum-resolved band model ===\n");

    for model in [AltermagnetBandModel::mnte(), AltermagnetBandModel::crsb()] {
        println!(
            "Material: {} altermagnet, a = {:.3e} m",
            model.symmetry, model.a_lattice
        );

        // Evaluate at a fixed radius inside the Brillouin zone (a*|k| = 0.4).
        let k_mag = 0.4 / model.a_lattice;
        println!("  spin splitting |d_up| - |d_down| around the zone (a|k| = 0.4):");
        for step in 0..5 {
            let phi = step as f64 * PI / 8.0;
            let kx = k_mag * phi.cos();
            let ky = k_mag * phi.sin();
            let (up, down) = model.spin_bands(kx, ky);
            let split = model.spin_splitting_at(kx, ky);
            println!(
                "    phi = {:>5.1} deg | E_up = [{:+.4}, {:+.4}]  E_dn = [{:+.4}, {:+.4}]  split = {:+.4} eV",
                phi.to_degrees(),
                up.lower,
                up.upper,
                down.lower,
                down.upper,
                split,
            );
        }

        let k_max = 0.8 / model.a_lattice;
        let polarization = model.net_spin_polarization(41, k_max)?;
        let hall = model.crystal_hall_conductivity(41, k_max)?;
        println!("  net spin polarization  = {polarization:+.3e} (compensated)");
        println!("  crystal Hall (no SOC)  = {hall:+.3e} (vanishes without SOC)\n");
    }

    // Dimensionless d-wave model with SOC switched on: the crystal Hall effect.
    println!("=== Crystal Hall effect with spin-orbit coupling ===");
    let soc_model = AltermagnetBandModel::new(
        AltermagneticSymmetry::DWave,
        1.0, // t_kin
        0.3, // delta_hyb
        0.8, // exchange_m (Neel splitting M)
        0.5, // t_am
        0.4, // lambda_soc (SOC on)
        0.2, // fermi_energy
        0.0, // neel_angle
        1.0, // a_lattice (dimensionless demo units)
    )?;
    let hall = soc_model.crystal_hall_conductivity(61, 2.5)?;

    let mut reversed = soc_model.clone();
    reversed.exchange_m = -soc_model.exchange_m;
    let hall_reversed = reversed.crystal_hall_conductivity(61, 2.5)?;

    println!("  sigma_xy (Neel +M) = {hall:+.6e}");
    println!("  sigma_xy (Neel -M) = {hall_reversed:+.6e}");
    println!("  -> the crystal Hall conductivity reverses with the Neel vector.");

    // Berry curvature of the two occupied spin channels at a representative k.
    let (kx, ky) = (0.7, 0.3);
    let berry_up =
        soc_model.berry_curvature(kx, ky, Spin::Up, spintronics::altermagnet::Band::Lower);
    let berry_dn =
        soc_model.berry_curvature(kx, ky, Spin::Down, spintronics::altermagnet::Band::Lower);
    println!(
        "  Berry curvature at k = ({kx}, {ky}) lower band: up = {berry_up:+.4e}, down = {berry_dn:+.4e}"
    );

    Ok(())
}
