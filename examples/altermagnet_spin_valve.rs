//! Altermagnet Spin Valve: Giant Magnetoresistance without Ferromagnetism
//!
//! **Difficulty**: ⭐⭐⭐ Advanced
//! **Category**: Altermagnetic spintronics
//! **Physics**: GMR-analog transport, zero-net-moment devices
//!
//! Conventional GMR spin valves need two ferromagnetic layers. This example
//! demonstrates an altermagnetic analog: two altermagnetic layers (CrSb or
//! MnTe), each with *zero net magnetization*, separated by a Cu spacer.
//! Rotating the relative Neel/crystal angle `theta` between the layers still
//! produces a strong, GMR-like magnetoresistance -- with no ferromagnetic
//! moment anywhere in the device.
//!
//! Reference: Smejkal et al., Phys. Rev. X 12, 040501 (2022) and 12, 031042
//! (2022); conventional GMR: Baibich et al., Phys. Rev. Lett. 61, 2472 (1988)

use std::f64::consts::PI;

use spintronics::altermagnet::AltermagnetSpinValve;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Altermagnet Spin Valve: GMR without Ferromagnetism ===\n");

    let crsb_valve = AltermagnetSpinValve::crsb_cu_crsb();
    let mnte_valve = AltermagnetSpinValve::mnte_based();

    println!(
        "Device 1: CrSb / Cu / CrSb (g-wave, n = {})",
        crsb_valve.harmonic_order()
    );
    println!(
        "Device 2: MnTe / Cu / MnTe (d-wave, n = {})\n",
        mnte_valve.harmonic_order()
    );

    // 1. Zero net magnetization -- the whole point of this device.
    println!("=== Zero Net Magnetization (both devices, always) ===");
    println!(
        "  CrSb valve net magnetization: {:.1} A/m",
        crsb_valve.net_magnetization()
    );
    println!(
        "  MnTe valve net magnetization: {:.1} A/m",
        mnte_valve.net_magnetization()
    );
    println!("  (No ferromagnetic layer is ever present in this device.)\n");

    // 2. GMR ratio vs relative angle.
    for (label, valve) in [("CrSb/Cu/CrSb", &crsb_valve), ("MnTe/Cu/MnTe", &mnte_valve)] {
        let n = valve.harmonic_order();
        let theta_swap = PI / f64::from(n);

        println!("=== {label}: GMR ratio vs relative angle (n = {n}) ===");
        println!("{:>14} {:>16}", "theta (deg)", "GMR ratio");
        println!("{}", "-".repeat(32));

        let n_points = 13;
        for i in 0..n_points {
            let theta = theta_swap * f64::from(i) / f64::from(n_points - 1);
            println!(
                "{:14.2} {:16.5}",
                theta.to_degrees(),
                valve.gmr_ratio(theta)
            );
        }
        println!(
            "  Channel-swap angle theta = pi/{n} = {:.2} deg -> maximal GMR = {:.4}\n",
            theta_swap.to_degrees(),
            valve.gmr_max
        );
    }

    // 3. Resistance and conductance at aligned vs. channel-swapped angles.
    println!("=== Resistance / Conductance (CrSb/Cu/CrSb, area = 100nm x 100nm) ===");
    let area = 1.0e-14; // m^2
    let n = crsb_valve.harmonic_order();
    let theta_swap = PI / f64::from(n);

    let r_aligned = crsb_valve.resistance(0.0, area)?;
    let r_swapped = crsb_valve.resistance(theta_swap, area)?;
    let g_aligned = crsb_valve.conductance(0.0, area)?;
    let g_swapped = crsb_valve.conductance(theta_swap, area)?;

    println!("  Aligned   (theta=0):      R = {r_aligned:.4e} Ohm, G = {g_aligned:.4e} S");
    println!("  Channel-swapped (theta=pi/n): R = {r_swapped:.4e} Ohm, G = {g_swapped:.4e} S");
    println!(
        "  Resistance ratio R_swap/R_aligned = {:.4} (matches 1 + gmr_max = {:.4})",
        r_swapped / r_aligned,
        1.0 + crsb_valve.gmr_max
    );
    println!(
        "  Sanity check R*G: aligned = {:.6}, swapped = {:.6} (should be 1.0)\n",
        r_aligned * g_aligned,
        r_swapped * g_swapped
    );

    // 4. Summary.
    println!("=== Summary ===");
    println!("This device reproduces the functional form of conventional GMR,");
    println!("(1 - cos(n*theta))/2, using the harmonic order n of the shared");
    println!("altermagnetic symmetry (d/g/i-wave -> n=2/4/6) in place of the");
    println!("usual ferromagnetic n=1 case -- while both layers remain");
    println!("exactly zero-net-moment altermagnets at every relative angle.");

    Ok(())
}
