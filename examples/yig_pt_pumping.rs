//! YIG/Pt Spin Pumping Simulation
//!
//! **Difficulty**: ⭐ Beginner
//! **Category**: Core Physics - Basic Spintronics
//! **Physics**: Spin pumping, Inverse Spin Hall Effect (ISHE)
//!
//! This example demonstrates the complete spin pumping + ISHE process:
//! 1. Magnetization dynamics in YIG under FMR conditions
//! 2. Spin current generation via spin pumping
//! 3. Electric field generation via inverse spin Hall effect in Pt
//!
//! This simulation reproduces the experimental setup from:
//! E. Saitoh et al., Appl. Phys. Lett. 88, 182509 (2006)

use spintronics::prelude::*;

fn main() {
    println!("=== YIG/Pt Spin Pumping + ISHE Simulation ===\n");

    // 1. Setup materials (YIG/Pt system)
    let yig = Ferromagnet::yig();
    println!("Ferromagnet: YIG");
    println!("  Gilbert damping α = {}", yig.alpha);
    println!("  Saturation magnetization M_s = {} A/m\n", yig.ms);

    let interface = SpinInterface::yig_pt();
    println!("Interface: YIG/Pt");
    println!("  Spin-mixing conductance g_r = {} Ω⁻¹m⁻²", interface.g_r);
    println!(
        "  Interface normal = ({:.1}, {:.1}, {:.1})\n",
        interface.normal.x, interface.normal.y, interface.normal.z
    );

    let pt_strip = InverseSpinHall::platinum();
    println!("Normal metal: Pt");
    println!("  Spin Hall angle θ_SH = {}", pt_strip.theta_sh);
    println!("  Resistivity ρ = {} Ω·m\n", pt_strip.rho);

    // 2. Initial state
    let mut m = Vector3::new(1.0, 0.0, 0.0); // Magnetization along x
    println!(
        "Initial magnetization: m = ({:.3}, {:.3}, {:.3})\n",
        m.x, m.y, m.z
    );

    // 3. External field (FMR condition: field perpendicular to initial m)
    let h_ext = Vector3::new(0.0, 0.0, 1.0); // Field along z
    let h_magnitude = 0.1; // Tesla
    let h_eff = h_ext * h_magnitude;
    println!(
        "External field: H = {:.3} T along ({:.1}, {:.1}, {:.1})\n",
        h_magnitude, h_ext.x, h_ext.y, h_ext.z
    );

    // 4. Time evolution parameters
    let dt = 1.0e-13; // 0.1 ps
    let total_time = 1.0e-9; // 1 ns
    let steps = (total_time / dt) as usize;
    let output_interval = steps / 20; // Output 20 snapshots

    println!("Time evolution:");
    println!("  Time step dt = {} ps", dt * 1e12);
    println!("  Total time = {} ns", total_time * 1e9);
    println!("  Total steps = {}\n", steps);

    println!("Time (ps)    |m_x|    |m_y|    |m_z|    |J_s| (J/m²)    E_field (V/m)");
    println!("------------------------------------------------------------------------");

    // 5. Main simulation loop
    for step in 0..steps {
        // Calculate LLG dynamics
        let dm_dt = calc_dm_dt(m, h_eff, GAMMA, yig.alpha);

        // Update magnetization (Euler method)
        m = (m + dm_dt * dt).normalize();

        // Calculate spin pumping current
        let js_vec = spin_pumping_current(&interface, m, dm_dt);

        // Calculate ISHE voltage
        // Spin current flows along interface normal (y direction)
        // Spin polarization is along js_vec direction
        let e_field = pt_strip.convert(interface.normal, js_vec);

        // Output at intervals
        if step % output_interval == 0 {
            let time_ps = (step as f64) * dt * 1e12;
            println!(
                "{:8.2}    {:6.3}  {:6.3}  {:6.3}  {:12.5e}  {:12.5e}",
                time_ps,
                m.x.abs(),
                m.y.abs(),
                m.z.abs(),
                js_vec.magnitude(),
                e_field.magnitude()
            );
        }
    }

    println!("\n=== Simulation Complete ===");

    // Final state analysis
    let dm_dt_final = calc_dm_dt(m, h_eff, GAMMA, yig.alpha);
    let js_final = spin_pumping_current(&interface, m, dm_dt_final);
    let e_final = pt_strip.convert(interface.normal, js_final);

    println!("\nFinal state:");
    println!("  Magnetization: m = ({:.3}, {:.3}, {:.3})", m.x, m.y, m.z);
    println!("  |m| = {:.6} (should be 1.0)", m.magnitude());
    println!(
        "  Spin current density: |J_s| = {:.5e} J/m²",
        js_final.magnitude()
    );
    println!("  Electric field: |E| = {:.5e} V/m", e_final.magnitude());

    // Estimate measurable voltage for a typical strip
    let strip_width = 1.0e-3; // 1 mm
    let voltage = e_final.magnitude() * strip_width;
    println!("\nEstimated voltage across {} mm strip:", strip_width * 1e3);
    println!("  V_ISHE = {:.3e} V = {:.3} μV", voltage, voltage * 1e6);

    println!("\nNote: This simulation uses simplified physics.");
    println!("For quantitative predictions, include:");
    println!("  - Anisotropy fields");
    println!("  - Demagnetization fields");
    println!("  - Spin diffusion in Pt layer");
    println!("  - Realistic FMR excitation conditions");
}
