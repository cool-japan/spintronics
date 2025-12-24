//! Magnon Propagation and Spin Pumping Detection
//!
//! This example demonstrates the complete magnon propagation process in YIG:
//! 1. RF excitation at one end (FMR condition)
//! 2. Magnon (spin wave) propagation through the chain
//! 3. Detection via spin pumping + ISHE at the other end
//!
//! This reproduces the experimental setup from Saitoh group's research on
//! magnon-driven spin currents and their electrical detection.

use spintronics::magnon::chain::{ChainParameters, SpinChain};
use spintronics::magnon::detector::SpinPumpingDetector;
use spintronics::magnon::solver::{MagnonSolver, RfExcitation};
use spintronics::vector3::Vector3;

fn main() {
    println!("=== Magnon Propagation in YIG Wire ===\n");

    // Physical parameters
    let params = ChainParameters::yig();
    let n_cells = 100; // 100 cells
    let total_length = n_cells as f64 * params.cell_size;

    println!("System Parameters:");
    println!("  Material: YIG (Yttrium Iron Garnet)");
    println!("  Gilbert damping α = {}", params.alpha);
    println!("  Saturation magnetization Ms = {} A/m", params.ms);
    println!("  Exchange stiffness A_ex = {} J/m", params.a_ex);
    println!("  Number of cells: {}", n_cells);
    println!("  Cell size: {} nm", params.cell_size * 1e9);
    println!("  Total length: {} nm\n", total_length * 1e9);

    // Create spin chain with small initial noise
    let mut chain = SpinChain::new_with_noise(n_cells, params.clone(), 0.01);

    // Time step (must be smaller than stability limit)
    let dt_max = params.max_stable_dt();
    let dt = dt_max * 0.5; // Use 50% of maximum for safety
    println!("Time Step:");
    println!("  Maximum stable dt: {} ps", dt_max * 1e12);
    println!("  Using dt: {} ps\n", dt * 1e12);

    // Setup RF excitation at the left end
    let rf_frequency = 10.0e9; // 10 GHz
    let rf_amplitude = 1000.0; // A/m
    let rf_excitation = RfExcitation::new(
        rf_frequency,
        rf_amplitude,
        Vector3::new(0.0, 1.0, 0.0), // Perpendicular to bias field
        (0, 5),                      // Excite first 5 cells
    );

    println!("RF Excitation:");
    println!("  Frequency: {} GHz", rf_frequency * 1e-9);
    println!("  Amplitude: {} A/m", rf_amplitude);
    println!("  Excitation region: cells 0-4\n");

    // Static bias field along x-axis
    let bias_field = Vector3::new(1.0e5, 0.0, 0.0); // 100 kA/m ≈ 0.126 T

    println!("Bias Field:");
    println!(
        "  Magnitude: {} kA/m (≈ {} mT)\n",
        bias_field.x * 1e-3,
        bias_field.x * 1.256
    );

    // Create solver
    let mut solver = MagnonSolver::new(dt, bias_field).with_rf_excitation(rf_excitation);

    // Setup detector at the right end
    let detector_position = n_cells - 1;
    let strip_width = 1.0e-3; // 1 mm
    let mut detector = SpinPumpingDetector::yig_pt(detector_position, strip_width);

    println!("Spin Pumping Detector:");
    println!("  Position: cell {} (right end)", detector_position);
    println!("  Material: Pt strip");
    println!("  Spin-mixing conductance g_r: {} Ω⁻¹m⁻²", detector.g_r);
    println!("  Spin Hall angle θ_SH: {}", detector.theta_sh);
    println!("  Strip width: {} mm\n", strip_width * 1e3);

    // Simulation parameters
    let total_time = 1.0e-9; // 1 nanosecond
    let total_steps = (total_time / dt) as usize;
    let output_interval = total_steps / 50; // 50 output points

    println!("Simulation:");
    println!("  Total time: {} ns", total_time * 1e9);
    println!("  Total steps: {}", total_steps);
    println!("  Output interval: every {} steps\n", output_interval);

    println!("Starting simulation...\n");
    println!("Time (ps)  | Spin Wave Pattern (y-component) | ISHE Voltage (μV)");
    println!("-----------|--------------------------------|------------------");

    let mut max_voltage = 0.0_f64;
    let mut avg_voltage = 0.0;
    let mut sample_count = 0;

    // Main simulation loop
    for step in 0..total_steps {
        // Evolve the system
        solver.step_heun(&mut chain);

        // Detect spin pumping signal
        let voltage = detector.detect(&chain, dt);

        // Track statistics
        if step > total_steps / 2 {
            // Only after transients
            max_voltage = max_voltage.max(voltage.abs());
            avg_voltage += voltage.abs();
            sample_count += 1;
        }

        // Visualization output
        if step % output_interval == 0 {
            let time_ps = solver.time * 1e12;

            // Create ASCII visualization of spin wave
            let wave_visual: String = chain
                .spins
                .iter()
                .step_by(2) // Sample every 2nd cell for display
                .map(|s| {
                    let amplitude = s.y * 50.0; // Amplify for visibility
                    if amplitude > 0.5 {
                        "▲"
                    } else if amplitude < -0.5 {
                        "▼"
                    } else {
                        "·"
                    }
                })
                .collect();

            println!(
                "{:9.2}  | {} | {:16.3e}",
                time_ps,
                wave_visual,
                voltage * 1e6 // Convert to μV
            );
        }
    }

    println!("\n=== Simulation Complete ===\n");

    // Final statistics
    avg_voltage /= sample_count as f64;

    println!("Results:");
    println!("  Maximum ISHE voltage: {:.3} μV", max_voltage * 1e6);
    println!("  Average ISHE voltage: {:.3} μV", avg_voltage * 1e6);

    // Analyze final state
    let final_avg_m = chain.average_magnetization();
    println!("\nFinal State:");
    println!(
        "  Average magnetization: ({:.3}, {:.3}, {:.3})",
        final_avg_m.x, final_avg_m.y, final_avg_m.z
    );

    // Check if magnon reached the detector
    let detector_m = chain.spins[detector_position];
    let detector_amplitude = (detector_m.y.powi(2) + detector_m.z.powi(2)).sqrt();
    println!("  Detector position amplitude: {:.3}", detector_amplitude);

    if detector_amplitude > 0.01 {
        println!("\n✓ Magnon successfully propagated to detector!");
    } else {
        println!("\n✗ Magnon did not reach detector (increase excitation or reduce damping)");
    }

    // Physical interpretation
    println!("\nPhysical Interpretation:");
    println!("  1. RF field at left end excites FMR (ferromagnetic resonance)");
    println!("  2. Magnons (quantized spin waves) propagate through YIG wire");
    println!("  3. At YIG/Pt interface, m × dm/dt generates spin current");
    println!("  4. Spin current converts to voltage via ISHE in Pt strip");
    println!("  5. Measured voltage ∝ magnon amplitude at detector position");

    println!("\nThis demonstrates Saitoh group's spin pumping mechanism:");
    println!("  - Pure spin current generation without charge flow");
    println!("  - Long-distance magnon transport in YIG");
    println!("  - Efficient spin-to-charge conversion in Pt");
}
