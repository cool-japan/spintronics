//! LLG Dynamics Visualization Example
//!
//! This example demonstrates how to export time-dependent magnetization dynamics:
//! - VTK: Time series for animation in ParaView
//! - CSV: Time evolution for plotting
//! - JSON: Complete simulation data with parameters
//!
//! Run with: cargo run --example llg_dynamics_visualization

use spintronics::prelude::*;
use spintronics::visualization::csv::write_time_series;
use spintronics::visualization::json::{JsonWriter, SimulationData};
use spintronics::visualization::vtk::VtkWriter;

fn main() -> std::io::Result<()> {
    println!("=== LLG Dynamics Visualization Example ===\n");

    // Simulation parameters
    let yig = Ferromagnet::yig();
    let dt = 1e-12; // 1 ps
    let total_time = 5e-9; // 5 ns
    let num_steps = (total_time / dt) as usize;
    let save_interval = 1000; // Save every 1000 steps

    println!("Material: YIG");
    println!("Gilbert damping: {}", yig.alpha);
    println!("Time step: {} s", dt);
    println!("Total time: {} s", total_time);
    println!("Save interval: {} steps\n", save_interval);

    // Initial state: tilted magnetization
    let mut m = Vector3::new(0.1, 0.0, 1.0).normalize();

    // External field (perpendicular)
    let h_ext = Vector3::new(0.0, 0.0, 1.0);
    let h_magnitude = 1000.0 / (GAMMA / (2.0 * std::f64::consts::PI)); // ~ 1000 Oe

    println!(
        "Initial magnetization: ({:.3}, {:.3}, {:.3})",
        m.x, m.y, m.z
    );
    println!(
        "External field: {:.1} Oe (z-direction)\n",
        h_magnitude * GAMMA / (2.0 * std::f64::consts::PI)
    );

    // Setup solver
    let solver = LlgSolver::new(yig.alpha, dt);

    // Data storage
    let mut times = Vec::new();
    let mut mx_data = Vec::new();
    let mut my_data = Vec::new();
    let mut mz_data = Vec::new();

    // VTK writer for time series (using simple 1D grid)
    let mut vtk_writer = VtkWriter::new("llg_dynamics");

    // JSON simulation data
    let mut sim_data = SimulationData::new("llg_precession", "0.1.0");
    sim_data.add_metadata("material", "YIG");
    sim_data.add_metadata("solver", "RK4");
    sim_data.add_parameter("alpha", yig.alpha);
    sim_data.add_parameter("dt", dt);
    sim_data.add_parameter("h_field", h_magnitude);
    sim_data.add_parameter("gamma", GAMMA);

    // Simulation loop
    println!("Running simulation...");

    for step in 0..num_steps {
        let t = step as f64 * dt;

        // Evolve with RK4 method
        m = solver.step_rk4(m, |_| h_ext * h_magnitude);

        // Save data
        if step % save_interval == 0 {
            times.push(t);
            mx_data.push(m.x);
            my_data.push(m.y);
            mz_data.push(m.z);

            // VTK snapshot (simple visualization as point cloud)
            let snapshot_spins = vec![m];
            vtk_writer.write_snapshot(&snapshot_spins, (1, 1, 1))?;

            // JSON time point
            sim_data.add_time_point(t, vec![m.x, m.y, m.z]);

            if step % (save_interval * 10) == 0 {
                println!(
                    "  t = {:.2} ns: m = ({:.3}, {:.3}, {:.3})",
                    t * 1e9,
                    m.x,
                    m.y,
                    m.z
                );
            }
        }
    }

    println!("\n=== Exporting Data ===\n");

    // === 1. CSV Export ===
    println!("1. Exporting CSV files...");

    // Time series data
    write_time_series(
        "llg_dynamics.csv",
        &times,
        &[&mx_data, &my_data, &mz_data],
        &["mx", "my", "mz"],
    )?;
    println!("   ✓ Saved: llg_dynamics.csv\n");

    // === 2. JSON Export ===
    println!("2. Exporting JSON...");
    JsonWriter::write("llg_dynamics.json", &sim_data)?;
    println!("   ✓ Saved: llg_dynamics.json\n");

    // === 3. VTK Time Series ===
    println!("3. VTK time series already saved");
    println!("   ✓ Saved: llg_dynamics_*.vtu ({} files)\n", times.len());

    // === Analysis ===
    println!("=== Analysis ===");

    // Calculate precession frequency
    if times.len() > 2 {
        // Find peaks in mx
        let mut peaks = Vec::new();
        for i in 1..(mx_data.len() - 1) {
            if mx_data[i] > mx_data[i - 1] && mx_data[i] > mx_data[i + 1] {
                peaks.push(times[i]);
            }
        }

        if peaks.len() >= 2 {
            let period = (peaks[peaks.len() - 1] - peaks[0]) / (peaks.len() - 1) as f64;
            let frequency = 1.0 / period;
            println!("\nPrecession frequency: {:.3} GHz", frequency * 1e-9);
            println!(
                "Expected (Larmor): {:.3} GHz",
                GAMMA * h_magnitude / (2.0 * std::f64::consts::PI) * 1e-9
            );
        }
    }

    // Final magnetization
    println!(
        "\nFinal magnetization: ({:.6}, {:.6}, {:.6})",
        m.x, m.y, m.z
    );
    println!("Magnitude: {:.6} (should be 1.0)", m.magnitude());

    println!("\n=== Visualization Instructions ===");
    println!("\nParaView (Animation):");
    println!("  1. Open llg_dynamics_*.vtu");
    println!("  2. Apply 'Glyph' filter");
    println!("  3. Press Play to see precession\n");

    println!("Python (Plotting):");
    println!("  import pandas as pd");
    println!("  import matplotlib.pyplot as plt");
    println!("  df = pd.read_csv('llg_dynamics.csv')");
    println!("  plt.plot(df['time']*1e9, df['mx'], label='mx')");
    println!("  plt.plot(df['time']*1e9, df['my'], label='my')");
    println!("  plt.plot(df['time']*1e9, df['mz'], label='mz')");
    println!("  plt.xlabel('Time (ns)')");
    println!("  plt.legend()");
    println!("  plt.show()");

    Ok(())
}
