//! Skyrmion Visualization Example
//!
//! This example demonstrates how to export skyrmion spin textures in multiple formats:
//! - VTK: For 3D visualization in ParaView
//! - CSV: For plotting cross-sections and profiles
//! - JSON: For structured data with metadata
//!
//! Run with: cargo run --example skyrmion_visualization

use std::f64::consts::PI;

use spintronics::prelude::*;
use spintronics::visualization::csv::CsvWriter;
use spintronics::visualization::json::{JsonWriter, SimulationData};
use spintronics::visualization::vtk::VtkWriter;

fn main() -> std::io::Result<()> {
    println!("=== Skyrmion Visualization Example ===\n");

    // Setup grid
    let nx = 64;
    let ny = 64;
    let nz = 1;
    let num_points = nx * ny * nz;

    println!("Grid size: {}x{}x{}", nx, ny, nz);
    println!("Total points: {}\n", num_points);

    // Create skyrmion texture
    let center_x = (nx / 2) as f64;
    let center_y = (ny / 2) as f64;
    let radius = 10.0;
    let chirality = 1;

    let mut spins = Vec::with_capacity(num_points);

    for j in 0..ny {
        for i in 0..nx {
            let x = i as f64 - center_x;
            let y = j as f64 - center_y;
            let r = (x * x + y * y).sqrt();

            let spin = if r < 0.1 {
                // Core: down
                Vector3::new(0.0, 0.0, -1.0)
            } else if r < radius * 2.0 {
                // Skyrmion profile
                let theta = PI * (1.0 - (-((r - radius) / 3.0).powi(2)).exp());
                let phi = (y.atan2(x) + chirality as f64 * PI / 2.0) % (2.0 * PI);

                Vector3::new(
                    theta.sin() * phi.cos(),
                    theta.sin() * phi.sin(),
                    theta.cos(),
                )
            } else {
                // Background: up
                Vector3::new(0.0, 0.0, 1.0)
            };

            spins.push(spin.normalize());
        }
    }

    println!("Skyrmion parameters:");
    println!("  Center: ({}, {})", center_x, center_y);
    println!("  Radius: {}", radius);
    println!("  Chirality: {}\n", chirality);

    // === 1. VTK Export for ParaView ===
    println!("1. Exporting to VTK format...");
    let mut vtk_writer = VtkWriter::new("skyrmion_3d");

    // Calculate topological charge density
    let mut topo_charge = vec![0.0; num_points];
    for j in 1..(ny - 1) {
        for i in 1..(nx - 1) {
            let idx = i + nx * j;
            let idx_right = (i + 1) + nx * j;
            let idx_up = i + nx * (j + 1);

            let m = &spins[idx];
            let m_dx = &spins[idx_right];
            let m_dy = &spins[idx_up];

            // Approximate topological charge density
            let dm_dx = (*m_dx - *m) * 1.0;
            let dm_dy = (*m_dy - *m) * 1.0;
            topo_charge[idx] = m.dot(&dm_dx.cross(&dm_dy)) / (4.0 * PI);
        }
    }

    vtk_writer.write_snapshot_with_scalar(
        &spins,
        &topo_charge,
        "TopologicalCharge",
        (nx, ny, nz),
    )?;
    println!("   ✓ Saved: skyrmion_3d_00000.vtu");
    println!("   Open in ParaView and apply 'Glyph' filter\n");

    // === 2. CSV Export for Plotting ===
    println!("2. Exporting to CSV format...");

    // Export cross-section along x-axis (y = center)
    let mut csv_writer = CsvWriter::new("skyrmion_cross_section.csv")?;
    csv_writer.write_header(&["x", "mx", "my", "mz", "magnitude"])?;

    let y_center = ny / 2;
    for i in 0..nx {
        let idx = i + nx * y_center;
        let spin = &spins[idx];
        csv_writer.write_row(&[i as f64, spin.x, spin.y, spin.z, spin.magnitude()])?;
    }
    println!("   ✓ Saved: skyrmion_cross_section.csv");

    // Export radial profile
    let mut radial_csv = CsvWriter::new("skyrmion_radial_profile.csv")?;
    radial_csv.write_header(&["radius", "mz", "in_plane_magnitude"])?;

    for i in 0..nx {
        let x = i as f64 - center_x;
        let idx = i + nx * y_center;
        let spin = &spins[idx];
        let in_plane = (spin.x * spin.x + spin.y * spin.y).sqrt();
        radial_csv.write_row(&[x.abs(), spin.z, in_plane])?;
    }
    println!("   ✓ Saved: skyrmion_radial_profile.csv");
    println!("   Plot with matplotlib/gnuplot\n");

    // === 3. JSON Export with Metadata ===
    println!("3. Exporting to JSON format...");
    let mut sim_data = SimulationData::new("skyrmion_texture", "0.1.0");

    // Add metadata
    sim_data.add_metadata("material", "Pt/CoFeB");
    sim_data.add_metadata("texture_type", "Neel_skyrmion");
    sim_data.add_metadata("chirality", &chirality.to_string());
    sim_data.add_metadata("description", "DMI-stabilized Neel-type skyrmion");

    // Add parameters
    sim_data.add_parameter("grid_nx", nx as f64);
    sim_data.add_parameter("grid_ny", ny as f64);
    sim_data.add_parameter("radius", radius);
    sim_data.add_parameter("center_x", center_x);
    sim_data.add_parameter("center_y", center_y);

    // Add vector snapshot
    sim_data.add_vector_snapshot(0.0, spins.clone());

    // Add scalar snapshot
    sim_data.add_scalar_snapshot(0.0, "topological_charge_density", topo_charge);

    // Calculate and add some analysis data
    let total_mz: f64 = spins.iter().map(|s| s.z).sum();
    let avg_mz = total_mz / num_points as f64;
    sim_data.add_parameter("average_mz", avg_mz);

    JsonWriter::write("skyrmion_data.json", &sim_data)?;
    println!("   ✓ Saved: skyrmion_data.json");
    println!("   Load with Python: json.load()\n");

    // === Summary ===
    println!("=== Export Complete ===");
    println!("\nVisualization Instructions:");
    println!("  ParaView: Open skyrmion_3d_00000.vtu");
    println!("            Apply Filters → Glyph (arrows)");
    println!("            Color by 'Sz' or 'TopologicalCharge'");
    println!("\n  Python:");
    println!("    import pandas as pd");
    println!("    df = pd.read_csv('skyrmion_cross_section.csv')");
    println!("    plt.plot(df['x'], df['mz'])");
    println!("\n  Analysis:");
    println!("    import json");
    println!("    with open('skyrmion_data.json') as f:");
    println!("        data = json.load(f)");
    println!("    print(data['parameters'])");

    Ok(())
}
