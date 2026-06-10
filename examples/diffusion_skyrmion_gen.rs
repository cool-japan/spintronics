//! Diffusion Model for Skyrmion Lattice Generation
//!
//! Trains a DDPM to learn the distribution of skyrmion spin textures,
//! then generates new textures and measures their topological charge.
//!
//! Run with: cargo run --example diffusion_skyrmion_gen --features autodiff

use spintronics::prelude::*;
use spintronics::texture::skyrmion::{Chirality, Helicity};

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Diffusion Model for Skyrmion Texture Generation ===\n");

    let nx = 16_usize;
    let ny = 16_usize;
    let domain_size = 200e-9_f64; // 200 nm domain
    let wall_width = 10e-9_f64;

    println!(
        "Building training dataset ({nx}×{ny} grid, {:.0}nm domain)...",
        domain_size * 1e9
    );

    // SkyrmionLattice::square(nx, ny, lattice_constant, skyrmion_radius, helicity, chirality)
    // Places skyrmions at grid positions (i*a, j*a).
    let spacing = 80e-9_f64; // lattice constant [m]
    let radius = 25e-9_f64; // skyrmion radius [m]

    // Build training textures from 2×2 skyrmion arrays with varying lattice constants
    let mut textures = Vec::new();
    for variant in 0..10_usize {
        // Slight lattice constant variation to increase dataset diversity
        let a = spacing * (1.0 + 0.05 * (variant as f64 - 5.0) / 5.0);
        let lattice =
            SkyrmionLattice::square(2, 2, a, radius, Helicity::Neel, Chirality::Clockwise);
        let texture =
            SpinTexture::from_skyrmion_lattice(&lattice, nx, ny, domain_size, wall_width)?;
        textures.push(texture);
    }

    println!(
        "Training dataset: {} textures of size {}×{}×3",
        textures.len(),
        nx,
        ny
    );

    // Build DDPM model: linear β-schedule, 100 steps, 2-layer MLP, hidden_dim=32
    let schedule = NoiseSchedule::linear(100, 1e-4, 0.02)?;
    let mut model = DiffusionModel::new(nx, ny, 32, 2, schedule)?;
    println!(
        "Model: 2-layer tanh MLP, hidden_dim=32, params={}",
        model.params.len()
    );

    println!("\nTraining for 30 epochs (Adam, lr=1e-3)...");
    let losses = model.train(&textures, 30, 1e-3)?;
    let first_loss = losses.first().copied().unwrap_or(0.0);
    let last_loss = losses.last().copied().unwrap_or(0.0);
    println!("  Initial loss: {:.4}", first_loss);
    println!("  Final loss:   {:.4}", last_loss);
    println!(
        "  Loss reduced: {}  (ratio = {:.2}×)",
        last_loss < first_loss,
        if last_loss > 0.0 {
            first_loss / last_loss
        } else {
            f64::INFINITY
        }
    );

    // Generate samples and measure topological charge
    println!("\nGenerating 3 spin texture samples via DDPM reverse diffusion...");
    let dx = domain_size / nx as f64;
    let dy = domain_size / ny as f64;

    for seed in [1_u64, 2, 3] {
        let sample = model.sample(nx, ny, seed)?;
        let q = model.topological_charge(&sample, dx, dy);
        println!("  Sample {seed}: topological charge Q = {q:.2}");
    }

    // Reference: compute topological charge of a known 2×2 skyrmion lattice
    // (each skyrmion contributes Q=±1, so 4 Clockwise → Q_total ≈ +4)
    let ref_lattice =
        SkyrmionLattice::square(2, 2, spacing, radius, Helicity::Neel, Chirality::Clockwise);
    let ref_tex =
        SpinTexture::from_skyrmion_lattice(&ref_lattice, nx, ny, domain_size, wall_width)?;
    let ref_q = model.topological_charge(&ref_tex, dx, dy);
    println!("\nReference 2×2 Clockwise Néel lattice: Q = {ref_q:.2} (should be ≈ +4)");
    println!(
        "  Integer topological charge (discrete): {}",
        ref_q.round() as i64
    );

    println!("\n=== Summary ===");
    println!("  Schedule: linear β ∈ [1e-4, 0.02], T=100 steps");
    println!(
        "  Architecture: MLP({} → 32 → {})",
        3 * nx * ny + 1,
        3 * nx * ny
    );
    println!("  Training: 30 epochs, Adam");
    println!(
        "  Loss improvement: {:.1}×",
        if last_loss > 0.0 {
            first_loss / last_loss
        } else {
            f64::INFINITY
        }
    );

    Ok(())
}
