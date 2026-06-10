//! Hybrid Quantum-Classical NN for Magnetic Hamiltonians
//!
//! Trains a neural network to reproduce the magnon spectrum of a 1D Heisenberg
//! ferromagnet via Bogoliubov transformation, with gradients from central
//! finite-difference (classical autodiff surrogate).
//!
//! Run with: cargo run --example variational_magnon_nn --features autodiff

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Hybrid Quantum-Classical NN for Magnon Hamiltonians ===\n");

    // Define the target Hamiltonian: 1D Heisenberg ferromagnet
    let n_modes = 8_usize;
    let j_exchange = 1.0_f64; // dimensionless exchange coupling
    let h_ext = 0.2_f64; // dimensionless external field
    let spin_s = 1.0_f64;

    let params = MagnonHamiltonianParams::new(n_modes, j_exchange, h_ext, spin_s)?;
    let target_freqs = params.reference_magnon_frequencies();

    println!("Target Hamiltonian: 1D Heisenberg chain, {n_modes} modes");
    println!("  J = {j_exchange:.2}, H_ext = {h_ext:.2}, S = {spin_s:.1}");
    println!("  Dispersion: ω_k = J·S·(1 − cos k) + H_ext");
    println!("  Target magnon frequencies:");
    for (i, f) in target_freqs.iter().enumerate() {
        let k_frac = (i + 1) as f64 / n_modes as f64;
        println!("    k[{i}] = {k_frac:.3}π  →  ω = {f:.4}");
    }

    // Build the quantum-classical optimizer
    // NN: [3] → 16 × 16 → [2] (Tanh hidden, Linear output)
    // Maps (k/π, J_norm, H_norm) → (A_k_raw, B_k_raw) → Bogoliubov ε_k
    let mut optimizer = QuantumClassicalOptimizer::new(
        params, 16,   // hidden dim
        2,    // depth (two hidden layers)
        0.01, // Adam learning rate
    );

    let (init_loss, _) = optimizer.compute_loss_and_gradients()?;
    println!("\nInitial loss: {:.4e}", init_loss);

    println!("Training for 50 steps (Adam + central finite-difference gradients)...");
    let result = optimizer.train(50)?;

    println!("Final loss:           {:.4e}", result.final_loss);
    println!(
        "Ground-state energy:  {:.4}  (Σ_k (ε_k − A_k)/2 ≤ 0)",
        result.ground_state_energy
    );
    println!(
        "Loss reduction:       {:.0}×",
        if result.final_loss > 0.0 {
            init_loss / result.final_loss
        } else {
            f64::INFINITY
        }
    );

    println!("\nComparison: target vs learned magnon frequencies:");
    println!(
        "  {:^5}  {:^10}  {:^10}  {:^10}",
        "k", "target ω", "learned ε_k", "rel. error"
    );
    let avg_err: f64 = result
        .target_frequencies
        .iter()
        .zip(result.final_magnon_frequencies.iter())
        .enumerate()
        .map(|(i, (&target, &learned))| {
            let rel_err = (learned - target).abs() / target.abs().max(1e-10);
            println!(
                "  k[{i}]  {:10.4}  {:11.4}  {:9.1}%",
                target,
                learned,
                rel_err * 100.0
            );
            rel_err
        })
        .sum::<f64>()
        / n_modes as f64;

    println!("\nMean relative error: {:.1}%", avg_err * 100.0);

    // Show loss history convergence
    println!("\nLoss curve (every 10 steps):");
    for (i, &loss) in result.loss_history.iter().enumerate() {
        if i % 10 == 0 || i == result.loss_history.len() - 1 {
            println!("  step {:3}: loss = {:.4e}", i + 1, loss);
        }
    }

    // Sanity check: ground state energy must be ≤ 0 for Bogoliubov transformation
    let gs_ok = result.ground_state_energy <= 0.0 + 1e-10;
    println!("\nBogoliubov stability (E_gs ≤ 0): {}", gs_ok);
    let freqs_ok = result.final_magnon_frequencies.iter().all(|&f| f > 0.0);
    println!("All magnon frequencies positive:    {}", freqs_ok);

    Ok(())
}
