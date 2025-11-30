//! Physical Reservoir Computing with Magnon Dynamics
//!
//! Demonstrates neuromorphic computing using spin wave non-linearity.
//! The magnon system acts as a non-linear, high-dimensional dynamical system
//! for computation, with only the readout layer trained.

use spintronics::ai::MagnonReservoir;
use spintronics::vector3::Vector3;

fn main() {
    println!("=== Physical Reservoir Computing with Magnons ===\n");

    // Create reservoir
    let n_spins = 30;
    let n_input = 5;
    let n_readout = 10;

    let mut reservoir = MagnonReservoir::uniform_distribution(
        n_spins, n_input, n_readout, 1.0e-20, // Exchange coupling
        0.01,    // Damping
    );

    println!("Reservoir Configuration:");
    println!("  Total spins: {}", n_spins);
    println!("  Input nodes: {}", n_input);
    println!("  Readout nodes: {}", n_readout);
    println!(
        "  State dimension: {} (3 components × {} nodes)",
        n_readout * 3,
        n_readout
    );
    println!();

    // Task 1: Constant function approximation
    println!("--- Task 1: Constant Function ---\n");

    let inputs_const = vec![
        vec![0.1, 0.2, 0.15, 0.3, 0.1],
        vec![0.3, 0.1, 0.25, 0.2, 0.15],
        vec![0.2, 0.3, 0.1, 0.15, 0.25],
        vec![0.15, 0.25, 0.3, 0.1, 0.2],
        vec![0.25, 0.15, 0.2, 0.3, 0.05],
    ];

    let targets_const = vec![1.0, 1.0, 1.0, 1.0, 1.0];

    println!("Training to output constant 1.0 regardless of input...");

    let h_ext = Vector3::new(0.0, 0.0, 1.0e5);
    let evolution_steps = 10;
    let dt = 1.0e-13;

    let mse_result = reservoir.train(&inputs_const, &targets_const, h_ext, evolution_steps, dt);

    match mse_result {
        Ok(mse) => {
            println!("  Training MSE: {:.6}", mse);
            println!();

            // Test the trained reservoir
            println!("Testing on training data:");
            for (i, input) in inputs_const.iter().enumerate() {
                reservoir.reset();
                let output = reservoir.process(input, h_ext, evolution_steps, dt);
                println!(
                    "  Input {:?} → Output: {:.4} (Target: {:.1})",
                    &input[..2],
                    output,
                    targets_const[i]
                );
            }
            println!();
        }
        Err(e) => {
            println!("  Training error: {}", e);
            println!();
        }
    }

    // Task 2: Linear function
    println!("--- Task 2: Linear Function ---\n");

    let mut reservoir2 = MagnonReservoir::uniform_distribution(25, 4, 8, 1.0e-20, 0.01);

    let inputs_linear = vec![
        vec![0.0, 0.0, 0.0, 0.0],
        vec![0.1, 0.0, 0.0, 0.0],
        vec![0.2, 0.0, 0.0, 0.0],
        vec![0.3, 0.0, 0.0, 0.0],
        vec![0.0, 0.1, 0.0, 0.0],
        vec![0.0, 0.2, 0.0, 0.0],
    ];

    // Target: sum of first two inputs
    let targets_linear: Vec<f64> = inputs_linear.iter().map(|inp| inp[0] + inp[1]).collect();

    println!("Training to compute f(x) = x[0] + x[1]...");

    let mse_result2 = reservoir2.train(&inputs_linear, &targets_linear, h_ext, 8, dt);

    match mse_result2 {
        Ok(mse) => {
            println!("  Training MSE: {:.6}", mse);
            println!();

            println!("Testing:");
            for (i, input) in inputs_linear.iter().enumerate() {
                reservoir2.reset();
                let output = reservoir2.process(input, h_ext, 8, dt);
                println!(
                    "  [{:.1}, {:.1}] → {:.4} (Target: {:.2})",
                    input[0], input[1], output, targets_linear[i]
                );
            }
            println!();
        }
        Err(e) => {
            println!("  Training error: {}", e);
            println!();
        }
    }

    // Task 3: Non-linear function (XOR-like)
    println!("--- Task 3: Non-linear Pattern Recognition ---\n");

    let mut reservoir3 = MagnonReservoir::uniform_distribution(40, 2, 12, 1.0e-20, 0.01);

    let inputs_nonlinear = vec![
        vec![0.0, 0.0], // Both low
        vec![0.0, 0.5], // First low, second high
        vec![0.5, 0.0], // First high, second low
        vec![0.5, 0.5], // Both high
        vec![0.1, 0.1], // Both slightly low
        vec![0.4, 0.4], // Both slightly high
    ];

    // XOR-like: high output when inputs differ
    let targets_nonlinear = vec![0.0, 1.0, 1.0, 0.0, 0.0, 0.0];

    println!("Training XOR-like function (output high when inputs differ)...");

    let mse_result3 = reservoir3.train(&inputs_nonlinear, &targets_nonlinear, h_ext, 12, dt);

    match mse_result3 {
        Ok(mse) => {
            println!("  Training MSE: {:.6}", mse);
            println!();

            println!("Testing:");
            for (i, input) in inputs_nonlinear.iter().enumerate() {
                reservoir3.reset();
                let output = reservoir3.process(input, h_ext, 12, dt);
                println!(
                    "  [{:.1}, {:.1}] → {:.4} (Target: {:.1})",
                    input[0], input[1], output, targets_nonlinear[i]
                );
            }
            println!();
        }
        Err(e) => {
            println!("  Training error: {}", e);
            println!();
        }
    }

    // Information about the reservoir
    println!("--- Reservoir Properties ---\n");

    reservoir.reset();
    let initial_state = reservoir.extract_state();
    println!("Initial state (all spins along +z):");
    println!("  State vector dimension: {}", initial_state.len());
    println!(
        "  First few components: [{:.3}, {:.3}, {:.3}, ...]",
        initial_state[0], initial_state[1], initial_state[2]
    );
    println!();

    // Evolve and show state change
    let test_input = vec![0.2, 0.3];
    reservoir.process(&test_input, h_ext, 10, dt);
    let evolved_state = reservoir.extract_state();

    println!(
        "After processing input [{:.1}, {:.1}]:",
        test_input[0], test_input[1]
    );
    println!(
        "  First few components: [{:.3}, {:.3}, {:.3}, ...]",
        evolved_state[0], evolved_state[1], evolved_state[2]
    );
    println!();

    println!("=== Summary ===");
    println!();
    println!("Physical Reservoir Computing leverages the non-linear dynamics of spin waves:");
    println!();
    println!("1. Input Layer:");
    println!("   - Scalar inputs encoded as magnetic field modulation");
    println!("   - Applied to specific spins in the chain");
    println!();
    println!("2. Reservoir Layer (physical substrate):");
    println!("   - Spin chain with LLG dynamics (non-linear)");
    println!("   - Exchange coupling creates spatial interactions");
    println!("   - Gilbert damping provides dissipation");
    println!("   - High-dimensional state space (N_spins × 3 components)");
    println!();
    println!("3. Readout Layer (trained):");
    println!("   - Linear combination of reservoir states");
    println!("   - Weights trained via gradient descent");
    println!("   - Only layer that requires training");
    println!();
    println!("Advantages:");
    println!("  - No need to train the complex reservoir dynamics");
    println!("  - Physical hardware naturally provides non-linearity");
    println!("  - Energy-efficient computation");
    println!("  - Potential for ultra-fast processing (GHz rates)");
    println!();
    println!("Applications:");
    println!("  - Time series prediction");
    println!("  - Pattern classification");
    println!("  - Signal processing");
    println!("  - Neuromorphic computing");
}
