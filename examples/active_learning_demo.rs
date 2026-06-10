//! Active Learning with Uncertainty Sampling
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Machine Learning / Sample-efficient learning
//! **Physics**: Query-efficient parameter fitting for expensive oracles
//!
//! ## Background
//!
//! When experimental data or DFT calls are expensive, naïve random sampling
//! wastes queries. Active learning picks the next query to maximize information
//! gain — typically where the model is most uncertain. The classic strategy:
//!
//!   1. Train an ensemble of MLPs from random subsets (bootstrap)
//!   2. For each candidate point, estimate uncertainty as ensemble std dev
//!   3. Query the point with the highest uncertainty
//!   4. Add the new (x, y) to the training set; retrain
//!
//! This example fits a "expensive" 1D function f(x) = sin(5x) · exp(-x²) using
//! three strategies on the same query budget, and compares final MSE.
//!
//! ## References
//! - Cohn, Atlas, Ladner, MLJ 15, 201 (1994) — query-by-committee
//! - Settles, "Active Learning Literature Survey", U. Wisconsin CS-TR-1648 (2009)
//! - Bernstein et al., npj Comput. Mater. 5, 99 (2019) — ML potentials with AL

use spintronics::autodiff::Activation;
use spintronics::prelude::*;

fn target_function(x: f64) -> f64 {
    (5.0 * x).sin() * (-x * x).exp()
}

fn evaluate_mse(learner: &ActiveLearner, test_x: &[f64], test_y: &[f64]) -> f64 {
    let mut acc = 0.0_f64;
    for (x, &y) in test_x.iter().zip(test_y.iter()) {
        let (mean, _std) = learner.predict_with_uncertainty(&[*x]).unwrap();
        let err = mean[0] - y;
        acc += err * err;
    }
    acc / test_x.len() as f64
}

fn run_strategy(
    strategy: QueryStrategy,
    label: &str,
) -> std::result::Result<(usize, f64), Box<dyn std::error::Error>> {
    println!("\n--- {label} ---\n");

    // Candidate pool: 200 points uniformly in [-2, 2]
    let candidate_pool: Vec<Vec<f64>> = (0..200)
        .map(|i| vec![-2.0 + 4.0 * (i as f64) / 199.0])
        .collect();

    let config = ActiveLearningConfig {
        n_committee: 5,
        n_initial_samples: 5,
        n_query_iterations: 15,
        n_inner_train: 80,
        lr: 1e-2,
        query_strategy: strategy,
    };

    let layer_sizes = [1, 16, 16, 1];
    let activations = [Activation::Tanh, Activation::Tanh, Activation::Linear];
    let mut learner = ActiveLearner::new(&layer_sizes, &activations, config, 42)?;

    // Oracle: evaluate target function
    let oracle = |x: &[f64]| -> Vec<f64> { vec![target_function(x[0])] };

    let result = learner.fit(oracle, &candidate_pool)?;

    // Build test grid (independent of training pool indices)
    let test_x: Vec<f64> = (0..101).map(|i| -2.0 + 4.0 * (i as f64) / 100.0).collect();
    let test_y: Vec<f64> = test_x.iter().copied().map(target_function).collect();
    let test_x_slices: Vec<Vec<f64>> = test_x.iter().map(|&x| vec![x]).collect();
    let mse = evaluate_mse(
        &learner,
        &test_x_slices.iter().map(|v| v[0]).collect::<Vec<f64>>(),
        &test_y,
    );

    println!("  Queries used:          {}", result.n_queries_used);
    println!("  Final training loss:   {:.4e}", result.final_loss);
    println!("  Independent-test MSE:  {mse:.4e}");
    println!(
        "  Queried indices:       {:?}",
        &result.queried_indices[..result.queried_indices.len().min(8)]
    );

    Ok((result.n_queries_used, mse))
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=============================================================");
    println!("  Active Learning: Uncertainty Sampling vs Random Baseline");
    println!("=============================================================");
    println!();
    println!("  Target: f(x) = sin(5x) · exp(-x²)  on  x ∈ [-2, 2]");
    println!("  Pool: 200 candidate points uniformly spaced");
    println!("  Query budget: 15 oracle calls + 5 initial random samples");
    println!("  Committee: 5 MLPs of [1 → 16 → 16 → 1]");

    let (n_rand, mse_rand) = run_strategy(QueryStrategy::RandomBaseline, "Random Baseline")?;
    let (n_qbc, mse_qbc) = run_strategy(QueryStrategy::QueryByCommittee, "Query By Committee")?;
    let (n_unc, mse_unc) =
        run_strategy(QueryStrategy::UncertaintySampling, "Uncertainty Sampling")?;

    println!("\n--- Comparison Summary ---\n");
    println!(
        "  {:>22}  {:>10}  {:>14}",
        "Strategy", "queries", "test MSE"
    );
    println!("  {}", "-".repeat(50));
    println!(
        "  {:>22}  {:>10}  {:>14.4e}",
        "RandomBaseline", n_rand, mse_rand
    );
    println!(
        "  {:>22}  {:>10}  {:>14.4e}",
        "QueryByCommittee", n_qbc, mse_qbc
    );
    println!(
        "  {:>22}  {:>10}  {:>14.4e}",
        "UncertaintySampling", n_unc, mse_unc
    );

    println!("\n=============================================================");
    println!("  Done. Active strategies typically achieve lower test MSE for");
    println!("  the same query budget by focusing on high-uncertainty regions.");
    println!("=============================================================\n");

    Ok(())
}
