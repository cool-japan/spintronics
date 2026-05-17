//! Bayesian Optimisation for Materials Parameter Search
//!
//! **Difficulty**: ⭐⭐⭐⭐
//! **Category**: Machine Learning / Sample-efficient search
//! **Physics**: Maximise a target observable over a 2D parameter pool
//!
//! ## Background
//!
//! Many spintronics figures of merit (skyrmion radius, switching speed,
//! spin-orbit torque efficiency) are produced by expensive DFT or
//! micromagnetic simulations. Random or grid search wastes evaluations.
//! **Bayesian Optimisation** with a Gaussian Process surrogate and an
//! Expected Improvement acquisition function balances exploration against
//! exploitation, typically reaching the optimum in 10× fewer queries than
//! random search.
//!
//! This demo:
//!   1. Defines a synthetic 2D oracle that is *deceptively* multimodal —
//!      `f(J, K) = exp(−10·((J − 0.7)² + (K − 0.4)²)) + 0.6·exp(−10·(J² + K²))`
//!      with a stronger peak at `(0.7, 0.4)`.
//!   2. Generates a 25 × 25 grid candidate pool on `[0, 1] × [0, 1]`.
//!   3. Runs a [`BayesianOptimizer`] with 5 random seed samples and 15 EI
//!      iterations to find the maximum.
//!   4. Reports the best `(J, K)` found and the distance to the true peak.
//!
//! ## References
//! - Brochu, Cora & de Freitas, "A Tutorial on Bayesian Optimization of
//!   Expensive Cost Functions, with Application to Active User Modeling and
//!   Hierarchical Reinforcement Learning", arXiv:1012.2599 (2010).
//! - Frazier, "A Tutorial on Bayesian Optimization", arXiv:1807.02811 (2018).

use spintronics::autodiff::{
    AcquisitionStrategy, BayesianOptConfig, BayesianOptimizer, GaussianProcess, GpConfig,
};

fn oracle(x: &[f64]) -> f64 {
    let j = x[0];
    let k = x[1];
    let peak_a = (-10.0 * ((j - 0.7).powi(2) + (k - 0.4).powi(2))).exp();
    let peak_b = 0.6 * (-10.0 * (j * j + k * k)).exp();
    peak_a + peak_b
}

fn main() {
    println!("Bayesian Optimisation — 2D Materials Parameter Search\n");

    // 1. Build the candidate pool (25 × 25 grid on [0, 1]²).
    let mut pool: Vec<Vec<f64>> = Vec::new();
    let resolution = 25;
    for i in 0..resolution {
        for j in 0..resolution {
            let xj = (i as f64) / ((resolution - 1) as f64);
            let xk = (j as f64) / ((resolution - 1) as f64);
            pool.push(vec![xj, xk]);
        }
    }
    println!("Candidate pool: {} grid points on [0,1]²", pool.len());

    // 2. Configure the optimiser.
    let mut cfg = BayesianOptConfig {
        n_initial_samples: 5,
        n_iterations: 15,
        gp_config: GaussianProcess::config_default(),
        acquisition: AcquisitionStrategy::ExpectedImprovement,
        maximize: true,
    };
    cfg.gp_config = GpConfig {
        length_scale: 0.18,
        signal_variance: 1.0,
        noise_variance: 1.0e-6,
        jitter: 1.0e-8,
    };
    let mut opt = BayesianOptimizer::new(cfg).expect("optimizer construction");

    // 3. Run the BO loop.
    println!("Running BO with 5 seed samples + 15 EI iterations...");
    let result = opt.optimize(oracle, &pool, 0x0BEE_FACE).expect("BO run");
    println!("Total evaluations: {}", result.n_evaluations);
    println!(
        "Best (J, K)       : ({:.4}, {:.4})",
        result.best_x[0], result.best_x[1]
    );
    println!("Best f            : {:.6}", result.best_y);

    // 4. Compare to the true optimum.
    let true_jk = [0.7_f64, 0.4_f64];
    let true_f = oracle(&true_jk);
    let dist =
        ((result.best_x[0] - true_jk[0]).powi(2) + (result.best_x[1] - true_jk[1]).powi(2)).sqrt();
    println!(
        "\nTrue optimum     : (J*, K*) = (0.7, 0.4), f* = {:.6}",
        true_f
    );
    println!("|x − x*|         : {:.4}", dist);
    println!(
        "Relative gap     : {:.4}",
        (true_f - result.best_y).abs() / true_f.abs().max(1.0e-12),
    );
}
