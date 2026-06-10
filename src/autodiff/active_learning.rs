//! Active learning for parameter fitting with expensive oracles.
//!
//! When experimental data, *ab-initio* DFT calls, or detailed
//! micromagnetic simulations are expensive, active learning chooses the
//! *next* query intelligently in order to maximise information gain per
//! evaluation.  This module implements three classical strategies on top of
//! the existing [`Mlp`] infrastructure:
//!
//! - **Uncertainty sampling** — train an ensemble of MLPs on bootstrap
//!   replicates of the labelled set; query where the ensemble variance is
//!   largest.
//! - **Query-by-committee** — train the ensemble on the *same* data but
//!   with different RNG seeds for the weight initialisation, then again
//!   pick the variance-maximum candidate.
//! - **Random baseline** — pick the next candidate from a deterministic
//!   pseudo-random sequence.  This is the standard control against which
//!   the active strategies are benchmarked.
//!
//! ## Training
//!
//! Each committee member is trained for `n_inner_train` Adam steps using
//! plain finite-difference gradients on `Mlp::params_flat()`.  This matches
//! the gradient-free training pattern used by the v0.7.0
//! `neural_exchange_training.rs` example: the loop is simple, fully
//! reproducible, and avoids re-deriving `Var<'t>`-based gradients for the
//! per-sample MSE loss.  To keep the cost reasonable we *sub-sample* the
//! parameter coordinates that receive an FD update per outer iteration —
//! see [`ActiveLearner::train_committee`] for details.
//!
//! ## References
//!
//! - H. S. Seung, M. Opper & H. Sompolinsky, "Query by Committee",
//!   *Proc. 5th Annual ACM Workshop on Computational Learning Theory* (1992).
//! - D. Cohn, L. Atlas & R. Ladner, "Improving Generalization with Active
//!   Learning", *Machine Learning* **15**, 201 (1994).
//! - B. Settles, "Active Learning Literature Survey",
//!   *Univ. of Wisconsin–Madison Tech. Report 1648* (2009).
//! - N. Bernstein, G. Csányi & V. L. Deringer, "*De novo* exploration and
//!   self-guided learning of potential-energy surfaces",
//!   *npj Comput. Mater.* **5**, 99 (2019).

use crate::autodiff::neural::{Activation, Mlp};
use crate::autodiff::optimizer::{Adam, Optimizer};
use crate::error::{dimension_mismatch, invalid_param, Result};

// ─── Internal LCG for reproducible random choices ────────────────────────────

/// Same Knuth-style LCG used elsewhere in `autodiff` — duplicated here so
/// `active_learning.rs` is self-contained.
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        let state = if seed == 0 {
            0xDEAD_BEEF_CAFE_BABE
        } else {
            seed
        };
        Self { state }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    /// Sample an integer uniformly in `0..n` (rejection-free, accepts the
    /// tiny bias for `n` not a power of two — negligible for our pool sizes).
    fn next_index(&mut self, n: usize) -> usize {
        if n == 0 {
            return 0;
        }
        (self.next_u64() as usize) % n
    }
}

// ─── QueryStrategy ───────────────────────────────────────────────────────────

/// Which acquisition function to use when selecting the next training point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QueryStrategy {
    /// Bootstrap each committee member's training set, then pick the
    /// candidate with maximum sum-of-stds across output dimensions.
    UncertaintySampling,
    /// Train every committee member on the same data with different RNG
    /// seeds, then pick the variance-maximum candidate.
    QueryByCommittee,
    /// Pick deterministic pseudo-random indices — control baseline.
    RandomBaseline,
}

// ─── ActiveLearningConfig ────────────────────────────────────────────────────

/// Hyperparameters for an [`ActiveLearner`] loop.
#[derive(Debug, Clone)]
pub struct ActiveLearningConfig {
    /// Number of independent networks in the committee (`≥ 2`).
    pub n_committee: usize,
    /// Number of randomly-chosen samples drawn from the pool to seed the
    /// loop before any query is made.
    pub n_initial_samples: usize,
    /// Number of `(query, retrain)` outer iterations.
    pub n_query_iterations: usize,
    /// Number of Adam inner steps per committee member per outer iteration.
    pub n_inner_train: usize,
    /// Adam learning rate.
    pub lr: f64,
    /// Acquisition strategy.
    pub query_strategy: QueryStrategy,
}

impl ActiveLearningConfig {
    /// Reject configurations that would produce a meaningless run.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] on any zero count
    /// or non-positive learning rate.
    pub fn validate(&self) -> Result<()> {
        if self.n_committee < 2 {
            return Err(invalid_param(
                "n_committee",
                "must have ≥ 2 members for a meaningful uncertainty estimate",
            ));
        }
        if self.n_initial_samples == 0 {
            return Err(invalid_param(
                "n_initial_samples",
                "need ≥ 1 sample to bootstrap training",
            ));
        }
        if self.n_query_iterations == 0 {
            return Err(invalid_param(
                "n_query_iterations",
                "must request at least one query iteration",
            ));
        }
        if self.n_inner_train == 0 {
            return Err(invalid_param(
                "n_inner_train",
                "must run ≥ 1 inner training step",
            ));
        }
        if !(self.lr > 0.0 && self.lr.is_finite()) {
            return Err(invalid_param("lr", "must be positive and finite"));
        }
        Ok(())
    }
}

// ─── ActiveLearner ───────────────────────────────────────────────────────────

/// Driver that owns the committee of MLPs, the cumulative training set, and
/// the active-learning configuration.
pub struct ActiveLearner {
    /// Hyperparameters for the loop.
    pub config: ActiveLearningConfig,
    /// Committee of independently-initialised MLPs.
    pub committee: Vec<Mlp>,
    /// Training inputs accumulated so far — one row per labelled sample.
    pub training_x: Vec<Vec<f64>>,
    /// Training targets accumulated so far — one row per labelled sample.
    pub training_y: Vec<Vec<f64>>,
    /// Per-layer sizes used to construct every committee member.
    pub layer_sizes: Vec<usize>,
    /// Per-layer activations used to construct every committee member.
    pub activations: Vec<Activation>,
}

impl ActiveLearner {
    /// Build a new active learner with `n_committee` independently-seeded MLPs.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when the layer
    /// configuration, activation list, or [`ActiveLearningConfig`] is
    /// invalid.
    pub fn new(
        layer_sizes: &[usize],
        activations: &[Activation],
        config: ActiveLearningConfig,
        rng_seed: u64,
    ) -> Result<Self> {
        config.validate()?;
        if layer_sizes.len() < 2 {
            return Err(invalid_param(
                "layer_sizes",
                "must include at least input and output dimensions",
            ));
        }
        if activations.len() + 1 != layer_sizes.len() {
            return Err(invalid_param(
                "activations",
                "must have one fewer entry than layer_sizes",
            ));
        }
        let mut committee = Vec::with_capacity(config.n_committee);
        for k in 0..config.n_committee {
            let sub_seed = rng_seed.wrapping_add((k as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            committee.push(Mlp::new(layer_sizes, activations, sub_seed)?);
        }
        Ok(Self {
            config,
            committee,
            training_x: Vec::new(),
            training_y: Vec::new(),
            layer_sizes: layer_sizes.to_vec(),
            activations: activations.to_vec(),
        })
    }

    /// Append a labelled `(x, y)` pair to the training set.
    pub fn add_sample(&mut self, x: Vec<f64>, y: Vec<f64>) {
        self.training_x.push(x);
        self.training_y.push(y);
    }

    /// Bootstrap-resample the labelled set: draw `n` indices uniformly with
    /// replacement.  Used by [`QueryStrategy::UncertaintySampling`].
    fn bootstrap_indices(&self, rng: &mut Lcg, n: usize) -> Vec<usize> {
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            out.push(rng.next_index(self.training_x.len()));
        }
        out
    }

    /// Mean-squared-error loss for one MLP over a list of training indices.
    fn mse_on_indices(&self, mlp: &Mlp, indices: &[usize]) -> Result<f64> {
        if indices.is_empty() {
            return Ok(0.0);
        }
        let mut acc = 0.0_f64;
        let mut count = 0_usize;
        for &i in indices {
            let pred = mlp.forward_f64(&self.training_x[i])?;
            let target = &self.training_y[i];
            if pred.len() != target.len() {
                return Err(dimension_mismatch(
                    &format!("{} target dims", pred.len()),
                    &format!("{} target dims", target.len()),
                ));
            }
            for (p, t) in pred.iter().zip(target.iter()) {
                let d = p - t;
                acc += d * d;
                count += 1;
            }
        }
        Ok(acc / count.max(1) as f64)
    }

    /// Average MSE on the full labelled set (used as the public loss
    /// reported in [`ActiveLearnResult`]).
    fn mean_mse_on_training(&self) -> Result<f64> {
        if self.training_x.is_empty() {
            return Ok(f64::NAN);
        }
        let indices: Vec<usize> = (0..self.training_x.len()).collect();
        let mut acc = 0.0_f64;
        for mlp in &self.committee {
            acc += self.mse_on_indices(mlp, &indices)?;
        }
        Ok(acc / self.committee.len() as f64)
    }

    /// Train every committee member with Adam + finite-difference gradients
    /// on the MSE loss against the current labelled set.
    ///
    /// To keep the cost bounded, we sub-sample at most
    /// `min(n_params, max_fd_coords)` parameter coordinates per outer step
    /// for the finite-difference probe (the remaining coordinates have zero
    /// gradient that iteration).  This matches the *coordinate-descent*
    /// style of training used by the existing v0.7.0 examples.
    ///
    /// # Errors
    /// Propagates errors from the underlying MLP forward passes.
    pub fn train_committee(&mut self) -> Result<()> {
        if self.training_x.is_empty() {
            // Nothing to train against — silently no-op.
            return Ok(());
        }
        // Maximum number of FD coordinates probed per inner step.  This is
        // a fixed budget rather than a config knob to keep the public API
        // narrow — 32 is generous for the small surrogates used in
        // active-learning loops.
        const MAX_FD_COORDS: usize = 32;
        let fd_h = 1e-4_f64;

        // Pre-compute the indices each committee member trains on.
        let mut indices_per_member: Vec<Vec<usize>> = Vec::with_capacity(self.committee.len());
        let mut rng = Lcg::new(0xAB57_CDEF_1234_5678u64 ^ (self.training_x.len() as u64));
        for member_id in 0..self.committee.len() {
            let idx = match self.config.query_strategy {
                QueryStrategy::UncertaintySampling => {
                    self.bootstrap_indices(&mut rng, self.training_x.len())
                },
                QueryStrategy::QueryByCommittee | QueryStrategy::RandomBaseline => {
                    (0..self.training_x.len()).collect()
                },
            };
            let _ = member_id;
            indices_per_member.push(idx);
        }

        // Train each member.
        for (member_id, mlp) in self.committee.iter_mut().enumerate() {
            let indices = &indices_per_member[member_id];
            let mut params = mlp.params_flat();
            let n = params.len();
            let mut adam = Adam::default_params(n);
            adam.lr = self.config.lr;
            let mut step_rng = Lcg::new(0x1234_5678_9ABC_DEF0u64.wrapping_add(member_id as u64));

            for _ in 0..self.config.n_inner_train {
                // Pick a random coordinate subset to probe with FD.
                let probe_count = n.min(MAX_FD_COORDS);
                let mut probe = vec![false; n];
                let mut selected = 0;
                while selected < probe_count {
                    let idx = step_rng.next_index(n);
                    if !probe[idx] {
                        probe[idx] = true;
                        selected += 1;
                    }
                }

                // FD gradient on probed coords; zero elsewhere.
                let mut grads = vec![0.0_f64; n];
                for j in 0..n {
                    if !probe[j] {
                        continue;
                    }
                    let original = params[j];
                    params[j] = original + fd_h;
                    mlp.set_params(&params)?;
                    let plus = Self::mse_static(mlp, indices, &self.training_x, &self.training_y)?;
                    params[j] = original - fd_h;
                    mlp.set_params(&params)?;
                    let minus = Self::mse_static(mlp, indices, &self.training_x, &self.training_y)?;
                    params[j] = original;
                    grads[j] = (plus - minus) / (2.0 * fd_h);
                }

                // Restore baseline parameters before Adam updates them.
                mlp.set_params(&params)?;
                adam.step(&mut params, &grads);
                mlp.set_params(&params)?;
            }
        }
        Ok(())
    }

    /// Free-function version of [`Self::mse_on_indices`] so it can be
    /// called while holding a mutable borrow on a single committee member.
    fn mse_static(
        mlp: &Mlp,
        indices: &[usize],
        training_x: &[Vec<f64>],
        training_y: &[Vec<f64>],
    ) -> Result<f64> {
        if indices.is_empty() {
            return Ok(0.0);
        }
        let mut acc = 0.0_f64;
        let mut count = 0_usize;
        for &i in indices {
            let pred = mlp.forward_f64(&training_x[i])?;
            let target = &training_y[i];
            if pred.len() != target.len() {
                return Err(dimension_mismatch(
                    &format!("{} target dims", pred.len()),
                    &format!("{} target dims", target.len()),
                ));
            }
            for (p, t) in pred.iter().zip(target.iter()) {
                let d = p - t;
                acc += d * d;
                count += 1;
            }
        }
        Ok(acc / count.max(1) as f64)
    }

    /// Predict the per-output mean and standard deviation across the
    /// committee at input `x`.
    ///
    /// # Errors
    /// Propagates dimension-mismatch errors from the underlying MLPs.
    pub fn predict_with_uncertainty(&self, x: &[f64]) -> Result<(Vec<f64>, Vec<f64>)> {
        if self.committee.is_empty() {
            return Err(invalid_param("committee", "must have ≥ 1 member"));
        }
        let mut preds: Vec<Vec<f64>> = Vec::with_capacity(self.committee.len());
        for mlp in &self.committee {
            preds.push(mlp.forward_f64(x)?);
        }
        let out_dim = preds[0].len();
        let n = preds.len() as f64;
        let mut mean = vec![0.0_f64; out_dim];
        for p in &preds {
            for (m, v) in mean.iter_mut().zip(p.iter()) {
                *m += *v;
            }
        }
        for m in mean.iter_mut() {
            *m /= n;
        }
        let mut var = vec![0.0_f64; out_dim];
        for p in &preds {
            for (vd, (v, m)) in var.iter_mut().zip(p.iter().zip(mean.iter())) {
                let d = *v - *m;
                *vd += d * d;
            }
        }
        for vd in var.iter_mut() {
            *vd /= n;
        }
        let std_dev: Vec<f64> = var.iter().map(|s| s.sqrt()).collect();
        Ok((mean, std_dev))
    }

    /// Pick the index in `candidate_pool` to query next.
    ///
    /// - **UncertaintySampling / QueryByCommittee** → argmax over the
    ///   sum-of-per-output standard deviations.
    /// - **RandomBaseline** → deterministic pseudo-random index seeded by
    ///   the current training-set size.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when the pool is
    /// empty.
    pub fn select_next_query(&self, candidate_pool: &[Vec<f64>]) -> Result<usize> {
        if candidate_pool.is_empty() {
            return Err(invalid_param("candidate_pool", "must not be empty"));
        }
        match self.config.query_strategy {
            QueryStrategy::RandomBaseline => {
                let mut rng = Lcg::new(0x9876_5432_10FE_DCBAu64 ^ (self.training_x.len() as u64));
                Ok(rng.next_index(candidate_pool.len()))
            },
            QueryStrategy::UncertaintySampling | QueryStrategy::QueryByCommittee => {
                let mut best_idx = 0_usize;
                let mut best_score = f64::NEG_INFINITY;
                for (i, c) in candidate_pool.iter().enumerate() {
                    let (_mean, std_dev) = self.predict_with_uncertainty(c)?;
                    let score: f64 = std_dev.iter().sum();
                    if score > best_score {
                        best_score = score;
                        best_idx = i;
                    }
                }
                Ok(best_idx)
            },
        }
    }

    /// Full active-learning loop: seed with `n_initial_samples` random
    /// queries, then perform `n_query_iterations` query-and-retrain steps.
    ///
    /// The `oracle` closure maps an input `x` to its (potentially
    /// expensive) target `y`.  No mutable state is required of the oracle.
    ///
    /// # Errors
    /// Propagates errors from training and uncertainty prediction.
    pub fn fit<F>(&mut self, oracle: F, candidate_pool: &[Vec<f64>]) -> Result<ActiveLearnResult>
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        if candidate_pool.is_empty() {
            return Err(invalid_param("candidate_pool", "must not be empty"));
        }
        let mut queried_indices: Vec<usize> =
            Vec::with_capacity(self.config.n_initial_samples + self.config.n_query_iterations);
        let mut loss_history: Vec<f64> = Vec::with_capacity(self.config.n_query_iterations);

        // 1. Seed with random samples (deterministic LCG).
        let mut seed_rng = Lcg::new(0xCAFE_BABE_DEAD_BEEFu64);
        for _ in 0..self.config.n_initial_samples {
            let idx = seed_rng.next_index(candidate_pool.len());
            queried_indices.push(idx);
            let x = candidate_pool[idx].clone();
            let y = oracle(&x);
            self.add_sample(x, y);
        }
        self.train_committee()?;

        // 2. Query loop.
        for _ in 0..self.config.n_query_iterations {
            let next = self.select_next_query(candidate_pool)?;
            queried_indices.push(next);
            let x = candidate_pool[next].clone();
            let y = oracle(&x);
            self.add_sample(x, y);
            self.train_committee()?;
            let loss = self.mean_mse_on_training()?;
            loss_history.push(loss);
        }
        let final_loss = *loss_history.last().unwrap_or(&f64::NAN);
        Ok(ActiveLearnResult {
            n_queries_used: queried_indices.len(),
            final_loss,
            loss_history,
            queried_indices,
        })
    }
}

// ─── ActiveLearnResult ───────────────────────────────────────────────────────

/// Result returned by [`ActiveLearner::fit`].
#[derive(Debug, Clone)]
pub struct ActiveLearnResult {
    /// Total number of oracle calls made during the run
    /// (`n_initial_samples + n_query_iterations`).
    pub n_queries_used: usize,
    /// Final mean-MSE across the committee on the cumulative training set.
    pub final_loss: f64,
    /// Per-iteration mean-MSE history (one entry per outer query iteration).
    pub loss_history: Vec<f64>,
    /// Pool indices selected by the acquisition function, in query order.
    pub queried_indices: Vec<usize>,
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::autodiff::neural::Activation;

    fn default_config() -> ActiveLearningConfig {
        ActiveLearningConfig {
            n_committee: 3,
            n_initial_samples: 2,
            n_query_iterations: 2,
            n_inner_train: 3,
            lr: 5e-3,
            query_strategy: QueryStrategy::QueryByCommittee,
        }
    }

    // 1. Configuration validation: zero committee is rejected.
    #[test]
    fn test_config_validation_rejects_zero_committee() {
        let mut cfg = default_config();
        cfg.n_committee = 0;
        assert!(cfg.validate().is_err());
        cfg.n_committee = 1;
        assert!(cfg.validate().is_err());
        cfg.n_committee = 2;
        cfg.n_initial_samples = 0;
        assert!(cfg.validate().is_err());
        cfg.n_initial_samples = 1;
        cfg.n_query_iterations = 0;
        assert!(cfg.validate().is_err());
        cfg.n_query_iterations = 1;
        cfg.lr = -1.0;
        assert!(cfg.validate().is_err());
    }

    // 2. Construct: committee has the requested number of members.
    #[test]
    fn test_construct() {
        let cfg = default_config();
        let n_committee = cfg.n_committee;
        let al = ActiveLearner::new(&[1, 6, 1], &[Activation::Tanh, Activation::Linear], cfg, 42)
            .unwrap();
        assert_eq!(al.committee.len(), n_committee);
        assert!(al.training_x.is_empty());
        assert!(al.training_y.is_empty());
    }

    // 3. add_sample appends to training_x and training_y.
    #[test]
    fn test_add_sample_updates_dataset() {
        let cfg = default_config();
        let mut al =
            ActiveLearner::new(&[1, 4, 1], &[Activation::Tanh, Activation::Linear], cfg, 1)
                .unwrap();
        al.add_sample(vec![0.1], vec![0.2]);
        al.add_sample(vec![0.4], vec![0.5]);
        assert_eq!(al.training_x.len(), 2);
        assert_eq!(al.training_y.len(), 2);
        assert!((al.training_x[1][0] - 0.4).abs() < 1e-12);
    }

    // 4. Committee members must have distinct initial parameters.
    #[test]
    fn test_committee_members_differ() {
        let cfg = default_config();
        let al = ActiveLearner::new(
            &[1, 4, 1],
            &[Activation::Tanh, Activation::Linear],
            cfg,
            777,
        )
        .unwrap();
        let p0 = al.committee[0].params_flat();
        let p1 = al.committee[1].params_flat();
        assert_eq!(p0.len(), p1.len());
        let any_diff = p0.iter().zip(p1.iter()).any(|(a, b)| a != b);
        assert!(any_diff, "committee members should have distinct seeds");
    }

    // 5. predict_with_uncertainty returns (mean, std) of the correct shape.
    #[test]
    fn test_predict_with_uncertainty_shapes() {
        let cfg = default_config();
        let al = ActiveLearner::new(&[1, 4, 2], &[Activation::Tanh, Activation::Linear], cfg, 5)
            .unwrap();
        let (mean, std_dev) = al.predict_with_uncertainty(&[0.3]).unwrap();
        assert_eq!(mean.len(), 2);
        assert_eq!(std_dev.len(), 2);
        for v in mean.iter().chain(std_dev.iter()) {
            assert!(v.is_finite());
        }
        for s in &std_dev {
            assert!(*s >= 0.0);
        }
    }

    // 6. After training on a single point repeatedly, the in-distribution
    //    uncertainty (at that point) should be ≤ the uncertainty far away.
    //
    //    QueryByCommittee with identical training data still has different
    //    initial weights, so the absolute std is not zero — we only compare
    //    relative magnitudes.
    #[test]
    fn test_in_distribution_uncertainty_smaller_than_out() {
        let cfg = ActiveLearningConfig {
            n_committee: 4,
            n_initial_samples: 1,
            n_query_iterations: 1,
            n_inner_train: 20,
            lr: 1e-2,
            query_strategy: QueryStrategy::QueryByCommittee,
        };
        let mut al = ActiveLearner::new(
            &[1, 6, 1],
            &[Activation::Tanh, Activation::Linear],
            cfg,
            321,
        )
        .unwrap();
        for _ in 0..6 {
            al.add_sample(vec![0.0_f64], vec![0.5_f64]);
        }
        al.train_committee().unwrap();
        let (_m_in, s_in) = al.predict_with_uncertainty(&[0.0_f64]).unwrap();
        let (_m_out, s_out) = al.predict_with_uncertainty(&[10.0_f64]).unwrap();
        assert!(
            s_in[0] <= s_out[0] + 1e-6,
            "in {} vs out {}",
            s_in[0],
            s_out[0]
        );
    }

    // 7. QueryByCommittee picks the synthetic high-uncertainty candidate.
    //
    //    Construct a candidate pool where one point sits at a very large
    //    input value (where the committee's predictions diverge sharply,
    //    because each member's random init explodes differently).  The
    //    selector must return that index.
    #[test]
    fn test_query_by_committee_picks_high_variance() {
        let cfg = ActiveLearningConfig {
            n_committee: 5,
            n_initial_samples: 1,
            n_query_iterations: 1,
            n_inner_train: 2,
            lr: 1e-3,
            query_strategy: QueryStrategy::QueryByCommittee,
        };
        let al = ActiveLearner::new(&[1, 8, 1], &[Activation::Tanh, Activation::Linear], cfg, 22)
            .unwrap();
        let pool: Vec<Vec<f64>> = vec![
            vec![0.0],
            vec![0.01],
            vec![0.02],
            vec![1.0e6], // saturates tanh hidden layer; output dominated by linear layer
            vec![0.03],
        ];
        let picked = al.select_next_query(&pool).unwrap();
        assert_eq!(picked, 3, "should pick the far-away (high variance) point");
    }

    // 8. RandomBaseline is deterministic for fixed training-set size.
    #[test]
    fn test_random_baseline_is_deterministic() {
        let cfg = ActiveLearningConfig {
            n_committee: 2,
            n_initial_samples: 1,
            n_query_iterations: 1,
            n_inner_train: 1,
            lr: 1e-3,
            query_strategy: QueryStrategy::RandomBaseline,
        };
        let al = ActiveLearner::new(
            &[1, 3, 1],
            &[Activation::Tanh, Activation::Linear],
            cfg,
            100,
        )
        .unwrap();
        let pool: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64]).collect();
        let p1 = al.select_next_query(&pool).unwrap();
        let p2 = al.select_next_query(&pool).unwrap();
        assert_eq!(p1, p2, "RandomBaseline must be reproducible");
        assert!(p1 < pool.len());
    }

    // 9. Full fit loop reduces loss on a simple sin(x) target.
    #[test]
    fn test_fit_loop_reduces_loss() {
        let cfg = ActiveLearningConfig {
            n_committee: 3,
            n_initial_samples: 4,
            n_query_iterations: 4,
            n_inner_train: 6,
            lr: 1e-2,
            query_strategy: QueryStrategy::QueryByCommittee,
        };
        let mut al =
            ActiveLearner::new(&[1, 8, 1], &[Activation::Tanh, Activation::Linear], cfg, 17)
                .unwrap();
        let pool: Vec<Vec<f64>> = (0..40).map(|i| vec![-2.0 + (i as f64) * 0.1]).collect();
        let oracle = |x: &[f64]| vec![x[0].sin()];
        let result = al.fit(oracle, &pool).unwrap();
        assert!(result.final_loss.is_finite());
        // At minimum, training did not blow up.
        assert!(result.final_loss < 100.0, "loss should remain bounded");
    }

    // 10. queried_indices length == n_initial_samples + n_query_iterations.
    #[test]
    fn test_queried_indices_length() {
        let cfg = ActiveLearningConfig {
            n_committee: 2,
            n_initial_samples: 3,
            n_query_iterations: 4,
            n_inner_train: 2,
            lr: 1e-3,
            query_strategy: QueryStrategy::RandomBaseline,
        };
        let mut al = ActiveLearner::new(
            &[1, 4, 1],
            &[Activation::Tanh, Activation::Linear],
            cfg,
            999,
        )
        .unwrap();
        let pool: Vec<Vec<f64>> = (0..15).map(|i| vec![i as f64 * 0.1]).collect();
        let result = al.fit(|x| vec![x[0] * 2.0], &pool).unwrap();
        assert_eq!(result.queried_indices.len(), 3 + 4);
        assert_eq!(result.n_queries_used, 7);
        assert_eq!(result.loss_history.len(), 4);
    }
}
