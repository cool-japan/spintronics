//! Bayesian optimization for materials parameter search.
//!
//! Uses a Gaussian Process (GP) surrogate with an RBF kernel to model an
//! expensive black-box function `f(x)` over a bounded parameter space.  Three
//! acquisition functions are provided to pick the next query point:
//!
//! - **Expected Improvement (EI)** balances exploration (high posterior
//!   variance) vs exploitation (predicted improvement over the current best).
//! - **Upper Confidence Bound (UCB)** picks `μ(x) + κ·σ(x)` with an
//!   exploration weight `κ`.
//! - **Posterior Variance** picks the most uncertain point — pure exploration.
//!
//! ## Use cases for spintronics materials discovery
//!
//! - Find the `(J, K, D)` combination that maximises a target (skyrmion
//!   radius, switching speed, anomalous Hall conductivity).
//! - Locate the phase transition between FM and spiral ground states.
//! - Sample-efficient hyper-parameter tuning for the v0.7.0
//!   [`crate::autodiff::neural::Mlp`] and v0.8.0
//!   [`crate::autodiff::equivariant::EquivariantMlp`] surrogates.
//!
//! ## Algorithm
//!
//! 1. **Kernel** — squared exponential (RBF):
//!    `k(x, x') = σ_f² · exp(−‖x − x'‖² / (2 ℓ²))`,
//!    optionally augmented by `σ_n² · δ_{xx'}` on the training diagonal.
//! 2. **Fit** — Cholesky factorise
//!    `K + (σ_n² + ε_jitter) I = L Lᵀ`,
//!    then solve `L Lᵀ α = y` (forward + back substitution) and cache
//!    `(L, α)`.
//! 3. **Predict** — posterior mean and variance at `x*`:
//!    `μ(x*) = k_*ᵀ α`,
//!    `σ²(x*) = k(x*, x*) − vᵀ v` where `L v = k_*`.
//! 4. **EI** — `EI(x*) = (μ − y_best) Φ(z) + σ φ(z)` with
//!    `z = (μ − y_best) / σ`, `Φ` the standard-normal CDF, `φ` the PDF.
//!    Implemented via Abramowitz–Stegun erf 7.1.26 (5-term polynomial,
//!    |error| < 1.5·10⁻⁷).
//! 5. **Selection** — grid search over a discrete candidate pool, picking the
//!    point with maximum acquisition.
//!
//! ## References
//!
//! - J. Mockus, *Bayesian Approach to Global Optimization* (1989).
//! - E. Brochu, V. M. Cora & N. de Freitas, "A Tutorial on Bayesian
//!   Optimization of Expensive Cost Functions, with Application to Active
//!   User Modeling and Hierarchical Reinforcement Learning",
//!   arXiv:1012.2599 (2010).
//! - P. I. Frazier, "A Tutorial on Bayesian Optimization", arXiv:1807.02811
//!   (2018).
//! - C. E. Rasmussen & C. K. I. Williams, *Gaussian Processes for Machine
//!   Learning*, MIT Press (2006).
//! - M. Abramowitz & I. A. Stegun, *Handbook of Mathematical Functions*,
//!   Dover (1965), §7.1.26.

use crate::error::{dimension_mismatch, invalid_param, numerical_error, Result};

// ─── Internal LCG for reproducible random samples ────────────────────────────

/// Same Knuth-style 64-bit LCG used elsewhere in `autodiff`, duplicated here
/// so this module stays self-contained.
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

    /// Uniform integer index in `0..n` (with negligible modulo bias for the
    /// candidate-pool sizes we encounter in BO).
    fn next_index(&mut self, n: usize) -> usize {
        if n == 0 {
            0
        } else {
            (self.next_u64() as usize) % n
        }
    }
}

// ─── erf / Φ / φ utilities ───────────────────────────────────────────────────

/// Abramowitz–Stegun 7.1.26 polynomial approximation of `erf(x)`.
///
/// Maximum absolute error on `|x| ≤ ∞` is `≤ 1.5·10⁻⁷`, which is more than
/// adequate for picking acquisition maxima.
fn erf_approx(x: f64) -> f64 {
    // Constants (AS 7.1.26).
    const A1: f64 = 0.254_829_592;
    const A2: f64 = -0.284_496_736;
    const A3: f64 = 1.421_413_741;
    const A4: f64 = -1.453_152_027;
    const A5: f64 = 1.061_405_429;
    const P: f64 = 0.327_591_1;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let ax = x.abs();
    let t = 1.0 / (1.0 + P * ax);
    let poly = (((A5 * t + A4) * t + A3) * t + A2) * t + A1;
    let y = 1.0 - poly * t * (-(ax * ax)).exp();
    sign * y
}

/// Standard-normal CDF `Φ(z) = ½ (1 + erf(z / √2))`.
fn normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf_approx(z / std::f64::consts::SQRT_2))
}

/// Standard-normal PDF `φ(z) = exp(−z²/2) / √(2π)`.
fn normal_pdf(z: f64) -> f64 {
    let inv_sqrt_two_pi = 1.0 / (2.0 * std::f64::consts::PI).sqrt();
    inv_sqrt_two_pi * (-0.5 * z * z).exp()
}

// ─── GpConfig ────────────────────────────────────────────────────────────────

/// Hyperparameters of the RBF Gaussian Process surrogate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpConfig {
    /// RBF length scale `ℓ`.  Larger values produce smoother surrogates.
    pub length_scale: f64,
    /// Signal variance `σ_f²`.
    pub signal_variance: f64,
    /// Noise variance `σ_n²` added to the diagonal of `K`.
    pub noise_variance: f64,
    /// Additional numerical jitter on the diagonal of `K` for Cholesky
    /// stability (`ε ≈ 10⁻⁶` works well for most problems).
    pub jitter: f64,
}

impl GpConfig {
    /// Validate that every entry is finite and positive (with the noise and
    /// jitter allowed to be zero only when summed they remain ≥ 0).
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when any field is
    /// not strictly positive (length / signal) or is negative (noise / jitter).
    pub fn validate(&self) -> Result<()> {
        if !(self.length_scale > 0.0 && self.length_scale.is_finite()) {
            return Err(invalid_param("length_scale", "must be positive and finite"));
        }
        if !(self.signal_variance > 0.0 && self.signal_variance.is_finite()) {
            return Err(invalid_param(
                "signal_variance",
                "must be positive and finite",
            ));
        }
        if !(self.noise_variance >= 0.0 && self.noise_variance.is_finite()) {
            return Err(invalid_param(
                "noise_variance",
                "must be non-negative and finite",
            ));
        }
        if !(self.jitter >= 0.0 && self.jitter.is_finite()) {
            return Err(invalid_param("jitter", "must be non-negative and finite"));
        }
        Ok(())
    }
}

// ─── GaussianProcess ─────────────────────────────────────────────────────────

/// Gaussian Process regressor over `ℝᵈ` with an RBF kernel.
///
/// After [`Self::fit`] the structure caches
///   - the lower-triangular Cholesky factor `L` of `K + σ_n² I + ε I`
///     (column-major, `n × n` flattened);
///   - `α = K⁻¹ y` (used directly in [`Self::predict_mean`]).
///
/// These caches survive across calls so that prediction is `O(n²)` after a
/// one-off `O(n³)` fit.
pub struct GaussianProcess {
    /// Kernel hyperparameters.
    pub config: GpConfig,
    /// Training inputs (one row per sample).
    pub x_train: Vec<Vec<f64>>,
    /// Training targets.
    pub y_train: Vec<f64>,
    /// `K⁻¹ y`, computed once at [`Self::fit`] time.
    alpha: Vec<f64>,
    /// Lower-triangular Cholesky factor (row-major, `n × n` flat).
    l_factor: Vec<f64>,
    /// Number of training rows (cached for convenience).
    n_train: usize,
}

impl GaussianProcess {
    /// Construct an *unfit* GP with the supplied hyperparameters.  Call
    /// [`Self::fit`] before any prediction.
    ///
    /// # Errors
    /// Propagates [`GpConfig::validate`].
    pub fn new(config: GpConfig) -> Result<Self> {
        config.validate()?;
        Ok(Self {
            config,
            x_train: Vec::new(),
            y_train: Vec::new(),
            alpha: Vec::new(),
            l_factor: Vec::new(),
            n_train: 0,
        })
    }

    /// Reasonable default hyperparameters for unit-scale inputs and outputs.
    pub fn config_default() -> GpConfig {
        GpConfig {
            length_scale: 1.0,
            signal_variance: 1.0,
            noise_variance: 1.0e-4,
            jitter: 1.0e-6,
        }
    }

    /// RBF kernel `k(a, b) = σ_f² · exp(−‖a − b‖² / (2 ℓ²))`.
    fn rbf_kernel(&self, a: &[f64], b: &[f64]) -> f64 {
        let mut sq = 0.0_f64;
        for (x, y) in a.iter().zip(b.iter()) {
            let d = x - y;
            sq += d * d;
        }
        let inv_2l2 = 1.0 / (2.0 * self.config.length_scale * self.config.length_scale);
        self.config.signal_variance * (-sq * inv_2l2).exp()
    }

    /// Fit the GP to a training set by Cholesky-factorising
    /// `K + (σ_n² + ε) I` and caching `K⁻¹ y`.
    ///
    /// If the Cholesky decomposition fails after one jitter retry (jitter
    /// scaled up by 10²), the method returns a
    /// [`crate::error::Error::NumericalError`].
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when `x_train` and
    /// `y_train` disagree in length, or a numerical error on Cholesky
    /// breakdown.
    pub fn fit(&mut self, x_train: Vec<Vec<f64>>, y_train: Vec<f64>) -> Result<()> {
        if x_train.len() != y_train.len() {
            return Err(dimension_mismatch(
                &format!("{} x rows", x_train.len()),
                &format!("{} y rows", y_train.len()),
            ));
        }
        if x_train.is_empty() {
            return Err(invalid_param("x_train", "must contain ≥ 1 training row"));
        }
        let n = x_train.len();
        // Sanity-check input width consistency.
        let d = x_train[0].len();
        for (i, row) in x_train.iter().enumerate() {
            if row.len() != d {
                return Err(dimension_mismatch(
                    &format!("{d} input dims (row 0)"),
                    &format!("{} input dims (row {i})", row.len()),
                ));
            }
        }
        self.x_train = x_train;
        self.y_train = y_train;
        self.n_train = n;

        // Build K + noise + jitter.
        let mut k_matrix = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                k_matrix[i * n + j] = self.rbf_kernel(&self.x_train[i], &self.x_train[j]);
            }
            k_matrix[i * n + i] += self.config.noise_variance + self.config.jitter;
        }

        // Attempt Cholesky; if it fails, bump the jitter once and retry.
        let l = match cholesky_lower(&k_matrix, n) {
            Some(l) => l,
            None => {
                let extra = self.config.jitter.max(1.0e-8) * 100.0;
                let mut k_retry = k_matrix.clone();
                for i in 0..n {
                    k_retry[i * n + i] += extra;
                }
                cholesky_lower(&k_retry, n).ok_or_else(|| {
                    numerical_error(
                        "GP kernel matrix is not positive definite even after jitter retry",
                    )
                })?
            },
        };
        self.l_factor = l;

        // Solve K α = y  ↔  L Lᵀ α = y.
        let z = forward_substitution(&self.l_factor, &self.y_train, n);
        self.alpha = backward_substitution_transposed(&self.l_factor, &z, n);
        Ok(())
    }

    /// Posterior mean at a single query point.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] if the GP has not
    /// been fitted, or [`crate::error::Error::DimensionMismatch`] if `x`
    /// width disagrees with the training inputs.
    pub fn predict_mean(&self, x: &[f64]) -> Result<f64> {
        if self.n_train == 0 {
            return Err(invalid_param("GaussianProcess", "must call fit() first"));
        }
        if x.len() != self.x_train[0].len() {
            return Err(dimension_mismatch(
                &format!("{} input dims", self.x_train[0].len()),
                &format!("{} input dims", x.len()),
            ));
        }
        let mut mu = 0.0_f64;
        for i in 0..self.n_train {
            mu += self.rbf_kernel(x, &self.x_train[i]) * self.alpha[i];
        }
        Ok(mu)
    }

    /// Posterior `(mean, standard_deviation)` at a single query point.
    ///
    /// Variance is clamped to zero from below for numerical safety (small
    /// negative values can appear from round-off in `kᵀ K⁻¹ k`).
    ///
    /// # Errors
    /// Returns the same errors as [`Self::predict_mean`].
    pub fn predict(&self, x: &[f64]) -> Result<(f64, f64)> {
        if self.n_train == 0 {
            return Err(invalid_param("GaussianProcess", "must call fit() first"));
        }
        if x.len() != self.x_train[0].len() {
            return Err(dimension_mismatch(
                &format!("{} input dims", self.x_train[0].len()),
                &format!("{} input dims", x.len()),
            ));
        }
        // k_star[i] = k(x, x_train[i]).
        let mut k_star = vec![0.0_f64; self.n_train];
        for (i, k_entry) in k_star.iter_mut().enumerate() {
            *k_entry = self.rbf_kernel(x, &self.x_train[i]);
        }
        // Mean.
        let mut mu = 0.0_f64;
        for (k, a) in k_star.iter().zip(self.alpha.iter()) {
            mu += k * a;
        }
        // Variance: k(x,x) − vᵀ v, with L v = k_star.
        let v = forward_substitution(&self.l_factor, &k_star, self.n_train);
        let mut vv = 0.0_f64;
        for vi in &v {
            vv += vi * vi;
        }
        let kxx = self.rbf_kernel(x, x);
        let var = (kxx - vv).max(0.0);
        Ok((mu, var.sqrt()))
    }
}

// ─── Cholesky / triangular solves ────────────────────────────────────────────

/// In-place lower-triangular Cholesky decomposition `A = L Lᵀ` of an `n × n`
/// symmetric positive-definite matrix stored row-major.  Returns `None` when
/// a non-positive pivot is encountered (i.e. the matrix is not positive
/// definite up to the supplied jitter).
fn cholesky_lower(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if !(sum > 0.0 && sum.is_finite()) {
                    return None;
                }
                l[i * n + j] = sum.sqrt();
            } else {
                let diag = l[j * n + j];
                if !(diag.abs() > 0.0 && diag.is_finite()) {
                    return None;
                }
                l[i * n + j] = sum / diag;
            }
        }
    }
    Some(l)
}

/// Solve `L y = b` for a lower-triangular matrix `L` stored row-major.
fn forward_substitution(l: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut y = vec![0.0_f64; n];
    for i in 0..n {
        let mut sum = b[i];
        for j in 0..i {
            sum -= l[i * n + j] * y[j];
        }
        y[i] = sum / l[i * n + i];
    }
    y
}

/// Solve `Lᵀ x = y` for a lower-triangular matrix `L` stored row-major.
fn backward_substitution_transposed(l: &[f64], y: &[f64], n: usize) -> Vec<f64> {
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = y[i];
        for j in (i + 1)..n {
            sum -= l[j * n + i] * x[j];
        }
        x[i] = sum / l[i * n + i];
    }
    x
}

// ─── AcquisitionStrategy ─────────────────────────────────────────────────────

/// Which acquisition function drives candidate selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AcquisitionStrategy {
    /// Expected Improvement (Mockus 1989).
    ExpectedImprovement,
    /// Upper Confidence Bound `μ(x) + κ·σ(x)`.  `kappa` is rounded to the
    /// nearest integer because it doubles as the variant payload.
    UpperConfidenceBound {
        /// Exploration weight `κ`.
        kappa: i32,
    },
    /// Pure exploration: pick the point with maximum posterior variance.
    PosteriorVariance,
}

// ─── BayesianOptConfig ───────────────────────────────────────────────────────

/// Hyperparameters that drive a [`BayesianOptimizer`] loop.
#[derive(Debug, Clone)]
pub struct BayesianOptConfig {
    /// How many random samples are drawn from the candidate pool before the
    /// first acquisition step.
    pub n_initial_samples: usize,
    /// How many `(propose, evaluate, refit)` iterations to run after seeding.
    pub n_iterations: usize,
    /// Hyperparameters of the GP surrogate.
    pub gp_config: GpConfig,
    /// Acquisition strategy.
    pub acquisition: AcquisitionStrategy,
    /// If `true`, optimisation maximises `f`; otherwise it minimises `f`.
    pub maximize: bool,
}

impl BayesianOptConfig {
    /// Validate that the counts are non-zero and the GP config is sane.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when
    /// `n_initial_samples == 0` or `n_iterations == 0`, or propagates
    /// [`GpConfig::validate`].
    pub fn validate(&self) -> Result<()> {
        if self.n_initial_samples == 0 {
            return Err(invalid_param(
                "n_initial_samples",
                "need ≥ 1 sample to seed the GP",
            ));
        }
        if self.n_iterations == 0 {
            return Err(invalid_param(
                "n_iterations",
                "must request ≥ 1 acquisition iteration",
            ));
        }
        self.gp_config.validate()?;
        Ok(())
    }
}

// ─── BayesianOptResult ───────────────────────────────────────────────────────

/// Final report of a [`BayesianOptimizer::optimize`] run.
#[derive(Debug, Clone)]
pub struct BayesianOptResult {
    /// Best `x` found (the argmax / argmin of `y` over `history_y`).
    pub best_x: Vec<f64>,
    /// Best `y` value found.
    pub best_y: f64,
    /// Total number of oracle calls.
    pub n_evaluations: usize,
    /// Every queried `x` in evaluation order.
    pub history_x: Vec<Vec<f64>>,
    /// Every observed `y` in evaluation order.
    pub history_y: Vec<f64>,
}

// ─── BayesianOptimizer ───────────────────────────────────────────────────────

/// Stateful Bayesian optimiser that owns the GP surrogate and the cumulative
/// evaluation history.
pub struct BayesianOptimizer {
    /// Configuration for this BO run.
    pub config: BayesianOptConfig,
    /// GP surrogate fit to the history.
    pub gp: GaussianProcess,
    /// Every queried `x` (one row per evaluation, in query order).
    pub history_x: Vec<Vec<f64>>,
    /// Every observed `y` (one entry per evaluation).
    pub history_y: Vec<f64>,
}

impl BayesianOptimizer {
    /// Build a fresh optimiser with no history.
    ///
    /// # Errors
    /// Propagates [`BayesianOptConfig::validate`].
    pub fn new(config: BayesianOptConfig) -> Result<Self> {
        config.validate()?;
        let gp = GaussianProcess::new(config.gp_config)?;
        Ok(Self {
            config,
            gp,
            history_x: Vec::new(),
            history_y: Vec::new(),
        })
    }

    /// Append a `(x, y)` evaluation to the history (used to seed the GP).
    pub fn add_sample(&mut self, x: Vec<f64>, y: f64) {
        self.history_x.push(x);
        self.history_y.push(y);
    }

    /// (Re-)fit the GP surrogate on the current evaluation history.
    ///
    /// # Errors
    /// Propagates errors from [`GaussianProcess::fit`], in particular
    /// numerical breakdown of the Cholesky factorisation.
    pub fn fit_gp(&mut self) -> Result<()> {
        // For maximisation we *can* fit on the raw values directly; for
        // minimisation EI is computed against y_best = min(history_y).  We
        // therefore fit on raw values regardless of `maximize`.
        self.gp.fit(self.history_x.clone(), self.history_y.clone())
    }

    /// Compute the acquisition score at a candidate point.  Higher score =
    /// more promising next query for every strategy variant.
    ///
    /// For Expected Improvement, the standard form (Mockus) is used.  When
    /// `maximize == true` improvement means `μ(x) − y_best`; when `false`
    /// the role of `y_best` is flipped and we maximise `y_best − μ(x)`.  We
    /// implement the latter by negating `mu` and the best value.
    ///
    /// # Errors
    /// Propagates errors from [`GaussianProcess::predict`].
    pub fn acquisition_score(&self, x: &[f64]) -> Result<f64> {
        let (mu, sigma) = self.gp.predict(x)?;
        match self.config.acquisition {
            AcquisitionStrategy::ExpectedImprovement => {
                let best_raw = self.current_best_y_for_score();
                // Effective sign: positive direction = "better".
                let (improvement_centre, effective_best) = if self.config.maximize {
                    (mu, best_raw)
                } else {
                    (-mu, -best_raw)
                };
                let diff = improvement_centre - effective_best;
                if sigma <= 0.0 {
                    return Ok(diff.max(0.0));
                }
                let z = diff / sigma;
                let ei = diff * normal_cdf(z) + sigma * normal_pdf(z);
                Ok(ei)
            },
            AcquisitionStrategy::UpperConfidenceBound { kappa } => {
                let k = kappa as f64;
                if self.config.maximize {
                    Ok(mu + k * sigma)
                } else {
                    Ok(-mu + k * sigma)
                }
            },
            AcquisitionStrategy::PosteriorVariance => Ok(sigma),
        }
    }

    /// Pick the index of the candidate with the maximum acquisition score.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when
    /// `candidate_pool` is empty, or propagates errors from
    /// [`Self::acquisition_score`].
    pub fn select_next_query(&self, candidate_pool: &[Vec<f64>]) -> Result<usize> {
        if candidate_pool.is_empty() {
            return Err(invalid_param(
                "candidate_pool",
                "must contain ≥ 1 candidate",
            ));
        }
        let mut best_idx = 0_usize;
        let mut best_score = f64::NEG_INFINITY;
        for (idx, x) in candidate_pool.iter().enumerate() {
            let s = self.acquisition_score(x)?;
            if s > best_score {
                best_score = s;
                best_idx = idx;
            }
        }
        Ok(best_idx)
    }

    /// Full BO loop: seed with `n_initial_samples` random draws from the
    /// pool, then iterate `n_iterations` rounds of `(propose, evaluate,
    /// refit)`.  The candidate pool is *discrete*; sampling without
    /// replacement is enforced so the same `x` is never queried twice.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when the pool is
    /// too small to fulfil the requested number of iterations, or propagates
    /// numerical errors from the GP fit.
    pub fn optimize<F>(
        &mut self,
        oracle: F,
        candidate_pool: &[Vec<f64>],
        rng_seed: u64,
    ) -> Result<BayesianOptResult>
    where
        F: Fn(&[f64]) -> f64,
    {
        let total_budget = self.config.n_initial_samples + self.config.n_iterations;
        if candidate_pool.len() < total_budget {
            return Err(invalid_param(
                "candidate_pool",
                &format!(
                    "needs ≥ {} candidates (got {})",
                    total_budget,
                    candidate_pool.len(),
                ),
            ));
        }

        let mut rng = Lcg::new(rng_seed);
        let mut visited = vec![false; candidate_pool.len()];

        // Seeding: random draws without replacement.
        for _ in 0..self.config.n_initial_samples {
            let mut pick = rng.next_index(candidate_pool.len());
            // Linear probe forward to find an unvisited slot.
            let mut probes = 0;
            while visited[pick] && probes < candidate_pool.len() {
                pick = (pick + 1) % candidate_pool.len();
                probes += 1;
            }
            visited[pick] = true;
            let x = candidate_pool[pick].clone();
            let y = oracle(&x);
            self.add_sample(x, y);
        }
        self.fit_gp()?;

        // BO loop.
        for _ in 0..self.config.n_iterations {
            // Build the *remaining* pool view; the index returned is into the
            // sub-pool, so we map it back to the global index for the
            // visited bitmap.
            let mut remaining: Vec<usize> = Vec::with_capacity(candidate_pool.len());
            for (i, &v) in visited.iter().enumerate() {
                if !v {
                    remaining.push(i);
                }
            }
            if remaining.is_empty() {
                break;
            }
            let sub_pool: Vec<Vec<f64>> = remaining
                .iter()
                .map(|&i| candidate_pool[i].clone())
                .collect();
            let local = self.select_next_query(&sub_pool)?;
            let global = remaining[local];
            visited[global] = true;
            let x = candidate_pool[global].clone();
            let y = oracle(&x);
            self.add_sample(x, y);
            self.fit_gp()?;
        }

        let (best_idx, best_y) = self.argbest();
        Ok(BayesianOptResult {
            best_x: self.history_x[best_idx].clone(),
            best_y,
            n_evaluations: self.history_y.len(),
            history_x: self.history_x.clone(),
            history_y: self.history_y.clone(),
        })
    }

    /// Internal helper: current best (index, value) in `history_y` under the
    /// chosen optimisation direction.
    fn argbest(&self) -> (usize, f64) {
        let mut best_idx = 0_usize;
        let mut best_val = if self.config.maximize {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
        for (i, &y) in self.history_y.iter().enumerate() {
            let better = if self.config.maximize {
                y > best_val
            } else {
                y < best_val
            };
            if better {
                best_val = y;
                best_idx = i;
            }
        }
        (best_idx, best_val)
    }

    /// Current best raw `y` value (without sign-flipping for direction).
    fn current_best_y_for_score(&self) -> f64 {
        let (_, best) = self.argbest();
        best
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // 1. GpConfig validation rejects bad inputs.
    #[test]
    fn test_gp_config_validation() {
        let good = GpConfig {
            length_scale: 1.0,
            signal_variance: 0.5,
            noise_variance: 1.0e-4,
            jitter: 1.0e-6,
        };
        assert!(good.validate().is_ok());
        let bad_ls = GpConfig {
            length_scale: 0.0,
            ..good
        };
        assert!(bad_ls.validate().is_err());
        let bad_sig = GpConfig {
            signal_variance: -0.1,
            ..good
        };
        assert!(bad_sig.validate().is_err());
        let bad_noise = GpConfig {
            noise_variance: -1.0,
            ..good
        };
        assert!(bad_noise.validate().is_err());
        let bad_jitter = GpConfig {
            jitter: f64::INFINITY,
            ..good
        };
        assert!(bad_jitter.validate().is_err());
        // GaussianProcess construction propagates the validation.
        assert!(GaussianProcess::new(bad_ls).is_err());
    }

    // 2. Fit on 5 points succeeds.
    #[test]
    fn test_gp_fit_5_points_succeeds() {
        let mut gp = GaussianProcess::new(GaussianProcess::config_default()).unwrap();
        let xs: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64 * 0.25]).collect();
        let ys: Vec<f64> = xs.iter().map(|x| x[0].sin()).collect();
        gp.fit(xs.clone(), ys.clone()).unwrap();
        assert_eq!(gp.n_train, 5);
        assert_eq!(gp.alpha.len(), 5);
        assert_eq!(gp.l_factor.len(), 25);
    }

    // 3. predict_mean at a training point recovers the training value within
    //    the noise floor.
    #[test]
    fn test_predict_at_training_point() {
        let mut gp = GaussianProcess::new(GpConfig {
            length_scale: 0.5,
            signal_variance: 1.0,
            noise_variance: 1.0e-8,
            jitter: 1.0e-10,
        })
        .unwrap();
        let xs: Vec<Vec<f64>> = (0..4)
            .map(|i| vec![(i as f64) * 0.3, (i as f64) * 0.1])
            .collect();
        let ys: Vec<f64> = vec![0.1, 0.4, -0.2, 0.9];
        gp.fit(xs.clone(), ys.clone()).unwrap();
        for (x, &y) in xs.iter().zip(ys.iter()) {
            let mu = gp.predict_mean(x).unwrap();
            assert!(approx(mu, y, 1.0e-4), "predict {mu} vs train {y}");
        }
    }

    // 4. Posterior variance is strictly positive at an unseen point.
    #[test]
    fn test_predict_variance_positive_at_unseen() {
        let mut gp = GaussianProcess::new(GaussianProcess::config_default()).unwrap();
        let xs: Vec<Vec<f64>> = vec![vec![0.0], vec![0.5], vec![1.0]];
        let ys: Vec<f64> = vec![0.0, 0.5, 1.0];
        gp.fit(xs, ys).unwrap();
        let (_mu, sigma) = gp.predict(&[2.0]).unwrap();
        assert!(sigma > 0.0, "expected positive sigma at unseen point");
    }

    // 5. RBF kernel symmetric: k(a, b) == k(b, a).
    #[test]
    fn test_rbf_kernel_symmetric() {
        let gp = GaussianProcess::new(GaussianProcess::config_default()).unwrap();
        let a = vec![0.3_f64, -0.2, 0.7];
        let b = vec![1.1_f64, 0.5, -0.4];
        let kab = gp.rbf_kernel(&a, &b);
        let kba = gp.rbf_kernel(&b, &a);
        assert!(approx(kab, kba, 1.0e-15));
    }

    // 6. RBF kernel attains its maximum at distance 0 (equal to σ_f²).
    #[test]
    fn test_rbf_kernel_max_at_zero_distance() {
        let cfg = GpConfig {
            length_scale: 0.8,
            signal_variance: 2.5,
            noise_variance: 0.0,
            jitter: 0.0,
        };
        let gp = GaussianProcess::new(cfg).unwrap();
        let a = vec![1.2_f64, -0.3];
        let kaa = gp.rbf_kernel(&a, &a);
        let kab = gp.rbf_kernel(&a, &[1.5_f64, 0.0]);
        assert!(approx(kaa, cfg.signal_variance, 1.0e-12));
        assert!(kaa > kab);
    }

    // 7. AcquisitionStrategy variants construct correctly.
    #[test]
    fn test_acquisition_strategy_variants() {
        let _ei = AcquisitionStrategy::ExpectedImprovement;
        let _ucb = AcquisitionStrategy::UpperConfidenceBound { kappa: 2 };
        let _pv = AcquisitionStrategy::PosteriorVariance;
        // Distinct discriminants.
        assert_ne!(
            AcquisitionStrategy::ExpectedImprovement,
            AcquisitionStrategy::PosteriorVariance,
        );
        assert_ne!(
            AcquisitionStrategy::UpperConfidenceBound { kappa: 1 },
            AcquisitionStrategy::UpperConfidenceBound { kappa: 2 },
        );
    }

    // 8. EI = 0 when sigma -> 0 at a point that cannot improve over the best.
    #[test]
    fn test_ei_zero_when_no_room_for_improvement() {
        let cfg = BayesianOptConfig {
            n_initial_samples: 1,
            n_iterations: 1,
            gp_config: GpConfig {
                length_scale: 0.5,
                signal_variance: 1.0,
                noise_variance: 1.0e-10,
                jitter: 1.0e-10,
            },
            acquisition: AcquisitionStrategy::ExpectedImprovement,
            maximize: true,
        };
        let mut opt = BayesianOptimizer::new(cfg).unwrap();
        // Single training point at x=0 with y=1.0 — this is the best seen.
        opt.add_sample(vec![0.0_f64], 1.0);
        opt.fit_gp().unwrap();
        // At the training point, sigma ≈ 0 and improvement ≈ 0.
        let ei = opt.acquisition_score(&[0.0_f64]).unwrap();
        assert!(
            ei.abs() < 1.0e-3,
            "EI at the best point should ≈ 0, got {ei}"
        );
    }

    // 9. BayesianOptimizer construct + config validation.
    #[test]
    fn test_bayesian_optimizer_construct() {
        let bad = BayesianOptConfig {
            n_initial_samples: 0,
            n_iterations: 5,
            gp_config: GaussianProcess::config_default(),
            acquisition: AcquisitionStrategy::ExpectedImprovement,
            maximize: false,
        };
        assert!(BayesianOptimizer::new(bad).is_err());
        let good = BayesianOptConfig {
            n_initial_samples: 2,
            n_iterations: 4,
            gp_config: GaussianProcess::config_default(),
            acquisition: AcquisitionStrategy::ExpectedImprovement,
            maximize: false,
        };
        let opt = BayesianOptimizer::new(good).unwrap();
        assert_eq!(opt.history_x.len(), 0);
        assert_eq!(opt.history_y.len(), 0);
    }

    // 10. optimize() reduces best_y on a quadratic when minimising.
    #[test]
    fn test_optimize_minimisation_reduces_best() {
        // f(x) = (x − 0.7)² over a 21-point grid in [0, 1].
        let pool: Vec<Vec<f64>> = (0..21).map(|i| vec![(i as f64) / 20.0]).collect();
        let cfg = BayesianOptConfig {
            n_initial_samples: 3,
            n_iterations: 8,
            gp_config: GpConfig {
                length_scale: 0.15,
                signal_variance: 1.0,
                noise_variance: 1.0e-6,
                jitter: 1.0e-8,
            },
            acquisition: AcquisitionStrategy::ExpectedImprovement,
            maximize: false,
        };
        let mut opt = BayesianOptimizer::new(cfg).unwrap();
        let oracle = |x: &[f64]| (x[0] - 0.7).powi(2);
        let result = opt.optimize(oracle, &pool, 0x00C0_FFEE_u64).unwrap();
        // Closest pool point to x=0.7 is 0.7 itself (i = 14), giving f = 0.
        // BO should find something within 0.2 of that minimum.
        assert!(
            result.best_y < 0.05,
            "BO did not converge: best_y = {}",
            result.best_y,
        );
        // The initial samples alone are unlikely to land on 0.7 every run,
        // so result.history_y[0] should usually be > result.best_y.
        let initial_best = result.history_y[..3]
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        assert!(
            result.best_y <= initial_best,
            "best_y {} should be ≤ initial best {}",
            result.best_y,
            initial_best,
        );
    }

    // 11. maximize=true flag actually drives EI to seek larger y.
    #[test]
    fn test_optimize_maximisation() {
        // f(x) = −(x − 0.3)² is maximised at x = 0.3.
        let pool: Vec<Vec<f64>> = (0..21).map(|i| vec![(i as f64) / 20.0]).collect();
        let cfg = BayesianOptConfig {
            n_initial_samples: 3,
            n_iterations: 8,
            gp_config: GpConfig {
                length_scale: 0.15,
                signal_variance: 1.0,
                noise_variance: 1.0e-6,
                jitter: 1.0e-8,
            },
            acquisition: AcquisitionStrategy::ExpectedImprovement,
            maximize: true,
        };
        let mut opt = BayesianOptimizer::new(cfg).unwrap();
        let oracle = |x: &[f64]| -(x[0] - 0.3).powi(2);
        let result = opt.optimize(oracle, &pool, 0xBEEF_u64).unwrap();
        // Maximum is 0; BO should land within 0.05.
        assert!(
            result.best_y > -0.05,
            "BO did not converge to maximum: best_y = {}",
            result.best_y,
        );
    }

    // 12. Seed reproducibility: same seed → same optimisation trace.
    #[test]
    fn test_seed_reproducibility() {
        let pool: Vec<Vec<f64>> = (0..15).map(|i| vec![(i as f64) / 14.0]).collect();
        let make_opt = || {
            BayesianOptimizer::new(BayesianOptConfig {
                n_initial_samples: 3,
                n_iterations: 4,
                gp_config: GpConfig {
                    length_scale: 0.2,
                    signal_variance: 1.0,
                    noise_variance: 1.0e-6,
                    jitter: 1.0e-8,
                },
                acquisition: AcquisitionStrategy::UpperConfidenceBound { kappa: 2 },
                maximize: false,
            })
            .unwrap()
        };
        let oracle = |x: &[f64]| (x[0] - 0.4).powi(2) + 0.1;
        let mut a = make_opt();
        let mut b = make_opt();
        let ra = a.optimize(oracle, &pool, 0xABCD_u64).unwrap();
        let rb = b.optimize(oracle, &pool, 0xABCD_u64).unwrap();
        assert_eq!(ra.history_y.len(), rb.history_y.len());
        for (x, y) in ra.history_y.iter().zip(rb.history_y.iter()) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
        assert_eq!(ra.best_y.to_bits(), rb.best_y.to_bits());
    }
}
