//! Gradient-based parameter optimizers for physics model fitting.
//!
//! # Included Optimizers
//!
//! | Name | Description |
//! |------|-------------|
//! | [`Sgd`]   | Stochastic Gradient Descent with optional momentum |
//! | [`Adam`]  | Adaptive Moment Estimation (Kingma & Ba, 2014) |
//! | [`LBfgs`] | Limited-memory BFGS (Nocedal, 1980) |
//!
//! # High-level fitting
//!
//! [`ParameterFitter`] wraps any of the above optimizers and provides a
//! convenient `fit` method that accepts a differentiable loss closure built
//! with the [`Tape`](super::tape::Tape) API.
//!
//! # References
//!
//! - D. P. Kingma & J. Ba, "Adam: A Method for Stochastic Optimization",
//!   *ICLR* (2015), arXiv:1412.6980
//! - J. Nocedal, "Updating Quasi-Newton Matrices with Limited Storage",
//!   *Math. Comp.* **35**, 773 (1980)
//! - D. C. Liu & J. Nocedal, "On the limited-memory BFGS method for
//!   large-scale optimization", *Math. Prog.* **45**, 503 (1989)

use crate::error::{invalid_param, Error as SpinError};

// ─── Optimizer trait ──────────────────────────────────────────────────────────

/// Common interface shared by all gradient-based optimizers.
pub trait Optimizer {
    /// Apply one gradient-descent step: update `params` using `grads`.
    ///
    /// Both slices must have the same length.
    fn step(&mut self, params: &mut [f64], grads: &[f64]);

    /// Human-readable optimizer name.
    fn name(&self) -> &'static str;

    /// Reset internal optimizer state (moments, velocities, history).
    ///
    /// Call this between independent fitting runs.
    fn reset(&mut self);
}

// ─── SGD ─────────────────────────────────────────────────────────────────────

/// Stochastic Gradient Descent with optional Polyak–Ruppert momentum.
///
/// Update rule:
/// ```text
/// v_{t+1} = momentum · v_t + g_t
/// p_{t+1} = p_t − lr · v_{t+1}
/// ```
///
/// When `momentum = 0`, this reduces to vanilla gradient descent:
/// `p_{t+1} = p_t − lr · g_t`.
///
/// # Reference
/// - Polyak (1964); Ruppert (1988); Qian (1999)
pub struct Sgd {
    /// Learning rate (step size).
    pub lr: f64,
    /// Momentum coefficient ∈ [0, 1).
    pub momentum: f64,
    velocity: Vec<f64>,
}

impl Sgd {
    /// Create an SGD optimizer with the given learning rate and momentum.
    ///
    /// # Errors
    /// Returns [`SpinError::InvalidParameter`] if `lr ≤ 0` or
    /// `momentum` is not in `[0, 1)`.
    pub fn new(lr: f64, momentum: f64) -> Result<Self, SpinError> {
        if lr <= 0.0 {
            return Err(invalid_param("lr", "learning rate must be positive"));
        }
        if !(0.0..1.0).contains(&momentum) {
            return Err(invalid_param("momentum", "must be in [0, 1)"));
        }
        Ok(Self {
            lr,
            momentum,
            velocity: Vec::new(),
        })
    }

    /// Create a vanilla SGD optimizer (no momentum).
    ///
    /// # Errors
    /// Returns [`SpinError::InvalidParameter`] if `lr ≤ 0`.
    pub fn vanilla(lr: f64) -> Result<Self, SpinError> {
        Self::new(lr, 0.0)
    }

    fn ensure_velocity(&mut self, n: usize) {
        if self.velocity.len() != n {
            self.velocity = vec![0.0; n];
        }
    }
}

impl Optimizer for Sgd {
    fn step(&mut self, params: &mut [f64], grads: &[f64]) {
        debug_assert_eq!(params.len(), grads.len());
        self.ensure_velocity(params.len());
        let lr = self.lr;
        let mu = self.momentum;
        for ((p, g), v) in params
            .iter_mut()
            .zip(grads.iter())
            .zip(self.velocity.iter_mut())
        {
            *v = mu * (*v) + (*g);
            *p -= lr * (*v);
        }
    }

    fn name(&self) -> &'static str {
        "SGD"
    }

    fn reset(&mut self) {
        for v in &mut self.velocity {
            *v = 0.0;
        }
    }
}

// ─── Adam ────────────────────────────────────────────────────────────────────

/// Adam (Adaptive Moment Estimation) optimizer.
///
/// Update rule at time step `t`:
/// ```text
/// m_t = β₁ · m_{t-1} + (1 − β₁) · g_t        (biased first moment)
/// v_t = β₂ · v_{t-1} + (1 − β₂) · g_t²        (biased second moment)
/// m̂_t = m_t / (1 − β₁^t)                      (bias-corrected)
/// v̂_t = v_t / (1 − β₂^t)                      (bias-corrected)
/// p_t = p_{t-1} − α · m̂_t / (√v̂_t + ε)
/// ```
///
/// # Reference
/// Kingma & Ba, *ICLR* 2015, arXiv:1412.6980
pub struct Adam {
    /// Learning rate α (step size).
    pub lr: f64,
    /// First moment decay rate β₁ ∈ (0, 1).
    pub beta1: f64,
    /// Second moment decay rate β₂ ∈ (0, 1).
    pub beta2: f64,
    /// Numerical stability constant ε.
    pub eps: f64,
    /// First moment vector `m`.
    m: Vec<f64>,
    /// Second moment vector `v`.
    v: Vec<f64>,
    /// Time step counter.
    t: usize,
}

impl Adam {
    /// Create an Adam optimizer with explicit hyperparameters and `n_params`
    /// parameters (pre-allocates moment vectors).
    ///
    /// # Errors
    /// Returns [`SpinError::InvalidParameter`] for out-of-range hyperparameters.
    pub fn new(
        lr: f64,
        beta1: f64,
        beta2: f64,
        eps: f64,
        n_params: usize,
    ) -> Result<Self, SpinError> {
        if lr <= 0.0 {
            return Err(invalid_param("lr", "learning rate must be positive"));
        }
        if !(0.0..1.0).contains(&beta1) {
            return Err(invalid_param("beta1", "must be in (0, 1)"));
        }
        if !(0.0..1.0).contains(&beta2) {
            return Err(invalid_param("beta2", "must be in (0, 1)"));
        }
        if eps <= 0.0 {
            return Err(invalid_param("eps", "must be positive"));
        }
        Ok(Self {
            lr,
            beta1,
            beta2,
            eps,
            m: vec![0.0; n_params],
            v: vec![0.0; n_params],
            t: 0,
        })
    }

    /// Create Adam with the default hyperparameters proposed in the paper:
    /// α = 1e−3, β₁ = 0.9, β₂ = 0.999, ε = 1e−8.
    pub fn default_params(n_params: usize) -> Self {
        Self {
            lr: 1e-3,
            beta1: 0.9,
            beta2: 0.999,
            eps: 1e-8,
            m: vec![0.0; n_params],
            v: vec![0.0; n_params],
            t: 0,
        }
    }

    fn ensure_state(&mut self, n: usize) {
        if self.m.len() != n {
            self.m = vec![0.0; n];
            self.v = vec![0.0; n];
            self.t = 0;
        }
    }
}

impl Optimizer for Adam {
    fn step(&mut self, params: &mut [f64], grads: &[f64]) {
        debug_assert_eq!(params.len(), grads.len());
        self.ensure_state(params.len());
        self.t += 1;
        let t = self.t as f64;
        let (b1, b2, lr, eps) = (self.beta1, self.beta2, self.lr, self.eps);
        let bc1 = 1.0 - b1.powf(t);
        let bc2 = 1.0 - b2.powf(t);
        for (i, (p, g)) in params.iter_mut().zip(grads.iter()).enumerate() {
            self.m[i] = b1 * self.m[i] + (1.0 - b1) * (*g);
            self.v[i] = b2 * self.v[i] + (1.0 - b2) * (*g) * (*g);
            let m_hat = self.m[i] / bc1;
            let v_hat = self.v[i] / bc2;
            *p -= lr * m_hat / (v_hat.sqrt() + eps);
        }
    }

    fn name(&self) -> &'static str {
        "Adam"
    }

    fn reset(&mut self) {
        for m in &mut self.m {
            *m = 0.0;
        }
        for v in &mut self.v {
            *v = 0.0;
        }
        self.t = 0;
    }
}

// ─── L-BFGS ──────────────────────────────────────────────────────────────────

/// Limited-memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) optimizer.
///
/// Maintains a sliding window of `history_size` curvature pairs
/// `(s_k, y_k)` where `s_k = x_{k+1} − x_k` and `y_k = ∇f_{k+1} − ∇f_k`.
///
/// The search direction is computed via the two-loop recursion of Nocedal
/// (1980) applied to the current gradient.  A simple Armijo-style
/// scaling with fixed `lr` is used as the line search.
///
/// # References
/// - J. Nocedal, "Updating Quasi-Newton Matrices", *Math. Comp.* **35** (1980)
/// - D. C. Liu & J. Nocedal, *Math. Prog.* **45**, 503 (1989)
pub struct LBfgs {
    /// Step size multiplier applied to the Hessian-scaled direction.
    pub lr: f64,
    /// Maximum number of (s, y) curvature pairs to retain.
    pub history_size: usize,
    /// Stored `s_k = x_{k+1} − x_k` vectors.
    s_history: Vec<Vec<f64>>,
    /// Stored `y_k = ∇f_{k+1} − ∇f_k` vectors.
    y_history: Vec<Vec<f64>>,
    /// `ρ_k = 1 / (y_k^T s_k)`.
    rho_history: Vec<f64>,
    /// Parameter vector from the previous step.
    prev_params: Vec<f64>,
    /// Gradient from the previous step.
    prev_grads: Vec<f64>,
    /// Number of steps taken (used to decide when history is valid).
    step_count: usize,
}

impl LBfgs {
    /// Create an L-BFGS optimizer with the given step size and history window.
    ///
    /// # Errors
    /// Returns [`SpinError::InvalidParameter`] if `lr ≤ 0` or `history_size == 0`.
    pub fn new(lr: f64, history_size: usize) -> Result<Self, SpinError> {
        if lr <= 0.0 {
            return Err(invalid_param("lr", "learning rate must be positive"));
        }
        if history_size == 0 {
            return Err(invalid_param("history_size", "must be at least 1"));
        }
        Ok(Self {
            lr,
            history_size,
            s_history: Vec::new(),
            y_history: Vec::new(),
            rho_history: Vec::new(),
            prev_params: Vec::new(),
            prev_grads: Vec::new(),
            step_count: 0,
        })
    }

    /// Compute `a^T b` (dot product of two equal-length slices).
    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(ai, bi)| ai * bi).sum()
    }

    /// Two-loop recursion — compute H_k * q where q = current gradient.
    fn two_loop_direction(&self, grads: &[f64]) -> Vec<f64> {
        let m = self.s_history.len();
        let mut q: Vec<f64> = grads.to_vec();
        let mut alphas = vec![0.0_f64; m];

        // First loop: newest pairs first
        for i in (0..m).rev() {
            let alpha_i = self.rho_history[i] * Self::dot(&self.s_history[i], &q);
            for (qj, yij) in q.iter_mut().zip(self.y_history[i].iter()) {
                *qj -= alpha_i * yij;
            }
            alphas[i] = alpha_i;
        }

        // Initial Hessian scaling: H_0 = (s_{m-1}^T y_{m-1}) / (y_{m-1}^T y_{m-1})
        let h0 = if m > 0 {
            let yty = Self::dot(&self.y_history[m - 1], &self.y_history[m - 1]);
            if yty.abs() < 1e-12 {
                1.0
            } else {
                Self::dot(&self.s_history[m - 1], &self.y_history[m - 1]) / yty
            }
        } else {
            1.0
        };

        let mut r: Vec<f64> = q.iter().map(|qi| h0 * qi).collect();

        // Second loop: oldest pairs first
        for (i, (s_k, y_k)) in self.s_history.iter().zip(self.y_history.iter()).enumerate() {
            let beta = self.rho_history[i] * Self::dot(y_k, &r);
            let delta = alphas[i] - beta;
            for (rj, sij) in r.iter_mut().zip(s_k.iter()) {
                *rj += delta * sij;
            }
        }

        // direction = -r
        r.iter().map(|ri| -ri).collect()
    }

    fn update_history(&mut self, new_params: &[f64], new_grads: &[f64]) {
        if self.prev_params.len() != new_params.len() {
            return; // first step — nothing to diff against
        }
        let s: Vec<f64> = new_params
            .iter()
            .zip(self.prev_params.iter())
            .map(|(np, pp)| np - pp)
            .collect();
        let y: Vec<f64> = new_grads
            .iter()
            .zip(self.prev_grads.iter())
            .map(|(ng, pg)| ng - pg)
            .collect();
        let sy = Self::dot(&s, &y);
        if sy.abs() < 1e-20 {
            return; // skip curvature pair — not positive definite
        }
        let rho = 1.0 / sy;
        if self.s_history.len() >= self.history_size {
            self.s_history.remove(0);
            self.y_history.remove(0);
            self.rho_history.remove(0);
        }
        self.s_history.push(s);
        self.y_history.push(y);
        self.rho_history.push(rho);
    }
}

impl Optimizer for LBfgs {
    fn step(&mut self, params: &mut [f64], grads: &[f64]) {
        debug_assert_eq!(params.len(), grads.len());

        // On first step there is no history; fall back to scaled gradient descent.
        if self.step_count > 0 {
            self.update_history(params, grads);
        }

        let direction = self.two_loop_direction(grads);
        let lr = self.lr;
        for (p, d) in params.iter_mut().zip(direction.iter()) {
            *p += lr * d;
        }

        self.prev_params = params.to_vec();
        self.prev_grads = grads.to_vec();
        self.step_count += 1;
    }

    fn name(&self) -> &'static str {
        "L-BFGS"
    }

    fn reset(&mut self) {
        self.s_history.clear();
        self.y_history.clear();
        self.rho_history.clear();
        self.prev_params.clear();
        self.prev_grads.clear();
        self.step_count = 0;
    }
}

// ─── OptimizerKind ───────────────────────────────────────────────────────────

/// Selector enum used by [`ParameterFitter`] to choose the optimizer.
#[derive(Clone, Debug)]
pub enum OptimizerKind {
    /// SGD with configurable momentum.
    Sgd {
        /// Momentum coefficient ∈ [0, 1).
        momentum: f64,
    },
    /// Adam with default hyperparameters.
    Adam,
    /// L-BFGS with configurable history window.
    LBfgs {
        /// Number of curvature pairs to store.
        history: usize,
    },
}

// ─── FitResult ───────────────────────────────────────────────────────────────

/// Result returned by [`ParameterFitter::fit`].
#[derive(Clone, Debug)]
pub struct FitResult {
    /// Parameter values at the end of optimization.
    pub final_params: Vec<f64>,
    /// Loss value at the end of optimization.
    pub final_loss: f64,
    /// Number of gradient steps taken.
    pub n_iterations: usize,
    /// `true` if the convergence criterion was met before `max_iter`.
    pub converged: bool,
    /// Loss value after each iteration (length == `n_iterations`).
    pub loss_history: Vec<f64>,
}

// ─── ParameterFitter ─────────────────────────────────────────────────────────

/// High-level fitting loop that repeatedly:
/// 1. Creates a fresh [`Tape`](super::tape::Tape),
/// 2. Evaluates the differentiable loss function,
/// 3. Runs `tape.backward`,
/// 4. Extracts gradients, and
/// 5. Calls the underlying optimizer's `step`.
///
/// A new tape is needed each iteration because gradients accumulate; reusing a
/// tape would double-count contributions from earlier passes.
pub struct ParameterFitter {
    /// Starting parameter values.
    pub initial_params: Vec<f64>,
    optimizer_kind: OptimizerKind,
    /// Learning rate passed to the underlying optimizer.
    pub lr: f64,
    /// Maximum number of gradient steps before forced termination.
    pub max_iter: usize,
    /// Convergence tolerance: stop when |loss_t − loss_{t-1}| < tol.
    pub tol: f64,
}

impl ParameterFitter {
    /// Create a fitter using the Adam optimizer.
    pub fn new_adam(initial_params: Vec<f64>, lr: f64, max_iter: usize, tol: f64) -> Self {
        Self {
            initial_params,
            optimizer_kind: OptimizerKind::Adam,
            lr,
            max_iter,
            tol,
        }
    }

    /// Create a fitter using SGD with `momentum`.
    pub fn new_sgd(
        initial_params: Vec<f64>,
        lr: f64,
        momentum: f64,
        max_iter: usize,
        tol: f64,
    ) -> Self {
        Self {
            initial_params,
            optimizer_kind: OptimizerKind::Sgd { momentum },
            lr,
            max_iter,
            tol,
        }
    }

    /// Create a fitter using L-BFGS with `history` curvature pairs.
    pub fn new_lbfgs(
        initial_params: Vec<f64>,
        lr: f64,
        history: usize,
        max_iter: usize,
        tol: f64,
    ) -> Self {
        Self {
            initial_params,
            optimizer_kind: OptimizerKind::LBfgs { history },
            lr,
            max_iter,
            tol,
        }
    }

    /// Run the optimization loop and return a [`FitResult`].
    ///
    /// `loss_fn` receives a fresh tape and a slice of [`Var`](super::tape::Var)
    /// leaf variables (one per parameter) and must return a scalar loss variable.
    pub fn fit<F>(&self, loss_fn: F) -> FitResult
    where
        F: for<'a> Fn(&'a super::tape::Tape, &[super::tape::Var<'a>]) -> super::tape::Var<'a>,
    {
        use super::tape::{Tape, Var};

        let n = self.initial_params.len();
        let mut params: Vec<f64> = self.initial_params.clone();
        let mut loss_history: Vec<f64> = Vec::with_capacity(self.max_iter);
        let mut prev_loss = f64::INFINITY;
        let mut converged = false;
        let mut n_iterations = 0;

        // Build optimizer
        let mut opt: Box<dyn Optimizer> = match &self.optimizer_kind {
            OptimizerKind::Sgd { momentum } => {
                Box::new(Sgd::new(self.lr, *momentum).expect("SGD params validated"))
            },
            OptimizerKind::Adam => {
                let mut a = Adam::default_params(n);
                a.lr = self.lr;
                Box::new(a)
            },
            OptimizerKind::LBfgs { history } => {
                Box::new(LBfgs::new(self.lr, *history).expect("LBfgs params validated"))
            },
        };

        for _ in 0..self.max_iter {
            // Fresh tape each iteration
            let tape = Tape::new();
            let leaves: Vec<Var<'_>> = params.iter().map(|&p| Var::leaf(&tape, p)).collect();
            let loss_var = loss_fn(&tape, &leaves);
            let loss_val = loss_var.value();
            tape.backward(loss_var);
            let grads: Vec<f64> = leaves.iter().map(|lv| lv.grad()).collect();

            loss_history.push(loss_val);
            n_iterations += 1;

            opt.step(&mut params, &grads);

            if (loss_val - prev_loss).abs() < self.tol {
                converged = true;
                break;
            }
            prev_loss = loss_val;
        }

        let final_loss = *loss_history.last().unwrap_or(&f64::NAN);
        FitResult {
            final_params: params,
            final_loss,
            n_iterations,
            converged,
            loss_history,
        }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // 1. SGD vanilla (momentum=0) on f(x) = x² → converges to x ≈ 0
    #[test]
    fn test_sgd_vanilla_quadratic() {
        let fitter = ParameterFitter::new_sgd(vec![2.0], 0.1, 0.0, 500, 1e-10);
        let result = fitter.fit(|_tape, leaves| leaves[0] * leaves[0]);
        assert!(
            result.final_params[0].abs() < 0.01,
            "SGD should converge near 0"
        );
    }

    // 2. SGD with momentum on f(x) = x²
    #[test]
    fn test_sgd_momentum_quadratic() {
        let fitter_vanilla = ParameterFitter::new_sgd(vec![3.0], 0.05, 0.0, 400, 1e-10);
        let fitter_momentum = ParameterFitter::new_sgd(vec![3.0], 0.05, 0.5, 200, 1e-10);
        let r_vanilla = fitter_vanilla.fit(|_tape, leaves| leaves[0] * leaves[0]);
        let r_momentum = fitter_momentum.fit(|_tape, leaves| leaves[0] * leaves[0]);
        // Momentum version should converge in fewer or comparable iterations
        assert!(r_momentum.n_iterations <= r_vanilla.n_iterations + 50);
        assert!(r_momentum.final_params[0].abs() < 0.05);
    }

    // 3. Adam on f(x, y) = x² + y² → converges to (0, 0)
    #[test]
    fn test_adam_2d_quadratic() {
        let fitter = ParameterFitter::new_adam(vec![3.0, -4.0], 0.1, 1000, 1e-10);
        let result = fitter.fit(|_tape, leaves| leaves[0] * leaves[0] + leaves[1] * leaves[1]);
        assert!(
            result.final_params[0].abs() < 0.01,
            "Adam: x should converge to 0"
        );
        assert!(
            result.final_params[1].abs() < 0.01,
            "Adam: y should converge to 0"
        );
    }

    // 4. L-BFGS on quadratic → converges in few steps
    #[test]
    fn test_lbfgs_quadratic() {
        let fitter = ParameterFitter::new_lbfgs(vec![5.0], 0.5, 10, 200, 1e-10);
        let result = fitter.fit(|_tape, leaves| leaves[0] * leaves[0]);
        assert!(
            result.final_params[0].abs() < 0.1,
            "L-BFGS should converge near 0"
        );
    }

    // 5. ParameterFitter with Adam on simple quadratic
    #[test]
    fn test_parameter_fitter_adam() {
        let fitter = ParameterFitter::new_adam(vec![10.0], 0.05, 2000, 1e-12);
        let result = fitter.fit(|_tape, leaves| {
            let diff = leaves[0] - 3.0_f64; // f = (x - 3)²
            diff * diff
        });
        assert!((result.final_params[0] - 3.0).abs() < 0.1);
    }

    // 6. Reset between fits works
    #[test]
    fn test_optimizer_reset() {
        let mut opt = Adam::default_params(2);
        let mut params = vec![1.0, 1.0];
        let grads = vec![0.5, 0.5];
        opt.step(&mut params, &grads);
        assert!(opt.t > 0);
        opt.reset();
        assert_eq!(opt.t, 0);
        assert!(opt.m.iter().all(|&m| m == 0.0));
        assert!(opt.v.iter().all(|&v| v == 0.0));
    }

    // 7. Gradient direction check: params should move toward minimum
    #[test]
    fn test_gradient_direction() {
        // f(x) = x² minimized at 0; from x=5, grad=10 > 0, so x should decrease
        let mut sgd = Sgd::vanilla(0.1).unwrap();
        let mut params = vec![5.0];
        let grads = vec![10.0]; // ∂(x²)/∂x at x=5
        sgd.step(&mut params, &grads);
        assert!(params[0] < 5.0, "params should decrease toward minimum");
    }

    // 8. Convergence tolerance: result.converged when tol is loose
    #[test]
    fn test_convergence_tolerance() {
        let fitter = ParameterFitter::new_adam(vec![1.0], 0.1, 10000, 1e-4);
        let result = fitter.fit(|_tape, leaves| leaves[0] * leaves[0]);
        assert!(result.converged, "should converge with loose tolerance");
    }

    // 9. max_iter respected even without convergence
    #[test]
    fn test_max_iter_limit() {
        let max = 17;
        let fitter = ParameterFitter::new_sgd(vec![100.0], 1e-8, 0.0, max, 1e-300);
        let result = fitter.fit(|_tape, leaves| leaves[0] * leaves[0]);
        assert_eq!(result.n_iterations, max);
        assert!(!result.converged);
    }

    // 10. Loss history has correct length
    #[test]
    fn test_loss_history_length() {
        let max = 50;
        let fitter = ParameterFitter::new_adam(vec![5.0], 0.1, max, 1e-300);
        let result = fitter.fit(|_tape, leaves| leaves[0] * leaves[0]);
        assert_eq!(result.loss_history.len(), result.n_iterations);
        // For a convex function, loss should generally decrease (not strictly enforced here,
        // but check first is larger than last for well-chosen lr)
        if result.loss_history.len() > 1 {
            assert!(
                result.loss_history[0] >= result.loss_history[result.loss_history.len() - 1],
                "loss should not increase overall on convex f"
            );
        }
    }

    // 11. SGD invalid parameters are rejected
    #[test]
    fn test_sgd_invalid_params() {
        assert!(Sgd::new(-0.01, 0.0).is_err(), "negative lr should error");
        assert!(Sgd::new(0.1, 1.0).is_err(), "momentum=1 should error");
        assert!(
            Sgd::new(0.1, -0.1).is_err(),
            "negative momentum should error"
        );
    }

    // 12. Adam invalid parameters are rejected
    #[test]
    fn test_adam_invalid_params() {
        assert!(
            Adam::new(0.0, 0.9, 0.999, 1e-8, 2).is_err(),
            "lr=0 should error"
        );
        assert!(
            Adam::new(0.001, 1.0, 0.999, 1e-8, 2).is_err(),
            "beta1=1 should error"
        );
    }

    // 13. LBfgs invalid parameters are rejected
    #[test]
    fn test_lbfgs_invalid_params() {
        assert!(LBfgs::new(0.0, 5).is_err(), "lr=0 should error");
        assert!(LBfgs::new(0.1, 0).is_err(), "history=0 should error");
    }

    // 14. Manually verify one Adam step with known values
    #[test]
    fn test_adam_one_step_manual() {
        let mut adam = Adam::new(0.001, 0.9, 0.999, 1e-8, 1).unwrap();
        let mut params = vec![1.0_f64];
        let grads = vec![2.0_f64];
        adam.step(&mut params, &grads);
        // m = 0.9*0 + 0.1*2 = 0.2; v = 0.999*0 + 0.001*4 = 0.004
        // m_hat = 0.2/0.1 = 2; v_hat = 0.004/0.001 = 4
        // p = 1.0 - 0.001 * 2 / (sqrt(4) + 1e-8) ≈ 1.0 - 0.001
        let expected = 1.0 - 0.001 * 2.0 / (4.0_f64.sqrt() + 1e-8);
        assert!((params[0] - expected).abs() < 1e-12);
    }

    // 15. Raw optimizer step: two-variable SGD sanity check
    #[test]
    fn test_sgd_two_var() {
        let mut sgd = Sgd::new(0.5, 0.0).unwrap();
        let mut params = vec![2.0, -3.0];
        let grads = vec![4.0, -6.0]; // ∂/∂x at (x,y)=(2,-3) for x²+y²
        sgd.step(&mut params, &grads);
        // p = p - lr*g: [2 - 2, -3 - (-3)] = [0, 0]
        assert!((params[0] - 0.0).abs() < 1e-14);
        assert!((params[1] - 0.0).abs() < 1e-14);
    }

    // 16. Tape interaction: verify fitter creates fresh tape each iteration
    #[test]
    fn test_fresh_tape_per_iter() {
        // If tapes were reused, grads would accumulate and explode
        let fitter = ParameterFitter::new_sgd(vec![1.0], 0.1, 0.0, 100, 1e-10);
        let result = fitter.fit(|_tape, leaves| leaves[0] * leaves[0]);
        // With correct fresh tape, should converge; with accumulation, would diverge
        assert!(result.final_params[0].abs() < 0.5);
    }

    // Extra: verify that ParameterFitter uses the tape correctly in a composite loss
    #[test]
    fn test_composite_loss_fitting() {
        // f(x) = (x - 2)^2 + (x - 2)^4 — minimum at x = 2
        // Use f64 subtraction: (x - 2.0) where 2.0 is a scalar constant
        let fitter = ParameterFitter::new_adam(vec![5.0], 0.05, 3000, 1e-12);
        let result = fitter.fit(|_tape, leaves| {
            let x = leaves[0];
            let diff = x - 2.0_f64; // Sub<f64> keeps gradient path through x only
            diff * diff + (diff * diff) * (diff * diff)
        });
        assert!((result.final_params[0] - 2.0).abs() < 0.2);
    }
}
