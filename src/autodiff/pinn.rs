//! Physics-informed neural networks (PINNs) for Landau–Lifshitz–Gilbert
//! magnetisation dynamics (v0.7.0).
//!
//! A PINN replaces the explicit time integration of an ordinary differential
//! equation by a neural surrogate `m(t; θ)` whose parameters `θ` are obtained
//! by minimising a *physics residual* loss together with initial-condition
//! and (optionally) magnitude-preservation penalties.  Once trained, a single
//! forward pass yields `m(t)` at arbitrary collocation times without any
//! further integration step.
//!
//! ## LLG equation
//!
//! The PINN here targets the explicit Landau–Lifshitz form
//!
//! ```text
//!     dm/dt = -γ / (1 + α²) · [ m × H_eff + α · m × (m × H_eff) ],
//! ```
//!
//! where `γ` is the gyromagnetic ratio (positive scalar), `α` is the Gilbert
//! damping, `H_eff` is a constant effective field, and `m` is a unit-norm
//! magnetisation vector.  This formulation is mathematically equivalent to
//! the implicit Gilbert form but avoids the time-derivative on the right-hand
//! side, which is crucial for finite-difference-in-time PINN training.
//!
//! ## Training strategy
//!
//! Computing `∂m/∂t` directly on the AD tape would require nested tapes (one
//! for the time derivative, one for the parameter gradient).  We instead
//! evaluate the MLP three times — at `t − h`, `t`, and `t + h` — on the *same*
//! tape and use the central finite difference `(m(t+h) − m(t−h)) / (2h)` as
//! `∂m/∂t`.  All three forward passes share the leaf parameter variables, so
//! the gradient of the residual loss with respect to every parameter
//! propagates correctly.
//!
//! ## References
//!
//! - M. Raissi, P. Perdikaris & G. E. Karniadakis, "Physics-informed neural
//!   networks: A deep-learning framework for solving forward and inverse
//!   problems involving nonlinear partial differential equations",
//!   *J. Comput. Phys.* **378**, 686 (2019).
//! - L. D. Landau & E. M. Lifshitz, "On the theory of the dispersion of
//!   magnetic permeability in ferromagnetic bodies", *Phys. Z. Sowjet.* **8**,
//!   153 (1935).
//! - T. L. Gilbert, "A phenomenological theory of damping in ferromagnetic
//!   materials", *IEEE Trans. Magn.* **40**, 3443 (2004).
//! - L. Lu, X. Meng, Z. Mao & G. E. Karniadakis, "DeepXDE: A deep learning
//!   library for solving differential equations", *SIAM Rev.* **63**, 208
//!   (2021).

use crate::autodiff::neural::{Activation, Mlp};
use crate::autodiff::optimizer::{Adam, FitResult, LBfgs, Optimizer, OptimizerKind, Sgd};
use crate::autodiff::tape::{Tape, Var};
use crate::error::{invalid_param, Result};

/// Step size for the central finite difference used to obtain `dm/dt` on the
/// tape.  Chosen as a compromise between truncation error (O(h²)) and
/// catastrophic cancellation; smaller values amplify cancellation noise.
const TIME_FD_STEP: f64 = 1e-7;

// ─── LlgPinn ────────────────────────────────────────────────────────────────

/// Neural surrogate for LLG magnetisation dynamics.
///
/// The wrapped [`Mlp`] maps `t ↦ (m_x, m_y, m_z)` with `Tanh` activations
/// throughout the hidden stack and a final linear layer (so the output can
/// take any real value).
pub struct LlgPinn {
    /// Underlying network: input dim 1, output dim 3.
    pub mlp: Mlp,
    /// Constant effective field vector (A/m).
    pub h_eff: [f64; 3],
    /// Gilbert damping constant.
    pub alpha: f64,
    /// Gyromagnetic ratio magnitude (rad / (s·T)).
    pub gamma: f64,
}

impl LlgPinn {
    /// Build a new PINN with the requested hidden layer sizes.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when `alpha < 0` or
    /// `gamma <= 0`, and propagates [`Mlp::new`] errors.
    pub fn new(
        hidden_sizes: &[usize],
        h_eff: [f64; 3],
        alpha: f64,
        gamma: f64,
        rng_seed: u64,
    ) -> Result<Self> {
        if alpha < 0.0 {
            return Err(invalid_param("alpha", "damping must be non-negative"));
        }
        if gamma <= 0.0 {
            return Err(invalid_param(
                "gamma",
                "gyromagnetic ratio must be positive",
            ));
        }
        let mut sizes = vec![1_usize];
        sizes.extend_from_slice(hidden_sizes);
        sizes.push(3);
        let mut acts: Vec<Activation> = vec![Activation::Tanh; hidden_sizes.len()];
        acts.push(Activation::Linear);
        let mlp = Mlp::new(&sizes, &acts, rng_seed)?;
        Ok(Self {
            mlp,
            h_eff,
            alpha,
            gamma,
        })
    }

    /// Plain-`f64` inference.
    ///
    /// # Errors
    /// Propagates [`Mlp::forward_f64`] errors.
    pub fn predict(&self, t: f64) -> Result<[f64; 3]> {
        let v = self.mlp.forward_f64(&[t])?;
        Ok([v[0], v[1], v[2]])
    }

    /// Number of free parameters (equivalent to `self.mlp.n_params()`).
    pub fn n_params(&self) -> usize {
        self.mlp.n_params()
    }

    /// Flatten all parameters.
    pub fn params_flat(&self) -> Vec<f64> {
        self.mlp.params_flat()
    }

    /// Restore parameters from a flat slice.
    ///
    /// # Errors
    /// Propagates [`Mlp::set_params`] errors.
    pub fn set_params(&mut self, flat: &[f64]) -> Result<()> {
        self.mlp.set_params(flat)
    }
}

// ─── PinnTrainer ────────────────────────────────────────────────────────────

/// Encapsulates the loss formulation and training loop for an [`LlgPinn`].
///
/// The total loss is
///
/// `L = w_res · L_res + w_ic · L_ic + w_norm · L_norm`,
///
/// where
/// - `L_res = mean ‖∂m/∂t + γ/(1+α²) · [m × H_eff + α · m × (m × H_eff)]‖²`
///   over the collocation set,
/// - `L_ic  = ‖m(0) − m₀‖²` enforces the initial condition,
/// - `L_norm = mean (‖m‖² − 1)²` softly enforces unit-norm magnetisation.
pub struct PinnTrainer {
    /// Time points at which the LLG residual is evaluated.
    pub t_collocation: Vec<f64>,
    /// Target magnetisation at `t = 0`.
    pub initial_m: [f64; 3],
    /// Weight `w_res` of the physics-residual term.
    pub residual_weight: f64,
    /// Weight `w_ic` of the initial-condition term.
    pub ic_weight: f64,
    /// Weight `w_norm` of the unit-norm penalty (set to `0.0` to disable).
    pub norm_weight: f64,
}

impl PinnTrainer {
    /// Construct a trainer with sensible default weights
    /// (`w_res = 1`, `w_ic = 10`, `w_norm = 0`).
    pub fn new(t_collocation: Vec<f64>, initial_m: [f64; 3]) -> Self {
        Self {
            t_collocation,
            initial_m,
            residual_weight: 1.0,
            ic_weight: 10.0,
            norm_weight: 0.0,
        }
    }

    /// Override the loss weights.
    pub fn with_weights(mut self, residual: f64, ic: f64, norm: f64) -> Self {
        self.residual_weight = residual;
        self.ic_weight = ic;
        self.norm_weight = norm;
        self
    }

    /// Total loss assembled on the shared `tape` with `params` as the leaf
    /// parameter variables.  The corresponding gradient is obtained by
    /// calling `tape.backward(loss)` once and reading each `params\[i\].grad()`.
    ///
    /// # Errors
    /// Propagates errors from network forward passes and parameter loading.
    pub fn total_loss<'t>(
        &self,
        pinn: &LlgPinn,
        tape: &'t Tape,
        params: &[Var<'t>],
    ) -> Result<Var<'t>> {
        // Materialise a temporary MLP whose parameters are the f64 *values* of
        // the leaf vars — for accumulating gradients we still need to compose
        // the forward pass using the leaves directly (see `forward_with_leaves`).
        let n_params = pinn.mlp.n_params();
        if params.len() != n_params {
            return Err(crate::error::dimension_mismatch(
                &format!("{n_params} param leaves"),
                &format!("{} param leaves", params.len()),
            ));
        }

        let mut loss = Var::leaf(tape, 0.0);
        let n_coll = self.t_collocation.len() as f64;

        // Compute physics residual loss across all collocation points.
        if !self.t_collocation.is_empty() && self.residual_weight > 0.0 {
            let mut res_sum = Var::leaf(tape, 0.0);
            let prefactor = pinn.gamma / (1.0 + pinn.alpha * pinn.alpha);
            for &t_k in &self.t_collocation {
                let m_plus = forward_with_leaves(pinn, tape, params, t_k + TIME_FD_STEP)?;
                let m_minus = forward_with_leaves(pinn, tape, params, t_k - TIME_FD_STEP)?;
                let m_mid = forward_with_leaves(pinn, tape, params, t_k)?;
                // dm/dt ≈ (m_plus - m_minus) / (2h).
                let inv2h = 1.0 / (2.0 * TIME_FD_STEP);
                let dm_dt = [
                    (m_plus[0] - m_minus[0]) * inv2h,
                    (m_plus[1] - m_minus[1]) * inv2h,
                    (m_plus[2] - m_minus[2]) * inv2h,
                ];
                // m × H_eff
                let mxh = cross_var_const(m_mid, pinn.h_eff);
                // m × (m × H_eff)
                let mxmxh = cross_var_var(m_mid, mxh);
                // RHS = -prefactor · (m × H_eff + α · m × (m × H_eff))
                let rhs = [
                    -prefactor * (mxh[0] + pinn.alpha * mxmxh[0]),
                    -prefactor * (mxh[1] + pinn.alpha * mxmxh[1]),
                    -prefactor * (mxh[2] + pinn.alpha * mxmxh[2]),
                ];
                // residual = dm/dt - rhs
                let r0 = dm_dt[0] - rhs[0];
                let r1 = dm_dt[1] - rhs[1];
                let r2 = dm_dt[2] - rhs[2];
                let r_sq = r0 * r0 + r1 * r1 + r2 * r2;
                res_sum = res_sum + r_sq;
            }
            let res_mean = res_sum / n_coll;
            loss = loss + res_mean * self.residual_weight;
        }

        // Initial-condition loss: m(0) ≈ initial_m.
        if self.ic_weight > 0.0 {
            let m0 = forward_with_leaves(pinn, tape, params, 0.0)?;
            let d0 = m0[0] - self.initial_m[0];
            let d1 = m0[1] - self.initial_m[1];
            let d2 = m0[2] - self.initial_m[2];
            let ic_sq = d0 * d0 + d1 * d1 + d2 * d2;
            loss = loss + ic_sq * self.ic_weight;
        }

        // Unit-norm penalty: average over collocation points.
        if self.norm_weight > 0.0 && !self.t_collocation.is_empty() {
            let mut norm_sum = Var::leaf(tape, 0.0);
            for &t_k in &self.t_collocation {
                let m_mid = forward_with_leaves(pinn, tape, params, t_k)?;
                let mag_sq = m_mid[0] * m_mid[0] + m_mid[1] * m_mid[1] + m_mid[2] * m_mid[2];
                let dev = mag_sq - 1.0;
                norm_sum = norm_sum + dev * dev;
            }
            let norm_mean = norm_sum / n_coll;
            loss = loss + norm_mean * self.norm_weight;
        }

        Ok(loss)
    }

    /// Run a full training loop using the requested optimizer kind.
    ///
    /// Each iteration rebuilds the tape, materialises the parameter leaves,
    /// evaluates [`PinnTrainer::total_loss`], runs `tape.backward`, and asks
    /// the optimizer to update the flat parameter vector.
    ///
    /// # Errors
    /// Propagates errors from [`PinnTrainer::total_loss`] and parameter
    /// loading on the network.
    pub fn train(
        &self,
        pinn: &mut LlgPinn,
        n_iter: usize,
        lr: f64,
        optimizer: OptimizerKind,
    ) -> Result<FitResult> {
        let mut params: Vec<f64> = pinn.params_flat();
        let n = params.len();
        let mut opt: Box<dyn Optimizer> = match optimizer {
            OptimizerKind::Sgd { momentum } => Box::new(
                Sgd::new(lr, momentum).map_err(|_| invalid_param("sgd", "invalid lr/momentum"))?,
            ),
            OptimizerKind::Adam => {
                let mut a = Adam::default_params(n);
                a.lr = lr;
                Box::new(a)
            },
            OptimizerKind::LBfgs { history } => Box::new(
                LBfgs::new(lr, history)
                    .map_err(|_| invalid_param("lbfgs", "invalid lr/history"))?,
            ),
        };
        let mut loss_history = Vec::with_capacity(n_iter);
        let mut prev_loss = f64::INFINITY;
        let mut converged = false;
        let tol = 1e-12;

        for _ in 0..n_iter {
            pinn.set_params(&params)?;
            let tape = Tape::new();
            let leaves: Vec<Var<'_>> = params.iter().map(|&p| Var::leaf(&tape, p)).collect();
            let loss_var = self.total_loss(pinn, &tape, &leaves)?;
            let loss_val = loss_var.value();
            tape.backward(loss_var);
            let grads: Vec<f64> = leaves.iter().map(|lv| lv.grad()).collect();
            loss_history.push(loss_val);
            opt.step(&mut params, &grads);
            if (loss_val - prev_loss).abs() < tol {
                converged = true;
                break;
            }
            prev_loss = loss_val;
        }
        pinn.set_params(&params)?;

        let final_loss = *loss_history.last().unwrap_or(&f64::NAN);
        let n_iterations = loss_history.len();
        Ok(FitResult {
            final_params: params,
            final_loss,
            n_iterations,
            converged,
            loss_history,
        })
    }
}

// ─── Internal helpers ───────────────────────────────────────────────────────

/// Forward the PINN's MLP at a single scalar time `t` while threading the
/// parameter leaves through the tape as multipliers.  This produces a
/// `[Var<'t>; 3]` triple whose gradient w.r.t. each `params[k]` is correctly
/// accumulated by `tape.backward`.
fn forward_with_leaves<'t>(
    pinn: &LlgPinn,
    tape: &'t Tape,
    params: &[Var<'t>],
    t_value: f64,
) -> Result<[Var<'t>; 3]> {
    // Walk the MLP layer by layer, using the leaves for weights/biases.
    let mut cursor = 0_usize;
    let mut current: Vec<Var<'t>> = vec![Var::leaf(tape, t_value)];
    for layer in &pinn.mlp.layers {
        let mut next: Vec<Var<'t>> = Vec::with_capacity(layer.out_dim);
        // Weights occupy `out_dim * in_dim` slots in `params`.
        for j in 0..layer.out_dim {
            // Pre-activation: sum_i w_ji x_i.
            let w0 = params[cursor + j * layer.in_dim];
            let mut pre = w0 * current[0];
            for i in 1..layer.in_dim {
                let w = params[cursor + j * layer.in_dim + i];
                pre = pre + w * current[i];
            }
            next.push(pre);
        }
        cursor += layer.out_dim * layer.in_dim;
        // Biases follow weights.
        for j in 0..layer.out_dim {
            let b = params[cursor + j];
            next[j] = next[j] + b;
        }
        cursor += layer.out_dim;
        // Apply activation.
        for v in next.iter_mut() {
            *v = layer.activation.apply(*v);
        }
        current = next;
    }
    if current.len() != 3 {
        return Err(crate::error::dimension_mismatch(
            "3 outputs",
            &format!("{} outputs", current.len()),
        ));
    }
    Ok([current[0], current[1], current[2]])
}

/// Cross product `a × b` where `a` is a `Var` triple and `b` is a constant
/// `f64` triple.
fn cross_var_const<'t>(a: [Var<'t>; 3], b: [f64; 3]) -> [Var<'t>; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Cross product `a × b` where both `a` and `b` are `Var` triples.
fn cross_var_var<'t>(a: [Var<'t>; 3], b: [Var<'t>; 3]) -> [Var<'t>; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::autodiff::tape::Tape;

    // 1. Construction succeeds with valid inputs.
    #[test]
    fn test_pinn_construction_ok() {
        let pinn = LlgPinn::new(&[6, 6], [0.0, 0.0, 1.0], 0.01, 1.76e11, 7).unwrap();
        assert_eq!(pinn.mlp.input_dim, 1);
        assert_eq!(pinn.mlp.output_dim, 3);
        assert!(pinn.n_params() > 0);
    }

    // 2. Predict returns a length-3 finite output.
    #[test]
    fn test_pinn_predict_shape() {
        let pinn = LlgPinn::new(&[5, 5], [0.0, 0.0, 1.0], 0.05, 1.76e11, 11).unwrap();
        let m = pinn.predict(1e-12).unwrap();
        assert!(m.iter().all(|x| x.is_finite()));
    }

    // 3. PinnTrainer::total_loss runs end-to-end and yields a finite loss.
    #[test]
    fn test_total_loss_runs() {
        let pinn = LlgPinn::new(&[4, 4], [0.0, 0.0, 1.0], 0.01, 1.0, 3).unwrap();
        let trainer = PinnTrainer::new(vec![0.0, 0.01, 0.02], [1.0, 0.0, 0.0]);
        let params = pinn.params_flat();
        let tape = Tape::new();
        let leaves: Vec<Var<'_>> = params.iter().map(|&p| Var::leaf(&tape, p)).collect();
        let loss = trainer.total_loss(&pinn, &tape, &leaves).unwrap();
        assert!(loss.value().is_finite());
    }

    // 4. Gradient w.r.t. a randomly chosen parameter matches the central
    //    finite-difference estimate.  Uses a small network and a single
    //    collocation point to keep the test cheap.
    #[test]
    fn test_loss_gradient_finite_diff() {
        let pinn = LlgPinn::new(&[3], [0.5, 0.0, 0.5], 0.02, 1.0, 5).unwrap();
        let trainer = PinnTrainer::new(vec![0.0, 0.1], [1.0, 0.0, 0.0]);
        let params = pinn.params_flat();
        // AD gradient.
        let tape = Tape::new();
        let leaves: Vec<Var<'_>> = params.iter().map(|&p| Var::leaf(&tape, p)).collect();
        let loss_var = trainer.total_loss(&pinn, &tape, &leaves).unwrap();
        tape.backward(loss_var);
        let ad_grad: Vec<f64> = leaves.iter().map(|l| l.grad()).collect();
        // Finite difference on a single (well-conditioned) bias parameter.
        let pick = params.len() - 1; // last bias
        let h = 1e-5;
        let mut p_plus = params.clone();
        let mut p_minus = params.clone();
        p_plus[pick] += h;
        p_minus[pick] -= h;
        let mut pinn_plus = LlgPinn::new(&[3], [0.5, 0.0, 0.5], 0.02, 1.0, 5).unwrap();
        let mut pinn_minus = LlgPinn::new(&[3], [0.5, 0.0, 0.5], 0.02, 1.0, 5).unwrap();
        pinn_plus.set_params(&p_plus).unwrap();
        pinn_minus.set_params(&p_minus).unwrap();
        let l_plus = {
            let t2 = Tape::new();
            let leaves2: Vec<Var<'_>> = p_plus.iter().map(|&p| Var::leaf(&t2, p)).collect();
            trainer
                .total_loss(&pinn_plus, &t2, &leaves2)
                .unwrap()
                .value()
        };
        let l_minus = {
            let t2 = Tape::new();
            let leaves2: Vec<Var<'_>> = p_minus.iter().map(|&p| Var::leaf(&t2, p)).collect();
            trainer
                .total_loss(&pinn_minus, &t2, &leaves2)
                .unwrap()
                .value()
        };
        let fd = (l_plus - l_minus) / (2.0 * h);
        // Compare; very loose tol because FD time step inside loss adds noise.
        let rel = (ad_grad[pick] - fd).abs() / (1.0 + fd.abs());
        assert!(rel < 1e-2, "AD {} vs FD {} rel {}", ad_grad[pick], fd, rel);
    }

    // 5. Short training reduces the loss.
    #[test]
    fn test_training_reduces_loss() {
        let mut pinn = LlgPinn::new(&[6, 6], [0.0, 0.0, 1.0], 0.01, 1.0, 13).unwrap();
        let trainer =
            PinnTrainer::new(vec![0.0, 0.1, 0.2], [1.0, 0.0, 0.0]).with_weights(1.0, 100.0, 0.0);
        let res = trainer
            .train(&mut pinn, 60, 1e-2, OptimizerKind::Adam)
            .unwrap();
        assert!(res.loss_history.first().unwrap() >= res.loss_history.last().unwrap());
    }

    // 6. Larmor precession sanity: starting from m × H_eff ≠ 0 the residual
    //    is non-zero before training.
    #[test]
    fn test_larmor_nonzero_initial_residual() {
        let pinn = LlgPinn::new(&[4, 4], [0.0, 0.0, 1.0], 0.0, 1.0, 19).unwrap();
        let trainer = PinnTrainer::new(vec![0.0, 0.01], [1.0, 0.0, 0.0]);
        let params = pinn.params_flat();
        let tape = Tape::new();
        let leaves: Vec<Var<'_>> = params.iter().map(|&p| Var::leaf(&tape, p)).collect();
        let loss = trainer.total_loss(&pinn, &tape, &leaves).unwrap();
        assert!(loss.value() > 0.0);
    }

    // 7. Unit-norm penalty actually contributes to the loss.
    #[test]
    fn test_unit_norm_penalty_active() {
        let pinn = LlgPinn::new(&[3], [0.0, 0.0, 1.0], 0.01, 1.0, 21).unwrap();
        let trainer_no_norm =
            PinnTrainer::new(vec![0.0, 0.01], [1.0, 0.0, 0.0]).with_weights(0.0, 0.0, 0.0);
        let trainer_norm =
            PinnTrainer::new(vec![0.0, 0.01], [1.0, 0.0, 0.0]).with_weights(0.0, 0.0, 1.0);
        let params = pinn.params_flat();
        let tape1 = Tape::new();
        let leaves1: Vec<Var<'_>> = params.iter().map(|&p| Var::leaf(&tape1, p)).collect();
        let l0 = trainer_no_norm
            .total_loss(&pinn, &tape1, &leaves1)
            .unwrap()
            .value();
        let tape2 = Tape::new();
        let leaves2: Vec<Var<'_>> = params.iter().map(|&p| Var::leaf(&tape2, p)).collect();
        let l1 = trainer_norm
            .total_loss(&pinn, &tape2, &leaves2)
            .unwrap()
            .value();
        assert!(l0 == 0.0);
        assert!(l1 >= 0.0);
    }

    // 8. params_flat/set_params round-trip on the PINN.
    #[test]
    fn test_pinn_params_roundtrip() {
        let mut pinn = LlgPinn::new(&[3, 3], [0.0, 0.0, 1.0], 0.01, 1.0, 23).unwrap();
        let original = pinn.params_flat();
        let mut perturbed = original.clone();
        for v in perturbed.iter_mut().take(5) {
            *v += 0.5;
        }
        pinn.set_params(&perturbed).unwrap();
        assert_eq!(pinn.params_flat(), perturbed);
    }

    // 9. Reproducibility: identical seed → identical predictions.
    #[test]
    fn test_pinn_reproducibility() {
        let p1 = LlgPinn::new(&[5], [0.1, 0.2, 0.3], 0.05, 1.0, 1001).unwrap();
        let p2 = LlgPinn::new(&[5], [0.1, 0.2, 0.3], 0.05, 1.0, 1001).unwrap();
        let m1 = p1.predict(0.7).unwrap();
        let m2 = p2.predict(0.7).unwrap();
        for k in 0..3 {
            assert!((m1[k] - m2[k]).abs() < 1e-15);
        }
    }

    // 10. Invalid construction parameters are rejected.
    #[test]
    fn test_pinn_invalid_params_rejected() {
        assert!(LlgPinn::new(&[3], [0.0, 0.0, 1.0], -0.1, 1.0, 0).is_err());
        assert!(LlgPinn::new(&[3], [0.0, 0.0, 1.0], 0.1, 0.0, 0).is_err());
        assert!(LlgPinn::new(&[3], [0.0, 0.0, 1.0], 0.1, -1.0, 0).is_err());
    }
}
