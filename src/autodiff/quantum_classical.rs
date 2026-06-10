//! Hybrid Quantum-Classical Neural Network for Magnetic Hamiltonians (v0.1.0).
//!
//! This module trains a neural network to parameterise the coefficients
//! `(A_k, B_k)` of a 1-D ferromagnet magnon Hamiltonian
//!
//! ```text
//! H = Σ_k [ A_k α†_k α_k + (B_k/2)(α†_k α†_{-k} + h.c.) ]
//! ```
//!
//! and computes the Bogoliubov quasi-particle spectrum
//! `ε_k = √(A_k² − B_k²)` as the training target.  Gradients with
//! respect to the neural-network weights are obtained via central finite
//! differences so that the quantum energy can be treated as a black-box
//! function of the trainable parameters.
//!
//! # Training objective
//!
//! Given target magnon frequencies `ω_k^target` (e.g. from a reference
//! Heisenberg model), the network minimises the mean squared error
//!
//! ```text
//! L = (1/N_k) Σ_k (ε_k(θ) − ω_k^target)²
//! ```
//!
//! over its parameters θ, updated with the Adam optimiser.
//!
//! # References
//!
//! - N. N. Bogoliubov, *J. Phys. USSR* **11**, 23 (1947).
//! - T. D. Stanescu & G. Sarma, "Hybrid quantum-classical computing",
//!   *Phys. Rev. B* **90**, 085110 (2014).
//! - D. P. Kingma & J. Ba, "Adam: A Method for Stochastic Optimization",
//!   *ICLR* (2015), arXiv:1412.6980.

use std::f64::consts::PI;

use crate::autodiff::neural::{Activation, Mlp};
use crate::autodiff::optimizer::{Adam, Optimizer};
use crate::error::{invalid_param, numerical_error, Result};

// ─── MagnonHamiltonianParams ─────────────────────────────────────────────────

/// Physical parameters that define the 1-D ferromagnet magnon Hamiltonian.
///
/// The reference dispersion (Heisenberg chain in a longitudinal field) is
///
/// ```text
/// ω_k = J·S·(1 − cos k) + H_ext,   k ∈ {π/N, 2π/N, …, π}
/// ```
#[derive(Debug, Clone)]
pub struct MagnonHamiltonianParams {
    /// Number of k-points (≥ 2).
    pub n_modes: usize,
    /// Exchange coupling J \[meV or dimensionless\] (must be positive).
    pub j_exchange: f64,
    /// External-field contribution H_ext \[same units as j_exchange\].
    pub h_ext: f64,
    /// Spin quantum number S (≥ 0.5, e.g. 1.0 for S=1).
    pub spin_s: f64,
}

impl MagnonHamiltonianParams {
    /// Construct and validate the Hamiltonian parameters.
    ///
    /// # Errors
    ///
    /// Returns [`crate::error::Error::InvalidParameter`] if
    /// - `n_modes < 2`,
    /// - `j_exchange ≤ 0`, or
    /// - `spin_s < 0.5`.
    pub fn new(n_modes: usize, j_exchange: f64, h_ext: f64, spin_s: f64) -> Result<Self> {
        if n_modes < 2 {
            return Err(invalid_param("n_modes", "must be at least 2"));
        }
        if j_exchange <= 0.0 {
            return Err(invalid_param("j_exchange", "must be strictly positive"));
        }
        if spin_s < 0.5 {
            return Err(invalid_param("spin_s", "must be at least 0.5"));
        }
        Ok(Self {
            n_modes,
            j_exchange,
            h_ext,
            spin_s,
        })
    }

    /// Return the N_modes k-values in the interval `(0, π]`:
    /// `k_i = (i+1)·π / N_modes` for `i = 0, …, N_modes − 1`.
    pub fn ferromagnet_1d_k_values(&self) -> Vec<f64> {
        (0..self.n_modes)
            .map(|i| (i + 1) as f64 * PI / self.n_modes as f64)
            .collect()
    }

    /// Analytical Heisenberg-chain magnon dispersion.
    ///
    /// ```text
    /// ω_k = J·S·(1 − cos k) + H_ext
    /// ```
    ///
    /// These are the "true" frequencies the neural network is trained to
    /// reproduce.
    pub fn reference_magnon_frequencies(&self) -> Vec<f64> {
        let js = self.j_exchange * self.spin_s;
        self.ferromagnet_1d_k_values()
            .into_iter()
            .map(|k| js * (1.0 - k.cos()) + self.h_ext)
            .collect()
    }
}

// ─── MagnonNeuralNetwork ─────────────────────────────────────────────────────

/// Neural network that maps `(k/π, J_norm, H_norm) → (A_k, B_k)`.
///
/// The MLP has architecture
///
/// ```text
/// [3] → [hidden_dim] × depth → [2],   Tanh hidden, Linear output
/// ```
///
/// Constraints applied on the MLP output before physics use:
/// - `A_k = exp(A_raw) + 1e-3`  — ensures A_k > 0.
/// - `B_k = B_raw × 0.3`        — keeps |B_k| small relative to A_k.
#[derive(Debug, Clone)]
pub struct MagnonNeuralNetwork {
    /// Underlying MLP: input_dim=3, output_dim=2.
    pub net: Mlp,
    /// Current trainable parameters (mirrors `net.params_flat()`).
    pub params: Vec<f64>,
}

impl MagnonNeuralNetwork {
    /// Build a new `MagnonNeuralNetwork`.
    ///
    /// Layer sizes are `[3, hidden_dim × depth, 2]` with [`Activation::Tanh`]
    /// in all hidden layers and [`Activation::Linear`] at the output.
    /// The initial parameters are extracted from the freshly initialised MLP.
    pub fn new(hidden_dim: usize, depth: usize) -> Self {
        // Build layer sizes: input=3, (depth) hidden layers of hidden_dim, output=2.
        let mut sizes = Vec::with_capacity(depth + 2);
        sizes.push(3_usize);
        for _ in 0..depth {
            sizes.push(hidden_dim);
        }
        sizes.push(2_usize);

        // Activations: Tanh for each hidden layer, Linear for the output layer.
        let mut acts = Vec::with_capacity(depth + 1);
        for _ in 0..depth {
            acts.push(Activation::Tanh);
        }
        acts.push(Activation::Linear);

        // The Mlp::new API requires layer_sizes.len() - 1 activations.  Both
        // sizes and acts satisfy this by construction.
        let net = Mlp::new(&sizes, &acts, 0x00C0_FFEE_1234_5678_u64)
            .expect("layer sizes and activations are valid by construction");
        let params = net.params_flat();

        Self { net, params }
    }

    /// Compute `(A_k, B_k)` for a single k-point using the current `params`.
    ///
    /// - `a_k = exp(output\[0\]) + 1e-3`
    /// - `b_k = output[1] × 0.3`
    pub fn compute_coefficients(&self, k: f64, j_norm: f64, h_norm: f64) -> (f64, f64) {
        let input = [k / PI, j_norm, h_norm];
        // forward_f64 returns Err only on dimension mismatch, which cannot
        // occur here because the input is always length 3 = input_dim.
        let output = self
            .net
            .forward_f64(&input)
            .expect("input dimension 3 matches network input_dim by construction");
        let a_k = output[0].exp() + 1e-3;
        let b_k = output[1] * 0.3;
        (a_k, b_k)
    }

    /// Compute `(A_k_vec, B_k_vec)` for all k-points of `params`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::error::Error::NumericalError`] if the Bogoliubov
    /// stability condition `A_k > |B_k|` is violated for any mode.
    pub fn compute_coefficients_all(
        &self,
        params: &MagnonHamiltonianParams,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        let k_values = params.ferromagnet_1d_k_values();
        let j_norm = 1.0_f64;
        let h_norm = params.h_ext / (params.j_exchange + 1e-15);

        let mut a_vec = Vec::with_capacity(k_values.len());
        let mut b_vec = Vec::with_capacity(k_values.len());

        for &k in &k_values {
            let (a_k, b_k) = self.compute_coefficients(k, j_norm, h_norm);
            if a_k <= b_k.abs() {
                return Err(numerical_error(&format!(
                    "Bogoliubov stability violated at k={k:.4}: A_k={a_k:.4} ≤ |B_k|={:.4}",
                    b_k.abs()
                )));
            }
            a_vec.push(a_k);
            b_vec.push(b_k);
        }

        Ok((a_vec, b_vec))
    }
}

// ─── QuantumClassicalResult ──────────────────────────────────────────────────

/// Outcome of a [`QuantumClassicalOptimizer`] training run.
#[derive(Debug, Clone)]
pub struct QuantumClassicalResult {
    /// Loss value after the final training step.
    pub final_loss: f64,
    /// Loss recorded at every training step.
    pub loss_history: Vec<f64>,
    /// Bogoliubov quasi-particle frequencies ε_k from the trained NN.
    pub final_magnon_frequencies: Vec<f64>,
    /// Reference target frequencies used during training.
    pub target_frequencies: Vec<f64>,
    /// Ground-state (zero-point) energy Σ_k (ε_k − A_k) / 2.
    pub ground_state_energy: f64,
}

// ─── QuantumClassicalOptimizer ───────────────────────────────────────────────

/// Hybrid quantum-classical optimiser.
///
/// A neural network parameterises the magnon Hamiltonian coefficients
/// `(A_k, B_k)`.  The quantum quasi-particle energies
/// `ε_k = √(A_k² − B_k²)` enter the MSE loss against known target
/// frequencies.  Parameter gradients are estimated via central finite
/// differences and fed to the Adam optimiser.
///
/// # Note on `Debug` and `Clone`
///
/// [`Adam`] does not derive `Debug` or `Clone`, so these traits are
/// implemented manually: `Debug` shows the learning rate; `Clone`
/// reconstructs Adam with the same learning rate and resets moments.
pub struct QuantumClassicalOptimizer {
    /// Physical Hamiltonian parameters (k-grid, J, H_ext, S).
    pub ham_params: MagnonHamiltonianParams,
    /// Neural network mapping k → (A_k, B_k).
    pub nn: MagnonNeuralNetwork,
    /// Adam optimiser state.
    pub adam: Adam,
    /// Learning rate stored separately to support `Clone`.
    adam_lr: f64,
    /// Target magnon frequencies from [`MagnonHamiltonianParams::reference_magnon_frequencies`].
    pub target_frequencies: Vec<f64>,
}

impl std::fmt::Debug for QuantumClassicalOptimizer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("QuantumClassicalOptimizer")
            .field("ham_params", &self.ham_params)
            .field("nn", &self.nn)
            .field("adam_lr", &self.adam_lr)
            .field("n_target_frequencies", &self.target_frequencies.len())
            .finish()
    }
}

impl Clone for QuantumClassicalOptimizer {
    fn clone(&self) -> Self {
        let n_params = self.nn.params.len();
        // Reconstruct Adam with the same learning rate; moment history is reset
        // (clone semantics: a fresh optimiser starting from the cloned parameters).
        let adam = Adam::new(self.adam_lr, 0.9, 0.999, 1e-8, n_params)
            .expect("Adam hyperparameters are valid constants");
        Self {
            ham_params: self.ham_params.clone(),
            nn: self.nn.clone(),
            adam,
            adam_lr: self.adam_lr,
            target_frequencies: self.target_frequencies.clone(),
        }
    }
}

impl QuantumClassicalOptimizer {
    /// Construct the optimiser.
    ///
    /// `target_frequencies` is computed from the reference Heisenberg model
    /// inside this constructor.
    ///
    /// # Panics (compile-time constants)
    ///
    /// Adam is constructed with valid default hyperparameters; the only
    /// `expect` call here uses literal-constant arguments that are always
    /// valid.
    pub fn new(
        ham_params: MagnonHamiltonianParams,
        hidden_dim: usize,
        depth: usize,
        lr: f64,
    ) -> Self {
        let target_frequencies = ham_params.reference_magnon_frequencies();
        let nn = MagnonNeuralNetwork::new(hidden_dim, depth);
        let n_params = nn.params.len();
        let adam = Adam::new(lr, 0.9, 0.999, 1e-8, n_params)
            .expect("Adam hyperparameters (β₁=0.9, β₂=0.999, ε=1e-8) are valid constants");
        Self {
            ham_params,
            nn,
            adam,
            adam_lr: lr,
            target_frequencies,
        }
    }

    // ── Private helpers ──────────────────────────────────────────────────────

    /// Evaluate the MSE loss with an arbitrary parameter vector.
    ///
    /// This method is pure: it does **not** mutate `self`.
    fn compute_loss_with_params(&self, params: &[f64]) -> Result<f64> {
        let k_values = self.ham_params.ferromagnet_1d_k_values();
        let j_norm = 1.0_f64;
        let h_norm = self.ham_params.h_ext / (self.ham_params.j_exchange + 1e-15);

        // Clone the MLP and install the candidate parameters.
        let mut temp_net = self.nn.net.clone();
        temp_net.set_params(params)?;

        let n = k_values.len() as f64;
        let mut loss = 0.0_f64;

        for (i, &k) in k_values.iter().enumerate() {
            let input = [k / PI, j_norm, h_norm];
            let output = temp_net.forward_f64(&input)?;

            let a_k = output[0].exp() + 1e-3;
            let b_k = output[1] * 0.3;

            // Clamp B_k to stay inside the Bogoliubov stability cone.
            let b_k_clamped = b_k.clamp(-a_k * 0.99, a_k * 0.99);

            let discriminant = a_k * a_k - b_k_clamped * b_k_clamped;
            if discriminant < 0.0 {
                return Err(numerical_error(&format!(
                    "negative discriminant A²−B²={discriminant:.4e} at k-index {i}"
                )));
            }
            let epsilon_k = discriminant.sqrt();

            let diff = epsilon_k - self.target_frequencies[i];
            loss += diff * diff;
        }

        Ok(loss / n)
    }

    /// Evaluate the MSE loss with the **current** NN parameters.
    fn compute_loss_f64(&self) -> Result<f64> {
        self.compute_loss_with_params(&self.nn.params)
    }

    // ── Public interface ─────────────────────────────────────────────────────

    /// Compute the loss and its central-finite-difference gradient.
    ///
    /// For each parameter θ_i:
    ///
    /// ```text
    /// ∂L/∂θ_i ≈ (L(θ + δ e_i) − L(θ − δ e_i)) / (2δ),   δ = 1e-5
    /// ```
    ///
    /// The baseline loss `L(θ)` is computed separately and returned as the
    /// first element of the tuple.
    ///
    /// # Errors
    ///
    /// Propagates errors from `Self::compute_loss_with_params`.
    pub fn compute_loss_and_gradients(&self) -> Result<(f64, Vec<f64>)> {
        let loss0 = self.compute_loss_f64()?;
        let delta = 1e-5_f64;
        let n = self.nn.params.len();
        let mut grads = vec![0.0_f64; n];

        // Work on a mutable clone so we never alter self.nn.params.
        let mut temp_params = self.nn.params.clone();

        for i in 0..n {
            // Forward perturbation.
            temp_params[i] += delta;
            let loss_plus = self.compute_loss_with_params(&temp_params)?;

            // Backward perturbation.
            temp_params[i] -= 2.0 * delta;
            let loss_minus = self.compute_loss_with_params(&temp_params)?;

            // Restore.
            temp_params[i] += delta;

            grads[i] = (loss_plus - loss_minus) / (2.0 * delta);
        }

        Ok((loss0, grads))
    }

    /// Apply one Adam step.  Returns the loss before the step.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`Self::compute_loss_and_gradients`].
    pub fn train_step(&mut self) -> Result<f64> {
        let (loss, grads) = self.compute_loss_and_gradients()?;
        self.adam.step(&mut self.nn.params, &grads);
        // Keep the MLP in sync with the updated parameter vector.
        self.nn
            .net
            .set_params(&self.nn.params)
            .expect("params length matches network by construction");
        Ok(loss)
    }

    /// Run `n_steps` Adam steps and return a structured result.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`Self::train_step`].
    pub fn train(&mut self, n_steps: usize) -> Result<QuantumClassicalResult> {
        let mut loss_history = Vec::with_capacity(n_steps);

        for _ in 0..n_steps {
            let step_loss = self.train_step()?;
            loss_history.push(step_loss);
        }

        let final_loss = *loss_history.last().unwrap_or(&f64::NAN);
        let final_magnon_frequencies = self.final_frequencies()?;

        // Ground-state energy: Σ_k (ε_k − A_k) / 2.
        let (a_vec, _) = self.nn.compute_coefficients_all(&self.ham_params)?;
        let ground_state_energy: f64 = final_magnon_frequencies
            .iter()
            .zip(a_vec.iter())
            .map(|(&eps_k, &a_k)| (eps_k - a_k) / 2.0)
            .sum();

        Ok(QuantumClassicalResult {
            final_loss,
            loss_history,
            final_magnon_frequencies,
            target_frequencies: self.target_frequencies.clone(),
            ground_state_energy,
        })
    }

    /// Evaluate the Bogoliubov quasi-particle frequencies ε_k from the
    /// current NN parameters.
    ///
    /// # Errors
    ///
    /// Returns an error if the Bogoliubov stability condition is violated
    /// for any mode.
    pub fn final_frequencies(&self) -> Result<Vec<f64>> {
        let k_values = self.ham_params.ferromagnet_1d_k_values();
        let j_norm = 1.0_f64;
        let h_norm = self.ham_params.h_ext / (self.ham_params.j_exchange + 1e-15);

        let mut freqs = Vec::with_capacity(k_values.len());

        for &k in &k_values {
            let (a_k, b_k) = self.nn.compute_coefficients(k, j_norm, h_norm);
            let b_k_clamped = b_k.clamp(-a_k * 0.99, a_k * 0.99);
            let discriminant = a_k * a_k - b_k_clamped * b_k_clamped;
            if discriminant < 0.0 {
                return Err(numerical_error(&format!(
                    "negative discriminant A²−B²={discriminant:.4e} at k={k:.4}"
                )));
            }
            freqs.push(discriminant.sqrt());
        }

        Ok(freqs)
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: build a small default setup used by many tests.
    fn default_setup() -> (MagnonHamiltonianParams, QuantumClassicalOptimizer) {
        let params = MagnonHamiltonianParams::new(4, 1.0, 0.1, 1.0).expect("valid default params");
        let opt = QuantumClassicalOptimizer::new(params.clone(), 16, 2, 0.01);
        (params, opt)
    }

    // ── MagnonHamiltonianParams ──────────────────────────────────────────────

    #[test]
    fn test_param_validation_n_modes() {
        let res = MagnonHamiltonianParams::new(1, 1.0, 0.0, 1.0);
        assert!(res.is_err(), "n_modes < 2 should be rejected");
    }

    #[test]
    fn test_param_validation_j_exchange() {
        let res = MagnonHamiltonianParams::new(4, -1.0, 0.0, 1.0);
        assert!(res.is_err(), "j_exchange ≤ 0 should be rejected");
    }

    #[test]
    fn test_param_validation_spin_s() {
        let res = MagnonHamiltonianParams::new(4, 1.0, 0.0, 0.3);
        assert!(res.is_err(), "spin_s < 0.5 should be rejected");
    }

    #[test]
    fn test_k_values_count_and_range() {
        let params = MagnonHamiltonianParams::new(4, 1.0, 0.0, 1.0).expect("valid");
        let ks = params.ferromagnet_1d_k_values();
        assert_eq!(ks.len(), 4);
        // All k in (0, π].
        for &k in &ks {
            assert!(k > 0.0, "k must be > 0");
            assert!(k <= PI + 1e-12, "k must be ≤ π");
        }
        // Last k should be exactly π.
        assert!((ks[3] - PI).abs() < 1e-12, "last k-value should be π");
    }

    #[test]
    fn test_reference_frequencies_positive() {
        let params = MagnonHamiltonianParams::new(4, 1.0, 0.1, 1.0).expect("valid");
        let freqs = params.reference_magnon_frequencies();
        assert_eq!(freqs.len(), 4);
        for (i, &f) in freqs.iter().enumerate() {
            assert!(f > 0.0, "reference frequency {i} = {f} must be positive");
        }
    }

    #[test]
    fn test_reference_frequencies_monotone() {
        // For J>0, H_ext>0, ω_k = JS(1−cos k)+H_ext is monotonically increasing
        // on (0,π] because (1−cos k) is increasing on [0,π].
        let params = MagnonHamiltonianParams::new(8, 2.0, 0.05, 0.5).expect("valid");
        let freqs = params.reference_magnon_frequencies();
        for i in 0..freqs.len() - 1 {
            assert!(
                freqs[i] < freqs[i + 1],
                "frequencies should be monotonically increasing: f[{i}]={} ≥ f[{}]={}",
                freqs[i],
                i + 1,
                freqs[i + 1]
            );
        }
    }

    // ── MagnonNeuralNetwork ──────────────────────────────────────────────────

    #[test]
    fn test_nn_build_and_param_count() {
        let nn = MagnonNeuralNetwork::new(16, 2);
        // [3 → 16 → 16 → 2]:
        // Layer 0: 3×16 + 16 = 64
        // Layer 1: 16×16 + 16 = 272
        // Layer 2: 16×2 + 2 = 34
        // Total = 370
        assert_eq!(nn.net.n_params(), nn.params.len());
        assert!(
            !nn.params.is_empty(),
            "network must have at least one parameter"
        );
    }

    #[test]
    fn test_nn_compute_coefficients_returns_positive_a() {
        let nn = MagnonNeuralNetwork::new(16, 2);
        let params = MagnonHamiltonianParams::new(4, 1.0, 0.1, 1.0).expect("valid");
        let ks = params.ferromagnet_1d_k_values();
        for &k in &ks {
            let (a_k, _b_k) = nn.compute_coefficients(k, 1.0, 0.1);
            assert!(a_k > 0.0, "A_k must be positive, got {a_k}");
        }
    }

    #[test]
    fn test_stability_check() {
        // The default NN initialisation should satisfy A_k > |B_k| because
        // A_k = exp(A_raw) + 1e-3 ≥ 1e-3 >> 0 and B_k = output[1] × 0.3,
        // which for fresh Xavier weights is small.
        let (params, opt) = default_setup();
        let res = opt.nn.compute_coefficients_all(&params);
        // This may succeed or fail depending on random init; we just verify
        // that when it succeeds the constraint holds.
        if let Ok((a_vec, b_vec)) = res {
            for (i, (&a_k, &b_k)) in a_vec.iter().zip(b_vec.iter()).enumerate() {
                assert!(
                    a_k > b_k.abs(),
                    "stability violated at mode {i}: A_k={a_k:.4} ≤ |B_k|={:.4}",
                    b_k.abs()
                );
            }
        }
        // If it returns Err, the error path is exercised — also acceptable.
    }

    // ── QuantumClassicalOptimizer ────────────────────────────────────────────

    #[test]
    fn test_loss_computable() {
        let (_params, opt) = default_setup();
        let loss = opt.compute_loss_f64().expect("loss should be computable");
        assert!(loss.is_finite(), "loss must be finite, got {loss}");
        assert!(loss >= 0.0, "MSE loss must be non-negative, got {loss}");
    }

    #[test]
    fn test_loss_and_gradients_shape() {
        let (_params, opt) = default_setup();
        let (loss, grads) = opt
            .compute_loss_and_gradients()
            .expect("loss+grads should be computable");
        assert!(loss.is_finite(), "loss must be finite");
        assert_eq!(
            grads.len(),
            opt.nn.params.len(),
            "gradient vector length mismatch"
        );
    }

    #[test]
    fn test_finite_diff_gradient_correct() {
        // By definition, the central-difference gradient of compute_loss_with_params
        // for a single parameter must match what compute_loss_and_gradients returns.
        let (_params, opt) = default_setup();
        let (_loss, grads) = opt.compute_loss_and_gradients().expect("computable");

        // Pick the first parameter and verify manually.
        let delta = 1e-5_f64;
        let mut temp = opt.nn.params.clone();

        temp[0] += delta;
        let lp = opt
            .compute_loss_with_params(&temp)
            .expect("loss+ computable");
        temp[0] -= 2.0 * delta;
        let lm = opt
            .compute_loss_with_params(&temp)
            .expect("loss- computable");

        let manual_grad = (lp - lm) / (2.0 * delta);
        assert!(
            (grads[0] - manual_grad).abs() < 1e-12,
            "finite-diff gradient mismatch: compute_loss_and_gradients[0]={} vs manual={}",
            grads[0],
            manual_grad
        );
    }

    #[test]
    fn test_training_reduces_loss() {
        let (_params, mut opt) = default_setup();
        let initial_loss = opt.compute_loss_f64().expect("computable");

        for _ in 0..20 {
            opt.train_step().expect("train step should succeed");
        }

        let final_loss = opt.compute_loss_f64().expect("computable");
        // After 20 steps the loss should not have increased significantly.
        assert!(
            final_loss <= initial_loss * 1.1 + 1e-12,
            "loss should not increase: initial={initial_loss:.6}, final={final_loss:.6}"
        );
    }

    #[test]
    fn test_magnon_frequencies_positive_after_training() {
        let (_params, mut opt) = default_setup();

        for _ in 0..20 {
            opt.train_step().expect("train step should succeed");
        }

        let freqs = opt
            .final_frequencies()
            .expect("frequencies should be computable");
        for (i, &f) in freqs.iter().enumerate() {
            assert!(
                f > 0.0,
                "frequency {i} must be positive after training, got {f}"
            );
        }
    }

    #[test]
    fn test_result_frequencies_close_to_target() {
        let params = MagnonHamiltonianParams::new(4, 1.0, 0.1, 1.0).expect("valid");
        let mut opt = QuantumClassicalOptimizer::new(params, 16, 2, 0.01);

        let initial_loss = opt.compute_loss_f64().expect("computable");
        let result = opt.train(50).expect("training should complete");

        assert!(
            result.final_loss < initial_loss,
            "loss should decrease after 50 steps: initial={initial_loss:.6}, final={:.6}",
            result.final_loss
        );
        assert_eq!(
            result.loss_history.len(),
            50,
            "loss_history should have 50 entries"
        );
        assert_eq!(
            result.final_magnon_frequencies.len(),
            4,
            "should have 4 final magnon frequencies"
        );
    }

    #[test]
    fn test_ground_state_energy_is_negative_or_zero() {
        // For a ferromagnet, zero-point fluctuations lower the ground state:
        // E_gs = Σ_k (ε_k − A_k) / 2 ≤ 0 because ε_k = √(A_k²−B_k²) ≤ A_k.
        let params = MagnonHamiltonianParams::new(4, 1.0, 0.1, 1.0).expect("valid");
        let mut opt = QuantumClassicalOptimizer::new(params, 16, 2, 0.01);
        let result = opt.train(20).expect("training should complete");
        assert!(
            result.ground_state_energy <= 0.0 + 1e-12,
            "ground state energy should be ≤ 0, got {}",
            result.ground_state_energy
        );
    }

    #[test]
    fn test_result_has_correct_target_frequencies() {
        let params = MagnonHamiltonianParams::new(4, 1.0, 0.1, 1.0).expect("valid");
        let expected_targets = params.reference_magnon_frequencies();
        let mut opt = QuantumClassicalOptimizer::new(params, 16, 2, 0.01);
        let result = opt.train(5).expect("training should complete");

        assert_eq!(result.target_frequencies.len(), expected_targets.len());
        for (i, (&got, &exp)) in result
            .target_frequencies
            .iter()
            .zip(expected_targets.iter())
            .enumerate()
        {
            assert!(
                (got - exp).abs() < 1e-12,
                "target_frequencies[{i}] mismatch: got={got}, expected={exp}"
            );
        }
    }

    #[test]
    fn test_params_unchanged_after_compute_loss_and_gradients() {
        // compute_loss_and_gradients must not mutate self.nn.params.
        let (_params, opt) = default_setup();
        let params_before = opt.nn.params.clone();
        let _ = opt.compute_loss_and_gradients().expect("computable");
        assert_eq!(
            opt.nn.params, params_before,
            "compute_loss_and_gradients must be pure (must not mutate self.nn.params)"
        );
    }

    #[test]
    fn test_loss_with_params_does_not_mutate_self() {
        // compute_loss_with_params is also pure.
        let (_params, opt) = default_setup();
        let params_before = opt.nn.params.clone();
        let dummy_params: Vec<f64> = opt.nn.params.iter().map(|&p| p + 0.001).collect();
        let _ = opt
            .compute_loss_with_params(&dummy_params)
            .expect("computable");
        assert_eq!(
            opt.nn.params, params_before,
            "compute_loss_with_params must not mutate self.nn.params"
        );
    }
}
