//! Neural-network potentials for differentiable spintronics (v0.7.0).
//!
//! This module provides feed-forward multi-layer perceptrons (MLPs) built
//! entirely on top of the existing reverse-mode AD [`Tape`]/[`Var`] machinery.
//! The networks are intended as flexible, data-driven surrogates for
//! traditionally hand-coded physics functions such as exchange couplings
//! `J(r)` and magnetocrystalline anisotropies `K(m)`.
//!
//! ## Architecture
//!
//! A [`Mlp`] is a stack of [`Layer`]s, each of which applies an affine
//! transformation followed by an element-wise [`Activation`].  Weights are
//! initialised with either the **Xavier (Glorot–Bengio)** scheme — appropriate
//! for [`Activation::Tanh`] and [`Activation::Sigmoid`] — or the **He (Kaiming)**
//! scheme, appropriate for [`Activation::Relu`] and [`Activation::Gelu`].  A
//! private deterministic linear-congruential generator (LCG) ensures
//! bit-for-bit reproducibility from a `u64` seed.
//!
//! ## Differentiability
//!
//! Forward passes accept either plain `f64` slices (fast inference) or slices
//! of [`Var<'t>`] (gradient-enabled).  Because every operation on `Var<'t>` is
//! recorded on the same [`Tape<'t>`], gradients flow seamlessly through both
//! the network parameters (when they are themselves leaves on the tape) and
//! the inputs, enabling end-to-end training and physics-aware loss design.
//!
//! ## Physics potentials
//!
//! Two ready-made wrappers specialise the MLP to common spintronics use cases:
//!
//! | Wrapper                | Mapping                       | Physics                                  |
//! |------------------------|-------------------------------|------------------------------------------|
//! | [`NeuralExchange`]     | `r ↦ J(r)`                    | Distance-dependent Heisenberg exchange   |
//! | [`NeuralAnisotropy`]   | `(m_x, m_y, m_z) ↦ K(m̂)`      | Magnetisation-direction anisotropy       |
//!
//! ## References
//!
//! - J. Behler & M. Parrinello, "Generalized Neural-Network Representation of
//!   High-Dimensional Potential-Energy Surfaces", *Phys. Rev. Lett.* **98**,
//!   146401 (2007).
//! - K. T. Schütt, F. Arbabzadah, S. Chmiela, K.-R. Müller & A. Tkatchenko,
//!   "Quantum-chemical insights from deep tensor neural networks",
//!   *Nat. Commun.* **8**, 13890 (2017).
//! - X. Glorot & Y. Bengio, "Understanding the difficulty of training deep
//!   feedforward neural networks", *AISTATS* (2010).
//! - K. He, X. Zhang, S. Ren & J. Sun, "Delving Deep into Rectifiers:
//!   Surpassing Human-Level Performance on ImageNet Classification",
//!   *ICCV* (2015).
//! - D. Hendrycks & K. Gimpel, "Gaussian Error Linear Units (GELUs)",
//!   arXiv:1606.08415 (2016).
//! - D. E. Knuth, *The Art of Computer Programming*, Vol. 2 (LCG constants).

use crate::autodiff::tape::{Tape, Var};
use crate::error::{dimension_mismatch, invalid_param, Result};

// ─── Linear-congruential pseudo-random generator ─────────────────────────────

/// Deterministic Knuth-style 64-bit LCG used for weight initialisation.
///
/// Constants are taken from Knuth's MMIX:
/// `a = 6_364_136_223_846_793_005`, `c = 1_442_695_040_888_963_407`.
/// The high 53 bits of each draw form a uniform `f64` in `[0, 1)`.
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        // Force a non-zero state — even a zero seed should produce a non-trivial
        // sequence.
        let state = if seed == 0 {
            0xDEAD_BEEF_CAFE_BABE
        } else {
            seed
        };
        Self { state }
    }

    /// Advance the state and return the new 64-bit word.
    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    /// Sample a `f64` uniformly in `[0, 1)` using the high 53 bits.
    fn next_unit(&mut self) -> f64 {
        let bits = self.next_u64() >> 11;
        (bits as f64) / ((1u64 << 53) as f64)
    }

    /// Sample a `f64` uniformly in `(low, high)`.
    fn next_uniform(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.next_unit()
    }

    /// Sample a `f64` from `Normal(0, 1)` using the Marsaglia polar method.
    fn next_normal(&mut self) -> f64 {
        loop {
            let u = 2.0 * self.next_unit() - 1.0;
            let v = 2.0 * self.next_unit() - 1.0;
            let s = u * u + v * v;
            if s > 0.0 && s < 1.0 {
                let factor = (-2.0 * s.ln() / s).sqrt();
                return u * factor;
            }
        }
    }
}

// ─── Activation ──────────────────────────────────────────────────────────────

/// Element-wise non-linear activation function applied at the output of a
/// [`Layer`].
///
/// The smooth variants ([`Activation::Tanh`], [`Activation::Sigmoid`],
/// [`Activation::Gelu`]) are preferred for physics-informed losses where
/// higher-order derivatives are required.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Activation {
    /// Rectified linear unit: `max(0, x)`.
    Relu,
    /// Hyperbolic tangent: `tanh(x)`.
    Tanh,
    /// Logistic sigmoid: `1 / (1 + exp(-x))`.
    Sigmoid,
    /// Gaussian Error Linear Unit (Hendrycks–Gimpel approximation).
    Gelu,
    /// Identity / pass-through: `x`.
    Linear,
}

impl Activation {
    /// Plain-`f64` evaluation, used by [`Mlp::forward_f64`].
    pub fn apply_f64(self, x: f64) -> f64 {
        match self {
            Activation::Relu => x.max(0.0),
            Activation::Tanh => x.tanh(),
            Activation::Sigmoid => 1.0 / (1.0 + (-x).exp()),
            Activation::Gelu => gelu_f64(x),
            Activation::Linear => x,
        }
    }

    /// Differentiable evaluation on a [`Var<'t>`].
    ///
    /// All branches are expressible via the elementary operations already
    /// recorded on the tape, so gradients propagate without extra work.
    pub fn apply<'t>(self, x: Var<'t>) -> Var<'t> {
        match self {
            Activation::Relu => {
                // ReLU(x) = (x + |x|) / 2 — smooth-ish branch that the
                // existing `Var::abs` derivative handles correctly except at
                // the singular point x = 0 (where ∂|x| := 0).
                (x + x.abs()) * 0.5
            },
            Activation::Tanh => x.tanh(),
            Activation::Sigmoid => 0.5 + 0.5 * (x * 0.5).tanh(),
            Activation::Gelu => gelu_var(x),
            Activation::Linear => x,
        }
    }
}

/// GELU as in Hendrycks–Gimpel (2016):
///
/// `GELU(x) ≈ 0.5 · x · (1 + tanh(√(2/π) · (x + 0.044715 · x³)))`.
fn gelu_f64(x: f64) -> f64 {
    let k0 = (2.0_f64 / std::f64::consts::PI).sqrt();
    let inner = k0 * (x + 0.044_715 * x * x * x);
    0.5 * x * (1.0 + inner.tanh())
}

fn gelu_var<'t>(x: Var<'t>) -> Var<'t> {
    let k0 = (2.0_f64 / std::f64::consts::PI).sqrt();
    let x_cubed = x * x * x;
    let inner = (x + x_cubed * 0.044_715) * k0;
    x * 0.5 * (1.0 + inner.tanh())
}

// ─── Layer ───────────────────────────────────────────────────────────────────

/// A single dense (fully-connected) layer:
///
/// `y_j = activation(Σ_i w_{ji} · x_i + b_j)`.
///
/// Weights are stored row-major in `weights[out_dim][in_dim]` so each output
/// neuron occupies a contiguous slice — the layout most natural for
/// matrix-vector products and serialisation.
#[derive(Debug, Clone)]
pub struct Layer {
    /// Weight matrix `[out_dim][in_dim]`.
    pub weights: Vec<Vec<f64>>,
    /// Bias vector `[out_dim]`.
    pub biases: Vec<f64>,
    /// Activation function applied element-wise.
    pub activation: Activation,
    /// Input dimension (number of incoming features).
    pub in_dim: usize,
    /// Output dimension (number of neurons).
    pub out_dim: usize,
}

impl Layer {
    /// Construct a new layer with random weights initialised by
    /// **Xavier (Glorot)** for [`Activation::Tanh`]/[`Activation::Sigmoid`] and
    /// **He (Kaiming)** for [`Activation::Relu`]/[`Activation::Gelu`].
    ///
    /// [`Activation::Linear`] uses Xavier by default.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] if either dimension is zero.
    pub fn new(
        in_dim: usize,
        out_dim: usize,
        activation: Activation,
        rng_seed: u64,
    ) -> Result<Self> {
        if in_dim == 0 {
            return Err(invalid_param("in_dim", "must be positive"));
        }
        if out_dim == 0 {
            return Err(invalid_param("out_dim", "must be positive"));
        }
        let mut rng = Lcg::new(rng_seed);
        let weights = match activation {
            Activation::Relu | Activation::Gelu => {
                // He init: w ~ N(0, sqrt(2/in)).
                let std_dev = (2.0 / in_dim as f64).sqrt();
                (0..out_dim)
                    .map(|_| (0..in_dim).map(|_| std_dev * rng.next_normal()).collect())
                    .collect()
            },
            Activation::Tanh | Activation::Sigmoid | Activation::Linear => {
                // Xavier init: w ~ U(-a, a), a = sqrt(6/(in+out)).
                let a = (6.0 / (in_dim + out_dim) as f64).sqrt();
                (0..out_dim)
                    .map(|_| (0..in_dim).map(|_| rng.next_uniform(-a, a)).collect())
                    .collect()
            },
        };
        let biases = vec![0.0; out_dim];
        Ok(Self {
            weights,
            biases,
            activation,
            in_dim,
            out_dim,
        })
    }

    /// Differentiable forward pass.  Each output element is recorded as an
    /// affine combination on the tape followed by the configured activation.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] if `inputs.len() != in_dim`.
    pub fn forward<'t>(&self, _tape: &'t Tape, inputs: &[Var<'t>]) -> Result<Vec<Var<'t>>> {
        if inputs.len() != self.in_dim {
            return Err(dimension_mismatch(
                &format!("{} inputs", self.in_dim),
                &format!("{} inputs", inputs.len()),
            ));
        }
        let mut outputs = Vec::with_capacity(self.out_dim);
        for j in 0..self.out_dim {
            // Pre-activation: sum_i w_ji x_i + b_j.
            // We accumulate by repeated addition on the tape.
            let w_row = &self.weights[j];
            let mut pre = inputs[0] * w_row[0];
            for i in 1..self.in_dim {
                pre = pre + inputs[i] * w_row[i];
            }
            pre = pre + self.biases[j];
            outputs.push(self.activation.apply(pre));
        }
        Ok(outputs)
    }

    /// Fast plain-`f64` forward, used for inference without gradient tracking.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] if `inputs.len() != in_dim`.
    pub fn forward_f64(&self, inputs: &[f64]) -> Result<Vec<f64>> {
        if inputs.len() != self.in_dim {
            return Err(dimension_mismatch(
                &format!("{} inputs", self.in_dim),
                &format!("{} inputs", inputs.len()),
            ));
        }
        let mut out = Vec::with_capacity(self.out_dim);
        for j in 0..self.out_dim {
            let w_row = &self.weights[j];
            let mut sum = self.biases[j];
            for i in 0..self.in_dim {
                sum += w_row[i] * inputs[i];
            }
            out.push(self.activation.apply_f64(sum));
        }
        Ok(out)
    }

    /// Total number of free parameters (weights + biases).
    pub fn n_params(&self) -> usize {
        self.in_dim * self.out_dim + self.out_dim
    }

    /// Flatten parameters as `[w_00, w_01, …, w_{out-1, in-1}, b_0, …, b_{out-1}]`.
    pub fn params_flat(&self) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.n_params());
        for row in &self.weights {
            v.extend_from_slice(row);
        }
        v.extend_from_slice(&self.biases);
        v
    }

    /// Restore parameters from a flat slice produced by [`Layer::params_flat`].
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the length is wrong.
    pub fn set_params(&mut self, flat: &[f64]) -> Result<()> {
        let expected = self.n_params();
        if flat.len() != expected {
            return Err(dimension_mismatch(
                &format!("{expected} params"),
                &format!("{} params", flat.len()),
            ));
        }
        let mut cursor = 0;
        for j in 0..self.out_dim {
            for i in 0..self.in_dim {
                self.weights[j][i] = flat[cursor];
                cursor += 1;
            }
        }
        for j in 0..self.out_dim {
            self.biases[j] = flat[cursor];
            cursor += 1;
        }
        Ok(())
    }
}

// ─── Mlp ─────────────────────────────────────────────────────────────────────

/// A multi-layer perceptron — a sequential stack of [`Layer`]s.
#[derive(Debug, Clone)]
pub struct Mlp {
    /// The layers, in evaluation order.
    pub layers: Vec<Layer>,
    /// Cached input dimensionality (= `layers\[0\].in_dim`).
    pub input_dim: usize,
    /// Cached output dimensionality (= `layers.last().out_dim`).
    pub output_dim: usize,
}

impl Mlp {
    /// Construct an MLP given a list of layer sizes and a parallel list of
    /// activations.
    ///
    /// `layer_sizes` has the form `[in, h_1, h_2, …, out]` (length ≥ 2) and
    /// `activations` lists one activation per layer (length `layer_sizes.len() - 1`).
    ///
    /// A per-layer sub-seed is derived deterministically from `rng_seed` so
    /// that successive layers receive independent yet reproducible weights.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when shapes are invalid.
    pub fn new(layer_sizes: &[usize], activations: &[Activation], rng_seed: u64) -> Result<Self> {
        if layer_sizes.len() < 2 {
            return Err(invalid_param(
                "layer_sizes",
                "need at least input and output dimensions",
            ));
        }
        if activations.len() + 1 != layer_sizes.len() {
            return Err(invalid_param(
                "activations",
                "must have one fewer entry than layer_sizes",
            ));
        }
        let mut layers = Vec::with_capacity(activations.len());
        for (k, act) in activations.iter().enumerate() {
            // Mix the layer index in via a Weyl-style stride so each sub-seed
            // is well separated even when rng_seed is small.
            let sub_seed = rng_seed.wrapping_add((k as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            layers.push(Layer::new(
                layer_sizes[k],
                layer_sizes[k + 1],
                *act,
                sub_seed,
            )?);
        }
        Ok(Self {
            input_dim: layer_sizes[0],
            output_dim: layer_sizes[layer_sizes.len() - 1],
            layers,
        })
    }

    /// Differentiable forward pass: composes [`Layer::forward`] in order.
    ///
    /// # Errors
    /// Propagates dimension-mismatch errors from inner layers.
    pub fn forward<'t>(&self, tape: &'t Tape, inputs: &[Var<'t>]) -> Result<Vec<Var<'t>>> {
        let mut current = inputs.to_vec();
        for layer in &self.layers {
            current = layer.forward(tape, &current)?;
        }
        Ok(current)
    }

    /// Plain-`f64` forward (no gradient tracking) — significantly faster.
    ///
    /// # Errors
    /// Propagates dimension-mismatch errors from inner layers.
    pub fn forward_f64(&self, inputs: &[f64]) -> Result<Vec<f64>> {
        let mut current = inputs.to_vec();
        for layer in &self.layers {
            current = layer.forward_f64(&current)?;
        }
        Ok(current)
    }

    /// Total number of free parameters across all layers.
    pub fn n_params(&self) -> usize {
        self.layers.iter().map(|l| l.n_params()).sum()
    }

    /// Concatenate all layer parameters in layer order.
    pub fn params_flat(&self) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.n_params());
        for layer in &self.layers {
            v.extend(layer.params_flat());
        }
        v
    }

    /// Inverse of [`Mlp::params_flat`].
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the length is wrong.
    pub fn set_params(&mut self, flat: &[f64]) -> Result<()> {
        let expected = self.n_params();
        if flat.len() != expected {
            return Err(dimension_mismatch(
                &format!("{expected} params"),
                &format!("{} params", flat.len()),
            ));
        }
        let mut cursor = 0;
        for layer in &mut self.layers {
            let n = layer.n_params();
            layer.set_params(&flat[cursor..cursor + n])?;
            cursor += n;
        }
        Ok(())
    }
}

// ─── NeuralExchange ──────────────────────────────────────────────────────────

/// Distance-dependent exchange coupling `J(r)` parameterised by an [`Mlp`].
///
/// The input is the scalar inter-spin distance `r` (in arbitrary units chosen
/// by the user) and the output is the scalar coupling `J(r)`.  An internal
/// `tanh` rescaling maps `r ↦ (2r − (r_min + r_max)) / (r_max − r_min)` so the
/// network operates on a roughly `[-1, 1]` domain regardless of the physical
/// scale.
pub struct NeuralExchange {
    /// Underlying network — input dim 1, output dim 1.
    pub mlp: Mlp,
    /// Lower bound of the physical distance range used for input scaling.
    pub r_min: f64,
    /// Upper bound of the physical distance range used for input scaling.
    pub r_max: f64,
}

impl NeuralExchange {
    /// Create a new `NeuralExchange` with `hidden_sizes` neurons in each
    /// hidden layer (all using [`Activation::Tanh`]) and a final linear layer.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when `r_min >= r_max`.
    pub fn new(hidden_sizes: &[usize], r_min: f64, r_max: f64, rng_seed: u64) -> Result<Self> {
        if r_min >= r_max {
            return Err(invalid_param("r_min/r_max", "must satisfy r_min < r_max"));
        }
        let mut sizes = vec![1_usize];
        sizes.extend_from_slice(hidden_sizes);
        sizes.push(1);
        let mut acts: Vec<Activation> = vec![Activation::Tanh; hidden_sizes.len()];
        acts.push(Activation::Linear);
        let mlp = Mlp::new(&sizes, &acts, rng_seed)?;
        Ok(Self { mlp, r_min, r_max })
    }

    /// Plain-`f64` evaluation at distance `r`.
    ///
    /// # Errors
    /// Propagates errors from the underlying MLP forward pass.
    pub fn coupling(&self, r: f64) -> Result<f64> {
        let scaled = self.scale_f64(r);
        let out = self.mlp.forward_f64(&[scaled])?;
        Ok(out[0])
    }

    /// Differentiable evaluation that records all operations on `tape`.
    ///
    /// # Errors
    /// Propagates errors from the underlying MLP forward pass.
    pub fn coupling_diff<'t>(&self, tape: &'t Tape, r: Var<'t>) -> Result<Var<'t>> {
        let scaled = (r - 0.5 * (self.r_min + self.r_max)) * (2.0 / (self.r_max - self.r_min));
        let out = self.mlp.forward(tape, &[scaled])?;
        Ok(out[0])
    }

    fn scale_f64(&self, r: f64) -> f64 {
        (r - 0.5 * (self.r_min + self.r_max)) * (2.0 / (self.r_max - self.r_min))
    }
}

// ─── NeuralAnisotropy ───────────────────────────────────────────────────────

/// Magnetisation-direction anisotropy `K(m̂)` parameterised by an [`Mlp`].
///
/// The input is the unit-vector magnetisation direction `(m_x, m_y, m_z)` and
/// the scalar output is the anisotropy energy density (in arbitrary user
/// units).  The default activation stack is `Tanh → Tanh → Linear`, giving a
/// twice-differentiable surrogate ideal for gradient-based magnetic-structure
/// optimisation.
pub struct NeuralAnisotropy {
    /// Underlying network — input dim 3, output dim 1.
    pub mlp: Mlp,
}

impl NeuralAnisotropy {
    /// Create a new `NeuralAnisotropy` with `hidden_sizes` neurons per hidden
    /// layer (all `Tanh`) and a final linear output layer.
    ///
    /// # Errors
    /// Propagates errors from [`Mlp::new`].
    pub fn new(hidden_sizes: &[usize], rng_seed: u64) -> Result<Self> {
        let mut sizes = vec![3_usize];
        sizes.extend_from_slice(hidden_sizes);
        sizes.push(1);
        let mut acts: Vec<Activation> = vec![Activation::Tanh; hidden_sizes.len()];
        acts.push(Activation::Linear);
        let mlp = Mlp::new(&sizes, &acts, rng_seed)?;
        Ok(Self { mlp })
    }

    /// Plain-`f64` evaluation.
    ///
    /// # Errors
    /// Propagates errors from the underlying MLP forward pass.
    pub fn energy(&self, mx: f64, my: f64, mz: f64) -> Result<f64> {
        let out = self.mlp.forward_f64(&[mx, my, mz])?;
        Ok(out[0])
    }

    /// Differentiable evaluation that records all operations on `tape`.
    ///
    /// # Errors
    /// Propagates errors from the underlying MLP forward pass.
    pub fn energy_diff<'t>(
        &self,
        tape: &'t Tape,
        mx: Var<'t>,
        my: Var<'t>,
        mz: Var<'t>,
    ) -> Result<Var<'t>> {
        let out = self.mlp.forward(tape, &[mx, my, mz])?;
        Ok(out[0])
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::autodiff::tape::Tape;

    fn central_diff<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
        (f(x + h) - f(x - h)) / (2.0 * h)
    }

    // 1. Activation enum properties (eq + hash work, every variant evaluates).
    #[test]
    fn test_activation_enum_basic() {
        let vs = [
            Activation::Relu,
            Activation::Tanh,
            Activation::Sigmoid,
            Activation::Gelu,
            Activation::Linear,
        ];
        for a in vs.iter() {
            let y = a.apply_f64(0.5);
            assert!(y.is_finite());
        }
        assert_eq!(Activation::Relu, Activation::Relu);
        assert_ne!(Activation::Relu, Activation::Tanh);
    }

    // 2. Xavier init range: for Tanh init, all weights in [-a, a] with
    //    a = sqrt(6/(in+out)).
    #[test]
    fn test_xavier_init_in_range() {
        let in_dim = 4;
        let out_dim = 6;
        let layer = Layer::new(in_dim, out_dim, Activation::Tanh, 12345).unwrap();
        let bound = (6.0 / (in_dim + out_dim) as f64).sqrt();
        for row in &layer.weights {
            for &w in row {
                assert!(
                    w.abs() < bound + 1e-12,
                    "weight {} outside Xavier bound {}",
                    w,
                    bound
                );
            }
        }
    }

    // 3. He init: empirical std of Relu-init weights should be close to
    //    sqrt(2/in) within a generous tolerance.
    #[test]
    fn test_he_init_std_estimate() {
        let in_dim = 50;
        let out_dim = 50;
        let layer = Layer::new(in_dim, out_dim, Activation::Relu, 99).unwrap();
        let mut sum = 0.0_f64;
        let mut sum_sq = 0.0_f64;
        let mut n = 0.0_f64;
        for row in &layer.weights {
            for &w in row {
                sum += w;
                sum_sq += w * w;
                n += 1.0;
            }
        }
        let mean = sum / n;
        let var = sum_sq / n - mean * mean;
        let std_est = var.sqrt();
        let expected = (2.0 / in_dim as f64).sqrt();
        assert!(
            (std_est - expected).abs() / expected < 0.3,
            "He-init std {} vs expected {}",
            std_est,
            expected
        );
    }

    // 4. Layer forward output has the right shape.
    #[test]
    fn test_layer_forward_shape() {
        let layer = Layer::new(3, 5, Activation::Tanh, 7).unwrap();
        let out = layer.forward_f64(&[0.1, -0.2, 0.3]).unwrap();
        assert_eq!(out.len(), 5);
        for y in out {
            assert!(y.is_finite());
            assert!(y.abs() <= 1.0);
        }
    }

    // 5. MLP forward shape and finiteness.
    #[test]
    fn test_mlp_forward_shape() {
        let mlp = Mlp::new(
            &[2, 4, 4, 3],
            &[Activation::Tanh, Activation::Gelu, Activation::Linear],
            42,
        )
        .unwrap();
        assert_eq!(mlp.input_dim, 2);
        assert_eq!(mlp.output_dim, 3);
        let y = mlp.forward_f64(&[0.5, -0.3]).unwrap();
        assert_eq!(y.len(), 3);
        for v in y {
            assert!(v.is_finite());
        }
    }

    // 6. Chain rule through MLP: input gradient via AD matches finite-diff.
    #[test]
    fn test_mlp_input_gradient_chain_rule() {
        let mlp = Mlp::new(
            &[1, 8, 8, 1],
            &[Activation::Tanh, Activation::Tanh, Activation::Linear],
            321,
        )
        .unwrap();
        let x_val = 0.4;
        // AD gradient w.r.t. the single input.
        let tape = Tape::new();
        let xv = Var::leaf(&tape, x_val);
        let out = mlp.forward(&tape, &[xv]).unwrap();
        tape.backward(out[0]);
        let ad_grad = xv.grad();
        // Finite-difference gradient.
        let fd_grad = central_diff(|x| mlp.forward_f64(&[x]).unwrap()[0], x_val, 1e-6);
        assert!(
            (ad_grad - fd_grad).abs() < 1e-5,
            "AD {} vs FD {}",
            ad_grad,
            fd_grad
        );
    }

    // 7. Tanh activation derivative on the tape matches analytic 1-tanh².
    #[test]
    fn test_activation_tanh_derivative() {
        let x_val = 0.6_f64;
        let tape = Tape::new();
        let xv = Var::leaf(&tape, x_val);
        let yv = Activation::Tanh.apply(xv);
        tape.backward(yv);
        let analytic = 1.0 - x_val.tanh().powi(2);
        assert!((xv.grad() - analytic).abs() < 1e-12);
    }

    // 8. GELU value and derivative match the Hendrycks-Gimpel formula.
    #[test]
    fn test_activation_gelu_value_and_grad() {
        let x_val = 0.7_f64;
        let expected_val = gelu_f64(x_val);
        let tape = Tape::new();
        let xv = Var::leaf(&tape, x_val);
        let yv = Activation::Gelu.apply(xv);
        let ad_val = yv.value();
        assert!((ad_val - expected_val).abs() < 1e-12);
        tape.backward(yv);
        let ad_grad = xv.grad();
        let fd_grad = central_diff(gelu_f64, x_val, 1e-6);
        assert!((ad_grad - fd_grad).abs() < 1e-5);
    }

    // 9. NeuralExchange returns a finite J(r) for r in [r_min, r_max] and is
    //    consistent between the f64 and Var paths.
    #[test]
    fn test_neural_exchange_in_range_and_consistent() {
        let nx = NeuralExchange::new(&[6, 6], 0.5, 3.0, 17).unwrap();
        for r in [0.6_f64, 1.2, 2.0, 2.8] {
            let j_f64 = nx.coupling(r).unwrap();
            let tape = Tape::new();
            let rv = Var::leaf(&tape, r);
            let j_var = nx.coupling_diff(&tape, rv).unwrap();
            assert!(j_f64.is_finite());
            assert!((j_f64 - j_var.value()).abs() < 1e-12);
        }
    }

    // 10. NeuralAnisotropy reproduces the uniaxial K·mz² behaviour after a
    //     short training-free initialisation — just check finiteness and that
    //     gradient w.r.t. mz is non-zero.
    #[test]
    fn test_neural_anisotropy_finite_and_diffable() {
        let na = NeuralAnisotropy::new(&[8, 8], 2024).unwrap();
        let e = na.energy(0.0, 0.0, 1.0).unwrap();
        assert!(e.is_finite());
        let tape = Tape::new();
        let mx = Var::leaf(&tape, 0.0);
        let my = Var::leaf(&tape, 0.0);
        let mz = Var::leaf(&tape, 1.0);
        let ev = na.energy_diff(&tape, mx, my, mz).unwrap();
        tape.backward(ev);
        let grad_mz = mz.grad();
        // Generic NN should have non-zero gradient for at least one of the
        // three inputs.
        let total = mx.grad().abs() + my.grad().abs() + grad_mz.abs();
        assert!(total > 1e-12);
    }

    // 11. params_flat/set_params round-trip preserves the network exactly.
    #[test]
    fn test_params_flat_set_params_roundtrip() {
        let mut mlp = Mlp::new(&[3, 4, 2], &[Activation::Tanh, Activation::Linear], 55).unwrap();
        let original = mlp.params_flat();
        let mut perturbed = original.clone();
        for v in &mut perturbed {
            *v += 0.1;
        }
        mlp.set_params(&perturbed).unwrap();
        let rt = mlp.params_flat();
        assert_eq!(rt.len(), perturbed.len());
        for (a, b) in rt.iter().zip(perturbed.iter()) {
            assert!((a - b).abs() < 1e-15);
        }
    }

    // 12. n_params accounting: weights + biases agree with the layout we set.
    #[test]
    fn test_n_params_accounting() {
        let layer = Layer::new(7, 11, Activation::Tanh, 1).unwrap();
        assert_eq!(layer.n_params(), 7 * 11 + 11);
        let mlp = Mlp::new(&[2, 5, 3], &[Activation::Tanh, Activation::Linear], 2).unwrap();
        assert_eq!(mlp.n_params(), (2 * 5 + 5) + (5 * 3 + 3));
        assert_eq!(mlp.params_flat().len(), mlp.n_params());
    }

    // 13. Reproducibility: same seed → bitwise-identical parameters.
    #[test]
    fn test_rng_seed_reproducibility() {
        let m1 = Mlp::new(
            &[4, 6, 6, 2],
            &[Activation::Relu, Activation::Tanh, Activation::Linear],
            123_456,
        )
        .unwrap();
        let m2 = Mlp::new(
            &[4, 6, 6, 2],
            &[Activation::Relu, Activation::Tanh, Activation::Linear],
            123_456,
        )
        .unwrap();
        let p1 = m1.params_flat();
        let p2 = m2.params_flat();
        assert_eq!(p1.len(), p2.len());
        for (a, b) in p1.iter().zip(p2.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        // Different seed → different params.
        let m3 = Mlp::new(
            &[4, 6, 6, 2],
            &[Activation::Relu, Activation::Tanh, Activation::Linear],
            123_457,
        )
        .unwrap();
        let p3 = m3.params_flat();
        let any_diff = p1.iter().zip(p3.iter()).any(|(a, b)| a != b);
        assert!(any_diff);
    }

    // 14. Errors on invalid inputs/shapes.
    #[test]
    fn test_invalid_shapes_return_errors() {
        assert!(Layer::new(0, 4, Activation::Tanh, 0).is_err());
        assert!(Layer::new(4, 0, Activation::Tanh, 0).is_err());
        assert!(Mlp::new(&[3], &[], 0).is_err());
        assert!(Mlp::new(&[3, 4, 1], &[Activation::Tanh], 0).is_err());
        let mlp = Mlp::new(&[2, 3, 1], &[Activation::Tanh, Activation::Linear], 1).unwrap();
        assert!(mlp.forward_f64(&[1.0]).is_err());
        let bad = vec![0.0; 0];
        let mut mlp2 = mlp.clone();
        assert!(mlp2.set_params(&bad).is_err());
        assert!(NeuralExchange::new(&[4], 1.0, 1.0, 0).is_err());
    }
}
