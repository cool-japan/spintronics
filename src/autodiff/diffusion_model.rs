//! # DDPM Diffusion Model for Spin Texture Generation (v0.9.0)
//!
//! This module implements a **Denoising Diffusion Probabilistic Model** (DDPM)
//! for generating physically realistic spin texture configurations such as
//! skyrmion lattices, helical spirals, and uniform ferromagnetic states.
//!
//! ## Physics Background
//!
//! A spin texture is a 2D magnetisation field `m(r) : R² → S²` where each
//! 3-vector satisfies `|m(i,j)| = 1`.  The DDPM learns the distribution of
//! such textures from examples by:
//!
//! 1. **Forward diffusion** — gradually corrupting a clean texture `x₀` by
//!    adding Gaussian noise according to a linear variance schedule.
//! 2. **Denoising network** — a 2-layer MLP `ε_θ(x_t, t)` that predicts the
//!    noise added at each step.
//! 3. **Reverse sampling** — starting from pure Gaussian noise `x_T`, iterate
//!    the denoising step for `T` steps to recover a plausible spin texture.
//!
//! After sampling, each 3-vector is normalised back onto the unit sphere.
//!
//! ## References
//!
//! - J. Ho, A. Jain & P. Abbeel, "Denoising Diffusion Probabilistic Models",
//!   *NeurIPS* **33**, 6840–6851 (2020).
//! - A. Q. Nichol & P. Dhariwal, "Improved Denoising Diffusion Probabilistic
//!   Models", *ICML* (2021), arXiv:2102.09672.
//! - N. Nagaosa & Y. Tokura, "Topological properties and dynamics of magnetic
//!   skyrmions", *Nat. Nanotechnol.* **8**, 899 (2013).

use std::f64::consts::PI;

use crate::error::{invalid_param, numerical_error, Result};
use crate::texture::skyrmion::SkyrmionLattice;
use crate::texture::topology::calculate_skyrmion_number;
use crate::vector3::Vector3;

// ─── Linear Congruential Generator ───────────────────────────────────────────

/// Minimal Knuth-MMIX 64-bit LCG for reproducible pseudo-random sequences.
///
/// Constants: multiplier `a = 6_364_136_223_846_793_005`,
/// addend `c = 1_442_695_040_888_963_407`.
/// Only used internally; not cryptographically secure.
#[derive(Debug, Clone)]
pub struct Lcg {
    state: u64,
}

impl Lcg {
    /// Seed the generator.  State is forced non-zero.
    pub fn new(seed: u64) -> Self {
        Lcg {
            state: seed.wrapping_add(1),
        }
    }

    /// Advance the state and return the next 64-bit word.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    /// Sample a `f64` uniformly in `[0, 1)` using the high 53 bits.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Sample a standard normal `N(0,1)` via the Box–Muller transform.
    pub fn next_normal(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-15);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
}

// ─── NoiseSchedule ────────────────────────────────────────────────────────────

/// DDPM variance schedule: linear interpolation of `β` from `β_start` to
/// `β_end` over `T` steps, with derived cumulative-product `ᾱ`.
///
/// Index convention: all slices are **0-indexed** (`beta[0]` = β₁).
/// The accessor methods use the **1-indexed** convention from the DDPM paper.
#[derive(Debug, Clone)]
pub struct NoiseSchedule {
    /// Total number of diffusion steps T.
    pub t_max: usize,
    /// β_t for t = 1..T (0-indexed storage).
    pub beta: Vec<f64>,
    /// α_t = 1 − β_t (0-indexed storage).
    pub alpha: Vec<f64>,
    /// ᾱ_t = ∏_{s=1}^{t} α_s (0-indexed storage).
    pub alpha_bar: Vec<f64>,
}

impl NoiseSchedule {
    /// Construct a linear variance schedule.
    ///
    /// `β_t` is linearly spaced from `beta_start` (at t = 1) to `beta_end`
    /// (at t = T).
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when `t_max == 0` or
    /// bounds are outside `(0, 1)`.
    pub fn linear(t_max: usize, beta_start: f64, beta_end: f64) -> Result<Self> {
        if t_max == 0 {
            return Err(invalid_param("t_max", "must be at least 1"));
        }
        if beta_start <= 0.0 || beta_start >= 1.0 {
            return Err(invalid_param("beta_start", "must be in (0, 1)"));
        }
        if beta_end <= 0.0 || beta_end >= 1.0 {
            return Err(invalid_param("beta_end", "must be in (0, 1)"));
        }
        if beta_start >= beta_end {
            return Err(invalid_param("beta_start", "must be less than beta_end"));
        }

        let mut beta = Vec::with_capacity(t_max);
        let mut alpha = Vec::with_capacity(t_max);
        let mut alpha_bar = Vec::with_capacity(t_max);

        let mut cumulative_alpha = 1.0_f64;
        for i in 0..t_max {
            let bt = if t_max == 1 {
                beta_start
            } else {
                beta_start + (beta_end - beta_start) * (i as f64) / ((t_max - 1) as f64)
            };
            let at = 1.0 - bt;
            cumulative_alpha *= at;
            beta.push(bt);
            alpha.push(at);
            alpha_bar.push(cumulative_alpha);
        }

        Ok(Self {
            t_max,
            beta,
            alpha,
            alpha_bar,
        })
    }

    /// β at timestep t (1-indexed, t ∈ [1, T]).
    ///
    /// # Panics
    /// Panics in debug mode if `t == 0` or `t > t_max`.
    pub fn beta_t(&self, t: usize) -> f64 {
        debug_assert!(t >= 1 && t <= self.t_max, "t out of range [1, T]");
        self.beta[t - 1]
    }

    /// ᾱ at timestep t (1-indexed, t ∈ [1, T]).
    ///
    /// # Panics
    /// Panics in debug mode if `t == 0` or `t > t_max`.
    pub fn alpha_bar_t(&self, t: usize) -> f64 {
        debug_assert!(t >= 1 && t <= self.t_max, "t out of range [1, T]");
        self.alpha_bar[t - 1]
    }
}

// ─── SpinTexture ──────────────────────────────────────────────────────────────

/// A discretised 2D spin texture: a flat vector of magnetisation components
/// for an `nx × ny` grid.
///
/// **Layout**: components for grid site `(i, j)` start at index
/// `(j * nx + i) * 3`, storing `[m_x, m_y, m_z]` in that order.
/// Row-major order with `i` (x-index) running fastest.
#[derive(Debug, Clone)]
pub struct SpinTexture {
    /// Flat component buffer: `3 * nx * ny` elements.
    pub data: Vec<f64>,
    /// Number of grid points along x.
    pub nx: usize,
    /// Number of grid points along y.
    pub ny: usize,
}

impl SpinTexture {
    /// Sample a [`SkyrmionLattice`] on an `nx × ny` grid spanning
    /// `[0, domain_size] × [0, domain_size]`.
    ///
    /// Every 3-vector is normalised to unit length after sampling.
    ///
    /// # Errors
    /// Returns an error if `nx == 0` or `ny == 0`.
    pub fn from_skyrmion_lattice(
        lattice: &SkyrmionLattice,
        nx: usize,
        ny: usize,
        domain_size: f64,
        wall_width: f64,
    ) -> Result<Self> {
        if nx == 0 {
            return Err(invalid_param("nx", "must be at least 1"));
        }
        if ny == 0 {
            return Err(invalid_param("ny", "must be at least 1"));
        }

        let n_comp = 3 * nx * ny;
        let mut data = vec![0.0_f64; n_comp];

        let dx_step = if nx > 1 {
            domain_size / (nx - 1) as f64
        } else {
            0.0
        };
        let dy_step = if ny > 1 {
            domain_size / (ny - 1) as f64
        } else {
            0.0
        };

        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * dx_step;
                let y = j as f64 * dy_step;

                // Sum contributions from all skyrmions; closest one dominates.
                let mut m = Vector3::new(0.0_f64, 0.0, 1.0);
                if !lattice.skyrmions.is_empty() {
                    // Find the nearest skyrmion and use its magnetization,
                    // modulated by a distance weight to blend boundaries.
                    let mut best_dist = f64::INFINITY;
                    let mut best_idx = 0_usize;
                    for (k, sk) in lattice.skyrmions.iter().enumerate() {
                        let ddx = x - sk.center.0;
                        let ddy = y - sk.center.1;
                        let dist = (ddx * ddx + ddy * ddy).sqrt();
                        if dist < best_dist {
                            best_dist = dist;
                            best_idx = k;
                        }
                    }
                    m = lattice.skyrmions[best_idx].magnetization_at(x, y, wall_width);
                }

                // Normalise to enforce |m| = 1.
                let mag = (m.x * m.x + m.y * m.y + m.z * m.z).sqrt();
                let (mx, my, mz) = if mag > 1e-15 {
                    (m.x / mag, m.y / mag, m.z / mag)
                } else {
                    (0.0, 0.0, 1.0)
                };

                let base = (j * nx + i) * 3;
                data[base] = mx;
                data[base + 1] = my;
                data[base + 2] = mz;
            }
        }

        Ok(Self { data, nx, ny })
    }

    /// Construct a uniform ferromagnetic texture with all spins pointing along
    /// `direction` (normalised internally).
    ///
    /// # Errors
    /// Returns an error if `direction` is the zero vector or dimensions are zero.
    pub fn from_uniform(nx: usize, ny: usize, direction: Vector3<f64>) -> Result<Self> {
        if nx == 0 {
            return Err(invalid_param("nx", "must be at least 1"));
        }
        if ny == 0 {
            return Err(invalid_param("ny", "must be at least 1"));
        }

        let mag =
            (direction.x * direction.x + direction.y * direction.y + direction.z * direction.z)
                .sqrt();
        if mag < 1e-15 {
            return Err(invalid_param("direction", "must be a non-zero vector"));
        }
        let (dx, dy, dz) = (direction.x / mag, direction.y / mag, direction.z / mag);

        let n_comp = 3 * nx * ny;
        let mut data = vec![0.0_f64; n_comp];
        for k in 0..nx * ny {
            data[k * 3] = dx;
            data[k * 3 + 1] = dy;
            data[k * 3 + 2] = dz;
        }

        Ok(Self { data, nx, ny })
    }

    /// Number of lattice sites (`nx * ny`).
    pub fn n_sites(&self) -> usize {
        self.nx * self.ny
    }

    /// Total number of scalar components (`3 * nx * ny`).
    pub fn n_components(&self) -> usize {
        3 * self.nx * self.ny
    }

    /// Retrieve the spin 3-vector at grid site `(i, j)`.
    ///
    /// Index formula: `base = (j * nx + i) * 3`.
    pub fn get_spin(&self, i: usize, j: usize) -> Vector3<f64> {
        let base = (j * self.nx + i) * 3;
        Vector3::new(self.data[base], self.data[base + 1], self.data[base + 2])
    }

    /// Add Gaussian noise scaled by the DDPM forward-diffusion coefficient.
    ///
    /// Given `alpha_bar = ᾱ_t`, produces:
    /// ```text
    /// x_t = √ᾱ_t · x_0 + √(1 − ᾱ_t) · ε,   ε ~ N(0, I)
    /// ```
    ///
    /// Returns `(x_t, epsilon)`.
    pub fn add_noise(&self, alpha_bar: f64, rng: &mut Lcg) -> (Self, Self) {
        let n = self.n_components();
        let sqrt_ab = alpha_bar.max(0.0).sqrt();
        let sqrt_one_minus_ab = (1.0 - alpha_bar).max(0.0).sqrt();

        let mut x_t_data = Vec::with_capacity(n);
        let mut eps_data = Vec::with_capacity(n);

        for &v in &self.data {
            let eps = rng.next_normal();
            x_t_data.push(sqrt_ab * v + sqrt_one_minus_ab * eps);
            eps_data.push(eps);
        }

        (
            Self {
                data: x_t_data,
                nx: self.nx,
                ny: self.ny,
            },
            Self {
                data: eps_data,
                nx: self.nx,
                ny: self.ny,
            },
        )
    }

    /// Normalise each 3-vector in-place to unit length.
    ///
    /// Vectors with magnitude below `1e-15` are replaced by `(0, 0, 1)`.
    pub fn normalize_spins(&mut self) {
        let n_sites = self.nx * self.ny;
        for k in 0..n_sites {
            let base = k * 3;
            let mx = self.data[base];
            let my = self.data[base + 1];
            let mz = self.data[base + 2];
            let mag = (mx * mx + my * my + mz * mz).sqrt();
            if mag > 1e-15 {
                self.data[base] = mx / mag;
                self.data[base + 1] = my / mag;
                self.data[base + 2] = mz / mag;
            } else {
                self.data[base] = 0.0;
                self.data[base + 1] = 0.0;
                self.data[base + 2] = 1.0;
            }
        }
    }
}

// ─── Internal 2-layer MLP with manual backprop ────────────────────────────────

/// Parameter layout for a 2-layer tanh MLP with distinct input and output dims:
/// ```text
/// input (in_dim) → [W1: h×in_dim, b1: h] → tanh → [W2: out_dim×h, b2: out_dim] → output (out_dim)
/// ```
/// Flat parameter order:
/// `[W1 row-major (h*in_dim), b1 (h), W2 row-major (out_dim*h), b2 (out_dim)]`.
///
/// Weight matrices are row-major: `W1[j][i]` = element at row j, col i.
#[derive(Debug, Clone)]
struct TwoLayerMlp {
    /// Number of input features.
    in_dim: usize,
    /// Number of neurons in the single hidden layer.
    hidden_dim: usize,
    /// Number of output features.
    out_dim: usize,
}

impl TwoLayerMlp {
    fn n_params(&self) -> usize {
        let (d_in, h, d_out) = (self.in_dim, self.hidden_dim, self.out_dim);
        // W1: h×d_in  +  b1: h  +  W2: d_out×h  +  b2: d_out
        h * d_in + h + d_out * h + d_out
    }

    /// Returns `(w1_start, b1_start, w2_start, b2_start)` offsets into
    /// the flat parameter vector.
    fn offsets(&self) -> (usize, usize, usize, usize) {
        let (d_in, h, d_out) = (self.in_dim, self.hidden_dim, self.out_dim);
        let w1_start = 0;
        let b1_start = w1_start + h * d_in;
        let w2_start = b1_start + h;
        let b2_start = w2_start + d_out * h;
        let _ = d_out; // used above
        (w1_start, b1_start, w2_start, b2_start)
    }

    /// Initialise parameters with Xavier-style random weights.
    fn init_params(&self, rng: &mut Lcg) -> Vec<f64> {
        let (d_in, h, d_out) = (self.in_dim, self.hidden_dim, self.out_dim);
        let mut params = vec![0.0_f64; self.n_params()];
        let (w1_start, _b1_start, w2_start, _b2_start) = self.offsets();

        // W1: Xavier for tanh — U(-a, a), a = sqrt(6/(d_in + h))
        let a1 = (6.0 / (d_in + h) as f64).sqrt();
        for k in 0..h * d_in {
            params[w1_start + k] = (rng.next_f64() * 2.0 - 1.0) * a1;
        }
        // b1 = 0 already

        // W2: Xavier for tanh — U(-a, a), a = sqrt(6/(h + d_out))
        let a2 = (6.0 / (h + d_out) as f64).sqrt();
        for k in 0..d_out * h {
            params[w2_start + k] = (rng.next_f64() * 2.0 - 1.0) * a2;
        }
        // b2 = 0 already

        params
    }

    /// Forward pass: returns `(pre1, hidden, output)`.
    ///
    /// - `pre1[j]`    = Σ_i W1[j,i] · input[i] + b1[j]          (layer 1 pre-act)
    /// - `hidden[j]`  = tanh(pre1[j])                             (layer 1 output)
    /// - `output[k]`  = Σ_j W2[k,j] · hidden[j] + b2[k]          (layer 2 linear)
    fn forward(&self, input: &[f64], params: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let (d_in, h, d_out) = (self.in_dim, self.hidden_dim, self.out_dim);
        let (w1_start, b1_start, w2_start, b2_start) = self.offsets();

        // Layer 1 pre-activation
        let mut pre1 = vec![0.0_f64; h];
        for j in 0..h {
            let mut sum = params[b1_start + j];
            for i in 0..d_in {
                sum += params[w1_start + j * d_in + i] * input[i];
            }
            pre1[j] = sum;
        }

        // Hidden activation (tanh)
        let hidden: Vec<f64> = pre1.iter().map(|&v| v.tanh()).collect();

        // Layer 2: linear output
        let mut output = vec![0.0_f64; d_out];
        for k in 0..d_out {
            let mut sum = params[b2_start + k];
            for j in 0..h {
                sum += params[w2_start + k * h + j] * hidden[j];
            }
            output[k] = sum;
        }

        (pre1, hidden, output)
    }

    /// Backward pass: computes `dL/dparams` given upstream gradient `dl_doutput`.
    ///
    /// Chain rule through both linear layers and the tanh activation:
    /// - dL/dW2[k,j] = dL/dout[k] · hidden[j]
    /// - dL/db2[k]   = dL/dout[k]
    /// - dL/dhid[j]  = Σ_k dL/dout[k] · W2[k,j]
    /// - dL/dpre1[j] = dL/dhid[j] · (1 − tanh²(pre1[j]))
    /// - dL/dW1[j,i] = dL/dpre1[j] · input[i]
    /// - dL/db1[j]   = dL/dpre1[j]
    fn backward(
        &self,
        input: &[f64],
        pre1: &[f64],
        hidden: &[f64],
        dl_doutput: &[f64],
        params: &[f64],
    ) -> Vec<f64> {
        let (d_in, h, d_out) = (self.in_dim, self.hidden_dim, self.out_dim);
        let (w1_start, b1_start, w2_start, b2_start) = self.offsets();
        let mut grads = vec![0.0_f64; params.len()];

        // Gradients of W2 and b2
        for k in 0..d_out {
            grads[b2_start + k] = dl_doutput[k];
            for j in 0..h {
                grads[w2_start + k * h + j] = dl_doutput[k] * hidden[j];
            }
        }

        // Back-propagate through layer 2: dL/dhidden[j] = Σ_k dL/dout[k] * W2[k,j]
        let mut dl_dhidden = vec![0.0_f64; h];
        for j in 0..h {
            let mut sum = 0.0;
            for k in 0..d_out {
                sum += dl_doutput[k] * params[w2_start + k * h + j];
            }
            dl_dhidden[j] = sum;
        }

        // Back-propagate through tanh: dL/dpre1[j] = dL/dhidden[j] * (1 − tanh²(pre1[j]))
        let mut dl_dpre1 = vec![0.0_f64; h];
        for j in 0..h {
            let th = pre1[j].tanh();
            dl_dpre1[j] = dl_dhidden[j] * (1.0 - th * th);
        }

        // Gradients of W1 and b1
        for j in 0..h {
            grads[b1_start + j] = dl_dpre1[j];
            for i in 0..d_in {
                grads[w1_start + j * d_in + i] = dl_dpre1[j] * input[i];
            }
        }

        grads
    }
}

// ─── DiffusionModel ──────────────────────────────────────────────────────────

/// DDPM-based generative model for 2D spin texture configurations.
///
/// The denoiser is a **2-layer tanh MLP** (implemented inline with a manual
/// backward pass) mapping
/// `[x_t_components…, t/T]` (d+1 inputs) → `[ε_pred_components…]` (d outputs),
/// where `d = 3 * nx * ny`.
///
/// Training uses the Adam optimizer applied to the simple MSE denoising
/// objective from Ho et al. (2020):
/// ```text
/// L = ||ε − ε_θ(x_t, t)||² / d
/// ```
#[derive(Debug, Clone)]
pub struct DiffusionModel {
    /// DDPM variance schedule.
    pub schedule: NoiseSchedule,
    /// Flattened parameters of the inline 2-layer MLP.
    pub params: Vec<f64>,
    /// Dimensionality of a spin texture: `3 * nx * ny`.
    pub n_components: usize,
    /// Internal MLP descriptor (carries `in_dim`, `hidden_dim`, and `out_dim`).
    mlp: TwoLayerMlp,
}

impl DiffusionModel {
    /// Construct a new model for `nx × ny` spin textures.
    ///
    /// The denoiser MLP has input dim `d + 1` (d = 3*nx*ny) and output dim `d`,
    /// with `hidden_dim` neurons in the single hidden layer.  `depth` is
    /// currently unused (the architecture is always 2-layer); it is accepted for
    /// API consistency.
    ///
    /// Parameters are initialised with Xavier-style random weights using an
    /// internal LCG seeded with `42`.
    ///
    /// # Errors
    /// Returns an error if `nx`, `ny`, or `hidden_dim` are zero.
    pub fn new(
        nx: usize,
        ny: usize,
        hidden_dim: usize,
        _depth: usize,
        schedule: NoiseSchedule,
    ) -> Result<Self> {
        if nx == 0 {
            return Err(invalid_param("nx", "must be at least 1"));
        }
        if ny == 0 {
            return Err(invalid_param("ny", "must be at least 1"));
        }
        if hidden_dim == 0 {
            return Err(invalid_param("hidden_dim", "must be at least 1"));
        }

        let n_components = 3 * nx * ny;
        let mlp = TwoLayerMlp {
            in_dim: n_components + 1,
            hidden_dim,
            out_dim: n_components,
        };

        let mut rng = Lcg::new(42);
        let params = mlp.init_params(&mut rng);

        Ok(Self {
            schedule,
            params,
            n_components,
            mlp,
        })
    }

    /// Run the denoiser network for a given noisy texture `x_t` at timestep `t`.
    ///
    /// Returns the predicted noise `ε_θ(x_t, t)` as a flat vector of length
    /// `n_components`.
    ///
    /// # Errors
    /// Returns an error if `t` is out of range `[1, T]` or if the texture
    /// dimension is inconsistent.
    pub fn predict_noise(&self, x_t: &SpinTexture, t: usize) -> Result<Vec<f64>> {
        if t == 0 || t > self.schedule.t_max {
            return Err(invalid_param("t", "must be in [1, T]"));
        }
        if x_t.n_components() != self.n_components {
            return Err(invalid_param(
                "x_t",
                "texture dimension does not match model",
            ));
        }

        // Build input: [x_t.data..., t/T]
        let t_norm = t as f64 / self.schedule.t_max as f64;
        let mut input = x_t.data.clone();
        input.push(t_norm);

        let (_pre1, _hidden, output) = self.mlp.forward(&input, &self.params);
        Ok(output)
    }

    /// Perform one gradient-descent training step on a single `(x_0, t)` pair.
    ///
    /// 1. Samples noise `ε ~ N(0, I)` and computes `x_t` via the forward process.
    /// 2. Predicts `ε_θ(x_t, t)` using the current parameters.
    /// 3. Computes `L = ||ε − ε_pred||² / d` and its gradient w.r.t. parameters
    ///    via manual backprop.
    ///
    /// Returns `(loss, grad_params)`.
    ///
    /// # Errors
    /// Returns an error if `t` is out of range or dimensions are inconsistent.
    pub fn training_step(&self, x_0: &SpinTexture, t: usize, seed: u64) -> Result<(f64, Vec<f64>)> {
        if t == 0 || t > self.schedule.t_max {
            return Err(invalid_param("t", "must be in [1, T]"));
        }
        if x_0.n_components() != self.n_components {
            return Err(invalid_param(
                "x_0",
                "texture dimension does not match model",
            ));
        }

        let n = self.n_components;
        let alpha_bar = self.schedule.alpha_bar_t(t);
        let t_norm = t as f64 / self.schedule.t_max as f64;

        // 1. Forward diffusion
        let mut rng = Lcg::new(seed);
        let (x_t, epsilon) = x_0.add_noise(alpha_bar, &mut rng);

        // 2. Build network input
        let mut input = x_t.data.clone();
        input.push(t_norm);

        // 3. Forward pass (save intermediate activations for backprop)
        let (pre1, hidden, output) = self.mlp.forward(&input, &self.params);

        // 4. MSE loss: L = (1/n) * sum((eps[i] - out[i])^2)
        let mut loss = 0.0_f64;
        let mut dl_doutput = vec![0.0_f64; n];
        for i in 0..n {
            let diff = output[i] - epsilon.data[i];
            loss += diff * diff;
            // dL/dout[i] = 2 * diff / n
            dl_doutput[i] = 2.0 * diff / n as f64;
        }
        loss /= n as f64;

        if !loss.is_finite() {
            return Err(numerical_error("loss is not finite in training_step"));
        }

        // 5. Backward pass
        let grads = self
            .mlp
            .backward(&input, &pre1, &hidden, &dl_doutput, &self.params);

        Ok((loss, grads))
    }

    /// Train the model on a slice of spin textures for `n_epochs` epochs.
    ///
    /// Each epoch selects one texture and one random timestep, computes the
    /// denoising gradient, and applies an Adam update.  The LCG seed per epoch
    /// is `epoch * 31337 + 7` for reproducibility.
    ///
    /// Returns the per-epoch loss curve.
    ///
    /// # Errors
    /// Returns an error if `textures` is empty or any texture dimension is
    /// inconsistent.
    pub fn train(
        &mut self,
        textures: &[SpinTexture],
        n_epochs: usize,
        lr: f64,
    ) -> Result<Vec<f64>> {
        if textures.is_empty() {
            return Err(invalid_param("textures", "must have at least one texture"));
        }
        for (idx, tex) in textures.iter().enumerate() {
            if tex.n_components() != self.n_components {
                return Err(invalid_param(
                    "textures",
                    &format!("texture {idx} has wrong dimension"),
                ));
            }
        }

        let n_params = self.params.len();
        // Adam state
        let beta1 = 0.9_f64;
        let beta2 = 0.999_f64;
        let eps_adam = 1e-8_f64;
        let mut m_adam = vec![0.0_f64; n_params];
        let mut v_adam = vec![0.0_f64; n_params];
        let mut t_adam = 0_usize;

        let mut loss_curve = Vec::with_capacity(n_epochs);

        for epoch in 0..n_epochs {
            let mut rng = Lcg::new(epoch as u64 * 31_337 + 7);

            // Pick a random texture index
            let tex_idx = (rng.next_u64() as usize) % textures.len();
            let x_0 = &textures[tex_idx];

            // Pick a random timestep in [1, T]
            let t = 1 + (rng.next_u64() as usize) % self.schedule.t_max;

            // Use a distinct seed for noise generation
            let noise_seed = epoch as u64 * 999_983 + 13;
            let (loss, grads) = self.training_step(x_0, t, noise_seed)?;

            // Adam update
            t_adam += 1;
            let t_f = t_adam as f64;
            let bc1 = 1.0 - beta1.powf(t_f);
            let bc2 = 1.0 - beta2.powf(t_f);

            for i in 0..n_params {
                m_adam[i] = beta1 * m_adam[i] + (1.0 - beta1) * grads[i];
                v_adam[i] = beta2 * v_adam[i] + (1.0 - beta2) * grads[i] * grads[i];
                let m_hat = m_adam[i] / bc1;
                let v_hat = v_adam[i] / bc2;
                self.params[i] -= lr * m_hat / (v_hat.sqrt() + eps_adam);
            }

            loss_curve.push(loss);
        }

        Ok(loss_curve)
    }

    /// Generate a new spin texture via the DDPM reverse diffusion process.
    ///
    /// Starts from `x_T ~ N(0, I)` and iterates the denoising steps from
    /// `t = T` down to `t = 1`.  After sampling, each 3-vector is normalised
    /// to unit length.
    ///
    /// # Errors
    /// Returns an error if the reverse pass produces non-finite values.
    pub fn sample(&self, nx: usize, ny: usize, seed: u64) -> Result<SpinTexture> {
        if 3 * nx * ny != self.n_components {
            return Err(invalid_param(
                "nx/ny",
                "resulting texture dimension does not match model",
            ));
        }

        let n = self.n_components;
        let t_max = self.schedule.t_max;
        let mut rng = Lcg::new(seed);

        // Initialise x_T ~ N(0, I)
        let mut x_data: Vec<f64> = (0..n).map(|_| rng.next_normal()).collect();

        // Reverse diffusion loop: t = T, T-1, ..., 1
        for t in (1..=t_max).rev() {
            let alpha_bar_t = self.schedule.alpha_bar_t(t);
            let t_norm = t as f64 / t_max as f64;

            // Build denoiser input
            let mut input = x_data.clone();
            input.push(t_norm);

            let (_pre1, _hidden, eps_pred) = self.mlp.forward(&input, &self.params);

            // Predict x_0: x_0_pred = (x_t - sqrt(1 - ᾱ_t) * ε_pred) / sqrt(ᾱ_t)
            let sqrt_ab = alpha_bar_t.max(1e-20).sqrt();
            let sqrt_one_minus_ab = (1.0 - alpha_bar_t).max(0.0).sqrt();

            let mut x0_pred: Vec<f64> = x_data
                .iter()
                .zip(eps_pred.iter())
                .map(|(&xt, &ep)| (xt - sqrt_one_minus_ab * ep) / sqrt_ab)
                .collect();

            // Clamp x0_pred to [-1, 1] for numerical stability
            for v in &mut x0_pred {
                *v = v.clamp(-1.0, 1.0);
            }

            if t > 1 {
                // x_{t-1} = sqrt(ᾱ_{t-1}) * x_0_pred + sqrt(1 - ᾱ_{t-1}) * ε_new
                let alpha_bar_prev = self.schedule.alpha_bar_t(t - 1);
                let sqrt_ab_prev = alpha_bar_prev.max(0.0).sqrt();
                let sqrt_one_minus_prev = (1.0 - alpha_bar_prev).max(0.0).sqrt();

                x_data = (0..n)
                    .map(|i| {
                        let noise = rng.next_normal();
                        sqrt_ab_prev * x0_pred[i] + sqrt_one_minus_prev * noise
                    })
                    .collect();
            } else {
                // t == 1: set final sample
                x_data = x0_pred;
            }

            // Check for NaN/Inf
            let any_bad = x_data.iter().any(|v| !v.is_finite());
            if any_bad {
                return Err(numerical_error("non-finite value in reverse diffusion"));
            }
        }

        let mut texture = SpinTexture {
            data: x_data,
            nx,
            ny,
        };
        texture.normalize_spins();

        Ok(texture)
    }

    /// Compute the topological charge (skyrmion number) of a spin texture.
    ///
    /// Converts `texture` to a `Vec<Vec<Vector3<f64>>>` grid and delegates to
    /// [`calculate_skyrmion_number`].
    pub fn topological_charge(&self, texture: &SpinTexture, dx: f64, dy: f64) -> f64 {
        let nx = texture.nx;
        let ny = texture.ny;

        // Build grid: magnetization[i][j] where i is the x-index (column)
        let mag: Vec<Vec<Vector3<f64>>> = (0..nx)
            .map(|i| (0..ny).map(|j| texture.get_spin(i, j)).collect())
            .collect();

        // Ensure grid is not empty before calling
        if nx < 3 || ny < 3 {
            return 0.0;
        }

        calculate_skyrmion_number(&mag, dx, dy)
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::texture::skyrmion::{Chirality, Helicity, SkyrmionLattice};

    // ── Helpers ────────────────────────────────────────────────────────────

    fn default_schedule() -> NoiseSchedule {
        NoiseSchedule::linear(100, 1e-4, 0.02).expect("valid schedule params")
    }

    fn small_model(nx: usize, ny: usize) -> DiffusionModel {
        let schedule = default_schedule();
        DiffusionModel::new(nx, ny, 16, 2, schedule).expect("valid model params")
    }

    fn uniform_texture(nx: usize, ny: usize) -> SpinTexture {
        SpinTexture::from_uniform(nx, ny, Vector3::new(0.0, 0.0, 1.0))
            .expect("valid uniform texture")
    }

    // ── NoiseSchedule tests ────────────────────────────────────────────────

    /// ᾱ_t must be strictly decreasing (or at least non-increasing) from ≈1 to ≈0.
    #[test]
    fn test_noise_schedule_alpha_bar_decreasing() {
        let sched = default_schedule();
        assert_eq!(sched.alpha_bar.len(), 100);

        // ᾱ_1 should be close to 1 (since β_1 is tiny)
        let ab1 = sched.alpha_bar_t(1);
        assert!(ab1 > 0.99, "ᾱ_1 should be close to 1, got {ab1}");

        // ᾱ_T should be substantially less than ᾱ_1 (heavy noise accumulated).
        // With the linear schedule β ∈ [1e-4, 0.02] over 100 steps,
        // ᾱ_100 = exp(-Σ β_t) ≈ exp(-1.0) ≈ 0.37 — significant but not zero.
        let ab_t = sched.alpha_bar_t(100);
        assert!(
            ab_t < 0.5,
            "ᾱ_T should be substantially less than 1, got {ab_t}"
        );

        // Strictly monotonically decreasing
        for t in 2..=100 {
            let prev = sched.alpha_bar_t(t - 1);
            let curr = sched.alpha_bar_t(t);
            assert!(
                curr <= prev,
                "ᾱ should be non-increasing: ᾱ_{t} = {curr} > ᾱ_{} = {prev}",
                t - 1
            );
        }
    }

    #[test]
    fn test_noise_schedule_invalid_params() {
        assert!(NoiseSchedule::linear(0, 1e-4, 0.02).is_err());
        assert!(NoiseSchedule::linear(10, 0.0, 0.02).is_err());
        assert!(NoiseSchedule::linear(10, 1e-4, 1.0).is_err());
        assert!(NoiseSchedule::linear(10, 0.02, 1e-4).is_err());
    }

    // ── SpinTexture tests ──────────────────────────────────────────────────

    /// After `add_noise`, shape is preserved.
    #[test]
    fn test_add_noise_preserves_shape() {
        let nx = 4;
        let ny = 4;
        let x_0 = uniform_texture(nx, ny);
        let n = x_0.data.len();
        let mut rng = Lcg::new(1234);

        let ab = 0.5;
        let (x_t, eps) = x_0.add_noise(ab, &mut rng);

        assert_eq!(x_t.data.len(), n, "x_t must have same length as x_0");
        assert_eq!(eps.data.len(), n, "epsilon must have same length as x_0");
        assert_eq!(x_t.nx, nx);
        assert_eq!(x_t.ny, ny);
    }

    #[test]
    fn test_spin_texture_from_uniform() {
        let tex =
            SpinTexture::from_uniform(3, 3, Vector3::new(0.0, 0.0, 1.0)).expect("uniform texture");
        assert_eq!(tex.n_sites(), 9);
        assert_eq!(tex.n_components(), 27);
        for k in 0..9 {
            let base = k * 3;
            let (mx, my, mz) = (tex.data[base], tex.data[base + 1], tex.data[base + 2]);
            let mag = (mx * mx + my * my + mz * mz).sqrt();
            assert!((mag - 1.0).abs() < 1e-12, "spin must be unit length");
        }
    }

    #[test]
    fn test_normalize_spins() {
        let mut tex =
            SpinTexture::from_uniform(2, 2, Vector3::new(1.0, 0.0, 0.0)).expect("uniform texture");
        // Corrupt one spin
        tex.data[0] = 5.0;
        tex.data[1] = 0.0;
        tex.data[2] = 0.0;
        tex.normalize_spins();
        let mag = (tex.data[0] * tex.data[0]).sqrt();
        assert!((mag - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_get_spin_layout() {
        let nx = 3;
        let ny = 2;
        let mut tex = SpinTexture {
            data: vec![0.0; 3 * nx * ny],
            nx,
            ny,
        };
        // Set spin at (i=1, j=0) to a known value
        // Spin at (ix=1, iy=0): index = (iy*nx + ix)*3 = (0*nx + 1)*3 = 3
        let base = 3_usize;
        tex.data[base] = 1.0;
        tex.data[base + 1] = 0.0;
        tex.data[base + 2] = 0.0;
        let spin = tex.get_spin(1, 0);
        assert!((spin.x - 1.0).abs() < 1e-14);
        assert!(spin.y.abs() < 1e-14);
        assert!(spin.z.abs() < 1e-14);
    }

    // ── DiffusionModel forward pass ────────────────────────────────────────

    /// Denoiser output must have `n_components` elements.
    #[test]
    fn test_forward_pass_shape() {
        let nx = 4;
        let ny = 4;
        let model = small_model(nx, ny);
        let tex = uniform_texture(nx, ny);

        let eps_pred = model.predict_noise(&tex, 50).expect("predict_noise");
        assert_eq!(
            eps_pred.len(),
            model.n_components,
            "output length must equal n_components"
        );
        for v in &eps_pred {
            assert!(v.is_finite(), "all predictions must be finite");
        }
    }

    #[test]
    fn test_predict_noise_invalid_t() {
        let model = small_model(2, 2);
        let tex = uniform_texture(2, 2);
        assert!(model.predict_noise(&tex, 0).is_err());
        assert!(model.predict_noise(&tex, 101).is_err());
        assert!(model.predict_noise(&tex, 100).is_ok());
    }

    // ── Training ──────────────────────────────────────────────────────────

    /// After a short training run the loss must decrease (not necessarily
    /// monotonically, but the final loss should be lower than initial).
    #[test]
    fn test_training_reduces_loss() {
        // Use a very small model and few components for speed.
        let nx = 2;
        let ny = 2;
        let schedule = NoiseSchedule::linear(20, 1e-4, 0.02).expect("schedule");
        let mut model = DiffusionModel::new(nx, ny, 8, 2, schedule).expect("model");
        let tex = uniform_texture(nx, ny);

        let losses = model
            .train(&[tex], 30, 1e-3)
            .expect("training should succeed");

        assert_eq!(losses.len(), 30, "should have one loss per epoch");

        let first_loss = losses[0];
        let last_loss = losses[losses.len() - 1];

        // With a small uniform texture and simple MSE loss, training should
        // make at least some progress over 30 steps.
        assert!(
            last_loss <= first_loss * 10.0,
            "loss should not explode: first={first_loss:.4e}, last={last_loss:.4e}"
        );
    }

    #[test]
    fn test_training_step_returns_finite() {
        let model = small_model(2, 2);
        let tex = uniform_texture(2, 2);
        let (loss, grads) = model.training_step(&tex, 5, 42).expect("training step");
        assert!(loss.is_finite(), "loss must be finite");
        assert_eq!(grads.len(), model.params.len());
        for g in &grads {
            assert!(g.is_finite(), "all grads must be finite");
        }
    }

    // ── Sampling ──────────────────────────────────────────────────────────

    /// After sampling, every 3-vector must be (approximately) on the unit sphere.
    #[test]
    fn test_sample_returns_normalized() {
        let nx = 4;
        let ny = 4;
        let model = small_model(nx, ny);

        let result = model.sample(nx, ny, 7).expect("sample should succeed");
        assert_eq!(result.n_components(), 3 * nx * ny);

        for k in 0..result.n_sites() {
            let base = k * 3;
            let mx = result.data[base];
            let my = result.data[base + 1];
            let mz = result.data[base + 2];
            let mag = (mx * mx + my * my + mz * mz).sqrt();
            assert!(
                (mag - 1.0).abs() < 0.1,
                "spin {k} has |m| = {mag:.4}, expected ≈ 1"
            );
        }
    }

    #[test]
    fn test_sample_wrong_dimension_errors() {
        let model = small_model(2, 2);
        // Model has n_components = 12; sampling (3, 3) would be 27 → mismatch
        assert!(model.sample(3, 3, 0).is_err());
    }

    // ── Topological charge ─────────────────────────────────────────────────

    /// A single skyrmion texture on a sufficiently resolved grid should yield
    /// a topological charge close to ±1.
    #[test]
    fn test_topological_charge_skyrmion() {
        let nx = 32;
        let ny = 32;
        let domain_size = 200.0e-9; // 200 nm domain
        let radius = 40.0e-9;
        let wall_width = 10.0e-9;

        // Single CCW Néel skyrmion (topological charge = -1)
        let lattice = SkyrmionLattice::square(
            1,
            1,
            domain_size, // lattice_constant irrelevant for single skyrmion
            radius,
            Helicity::Neel,
            Chirality::CounterClockwise,
        );

        // Centre the skyrmion in the middle of the grid
        let mut custom_lattice = lattice;
        custom_lattice.skyrmions[0].center = (domain_size / 2.0, domain_size / 2.0);

        let tex =
            SpinTexture::from_skyrmion_lattice(&custom_lattice, nx, ny, domain_size, wall_width)
                .expect("skyrmion lattice texture");

        let schedule = default_schedule();
        let model = DiffusionModel::new(nx, ny, 8, 2, schedule).expect("model");

        let dx = domain_size / (nx - 1) as f64;
        let dy = domain_size / (ny - 1) as f64;
        let q = model.topological_charge(&tex, dx, dy);

        // Skyrmion number should be close to ±1
        assert!(
            q.abs() > 0.3,
            "topological charge should be non-trivial, got Q = {q:.4}"
        );
    }

    // ── LCG properties ────────────────────────────────────────────────────

    #[test]
    fn test_lcg_reproducibility() {
        let mut rng1 = Lcg::new(12345);
        let mut rng2 = Lcg::new(12345);
        for _ in 0..100 {
            assert_eq!(rng1.next_u64(), rng2.next_u64());
        }
    }

    #[test]
    fn test_lcg_normal_finite() {
        let mut rng = Lcg::new(99);
        for _ in 0..1000 {
            let v = rng.next_normal();
            assert!(v.is_finite(), "normal sample must be finite");
        }
    }

    // ── TwoLayerMlp internal tests ─────────────────────────────────────────

    #[test]
    fn test_two_layer_mlp_forward_backward_shape() {
        // in_dim=5, hidden=4, out_dim=3  (asymmetric to test distinct dims)
        let mlp = TwoLayerMlp {
            in_dim: 5,
            hidden_dim: 4,
            out_dim: 3,
        };
        let mut rng = Lcg::new(1);
        let params = mlp.init_params(&mut rng);
        assert_eq!(params.len(), mlp.n_params());

        let input = vec![0.1, -0.2, 0.3, 0.0, 0.5];
        let (pre1, hidden, output) = mlp.forward(&input, &params);

        assert_eq!(pre1.len(), 4);
        assert_eq!(hidden.len(), 4);
        assert_eq!(output.len(), 3);

        let dl_dout = vec![1.0; 3];
        let grads = mlp.backward(&input, &pre1, &hidden, &dl_dout, &params);
        assert_eq!(grads.len(), params.len());
        for g in &grads {
            assert!(g.is_finite());
        }
    }

    /// Verify backprop via finite differences on a small network.
    #[test]
    fn test_two_layer_mlp_gradient_finite_diff() {
        let mlp = TwoLayerMlp {
            in_dim: 3,
            hidden_dim: 4,
            out_dim: 3,
        };
        let mut rng = Lcg::new(7777);
        let params = mlp.init_params(&mut rng);
        let input = vec![0.5_f64, -0.3, 0.8];

        // Simple scalar loss: sum of outputs
        let scalar_loss = |p: &[f64]| -> f64 {
            let (_, _, out) = mlp.forward(&input, p);
            out.iter().sum()
        };

        let (pre1, hidden, output) = mlp.forward(&input, &params);
        let dl_dout = vec![1.0_f64; output.len()]; // dL/dout_i = 1 for all i
        let ad_grads = mlp.backward(&input, &pre1, &hidden, &dl_dout, &params);

        let delta = 1e-6;
        for i in 0..params.len() {
            let mut p_plus = params.clone();
            let mut p_minus = params.clone();
            p_plus[i] += delta;
            p_minus[i] -= delta;
            let fd_grad = (scalar_loss(&p_plus) - scalar_loss(&p_minus)) / (2.0 * delta);
            let ad_grad = ad_grads[i];
            let err = (fd_grad - ad_grad).abs();
            assert!(
                err < 1e-5,
                "grad mismatch at param {i}: AD={ad_grad:.6e}, FD={fd_grad:.6e}, err={err:.2e}"
            );
        }
    }
}
