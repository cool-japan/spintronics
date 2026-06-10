//! O(3)-equivariant neural network layers for spin Hamiltonians.
//!
//! Physical spin energies `E(S_1, …, S_N)` are *invariant* under simultaneous
//! rotation of all spin inputs (SO(3) symmetry).  Standard MLPs must *learn*
//! this symmetry from data, wasting capacity; equivariant layers bake it in
//! by construction.
//!
//! This module implements a **Cartesian-tensor formulation** (simpler than
//! the full Clebsch–Gordan / E(3)-equivariant machinery):
//!   - Each layer has scalar channels (`s_i ∈ ℝ`) and vector channels
//!     (`v_i ∈ ℝ³`).
//!   - Linear transforms mix scalars with `|v|` magnitudes (rotation-invariant)
//!     and vectors with linear combinations of vectors plus scalar-vector
//!     gating (rotation-equivariant).
//!   - Gated activation: `σ(s) · v` keeps vectors equivariant; `tanh(s)` on
//!     scalars is non-linear yet rotation-invariant.
//!   - Final scalar output → rotation-invariant energy.
//!   - Final vector output → rotation-equivariant field.
//!
//! ## Why does this work?
//!
//! Under a rotation `R ∈ SO(3)`, vectors transform as `v ↦ R v`, while scalars
//! and inner products `v · w` are invariant.  The forward equations are
//! constructed exclusively from operations that respect this:
//!
//! ```text
//! s_out[i] = Σ_j w_ss[i][j] · s[j]  +  Σ_k w_vs[i][k] · ‖v[k]‖  +  bias_s[i]
//! v_out[i] = Σ_k w_vv[i][k] · v[k]  +  Σ_j w_sv[i][j] · s[j] · v[j]     (paired)
//! ```
//!
//! The first sum on `v_out` is linear in the input vectors (equivariant),
//! and the second multiplies each input vector by an input scalar — still
//! linear in the *vectors*, with a scalar coefficient that is rotation
//! invariant.  No `bias` is added to vectors because a constant offset would
//! break equivariance.
//!
//! ## Reference
//!
//! K. T. Schütt, O. Unke & M. Gastegger, "Equivariant Message Passing for the
//! Prediction of Tensorial Properties and Molecular Spectra",
//! *ICML* (2021), arXiv:2102.03150.

use crate::error::{dimension_mismatch, invalid_param, Result};
use crate::vector3::Vector3;

// ─── Linear-congruential pseudo-random generator ─────────────────────────────

/// Deterministic Knuth-style 64-bit LCG used for weight initialisation.
///
/// Constants are taken from MMIX:
/// `a = 6_364_136_223_846_793_005`, `c = 1_442_695_040_888_963_407`.
/// The high 53 bits of each draw form a uniform `f64` in `[0, 1)`.
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

    fn next_unit(&mut self) -> f64 {
        let bits = self.next_u64() >> 11;
        (bits as f64) / ((1u64 << 53) as f64)
    }

    fn next_uniform(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.next_unit()
    }

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

// ─── EquivariantConfig ───────────────────────────────────────────────────────

/// Shape descriptor for an [`EquivariantLinear`] layer.
///
/// A layer maps `(n_scalar_in scalars, n_vector_in vectors)` to
/// `(n_scalar_out scalars, n_vector_out vectors)`.  For the scalar-vector
/// mixing term (the `w_sv · s · v` channel) we *require*
/// `n_scalar_in == n_vector_in`; when they differ, the cross term is
/// silently zeroed and the layer falls back to pure linear vector mixing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EquivariantConfig {
    /// Number of input scalar channels.
    pub n_scalar_in: usize,
    /// Number of input vector channels.
    pub n_vector_in: usize,
    /// Number of output scalar channels.
    pub n_scalar_out: usize,
    /// Number of output vector channels.
    pub n_vector_out: usize,
}

impl EquivariantConfig {
    /// Construct a new shape descriptor.
    pub const fn new(
        n_scalar_in: usize,
        n_vector_in: usize,
        n_scalar_out: usize,
        n_vector_out: usize,
    ) -> Self {
        Self {
            n_scalar_in,
            n_vector_in,
            n_scalar_out,
            n_vector_out,
        }
    }

    /// `true` when the scalar-vector mixing block can be used (matching dims).
    pub const fn has_paired_sv(&self) -> bool {
        self.n_scalar_in == self.n_vector_in && self.n_scalar_in > 0
    }

    /// Validate that at least one of the input/output dimensions is non-zero.
    ///
    /// We tolerate `n_scalar_in == 0` (energy network with vector inputs only)
    /// and `n_vector_out == 0` (pure invariant readout).
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when both input
    /// dimensions are zero, or both output dimensions are zero.
    pub fn validate(&self) -> Result<()> {
        if self.n_scalar_in == 0 && self.n_vector_in == 0 {
            return Err(invalid_param(
                "EquivariantConfig",
                "must have at least one input channel (scalar or vector)",
            ));
        }
        if self.n_scalar_out == 0 && self.n_vector_out == 0 {
            return Err(invalid_param(
                "EquivariantConfig",
                "must have at least one output channel (scalar or vector)",
            ));
        }
        Ok(())
    }
}

// ─── EquivariantLinear ───────────────────────────────────────────────────────

/// Single-layer Cartesian-equivariant linear map.
///
/// All weight matrices are stored row-major: `w[output_row][input_col]`.
///
/// | Weight | Shape                            | Role                                   |
/// |--------|----------------------------------|----------------------------------------|
/// | `w_ss` | `[n_scalar_out][n_scalar_in]`    | scalar → scalar                        |
/// | `w_vs` | `[n_scalar_out][n_vector_in]`    | ‖vector‖ → scalar (invariant)          |
/// | `w_sv` | `[n_vector_out][n_scalar_in]`    | scalar × paired vector → vector        |
/// | `w_vv` | `[n_vector_out][n_vector_in]`    | vector → vector                        |
/// | `bias_s` | `[n_scalar_out]`               | scalar bias (vectors must not have one) |
///
/// The `w_sv` block is exercised only when
/// [`EquivariantConfig::has_paired_sv`] returns `true`; otherwise it is held
/// at zero to keep the layer well-defined.
pub struct EquivariantLinear {
    /// Shape descriptor.
    pub config: EquivariantConfig,
    /// Scalar → scalar weights `[n_scalar_out][n_scalar_in]`.
    pub w_ss: Vec<Vec<f64>>,
    /// Vector-magnitude → scalar weights `[n_scalar_out][n_vector_in]`.
    pub w_vs: Vec<Vec<f64>>,
    /// Scalar → vector gating weights `[n_vector_out][n_scalar_in]`
    /// (multiplies the paired input vector).
    pub w_sv: Vec<Vec<f64>>,
    /// Vector → vector weights `[n_vector_out][n_vector_in]`.
    pub w_vv: Vec<Vec<f64>>,
    /// Scalar bias `[n_scalar_out]`.
    pub bias_s: Vec<f64>,
}

impl EquivariantLinear {
    /// Construct an equivariant linear layer with randomly-initialised weights.
    ///
    /// Initialisation:
    ///   - `w_ss`, `w_vv`: Xavier-uniform with `a = √(6/(in+out))` — appropriate
    ///     for `tanh`-style gated activations downstream.
    ///   - `w_vs`, `w_sv`: He-normal with `σ = √(2/in)` — these feed into the
    ///     possibly-unbounded magnitude channel and into the gated vector
    ///     channel, so He keeps activations from collapsing.
    ///   - `bias_s`: zero.
    ///
    /// # Errors
    /// Propagates [`EquivariantConfig::validate`].
    pub fn new(config: EquivariantConfig, rng_seed: u64) -> Result<Self> {
        config.validate()?;
        let mut rng = Lcg::new(rng_seed);

        // w_ss: Xavier(in=n_s_in, out=n_s_out)
        let w_ss = if config.n_scalar_in > 0 && config.n_scalar_out > 0 {
            let a = (6.0 / (config.n_scalar_in + config.n_scalar_out) as f64).sqrt();
            (0..config.n_scalar_out)
                .map(|_| {
                    (0..config.n_scalar_in)
                        .map(|_| rng.next_uniform(-a, a))
                        .collect()
                })
                .collect()
        } else {
            vec![vec![0.0; config.n_scalar_in]; config.n_scalar_out]
        };

        // w_vs: He(in=n_v_in, out=n_s_out)
        let w_vs = if config.n_vector_in > 0 && config.n_scalar_out > 0 {
            let std_dev = (2.0 / config.n_vector_in as f64).sqrt();
            (0..config.n_scalar_out)
                .map(|_| {
                    (0..config.n_vector_in)
                        .map(|_| std_dev * rng.next_normal())
                        .collect()
                })
                .collect()
        } else {
            vec![vec![0.0; config.n_vector_in]; config.n_scalar_out]
        };

        // w_sv: He(in=n_s_in, out=n_v_out) — only used when paired
        let w_sv = if config.has_paired_sv() && config.n_vector_out > 0 {
            let std_dev = (2.0 / config.n_scalar_in as f64).sqrt();
            (0..config.n_vector_out)
                .map(|_| {
                    (0..config.n_scalar_in)
                        .map(|_| std_dev * rng.next_normal())
                        .collect()
                })
                .collect()
        } else {
            vec![vec![0.0; config.n_scalar_in]; config.n_vector_out]
        };

        // w_vv: Xavier(in=n_v_in, out=n_v_out)
        let w_vv = if config.n_vector_in > 0 && config.n_vector_out > 0 {
            let a = (6.0 / (config.n_vector_in + config.n_vector_out) as f64).sqrt();
            (0..config.n_vector_out)
                .map(|_| {
                    (0..config.n_vector_in)
                        .map(|_| rng.next_uniform(-a, a))
                        .collect()
                })
                .collect()
        } else {
            vec![vec![0.0; config.n_vector_in]; config.n_vector_out]
        };

        let bias_s = vec![0.0_f64; config.n_scalar_out];

        Ok(Self {
            config,
            w_ss,
            w_vs,
            w_sv,
            w_vv,
            bias_s,
        })
    }

    /// Forward pass — `(scalars, vectors) ↦ (out_scalars, out_vectors)`.
    ///
    /// Mathematically, for each output index `i`:
    ///
    /// ```text
    /// s_out[i] = Σ_j w_ss[i][j] · s[j] + Σ_k w_vs[i][k] · ‖v[k]‖ + bias_s[i]
    /// v_out[i] = Σ_k w_vv[i][k] · v[k] + Σ_j w_sv[i][j] · s[j] · v[j]
    /// ```
    ///
    /// where the second term in `v_out` is only evaluated when the layer's
    /// scalar and vector input dimensions match.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the input
    /// dimensions do not match [`EquivariantConfig`].
    pub fn forward(
        &self,
        scalars: &[f64],
        vectors: &[Vector3<f64>],
    ) -> Result<(Vec<f64>, Vec<Vector3<f64>>)> {
        if scalars.len() != self.config.n_scalar_in {
            return Err(dimension_mismatch(
                &format!("{} scalar inputs", self.config.n_scalar_in),
                &format!("{} scalar inputs", scalars.len()),
            ));
        }
        if vectors.len() != self.config.n_vector_in {
            return Err(dimension_mismatch(
                &format!("{} vector inputs", self.config.n_vector_in),
                &format!("{} vector inputs", vectors.len()),
            ));
        }

        // Pre-compute vector magnitudes (rotation-invariant scalars).
        let mut mags = Vec::with_capacity(vectors.len());
        for v in vectors {
            mags.push(v.magnitude());
        }

        // Scalar output: w_ss · s + w_vs · |v| + bias_s.
        let mut out_scalars = vec![0.0_f64; self.config.n_scalar_out];
        for (i, out_s) in out_scalars.iter_mut().enumerate() {
            let mut acc = self.bias_s[i];
            for (j, &sj) in scalars.iter().enumerate() {
                acc += self.w_ss[i][j] * sj;
            }
            for (k, &mk) in mags.iter().enumerate() {
                acc += self.w_vs[i][k] * mk;
            }
            *out_s = acc;
        }

        // Vector output: w_vv · v  + (optional) Σ_j w_sv[i][j] · s[j] · v[j].
        let mut out_vectors = vec![Vector3::<f64>::zero(); self.config.n_vector_out];
        let paired = self.config.has_paired_sv();
        for (i, out_v) in out_vectors.iter_mut().enumerate() {
            let mut acc = Vector3::<f64>::zero();
            for (k, vk) in vectors.iter().enumerate() {
                acc = acc + (*vk * self.w_vv[i][k]);
            }
            if paired {
                for (j, &sj) in scalars.iter().enumerate() {
                    // Pair index j ↔ vector j (n_scalar_in == n_vector_in).
                    let coef = self.w_sv[i][j] * sj;
                    acc = acc + vectors[j] * coef;
                }
            }
            *out_v = acc;
        }

        Ok((out_scalars, out_vectors))
    }

    /// Total number of free parameters (sum of all weight matrices + biases).
    pub fn n_params(&self) -> usize {
        let c = &self.config;
        c.n_scalar_out * c.n_scalar_in        // w_ss
            + c.n_scalar_out * c.n_vector_in  // w_vs
            + c.n_vector_out * c.n_scalar_in  // w_sv
            + c.n_vector_out * c.n_vector_in  // w_vv
            + c.n_scalar_out // bias_s
    }

    /// Flatten all parameters in a stable, reversible order:
    /// `w_ss → w_vs → w_sv → w_vv → bias_s`.
    pub fn params_flat(&self) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.n_params());
        for row in &self.w_ss {
            v.extend_from_slice(row);
        }
        for row in &self.w_vs {
            v.extend_from_slice(row);
        }
        for row in &self.w_sv {
            v.extend_from_slice(row);
        }
        for row in &self.w_vv {
            v.extend_from_slice(row);
        }
        v.extend_from_slice(&self.bias_s);
        v
    }

    /// Inverse of [`EquivariantLinear::params_flat`].
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the slice
    /// length does not equal [`EquivariantLinear::n_params`].
    pub fn set_params(&mut self, flat: &[f64]) -> Result<()> {
        let expected = self.n_params();
        if flat.len() != expected {
            return Err(dimension_mismatch(
                &format!("{expected} params"),
                &format!("{} params", flat.len()),
            ));
        }
        let c = self.config;
        let mut cur = 0_usize;

        for i in 0..c.n_scalar_out {
            for j in 0..c.n_scalar_in {
                self.w_ss[i][j] = flat[cur];
                cur += 1;
            }
        }
        for i in 0..c.n_scalar_out {
            for k in 0..c.n_vector_in {
                self.w_vs[i][k] = flat[cur];
                cur += 1;
            }
        }
        for i in 0..c.n_vector_out {
            for j in 0..c.n_scalar_in {
                self.w_sv[i][j] = flat[cur];
                cur += 1;
            }
        }
        for i in 0..c.n_vector_out {
            for k in 0..c.n_vector_in {
                self.w_vv[i][k] = flat[cur];
                cur += 1;
            }
        }
        for i in 0..c.n_scalar_out {
            self.bias_s[i] = flat[cur];
            cur += 1;
        }
        Ok(())
    }
}

// ─── EquivariantMlp ──────────────────────────────────────────────────────────

/// Multi-layer perceptron composed of [`EquivariantLinear`] blocks with
/// equivariance-preserving gated activations between them.
///
/// Between two linear layers we apply:
///   - `tanh(s)` element-wise to scalars (still invariant under rotations
///     because each scalar is itself invariant);
///   - `v ← v` unchanged for vectors (the gated multiplication that
///     gives this architecture its expressivity is already absorbed into
///     the linear layers via the `w_sv` block).
///
/// This separation makes the proof of equivariance straightforward: every
/// operation that touches a vector is either linear in that vector, or
/// scales it by a rotation-invariant scalar.  The composition therefore
/// preserves the equivariance property.
pub struct EquivariantMlp {
    /// Layers, in evaluation order.
    pub layers: Vec<EquivariantLinear>,
    /// Cached input shape (= `layers\[0\].config`).
    pub input_config: EquivariantConfig,
    /// Cached output shape (= `layers.last().config`).
    pub output_config: EquivariantConfig,
}

impl EquivariantMlp {
    /// Build a stacked equivariant MLP from a list of layer configurations.
    ///
    /// `layer_configs[k+1].n_scalar_in` must equal `layer_configs[k].n_scalar_out`
    /// (and similarly for vector channels) so successive layers compose.
    ///
    /// A per-layer sub-seed is derived from `rng_seed` via a Weyl-style stride
    /// to keep layer-to-layer weights independent.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when `layer_configs`
    /// is empty or when adjacent layers' shapes are incompatible.
    pub fn new(layer_configs: &[EquivariantConfig], rng_seed: u64) -> Result<Self> {
        if layer_configs.is_empty() {
            return Err(invalid_param("layer_configs", "must contain ≥ 1 layer"));
        }
        for window in layer_configs.windows(2) {
            let prev = window[0];
            let next = window[1];
            if prev.n_scalar_out != next.n_scalar_in {
                return Err(invalid_param(
                    "layer_configs",
                    "scalar output of layer k must equal scalar input of layer k+1",
                ));
            }
            if prev.n_vector_out != next.n_vector_in {
                return Err(invalid_param(
                    "layer_configs",
                    "vector output of layer k must equal vector input of layer k+1",
                ));
            }
        }
        let mut layers = Vec::with_capacity(layer_configs.len());
        for (k, cfg) in layer_configs.iter().enumerate() {
            let sub_seed = rng_seed.wrapping_add((k as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            layers.push(EquivariantLinear::new(*cfg, sub_seed)?);
        }
        Ok(Self {
            input_config: layer_configs[0],
            output_config: *layer_configs.last().expect("non-empty checked above"),
            layers,
        })
    }

    /// Differentiability-free forward pass through every layer.  A `tanh`
    /// non-linearity is applied to scalars *between* (not after) layers; the
    /// final layer's outputs are returned linearly so they can serve as
    /// energies, fields, or further downstream inputs.
    ///
    /// # Errors
    /// Propagates dimension-mismatch errors from inner layers.
    pub fn forward(
        &self,
        scalars: &[f64],
        vectors: &[Vector3<f64>],
    ) -> Result<(Vec<f64>, Vec<Vector3<f64>>)> {
        let mut s = scalars.to_vec();
        let mut v = vectors.to_vec();
        let last_idx = self.layers.len() - 1;
        for (k, layer) in self.layers.iter().enumerate() {
            let (out_s, out_v) = layer.forward(&s, &v)?;
            s = out_s;
            v = out_v;
            if k != last_idx {
                // Element-wise tanh on scalars; vectors pass through unchanged.
                for sj in s.iter_mut() {
                    *sj = sj.tanh();
                }
            }
        }
        Ok((s, v))
    }

    /// Convenience: compute the rotation-invariant scalar energy of a spin
    /// configuration `spins`.
    ///
    /// The network is expected to have `input_config.n_scalar_in == 0`,
    /// `input_config.n_vector_in == spins.len()`, and
    /// `output_config.n_scalar_out ≥ 1`.  The first output scalar is taken
    /// as the energy (additional outputs are ignored).
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the shapes do
    /// not match, or [`crate::error::Error::InvalidParameter`] when the
    /// network exposes no scalar output.
    pub fn energy(&self, spins: &[Vector3<f64>]) -> Result<f64> {
        if self.output_config.n_scalar_out == 0 {
            return Err(invalid_param(
                "EquivariantMlp",
                "network must expose ≥ 1 scalar output to act as an energy",
            ));
        }
        let scalars: Vec<f64> = Vec::new();
        let (out_s, _out_v) = self.forward(&scalars, spins)?;
        Ok(out_s[0])
    }

    /// Total number of free parameters across all layers.
    pub fn n_params(&self) -> usize {
        self.layers.iter().map(|l| l.n_params()).sum()
    }

    /// Concatenate every layer's flat parameter vector in layer order.
    pub fn params_flat(&self) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.n_params());
        for layer in &self.layers {
            v.extend(layer.params_flat());
        }
        v
    }

    /// Inverse of [`EquivariantMlp::params_flat`].
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the length
    /// is wrong.
    pub fn set_params(&mut self, flat: &[f64]) -> Result<()> {
        let expected = self.n_params();
        if flat.len() != expected {
            return Err(dimension_mismatch(
                &format!("{expected} params"),
                &format!("{} params", flat.len()),
            ));
        }
        let mut cursor = 0_usize;
        for layer in &mut self.layers {
            let n = layer.n_params();
            layer.set_params(&flat[cursor..cursor + n])?;
            cursor += n;
        }
        Ok(())
    }
}

// ─── Rotation utilities ──────────────────────────────────────────────────────

/// Apply a 3×3 rotation matrix to a [`Vector3<f64>`]: `R · v`.
///
/// The matrix is interpreted as row-major: `r\[i\][j]` is element `R_{ij}`.
pub fn rotate_vector(r: &[[f64; 3]; 3], v: Vector3<f64>) -> Vector3<f64> {
    Vector3::new(
        r[0][0] * v.x + r[0][1] * v.y + r[0][2] * v.z,
        r[1][0] * v.x + r[1][1] * v.y + r[1][2] * v.z,
        r[2][0] * v.x + r[2][1] * v.y + r[2][2] * v.z,
    )
}

/// Sample a (close-to) uniformly-distributed rotation matrix on SO(3).
///
/// The construction uses a Rodrigues rotation around a Marsaglia-normal
/// random axis with a random angle in `[0, π]`.  This is *not* the
/// mathematically optimal Haar measure (which would use Shoemake's algorithm
/// on quaternions), but it suffices for the equivariance tests in this module:
/// we only need a rotation matrix `R` with `R · Rᵀ = I` and `det(R) = +1`.
pub fn random_so3(rng_seed: u64) -> [[f64; 3]; 3] {
    let mut rng = Lcg::new(rng_seed);

    // Sample a uniform unit vector via three independent normals.
    let mut axis_x = rng.next_normal();
    let mut axis_y = rng.next_normal();
    let mut axis_z = rng.next_normal();
    let norm = (axis_x * axis_x + axis_y * axis_y + axis_z * axis_z).sqrt();
    if norm < 1e-14 {
        // Astronomical odds; fall back to the canonical z-axis.
        axis_x = 0.0;
        axis_y = 0.0;
        axis_z = 1.0;
    } else {
        axis_x /= norm;
        axis_y /= norm;
        axis_z /= norm;
    }
    let theta = std::f64::consts::PI * rng.next_unit();
    let c = theta.cos();
    let s = theta.sin();
    let cm1 = 1.0 - c;
    let (ax, ay, az) = (axis_x, axis_y, axis_z);

    // Rodrigues: R = I + sin(θ)·K + (1−cos(θ))·K²,
    // where K is the skew-symmetric matrix of `axis`.
    [
        [
            c + ax * ax * cm1,
            ax * ay * cm1 - az * s,
            ax * az * cm1 + ay * s,
        ],
        [
            ay * ax * cm1 + az * s,
            c + ay * ay * cm1,
            ay * az * cm1 - ax * s,
        ],
        [
            az * ax * cm1 - ay * s,
            az * ay * cm1 + ax * s,
            c + az * az * cm1,
        ],
    ]
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn vec_approx(a: Vector3<f64>, b: Vector3<f64>, tol: f64) -> bool {
        approx(a.x, b.x, tol) && approx(a.y, b.y, tol) && approx(a.z, b.z, tol)
    }

    // 1. Construct + config validation: zero-everything inputs/outputs error.
    #[test]
    fn test_config_validation() {
        let good = EquivariantConfig::new(0, 3, 4, 0);
        assert!(good.validate().is_ok());
        let no_input = EquivariantConfig::new(0, 0, 1, 1);
        assert!(no_input.validate().is_err());
        let no_output = EquivariantConfig::new(2, 2, 0, 0);
        assert!(no_output.validate().is_err());
        assert!(EquivariantLinear::new(no_input, 0).is_err());
    }

    // 2. n_params() agrees with the per-block accounting.
    #[test]
    fn test_n_params_accounting() {
        let cfg = EquivariantConfig::new(3, 3, 4, 2);
        let layer = EquivariantLinear::new(cfg, 13).unwrap();
        let expected = 4 * 3 + 4 * 3 + 2 * 3 + 2 * 3 + 4; // ss + vs + sv + vv + bias_s
        assert_eq!(layer.n_params(), expected);
        assert_eq!(layer.params_flat().len(), expected);
    }

    // 3. params_flat / set_params round-trip preserves every weight bit-for-bit.
    #[test]
    fn test_params_roundtrip() {
        let cfg = EquivariantConfig::new(2, 2, 3, 3);
        let mut layer = EquivariantLinear::new(cfg, 99).unwrap();
        let original = layer.params_flat();
        let mut perturbed = original.clone();
        for (i, v) in perturbed.iter_mut().enumerate() {
            *v += 0.01 * (i as f64);
        }
        layer.set_params(&perturbed).unwrap();
        let rt = layer.params_flat();
        for (a, b) in rt.iter().zip(perturbed.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    // 4. Rotation invariance of the scalar energy output.
    //
    // We construct a 3-spin → 1-scalar network: rotate all input spins by R
    // and confirm the scalar energy is preserved to machine precision.
    #[test]
    fn test_rotation_invariance_of_energy() {
        let layer1 = EquivariantConfig::new(0, 3, 4, 4);
        let layer2 = EquivariantConfig::new(4, 4, 1, 0);
        let mlp = EquivariantMlp::new(&[layer1, layer2], 42).unwrap();
        let spins = vec![
            Vector3::new(0.5, -0.3, 0.7),
            Vector3::new(-0.2, 0.8, 0.1),
            Vector3::new(0.9, 0.1, -0.4),
        ];
        let r = random_so3(7);
        let rotated: Vec<Vector3<f64>> = spins.iter().map(|s| rotate_vector(&r, *s)).collect();
        let e0 = mlp.energy(&spins).unwrap();
        let er = mlp.energy(&rotated).unwrap();
        assert!(
            approx(e0, er, 1e-10),
            "energy not invariant: {} vs {}",
            e0,
            er
        );
    }

    // 5. Rotation equivariance of the vector output.
    //
    // For a network with vector output, the output should rotate by exactly R
    // when the inputs are rotated by R.
    #[test]
    fn test_rotation_equivariance_of_vector_output() {
        let cfg = EquivariantConfig::new(0, 3, 0, 2);
        let layer = EquivariantLinear::new(cfg, 11).unwrap();
        let spins = vec![
            Vector3::new(0.3, 0.4, -0.5),
            Vector3::new(-0.6, 0.2, 0.8),
            Vector3::new(0.1, -0.9, 0.4),
        ];
        let r = random_so3(123);
        let rotated_in: Vec<Vector3<f64>> = spins.iter().map(|s| rotate_vector(&r, *s)).collect();
        let (_s0, v0) = layer.forward(&[], &spins).unwrap();
        let (_sr, vr) = layer.forward(&[], &rotated_in).unwrap();
        for (orig, after) in v0.iter().zip(vr.iter()) {
            let expected = rotate_vector(&r, *orig);
            assert!(
                vec_approx(expected, *after, 1e-10),
                "vector equivariance broken: expected R·v = ({},{},{}) got ({},{},{})",
                expected.x,
                expected.y,
                expected.z,
                after.x,
                after.y,
                after.z,
            );
        }
    }

    // 6. Identity at zero weights: a layer whose every weight is zero outputs
    //    zeros (plus the zero bias).
    #[test]
    fn test_zero_weights_zero_output() {
        let cfg = EquivariantConfig::new(0, 2, 1, 2);
        let mut layer = EquivariantLinear::new(cfg, 0).unwrap();
        let zeros = vec![0.0_f64; layer.n_params()];
        layer.set_params(&zeros).unwrap();
        let (out_s, out_v) = layer
            .forward(
                &[],
                &[Vector3::new(1.0, 2.0, 3.0), Vector3::new(0.5, 0.5, 0.5)],
            )
            .unwrap();
        assert_eq!(out_s, vec![0.0]);
        for v in &out_v {
            assert!(vec_approx(*v, Vector3::zero(), 1e-15));
        }
    }

    // 7. Ferromagnetically-aligned spins yield a different energy from a
    //    disordered configuration (after random init, the network is generic
    //    enough that the two energies do not happen to coincide).
    #[test]
    fn test_fm_vs_disordered_energy_differs() {
        let layer1 = EquivariantConfig::new(0, 4, 5, 4);
        let layer2 = EquivariantConfig::new(5, 4, 1, 0);
        let mlp = EquivariantMlp::new(&[layer1, layer2], 2024).unwrap();
        let fm = vec![Vector3::unit_z(); 4];
        let disordered = vec![
            Vector3::unit_z(),
            Vector3::unit_x(),
            Vector3::unit_y(),
            Vector3::new(0.0, 0.0, -1.0),
        ];
        let e_fm = mlp.energy(&fm).unwrap();
        let e_dis = mlp.energy(&disordered).unwrap();
        assert!(e_fm.is_finite() && e_dis.is_finite());
        assert!(
            (e_fm - e_dis).abs() > 1e-8,
            "FM and disordered should differ, got {} vs {}",
            e_fm,
            e_dis
        );
    }

    // 8. Multi-layer composition preserves invariance/equivariance through
    //    several stacks.
    #[test]
    fn test_multi_layer_invariance() {
        let cfgs = [
            EquivariantConfig::new(0, 3, 5, 5),
            EquivariantConfig::new(5, 5, 4, 4),
            EquivariantConfig::new(4, 4, 1, 1),
        ];
        let mlp = EquivariantMlp::new(&cfgs, 314).unwrap();
        let spins = vec![
            Vector3::new(0.2, 0.3, 0.6),
            Vector3::new(-0.5, 0.4, 0.1),
            Vector3::new(0.7, -0.2, 0.0),
        ];
        let r = random_so3(2718);
        let rotated: Vec<Vector3<f64>> = spins.iter().map(|s| rotate_vector(&r, *s)).collect();
        let (s0, v0) = mlp.forward(&[], &spins).unwrap();
        let (sr, vr) = mlp.forward(&[], &rotated).unwrap();
        assert!(approx(s0[0], sr[0], 1e-10));
        let expected_v = rotate_vector(&r, v0[0]);
        assert!(vec_approx(expected_v, vr[0], 1e-10));
    }

    // 9. rotate_vector preserves vector magnitude (a basic sanity check on the
    //    rotation matrix structure).
    #[test]
    fn test_rotate_vector_preserves_magnitude() {
        let r = random_so3(8);
        let v = Vector3::new(0.7, -0.3, 1.2);
        let m0 = v.magnitude();
        let mr = rotate_vector(&r, v).magnitude();
        assert!(approx(m0, mr, 1e-12));
    }

    // 10. random_so3 produces an orthonormal matrix (R · Rᵀ = I).
    #[test]
    fn test_random_so3_is_orthonormal() {
        let r = random_so3(31415);
        // Compute R · Rᵀ.
        let mut rrt = [[0.0_f64; 3]; 3];
        for (i, row_i) in r.iter().enumerate() {
            for (j, row_j) in r.iter().enumerate() {
                let acc: f64 = row_i.iter().zip(row_j.iter()).map(|(a, b)| a * b).sum();
                rrt[i][j] = acc;
            }
        }
        for (i, row) in rrt.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approx(val, expected, 1e-12),
                    "RRᵀ[{}][{}] = {} ≠ {}",
                    i,
                    j,
                    val,
                    expected
                );
            }
        }
    }

    // 11. Seed reproducibility: same seed → bitwise-identical parameters and
    //     bitwise-identical SO(3) matrices.
    #[test]
    fn test_seed_reproducibility() {
        let cfg = EquivariantConfig::new(2, 2, 3, 3);
        let l1 = EquivariantLinear::new(cfg, 555).unwrap();
        let l2 = EquivariantLinear::new(cfg, 555).unwrap();
        let p1 = l1.params_flat();
        let p2 = l2.params_flat();
        for (a, b) in p1.iter().zip(p2.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        let r1 = random_so3(777);
        let r2 = random_so3(777);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(r1[i][j].to_bits(), r2[i][j].to_bits());
            }
        }
    }

    // 12. Gated activation (tanh on scalars between layers) preserves
    //     equivariance of the vector output: rotate inputs, get rotated outputs.
    //     This is tested through a deeper (3-layer) MLP, complementing test 8.
    #[test]
    fn test_gated_activation_preserves_equivariance() {
        let cfgs = [
            EquivariantConfig::new(2, 2, 4, 4),
            EquivariantConfig::new(4, 4, 3, 3),
            EquivariantConfig::new(3, 3, 1, 2),
        ];
        let mlp = EquivariantMlp::new(&cfgs, 654321).unwrap();
        let scalars = vec![0.4_f64, -0.2];
        let spins = vec![Vector3::new(0.1, 0.5, -0.2), Vector3::new(0.3, -0.6, 0.7)];
        let r = random_so3(987);
        let rotated_spins: Vec<Vector3<f64>> =
            spins.iter().map(|s| rotate_vector(&r, *s)).collect();
        let (s0, v0) = mlp.forward(&scalars, &spins).unwrap();
        // Scalars are rotation invariants — pass them unchanged on both sides.
        let (sr, vr) = mlp.forward(&scalars, &rotated_spins).unwrap();
        assert!(approx(s0[0], sr[0], 1e-10));
        for (o, n) in v0.iter().zip(vr.iter()) {
            let expected = rotate_vector(&r, *o);
            assert!(vec_approx(expected, *n, 1e-10));
        }
    }

    // 13. Large input batch (many spins) runs to completion and yields finite,
    //     non-trivial outputs.
    #[test]
    fn test_large_input_runs() {
        let n_spins = 64;
        let cfgs = [
            EquivariantConfig::new(0, n_spins, 8, n_spins),
            EquivariantConfig::new(8, n_spins, 1, 0),
        ];
        let mlp = EquivariantMlp::new(&cfgs, 4242).unwrap();
        let mut spins = Vec::with_capacity(n_spins);
        for i in 0..n_spins {
            let phase = (i as f64) * 0.15;
            spins.push(Vector3::new(phase.cos(), phase.sin(), 0.0));
        }
        let e = mlp.energy(&spins).unwrap();
        assert!(e.is_finite(), "energy must be finite for many-spin input");
    }
}
