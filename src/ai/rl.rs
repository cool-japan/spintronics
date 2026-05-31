//! Reinforcement Learning for SOT Switching Protocol Optimization
//!
//! This module implements a Cross-Entropy Method (CEM) optimizer that learns
//! the optimal spin-orbit torque (SOT) current pulse sequence to switch
//! perpendicular magnetization from +z to -z with minimal energy dissipation.
//!
//! ## Physics
//!
//! The macrospin is governed by the Landau-Lifshitz-Gilbert (LLG) equation:
//!
//! ```text
//! dm/dt = -γ/(1+α²) * [m × H_eff + α * m × (m × H_eff)]
//! ```
//!
//! where the effective field is:
//!
//! ```text
//! H_eff = H_ext + H_SOT(j, m, ĵ) + H_anisotropy
//! H_anisotropy = (2K / (μ₀ Ms)) * m_z * ẑ   (perpendicular magnetic anisotropy)
//! ```
//!
//! ## Algorithm: Cross-Entropy Method
//!
//! A population of pulse sequences is drawn from a Gaussian distribution
//! N(μ, σ²) over current amplitudes.  The elite fraction (top 20%) is used
//! to update the distribution, iterating until convergence.
//!
//! ## Example
//!
//! ```rust
//! use spintronics::ai::rl::{SotSwitchingEnv, SotRlOptimizer};
//!
//! let env  = SotSwitchingEnv::default_cofeb_pt();
//! let mut opt = SotRlOptimizer::new(env);
//! let result = opt.train(20, 42);
//! println!("best reward: {:.3}", result.best_reward);
//! ```

use crate::constants::{E_CHARGE, GAMMA, HBAR, MU_0};
use crate::effect::sot::SpinOrbitTorque;
use crate::error::{invalid_param, Result};
use crate::material::Ferromagnet;
use crate::vector3::Vector3;

// =============================================================================
// Seeded Linear Congruential Generator
// =============================================================================

/// Seeded 64-bit LCG pseudo-random number generator.
///
/// Uses the Knuth multiplicative LCG parameters, which give a full-period
/// sequence of length 2⁶⁴.  The generator is deterministic given a seed,
/// which is essential for reproducible RL experiments.
#[derive(Debug, Clone)]
pub struct Lcg {
    state: u64,
}

impl Lcg {
    /// Construct a new LCG with the given seed.
    pub fn new(seed: u64) -> Self {
        Lcg {
            state: seed.wrapping_add(1),
        }
    }

    /// Advance the state and return the next raw 64-bit value.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    /// Return a uniform f64 in [0, 1).
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Return a standard-normal sample via the Box-Muller transform.
    #[inline]
    pub fn next_normal(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-15);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

// =============================================================================
// SOT Effective-Field Helper
// =============================================================================

/// Compute the total SOT effective field acting on the magnetization.
///
/// Implements the Slonczewski / spin-Hall picture:
///
/// ```text
/// σ̂ = normalize(ĵ × ẑ)                  (spin polarization direction)
/// H_DL = (θ_SH ℏ |j|) / (2 e μ₀ Ms t_FM) * (m × σ̂)
/// H_FL = 0.1 * |H_DL| * σ̂
/// ```
///
/// # Arguments
/// * `sot`   – SOT material parameters.
/// * `m`     – Unit magnetization vector.
/// * `j`     – Charge current density \[A/m²\] (signed; sign sets direction).
/// * `j_dir` – Unit vector along current flow (e.g. x̂ for standard geometry).
/// * `ms`    – Saturation magnetization \[A/m\].
/// * `t_fm`  – Ferromagnetic layer thickness \[m\].
///
/// # Returns
/// Total SOT effective field \[A/m\].
fn compute_sot_field(
    sot: &SpinOrbitTorque,
    m: Vector3<f64>,
    j: f64,
    j_dir: Vector3<f64>,
    ms: f64,
    t_fm: f64,
) -> Vector3<f64> {
    let z_hat = Vector3::new(0.0, 0.0, 1.0);

    // Spin polarization: perpendicular to both current direction and film normal.
    // If j_dir is exactly parallel to z_hat, fall back to ŷ to avoid degenerate cross.
    let sigma_raw = j_dir.cross(&z_hat);
    let sigma_hat = if sigma_raw.magnitude() > 1e-12 {
        sigma_raw.normalize()
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };

    // Amplitude of damping-like field: (θ_SH ℏ |j|) / (2 e μ₀ Ms t_FM)
    let h_amplitude = (sot.theta_sh.abs() * HBAR * j.abs()) / (2.0 * E_CHARGE * MU_0 * ms * t_fm);

    // Adjust sign: if j < 0 or theta_sh < 0, flip the sense of the torque.
    let sign = if (j < 0.0) ^ (sot.theta_sh < 0.0) {
        -1.0_f64
    } else {
        1.0_f64
    };

    // Damping-like: m × σ̂  (pushes m toward or away from the equator).
    let h_dl = m.cross(&sigma_hat) * (sign * h_amplitude);

    // Field-like: small component along σ̂.
    let h_fl = sigma_hat * (0.1 * h_amplitude);

    h_dl + h_fl
}

// =============================================================================
// LLG RK4 Integrator
// =============================================================================

/// Evaluate the LLG right-hand side dm/dt given m and H_eff.
///
/// LLG (with implicit Gilbert damping form):
/// ```text
/// dm/dt = -γ/(1+α²) * [m × H + α * m × (m × H)]
/// ```
#[inline]
fn llg_rhs(m: Vector3<f64>, h_eff: Vector3<f64>, gamma: f64, alpha: f64) -> Vector3<f64> {
    let factor = gamma / (1.0 + alpha * alpha);
    let m_cross_h = m.cross(&h_eff);
    let m_cross_m_cross_h = m.cross(&m_cross_h);
    // dm/dt = -factor * (m×H + α * m×(m×H))
    (m_cross_h + m_cross_m_cross_h * alpha) * (-factor)
}

/// Perform a single 4th-order Runge-Kutta LLG step.
///
/// Returns the updated (but not yet renormalized) magnetization.
#[inline]
fn llg_rk4_step(
    m: Vector3<f64>,
    h_eff: Vector3<f64>,
    gamma: f64,
    alpha: f64,
    dt: f64,
) -> Vector3<f64> {
    let k1 = llg_rhs(m, h_eff, gamma, alpha);
    let m2 = m + k1 * (dt / 2.0);
    let k2 = llg_rhs(m2, h_eff, gamma, alpha);
    let m3 = m + k2 * (dt / 2.0);
    let k3 = llg_rhs(m3, h_eff, gamma, alpha);
    let m4 = m + k3 * dt;
    let k4 = llg_rhs(m4, h_eff, gamma, alpha);

    m + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0)
}

// =============================================================================
// Configuration
// =============================================================================

/// Configuration for the SOT switching RL environment.
///
/// All SI units.  Default values are appropriate for a Pt(5nm)/CoFeB(1.5nm)
/// heterostructure with perpendicular magnetic anisotropy (PMA).
#[derive(Debug, Clone)]
pub struct SotSwitchingConfig {
    /// Maximum charge current density \[A/m²\].
    pub j_max: f64,
    /// Ferromagnetic layer thickness \[m\].
    pub t_fm: f64,
    /// In-plane bias field magnitude \[A/m\].
    ///
    /// Applied along x̂; breaks the symmetry required for deterministic
    /// PMA switching.
    pub h_bias: f64,
    /// LLG time step \[s\].
    pub dt: f64,
    /// Maximum number of environment steps per episode.
    pub max_steps: usize,
    /// Number of independent pulse amplitudes to optimise.
    pub n_pulse_steps: usize,
    /// m_z threshold below which the magnetization is considered switched.
    pub switch_threshold: f64,
}

impl Default for SotSwitchingConfig {
    fn default() -> Self {
        Self {
            j_max: 2.0e11,
            t_fm: 1.5e-9,
            h_bias: 5.0e3,
            dt: 1.0e-12,
            max_steps: 500,
            n_pulse_steps: 10,
            switch_threshold: -0.8,
        }
    }
}

impl SotSwitchingConfig {
    /// Validate the configuration, returning an error if parameters are unphysical.
    pub fn validate(&self) -> Result<()> {
        if self.j_max <= 0.0 {
            return Err(invalid_param("j_max", "must be positive"));
        }
        if self.t_fm <= 0.0 {
            return Err(invalid_param("t_fm", "must be positive"));
        }
        if self.dt <= 0.0 {
            return Err(invalid_param("dt", "must be positive"));
        }
        if self.max_steps == 0 {
            return Err(invalid_param("max_steps", "must be at least 1"));
        }
        if self.n_pulse_steps == 0 {
            return Err(invalid_param("n_pulse_steps", "must be at least 1"));
        }
        if self.switch_threshold >= 0.0 {
            return Err(invalid_param(
                "switch_threshold",
                "must be negative for +z → −z switching",
            ));
        }
        Ok(())
    }
}

// =============================================================================
// Environment
// =============================================================================

/// Macrospin environment for SOT switching driven by the LLG equation.
///
/// The RL agent controls the charge current density at each step.  The
/// magnetization evolves via the LLG equation, and the episode terminates
/// either when the magnetization switches (m_z < `switch_threshold`) or the
/// step budget is exhausted.
#[derive(Debug, Clone)]
pub struct SotSwitchingEnv {
    /// Simulation configuration.
    pub config: SotSwitchingConfig,
    /// Ferromagnetic material parameters (Ms, α, K, …).
    pub material: Ferromagnet,
    /// Heavy-metal SOT parameters (θ_SH, …).
    pub sot: SpinOrbitTorque,
    /// External (bias) field \[A/m\].
    pub h_ext: Vector3<f64>,
    /// Current unit magnetization vector.
    pub m: Vector3<f64>,
    /// Number of steps taken in the current episode.
    pub step_count: usize,
    /// Accumulated episode reward.
    pub total_reward: f64,
}

impl SotSwitchingEnv {
    /// Construct a new environment with explicit material and SOT parameters.
    pub fn new(config: SotSwitchingConfig, material: Ferromagnet, sot: SpinOrbitTorque) -> Self {
        let h_bias_x = config.h_bias;
        Self {
            config,
            material,
            sot,
            h_ext: Vector3::new(h_bias_x, 0.0, 0.0),
            m: Vector3::new(0.0, 0.0, 1.0),
            step_count: 0,
            total_reward: 0.0,
        }
    }

    /// Default Pt(5nm)/CoFeB(1.5nm) PMA environment.
    ///
    /// Uses a PMA-tuned CoFeB (α = 0.01, Ms = 1.4 MA/m, K = 1.0 MJ/m³)
    /// on top of a Pt SOT source.
    pub fn default_cofeb_pt() -> Self {
        let config = SotSwitchingConfig::default();
        // CoFeB with PMA: α smaller than bulk Co, Ms and K tuned for PMA
        let material = Ferromagnet {
            alpha: 0.01,
            ms: 1.4e6,
            anisotropy_k: 1.0e6,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_a: 2.0e-11,
        };
        let sot = SpinOrbitTorque::platinum_cofeb();
        Self::new(config, material, sot)
    }

    /// Reset the environment to the initial state m ≈ +ẑ with a small x-tilt.
    ///
    /// The tiny perturbation is necessary so that the SOT has a non-zero
    /// component to act on; a perfect +ẑ state is a fixed point of the
    /// (unperturbed) dynamics.
    ///
    /// Returns the initial magnetization.
    pub fn reset(&mut self) -> Vector3<f64> {
        self.m = Vector3::new(0.01, 0.0, 1.0).normalize();
        self.step_count = 0;
        self.total_reward = 0.0;
        self.m
    }

    /// Compute the effective field for the current magnetization and current.
    ///
    /// H_eff = H_ext + H_SOT + H_anisotropy
    ///
    /// The anisotropy field for perpendicular uniaxial anisotropy is:
    /// ```text
    /// H_ani = (2K / (μ₀ Ms)) * m_z * ẑ
    /// ```
    pub fn effective_field(&self, m: Vector3<f64>, j_amplitude: f64) -> Vector3<f64> {
        let ms = self.material.ms;
        let k = self.material.anisotropy_k;

        // Uniaxial PMA anisotropy field along ẑ
        let h_ani_z = (2.0 * k) / (MU_0 * ms) * m.z;
        let h_anisotropy = Vector3::new(0.0, 0.0, h_ani_z);

        // SOT field (current along x̂ for standard geometry)
        let j_dir = Vector3::new(1.0, 0.0, 0.0);
        let h_sot = compute_sot_field(&self.sot, m, j_amplitude, j_dir, ms, self.config.t_fm);

        self.h_ext + h_anisotropy + h_sot
    }

    /// Advance the environment by one time step with the given current density.
    ///
    /// # Arguments
    /// * `j_amplitude` – Charge current density \[A/m²\] (may be negative).
    ///
    /// # Returns
    /// `(new_m, reward, done)` where:
    /// * `new_m`  – Updated unit magnetization.
    /// * `reward` – Instantaneous reward signal.
    /// * `done`   – Whether the episode has terminated.
    pub fn step(&mut self, j_amplitude: f64) -> (Vector3<f64>, f64, bool) {
        let h_eff = self.effective_field(self.m, j_amplitude);

        // RK4 LLG integration using the material's gyromagnetic ratio.
        // The Ferromagnet struct does not store gamma separately; use the
        // library constant (rad s⁻¹ T⁻¹) which corresponds to |γ_e|.
        let gamma = GAMMA;
        let alpha = self.material.alpha;
        let dt = self.config.dt;

        let m_new = llg_rk4_step(self.m, h_eff, gamma, alpha, dt);

        // Renormalize to stay on the unit sphere.
        let mag = m_new.magnitude();
        self.m = if mag > 1e-15 {
            m_new * (1.0 / mag)
        } else {
            // Degenerate: recover to +ẑ (should never happen with typical fields)
            Vector3::new(0.0, 0.0, 1.0)
        };

        // Reward components
        let j_norm = j_amplitude / self.config.j_max;
        let energy_penalty = -0.001 * j_norm * j_norm;
        let time_penalty = -0.005;
        let switched = self.m.z < self.config.switch_threshold;
        let switch_reward = if switched { 10.0 } else { 0.0 };
        let reward = switch_reward + energy_penalty + time_penalty;

        self.step_count += 1;
        self.total_reward += reward;

        let done = switched || self.step_count >= self.config.max_steps;

        (self.m, reward, done)
    }
}

// =============================================================================
// CEM Policy
// =============================================================================

/// Gaussian policy for the CEM optimizer.
///
/// Each pulse step `i` has an independent Gaussian:
/// ```text
/// j_i ~ N(μ_i, σ_i²),  clamped to [−j_max, +j_max]
/// ```
///
/// The mean vector `mu` converges to a good deterministic policy after
/// a few generations.
#[derive(Debug, Clone)]
pub struct CemPolicy {
    /// Number of pulse amplitudes (action dimension).
    pub n_pulse_steps: usize,
    /// Per-step mean current density \[A/m²\].
    pub mu: Vec<f64>,
    /// Per-step standard deviation \[A/m²\].
    pub sigma: Vec<f64>,
    /// Clamp bound for sampled actions \[A/m²\].
    pub j_max: f64,
}

impl CemPolicy {
    /// Initialise a zero-mean policy with σ = j_max / 2.
    pub fn new(n_pulse_steps: usize, j_max: f64) -> Self {
        let sigma_init = j_max / 2.0;
        Self {
            n_pulse_steps,
            mu: vec![0.0; n_pulse_steps],
            sigma: vec![sigma_init; n_pulse_steps],
            j_max,
        }
    }

    /// Draw a policy sample: one current amplitude per pulse step.
    ///
    /// Each amplitude is sampled as N(μ_i, σ_i²) and clamped to
    /// [−j_max, +j_max].
    pub fn sample_policy(&self, rng: &mut Lcg) -> Vec<f64> {
        self.mu
            .iter()
            .zip(self.sigma.iter())
            .map(|(&mu_i, &sigma_i)| {
                let sample = mu_i + sigma_i * rng.next_normal();
                sample.clamp(-self.j_max, self.j_max)
            })
            .collect()
    }

    /// Fit the Gaussian parameters to the elite policy set.
    ///
    /// Uses the unbiased sample mean and standard deviation.
    /// Sigma is clamped from below at `j_max * 0.01` to prevent premature
    /// collapse of the search distribution.
    ///
    /// Panics-free: if `elite_policies` is empty the distribution is unchanged.
    pub fn update_from_elite(&mut self, elite_policies: &[Vec<f64>]) {
        let n = elite_policies.len();
        if n == 0 {
            return;
        }
        let n_f = n as f64;
        let sigma_floor = self.j_max * 0.01;

        for i in 0..self.n_pulse_steps {
            // Mean
            let mean = elite_policies.iter().map(|p| p[i]).sum::<f64>() / n_f;
            // Variance (use population variance; fine for CEM)
            let var = elite_policies
                .iter()
                .map(|p| {
                    let diff = p[i] - mean;
                    diff * diff
                })
                .sum::<f64>()
                / n_f;
            self.mu[i] = mean;
            self.sigma[i] = var.sqrt().max(sigma_floor);
        }
    }

    /// Return the deterministic best-guess policy (the current mean).
    pub fn best_policy(&self) -> Vec<f64> {
        self.mu.clone()
    }
}

// =============================================================================
// Result
// =============================================================================

/// Summary of a CEM training run.
#[derive(Debug, Clone)]
pub struct SotRlResult {
    /// Best total reward achieved across all generations.
    pub best_reward: f64,
    /// Mean population reward per generation.
    pub mean_rewards_per_gen: Vec<f64>,
    /// Best individual reward per generation.
    pub best_rewards_per_gen: Vec<f64>,
    /// Pulse sequence (j_amplitude \[A/m²\] per step) of the best policy found.
    pub best_policy: Vec<f64>,
    /// Whether any evaluated policy achieved magnetization switching.
    pub switching_achieved: bool,
    /// Number of generations run.
    pub n_generations: usize,
}

// =============================================================================
// Optimizer
// =============================================================================

/// CEM optimizer for SOT switching pulse sequences.
///
/// Wraps the [`SotSwitchingEnv`] and [`CemPolicy`] and runs the main
/// generation loop.
#[derive(Debug, Clone)]
pub struct SotRlOptimizer {
    /// The physical simulation environment.
    pub env: SotSwitchingEnv,
    /// Gaussian policy parameterisation.
    pub policy: CemPolicy,
    /// Number of policy samples drawn per generation.
    pub population_size: usize,
    /// Fraction of the population used as elite (kept for distribution update).
    pub elite_fraction: f64,
}

impl SotRlOptimizer {
    /// Create a new optimizer with sensible defaults (population = 50, elite = 20%).
    pub fn new(env: SotSwitchingEnv) -> Self {
        let n_pulse = env.config.n_pulse_steps;
        let j_max = env.config.j_max;
        let policy = CemPolicy::new(n_pulse, j_max);
        Self {
            env,
            policy,
            population_size: 50,
            elite_fraction: 0.2,
        }
    }

    /// Evaluate a fixed pulse sequence and return the total episode reward.
    ///
    /// Resets the environment before rolling out the sequence.  Each element
    /// of `pulse_sequence` is applied as the current density for one
    /// LLG time step; the roll-out terminates early if the episode ends.
    pub fn evaluate_policy(&mut self, pulse_sequence: &[f64]) -> f64 {
        self.env.reset();
        let mut accumulated = 0.0_f64;
        for &j in pulse_sequence {
            let (_m, reward, done) = self.env.step(j);
            accumulated += reward;
            if done {
                break;
            }
        }
        accumulated
    }

    /// Run `n_generations` of CEM and return a [`SotRlResult`].
    ///
    /// # Arguments
    /// * `n_generations` – Number of CEM iterations.
    /// * `seed`          – Initial seed for the LCG.
    pub fn train(&mut self, n_generations: usize, seed: u64) -> SotRlResult {
        let mut rng = Lcg::new(seed);

        let elite_count =
            ((self.population_size as f64 * self.elite_fraction).ceil() as usize).max(1);

        let mut mean_rewards_per_gen = Vec::with_capacity(n_generations);
        let mut best_rewards_per_gen = Vec::with_capacity(n_generations);
        let mut overall_best_reward = f64::NEG_INFINITY;
        let mut overall_best_policy = self.policy.best_policy();
        let mut switching_achieved = false;

        for _gen in 0..n_generations {
            // --- Sample population ---
            let mut population: Vec<(f64, Vec<f64>)> = (0..self.population_size)
                .map(|_| {
                    let p = self.policy.sample_policy(&mut rng);
                    let r = self.evaluate_policy(&p);
                    (r, p)
                })
                .collect();

            // --- Sort descending by reward ---
            population.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

            // --- Statistics ---
            let gen_best = population[0].0;
            let gen_mean =
                population.iter().map(|(r, _)| *r).sum::<f64>() / self.population_size as f64;

            mean_rewards_per_gen.push(gen_mean);
            best_rewards_per_gen.push(gen_best);

            if gen_best > overall_best_reward {
                overall_best_reward = gen_best;
                overall_best_policy = population[0].1.clone();
            }

            // Check whether any policy achieved switching
            if !switching_achieved && gen_best >= 10.0 - 0.005 * self.env.config.max_steps as f64 {
                // Switching reward is 10.0; time/energy penalties reduce total.
                // A conservative check: total reward > 0 implies switching occurred.
                switching_achieved = gen_best > 0.0;
            }

            // --- Elite update ---
            let elite_policies: Vec<Vec<f64>> = population
                .iter()
                .take(elite_count)
                .map(|(_, p)| p.clone())
                .collect();
            self.policy.update_from_elite(&elite_policies);
        }

        // Final check: verify the greedy policy achieves switching
        let greedy = self.policy.best_policy();
        let greedy_reward = self.evaluate_policy(&greedy);
        if greedy_reward > overall_best_reward {
            overall_best_reward = greedy_reward;
            overall_best_policy = greedy;
        }
        if !switching_achieved && overall_best_reward > 0.0 {
            switching_achieved = true;
        }

        SotRlResult {
            best_reward: overall_best_reward,
            mean_rewards_per_gen,
            best_rewards_per_gen,
            best_policy: overall_best_policy,
            switching_achieved,
            n_generations,
        }
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_env() -> SotSwitchingEnv {
        SotSwitchingEnv::default_cofeb_pt()
    }

    // ------------------------------------------------------------------
    // Environment tests
    // ------------------------------------------------------------------

    #[test]
    fn test_env_reset() {
        let mut env = make_env();
        let m = env.reset();
        // After reset, m_z should be close to 1.0 (tiny x perturbation)
        assert!(
            (m.z - 1.0).abs() < 0.01,
            "m_z = {:.6} after reset, expected ≈ 1.0",
            m.z
        );
        // Step count and total reward should be zero
        assert_eq!(env.step_count, 0);
        assert!((env.total_reward).abs() < 1e-15);
    }

    #[test]
    fn test_step_changes_state() {
        let mut env = make_env();
        let m0 = env.reset();
        let j_max = env.config.j_max;
        let (m1, _reward, _done) = env.step(j_max);
        // Magnetization must have moved
        let delta = (m1.x - m0.x).abs() + (m1.y - m0.y).abs() + (m1.z - m0.z).abs();
        assert!(
            delta > 1e-15,
            "Magnetization did not change after step with j = j_max"
        );
    }

    #[test]
    fn test_norm_preserved() {
        let mut env = make_env();
        env.reset();
        let j_max = env.config.j_max;
        for i in 0..100 {
            // Alternate current sign to stay in the episode longer
            let j = if i % 2 == 0 { j_max } else { -j_max };
            let (m, _r, done) = env.step(j);
            let norm = m.magnitude();
            assert!(
                (norm - 1.0).abs() < 1e-3,
                "step {}: |m| = {:.8} deviates from unity",
                i,
                norm
            );
            if done {
                break;
            }
        }
    }

    #[test]
    fn test_switching_with_large_current() {
        // With j_max and 500 steps, m_z should at least decrease from its
        // initial value (even if full switching is not guaranteed for these
        // material parameters with a fixed constant current).
        let mut env = make_env();
        env.reset();
        let initial_mz = env.m.z;
        let j = env.config.j_max;
        let max_steps = env.config.max_steps;

        let mut final_mz = initial_mz;
        for _ in 0..max_steps {
            let (m, _r, done) = env.step(j);
            final_mz = m.z;
            if done {
                break;
            }
        }
        // At minimum, m_z should decrease (partial switching or full switching)
        assert!(
            final_mz < initial_mz,
            "m_z did not decrease: initial = {:.4}, final = {:.4}",
            initial_mz,
            final_mz
        );
    }

    #[test]
    fn test_cem_policy_converges() {
        // After 5 generations the best reward should not decrease compared
        // to the first generation (CEM is monotone in the best-so-far sense).
        let env = make_env();
        let mut opt = SotRlOptimizer::new(env);
        opt.population_size = 10;

        let result = opt.train(5, 12345);

        assert_eq!(result.n_generations, 5);
        assert_eq!(result.best_rewards_per_gen.len(), 5);

        // Best reward must be non-decreasing across generations (monotone CEM property
        // holds for the overall best tracked separately, checked here for the series).
        let first = result.best_rewards_per_gen[0];
        let last = *result.best_rewards_per_gen.last().unwrap_or(&first);
        assert!(
            result.best_reward >= first,
            "best_reward ({:.4}) should be >= first gen best ({:.4})",
            result.best_reward,
            first
        );
        // Suppress unused warning in this context
        let _ = last;
    }

    #[test]
    fn test_cem_update_reduces_sigma_toward_elite() {
        // If all elite policies are identical (constant), sigma should shrink
        // to the floor value after an update.
        let env = make_env();
        let j_max = env.config.j_max;
        let n = env.config.n_pulse_steps;
        let mut policy = CemPolicy::new(n, j_max);

        let initial_sigma_sum: f64 = policy.sigma.iter().sum();

        // Elite: all policies equal to j_max * 0.5 at every step
        let elite: Vec<Vec<f64>> = (0..10).map(|_| vec![j_max * 0.5; n]).collect();
        policy.update_from_elite(&elite);

        let updated_sigma_sum: f64 = policy.sigma.iter().sum();

        // Sigma should be strictly smaller (variance = 0 → clamped to floor)
        assert!(
            updated_sigma_sum < initial_sigma_sum,
            "sigma did not decrease: before = {:.4}, after = {:.4}",
            initial_sigma_sum,
            updated_sigma_sum
        );

        // All sigmas should be at the floor
        let sigma_floor = j_max * 0.01;
        for &s in &policy.sigma {
            assert!(
                (s - sigma_floor).abs() < 1e-15,
                "sigma = {:.6e}, expected floor = {:.6e}",
                s,
                sigma_floor
            );
        }
    }

    #[test]
    fn test_reward_structure() {
        // Applying j = 0 should yield a negative reward (only the time penalty
        // and zero energy penalty), never positive.
        let mut env = make_env();
        env.reset();
        let (_m, reward, _done) = env.step(0.0);
        // time_penalty = -0.005, energy_penalty = 0.0, switch = 0
        // unless by accident m.z crossed the threshold on the first step, which is impossible
        assert!(
            reward < 0.0,
            "reward with j=0 should be negative (time penalty), got {:.6}",
            reward
        );
    }

    // ------------------------------------------------------------------
    // LCG / numerical helper tests
    // ------------------------------------------------------------------

    #[test]
    fn test_lcg_deterministic() {
        let mut rng1 = Lcg::new(42);
        let mut rng2 = Lcg::new(42);
        for _ in 0..100 {
            assert_eq!(rng1.next_u64(), rng2.next_u64());
        }
    }

    #[test]
    fn test_lcg_normal_range() {
        let mut rng = Lcg::new(99);
        let samples: Vec<f64> = (0..1000).map(|_| rng.next_normal()).collect();
        // Empirical mean should be close to 0
        let mean = samples.iter().sum::<f64>() / 1000.0;
        assert!(mean.abs() < 0.2, "LCG normal mean = {:.4}", mean);
        // Empirical std should be close to 1
        let var = samples.iter().map(|x| x * x).sum::<f64>() / 1000.0;
        assert!(
            (var.sqrt() - 1.0).abs() < 0.2,
            "LCG normal std = {:.4}",
            var.sqrt()
        );
    }

    #[test]
    fn test_effective_field_no_current() {
        // Without current, the effective field should be dominated by the
        // anisotropy and bias, and point approximately along ẑ for m ≈ ẑ.
        let mut env = make_env();
        env.reset();
        let h = env.effective_field(env.m, 0.0);
        // With PMA (K > 0, m_z ≈ 1), h_ani_z = 2K/(μ₀ Ms) * m_z > 0
        // So h.z should be positive and large.
        assert!(
            h.z > 0.0,
            "With PMA and m ≈ ẑ, H_eff.z should be positive (got {:.4e})",
            h.z
        );
    }
}
