//! Magnetic structure optimisation via differentiable energy minimisation
//! (v0.7.0).
//!
//! Given a Heisenberg-style spin Hamiltonian, this module finds magnetic
//! ground states by gradient descent in spherical coordinates `(θ_i, φ_i)`
//! per site.  Choosing spherical coordinates rather than Cartesian
//! components bakes the unit-norm constraint into the parameterisation —
//! there is no need for Lagrange multipliers or projection steps.
//!
//! ## Energy functional
//!
//! ```text
//!     H = − J Σ_<ij> m_i · m_j
//!         − K_u Σ_i (m_i · ẑ)²
//!         − μ_0 M_s Σ_i m_i · H_ext
//! ```
//!
//! All three terms are recorded on the AD tape, so the gradient of the total
//! energy with respect to every `(θ_i, φ_i)` is obtained in a single backward
//! pass.  Any [`OptimizerKind`] from [`crate::autodiff::optimizer`] can be
//! plugged in.
//!
//! ## Pre-built bond topologies
//!
//! Two convenience constructors generate common 1-D bond lists:
//!
//! | Constructor               | Topology |
//! |--------------------------|----------|
//! | [`EnergyFunctional::with_chain_bonds`] | Open nearest-neighbour chain `0–1–2–…–(N−1)`            |
//! | [`EnergyFunctional::with_ring_bonds`]  | Periodic ring, with the extra bond `(N−1, 0)`           |
//!
//! Arbitrary user-defined bond lists are also accepted via
//! [`EnergyFunctional::with_bonds`].
//!
//! ## References
//!
//! - W. Heisenberg, "Zur Theorie des Ferromagnetismus", *Z. Phys.* **49**,
//!   619 (1928).
//! - L. Néel, "Propriétés magnétiques des ferrites; ferrimagnétisme et
//!   antiferromagnétisme", *Ann. Phys. (Paris)* **3**, 137 (1948).
//! - U. Nowak, "Classical Spin Models", in *Handbook of Magnetism and
//!   Advanced Magnetic Materials*, Wiley (2007).

use crate::autodiff::optimizer::{Adam, LBfgs, Optimizer, OptimizerKind, Sgd};
use crate::autodiff::tape::{Tape, Var};
use crate::constants::MU_0;
use crate::error::{dimension_mismatch, invalid_param, Result};
use crate::vector3::Vector3;

// ─── SpinConfig ─────────────────────────────────────────────────────────────

/// A configuration of `N` classical Heisenberg spins in spherical coordinates.
///
/// `theta_i ∈ [0, π]` is the polar angle measured from the `+ẑ` axis and
/// `phi_i ∈ [0, 2π)` is the azimuthal angle measured from `+x̂`.  The
/// corresponding unit vector is
///
/// `m_i = (sin θ_i cos φ_i, sin θ_i sin φ_i, cos θ_i)`.
#[derive(Debug, Clone)]
pub struct SpinConfig {
    /// Polar angles, length = N.
    pub theta: Vec<f64>,
    /// Azimuthal angles, length = N.
    pub phi: Vec<f64>,
}

impl SpinConfig {
    /// Create a configuration of `n` spins all pointing along `+ẑ`.
    pub fn new(n: usize) -> Self {
        Self {
            theta: vec![0.0; n],
            phi: vec![0.0; n],
        }
    }

    /// Random orientation per spin using the same LCG scheme as
    /// [`crate::autodiff::neural`] so seeds remain reproducible.
    pub fn random(n: usize, rng_seed: u64) -> Self {
        let mut state: u64 = if rng_seed == 0 {
            0xDEAD_BEEF_CAFE_BABE
        } else {
            rng_seed
        };
        let mut draw_unit = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let bits = state >> 11;
            (bits as f64) / ((1u64 << 53) as f64)
        };
        let mut theta = Vec::with_capacity(n);
        let mut phi = Vec::with_capacity(n);
        for _ in 0..n {
            // Sample θ via the inverse-CDF for uniform surface measure:
            // cos θ ~ Uniform(-1, 1)  ⇒  θ = acos(1 - 2u)
            let u = draw_unit();
            let v = draw_unit();
            theta.push((1.0 - 2.0 * u).acos());
            phi.push(2.0 * std::f64::consts::PI * v);
        }
        Self { theta, phi }
    }

    /// Convert each `(θ_i, φ_i)` to a Cartesian unit vector.
    pub fn to_vectors(&self) -> Vec<Vector3<f64>> {
        self.theta
            .iter()
            .zip(self.phi.iter())
            .map(|(&t, &p)| Vector3::new(t.sin() * p.cos(), t.sin() * p.sin(), t.cos()))
            .collect()
    }

    /// Inverse of [`SpinConfig::to_vectors`].  Each input vector is treated
    /// as a direction; its norm is ignored.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] if any vector has
    /// numerically vanishing magnitude.
    pub fn from_vectors(v: &[Vector3<f64>]) -> Result<Self> {
        let mut theta = Vec::with_capacity(v.len());
        let mut phi = Vec::with_capacity(v.len());
        for vi in v {
            let mag = vi.magnitude();
            if mag < 1e-12 {
                return Err(invalid_param(
                    "spin vector",
                    "zero-magnitude vector cannot be converted to (θ, φ)",
                ));
            }
            let nz = vi.z / mag;
            let nz_clamped = nz.clamp(-1.0, 1.0);
            theta.push(nz_clamped.acos());
            phi.push(vi.y.atan2(vi.x));
        }
        Ok(Self { theta, phi })
    }

    /// Flatten as `[θ_0, φ_0, θ_1, φ_1, …]`.
    pub fn params_flat(&self) -> Vec<f64> {
        let n = self.theta.len();
        let mut v = Vec::with_capacity(2 * n);
        for i in 0..n {
            v.push(self.theta[i]);
            v.push(self.phi[i]);
        }
        v
    }

    /// Inverse of [`SpinConfig::params_flat`].
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the length is wrong.
    pub fn set_params(&mut self, flat: &[f64]) -> Result<()> {
        let n = self.theta.len();
        if flat.len() != 2 * n {
            return Err(dimension_mismatch(
                &format!("{} params", 2 * n),
                &format!("{} params", flat.len()),
            ));
        }
        for i in 0..n {
            self.theta[i] = flat[2 * i];
            self.phi[i] = flat[2 * i + 1];
        }
        Ok(())
    }

    /// Number of spins in the configuration.
    pub fn n_spins(&self) -> usize {
        self.theta.len()
    }
}

// ─── EnergyFunctional ───────────────────────────────────────────────────────

/// Heisenberg-style energy functional whose total energy is differentiable
/// w.r.t. every spherical angle.
#[derive(Debug, Clone)]
pub struct EnergyFunctional {
    /// Heisenberg exchange `J` per bond (J/atom). Positive ⇒ ferromagnetic.
    pub exchange_j: f64,
    /// Uniaxial anisotropy constant `K_u` along ẑ (J/atom).
    pub anisotropy_k: f64,
    /// External field vector (A/m).
    pub h_ext: [f64; 3],
    /// Saturation magnetisation `M_s` (A/m).
    pub ms: f64,
    /// Exchange bond list `(i, j)` with `i < j`.
    pub bonds: Vec<(usize, usize)>,
}

impl EnergyFunctional {
    /// Construct an empty (no-bond) energy functional.  Use one of the
    /// `with_*` builders to populate the bond list.
    pub fn new(exchange_j: f64, anisotropy_k: f64, h_ext: [f64; 3], ms: f64) -> Self {
        Self {
            exchange_j,
            anisotropy_k,
            h_ext,
            ms,
            bonds: Vec::new(),
        }
    }

    /// Replace the bond list with `bonds`.
    pub fn with_bonds(mut self, bonds: Vec<(usize, usize)>) -> Self {
        self.bonds = bonds;
        self
    }

    /// Add open 1-D chain bonds `0–1–2–…–(N−1)`.
    pub fn with_chain_bonds(mut self, n_spins: usize) -> Self {
        let mut bonds = Vec::with_capacity(n_spins.saturating_sub(1));
        for i in 0..n_spins.saturating_sub(1) {
            bonds.push((i, i + 1));
        }
        self.bonds = bonds;
        self
    }

    /// Add periodic 1-D ring bonds (chain + closing bond `(N−1, 0)`).
    pub fn with_ring_bonds(mut self, n_spins: usize) -> Self {
        let mut bonds = Vec::with_capacity(n_spins);
        for i in 0..n_spins.saturating_sub(1) {
            bonds.push((i, i + 1));
        }
        if n_spins >= 2 {
            bonds.push((n_spins - 1, 0));
        }
        self.bonds = bonds;
        self
    }

    /// Total energy as a [`Var<'t>`] on the tape.  `theta_phi` is the flat
    /// parameter vector `[θ_0, φ_0, θ_1, φ_1, …]`.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the length is
    /// not even or any bond references a missing spin.
    pub fn energy<'t>(&self, tape: &'t Tape, theta_phi: &[Var<'t>]) -> Result<Var<'t>> {
        if theta_phi.len() % 2 != 0 {
            return Err(invalid_param(
                "theta_phi",
                "length must be even ([θ_0, φ_0, θ_1, φ_1, …])",
            ));
        }
        let n = theta_phi.len() / 2;
        // Precompute per-spin (m_x, m_y, m_z) as Var<'t>.
        let mut mx = Vec::with_capacity(n);
        let mut my = Vec::with_capacity(n);
        let mut mz = Vec::with_capacity(n);
        for i in 0..n {
            let t = theta_phi[2 * i];
            let p = theta_phi[2 * i + 1];
            let sin_t = t.sin();
            mx.push(sin_t * p.cos());
            my.push(sin_t * p.sin());
            mz.push(t.cos());
        }
        let mut energy = Var::leaf(tape, 0.0);

        // Exchange:  - J Σ_<ij> m_i · m_j
        if self.exchange_j != 0.0 && !self.bonds.is_empty() {
            let mut ex_sum = Var::leaf(tape, 0.0);
            for &(i, j) in &self.bonds {
                if i >= n || j >= n {
                    return Err(dimension_mismatch(
                        &format!("indices < {n}"),
                        &format!("bond ({i}, {j})"),
                    ));
                }
                let dot = mx[i] * mx[j] + my[i] * my[j] + mz[i] * mz[j];
                ex_sum = ex_sum + dot;
            }
            energy = energy + ex_sum * (-self.exchange_j);
        }

        // Uniaxial anisotropy:  - K_u Σ (m_i^z)²
        if self.anisotropy_k != 0.0 {
            let mut aniso_sum = Var::leaf(tape, 0.0);
            for mz_i in mz.iter().take(n) {
                aniso_sum = aniso_sum + *mz_i * *mz_i;
            }
            energy = energy + aniso_sum * (-self.anisotropy_k);
        }

        // Zeeman:  - μ_0 M_s Σ m_i · H_ext
        let zeeman_coeff = -MU_0 * self.ms;
        if zeeman_coeff != 0.0
            && (self.h_ext[0] != 0.0 || self.h_ext[1] != 0.0 || self.h_ext[2] != 0.0)
        {
            let mut zee_sum = Var::leaf(tape, 0.0);
            for i in 0..n {
                let dot = mx[i] * self.h_ext[0] + my[i] * self.h_ext[1] + mz[i] * self.h_ext[2];
                zee_sum = zee_sum + dot;
            }
            energy = energy + zee_sum * zeeman_coeff;
        }

        Ok(energy)
    }

    /// Plain-`f64` evaluation of the total energy (no tape, no gradients).
    pub fn energy_f64(&self, config: &SpinConfig) -> f64 {
        let n = config.n_spins();
        let vs: Vec<Vector3<f64>> = config.to_vectors();
        let mut e = 0.0_f64;
        if self.exchange_j != 0.0 {
            let mut s = 0.0_f64;
            for &(i, j) in &self.bonds {
                if i < n && j < n {
                    s += vs[i].dot(&vs[j]);
                }
            }
            e -= self.exchange_j * s;
        }
        if self.anisotropy_k != 0.0 {
            let mut s = 0.0_f64;
            for vi in &vs {
                s += vi.z * vi.z;
            }
            e -= self.anisotropy_k * s;
        }
        let zee_coeff = -MU_0 * self.ms;
        if zee_coeff != 0.0 {
            let mut s = 0.0_f64;
            for vi in &vs {
                s += vi.x * self.h_ext[0] + vi.y * self.h_ext[1] + vi.z * self.h_ext[2];
            }
            e += zee_coeff * s;
        }
        e
    }
}

// ─── StructureOptResult ─────────────────────────────────────────────────────

/// Result of [`MagneticStructureOptimizer::optimize`].
#[derive(Debug, Clone)]
pub struct StructureOptResult {
    /// Final spin configuration after optimisation.
    pub final_config: SpinConfig,
    /// Final total energy.
    pub final_energy: f64,
    /// Number of gradient steps actually taken.
    pub n_iterations: usize,
    /// Whether the convergence criterion was met.
    pub converged: bool,
    /// Per-iteration total-energy history.
    pub energy_history: Vec<f64>,
}

// ─── MagneticStructureOptimizer ─────────────────────────────────────────────

/// High-level driver that minimises an [`EnergyFunctional`] over a
/// [`SpinConfig`].
pub struct MagneticStructureOptimizer {
    /// Starting spin configuration (will be cloned on each `optimize` call).
    pub initial: SpinConfig,
    /// Energy functional to minimise.
    pub energy: EnergyFunctional,
    /// Choice of inner gradient-based optimiser.
    pub optimizer_kind: OptimizerKind,
    /// Learning rate forwarded to the optimiser.
    pub lr: f64,
    /// Maximum number of gradient steps.
    pub max_iter: usize,
    /// Convergence tolerance on |ΔE|.
    pub tol: f64,
}

impl MagneticStructureOptimizer {
    /// Construct an optimiser with sensible defaults: Adam, `lr = 0.05`,
    /// `max_iter = 500`, `tol = 1e-12`.
    pub fn new(initial: SpinConfig, energy: EnergyFunctional) -> Self {
        Self {
            initial,
            energy,
            optimizer_kind: OptimizerKind::Adam,
            lr: 0.05,
            max_iter: 500,
            tol: 1e-12,
        }
    }

    /// Override the optimiser kind and learning rate.
    pub fn with_optimizer(mut self, kind: OptimizerKind, lr: f64) -> Self {
        self.optimizer_kind = kind;
        self.lr = lr;
        self
    }

    /// Override `max_iter`.
    pub fn with_max_iter(mut self, n: usize) -> Self {
        self.max_iter = n;
        self
    }

    /// Override `tol`.
    pub fn with_tol(mut self, tol: f64) -> Self {
        self.tol = tol;
        self
    }

    /// Run the optimisation loop.
    ///
    /// # Errors
    /// Propagates errors from [`EnergyFunctional::energy`].
    pub fn optimize(&mut self) -> Result<StructureOptResult> {
        let n = self.initial.n_spins();
        let mut params = self.initial.params_flat();
        let n_p = params.len();
        let mut opt: Box<dyn Optimizer> = match self.optimizer_kind {
            OptimizerKind::Sgd { momentum } => Box::new(
                Sgd::new(self.lr, momentum)
                    .map_err(|_| invalid_param("sgd", "invalid lr/momentum"))?,
            ),
            OptimizerKind::Adam => {
                let mut a = Adam::default_params(n_p);
                a.lr = self.lr;
                Box::new(a)
            },
            OptimizerKind::LBfgs { history } => Box::new(
                LBfgs::new(self.lr, history)
                    .map_err(|_| invalid_param("lbfgs", "invalid lr/history"))?,
            ),
        };
        let mut energy_history = Vec::with_capacity(self.max_iter);
        let mut prev_e = f64::INFINITY;
        let mut converged = false;

        for _ in 0..self.max_iter {
            let tape = Tape::new();
            let leaves: Vec<Var<'_>> = params.iter().map(|&p| Var::leaf(&tape, p)).collect();
            let energy_var = self.energy.energy(&tape, &leaves)?;
            let energy_val = energy_var.value();
            tape.backward(energy_var);
            let grads: Vec<f64> = leaves.iter().map(|l| l.grad()).collect();
            energy_history.push(energy_val);
            opt.step(&mut params, &grads);
            if (energy_val - prev_e).abs() < self.tol {
                converged = true;
                break;
            }
            prev_e = energy_val;
        }

        let mut final_config = SpinConfig::new(n);
        final_config.set_params(&params)?;
        let final_energy = *energy_history.last().unwrap_or(&f64::NAN);
        let n_iterations = energy_history.len();
        Ok(StructureOptResult {
            final_config,
            final_energy,
            n_iterations,
            converged,
            energy_history,
        })
    }
}

// ─── Convenience preset functions ───────────────────────────────────────────

/// Find the ferromagnetic ground state of an `n_spins`-site Heisenberg chain
/// with positive exchange and a field along `+ẑ`.
///
/// # Errors
/// Propagates errors from the inner optimiser.
pub fn find_fm_ground_state(
    n_spins: usize,
    exchange_j: f64,
    h_ext_z: f64,
    ms: f64,
) -> Result<SpinConfig> {
    if n_spins == 0 {
        return Err(invalid_param("n_spins", "must be positive"));
    }
    let initial = SpinConfig::random(n_spins, 0xABCD_1234);
    let energy =
        EnergyFunctional::new(exchange_j, 0.0, [0.0, 0.0, h_ext_z], ms).with_chain_bonds(n_spins);
    let mut opt = MagneticStructureOptimizer::new(initial, energy)
        .with_optimizer(OptimizerKind::Adam, 0.05)
        .with_max_iter(800)
        .with_tol(1e-14);
    let res = opt.optimize()?;
    Ok(res.final_config)
}

/// Find the antiferromagnetic Néel ground state of an `n_spins`-site
/// Heisenberg chain with **negative** exchange (`exchange_j_afm > 0` in this
/// convention: the value supplied is the magnitude and is converted to a
/// negative `J` internally so the resulting Hamiltonian favours alternation).
///
/// # Errors
/// Propagates errors from the inner optimiser.
pub fn find_afm_ground_state(n_spins: usize, exchange_j_afm: f64) -> Result<SpinConfig> {
    if n_spins == 0 {
        return Err(invalid_param("n_spins", "must be positive"));
    }
    if exchange_j_afm <= 0.0 {
        return Err(invalid_param(
            "exchange_j_afm",
            "magnitude must be positive",
        ));
    }
    // Seed the initial configuration with a Néel pattern so the optimiser
    // locks into the global minimum rather than a domain-walled local one.
    let mut initial = SpinConfig::new(n_spins);
    for i in 0..n_spins {
        initial.theta[i] = if i % 2 == 0 {
            0.0
        } else {
            std::f64::consts::PI
        };
        initial.phi[i] = 0.0;
    }
    let energy =
        EnergyFunctional::new(-exchange_j_afm, 0.0, [0.0, 0.0, 0.0], 0.0).with_chain_bonds(n_spins);
    let mut opt = MagneticStructureOptimizer::new(initial, energy)
        .with_optimizer(OptimizerKind::Adam, 0.02)
        .with_max_iter(400)
        .with_tol(1e-14);
    let res = opt.optimize()?;
    Ok(res.final_config)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::autodiff::tape::Tape;

    // 1. SpinConfig::to_vectors then from_vectors is the identity (up to atan2
    //    sign conventions).
    #[test]
    fn test_spin_config_roundtrip_vectors() {
        let cfg = SpinConfig::random(8, 99);
        let vs = cfg.to_vectors();
        let cfg2 = SpinConfig::from_vectors(&vs).unwrap();
        let vs2 = cfg2.to_vectors();
        for (v1, v2) in vs.iter().zip(vs2.iter()).take(cfg.n_spins()) {
            assert!((v1.x - v2.x).abs() < 1e-10);
            assert!((v1.y - v2.y).abs() < 1e-10);
            assert!((v1.z - v2.z).abs() < 1e-10);
        }
    }

    // 2. from_vectors rejects zero-magnitude inputs.
    #[test]
    fn test_from_vectors_zero_rejected() {
        let v = vec![Vector3::new(0.0_f64, 0.0, 0.0)];
        assert!(SpinConfig::from_vectors(&v).is_err());
    }

    // 3. Heisenberg sign convention: two parallel spins with J > 0 give a
    //    negative contribution (-J), antiparallel give +J.
    #[test]
    fn test_heisenberg_sign() {
        // Two spins, both along +z.
        let cfg_parallel = SpinConfig::new(2);
        let mut cfg_anti = SpinConfig::new(2);
        cfg_anti.theta[1] = std::f64::consts::PI;
        let ef = EnergyFunctional::new(1.0, 0.0, [0.0, 0.0, 0.0], 0.0).with_chain_bonds(2);
        let e_para = ef.energy_f64(&cfg_parallel);
        let e_anti = ef.energy_f64(&cfg_anti);
        // Parallel: m_i·m_j = +1 → energy = -1.  Antiparallel → +1.
        assert!((e_para - (-1.0)).abs() < 1e-12);
        assert!((e_anti - 1.0).abs() < 1e-12);
    }

    // 4. Finite-difference vs AD gradient on the total energy for a small
    //    random configuration.
    #[test]
    fn test_energy_gradient_finite_diff() {
        let cfg = SpinConfig::random(3, 17);
        let ef = EnergyFunctional::new(0.5, 0.1, [0.0, 0.0, 1e3], 1e5).with_chain_bonds(3);
        let params = cfg.params_flat();
        // AD gradient.
        let tape = Tape::new();
        let leaves: Vec<Var<'_>> = params.iter().map(|&p| Var::leaf(&tape, p)).collect();
        let e_var = ef.energy(&tape, &leaves).unwrap();
        tape.backward(e_var);
        let ad: Vec<f64> = leaves.iter().map(|l| l.grad()).collect();
        // FD gradient.
        let h = 1e-6;
        for k in 0..params.len() {
            let mut p_plus = params.clone();
            let mut p_minus = params.clone();
            p_plus[k] += h;
            p_minus[k] -= h;
            let mut cfg_p = cfg.clone();
            let mut cfg_m = cfg.clone();
            cfg_p.set_params(&p_plus).unwrap();
            cfg_m.set_params(&p_minus).unwrap();
            let fd = (ef.energy_f64(&cfg_p) - ef.energy_f64(&cfg_m)) / (2.0 * h);
            assert!(
                (ad[k] - fd).abs() < 1e-4,
                "k={} AD {} vs FD {}",
                k,
                ad[k],
                fd
            );
        }
    }

    // 5. FM ground state: after optimisation, all spins are nearly parallel.
    #[test]
    fn test_fm_ground_state_parallel() {
        let cfg = find_fm_ground_state(4, 1.0, 1e6, 1e5).unwrap();
        let vs = cfg.to_vectors();
        // Compare every pair of spins via dot product; should all be ≈ +1.
        for i in 1..vs.len() {
            let dot = vs[0].dot(&vs[i]);
            assert!(
                dot > 0.95,
                "spin {} not parallel to spin 0 (dot = {})",
                i,
                dot
            );
        }
    }

    // 6. AFM ground state alternates: consecutive dot products ≈ -1.
    #[test]
    fn test_afm_ground_state_alternates() {
        let cfg = find_afm_ground_state(4, 1.0).unwrap();
        let vs = cfg.to_vectors();
        for i in 0..vs.len() - 1 {
            let dot = vs[i].dot(&vs[i + 1]);
            assert!(
                dot < -0.95,
                "bond ({}, {}) not antiparallel (dot = {})",
                i,
                i + 1,
                dot
            );
        }
    }

    // 7. Uniaxial anisotropy aligns m with ẑ when started slightly off the
    //    saddle at θ = π/2 (the in-plane saddle is unstable but has zero
    //    gradient, so a tiny perturbation is needed to break the symmetry).
    #[test]
    fn test_uniaxial_anisotropy_aligns_z() {
        let mut initial = SpinConfig::new(1);
        // Perturb away from the unstable saddle at θ = π/2.
        initial.theta[0] = std::f64::consts::FRAC_PI_2 - 0.1;
        let ef = EnergyFunctional::new(0.0, 1.0, [0.0, 0.0, 0.0], 0.0);
        let mut opt = MagneticStructureOptimizer::new(initial, ef)
            .with_optimizer(OptimizerKind::Adam, 0.05)
            .with_max_iter(2000)
            .with_tol(1e-15);
        let res = opt.optimize().unwrap();
        let v = res.final_config.to_vectors()[0];
        // Either +z or -z is a valid minimum (degenerate).
        assert!(v.z.abs() > 0.95, "|m_z| = {} not near 1", v.z.abs());
    }

    // 8. Zeeman energy minimised when m ∥ H.
    #[test]
    fn test_zeeman_minimised_when_parallel() {
        let mut initial = SpinConfig::new(1);
        initial.theta[0] = 0.0; // along +z
        let mut initial_anti = SpinConfig::new(1);
        initial_anti.theta[0] = std::f64::consts::PI; // along -z
        let ef = EnergyFunctional::new(0.0, 0.0, [0.0, 0.0, 1e6], 1e5);
        let e_para = ef.energy_f64(&initial);
        let e_anti = ef.energy_f64(&initial_anti);
        assert!(e_para < e_anti);
    }

    // 9. Ring topology has exactly `N` bonds (one extra closing bond).
    #[test]
    fn test_ring_topology_bond_count() {
        let ef = EnergyFunctional::new(1.0, 0.0, [0.0, 0.0, 0.0], 0.0).with_ring_bonds(5);
        assert_eq!(ef.bonds.len(), 5);
        let chain = EnergyFunctional::new(1.0, 0.0, [0.0, 0.0, 0.0], 0.0).with_chain_bonds(5);
        assert_eq!(chain.bonds.len(), 4);
    }

    // 10. Convergence within tolerance terminates early.
    #[test]
    fn test_convergence_within_tol() {
        // Trivial: a single spin in a field — should converge in very few
        // iterations and report converged = true.
        let mut initial = SpinConfig::new(1);
        initial.theta[0] = 0.0;
        let ef = EnergyFunctional::new(0.0, 1.0, [0.0, 0.0, 1e3], 1e3);
        let mut opt = MagneticStructureOptimizer::new(initial, ef)
            .with_optimizer(OptimizerKind::Adam, 0.05)
            .with_max_iter(1000)
            .with_tol(1e-10);
        let res = opt.optimize().unwrap();
        assert!(res.converged);
        assert!(res.n_iterations < 1000);
    }

    // 11. params_flat/set_params round-trip preserves SpinConfig.
    #[test]
    fn test_spin_config_params_roundtrip() {
        let cfg = SpinConfig::random(5, 31);
        let flat = cfg.params_flat();
        let mut cfg2 = SpinConfig::new(5);
        cfg2.set_params(&flat).unwrap();
        for i in 0..5 {
            assert!((cfg.theta[i] - cfg2.theta[i]).abs() < 1e-15);
            assert!((cfg.phi[i] - cfg2.phi[i]).abs() < 1e-15);
        }
    }
}
