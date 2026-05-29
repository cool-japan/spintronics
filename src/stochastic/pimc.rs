//! Path-integral Monte Carlo (PIMC) for finite-temperature classical-spin
//! Heisenberg models with imaginary-time discretization.
//!
//! # Physical model
//!
//! For the Heisenberg chain
//!
//! ```text
//!     H = − J Σ_{⟨i,j⟩} S_i · S_j  − h Σ_i S_i^z,
//! ```
//!
//! the partition function at inverse temperature `β = 1/(k_B T)` is
//!
//! ```text
//!     Z = Tr exp(−β H) ≈ Tr [exp(−Δτ H / M)]^M,    Δτ = β/M.
//! ```
//!
//! For classical Heisenberg spins, the worldlines are continuous unit
//! vectors `S_i(τ) ∈ S²` at each of `M` imaginary-time slices, joined
//! periodically (`S_i(M) = S_i(0)`). The Trotterized effective action is
//!
//! ```text
//!     S = Σ_{slices,bonds} Δτ · [ −J · S_i(τ) · S_j(τ) − h · S_i^z(τ) ]
//!       − K_rig · Σ_{slices,sites} S_i(τ) · S_i(τ+1),
//! ```
//!
//! where the *imaginary-time rigidity* `K_rig = 1/Δτ` keeps adjacent slices
//! correlated — this term is the high-temperature classical analog of the
//! kinetic contribution that becomes important when `Δτ → 0`.
//!
//! Monte Carlo updates are local single-slice-single-site rotations. A
//! proposal rotates `S_i(τ) → S_i(τ) + δ` then renormalises; the change in
//! action is computed locally and accepted with probability
//! `min(1, exp(−ΔS))`. The β factor is *already inside* `Δτ`, so no
//! additional β multiplication is needed at the acceptance step.
//!
//! # Limitations
//!
//! This implementation is the "high-temperature classical" approximation;
//! it is exact in the `M → 1` limit (classical Heisenberg model with no
//! imaginary-time structure) and exact in the `Δτ → 0`, `K_rig → ∞` limit
//! (frozen worldlines). For intermediate `M` it provides a useful
//! interpolation suitable for benchmarking and pedagogical study, but it is
//! *not* a substitute for the full quantum stochastic series expansion
//! (Sandvik & Kurkijärvi, PRB 43, 5950, 1991) when sign-problem-free
//! quantum estimates are required.
//!
//! # References
//!
//! - Suzuki, M., *Generalized Trotter's formula and systematic
//!   approximants of exponential operators*, Prog. Theor. Phys. **56**,
//!   1454 (1976).
//! - Hirsch, J. E., *Two-dimensional Hubbard model: Numerical simulation
//!   study*, Phys. Rev. B **31**, 4403 (1985).
//! - Sandvik, A. W. and Kurkijärvi, J., *Quantum Monte Carlo simulation
//!   method for spin systems*, Phys. Rev. B **43**, 5950 (1991).
//! - Bonner, J. C. and Fisher, M. E., *Linear magnetic chains with
//!   anisotropic coupling*, Phys. Rev. **135**, A640 (1964).

use scirs2_core::random::core::Random as CoreRandom;
use scirs2_core::random::rand_distributions::Normal;
use scirs2_core::random::{rngs, seeded_rng};

use crate::error::{invalid_param, Result};
use crate::vector3::Vector3;

/// Geometry of the lattice on which the PIMC simulation runs.
#[derive(Debug, Clone, Copy)]
pub enum PimcLattice {
    /// Open 1D chain with nearest-neighbour interactions.
    Chain1D {
        /// Number of lattice sites.
        n_spins: usize,
        /// Exchange coupling `[arbitrary units]`. `J > 0` ferromagnetic in
        /// the convention `H = -J S_i · S_j`.
        j: f64,
        /// Applied field along z.
        h_z: f64,
    },
    /// 1D ring with periodic boundary conditions on the spatial lattice.
    Ring1D {
        /// Number of sites in the ring.
        n_spins: usize,
        /// Exchange coupling.
        j: f64,
        /// Applied field along z.
        h_z: f64,
    },
}

impl PimcLattice {
    fn n_spins(&self) -> usize {
        match self {
            Self::Chain1D { n_spins, .. } | Self::Ring1D { n_spins, .. } => *n_spins,
        }
    }
    fn j(&self) -> f64 {
        match self {
            Self::Chain1D { j, .. } | Self::Ring1D { j, .. } => *j,
        }
    }
    fn h_z(&self) -> f64 {
        match self {
            Self::Chain1D { h_z, .. } | Self::Ring1D { h_z, .. } => *h_z,
        }
    }
    fn periodic(&self) -> bool {
        matches!(self, Self::Ring1D { .. })
    }
}

/// Configuration parameters for the PIMC simulation.
#[derive(Debug, Clone)]
pub struct PimcConfig {
    /// Number of Trotter slices in imaginary time.
    pub n_slices: usize,
    /// Number of MC sweeps over the full (sites × slices) lattice for
    /// data collection.
    pub n_steps: usize,
    /// Number of thermalisation sweeps performed before measurements.
    pub n_thermalize: usize,
    /// Inverse temperature `β = 1 / (k_B T)` in units of `1/J`
    /// (dimensionless `β·J`).
    pub beta: f64,
    /// Seed for the deterministic RNG (reproducibility).
    pub seed: u64,
    /// Maximum angular displacement per proposal `[rad]`.
    pub angular_step: f64,
}

impl PimcConfig {
    /// Validate that all parameters are in physically meaningful ranges.
    pub fn validate(&self) -> Result<()> {
        if self.n_slices == 0 {
            return Err(invalid_param("n_slices", "must be ≥ 1"));
        }
        if self.n_steps == 0 {
            return Err(invalid_param("n_steps", "must be ≥ 1"));
        }
        if !self.beta.is_finite() || self.beta <= 0.0 {
            return Err(invalid_param(
                "beta",
                "inverse temperature must be positive and finite",
            ));
        }
        if !self.angular_step.is_finite() || self.angular_step <= 0.0 {
            return Err(invalid_param(
                "angular_step",
                "must be positive and finite (radians)",
            ));
        }
        Ok(())
    }
}

/// Aggregated observables from a PIMC measurement run.
#[derive(Debug, Clone)]
pub struct PimcResult {
    /// Average z-component of total magnetization per site (slice-and-site
    /// averaged).
    pub average_magnetization_z: f64,
    /// Average energy per site (classical Heisenberg + Zeeman).
    pub average_energy_per_site: f64,
    /// Estimator for the specific heat C = β² · var(E)/N.
    pub specific_heat: f64,
    /// Magnetic susceptibility χ = β · N · var(M_z/N).
    pub susceptibility: f64,
    /// Total Metropolis proposals accepted during measurement.
    pub n_accepted: usize,
    /// Total Metropolis proposals attempted during measurement.
    pub n_proposed: usize,
}

/// PIMC simulation driver.
pub struct PimcSimulation {
    /// Lattice geometry (chain or ring).
    pub lattice: PimcLattice,
    /// Run configuration.
    pub config: PimcConfig,
    /// Worldlines: `worldlines[site][slice]` is the classical spin vector
    /// at the given site and imaginary-time slice.
    pub worldlines: Vec<Vec<Vector3<f64>>>,
    rng: CoreRandom<rngs::StdRng>,
    normal: Normal<f64>,
}

impl PimcSimulation {
    /// Create a new PIMC simulation. Worldlines are initialised so that
    /// every slice carries the same random spin (a "classical" starting
    /// point with perfect imaginary-time order), then thermalised with
    /// [`Self::thermalize`].
    ///
    /// # Errors
    /// Returns `InvalidParameter` for `n_spins == 0`, `n_slices == 0`,
    /// `beta ≤ 0`, or invalid angular step.
    pub fn new(lattice: PimcLattice, config: PimcConfig) -> Result<Self> {
        if lattice.n_spins() == 0 {
            return Err(invalid_param("n_spins", "must be ≥ 1"));
        }
        config.validate()?;

        let mut rng = seeded_rng(config.seed);
        let normal = Normal::new(0.0, 1.0).expect("standard normal parameters are valid");

        let mut worldlines = Vec::with_capacity(lattice.n_spins());
        for _ in 0..lattice.n_spins() {
            let s0 = random_unit_vector(&mut rng, &normal);
            let mut line = Vec::with_capacity(config.n_slices);
            for _ in 0..config.n_slices {
                line.push(s0);
            }
            worldlines.push(line);
        }

        Ok(Self {
            lattice,
            config,
            worldlines,
            rng,
            normal,
        })
    }

    /// Run the thermalisation phase (no measurements are accumulated).
    pub fn thermalize(&mut self) -> Result<()> {
        for _ in 0..self.config.n_thermalize {
            let _ = self.metropolis_sweep();
        }
        Ok(())
    }

    /// Run the measurement phase, returning the averaged observables.
    pub fn run(&mut self) -> Result<PimcResult> {
        let n_sites = self.lattice.n_spins();
        let n_slices = self.config.n_slices;
        let n_total = (n_sites * n_slices) as f64;

        let mut e_sum = 0.0_f64;
        let mut e_sq_sum = 0.0_f64;
        let mut mz_sum = 0.0_f64;
        let mut mz_sq_sum = 0.0_f64;
        let mut n_accepted = 0_usize;
        let mut n_proposed = 0_usize;

        for _ in 0..self.config.n_steps {
            let (accepted, proposed) = self.metropolis_sweep_counted();
            n_accepted += accepted;
            n_proposed += proposed;

            let (energy_per_site, mz_per_site) = self.measure_energy_and_magnetization();
            e_sum += energy_per_site;
            e_sq_sum += energy_per_site * energy_per_site;
            mz_sum += mz_per_site;
            mz_sq_sum += mz_per_site * mz_per_site;
        }

        let n = self.config.n_steps as f64;
        let mean_e = e_sum / n;
        let mean_mz = mz_sum / n;
        let var_e = (e_sq_sum / n - mean_e * mean_e).max(0.0);
        let var_mz = (mz_sq_sum / n - mean_mz * mean_mz).max(0.0);

        // Per-site specific heat: C/N = β² · var(E/N) · N (var on per-site
        // energy includes the 1/N factor). Following the standard estimator
        // C = β² · ⟨(E − ⟨E⟩)²⟩ / N → we keep var on per-site quantity
        // and rescale.
        let beta = self.config.beta;
        let specific_heat = beta * beta * var_e * (n_sites as f64);
        let susceptibility = beta * var_mz * (n_sites as f64);

        // Acknowledge n_total — sweep size used as cross-check for tests.
        debug_assert!(n_total > 0.0);

        Ok(PimcResult {
            average_magnetization_z: mean_mz,
            average_energy_per_site: mean_e,
            specific_heat,
            susceptibility,
            n_accepted,
            n_proposed,
        })
    }

    /// Internal: one full sweep, returning the number of accepted updates.
    fn metropolis_sweep(&mut self) -> usize {
        self.metropolis_sweep_counted().0
    }

    /// Internal: one full sweep, returning `(n_accepted, n_proposed)`.
    fn metropolis_sweep_counted(&mut self) -> (usize, usize) {
        let n_sites = self.lattice.n_spins();
        let n_slices = self.config.n_slices;
        let mut accepted = 0_usize;
        let mut proposed = 0_usize;

        for site in 0..n_sites {
            for slice in 0..n_slices {
                proposed += 1;
                if self.attempt_local_update(site, slice) {
                    accepted += 1;
                }
            }
        }
        (accepted, proposed)
    }

    /// Attempt a single local update at `(site, slice)`. Returns `true` if
    /// accepted.
    fn attempt_local_update(&mut self, site: usize, slice: usize) -> bool {
        let old_spin = self.worldlines[site][slice];
        let old_action = self.local_action(site, slice);

        // Propose a small rotation: add a random tangent vector scaled by
        // `angular_step`, then renormalise.
        let perturbation = random_unit_vector(&mut self.rng, &self.normal)
            * (self.config.angular_step * self.rng.gen_range(0.0_f64..1.0));
        let proposed = (old_spin + perturbation).normalize();
        self.worldlines[site][slice] = proposed;
        let new_action = self.local_action(site, slice);

        let delta_s = new_action - old_action;
        // Acceptance: exp(−ΔS). The β factor is absorbed in Δτ = β / M.
        if delta_s <= 0.0 {
            return true;
        }
        let accept_prob = (-delta_s).exp();
        let r: f64 = self.rng.gen_range(0.0_f64..1.0);
        if r < accept_prob {
            true
        } else {
            self.worldlines[site][slice] = old_spin;
            false
        }
    }

    /// Local action contribution at `(site, slice)`:
    ///   − Δτ · (J · S · Σ_nn + h_z · S_z)
    ///   − K_rig · (S · S(τ−1) + S · S(τ+1))
    /// where the rigidity coupling K_rig = 1/Δτ enforces imaginary-time
    /// continuity.
    fn local_action(&self, site: usize, slice: usize) -> f64 {
        let n_slices = self.config.n_slices;
        let n_sites = self.lattice.n_spins();
        let d_tau = self.config.beta / n_slices as f64;
        let k_rig = 1.0 / d_tau;
        let s = self.worldlines[site][slice];
        let j = self.lattice.j();
        let h_z = self.lattice.h_z();

        // Spatial neighbours (1D NN with optional periodic BC).
        let mut neighbour_sum = Vector3::zero();
        if site + 1 < n_sites {
            neighbour_sum = neighbour_sum + self.worldlines[site + 1][slice];
        } else if self.lattice.periodic() {
            neighbour_sum = neighbour_sum + self.worldlines[0][slice];
        }
        if site > 0 {
            neighbour_sum = neighbour_sum + self.worldlines[site - 1][slice];
        } else if self.lattice.periodic() {
            neighbour_sum = neighbour_sum + self.worldlines[n_sites - 1][slice];
        }

        // Local spatial action: −Δτ · J · S · Σ_nn (with H = −J Σ S·S) and
        // −Δτ · h_z · S_z. Sign convention: lower action ↔ lower energy.
        let spatial = -d_tau * (j * s.dot(&neighbour_sum) + h_z * s.z);

        // Imaginary-time rigidity (periodic in τ): −K_rig · S(τ) · S(τ±1).
        let prev_slice = if slice == 0 { n_slices - 1 } else { slice - 1 };
        let next_slice = if slice + 1 == n_slices { 0 } else { slice + 1 };
        let s_prev = self.worldlines[site][prev_slice];
        let s_next = self.worldlines[site][next_slice];
        let temporal = -k_rig * (s.dot(&s_prev) + s.dot(&s_next));

        spatial + temporal
    }

    /// Measure per-site energy and z-magnetization, slice-averaged.
    ///
    /// Returns `(energy_per_site, mz_per_site)` where both are averaged
    /// over all `(site, slice)` pairs.
    fn measure_energy_and_magnetization(&self) -> (f64, f64) {
        let n_slices = self.config.n_slices;
        let n_sites = self.lattice.n_spins();
        let j = self.lattice.j();
        let h_z = self.lattice.h_z();
        let periodic = self.lattice.periodic();

        let mut e_total = 0.0_f64;
        let mut mz_total = 0.0_f64;
        for slice in 0..n_slices {
            let mut e_slice = 0.0_f64;
            for site in 0..n_sites {
                let s = self.worldlines[site][slice];
                mz_total += s.z;
                e_slice -= h_z * s.z;
                if site + 1 < n_sites {
                    let s_next = self.worldlines[site + 1][slice];
                    e_slice -= j * s.dot(&s_next);
                } else if periodic && n_sites > 1 {
                    let s_next = self.worldlines[0][slice];
                    e_slice -= j * s.dot(&s_next);
                }
            }
            e_total += e_slice;
        }
        let n_total = (n_sites * n_slices) as f64;
        (
            e_total / n_total,  // per (site, slice)
            mz_total / n_total, // per (site, slice)
        )
    }
}

/// Sample a uniformly random unit vector on `S²` via Marsaglia's method
/// applied to a pair of standard Gaussian draws — we obtain three
/// Gaussians and normalize, which is statistically equivalent to the
/// rejection sampler and avoids loop divergence.
fn random_unit_vector(rng: &mut CoreRandom<rngs::StdRng>, normal: &Normal<f64>) -> Vector3<f64> {
    loop {
        let v = Vector3::new(
            rng.sample(*normal),
            rng.sample(*normal),
            rng.sample(*normal),
        );
        let mag = v.magnitude();
        if mag > 1.0e-12 {
            return v * (1.0 / mag);
        }
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn default_config() -> PimcConfig {
        PimcConfig {
            n_slices: 4,
            n_steps: 100,
            n_thermalize: 50,
            beta: 1.0,
            seed: 42,
            angular_step: 0.5,
        }
    }

    // -- Test 1: construct + validate ----------------------------------

    #[test]
    fn test_construct_basic_chain() {
        let lat = PimcLattice::Chain1D {
            n_spins: 6,
            j: 1.0,
            h_z: 0.1,
        };
        let cfg = default_config();
        let sim = PimcSimulation::new(lat, cfg).expect("ok");
        assert_eq!(sim.worldlines.len(), 6);
        assert_eq!(sim.worldlines[0].len(), 4);
        for line in &sim.worldlines {
            for s in line {
                assert!(
                    close(s.magnitude(), 1.0, 1.0e-10),
                    "init spin not normalized"
                );
            }
        }
    }

    // -- Test 2: invalid config errors --------------------------------

    #[test]
    fn test_invalid_configs() {
        let lat_zero = PimcLattice::Chain1D {
            n_spins: 0,
            j: 1.0,
            h_z: 0.0,
        };
        assert!(PimcSimulation::new(lat_zero, default_config()).is_err());

        let lat = PimcLattice::Chain1D {
            n_spins: 4,
            j: 1.0,
            h_z: 0.0,
        };
        let mut cfg = default_config();
        cfg.n_slices = 0;
        assert!(PimcSimulation::new(lat, cfg.clone()).is_err());

        let mut cfg = default_config();
        cfg.beta = 0.0;
        assert!(PimcSimulation::new(lat, cfg.clone()).is_err());

        let mut cfg = default_config();
        cfg.beta = -1.0;
        assert!(PimcSimulation::new(lat, cfg.clone()).is_err());

        let mut cfg = default_config();
        cfg.angular_step = 0.0;
        assert!(PimcSimulation::new(lat, cfg.clone()).is_err());

        let mut cfg = default_config();
        cfg.n_steps = 0;
        assert!(PimcSimulation::new(lat, cfg).is_err());
    }

    // -- Test 3: thermalisation runs to completion ---------------------

    #[test]
    fn test_thermalisation_runs() {
        let lat = PimcLattice::Chain1D {
            n_spins: 4,
            j: 1.0,
            h_z: 0.0,
        };
        let cfg = PimcConfig {
            n_slices: 3,
            n_steps: 50,
            n_thermalize: 30,
            beta: 1.0,
            seed: 7,
            angular_step: 0.5,
        };
        let mut sim = PimcSimulation::new(lat, cfg).expect("ok");
        sim.thermalize().expect("thermalize ok");

        // Worldlines still on the unit sphere after thermalisation.
        for line in &sim.worldlines {
            for s in line {
                assert!(close(s.magnitude(), 1.0, 1.0e-10));
            }
        }
    }

    // -- Test 4: FM chain at very high T → ⟨M_z⟩ ≈ 0 ------------------

    #[test]
    fn test_high_temperature_zero_magnetization() {
        // β = 0.01 (i.e. T very large) with small field → magnetisation
        // should be very close to 0 with statistical fluctuations.
        let lat = PimcLattice::Chain1D {
            n_spins: 8,
            j: 1.0,
            h_z: 0.05,
        };
        let cfg = PimcConfig {
            n_slices: 4,
            n_steps: 600,
            n_thermalize: 200,
            beta: 0.01,
            seed: 99,
            angular_step: 1.5, // large step at high T
        };
        let mut sim = PimcSimulation::new(lat, cfg).expect("ok");
        sim.thermalize().expect("ok");
        let res = sim.run().expect("ok");
        assert!(
            res.average_magnetization_z.abs() < 0.5,
            "⟨M_z⟩ = {} at high T should be ≈ 0 (within statistics)",
            res.average_magnetization_z
        );
    }

    // -- Test 5: FM chain at very low T with field → ⟨M_z⟩ > 0 --------

    #[test]
    fn test_low_temperature_polarization() {
        // β large, J > 0 ferromagnetic, strong field along z. Expect
        // alignment along z.
        let lat = PimcLattice::Chain1D {
            n_spins: 6,
            j: 1.0,
            h_z: 2.0,
        };
        let cfg = PimcConfig {
            n_slices: 4,
            n_steps: 400,
            n_thermalize: 200,
            beta: 5.0,
            seed: 11,
            angular_step: 0.4,
        };
        let mut sim = PimcSimulation::new(lat, cfg).expect("ok");
        sim.thermalize().expect("ok");
        let res = sim.run().expect("ok");
        assert!(
            res.average_magnetization_z > 0.3,
            "⟨M_z⟩ = {} should be positive and substantial at low T with h>0",
            res.average_magnetization_z
        );
    }

    // -- Test 6: AFM chain (J < 0) → no net z-magnetization at h=0 -----

    #[test]
    fn test_afm_no_net_magnetization() {
        // Antiferromagnetic J < 0 on a ring, no field: total magnetization
        // should be near zero by symmetry.
        let lat = PimcLattice::Ring1D {
            n_spins: 8,
            j: -1.0,
            h_z: 0.0,
        };
        let cfg = PimcConfig {
            n_slices: 4,
            n_steps: 400,
            n_thermalize: 200,
            beta: 2.0,
            seed: 23,
            angular_step: 0.5,
        };
        let mut sim = PimcSimulation::new(lat, cfg).expect("ok");
        sim.thermalize().expect("ok");
        let res = sim.run().expect("ok");
        assert!(
            res.average_magnetization_z.abs() < 0.5,
            "AFM ⟨M_z⟩ = {} should be small at h=0",
            res.average_magnetization_z
        );
    }

    // -- Test 7: susceptibility positive -----------------------------

    #[test]
    fn test_susceptibility_positive() {
        let lat = PimcLattice::Chain1D {
            n_spins: 4,
            j: 1.0,
            h_z: 0.1,
        };
        let cfg = PimcConfig {
            n_slices: 3,
            n_steps: 300,
            n_thermalize: 100,
            beta: 1.0,
            seed: 5,
            angular_step: 0.5,
        };
        let mut sim = PimcSimulation::new(lat, cfg).expect("ok");
        sim.thermalize().expect("ok");
        let res = sim.run().expect("ok");
        assert!(
            res.susceptibility >= 0.0,
            "χ must be ≥ 0, got {}",
            res.susceptibility
        );
    }

    // -- Test 8: reproducibility — same seed → identical results -------

    #[test]
    fn test_reproducibility_same_seed() {
        let lat = PimcLattice::Chain1D {
            n_spins: 4,
            j: 1.0,
            h_z: 0.1,
        };
        let cfg = PimcConfig {
            n_slices: 3,
            n_steps: 50,
            n_thermalize: 20,
            beta: 1.0,
            seed: 12345,
            angular_step: 0.5,
        };
        let mut a = PimcSimulation::new(lat, cfg.clone()).expect("ok");
        let mut b = PimcSimulation::new(lat, cfg).expect("ok");
        a.thermalize().expect("ok");
        b.thermalize().expect("ok");
        let ra = a.run().expect("ok");
        let rb = b.run().expect("ok");
        assert!(close(
            ra.average_magnetization_z,
            rb.average_magnetization_z,
            1.0e-12
        ));
        assert!(close(
            ra.average_energy_per_site,
            rb.average_energy_per_site,
            1.0e-12
        ));
        assert_eq!(ra.n_accepted, rb.n_accepted);
        assert_eq!(ra.n_proposed, rb.n_proposed);
    }

    // -- Test 9: acceptance ratio in a reasonable window ---------------

    #[test]
    fn test_acceptance_ratio_in_window() {
        let lat = PimcLattice::Chain1D {
            n_spins: 4,
            j: 1.0,
            h_z: 0.0,
        };
        let cfg = PimcConfig {
            n_slices: 3,
            n_steps: 400,
            n_thermalize: 100,
            beta: 1.0,
            seed: 33,
            angular_step: 0.3, // tuned for a healthy acceptance
        };
        let mut sim = PimcSimulation::new(lat, cfg).expect("ok");
        sim.thermalize().expect("ok");
        let res = sim.run().expect("ok");
        assert!(res.n_proposed > 0);
        let ratio = res.n_accepted as f64 / res.n_proposed as f64;
        assert!(
            ratio > 0.05 && ratio < 0.99,
            "acceptance ratio {} outside healthy range",
            ratio
        );
    }

    // -- Test 10: n_accepted ≤ n_proposed -----------------------------

    #[test]
    fn test_accepted_le_proposed() {
        let lat = PimcLattice::Ring1D {
            n_spins: 5,
            j: 1.0,
            h_z: 0.05,
        };
        let cfg = PimcConfig {
            n_slices: 4,
            n_steps: 100,
            n_thermalize: 50,
            beta: 0.5,
            seed: 4,
            angular_step: 0.5,
        };
        let mut sim = PimcSimulation::new(lat, cfg).expect("ok");
        sim.thermalize().expect("ok");
        let res = sim.run().expect("ok");
        assert!(res.n_accepted <= res.n_proposed);
        assert_eq!(res.n_proposed, 100 * 5 * 4);
    }

    // -- Test 11: specific heat ≥ 0 -----------------------------------

    #[test]
    fn test_specific_heat_nonnegative() {
        let lat = PimcLattice::Chain1D {
            n_spins: 4,
            j: 1.0,
            h_z: 0.0,
        };
        let cfg = PimcConfig {
            n_slices: 4,
            n_steps: 400,
            n_thermalize: 150,
            beta: 1.5,
            seed: 19,
            angular_step: 0.4,
        };
        let mut sim = PimcSimulation::new(lat, cfg).expect("ok");
        sim.thermalize().expect("ok");
        let res = sim.run().expect("ok");
        assert!(
            res.specific_heat >= 0.0,
            "C = {} must be ≥ 0",
            res.specific_heat
        );
    }

    // -- Bonus: random_unit_vector samples to unit length -------------

    #[test]
    fn test_random_unit_vector_unit_length() {
        let mut rng = seeded_rng(42);
        let normal = Normal::new(0.0, 1.0).expect("ok");
        for _ in 0..50 {
            let v = random_unit_vector(&mut rng, &normal);
            assert!(close(v.magnitude(), 1.0, 1.0e-12));
        }
    }
}
