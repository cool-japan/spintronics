//! Spin ice physics
//!
//! This module implements the physics of spin ice materials, including:
//! - Ice rules ("2-in, 2-out" constraint on tetrahedra)
//! - Pauling residual entropy
//! - Emergent magnetic monopoles
//! - Constrained Monte Carlo simulation
//!
//! # Physics Background
//!
//! Spin ice materials (e.g., Dy₂Ti₂O₇, Ho₂Ti₂O₇) are frustrated magnets
//! on the pyrochlore lattice where large Ising-like spins point along
//! local \[111\] directions. The ground state manifold is degenerate with
//! the "ice rule" constraint: 2 spins point in and 2 spins point out
//! of each tetrahedron.
//!
//! This degeneracy leads to a residual entropy at T → 0 given by
//! Pauling's formula: S_P = (R/2) ln(3/2) per mole of spins,
//! identical to the entropy of water ice.
//!
//! Excitations above the ice manifold are magnetic monopoles: point-like
//! sources of H-field that interact via a Coulomb potential.
//!
//! # References
//!
//! - S.T. Bramwell & M.J.P. Gingras, "Spin Ice State in Frustrated
//!   Magnetic Pyrochlore Materials", Science 294, 1495-1501 (2001)
//! - C. Castelnovo, R. Moessner, S.L. Sondhi, "Magnetic monopoles in
//!   spin ice", Nature 451, 42-45 (2008)
//! - A.P. Ramirez et al., "Zero-point entropy in 'spin ice'",
//!   Nature 399, 333-335 (1999)

use super::lattice::{FrustratedLattice, Xorshift64};
use crate::constants::{KB, MU_0, NA};
use crate::error::{Error, Result};
use crate::vector3::Vector3;

/// Gas constant R = N_A * k_B [J/(mol·K)]
const R_GAS: f64 = NA * KB;

/// Parameters for a spin ice material
#[derive(Debug, Clone)]
pub struct SpinIceParams {
    /// Nearest-neighbor exchange coupling J_nn \[K\] (in temperature units)
    pub j_nn: f64,
    /// Dipolar coupling constant D_nn \[K\]
    pub d_nn: f64,
    /// Effective nearest-neighbor coupling J_eff = J_nn + D_nn \[K\]
    pub j_eff: f64,
    /// Magnetic moment per ion \[μ_B\]
    pub moment: f64,
    /// Nearest-neighbor distance \[m\]
    pub nn_distance: f64,
    /// Material name
    pub name: String,
}

impl SpinIceParams {
    /// Parameters for Dy₂Ti₂O₇
    ///
    /// Dysprosium titanate is the prototypical spin ice material.
    /// J_nn ≈ -3.72 K, D_nn ≈ 2.35 K, moment ≈ 10 μ_B
    pub fn dy2ti2o7() -> Self {
        Self {
            j_nn: -3.72,
            d_nn: 2.35,
            j_eff: -3.72 + 2.35,
            moment: 10.0,
            nn_distance: 3.54e-10, // a/√8 with a = 10.12 Å
            name: "Dy2Ti2O7".to_string(),
        }
    }

    /// Parameters for Ho₂Ti₂O₇
    ///
    /// Holmium titanate spin ice.
    /// J_nn ≈ -1.56 K, D_nn ≈ 2.35 K, moment ≈ 10 μ_B
    pub fn ho2ti2o7() -> Self {
        Self {
            j_nn: -1.56,
            d_nn: 2.35,
            j_eff: -1.56 + 2.35,
            moment: 10.0,
            nn_distance: 3.54e-10,
            name: "Ho2Ti2O7".to_string(),
        }
    }
}

/// Represents a tetrahedron in the pyrochlore lattice
#[derive(Debug, Clone)]
pub struct Tetrahedron {
    /// Indices of the 4 spins in this tetrahedron
    pub site_indices: [usize; 4],
    /// Local Ising axes for each site (pointing toward/away from center)
    pub ising_axes: [Vector3<f64>; 4],
}

/// Spin ice system on a pyrochlore lattice
#[derive(Debug, Clone)]
pub struct SpinIce {
    /// Underlying pyrochlore lattice
    pub lattice: FrustratedLattice,
    /// Material parameters
    pub params: SpinIceParams,
    /// Tetrahedra (groups of 4 corner-sharing sites)
    pub tetrahedra: Vec<Tetrahedron>,
    /// Ising spin values (+1 or -1 along local axis)
    pub ising_spins: Vec<f64>,
}

impl SpinIce {
    /// Create a new spin ice system
    ///
    /// # Arguments
    ///
    /// * `nx` - System size along x
    /// * `ny` - System size along y
    /// * `params` - Spin ice material parameters
    ///
    /// # Errors
    ///
    /// Returns error if lattice creation fails.
    pub fn new(nx: usize, ny: usize, params: SpinIceParams) -> Result<Self> {
        let coupling_j = params.j_eff.abs() * KB; // Convert from K to J
        let lattice =
            FrustratedLattice::pyrochlore(nx, ny, coupling_j, params.nn_distance * 8.0_f64.sqrt())?;

        let n_sites = lattice.num_sites();

        // Define Ising axes for each sublattice
        let ising_axes = [
            Vector3::new(1.0, 1.0, 1.0).normalize(),
            Vector3::new(1.0, -1.0, -1.0).normalize(),
            Vector3::new(-1.0, 1.0, -1.0).normalize(),
            Vector3::new(-1.0, -1.0, 1.0).normalize(),
        ];

        // Build tetrahedra from the lattice structure
        // Each unit cell in pyrochlore contains one "up" tetrahedron
        let nz = nx.min(ny);
        let n_cells = nx * ny * nz;
        let mut tetrahedra = Vec::with_capacity(n_cells);

        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let base = 4 * ((iz * ny + iy) * nx + ix);
                    let tet = Tetrahedron {
                        site_indices: [base, base + 1, base + 2, base + 3],
                        ising_axes,
                    };
                    tetrahedra.push(tet);
                }
            }
        }

        // Initialize all Ising spins to +1 (not satisfying ice rules yet)
        let ising_spins = vec![1.0; n_sites];

        let mut ice = Self {
            lattice,
            params,
            tetrahedra,
            ising_spins,
        };

        // Set the actual spin vectors from Ising values
        ice.update_spin_vectors();

        Ok(ice)
    }

    /// Update the 3D spin vectors from the Ising spin values
    fn update_spin_vectors(&mut self) {
        let ising_axes = [
            Vector3::new(1.0, 1.0, 1.0).normalize(),
            Vector3::new(1.0, -1.0, -1.0).normalize(),
            Vector3::new(-1.0, 1.0, -1.0).normalize(),
            Vector3::new(-1.0, -1.0, 1.0).normalize(),
        ];

        for (i, spin) in self.lattice.spins.iter_mut().enumerate() {
            let sub = i % 4;
            *spin = ising_axes[sub] * self.ising_spins[i];
        }
    }

    /// Check if a tetrahedron satisfies the ice rule (2-in, 2-out)
    ///
    /// The "in" or "out" direction is defined relative to the tetrahedron center.
    /// For the "up" tetrahedron, spins along +\[111\] directions point "out".
    pub fn check_ice_rule(&self, tet_index: usize) -> bool {
        if tet_index >= self.tetrahedra.len() {
            return false;
        }

        let tet = &self.tetrahedra[tet_index];
        let mut sum = 0.0;
        for &idx in &tet.site_indices {
            if idx < self.ising_spins.len() {
                sum += self.ising_spins[idx];
            }
        }
        // Ice rule: 2-in, 2-out means the sum of Ising spins = 0
        // (2 × (+1) + 2 × (-1) = 0)
        sum.abs() < 0.5
    }

    /// Count the number of tetrahedra violating the ice rule
    pub fn count_ice_rule_violations(&self) -> usize {
        (0..self.tetrahedra.len())
            .filter(|&i| !self.check_ice_rule(i))
            .count()
    }

    /// Calculate the fraction of tetrahedra satisfying the ice rule
    pub fn ice_rule_fraction(&self) -> f64 {
        if self.tetrahedra.is_empty() {
            return 0.0;
        }
        let satisfied = self.tetrahedra.len() - self.count_ice_rule_violations();
        satisfied as f64 / self.tetrahedra.len() as f64
    }

    /// Initialize the system to a state satisfying ice rules
    ///
    /// Uses a greedy algorithm: go tetrahedron by tetrahedron and flip spins
    /// to satisfy 2-in, 2-out on each.
    ///
    /// # Arguments
    ///
    /// * `seed` - PRNG seed
    ///
    /// # Errors
    ///
    /// Returns error if seed is zero.
    pub fn initialize_ice_state(&mut self, seed: u64) -> Result<()> {
        let mut rng = Xorshift64::new(seed)?;

        // First randomize all spins to ±1
        for spin in self.ising_spins.iter_mut() {
            *spin = if rng.next_f64() < 0.5 { 1.0 } else { -1.0 };
        }

        // Now fix tetrahedra one by one
        for tet_idx in 0..self.tetrahedra.len() {
            let tet = &self.tetrahedra[tet_idx];
            let indices = tet.site_indices;

            let sum: f64 = indices.iter().map(|&i| self.ising_spins[i]).sum();

            // Need sum = 0 for ice rule
            if (sum - 0.0).abs() > 0.5 {
                // Flip spins to get closer to zero
                let n_to_flip = (sum.abs() / 2.0).round() as usize;
                let mut flipped = 0;
                for &idx in &indices {
                    if flipped >= n_to_flip {
                        break;
                    }
                    if (sum > 0.0 && self.ising_spins[idx] > 0.0)
                        || (sum < 0.0 && self.ising_spins[idx] < 0.0)
                    {
                        self.ising_spins[idx] = -self.ising_spins[idx];
                        flipped += 1;
                    }
                }
            }
        }

        self.update_spin_vectors();
        Ok(())
    }

    /// Run Monte Carlo with ice-rule-preserving loop moves
    ///
    /// Uses single spin flips with Metropolis criterion, but only accepts
    /// moves that do not create ice rule violations on any affected tetrahedron.
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature \[K\]
    /// * `n_sweeps` - Number of sweeps
    /// * `seed` - PRNG seed
    ///
    /// # Errors
    ///
    /// Returns error if temperature is negative or seed is zero.
    pub fn constrained_mc(&mut self, temperature: f64, n_sweeps: usize, seed: u64) -> Result<f64> {
        if temperature < 0.0 {
            return Err(Error::InvalidParameter {
                param: "temperature".to_string(),
                reason: "temperature must be non-negative".to_string(),
            });
        }

        let mut rng = Xorshift64::new(seed)?;
        let n_sites = self.lattice.num_sites();
        let beta = if temperature > 1e-30 {
            1.0 / (KB * temperature)
        } else {
            f64::INFINITY
        };

        for _sweep in 0..n_sweeps {
            for _step in 0..n_sites {
                let site = (rng.next_u64() as usize) % n_sites;
                let old_ising = self.ising_spins[site];

                // Flip the Ising spin
                self.ising_spins[site] = -old_ising;

                // Check ice rules on all tetrahedra containing this site
                let mut violates = false;
                for (tet_idx, tet) in self.tetrahedra.iter().enumerate() {
                    if tet.site_indices.contains(&site) && !self.check_ice_rule(tet_idx) {
                        violates = true;
                        break;
                    }
                }

                if violates {
                    // Reject: restore
                    self.ising_spins[site] = old_ising;
                    continue;
                }

                // Calculate energy change (simplified nearest-neighbor model)
                // Flipping spin i changes energy by -2 * J_eff * σ_i * Σ_j σ_j
                let mut neighbor_sum = 0.0;
                if site < self.lattice.neighbors.len() {
                    for &j in &self.lattice.neighbors[site] {
                        if j < self.ising_spins.len() {
                            neighbor_sum += self.ising_spins[j];
                        }
                    }
                }
                let delta_e = -2.0 * self.params.j_eff * KB * (-old_ising) * neighbor_sum;

                // Metropolis criterion
                if delta_e > 0.0 {
                    if beta.is_infinite() {
                        self.ising_spins[site] = old_ising;
                    } else {
                        let acceptance = (-beta * delta_e).exp();
                        if rng.next_f64() >= acceptance {
                            self.ising_spins[site] = old_ising;
                        }
                    }
                }
            }
        }

        self.update_spin_vectors();
        Ok(self.lattice.total_energy() / n_sites as f64)
    }

    /// Calculate total magnetization of the spin ice
    pub fn total_magnetization(&self) -> Vector3<f64> {
        self.lattice.average_magnetization()
    }
}

/// Calculate the Pauling residual entropy per mole of spins
///
/// S_P = (R/2) ln(3/2)
///
/// This is the entropy of the ice-rule-obeying ground state manifold,
/// arising from the macroscopic degeneracy of 2-in/2-out configurations
/// on the pyrochlore lattice.
///
/// # Returns
///
/// Pauling entropy in J/(mol·K)
pub fn pauling_entropy() -> f64 {
    0.5 * R_GAS * (1.5_f64).ln()
}

/// Calculate the Pauling entropy per spin in units of k_B
///
/// s_P = (1/2) ln(3/2) per spin
pub fn pauling_entropy_per_spin_kb() -> f64 {
    0.5 * (1.5_f64).ln()
}

/// Calculate the magnetic monopole creation energy
///
/// When an ice rule is violated (3-in/1-out or 1-in/3-out), it creates
/// a monopole-antimonopole pair. The energy cost is:
///
/// E_monopole = 2 * J_eff * μ² / (r_nn)
///
/// where J_eff is the effective nearest-neighbor coupling.
///
/// In practice, for Dy₂Ti₂O₇, the monopole creation energy is
/// approximately 4.35 K (≈ 6e-23 J).
///
/// # Arguments
///
/// * `j_eff` - Effective coupling \[K\]
/// * `moment_mu_b` - Magnetic moment in Bohr magnetons
/// * `nn_distance` - Nearest-neighbor distance \[m\]
///
/// # Returns
///
/// Monopole creation energy \[J\]
///
/// # Errors
///
/// Returns error if nn_distance is not positive.
pub fn monopole_creation_energy(j_eff: f64, moment_mu_b: f64, nn_distance: f64) -> Result<f64> {
    if nn_distance <= 0.0 {
        return Err(Error::InvalidParameter {
            param: "nn_distance".to_string(),
            reason: "nearest-neighbor distance must be positive".to_string(),
        });
    }

    // Energy to create a monopole pair from the ice manifold
    // Each violated tetrahedron costs additional J_eff per misaligned bond
    // A 3-in/1-out defect has 2 more "wrong" bonds than 2-in/2-out
    let mu = moment_mu_b * crate::constants::MU_B;
    let energy = 2.0 * j_eff.abs() * KB
        + MU_0 * mu * mu / (4.0 * std::f64::consts::PI * nn_distance.powi(3));

    Ok(energy)
}

/// Calculate the Coulomb interaction between magnetic monopoles
///
/// Magnetic monopoles in spin ice interact via a magnetic Coulomb law:
///
/// V(r) = (μ₀/4π) * Q_m² / r
///
/// where Q_m = 2μ/a_d is the magnetic charge (a_d = diamond lattice constant).
///
/// # Arguments
///
/// * `separation` - Distance between monopoles \[m\]
/// * `moment_mu_b` - Magnetic moment per ion \[μ_B\]
/// * `diamond_lattice_const` - Diamond sublattice constant \[m\]
///
/// # Returns
///
/// Interaction energy \[J\]
///
/// # Errors
///
/// Returns error if separation or lattice constant is not positive.
pub fn monopole_coulomb_interaction(
    separation: f64,
    moment_mu_b: f64,
    diamond_lattice_const: f64,
) -> Result<f64> {
    if separation <= 0.0 {
        return Err(Error::InvalidParameter {
            param: "separation".to_string(),
            reason: "monopole separation must be positive".to_string(),
        });
    }
    if diamond_lattice_const <= 0.0 {
        return Err(Error::InvalidParameter {
            param: "diamond_lattice_const".to_string(),
            reason: "lattice constant must be positive".to_string(),
        });
    }

    let mu = moment_mu_b * crate::constants::MU_B;
    let q_m = 2.0 * mu / diamond_lattice_const;

    // V = (μ₀/4π) * Q_m² / r
    let energy = MU_0 * q_m * q_m / (4.0 * std::f64::consts::PI * separation);

    Ok(energy)
}

/// Calculate the density of monopoles at a given temperature
///
/// In the dilute monopole gas approximation:
/// n_m ∝ exp(-E_create / (2 k_B T))
///
/// Returns the Boltzmann factor (dimensionless) for monopole density.
///
/// # Arguments
///
/// * `temperature` - Temperature \[K\]
/// * `creation_energy` - Monopole creation energy \[J\]
///
/// # Errors
///
/// Returns error if temperature is not positive.
pub fn monopole_density_factor(temperature: f64, creation_energy: f64) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(Error::InvalidParameter {
            param: "temperature".to_string(),
            reason: "temperature must be positive for thermal activation".to_string(),
        });
    }

    // Each monopole costs E_create/2 (pair creation splits cost)
    let exponent = -creation_energy / (2.0 * KB * temperature);
    Ok(exponent.exp())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frustrated::lattice::LatticeType;

    #[test]
    fn test_pauling_entropy() {
        let s_p = pauling_entropy();
        // S_P = (R/2) ln(3/2) ≈ 0.5 * 8.314 * 0.4055 ≈ 1.685 J/(mol·K)
        let expected = 0.5 * R_GAS * (1.5_f64).ln();
        assert!(
            (s_p - expected).abs() < 1e-10,
            "Pauling entropy = {}, expected {}",
            s_p,
            expected
        );
        // Numerical check
        assert!(
            (s_p - 1.6854).abs() < 0.01,
            "Pauling entropy = {} J/(mol·K), expected ~1.685",
            s_p
        );
    }

    #[test]
    fn test_pauling_entropy_per_spin() {
        let s_per_spin = pauling_entropy_per_spin_kb();
        let expected = 0.5 * (1.5_f64).ln();
        assert!(
            (s_per_spin - expected).abs() < 1e-15,
            "Pauling entropy per spin = {}, expected {}",
            s_per_spin,
            expected
        );
    }

    #[test]
    fn test_spin_ice_creation() {
        let params = SpinIceParams::dy2ti2o7();
        let ice = SpinIce::new(3, 3, params).expect("failed to create spin ice");
        assert_eq!(ice.lattice.lattice_type, LatticeType::Pyrochlore);
        assert!(ice.lattice.num_sites() > 0);
    }

    #[test]
    fn test_ice_rules_after_initialization() {
        let params = SpinIceParams::dy2ti2o7();
        let mut ice = SpinIce::new(3, 3, params).expect("failed to create spin ice");
        ice.initialize_ice_state(42)
            .expect("failed to init ice state");

        // After initialization, most tetrahedra should satisfy ice rules
        let fraction = ice.ice_rule_fraction();
        assert!(
            fraction > 0.5,
            "ice rule fraction = {}, expected > 0.5 after initialization",
            fraction
        );
    }

    #[test]
    fn test_monopole_energy_positive() {
        let params = SpinIceParams::dy2ti2o7();
        let energy = monopole_creation_energy(params.j_eff, params.moment, params.nn_distance)
            .expect("failed to compute monopole energy");
        assert!(
            energy > 0.0,
            "monopole creation energy = {}, must be positive",
            energy
        );
    }

    #[test]
    fn test_monopole_coulomb_positive() {
        let separation = 5e-10; // 5 Angstrom
        let a_d = 4.3e-10; // diamond lattice constant
        let energy = monopole_coulomb_interaction(separation, 10.0, a_d)
            .expect("failed to compute Coulomb energy");
        assert!(
            energy > 0.0,
            "Coulomb energy = {}, must be positive (repulsive for like charges)",
            energy
        );
    }

    #[test]
    fn test_monopole_coulomb_decreases_with_distance() {
        let a_d = 4.3e-10;
        let e1 = monopole_coulomb_interaction(5e-10, 10.0, a_d).expect("close distance");
        let e2 = monopole_coulomb_interaction(10e-10, 10.0, a_d).expect("far distance");
        assert!(
            e1 > e2,
            "Coulomb energy should decrease with distance: {} > {}",
            e1,
            e2
        );
    }

    #[test]
    fn test_monopole_density_decreases_with_lower_temp() {
        let e_create = 6e-23; // ~4.3 K creation energy
        let n_high = monopole_density_factor(10.0, e_create).expect("high T");
        let n_low = monopole_density_factor(1.0, e_create).expect("low T");
        assert!(
            n_high > n_low,
            "monopole density should increase with temperature"
        );
    }

    #[test]
    fn test_spin_ice_params_dy2ti2o7() {
        let p = SpinIceParams::dy2ti2o7();
        assert!((p.moment - 10.0).abs() < 0.1);
        assert!(p.j_nn < 0.0); // ferromagnetic exchange
        assert!(p.d_nn > 0.0); // dipolar
        assert!((p.j_eff - (p.j_nn + p.d_nn)).abs() < 1e-10);
    }

    #[test]
    fn test_spin_ice_params_ho2ti2o7() {
        let p = SpinIceParams::ho2ti2o7();
        assert!((p.moment - 10.0).abs() < 0.1);
        assert!(p.j_eff > 0.0); // Ho2Ti2O7 has J_eff > 0
    }

    #[test]
    fn test_monopole_energy_invalid_distance() {
        let result = monopole_creation_energy(1.0, 10.0, 0.0);
        assert!(result.is_err());
        let result = monopole_creation_energy(1.0, 10.0, -1.0);
        assert!(result.is_err());
    }
}
