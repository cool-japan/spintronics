//! Spin chain structure for magnon propagation simulations
//!
//! Represents a one-dimensional chain of magnetic moments (spins) that can
//! support magnon (spin wave) propagation through exchange coupling.

use crate::constants::GAMMA;
use crate::vector3::Vector3;

/// Physical parameters for the spin chain
#[derive(Debug, Clone)]
pub struct ChainParameters {
    /// Exchange stiffness constant \[J/m\]
    /// YIG: ~3.7e-12, Permalloy: ~1.3e-11
    pub a_ex: f64,

    /// Saturation magnetization \[A/m\]
    pub ms: f64,

    /// Gilbert damping constant
    pub alpha: f64,

    /// Cell size (discretization length) \[m\]
    pub cell_size: f64,

    /// Vacuum permeability [H/m]
    pub mu0: f64,
}

impl Default for ChainParameters {
    fn default() -> Self {
        Self {
            a_ex: 3.7e-12,     // YIG
            ms: 1.4e5,         // YIG
            alpha: 0.0001,     // YIG (ultra-low damping)
            cell_size: 5.0e-9, // 5 nm
            mu0: 1.256637e-6,
        }
    }
}

impl ChainParameters {
    /// Create parameters for YIG
    pub fn yig() -> Self {
        Self::default()
    }

    /// Create parameters for Permalloy
    pub fn permalloy() -> Self {
        Self {
            a_ex: 1.3e-11,
            ms: 8.0e5,
            alpha: 0.01,
            ..Self::default()
        }
    }

    /// Calculate maximum recommended time step
    ///
    /// For LLG simulations with exchange coupling, empirically:
    /// dt ~ 0.1 ps is typically stable for nm-scale discretization
    ///
    /// The actual stability limit depends on the numerical method:
    /// - Euler: most restrictive
    /// - Heun: ~2x more stable
    /// - RK4: ~4x more stable
    ///
    /// This returns a conservative estimate suitable for Euler/Heun methods.
    pub fn max_stable_dt(&self) -> f64 {
        // Base time step: 0.1 ps
        let base_dt = 1.0e-13;

        // Scale inversely with cell size squared (finer discretization needs smaller dt)
        let reference_size = 5.0e-9; // 5 nm reference
        let scaling = (self.cell_size / reference_size).powi(2);

        base_dt * scaling
    }
}

/// One-dimensional chain of spins
///
/// Each spin is represented as a normalized magnetization vector.
/// The chain supports exchange coupling between nearest neighbors.
#[derive(Debug, Clone)]
pub struct SpinChain {
    /// Normalized magnetization vectors for each cell
    pub spins: Vec<Vector3<f64>>,

    /// Number of cells in the chain
    pub n_cells: usize,

    /// Physical parameters
    pub params: ChainParameters,
}

impl SpinChain {
    /// Create a new spin chain with all spins aligned along x-axis
    ///
    /// # Arguments
    /// * `n_cells` - Number of discrete cells in the chain
    /// * `params` - Physical parameters for the material
    pub fn new(n_cells: usize, params: ChainParameters) -> Self {
        let spins = vec![Vector3::new(1.0, 0.0, 0.0); n_cells];
        Self {
            spins,
            n_cells,
            params,
        }
    }

    /// Create a chain with initial noise to break symmetry
    ///
    /// This is useful for studying realistic dynamics where perfect
    /// alignment never occurs in practice.
    pub fn new_with_noise(n_cells: usize, params: ChainParameters, noise_amplitude: f64) -> Self {
        let mut spins = Vec::with_capacity(n_cells);
        for i in 0..n_cells {
            let noise = (i as f64 * 0.1).sin() * noise_amplitude;
            spins.push(Vector3::new(1.0, noise, 0.0).normalize());
        }
        Self {
            spins,
            n_cells,
            params,
        }
    }

    /// Calculate the effective magnetic field at a given cell
    ///
    /// H_eff = H_ext + H_exchange + H_anisotropy
    ///
    /// # Arguments
    /// * `idx` - Cell index
    /// * `h_ext` - External magnetic field \[A/m\]
    ///
    /// # Returns
    /// Effective field vector \[A/m\]
    pub fn calc_effective_field(&self, idx: usize, h_ext: Vector3<f64>) -> Vector3<f64> {
        // Start with external field
        let mut h_eff = h_ext;

        // Add exchange field
        h_eff = h_eff + self.calc_exchange_field(idx);

        // Could add anisotropy field here if needed
        // h_eff = h_eff + self.calc_anisotropy_field(idx);

        h_eff
    }

    /// Calculate exchange field using finite differences
    ///
    /// H_ex = (2 A_ex / (μ0 Ms)) ∇²m
    ///
    /// Discretized: ∇²m ≈ (m[i+1] - 2m[i] + m[i-1]) / dx²
    fn calc_exchange_field(&self, idx: usize) -> Vector3<f64> {
        let m_i = self.spins[idx];

        // Neumann boundary conditions: dm/dx = 0 at boundaries
        let m_prev = if idx == 0 { m_i } else { self.spins[idx - 1] };

        let m_next = if idx == self.n_cells - 1 {
            m_i
        } else {
            self.spins[idx + 1]
        };

        // Laplacian
        let laplacian = (m_next + m_prev - m_i * 2.0) * (1.0 / self.params.cell_size.powi(2));

        // Exchange field coefficient
        let coeff_ex = 2.0 * self.params.a_ex / (self.params.mu0 * self.params.ms);

        laplacian * coeff_ex
    }

    /// Calculate LLG torque: dm/dt = -γ/(1+α²) [m × H_eff + α m × (m × H_eff)]
    ///
    /// # Arguments
    /// * `idx` - Cell index
    /// * `h_eff` - Effective field \[A/m\]
    ///
    /// # Returns
    /// Time derivative of magnetization [1/s]
    pub fn calc_llg_torque(&self, idx: usize, h_eff: Vector3<f64>) -> Vector3<f64> {
        let m = self.spins[idx];
        let alpha = self.params.alpha;

        // m × H_eff
        let m_cross_h = m.cross(&h_eff);

        // Precession term
        let term1 = m_cross_h;

        // Damping term: m × (m × H_eff)
        let term2 = m.cross(&m_cross_h);

        // Prefactor: -γ/(1+α²)
        let prefactor = -GAMMA / (1.0 + alpha * alpha);

        (term1 + term2 * alpha) * prefactor
    }

    /// Get total length of the chain \[m\]
    pub fn length(&self) -> f64 {
        self.n_cells as f64 * self.params.cell_size
    }

    /// Get average magnetization across the chain
    pub fn average_magnetization(&self) -> Vector3<f64> {
        let sum = self
            .spins
            .iter()
            .fold(Vector3::new(0.0, 0.0, 0.0), |acc, &m| acc + m);
        sum * (1.0 / self.n_cells as f64)
    }

    /// Reset all spins to align along +z direction
    ///
    /// Useful for initializing reservoir states
    pub fn reset_to_z(&mut self) {
        for i in 0..self.n_cells {
            self.spins[i] = Vector3::new(0.0, 0.0, 1.0);
        }
    }

    /// Evolve the entire chain by one time step using Heun's method
    ///
    /// Simple interface for reservoir computing without needing a separate Solver
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[A/m\]
    /// * `dt` - Time step \[s\]
    #[allow(clippy::needless_range_loop)]
    pub fn evolve_heun(&mut self, h_ext: Vector3<f64>, dt: f64) {
        // First stage: Euler predictor step
        let mut k1 = Vec::with_capacity(self.n_cells);
        for i in 0..self.n_cells {
            let h_eff = self.calc_effective_field(i, h_ext);
            let torque = self.calc_llg_torque(i, h_eff);
            k1.push(torque);
        }

        // Temporarily advance spins for second stage
        let original_spins = self.spins.clone();
        for i in 0..self.n_cells {
            self.spins[i] = (self.spins[i] + k1[i] * dt).normalize();
        }

        // Second stage: evaluate at predicted point
        let mut k2 = Vec::with_capacity(self.n_cells);
        for i in 0..self.n_cells {
            let h_eff = self.calc_effective_field(i, h_ext);
            let torque = self.calc_llg_torque(i, h_eff);
            k2.push(torque);
        }

        // Final update: average of k1 and k2
        for i in 0..self.n_cells {
            let dm_dt = (k1[i] + k2[i]) * 0.5;
            self.spins[i] = (original_spins[i] + dm_dt * dt).normalize();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chain_creation() {
        let chain = SpinChain::new(100, ChainParameters::default());
        assert_eq!(chain.n_cells, 100);
        assert_eq!(chain.spins.len(), 100);
    }

    #[test]
    fn test_spin_normalization() {
        let chain = SpinChain::new_with_noise(50, ChainParameters::default(), 0.1);
        for spin in &chain.spins {
            let mag = spin.magnitude();
            assert!((mag - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_exchange_field_uniform() {
        let chain = SpinChain::new(10, ChainParameters::default());
        // For uniform magnetization, exchange field should be zero
        let h_ex = chain.calc_exchange_field(5);
        assert!(h_ex.magnitude() < 1e-10);
    }

    #[test]
    fn test_max_stable_dt() {
        let params = ChainParameters::yig();
        let dt_max = params.max_stable_dt();
        assert!(dt_max > 0.0);
        // For 5nm cells, dt should be around 0.1 ps
        assert!((dt_max - 1.0e-13).abs() < 1.0e-14);

        // Larger cells should allow larger dt
        let params_large = ChainParameters {
            cell_size: 10.0e-9,
            ..params
        };
        assert!(params_large.max_stable_dt() > dt_max);
    }

    #[test]
    fn test_llg_perpendicular_to_m() {
        let chain = SpinChain::new(10, ChainParameters::default());
        let h_ext = Vector3::new(0.0, 0.0, 1000.0);
        let h_eff = chain.calc_effective_field(5, h_ext);
        let torque = chain.calc_llg_torque(5, h_eff);

        // Torque should be perpendicular to m
        assert!(torque.dot(&chain.spins[5]).abs() < 1e-6);
    }
}
