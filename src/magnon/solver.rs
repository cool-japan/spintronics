//! Time evolution solver for magnon dynamics
//!
//! Implements various numerical methods for solving the LLG equation
//! in a spin chain, including Euler, Heun, and Runge-Kutta methods.

use super::chain::SpinChain;
use crate::vector3::Vector3;

/// RF (radio-frequency) excitation field for spin wave generation
#[derive(Debug, Clone)]
pub struct RfExcitation {
    /// Frequency \[Hz\]
    pub frequency: f64,

    /// Amplitude \[A/m\]
    pub amplitude: f64,

    /// Direction of the oscillating field
    pub direction: Vector3<f64>,

    /// Range of cells to apply excitation (start, end)
    pub excitation_range: (usize, usize),
}

impl RfExcitation {
    /// Create a new RF excitation
    ///
    /// # Arguments
    /// * `frequency` - RF frequency in Hz
    /// * `amplitude` - Field amplitude in A/m
    /// * `direction` - Direction of oscillating field (will be normalized)
    /// * `excitation_range` - (start, end) indices of cells to excite
    pub fn new(
        frequency: f64,
        amplitude: f64,
        direction: Vector3<f64>,
        excitation_range: (usize, usize),
    ) -> Self {
        Self {
            frequency,
            amplitude,
            direction: direction.normalize(),
            excitation_range,
        }
    }

    /// Calculate RF field at given time
    pub fn field_at(&self, t: f64) -> Vector3<f64> {
        let omega = 2.0 * std::f64::consts::PI * self.frequency;
        self.direction * (self.amplitude * (omega * t).sin())
    }

    /// Check if a given cell index is in the excitation range
    pub fn applies_to(&self, idx: usize) -> bool {
        idx >= self.excitation_range.0 && idx < self.excitation_range.1
    }
}

/// Magnon propagation solver
pub struct MagnonSolver {
    /// Static bias field \[A/m\]
    pub bias_field: Vector3<f64>,

    /// Optional RF excitation
    pub rf_excitation: Option<RfExcitation>,

    /// Time step \[s\]
    pub dt: f64,

    /// Current simulation time \[s\]
    pub time: f64,

    /// Preallocated workspace buffers to avoid repeated allocations
    workspace_k1: Vec<Vector3<f64>>,
    workspace_k2: Vec<Vector3<f64>>,
    workspace_k3: Vec<Vector3<f64>>,
    workspace_k4: Vec<Vector3<f64>>,
    workspace_spins: Vec<Vector3<f64>>,
}

impl MagnonSolver {
    /// Create a new magnon solver
    ///
    /// # Arguments
    /// * `dt` - Time step in seconds
    /// * `bias_field` - Static external field in A/m
    pub fn new(dt: f64, bias_field: Vector3<f64>) -> Self {
        Self {
            bias_field,
            rf_excitation: None,
            dt,
            time: 0.0,
            workspace_k1: Vec::new(),
            workspace_k2: Vec::new(),
            workspace_k3: Vec::new(),
            workspace_k4: Vec::new(),
            workspace_spins: Vec::new(),
        }
    }

    /// Ensure workspace buffers have correct size
    fn ensure_workspace_size(&mut self, n: usize) {
        let zero = Vector3::new(0.0, 0.0, 0.0);
        if self.workspace_k1.len() != n {
            self.workspace_k1.resize(n, zero);
            self.workspace_k2.resize(n, zero);
            self.workspace_k3.resize(n, zero);
            self.workspace_k4.resize(n, zero);
            self.workspace_spins.resize(n, zero);
        }
    }

    /// Set RF excitation
    pub fn with_rf_excitation(mut self, rf: RfExcitation) -> Self {
        self.rf_excitation = Some(rf);
        self
    }

    /// Get total external field for a given cell at current time
    fn get_external_field(&self, idx: usize) -> Vector3<f64> {
        let mut h_total = self.bias_field;

        // Add RF field if applicable
        if let Some(ref rf) = self.rf_excitation {
            if rf.applies_to(idx) {
                h_total = h_total + rf.field_at(self.time);
            }
        }

        h_total
    }

    /// Evolve the spin chain by one time step using Forward Euler method
    ///
    /// This is a simple first-order method. For better accuracy, use
    /// `step_heun` or `step_rk4`.
    #[allow(clippy::needless_range_loop)]
    pub fn step_euler(&mut self, chain: &mut SpinChain) {
        let n = chain.n_cells;
        self.ensure_workspace_size(n);

        // Calculate torques for all cells (reuse workspace_k1)
        for i in 0..n {
            let h_ext = self.get_external_field(i);
            let h_eff = chain.calc_effective_field(i, h_ext);
            self.workspace_k1[i] = chain.calc_llg_torque(i, h_eff);
        }

        // Update spins
        for i in 0..n {
            chain.spins[i] = (chain.spins[i] + self.workspace_k1[i] * self.dt).normalize();
        }

        self.time += self.dt;
    }

    /// Evolve using Heun's method (2nd order Runge-Kutta)
    ///
    /// More accurate than Euler with only modest computational overhead.
    #[allow(clippy::needless_range_loop)]
    pub fn step_heun(&mut self, chain: &mut SpinChain) {
        let n = chain.n_cells;
        self.ensure_workspace_size(n);

        // Store original spins in workspace buffer
        self.workspace_spins.copy_from_slice(&chain.spins);

        // Step 1: Calculate k1 (torque at current state)
        for i in 0..n {
            let h_ext = self.get_external_field(i);
            let h_eff = chain.calc_effective_field(i, h_ext);
            self.workspace_k1[i] = chain.calc_llg_torque(i, h_eff);
        }

        // Step 2: Predict using k1
        for i in 0..n {
            chain.spins[i] = (self.workspace_spins[i] + self.workspace_k1[i] * self.dt).normalize();
        }

        // Step 3: Calculate k2 (torque at predicted state)
        // Update time for RF field
        self.time += self.dt;

        for i in 0..n {
            let h_ext = self.get_external_field(i);
            let h_eff = chain.calc_effective_field(i, h_ext);
            self.workspace_k2[i] = chain.calc_llg_torque(i, h_eff);
        }

        // Step 4: Average k1 and k2 for final update
        for i in 0..n {
            let avg_torque = (self.workspace_k1[i] + self.workspace_k2[i]) * 0.5;
            chain.spins[i] = (self.workspace_spins[i] + avg_torque * self.dt).normalize();
        }

        // Time is already updated
    }

    /// Evolve using 4th order Runge-Kutta method
    ///
    /// Most accurate method, recommended for production simulations.
    #[allow(dead_code)]
    #[allow(clippy::needless_range_loop)]
    pub fn step_rk4(&mut self, chain: &mut SpinChain) {
        let n = chain.n_cells;
        self.ensure_workspace_size(n);

        // Store original spins and time
        self.workspace_spins.copy_from_slice(&chain.spins);
        let original_time = self.time;

        // k1
        for i in 0..n {
            let h_ext = self.get_external_field(i);
            let h_eff = chain.calc_effective_field(i, h_ext);
            self.workspace_k1[i] = chain.calc_llg_torque(i, h_eff);
        }

        // k2
        self.time = original_time + 0.5 * self.dt;
        for i in 0..n {
            chain.spins[i] =
                (self.workspace_spins[i] + self.workspace_k1[i] * (0.5 * self.dt)).normalize();
        }
        for i in 0..n {
            let h_ext = self.get_external_field(i);
            let h_eff = chain.calc_effective_field(i, h_ext);
            self.workspace_k2[i] = chain.calc_llg_torque(i, h_eff);
        }

        // k3
        for i in 0..n {
            chain.spins[i] =
                (self.workspace_spins[i] + self.workspace_k2[i] * (0.5 * self.dt)).normalize();
        }
        for i in 0..n {
            let h_ext = self.get_external_field(i);
            let h_eff = chain.calc_effective_field(i, h_ext);
            self.workspace_k3[i] = chain.calc_llg_torque(i, h_eff);
        }

        // k4
        self.time = original_time + self.dt;
        for i in 0..n {
            chain.spins[i] = (self.workspace_spins[i] + self.workspace_k3[i] * self.dt).normalize();
        }
        for i in 0..n {
            let h_ext = self.get_external_field(i);
            let h_eff = chain.calc_effective_field(i, h_ext);
            self.workspace_k4[i] = chain.calc_llg_torque(i, h_eff);
        }

        // Final update: m_new = m_old + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        for i in 0..n {
            let increment = (self.workspace_k1[i]
                + self.workspace_k2[i] * 2.0
                + self.workspace_k3[i] * 2.0
                + self.workspace_k4[i])
                * (self.dt / 6.0);
            chain.spins[i] = (self.workspace_spins[i] + increment).normalize();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::magnon::chain::ChainParameters;

    #[test]
    fn test_rf_excitation() {
        let rf = RfExcitation::new(10.0e9, 1000.0, Vector3::new(0.0, 1.0, 0.0), (0, 5));

        assert!(rf.applies_to(0));
        assert!(rf.applies_to(4));
        assert!(!rf.applies_to(5));
        assert!(!rf.applies_to(10));

        let field = rf.field_at(0.0);
        assert!(field.magnitude() < 1e-10); // sin(0) = 0
    }

    #[test]
    fn test_solver_time_advance() {
        let mut chain = SpinChain::new(10, ChainParameters::default());
        let mut solver = MagnonSolver::new(1.0e-13, Vector3::new(1000.0, 0.0, 0.0));

        assert_eq!(solver.time, 0.0);
        solver.step_euler(&mut chain);
        assert!((solver.time - 1.0e-13).abs() < 1e-20);
    }

    #[test]
    fn test_magnetization_conservation() {
        let mut chain = SpinChain::new(10, ChainParameters::default());
        let mut solver = MagnonSolver::new(1.0e-14, Vector3::new(10000.0, 0.0, 0.0));

        for _ in 0..100 {
            solver.step_euler(&mut chain);
        }

        // All spins should remain normalized
        for spin in &chain.spins {
            assert!((spin.magnitude() - 1.0).abs() < 1e-8);
        }
    }
}
