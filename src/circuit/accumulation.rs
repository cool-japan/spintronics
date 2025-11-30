//! Spin accumulation and relaxation
//!
//! Models the build-up and decay of spin polarization in non-magnetic conductors.

/// Spin relaxation mechanisms
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RelaxationMechanism {
    /// Elliott-Yafet (momentum scattering)
    ElliottYafet,
    /// D'yakonov-Perel (spin precession)
    DyakonovPerel,
    /// Bir-Aronov-Pikus (electron-hole exchange)
    BirAronovPikus,
}

/// Spin relaxation time parameters
#[derive(Debug, Clone)]
pub struct RelaxationTime {
    /// Spin relaxation time \[s\]
    pub tau_s: f64,

    /// Dominant relaxation mechanism
    pub mechanism: RelaxationMechanism,
}

impl Default for RelaxationTime {
    fn default() -> Self {
        Self {
            tau_s: 100.0e-12, // 100 ps (typical for Cu)
            mechanism: RelaxationMechanism::ElliottYafet,
        }
    }
}

impl RelaxationTime {
    /// Create relaxation time for copper
    pub fn copper() -> Self {
        Self {
            tau_s: 100.0e-12,
            mechanism: RelaxationMechanism::ElliottYafet,
        }
    }

    /// Create relaxation time for aluminum
    pub fn aluminum() -> Self {
        Self {
            tau_s: 200.0e-12,
            mechanism: RelaxationMechanism::ElliottYafet,
        }
    }

    /// Create relaxation time for silicon
    pub fn silicon() -> Self {
        Self {
            tau_s: 1.0e-9, // 1 ns (longer in semiconductors)
            mechanism: RelaxationMechanism::DyakonovPerel,
        }
    }

    /// Calculate spin diffusion length
    ///
    /// λ_s = √(D × τ_s)
    ///
    /// where D is the diffusion constant
    pub fn spin_diffusion_length(&self, diffusion_const: f64) -> f64 {
        (diffusion_const * self.tau_s).sqrt()
    }
}

/// Spin accumulation in a conductor
#[derive(Debug, Clone)]
pub struct SpinAccumulation {
    /// Spin chemical potential \[J\]
    pub mu_s: f64,

    /// Charge chemical potential \[J\]
    pub mu_c: f64,

    /// Position \[m\]
    pub position: f64,

    /// Relaxation time
    pub relaxation: RelaxationTime,
}

impl Default for SpinAccumulation {
    fn default() -> Self {
        Self {
            mu_s: 0.0,
            mu_c: 0.0,
            position: 0.0,
            relaxation: RelaxationTime::default(),
        }
    }
}

impl SpinAccumulation {
    /// Create a new spin accumulation
    pub fn new(position: f64, relaxation: RelaxationTime) -> Self {
        Self {
            mu_s: 0.0,
            mu_c: 0.0,
            position,
            relaxation,
        }
    }

    /// Set spin accumulation from injected spin current
    ///
    /// μ_s = (2 R_s / A) × I_s
    ///
    /// where R_s is spin resistance and A is area
    pub fn from_spin_current(&mut self, spin_current: f64, spin_resistance: f64, area: f64) {
        self.mu_s = (2.0 * spin_resistance / area) * spin_current;
    }

    /// Calculate spin polarization
    ///
    /// P_s = μ_s / (2 k_B T)
    ///
    /// For low energies, can use: P_s ≈ μ_s / E_F
    pub fn polarization(&self, fermi_energy: f64) -> f64 {
        if fermi_energy > 0.0 {
            self.mu_s / fermi_energy
        } else {
            0.0
        }
    }

    /// Time evolution of spin accumulation
    ///
    /// dμ_s/dt = -μ_s / τ_s + source
    ///
    /// # Arguments
    /// * `source` - Source term [J/s]
    /// * `dt` - Time step \[s\]
    pub fn evolve(&mut self, source: f64, dt: f64) {
        let decay = -self.mu_s / self.relaxation.tau_s;
        let dmu_dt = decay + source;
        self.mu_s += dmu_dt * dt;
    }

    /// Steady-state spin accumulation with continuous injection
    ///
    /// μ_s = source × τ_s
    pub fn steady_state(&mut self, source: f64) {
        self.mu_s = source * self.relaxation.tau_s;
    }

    /// Calculate spin current from spatial gradient
    ///
    /// I_s = -σ_s ∇μ_s
    ///
    /// where σ_s is spin conductivity
    pub fn spin_current_from_gradient(&self, gradient: f64, spin_conductivity: f64) -> f64 {
        -spin_conductivity * gradient
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_relaxation_time() {
        let cu = RelaxationTime::copper();
        let al = RelaxationTime::aluminum();
        let si = RelaxationTime::silicon();

        assert!(cu.tau_s < al.tau_s);
        assert!(al.tau_s < si.tau_s);
    }

    #[test]
    fn test_spin_diffusion_length() {
        let relaxation = RelaxationTime::copper();
        let d = 1.0e-3; // Diffusion constant for Cu

        let lambda_s = relaxation.spin_diffusion_length(d);
        assert!(lambda_s > 0.0);
        assert!(lambda_s < 1.0e-3); // Should be μm scale
    }

    #[test]
    fn test_spin_accumulation_evolution() {
        let mut acc = SpinAccumulation::new(0.0, RelaxationTime::copper());

        // Apply constant source
        let source = 1.0e-15; // J/s
        let dt = 1.0e-12; // 1 ps

        for _ in 0..1000 {
            acc.evolve(source, dt);
        }

        // Should approach steady state
        let expected_steady = source * acc.relaxation.tau_s;
        assert!((acc.mu_s - expected_steady).abs() / expected_steady < 0.1);
    }

    #[test]
    fn test_steady_state() {
        let mut acc = SpinAccumulation::new(0.0, RelaxationTime::aluminum());
        let source = 1.0e-15;

        acc.steady_state(source);

        assert!((acc.mu_s - source * acc.relaxation.tau_s).abs() < 1e-20);
    }
}
