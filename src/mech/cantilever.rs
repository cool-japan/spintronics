//! Nanomechanical cantilever resonator
//!
//! Models a vibrating cantilever used for detecting spin dynamics

/// Cantilever vibration mode
#[derive(Debug, Clone, Copy)]
pub enum CantileverMode {
    /// Fundamental mode
    Fundamental,
    /// First overtone
    FirstOvertone,
    /// Torsional mode
    Torsional,
}

/// Nanomechanical cantilever
#[derive(Debug, Clone)]
pub struct Cantilever {
    /// Resonance frequency \[Hz\]
    pub f0: f64,

    /// Quality factor (dimensionless)
    pub q_factor: f64,

    /// Effective mass \[kg\]
    pub mass: f64,

    /// Spring constant [N/m]
    pub k_spring: f64,

    /// Current displacement \[m\]
    pub displacement: f64,

    /// Current velocity \[m/s\]
    pub velocity: f64,
}

impl Cantilever {
    /// Create a new cantilever
    pub fn new(f0: f64, q_factor: f64, mass: f64) -> Self {
        let omega = 2.0 * std::f64::consts::PI * f0;
        let k_spring = mass * omega * omega;

        Self {
            f0,
            q_factor,
            mass,
            k_spring,
            displacement: 0.0,
            velocity: 0.0,
        }
    }

    /// Typical AFM cantilever parameters
    pub fn afm_cantilever() -> Self {
        Self::new(
            300.0e3, // 300 kHz
            10000.0, // High Q
            1.0e-12, // 1 pg
        )
    }

    /// Evolve cantilever motion with external force
    ///
    /// Equation of motion:
    /// m d²x/dt² + (m ω₀ / Q) dx/dt + k x = F(t)
    pub fn evolve(&mut self, force: f64, dt: f64) {
        let omega0 = 2.0 * std::f64::consts::PI * self.f0;
        let damping = self.mass * omega0 / self.q_factor;

        // Current acceleration
        let spring_force = -self.k_spring * self.displacement;
        let damping_force = -damping * self.velocity;
        let total_force = spring_force + damping_force + force;
        let a1 = total_force / self.mass;

        // Update position with current velocity
        self.displacement += self.velocity * dt + 0.5 * a1 * dt * dt;

        // Calculate new acceleration at updated position
        let spring_force_new = -self.k_spring * self.displacement;
        let damping_force_new = -damping * self.velocity; // Use old velocity for stability
        let total_force_new = spring_force_new + damping_force_new + force;
        let a2 = total_force_new / self.mass;

        // Update velocity using average acceleration (velocity Verlet)
        self.velocity += 0.5 * (a1 + a2) * dt;
    }

    /// Get current kinetic energy
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.velocity * self.velocity
    }

    /// Get current potential energy
    pub fn potential_energy(&self) -> f64 {
        0.5 * self.k_spring * self.displacement * self.displacement
    }

    /// Total energy
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cantilever_creation() {
        let cant = Cantilever::afm_cantilever();
        assert!(cant.f0 > 0.0);
        assert!(cant.q_factor > 0.0);
    }

    #[test]
    fn test_energy_conservation() {
        let mut cant = Cantilever::afm_cantilever();
        cant.displacement = 1.0e-9; // 1 nm initial displacement

        let e0 = cant.total_energy();

        // Evolve without external force (should conserve energy in Q→∞ limit)
        for _ in 0..100 {
            cant.evolve(0.0, 1.0e-8);
        }

        // Energy should decrease only due to damping
        assert!(cant.total_energy() < e0);
        assert!(cant.total_energy() > 0.0);
    }
}
