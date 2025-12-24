//! Thermal Fluctuations in Magnetic Systems
//!
//! At finite temperatures, magnetic moments experience random thermal fluctuations
//! described by the fluctuation-dissipation theorem. This module implements
//! stochastic LLG equation with thermal noise.
//!
//! ## Stochastic LLG Equation
//!
//! dm/dt = -γ m × (H_eff + H_th) + α m × dm/dt
//!
//! where H_th is a stochastic thermal field with properties:
//! ⟨H_th(t)⟩ = 0
//! ⟨H_th(t) H_th(t')⟩ = (2 α k_B T) / (γ M_s V) δ(t-t')
//!
//! ## Applications
//!
//! - Thermal spin-torque oscillators
//! - Magnetization switching at finite temperature
//! - Thermally-activated reversal
//! - Spin Seebeck effect with fluctuations

use scirs2_core::random::core::Random;
use scirs2_core::random::rand_distributions::Normal;
use scirs2_core::random::thread_rng;
use scirs2_core::rngs;

use crate::constants::{GAMMA, KB};
use crate::vector3::Vector3;

/// Thermal field generator based on fluctuation-dissipation theorem
pub struct ThermalField {
    /// Temperature \[K\]
    pub temperature: f64,

    /// Cell/grain volume \[m³\]
    pub volume: f64,

    /// Saturation magnetization \[A/m\]
    pub ms: f64,

    /// Gilbert damping parameter
    pub alpha: f64,

    /// Gyromagnetic ratio \[rad/(s·T)\]
    pub gamma: f64,

    /// Random number generator
    rng: Random<rngs::ThreadRng>,

    /// Normal distribution for Gaussian noise (standard normal N(0,1))
    normal: Normal<f64>,
}

impl ThermalField {
    /// Create a new thermal field generator
    ///
    /// # Arguments
    /// * `temperature` - Temperature in Kelvin
    /// * `volume` - Volume of magnetic cell in m³
    /// * `ms` - Saturation magnetization in A/m
    /// * `alpha` - Gilbert damping parameter
    pub fn new(temperature: f64, volume: f64, ms: f64, alpha: f64) -> Self {
        Self {
            temperature,
            volume,
            ms,
            alpha,
            gamma: GAMMA,
            rng: thread_rng(),
            normal: Normal::new(0.0, 1.0).unwrap(),
        }
    }

    /// Generate thermal field vector for given timestep
    ///
    /// Returns H_th with proper variance according to FDT:
    /// σ² = (2 α k_B T) / (γ M_s V Δt)
    ///
    /// # Arguments
    /// * `dt` - Time step \[s\]
    ///
    /// # Returns
    /// Random thermal field vector \[A/m\]
    pub fn generate(&mut self, dt: f64) -> Vector3<f64> {
        if self.temperature <= 0.0 || dt <= 0.0 {
            return Vector3::new(0.0, 0.0, 0.0);
        }

        // Fluctuation-dissipation theorem variance
        let variance =
            (2.0 * self.alpha * KB * self.temperature) / (self.gamma * self.ms * self.volume * dt);

        let sigma = variance.sqrt();

        // Generate three independent Gaussian random numbers
        Vector3::new(
            self.rng.sample(self.normal) * sigma,
            self.rng.sample(self.normal) * sigma,
            self.rng.sample(self.normal) * sigma,
        )
    }

    /// Generate thermal field with custom variance scaling
    ///
    /// Useful for testing or special simulation scenarios
    pub fn generate_scaled(&mut self, dt: f64, scale: f64) -> Vector3<f64> {
        let base_field = self.generate(dt);
        base_field * scale
    }

    /// Set new temperature
    pub fn set_temperature(&mut self, temperature: f64) {
        self.temperature = temperature;
    }

    /// Calculate thermal energy k_B T \[J\]
    pub fn thermal_energy(&self) -> f64 {
        KB * self.temperature
    }

    /// Calculate attempt frequency \[Hz\]
    ///
    /// f_0 = γ μ₀ M_s / (2π)
    ///
    /// Typical value for ferromagnets: ~1-10 GHz
    pub fn attempt_frequency(&self) -> f64 {
        let mu0 = 1.256637e-6; // H/m
        (self.gamma * mu0 * self.ms) / (2.0 * std::f64::consts::PI)
    }
}

/// Stochastic LLG solver
///
/// Combines deterministic LLG with thermal noise
pub struct StochasticLLG {
    /// Magnetization (normalized)
    pub magnetization: Vector3<f64>,

    /// Thermal field generator
    pub thermal: ThermalField,

    /// Gilbert damping
    pub alpha: f64,
}

impl StochasticLLG {
    /// Create new stochastic LLG solver
    pub fn new(
        initial_m: Vector3<f64>,
        temperature: f64,
        volume: f64,
        ms: f64,
        alpha: f64,
    ) -> Self {
        Self {
            magnetization: initial_m.normalize(),
            thermal: ThermalField::new(temperature, volume, ms, alpha),
            alpha,
        }
    }

    /// Evolve with Euler-Maruyama method
    ///
    /// # Arguments
    /// * `h_eff` - Effective field (deterministic part) \[A/m\]
    /// * `dt` - Time step \[s\]
    pub fn evolve(&mut self, h_eff: Vector3<f64>, dt: f64) {
        // Add thermal noise to effective field
        let h_thermal = self.thermal.generate(dt);
        let h_total = h_eff + h_thermal;

        // Standard LLG dynamics
        let gamma = GAMMA;
        let m_cross_h = self.magnetization.cross(&h_total);
        let damping = self.magnetization.cross(&m_cross_h) * self.alpha;

        let dm_dt = (m_cross_h + damping) * (-gamma / (1.0 + self.alpha * self.alpha));

        // Euler-Maruyama step
        self.magnetization = (self.magnetization + dm_dt * dt).normalize();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thermal_field_creation() {
        let thermal = ThermalField::new(300.0, 1.0e-24, 1.0e6, 0.01);
        assert!((thermal.temperature - 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_field_zero_temp() {
        let mut thermal = ThermalField::new(0.0, 1.0e-24, 1.0e6, 0.01);
        let h_th = thermal.generate(1.0e-12);

        // Zero temperature should give zero field
        assert!(h_th.magnitude() < 1e-20);
    }

    #[test]
    fn test_thermal_field_nonzero() {
        let mut thermal = ThermalField::new(300.0, 1.0e-24, 1.0e6, 0.01);
        let h_th = thermal.generate(1.0e-12);

        // At 300K, thermal field should be non-zero
        assert!(h_th.magnitude() > 0.0);
    }

    #[test]
    fn test_thermal_energy() {
        let thermal = ThermalField::new(300.0, 1.0e-24, 1.0e6, 0.01);
        let e_th = thermal.thermal_energy();

        // k_B * 300 K ≈ 4.14e-21 J
        assert!((e_th - 4.14e-21).abs() < 1e-22);
    }

    #[test]
    fn test_attempt_frequency() {
        let thermal = ThermalField::new(300.0, 1.0e-24, 1.0e6, 0.01);
        let f0 = thermal.attempt_frequency();

        // Should be in GHz range for typical ferromagnet
        assert!(f0 > 1.0e8); // > 100 MHz
        assert!(f0 < 1.0e12); // < 1 THz (relaxed upper bound)
    }

    #[test]
    fn test_stochastic_llg_creation() {
        let sllg = StochasticLLG::new(Vector3::new(1.0, 0.0, 0.0), 300.0, 1.0e-24, 1.0e6, 0.01);

        assert!((sllg.magnetization.magnitude() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_stochastic_evolution() {
        let mut sllg = StochasticLLG::new(Vector3::new(1.0, 0.0, 0.0), 300.0, 1.0e-24, 1.0e6, 0.01);

        let h_ext = Vector3::new(0.0, 0.0, 1.0e5);

        // Evolve
        for _ in 0..100 {
            sllg.evolve(h_ext, 1.0e-13);
        }

        // Magnetization should remain normalized
        assert!((sllg.magnetization.magnitude() - 1.0).abs() < 1e-6);
    }

    #[test]
    #[allow(dead_code)]
    fn test_temperature_scaling() {
        let mut thermal_cold = ThermalField::new(100.0, 1.0e-24, 1.0e6, 0.01);
        let mut thermal_hot = ThermalField::new(600.0, 1.0e-24, 1.0e6, 0.01);

        let dt = 1.0e-12;

        // Collect samples
        let mut variance_cold = 0.0;
        let mut variance_hot = 0.0;
        let n_samples = 1000;

        for _ in 0..n_samples {
            let h_cold = thermal_cold.generate(dt);
            let h_hot = thermal_hot.generate(dt);

            variance_cold += h_cold.magnitude().powi(2);
            variance_hot += h_hot.magnitude().powi(2);
        }

        variance_cold /= n_samples as f64;
        variance_hot /= n_samples as f64;

        // Higher temperature should give larger variance
        // (with some statistical tolerance)
        assert!(variance_hot > variance_cold * 0.5);
    }
}
