//! Spin-dependent resistors and conductances
//!
//! Models electrical resistance with spin-dependent scattering,
//! implementing the two-current model (Mott model).

/// Spin resistor with spin-dependent channels
///
/// In the two-current model, spin-up and spin-down electrons
/// experience different scattering rates, leading to spin accumulation.
#[derive(Debug, Clone)]
pub struct SpinResistor {
    /// Resistance for spin-up channel \[Ω\]
    pub r_up: f64,

    /// Resistance for spin-down channel \[Ω\]
    pub r_down: f64,

    /// Length \[m\]
    pub length: f64,

    /// Cross-sectional area \[m²\]
    pub area: f64,
}

impl SpinResistor {
    /// Create a new spin resistor
    pub fn new(r_up: f64, r_down: f64, length: f64, area: f64) -> Self {
        Self {
            r_up,
            r_down,
            length,
            area,
        }
    }

    /// Create from bulk resistivity and spin asymmetry
    ///
    /// # Arguments
    /// * `rho` - Average resistivity [Ω·m]
    /// * `beta` - Spin asymmetry parameter (0 to 1)
    /// * `length` - Length \[m\]
    /// * `area` - Cross-sectional area \[m²\]
    pub fn from_bulk(rho: f64, beta: f64, length: f64, area: f64) -> Self {
        let r_avg = rho * length / area;

        // R_↑ = R_avg × (1 - β)
        // R_↓ = R_avg × (1 + β)
        let r_up = r_avg * (1.0 - beta);
        let r_down = r_avg * (1.0 + beta);

        Self {
            r_up,
            r_down,
            length,
            area,
        }
    }

    /// Calculate total resistance (parallel combination)
    pub fn total_resistance(&self) -> f64 {
        1.0 / (1.0 / self.r_up + 1.0 / self.r_down)
    }

    /// Calculate spin polarization of current
    ///
    /// P = (I_↑ - I_↓) / (I_↑ + I_↓) = (R_↓ - R_↑) / (R_↓ + R_↑)
    pub fn current_polarization(&self) -> f64 {
        (self.r_down - self.r_up) / (self.r_down + self.r_up)
    }

    /// Calculate spin accumulation for given total current
    ///
    /// μ_s = I × R_s
    ///
    /// where R_s is the spin resistance
    pub fn spin_accumulation(&self, current: f64) -> f64 {
        let r_s = (self.r_up - self.r_down).abs() / 2.0;
        current * r_s
    }
}

/// Spin conductance (inverse of resistance)
#[derive(Debug, Clone)]
pub struct SpinConductance {
    /// Conductance for spin-up channel \[S\]
    pub g_up: f64,

    /// Conductance for spin-down channel \[S\]
    pub g_down: f64,
}

impl SpinConductance {
    /// Create from spin resistor
    pub fn from_resistor(resistor: &SpinResistor) -> Self {
        Self {
            g_up: 1.0 / resistor.r_up,
            g_down: 1.0 / resistor.r_down,
        }
    }

    /// Total conductance
    pub fn total(&self) -> f64 {
        self.g_up + self.g_down
    }

    /// Spin-dependent conductance
    pub fn spin(&self) -> f64 {
        self.g_up - self.g_down
    }

    /// Polarization
    pub fn polarization(&self) -> f64 {
        (self.g_up - self.g_down) / (self.g_up + self.g_down)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spin_resistor_creation() {
        let resistor = SpinResistor::new(100.0, 200.0, 1.0e-6, 1.0e-12);
        assert_eq!(resistor.r_up, 100.0);
        assert_eq!(resistor.r_down, 200.0);
    }

    #[test]
    fn test_from_bulk() {
        let resistor = SpinResistor::from_bulk(1.0e-7, 0.5, 1.0e-6, 1.0e-12);
        assert!(resistor.r_down > resistor.r_up);
    }

    #[test]
    fn test_total_resistance() {
        let resistor = SpinResistor::new(100.0, 200.0, 1.0e-6, 1.0e-12);
        let r_total = resistor.total_resistance();

        // Should be harmonic mean
        let expected = 1.0 / (1.0 / 100.0 + 1.0 / 200.0);
        assert!((r_total - expected).abs() < 1e-10);
    }

    #[test]
    fn test_current_polarization() {
        let resistor = SpinResistor::new(100.0, 200.0, 1.0e-6, 1.0e-12);
        let p = resistor.current_polarization();

        // P = (200-100)/(200+100) = 1/3
        assert!((p - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_conductance_conversion() {
        let resistor = SpinResistor::new(100.0, 200.0, 1.0e-6, 1.0e-12);
        let conductance = SpinConductance::from_resistor(&resistor);

        assert!((conductance.g_up - 0.01).abs() < 1e-10);
        assert!((conductance.g_down - 0.005).abs() < 1e-10);
    }
}
