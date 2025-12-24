//! Spin Peltier Effect
//!
//! The spin Peltier effect is the reciprocal of the spin Seebeck effect:
//! a spin current induces a heat current at a ferromagnet/normal-metal interface.
//!
//! Q = Π_s × J_s
//!
//! where Π_s is the spin Peltier coefficient and J_s is the spin current.

use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Spin Peltier Effect
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpinPeltier {
    /// Spin Peltier coefficient [W·m²/J]
    ///
    /// Related to spin Seebeck coefficient by Onsager reciprocity:
    /// Π_s = T × S_s
    pub pi_s: f64,

    /// Temperature \[K\]
    pub temperature: f64,

    /// Interface area \[m²\]
    pub area: f64,
}

impl Default for SpinPeltier {
    fn default() -> Self {
        Self {
            pi_s: 1.0e-6, // Typical value
            temperature: 300.0,
            area: 1.0e-12, // 1 μm²
        }
    }
}

impl SpinPeltier {
    /// Create spin Peltier properties for YIG/Pt interface
    pub fn yig_pt(temperature: f64) -> Self {
        // Using Onsager relation: Π_s = T × S_s
        let seebeck_spin = 1.0e-6; // Spin Seebeck coefficient
        Self {
            pi_s: temperature * seebeck_spin,
            temperature,
            area: 1.0e-12,
        }
    }

    /// Calculate heat current from spin current
    ///
    /// # Arguments
    /// * `js_magnitude` - Magnitude of spin current density \[J/m²\]
    ///
    /// # Returns
    /// Heat current \[W\]
    #[inline]
    pub fn heat_current(&self, js_magnitude: f64) -> f64 {
        self.pi_s * js_magnitude * self.area
    }

    /// Calculate temperature change rate
    ///
    /// # Arguments
    /// * `js_magnitude` - Spin current density \[J/m²\]
    /// * `heat_capacity` - Volumetric heat capacity [J/(m³·K)]
    /// * `volume` - Volume \[m³\]
    ///
    /// # Returns
    /// dT/dt [K/s]
    #[inline]
    pub fn temperature_change_rate(
        &self,
        js_magnitude: f64,
        heat_capacity: f64,
        volume: f64,
    ) -> f64 {
        let q = self.heat_current(js_magnitude);
        q / (heat_capacity * volume)
    }

    /// Update Peltier coefficient for new temperature (Onsager relation)
    pub fn update_temperature(&mut self, new_temperature: f64) {
        let seebeck_spin = self.pi_s / self.temperature;
        self.temperature = new_temperature;
        self.pi_s = new_temperature * seebeck_spin;
    }

    /// Builder method to set spin Peltier coefficient
    pub fn with_pi_s(mut self, pi_s: f64) -> Self {
        self.pi_s = pi_s;
        self
    }

    /// Builder method to set temperature
    pub fn with_temperature(mut self, temperature: f64) -> Self {
        self.temperature = temperature;
        self
    }

    /// Builder method to set interface area
    pub fn with_area(mut self, area: f64) -> Self {
        self.area = area;
        self
    }
}

impl fmt::Display for SpinPeltier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SpinPeltier: Π_s={:.2e} W·m²/J, T={:.0} K, A={:.2e} m²",
            self.pi_s, self.temperature, self.area
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_heat_current_proportional_to_spin_current() {
        let peltier = SpinPeltier::default();

        let q1 = peltier.heat_current(1.0e10);
        let q2 = peltier.heat_current(2.0e10);

        assert!((q2 / q1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_zero_spin_current() {
        let peltier = SpinPeltier::default();
        let q = peltier.heat_current(0.0);
        assert!(q.abs() < 1e-30);
    }

    #[test]
    fn test_onsager_relation() {
        let t1 = 300.0;
        let t2 = 400.0;
        let peltier1 = SpinPeltier::yig_pt(t1);
        let peltier2 = SpinPeltier::yig_pt(t2);

        // Π_s should scale with temperature
        assert!((peltier2.pi_s / peltier1.pi_s - t2 / t1).abs() < 1e-10);
    }

    #[test]
    fn test_temperature_update() {
        let mut peltier = SpinPeltier::yig_pt(300.0);
        let original_ratio = peltier.pi_s / peltier.temperature;

        peltier.update_temperature(400.0);

        let new_ratio = peltier.pi_s / peltier.temperature;
        assert!((new_ratio - original_ratio).abs() < 1e-15);
    }
}
