//! Thermal magnon transport
//!
//! Magnons (spin waves) can carry both angular momentum and heat.
//! This module implements thermal transport by magnons in magnetic insulators.

use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::constants::KB;

/// Magnon thermal conductivity
///
/// In magnetic insulators like YIG, magnons are the primary heat carriers
/// at low temperatures.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MagnonThermalConductivity {
    /// Magnon thermal conductivity \[W/(m·K)\]
    pub kappa_magnon: f64,

    /// Temperature \[K\]
    pub temperature: f64,

    /// Magnon mean free path \[m\]
    pub mean_free_path: f64,
}

impl Default for MagnonThermalConductivity {
    fn default() -> Self {
        Self {
            kappa_magnon: 1.0, // W/(m·K) for YIG at room temperature
            temperature: 300.0,
            mean_free_path: 1.0e-6,
        }
    }
}

impl MagnonThermalConductivity {
    /// Calculate thermal conductivity at given temperature
    ///
    /// For magnons, κ(T) typically follows T² or T³ behavior at low T
    #[inline]
    pub fn conductivity_at_temperature(&self, temp: f64) -> f64 {
        if temp < 1.0 {
            return 0.0;
        }

        // Simplified model: κ ∝ T² at low T, saturates at high T
        let t_ratio = temp / self.temperature;

        if temp < 50.0 {
            // Low temperature: T² dependence
            self.kappa_magnon * t_ratio.powi(2)
        } else {
            // High temperature: approaches constant
            self.kappa_magnon * (1.0 - (-t_ratio).exp())
        }
    }
}

/// Thermal magnon transport calculator
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ThermalMagnonTransport {
    /// Magnon thermal conductivity calculator
    pub conductivity: MagnonThermalConductivity,

    /// Spin Seebeck coefficient \[V/K\]
    pub seebeck_coefficient: f64,
}

impl Default for ThermalMagnonTransport {
    fn default() -> Self {
        Self {
            conductivity: MagnonThermalConductivity::default(),
            seebeck_coefficient: 1.0e-6, // 1 μV/K
        }
    }
}

impl ThermalMagnonTransport {
    /// Create transport properties for YIG
    pub fn yig() -> Self {
        Self {
            conductivity: MagnonThermalConductivity {
                kappa_magnon: 1.0,
                temperature: 300.0,
                mean_free_path: 10.0e-6, // 10 μm (YIG has long mean free path)
            },
            seebeck_coefficient: 1.0e-6,
        }
    }

    /// Calculate heat flux from temperature gradient
    ///
    /// Q = -κ ∇T
    ///
    /// # Arguments
    /// * `grad_t_magnitude` - Magnitude of temperature gradient \[K/m\]
    /// * `temperature` - Temperature \[K\]
    ///
    /// # Returns
    /// Heat flux \[W/m²\]
    #[inline]
    pub fn heat_flux(&self, grad_t_magnitude: f64, temperature: f64) -> f64 {
        let kappa = self.conductivity.conductivity_at_temperature(temperature);
        kappa * grad_t_magnitude
    }

    /// Calculate magnon chemical potential from temperature gradient
    ///
    /// μ_m = S_m × ∇T
    ///
    /// where S_m is the magnon Seebeck coefficient
    #[inline]
    pub fn magnon_chemical_potential(&self, grad_t_magnitude: f64) -> f64 {
        self.seebeck_coefficient * grad_t_magnitude
    }

    /// Calculate thermal magnon accumulation
    ///
    /// Δμ_m / k_B T gives the effective magnon population change
    #[inline]
    pub fn thermal_magnon_accumulation(&self, grad_t_magnitude: f64, temperature: f64) -> f64 {
        if temperature < 1.0 {
            return 0.0;
        }

        let mu_m = self.magnon_chemical_potential(grad_t_magnitude);
        mu_m / (KB * temperature)
    }
}

impl fmt::Display for MagnonThermalConductivity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "MagnonThermalConductivity: κ={:.2} W/(m·K), T={:.0} K, λ={:.2} μm",
            self.kappa_magnon,
            self.temperature,
            self.mean_free_path * 1e6
        )
    }
}

impl fmt::Display for ThermalMagnonTransport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "ThermalMagnonTransport: S={:.2e} V/K, {}",
            self.seebeck_coefficient, self.conductivity
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thermal_conductivity_temperature_dependence() {
        let kappa = MagnonThermalConductivity::default();

        let k_low = kappa.conductivity_at_temperature(10.0);
        let k_mid = kappa.conductivity_at_temperature(100.0);
        let k_high = kappa.conductivity_at_temperature(300.0);

        // Should increase with temperature
        assert!(k_low < k_mid);
        assert!(k_mid < k_high);
    }

    #[test]
    fn test_heat_flux_proportional_to_gradient() {
        let transport = ThermalMagnonTransport::yig();

        let q1 = transport.heat_flux(1000.0, 300.0);
        let q2 = transport.heat_flux(2000.0, 300.0);

        // Double gradient should double flux
        assert!((q2 / q1 - 2.0).abs() < 0.1);
    }

    #[test]
    fn test_zero_gradient_zero_flux() {
        let transport = ThermalMagnonTransport::default();
        let q = transport.heat_flux(0.0, 300.0);
        assert!(q.abs() < 1e-30);
    }

    #[test]
    fn test_magnon_accumulation() {
        let transport = ThermalMagnonTransport::yig();
        let accumulation = transport.thermal_magnon_accumulation(1.0e6, 300.0);
        assert!(accumulation > 0.0);
    }
}
