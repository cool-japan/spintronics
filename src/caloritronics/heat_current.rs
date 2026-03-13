//! Unified Heat Current Calculator for Spin Caloritronics
//!
//! This module provides a comprehensive calculator for heat currents arising
//! from multiple physical mechanisms:
//!
//! 1. **Fourier conduction**: j_Q = −κ·∇T (classical thermal conduction)
//! 2. **Peltier effect**: j_Q = Π·j_c (heat drag by charge carriers)
//! 3. **Spin Peltier effect**: j_Q = Π_s·j_s (heat drag by spin current)
//!
//! In spin caloritronics, all three mechanisms can operate simultaneously and
//! their interference is important for device design.
//!
//! ## Physical Background
//!
//! The total heat current in a spintronic device is:
//!
//!   j_Q = −κ·∇T + Π·j_c + Π_s·j_s
//!
//! where:
//! - κ \[W/(m·K)\] is the total thermal conductivity (electron + phonon + magnon)
//! - Π = T·S_e \[V\] is the (charge) Peltier coefficient
//! - Π_s = T·S_s \[V\] is the spin Peltier coefficient
//!
//! The spin Peltier effect was first observed in:
//! J. Flipse et al., "Observation of the spin Peltier effect for magnetic insulators",
//! Phys. Rev. Lett. 113, 027601 (2014)

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::vector3::Vector3;

/// Unified heat current calculator combining Fourier, Peltier, and spin Peltier mechanisms.
///
/// This calculator composes the three distinct heat-current contributions that can
/// coexist in a spin caloritronic heterostructure. Each mechanism is physically
/// independent and additive in the linear-response regime.
///
/// # Example
///
/// ```rust
/// use spintronics::caloritronics::HeatCurrentCalculator;
/// use spintronics::Vector3;
///
/// // YIG/Pt system at 300 K
/// let calc = HeatCurrentCalculator::new(
///     46.0,    // thermal conductivity κ [W/(m·K)]
///     -1.5e-3, // Peltier Π = T·S_e = 300·(-5e-6) [V]
///     0.3,     // spin Peltier Π_s = T·S_s = 300·1e-3 [V]
/// );
///
/// let grad_t  = Vector3::new(1000.0, 0.0, 0.0);  // 1000 K/m
/// let j_charge = Vector3::new(1.0e6, 0.0, 0.0);  // 1 MA/m²
/// let j_spin   = Vector3::new(0.0, 0.0, 1.0e3);  // small spin current
///
/// let j_q_total = calc.total_heat_current(&grad_t, &j_charge, &j_spin);
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HeatCurrentCalculator {
    /// Thermal conductivity κ [W/(m·K)]
    ///
    /// Includes electronic, phononic, and magnonic contributions.
    pub thermal_conductivity: f64,

    /// (Charge) Peltier coefficient Π \[V\]
    ///
    /// Related to Seebeck coefficient by Kelvin relation: Π = T·S_e
    pub peltier_coeff: f64,

    /// Spin Peltier coefficient Π_s \[V\]
    ///
    /// Reciprocal of the spin Seebeck effect via Kelvin relation: Π_s = T·S_s.
    /// Represents heat transported per unit spin current.
    pub spin_peltier_coeff: f64,
}

impl HeatCurrentCalculator {
    /// Create a new heat current calculator.
    ///
    /// # Arguments
    /// * `thermal_conductivity` - κ [W/(m·K)]
    /// * `peltier_coeff` - Π = T·S_e \[V\]
    /// * `spin_peltier_coeff` - Π_s = T·S_s \[V\]
    pub fn new(thermal_conductivity: f64, peltier_coeff: f64, spin_peltier_coeff: f64) -> Self {
        Self {
            thermal_conductivity,
            peltier_coeff,
            spin_peltier_coeff,
        }
    }

    /// Compute the Fourier (diffusive) heat current.
    ///
    /// Fourier's law of heat conduction:
    ///   j_Q^Fourier = −κ · ∇T   [W/m²]
    ///
    /// The minus sign reflects that heat flows from hot to cold (down the gradient).
    ///
    /// # Arguments
    /// * `grad_t` - Temperature gradient [K/m]
    #[inline]
    pub fn fourier_current(&self, grad_t: &Vector3<f64>) -> Vector3<f64> {
        Vector3::new(
            -self.thermal_conductivity * grad_t.x,
            -self.thermal_conductivity * grad_t.y,
            -self.thermal_conductivity * grad_t.z,
        )
    }

    /// Compute the (charge) Peltier heat current.
    ///
    /// The Peltier effect causes charge carriers to transport heat as they move:
    ///   j_Q^Peltier = Π · j_c   [W/m²]
    ///
    /// This is the reciprocal of the Seebeck effect (by Kelvin relation).
    ///
    /// # Arguments
    /// * `j_charge` - Charge current density [A/m²]
    #[inline]
    pub fn peltier_current(&self, j_charge: &Vector3<f64>) -> Vector3<f64> {
        Vector3::new(
            self.peltier_coeff * j_charge.x,
            self.peltier_coeff * j_charge.y,
            self.peltier_coeff * j_charge.z,
        )
    }

    /// Compute the spin Peltier heat current.
    ///
    /// Spin carriers (magnons, conduction electrons with net spin polarization)
    /// transport heat proportional to their spin current:
    ///   j_Q^sPeltier = Π_s · j_s   [W/m²]
    ///
    /// This is the reciprocal of the spin Seebeck effect (Kelvin relation).
    ///
    /// # Arguments
    /// * `j_spin` - Spin current density [A/m²]
    #[inline]
    pub fn spin_peltier_current(&self, j_spin: &Vector3<f64>) -> Vector3<f64> {
        Vector3::new(
            self.spin_peltier_coeff * j_spin.x,
            self.spin_peltier_coeff * j_spin.y,
            self.spin_peltier_coeff * j_spin.z,
        )
    }

    /// Compute the total heat current from all three mechanisms.
    ///
    /// Superposition of Fourier, Peltier, and spin Peltier contributions:
    ///   j_Q = −κ·∇T + Π·j_c + Π_s·j_s   [W/m²]
    ///
    /// # Arguments
    /// * `grad_t` - Temperature gradient [K/m]
    /// * `j_charge` - Charge current density [A/m²]
    /// * `j_spin` - Spin current density [A/m²]
    pub fn total_heat_current(
        &self,
        grad_t: &Vector3<f64>,
        j_charge: &Vector3<f64>,
        j_spin: &Vector3<f64>,
    ) -> Vector3<f64> {
        let fourier = self.fourier_current(grad_t);
        let peltier = self.peltier_current(j_charge);
        let spin_peltier = self.spin_peltier_current(j_spin);
        Vector3::new(
            fourier.x + peltier.x + spin_peltier.x,
            fourier.y + peltier.y + spin_peltier.y,
            fourier.z + peltier.z + spin_peltier.z,
        )
    }
}
