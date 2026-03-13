//! Orbital Hall Effect (OHE) Module
//!
//! This module implements the Orbital Hall Effect, where charge current generates
//! a transverse orbital angular momentum (OAM) current, analogous to the Spin Hall
//! Effect but involving orbital degrees of freedom instead of spin.
//!
//! ## Physics Background
//!
//! The Orbital Hall Effect arises from the coupling between electron orbital motion
//! and crystal momentum. Unlike the Spin Hall Effect which requires strong spin-orbit
//! coupling (SOC), the OHE can be giant even in light metals with weak SOC.
//!
//! Key insight: In materials like Cr and Ti, the orbital Hall conductivity can be
//! orders of magnitude larger than the spin Hall conductivity of heavy metals like Pt,
//! despite having much weaker SOC. This is because the OHE originates from orbital
//! texture in momentum space, which does not require SOC.
//!
//! ## Orbital-to-Spin Conversion
//!
//! At interfaces with heavy metals or materials with strong SOC, orbital currents
//! can be converted into spin currents via the orbital-to-spin conversion mechanism:
//!
//! J_s = eta_LS * J_L
//!
//! where eta_LS is the conversion efficiency that depends on the SOC strength of
//! the adjacent material.
//!
//! ## Key References
//!
//! - H. Kontani et al., "Giant Orbital Hall Effect in Transition Metals",
//!   Phys. Rev. Lett. 102, 016601 (2009)
//! - D. Go et al., "Intrinsic Spin and Orbital Hall Effects from Orbital Texture",
//!   Phys. Rev. Lett. 121, 086602 (2018)
//! - S. Lee et al., "Efficient conversion of orbital Hall current to spin current
//!   for spin-orbit torque switching", Commun. Phys. 4, 234 (2021)
//! - D. Jo et al., "Gigantic intrinsic orbital Hall effects in weakly spin-orbit
//!   coupled metals", Phys. Rev. B 98, 214405 (2018)

use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::constants::{E_CHARGE, HBAR};
use crate::error::{Error, Result};
use crate::vector3::Vector3;

// =============================================================================
// Material Database
// =============================================================================

/// Orbital Hall material parameters
///
/// Contains the physical properties relevant to the Orbital Hall Effect
/// for a given material. The orbital Hall conductivity sigma_oh is given
/// in units of (hbar/e) (Ohm*cm)^-1.
///
/// # Material Categories
///
/// - **Light metals** (Cr, Ti, V): Giant OHE despite weak SOC
/// - **Noble metals** (Cu, Al): Moderate OHE
/// - **Heavy metals** (Pt, W): OHE coexists with strong SHE
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OrbitalHallMaterial {
    /// Material name
    pub name: &'static str,
    /// Orbital Hall conductivity in units of (hbar/e) (Ohm*cm)^-1
    pub sigma_oh: f64,
    /// Electrical resistivity \[Ohm*m\]
    pub resistivity: f64,
    /// Spin-orbit coupling strength (relative, dimensionless)
    /// Normalized such that Pt = 1.0
    pub soc_strength: f64,
    /// Orbital Hall angle theta_OH = sigma_OH / sigma_charge (dimensionless)
    pub orbital_hall_angle: f64,
}

impl OrbitalHallMaterial {
    /// Create a new OrbitalHallMaterial with parameter validation
    ///
    /// # Arguments
    /// * `name` - Material identifier
    /// * `sigma_oh` - Orbital Hall conductivity [(hbar/e)(Ohm*cm)^-1]
    /// * `resistivity` - Electrical resistivity \[Ohm*m\]
    /// * `soc_strength` - Relative SOC strength (dimensionless)
    /// * `orbital_hall_angle` - Orbital Hall angle (dimensionless)
    ///
    /// # Errors
    /// Returns error if resistivity is non-positive or soc_strength is negative
    pub fn new(
        name: &'static str,
        sigma_oh: f64,
        resistivity: f64,
        soc_strength: f64,
        orbital_hall_angle: f64,
    ) -> Result<Self> {
        if resistivity <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "resistivity".to_string(),
                reason: "must be positive".to_string(),
            });
        }
        if soc_strength < 0.0 {
            return Err(Error::InvalidParameter {
                param: "soc_strength".to_string(),
                reason: "must be non-negative".to_string(),
            });
        }
        Ok(Self {
            name,
            sigma_oh,
            resistivity,
            soc_strength,
            orbital_hall_angle,
        })
    }

    /// Chromium: Giant OHE in a light 3d transition metal
    ///
    /// Cr exhibits one of the largest known orbital Hall conductivities
    /// (~5000 (hbar/e)(Ohm*cm)^-1) despite having very weak spin-orbit coupling.
    /// This demonstrates that the OHE originates from orbital texture, not SOC.
    ///
    /// Reference: Jo et al., Phys. Rev. B 98, 214405 (2018)
    pub fn chromium() -> Self {
        Self {
            name: "Cr",
            sigma_oh: 5000.0,
            resistivity: 1.29e-7,
            soc_strength: 0.04,
            orbital_hall_angle: 0.065,
        }
    }

    /// Titanium: Large OHE in a light metal
    ///
    /// Ti has a large orbital Hall conductivity (~3000 (hbar/e)(Ohm*cm)^-1)
    /// with very weak SOC, making it attractive for orbitronic applications.
    pub fn titanium() -> Self {
        Self {
            name: "Ti",
            sigma_oh: 3000.0,
            resistivity: 4.20e-7,
            soc_strength: 0.02,
            orbital_hall_angle: 0.126,
        }
    }

    /// Vanadium: Strong OHE in a 3d metal
    ///
    /// V shows significant orbital Hall conductivity arising from
    /// its d-orbital texture in momentum space.
    pub fn vanadium() -> Self {
        Self {
            name: "V",
            sigma_oh: 2500.0,
            resistivity: 1.97e-7,
            soc_strength: 0.03,
            orbital_hall_angle: 0.049,
        }
    }

    /// Copper: Moderate OHE in a noble metal
    ///
    /// Cu has a moderate orbital Hall conductivity but excellent
    /// electrical conductivity, making it useful as an orbital current source
    /// in combination with strong-SOC conversion layers.
    pub fn copper() -> Self {
        Self {
            name: "Cu",
            sigma_oh: 1500.0,
            resistivity: 1.68e-8,
            soc_strength: 0.01,
            orbital_hall_angle: 0.003,
        }
    }

    /// Aluminum: Light-metal OHE source
    ///
    /// Al provides orbital currents with negligible SOC,
    /// serving as a testbed for pure orbital transport studies.
    pub fn aluminum() -> Self {
        Self {
            name: "Al",
            sigma_oh: 800.0,
            resistivity: 2.65e-8,
            soc_strength: 0.005,
            orbital_hall_angle: 0.002,
        }
    }

    /// Platinum: Heavy metal with both OHE and SHE
    ///
    /// In Pt, both orbital and spin Hall effects coexist.
    /// The spin Hall conductivity is ~2000 (hbar/e)(Ohm*cm)^-1,
    /// while the orbital Hall conductivity adds on top of that.
    /// The strong SOC in Pt efficiently converts orbital to spin current internally.
    pub fn platinum() -> Self {
        Self {
            name: "Pt",
            sigma_oh: 4000.0,
            resistivity: 1.06e-7,
            soc_strength: 1.0,
            orbital_hall_angle: 0.042,
        }
    }

    /// Tungsten: Heavy metal with large OHE and negative SHE
    ///
    /// W (beta-phase) has both large OHE and large negative SHE.
    /// The interplay between orbital and spin channels leads to
    /// complex torque symmetries.
    pub fn tungsten() -> Self {
        Self {
            name: "W",
            sigma_oh: 3500.0,
            resistivity: 2.0e-7,
            soc_strength: 0.8,
            orbital_hall_angle: 0.070,
        }
    }

    /// Manganese: Antiferromagnetic metal with OHE
    ///
    /// Mn exhibits significant OHE. In antiferromagnetic ordering,
    /// the orbital current may couple to the Neel order parameter.
    pub fn manganese() -> Self {
        Self {
            name: "Mn",
            sigma_oh: 2000.0,
            resistivity: 1.44e-6,
            soc_strength: 0.05,
            orbital_hall_angle: 0.288,
        }
    }

    /// Get the charge conductivity [1/(Ohm*m)]
    pub fn charge_conductivity(&self) -> f64 {
        1.0 / self.resistivity
    }

    /// Check if this material is a light metal (weak SOC)
    ///
    /// Light metals are defined as having SOC strength < 0.1 relative to Pt.
    /// These are particularly interesting for orbitronics as they demonstrate
    /// that giant OHE does not require strong SOC.
    pub fn is_light_metal(&self) -> bool {
        self.soc_strength < 0.1
    }

    /// Compute the ratio of orbital Hall conductivity to Pt's spin Hall conductivity
    ///
    /// This ratio quantifies how much larger the OHE can be compared to
    /// the conventional SHE in Pt (sigma_SH_Pt ~ 2000 (hbar/e)(Ohm*cm)^-1).
    /// Values > 1.0 indicate the OHE exceeds Pt's SHE.
    pub fn ohe_to_pt_she_ratio(&self) -> f64 {
        let sigma_sh_pt = 2000.0; // Pt spin Hall conductivity
        self.sigma_oh / sigma_sh_pt
    }
}

impl fmt::Display for OrbitalHallMaterial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "OrbitalHallMaterial({}: sigma_OH={:.0} (hbar/e)(Ohm*cm)^-1, \
             rho={:.2e} Ohm*m, SOC={:.3}, theta_OH={:.4})",
            self.name, self.sigma_oh, self.resistivity, self.soc_strength, self.orbital_hall_angle
        )
    }
}

// =============================================================================
// Orbital-to-Spin Converter
// =============================================================================

/// Orbital-to-spin conversion at interfaces
///
/// When an orbital current reaches an interface with a material possessing
/// strong spin-orbit coupling, the orbital angular momentum is partially
/// converted into spin angular momentum. The efficiency of this conversion
/// depends on the SOC strength and interfacial properties.
///
/// The conversion follows: J_s = eta_LS * J_L
///
/// where:
/// - J_s is the resulting spin current density
/// - eta_LS is the orbital-to-spin conversion efficiency (-1 <= eta_LS <= 1)
/// - J_L is the incoming orbital current density
///
/// # Interface Types
///
/// - **Light metal / Heavy metal**: e.g., Cu/Pt, Ti/Pt - efficient conversion
/// - **Light metal / Oxide**: e.g., Cu/BiOx - conversion via Rashba SOC
/// - **Light metal / Topological insulator**: e.g., Cu/Bi2Se3 - surface SOC
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OrbitalToSpinConverter {
    /// Conversion efficiency eta_LS (dimensionless, -1 to 1)
    pub conversion_efficiency: f64,
    /// Description of the interface type
    pub interface_type: &'static str,
}

impl OrbitalToSpinConverter {
    /// Create a new converter with validation
    ///
    /// # Arguments
    /// * `conversion_efficiency` - eta_LS, must be in [-1, 1]
    /// * `interface_type` - Description of the interface
    ///
    /// # Errors
    /// Returns error if efficiency is outside [-1, 1]
    pub fn new(conversion_efficiency: f64, interface_type: &'static str) -> Result<Self> {
        if conversion_efficiency.abs() > 1.0 {
            return Err(Error::InvalidParameter {
                param: "conversion_efficiency".to_string(),
                reason: "must be between -1.0 and 1.0".to_string(),
            });
        }
        Ok(Self {
            conversion_efficiency,
            interface_type,
        })
    }

    /// Cu/Pt interface: efficient orbital-to-spin conversion
    ///
    /// Pt's strong SOC provides high conversion efficiency (~0.6).
    /// This is one of the most studied L-S conversion interfaces.
    pub fn cu_pt() -> Self {
        Self {
            conversion_efficiency: 0.6,
            interface_type: "Cu/Pt (heavy metal SOC)",
        }
    }

    /// Ti/Pt interface: light metal source with heavy metal converter
    ///
    /// Ti generates large orbital current, Pt converts it to spin.
    /// The combined system can produce torques exceeding those from
    /// Pt alone via SHE.
    pub fn ti_pt() -> Self {
        Self {
            conversion_efficiency: 0.55,
            interface_type: "Ti/Pt (light/heavy metal)",
        }
    }

    /// Cr/Pt interface: giant OHE source with efficient conversion
    pub fn cr_pt() -> Self {
        Self {
            conversion_efficiency: 0.5,
            interface_type: "Cr/Pt (giant OHE/heavy metal)",
        }
    }

    /// Cu/BiOx interface: oxide Rashba conversion
    ///
    /// The Bi-oxide interface provides Rashba-type SOC for conversion.
    /// Lower efficiency than metallic interfaces but interesting for
    /// oxide-based devices.
    pub fn cu_biox() -> Self {
        Self {
            conversion_efficiency: 0.25,
            interface_type: "Cu/BiOx (Rashba SOC)",
        }
    }

    /// Cu/Bi2Se3 interface: topological insulator conversion
    ///
    /// The topological surface states of Bi2Se3 provide efficient
    /// orbital-to-spin conversion via surface SOC.
    pub fn cu_bi2se3() -> Self {
        Self {
            conversion_efficiency: 0.7,
            interface_type: "Cu/Bi2Se3 (TI surface SOC)",
        }
    }

    /// Weak conversion (no significant SOC at interface)
    pub fn weak() -> Self {
        Self {
            conversion_efficiency: 0.05,
            interface_type: "weak SOC interface",
        }
    }

    /// Convert orbital current density to spin current density
    ///
    /// Applies the conversion: J_s = eta_LS * J_L
    ///
    /// # Arguments
    /// * `j_orbital` - Orbital current density vector [A/m^2] (in units of hbar)
    ///
    /// # Returns
    /// Spin current density vector [A/m^2] (in units of hbar)
    #[inline]
    pub fn convert_to_spin(&self, j_orbital: Vector3<f64>) -> Vector3<f64> {
        j_orbital * self.conversion_efficiency
    }
}

impl fmt::Display for OrbitalToSpinConverter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "OrbitalToSpinConverter(eta_LS={:.3}, type={})",
            self.conversion_efficiency, self.interface_type
        )
    }
}

// =============================================================================
// Orbital Hall Effect Calculator
// =============================================================================

/// Orbital Hall Effect calculator
///
/// Computes orbital current densities and related quantities from the
/// Orbital Hall Effect. The OHE generates a transverse orbital angular
/// momentum current from a longitudinal charge current:
///
/// J_L = sigma_OH * (E x z_hat)
///
/// where E is the applied electric field and z_hat is the interface normal.
///
/// # Example
///
/// ```
/// use spintronics::orbitronics::orbital_hall::{
///     OrbitalHallEffect, OrbitalHallMaterial, OrbitalToSpinConverter,
/// };
/// use spintronics::vector3::Vector3;
///
/// let cr = OrbitalHallMaterial::chromium();
/// let converter = OrbitalToSpinConverter::cr_pt();
/// let ohe = OrbitalHallEffect::new(cr, converter);
///
/// // Apply electric field along x
/// let e_field = Vector3::new(1.0e5, 0.0, 0.0); // V/m
/// let j_orbital = ohe.orbital_current_density(e_field);
///
/// // Orbital current flows along y (transverse)
/// assert!(j_orbital.y.abs() > j_orbital.x.abs());
/// ```
#[derive(Debug, Clone)]
pub struct OrbitalHallEffect {
    /// Material providing the orbital Hall effect
    pub material: OrbitalHallMaterial,
    /// Interface converter for orbital-to-spin conversion
    pub converter: OrbitalToSpinConverter,
}

impl OrbitalHallEffect {
    /// Create a new OrbitalHallEffect system
    pub fn new(material: OrbitalHallMaterial, converter: OrbitalToSpinConverter) -> Self {
        Self {
            material,
            converter,
        }
    }

    /// Calculate the orbital current density from the Orbital Hall Effect
    ///
    /// J_L = sigma_OH * (E x z_hat)
    ///
    /// The orbital current flows transverse to both the electric field and
    /// the interface normal (z-direction). The orbital current density is
    /// in units of (hbar / (2e)) * [A/m^2].
    ///
    /// # Arguments
    /// * `e_field` - Applied electric field vector \[V/m\]
    ///
    /// # Returns
    /// Orbital current density vector in units of (hbar/(2e)) * [A/m^2]
    ///
    /// # Note
    /// sigma_oh is in (hbar/e)(Ohm*cm)^-1 so we convert:
    /// J_L [A/m^2 * hbar/(2e)] = sigma_oh [hbar/e / (Ohm*cm)] * E \[V/m\]
    /// Need unit conversion: 1/(Ohm*cm) = 100/(Ohm*m)
    #[inline]
    pub fn orbital_current_density(&self, e_field: Vector3<f64>) -> Vector3<f64> {
        let z_hat = Vector3::unit_z();
        let cross = e_field.cross(&z_hat);
        // Convert sigma_oh from (hbar/e)/(Ohm*cm) to SI: multiply by 100
        // J_L = sigma_oh * (E x z_hat)
        cross * (self.material.sigma_oh * 100.0)
    }

    /// Calculate orbital current density from charge current density
    ///
    /// Uses J_L = theta_OH * J_c x z_hat where theta_OH = sigma_OH / sigma_charge
    ///
    /// # Arguments
    /// * `j_charge` - Charge current density magnitude [A/m^2]
    /// * `current_direction` - Unit vector of current flow direction
    ///
    /// # Returns
    /// Orbital angular momentum current density vector
    #[inline]
    pub fn orbital_current_from_charge(
        &self,
        j_charge: f64,
        current_direction: Vector3<f64>,
    ) -> Vector3<f64> {
        let z_hat = Vector3::unit_z();
        let j_vec = current_direction * j_charge;
        let cross = j_vec.cross(&z_hat);
        // Use orbital Hall angle for direct current-to-orbital conversion
        // J_L = theta_OH * (J_c x z_hat) in units of hbar/(2e)
        cross * self.material.orbital_hall_angle
    }

    /// Calculate the spin current after orbital-to-spin conversion
    ///
    /// Two-step process:
    /// 1. OHE generates orbital current: J_L = sigma_OH * (E x z_hat)
    /// 2. Interface converts to spin: J_s = eta_LS * J_L
    ///
    /// # Arguments
    /// * `e_field` - Applied electric field vector \[V/m\]
    ///
    /// # Returns
    /// Effective spin current density vector
    pub fn effective_spin_current(&self, e_field: Vector3<f64>) -> Vector3<f64> {
        let j_orbital = self.orbital_current_density(e_field);
        self.converter.convert_to_spin(j_orbital)
    }

    /// Calculate the effective spin Hall conductivity from OHE + conversion
    ///
    /// sigma_SH_eff = eta_LS * sigma_OH
    ///
    /// This effective conductivity can exceed the intrinsic SHE of heavy metals
    /// when the OHE source material has large sigma_OH and the conversion
    /// efficiency is reasonable.
    pub fn effective_spin_hall_conductivity(&self) -> f64 {
        self.converter.conversion_efficiency * self.material.sigma_oh
    }

    /// Compare effective OHE-based spin current generation with direct SHE in Pt
    ///
    /// Returns the ratio sigma_SH_eff(OHE) / sigma_SH(Pt).
    /// Values > 1 indicate the OHE pathway is more efficient.
    pub fn enhancement_over_pt_she(&self) -> f64 {
        let sigma_sh_pt = 2000.0; // Pt SHE conductivity [(hbar/e)(Ohm*cm)^-1]
        self.effective_spin_hall_conductivity() / sigma_sh_pt
    }

    /// Compute the orbital current magnitude for a given charge current
    ///
    /// Returns |J_L| = |theta_OH| * |J_c| in natural units (hbar/(2e) * A/m^2)
    ///
    /// # Arguments
    /// * `j_charge` - Charge current density magnitude [A/m^2]
    pub fn orbital_current_magnitude(&self, j_charge: f64) -> f64 {
        self.material.orbital_hall_angle.abs() * j_charge.abs()
    }

    /// Estimate the orbital accumulation at the interface
    ///
    /// The orbital accumulation is analogous to spin accumulation and is given by:
    /// L_acc = (sigma_OH * E * lambda_L) / sigma_charge
    ///
    /// where lambda_L is the orbital diffusion length (assumed ~10 nm for metals).
    ///
    /// # Arguments
    /// * `e_field_magnitude` - Electric field magnitude \[V/m\]
    /// * `orbital_diffusion_length` - Orbital diffusion length \[m\] (typically ~10 nm)
    ///
    /// # Returns
    /// Orbital accumulation in units of hbar
    pub fn orbital_accumulation(
        &self,
        e_field_magnitude: f64,
        orbital_diffusion_length: f64,
    ) -> f64 {
        // sigma_oh in (hbar/e)/(Ohm*cm), convert to SI
        let sigma_oh_si = self.material.sigma_oh * 100.0;
        let sigma_charge = self.material.charge_conductivity();
        sigma_oh_si * e_field_magnitude * orbital_diffusion_length / sigma_charge
    }
}

impl fmt::Display for OrbitalHallEffect {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "OrbitalHallEffect(material={}, converter={})",
            self.material.name, self.converter.interface_type
        )
    }
}

// =============================================================================
// Predefined OHE Systems
// =============================================================================

/// Create a Cr/Pt OHE system (giant orbital Hall source + efficient converter)
pub fn cr_pt_system() -> OrbitalHallEffect {
    OrbitalHallEffect::new(
        OrbitalHallMaterial::chromium(),
        OrbitalToSpinConverter::cr_pt(),
    )
}

/// Create a Ti/Pt OHE system (large orbital Hall + efficient converter)
pub fn ti_pt_system() -> OrbitalHallEffect {
    OrbitalHallEffect::new(
        OrbitalHallMaterial::titanium(),
        OrbitalToSpinConverter::ti_pt(),
    )
}

/// Create a Cu/Pt OHE system (moderate orbital Hall + efficient converter)
pub fn cu_pt_system() -> OrbitalHallEffect {
    OrbitalHallEffect::new(
        OrbitalHallMaterial::copper(),
        OrbitalToSpinConverter::cu_pt(),
    )
}

// =============================================================================
// Utility: Material comparison
// =============================================================================

/// Compare orbital Hall conductivities across a standard set of materials
///
/// Returns a sorted list of (material_name, sigma_OH) pairs in descending order.
pub fn material_ohe_ranking() -> Vec<(&'static str, f64)> {
    let materials = [
        OrbitalHallMaterial::chromium(),
        OrbitalHallMaterial::titanium(),
        OrbitalHallMaterial::vanadium(),
        OrbitalHallMaterial::copper(),
        OrbitalHallMaterial::aluminum(),
        OrbitalHallMaterial::platinum(),
        OrbitalHallMaterial::tungsten(),
        OrbitalHallMaterial::manganese(),
    ];
    let mut ranking: Vec<(&'static str, f64)> =
        materials.iter().map(|m| (m.name, m.sigma_oh)).collect();
    ranking.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    ranking
}

/// Convert orbital Hall conductivity units
///
/// Converts from (hbar/e)(Ohm*cm)^-1 to SI: (hbar/(2e)) * S/m
///
/// # Arguments
/// * `sigma_oh_cgs` - Orbital Hall conductivity in (hbar/e)(Ohm*cm)^-1
///
/// # Returns
/// Orbital Hall conductivity in SI units [(hbar/(2e)) * S/m]
pub fn sigma_oh_to_si(sigma_oh_cgs: f64) -> f64 {
    // 1/(Ohm*cm) = 100 S/m, and we include the hbar/(2e) prefactor
    sigma_oh_cgs * 100.0 * HBAR / (2.0 * E_CHARGE)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cr_giant_ohe() {
        let cr = OrbitalHallMaterial::chromium();
        // Cr should have giant OHE > 4000 (hbar/e)(Ohm*cm)^-1
        assert!(
            cr.sigma_oh > 4000.0,
            "Cr sigma_OH should be giant, got {}",
            cr.sigma_oh
        );
        // Cr is a light metal with weak SOC
        assert!(
            cr.is_light_metal(),
            "Cr should be classified as light metal"
        );
        // Cr's OHE should exceed Pt's SHE
        assert!(
            cr.ohe_to_pt_she_ratio() > 1.0,
            "Cr OHE should exceed Pt SHE"
        );
    }

    #[test]
    fn test_ti_large_ohe() {
        let ti = OrbitalHallMaterial::titanium();
        // Ti has large OHE
        assert!(
            ti.sigma_oh > 2000.0,
            "Ti sigma_OH should be large, got {}",
            ti.sigma_oh
        );
        // Ti is a light metal
        assert!(
            ti.is_light_metal(),
            "Ti should be classified as light metal"
        );
    }

    #[test]
    fn test_ohe_sign_cr_ti() {
        // Both Cr and Ti should have positive orbital Hall conductivity
        // based on theoretical predictions (Kontani et al.)
        let cr = OrbitalHallMaterial::chromium();
        let ti = OrbitalHallMaterial::titanium();
        assert!(cr.sigma_oh > 0.0, "Cr sigma_OH should be positive");
        assert!(ti.sigma_oh > 0.0, "Ti sigma_OH should be positive");
    }

    #[test]
    fn test_giant_ohe_in_light_metals_vs_pt_she() {
        // Key physics: light metals can have OHE >> Pt SHE despite weaker SOC
        let cr = OrbitalHallMaterial::chromium();
        let pt = OrbitalHallMaterial::platinum();
        let sigma_sh_pt = 2000.0; // Pt spin Hall conductivity

        // Cr OHE should be much larger than Pt SHE
        assert!(
            cr.sigma_oh > sigma_sh_pt * 2.0,
            "Cr OHE ({}) should be >> Pt SHE ({})",
            cr.sigma_oh,
            sigma_sh_pt
        );
        // But Cr has much weaker SOC than Pt
        assert!(
            cr.soc_strength < pt.soc_strength * 0.1,
            "Cr SOC ({}) should be << Pt SOC ({})",
            cr.soc_strength,
            pt.soc_strength
        );
    }

    #[test]
    fn test_orbital_to_spin_conversion_efficiency_bounds() {
        // Conversion efficiency must be bounded by [-1, 1]
        let valid = OrbitalToSpinConverter::new(0.5, "test").expect("valid converter");
        assert!((valid.conversion_efficiency - 0.5).abs() < 1e-10);

        let neg_valid =
            OrbitalToSpinConverter::new(-0.8, "test negative").expect("valid negative converter");
        assert!((neg_valid.conversion_efficiency - (-0.8)).abs() < 1e-10);

        // Out of bounds should fail
        let too_large = OrbitalToSpinConverter::new(1.5, "invalid");
        assert!(too_large.is_err(), "efficiency > 1 should be rejected");
        let too_negative = OrbitalToSpinConverter::new(-1.5, "invalid");
        assert!(too_negative.is_err(), "efficiency < -1 should be rejected");
    }

    #[test]
    fn test_orbital_current_density_direction() {
        let ohe = cr_pt_system();
        // E along x => J_L along -y (E x z = x x z = -y)
        // Actually: x cross z = (0, -1, 0) * (-1) ... let's compute:
        // e_x cross e_z = (0*0 - 0*1, 0*0 - 1*0, 1*1 - 0*0) = (0, 0, 1)
        // Wait: cross product of (1,0,0) x (0,0,1) = (0*1-0*0, 0*0-1*1, 1*0-0*0) = (0,-1,0)
        let e_field = Vector3::new(1.0e5, 0.0, 0.0);
        let j_l = ohe.orbital_current_density(e_field);

        // E x z_hat = (1,0,0) x (0,0,1) = (0,-1,0)
        // So j_l should be along -y
        assert!(
            j_l.y < 0.0,
            "Orbital current should be along -y for E along x, got j_l.y={}",
            j_l.y
        );
        assert!(j_l.x.abs() < 1e-10, "j_l.x should be zero");
        assert!(j_l.z.abs() < 1e-10, "j_l.z should be zero");
    }

    #[test]
    fn test_orbital_current_density_magnitude() {
        let ohe = cr_pt_system();
        let e_field = Vector3::new(1.0e5, 0.0, 0.0);
        let j_l = ohe.orbital_current_density(e_field);
        let mag = j_l.magnitude();

        // sigma_oh = 5000, E = 1e5, factor 100 for unit conversion
        // |J_L| = 5000 * 100 * 1e5 = 5e10
        let expected = 5000.0 * 100.0 * 1.0e5;
        let rel_err = ((mag - expected) / expected).abs();
        assert!(
            rel_err < 1e-10,
            "Orbital current magnitude mismatch: got {}, expected {}",
            mag,
            expected
        );
    }

    #[test]
    fn test_effective_spin_current_conversion() {
        let ohe = cr_pt_system();
        let e_field = Vector3::new(1.0e5, 0.0, 0.0);

        let j_orbital = ohe.orbital_current_density(e_field);
        let j_spin = ohe.effective_spin_current(e_field);

        // J_s = eta_LS * J_L, so |J_s| = |eta_LS| * |J_L|
        let eta = ohe.converter.conversion_efficiency;
        let expected_ratio = eta.abs();
        let actual_ratio = j_spin.magnitude() / j_orbital.magnitude();
        assert!(
            (actual_ratio - expected_ratio).abs() < 1e-10,
            "Spin/orbital ratio mismatch"
        );
    }

    #[test]
    fn test_light_metal_enhancement_factor() {
        // The Cr/Pt system should show enhancement over bare Pt SHE
        let cr_system = cr_pt_system();
        let enhancement = cr_system.enhancement_over_pt_she();

        // sigma_SH_eff = eta_LS(0.5) * sigma_OH(5000) = 2500
        // Enhancement = 2500 / 2000 = 1.25
        assert!(
            enhancement > 1.0,
            "Cr/Pt OHE should enhance over Pt SHE, got {}",
            enhancement
        );
    }

    #[test]
    fn test_material_parameter_validation() {
        // Valid material
        let mat = OrbitalHallMaterial::new("Test", 1000.0, 1e-7, 0.5, 0.01);
        assert!(mat.is_ok(), "Valid material should be accepted");

        // Invalid: zero resistivity
        let bad_rho = OrbitalHallMaterial::new("Bad", 1000.0, 0.0, 0.5, 0.01);
        assert!(bad_rho.is_err(), "Zero resistivity should be rejected");

        // Invalid: negative resistivity
        let neg_rho = OrbitalHallMaterial::new("Bad", 1000.0, -1.0, 0.5, 0.01);
        assert!(neg_rho.is_err(), "Negative resistivity should be rejected");

        // Invalid: negative SOC
        let neg_soc = OrbitalHallMaterial::new("Bad", 1000.0, 1e-7, -0.1, 0.01);
        assert!(neg_soc.is_err(), "Negative SOC should be rejected");
    }

    #[test]
    fn test_material_ohe_ranking() {
        let ranking = material_ohe_ranking();
        assert!(!ranking.is_empty(), "Ranking should not be empty");

        // First entry should have the highest sigma_OH
        assert_eq!(
            ranking[0].0, "Cr",
            "Cr should have the highest OHE conductivity"
        );

        // Values should be in descending order
        for i in 1..ranking.len() {
            assert!(
                ranking[i - 1].1 >= ranking[i].1,
                "Ranking should be descending: {} ({}) >= {} ({})",
                ranking[i - 1].0,
                ranking[i - 1].1,
                ranking[i].0,
                ranking[i].1
            );
        }
    }

    #[test]
    fn test_orbital_accumulation() {
        let ohe = cr_pt_system();
        let e_mag = 1.0e5; // V/m
        let lambda_l = 10.0e-9; // 10 nm orbital diffusion length

        let acc = ohe.orbital_accumulation(e_mag, lambda_l);
        // Accumulation should be positive and finite
        assert!(acc > 0.0, "Orbital accumulation should be positive");
        assert!(acc.is_finite(), "Orbital accumulation should be finite");
    }

    #[test]
    fn test_current_direction_dependence() {
        let ohe = cr_pt_system();

        // E along x => J_L along -y
        let e_x = Vector3::new(1.0e5, 0.0, 0.0);
        let j_x = ohe.orbital_current_density(e_x);
        assert!(j_x.y < 0.0, "E_x should give J_L in -y");

        // E along y => J_L along +x
        // (0,1,0) x (0,0,1) = (1*1-0*0, 0*0-0*1, 0*0-1*0) = (1,0,0)
        let e_y = Vector3::new(0.0, 1.0e5, 0.0);
        let j_y = ohe.orbital_current_density(e_y);
        assert!(j_y.x > 0.0, "E_y should give J_L in +x");

        // E along z => J_L = 0 (z x z = 0)
        let e_z = Vector3::new(0.0, 0.0, 1.0e5);
        let j_z = ohe.orbital_current_density(e_z);
        assert!(
            j_z.magnitude() < 1e-5,
            "E_z should give zero orbital current"
        );
    }

    #[test]
    fn test_sigma_oh_unit_conversion() {
        let sigma_cgs = 5000.0; // (hbar/e)(Ohm*cm)^-1
        let sigma_si = sigma_oh_to_si(sigma_cgs);
        // Should be positive and finite
        assert!(sigma_si > 0.0, "SI sigma should be positive");
        assert!(sigma_si.is_finite(), "SI sigma should be finite");
        // Check approximate magnitude: 5000 * 100 * hbar/(2e) ~ 5e5 * 3.3e-16 ~ 1.6e-10
        let approx = 5000.0 * 100.0 * HBAR / (2.0 * E_CHARGE);
        let rel_err = ((sigma_si - approx) / approx).abs();
        assert!(rel_err < 1e-10, "Unit conversion mismatch");
    }

    #[test]
    fn test_display_formats() {
        let cr = OrbitalHallMaterial::chromium();
        let display = format!("{}", cr);
        assert!(
            display.contains("Cr"),
            "Display should contain material name"
        );

        let conv = OrbitalToSpinConverter::cu_pt();
        let display = format!("{}", conv);
        assert!(
            display.contains("Cu/Pt"),
            "Display should contain interface type"
        );

        let ohe = cr_pt_system();
        let display = format!("{}", ohe);
        assert!(
            display.contains("Cr"),
            "OHE display should contain material"
        );
    }
}
