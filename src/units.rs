//! Unit consistency checks for physical quantities
//!
//! This module provides runtime checks to ensure physical quantities
//! have the correct units and magnitudes expected in spintronics calculations.
//!
//! # Example
//! ```
//! use spintronics::units::*;
//!
//! // Check if a value is a valid magnetization (A/m)
//! let ms = 1.4e5; // YIG saturation magnetization
//! assert!(is_valid_magnetization(ms));
//!
//! // Check if a value is a valid damping parameter
//! let alpha = 0.01;
//! assert!(is_valid_damping(alpha));
//! ```

use crate::constants::*;

/// Check if a magnetization value is physically reasonable \[A/m\]
///
/// Valid range: 1e3 to 1e7 A/m
/// - Lower bound: weak ferrimagnets
/// - Upper bound: strong ferromagnets like Fe, Co
pub fn is_valid_magnetization(ms: f64) -> bool {
    (1e3..=1e7).contains(&ms)
}

/// Check if a damping parameter is physically reasonable (dimensionless)
///
/// Valid range: 1e-5 to 1.0
/// - Lower bound: ultra-low damping materials (YIG)
/// - Upper bound: heavily damped or critical systems
pub fn is_valid_damping(alpha: f64) -> bool {
    (1e-5..=1.0).contains(&alpha)
}

/// Check if an exchange stiffness is physically reasonable \[J/m\]
///
/// Valid range: 1e-12 to 1e-10 J/m
/// - Typical ferromagnets: 1-50 pJ/m
pub fn is_valid_exchange_stiffness(a_ex: f64) -> bool {
    (1e-12..=1e-10).contains(&a_ex)
}

/// Check if a spin Hall angle is physically reasonable (dimensionless)
///
/// Valid range: -1.0 to 1.0
/// - Pt: ~0.07-0.15
/// - Ta: ~-0.12 to -0.15
/// - W: ~-0.30 to -0.50
pub fn is_valid_spin_hall_angle(theta_sh: f64) -> bool {
    theta_sh.abs() <= 1.0
}

/// Check if a temperature is physically reasonable \[K\]
///
/// Valid range: 0 to 2000 K
/// - Lower bound: absolute zero (exclusive)
/// - Upper bound: above Curie temperature of most materials
pub fn is_valid_temperature(temp: f64) -> bool {
    temp > 0.0 && (0.0..=2000.0).contains(&temp)
}

/// Check if a magnetic field is physically reasonable \[T\]
///
/// Valid range: 0 to 100 T
/// - Typical lab fields: 0-10 T
/// - Pulsed fields: up to ~100 T
pub fn is_valid_magnetic_field(h: f64) -> bool {
    (0.0..=100.0).contains(&h)
}

/// Check if a thickness is physically reasonable \[m\]
///
/// Valid range: 0.1 nm to 1 mm
pub fn is_valid_thickness(thickness: f64) -> bool {
    (1e-10..=1e-3).contains(&thickness)
}

/// Check if a resistivity is physically reasonable \[Ω·m\]
///
/// Valid range: 1e-9 to 1e-3 Ω·m
/// - Metals: 1e-8 to 1e-6 Ω·m
/// - Semiconductors: up to 1e-3 Ω·m
pub fn is_valid_resistivity(rho: f64) -> bool {
    (1e-9..=1e-3).contains(&rho)
}

/// Check if a spin diffusion length is physically reasonable \[m\]
///
/// Valid range: 0.1 nm to 1 μm
pub fn is_valid_spin_diffusion_length(lambda_sd: f64) -> bool {
    (1e-10..=1e-6).contains(&lambda_sd)
}

/// Check if a DMI constant is physically reasonable \[J/m²\]
///
/// Valid range: 1e-5 to 1e-2 J/m²
/// - Typical values: 0.1-10 mJ/m²
pub fn is_valid_dmi_constant(d: f64) -> bool {
    (1e-5..=1e-2).contains(&d.abs())
}

/// Check if a gyromagnetic ratio is close to the electron value \[rad/(s·T)\]
///
/// Checks if within 10% of the standard electron gyromagnetic ratio
pub fn is_valid_gyromagnetic_ratio(gamma: f64) -> bool {
    let rel_error = ((gamma - GAMMA) / GAMMA).abs();
    rel_error < 0.1
}

/// Validate energy scale is reasonable \[J\]
///
/// Valid range: 1e-25 to 1e-15 J (approximately 1 μeV to 1 eV)
pub fn is_valid_energy(energy: f64) -> bool {
    (1e-25..=1e-15).contains(&energy.abs())
}

/// Check if a current density is physically reasonable \[A/m²\]
///
/// Valid range: 1e6 to 1e12 A/m²
/// - Lower bound: minimum for spin-torque effects
/// - Upper bound: near breakdown limit
pub fn is_valid_current_density(j: f64) -> bool {
    (1e6..=1e12).contains(&j.abs())
}

/// Check if a voltage is physically reasonable \[V\]
///
/// Valid range: 1 nV to 1000 V
pub fn is_valid_voltage(v: f64) -> bool {
    (1e-9..=1e3).contains(&v.abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_magnetization_validation() {
        // YIG
        assert!(is_valid_magnetization(1.4e5));
        // Permalloy
        assert!(is_valid_magnetization(8.6e5));
        // Iron
        assert!(is_valid_magnetization(1.7e6));

        // Too small
        assert!(!is_valid_magnetization(100.0));
        // Too large
        assert!(!is_valid_magnetization(1e8));
        // Negative
        assert!(!is_valid_magnetization(-1e5));
    }

    #[test]
    fn test_damping_validation() {
        // YIG (ultra-low)
        assert!(is_valid_damping(1e-4));
        // Permalloy (typical)
        assert!(is_valid_damping(0.008));
        // Large damping
        assert!(is_valid_damping(0.5));

        // Too small
        assert!(!is_valid_damping(1e-6));
        // Too large
        assert!(!is_valid_damping(1.5));
        // Negative
        assert!(!is_valid_damping(-0.01));
    }

    #[test]
    fn test_spin_hall_angle_validation() {
        // Pt
        assert!(is_valid_spin_hall_angle(0.1));
        // Ta (negative)
        assert!(is_valid_spin_hall_angle(-0.12));
        // W (large negative)
        assert!(is_valid_spin_hall_angle(-0.4));

        // Too large
        assert!(!is_valid_spin_hall_angle(1.5));
        assert!(!is_valid_spin_hall_angle(-1.5));
    }

    #[test]
    fn test_temperature_validation() {
        // Room temperature
        assert!(is_valid_temperature(300.0));
        // Liquid nitrogen
        assert!(is_valid_temperature(77.0));
        // High temperature
        assert!(is_valid_temperature(1500.0));

        // Absolute zero
        assert!(!is_valid_temperature(0.0));
        // Negative
        assert!(!is_valid_temperature(-10.0));
        // Too high
        assert!(!is_valid_temperature(3000.0));
    }

    #[test]
    fn test_gyromagnetic_ratio_validation() {
        // Standard electron value
        assert!(is_valid_gyromagnetic_ratio(GAMMA));
        // Within 5% tolerance
        assert!(is_valid_gyromagnetic_ratio(GAMMA * 1.05));

        // Outside 10% tolerance
        assert!(!is_valid_gyromagnetic_ratio(GAMMA * 1.2));
        assert!(!is_valid_gyromagnetic_ratio(GAMMA * 0.8));
    }

    #[test]
    fn test_thickness_validation() {
        // 1 nm
        assert!(is_valid_thickness(1e-9));
        // 10 nm
        assert!(is_valid_thickness(10e-9));
        // 1 μm
        assert!(is_valid_thickness(1e-6));

        // Too thin (atomic scale)
        assert!(!is_valid_thickness(1e-11));
        // Too thick
        assert!(!is_valid_thickness(1e-2));
    }
}
