//! Physical constants for spintronics simulations
//!
//! This module defines fundamental physical constants used in spin dynamics
//! and spintronics calculations, including constants specific to magnetic
//! phenomena discovered by Prof. Saitoh's research group.
//!
//! # Constants Overview
//!
//! ## Fundamental Constants
//! - `HBAR`: Reduced Planck constant
//! - `E_CHARGE`: Elementary charge
//! - `KB`: Boltzmann constant
//! - `C_LIGHT`: Speed of light in vacuum
//! - `NA`: Avogadro constant
//!
//! ## Electromagnetic Constants
//! - `MU_0`: Vacuum permeability
//! - `EPSILON_0`: Vacuum permittivity
//! - `ALPHA_FS`: Fine structure constant
//!
//! ## Magnetic Constants
//! - `GAMMA`: Gyromagnetic ratio for electron
//! - `MU_B`: Bohr magneton
//! - `G_LANDE`: Landé g-factor for free electron
//!
//! ## Particle Constants
//! - `ME`: Electron mass
//! - `MP`: Proton mass

// =============================================================================
// Fundamental Constants
// =============================================================================

/// Reduced Planck constant (Dirac constant) [J·s]
///
/// ℏ = h / (2π) where h is Planck's constant
pub const HBAR: f64 = 1.054_571_817e-34;

/// Planck constant [J·s]
///
/// The quantum of action
pub const H_PLANCK: f64 = 6.626_070_15e-34;

/// Elementary charge \[C\]
///
/// The electric charge carried by a single proton
pub const E_CHARGE: f64 = 1.602_176_634e-19;

/// Boltzmann constant [J/K]
///
/// Relates temperature to energy
pub const KB: f64 = 1.380_649e-23;

/// Speed of light in vacuum [m/s]
///
/// Exact value by definition since 2019
pub const C_LIGHT: f64 = 299_792_458.0;

/// Avogadro constant [mol⁻¹]
///
/// Number of particles in one mole
pub const NA: f64 = 6.022_140_76e23;

// =============================================================================
// Electromagnetic Constants
// =============================================================================

/// Vacuum permeability [H/m] or [N/A²]
///
/// μ₀ = 4π × 10⁻⁷ (exact before 2019, now derived)
pub const MU_0: f64 = 1.256_637_062_12e-6;

/// Vacuum permittivity [F/m]
///
/// ε₀ = 1 / (μ₀ c²)
pub const EPSILON_0: f64 = 8.854_187_812_8e-12;

/// Fine structure constant (dimensionless)
///
/// α = e² / (4π ε₀ ℏ c) ≈ 1/137
pub const ALPHA_FS: f64 = 7.297_352_569_3e-3;

// =============================================================================
// Magnetic Constants
// =============================================================================

/// Gyromagnetic ratio for electron \[rad/(s·T)\]
///
/// γ = g μ_B / ℏ where g ≈ 2.002319 for free electron
/// This is the magnitude |γ|; note that γ < 0 for electrons
pub const GAMMA: f64 = 1.760_859_630_23e11;

/// Bohr magneton [J/T]
///
/// μ_B = eℏ / (2 m_e)
pub const MU_B: f64 = 9.274_010_078_3e-24;

/// Landé g-factor for free electron (dimensionless)
///
/// Deviation from 2 is due to QED corrections
pub const G_LANDE: f64 = 2.002_319_304_362_56;

/// Nuclear magneton [J/T]
///
/// μ_N = eℏ / (2 m_p)
pub const MU_N: f64 = 5.050_783_746_1e-27;

// =============================================================================
// Particle Constants
// =============================================================================

/// Electron mass \[kg\]
pub const ME: f64 = 9.109_383_701_5e-31;

/// Proton mass \[kg\]
pub const MP: f64 = 1.672_621_923_7e-27;

/// Electron charge-to-mass ratio [C/kg]
///
/// e / m_e
pub const E_OVER_ME: f64 = 1.758_820_010_8e11;

// =============================================================================
// Derived/Useful Constants for Spintronics
// =============================================================================

/// Spin quantum ℏ/2 [J·s]
///
/// Magnitude of electron spin angular momentum along any axis
pub const SPIN_QUANTUM: f64 = HBAR / 2.0;

/// Thermal voltage at 300 K \[V\]
///
/// V_T = k_B T / e at room temperature
pub const THERMAL_VOLTAGE_300K: f64 = 0.025_85;

/// Magnetic flux quantum \[Wb\]
///
/// Φ₀ = h / (2e) - fundamental unit of magnetic flux in superconductors
pub const FLUX_QUANTUM: f64 = 2.067_833_848e-15;

/// Conductance quantum \[S\]
///
/// G₀ = 2e² / h - fundamental unit of conductance
pub const CONDUCTANCE_QUANTUM: f64 = 7.748_091_729e-5;

/// Resistance quantum \[Ω\]
///
/// R_K = h / e² - von Klitzing constant
pub const RESISTANCE_QUANTUM: f64 = 25_812.807_45;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::assertions_on_constants)]
    fn test_constants_positive() {
        assert!(HBAR > 0.0);
        assert!(H_PLANCK > 0.0);
        assert!(GAMMA > 0.0);
        assert!(E_CHARGE > 0.0);
        assert!(MU_B > 0.0);
        assert!(KB > 0.0);
        assert!(C_LIGHT > 0.0);
        assert!(MU_0 > 0.0);
        assert!(EPSILON_0 > 0.0);
        assert!(ME > 0.0);
        assert!(MP > 0.0);
    }

    #[test]
    fn test_derived_constants() {
        // Check c² = 1 / (μ₀ ε₀)
        let c_squared = 1.0 / (MU_0 * EPSILON_0);
        let relative_error = ((c_squared - C_LIGHT * C_LIGHT) / (C_LIGHT * C_LIGHT)).abs();
        assert!(relative_error < 1e-8);

        // Check μ_B = eℏ / (2 m_e)
        let mu_b_calc = E_CHARGE * HBAR / (2.0 * ME);
        let relative_error = ((mu_b_calc - MU_B) / MU_B).abs();
        assert!(relative_error < 1e-8);

        // Check γ = g μ_B / ℏ
        let gamma_calc = G_LANDE * MU_B / HBAR;
        let relative_error = ((gamma_calc - GAMMA) / GAMMA).abs();
        assert!(relative_error < 1e-8);
    }

    #[test]
    fn test_fine_structure() {
        // α ≈ 1/137
        let alpha_inverse = 1.0 / ALPHA_FS;
        assert!((alpha_inverse - 137.036).abs() < 0.001);
    }
}
