//! Physical constants for spintronics simulations
//!
//! This module defines fundamental physical constants used in spin dynamics
//! and spintronics calculations, including constants specific to magnetic
//! phenomena discovered by Prof. Saitoh's research group.

/// Reduced Planck constant (Dirac constant) [J·s]
pub const HBAR: f64 = 1.054_571_8e-34;

/// Gyromagnetic ratio for electron \[rad/(s·T)\]
/// γ = g * μ_B / ℏ where g ≈ 2 for free electron
pub const GAMMA: f64 = 1.760_859_6e11;

/// Elementary charge \[C\]
pub const E_CHARGE: f64 = 1.602_176_6e-19;

/// Bohr magneton [J/T]
pub const MU_B: f64 = 9.274_009_9e-24;

/// Boltzmann constant [J/K]
pub const KB: f64 = 1.380_649e-23;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::assertions_on_constants)]
    fn test_constants_positive() {
        assert!(HBAR > 0.0);
        assert!(GAMMA > 0.0);
        assert!(E_CHARGE > 0.0);
        assert!(MU_B > 0.0);
        assert!(KB > 0.0);
    }
}
