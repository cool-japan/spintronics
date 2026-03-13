//! Magnon propagation and spin wave dynamics
//!
//! This module simulates magnon (spin wave) propagation in magnetic materials,
//! including the excitation, propagation, and detection of spin waves via
//! spin pumping and inverse spin Hall effect.

pub mod bec;
pub mod chain;
pub mod detector;
#[cfg(feature = "scirs2")]
pub mod parallel;
pub mod solver;
pub mod spin_density_wave;

pub use bec::{
    bec_temperature, classical_distribution, critical_density, magnon_distribution,
    magnon_supercurrent, supercurrent_profile, GrossPitaevskiiSolver, MagnonCondensate,
    ParametricPumping,
};
pub use chain::SpinChain;
pub use detector::SpinPumpingDetector;
#[cfg(feature = "scirs2")]
pub use parallel::MultiDomainSystem;
pub use solver::MagnonSolver;
pub use spin_density_wave::{
    condensation_energy, elastic_energy, enhanced_susceptibility, resistivity_anomaly,
    total_sdw_energy, ChromiumSdw, SdwGapSolver, SdwPolarization, SpinDensityWave,
};
