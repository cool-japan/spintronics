//! Magnon propagation and spin wave dynamics
//!
//! This module simulates magnon (spin wave) propagation in magnetic materials,
//! including the excitation, propagation, and detection of spin waves via
//! spin pumping and inverse spin Hall effect.

pub mod bec;
pub mod chain;
pub mod detector;
pub mod nonlinear;
#[cfg(feature = "scirs2")]
pub mod parallel;
pub mod solver;
#[cfg(feature = "scirs2")]
pub mod spectral;
pub mod spin_density_wave;

pub use bec::{
    bec_temperature, classical_distribution, critical_density, magnon_distribution,
    magnon_supercurrent, supercurrent_profile, GrossPitaevskiiSolver, MagnonCondensate,
    ParametricPumping,
};
pub use chain::SpinChain;
pub use detector::SpinPumpingDetector;
pub use nonlinear::{
    four_magnon_relaxation_rate, magnon_magnon_interaction_energy,
    suhl_spin_wave_instability_power, FourMagnonScattering, NonlinearFmrLinewidth,
    ParametricAmplification,
};
#[cfg(feature = "scirs2")]
pub use parallel::MultiDomainSystem;
pub use solver::MagnonSolver;
#[cfg(feature = "scirs2")]
pub use spectral::SpectralMagnonSolver;
pub use spin_density_wave::{
    condensation_energy, elastic_energy, enhanced_susceptibility, resistivity_anomaly,
    total_sdw_energy, ChromiumSdw, SdwGapSolver, SdwPolarization, SdwRelaxationDynamics,
    SpinDensityWave,
};
