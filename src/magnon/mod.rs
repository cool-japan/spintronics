//! Magnon propagation and spin wave dynamics
//!
//! This module simulates magnon (spin wave) propagation in magnetic materials,
//! including the excitation, propagation, and detection of spin waves via
//! spin pumping and inverse spin Hall effect.

pub mod chain;
pub mod detector;
#[cfg(feature = "scirs2")]
pub mod parallel;
pub mod solver;

pub use chain::SpinChain;
pub use detector::SpinPumpingDetector;
#[cfg(feature = "scirs2")]
pub use parallel::MultiDomainSystem;
pub use solver::MagnonSolver;
