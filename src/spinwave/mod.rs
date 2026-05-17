//! Spin wave theory module
//!
//! This module provides a comprehensive framework for spin wave (magnon) physics
//! in ferromagnetic thin films and nanostructures. It covers:
//!
//! - **Dispersion relations**: From the simple Kittel formula through full
//!   Kalinikos-Slavin dipole-exchange theory
//! - **Mode classification**: Damon-Eshbach surface modes, backward volume
//!   magnetostatic waves (BVMSW), and forward volume magnetostatic waves (FVMSW)
//! - **Quantized modes**: Standing spin waves in confined geometries including
//!   stripes, disks, and rectangular elements
//!
//! # Physical Background
//!
//! Spin waves are collective excitations of the magnetic order parameter in a
//! ferromagnet. They are the magnetic analog of phonons in a crystal lattice.
//! The quantum of a spin wave is called a magnon.
//!
//! The dispersion relation omega(k) depends on the interplay of:
//! - **Exchange interaction**: Dominant at short wavelengths (large k), gives
//!   quadratic dispersion omega ~ k^2
//! - **Dipolar interaction**: Dominant at long wavelengths (small k), gives
//!   anisotropic dispersion depending on the angle between k and M
//! - **Zeeman energy**: Sets the uniform mode (FMR) frequency via the Kittel formula
//!
//! # Module Organization
//!
//! - [`dispersion`]: Core dispersion relations (Kittel, exchange, Kalinikos-Slavin)
//! - [`modes`]: Magnetostatic mode types (DE, BVMSW, FVMSW)
//! - [`quantization`]: Quantized modes in finite geometries
//! - [`nanodisk`]: Confined modes in magnetic nanodisks (v0.7.0)
//! - [`magnonic_crystal`]: 1D/2D magnonic crystals with plane-wave band structure (v0.7.0)
//! - [`semi_infinite_de`]: Semi-infinite Damon-Eshbach surface modes (v0.7.0)
//!
//! # References
//!
//! - C. Kittel, "On the Theory of Ferromagnetic Resonance Absorption",
//!   Phys. Rev. 73, 155 (1948)
//! - R. W. Damon and J. R. Eshbach, "Magnetostatic modes of a ferromagnet slab",
//!   J. Phys. Chem. Solids 19, 308 (1961)
//! - B. A. Kalinikos and A. N. Slavin, "Theory of dipole-exchange spin wave
//!   spectrum for ferromagnetic films with mixed exchange boundary conditions",
//!   J. Phys. C 19, 7013 (1986)

pub mod bvmsw;
pub mod damon_eshbach;
pub mod dispersion;
pub mod magnonic_crystal;
pub mod modes;
pub mod nanodisk;
pub mod quantization;
pub mod semi_infinite_de;
pub mod surface;

pub use bvmsw::BackwardVolumeMSW;
pub use damon_eshbach::DamonEshbachDetailed;
pub use dispersion::SpinWaveDispersion;
pub use magnonic_crystal::{MagnonicCrystal1D, MagnonicCrystal2D};
pub use modes::{SpinWaveMode, SpinWaveModeCalculator};
pub use nanodisk::NanodiskSpinWaves;
pub use quantization::{NanostructureGeometry, QuantizedModes};
pub use semi_infinite_de::SemiInfiniteDamonEshbach;
pub use surface::SurfaceSpinWave;
