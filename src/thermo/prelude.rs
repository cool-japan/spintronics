//! Thermo module prelude
//!
//! Convenient imports for thermoelectric and thermomagnetic effects.
//!
//! ```rust
//! use spintronics::thermo::prelude::*;
//! ```

// Core thermal effects
pub use super::ane::AnomalousNernst;
// Magnon thermal transport
pub use super::magnon_thermal::{MagnonThermalConductivity, ThermalMagnonTransport};
// Thermal multilayers
pub use super::multilayer::{Layer, MultilayerStack, ThermalBoundary};
pub use super::peltier::SpinPeltier;

// Type aliases for convenience
/// Anomalous Nernst Effect (ANE)
pub type ANE = AnomalousNernst;
/// Spin Peltier Effect (SPE)
pub type SPE = SpinPeltier;
