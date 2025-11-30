//! Thermoelectric and thermomagnetic effects in spintronics
//!
//! This module implements thermal spin transport phenomena including:
//! - Spin Seebeck Effect (SSE): Temperature gradient → spin current
//! - Anomalous Nernst Effect (ANE): Temperature gradient → transverse voltage
//! - Spin Peltier Effect: Spin current → heat current
//! - Thermal magnon transport in magnetic insulators
//! - Multilayer thermal transport

pub mod ane;
pub mod magnon_thermal;
pub mod multilayer;
pub mod peltier;

pub use ane::AnomalousNernst;
pub use magnon_thermal::{MagnonThermalConductivity, ThermalMagnonTransport};
pub use multilayer::{Layer, MultilayerStack, ThermalBoundary};
pub use peltier::SpinPeltier;
