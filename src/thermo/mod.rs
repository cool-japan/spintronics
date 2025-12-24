//! Thermoelectric and thermomagnetic effects in spintronics
//!
//! This module implements thermal spin transport phenomena including:
//! - **Spin Seebeck Effect (SSE)**: Temperature gradient → spin current
//! - **Anomalous Nernst Effect (ANE)**: Temperature gradient → transverse voltage
//! - **Spin Peltier Effect**: Spin current → heat current (reciprocal of SSE)
//! - **Thermal magnon transport**: Heat carried by magnons in magnetic insulators
//! - **Multilayer thermal transport**: Heat flow in FM/NM/FM stacks
//!
//! # Quick Start
//!
//! ```rust
//! use spintronics::thermo::prelude::*;
//!
//! let ane = AnomalousNernst::cofeb();
//! let peltier = SpinPeltier::yig_pt(300.0);
//! ```

pub mod ane;
pub mod magnon_thermal;
pub mod multilayer;
pub mod peltier;
pub mod prelude;

pub use ane::AnomalousNernst;
pub use magnon_thermal::{MagnonThermalConductivity, ThermalMagnonTransport};
pub use multilayer::{Layer, MultilayerStack, ThermalBoundary};
pub use peltier::SpinPeltier;
