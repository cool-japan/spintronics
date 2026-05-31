//! Spin-charge conversion effects
//!
//! This module implements various spin-charge conversion phenomena:
//!
//! - **Exchange Bias**: FM/AFM interface coupling → hysteresis loop shift
//! - **Inverse Spin Hall Effect (ISHE)**: Spin current → charge current
//! - **Rashba Effect**: 2DEG spin splitting and Edelstein effect
//! - **Spin-Orbit Torque (SOT)**: Current-driven magnetization switching
//! - **Spin Nernst Effect (SNE)**: Thermal gradient → transverse spin current
//! - **Spin Seebeck Effect (SSE)**: Thermal generation of spin current
//! - **Topological Hall Effect (THE)**: Berry phase from skyrmion textures
//!
//! # Quick Start
//!
//! ```rust
//! use spintronics::effect::prelude::*;
//!
//! let ishe = InverseSpinHall::platinum();
//! let sot = SpinOrbitTorque::platinum_cofeb();
//! ```

pub mod exchange_bias;
pub mod ishe;
pub mod optical_switching;
pub mod prelude;
pub mod rashba;
pub mod smr;
pub mod sot;
pub mod spin_nernst;
pub mod sse;
pub mod stno;
pub mod topological_hall;

pub use exchange_bias::{
    ExchangeBias, LoopShiftResult, J_EB_BFCO_COFE, J_EB_COO_CO, J_EB_FEMN_NIFE, J_EB_IRMN_CO,
};
pub use ishe::InverseSpinHall;
pub use optical_switching::{
    CircularHelicity, LaserPulseParams, OpticalMagneticMaterial, OpticalSwitchResult,
    OpticalSwitching,
};
pub use rashba::RashbaSystem;
pub use smr::{SpinHallMagnetoresistance, UnidirectionalSmr};
pub use sot::SpinOrbitTorque;
pub use spin_nernst::SpinNernst;
pub use sse::SpinSeebeck;
pub use stno::{SpinTorqueOscillator, SpinTorqueOscillatorConfig};
pub use topological_hall::TopologicalHall;
