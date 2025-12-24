//! Spin-charge conversion effects
//!
//! This module implements various spin-charge conversion phenomena:
//!
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

pub mod ishe;
pub mod prelude;
pub mod rashba;
pub mod sot;
pub mod spin_nernst;
pub mod sse;
pub mod topological_hall;

pub use ishe::InverseSpinHall;
pub use rashba::RashbaSystem;
pub use sot::SpinOrbitTorque;
pub use spin_nernst::SpinNernst;
pub use sse::SpinSeebeck;
pub use topological_hall::TopologicalHall;
