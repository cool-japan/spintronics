//! Effect module prelude
//!
//! Convenient imports for spin-charge conversion effects.
//!
//! ```rust
//! use spintronics::effect::prelude::*;
//! ```

pub use super::{
    InverseSpinHall, RashbaSystem, SpinNernst, SpinOrbitTorque, SpinSeebeck, TopologicalHall,
};

// Type aliases for convenience (common abbreviations in literature)
/// Inverse Spin Hall Effect (ISHE)
pub type ISHE = InverseSpinHall;
/// Spin-Orbit Torque (SOT)
pub type SOT = SpinOrbitTorque;
/// Spin Nernst Effect (SNE)
pub type SNE = SpinNernst;
/// Spin Seebeck Effect (SSE)
pub type SSE = SpinSeebeck;
/// Topological Hall Effect (THE)
pub type THE = TopologicalHall;
/// Rashba 2D electron gas
pub type Rashba2DEG = RashbaSystem;
