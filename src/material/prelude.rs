//! Material module prelude
//!
//! Convenient imports for magnetic and topological materials.
//!
//! ```rust
//! use spintronics::material::prelude::*;
//! ```

// Core material types
// Antiferromagnetic materials
pub use super::antiferromagnet::{AfmStructure, Antiferromagnet};
// Topological materials
pub use super::topological::{surface_spin_texture, TopologicalClass, TopologicalInsulator};
// Material traits for generic programming
pub use super::traits::{
    InterfaceMaterial, MagneticMaterial, SpinChargeConverter, TemperatureDependent,
    TopologicalMaterial,
};
pub use super::weyl::{MagneticState, WeylSemimetal, WeylType};
pub use super::{
    Ferromagnet, Magnetic2D, MagneticMultilayer, MagneticOrdering, SpacerLayer, SpinInterface,
    ThermalFerromagnet,
};

// Type aliases for common materials
/// Ferromagnet (FM)
pub type FM = Ferromagnet;
/// Topological Insulator (TI)
pub type TI = TopologicalInsulator;
/// Weyl Semimetal (WSM)
pub type WSM = WeylSemimetal;
/// Antiferromagnet (AFM)
pub type AFM = Antiferromagnet;
