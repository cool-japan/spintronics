//! Texture module prelude
//!
//! Convenient imports for magnetic textures and topological structures.
//!
//! ```rust
//! use spintronics::texture::prelude::*;
//! ```

// Skyrmions and lattices
// DMI (Dzyaloshinskii-Moriya Interaction)
pub use super::dmi::{DmiParameters, DmiType};
// Domain walls
pub use super::domain_wall::{DomainWall, WallType};
pub use super::skyrmion::{Chirality, Helicity, LatticeType, Skyrmion, SkyrmionLattice};
// Topological charge calculations
pub use super::topology::{calculate_skyrmion_number, TopologicalCharge};

// Type aliases for convenience
/// Dzyaloshinskii-Moriya Interaction (DMI)
pub type DMI = DmiParameters;
/// Domain Wall (DW)
pub type DW = DomainWall;
/// Skyrmion (Sk)
pub type Sk = Skyrmion;
