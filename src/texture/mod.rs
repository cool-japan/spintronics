//! Magnetic textures and topological structures
//!
//! This module implements:
//! - Skyrmions: Topologically protected spin textures
//! - Domain walls: Interfaces between magnetic domains
//! - Vortices: Circular spin configurations
//! - DMI (Dzyaloshinskii-Moriya Interaction): Chiral magnetic interactions
//! - Topological charge calculations

pub mod dmi;
pub mod domain_wall;
pub mod skyrmion;
pub mod topology;

pub use dmi::{DmiParameters, DmiType};
pub use domain_wall::{DomainWall, WallType};
pub use skyrmion::{SkymionLattice, Skyrmion};
pub use topology::{calculate_skyrmion_number, TopologicalCharge};
