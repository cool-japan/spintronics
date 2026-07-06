//! Magnetic textures and topological structures
//!
//! This module implements:
//! - **Skyrmions**: Topologically protected spin textures (Néel and Bloch types)
//! - **Hopfions**: 3D topological solitons characterized by the Hopf invariant
//! - **Domain walls**: Interfaces between magnetic domains
//! - **Vortices**: Circular spin configurations
//! - **DMI**: Dzyaloshinskii-Moriya Interaction (interface and bulk)
//! - **Topological charge**: Skyrmion number calculations
//!
//! # Quick Start
//!
//! ```rust
//! use spintronics::texture::prelude::*;
//!
//! let dmi = DmiParameters::pt_co();
//! let skyrmion = Skyrmion::new(
//!     (0.0, 0.0),
//!     50.0e-9,
//!     Helicity::Neel,
//!     Chirality::CounterClockwise,
//! );
//! ```

pub mod dmi;
pub mod domain_wall;
pub mod dw_dynamics;
pub mod hopfion;
pub mod hopfion_dynamics;
pub mod hopfion_stability_modes;
pub mod prelude;
pub mod skyrmion;
pub mod topology;

pub use dmi::{DmiParameters, DmiType};
pub use domain_wall::{DomainWall, WallType};
pub use dw_dynamics::{DwMaterial, DwSotDynamics, DwSttDynamics, WalkerBreakdown};
pub use hopfion::{
    HopfFibration, HopfInvariant, Hopfion, HopfionEnergy, HopfionEnergyParams, HopfionStability,
};
pub use hopfion_dynamics::{HopfionDynamicsConfig, HopfionDynamicsResult, HopfionDynamicsSolver};
pub use hopfion_stability_modes::{
    CollectiveMode, CollectiveStepSizes, EigenmodeSpectrum, HopfionEigenmodeStability,
};
pub use skyrmion::{Chirality, Helicity, Skyrmion, SkyrmionLattice};
pub use topology::{calculate_skyrmion_number, TopologicalCharge};
