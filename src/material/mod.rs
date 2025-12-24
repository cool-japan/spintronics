//! Material properties for spintronics simulations
//!
//! This module provides material parameters for various spintronic systems:
//!
//! - **Ferromagnets**: YIG, Permalloy, CoFeB, CoFe, Fe, Co, Ni
//! - **Antiferromagnets**: NiO, MnF₂, FeF₂, Cr₂O₃
//! - **Topological Insulators**: Bi₂Se₃, Bi₂Te₃, Sb₂Te₃
//! - **Weyl Semimetals**: TaAs, NbAs, MoTe₂
//! - **2D Magnets**: CrI₃, Fe₃GeTe₂, MnBi₂Te₄
//! - **Interfaces**: Spin mixing conductance for various FM/NM interfaces
//!
//! # Quick Start
//!
//! ```rust
//! use spintronics::material::prelude::*;
//!
//! let yig = Ferromagnet::yig();
//! let interface = SpinInterface::yig_pt();
//! let ti = TopologicalInsulator::bi2se3();
//! ```

pub mod antiferromagnet;
pub mod ferromagnet;
pub mod interface;
pub mod magnetic_2d;
pub mod multilayer;
pub mod prelude;
pub mod temperature;
pub mod topological;
pub mod traits;
pub mod weyl;

pub use antiferromagnet::{AfmStructure, Antiferromagnet};
pub use ferromagnet::Ferromagnet;
pub use interface::SpinInterface;
pub use magnetic_2d::{Magnetic2D, MagneticOrdering};
pub use multilayer::{MagneticMultilayer, MultilayerType, SpacerLayer};
pub use temperature::ThermalFerromagnet;
pub use topological::{surface_spin_texture, TopologicalClass, TopologicalInsulator};
pub use traits::{
    InterfaceMaterial, MagneticMaterial, SpinChargeConverter, TemperatureDependent,
    TopologicalMaterial,
};
pub use weyl::{MagneticState, WeylSemimetal, WeylType};
