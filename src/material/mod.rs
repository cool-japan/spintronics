//! Material properties for spintronics simulations

pub mod antiferromagnet;
pub mod ferromagnet;
pub mod interface;
pub mod magnetic_2d;
pub mod multilayer;
pub mod temperature;
pub mod topological;
pub mod weyl;

pub use antiferromagnet::{AfmStructure, Antiferromagnet};
pub use ferromagnet::Ferromagnet;
pub use interface::SpinInterface;
pub use magnetic_2d::{Magnetic2D, MagneticOrdering};
pub use multilayer::{MagneticMultilayer, MultilayerType, SpacerLayer};
pub use temperature::ThermalFerromagnet;
pub use topological::{surface_spin_texture, TopologicalClass, TopologicalInsulator};
pub use weyl::{MagneticState, WeylSemimetal, WeylType};
