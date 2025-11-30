//! Commonly used types and functions
//!
//! This module re-exports the most frequently used items from the library
//! for convenient access via `use spintronics::prelude::*;`

pub use crate::constants::{E_CHARGE, GAMMA, HBAR, KB, MU_B};

pub use crate::vector3::Vector3;

pub use crate::material::{
    AfmStructure, Antiferromagnet, Ferromagnet, Magnetic2D, MagneticMultilayer, MagneticOrdering,
    MagneticState, MultilayerType, SpacerLayer, SpinInterface, ThermalFerromagnet,
    TopologicalClass, TopologicalInsulator, WeylSemimetal, WeylType,
};

pub use crate::dynamics::{calc_dm_dt, LlgSolver};

pub use crate::transport::{spin_pumping_current, SpinDiffusion};

pub use crate::effect::{
    InverseSpinHall, RashbaSystem, SpinNernst, SpinOrbitTorque, SpinSeebeck, TopologicalHall,
};

#[cfg(not(target_arch = "wasm32"))]
pub use crate::magnon::{MagnonSolver, MultiDomainSystem, SpinChain, SpinPumpingDetector};

pub use crate::texture::{DmiParameters, DmiType};

pub use crate::visualization::{CsvWriter, JsonWriter, SimulationData, VtkWriter};

pub use crate::io::{OvfData, OvfFormat, OvfReader, OvfWriter};

pub use crate::error::{Error, Result};
