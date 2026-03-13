//! Commonly used types and functions
//!
//! This module re-exports the most frequently used items from the library
//! for convenient access via `use spintronics::prelude::*;`

// Altermagnet types
pub use crate::altermagnet::{Altermagnet, AltermagnetTransport, AltermagneticSymmetry};
// Builder
pub use crate::builder::{Simulation, SimulationBuilder, SimulationResult, SolverKind};
// Caloritronics
pub use crate::caloritronics::{
    AllCurrents, CaloritronicsResult, HeatCurrentCalculator, OnsagerMatrix,
    SpinCaloritronicsMaterial,
};
// Physical constants
pub use crate::constants::{
    // Electromagnetic
    ALPHA_FS,
    // Derived/Spintronics
    CONDUCTANCE_QUANTUM,
    // Fundamental
    C_LIGHT,
    EPSILON_0,
    E_CHARGE,
    // Particle
    E_OVER_ME,
    FLUX_QUANTUM,
    GAMMA,
    // Magnetic
    G_LANDE,
    HBAR,
    H_PLANCK,
    KB,
    ME,
    MP,
    MU_0,
    MU_B,
    MU_N,
    NA,
    RESISTANCE_QUANTUM,
    SPIN_QUANTUM,
    THERMAL_VOLTAGE_300K,
};
// Dynamics
pub use crate::dynamics::{calc_dm_dt, LlbMaterial, LlbResult, LlbSolver, LlgSolver};
// Effects
pub use crate::effect::{
    InverseSpinHall, RashbaSystem, SpinNernst, SpinOrbitTorque, SpinSeebeck, TopologicalHall,
};
// Core types
pub use crate::error::{Error, Result};
// Frustrated magnets
pub use crate::frustrated::{
    frustration_parameter, kagome_magnon_bands, pauling_entropy, FrustratedLattice,
    KagomeMagnonConfig, LatticeType, SpinIce, SpinIceParams,
};
// Visualization and I/O
pub use crate::io::{OvfData, OvfFormat, OvfReader, OvfWriter};
// Magnon physics (not available on WASM)
#[cfg(all(not(target_arch = "wasm32"), feature = "scirs2"))]
pub use crate::magnon::MultiDomainSystem;
#[cfg(not(target_arch = "wasm32"))]
pub use crate::magnon::{MagnonSolver, SpinChain, SpinPumpingDetector};
// Material types
pub use crate::material::{
    AfmStructure, Antiferromagnet, Ferromagnet, Magnetic2D, MagneticMultilayer, MagneticOrdering,
    MagneticState, MultilayerType, SpacerLayer, SpinInterface, ThermalFerromagnet,
    TopologicalClass, TopologicalInsulator, WeylSemimetal, WeylType,
};
// Material traits (v0.2.0)
pub use crate::material::{
    InterfaceMaterial, MagneticMaterial, SpinChargeConverter, TemperatureDependent,
    TopologicalMaterial,
};
// Memory management (v0.2.0)
pub use crate::memory::{
    get_f64_vec, get_spin_array, put_f64_vec, put_spin_array, HeunWorkspace, Rk4Workspace,
    SpinArrayPool, VectorPool,
};
// Orbitronics
pub use crate::orbitronics::{
    OrbitalHallEffect, OrbitalHallMaterial, OrbitalRashba, OrbitalToSpinConverter, OrbitalTorque,
};
// Spin wave theory
pub use crate::spinwave::{
    NanostructureGeometry, QuantizedModes, SpinWaveDispersion, SpinWaveMode, SpinWaveModeCalculator,
};
// Magnetic textures
pub use crate::texture::{
    calculate_skyrmion_number, Chirality, DmiParameters, DmiType, DomainWall, Helicity, Skyrmion,
    SkyrmionLattice, TopologicalCharge, WallType,
};
// Hopfion dynamics
pub use crate::texture::{HopfionDynamicsConfig, HopfionDynamicsResult, HopfionDynamicsSolver};
// Thermal effects
pub use crate::thermo::{AnomalousNernst, SpinPeltier};
// Transport
pub use crate::transport::{spin_pumping_current, SpinDiffusion};
// Unit validation utilities
pub use crate::units::{
    is_valid_current_density, is_valid_damping, is_valid_dmi_constant, is_valid_energy,
    is_valid_exchange_stiffness, is_valid_gyromagnetic_ratio, is_valid_magnetic_field,
    is_valid_magnetization, is_valid_resistivity, is_valid_spin_diffusion_length,
    is_valid_spin_hall_angle, is_valid_temperature, is_valid_thickness, is_valid_voltage,
};
pub use crate::vector3::Vector3;
pub use crate::visualization::{
    CsvWriter, Hdf5Reader, Hdf5Writer, JsonWriter, SimulationData, VtkWriter,
};
