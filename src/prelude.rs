//! Commonly used types and functions
//!
//! This module re-exports the most frequently used items from the library
//! for convenient access via `use spintronics::prelude::*;`

// Altermagnet types
pub use crate::altermagnet::{Altermagnet, AltermagnetTransport, AltermagneticSymmetry};
// Math primitives (v0.4.0)
pub use crate::math::{CMatrix, Complex};
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
pub use crate::magnon::{
    four_magnon_relaxation_rate, magnon_magnon_interaction_energy,
    suhl_spin_wave_instability_power, FourMagnonScattering, MagnonSolver, NonlinearFmrLinewidth,
    ParametricAmplification, SpinChain, SpinPumpingDetector,
};
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
// Quantum magnonics (v0.4.0)
pub use crate::quantum::{BogoliubovTransform, HolsteinPrimakoff, HpOrder, ZeroPointFluctuations};
// NEGF non-equilibrium transport (v0.4.0)
pub use crate::negf::{
    BoundaryCondition, GreenFunction, Hamiltonian1D, KeldyshSolver, LeadSelfEnergy, SanchoRubio,
    ShotNoise, SpinAccumulation1D, TransportCalculator,
};
// Topological magnon bands (v0.4.0)
// Note: LatticeType from topomagnon is aliased as TopoLatticeType to avoid clash with frustrated::LatticeType
pub use crate::topomagnon::band_model::LatticeType as TopoLatticeType;
pub use crate::topomagnon::{
    BerryCurvature, ChernNumber, EdgeMode, EdgeModes, EdgeSide, KaneMeleModel, MagnonBandModel,
    MagnonHallConductivity,
};
// Advanced topology (v0.6.0)
pub use crate::topomagnon::{
    AxionElectrodynamics, AxionMagnonPhoton, BbhModel, BreathingKagomeModel, CornerStateSolver,
    HotiLattice, MagnonBandModel3D, WilsonLoop,
};
// Cavity extensions (v0.4.0)
pub use crate::cavity::{
    Branch, BrillouinScattering, MagnonPolariton, MagnonicFrequencyComb, MicrowaveToOptical,
    MultiModePolariton, OptomagnonicCoupling, TavisCummings,
};
// Random anisotropy disorder (v0.4.0)
pub use crate::material::{RandomAnisotropy, RandomAnisotropyDistribution};
// Spin wave theory
pub use crate::spinwave::{
    NanostructureGeometry, QuantizedModes, SpinWaveDispersion, SpinWaveMode, SpinWaveModeCalculator,
};
// Spin wave extensions (v0.6.0)
#[cfg(all(feature = "scirs2", not(target_arch = "wasm32")))]
pub use crate::magnon::spectral::SpectralMagnonSolver;
pub use crate::spinwave::{BackwardVolumeMSW, DamonEshbachDetailed, SurfaceSpinWave};
// Magnetic textures
pub use crate::texture::{
    calculate_skyrmion_number, Chirality, DmiParameters, DmiType, DomainWall, Helicity, Skyrmion,
    SkyrmionLattice, TopologicalCharge, WallType,
};
// Hopfion dynamics
pub use crate::texture::{HopfionDynamicsConfig, HopfionDynamicsResult, HopfionDynamicsSolver};
// Noncollinear magnetism (v0.5.0)
pub use crate::noncollinear::{
    ExchangeInteraction, LuttingerTisza, SpinSpiral, SpiralChirality, SpiralType,
};
// Multiferroic / magnetoelectric coupling (v0.5.0)
pub use crate::multiferroic::{
    dzyaloshinskii_moriya_polarization, exchange_striction_polarization, toroidal_moment,
    InverseMagnetoelectric, KnbMechanism, MagnetoelectricTensor, MultiferroicType,
};
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
// ML autodiff (v0.6.0)
#[cfg(feature = "autodiff")]
pub use crate::autodiff::{Adam, FitResult, LBfgs, OptimizerKind, ParameterFitter, Sgd, Tape, Var};
pub use crate::vector3::Vector3;
pub use crate::visualization::{
    CsvWriter, Hdf5Reader, Hdf5Writer, JsonWriter, SimulationData, VtkWriter,
};
// Data export formats (v0.6.0)
#[cfg(feature = "netcdf")]
pub use crate::visualization::netcdf::{NetCdfReader, NetCdfWriter};
#[cfg(feature = "vti")]
pub use crate::visualization::vti::VtiWriter;
#[cfg(feature = "xdmf")]
pub use crate::visualization::xdmf::{XdmfTimeStep, XdmfWriter};
#[cfg(feature = "zarr")]
pub use crate::visualization::zarr::{ZarrDtype, ZarrStore};
