//! Commonly used types and functions
//!
//! This module re-exports the most frequently used items from the library
//! for convenient access via `use spintronics::prelude::*;`

// Altermagnet types
pub use crate::altermagnet::{Altermagnet, AltermagnetTransport, AltermagneticSymmetry};
// ML enhancements (v0.7.0)
#[cfg(feature = "autodiff")]
pub use crate::autodiff::{
    find_afm_ground_state, find_fm_ground_state, Activation, EnergyFunctional, Layer, LlgPinn,
    MagneticStructureOptimizer, Mlp, NeuralAnisotropy, NeuralExchange, PinnTrainer, SpinConfig,
    StructureOptResult,
};
// Advanced ML Phase 3 (v0.8.0)
#[cfg(feature = "autodiff")]
pub use crate::autodiff::{
    random_so3, rotate_vector, ActiveLearnResult, ActiveLearner, ActiveLearningConfig,
    EquivariantConfig, EquivariantLinear, EquivariantMlp, QueryStrategy,
};
// Advanced ML Phase 4 (v0.9.0)
#[cfg(feature = "autodiff")]
pub use crate::autodiff::{
    AcquisitionStrategy, BayesianOptConfig, BayesianOptResult, BayesianOptimizer, GaussianProcess,
    GpConfig, GraphMessagePassingLayer, GraphMlp, LatticeGraph, NodeFeatures,
};
// ML autodiff (v0.6.0)
#[cfg(feature = "autodiff")]
pub use crate::autodiff::{Adam, FitResult, LBfgs, OptimizerKind, ParameterFitter, Sgd, Tape, Var};
// Builder
pub use crate::builder::{Simulation, SimulationBuilder, SimulationResult, SolverKind};
// Caloritronics
pub use crate::caloritronics::{
    AllCurrents, CaloritronicsResult, HeatCurrentCalculator, OnsagerMatrix,
    SpinCaloritronicsMaterial,
};
// Cavity extensions (v0.4.0)
pub use crate::cavity::{
    Branch, BrillouinScattering, MagnonPolariton, MagnonicFrequencyComb, MicrowaveToOptical,
    MultiModePolariton, OptomagnonicCoupling, TavisCummings,
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
// Stiff/diffusion integrators (v0.7.0).
// Note: BoundaryCondition is aliased as DiffusionBoundary to avoid clash with negf::BoundaryCondition.
pub use crate::dynamics::integrators::BoundaryCondition as DiffusionBoundary;
pub use crate::dynamics::integrators::{
    CrankNicolsonDiffusion, ImplicitMidpointNewton, SpinDiffusionCrankNicolson,
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
// GPU acceleration (v0.9.0) — Device trait always available; CudaDevice gated.
#[cfg(feature = "cuda")]
pub use crate::gpu::CudaDevice;
pub use crate::gpu::{available_devices, select_best_device, CpuDevice, Device};
// Visualization and I/O
pub use crate::io::{OvfData, OvfFormat, OvfReader, OvfWriter};
// Spin wave extensions (v0.6.0)
#[cfg(all(feature = "scirs2", not(target_arch = "wasm32")))]
pub use crate::magnon::spectral::SpectralMagnonSolver;
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
// Random anisotropy disorder (v0.4.0)
pub use crate::material::{RandomAnisotropy, RandomAnisotropyDistribution};
// Math primitives (v0.4.0)
pub use crate::math::{CMatrix, Complex};
// Memory management (v0.2.0)
pub use crate::memory::{
    get_f64_vec, get_spin_array, put_f64_vec, put_spin_array, HeunWorkspace, Rk4Workspace,
    SpinArrayPool, VectorPool,
};
// Multiferroic / magnetoelectric coupling (v0.5.0)
pub use crate::multiferroic::{
    dzyaloshinskii_moriya_polarization, exchange_striction_polarization, toroidal_moment,
    InverseMagnetoelectric, KnbMechanism, MagnetoelectricTensor, MultiferroicType,
};
// NEGF non-equilibrium transport (v0.4.0)
pub use crate::negf::{
    BoundaryCondition, GreenFunction, Hamiltonian1D, KeldyshSolver, LeadSelfEnergy, SanchoRubio,
    ShotNoise, SpinAccumulation1D, TransportCalculator,
};
// Noncollinear magnetism (v0.5.0)
pub use crate::noncollinear::{
    ExchangeInteraction, LuttingerTisza, SpinSpiral, SpiralChirality, SpiralType,
};
// Orbitronics
pub use crate::orbitronics::{
    OrbitalHallEffect, OrbitalHallMaterial, OrbitalRashba, OrbitalToSpinConverter, OrbitalTorque,
};
// Quantum magnonics (v0.4.0)
pub use crate::quantum::{BogoliubovTransform, HolsteinPrimakoff, HpOrder, ZeroPointFluctuations};
pub use crate::spinwave::{BackwardVolumeMSW, DamonEshbachDetailed, SurfaceSpinWave};
// Advanced spin wave models (v0.7.0)
pub use crate::spinwave::{
    MagnonicCrystal1D, MagnonicCrystal2D, NanodiskSpinWaves, SemiInfiniteDamonEshbach,
};
// Spin wave theory
pub use crate::spinwave::{
    NanostructureGeometry, QuantizedModes, SpinWaveDispersion, SpinWaveMode, SpinWaveModeCalculator,
};
// Stochastic methods (v0.8.0)
#[cfg(feature = "scirs2")]
pub use crate::stochastic::{
    HeunAdaptive, ImplicitMilstein, PimcConfig, PimcLattice, PimcResult, PimcSimulation,
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
// Topological magnon bands (v0.4.0)
// Note: LatticeType from topomagnon is aliased as TopoLatticeType to avoid clash with frustrated::LatticeType
pub use crate::topomagnon::band_model::LatticeType as TopoLatticeType;
// Advanced topology (v0.6.0)
pub use crate::topomagnon::{
    AxionElectrodynamics, AxionMagnonPhoton, BbhModel, BreathingKagomeModel, CornerStateSolver,
    HotiLattice, MagnonBandModel3D, WilsonLoop,
};
pub use crate::topomagnon::{
    BerryCurvature, ChernNumber, EdgeMode, EdgeModes, EdgeSide, KaneMeleModel, MagnonBandModel,
    MagnonHallConductivity,
};
// Transport
pub use crate::transport::{spin_pumping_current, SpinDiffusion};
// Unit validation utilities
pub use crate::units::{
    is_valid_current_density, is_valid_damping, is_valid_dmi_constant, is_valid_energy,
    is_valid_exchange_stiffness, is_valid_gyromagnetic_ratio, is_valid_magnetic_field,
    is_valid_magnetization, is_valid_resistivity, is_valid_spin_diffusion_length,
    is_valid_spin_hall_angle, is_valid_temperature, is_valid_thickness, is_valid_voltage,
};
// More experimental validations (v0.9.0)
pub use crate::validation::experimental::boona_2014::Boona2014Validation;
// Experimental validation (v0.7.0)
pub use crate::validation::experimental::demidov_2006::Demidov2006Validation;
pub use crate::validation::experimental::garello_2013::Garello2013Validation;
// More experimental validations (v0.8.0)
pub use crate::validation::experimental::liu_2012::Liu2012Validation;
pub use crate::validation::experimental::mosendz_2010::Mosendz2010Validation;
pub use crate::validation::experimental::saitoh_2006::Saitoh2006Validation;
pub use crate::validation::experimental::uchida_2008::Uchida2008Validation;
pub use crate::validation::experimental::ValidationResult as ExperimentalValidationResult;
pub use crate::vector3::Vector3;
// Data export formats (v0.6.0)
#[cfg(feature = "netcdf")]
pub use crate::visualization::netcdf::{NetCdfReader, NetCdfWriter};
#[cfg(feature = "vti")]
pub use crate::visualization::vti::VtiWriter;
#[cfg(feature = "xdmf")]
pub use crate::visualization::xdmf::{XdmfTimeStep, XdmfWriter};
#[cfg(feature = "zarr")]
pub use crate::visualization::zarr::{ZarrDtype, ZarrStore};
pub use crate::visualization::{
    CsvWriter, Hdf5Reader, Hdf5Writer, JsonWriter, SimulationData, VtkWriter,
};
