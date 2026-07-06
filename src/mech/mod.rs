//! Nanomechanical spintronics
//!
//! Coupling between magnetic moments and mechanical motion at the nanoscale.
//!
//! ## Key Effects
//!
//! 1. **Barnett Effect**: Rotation → Magnetization
//!    - Mechanical rotation creates magnetization
//!    - H_eff = -Ω / γ
//!
//! 2. **Einstein-de Haas Effect**: Magnetization → Rotation
//!    - Change in magnetization causes mechanical rotation
//!    - Conservation of angular momentum
//!
//! 3. **Coupled Dynamics**: LLG + Newton's equations
//!    - Spin and mechanical degrees of freedom exchange angular momentum
//!
//! ## Applications
//! - Nanomechanical resonators as spin detectors
//! - Optomechanical spin control
//! - Angular momentum transfer experiments

pub mod barnett_effect;
pub mod cantilever;
pub mod coupled_dynamics;
pub mod einstein_de_haas;
pub mod magnetoelastic;
pub mod saw;
pub mod strain_driven_dynamics;
pub mod straintronics;

pub use barnett_effect::BarnettMagnetization;
pub use cantilever::{Cantilever, CantileverMode};
pub use coupled_dynamics::SpinMechanicalCoupling;
pub use einstein_de_haas::EinsteinDeHaas;
pub use magnetoelastic::{
    MagnetoelasticMaterial, PiezoelectricSubstrate, StrainTensor, StressTensor,
};
pub use saw::{
    MagnetoelasticMaterial as SawMagnetoelastic, PiezoSubstrate, SawMagnetoacoustics, SawSource,
    SawSpinWaveExcitation,
};
pub use strain_driven_dynamics::{
    PiezoAcStrainDrive, SawStrainDrive, StrainDrive, StrainDrivenLlgDriver, StrainDrivenStepSample,
    StrainDrivenTrajectory, OPTIMAL_COUPLING_ANGLE_RAD,
};
pub use straintronics::{MultiferroicSwitchingResult, StraintronicDevice, VcmaParameters};
