//! Fluid spintronics and spin-vorticity coupling
//!
//! This module implements the groundbreaking discovery by Prof. Saitoh's group:
//! **spin current generation in liquid metals via vorticity**.
//!
//! Key papers:
//! - R. Takahashi et al., "Spin hydrodynamic generation",
//!   Nature Physics 12, 52-56 (2016)
//! - D. Kobayashi et al., "Spin current generation using a surface acoustic wave
//!   generated via spin-rotation coupling", Phys. Rev. Lett. 119, 077202 (2017)
//!
//! ## Physical Principle
//!
//! Fluid vorticity (rotation) couples to electron spin via the Barnett effect:
//! - Vorticity ω = ∇ × v creates an effective magnetic field
//! - H_Barnett = ω / γ
//! - This polarizes spins in rotating fluids (mercury, galinstan, etc.)

pub mod barnett;
pub mod navier_stokes;
pub mod vorticity;

pub use barnett::BarnettField;
pub use navier_stokes::{FluidField, NavierStokes};
pub use vorticity::VorticityCalculator;
