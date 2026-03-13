//! Higher-order time integrators for spin dynamics
//!
//! This module provides a suite of numerical integration methods suitable for
//! magnetization dynamics and general ODE systems expressed as systems of `Vector3`
//! equations. The integrators range from embedded adaptive Runge-Kutta pairs
//! (Dormand-Prince 5(4) and 8(7)) to structure-preserving symplectic methods
//! (Velocity Verlet, Yoshida 4th order, Forest-Ruth) and a semi-implicit solver
//! for stiff problems.
//!
//! # Design
//!
//! All integrators implement the [`Integrator`] trait, which provides a uniform
//! `step` interface returning an [`IntegratorOutput`]. The [`AdaptiveIntegrator`]
//! wrapper adds automatic step size control on top of any integrator that provides
//! error estimates.
//!
//! # Physical Motivation
//!
//! Spin dynamics governed by the Landau-Lifshitz-Gilbert equation exhibit both
//! precessional (energy-conserving) and dissipative behaviour. Symplectic methods
//! excel at long-time energy conservation in the undamped limit, while adaptive
//! Runge-Kutta methods efficiently handle varying timescales that arise from
//! exchange coupling, anisotropy, and applied field pulses.

mod adaptive;
mod dormand_prince;
mod rhs_fn;
mod semi_implicit;
mod symplectic;

#[cfg(test)]
mod tests;

pub use adaptive::AdaptiveIntegrator;
pub use dormand_prince::{DormandPrince45, DormandPrince87};
pub use rhs_fn::{Integrator, IntegratorOutput, RhsFn};
pub use semi_implicit::SemiImplicit;
pub use symplectic::{ForestRuth, VelocityVerlet, Yoshida4};
