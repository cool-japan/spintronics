//! Type-state markers for compile-time builder validation
//!
//! These zero-sized types serve as compile-time markers in the
//! [`SimulationBuilder`](super::SimulationBuilder) type-state pattern.
//! By encoding the presence or absence of required configuration in the type
//! system, we guarantee at compile time that a simulation cannot be built
//! without all mandatory parameters (material, external field, solver).
//!
//! # Example
//!
//! ```rust
//! use spintronics::builder::{SimulationBuilder, states::*};
//!
//! // Initially, the builder is in the unconfigured state:
//! //   SimulationBuilder<NoMaterial, NoField, NoSolver>
//! let builder = SimulationBuilder::new();
//!
//! // After setting material, the M type parameter becomes WithMaterial.
//! // The .build() method is only available when all three type parameters
//! // are in their "With*" state.
//! ```

/// Marker: no material has been set yet.
///
/// The builder starts in this state for the material type parameter.
/// Call [`.material()`](super::SimulationBuilder::material) to transition
/// to [`WithMaterial`].
pub struct NoMaterial;

/// Marker: a [`Ferromagnet`](crate::material::Ferromagnet) material has been configured.
pub struct WithMaterial;

/// Marker: no external field has been set yet.
///
/// The builder starts in this state for the field type parameter.
/// Call [`.external_field()`](super::SimulationBuilder::external_field) to transition
/// to [`WithField`].
pub struct NoField;

/// Marker: an external magnetic field has been configured.
pub struct WithField;

/// Marker: no solver has been set yet.
///
/// The builder starts in this state for the solver type parameter.
/// Call [`.solver_rk4()`](super::SimulationBuilder::solver_rk4) to transition
/// to [`WithSolver`].
pub struct NoSolver;

/// Marker: a solver (e.g., RK4) has been configured.
pub struct WithSolver;
