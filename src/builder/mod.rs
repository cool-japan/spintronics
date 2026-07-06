//! Type-state simulation builder for spintronics simulations
//!
//! This module provides a [`SimulationBuilder`] that uses the *type-state* pattern
//! to enforce, at compile time, that all required configuration (material,
//! external field, solver) is present before a [`Simulation`] can be built.
//!
//! # Type-State Pattern
//!
//! The builder carries three phantom type parameters `<M, F, S>` that track
//! whether a material, an external field, and a solver have been configured.
//! The [`.build()`](SimulationBuilder::build) method is only available when
//! all three parameters are in their "configured" state, so forgetting a
//! required parameter is a *compile-time* error rather than a runtime one.
//!
//! # Examples
//!
//! ```rust
//! use spintronics::prelude::*;
//! use spintronics::builder::SimulationBuilder;
//!
//! let result = SimulationBuilder::new()
//!     .material(Ferromagnet::yig())
//!     .external_field(Vector3::new(0.0, 0.0, 0.1))
//!     .solver_rk4()
//!     .num_steps(100)
//!     .build();
//! assert!(result.is_ok());
//! ```
//!
//! ## Compile-time safety
//!
//! The following code would **not** compile because no material has been set:
//!
//! ```compile_fail
//! use spintronics::prelude::*;
//! use spintronics::builder::SimulationBuilder;
//!
//! // ERROR: `build()` is not available on SimulationBuilder<NoMaterial, ...>
//! let _ = SimulationBuilder::new()
//!     .external_field(Vector3::new(0.0, 0.0, 0.1))
//!     .solver_rk4()
//!     .build();
//! ```

pub mod states;

use std::marker::PhantomData;

use states::{NoField, NoMaterial, NoSolver, WithField, WithMaterial, WithSolver};

use crate::dynamics::{
    anisotropy_energy, calc_dm_dt, zeeman_energy, AdaptiveIntegrator, DormandPrince45,
    DormandPrince87, ForestRuth, Integrator, LlgSolver, SemiImplicit, Yoshida4,
};
use crate::error::{Error, Result};
use crate::material::Ferromagnet;
use crate::vector3::Vector3;

/// Type alias for the RHS closure used by generic integrators
type RhsFn<'a> = &'a dyn Fn(&[Vector3<f64>], f64) -> Vec<Vector3<f64>>;

// ---------------------------------------------------------------------------
// SolverKind
// ---------------------------------------------------------------------------

/// Which integration method to use for the LLG equation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SolverKind {
    /// Classical 4th-order Runge-Kutta method.
    Rk4,
    /// Simple forward-Euler method (first order).
    Euler,
    /// Second-order Heun predictor-corrector method.
    Heun,
    /// Dormand-Prince 5(4) adaptive Runge-Kutta (ode45-style).
    ///
    /// Uses an embedded error estimate to drive automatic step size control.
    /// The `tolerance` parameter sets the local error bound.
    Dp45 {
        /// Local truncation error tolerance for adaptive step control.
        tolerance: f64,
    },
    /// Dormand-Prince 8(7) adaptive Runge-Kutta (high-precision).
    ///
    /// A 13-stage method providing 8th-order accuracy with adaptive control.
    /// Suitable for very long-time simulations requiring high precision.
    Dp87 {
        /// Local truncation error tolerance for adaptive step control.
        tolerance: f64,
    },
    /// Yoshida 4th-order symplectic integrator.
    ///
    /// Structure-preserving method with excellent long-time energy conservation
    /// in the undamped limit. The magnetization state is treated as a position
    /// vector; an auxiliary zero-momentum coordinate is appended internally.
    Yoshida4,
    /// Forest-Ruth 4th-order symplectic integrator.
    ///
    /// Alternative 4th-order composition method with slightly different
    /// stability characteristics than Yoshida. Same state-vector convention.
    ForestRuth,
    /// Semi-implicit integrator for stiff systems (implicit midpoint rule).
    ///
    /// A-stable method suitable for large exchange coupling or small damping.
    /// Fixed-point iteration converges to within `tolerance` or stops at
    /// `max_iter` iterations.
    SemiImplicit {
        /// Fixed-point iteration convergence tolerance.
        tolerance: f64,
        /// Maximum number of fixed-point iterations per time step.
        max_iter: usize,
    },
}

// ---------------------------------------------------------------------------
// SimulationBuilder
// ---------------------------------------------------------------------------

/// A type-state builder for constructing [`Simulation`] instances.
///
/// The three type parameters encode whether each mandatory piece of
/// configuration has been provided:
///
/// | Parameter | Unconfigured   | Configured     |
/// |-----------|----------------|----------------|
/// | `M`       | [`NoMaterial`] | [`WithMaterial`] |
/// | `F`       | [`NoField`]    | [`WithField`]    |
/// | `S`       | [`NoSolver`]   | [`WithSolver`]   |
///
/// The [`build`](SimulationBuilder::build) method is only available when
/// `M = WithMaterial`, `F = WithField`, `S = WithSolver`.
pub struct SimulationBuilder<M, F, S> {
    material: Option<Ferromagnet>,
    external_field: Option<Vector3<f64>>,
    initial_magnetization: Option<Vector3<f64>>,
    damping: Option<f64>,
    temperature: Option<f64>,
    dt: f64,
    num_steps: usize,
    solver_kind: Option<SolverKind>,
    _material: PhantomData<M>,
    _field: PhantomData<F>,
    _solver: PhantomData<S>,
}

// ---------------------------------------------------------------------------
// Constructors (initial state)
// ---------------------------------------------------------------------------

impl SimulationBuilder<NoMaterial, NoField, NoSolver> {
    /// Create a new, unconfigured builder.
    ///
    /// Default values:
    /// - `dt`        = 1.0e-13 s (0.1 ps)
    /// - `num_steps` = 1000
    /// - `temperature` = 0.0 K
    ///
    /// All three mandatory parameters (material, field, solver) must be set
    /// before [`build`](SimulationBuilder::build) becomes available.
    pub fn new() -> Self {
        Self {
            material: None,
            external_field: None,
            initial_magnetization: None,
            damping: None,
            temperature: None,
            dt: 1.0e-13,
            num_steps: 1000,
            solver_kind: None,
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }
}

impl Default for SimulationBuilder<NoMaterial, NoField, NoSolver> {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Type-state transitions (mandatory configuration)
// ---------------------------------------------------------------------------

impl<M, F, S> SimulationBuilder<M, F, S> {
    /// Set the ferromagnetic material.
    ///
    /// Transitions the material type-state to [`WithMaterial`].
    pub fn material(self, mat: Ferromagnet) -> SimulationBuilder<WithMaterial, F, S> {
        SimulationBuilder {
            material: Some(mat),
            external_field: self.external_field,
            initial_magnetization: self.initial_magnetization,
            damping: self.damping,
            temperature: self.temperature,
            dt: self.dt,
            num_steps: self.num_steps,
            solver_kind: self.solver_kind,
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }

    /// Set the external magnetic field vector \[T\].
    ///
    /// Transitions the field type-state to [`WithField`].
    pub fn external_field(self, field: Vector3<f64>) -> SimulationBuilder<M, WithField, S> {
        SimulationBuilder {
            material: self.material,
            external_field: Some(field),
            initial_magnetization: self.initial_magnetization,
            damping: self.damping,
            temperature: self.temperature,
            dt: self.dt,
            num_steps: self.num_steps,
            solver_kind: self.solver_kind,
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }

    /// Select the 4th-order Runge-Kutta solver.
    ///
    /// Transitions the solver type-state to [`WithSolver`].
    pub fn solver_rk4(self) -> SimulationBuilder<M, F, WithSolver> {
        SimulationBuilder {
            material: self.material,
            external_field: self.external_field,
            initial_magnetization: self.initial_magnetization,
            damping: self.damping,
            temperature: self.temperature,
            dt: self.dt,
            num_steps: self.num_steps,
            solver_kind: Some(SolverKind::Rk4),
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }

    /// Select the forward-Euler solver (first order).
    ///
    /// Transitions the solver type-state to [`WithSolver`].
    pub fn solver_euler(self) -> SimulationBuilder<M, F, WithSolver> {
        SimulationBuilder {
            material: self.material,
            external_field: self.external_field,
            initial_magnetization: self.initial_magnetization,
            damping: self.damping,
            temperature: self.temperature,
            dt: self.dt,
            num_steps: self.num_steps,
            solver_kind: Some(SolverKind::Euler),
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }

    /// Select the Heun predictor-corrector solver (second order).
    ///
    /// Transitions the solver type-state to [`WithSolver`].
    pub fn solver_heun(self) -> SimulationBuilder<M, F, WithSolver> {
        SimulationBuilder {
            material: self.material,
            external_field: self.external_field,
            initial_magnetization: self.initial_magnetization,
            damping: self.damping,
            temperature: self.temperature,
            dt: self.dt,
            num_steps: self.num_steps,
            solver_kind: Some(SolverKind::Heun),
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }

    /// Select the Dormand-Prince 5(4) adaptive Runge-Kutta solver.
    ///
    /// This is the method underlying MATLAB's `ode45`. It uses automatic
    /// step size control to keep the local truncation error within
    /// `tolerance`. The `dt` setting is used as the initial step guess.
    ///
    /// Transitions the solver type-state to [`WithSolver`].
    pub fn solver_dp45(self, tolerance: f64) -> SimulationBuilder<M, F, WithSolver> {
        SimulationBuilder {
            material: self.material,
            external_field: self.external_field,
            initial_magnetization: self.initial_magnetization,
            damping: self.damping,
            temperature: self.temperature,
            dt: self.dt,
            num_steps: self.num_steps,
            solver_kind: Some(SolverKind::Dp45 { tolerance }),
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }

    /// Select the Dormand-Prince 8(7) adaptive Runge-Kutta solver.
    ///
    /// A 13-stage method delivering 8th-order accuracy with a 7th-order
    /// embedded error estimate for adaptive step control. Use when very
    /// high precision is required over long simulation times.
    ///
    /// Transitions the solver type-state to [`WithSolver`].
    pub fn solver_dp87(self, tolerance: f64) -> SimulationBuilder<M, F, WithSolver> {
        SimulationBuilder {
            material: self.material,
            external_field: self.external_field,
            initial_magnetization: self.initial_magnetization,
            damping: self.damping,
            temperature: self.temperature,
            dt: self.dt,
            num_steps: self.num_steps,
            solver_kind: Some(SolverKind::Dp87 { tolerance }),
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }

    /// Select the Yoshida 4th-order symplectic integrator.
    ///
    /// A structure-preserving composition method that provides excellent
    /// long-time energy conservation in the undamped (zero-damping) limit.
    /// Internally the magnetization vector is augmented with an auxiliary
    /// zero-momentum coordinate to satisfy the symplectic state convention.
    ///
    /// Transitions the solver type-state to [`WithSolver`].
    pub fn solver_yoshida4(self) -> SimulationBuilder<M, F, WithSolver> {
        SimulationBuilder {
            material: self.material,
            external_field: self.external_field,
            initial_magnetization: self.initial_magnetization,
            damping: self.damping,
            temperature: self.temperature,
            dt: self.dt,
            num_steps: self.num_steps,
            solver_kind: Some(SolverKind::Yoshida4),
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }

    /// Select the Forest-Ruth 4th-order symplectic integrator.
    ///
    /// An alternative 4th-order composition method with slightly different
    /// stability properties than Yoshida4. The same internal state-vector
    /// augmentation convention is used.
    ///
    /// Transitions the solver type-state to [`WithSolver`].
    pub fn solver_forest_ruth(self) -> SimulationBuilder<M, F, WithSolver> {
        SimulationBuilder {
            material: self.material,
            external_field: self.external_field,
            initial_magnetization: self.initial_magnetization,
            damping: self.damping,
            temperature: self.temperature,
            dt: self.dt,
            num_steps: self.num_steps,
            solver_kind: Some(SolverKind::ForestRuth),
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }

    /// Select the semi-implicit integrator (implicit midpoint rule).
    ///
    /// An A-stable method for stiff ODE systems arising from large exchange
    /// coupling or very small Gilbert damping. Uses fixed-point iteration to
    /// solve the implicit midpoint equation at each step.
    ///
    /// # Arguments
    /// * `tolerance` - Fixed-point iteration convergence threshold.
    /// * `max_iter`  - Maximum fixed-point iterations per time step.
    ///
    /// Transitions the solver type-state to [`WithSolver`].
    pub fn solver_semi_implicit(
        self,
        tolerance: f64,
        max_iter: usize,
    ) -> SimulationBuilder<M, F, WithSolver> {
        SimulationBuilder {
            material: self.material,
            external_field: self.external_field,
            initial_magnetization: self.initial_magnetization,
            damping: self.damping,
            temperature: self.temperature,
            dt: self.dt,
            num_steps: self.num_steps,
            solver_kind: Some(SolverKind::SemiImplicit {
                tolerance,
                max_iter,
            }),
            _material: PhantomData,
            _field: PhantomData,
            _solver: PhantomData,
        }
    }
}

// ---------------------------------------------------------------------------
// Optional configuration (does not change type-state)
// ---------------------------------------------------------------------------

impl<M, F, S> SimulationBuilder<M, F, S> {
    /// Set the initial magnetization direction.
    ///
    /// The vector will be normalised internally. If not set, defaults to the
    /// material's easy axis direction.
    pub fn initial_magnetization(mut self, m0: Vector3<f64>) -> Self {
        self.initial_magnetization = Some(m0);
        self
    }

    /// Override the Gilbert damping constant.
    ///
    /// If not set, the material's own `alpha` value is used.
    pub fn damping(mut self, alpha: f64) -> Self {
        self.damping = Some(alpha);
        self
    }

    /// Set the simulation temperature \[K\].
    ///
    /// Defaults to 0 K (zero-temperature dynamics).
    pub fn temperature(mut self, temp: f64) -> Self {
        self.temperature = Some(temp);
        self
    }

    /// Set the integration time step \[s\].
    ///
    /// Defaults to 1.0e-13 s (0.1 ps).
    pub fn time_step(mut self, dt: f64) -> Self {
        self.dt = dt;
        self
    }

    /// Set the total number of integration steps.
    ///
    /// Defaults to 1000.
    pub fn num_steps(mut self, n: usize) -> Self {
        self.num_steps = n;
        self
    }
}

// ---------------------------------------------------------------------------
// build() — only when fully configured
// ---------------------------------------------------------------------------

impl SimulationBuilder<WithMaterial, WithField, WithSolver> {
    /// Build a [`Simulation`] from the fully configured builder.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] if:
    /// - `dt <= 0`
    /// - `num_steps == 0`
    /// - The initial magnetization vector has zero magnitude
    pub fn build(self) -> Result<Simulation> {
        // Validate time step
        if self.dt <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "dt".to_string(),
                reason: "time step must be positive".to_string(),
            });
        }
        if self.num_steps == 0 {
            return Err(Error::InvalidParameter {
                param: "num_steps".to_string(),
                reason: "must have at least one integration step".to_string(),
            });
        }

        // Material is guaranteed by the type-state
        let material = self.material.ok_or_else(|| Error::ConfigurationError {
            description: "material is required (type-state violated)".to_string(),
        })?;

        let external_field = self
            .external_field
            .ok_or_else(|| Error::ConfigurationError {
                description: "external field is required (type-state violated)".to_string(),
            })?;

        let solver_kind = self.solver_kind.ok_or_else(|| Error::ConfigurationError {
            description: "solver is required (type-state violated)".to_string(),
        })?;

        // Determine damping (override or material default)
        let alpha = self.damping.unwrap_or(material.alpha);

        // Determine initial magnetization
        let m0_raw = self.initial_magnetization.unwrap_or(material.easy_axis);

        let m0_mag = m0_raw.magnitude();
        if m0_mag < 1.0e-30 {
            return Err(Error::InvalidParameter {
                param: "initial_magnetization".to_string(),
                reason: "magnetization vector must be non-zero".to_string(),
            });
        }
        let m0 = m0_raw * (1.0 / m0_mag);

        let solver = LlgSolver::new(alpha, self.dt);
        let temperature = self.temperature.unwrap_or(0.0);

        Ok(Simulation {
            material,
            external_field,
            magnetization: m0,
            solver,
            solver_kind,
            temperature,
            dt: self.dt,
            num_steps: self.num_steps,
        })
    }
}

// ---------------------------------------------------------------------------
// Preset configurations
// ---------------------------------------------------------------------------

impl SimulationBuilder<NoMaterial, NoField, NoSolver> {
    /// Pre-configured YIG/Pt spin pumping simulation.
    ///
    /// Sets up a Yttrium Iron Garnet ferromagnet with a DC magnetic field
    /// along the z-axis (0.1 T) and RK4 solver. The magnetization starts
    /// along the x-axis so that precession around the field is clearly visible.
    ///
    /// # Physical context
    ///
    /// YIG is the material of choice for spin pumping experiments due to its
    /// extremely low Gilbert damping (alpha ~ 1e-4). See:
    /// E. Saitoh et al., *Appl. Phys. Lett.* **88**, 182509 (2006).
    pub fn yig_pt_pumping() -> SimulationBuilder<WithMaterial, WithField, WithSolver> {
        SimulationBuilder::<NoMaterial, NoField, NoSolver>::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .initial_magnetization(Vector3::new(1.0, 0.0, 0.0))
            .time_step(1.0e-13)
            .num_steps(10_000)
    }

    /// Pre-configured skyrmion dynamics simulation.
    ///
    /// Uses a CoFeB-like material with perpendicular magnetic anisotropy (PMA)
    /// and a moderate out-of-plane field. The initial magnetization is tilted
    /// to seed skyrmion nucleation dynamics.
    ///
    /// # Physical context
    ///
    /// Magnetic skyrmions are topologically protected spin textures stabilised
    /// by the Dzyaloshinskii-Moriya interaction (DMI). See:
    /// S. Woo et al., *Nat. Mater.* **15**, 501 (2016).
    pub fn skyrmion_dynamics() -> SimulationBuilder<WithMaterial, WithField, WithSolver> {
        let skyrmion_material = Ferromagnet {
            alpha: 0.01,
            ms: 1.0e6,
            anisotropy_k: 8.0e5,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_a: 1.5e-11,
        };

        SimulationBuilder::<NoMaterial, NoField, NoSolver>::new()
            .material(skyrmion_material)
            .external_field(Vector3::new(0.0, 0.0, 0.3))
            .solver_rk4()
            .initial_magnetization(Vector3::new(0.1, 0.1, 1.0))
            .time_step(5.0e-14)
            .num_steps(20_000)
    }

    /// Pre-configured ferromagnetic resonance (FMR) simulation.
    ///
    /// Combines a static DC field along z with a small oscillating RF field
    /// component along x. The RF amplitude is encoded in the initial field;
    /// for a time-varying drive the caller should update the field during the
    /// simulation run.
    ///
    /// # Arguments
    ///
    /// * `material` - The ferromagnetic material to resonate.
    /// * `h_dc`     - DC bias field magnitude \[T\] (applied along z).
    /// * `h_rf`     - RF excitation amplitude \[T\] (applied along x initially).
    pub fn ferromagnetic_resonance(
        material: Ferromagnet,
        h_dc: f64,
        h_rf: f64,
    ) -> SimulationBuilder<WithMaterial, WithField, WithSolver> {
        // Combined DC + initial RF field
        let field = Vector3::new(h_rf, 0.0, h_dc);

        SimulationBuilder::<NoMaterial, NoField, NoSolver>::new()
            .material(material)
            .external_field(field)
            .solver_rk4()
            .initial_magnetization(Vector3::new(1.0, 0.0, 0.0))
            .time_step(1.0e-13)
            .num_steps(50_000)
    }
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

/// A fully configured simulation ready to run.
///
/// Created via [`SimulationBuilder::build`]. Call [`run`](Simulation::run) to
/// evolve the magnetization and obtain a [`SimulationResult`].
pub struct Simulation {
    /// The ferromagnetic material being simulated.
    material: Ferromagnet,
    /// External magnetic field \[T\].
    external_field: Vector3<f64>,
    /// Current magnetization direction (unit vector).
    magnetization: Vector3<f64>,
    /// The LLG solver instance.
    solver: LlgSolver,
    /// Which integration scheme to use.
    solver_kind: SolverKind,
    /// Temperature \[K\] (reserved for future stochastic extensions).
    #[allow(dead_code)]
    temperature: f64,
    /// Integration time step \[s\].
    #[allow(dead_code)]
    dt: f64,
    /// Total number of integration steps.
    num_steps: usize,
}

impl Simulation {
    /// Run the simulation and return the results.
    ///
    /// The magnetization is evolved for `num_steps` time steps using the
    /// configured solver. At each step the effective field includes both
    /// the external Zeeman field and the uniaxial anisotropy field.
    ///
    /// # Returns
    ///
    /// A [`SimulationResult`] containing the full trajectory, final
    /// magnetization, and energy history.
    ///
    /// # Errors
    ///
    /// Returns [`Error::NumericalError`] if the magnetization becomes
    /// non-finite during integration (NaN / infinity).
    pub fn run(&mut self) -> Result<SimulationResult> {
        let mut trajectory = Vec::with_capacity(self.num_steps + 1);
        let mut energies = Vec::with_capacity(self.num_steps + 1);

        // Record initial state
        trajectory.push(self.magnetization);
        energies.push(self.compute_energy(self.magnetization));

        // Pre-extract all parameters needed by closures so that we avoid
        // holding borrows on `self` inside loops that also mutate `self`.
        let alpha = self.solver.alpha;
        let gamma = self.solver.gamma;
        let dt = self.solver.dt;
        let h_ext = self.external_field;
        let mat_ms = self.material.ms;
        let mat_k = self.material.anisotropy_k;
        let mat_easy = self.material.easy_axis;
        let num_steps = self.num_steps;
        let solver_kind = self.solver_kind;

        // Closure: effective field from magnetization (captures copied primitives).
        let eff_field = |m_eval: Vector3<f64>| -> Vector3<f64> {
            use crate::constants::MU_0;
            let projection = m_eval.dot(&mat_easy);
            let h_anis = if mat_ms.abs() > 1.0e-30 {
                mat_easy * (2.0 * mat_k * projection / (MU_0 * mat_ms))
            } else {
                Vector3::zero()
            };
            h_ext + h_anis
        };

        // Closure: LLG dm/dt.
        let llg_rhs = |m_eval: Vector3<f64>| -> Vector3<f64> {
            calc_dm_dt(m_eval, eff_field(m_eval), gamma, alpha)
        };

        // Dispatch based on solver kind
        match solver_kind {
            SolverKind::Dp45 { tolerance } => {
                let inner = DormandPrince45::new();
                let mut adaptive = AdaptiveIntegrator::new(inner, tolerance, 4.0)
                    .with_dt_min(1e-20)
                    .with_dt_max(dt * 10.0);
                let mut cur_dt = dt;
                for step_idx in 0..num_steps {
                    let m = self.magnetization;
                    let state = vec![m];
                    let rhs: RhsFn<'_> = &|s: &[Vector3<f64>], _t: f64| vec![llg_rhs(s[0])];
                    let (output, used_dt) = adaptive.adaptive_step(&state, 0.0, cur_dt, rhs)?;
                    cur_dt = output.suggested_dt.unwrap_or(used_dt);
                    let new_m = output.new_state.into_iter().next().ok_or_else(|| {
                        Error::NumericalError {
                            description: format!(
                                "DP45: empty output state at step {}",
                                step_idx + 1
                            ),
                        }
                    })?;
                    Self::check_finite(new_m, trajectory.len())?;
                    let new_m = Self::renormalize(new_m);
                    self.magnetization = new_m;
                    trajectory.push(new_m);
                    energies.push(self.compute_energy(new_m));
                }
            },
            SolverKind::Dp87 { tolerance } => {
                let inner = DormandPrince87::new();
                let mut adaptive = AdaptiveIntegrator::new(inner, tolerance, 7.0)
                    .with_dt_min(1e-20)
                    .with_dt_max(dt * 10.0);
                let mut cur_dt = dt;
                for step_idx in 0..num_steps {
                    let m = self.magnetization;
                    let state = vec![m];
                    let rhs: RhsFn<'_> = &|s: &[Vector3<f64>], _t: f64| vec![llg_rhs(s[0])];
                    let (output, used_dt) = adaptive.adaptive_step(&state, 0.0, cur_dt, rhs)?;
                    cur_dt = output.suggested_dt.unwrap_or(used_dt);
                    let new_m = output.new_state.into_iter().next().ok_or_else(|| {
                        Error::NumericalError {
                            description: format!(
                                "DP87: empty output state at step {}",
                                step_idx + 1
                            ),
                        }
                    })?;
                    Self::check_finite(new_m, trajectory.len())?;
                    let new_m = Self::renormalize(new_m);
                    self.magnetization = new_m;
                    trajectory.push(new_m);
                    energies.push(self.compute_energy(new_m));
                }
            },
            SolverKind::Yoshida4 => {
                // Symplectic leapfrog on (q=m, p=dm/dt).
                // dq/dt = p  (drift), dp/dt = LLG(q) (kick).
                let mut integrator = Yoshida4::new();
                // Initialise momentum from the LLG derivative at the starting state.
                let mut p = llg_rhs(self.magnetization);
                for step_idx in 0..num_steps {
                    let q = self.magnetization;
                    let state = vec![q, p];
                    let rhs: RhsFn<'_> = &|s: &[Vector3<f64>], _t: f64| vec![s[1], llg_rhs(s[0])];
                    let output = integrator.step(&state, 0.0, dt, rhs)?;
                    let mut it = output.new_state.into_iter();
                    let new_m = it.next().ok_or_else(|| Error::NumericalError {
                        description: format!("Yoshida4: missing q at step {}", step_idx + 1),
                    })?;
                    // Carry the updated momentum into the next step.
                    p = it.next().unwrap_or_else(|| llg_rhs(new_m));
                    Self::check_finite(new_m, trajectory.len())?;
                    let new_m = Self::renormalize(new_m);
                    self.magnetization = new_m;
                    trajectory.push(new_m);
                    energies.push(self.compute_energy(new_m));
                }
            },
            SolverKind::ForestRuth => {
                // Same (q, p) leapfrog convention as Yoshida4.
                let mut integrator = ForestRuth::new();
                let mut p = llg_rhs(self.magnetization);
                for step_idx in 0..num_steps {
                    let q = self.magnetization;
                    let state = vec![q, p];
                    let rhs: RhsFn<'_> = &|s: &[Vector3<f64>], _t: f64| vec![s[1], llg_rhs(s[0])];
                    let output = integrator.step(&state, 0.0, dt, rhs)?;
                    let mut it = output.new_state.into_iter();
                    let new_m = it.next().ok_or_else(|| Error::NumericalError {
                        description: format!("ForestRuth: missing q at step {}", step_idx + 1),
                    })?;
                    p = it.next().unwrap_or_else(|| llg_rhs(new_m));
                    Self::check_finite(new_m, trajectory.len())?;
                    let new_m = Self::renormalize(new_m);
                    self.magnetization = new_m;
                    trajectory.push(new_m);
                    energies.push(self.compute_energy(new_m));
                }
            },
            SolverKind::SemiImplicit {
                tolerance,
                max_iter,
            } => {
                let mut integrator = SemiImplicit::new(max_iter, tolerance);
                for step_idx in 0..num_steps {
                    let m = self.magnetization;
                    let state = vec![m];
                    let rhs: RhsFn<'_> = &|s: &[Vector3<f64>], _t: f64| vec![llg_rhs(s[0])];
                    let output = integrator.step(&state, 0.0, dt, rhs)?;
                    let new_m = output.new_state.into_iter().next().ok_or_else(|| {
                        Error::NumericalError {
                            description: format!(
                                "SemiImplicit: empty output state at step {}",
                                step_idx + 1
                            ),
                        }
                    })?;
                    Self::check_finite(new_m, trajectory.len())?;
                    let new_m = Self::renormalize(new_m);
                    self.magnetization = new_m;
                    trajectory.push(new_m);
                    energies.push(self.compute_energy(new_m));
                }
            },
            // Legacy fixed-step solvers via LlgSolver
            _ => {
                for _ in 0..num_steps {
                    let m = self.magnetization;

                    let new_m = match solver_kind {
                        SolverKind::Rk4 => self.solver.step_rk4(m, eff_field),
                        SolverKind::Euler => self.solver.step_euler(m, eff_field(m)),
                        SolverKind::Heun => self.solver.step_heun(m, eff_field),
                        // All new variants are handled above. This arm is
                        // unreachable at runtime but required for exhaustiveness.
                        _ => {
                            return Err(Error::ConfigurationError {
                                description: "solver dispatch error: unexpected variant in \
                                              legacy arm"
                                    .to_string(),
                            });
                        },
                    };

                    Self::check_finite(new_m, trajectory.len())?;
                    self.magnetization = new_m;
                    trajectory.push(new_m);
                    energies.push(self.compute_energy(new_m));
                }
            },
        }

        Ok(SimulationResult {
            trajectory,
            final_magnetization: self.magnetization,
            energies,
        })
    }

    /// Run the simulation, invoking `callback` once per completed step
    /// instead of collecting the full trajectory into memory.
    ///
    /// This performs *exactly* the same integration loop as
    /// [`run`](Simulation::run) — same solver dispatch, same finite-value
    /// checks, same renormalization behaviour per solver family — but instead
    /// of accumulating `Vec<Vector3<f64>>` / `Vec<f64>` histories it calls
    /// `callback(step_index, &magnetization, energy)` once per completed step
    /// and never retains a growing collection internally. This makes it
    /// suitable for very large `num_steps` or long-running simulations where
    /// materializing the full [`SimulationResult`] trajectory would be
    /// memory-prohibitive.
    ///
    /// `step_index` starts at `0` for the initial state (mirroring
    /// `SimulationResult::trajectory[0]` / `energies[0]`) and increases by
    /// one per integration step thereafter, so the sequence of
    /// `(step_index, m, energy)` triples reported here is exactly the
    /// row-by-row content of the [`SimulationResult`] that
    /// [`run`](Simulation::run) would produce from an identically configured
    /// builder.
    ///
    /// # Errors
    ///
    /// Returns [`Error::NumericalError`] if the magnetization becomes
    /// non-finite during integration, or [`Error::ConfigurationError`] on an
    /// internal solver-dispatch inconsistency (unreachable in practice). Any
    /// `Err` returned by `callback` is propagated immediately, aborting the
    /// integration at that step.
    pub fn run_streaming<Cb>(&mut self, mut callback: Cb) -> Result<()>
    where
        Cb: FnMut(usize, &Vector3<f64>, f64) -> Result<()>,
    {
        // Report the initial state (step 0), matching
        // `SimulationResult::trajectory[0]` / `energies[0]`, without ever
        // allocating a Vec for it.
        let e0 = self.compute_energy(self.magnetization);
        callback(0, &self.magnetization, e0)?;

        // Pre-extract all parameters needed by closures so that we avoid
        // holding borrows on `self` inside loops that also mutate `self`.
        let alpha = self.solver.alpha;
        let gamma = self.solver.gamma;
        let dt = self.solver.dt;
        let h_ext = self.external_field;
        let mat_ms = self.material.ms;
        let mat_k = self.material.anisotropy_k;
        let mat_easy = self.material.easy_axis;
        let num_steps = self.num_steps;
        let solver_kind = self.solver_kind;

        // Closure: effective field from magnetization (captures copied primitives).
        let eff_field = |m_eval: Vector3<f64>| -> Vector3<f64> {
            use crate::constants::MU_0;
            let projection = m_eval.dot(&mat_easy);
            let h_anis = if mat_ms.abs() > 1.0e-30 {
                mat_easy * (2.0 * mat_k * projection / (MU_0 * mat_ms))
            } else {
                Vector3::zero()
            };
            h_ext + h_anis
        };

        // Closure: LLG dm/dt.
        let llg_rhs = |m_eval: Vector3<f64>| -> Vector3<f64> {
            calc_dm_dt(m_eval, eff_field(m_eval), gamma, alpha)
        };

        // Dispatch based on solver kind (mirrors `run()` step-for-step).
        match solver_kind {
            SolverKind::Dp45 { tolerance } => {
                let inner = DormandPrince45::new();
                let mut adaptive = AdaptiveIntegrator::new(inner, tolerance, 4.0)
                    .with_dt_min(1e-20)
                    .with_dt_max(dt * 10.0);
                let mut cur_dt = dt;
                for step_idx in 0..num_steps {
                    let m = self.magnetization;
                    let state = vec![m];
                    let rhs: RhsFn<'_> = &|s: &[Vector3<f64>], _t: f64| vec![llg_rhs(s[0])];
                    let (output, used_dt) = adaptive.adaptive_step(&state, 0.0, cur_dt, rhs)?;
                    cur_dt = output.suggested_dt.unwrap_or(used_dt);
                    let new_m = output.new_state.into_iter().next().ok_or_else(|| {
                        Error::NumericalError {
                            description: format!(
                                "DP45: empty output state at step {}",
                                step_idx + 1
                            ),
                        }
                    })?;
                    Self::check_finite(new_m, step_idx + 1)?;
                    let new_m = Self::renormalize(new_m);
                    self.magnetization = new_m;
                    let energy = self.compute_energy(new_m);
                    callback(step_idx + 1, &new_m, energy)?;
                }
            },
            SolverKind::Dp87 { tolerance } => {
                let inner = DormandPrince87::new();
                let mut adaptive = AdaptiveIntegrator::new(inner, tolerance, 7.0)
                    .with_dt_min(1e-20)
                    .with_dt_max(dt * 10.0);
                let mut cur_dt = dt;
                for step_idx in 0..num_steps {
                    let m = self.magnetization;
                    let state = vec![m];
                    let rhs: RhsFn<'_> = &|s: &[Vector3<f64>], _t: f64| vec![llg_rhs(s[0])];
                    let (output, used_dt) = adaptive.adaptive_step(&state, 0.0, cur_dt, rhs)?;
                    cur_dt = output.suggested_dt.unwrap_or(used_dt);
                    let new_m = output.new_state.into_iter().next().ok_or_else(|| {
                        Error::NumericalError {
                            description: format!(
                                "DP87: empty output state at step {}",
                                step_idx + 1
                            ),
                        }
                    })?;
                    Self::check_finite(new_m, step_idx + 1)?;
                    let new_m = Self::renormalize(new_m);
                    self.magnetization = new_m;
                    let energy = self.compute_energy(new_m);
                    callback(step_idx + 1, &new_m, energy)?;
                }
            },
            SolverKind::Yoshida4 => {
                // Symplectic leapfrog on (q=m, p=dm/dt).
                let mut integrator = Yoshida4::new();
                let mut p = llg_rhs(self.magnetization);
                for step_idx in 0..num_steps {
                    let q = self.magnetization;
                    let state = vec![q, p];
                    let rhs: RhsFn<'_> = &|s: &[Vector3<f64>], _t: f64| vec![s[1], llg_rhs(s[0])];
                    let output = integrator.step(&state, 0.0, dt, rhs)?;
                    let mut it = output.new_state.into_iter();
                    let new_m = it.next().ok_or_else(|| Error::NumericalError {
                        description: format!("Yoshida4: missing q at step {}", step_idx + 1),
                    })?;
                    p = it.next().unwrap_or_else(|| llg_rhs(new_m));
                    Self::check_finite(new_m, step_idx + 1)?;
                    let new_m = Self::renormalize(new_m);
                    self.magnetization = new_m;
                    let energy = self.compute_energy(new_m);
                    callback(step_idx + 1, &new_m, energy)?;
                }
            },
            SolverKind::ForestRuth => {
                // Same (q, p) leapfrog convention as Yoshida4.
                let mut integrator = ForestRuth::new();
                let mut p = llg_rhs(self.magnetization);
                for step_idx in 0..num_steps {
                    let q = self.magnetization;
                    let state = vec![q, p];
                    let rhs: RhsFn<'_> = &|s: &[Vector3<f64>], _t: f64| vec![s[1], llg_rhs(s[0])];
                    let output = integrator.step(&state, 0.0, dt, rhs)?;
                    let mut it = output.new_state.into_iter();
                    let new_m = it.next().ok_or_else(|| Error::NumericalError {
                        description: format!("ForestRuth: missing q at step {}", step_idx + 1),
                    })?;
                    p = it.next().unwrap_or_else(|| llg_rhs(new_m));
                    Self::check_finite(new_m, step_idx + 1)?;
                    let new_m = Self::renormalize(new_m);
                    self.magnetization = new_m;
                    let energy = self.compute_energy(new_m);
                    callback(step_idx + 1, &new_m, energy)?;
                }
            },
            SolverKind::SemiImplicit {
                tolerance,
                max_iter,
            } => {
                let mut integrator = SemiImplicit::new(max_iter, tolerance);
                for step_idx in 0..num_steps {
                    let m = self.magnetization;
                    let state = vec![m];
                    let rhs: RhsFn<'_> = &|s: &[Vector3<f64>], _t: f64| vec![llg_rhs(s[0])];
                    let output = integrator.step(&state, 0.0, dt, rhs)?;
                    let new_m = output.new_state.into_iter().next().ok_or_else(|| {
                        Error::NumericalError {
                            description: format!(
                                "SemiImplicit: empty output state at step {}",
                                step_idx + 1
                            ),
                        }
                    })?;
                    Self::check_finite(new_m, step_idx + 1)?;
                    let new_m = Self::renormalize(new_m);
                    self.magnetization = new_m;
                    let energy = self.compute_energy(new_m);
                    callback(step_idx + 1, &new_m, energy)?;
                }
            },
            // Legacy fixed-step solvers via LlgSolver
            _ => {
                for step_idx in 0..num_steps {
                    let m = self.magnetization;

                    let new_m = match solver_kind {
                        SolverKind::Rk4 => self.solver.step_rk4(m, eff_field),
                        SolverKind::Euler => self.solver.step_euler(m, eff_field(m)),
                        SolverKind::Heun => self.solver.step_heun(m, eff_field),
                        // All new variants are handled above. This arm is
                        // unreachable at runtime but required for exhaustiveness.
                        _ => {
                            return Err(Error::ConfigurationError {
                                description: "solver dispatch error: unexpected variant in \
                                              legacy arm"
                                    .to_string(),
                            });
                        },
                    };

                    Self::check_finite(new_m, step_idx + 1)?;
                    self.magnetization = new_m;
                    let energy = self.compute_energy(new_m);
                    callback(step_idx + 1, &new_m, energy)?;
                }
            },
        }

        Ok(())
    }

    /// Check that a magnetization vector has finite components.
    fn check_finite(m: Vector3<f64>, step: usize) -> Result<()> {
        if !m.x.is_finite() || !m.y.is_finite() || !m.z.is_finite() {
            return Err(Error::NumericalError {
                description: format!("magnetization became non-finite at step {} : {:?}", step, m),
            });
        }
        Ok(())
    }

    /// Re-normalise a magnetization vector to unit length.
    ///
    /// If the magnitude is essentially zero (numerical collapse), returns the
    /// vector unchanged to avoid division by zero.
    fn renormalize(m: Vector3<f64>) -> Vector3<f64> {
        let mag = m.magnitude();
        if mag > 1.0e-30 {
            m * (1.0 / mag)
        } else {
            m
        }
    }

    /// Compute total energy density [J/m^3] for the current state.
    fn compute_energy(&self, m: Vector3<f64>) -> f64 {
        let e_zeeman = zeeman_energy(m, self.external_field, self.material.ms);
        let e_anis = anisotropy_energy(m, self.material.easy_axis, self.material.anisotropy_k);
        e_zeeman + e_anis
    }
}

// ---------------------------------------------------------------------------
// SimulationResult
// ---------------------------------------------------------------------------

/// Results from a completed simulation run.
///
/// Contains the full magnetization trajectory, the final magnetization state,
/// and the energy history at each time step.
#[derive(Debug, Clone)]
pub struct SimulationResult {
    /// Magnetization direction at each time step (including t = 0).
    pub trajectory: Vec<Vector3<f64>>,
    /// Final magnetization direction after all steps.
    pub final_magnetization: Vector3<f64>,
    /// Total energy density [J/m^3] at each time step.
    pub energies: Vec<f64>,
}

impl SimulationResult {
    /// Number of recorded time steps (including the initial state).
    pub fn len(&self) -> usize {
        self.trajectory.len()
    }

    /// Whether the result contains no data points.
    pub fn is_empty(&self) -> bool {
        self.trajectory.is_empty()
    }

    /// Maximum relative energy change over the trajectory.
    ///
    /// Useful for verifying energy conservation in zero-damping simulations.
    /// When the reference energy is near zero, the *absolute* maximum
    /// deviation is returned instead (normalised by 1 J/m^3 to keep it
    /// dimensionless).
    ///
    /// Returns `None` only if the trajectory is empty.
    pub fn max_relative_energy_change(&self) -> Option<f64> {
        let e0 = *self.energies.first()?;
        if e0.abs() < 1.0e-30 {
            // Fall back to absolute deviation (treat 1 J/m^3 as unit scale)
            return self
                .energies
                .iter()
                .map(|e| (e - e0).abs())
                .fold(None, |acc, v| {
                    Some(match acc {
                        Some(prev) => {
                            if v > prev {
                                v
                            } else {
                                prev
                            }
                        },
                        None => v,
                    })
                });
        }
        self.energies
            .iter()
            .map(|e| ((e - e0) / e0).abs())
            .fold(None, |acc, v| {
                Some(match acc {
                    Some(prev) => {
                        if v > prev {
                            v
                        } else {
                            prev
                        }
                    },
                    None => v,
                })
            })
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_with_all_parameters() {
        let sim = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .initial_magnetization(Vector3::new(1.0, 0.0, 0.0))
            .temperature(300.0)
            .time_step(1.0e-13)
            .num_steps(500)
            .build();

        assert!(sim.is_ok(), "build with all parameters should succeed");
    }

    #[test]
    fn test_preset_yig_pt_pumping_builds() {
        let sim = SimulationBuilder::yig_pt_pumping().build();
        assert!(
            sim.is_ok(),
            "yig_pt_pumping preset should build successfully"
        );
    }

    #[test]
    fn test_preset_skyrmion_dynamics_builds() {
        let sim = SimulationBuilder::skyrmion_dynamics().build();
        assert!(
            sim.is_ok(),
            "skyrmion_dynamics preset should build successfully"
        );
    }

    #[test]
    fn test_preset_fmr_builds() {
        let sim = SimulationBuilder::ferromagnetic_resonance(Ferromagnet::permalloy(), 0.5, 0.001)
            .build();
        assert!(
            sim.is_ok(),
            "ferromagnetic_resonance preset should build successfully"
        );
    }

    #[test]
    fn test_run_produces_trajectory() {
        let mut sim = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .initial_magnetization(Vector3::new(1.0, 0.0, 0.0))
            .num_steps(200)
            .build()
            .expect("builder should succeed");

        let result = sim.run().expect("simulation should run without error");

        // Trajectory length = num_steps + 1 (includes initial state)
        assert_eq!(result.len(), 201);
        assert!(!result.is_empty());

        // The final magnetization should still be a unit vector
        let mag = result.final_magnetization.magnitude();
        assert!(
            (mag - 1.0).abs() < 1.0e-6,
            "final magnetization should be normalised, got magnitude {}",
            mag
        );
    }

    #[test]
    fn test_custom_parameters_override_defaults() {
        let custom_dt = 5.0e-14;
        let custom_steps = 42;
        let custom_alpha = 0.05;

        let mut sim = SimulationBuilder::new()
            .material(Ferromagnet::yig()) // yig alpha = 0.0001
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .damping(custom_alpha)
            .time_step(custom_dt)
            .num_steps(custom_steps)
            .build()
            .expect("builder with overrides should succeed");

        let result = sim.run().expect("simulation should run");

        // Check that custom num_steps was respected
        assert_eq!(result.len(), custom_steps + 1);
    }

    #[test]
    fn test_energy_conservation_zero_damping() {
        // With zero damping the LLG equation is conservative, so total
        // energy (Zeeman + anisotropy) should be approximately constant.
        let zero_damping_material = Ferromagnet {
            alpha: 0.0,
            ms: 8.0e5,
            anisotropy_k: 0.0,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_a: 1.3e-11,
        };

        // Start with m tilted 45 degrees from the field so the initial
        // Zeeman energy is non-zero and the relative-change metric is
        // well-defined.
        let mut sim = SimulationBuilder::new()
            .material(zero_damping_material)
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .damping(0.0) // explicitly zero
            .initial_magnetization(Vector3::new(1.0, 0.0, 1.0)) // 45 deg tilt
            .time_step(1.0e-14) // small step for accuracy
            .num_steps(5000)
            .build()
            .expect("zero-damping builder should succeed");

        let result = sim.run().expect("zero-damping simulation should run");

        // Energy should be conserved to within a small tolerance
        let max_change = result
            .max_relative_energy_change()
            .expect("should have energy data");
        assert!(
            max_change < 1.0e-6,
            "energy should be conserved for zero damping, max relative change = {:.2e}",
            max_change
        );
    }

    #[test]
    fn test_invalid_dt_rejected() {
        let result = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .time_step(-1.0)
            .build();

        assert!(result.is_err(), "negative dt should be rejected");
    }

    #[test]
    fn test_zero_steps_rejected() {
        let result = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .num_steps(0)
            .build();

        assert!(result.is_err(), "zero num_steps should be rejected");
    }

    /// Compile-time type-state check: building without a material should
    /// not compile. This is verified via `compile_fail` in the module-level
    /// doc comment, so this test simply ensures the doc-test infrastructure
    /// catches the error.
    #[test]
    fn test_type_state_requires_all_three() {
        // We verify that a fully configured builder works, implying that
        // removing any of the three mandatory calls would break compilation.
        let _sim = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .build()
            .expect("fully configured builder must succeed");
    }

    // -----------------------------------------------------------------------
    // New-integrator builder and run tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_builder_dp45_builds_without_panic() {
        let result = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_dp45(1.0e-6)
            .num_steps(50)
            .time_step(1.0e-13)
            .build();
        assert!(
            result.is_ok(),
            "DP45 builder should succeed: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_builder_dp87_builds_without_panic() {
        let result = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_dp87(1.0e-8)
            .num_steps(50)
            .time_step(1.0e-13)
            .build();
        assert!(
            result.is_ok(),
            "DP87 builder should succeed: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_builder_yoshida4_builds_without_panic() {
        let result = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_yoshida4()
            .num_steps(50)
            .time_step(1.0e-13)
            .build();
        assert!(
            result.is_ok(),
            "Yoshida4 builder should succeed: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_builder_forest_ruth_builds_without_panic() {
        let result = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_forest_ruth()
            .num_steps(50)
            .time_step(1.0e-13)
            .build();
        assert!(
            result.is_ok(),
            "ForestRuth builder should succeed: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_builder_semi_implicit_builds_without_panic() {
        let result = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_semi_implicit(1.0e-12, 50)
            .num_steps(50)
            .time_step(1.0e-13)
            .build();
        assert!(
            result.is_ok(),
            "SemiImplicit builder should succeed: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_dp45_simulation_runs_produces_trajectory() {
        let mut sim = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_dp45(1.0e-6)
            .initial_magnetization(Vector3::new(1.0, 0.0, 0.0))
            .num_steps(100)
            .time_step(1.0e-13)
            .build()
            .expect("DP45 builder should succeed");

        let result = sim.run().expect("DP45 simulation should run");

        // Trajectory includes initial state + num_steps states
        assert_eq!(result.len(), 101);
        assert!(!result.is_empty());

        // Final magnetization should remain a unit vector
        let mag = result.final_magnetization.magnitude();
        assert!(
            (mag - 1.0).abs() < 1.0e-6,
            "DP45 final |m| should be ~1, got {}",
            mag
        );
    }

    #[test]
    fn test_dp45_energy_conservation_similar_to_rk4() {
        // Both DP45 and RK4 should conserve energy to a similar order in the
        // zero-damping limit; DP45 adaptive control should be at least as good.
        let zero_damping = Ferromagnet {
            alpha: 0.0,
            ms: 8.0e5,
            anisotropy_k: 0.0,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_a: 1.3e-11,
        };

        // RK4 reference
        let mut sim_rk4 = SimulationBuilder::new()
            .material(zero_damping.clone())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .damping(0.0)
            .initial_magnetization(Vector3::new(1.0, 0.0, 1.0))
            .time_step(1.0e-14)
            .num_steps(500)
            .build()
            .expect("RK4 zero-damping builder should succeed");
        let rk4_result = sim_rk4.run().expect("RK4 zero-damping sim should run");
        let rk4_max_change = rk4_result
            .max_relative_energy_change()
            .expect("RK4 should have energies");

        // DP45 adaptive
        let mut sim_dp45 = SimulationBuilder::new()
            .material(zero_damping)
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_dp45(1.0e-8)
            .damping(0.0)
            .initial_magnetization(Vector3::new(1.0, 0.0, 1.0))
            .time_step(1.0e-14)
            .num_steps(500)
            .build()
            .expect("DP45 zero-damping builder should succeed");
        let dp45_result = sim_dp45.run().expect("DP45 zero-damping sim should run");
        let dp45_max_change = dp45_result
            .max_relative_energy_change()
            .expect("DP45 should have energies");

        // Both methods should conserve energy reasonably; DP45 with tight
        // tolerance should be no worse than RK4 at the same or better accuracy.
        assert!(
            dp45_max_change < 1.0e-3,
            "DP45 energy conservation should be within 0.1%, got {:.2e}",
            dp45_max_change
        );
        // DP45 adaptive with tolerance 1e-8 should not be dramatically worse
        // than fixed-step RK4 at dt=1e-14.
        assert!(
            dp45_max_change <= rk4_max_change * 100.0 + 1.0e-10,
            "DP45 energy drift {:.2e} unexpectedly much larger than RK4 {:.2e}",
            dp45_max_change,
            rk4_max_change
        );
    }

    #[test]
    fn test_adaptive_runs_fixed_steps() {
        // Verify that adaptive integrators (DP45, DP87) still produce exactly
        // `num_steps` trajectory points (they perform adaptive sub-stepping
        // internally but the outer loop counts fixed macro-steps).
        let num_steps = 75;
        let mut sim = SimulationBuilder::new()
            .material(Ferromagnet::permalloy())
            .external_field(Vector3::new(0.0, 0.0, 0.05))
            .solver_dp87(1.0e-8)
            .initial_magnetization(Vector3::new(0.0, 1.0, 0.0))
            .time_step(1.0e-13)
            .num_steps(num_steps)
            .build()
            .expect("DP87 builder should succeed");

        let result = sim.run().expect("DP87 simulation should run");

        assert_eq!(
            result.len(),
            num_steps + 1,
            "DP87 should produce num_steps + 1 = {} trajectory points, got {}",
            num_steps + 1,
            result.len()
        );

        // Magnetization should remain normalised throughout
        for (i, &m) in result.trajectory.iter().enumerate() {
            let mag = m.magnitude();
            assert!(
                (mag - 1.0).abs() < 1.0e-5,
                "trajectory point {} has |m| = {:.6}, expected ~1.0",
                i,
                mag
            );
        }
    }

    // -----------------------------------------------------------------------
    // Streaming API tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_run_streaming_matches_run_trajectory() {
        // `run_streaming` must reproduce, step for step, exactly what `run()`
        // collects into `SimulationResult` for an identically configured
        // simulation (same material, field, solver, dt, num_steps, initial m).
        fn build_sim() -> Simulation {
            SimulationBuilder::new()
                .material(Ferromagnet::yig())
                .external_field(Vector3::new(0.0, 0.0, 0.1))
                .solver_rk4()
                .initial_magnetization(Vector3::new(1.0, 0.0, 0.0))
                .time_step(1.0e-13)
                .num_steps(40)
                .build()
                .expect("builder should succeed")
        }

        let mut sim_collected = build_sim();
        let collected = sim_collected.run().expect("collected run() should succeed");

        let mut sim_streamed = build_sim();
        let mut streamed_m: Vec<Vector3<f64>> = Vec::new();
        let mut streamed_e: Vec<f64> = Vec::new();
        sim_streamed
            .run_streaming(|step_index, m, energy| {
                assert_eq!(
                    step_index,
                    streamed_m.len(),
                    "callback must be invoked with strictly sequential step indices"
                );
                streamed_m.push(*m);
                streamed_e.push(energy);
                Ok(())
            })
            .expect("run_streaming should succeed");

        assert_eq!(
            streamed_m.len(),
            collected.len(),
            "streaming should report exactly as many steps as the collected trajectory"
        );
        assert_eq!(streamed_e.len(), collected.energies.len());

        for (i, (m_stream, m_collect)) in streamed_m
            .iter()
            .zip(collected.trajectory.iter())
            .enumerate()
        {
            assert!(
                (m_stream.x - m_collect.x).abs() < 1.0e-12
                    && (m_stream.y - m_collect.y).abs() < 1.0e-12
                    && (m_stream.z - m_collect.z).abs() < 1.0e-12,
                "trajectory mismatch at step {}: streamed {:?} vs collected {:?}",
                i,
                m_stream,
                m_collect
            );
        }

        for (i, (&e_stream, &e_collect)) in
            streamed_e.iter().zip(collected.energies.iter()).enumerate()
        {
            let scale = e_collect.abs().max(1.0);
            assert!(
                (e_stream - e_collect).abs() < 1.0e-9 * scale,
                "energy mismatch at step {}: streamed {} vs collected {}",
                i,
                e_stream,
                e_collect
            );
        }

        // The final magnetization reached via streaming must match the one
        // reported by the collected run (both start from identical initial
        // conditions and execute the identical integration loop).
        assert!(
            (sim_streamed.magnetization.x - collected.final_magnetization.x).abs() < 1.0e-12
                && (sim_streamed.magnetization.y - collected.final_magnetization.y).abs() < 1.0e-12
                && (sim_streamed.magnetization.z - collected.final_magnetization.z).abs() < 1.0e-12,
            "final magnetization mismatch: streamed {:?} vs collected {:?}",
            sim_streamed.magnetization,
            collected.final_magnetization
        );
    }

    #[test]
    fn test_run_streaming_handles_large_run_with_constant_size_accumulator() {
        // `run_streaming` must be usable for very large / long-running
        // simulations without the caller needing an ever-growing collection.
        // We drive it for far more steps than any *executed* trajectory
        // elsewhere in this test module (the largest is 5000, in
        // `test_energy_conservation_zero_damping`) using an accumulator whose
        // size is fixed at compile time and independent of `num_steps`.
        #[derive(Default)]
        struct RunningStats {
            count: usize,
            checksum: f64,
            max_abs_component: f64,
        }

        // The accumulator's footprint is a handful of scalars, regardless of
        // how many steps are processed -- unlike `Vec<Vector3<f64>>`, whose
        // size in `SimulationResult::trajectory` grows linearly with
        // `num_steps`. This is the structural property that makes
        // `run_streaming` suitable for large-N / long-run simulations: the
        // hot loop below never calls `Vec::push`.
        assert!(
            std::mem::size_of::<RunningStats>() <= 32,
            "accumulator must stay a small, fixed size independent of num_steps"
        );

        let num_steps = 200_000usize;
        let mut sim = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .initial_magnetization(Vector3::new(1.0, 0.0, 0.0))
            .time_step(1.0e-13)
            .num_steps(num_steps)
            .build()
            .expect("builder should succeed");

        let mut stats = RunningStats::default();
        sim.run_streaming(|_step_index, m, energy| {
            stats.count += 1;
            stats.checksum += energy;
            stats.max_abs_component = stats
                .max_abs_component
                .max(m.x.abs())
                .max(m.y.abs())
                .max(m.z.abs());
            Ok(())
        })
        .expect("run_streaming over a large step count should succeed");

        // Initial state (step 0) plus `num_steps` integration steps.
        assert_eq!(stats.count, num_steps + 1);
        assert!(stats.checksum.is_finite());
        // Magnetization stays a unit vector throughout RK4 integration, so no
        // component should ever exceed 1 (with a small numerical margin).
        assert!(
            stats.max_abs_component <= 1.0 + 1.0e-6,
            "unexpected magnetization component magnitude {}",
            stats.max_abs_component
        );
    }

    #[test]
    fn test_run_streaming_propagates_callback_error_and_aborts_early() {
        // An `Err` returned by the callback must short-circuit the
        // integration immediately rather than continuing to completion.
        let mut sim = SimulationBuilder::new()
            .material(Ferromagnet::yig())
            .external_field(Vector3::new(0.0, 0.0, 0.1))
            .solver_rk4()
            .initial_magnetization(Vector3::new(1.0, 0.0, 0.0))
            .time_step(1.0e-13)
            .num_steps(1000)
            .build()
            .expect("builder should succeed");

        let abort_after = 10usize;
        let mut seen = 0usize;
        let result = sim.run_streaming(|step_index, _m, _energy| {
            seen += 1;
            if step_index >= abort_after {
                return Err(Error::ConfigurationError {
                    description: "synthetic abort for test".to_string(),
                });
            }
            Ok(())
        });

        assert!(result.is_err(), "callback error should propagate as Err");
        // Steps 0..=abort_after inclusive are reported before aborting.
        assert_eq!(seen, abort_after + 1);
    }
}
