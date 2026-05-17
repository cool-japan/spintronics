//! # Reverse-mode Automatic Differentiation (v0.6.0)
//!
//! This module provides a **tape-based reverse-mode AD engine** (also called
//! *back-propagation* or the *adjoint method*) together with a suite of
//! gradient-based optimizers and differentiable spintronics physics functions.
//!
//! ## How it works
//!
//! Every arithmetic or transcendental operation on a [`Var`] node records
//! itself on a shared [`Tape`] (the *Wengert list*).  After the forward
//! evaluation of a scalar loss `z`, calling [`Tape::backward`] propagates
//! partial derivatives from `z` back through the recorded operations in
//! reverse order, accumulating into the `.grad()` field of each leaf variable.
//!
//! This is the same algorithm that powers modern deep-learning frameworks
//! (PyTorch `autograd`, JAX, etc.), here implemented from scratch in safe Rust
//! without any external dependencies.
//!
//! ## Quick-start
//!
//! ```rust
//! use spintronics::autodiff::{Tape, Var};
//!
//! let tape = Tape::new();
//! let x = Var::leaf(&tape, 3.0);
//! let y = Var::leaf(&tape, 4.0);
//! let z = (x * y) + x.sin();   // z = x·y + sin(x)
//! tape.backward(z);
//! // dz/dx = y + cos(x), dz/dy = x
//! let dzdx = x.grad();
//! let dzdy = y.grad();
//! assert!((dzdx - (4.0 + 3.0_f64.cos())).abs() < 1e-12);
//! assert!((dzdy - 3.0).abs() < 1e-12);
//! ```
//!
//! ## Optimizers
//!
//! | Type | Description |
//! |------|-------------|
//! | [`Sgd`]   | Stochastic Gradient Descent (with optional momentum) |
//! | [`Adam`]  | Adaptive Moment Estimation — Kingma & Ba (2015) |
//! | [`LBfgs`] | Limited-memory BFGS — Nocedal (1980) |
//!
//! [`ParameterFitter`] wraps any optimizer and exposes a single [`ParameterFitter::fit`]
//! method for end-to-end physics model fitting.
//!
//! ## Differentiable physics functions
//!
//! | Function | Physics |
//! |----------|---------|
//! | [`kittel_frequency_diff`]  | FMR resonance frequency |
//! | [`zeeman_energy_diff`]     | Zeeman coupling |
//! | [`exchange_energy_diff`]   | Heisenberg nearest-neighbour exchange |
//! | [`dmi_energy_diff`]        | Interfacial DMI |
//! | [`anisotropy_energy_diff`] | Uniaxial anisotropy |
//!
//! ## References
//!
//! - A. Griewank & A. Walther, *Evaluating Derivatives: Principles and
//!   Techniques of Algorithmic Differentiation*, 2nd ed., SIAM (2008).
//! - A. G. Baydin, B. A. Pearlmutter, A. A. Radul & J. M. Siskind,
//!   "Automatic Differentiation in Machine Learning: a Survey",
//!   *J. Mach. Learn. Res.* **18**, 1–43 (2018).
//! - S. Linnainmaa, "Taylor expansion of the accumulated rounding error",
//!   *BIT Numer. Math.* **16**, 146–160 (1976).
//! - D. P. Kingma & J. Ba, "Adam: A Method for Stochastic Optimization",
//!   *ICLR* (2015), arXiv:1412.6980.
//! - J. Nocedal, "Updating Quasi-Newton Matrices with Limited Storage",
//!   *Math. Comp.* **35**, 773–782 (1980).

pub mod optimizer;
pub mod physics_fns;
pub mod tape;

// ─── Re-exports ───────────────────────────────────────────────────────────────

pub use optimizer::{Adam, FitResult, LBfgs, Optimizer, OptimizerKind, ParameterFitter, Sgd};
pub use physics_fns::{
    anisotropy_energy_diff, dmi_energy_diff, exchange_energy_diff, kittel_frequency_diff,
    llg_torque_norm_diff, zeeman_energy_diff,
};
pub use tape::{check_gradient, finite_diff_grad, Tape, Var};
