//! Physical Validity Checks and Experimental Validation
//!
//! Two-pronged validation:
//! 1. **Parameter checks** — runtime assertions that physical parameters
//!    satisfy basic consistency requirements (positive temperatures, normalized
//!    magnetization, damping ∈ [0,1), etc.). Re-exported at module root for
//!    backward compatibility.
//! 2. **Experimental validation** — comparison of simulation results against
//!    landmark experimental papers, with reference data tables and quantitative
//!    relative-error metrics.
//!
//! ## Examples
//!
//! ```
//! use spintronics::validation::{check_damping, check_magnetization};
//!
//! assert!(check_damping(0.01).is_ok());
//! assert!(check_magnetization(1.4e5).is_ok());
//! ```
//!
//! ## Module Layout
//!
//! - [`parameter_checks`] — `check_*` functions for runtime parameter validation
//!   (re-exported at the module root for backward compatibility).
//! - [`experimental`] — landmark-paper validations:
//!   - [`experimental::demidov_2006`] — Damon-Eshbach BLS in YIG films
//!     (Demidov *et al.*, PRL **96**, 097202).
//!   - [`experimental::saitoh_2006`] — Inverse Spin Hall Effect in Pt/Py
//!     (Saitoh *et al.*, APL **88**, 182509).
//!   - [`experimental::uchida_2008`] — Longitudinal Spin Seebeck Effect in Pt/YIG
//!     (Uchida *et al.*, Nature **455**, 778).
//!   - [`experimental::mosendz_2010`] — Quantitative spin pumping and refined Pt
//!     spin Hall angle in Py/Pt (Mosendz *et al.*, PRL **104**, 046601).
//!   - [`experimental::liu_2012`] — SOT switching with β-Ta in Ta/CoFeB/MgO
//!     (Liu *et al.*, Science **336**, 555).

pub mod experimental;
pub mod parameter_checks;
pub mod standard_problems;

// Backward-compatibility re-export: every `check_*` (and any other public item)
// from `parameter_checks` is reachable directly from `crate::validation::*`.
pub use parameter_checks::*;
