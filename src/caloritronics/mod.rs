//! Spin Caloritronics: Unified module for thermoelectric and spin-caloritronic effects.
//!
//! This module combines Onsager transport theory with the spin Seebeck effect,
//! spin Peltier effect, anomalous Nernst effect, and spin Nernst effect into a
//! single unified framework.
//!
//! ## Architecture
//!
//! ```text
//! caloritronics/
//! ├── onsager.rs        – Onsager transport matrix (charge, spin, heat coupling)
//! ├── heat_current.rs   – Heat current calculator (Fourier + Peltier + spin Peltier)
//! └── cross_effects.rs  – Unified SpinCaloritronicsMaterial + CaloritronicsResult
//! ```
//!
//! ## Quick Start
//!
//! ```rust
//! use spintronics::caloritronics::{SpinCaloritronicsMaterial, OnsagerMatrix};
//! use spintronics::Vector3;
//!
//! // Create a YIG/Pt material at 300 K
//! let mat = SpinCaloritronicsMaterial::yig_pt(300.0);
//!
//! // Apply a temperature gradient of 1000 K/m along x
//! let grad_t = Vector3::new(1000.0, 0.0, 0.0);
//! let j_spin  = Vector3::zero();
//!
//! // Compute all caloritronic observables
//! let result = mat.compute_all(&grad_t, &j_spin).expect("computation failed");
//! println!("SSE current: {:?}", result.spin_seebeck_current);
//! println!("Nernst voltage: {} V/m", result.nernst_voltage);
//! ```
//!
//! ## Physical Background
//!
//! Spin caloritronics studies the interplay between heat currents and spin currents,
//! analogous to classical spintronics (charge + spin) and thermoelectrics (charge + heat).
//! The Onsager reciprocal relations provide the unifying theoretical framework:
//!
//! - **Spin Seebeck effect** (SSE): ∇T → spin current; coefficient S_s
//! - **Spin Peltier effect** (SPE): j_s → heat current; coefficient Π_s = T·S_s
//! - **Anomalous Nernst effect** (ANE): ∇T → transverse charge voltage; via Hall angle θ_H
//! - **Spin Nernst effect** (SNE): ∇T → transverse spin current; via spin Nernst angle
//!
//! ## References
//!
//! - G. E. W. Bauer, E. Saitoh, B. J. van Wees, "Spin caloritronics",
//!   *Nat. Mater.* **11**, 391–399 (2012)
//! - K. Uchida et al., "Observation of longitudinal spin-Seebeck effect in magnetic
//!   insulators", *Appl. Phys. Lett.* **97**, 172505 (2010)
//! - J. Flipse et al., "Observation of the spin Peltier effect for magnetic insulators",
//!   *Phys. Rev. Lett.* **113**, 027601 (2014)

pub mod cross_effects;
pub mod heat_current;
pub mod onsager;

pub use cross_effects::{CaloritronicsResult, SpinCaloritronicsMaterial};
pub use heat_current::HeatCurrentCalculator;
pub use onsager::{AllCurrents, OnsagerMatrix};
