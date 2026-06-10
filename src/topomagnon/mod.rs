//! Topological magnon bands — bosonic analogues of electronic Chern insulators.
//!
//! This module implements the theory of topological magnons: spin-wave excitations
//! in magnetic insulators that carry non-trivial Berry curvature and topological
//! Chern numbers, leading to the magnon Hall effect and robust edge states.
//!
//! # Theoretical Background
//!
//! In ferromagnetic insulators the low-energy excitations are magnons (quantised
//! spin waves). In materials with broken time-reversal symmetry — either by an
//! external field or by the Dzyaloshinskii-Moriya interaction (DMI) — the magnon
//! bands acquire Berry curvature and, when integrated over the Brillouin zone,
//! a topological Chern number. This is the bosonic analogue of the quantum
//! anomalous Hall effect first proposed for electrons by Haldane.
//!
//! ## Key References
//!
//! - F. D. M. Haldane, "Model for a quantum Hall effect without Landau levels:
//!   Condensed-matter realization of the 'parity anomaly'",
//!   *Phys. Rev. Lett.* **61**, 2015 (1988).
//!   The original Haldane model on the honeycomb lattice with imaginary next-nearest-
//!   neighbour hopping; the bosonic analogue drives topological magnon physics.
//!
//! - T. Fukui, Y. Hatsugai, H. Suzuki, "Chern numbers in discretized Brillouin zone:
//!   Efficient method of computing (spin) Hall conductances",
//!   *J. Phys. Soc. Jpn.* **74**, 1674 (2005).
//!   The lattice Chern-number algorithm used in [`chern_number`].
//!
//! - R. Matsumoto, S. Murakami, "Theoretical prediction of a rotating magnon wave
//!   packet in ferromagnets",
//!   *Phys. Rev. Lett.* **106**, 197202 (2011).
//!   Thermal (magnon) Hall conductivity formula implemented in [`magnon_hall`].
//!
//! - Y. Onose, T. Ideue, H. Katsura, Y. Shiomi, N. Nagaosa, Y. Tokura,
//!   "Observation of the magnon Hall effect",
//!   *Science* **329**, 297 (2010).
//!   First experimental observation of magnon Hall transport in Lu₂V₂O₇.
//!
//! - S. A. Owerre, "A first theoretical realization of honeycomb topological magnon
//!   insulator",
//!   *J. Phys.: Condens. Matter* **28**, 386001 (2016).
//!   Topological magnons on the honeycomb lattice with DMI-driven Haldane mass term.
//!
//! # Module Structure
//!
//! | Sub-module | Contents |
//! |---|---|
//! | [`band_model`] | Bloch Hamiltonians for honeycomb, kagome, and square-DMI magnon lattices |
//! | [`berry_curvature`] | Berry curvature Ω_n(k) via the sum-over-states formula |
//! | [`chern_number`] | Lattice Chern numbers via Fukui-Hatsugai and Wilson-loop methods |
//! | [`edge_modes`] | Strip-geometry Hamiltonians and topological edge modes |
//! | [`magnon_hall`] | Matsumoto-Murakami thermal Hall conductivity |
//!
//! # Quick Start
//!
//! ```rust
//! use spintronics::topomagnon::{MagnonBandModel, ChernNumber};
//!
//! // Build a honeycomb Haldane-magnon model with DMI
//! let model = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.2, 0.1).unwrap();
//!
//! // Compute the Chern number of the lower band on a 20×20 k-grid
//! let chern = ChernNumber::new(&model);
//! let c = chern.compute(0, 20, 20).unwrap();
//! assert_eq!(c.abs(), 1);
//! ```

pub mod axion;
pub mod band_model;
pub mod band_model_3d;
pub mod berry_curvature;
pub mod chern_number;
pub mod edge_modes;
pub mod hoti;
pub mod magnon_hall;
pub mod qsh;
pub mod wilson;

pub use axion::{AxionElectrodynamics, AxionMagnonPhoton};
pub use band_model::{LatticeType, MagnonBandModel};
pub use band_model_3d::{LatticeType3D, MagnonBandModel3D};
pub use berry_curvature::BerryCurvature;
pub use chern_number::ChernNumber;
pub use edge_modes::{EdgeMode, EdgeModes, EdgeSide};
pub use hoti::{BbhModel, BreathingKagomeModel, CornerStateSolver, HotiLattice};
pub use magnon_hall::MagnonHallConductivity;
pub use qsh::KaneMeleModel;
pub use wilson::WilsonLoop;
