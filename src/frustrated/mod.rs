//! Frustrated magnets module
//!
//! This module implements physics of geometrically frustrated magnetic systems,
//! where the lattice geometry prevents simultaneous minimization of all
//! pairwise interactions.
//!
//! # Submodules
//!
//! - [`lattice`] - Frustrated lattice geometries (triangular, kagome, pyrochlore)
//! - [`spin_ice`] - Spin ice physics (ice rules, Pauling entropy, monopoles)
//! - [`kagome`] - Kagome lattice specific physics (magnon bands, spin liquid)
//!
//! # Overview
//!
//! Geometric frustration is a fundamental concept in condensed matter physics.
//! When antiferromagnetic interactions exist on lattices with triangular motifs,
//! the system cannot find a unique ground state, leading to macroscopic
//! degeneracy and exotic phenomena:
//!
//! - **Spin ice**: Residual entropy and emergent magnetic monopoles
//! - **Spin liquids**: No long-range order down to T = 0
//! - **Flat magnon bands**: Localized excitations from frustration
//! - **Dirac magnons**: Topological band crossings
//!
//! # Quick Start
//!
//! ```rust
//! use spintronics::frustrated::lattice::{FrustratedLattice, frustration_parameter};
//! use spintronics::frustrated::spin_ice::{pauling_entropy, SpinIceParams, SpinIce};
//! use spintronics::frustrated::kagome::{kagome_magnon_bands, KagomeMagnonConfig};
//!
//! // Create a triangular antiferromagnet
//! let mut lattice = FrustratedLattice::triangular(6, 6, 1e-21, 3e-10)
//!     .expect("lattice creation failed");
//! lattice.set_120_degree_order();
//!
//! // Calculate frustration parameter for a known material
//! let f = frustration_parameter(-500.0, 3.5).expect("valid parameters");
//! assert!(f > 100.0); // strongly frustrated
//!
//! // Pauling entropy of spin ice
//! let s_p = pauling_entropy();
//! assert!(s_p > 0.0);
//!
//! // Kagome magnon flat band
//! let config = KagomeMagnonConfig::default();
//! let bands = kagome_magnon_bands(1.0, 0.0, &config);
//! assert!(bands.bands[0] < 0.3); // flat band near zero
//! ```

pub mod kagome;
pub mod lattice;
pub mod spin_ice;

// Re-exports for convenience
pub use kagome::{
    dirac_magnon_energy, is_spin_liquid_candidate, kagome_band_path, kagome_magnon_bands,
    neel_120_order_parameter, spin_correlation_function, KagomeMagnonBands, KagomeMagnonConfig,
};
pub use lattice::{
    curie_weiss_temperature, frustration_parameter, FrustratedLattice, LatticeType, Xorshift64,
};
pub use spin_ice::{
    monopole_coulomb_interaction, monopole_creation_energy, monopole_density_factor,
    pauling_entropy, pauling_entropy_per_spin_kb, SpinIce, SpinIceParams, Tetrahedron,
};
