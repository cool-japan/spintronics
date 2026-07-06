//! Orbitronics Module
//!
//! This module implements orbital angular momentum (OAM) transport phenomena,
//! a rapidly growing field in condensed matter physics that studies how orbital
//! currents can be generated, manipulated, and converted to achieve functionality
//! beyond conventional spin-based spintronics.
//!
//! ## Overview
//!
//! Orbitronics exploits the orbital angular momentum of electrons (as opposed to
//! their spin) for information processing and magnetization control. The key
//! insight is that orbital currents can be generated even in light metals with
//! weak spin-orbit coupling (SOC), through the Orbital Hall Effect (OHE).
//!
//! ## Key Physics
//!
//! ### Orbital Hall Effect (OHE)
//! Analogous to the Spin Hall Effect, the OHE generates transverse orbital
//! angular momentum current from longitudinal charge current. Crucially, the
//! OHE does not require strong SOC and can be giant in light 3d transition
//! metals like Cr and Ti.
//!
//! ### Orbital-to-Spin Conversion
//! At interfaces with materials possessing strong SOC (e.g., Pt, Bi), orbital
//! currents can be efficiently converted to spin currents, enabling
//! magnetization manipulation.
//!
//! ### Orbital Torques
//! Orbital currents can exert torques on ferromagnetic layers, providing an
//! alternative pathway to spin-orbit torques for current-driven magnetization
//! switching.
//!
//! ### Orbital Rashba Effect
//! At interfaces lacking inversion symmetry, orbital angular momentum splits
//! in momentum space, leading to the orbital Edelstein effect.
//!
//! ### Local d-orbital moments
//! A separate, loosely-coupled sub-theme: local (atomic-like) d-orbital
//! magnetic moments of 3d transition-metal ions in a static ligand field —
//! crystal-field splitting, Hund's-rule high-spin/low-spin states, orbital
//! quenching, and effective moments. See [`crystal_field`] and
//! [`d_orbital_moment`].
//!
//! ## Modules
//!
//! - [`orbital_hall`]: Orbital Hall Effect, material database, orbital-to-spin conversion
//! - [`orbital_torque`]: Orbital torques, orbital Rashba effect, switching efficiency
//! - [`crystal_field`]: Real d-orbital basis, orbital angular momentum operators,
//!   crystal-field and spin-orbit Hamiltonians
//! - [`d_orbital_moment`]: Hund's-rule free-ion terms, high-spin/low-spin,
//!   orbital quenching, effective moments, preset ions
//!
//! ## References
//!
//! - H. Kontani et al., "Giant Orbital Hall Effect in Transition Metals",
//!   Phys. Rev. Lett. 102, 016601 (2009)
//! - D. Go et al., "Intrinsic Spin and Orbital Hall Effects from Orbital Texture",
//!   Phys. Rev. Lett. 121, 086602 (2018)
//! - D. Go and H.-W. Lee, "Orbital torque: Torque generation by orbital
//!   current injection", Phys. Rev. Research 2, 013177 (2020)
//! - S. Lee et al., "Efficient conversion of orbital Hall current to spin
//!   current for spin-orbit torque switching", Commun. Phys. 4, 234 (2021)
//! - C. J. Ballhausen, "Introduction to Ligand Field Theory" (McGraw-Hill, 1962)
//! - J. S. Griffith, "The Theory of Transition-Metal Ions" (Cambridge, 1961)

pub mod crystal_field;
pub mod d_orbital_moment;
pub mod orbital_hall;
pub mod orbital_torque;

// Re-exports for convenient access
pub use crystal_field::{
    crystal_field_hamiltonian, orbital_angular_momentum_operators, orbital_l_squared,
    soc_hamiltonian, CrystalFieldEnvironment, DOrbital, OrbitalSymmetry,
};
pub use d_orbital_moment::{free_ion_term, CrystalFieldModel, FreeIonTerm, GroundTermSymmetry};
pub use orbital_hall::{
    cr_pt_system, cu_pt_system, material_ohe_ranking, sigma_oh_to_si, ti_pt_system,
    OrbitalHallEffect, OrbitalHallMaterial, OrbitalToSpinConverter,
};
pub use orbital_torque::{cr_pt_cofeb_system, ti_pt_cofeb_system, OrbitalRashba, OrbitalTorque};
