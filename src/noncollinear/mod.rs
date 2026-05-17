//! Noncollinear magnetic order — spiral structures and the Luttinger-Tisza method
//!
//! This module implements the physics of noncollinear magnetic ordering, where
//! the spin directions vary periodically in space without being simply parallel
//! or antiparallel as in collinear ferro- and antiferromagnets. Such structures
//! arise from competing exchange interactions (frustration), Dzyaloshinskii-Moriya
//! interaction (DMI), single-ion anisotropy, and long-range dipolar forces.
//!
//! # Physics Background
//!
//! ## Spin Spirals and Helical Structures
//!
//! A spin spiral is characterized by a propagation (ordering) wavevector **q** such
//! that the spin at position **r** rotates as the observer moves along **q**. The
//! two canonical geometries are:
//!
//! - **Cycloidal spiral**: the rotation plane contains **q**. Found in orthorhombic
//!   manganites such as TbMnO₃, which exhibit strong multiferroic coupling.
//! - **Helical spiral**: the rotation plane is perpendicular to **q** (screw-like).
//!   Found in MnSi, FeGe, and Cu₂OSeO₃.
//! - **Conical spiral**: a cone of spins winding about a net ferromagnetic axis.
//! - **Fan structure**: a partial cone that fans symmetrically about an easy axis.
//!
//! ## Multiferroicity via the KNB Mechanism
//!
//! Katsura, Nagaosa, and Balatsky (2005) showed that a cycloidal spin spiral
//! generates an electric polarization
//!
//! **P** ~ **e**_{ij} × (**S**_i × **S**_j)
//!
//! where **e**_{ij} is the bond direction. Equivalently, for a spiral with
//! propagation vector **q** and rotation axis **ê**:
//!
//! **P** ∝ **q** × **ê**
//!
//! This mechanism explains the large magnetically-induced polarization in TbMnO₃,
//! GdMnO₃, and related frustrated magnets (Mostovoy 2006).
//!
//! ## Luttinger-Tisza Method
//!
//! The Luttinger-Tisza (LT) method (1946) finds the **classical** ground state
//! of the Heisenberg model H = Σ_{ij} J_{ij} **S**_i · **S**_j by working in
//! Fourier space. For a Bravais lattice, the exchange matrix in reciprocal space is
//!
//! J(**q**) = Σ_{δ} J_{0δ} exp(i **q** · **δ**)
//!
//! The ground-state ordering wavevector **q*** minimizes J(**q***), which for
//! real (centrosymmetric) systems reduces to minimizing
//!
//! J(**q**) = Σ_{δ} J_{0δ} cos(**q** · **δ**)
//!
//! The method correctly predicts spiral order in frustrated J₁-J₂ chains and
//! many 2D and 3D frustrated magnets.
//!
//! # Submodules
//!
//! - [`spiral`] — Spin-spiral structures (cycloidal, helical, conical, fan)
//! - [`luttinger_tisza`] — LT classical ground-state search for Heisenberg models
//!
//! # Key References
//!
//! - A. Yoshimori, "A New Type of Antiferromagnetic Structure in the Rutile Type Crystal",
//!   *J. Phys. Soc. Jpn.* **14**, 807–821 (1959) — first prediction of helical spin structures
//! - J. Villain, "La structure des substances magnétiques",
//!   *J. Phys. Chem. Solids* **11**, 303–309 (1959) — frustrated spin systems and spirals
//! - J. M. Luttinger and L. Tisza, "Theory of Dipole Interaction in Crystals",
//!   *Phys. Rev.* **70**, 954–964 (1946) — Fourier method for classical ground states
//! - M. Mostovoy, "Ferroelectricity in Spiral Magnets",
//!   *Phys. Rev. Lett.* **96**, 067601 (2006) — cycloidal magnetism and ferroelectricity
//!
//! # Quick Start
//!
//! ```rust
//! use spintronics::noncollinear::{SpinSpiral, LuttingerTisza, ExchangeInteraction};
//! use spintronics::vector3::Vector3;
//! use std::f64::consts::PI;
//!
//! // TbMnO3 cycloidal spiral (multiferroic material)
//! let tbmno3 = SpinSpiral::terbium_manganese_oxide();
//! let p_dir = tbmno3.electric_polarization_direction();
//! assert!(p_dir.is_some()); // cycloidal → finite polarization
//!
//! // J1-J2 chain: frustrated → spiral ground state
//! let a = 3.0e-10_f64;
//! let lt = LuttingerTisza::j1j2_chain(1.0, 0.5, a, 0.5);
//! assert!(lt.is_spiral_ground_state(1e-3));
//! ```

pub mod luttinger_tisza;
pub mod spiral;

pub use luttinger_tisza::{ExchangeInteraction, LuttingerTisza};
pub use spiral::{SpinSpiral, SpiralChirality, SpiralType};
