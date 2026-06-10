//! Multiferroic Materials and Magnetoelectric Coupling
//!
//! This module implements the physics of multiferroic materials, which simultaneously
//! exhibit two or more ferroic orders (ferroelectricity, ferromagnetism, and/or
//! ferroelasticity), and the magnetoelectric (ME) coupling between them.
//!
//! ## Overview
//!
//! Multiferroics represent one of the most actively studied frontiers in condensed
//! matter physics due to their potential for electrically controlled magnetic memory
//! and magnetically tunable ferroelectric devices.  The simultaneous ordering of
//! electric and magnetic degrees of freedom, and their mutual coupling, enables
//! unprecedented cross-control of material properties.
//!
//! ## Classification of Multiferroics
//!
//! ### Type-I Multiferroics
//! Ferroelectricity and magnetism arise from **different** microscopic mechanisms
//! and are largely independent order parameters.  The strong spontaneous polarization
//! (P_s ~ μC/cm²) makes them technologically attractive, but the weak ME coupling
//! limits cross-control.
//!
//! **Canonical example:** BiFeO₃ (BFO)
//! - Ferroelectric from Bi³⁺ lone-pair (6s²) stereochemical activity
//! - G-type antiferromagnetism from Fe³⁺ super-exchange (T_N = 643 K)
//! - Weak linear ME coupling due to broken inversion + time-reversal
//!
//! ### Type-II Multiferroics (Spin-induced ferroelectrics)
//! Ferroelectricity is **caused by** magnetic ordering — typically a non-collinear
//! or spiral spin structure that breaks inversion symmetry through the spin-current
//! (Dzyaloshinskii-Moriya) or exchange-striction mechanism.  The coupling is
//! intrinsically strong but P_s is small (nC/cm²).
//!
//! **Canonical example:** TbMnO₃
//! - Cycloidal Mn³⁺ spin spiral induces P via the KNB mechanism
//! - T_FE = 27 K, T_N = 41 K
//!
//! ### Type-III Multiferroics
//! Lone-pair polarization combined with spiral magnetic order; intermediate
//! behavior between Type-I and Type-II.
//!
//! ## Magnetoelectric Effect
//!
//! The linear magnetoelectric effect is described by the coupling tensor α_ij:
//!
//! ```text
//! P_i = ε₀ χ_e E_i + α_ij H_j     (electric polarization from magnetic field)
//! M_i = μ₀ χ_m H_i + α_ij E_j     (magnetization from electric field, via α^T)
//! ```
//!
//! The thermodynamic bound on the ME coupling (Brown-Hornreich-Shtrikman):
//!
//! ```text
//! α_ij² < ε₀ χ_e,ii · μ₀ χ_m,jj
//! ```
//!
//! ## Spin-Current Mechanism (KNB)
//!
//! For cycloidal spin spirals the Katsura-Nagaosa-Balatsky (KNB) mechanism
//! gives a microscopic expression for the induced polarization:
//!
//! ```text
//! P ∝ e_ij × (S_i × S_j)
//! ```
//!
//! where e_ij is the unit vector connecting sites i and j.  This vanishes for
//! helical spirals (q ∥ rotation_axis) and is maximal for cycloidal spirals.
//!
//! ## Key References
//!
//! - S.-W. Cheong & M. Mostovoy, "Multiferroics: a magnetic twist for ferroelectricity",
//!   Nat. Mater. **6**, 13–20 (2007)
//! - M. Mostovoy, "Ferroelectricity in Spiral Magnets",
//!   Phys. Rev. Lett. **96**, 067601 (2006)
//! - H. Katsura, N. Nagaosa & A. V. Balatsky, "Spin Current and Magnetoelectric Effect
//!   in Noncollinear Magnets", Phys. Rev. Lett. **95**, 057205 (2005)
//! - N. A. Hill, "Why Are There so Few Magnetic Ferroelectrics?",
//!   J. Phys. Chem. B **104**, 6694–6709 (2000)
//! - N. A. Spaldin & M. Fiebig, "The Renaissance of Magnetoelectric Multiferroics",
//!   Science **309**, 391–392 (2005)
//!
//! ## Example
//!
//! ```rust
//! use spintronics::multiferroic::{
//!     MagnetoelectricTensor, KnbMechanism, dzyaloshinskii_moriya_polarization,
//! };
//! use spintronics::vector3::Vector3;
//!
//! // BiFeO3 — prototypical Type-I multiferroic
//! let bfo = MagnetoelectricTensor::bife_o3();
//! assert!(bfo.is_linear_magnetoelectric());
//!
//! // Cr2O3 — prototypical linear magnetoelectric
//! let cr2o3 = MagnetoelectricTensor::cr2_o3();
//! let h = Vector3::new(0.0, 0.0, 1e6); // H = 1 MA/m along z
//! let p = cr2o3.electric_polarization_from_field(&h);
//! assert!(p.magnitude() > 0.0);
//!
//! // KNB mechanism for TbMnO3-like cycloidal spiral
//! let knb = KnbMechanism::tbmno3();
//! let s_i = Vector3::new(1.0, 0.0, 0.0);
//! let s_j = Vector3::new(0.0, 1.0, 0.0);
//! let e_ij = Vector3::new(1.0, 0.0, 0.0);
//! let p_knb = knb.polarization_from_spiral(&s_i, &s_j, &e_ij);
//! assert!(p_knb.magnitude() > 0.0);
//!
//! // DM-based spin-current polarization (free function).
//! // Use a bond vector NOT parallel to (s_i × s_j) = (0,0,1) so that
//! // P_DM = ê_ij × (s_i × s_j) is non-zero.
//! let bond = Vector3::new(1.0, 0.0, 0.0);
//! let p_dm = dzyaloshinskii_moriya_polarization(&s_i, &s_j, &bond);
//! assert!(p_dm.magnitude() > 0.0);
//! ```

pub mod magnetoelectric;
pub mod spin_polarization;

pub use magnetoelectric::{
    dzyaloshinskii_moriya_polarization, exchange_striction_polarization, toroidal_moment,
    MagnetoelectricTensor, MultiferroicType,
};
pub use spin_polarization::{
    magnon_drag_contribution, spin_current_from_spiral, InverseMagnetoelectric, KnbMechanism,
};
