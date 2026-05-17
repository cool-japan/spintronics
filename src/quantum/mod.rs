//! Quantum magnonics — Holstein-Primakoff, Bogoliubov, and zero-point fluctuations.
//!
//! This module implements the key quantum-mechanical transformations used in the
//! magnon physics of ferro- and antiferromagnets: the Holstein-Primakoff mapping,
//! the Bogoliubov canonical transformation for magnon Hamiltonians, and zero-point
//! quantum fluctuation amplitudes for confined spin-wave modes.
//!
//! # Key References
//!
//! - T. Holstein, H. Primakoff, "Field Dependence of the Intrinsic Domain
//!   Magnetization of a Ferromagnet", *Phys. Rev.* **58**, 1098 (1940).
//!   The original HP transformation mapping spin operators to bosonic creation/
//!   annihilation operators.
//!
//! - N. N. Bogoliubov, "On the Theory of Superfluidity", *J. Phys. USSR* **11**, 23
//!   (1947). The canonical Bogoliubov transformation diagonalizing quadratic Bose
//!   Hamiltonians via a squeeze transformation (cosh/sinh rotation in Fock space).
//!
//! - A. Kamra, W. Belzig, "Super-Poissonian Shot Noise of Squeezed-Magnon Mediated
//!   Spin Transport", *Phys. Rev. Lett.* **116**, 146601 (2016). Demonstrates that
//!   the Bogoliubov vacuum of an antiferromagnet is a *squeezed* magnon state with
//!   non-trivial quantum correlations, enhancing spin noise beyond the Poisson limit.
//!
//! # Module Organization
//!
//! - [`holstein_primakoff`]: Mapping of spin-`S` operators to bosonic operators in
//!   powers of 1/S. Handles both the linear (leading-order) approximation and the
//!   quadratic 1/S correction.
//! - [`bogoliubov`]: Bogoliubov (squeeze) transformation for diagonalizing bilinear
//!   magnon Hamiltonians that include pairing terms. Includes a convenience
//!   constructor for the AFM square lattice.
//! - [`zero_point`]: Zero-point fluctuation amplitudes and Casimir-like energies for
//!   quantized spin-wave modes in confined nanostructures.
//!
//! # Example
//!
//! ```rust
//! use spintronics::quantum::{HolsteinPrimakoff, HpOrder, BogoliubovTransform};
//!
//! // Build Holstein-Primakoff for spin-1 ferromagnet
//! let hp = HolsteinPrimakoff::ferromagnet_default();
//! assert!((hp.transform_sz(0.0) - 1.0).abs() < 1e-12);
//!
//! // Diagonalize a simple pairing Hamiltonian A_k = 10, B_k = 3
//! let a = vec![10.0];
//! let b = vec![3.0];
//! let bogo = BogoliubovTransform::from_hamiltonian(&a, &b)
//!     .expect("stable Hamiltonian");
//! let omega = bogo.magnon_frequency(0).expect("valid index");
//! let expected = (10.0_f64 * 10.0 - 3.0 * 3.0_f64).sqrt();
//! assert!((omega - expected).abs() < 1e-10);
//! ```

pub mod bogoliubov;
pub mod holstein_primakoff;
pub mod zero_point;

pub use bogoliubov::BogoliubovTransform;
pub use holstein_primakoff::{HolsteinPrimakoff, HpOrder};
pub use zero_point::ZeroPointFluctuations;
