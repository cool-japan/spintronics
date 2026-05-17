//! Mathematical primitives for spintronics quantum calculations.
//!
//! This module provides fundamental mathematical types used across the library's
//! quantum, topological, and transport physics modules:
//!
//! - [`Complex`] — a lightweight complex number (`re + i·im`) with constants `ZERO`, `ONE`, `I`
//!   and operations `add`, `sub`, `mul`, `div`, `scale`, `exp`, `pow_n`, `conj`, `neg`.
//! - [`CMatrix`] — a small dense N×N complex matrix with Gauss-Jordan inversion and
//!   Jacobi Hermitian eigendecomposition (intended for N ≤ 64).
//!
//! # Design Philosophy
//!
//! These types are purpose-built for physics code in this library. They use explicit
//! method-call syntax (`a.mul(&b)`) rather than operator overloading, to make complex
//! arithmetic in physics equations unambiguous and easy to read alongside real-valued
//! operations.
//!
//! No external numerical linear algebra dependencies are needed — the implementations
//! are self-contained and correct for the small matrix sizes encountered in:
//!
//! - Magnon band Hamiltonians (2×2 honeycomb, 3×3 kagome, up to 60×60 edge strips)
//! - NEGF transport (1D tight-binding ≤ 32 sites → 32×32 matrices)
//! - Bogoliubov / Holstein-Primakoff transformations (analytically small)
//!
//! # Example
//!
//! ```rust
//! use spintronics::math::{Complex, CMatrix};
//!
//! // Complex arithmetic
//! let z = Complex::from_polar(1.0, std::f64::consts::FRAC_PI_2);
//! assert!((z.im - 1.0).abs() < 1e-14);
//!
//! // Small matrix inversion
//! let a = CMatrix::from_rows(vec![
//!     vec![Complex::new(2.0, 0.0), Complex::ZERO],
//!     vec![Complex::ZERO, Complex::new(3.0, 0.0)],
//! ]).unwrap();
//! let inv = a.inverse().unwrap();
//! assert!((inv.get(0, 0).re - 0.5).abs() < 1e-12);
//! ```

pub mod complex;
pub mod matrix;

pub use complex::Complex;
pub use matrix::CMatrix;
