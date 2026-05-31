//! NIST muMAG standard micromagnetic benchmark problems.
//!
//! These problems provide well-defined, community-accepted test cases for
//! micromagnetics codes, enabling quantitative comparison between different
//! implementations (OOMMF, mumax3, Fidimag, this library, etc.).
//!
//! ## Standard Problems
//!
//! | Problem | Physics                                          | Status          |
//! |---------|--------------------------------------------------|-----------------|
//! | SP#2    | Hysteresis loop in a thin film                   | Not implemented |
//! | SP#3    | Single-domain limit: cube flower↔vortex crossover | Implemented     |
//! | SP#4    | Dynamic switching in a thin film strip            | Not implemented |
//!
//! ## Overview of SP#3
//!
//! SP#3 asks: for a Permalloy cube of edge length L, which equilibrium remanence
//! state is lower in energy — the near-uniform "flower" state or the
//! flux-closure "vortex" state?
//!
//! A competition between exchange energy (which favours uniform alignment) and
//! magnetostatic energy (which favours flux closure) determines a critical length
//! scale L_c. Below L_c the flower state wins; above L_c the vortex wins.
//!
//! The muMAG reference value is L_c ≈ 8.47 l_ex where
//! `l_ex = sqrt(A / (0.5 μ₀ Ms²))` is the exchange length (~5.7 nm for Permalloy).
//!
//! # References
//! - muMAG Standard Problems: <https://www.ctcms.nist.gov/~rdm/mumag.org.html>
//! - R. D. McMichael and M. J. Donahue, IEEE Trans. Magn. 33, 4167 (1997)
//! - M. J. Donahue and D. G. Porter, OOMMF User's Guide, NIST Interagency Report
//!   6376, National Institute of Standards and Technology, Gaithersburg, MD (1999)

pub mod sp3;

pub use sp3::{Sp3Config, Sp3Result, StableState, StandardProblem3};
