//! Non-equilibrium Green's functions (NEGF) for mesoscopic spintronics.
//!
//! This module implements the Keldysh–Landauer formalism for coherent quantum
//! transport through 1D tight-binding chains coupled to fermionic reservoirs.
//! It covers:
//!
//! - **Green's functions** (retarded, advanced, lesser, greater) for finite chains
//!   with wide-band-limit self-energies.
//! - **Landauer transmission** T(E) and Büttiker current I(V).
//! - **Keldysh non-equilibrium occupation** densities and per-site carrier filling.
//! - **Shot noise** and Fano factors in the Blanter–Büttiker framework.
//! - **Spin accumulation** and 1D spin diffusion with explicit and implicit solvers.
//!
//! # Physical References
//!
//! - L. V. Keldysh, Sov. Phys. JETP **20**, 1018 (1965) — non-equilibrium formalism
//! - S. Datta, *Electronic Transport in Mesoscopic Systems* (Cambridge, 1995) — NEGF textbook
//! - M. Büttiker, Phys. Rev. Lett. **65**, 2901 (1990) — shot noise / Landauer approach
//! - M. P. López Sancho, J. M. López Sancho, J. Rubio,
//!   J. Phys. F **15**, 851 (1985) — iterative surface Green's function
//! - Ya. M. Blanter, M. Büttiker, Phys. Rep. **336**, 1 (2000) — shot noise review
//!
//! # Quick Start
//!
//! ```rust
//! use spintronics::negf::{Hamiltonian1D, LeadSelfEnergy, GreenFunction, TransportCalculator};
//!
//! // 5-site tight-binding chain, hopping t = 1 eV, on-site ε = 0
//! let h = Hamiltonian1D::from_uniform(5, 0.0, 1.0).expect("valid parameters");
//! let sl = LeadSelfEnergy::new(0.5, 0.0).expect("valid");
//! let sr = LeadSelfEnergy::new(0.5, 0.0).expect("valid");
//! let gf = GreenFunction::new(h, sl, sr, 1e-3).expect("valid");
//! let tc = TransportCalculator::new(gf, -3.0, 3.0, 50).expect("valid");
//!
//! // Compute transmission at band centre
//! let t_at_ef = tc.transmission(0.0).expect("no singular matrix");
//! assert!(t_at_ef >= 0.0 && t_at_ef <= 1.0 + 1e-8);
//! ```

pub mod accumulation;
pub mod green_function;
pub mod keldysh;
pub mod shot_noise;

pub use accumulation::{BoundaryCondition, SpinAccumulation1D};
pub use green_function::{
    GreenFunction, Hamiltonian1D, LeadSelfEnergy, SanchoRubio, TransportCalculator,
};
pub use keldysh::KeldyshSolver;
pub use shot_noise::ShotNoise;
