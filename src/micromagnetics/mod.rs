//! Finite-difference micromagnetics solver
//!
//! Implements a structured-grid FD micromagnetics solver with:
//! - Newell analytic demagnetization tensor (diagonal components)
//! - Full effective field: exchange (FD Laplacian) + demag + Zeeman + anisotropy
//! - LLG time integration (RK4)
//!
//! The demagnetization tensor is computed using the analytic Newell (1993) formulas,
//! which give the exact dipolar interaction between uniformly magnetized rectangular
//! prisms. Only the diagonal components (N_xx, N_yy, N_zz) are retained; the
//! off-diagonal terms vanish for cubic/tetragonal grids aligned with the coordinate
//! axes and represent small corrections for other geometries.
//!
//! # References
//! - A. J. Newell, W. Williams, D. J. Dunlop, J. Geophys. Res. 98, 9551 (1993)
//! - D. V. Berkov, J. Magn. Magn. Mater. 186, 199 (1998)
//! - muMAG Standard Problems: <https://www.ctcms.nist.gov/~rdm/mumag.org.html>

pub mod demag;
pub mod grid;

pub use demag::{DemagField, NewellTensor};
pub use grid::{GridConfig, LlgResult, MicromagneticGrid};
