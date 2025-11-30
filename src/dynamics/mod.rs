//! Magnetization dynamics solvers

pub mod llg;

pub use llg::{calc_dm_dt, LlgSolver};
