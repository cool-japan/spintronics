//! Magnetization dynamics solvers

pub mod llg;

pub use llg::{anisotropy_energy, calc_dm_dt, exchange_energy, zeeman_energy, LlgSolver};
