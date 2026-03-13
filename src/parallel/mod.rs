//! Enhanced parallel computing for spintronics simulations
//!
//! Provides domain decomposition and parallel parameter sweeps
//! using rayon for thread-level parallelism.

pub mod decomposition;
pub mod sweep;

pub use decomposition::{Domain, DomainDecomposition, DomainDecomposition2D, Tile2D};
pub use sweep::{
    field_sweep, parallel_sweep, parallel_sweep_2d, parallel_sweep_with_progress,
    temperature_sweep, FieldSweepResult, ParameterSweep, SweepPoint2D, TemperatureSweepResult,
};
