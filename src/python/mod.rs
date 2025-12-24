//! Python bindings for the spintronics library using PyO3
//!
//! This module provides Python bindings for the core spintronics functionality,
//! enabling Python users to access high-performance Rust simulations.
//!
//! ## Usage from Python
//!
//! ```python
//! import spintronics
//!
//! # Create a 3D vector
//! v = spintronics.Vector3(1.0, 0.0, 0.0)
//! print(f"Magnitude: {v.magnitude()}")
//!
//! # Create materials
//! yig = spintronics.Ferromagnet.yig()
//! pt = spintronics.InverseSpinHall.platinum()
//!
//! # Run LLG simulation
//! sim = spintronics.LlgSimulator(yig)
//! sim.set_external_field(0.0, 0.0, 0.1)
//! result = sim.evolve(1e-9, 1000)
//! ```

use pyo3::prelude::*;

mod dynamics;
mod effects;
mod materials;
mod simulation;
mod vector;

pub use dynamics::PyLlgSimulator;
pub use effects::PyInverseSpinHall;
pub use materials::{PyFerromagnet, PySpinInterface};
pub use simulation::PySpinPumpingSimulation;
pub use vector::PyVector3;

/// Python module for spintronics simulations
#[pymodule]
pub fn spintronics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Core types
    m.add_class::<PyVector3>()?;

    // Materials
    m.add_class::<PyFerromagnet>()?;
    m.add_class::<PySpinInterface>()?;

    // Effects
    m.add_class::<PyInverseSpinHall>()?;

    // Dynamics
    m.add_class::<PyLlgSimulator>()?;

    // High-level simulations
    m.add_class::<PySpinPumpingSimulation>()?;

    // Physical constants
    m.add("HBAR", crate::constants::HBAR)?;
    m.add("GAMMA", crate::constants::GAMMA)?;
    m.add("E_CHARGE", crate::constants::E_CHARGE)?;
    m.add("MU_B", crate::constants::MU_B)?;
    m.add("KB", crate::constants::KB)?;

    Ok(())
}
