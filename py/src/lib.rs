//! Python bindings entry point for the spintronics library.
//!
//! This crate provides the `cdylib` entry point for maturin-based PyPI packaging.
//! All actual binding implementations live in the main `spintronics` crate under
//! the `python` feature flag (`src/python/`).
//!
//! ## Building
//!
//! ```bash
//! cd py
//! maturin develop --features pyo3/extension-module
//! ```
//!
//! ## Usage from Python
//!
//! ```python
//! import spintronics
//!
//! yig = spintronics.Ferromagnet.yig()
//! sim = spintronics.LlgSimulator(yig)
//! sim.set_external_field(0.0, 0.0, 0.1)
//! trajectory = sim.evolve(1e-9, 1000)
//! ```

use pyo3::prelude::*;

/// Register the spintronics Python module.
///
/// Delegates all class and constant registrations to the main crate's
/// `python::spintronics` function, which is defined in `src/python/mod.rs`.
/// The `::spintronics` path explicitly refers to the external crate (as opposed
/// to the pymodule function of the same name being declared below).
#[pymodule]
fn spintronics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    ::spintronics::python::spintronics(m)
}
