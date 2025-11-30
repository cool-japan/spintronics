//! Visualization and Data Export
//!
//! This module provides tools for exporting simulation data to various formats
//! for visualization and analysis.
//!
//! ## Available Exporters
//!
//! - **VTK (ParaView)**: 3D vector field visualization, time series animation
//! - **CSV**: Simple tabular data for plotting and analysis
//! - **JSON**: Structured data export with metadata
//!
//! ## Example Usage
//!
//! ### VTK Export for ParaView
//!
//! ```rust
//! use spintronics::visualization::vtk::VtkWriter;
//! use spintronics::Vector3;
//!
//! let mut writer = VtkWriter::new("skyrmion_dynamics");
//!
//! // Simulate and export
//! let spins = vec![
//!     Vector3::new(1.0, 0.0, 0.0),
//!     Vector3::new(0.0, 1.0, 0.0),
//!     Vector3::new(0.0, 0.0, 1.0),
//! ];
//!
//! writer.write_snapshot(&spins, (3, 1, 1)).unwrap();
//! ```
//!
//! ### CSV Export for Plotting
//!
//! ```rust
//! use spintronics::visualization::csv::CsvWriter;
//!
//! let mut writer = CsvWriter::new("magnetization.csv").unwrap();
//! writer.write_header(&["time", "mx", "my", "mz"]).unwrap();
//! writer.write_row(&[0.0, 1.0, 0.0, 0.0]).unwrap();
//! ```
//!
//! ### JSON Export with Metadata
//!
//! ```rust
//! use spintronics::visualization::json::{JsonWriter, SimulationData};
//!
//! let mut data = SimulationData::new("LLG_dynamics", "0.1.0");
//! data.add_metadata("material", "YIG");
//! data.add_parameter("alpha", 0.001);
//! data.add_time_point(0.0, vec![1.0, 0.0, 0.0]);
//!
//! JsonWriter::write("simulation.json", &data).unwrap();
//! ```

pub mod csv;
pub mod json;
pub mod vtk;

pub use csv::CsvWriter;
pub use json::{JsonWriter, SimulationData};
pub use vtk::VtkWriter;
