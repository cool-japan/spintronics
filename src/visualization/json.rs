//! JSON Export for Structured Data with Metadata
//!
//! This module exports simulation data to JSON format, ideal for:
//! - Structured data with metadata (material properties, simulation parameters)
//! - Web applications and interactive visualizations
//! - Data interchange between different analysis tools
//! - Long-term archival with self-describing format
//!
//! ## Features
//!
//! - **Metadata preservation**: Material properties, simulation parameters
//! - **Nested structures**: Complex hierarchical data
//! - **Human-readable**: Easy to inspect and debug
//! - **Language-agnostic**: Easy to parse in Python, JavaScript, etc.
//!
//! ## Example Usage
//!
//! ```rust
//! use spintronics::visualization::json::{JsonWriter, SimulationData};
//! use spintronics::Vector3;
//!
//! let mut data = SimulationData::new("LLG_dynamics", "0.1.0");
//! data.add_metadata("material", "YIG");
//! data.add_metadata("temperature", "300K");
//! data.add_parameter("alpha", 0.001);
//! data.add_parameter("dt", 1e-12);
//!
//! // Add time series
//! data.time_series.push((0.0, vec![1.0, 0.0, 0.0]));
//! data.time_series.push((1e-12, vec![0.99, 0.01, 0.0]));
//!
//! JsonWriter::write("simulation.json", &data).unwrap();
//! ```

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Result, Write};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::vector3::Vector3;

/// Simulation data container with metadata
///
/// With the `serde` feature enabled, this struct can be serialized/deserialized
/// directly using serde_json for more flexible data handling.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SimulationData {
    /// Simulation type/name
    pub name: String,

    /// Version identifier
    pub version: String,

    /// String metadata (material name, notes, etc.)
    pub metadata: HashMap<String, String>,

    /// Numeric parameters (alpha, temperature, field, etc.)
    pub parameters: HashMap<String, f64>,

    /// Time series data: (time, values)
    pub time_series: Vec<(f64, Vec<f64>)>,

    /// Vector field snapshots: (time, vectors)
    pub vector_snapshots: Vec<(f64, Vec<Vector3<f64>>)>,

    /// Scalar field snapshots: (time, field_name, values)
    pub scalar_snapshots: Vec<(f64, String, Vec<f64>)>,
}

impl SimulationData {
    /// Create new simulation data container
    pub fn new(name: &str, version: &str) -> Self {
        Self {
            name: name.to_string(),
            version: version.to_string(),
            metadata: HashMap::new(),
            parameters: HashMap::new(),
            time_series: Vec::new(),
            vector_snapshots: Vec::new(),
            scalar_snapshots: Vec::new(),
        }
    }

    /// Add string metadata
    pub fn add_metadata(&mut self, key: &str, value: &str) {
        self.metadata.insert(key.to_string(), value.to_string());
    }

    /// Add numeric parameter
    pub fn add_parameter(&mut self, key: &str, value: f64) {
        self.parameters.insert(key.to_string(), value);
    }

    /// Add time series data point
    pub fn add_time_point(&mut self, time: f64, values: Vec<f64>) {
        self.time_series.push((time, values));
    }

    /// Add vector field snapshot
    pub fn add_vector_snapshot(&mut self, time: f64, vectors: Vec<Vector3<f64>>) {
        self.vector_snapshots.push((time, vectors));
    }

    /// Add scalar field snapshot
    pub fn add_scalar_snapshot(&mut self, time: f64, field_name: &str, values: Vec<f64>) {
        self.scalar_snapshots
            .push((time, field_name.to_string(), values));
    }
}

/// JSON writer for simulation data
pub struct JsonWriter;

impl JsonWriter {
    /// Write simulation data to JSON file
    ///
    /// # Arguments
    /// * `filename` - Output JSON filename
    /// * `data` - Simulation data to export
    pub fn write(filename: &str, data: &SimulationData) -> Result<()> {
        let file = File::create(filename)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "{{")?;
        writeln!(writer, "  \"name\": \"{}\",", data.name)?;
        writeln!(writer, "  \"version\": \"{}\",", data.version)?;

        // Metadata
        writeln!(writer, "  \"metadata\": {{")?;
        let metadata_entries: Vec<_> = data.metadata.iter().collect();
        for (i, (key, value)) in metadata_entries.iter().enumerate() {
            let comma = if i < metadata_entries.len() - 1 {
                ","
            } else {
                ""
            };
            writeln!(writer, "    \"{}\": \"{}\"{}", key, value, comma)?;
        }
        writeln!(writer, "  }},")?;

        // Parameters
        writeln!(writer, "  \"parameters\": {{")?;
        let param_entries: Vec<_> = data.parameters.iter().collect();
        for (i, (key, value)) in param_entries.iter().enumerate() {
            let comma = if i < param_entries.len() - 1 { "," } else { "" };
            writeln!(writer, "    \"{}\": {}{}", key, value, comma)?;
        }
        writeln!(writer, "  }},")?;

        // Time series
        writeln!(writer, "  \"time_series\": [")?;
        for (i, (time, values)) in data.time_series.iter().enumerate() {
            let comma = if i < data.time_series.len() - 1 {
                ","
            } else {
                ""
            };
            let values_str = values
                .iter()
                .map(|v| format!("{:.10e}", v))
                .collect::<Vec<_>>()
                .join(", ");
            writeln!(
                writer,
                "    {{\"time\": {:.10e}, \"values\": [{}]}}{}",
                time, values_str, comma
            )?;
        }
        writeln!(writer, "  ],")?;

        // Vector snapshots
        writeln!(writer, "  \"vector_snapshots\": [")?;
        for (i, (time, vectors)) in data.vector_snapshots.iter().enumerate() {
            let comma = if i < data.vector_snapshots.len() - 1 {
                ","
            } else {
                ""
            };
            writeln!(writer, "    {{\"time\": {:.10e}, \"vectors\": [", time)?;
            for (j, v) in vectors.iter().enumerate() {
                let vec_comma = if j < vectors.len() - 1 { "," } else { "" };
                writeln!(
                    writer,
                    "      [{:.10e}, {:.10e}, {:.10e}]{}",
                    v.x, v.y, v.z, vec_comma
                )?;
            }
            writeln!(writer, "    ]}}{}", comma)?;
        }
        writeln!(writer, "  ],")?;

        // Scalar snapshots
        writeln!(writer, "  \"scalar_snapshots\": [")?;
        for (i, (time, name, values)) in data.scalar_snapshots.iter().enumerate() {
            let comma = if i < data.scalar_snapshots.len() - 1 {
                ","
            } else {
                ""
            };
            let values_str = values
                .iter()
                .map(|v| format!("{:.10e}", v))
                .collect::<Vec<_>>()
                .join(", ");
            writeln!(
                writer,
                "    {{\"time\": {:.10e}, \"field\": \"{}\", \"values\": [{}]}}{}",
                time, name, values_str, comma
            )?;
        }
        writeln!(writer, "  ]")?;

        writeln!(writer, "}}")?;

        writer.flush()?;
        Ok(())
    }

    /// Write minimal JSON (only time series, no metadata)
    ///
    /// Useful for quick exports without full metadata
    pub fn write_simple(filename: &str, times: &[f64], data: &[Vec<f64>]) -> Result<()> {
        let file = File::create(filename)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "{{")?;
        writeln!(writer, "  \"data\": [")?;

        for (i, (time, values)) in times.iter().zip(data.iter()).enumerate() {
            let comma = if i < times.len() - 1 { "," } else { "" };
            let values_str = values
                .iter()
                .map(|v| format!("{:.10e}", v))
                .collect::<Vec<_>>()
                .join(", ");
            writeln!(
                writer,
                "    {{\"t\": {:.10e}, \"v\": [{}]}}{}",
                time, values_str, comma
            )?;
        }

        writeln!(writer, "  ]")?;
        writeln!(writer, "}}")?;

        writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simulation_data_creation() {
        let data = SimulationData::new("test", "1.0");
        assert_eq!(data.name, "test");
        assert_eq!(data.version, "1.0");
    }

    #[test]
    fn test_add_metadata() {
        let mut data = SimulationData::new("test", "1.0");
        data.add_metadata("material", "YIG");
        assert_eq!(data.metadata.get("material"), Some(&"YIG".to_string()));
    }

    #[test]
    fn test_add_parameter() {
        let mut data = SimulationData::new("test", "1.0");
        data.add_parameter("alpha", 0.001);
        assert_eq!(data.parameters.get("alpha"), Some(&0.001));
    }

    #[test]
    fn test_add_time_point() {
        let mut data = SimulationData::new("test", "1.0");
        data.add_time_point(0.0, vec![1.0, 0.0, 0.0]);
        assert_eq!(data.time_series.len(), 1);
    }

    #[test]
    fn test_add_vector_snapshot() {
        let mut data = SimulationData::new("test", "1.0");
        let vectors = vec![Vector3::new(1.0, 0.0, 0.0)];
        data.add_vector_snapshot(0.0, vectors);
        assert_eq!(data.vector_snapshots.len(), 1);
    }

    #[test]
    fn test_add_scalar_snapshot() {
        let mut data = SimulationData::new("test", "1.0");
        data.add_scalar_snapshot(0.0, "energy", vec![1.5, 2.3]);
        assert_eq!(data.scalar_snapshots.len(), 1);
    }

    #[test]
    fn test_json_write() {
        let mut data = SimulationData::new("llg_test", "0.1.0");
        data.add_metadata("material", "YIG");
        data.add_parameter("alpha", 0.001);
        data.add_time_point(0.0, vec![1.0, 0.0, 0.0]);

        let result = JsonWriter::write("/tmp/test_sim.json", &data);
        assert!(result.is_ok());
    }

    #[test]
    fn test_json_write_simple() {
        let times = vec![0.0, 1.0, 2.0];
        let data = vec![vec![1.0, 0.0], vec![0.5, 0.5], vec![0.0, 1.0]];

        let result = JsonWriter::write_simple("/tmp/test_simple.json", &times, &data);
        assert!(result.is_ok());
    }

    #[test]
    fn test_full_simulation_export() {
        let mut data = SimulationData::new("skyrmion_dynamics", "0.1.0");

        // Metadata
        data.add_metadata("material", "Pt/CoFeB");
        data.add_metadata("notes", "DMI-stabilized skyrmion");

        // Parameters
        data.add_parameter("alpha", 0.3);
        data.add_parameter("temperature", 300.0);
        data.add_parameter("dmi_d", 3.5e-3);

        // Time series
        data.add_time_point(0.0, vec![1.0, 0.0, 0.0]);
        data.add_time_point(1e-12, vec![0.99, 0.01, 0.0]);

        // Vector snapshot
        let vectors = vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)];
        data.add_vector_snapshot(0.0, vectors);

        // Scalar snapshot
        data.add_scalar_snapshot(0.0, "energy_density", vec![1.5, 2.3]);

        let result = JsonWriter::write("/tmp/test_full.json", &data);
        assert!(result.is_ok());
    }
}
