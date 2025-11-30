//! CSV Export for Time Series and Spatial Data
//!
//! This module exports simulation data to CSV (Comma-Separated Values) format
//! for analysis and plotting in tools like Python (matplotlib, pandas), MATLAB,
//! Excel, Origin, and gnuplot.
//!
//! ## Features
//!
//! - **Time series**: Magnetization dynamics, current evolution
//! - **Spatial profiles**: Spin configurations along a line or plane
//! - **Scalar fields**: Energy density, topological charge
//! - **Vector fields**: Spin vectors with components
//!
//! ## Example Usage
//!
//! ```rust
//! use spintronics::visualization::csv::CsvWriter;
//! use spintronics::Vector3;
//!
//! let mut writer = CsvWriter::new("magnetization_dynamics.csv").unwrap();
//!
//! // Write header
//! writer.write_header(&["time", "mx", "my", "mz"]).unwrap();
//!
//! // Write time series data
//! for t in 0..100 {
//!     let time = t as f64 * 0.01;
//!     let m = Vector3::new(time.cos(), time.sin(), 0.0);
//!     writer.write_row(&[time, m.x, m.y, m.z]).unwrap();
//! }
//! ```

use crate::vector3::Vector3;
use std::fs::File;
use std::io::{BufWriter, Result, Write};

/// CSV file writer for time series and spatial data
pub struct CsvWriter {
    writer: BufWriter<File>,
}

impl CsvWriter {
    /// Create a new CSV writer
    ///
    /// # Arguments
    /// * `filename` - Output CSV filename (e.g., "data.csv")
    pub fn new(filename: &str) -> Result<Self> {
        let file = File::create(filename)?;
        let writer = BufWriter::new(file);
        Ok(Self { writer })
    }

    /// Write header row with column names
    ///
    /// # Arguments
    /// * `columns` - Array of column names
    pub fn write_header(&mut self, columns: &[&str]) -> Result<()> {
        writeln!(self.writer, "{}", columns.join(","))
    }

    /// Write a data row
    ///
    /// # Arguments
    /// * `values` - Array of numeric values
    pub fn write_row(&mut self, values: &[f64]) -> Result<()> {
        let line = values
            .iter()
            .map(|v| format!("{:.10e}", v))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(self.writer, "{}", line)
    }

    /// Write vector field data
    ///
    /// # Arguments
    /// * `vectors` - Array of 3D vectors
    /// * `include_magnitude` - If true, add magnitude column
    pub fn write_vectors(
        &mut self,
        vectors: &[Vector3<f64>],
        include_magnitude: bool,
    ) -> Result<()> {
        if include_magnitude {
            for v in vectors {
                writeln!(
                    self.writer,
                    "{:.10e},{:.10e},{:.10e},{:.10e}",
                    v.x,
                    v.y,
                    v.z,
                    v.magnitude()
                )?;
            }
        } else {
            for v in vectors {
                writeln!(self.writer, "{:.10e},{:.10e},{:.10e}", v.x, v.y, v.z)?;
            }
        }
        Ok(())
    }

    /// Flush buffered data to disk
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()
    }
}

/// Write time series data to CSV
///
/// Convenience function for simple time series export
///
/// # Arguments
/// * `filename` - Output CSV filename
/// * `times` - Time points
/// * `data` - Data arrays (each array is one column)
/// * `column_names` - Names for each data column
pub fn write_time_series(
    filename: &str,
    times: &[f64],
    data: &[&[f64]],
    column_names: &[&str],
) -> Result<()> {
    let mut writer = CsvWriter::new(filename)?;

    // Write header
    let mut header = vec!["time"];
    header.extend_from_slice(column_names);
    writer.write_header(&header)?;

    // Write data rows
    for (i, &t) in times.iter().enumerate() {
        let mut row = vec![t];
        for dataset in data {
            if i < dataset.len() {
                row.push(dataset[i]);
            }
        }
        writer.write_row(&row)?;
    }

    writer.flush()?;
    Ok(())
}

/// Write 2D spatial data (grid) to CSV
///
/// Exports a 2D grid with coordinates and values
///
/// # Arguments
/// * `filename` - Output CSV filename
/// * `data` - 2D array (flattened) of scalar values
/// * `dims` - Grid dimensions (nx, ny)
pub fn write_2d_grid(filename: &str, data: &[f64], dims: (usize, usize)) -> Result<()> {
    let (nx, ny) = dims;

    if data.len() != nx * ny {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("Expected {} points, got {}", nx * ny, data.len()),
        ));
    }

    let mut writer = CsvWriter::new(filename)?;
    writer.write_header(&["x", "y", "value"])?;

    for j in 0..ny {
        for i in 0..nx {
            let idx = i + nx * j;
            writer.write_row(&[i as f64, j as f64, data[idx]])?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Write vector field to CSV with spatial coordinates
///
/// Exports a 3D vector field on a regular grid
///
/// # Arguments
/// * `filename` - Output CSV filename
/// * `vectors` - Array of 3D vectors (flattened from grid)
/// * `dims` - Grid dimensions (nx, ny, nz)
pub fn write_vector_field(
    filename: &str,
    vectors: &[Vector3<f64>],
    dims: (usize, usize, usize),
) -> Result<()> {
    let (nx, ny, nz) = dims;

    if vectors.len() != nx * ny * nz {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("Expected {} vectors, got {}", nx * ny * nz, vectors.len()),
        ));
    }

    let mut writer = CsvWriter::new(filename)?;
    writer.write_header(&["x", "y", "z", "vx", "vy", "vz", "magnitude"])?;

    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let idx = i + nx * j + nx * ny * k;
                let v = &vectors[idx];
                writer.write_row(&[i as f64, j as f64, k as f64, v.x, v.y, v.z, v.magnitude()])?;
            }
        }
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_csv_writer_creation() {
        let writer = CsvWriter::new("/tmp/test_csv.csv");
        assert!(writer.is_ok());
    }

    #[test]
    fn test_write_header() {
        let mut writer = CsvWriter::new("/tmp/test_header.csv").unwrap();
        let result = writer.write_header(&["time", "mx", "my", "mz"]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_write_row() {
        let mut writer = CsvWriter::new("/tmp/test_row.csv").unwrap();
        writer.write_header(&["x", "y", "z"]).unwrap();
        let result = writer.write_row(&[1.0, 2.0, 3.0]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_write_vectors() {
        let mut writer = CsvWriter::new("/tmp/test_vectors.csv").unwrap();
        let vectors = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let result = writer.write_vectors(&vectors, true);
        assert!(result.is_ok());
    }

    #[test]
    fn test_time_series() {
        let times = vec![0.0, 1.0, 2.0, 3.0];
        let mx = vec![1.0, 0.5, 0.0, -0.5];
        let my = vec![0.0, 0.5, 1.0, 0.5];

        let result = write_time_series(
            "/tmp/test_time_series.csv",
            &times,
            &[&mx, &my],
            &["mx", "my"],
        );
        assert!(result.is_ok());
    }

    #[test]
    fn test_2d_grid() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let result = write_2d_grid("/tmp/test_grid.csv", &data, (2, 2));
        assert!(result.is_ok());
    }

    #[test]
    fn test_2d_grid_size_mismatch() {
        let data = vec![1.0, 2.0, 3.0];
        let result = write_2d_grid("/tmp/test_error.csv", &data, (2, 2));
        assert!(result.is_err());
    }

    #[test]
    fn test_vector_field_export() {
        let vectors = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 0.0),
        ];
        let result = write_vector_field("/tmp/test_vfield.csv", &vectors, (2, 2, 1));
        assert!(result.is_ok());
    }

    #[test]
    fn test_vector_field_size_mismatch() {
        let vectors = vec![Vector3::new(1.0, 0.0, 0.0)];
        let result = write_vector_field("/tmp/test_vfield_error.csv", &vectors, (2, 2, 1));
        assert!(result.is_err());
    }
}
