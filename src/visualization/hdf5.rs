//! HDF5 Export for Large-Scale Simulation Data
//!
//! This module provides HDF5 file format support for storing and retrieving
//! large simulation datasets efficiently. HDF5 (Hierarchical Data Format 5)
//! is ideal for:
//!
//! - **Large datasets**: Efficient storage of multi-GB simulation data
//! - **Hierarchical organization**: Group related data (parameters, time series, snapshots)
//! - **Compression**: Built-in compression for reduced file sizes
//! - **Cross-platform**: Read by Python (h5py), MATLAB, Julia, C/C++
//! - **Partial I/O**: Read/write subsets without loading entire file
//!
//! ## Example Usage
//!
//! ```rust,ignore
//! use spintronics::visualization::hdf5::Hdf5Writer;
//! use spintronics::Vector3;
//!
//! let mut writer = Hdf5Writer::create("simulation.h5").unwrap();
//!
//! // Write metadata
//! writer.write_attribute("material", "YIG").unwrap();
//! writer.write_attribute("version", "0.2.0").unwrap();
//!
//! // Write parameters
//! writer.write_scalar("parameters/alpha", 0.001).unwrap();
//! writer.write_scalar("parameters/temperature", 300.0).unwrap();
//!
//! // Write time series
//! let times = vec![0.0, 1e-12, 2e-12, 3e-12];
//! writer.write_array("time", &times).unwrap();
//!
//! // Write magnetization trajectory
//! let mx = vec![1.0, 0.99, 0.97, 0.95];
//! writer.write_array("magnetization/mx", &mx).unwrap();
//!
//! // Write vector field snapshot
//! let spins: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
//! writer.write_vector_field("snapshots/t0/spins", &spins).unwrap();
//! ```
//!
//! ## Reading in Python
//!
//! ```python
//! import h5py
//! import numpy as np
//!
//! with h5py.File('simulation.h5', 'r') as f:
//!     # Read metadata
//!     material = f.attrs['material']
//!
//!     # Read parameters
//!     alpha = f['parameters/alpha'][()]
//!
//!     # Read time series
//!     times = f['time'][:]
//!     mx = f['magnetization/mx'][:]
//!
//!     # Read vector field
//!     spins = f['snapshots/t0/spins'][:]
//! ```

#[cfg(not(feature = "hdf5"))]
use std::io;

#[cfg(feature = "hdf5")]
use hdf5::{File, Result as Hdf5Result};

#[cfg(feature = "hdf5")]
use crate::vector3::Vector3;

/// HDF5 file writer for simulation data
///
/// Provides a hierarchical interface for storing simulation data
/// with metadata, parameters, and array datasets.
#[cfg(feature = "hdf5")]
pub struct Hdf5Writer {
    file: File,
}

#[cfg(feature = "hdf5")]
impl Hdf5Writer {
    /// Create a new HDF5 file for writing
    ///
    /// # Arguments
    /// * `filename` - Path to the HDF5 file to create
    ///
    /// # Returns
    /// A new Hdf5Writer instance
    pub fn create(filename: &str) -> Hdf5Result<Self> {
        let file = File::create(filename)?;
        Ok(Self { file })
    }

    /// Open an existing HDF5 file for appending
    pub fn open(filename: &str) -> Hdf5Result<Self> {
        let file = File::open_rw(filename)?;
        Ok(Self { file })
    }

    /// Write a string attribute to the root group
    pub fn write_attribute(&self, name: &str, value: &str) -> Hdf5Result<()> {
        use std::str::FromStr;
        let attr = self
            .file
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create(name)?;
        let unicode_value = hdf5::types::VarLenUnicode::from_str(value)
            .map_err(|_| hdf5::Error::Internal("String conversion error".to_string()))?;
        attr.write_scalar(&unicode_value)?;
        Ok(())
    }

    /// Write a scalar value to a dataset
    pub fn write_scalar(&self, path: &str, value: f64) -> Hdf5Result<()> {
        // Create parent groups if needed
        self.ensure_groups(path)?;

        let dataset = self.file.new_dataset::<f64>().create(path)?;
        dataset.write_scalar(&value)?;
        Ok(())
    }

    /// Write a 1D array to a dataset
    pub fn write_array(&self, path: &str, data: &[f64]) -> Hdf5Result<()> {
        self.ensure_groups(path)?;

        let dataset = self
            .file
            .new_dataset::<f64>()
            .shape([data.len()])
            .create(path)?;
        dataset.write(data)?;
        Ok(())
    }

    /// Write a 2D array to a dataset
    pub fn write_array_2d(&self, path: &str, data: &[Vec<f64>]) -> Hdf5Result<()> {
        if data.is_empty() {
            return Ok(());
        }

        self.ensure_groups(path)?;

        let rows = data.len();
        let cols = data[0].len();

        // Flatten the data
        let flat: Vec<f64> = data.iter().flat_map(|row| row.iter().copied()).collect();

        let dataset = self
            .file
            .new_dataset::<f64>()
            .shape([rows, cols])
            .create(path)?;
        dataset.write(&flat)?;
        Ok(())
    }

    /// Write a vector field (array of Vector3) to a dataset
    pub fn write_vector_field(&self, path: &str, vectors: &[Vector3<f64>]) -> Hdf5Result<()> {
        self.ensure_groups(path)?;

        // Convert to [N, 3] array
        let data: Vec<f64> = vectors.iter().flat_map(|v| [v.x, v.y, v.z]).collect();

        let dataset = self
            .file
            .new_dataset::<f64>()
            .shape([vectors.len(), 3])
            .create(path)?;
        dataset.write(&data)?;
        Ok(())
    }

    /// Write a time series with associated data
    pub fn write_time_series(
        &self,
        group_path: &str,
        times: &[f64],
        data: &[(String, Vec<f64>)],
    ) -> Hdf5Result<()> {
        // Create the group
        self.file.create_group(group_path)?;

        // Write time array
        self.write_array(&format!("{}/time", group_path), times)?;

        // Write each data column
        for (name, values) in data {
            self.write_array(&format!("{}/{}", group_path, name), values)?;
        }

        Ok(())
    }

    /// Ensure parent groups exist for a dataset path
    fn ensure_groups(&self, path: &str) -> Hdf5Result<()> {
        let parts: Vec<&str> = path.split('/').collect();
        if parts.len() <= 1 {
            return Ok(());
        }

        let mut current_path = String::new();
        for part in &parts[..parts.len() - 1] {
            if current_path.is_empty() {
                current_path = part.to_string();
            } else {
                current_path = format!("{}/{}", current_path, part);
            }

            // Try to create group, ignore if it already exists
            let _ = self.file.create_group(&current_path);
        }

        Ok(())
    }
}

/// HDF5 file reader for simulation data
#[cfg(feature = "hdf5")]
pub struct Hdf5Reader {
    file: File,
}

#[cfg(feature = "hdf5")]
impl Hdf5Reader {
    /// Open an HDF5 file for reading
    pub fn open(filename: &str) -> Hdf5Result<Self> {
        let file = File::open(filename)?;
        Ok(Self { file })
    }

    /// Read a string attribute from the root group
    pub fn read_attribute(&self, name: &str) -> Hdf5Result<String> {
        let attr = self.file.attr(name)?;
        let value: hdf5::types::VarLenUnicode = attr.read_scalar()?;
        Ok(value.to_string())
    }

    /// Read a scalar value from a dataset
    pub fn read_scalar(&self, path: &str) -> Hdf5Result<f64> {
        let dataset = self.file.dataset(path)?;
        let value: f64 = dataset.read_scalar()?;
        Ok(value)
    }

    /// Read a 1D array from a dataset
    pub fn read_array(&self, path: &str) -> Hdf5Result<Vec<f64>> {
        let dataset = self.file.dataset(path)?;
        let data: Vec<f64> = dataset.read_raw()?;
        Ok(data)
    }

    /// Read a vector field from a dataset
    pub fn read_vector_field(&self, path: &str) -> Hdf5Result<Vec<Vector3<f64>>> {
        let dataset = self.file.dataset(path)?;
        let shape = dataset.shape();

        if shape.len() != 2 || shape[1] != 3 {
            return Err(hdf5::Error::from("Invalid vector field shape"));
        }

        let data: Vec<f64> = dataset.read_raw()?;
        let vectors: Vec<Vector3<f64>> = data
            .chunks(3)
            .map(|chunk| Vector3::new(chunk[0], chunk[1], chunk[2]))
            .collect();

        Ok(vectors)
    }

    /// List all datasets in a group
    pub fn list_datasets(&self, group_path: &str) -> Hdf5Result<Vec<String>> {
        let group = if group_path.is_empty() || group_path == "/" {
            self.file.as_group()?
        } else {
            self.file.group(group_path)?
        };

        let names: Vec<String> = group.member_names()?;
        Ok(names)
    }
}

/// Stub implementation when HDF5 feature is not enabled
#[cfg(not(feature = "hdf5"))]
pub struct Hdf5Writer;

#[cfg(not(feature = "hdf5"))]
impl Hdf5Writer {
    /// Create a new HDF5 file (stub - returns error)
    pub fn create(_filename: &str) -> io::Result<Self> {
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "HDF5 support not enabled. Compile with --features hdf5",
        ))
    }
}

/// Stub implementation when HDF5 feature is not enabled
#[cfg(not(feature = "hdf5"))]
pub struct Hdf5Reader;

#[cfg(not(feature = "hdf5"))]
impl Hdf5Reader {
    /// Open an HDF5 file (stub - returns error)
    pub fn open(_filename: &str) -> io::Result<Self> {
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "HDF5 support not enabled. Compile with --features hdf5",
        ))
    }
}

#[cfg(test)]
mod tests {
    #[allow(unused_imports)]
    use super::*;

    #[test]
    #[cfg(not(feature = "hdf5"))]
    fn test_hdf5_disabled() {
        let result = Hdf5Writer::create("/tmp/test.h5");
        assert!(result.is_err());
    }
}
