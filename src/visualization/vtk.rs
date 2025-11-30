//! VTK Export for ParaView Visualization
//!
//! This module exports simulation results to VTK (Visualization Toolkit) format,
//! specifically VTU (VTK Unstructured Grid XML) files that can be opened in
//! ParaView, VisIt, and other scientific visualization software.
//!
//! ## Visualization Capabilities
//!
//! - **Spin configurations**: Vector field visualization with glyphs (arrows)
//! - **Magnetization textures**: Skyrmions, domain walls, vortices
//! - **Time series**: Animation of spin dynamics
//! - **Scalar fields**: Energy density, topological charge density
//! - **Multi-field visualization**: Combine spin, temperature, current, etc.
//!
//! ## ParaView Workflow
//!
//! 1. Export simulation data using `VtkWriter`
//! 2. Open .vtu files in ParaView
//! 3. Apply "Glyph" filter to visualize spin vectors as arrows
//! 4. Use "Slice" or "Volume" rendering for 3D data
//! 5. Create animations from time series
//!
//! ## Example Usage
//!
//! ```rust
//! use spintronics::visualization::vtk::VtkWriter;
//! use spintronics::Vector3;
//!
//! let mut writer = VtkWriter::new("spin_dynamics");
//!
//! // Create a simple spin configuration
//! let spins = vec![
//!     Vector3::new(1.0, 0.0, 0.0),
//!     Vector3::new(0.0, 1.0, 0.0),
//!     Vector3::new(0.0, 0.0, 1.0),
//! ];
//!
//! // Export snapshot
//! writer.write_snapshot(&spins, (3, 1, 1)).unwrap();
//! ```

use crate::vector3::Vector3;
use std::fs::File;
use std::io::{BufWriter, Result, Write};

/// VTK file writer for exporting spin configurations to ParaView
///
/// Generates VTU (VTK Unstructured Grid XML) files with proper formatting
/// for scientific visualization software.
pub struct VtkWriter {
    /// Base filename (without extension)
    filename_base: String,

    /// Current time step counter
    step_counter: usize,

    /// Include cell data (for volume rendering)
    include_cells: bool,
}

impl VtkWriter {
    /// Create a new VTK writer
    ///
    /// # Arguments
    /// * `filename_base` - Base name for output files (e.g., "skyrmion")
    ///
    /// Files will be named: `{filename_base}_{step:05}.vtu`
    pub fn new(filename_base: &str) -> Self {
        Self {
            filename_base: filename_base.to_string(),
            step_counter: 0,
            include_cells: false,
        }
    }

    /// Enable cell data export for volume rendering
    pub fn with_cells(mut self) -> Self {
        self.include_cells = true;
        self
    }

    /// Reset step counter
    pub fn reset_counter(&mut self) {
        self.step_counter = 0;
    }

    /// Write a snapshot of spin configuration
    ///
    /// # Arguments
    /// * `spins` - Array of spin vectors (flattened from 3D grid)
    /// * `dims` - Grid dimensions (nx, ny, nz)
    ///
    /// # Grid Layout
    /// Spins are indexed as: `spins[i + nx*j + nx*ny*k]` for point (i,j,k)
    pub fn write_snapshot(
        &mut self,
        spins: &[Vector3<f64>],
        dims: (usize, usize, usize),
    ) -> Result<()> {
        let (nx, ny, nz) = dims;
        let num_points = nx * ny * nz;

        if spins.len() != num_points {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Expected {} spins, got {}", num_points, spins.len()),
            ));
        }

        let filename = format!("{}_{:05}.vtu", self.filename_base, self.step_counter);
        let file = File::create(&filename)?;
        let mut writer = BufWriter::new(file);

        self.write_header(&mut writer)?;
        self.write_point_data(&mut writer, spins)?;
        self.write_points(&mut writer, dims)?;

        if self.include_cells {
            self.write_cells(&mut writer, dims)?;
        } else {
            self.write_empty_cells(&mut writer)?;
        }

        self.write_footer(&mut writer)?;

        self.step_counter += 1;
        Ok(())
    }

    /// Write snapshot with additional scalar field
    ///
    /// Useful for visualizing energy density, topological charge, etc.
    pub fn write_snapshot_with_scalar(
        &mut self,
        spins: &[Vector3<f64>],
        scalar_data: &[f64],
        scalar_name: &str,
        dims: (usize, usize, usize),
    ) -> Result<()> {
        let (nx, ny, nz) = dims;
        let num_points = nx * ny * nz;

        if spins.len() != num_points || scalar_data.len() != num_points {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Size mismatch in data arrays",
            ));
        }

        let filename = format!("{}_{:05}.vtu", self.filename_base, self.step_counter);
        let file = File::create(&filename)?;
        let mut writer = BufWriter::new(file);

        self.write_header(&mut writer)?;

        // Write point data with additional scalar field
        writeln!(
            writer,
            "      <PointData Vectors=\"Spin\" Scalars=\"{}\">\n",
            scalar_name
        )?;
        self.write_vector_array(&mut writer, spins, "Spin")?;
        self.write_scalar_component(&mut writer, spins, "Sz")?;
        self.write_scalar_array(&mut writer, scalar_data, scalar_name)?;
        writeln!(writer, "      </PointData>")?;

        self.write_points(&mut writer, dims)?;

        if self.include_cells {
            self.write_cells(&mut writer, dims)?;
        } else {
            self.write_empty_cells(&mut writer)?;
        }

        self.write_footer(&mut writer)?;

        self.step_counter += 1;
        Ok(())
    }

    fn write_header(&self, writer: &mut BufWriter<File>) -> Result<()> {
        writeln!(writer, "<?xml version=\"1.0\"?>")?;
        writeln!(
            writer,
            "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">"
        )?;
        writeln!(writer, "  <UnstructuredGrid>")?;
        Ok(())
    }

    fn write_point_data(&self, writer: &mut BufWriter<File>, spins: &[Vector3<f64>]) -> Result<()> {
        let num_points = spins.len();
        let num_cells = if self.include_cells {
            num_points.saturating_sub(1)
        } else {
            0
        };

        writeln!(
            writer,
            "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
            num_points, num_cells
        )?;

        writeln!(writer, "      <PointData Vectors=\"Spin\" Scalars=\"Sz\">")?;

        // Spin vectors
        self.write_vector_array(writer, spins, "Spin")?;

        // Z-component for color mapping
        self.write_scalar_component(writer, spins, "Sz")?;

        // Magnetization magnitude
        self.write_magnitude_array(writer, spins)?;

        writeln!(writer, "      </PointData>")?;
        Ok(())
    }

    fn write_vector_array(
        &self,
        writer: &mut BufWriter<File>,
        vectors: &[Vector3<f64>],
        name: &str,
    ) -> Result<()> {
        writeln!(
            writer,
            "        <DataArray type=\"Float32\" Name=\"{}\" NumberOfComponents=\"3\" format=\"ascii\">",
            name
        )?;
        for v in vectors {
            write!(writer, "{:.6} {:.6} {:.6} ", v.x, v.y, v.z)?;
        }
        writeln!(writer, "\n        </DataArray>")?;
        Ok(())
    }

    fn write_scalar_component(
        &self,
        writer: &mut BufWriter<File>,
        vectors: &[Vector3<f64>],
        name: &str,
    ) -> Result<()> {
        writeln!(
            writer,
            "        <DataArray type=\"Float32\" Name=\"{}\" NumberOfComponents=\"1\" format=\"ascii\">",
            name
        )?;
        for v in vectors {
            write!(writer, "{:.6} ", v.z)?;
        }
        writeln!(writer, "\n        </DataArray>")?;
        Ok(())
    }

    fn write_scalar_array(
        &self,
        writer: &mut BufWriter<File>,
        scalars: &[f64],
        name: &str,
    ) -> Result<()> {
        writeln!(
            writer,
            "        <DataArray type=\"Float32\" Name=\"{}\" NumberOfComponents=\"1\" format=\"ascii\">",
            name
        )?;
        for s in scalars {
            write!(writer, "{:.6} ", s)?;
        }
        writeln!(writer, "\n        </DataArray>")?;
        Ok(())
    }

    fn write_magnitude_array(
        &self,
        writer: &mut BufWriter<File>,
        vectors: &[Vector3<f64>],
    ) -> Result<()> {
        writeln!(
            writer,
            "        <DataArray type=\"Float32\" Name=\"Magnitude\" NumberOfComponents=\"1\" format=\"ascii\">"
        )?;
        for v in vectors {
            write!(writer, "{:.6} ", v.magnitude())?;
        }
        writeln!(writer, "\n        </DataArray>")?;
        Ok(())
    }

    fn write_points(
        &self,
        writer: &mut BufWriter<File>,
        dims: (usize, usize, usize),
    ) -> Result<()> {
        let (nx, ny, nz) = dims;

        writeln!(writer, "      <Points>")?;
        writeln!(
            writer,
            "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"
        )?;

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    write!(writer, "{} {} {} ", i, j, k)?;
                }
            }
        }

        writeln!(writer, "\n        </DataArray>")?;
        writeln!(writer, "      </Points>")?;
        Ok(())
    }

    fn write_empty_cells(&self, writer: &mut BufWriter<File>) -> Result<()> {
        writeln!(writer, "      <Cells>")?;
        writeln!(
            writer,
            "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n        </DataArray>"
        )?;
        writeln!(
            writer,
            "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n        </DataArray>"
        )?;
        writeln!(
            writer,
            "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n        </DataArray>"
        )?;
        writeln!(writer, "      </Cells>")?;
        Ok(())
    }

    fn write_cells(&self, writer: &mut BufWriter<File>, dims: (usize, usize, usize)) -> Result<()> {
        let (nx, ny, nz) = dims;

        writeln!(writer, "      <Cells>")?;

        // Connectivity for voxel cells
        writeln!(
            writer,
            "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
        )?;

        for k in 0..(nz - 1) {
            for j in 0..(ny - 1) {
                for i in 0..(nx - 1) {
                    let idx = |di: usize, dj: usize, dk: usize| {
                        (i + di) + nx * (j + dj) + nx * ny * (k + dk)
                    };

                    // VTK voxel: 8 vertices in specific order
                    write!(
                        writer,
                        "{} {} {} {} {} {} {} {} ",
                        idx(0, 0, 0),
                        idx(1, 0, 0),
                        idx(0, 1, 0),
                        idx(1, 1, 0),
                        idx(0, 0, 1),
                        idx(1, 0, 1),
                        idx(0, 1, 1),
                        idx(1, 1, 1)
                    )?;
                }
            }
        }
        writeln!(writer, "\n        </DataArray>")?;

        // Offsets
        writeln!(
            writer,
            "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
        )?;
        let num_cells = (nx - 1) * (ny - 1) * (nz - 1);
        for i in 1..=num_cells {
            write!(writer, "{} ", i * 8)?;
        }
        writeln!(writer, "\n        </DataArray>")?;

        // Types (11 = VTK_VOXEL)
        writeln!(
            writer,
            "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
        )?;
        for _ in 0..num_cells {
            write!(writer, "11 ")?;
        }
        writeln!(writer, "\n        </DataArray>")?;

        writeln!(writer, "      </Cells>")?;
        Ok(())
    }

    fn write_footer(&self, writer: &mut BufWriter<File>) -> Result<()> {
        writeln!(writer, "    </Piece>")?;
        writeln!(writer, "  </UnstructuredGrid>")?;
        writeln!(writer, "</VTKFile>")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vtk_writer_creation() {
        let writer = VtkWriter::new("test_output");
        assert_eq!(writer.step_counter, 0);
        assert_eq!(writer.filename_base, "test_output");
    }

    #[test]
    fn test_write_simple_snapshot() {
        let mut writer = VtkWriter::new("/tmp/test_spin");

        let spins = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(-1.0, 0.0, 0.0),
        ];

        let result = writer.write_snapshot(&spins, (2, 2, 1));
        assert!(result.is_ok());
        assert_eq!(writer.step_counter, 1);
    }

    #[test]
    fn test_write_multiple_snapshots() {
        let mut writer = VtkWriter::new("/tmp/test_series");

        let spins = vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)];

        for _ in 0..3 {
            writer.write_snapshot(&spins, (2, 1, 1)).unwrap();
        }

        assert_eq!(writer.step_counter, 3);
    }

    #[test]
    fn test_size_mismatch_error() {
        let mut writer = VtkWriter::new("/tmp/test_error");

        let spins = vec![Vector3::new(1.0, 0.0, 0.0)];

        // Wrong dimensions (expects 4 spins)
        let result = writer.write_snapshot(&spins, (2, 2, 1));
        assert!(result.is_err());
    }

    #[test]
    fn test_with_scalar_field() {
        let mut writer = VtkWriter::new("/tmp/test_scalar");

        let spins = vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)];

        let energy = vec![1.5, 2.3];

        let result = writer.write_snapshot_with_scalar(&spins, &energy, "Energy", (2, 1, 1));
        assert!(result.is_ok());
    }

    #[test]
    fn test_reset_counter() {
        let mut writer = VtkWriter::new("/tmp/test_reset");

        let spins = vec![Vector3::new(1.0, 0.0, 0.0)];

        writer.write_snapshot(&spins, (1, 1, 1)).unwrap();
        assert_eq!(writer.step_counter, 1);

        writer.reset_counter();
        assert_eq!(writer.step_counter, 0);
    }

    #[test]
    fn test_with_cells() {
        let writer = VtkWriter::new("/tmp/test_cells").with_cells();
        assert!(writer.include_cells);
    }
}
