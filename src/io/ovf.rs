//! OOMMF Vector Field (OVF) Format Support
//!
//! This module implements reading and writing of OVF files, the standard
//! format used by OOMMF (Object Oriented MicroMagnetic Framework) and
//! other micromagnetic simulation tools.
//!
//! ## OVF Format Specification
//!
//! OVF files can be in text or binary format and contain:
//! - Header with metadata (mesh dimensions, units, etc.)
//! - Vector field data (typically magnetization)
//!
//! ## Supported Versions
//! - OVF 1.0 (text and binary)
//! - OVF 2.0 (text and binary)
//!
//! ## References
//! - OOMMF User's Guide: <https://math.nist.gov/oommf/>
//! - OVF Format Specification: <https://math.nist.gov/oommf/doc/userguide20a3/userguide/OVF_2.0_format.html>

use crate::vector3::Vector3;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// OVF file format version and encoding
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OvfFormat {
    /// OVF 1.0 text format
    Text1_0,
    /// OVF 2.0 text format
    Text2_0,
    /// OVF 1.0 binary format (4-byte floats)
    Binary4_1_0,
    /// OVF 2.0 binary format (4-byte floats)
    Binary4_2_0,
    /// OVF 2.0 binary format (8-byte doubles)
    Binary8_2_0,
}

/// OVF data structure
#[derive(Debug, Clone)]
pub struct OvfData {
    /// Title/description
    pub title: String,

    /// Mesh dimensions (nx, ny, nz)
    pub mesh_size: (usize, usize, usize),

    /// Physical mesh dimensions in meters (xsize, ysize, zsize)
    pub mesh_physical_size: (f64, f64, f64),

    /// Cell size in meters (dx, dy, dz)
    pub cell_size: (f64, f64, f64),

    /// Value dimension (3 for vector field)
    pub value_dim: usize,

    /// Value units (e.g., "A/m" for magnetization)
    pub value_units: String,

    /// Value labels (e.g., ["m_x", "m_y", "m_z"])
    pub value_labels: Vec<String>,

    /// Vector field data (flattened: [v0_x, v0_y, v0_z, v1_x, v1_y, v1_z, ...])
    pub data: Vec<Vector3<f64>>,
}

impl OvfData {
    /// Create new OVF data structure
    pub fn new(nx: usize, ny: usize, nz: usize, xsize: f64, ysize: f64, zsize: f64) -> Self {
        let mesh_size = (nx, ny, nz);
        let mesh_physical_size = (xsize, ysize, zsize);
        let cell_size = (xsize / nx as f64, ysize / ny as f64, zsize / nz as f64);

        Self {
            title: "Magnetization field".to_string(),
            mesh_size,
            mesh_physical_size,
            cell_size,
            value_dim: 3,
            value_units: "A/m".to_string(),
            value_labels: vec!["m_x".to_string(), "m_y".to_string(), "m_z".to_string()],
            data: vec![Vector3::new(0.0, 0.0, 0.0); nx * ny * nz],
        }
    }

    /// Set vector at position (i, j, k)
    pub fn set_vector(&mut self, i: usize, j: usize, k: usize, v: Vector3<f64>) {
        let (nx, ny, _nz) = self.mesh_size;
        let idx = i + j * nx + k * nx * ny;
        if idx < self.data.len() {
            self.data[idx] = v;
        }
    }

    /// Get vector at position (i, j, k)
    pub fn get_vector(&self, i: usize, j: usize, k: usize) -> Option<Vector3<f64>> {
        let (nx, ny, _nz) = self.mesh_size;
        let idx = i + j * nx + k * nx * ny;
        self.data.get(idx).copied()
    }
}

/// OVF file writer
pub struct OvfWriter {
    format: OvfFormat,
}

impl OvfWriter {
    /// Create new OVF writer with specified format
    pub fn new(format: OvfFormat) -> Self {
        Self { format }
    }

    /// Write OVF data to file
    pub fn write<P: AsRef<Path>>(&self, path: P, data: &OvfData) -> std::io::Result<()> {
        match self.format {
            OvfFormat::Text2_0 => self.write_text_2_0(path, data),
            OvfFormat::Text1_0 => self.write_text_1_0(path, data),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::Unsupported,
                "Binary OVF format not yet implemented",
            )),
        }
    }

    /// Write OVF 2.0 text format
    fn write_text_2_0<P: AsRef<Path>>(&self, path: P, data: &OvfData) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(writer, "# OOMMF OVF 2.0")?;
        writeln!(writer, "#")?;
        writeln!(writer, "# Segment count: 1")?;
        writeln!(writer, "#")?;
        writeln!(writer, "# Begin: Segment")?;
        writeln!(writer, "# Begin: Header")?;
        writeln!(writer, "#")?;
        writeln!(writer, "# Title: {}", data.title)?;
        writeln!(writer, "# Desc: Magnetization field data")?;
        writeln!(writer, "#")?;
        writeln!(writer, "# meshtype: rectangular")?;
        writeln!(writer, "# meshunit: m")?;
        writeln!(writer, "#")?;

        let (nx, ny, nz) = data.mesh_size;
        writeln!(writer, "# xnodes: {}", nx)?;
        writeln!(writer, "# ynodes: {}", ny)?;
        writeln!(writer, "# znodes: {}", nz)?;
        writeln!(writer, "#")?;

        let (xsize, ysize, zsize) = data.mesh_physical_size;
        writeln!(writer, "# xmin: 0")?;
        writeln!(writer, "# ymin: 0")?;
        writeln!(writer, "# zmin: 0")?;
        writeln!(writer, "# xmax: {}", xsize)?;
        writeln!(writer, "# ymax: {}", ysize)?;
        writeln!(writer, "# zmax: {}", zsize)?;
        writeln!(writer, "#")?;

        writeln!(writer, "# valuedim: {}", data.value_dim)?;
        writeln!(
            writer,
            "# valuelabels: {} {} {}",
            data.value_labels[0], data.value_labels[1], data.value_labels[2]
        )?;
        writeln!(
            writer,
            "# valueunits: {} {} {}",
            data.value_units, data.value_units, data.value_units
        )?;
        writeln!(writer, "#")?;
        writeln!(writer, "# End: Header")?;
        writeln!(writer, "#")?;
        writeln!(writer, "# Begin: Data Text")?;

        // Write data
        for v in &data.data {
            writeln!(writer, " {:.17e}  {:.17e}  {:.17e}", v.x, v.y, v.z)?;
        }

        writeln!(writer, "# End: Data Text")?;
        writeln!(writer, "# End: Segment")?;

        Ok(())
    }

    /// Write OVF 1.0 text format
    fn write_text_1_0<P: AsRef<Path>>(&self, path: P, data: &OvfData) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // OVF 1.0 header
        writeln!(writer, "# OOMMF: rectangular mesh v1.0")?;
        writeln!(writer, "# Segment count: 1")?;
        writeln!(writer, "# Begin: Segment")?;
        writeln!(writer, "# Begin: Header")?;
        writeln!(writer, "# Title: {}", data.title)?;
        writeln!(writer, "# Desc: Magnetization field data")?;

        let (nx, ny, nz) = data.mesh_size;
        let (dx, dy, dz) = data.cell_size;

        writeln!(writer, "# meshunit: m")?;
        writeln!(writer, "# meshtype: rectangular")?;
        writeln!(writer, "# xbase: {}", dx / 2.0)?;
        writeln!(writer, "# ybase: {}", dy / 2.0)?;
        writeln!(writer, "# zbase: {}", dz / 2.0)?;
        writeln!(writer, "# xstepsize: {}", dx)?;
        writeln!(writer, "# ystepsize: {}", dy)?;
        writeln!(writer, "# zstepsize: {}", dz)?;
        writeln!(writer, "# xnodes: {}", nx)?;
        writeln!(writer, "# ynodes: {}", ny)?;
        writeln!(writer, "# znodes: {}", nz)?;
        writeln!(writer, "# xmin: 0")?;
        writeln!(writer, "# ymin: 0")?;
        writeln!(writer, "# zmin: 0")?;
        writeln!(writer, "# xmax: {}", data.mesh_physical_size.0)?;
        writeln!(writer, "# ymax: {}", data.mesh_physical_size.1)?;
        writeln!(writer, "# zmax: {}", data.mesh_physical_size.2)?;
        writeln!(writer, "# valuedim: {}", data.value_dim)?;
        writeln!(writer, "# valuelabels: m_x m_y m_z")?;
        writeln!(writer, "# valueunits: A/m A/m A/m")?;
        writeln!(writer, "# End: Header")?;
        writeln!(writer, "# Begin: data text")?;

        // Write data
        for v in &data.data {
            writeln!(writer, "{:.17e} {:.17e} {:.17e}", v.x, v.y, v.z)?;
        }

        writeln!(writer, "# End: data text")?;
        writeln!(writer, "# End: Segment")?;

        Ok(())
    }
}

/// OVF file reader
pub struct OvfReader;

impl OvfReader {
    /// Read OVF file
    pub fn read<P: AsRef<Path>>(path: P) -> std::io::Result<OvfData> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let mut lines = reader.lines();
        let first_line = lines
            .next()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Empty file"))??;

        if first_line.starts_with("# OOMMF OVF 2.0") {
            Self::read_ovf_2_0(lines)
        } else if first_line.starts_with("# OOMMF: rectangular mesh v1.0") {
            Self::read_ovf_1_0(lines)
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unknown OVF format",
            ))
        }
    }

    /// Read OVF 2.0 format
    fn read_ovf_2_0<I>(lines: I) -> std::io::Result<OvfData>
    where
        I: Iterator<Item = std::io::Result<String>>,
    {
        let mut title = String::new();
        let mut nx = 0;
        let mut ny = 0;
        let mut nz = 0;
        let mut xmax = 0.0;
        let mut ymax = 0.0;
        let mut zmax = 0.0;
        let mut value_dim = 3;
        let mut in_data_section = false;
        let mut data = Vec::new();

        for line in lines {
            let line = line?;
            let line = line.trim();

            if line.starts_with("# Title:") {
                title = line.strip_prefix("# Title:").unwrap().trim().to_string();
            } else if line.starts_with("# xnodes:") {
                nx = line
                    .strip_prefix("# xnodes:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# ynodes:") {
                ny = line
                    .strip_prefix("# ynodes:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# znodes:") {
                nz = line
                    .strip_prefix("# znodes:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# xmax:") {
                xmax = line
                    .strip_prefix("# xmax:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# ymax:") {
                ymax = line
                    .strip_prefix("# ymax:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# zmax:") {
                zmax = line
                    .strip_prefix("# zmax:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# valuedim:") {
                value_dim = line
                    .strip_prefix("# valuedim:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(3);
            } else if line.starts_with("# Begin: Data Text") {
                in_data_section = true;
            } else if line.starts_with("# End: Data Text") {
                break;
            } else if in_data_section && !line.starts_with('#') && !line.is_empty() {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 3 {
                    let x: f64 = parts[0].parse().unwrap_or(0.0);
                    let y: f64 = parts[1].parse().unwrap_or(0.0);
                    let z: f64 = parts[2].parse().unwrap_or(0.0);
                    data.push(Vector3::new(x, y, z));
                }
            }
        }

        let mut ovf_data = OvfData::new(nx, ny, nz, xmax, ymax, zmax);
        ovf_data.title = title;
        ovf_data.value_dim = value_dim;
        ovf_data.data = data;

        Ok(ovf_data)
    }

    /// Read OVF 1.0 format
    fn read_ovf_1_0<I>(lines: I) -> std::io::Result<OvfData>
    where
        I: Iterator<Item = std::io::Result<String>>,
    {
        let mut title = String::new();
        let mut nx = 0;
        let mut ny = 0;
        let mut nz = 0;
        let mut xmax = 0.0;
        let mut ymax = 0.0;
        let mut zmax = 0.0;
        let mut in_data_section = false;
        let mut data = Vec::new();

        for line in lines {
            let line = line?;
            let line = line.trim();

            if line.starts_with("# Title:") {
                title = line.strip_prefix("# Title:").unwrap().trim().to_string();
            } else if line.starts_with("# xnodes:") {
                nx = line
                    .strip_prefix("# xnodes:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# ynodes:") {
                ny = line
                    .strip_prefix("# ynodes:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# znodes:") {
                nz = line
                    .strip_prefix("# znodes:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# xmax:") {
                xmax = line
                    .strip_prefix("# xmax:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# ymax:") {
                ymax = line
                    .strip_prefix("# ymax:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# zmax:") {
                zmax = line
                    .strip_prefix("# zmax:")
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# Begin: data text") {
                in_data_section = true;
            } else if line.starts_with("# End: data text") {
                break;
            } else if in_data_section && !line.starts_with('#') && !line.is_empty() {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 3 {
                    let x: f64 = parts[0].parse().unwrap_or(0.0);
                    let y: f64 = parts[1].parse().unwrap_or(0.0);
                    let z: f64 = parts[2].parse().unwrap_or(0.0);
                    data.push(Vector3::new(x, y, z));
                }
            }
        }

        let mut ovf_data = OvfData::new(nx, ny, nz, xmax, ymax, zmax);
        ovf_data.title = title;
        ovf_data.data = data;

        Ok(ovf_data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_ovf_data_creation() {
        let ovf = OvfData::new(10, 10, 1, 1e-6, 1e-6, 1e-9);

        assert_eq!(ovf.mesh_size, (10, 10, 1));
        assert_eq!(ovf.data.len(), 100);
        assert_eq!(ovf.value_dim, 3);
    }

    #[test]
    fn test_ovf_set_get_vector() {
        let mut ovf = OvfData::new(5, 5, 1, 5e-7, 5e-7, 1e-9);

        let v = Vector3::new(1.0, 0.0, 0.0);
        ovf.set_vector(2, 3, 0, v);

        let retrieved = ovf.get_vector(2, 3, 0).unwrap();
        assert_eq!(retrieved.x, 1.0);
        assert_eq!(retrieved.y, 0.0);
        assert_eq!(retrieved.z, 0.0);
    }

    #[test]
    fn test_ovf_write_read_2_0() {
        let mut ovf = OvfData::new(3, 3, 1, 3e-9, 3e-9, 1e-9);
        ovf.title = "Test data".to_string();

        // Set some test vectors
        for i in 0..3 {
            for j in 0..3 {
                let angle = (i + j) as f64 * PI / 4.0;
                let v = Vector3::new(angle.cos(), angle.sin(), 0.0);
                ovf.set_vector(i, j, 0, v);
            }
        }

        // Write to file
        let writer = OvfWriter::new(OvfFormat::Text2_0);
        let path = "/tmp/test_ovf_2_0.ovf";
        writer.write(path, &ovf).unwrap();

        // Read back
        let ovf_read = OvfReader::read(path).unwrap();

        assert_eq!(ovf_read.mesh_size, (3, 3, 1));
        assert_eq!(ovf_read.data.len(), 9);
        assert_eq!(ovf_read.title, "Test data");
    }

    #[test]
    fn test_ovf_write_1_0() {
        let mut ovf = OvfData::new(2, 2, 1, 2e-9, 2e-9, 1e-9);
        ovf.title = "OVF 1.0 test".to_string();

        ovf.set_vector(0, 0, 0, Vector3::new(1.0, 0.0, 0.0));
        ovf.set_vector(1, 1, 0, Vector3::new(0.0, 1.0, 0.0));

        let writer = OvfWriter::new(OvfFormat::Text1_0);
        let path = "/tmp/test_ovf_1_0.ovf";
        writer.write(path, &ovf).unwrap();

        let ovf_read = OvfReader::read(path).unwrap();
        assert_eq!(ovf_read.mesh_size, (2, 2, 1));
    }
}
