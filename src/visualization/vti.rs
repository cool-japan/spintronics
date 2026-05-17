//! VTK ImageData (VTI) Export — binary appended data with base64 encoding
//!
//! Writes `.vti` files in VTK XML ImageData format. Binary data is encoded
//! in base64 (little-endian f32, with a 4-byte uint32 byte-count header).
//! This format is directly loadable in ParaView and VisIt.
//!
//! ## Format
//!
//! - Regular structured grids with uniform or per-axis spacing
//! - Each DataArray stores a 4-byte LE uint32 byte count, then the f32 values
//! - Multiple fields (vector + scalar) in one file via `write_multi_field`
//!
//! ## Example
//!
//! ```rust
//! use spintronics::visualization::vti::VtiWriter;
//! use spintronics::Vector3;
//!
//! let writer = VtiWriter::with_uniform_spacing(4, 4, 1, 1e-9);
//! let spins = vec![Vector3::new(1.0, 0.0, 0.0); 16];
//! // writer.write_vector_field(path, "spin", &spins).unwrap();
//! ```

#![cfg(feature = "vti")]

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use base64::prelude::BASE64_STANDARD;
use base64::Engine;

use crate::error::{invalid_param, Error};
use crate::vector3::Vector3;

/// Typed field data: either a vector field or a scalar field.
pub enum FieldData {
    /// 3-component vector field (stored as f32 triplets)
    Vector(Vec<Vector3<f64>>),
    /// Single-component scalar field
    Scalar(Vec<f64>),
}

/// VTK ImageData writer (`.vti`).
///
/// Generates XML-format VTK ImageData files with binary appended data
/// encoded in base64. Each data array contains a 4-byte LE uint32
/// byte count followed by the f32 data values.
pub struct VtiWriter {
    /// Number of points in the x direction
    pub nx: usize,
    /// Number of points in the y direction
    pub ny: usize,
    /// Number of points in the z direction
    pub nz: usize,
    /// Origin (x0, y0, z0) in physical units
    pub origin: [f64; 3],
    /// Grid spacing (dx, dy, dz) in physical units
    pub spacing: [f64; 3],
}

impl VtiWriter {
    /// Create a new `VtiWriter` with explicit origin and per-axis spacing.
    ///
    /// Returns an error if any dimension is zero or any spacing is
    /// non-positive.
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        origin: [f64; 3],
        spacing: [f64; 3],
    ) -> Result<Self, Error> {
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(invalid_param("dims", "all dimensions must be >= 1"));
        }
        if spacing.iter().any(|&s| s <= 0.0) {
            return Err(invalid_param(
                "spacing",
                "all spacing values must be positive",
            ));
        }
        Ok(Self {
            nx,
            ny,
            nz,
            origin,
            spacing,
        })
    }

    /// Create a `VtiWriter` with origin at zero and uniform isotropic spacing.
    pub fn with_uniform_spacing(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            origin: [0.0, 0.0, 0.0],
            spacing: [dx, dx, dx],
        }
    }

    /// Write a vector field to a `.vti` file.
    ///
    /// `data` must have exactly `nx * ny * nz` elements.
    pub fn write_vector_field(
        &self,
        path: &Path,
        name: &str,
        data: &[Vector3<f64>],
    ) -> Result<(), Error> {
        let expected = self.nx * self.ny * self.nz;
        if data.len() != expected {
            return Err(invalid_param(
                "data",
                &format!("expected {} points, got {}", expected, data.len()),
            ));
        }
        let mut fields = HashMap::new();
        let owned: Vec<Vector3<f64>> = data.to_vec();
        fields.insert(name, FieldData::Vector(owned));
        self.write_multi_field(path, &fields)
    }

    /// Write a scalar field to a `.vti` file.
    ///
    /// `data` must have exactly `nx * ny * nz` elements.
    pub fn write_scalar_field(&self, path: &Path, name: &str, data: &[f64]) -> Result<(), Error> {
        let expected = self.nx * self.ny * self.nz;
        if data.len() != expected {
            return Err(invalid_param(
                "data",
                &format!("expected {} points, got {}", expected, data.len()),
            ));
        }
        let mut fields = HashMap::new();
        let owned: Vec<f64> = data.to_vec();
        fields.insert(name, FieldData::Scalar(owned));
        self.write_multi_field(path, &fields)
    }

    /// Write multiple named fields (vector and/or scalar) to a single `.vti` file.
    pub fn write_multi_field(
        &self,
        path: &Path,
        fields: &HashMap<&str, FieldData>,
    ) -> Result<(), Error> {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let [ox, oy, oz] = self.origin;
        let [dx, dy, dz] = self.spacing;

        // Determine default scalar/vector names for the PointData element
        let first_scalar = fields
            .iter()
            .find(|(_, v)| matches!(v, FieldData::Scalar(_)))
            .map(|(k, _)| *k);
        let first_vector = fields
            .iter()
            .find(|(_, v)| matches!(v, FieldData::Vector(_)))
            .map(|(k, _)| *k);

        let pointdata_attrs = match (first_scalar, first_vector) {
            (Some(s), Some(v)) => format!(" Scalars=\"{}\" Vectors=\"{}\"", s, v),
            (Some(s), None) => format!(" Scalars=\"{}\"", s),
            (None, Some(v)) => format!(" Vectors=\"{}\"", v),
            (None, None) => String::new(),
        };

        writeln!(w, "<?xml version=\"1.0\"?>")?;
        writeln!(
            w,
            "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"
        )?;
        writeln!(
            w,
            "  <ImageData WholeExtent=\"0 {} 0 {} 0 {}\" \
             Origin=\"{} {} {}\" Spacing=\"{} {} {}\">",
            nx - 1,
            ny - 1,
            nz - 1,
            ox,
            oy,
            oz,
            dx,
            dy,
            dz
        )?;
        writeln!(
            w,
            "    <Piece Extent=\"0 {} 0 {} 0 {}\">",
            nx - 1,
            ny - 1,
            nz - 1
        )?;
        writeln!(w, "      <PointData{}>", pointdata_attrs)?;

        // Write each field as a DataArray
        // Sort keys for deterministic output order
        let mut sorted_keys: Vec<&&str> = fields.keys().collect();
        sorted_keys.sort();
        for key in sorted_keys {
            let field = &fields[key];
            self.write_data_array(&mut w, key, field)?;
        }

        writeln!(w, "      </PointData>")?;
        writeln!(w, "      <CellData/>")?;
        writeln!(w, "    </Piece>")?;
        writeln!(w, "  </ImageData>")?;
        writeln!(w, "</VTKFile>")?;

        Ok(())
    }

    /// Encode one field as a `<DataArray>` with base64-encoded binary payload.
    ///
    /// Binary layout: `[byte_count: u32 LE][f32 values: LE f32 each]`
    fn write_data_array(
        &self,
        w: &mut BufWriter<File>,
        name: &str,
        field: &FieldData,
    ) -> Result<(), Error> {
        let (ncomp, encoded) = match field {
            FieldData::Vector(vecs) => {
                let encoded = encode_vector_field_base64(vecs);
                (3usize, encoded)
            },
            FieldData::Scalar(vals) => {
                let encoded = encode_scalar_field_base64(vals);
                (1usize, encoded)
            },
        };

        writeln!(
            w,
            "        <DataArray type=\"Float32\" Name=\"{}\" \
             NumberOfComponents=\"{}\" format=\"binary\">",
            name, ncomp
        )?;
        writeln!(w, "          {}", encoded)?;
        writeln!(w, "        </DataArray>")?;
        Ok(())
    }
}

/// Encode a vector field as base64 binary: u32 LE byte count + f32 LE components.
fn encode_vector_field_base64(vecs: &[Vector3<f64>]) -> String {
    let n_floats = vecs.len() * 3;
    let byte_count = (n_floats * 4) as u32;

    let mut buf: Vec<u8> = Vec::with_capacity(4 + n_floats * 4);
    buf.extend_from_slice(&byte_count.to_le_bytes());
    for v in vecs {
        buf.extend_from_slice(&(v.x as f32).to_bits().to_le_bytes());
        buf.extend_from_slice(&(v.y as f32).to_bits().to_le_bytes());
        buf.extend_from_slice(&(v.z as f32).to_bits().to_le_bytes());
    }
    BASE64_STANDARD.encode(&buf)
}

/// Encode a scalar field as base64 binary: u32 LE byte count + f32 LE values.
fn encode_scalar_field_base64(vals: &[f64]) -> String {
    let n_floats = vals.len();
    let byte_count = (n_floats * 4) as u32;

    let mut buf: Vec<u8> = Vec::with_capacity(4 + n_floats * 4);
    buf.extend_from_slice(&byte_count.to_le_bytes());
    for &v in vals {
        buf.extend_from_slice(&(v as f32).to_bits().to_le_bytes());
    }
    BASE64_STANDARD.encode(&buf)
}

/// Decode a base64-encoded VTI binary DataArray back to f32 values.
///
/// Returns `(byte_count_header, values)`.
pub fn decode_vti_base64(encoded: &str) -> Result<(u32, Vec<f32>), Error> {
    let raw = BASE64_STANDARD
        .decode(encoded.trim())
        .map_err(|e| invalid_param("base64", &e.to_string()))?;

    if raw.len() < 4 {
        return Err(invalid_param("base64", "buffer too short for header"));
    }

    let byte_count = u32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]]);
    let expected_bytes = byte_count as usize;

    if raw.len() < 4 + expected_bytes {
        return Err(invalid_param(
            "base64",
            "buffer shorter than declared byte count",
        ));
    }

    let n_floats = expected_bytes / 4;
    let mut values = Vec::with_capacity(n_floats);
    for i in 0..n_floats {
        let offset = 4 + i * 4;
        let bits = u32::from_le_bytes([
            raw[offset],
            raw[offset + 1],
            raw[offset + 2],
            raw[offset + 3],
        ]);
        values.push(f32::from_bits(bits));
    }
    Ok((byte_count, values))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn temp_path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(name)
    }

    // Read back the contents of a VTI file as a string
    fn read_file(path: &Path) -> String {
        fs::read_to_string(path).expect("should read temp file")
    }

    #[test]
    fn test_with_uniform_spacing() {
        let w = VtiWriter::with_uniform_spacing(4, 4, 2, 1e-9);
        assert_eq!(w.nx, 4);
        assert_eq!(w.ny, 4);
        assert_eq!(w.nz, 2);
        assert_eq!(w.spacing, [1e-9, 1e-9, 1e-9]);
        assert_eq!(w.origin, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_new_invalid_dims() {
        let result = VtiWriter::new(0, 2, 2, [0.0; 3], [1.0; 3]);
        assert!(result.is_err());
    }

    #[test]
    fn test_new_invalid_spacing() {
        let result = VtiWriter::new(2, 2, 2, [0.0; 3], [1.0, -1.0, 1.0]);
        assert!(result.is_err());
    }

    #[test]
    fn test_vector_field_size_mismatch() {
        let w = VtiWriter::with_uniform_spacing(2, 2, 2, 1.0);
        let path = temp_path("vti_size_mismatch.vti");
        let data = vec![Vector3::new(1.0, 0.0, 0.0); 3]; // wrong size
        let result = w.write_vector_field(&path, "spin", &data);
        assert!(result.is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_scalar_field_size_mismatch() {
        let w = VtiWriter::with_uniform_spacing(2, 2, 2, 1.0);
        let path = temp_path("vti_scalar_mismatch.vti");
        let data = vec![0.0f64; 3];
        let result = w.write_scalar_field(&path, "energy", &data);
        assert!(result.is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_write_vector_field_creates_file() {
        let w = VtiWriter::with_uniform_spacing(2, 2, 2, 1e-9);
        let path = temp_path("vti_vector_test.vti");
        let data: Vec<Vector3<f64>> = (0..8)
            .map(|i| Vector3::new(i as f64 * 0.1, 0.0, 1.0))
            .collect();
        w.write_vector_field(&path, "magnetization", &data)
            .expect("write should succeed");
        assert!(path.exists());
        let content = read_file(&path);
        assert!(content.contains("ImageData"));
        assert!(content.contains("magnetization"));
        assert!(content.contains("Float32"));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_write_scalar_field_creates_file() {
        let w = VtiWriter::with_uniform_spacing(3, 3, 1, 0.5e-9);
        let path = temp_path("vti_scalar_test.vti");
        let data: Vec<f64> = (0..9).map(|i| i as f64).collect();
        w.write_scalar_field(&path, "energy", &data)
            .expect("write should succeed");
        assert!(path.exists());
        let content = read_file(&path);
        assert!(content.contains("energy"));
        assert!(content.contains("NumberOfComponents=\"1\""));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_binary_roundtrip_vector() {
        let w = VtiWriter::with_uniform_spacing(2, 2, 2, 1.0);
        let path = temp_path("vti_roundtrip_vec.vti");

        let spins: Vec<Vector3<f64>> = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.5, 0.5, 0.0),
            Vector3::new(0.0, 0.5, 0.5),
            Vector3::new(0.5, 0.0, 0.5),
            Vector3::new(0.0, 0.0, -1.0),
        ];

        w.write_vector_field(&path, "spin", &spins)
            .expect("write should succeed");

        // Parse the file and extract the base64 DataArray
        let content = read_file(&path);
        let start = content
            .find("<DataArray")
            .expect("DataArray tag must exist");
        let end = content
            .find("</DataArray>")
            .expect("DataArray close must exist");
        let section = &content[start..end];

        // Extract base64 content between the closing '>' of the open tag and '</DataArray>'
        let tag_end = section.find('>').expect("closing >") + 1;
        let encoded = section[tag_end..].trim();

        let (byte_count, floats) = decode_vti_base64(encoded).expect("decode should succeed");

        // byte_count must match: 8 vectors * 3 components * 4 bytes/f32 = 96
        assert_eq!(byte_count, 96);
        assert_eq!(floats.len(), 24);

        // Check x-component of first spin
        assert!((floats[0] - 1.0f32).abs() < 1e-6);
        // Check y-component of second spin
        assert!((floats[4] - 1.0f32).abs() < 1e-6);
        // Check z-component of third spin
        assert!((floats[8] - 1.0f32).abs() < 1e-6);

        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_binary_header_byte_count_scalar() {
        let w = VtiWriter::with_uniform_spacing(2, 2, 1, 1.0);
        let path = temp_path("vti_scalar_bytecount.vti");
        let data = vec![1.0f64, 2.0, 3.0, 4.0];

        w.write_scalar_field(&path, "energy", &data)
            .expect("write should succeed");

        let content = read_file(&path);
        let tag_end = content.find('>').unwrap();
        // Find the binary DataArray content
        let da_start = content.find("<DataArray").unwrap();
        let da_section = &content[da_start..];
        let inner_start = da_section.find('>').unwrap() + 1;
        let inner_end = da_section.find("</DataArray>").unwrap();
        let encoded = da_section[inner_start..inner_end].trim();

        let (byte_count, floats) = decode_vti_base64(encoded).expect("decode must succeed");
        // 4 scalars * 4 bytes = 16
        assert_eq!(byte_count, 16);
        assert_eq!(floats.len(), 4);
        assert!((floats[0] - 1.0f32).abs() < 1e-6);
        assert!((floats[3] - 4.0f32).abs() < 1e-6);

        let _ = fs::remove_file(&path);
        let _ = tag_end; // suppress unused warning
    }

    #[test]
    fn test_xml_structure_extent() {
        let w = VtiWriter::with_uniform_spacing(5, 4, 3, 2.0);
        let path = temp_path("vti_extent.vti");
        let n = 5 * 4 * 3;
        let data: Vec<Vector3<f64>> = (0..n).map(|_| Vector3::new(0.0, 0.0, 1.0)).collect();

        w.write_vector_field(&path, "m", &data)
            .expect("write should succeed");

        let content = read_file(&path);
        // WholeExtent should be "0 4 0 3 0 2"
        assert!(content.contains("WholeExtent=\"0 4 0 3 0 2\""));
        assert!(content.contains("Spacing=\"2 2 2\""));

        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_multi_field() {
        let w = VtiWriter::with_uniform_spacing(2, 2, 1, 1.0);
        let path = temp_path("vti_multi.vti");
        let mut fields = HashMap::new();
        let vecs: Vec<Vector3<f64>> = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(-1.0, 0.0, 0.0),
        ];
        let scalars = vec![0.1f64, 0.2, 0.3, 0.4];
        fields.insert("spin", FieldData::Vector(vecs));
        fields.insert("energy", FieldData::Scalar(scalars));
        w.write_multi_field(&path, &fields)
            .expect("write should succeed");
        let content = read_file(&path);
        assert!(content.contains("spin"));
        assert!(content.contains("energy"));
        assert!(content.contains("NumberOfComponents=\"3\""));
        assert!(content.contains("NumberOfComponents=\"1\""));
        let _ = fs::remove_file(&path);
    }
}
