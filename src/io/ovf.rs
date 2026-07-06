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

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::vector3::Vector3;

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
            OvfFormat::Binary4_2_0 => self.write_binary_2_0(path, data, BinaryWidth::F32),
            OvfFormat::Binary8_2_0 => self.write_binary_2_0(path, data, BinaryWidth::F64),
            OvfFormat::Binary4_1_0 => self.write_binary_1_0(path, data, BinaryWidth::F32),
        }
    }

    /// Write the OVF 2.0 text header.
    ///
    /// Emits every line up to (and including) the metadata block, stopping just
    /// before the `# Begin: Data ...` marker. Both the text and binary OVF 2.0
    /// writers share this helper so that their headers are byte-identical.
    fn write_header_2_0<W: Write>(&self, writer: &mut W, data: &OvfData) -> std::io::Result<()> {
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

        Ok(())
    }

    /// Write the OVF 1.0 text header.
    ///
    /// Emits every line up to (and including) the metadata block, stopping just
    /// before the `# Begin: data ...` marker. Both the text and binary OVF 1.0
    /// writers share this helper so that their headers are byte-identical.
    fn write_header_1_0<W: Write>(&self, writer: &mut W, data: &OvfData) -> std::io::Result<()> {
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

        Ok(())
    }

    /// Write OVF 2.0 text format
    fn write_text_2_0<P: AsRef<Path>>(&self, path: P, data: &OvfData) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        self.write_header_2_0(&mut writer, data)?;
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

        self.write_header_1_0(&mut writer, data)?;
        writeln!(writer, "# Begin: data text")?;

        // Write data
        for v in &data.data {
            writeln!(writer, "{:.17e} {:.17e} {:.17e}", v.x, v.y, v.z)?;
        }

        writeln!(writer, "# End: data text")?;
        writeln!(writer, "# End: Segment")?;

        Ok(())
    }

    /// Write OVF 2.0 binary format (little-endian per the OVF 2.0 specification).
    ///
    /// The text header is identical to [`Self::write_text_2_0`]; only the data
    /// section differs. The binary block consists of a control/check value (used
    /// by readers to verify byte order) followed by `nx*ny*nz*value_dim` field
    /// floats in node order (`idx = i + j*nx + k*nx*ny`).
    fn write_binary_2_0<P: AsRef<Path>>(
        &self,
        path: P,
        data: &OvfData,
        width: BinaryWidth,
    ) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        self.write_header_2_0(&mut writer, data)?;
        let tag = width.bits();
        writeln!(writer, "# Begin: Data Binary {}", tag)?;
        // Flush the buffered text so the raw binary block follows the newline
        // terminating the begin marker exactly.
        writer.flush()?;

        write_binary_block(&mut writer, data, width, Endianness::Little)?;

        writeln!(writer, "# End: Data Binary {}", tag)?;
        writeln!(writer, "# End: Segment")?;

        Ok(())
    }

    /// Write OVF 1.0 binary format (big-endian per the OVF 1.0 specification).
    ///
    /// Structurally identical to [`Self::write_binary_2_0`] but emits the OVF 1.0
    /// header and serializes the binary block in big-endian (MSB-first) order.
    fn write_binary_1_0<P: AsRef<Path>>(
        &self,
        path: P,
        data: &OvfData,
        width: BinaryWidth,
    ) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        self.write_header_1_0(&mut writer, data)?;
        let tag = width.bits();
        writeln!(writer, "# Begin: data binary {}", tag)?;
        writer.flush()?;

        write_binary_block(&mut writer, data, width, Endianness::Big)?;

        writeln!(writer, "# End: data binary {}", tag)?;
        writeln!(writer, "# End: Segment")?;

        Ok(())
    }
}

/// Width of a binary OVF floating-point value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BinaryWidth {
    /// Single precision (4-byte `f32`).
    F32,
    /// Double precision (8-byte `f64`).
    F64,
}

impl BinaryWidth {
    /// Number of bytes occupied by one value of this width.
    fn bytes(self) -> usize {
        match self {
            BinaryWidth::F32 => 4,
            BinaryWidth::F64 => 8,
        }
    }

    /// Tag used in the `Data Binary N` markers (4 or 8).
    fn bits(self) -> usize {
        match self {
            BinaryWidth::F32 => 4,
            BinaryWidth::F64 => 8,
        }
    }
}

/// Byte order of the binary data section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Endianness {
    /// Most-significant byte first (OVF 1.0).
    Big,
    /// Least-significant byte first (OVF 2.0).
    Little,
}

/// Control/check value for `Data Binary 4` (exactly representable as `f32`).
const OVF_CONTROL_F32: f32 = 1234567.0_f32;

/// Control/check value for `Data Binary 8` (exactly representable as `f64`).
const OVF_CONTROL_F64: f64 = 123456789012345.0_f64;

/// Serialize the control value followed by the flattened field data.
///
/// The control value lets readers detect a byte-order mismatch: it is written
/// in the same endianness as the field data, so a reader interpreting the bytes
/// with the wrong order will not recover the expected constant.
fn write_binary_block<W: Write>(
    writer: &mut W,
    data: &OvfData,
    width: BinaryWidth,
    endian: Endianness,
) -> std::io::Result<()> {
    match width {
        BinaryWidth::F32 => {
            let control = OVF_CONTROL_F32;
            writer.write_all(&to_bytes_f32(control, endian))?;
            for v in &data.data {
                writer.write_all(&to_bytes_f32(v.x as f32, endian))?;
                writer.write_all(&to_bytes_f32(v.y as f32, endian))?;
                writer.write_all(&to_bytes_f32(v.z as f32, endian))?;
            }
        },
        BinaryWidth::F64 => {
            let control = OVF_CONTROL_F64;
            writer.write_all(&to_bytes_f64(control, endian))?;
            for v in &data.data {
                writer.write_all(&to_bytes_f64(v.x, endian))?;
                writer.write_all(&to_bytes_f64(v.y, endian))?;
                writer.write_all(&to_bytes_f64(v.z, endian))?;
            }
        },
    }
    // Terminate the raw binary block with a newline before the end marker.
    writer.write_all(b"\n")?;
    Ok(())
}

/// Encode an `f32` in the requested byte order.
fn to_bytes_f32(value: f32, endian: Endianness) -> [u8; 4] {
    match endian {
        Endianness::Big => value.to_be_bytes(),
        Endianness::Little => value.to_le_bytes(),
    }
}

/// Encode an `f64` in the requested byte order.
fn to_bytes_f64(value: f64, endian: Endianness) -> [u8; 8] {
    match endian {
        Endianness::Big => value.to_be_bytes(),
        Endianness::Little => value.to_le_bytes(),
    }
}

/// Decode an `f32` from a 4-byte slice in the requested byte order.
fn read_f32(bytes: &[u8], endian: Endianness) -> f32 {
    let mut buf = [0u8; 4];
    buf.copy_from_slice(&bytes[..4]);
    match endian {
        Endianness::Big => f32::from_be_bytes(buf),
        Endianness::Little => f32::from_le_bytes(buf),
    }
}

/// Decode an `f64` from an 8-byte slice in the requested byte order.
fn read_f64(bytes: &[u8], endian: Endianness) -> f64 {
    let mut buf = [0u8; 8];
    buf.copy_from_slice(&bytes[..8]);
    match endian {
        Endianness::Big => f64::from_be_bytes(buf),
        Endianness::Little => f64::from_le_bytes(buf),
    }
}

/// OVF specification version, used to select the binary byte order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OvfVersion {
    /// OVF 1.0 (binary data is big-endian).
    V1_0,
    /// OVF 2.0 (binary data is little-endian).
    V2_0,
}

/// Encoding of the OVF data section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DataEncoding {
    /// ASCII text records.
    Text,
    /// Raw binary records of the given floating-point width.
    Binary(BinaryWidth),
}

/// Ensure that `len` bytes are available in `raw` starting at `offset`.
///
/// Returns an [`std::io::ErrorKind::InvalidData`] error on a truncated file so
/// that decoding never reads past the end of the buffer.
fn ensure_available(raw: &[u8], offset: usize, len: usize) -> std::io::Result<()> {
    let fits = match offset.checked_add(len) {
        Some(end) => end <= raw.len(),
        None => false,
    };
    if !fits {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Truncated OVF binary data section",
        ));
    }
    Ok(())
}

/// OVF file reader
pub struct OvfReader;

impl OvfReader {
    /// Read an OVF file, auto-detecting the version (1.0/2.0) and the data
    /// encoding (text or binary).
    ///
    /// The encoding is determined from the `# Begin: Data ...` marker in the
    /// (always textual) header, matched case-insensitively. Text files are
    /// dispatched to the line-based reader; binary files are dispatched to the
    /// byte-oriented reader, which selects big-endian for OVF 1.0 and
    /// little-endian for OVF 2.0 per the respective specifications.
    pub fn read<P: AsRef<Path>>(path: P) -> std::io::Result<OvfData> {
        let path = path.as_ref();
        let (version, encoding) = Self::detect_format(path)?;

        match encoding {
            DataEncoding::Text => Self::read_text(path),
            DataEncoding::Binary(width) => Self::read_binary(path, version, width),
        }
    }

    /// Detect the OVF version and data encoding by scanning the header lines.
    ///
    /// Reads the file with a [`BufReader`], inspecting lines until the
    /// `# Begin: Data ...` marker is found. The first line yields the version;
    /// the begin marker yields the encoding (text, binary 4, or binary 8).
    fn detect_format(path: &Path) -> std::io::Result<(OvfVersion, DataEncoding)> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let mut version: Option<OvfVersion> = None;
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();

            if version.is_none() {
                if trimmed.starts_with("# OOMMF OVF 2.0") {
                    version = Some(OvfVersion::V2_0);
                } else if trimmed.starts_with("# OOMMF: rectangular mesh v1.0") {
                    version = Some(OvfVersion::V1_0);
                }
            }

            let lower = trimmed.to_ascii_lowercase();
            if let Some(rest) = lower.strip_prefix("# begin: ") {
                if let Some(kind) = rest.strip_prefix("data ") {
                    let kind = kind.trim();
                    let encoding = if kind.starts_with("text") {
                        DataEncoding::Text
                    } else if kind.starts_with("binary 4") {
                        DataEncoding::Binary(BinaryWidth::F32)
                    } else if kind.starts_with("binary 8") {
                        DataEncoding::Binary(BinaryWidth::F64)
                    } else {
                        return Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "Unknown OVF data encoding",
                        ));
                    };
                    let version = version.ok_or_else(|| {
                        std::io::Error::new(std::io::ErrorKind::InvalidData, "Unknown OVF format")
                    })?;
                    return Ok((version, encoding));
                }
            }
        }

        Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Missing OVF data section",
        ))
    }

    /// Read a text-encoded OVF file using the line-based parsers.
    fn read_text(path: &Path) -> std::io::Result<OvfData> {
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

    /// Read a binary-encoded OVF file.
    ///
    /// The whole file is loaded into memory; the textual header is parsed line
    /// by line for mesh metadata, then the raw binary block is decoded starting
    /// at the byte immediately after the newline terminating the
    /// `# Begin: Data Binary N` marker. The control value is validated against
    /// the width-specific check constant to catch byte-order mismatches and
    /// corruption.
    fn read_binary(
        path: &Path,
        version: OvfVersion,
        width: BinaryWidth,
    ) -> std::io::Result<OvfData> {
        let endian = match version {
            OvfVersion::V1_0 => Endianness::Big,
            OvfVersion::V2_0 => Endianness::Little,
        };

        let raw = std::fs::read(path)?;

        // Parse the textual header, accumulating metadata and locating the byte
        // offset of the data block (just past the begin-marker newline).
        let mut title = String::new();
        let mut nx = 0usize;
        let mut ny = 0usize;
        let mut nz = 0usize;
        let mut xmax = 0.0f64;
        let mut ymax = 0.0f64;
        let mut zmax = 0.0f64;
        let mut value_dim = 3usize;
        let mut data_offset: Option<usize> = None;

        let mut cursor = 0usize;
        while cursor < raw.len() {
            // Locate the next newline; the header is guaranteed to terminate the
            // begin marker with one before the binary block starts.
            let newline = match raw[cursor..].iter().position(|&b| b == b'\n') {
                Some(rel) => cursor + rel,
                None => break,
            };
            let line = String::from_utf8_lossy(&raw[cursor..newline]);
            let line = line.trim();

            if line.starts_with("# Title:") {
                if let Some(rest) = line.strip_prefix("# Title:") {
                    title = rest.trim().to_string();
                }
            } else if line.starts_with("# xnodes:") {
                nx = line
                    .strip_prefix("# xnodes:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# ynodes:") {
                ny = line
                    .strip_prefix("# ynodes:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# znodes:") {
                nz = line
                    .strip_prefix("# znodes:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# xmax:") {
                xmax = line
                    .strip_prefix("# xmax:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# ymax:") {
                ymax = line
                    .strip_prefix("# ymax:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# zmax:") {
                zmax = line
                    .strip_prefix("# zmax:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# valuedim:") {
                value_dim = line
                    .strip_prefix("# valuedim:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(3);
            } else if line
                .to_ascii_lowercase()
                .starts_with("# begin: data binary")
            {
                // The binary block begins immediately after this line's newline.
                data_offset = Some(newline + 1);
                break;
            }

            cursor = newline + 1;
        }

        let mut offset = data_offset.ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Missing OVF binary data section",
            )
        })?;

        let value_bytes = width.bytes();
        let comp_count = value_dim.max(1);
        let node_count = nx * ny * nz;

        // Validate the control value, which is written ahead of the field data.
        let control_ok = match width {
            BinaryWidth::F32 => {
                ensure_available(&raw, offset, value_bytes)?;
                let control = read_f32(&raw[offset..], endian);
                control == OVF_CONTROL_F32
            },
            BinaryWidth::F64 => {
                ensure_available(&raw, offset, value_bytes)?;
                let control = read_f64(&raw[offset..], endian);
                control == OVF_CONTROL_F64
            },
        };
        if !control_ok {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "OVF binary control number mismatch (wrong byte order or corrupt file)",
            ));
        }
        offset += value_bytes;

        // Decode `node_count` records of `comp_count` components each, keeping
        // the first three for the stored `Vector3` (value_dim is 3 here).
        let mut data = Vec::with_capacity(node_count);
        for _ in 0..node_count {
            let mut comps = [0.0f64; 3];
            for (c, comp) in comps.iter_mut().enumerate().take(comp_count.min(3)) {
                ensure_available(&raw, offset + c * value_bytes, value_bytes)?;
                *comp = match width {
                    BinaryWidth::F32 => read_f32(&raw[offset + c * value_bytes..], endian) as f64,
                    BinaryWidth::F64 => read_f64(&raw[offset + c * value_bytes..], endian),
                };
            }
            // Skip any trailing components beyond the three we retain.
            ensure_available(&raw, offset, comp_count * value_bytes)?;
            offset += comp_count * value_bytes;
            data.push(Vector3::new(comps[0], comps[1], comps[2]));
        }

        let mut ovf_data = OvfData::new(nx, ny, nz, xmax, ymax, zmax);
        ovf_data.title = title;
        ovf_data.value_dim = value_dim;
        ovf_data.data = data;

        Ok(ovf_data)
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
                if let Some(rest) = line.strip_prefix("# Title:") {
                    title = rest.trim().to_string();
                }
            } else if line.starts_with("# xnodes:") {
                nx = line
                    .strip_prefix("# xnodes:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# ynodes:") {
                ny = line
                    .strip_prefix("# ynodes:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# znodes:") {
                nz = line
                    .strip_prefix("# znodes:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# xmax:") {
                xmax = line
                    .strip_prefix("# xmax:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# ymax:") {
                ymax = line
                    .strip_prefix("# ymax:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# zmax:") {
                zmax = line
                    .strip_prefix("# zmax:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# valuedim:") {
                value_dim = line
                    .strip_prefix("# valuedim:")
                    .unwrap_or("")
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
                if let Some(rest) = line.strip_prefix("# Title:") {
                    title = rest.trim().to_string();
                }
            } else if line.starts_with("# xnodes:") {
                nx = line
                    .strip_prefix("# xnodes:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# ynodes:") {
                ny = line
                    .strip_prefix("# ynodes:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# znodes:") {
                nz = line
                    .strip_prefix("# znodes:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0);
            } else if line.starts_with("# xmax:") {
                xmax = line
                    .strip_prefix("# xmax:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# ymax:") {
                ymax = line
                    .strip_prefix("# ymax:")
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            } else if line.starts_with("# zmax:") {
                zmax = line
                    .strip_prefix("# zmax:")
                    .unwrap_or("")
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
    use std::f64::consts::PI;

    use super::*;

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

        let retrieved = ovf
            .get_vector(2, 3, 0)
            .expect("vector at (2,3,0) should exist");
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
        let mut path = std::env::temp_dir();
        path.push("test_ovf_2_0.ovf");
        writer
            .write(&path, &ovf)
            .expect("OVF 2.0 write should succeed");

        // Read back
        let ovf_read = OvfReader::read(&path).expect("OVF 2.0 read should succeed");

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
        let mut path = std::env::temp_dir();
        path.push("test_ovf_1_0.ovf");
        writer
            .write(&path, &ovf)
            .expect("OVF 1.0 write should succeed");

        let ovf_read = OvfReader::read(&path).expect("OVF 1.0 read should succeed");
        assert_eq!(ovf_read.mesh_size, (2, 2, 1));
    }

    /// Build a distinctive 3x2x2 field with varied components, including
    /// negatives and small magnitudes, to exercise binary serialization.
    fn make_test_field() -> OvfData {
        let mut ovf = OvfData::new(3, 2, 2, 3e-9, 2e-9, 2e-9);
        ovf.title = "Binary OVF test".to_string();

        let (nx, ny, nz) = ovf.mesh_size;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = (i + j * nx + k * nx * ny) as f64;
                    let v = Vector3::new(
                        (idx + 1.0) * 1.0e3,
                        -(idx + 0.5) * 2.5,
                        (idx - 4.0) * 1.0e-7,
                    );
                    ovf.set_vector(i, j, k, v);
                }
            }
        }
        ovf
    }

    #[test]
    fn test_ovf_binary_8_2_0_roundtrip() {
        let ovf = make_test_field();

        let writer = OvfWriter::new(OvfFormat::Binary8_2_0);
        let mut path = std::env::temp_dir();
        path.push("test_ovf_binary_8_2_0.ovf");
        writer
            .write(&path, &ovf)
            .expect("OVF binary 8 write should succeed");

        let read = OvfReader::read(&path).expect("OVF binary 8 read should succeed");

        assert_eq!(read.mesh_size, ovf.mesh_size);
        assert_eq!(read.title, ovf.title);
        assert_eq!(read.value_dim, ovf.value_dim);
        assert_eq!(read.data.len(), ovf.data.len());

        // Binary 8 is lossless: exact f64 equality must hold.
        for (a, b) in read.data.iter().zip(ovf.data.iter()) {
            assert_eq!(a.x, b.x);
            assert_eq!(a.y, b.y);
            assert_eq!(a.z, b.z);
        }
    }

    #[test]
    fn test_ovf_binary_4_2_0_roundtrip() {
        let ovf = make_test_field();

        let writer = OvfWriter::new(OvfFormat::Binary4_2_0);
        let mut path = std::env::temp_dir();
        path.push("test_ovf_binary_4_2_0.ovf");
        writer
            .write(&path, &ovf)
            .expect("OVF binary 4 write should succeed");

        let read = OvfReader::read(&path).expect("OVF binary 4 read should succeed");

        assert_eq!(read.mesh_size, ovf.mesh_size);
        assert_eq!(read.title, ovf.title);
        assert_eq!(read.value_dim, ovf.value_dim);

        // Binary 4 truncates to f32: compare against the original cast through f32.
        for (a, b) in read.data.iter().zip(ovf.data.iter()) {
            assert_eq!(a.x, (b.x as f32) as f64);
            assert_eq!(a.y, (b.y as f32) as f64);
            assert_eq!(a.z, (b.z as f32) as f64);
        }
    }

    #[test]
    fn test_ovf_binary_4_1_0_roundtrip_big_endian() {
        let ovf = make_test_field();

        let writer = OvfWriter::new(OvfFormat::Binary4_1_0);
        let mut path = std::env::temp_dir();
        path.push("test_ovf_binary_4_1_0.ovf");
        writer
            .write(&path, &ovf)
            .expect("OVF binary 4 (1.0, big-endian) write should succeed");

        // The reader must select big-endian from the OVF 1.0 version line.
        let read = OvfReader::read(&path).expect("OVF binary 4 (1.0) read should succeed");

        assert_eq!(read.mesh_size, ovf.mesh_size);
        assert_eq!(read.data.len(), ovf.data.len());

        for (a, b) in read.data.iter().zip(ovf.data.iter()) {
            let rel = |x: f64, y: f64| (x - y).abs() <= 1e-6 * y.abs().max(1.0);
            assert!(rel(a.x, b.x), "x mismatch: {} vs {}", a.x, b.x);
            assert!(rel(a.y, b.y), "y mismatch: {} vs {}", a.y, b.y);
            assert!(rel(a.z, b.z), "z mismatch: {} vs {}", a.z, b.z);
        }
    }

    #[test]
    fn test_ovf_binary_control_value_mismatch_errors() {
        use std::io::Write as _;

        // Construct a syntactically valid OVF 2.0 binary header but seed the
        // control value with the wrong number; the reader must reject it.
        let mut path = std::env::temp_dir();
        path.push("test_ovf_binary_bad_control.ovf");

        let nx = 1usize;
        let ny = 1usize;
        let nz = 1usize;

        {
            let file = File::create(&path).expect("create bad-control file");
            let mut writer = BufWriter::new(file);
            writeln!(writer, "# OOMMF OVF 2.0").expect("write header");
            writeln!(writer, "# Segment count: 1").expect("write header");
            writeln!(writer, "# Begin: Segment").expect("write header");
            writeln!(writer, "# Begin: Header").expect("write header");
            writeln!(writer, "# Title: bad control").expect("write header");
            writeln!(writer, "# xnodes: {}", nx).expect("write header");
            writeln!(writer, "# ynodes: {}", ny).expect("write header");
            writeln!(writer, "# znodes: {}", nz).expect("write header");
            writeln!(writer, "# xmax: 1e-9").expect("write header");
            writeln!(writer, "# ymax: 1e-9").expect("write header");
            writeln!(writer, "# zmax: 1e-9").expect("write header");
            writeln!(writer, "# valuedim: 3").expect("write header");
            writeln!(writer, "# End: Header").expect("write header");
            writeln!(writer, "# Begin: Data Binary 4").expect("write marker");
            writer.flush().expect("flush header");

            // Deliberately wrong control value (little-endian, as OVF 2.0 expects).
            writer
                .write_all(&9999.0_f32.to_le_bytes())
                .expect("write bad control");
            // One node of data so the field block is otherwise well-formed.
            for _ in 0..(nx * ny * nz * 3) {
                writer
                    .write_all(&0.0_f32.to_le_bytes())
                    .expect("write data");
            }
            writer.write_all(b"\n").expect("write block terminator");
            writeln!(writer, "# End: Data Binary 4").expect("write end marker");
            writeln!(writer, "# End: Segment").expect("write end marker");
        }

        let err = OvfReader::read(&path).expect_err("bad control value must error");
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    }

    #[test]
    fn test_ovf_binary_2_0_header_matches_text() {
        let ovf = make_test_field();

        // Write the same data as binary-2.0 and text-2.0 and confirm both yield
        // identical metadata round-trips.
        let mut bin_path = std::env::temp_dir();
        bin_path.push("test_ovf_binary_2_0_meta.ovf");
        OvfWriter::new(OvfFormat::Binary4_2_0)
            .write(&bin_path, &ovf)
            .expect("binary write should succeed");

        let mut text_path = std::env::temp_dir();
        text_path.push("test_ovf_text_2_0_meta.ovf");
        OvfWriter::new(OvfFormat::Text2_0)
            .write(&text_path, &ovf)
            .expect("text write should succeed");

        let bin_read = OvfReader::read(&bin_path).expect("binary read should succeed");
        let text_read = OvfReader::read(&text_path).expect("text read should succeed");

        assert_eq!(bin_read.mesh_size, text_read.mesh_size);
        assert_eq!(bin_read.title, text_read.title);
        assert_eq!(bin_read.value_dim, text_read.value_dim);
    }
}
