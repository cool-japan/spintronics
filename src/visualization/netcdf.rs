//! Pure-Rust NetCDF3 Classic Format Writer/Reader
//!
//! Implements the NetCDF3 Classic (CDF-1) binary format as defined in
//! the UNIDATA Network Common Data Form specification. The format uses
//! XDR (big-endian) encoding for all integers, with padding to 4-byte
//! boundaries.
//!
//! ## NetCDF3 Classic binary layout
//!
//! ```text
//! magic       : b"CDF\x01"
//! numrecs     : i32  (0 for fixed-size files)
//! dim_list    : tag=0x0000000A, ndims i32, [name_len i32, name bytes (padded), size i32]...
//! att_list    : tag=0x0000000C, natts i32, [name, type=NC_CHAR, len, chars (padded)]...
//! var_list    : tag=0x0000000E, nvars i32, [name, ndimids, atts, type, vsize, begin]...
//! data        : each variable's values packed in XDR order (big-endian)
//! ```
//!
//! ## CF conventions
//!
//! `write_vector_field_cf` creates `{name}_x`, `{name}_y`, `{name}_z`
//! variables with CF-standard `units = "1"` attributes, ready for
//! import into NCO, Xarray, CDO, etc.
//!
//! ## Example
//!
//! ```rust
//! use spintronics::visualization::netcdf::{NetCdfWriter, NetCdfData};
//! let mut nc = NetCdfWriter::new();
//! nc.add_dimension("x", 10);
//! nc.add_variable("temperature", vec!["x"], "K", "Temperature",
//!     NetCdfData::Float64((0..10).map(|i| i as f64 * 1.5).collect()));
//! // nc.write(path).unwrap();
//! ```

#![cfg(feature = "netcdf")]

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use crate::error::{invalid_param, Error};
use crate::vector3::Vector3;

// ---------------------------------------------------------------------------
// NetCDF3 binary constants (XDR tags)
// ---------------------------------------------------------------------------
const NC_DIMENSION: u32 = 0x0000_000A;
const NC_ATTRIBUTE: u32 = 0x0000_000C;
const NC_VARIABLE: u32 = 0x0000_000E;

// NC type codes
const NC_CHAR: u32 = 2;
const NC_INT: u32 = 4;
const NC_FLOAT: u32 = 5;
const NC_DOUBLE: u32 = 6;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A named dimension with a fixed size.
pub struct NetCdfDimension {
    /// Dimension name (e.g., "x", "time", "points")
    pub name: String,
    /// Number of elements along this dimension
    pub size: usize,
}

/// Variable data storage variants.
pub enum NetCdfData {
    /// 64-bit IEEE floats (NC_DOUBLE)
    Float64(Vec<f64>),
    /// 32-bit IEEE floats (NC_FLOAT)
    Float32(Vec<f32>),
    /// 32-bit signed integers (NC_INT)
    Int32(Vec<i32>),
}

impl NetCdfData {
    fn len(&self) -> usize {
        match self {
            NetCdfData::Float64(v) => v.len(),
            NetCdfData::Float32(v) => v.len(),
            NetCdfData::Int32(v) => v.len(),
        }
    }

    fn type_code(&self) -> u32 {
        match self {
            NetCdfData::Float64(_) => NC_DOUBLE,
            NetCdfData::Float32(_) => NC_FLOAT,
            NetCdfData::Int32(_) => NC_INT,
        }
    }

    /// Byte size of one element for this type
    fn elem_bytes(&self) -> usize {
        match self {
            NetCdfData::Float64(_) => 8,
            NetCdfData::Float32(_) => 4,
            NetCdfData::Int32(_) => 4,
        }
    }
}

/// A named variable with dimension references, CF attributes, and data.
pub struct NetCdfVariable {
    /// Variable name
    pub name: String,
    /// Ordered dimension names
    pub dim_names: Vec<String>,
    /// CF units string (e.g., "m/s", "K", "1")
    pub units: String,
    /// Longer descriptive name
    pub long_name: String,
    /// The actual data values
    pub data: NetCdfData,
}

/// NetCDF3 Classic writer.
///
/// Accumulate dimensions, global attributes, and variables then call
/// `write` to produce a standards-compliant `.nc` file.
pub struct NetCdfWriter {
    /// Ordered list of dimensions
    pub dimensions: Vec<NetCdfDimension>,
    /// Global attributes (written as NC_CHAR strings)
    pub global_attributes: HashMap<String, String>,
    /// Ordered list of variables
    pub variables: Vec<NetCdfVariable>,
}

impl NetCdfWriter {
    /// Create an empty writer.
    pub fn new() -> Self {
        Self {
            dimensions: Vec::new(),
            global_attributes: HashMap::new(),
            variables: Vec::new(),
        }
    }

    /// Add a named fixed-size dimension.
    pub fn add_dimension(&mut self, name: &str, size: usize) {
        self.dimensions.push(NetCdfDimension {
            name: name.to_string(),
            size,
        });
    }

    /// Add a global string attribute.
    pub fn add_global_attribute(&mut self, name: &str, value: &str) {
        self.global_attributes
            .insert(name.to_string(), value.to_string());
    }

    /// Add a variable.
    pub fn add_variable(
        &mut self,
        name: &str,
        dim_names: Vec<&str>,
        units: &str,
        long_name: &str,
        data: NetCdfData,
    ) {
        self.variables.push(NetCdfVariable {
            name: name.to_string(),
            dim_names: dim_names.iter().map(|s| s.to_string()).collect(),
            units: units.to_string(),
            long_name: long_name.to_string(),
            data,
        });
    }

    /// Add three variables `{name}_x`, `{name}_y`, `{name}_z` for a vector field
    /// using CF-convention `units = "1"` and descriptive long names.
    pub fn write_vector_field_cf(
        &mut self,
        name: &str,
        data: &[Vector3<f64>],
        dim_names: Vec<&str>,
    ) {
        let xs: Vec<f64> = data.iter().map(|v| v.x).collect();
        let ys: Vec<f64> = data.iter().map(|v| v.y).collect();
        let zs: Vec<f64> = data.iter().map(|v| v.z).collect();

        self.add_variable(
            &format!("{}_x", name),
            dim_names.clone(),
            "1",
            &format!("{} x-component", name),
            NetCdfData::Float64(xs),
        );
        self.add_variable(
            &format!("{}_y", name),
            dim_names.clone(),
            "1",
            &format!("{} y-component", name),
            NetCdfData::Float64(ys),
        );
        self.add_variable(
            &format!("{}_z", name),
            dim_names,
            "1",
            &format!("{} z-component", name),
            NetCdfData::Float64(zs),
        );
    }

    /// Write the NetCDF3 Classic binary file to `path`.
    pub fn write(&self, path: &Path) -> Result<(), Error> {
        // Build dim-name → index map
        let dim_map: HashMap<&str, usize> = self
            .dimensions
            .iter()
            .enumerate()
            .map(|(i, d)| (d.name.as_str(), i))
            .collect();

        // Validate variables reference known dimensions
        for var in &self.variables {
            for dn in &var.dim_names {
                if !dim_map.contains_key(dn.as_str()) {
                    return Err(invalid_param(
                        "dimension",
                        &format!(
                            "variable '{}' references unknown dimension '{}'",
                            var.name, dn
                        ),
                    ));
                }
            }
        }

        // --- Layout: compute variable sizes and file begin offsets ---
        // Header size computation (approximate first; we'll fix begin offsets in a second pass)
        let header_size =
            compute_header_size(&self.dimensions, &self.global_attributes, &self.variables);

        // begin offset for first variable = header_size; subsequent ones follow
        let mut begins: Vec<u64> = Vec::with_capacity(self.variables.len());
        let mut offset = header_size as u64;
        for var in &self.variables {
            begins.push(offset);
            let vsize = padded_data_size(var.data.len(), var.data.elem_bytes());
            offset += vsize as u64;
        }

        // --- Write file ---
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        // Magic + numrecs
        w.write_all(b"CDF\x01")?;
        write_i32(&mut w, 0)?; // numrecs = 0 (no record variables)

        // Dim list
        write_dim_list(&mut w, &self.dimensions)?;

        // Global att list
        write_att_list(&mut w, &self.global_attributes)?;

        // Var list
        write_var_list(&mut w, &self.variables, &dim_map, &begins)?;

        // Data section
        for var in &self.variables {
            write_variable_data(&mut w, &var.data)?;
        }

        w.flush()?;
        Ok(())
    }
}

impl Default for NetCdfWriter {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// NetCDF3 reader (for round-trip tests)
// ---------------------------------------------------------------------------

/// Minimal NetCDF3 Classic reader.
///
/// Only reads the magic bytes and variable data; enough for test verification.
pub struct NetCdfReader;

impl NetCdfReader {
    /// Read a Float64 variable by name from a NetCDF3 file.
    ///
    /// Returns the flat data array for the named variable.
    pub fn read_variable_f64(path: &Path, var_name: &str) -> Result<Vec<f64>, Error> {
        let mut file = File::open(path)?;

        // Verify magic
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if &magic != b"CDF\x01" {
            return Err(invalid_param("magic", "not a NetCDF3 Classic file"));
        }

        // numrecs
        let _numrecs = read_be_i32(&mut file)?;

        // Parse dim_list
        let dims = parse_dim_list(&mut file)?;

        // Parse att_list (global)
        let _global_atts = parse_att_list(&mut file)?;

        // Parse var_list; find our variable
        let vars = parse_var_list(&mut file)?;

        let var = vars.iter().find(|v| v.name == var_name).ok_or_else(|| {
            invalid_param("var_name", &format!("variable '{}' not found", var_name))
        })?;

        if var.type_code != NC_DOUBLE {
            return Err(invalid_param("type", "variable is not Float64 (NC_DOUBLE)"));
        }

        // Calculate number of elements from dimensions
        let n_elems: usize = if var.dimids.is_empty() {
            1
        } else {
            var.dimids
                .iter()
                .map(|&id| dims.get(id).map(|d| d.size).unwrap_or(1))
                .product()
        };

        // Seek to begin
        file.seek(SeekFrom::Start(var.begin))?;

        let mut values = Vec::with_capacity(n_elems);
        for _ in 0..n_elems {
            let mut buf = [0u8; 8];
            file.read_exact(&mut buf)?;
            // XDR / big-endian IEEE 754 double
            values.push(f64::from_be_bytes(buf));
        }

        Ok(values)
    }

    /// Read and verify the magic bytes of a NetCDF3 file.
    pub fn verify_magic(path: &Path) -> Result<(), Error> {
        let mut file = File::open(path)?;
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if &magic != b"CDF\x01" {
            return Err(invalid_param("magic", "not a NetCDF3 Classic file"));
        }
        Ok(())
    }

    /// Return the names of all variables in a NetCDF3 file.
    pub fn list_variables(path: &Path) -> Result<Vec<String>, Error> {
        let mut file = File::open(path)?;
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if &magic != b"CDF\x01" {
            return Err(invalid_param("magic", "not a NetCDF3 Classic file"));
        }
        let _numrecs = read_be_i32(&mut file)?;
        let _dims = parse_dim_list(&mut file)?;
        let _atts = parse_att_list(&mut file)?;
        let vars = parse_var_list(&mut file)?;
        Ok(vars.into_iter().map(|v| v.name).collect())
    }
}

// ---------------------------------------------------------------------------
// Internal write helpers
// ---------------------------------------------------------------------------

fn write_i32(w: &mut impl Write, v: i32) -> Result<(), Error> {
    w.write_all(&v.to_be_bytes())?;
    Ok(())
}

fn write_u32(w: &mut impl Write, v: u32) -> Result<(), Error> {
    w.write_all(&v.to_be_bytes())?;
    Ok(())
}

/// Write a NetCDF3 padded string: i32 length, then bytes, then zero-padding to 4-byte boundary.
fn write_nc_string(w: &mut impl Write, s: &str) -> Result<(), Error> {
    let bytes = s.as_bytes();
    write_i32(w, bytes.len() as i32)?;
    w.write_all(bytes)?;
    let pad = (4 - (bytes.len() % 4)) % 4;
    w.write_all(&[0u8; 4][..pad])?;
    Ok(())
}

fn write_dim_list(w: &mut impl Write, dims: &[NetCdfDimension]) -> Result<(), Error> {
    if dims.is_empty() {
        write_u32(w, 0)?; // absent tag
        write_i32(w, 0)?;
        return Ok(());
    }
    write_u32(w, NC_DIMENSION)?;
    write_i32(w, dims.len() as i32)?;
    for dim in dims {
        write_nc_string(w, &dim.name)?;
        write_i32(w, dim.size as i32)?;
    }
    Ok(())
}

fn write_att_list(w: &mut impl Write, atts: &HashMap<String, String>) -> Result<(), Error> {
    if atts.is_empty() {
        write_u32(w, 0)?; // absent tag
        write_i32(w, 0)?;
        return Ok(());
    }
    write_u32(w, NC_ATTRIBUTE)?;
    write_i32(w, atts.len() as i32)?;
    // Sort for determinism
    let mut pairs: Vec<(&String, &String)> = atts.iter().collect();
    pairs.sort_by_key(|(k, _)| k.as_str());
    for (name, value) in pairs {
        write_nc_string(w, name)?;
        write_u32(w, NC_CHAR)?;
        let bytes = value.as_bytes();
        write_i32(w, bytes.len() as i32)?;
        w.write_all(bytes)?;
        let pad = (4 - (bytes.len() % 4)) % 4;
        w.write_all(&[0u8; 4][..pad])?;
    }
    Ok(())
}

fn write_var_list(
    w: &mut impl Write,
    vars: &[NetCdfVariable],
    dim_map: &HashMap<&str, usize>,
    begins: &[u64],
) -> Result<(), Error> {
    if vars.is_empty() {
        write_u32(w, 0)?;
        write_i32(w, 0)?;
        return Ok(());
    }
    write_u32(w, NC_VARIABLE)?;
    write_i32(w, vars.len() as i32)?;
    for (idx, var) in vars.iter().enumerate() {
        write_nc_string(w, &var.name)?;
        // dim ids
        write_i32(w, var.dim_names.len() as i32)?;
        for dn in &var.dim_names {
            let did = *dim_map.get(dn.as_str()).unwrap_or(&0);
            write_i32(w, did as i32)?;
        }
        // per-variable att list (units, long_name)
        let mut var_atts: HashMap<String, String> = HashMap::new();
        if !var.units.is_empty() {
            var_atts.insert("units".to_string(), var.units.clone());
        }
        if !var.long_name.is_empty() {
            var_atts.insert("long_name".to_string(), var.long_name.clone());
        }
        write_att_list(w, &var_atts)?;
        // nc_type
        write_u32(w, var.data.type_code())?;
        // vsize: padded data size in bytes
        let vsize = padded_data_size(var.data.len(), var.data.elem_bytes());
        write_i32(w, vsize as i32)?;
        // begin: file offset (stored as i64 in CDF-1 as two i32 words; we use i32 for simple files)
        // NetCDF3 classic uses i32 for begin (offset < 2^31)
        write_i32(w, begins[idx] as i32)?;
    }
    Ok(())
}

fn write_variable_data(w: &mut impl Write, data: &NetCdfData) -> Result<(), Error> {
    match data {
        NetCdfData::Float64(vals) => {
            for &v in vals {
                w.write_all(&v.to_be_bytes())?;
            }
            // Padding (each f64 is 8 bytes; pad to 4-byte boundary — always aligned)
            let raw_bytes = vals.len() * 8;
            let pad = (4 - (raw_bytes % 4)) % 4;
            w.write_all(&[0u8; 4][..pad])?;
        },
        NetCdfData::Float32(vals) => {
            for &v in vals {
                w.write_all(&v.to_bits().to_be_bytes())?;
            }
            let raw_bytes = vals.len() * 4;
            let pad = (4 - (raw_bytes % 4)) % 4;
            w.write_all(&[0u8; 4][..pad])?;
        },
        NetCdfData::Int32(vals) => {
            for &v in vals {
                w.write_all(&v.to_be_bytes())?;
            }
            let raw_bytes = vals.len() * 4;
            let pad = (4 - (raw_bytes % 4)) % 4;
            w.write_all(&[0u8; 4][..pad])?;
        },
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Size / layout computations
// ---------------------------------------------------------------------------

fn nc_string_size(s: &str) -> usize {
    4 + s.len() + (4 - (s.len() % 4)) % 4
}

fn att_entry_size(name: &str, value: &str) -> usize {
    let bytes = value.len();
    nc_string_size(name) + 4 + 4 + bytes + (4 - (bytes % 4)) % 4
}

fn var_att_list_size(units: &str, long_name: &str) -> usize {
    let mut atts_non_empty = 0;
    let mut size = 8; // tag + count
    if !units.is_empty() {
        size += att_entry_size("units", units);
        atts_non_empty += 1;
    }
    if !long_name.is_empty() {
        size += att_entry_size("long_name", long_name);
        atts_non_empty += 1;
    }
    if atts_non_empty == 0 {
        size = 8;
    }
    size
}

fn var_header_size(var: &NetCdfVariable) -> usize {
    nc_string_size(&var.name)
        + 4 // ndims
        + var.dim_names.len() * 4
        + var_att_list_size(&var.units, &var.long_name)
        + 4 // nc_type
        + 4 // vsize
        + 4 // begin (i32)
}

fn global_att_list_size(atts: &HashMap<String, String>) -> usize {
    if atts.is_empty() {
        return 8;
    }
    let mut size = 8; // tag + count
    for (k, v) in atts {
        size += att_entry_size(k, v);
    }
    size
}

fn padded_data_size(n_elems: usize, elem_bytes: usize) -> usize {
    let raw = n_elems * elem_bytes;
    raw + (4 - (raw % 4)) % 4
}

fn compute_header_size(
    dims: &[NetCdfDimension],
    global_atts: &HashMap<String, String>,
    vars: &[NetCdfVariable],
) -> usize {
    // magic(4) + numrecs(4)
    let mut size = 8;

    // dim_list
    if dims.is_empty() {
        size += 8;
    } else {
        size += 8; // tag + count
        for d in dims {
            size += nc_string_size(&d.name) + 4; // + size i32
        }
    }

    // global att_list
    size += global_att_list_size(global_atts);

    // var_list
    if vars.is_empty() {
        size += 8;
    } else {
        size += 8; // tag + count
        for v in vars {
            size += var_header_size(v);
        }
    }

    size
}

// ---------------------------------------------------------------------------
// Internal read helpers (for NetCdfReader)
// ---------------------------------------------------------------------------

fn read_be_i32(r: &mut impl Read) -> Result<i32, Error> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(i32::from_be_bytes(buf))
}

fn read_be_u32(r: &mut impl Read) -> Result<u32, Error> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(u32::from_be_bytes(buf))
}

fn read_nc_string(r: &mut impl Read) -> Result<String, Error> {
    let len = read_be_i32(r)? as usize;
    let mut bytes = vec![0u8; len];
    r.read_exact(&mut bytes)?;
    let pad = (4 - (len % 4)) % 4;
    let mut padding = vec![0u8; pad];
    r.read_exact(&mut padding)?;
    Ok(String::from_utf8_lossy(&bytes).into_owned())
}

struct ParsedDim {
    _name: String,
    size: usize,
}

fn parse_dim_list(r: &mut impl Read) -> Result<Vec<ParsedDim>, Error> {
    let tag = read_be_u32(r)?;
    let ndims = read_be_i32(r)? as usize;
    if tag == 0 {
        return Ok(Vec::new());
    }
    if tag != NC_DIMENSION {
        return Err(invalid_param("tag", "expected NC_DIMENSION in dim_list"));
    }
    let mut dims = Vec::with_capacity(ndims);
    for _ in 0..ndims {
        let name = read_nc_string(r)?;
        let size = read_be_i32(r)? as usize;
        dims.push(ParsedDim { _name: name, size });
    }
    Ok(dims)
}

struct ParsedAtt {
    _name: String,
}

fn parse_att_list(r: &mut impl Read) -> Result<Vec<ParsedAtt>, Error> {
    let tag = read_be_u32(r)?;
    let natts = read_be_i32(r)? as usize;
    if tag == 0 {
        return Ok(Vec::new());
    }
    if tag != NC_ATTRIBUTE {
        return Err(invalid_param("tag", "expected NC_ATTRIBUTE in att_list"));
    }
    let mut atts = Vec::with_capacity(natts);
    for _ in 0..natts {
        let name = read_nc_string(r)?;
        let atype = read_be_u32(r)?;
        let len = read_be_i32(r)? as usize;
        let elem_bytes = match atype {
            NC_CHAR => 1,
            NC_INT => 4,
            NC_FLOAT => 4,
            NC_DOUBLE => 8,
            _ => 1,
        };
        let raw_bytes = len * elem_bytes;
        let mut buf = vec![0u8; raw_bytes];
        r.read_exact(&mut buf)?;
        let pad = (4 - (raw_bytes % 4)) % 4;
        let mut pbuf = vec![0u8; pad];
        r.read_exact(&mut pbuf)?;
        atts.push(ParsedAtt { _name: name });
    }
    Ok(atts)
}

struct ParsedVar {
    name: String,
    dimids: Vec<usize>,
    type_code: u32,
    begin: u64,
}

fn parse_var_list(r: &mut impl Read) -> Result<Vec<ParsedVar>, Error> {
    let tag = read_be_u32(r)?;
    let nvars = read_be_i32(r)? as usize;
    if tag == 0 {
        return Ok(Vec::new());
    }
    if tag != NC_VARIABLE {
        return Err(invalid_param("tag", "expected NC_VARIABLE in var_list"));
    }
    let mut vars = Vec::with_capacity(nvars);
    for _ in 0..nvars {
        let name = read_nc_string(r)?;
        let ndimids = read_be_i32(r)? as usize;
        let mut dimids = Vec::with_capacity(ndimids);
        for _ in 0..ndimids {
            dimids.push(read_be_i32(r)? as usize);
        }
        // skip var att list
        let _vatts = parse_att_list(r)?;
        let nc_type = read_be_u32(r)?;
        let _vsize = read_be_i32(r)?;
        let begin = read_be_i32(r)? as u64;
        vars.push(ParsedVar {
            name,
            dimids,
            type_code: nc_type,
            begin,
        });
    }
    Ok(vars)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn temp_path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(name)
    }

    #[test]
    fn test_verify_magic() {
        let mut nc = NetCdfWriter::new();
        nc.add_dimension("x", 4);
        nc.add_variable(
            "v",
            vec!["x"],
            "1",
            "test",
            NetCdfData::Float64(vec![1.0; 4]),
        );
        let path = temp_path("nc_magic.nc");
        nc.write(&path).expect("write should succeed");
        NetCdfReader::verify_magic(&path).expect("magic should be valid");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_roundtrip_f64_variable() {
        let values: Vec<f64> = (0..10).map(|i| i as f64 * 1.5).collect();
        let mut nc = NetCdfWriter::new();
        nc.add_dimension("points", 10);
        nc.add_variable(
            "temperature",
            vec!["points"],
            "K",
            "Air temperature",
            NetCdfData::Float64(values.clone()),
        );
        let path = temp_path("nc_roundtrip_f64.nc");
        nc.write(&path).expect("write should succeed");

        let read_back =
            NetCdfReader::read_variable_f64(&path, "temperature").expect("read should succeed");
        assert_eq!(read_back.len(), values.len());
        for (a, b) in read_back.iter().zip(values.iter()) {
            assert!((a - b).abs() < 1e-12, "value mismatch: {} vs {}", a, b);
        }
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_vector_field_cf_creates_xyz_vars() {
        let data: Vec<Vector3<f64>> = vec![
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(4.0, 5.0, 6.0),
            Vector3::new(7.0, 8.0, 9.0),
        ];
        let mut nc = NetCdfWriter::new();
        nc.add_dimension("n", 3);
        nc.write_vector_field_cf("spin", &data, vec!["n"]);

        let path = temp_path("nc_vector_cf.nc");
        nc.write(&path).expect("write should succeed");

        let vars = NetCdfReader::list_variables(&path).expect("list should succeed");
        assert!(vars.contains(&"spin_x".to_string()));
        assert!(vars.contains(&"spin_y".to_string()));
        assert!(vars.contains(&"spin_z".to_string()));

        // Verify x component values
        let xs = NetCdfReader::read_variable_f64(&path, "spin_x").expect("read x");
        assert!((xs[0] - 1.0).abs() < 1e-12);
        assert!((xs[1] - 4.0).abs() < 1e-12);
        assert!((xs[2] - 7.0).abs() < 1e-12);

        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_unknown_dimension_error() {
        let mut nc = NetCdfWriter::new();
        nc.add_variable(
            "v",
            vec!["nonexistent_dim"],
            "1",
            "test",
            NetCdfData::Float64(vec![]),
        );
        let path = temp_path("nc_bad_dim.nc");
        let result = nc.write(&path);
        assert!(result.is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_global_attributes_written() {
        let mut nc = NetCdfWriter::new();
        nc.add_global_attribute("Conventions", "CF-1.8");
        nc.add_global_attribute("institution", "COOLJAPAN");
        nc.add_dimension("x", 2);
        nc.add_variable(
            "v",
            vec!["x"],
            "1",
            "test",
            NetCdfData::Float64(vec![1.0, 2.0]),
        );
        let path = temp_path("nc_global_atts.nc");
        nc.write(&path).expect("write should succeed");
        // File must at least be created and have valid magic
        NetCdfReader::verify_magic(&path).expect("magic ok");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_roundtrip_multiple_variables() {
        let mut nc = NetCdfWriter::new();
        nc.add_dimension("time", 5);
        let temp_vals: Vec<f64> = (0..5).map(|i| 300.0 + i as f64).collect();
        let pressure_vals: Vec<f64> = (0..5).map(|i| 1013.0 - i as f64 * 0.1).collect();
        nc.add_variable(
            "temperature",
            vec!["time"],
            "K",
            "Temperature",
            NetCdfData::Float64(temp_vals.clone()),
        );
        nc.add_variable(
            "pressure",
            vec!["time"],
            "hPa",
            "Pressure",
            NetCdfData::Float64(pressure_vals.clone()),
        );
        let path = temp_path("nc_multi_var.nc");
        nc.write(&path).expect("write should succeed");

        let t = NetCdfReader::read_variable_f64(&path, "temperature").unwrap();
        let p = NetCdfReader::read_variable_f64(&path, "pressure").unwrap();
        for (a, b) in t.iter().zip(temp_vals.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
        for (a, b) in p.iter().zip(pressure_vals.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_list_variables() {
        let mut nc = NetCdfWriter::new();
        nc.add_dimension("n", 3);
        nc.add_variable(
            "alpha",
            vec!["n"],
            "1",
            "",
            NetCdfData::Float64(vec![1.0, 2.0, 3.0]),
        );
        nc.add_variable(
            "beta",
            vec!["n"],
            "1",
            "",
            NetCdfData::Float64(vec![4.0, 5.0, 6.0]),
        );
        let path = temp_path("nc_list_vars.nc");
        nc.write(&path).expect("write should succeed");

        let vars = NetCdfReader::list_variables(&path).unwrap();
        assert!(vars.contains(&"alpha".to_string()));
        assert!(vars.contains(&"beta".to_string()));
        assert_eq!(vars.len(), 2);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_default_impl() {
        let nc = NetCdfWriter::default();
        assert!(nc.dimensions.is_empty());
        assert!(nc.variables.is_empty());
    }
}
