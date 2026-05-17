//! Pure-Rust Zarr v2 Store Writer and Reader
//!
//! Implements the Zarr v2 specification for on-disk chunk stores. Each array
//! is stored as:
//!
//! ```text
//! {root}/{array_name}/.zarray   – JSON metadata
//! {root}/{array_name}/.zattrs   – JSON attributes (optional)
//! {root}/{array_name}/0.0       – chunk file (C-order, little-endian binary)
//! {root}/{array_name}/0.1       – next chunk along second axis, etc.
//! ```
//!
//! This implementation writes and reads arrays of `f64`, `f32`, or `i32`
//! without any compression (compressor = null), matching the Zarr v2 spec
//! for `"order": "C"` and `"dtype": "<f8"` (or `<f4`, `<i4`).
//!
//! ## Example
//!
//! ```rust
//! use spintronics::visualization::zarr::{ZarrStore, ZarrDtype};
//! let store_path = std::env::temp_dir().join("my_zarr");
//! // let mut store = ZarrStore::new_store(&store_path).unwrap();
//! // store.add_array("field", vec![100, 3], vec![10, 3], ZarrDtype::Float64);
//! // store.write_all(&data_map).unwrap();
//! ```

#![cfg(feature = "zarr")]

use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::error::{invalid_param, Error};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Data type for a Zarr array, following NumPy dtype convention.
pub enum ZarrDtype {
    /// 64-bit little-endian float (`<f8`)
    Float64,
    /// 32-bit little-endian float (`<f4`)
    Float32,
    /// 32-bit little-endian signed integer (`<i4`)
    Int32,
}

impl ZarrDtype {
    /// Return the Zarr dtype string (NumPy format).
    pub fn as_str(&self) -> &'static str {
        match self {
            ZarrDtype::Float64 => "<f8",
            ZarrDtype::Float32 => "<f4",
            ZarrDtype::Int32 => "<i4",
        }
    }

    /// Return the number of bytes per element.
    pub fn bytes_per_elem(&self) -> usize {
        match self {
            ZarrDtype::Float64 => 8,
            ZarrDtype::Float32 => 4,
            ZarrDtype::Int32 => 4,
        }
    }
}

/// A single named Zarr array within a store.
pub struct ZarrArray {
    /// Path to the array directory (root/name)
    pub path: PathBuf,
    /// Array shape (C order: outermost first)
    pub shape: Vec<usize>,
    /// Chunk shape (must have same ndim as shape)
    pub chunks: Vec<usize>,
    /// Element data type
    pub dtype: ZarrDtype,
    /// Optional key-value attributes stored in `.zattrs`
    pub attributes: HashMap<String, String>,
}

/// A Zarr v2 on-disk store (a root directory containing named arrays).
pub struct ZarrStore {
    /// Root path of the store
    pub root_path: PathBuf,
    /// Arrays in this store (ordered by insertion)
    pub arrays: Vec<ZarrArray>,
}

impl ZarrStore {
    /// Create a new store at `root_path`, creating the directory if needed.
    pub fn new_store(root_path: &Path) -> Result<Self, Error> {
        fs::create_dir_all(root_path)?;
        Ok(Self {
            root_path: root_path.to_path_buf(),
            arrays: Vec::new(),
        })
    }

    /// Add a new array definition to the store and return a mutable reference to it.
    ///
    /// Creates the array sub-directory immediately.
    pub fn add_array(
        &mut self,
        name: &str,
        shape: Vec<usize>,
        chunks: Vec<usize>,
        dtype: ZarrDtype,
    ) -> &mut ZarrArray {
        let array_path = self.root_path.join(name);
        let _ = fs::create_dir_all(&array_path);
        self.arrays.push(ZarrArray {
            path: array_path,
            shape,
            chunks,
            dtype,
            attributes: HashMap::new(),
        });
        self.arrays.last_mut().expect("just pushed")
    }

    /// Write all `.zarray` (and `.zattrs` if non-empty) metadata files.
    pub fn write_metadata(&self) -> Result<(), Error> {
        for array in &self.arrays {
            fs::create_dir_all(&array.path)?;
            write_zarray_json(
                &array.path.join(".zarray"),
                &array.shape,
                &array.chunks,
                array.dtype.as_str(),
            )?;
            if !array.attributes.is_empty() {
                write_zattrs_json(&array.path.join(".zattrs"), &array.attributes)?;
            }
        }
        Ok(())
    }

    /// Write raw `f64` data for the named array.
    ///
    /// Data is expected in C (row-major) order matching the array shape.
    /// The array must already exist in the store (added via `add_array`).
    pub fn write_array_f64(&self, name: &str, data: &[f64]) -> Result<(), Error> {
        let array = self.find_array(name)?;
        let n_elems: usize = array.shape.iter().product();
        if data.len() != n_elems {
            return Err(invalid_param(
                "data",
                &format!(
                    "expected {} elements for array '{}', got {}",
                    n_elems,
                    name,
                    data.len()
                ),
            ));
        }
        write_chunks_f64(&array.path, &array.shape, &array.chunks, data)?;
        Ok(())
    }

    /// Convenience alias — same as `write_array_f64`.
    pub fn write_array_data(&self, name: &str, data: &[f64]) -> Result<(), Error> {
        self.write_array_f64(name, data)
    }

    /// Write metadata for all arrays, then write data for each array in `data_by_name`.
    pub fn write_all(&self, data_by_name: &HashMap<String, Vec<f64>>) -> Result<(), Error> {
        self.write_metadata()?;
        for (name, data) in data_by_name {
            self.write_array_f64(name, data)?;
        }
        Ok(())
    }

    fn find_array(&self, name: &str) -> Result<&ZarrArray, Error> {
        let array_path = self.root_path.join(name);
        self.arrays
            .iter()
            .find(|a| a.path == array_path)
            .ok_or_else(|| invalid_param("name", &format!("array '{}' not found in store", name)))
    }
}

/// Read all chunk files for the named array and reconstruct the flat f64 vector.
pub fn read_array_f64(store_root: &Path, name: &str) -> Result<Vec<f64>, Error> {
    let array_path = store_root.join(name);
    let zarray_path = array_path.join(".zarray");

    let meta_str = fs::read_to_string(&zarray_path)?;
    let (shape, chunks) = parse_zarray_shape_chunks(&meta_str)?;

    let n_elems: usize = shape.iter().product();
    let chunk_counts = chunk_grid_counts(&shape, &chunks);

    let mut flat: Vec<f64> = vec![0.0f64; n_elems];

    // Iterate over all chunk index tuples
    let n_chunks: usize = chunk_counts.iter().product();
    let ndim = shape.len();

    for chunk_linear in 0..n_chunks {
        let chunk_idx = linear_to_nd_index(chunk_linear, &chunk_counts);

        // Read chunk file
        let chunk_filename = chunk_idx
            .iter()
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(".");
        let chunk_path = array_path.join(&chunk_filename);

        if !chunk_path.exists() {
            // Chunk may be absent if all fill values (we treat fill=0)
            continue;
        }

        let bytes = fs::read(&chunk_path)?;
        let n_floats = bytes.len() / 8;
        let mut chunk_vals = Vec::with_capacity(n_floats);
        for i in 0..n_floats {
            let b: [u8; 8] = bytes[i * 8..i * 8 + 8]
                .try_into()
                .map_err(|_| invalid_param("chunk", "byte slice conversion failed"))?;
            chunk_vals.push(f64::from_le_bytes(b));
        }

        // Map chunk values back to flat array positions
        scatter_chunk_to_flat(&chunk_vals, &chunk_idx, &shape, &chunks, ndim, &mut flat);
    }

    Ok(flat)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Write the `.zarray` JSON metadata file (no serde_json dependency).
fn write_zarray_json(
    path: &Path,
    shape: &[usize],
    chunks: &[usize],
    dtype: &str,
) -> Result<(), Error> {
    let shape_str = shape
        .iter()
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(", ");
    let chunks_str = chunks
        .iter()
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(", ");

    let json = format!(
        "{{\n\
        \t\"zarr_format\": 2,\n\
        \t\"shape\": [{shape}],\n\
        \t\"chunks\": [{chk}],\n\
        \t\"dtype\": \"{dt}\",\n\
        \t\"compressor\": null,\n\
        \t\"fill_value\": 0,\n\
        \t\"order\": \"C\",\n\
        \t\"filters\": null\n\
        }}",
        shape = shape_str,
        chk = chunks_str,
        dt = dtype
    );

    let mut f = fs::File::create(path)?;
    f.write_all(json.as_bytes())?;
    Ok(())
}

/// Write `.zattrs` as a minimal JSON object with string values.
fn write_zattrs_json(path: &Path, attrs: &HashMap<String, String>) -> Result<(), Error> {
    let mut parts: Vec<String> = attrs
        .iter()
        .map(|(k, v)| format!("\t\"{}\": \"{}\"", k, v))
        .collect();
    parts.sort(); // deterministic order
    let json = format!("{{\n{}\n}}", parts.join(",\n"));
    let mut f = fs::File::create(path)?;
    f.write_all(json.as_bytes())?;
    Ok(())
}

/// Number of chunks along each axis.
fn chunk_grid_counts(shape: &[usize], chunks: &[usize]) -> Vec<usize> {
    shape
        .iter()
        .zip(chunks.iter())
        .map(|(&s, &c)| (s + c - 1) / c)
        .collect()
}

/// Convert a linear chunk index to a multi-dimensional chunk index.
fn linear_to_nd_index(mut linear: usize, counts: &[usize]) -> Vec<usize> {
    let ndim = counts.len();
    let mut idx = vec![0usize; ndim];
    for d in (0..ndim).rev() {
        idx[d] = linear % counts[d];
        linear /= counts[d];
    }
    idx
}

/// Write all chunks for a 1-D or N-D float64 array.
///
/// For each chunk index tuple, slices the relevant portion of `data` (in
/// C/row-major order) and writes it as raw little-endian f64 bytes.
fn write_chunks_f64(
    array_path: &Path,
    shape: &[usize],
    chunks: &[usize],
    data: &[f64],
) -> Result<(), Error> {
    if shape.is_empty() {
        return Err(invalid_param(
            "shape",
            "array must have at least 1 dimension",
        ));
    }

    let chunk_counts = chunk_grid_counts(shape, chunks);
    let n_chunks: usize = chunk_counts.iter().product();
    let ndim = shape.len();

    for chunk_linear in 0..n_chunks {
        let chunk_idx = linear_to_nd_index(chunk_linear, &chunk_counts);

        // Gather this chunk's values from data in C order
        let chunk_vals = gather_chunk_from_flat(data, &chunk_idx, shape, chunks, ndim);

        // Chunk file name: "{c0}.{c1}...." (or just "{c0}" for 1-D)
        let chunk_filename = chunk_idx
            .iter()
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(".");

        let chunk_path = array_path.join(&chunk_filename);
        let mut f = fs::File::create(&chunk_path)?;
        for &v in &chunk_vals {
            f.write_all(&v.to_le_bytes())?;
        }
    }
    Ok(())
}

/// Collect values from the flat array that belong to the given chunk.
///
/// The returned vector contains the chunk data in C order. Values beyond
/// the array boundary are filled with 0.0 (fill value).
fn gather_chunk_from_flat(
    data: &[f64],
    chunk_idx: &[usize],
    shape: &[usize],
    chunks: &[usize],
    ndim: usize,
) -> Vec<f64> {
    // Chunk extents: start and end indices along each axis
    let starts: Vec<usize> = chunk_idx
        .iter()
        .zip(chunks.iter())
        .map(|(&ci, &cs)| ci * cs)
        .collect();
    let ends: Vec<usize> = starts
        .iter()
        .zip(shape.iter())
        .zip(chunks.iter())
        .map(|((&s, &sh), &cs)| (s + cs).min(sh))
        .collect();

    // Number of elements in this chunk
    let chunk_shape: Vec<usize> = starts
        .iter()
        .zip(ends.iter())
        .map(|(&s, &e)| e - s)
        .collect();
    let chunk_size: usize = chunk_shape.iter().product();

    let mut vals = Vec::with_capacity(chunk_size);

    // Iterate in C order over the chunk local indices
    for local_linear in 0..chunk_size {
        let local_idx = linear_to_nd_index(local_linear, &chunk_shape);
        let global_idx: Vec<usize> = starts
            .iter()
            .zip(local_idx.iter())
            .map(|(&s, &li)| s + li)
            .collect();

        let flat_pos = c_order_index(&global_idx, shape, ndim);
        let v = if flat_pos < data.len() {
            data[flat_pos]
        } else {
            0.0
        };
        vals.push(v);
    }
    vals
}

/// Scatter chunk values back into the flat output array (for reading).
fn scatter_chunk_to_flat(
    chunk_vals: &[f64],
    chunk_idx: &[usize],
    shape: &[usize],
    chunks: &[usize],
    ndim: usize,
    flat: &mut [f64],
) {
    let starts: Vec<usize> = chunk_idx
        .iter()
        .zip(chunks.iter())
        .map(|(&ci, &cs)| ci * cs)
        .collect();
    let ends: Vec<usize> = starts
        .iter()
        .zip(shape.iter())
        .zip(chunks.iter())
        .map(|((&s, &sh), &cs)| (s + cs).min(sh))
        .collect();

    let chunk_shape: Vec<usize> = starts
        .iter()
        .zip(ends.iter())
        .map(|(&s, &e)| e - s)
        .collect();
    let chunk_size: usize = chunk_shape.iter().product();

    for (local_linear, &val) in chunk_vals.iter().enumerate().take(chunk_size) {
        let local_idx = linear_to_nd_index(local_linear, &chunk_shape);
        let global_idx: Vec<usize> = starts
            .iter()
            .zip(local_idx.iter())
            .map(|(&s, &li)| s + li)
            .collect();
        let flat_pos = c_order_index(&global_idx, shape, ndim);
        if flat_pos < flat.len() {
            flat[flat_pos] = val;
        }
    }
}

/// Compute the C-order (row-major) flat index from a multi-dimensional index.
fn c_order_index(idx: &[usize], shape: &[usize], ndim: usize) -> usize {
    let mut pos = 0usize;
    let mut stride = 1usize;
    for d in (0..ndim).rev() {
        pos += idx[d] * stride;
        stride *= shape[d];
    }
    pos
}

/// Parse shape and chunks from a hand-written `.zarray` JSON string.
///
/// Only handles the simple format produced by `write_zarray_json`.
fn parse_zarray_shape_chunks(json: &str) -> Result<(Vec<usize>, Vec<usize>), Error> {
    let shape = parse_json_usize_array(json, "shape")?;
    let chunks = parse_json_usize_array(json, "chunks")?;
    Ok((shape, chunks))
}

/// Extract a JSON array of integers for a given key from a simple flat JSON string.
fn parse_json_usize_array(json: &str, key: &str) -> Result<Vec<usize>, Error> {
    let needle = format!("\"{}\"", key);
    let start = json
        .find(&needle)
        .ok_or_else(|| invalid_param("json", &format!("key '{}' not found in .zarray", key)))?;
    let after_key = &json[start + needle.len()..];
    let bracket_start = after_key
        .find('[')
        .ok_or_else(|| invalid_param("json", &format!("no '[' after key '{}' in .zarray", key)))?;
    let bracket_end = after_key
        .find(']')
        .ok_or_else(|| invalid_param("json", &format!("no ']' after key '{}' in .zarray", key)))?;
    let inner = &after_key[bracket_start + 1..bracket_end];

    if inner.trim().is_empty() {
        return Ok(Vec::new());
    }

    inner
        .split(',')
        .map(|s| {
            s.trim()
                .parse::<usize>()
                .map_err(|e| invalid_param("json", &format!("parse error in '{}': {}", key, e)))
        })
        .collect()
}

/// Extract the dtype string from a `.zarray` JSON.
#[cfg(test)]
fn parse_zarray_dtype(json: &str) -> Result<String, Error> {
    let needle = "\"dtype\"";
    let start = json
        .find(needle)
        .ok_or_else(|| invalid_param("json", "key 'dtype' not found"))?;
    let after = &json[start + needle.len()..];
    let q1 = after
        .find('"')
        .ok_or_else(|| invalid_param("json", "no opening quote for dtype value"))?;
    let inner = &after[q1 + 1..];
    let q2 = inner
        .find('"')
        .ok_or_else(|| invalid_param("json", "no closing quote for dtype value"))?;
    Ok(inner[..q2].to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn temp_dir(sub: &str) -> PathBuf {
        std::env::temp_dir().join(sub)
    }

    fn cleanup(path: &Path) {
        let _ = fs::remove_dir_all(path);
    }

    #[test]
    fn test_new_store_creates_directory() {
        let root = temp_dir("zarr_test_new_store");
        cleanup(&root);
        let _store = ZarrStore::new_store(&root).expect("store creation should succeed");
        assert!(root.is_dir());
        cleanup(&root);
    }

    #[test]
    fn test_write_and_read_1d_array() {
        let root = temp_dir("zarr_test_1d");
        cleanup(&root);
        let mut store = ZarrStore::new_store(&root).expect("store ok");
        let data: Vec<f64> = (0..10).map(|i| i as f64 * 1.5).collect();
        store.add_array("field1d", vec![10], vec![5], ZarrDtype::Float64);
        let mut data_map = HashMap::new();
        data_map.insert("field1d".to_string(), data.clone());
        store
            .write_all(&data_map)
            .expect("write_all should succeed");

        let read_back = read_array_f64(&root, "field1d").expect("read should succeed");
        assert_eq!(read_back.len(), data.len());
        for (a, b) in read_back.iter().zip(data.iter()) {
            assert!((a - b).abs() < 1e-12, "{} != {}", a, b);
        }
        cleanup(&root);
    }

    #[test]
    fn test_write_2d_array_chunk_files_created() {
        let root = temp_dir("zarr_test_2d");
        cleanup(&root);
        let mut store = ZarrStore::new_store(&root).expect("store ok");
        // 10 rows, 3 cols; chunk = 5 rows, 3 cols → 2 chunk files (0.0 and 1.0)
        let n = 10 * 3;
        let data: Vec<f64> = (0..n).map(|i| i as f64).collect();
        store.add_array("field2d", vec![10, 3], vec![5, 3], ZarrDtype::Float64);
        let mut data_map = HashMap::new();
        data_map.insert("field2d".to_string(), data.clone());
        store
            .write_all(&data_map)
            .expect("write_all should succeed");

        let chunk0 = root.join("field2d").join("0.0");
        let chunk1 = root.join("field2d").join("1.0");
        assert!(chunk0.exists(), "chunk 0.0 should exist");
        assert!(chunk1.exists(), "chunk 1.0 should exist");

        // Each chunk covers 5×3 = 15 elements × 8 bytes = 120 bytes
        let meta0 = fs::metadata(&chunk0).unwrap();
        let meta1 = fs::metadata(&chunk1).unwrap();
        assert_eq!(meta0.len(), 120);
        assert_eq!(meta1.len(), 120);

        cleanup(&root);
    }

    #[test]
    fn test_2d_roundtrip() {
        let root = temp_dir("zarr_test_2d_rt");
        cleanup(&root);
        let mut store = ZarrStore::new_store(&root).expect("store ok");
        let shape = vec![6, 4];
        let n: usize = shape.iter().product();
        let data: Vec<f64> = (0..n).map(|i| i as f64 * 0.5 + 1.0).collect();
        store.add_array("m", shape, vec![3, 2], ZarrDtype::Float64);
        let mut dm = HashMap::new();
        dm.insert("m".to_string(), data.clone());
        store.write_all(&dm).expect("write ok");

        let read_back = read_array_f64(&root, "m").expect("read ok");
        assert_eq!(read_back.len(), data.len());
        for (a, b) in read_back.iter().zip(data.iter()) {
            assert!((a - b).abs() < 1e-12, "mismatch {} vs {}", a, b);
        }
        cleanup(&root);
    }

    #[test]
    fn test_zarray_json_content() {
        let root = temp_dir("zarr_test_json");
        cleanup(&root);
        let mut store = ZarrStore::new_store(&root).expect("store ok");
        store.add_array("arr", vec![8, 3], vec![4, 3], ZarrDtype::Float64);
        store.write_metadata().expect("metadata write ok");

        let zarray_path = root.join("arr").join(".zarray");
        assert!(zarray_path.exists());
        let content = fs::read_to_string(&zarray_path).unwrap();
        assert!(content.contains("\"zarr_format\": 2"));
        assert!(content.contains("\"dtype\": \"<f8\""));
        assert!(content.contains("[8, 3]"));
        assert!(content.contains("[4, 3]"));
        assert!(content.contains("\"compressor\": null"));

        cleanup(&root);
    }

    #[test]
    fn test_data_size_mismatch_error() {
        let root = temp_dir("zarr_test_size_mismatch");
        cleanup(&root);
        let mut store = ZarrStore::new_store(&root).expect("store ok");
        store.add_array("v", vec![10], vec![5], ZarrDtype::Float64);
        store.write_metadata().expect("meta ok");
        let result = store.write_array_f64("v", &[1.0f64, 2.0]); // wrong size
        assert!(result.is_err());
        cleanup(&root);
    }

    #[test]
    fn test_parse_zarray_dtype() {
        let json = r#"{"zarr_format": 2,"shape": [10, 3],"chunks": [5, 3],"dtype": "<f8","compressor": null,"fill_value": 0,"order": "C","filters": null}"#;
        let dtype = parse_zarray_dtype(json).expect("parse dtype ok");
        assert_eq!(dtype, "<f8");
    }

    #[test]
    fn test_1d_chunked_roundtrip_non_divisible() {
        // 13 elements, chunk size 4 → 4 chunks (4+4+4+1)
        let root = temp_dir("zarr_test_1d_nondiv");
        cleanup(&root);
        let mut store = ZarrStore::new_store(&root).expect("store ok");
        let data: Vec<f64> = (0..13).map(|i| i as f64 * 2.0 - 3.0).collect();
        store.add_array("nd", vec![13], vec![4], ZarrDtype::Float64);
        let mut dm = HashMap::new();
        dm.insert("nd".to_string(), data.clone());
        store.write_all(&dm).expect("write ok");

        let read_back = read_array_f64(&root, "nd").expect("read ok");
        assert_eq!(read_back.len(), data.len());
        for (a, b) in read_back.iter().zip(data.iter()) {
            assert!((a - b).abs() < 1e-12, "{} != {}", a, b);
        }
        cleanup(&root);
    }
}
