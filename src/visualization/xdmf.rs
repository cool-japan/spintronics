//! XDMF Writer — XML descriptor with raw binary heavy data
//!
//! XDMF (eXtensible Data Model and Format) uses XML to describe
//! structured or unstructured grids and their associated data. Heavy
//! data (the actual floating-point arrays) is written to a companion
//! `.bin` file in raw little-endian f64 format. The XDMF XML points
//! to this binary file using seek offsets.
//!
//! This implementation supports:
//! - Uniform 3D co-rectilinear meshes (ORIGIN_DXDYDZ topology)
//! - Temporal collections (multiple time steps)
//! - Vector and scalar node-centred attributes
//! - `Format="Binary"` references with explicit seek offsets

#![cfg(feature = "xdmf")]

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::error::{invalid_param, Error};
use crate::vector3::Vector3;

/// A snapshot at one simulation time.
pub struct XdmfTimeStep {
    /// Physical time value
    pub time: f64,
    /// Named fields at this time step
    pub fields: HashMap<String, XdmfFieldData>,
}

/// Field data variants for XDMF output.
pub enum XdmfFieldData {
    /// 3-component vector field stored as f64
    Vector(Vec<Vector3<f64>>),
    /// Single-component scalar field stored as f64
    Scalar(Vec<f64>),
}

/// Type of the XDMF grid topology to emit.
pub enum XdmfGridType {
    /// Regular 3-D co-rectilinear mesh (ORIGIN_DXDYDZ)
    Uniform3D,
    /// Placeholder for future unstructured support
    Unstructured,
}

/// XDMF writer for uniform 3-D structured grids with temporal data.
///
/// Call `add_time_step` for each simulation frame, then `write` to
/// produce the `.xmf` descriptor and the companion `.bin` binary file.
pub struct XdmfWriter {
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
    /// Ordered sequence of time steps
    pub time_steps: Vec<XdmfTimeStep>,
}

impl XdmfWriter {
    /// Create a new writer for a uniform 3-D co-rectilinear mesh.
    pub fn new_uniform_grid(
        nx: usize,
        ny: usize,
        nz: usize,
        origin: [f64; 3],
        spacing: [f64; 3],
    ) -> Self {
        Self {
            nx,
            ny,
            nz,
            origin,
            spacing,
            time_steps: Vec::new(),
        }
    }

    /// Append a time step with named field data.
    pub fn add_time_step(&mut self, time: f64, fields: HashMap<String, XdmfFieldData>) {
        self.time_steps.push(XdmfTimeStep { time, fields });
    }

    /// Write the XDMF descriptor XML to `xmf_path` and binary payload to `bin_path`.
    ///
    /// Binary layout: fields are written in the order they appear (sorted by
    /// name within each time step, time steps in order). Each field occupies
    /// `n_points * n_components * 8` bytes of raw little-endian f64.
    pub fn write(&self, xmf_path: &Path, bin_path: &Path) -> Result<(), Error> {
        if self.nx == 0 || self.ny == 0 || self.nz == 0 {
            return Err(invalid_param("dims", "all grid dimensions must be >= 1"));
        }
        let n_points = self.nx * self.ny * self.nz;

        // --- Binary pass: write all f64 data and record byte offsets ---
        // offsets[t_idx][field_name] = byte offset in bin file
        let mut offsets: Vec<HashMap<String, u64>> = Vec::with_capacity(self.time_steps.len());
        {
            let bin_file = File::create(bin_path)?;
            let mut bin_writer = BufWriter::new(bin_file);
            let mut cursor: u64 = 0;

            for step in &self.time_steps {
                let mut step_offsets: HashMap<String, u64> = HashMap::new();
                let mut sorted_names: Vec<&String> = step.fields.keys().collect();
                sorted_names.sort();

                for name in sorted_names {
                    step_offsets.insert(name.clone(), cursor);
                    let field = &step.fields[name];
                    cursor += write_field_binary(&mut bin_writer, field, n_points)?;
                }
                offsets.push(step_offsets);
            }
            bin_writer.flush()?;
        }

        // --- XML pass ---
        let xmf_file = File::create(xmf_path)?;
        let mut w = BufWriter::new(xmf_file);

        // Use a relative path for the binary reference inside the XML so the
        // pair of files remains portable.
        let bin_name = bin_path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("data.bin");

        writeln!(w, "<?xml version=\"1.0\" ?>")?;
        writeln!(w, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>")?;
        writeln!(w, "<Xdmf Version=\"3.0\">")?;
        writeln!(w, "  <Domain>")?;
        writeln!(
            w,
            "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">"
        )?;

        for (t_idx, step) in self.time_steps.iter().enumerate() {
            let step_offsets = &offsets[t_idx];
            self.write_time_step_xml(&mut w, t_idx, step, step_offsets, bin_name, n_points)?;
        }

        writeln!(w, "    </Grid>")?;
        writeln!(w, "  </Domain>")?;
        writeln!(w, "</Xdmf>")?;
        w.flush()?;

        Ok(())
    }

    fn write_time_step_xml(
        &self,
        w: &mut BufWriter<File>,
        t_idx: usize,
        step: &XdmfTimeStep,
        step_offsets: &HashMap<String, u64>,
        bin_name: &str,
        n_points: usize,
    ) -> Result<(), Error> {
        let [ox, oy, oz] = self.origin;
        let [dx, dy, dz] = self.spacing;
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);

        writeln!(
            w,
            "      <Grid Name=\"Mesh_t{}\" GridType=\"Uniform\">",
            t_idx
        )?;
        writeln!(w, "        <Time Value=\"{}\"/>", step.time)?;
        writeln!(
            w,
            "        <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"{} {} {}\"/>",
            nz, ny, nx
        )?;
        writeln!(w, "        <Geometry GeometryType=\"ORIGIN_DXDYDZ\">")?;
        writeln!(
            w,
            "          <DataItem Format=\"XML\" Dimensions=\"3\">{} {} {}</DataItem>",
            oz, oy, ox
        )?;
        writeln!(
            w,
            "          <DataItem Format=\"XML\" Dimensions=\"3\">{} {} {}</DataItem>",
            dz, dy, dx
        )?;
        writeln!(w, "        </Geometry>")?;

        // Attributes sorted for determinism
        let mut sorted_names: Vec<&String> = step.fields.keys().collect();
        sorted_names.sort();

        for field_name in sorted_names {
            let field = &step.fields[field_name];
            let offset = step_offsets.get(field_name).copied().unwrap_or(0);

            match field {
                XdmfFieldData::Vector(_) => {
                    writeln!(
                        w,
                        "        <Attribute Name=\"{}\" AttributeType=\"Vector\" Center=\"Node\">",
                        field_name
                    )?;
                    writeln!(
                        w,
                        "          <DataItem Format=\"Binary\" DataType=\"Float\" \
                         Precision=\"8\" Dimensions=\"{} 3\" Endian=\"Little\" Seek=\"{}\">",
                        n_points, offset
                    )?;
                    writeln!(w, "            {}", bin_name)?;
                    writeln!(w, "          </DataItem>")?;
                    writeln!(w, "        </Attribute>")?;
                },
                XdmfFieldData::Scalar(_) => {
                    writeln!(
                        w,
                        "        <Attribute Name=\"{}\" AttributeType=\"Scalar\" Center=\"Node\">",
                        field_name
                    )?;
                    writeln!(
                        w,
                        "          <DataItem Format=\"Binary\" DataType=\"Float\" \
                         Precision=\"8\" Dimensions=\"{}\" Endian=\"Little\" Seek=\"{}\">",
                        n_points, offset
                    )?;
                    writeln!(w, "            {}", bin_name)?;
                    writeln!(w, "          </DataItem>")?;
                    writeln!(w, "        </Attribute>")?;
                },
            }
        }

        writeln!(w, "      </Grid>")?;
        Ok(())
    }
}

/// Write one field to the binary writer; return bytes written.
fn write_field_binary(
    w: &mut BufWriter<File>,
    field: &XdmfFieldData,
    n_points: usize,
) -> Result<u64, Error> {
    match field {
        XdmfFieldData::Vector(vecs) => {
            if vecs.len() != n_points {
                return Err(invalid_param(
                    "field",
                    &format!("expected {} vectors, got {}", n_points, vecs.len()),
                ));
            }
            for v in vecs {
                w.write_all(&v.x.to_le_bytes())?;
                w.write_all(&v.y.to_le_bytes())?;
                w.write_all(&v.z.to_le_bytes())?;
            }
            Ok((n_points * 3 * 8) as u64)
        },
        XdmfFieldData::Scalar(vals) => {
            if vals.len() != n_points {
                return Err(invalid_param(
                    "field",
                    &format!("expected {} scalars, got {}", n_points, vals.len()),
                ));
            }
            for &v in vals {
                w.write_all(&v.to_le_bytes())?;
            }
            Ok((n_points * 8) as u64)
        },
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    fn temp_path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(name)
    }

    #[test]
    fn test_new_uniform_grid() {
        let w = XdmfWriter::new_uniform_grid(4, 4, 2, [0.0; 3], [1e-9; 3]);
        assert_eq!(w.nx, 4);
        assert_eq!(w.ny, 4);
        assert_eq!(w.nz, 2);
        assert!(w.time_steps.is_empty());
    }

    #[test]
    fn test_write_single_vector_step() {
        let mut w = XdmfWriter::new_uniform_grid(2, 2, 2, [0.0; 3], [1.0; 3]);
        let n = 8;
        let vecs: Vec<Vector3<f64>> = (0..n)
            .map(|i| Vector3::new(i as f64 * 0.1, 0.0, 1.0))
            .collect();
        let mut fields = HashMap::new();
        fields.insert("magnetization".to_string(), XdmfFieldData::Vector(vecs));
        w.add_time_step(0.0, fields);

        let xmf_path = temp_path("xdmf_single.xmf");
        let bin_path = temp_path("xdmf_single.bin");

        w.write(&xmf_path, &bin_path).expect("write should succeed");
        assert!(xmf_path.exists());
        assert!(bin_path.exists());

        let xml = fs::read_to_string(&xmf_path).unwrap();
        assert!(xml.contains("Xdmf"));
        assert!(xml.contains("TimeSeries"));
        assert!(xml.contains("magnetization"));
        assert!(xml.contains("Vector"));

        let _ = fs::remove_file(&xmf_path);
        let _ = fs::remove_file(&bin_path);
    }

    #[test]
    fn test_binary_file_size_vector() {
        let n_points = 2 * 3 * 4; // 24
        let mut w = XdmfWriter::new_uniform_grid(2, 3, 4, [0.0; 3], [1.0; 3]);
        let vecs: Vec<Vector3<f64>> = (0..n_points).map(|_| Vector3::new(0.0, 0.0, 1.0)).collect();
        let mut fields = HashMap::new();
        fields.insert("spin".to_string(), XdmfFieldData::Vector(vecs));
        w.add_time_step(0.0, fields);

        let xmf_path = temp_path("xdmf_binsize_vec.xmf");
        let bin_path = temp_path("xdmf_binsize_vec.bin");
        w.write(&xmf_path, &bin_path).expect("write should succeed");

        let meta = fs::metadata(&bin_path).unwrap();
        // n_points * 3 components * 8 bytes each
        assert_eq!(meta.len(), (n_points * 3 * 8) as u64);

        let _ = fs::remove_file(&xmf_path);
        let _ = fs::remove_file(&bin_path);
    }

    #[test]
    fn test_binary_file_size_scalar() {
        let n_points = 3 * 3 * 3; // 27
        let mut w = XdmfWriter::new_uniform_grid(3, 3, 3, [0.0; 3], [1.0; 3]);
        let vals: Vec<f64> = (0..n_points).map(|i| i as f64).collect();
        let mut fields = HashMap::new();
        fields.insert("energy".to_string(), XdmfFieldData::Scalar(vals));
        w.add_time_step(1.0, fields);

        let xmf_path = temp_path("xdmf_binsize_scalar.xmf");
        let bin_path = temp_path("xdmf_binsize_scalar.bin");
        w.write(&xmf_path, &bin_path).expect("write should succeed");

        let meta = fs::metadata(&bin_path).unwrap();
        assert_eq!(meta.len(), (n_points * 8) as u64);

        let _ = fs::remove_file(&xmf_path);
        let _ = fs::remove_file(&bin_path);
    }

    #[test]
    fn test_multiple_time_steps() {
        let mut w = XdmfWriter::new_uniform_grid(2, 2, 1, [0.0; 3], [1.0; 3]);
        let n = 4;
        for step_idx in 0..3 {
            let vecs: Vec<Vector3<f64>> = (0..n)
                .map(|_| Vector3::new(step_idx as f64, 0.0, 1.0))
                .collect();
            let mut fields = HashMap::new();
            fields.insert("m".to_string(), XdmfFieldData::Vector(vecs));
            w.add_time_step(step_idx as f64 * 1e-12, fields);
        }

        let xmf_path = temp_path("xdmf_multi_step.xmf");
        let bin_path = temp_path("xdmf_multi_step.bin");
        w.write(&xmf_path, &bin_path).expect("write should succeed");

        let xml = fs::read_to_string(&xmf_path).unwrap();
        assert!(xml.contains("Mesh_t0"));
        assert!(xml.contains("Mesh_t1"));
        assert!(xml.contains("Mesh_t2"));

        // Binary should contain 3 steps * 4 points * 3 components * 8 bytes
        let meta = fs::metadata(&bin_path).unwrap();
        assert_eq!(meta.len(), (3 * 4 * 3 * 8) as u64);

        let _ = fs::remove_file(&xmf_path);
        let _ = fs::remove_file(&bin_path);
    }

    #[test]
    fn test_seek_offsets_in_xml() {
        let mut w = XdmfWriter::new_uniform_grid(2, 2, 1, [0.0; 3], [1.0; 3]);
        let n = 4;

        // Two fields per time step: scalar then vector (sorted alphabetically)
        let mut fields = HashMap::new();
        fields.insert(
            "energy".to_string(),
            XdmfFieldData::Scalar((0..n).map(|i| i as f64).collect()),
        );
        fields.insert(
            "spin".to_string(),
            XdmfFieldData::Vector((0..n).map(|_| Vector3::new(0.0, 0.0, 1.0)).collect()),
        );
        w.add_time_step(0.0, fields);

        let xmf_path = temp_path("xdmf_offsets.xmf");
        let bin_path = temp_path("xdmf_offsets.bin");
        w.write(&xmf_path, &bin_path).expect("write should succeed");

        let xml = fs::read_to_string(&xmf_path).unwrap();
        // energy comes first (alphabetical), offset 0
        assert!(xml.contains("Seek=\"0\""));
        // spin comes after energy: offset = n*1*8 = 32
        assert!(xml.contains(&format!("Seek=\"{}\"", n * 8)));

        let _ = fs::remove_file(&xmf_path);
        let _ = fs::remove_file(&bin_path);
    }

    #[test]
    fn test_xml_topology_dimensions() {
        let mut w = XdmfWriter::new_uniform_grid(5, 4, 3, [0.0; 3], [1.0; 3]);
        let n = 5 * 4 * 3;
        let mut fields = HashMap::new();
        fields.insert(
            "m".to_string(),
            XdmfFieldData::Vector((0..n).map(|_| Vector3::new(0.0, 0.0, 1.0)).collect()),
        );
        w.add_time_step(0.0, fields);

        let xmf_path = temp_path("xdmf_topo.xmf");
        let bin_path = temp_path("xdmf_topo.bin");
        w.write(&xmf_path, &bin_path).expect("write should succeed");

        let xml = fs::read_to_string(&xmf_path).unwrap();
        // Dimensions are nz ny nx
        assert!(xml.contains("Dimensions=\"3 4 5\""));

        let _ = fs::remove_file(&xmf_path);
        let _ = fs::remove_file(&bin_path);
    }
}
