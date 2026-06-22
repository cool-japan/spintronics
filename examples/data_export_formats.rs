//! 3D Magnetic Texture Export: VTK, VTI, NetCDF, and Zarr Formats
//!
//! **Difficulty**: ⭐⭐
//! **Category**: Data Export and Visualisation
//! **Physics**: Skyrmion-like magnetic textures, 3D spin configurations
//!
//! ## Background
//!
//! Scientific data portability requires formats understood by different analysis
//! pipelines.  This example writes a single skyrmion-like spin texture on an
//! 8×8×4 grid to four file formats:
//!
//! - **VTU** — VTK UnstructuredGrid XML (existing, ParaView)
//! - **VTI** — VTK ImageData XML with base64 binary (new in v0.6.0, ParaView)
//! - **NetCDF3 Classic** — pure-Rust XDR binary (new in v0.6.0, Xarray/NCO)
//! - **Zarr v2** — pure-Rust chunk store (new in v0.6.0, Zarr/Python)
//!
//! The magnetic texture used is a Bloch-type skyrmion profile:
//!
//! ```text
//! m_z  = cos(r/R)
//! m_xy = sin(r/R) · (ŷ, -x̂) / r   (Bloch-type rotation)
//! ```
//!
//! where `R = 4` (in lattice units) and `r` is the radial distance from the
//! centre of each z-layer.
//!
//! ## References
//!
//! - VTK File Formats for Visualization Toolkit, Kitware Inc.
//! - XDMF specification v3.0
//! - NetCDF CF Conventions 1.10
//! - Zarr v2 Specification

use std::fs;
use std::path::PathBuf;

#[cfg(feature = "vti")]
use spintronics::prelude::VtiWriter;
use spintronics::prelude::*;
#[cfg(feature = "netcdf")]
use spintronics::prelude::{NetCdfReader, NetCdfWriter};
#[cfg(feature = "zarr")]
use spintronics::prelude::{ZarrDtype, ZarrStore};
#[cfg(feature = "netcdf")]
use spintronics::visualization::netcdf::NetCdfData;
#[cfg(feature = "zarr")]
use spintronics::visualization::zarr::read_array_f64;

// ──────────────────────────────────────────────────────────────────────────────
// Grid dimensions
// ──────────────────────────────────────────────────────────────────────────────
const NX: usize = 8;
const NY: usize = 8;
const NZ: usize = 4;
const N_TOTAL: usize = NX * NY * NZ;

/// Lattice constant (physical spacing between grid points).
const A_LAT: f64 = 1.0e-9; // 1 nm

/// Skyrmion radius in lattice units.
const SKYRMION_R: f64 = 4.0;

// ──────────────────────────────────────────────────────────────────────────────
// Build the skyrmion-like texture
// ──────────────────────────────────────────────────────────────────────────────
fn build_skyrmion_texture() -> Vec<Vector3<f64>> {
    // Centre of each z-layer in fractional lattice units.
    let cx = (NX as f64 - 1.0) * 0.5;
    let cy = (NY as f64 - 1.0) * 0.5;

    let mut spins = Vec::with_capacity(N_TOTAL);
    for iz in 0..NZ {
        for iy in 0..NY {
            for ix in 0..NX {
                let dx = ix as f64 - cx;
                let dy = iy as f64 - cy;
                let r = (dx * dx + dy * dy).sqrt();

                // Skyrmion profile
                let theta_sky = r / SKYRMION_R; // polar angle w.r.t. z-axis
                let m_z = theta_sky.cos();

                // In-plane component: Bloch-type rotation (perpendicular to r)
                let m_xy_norm = theta_sky.sin();
                let (m_x, m_y) = if r < 1e-10 {
                    // At the core: spin points along +z
                    (0.0, 0.0)
                } else {
                    // Bloch chirality: rotate the r-hat by 90° in-plane
                    let r_hat_x = dx / r;
                    let r_hat_y = dy / r;
                    // Bloch: phi rotated 90° → (-r_hat_y, r_hat_x)
                    (m_xy_norm * (-r_hat_y), m_xy_norm * r_hat_x)
                };

                // Clamp to unit sphere (floating-point safety)
                let norm = (m_x * m_x + m_y * m_y + m_z * m_z).sqrt().max(1e-15);
                spins.push(Vector3::new(m_x / norm, m_y / norm, m_z / norm));
                let _ = iz; // iz does not alter the in-plane profile for this simple model
            }
        }
    }
    spins
}

// ──────────────────────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────────────────────
fn file_size_bytes(path: &PathBuf) -> u64 {
    fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}

// Used only by the vti / netcdf / zarr feature blocks below; suppress the
// dead-code lint solely in the configuration where none of them are enabled.
#[cfg_attr(
    not(any(feature = "vti", feature = "netcdf", feature = "zarr")),
    allow(dead_code)
)]
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(name)
}

// ──────────────────────────────────────────────────────────────────────────────
// main
// ──────────────────────────────────────────────────────────────────────────────
fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    // ─────────────────────────────────────────────────────────────────────────
    // Section 1: Generate the 3D magnetic texture
    // ─────────────────────────────────────────────────────────────────────────
    println!("======================================================");
    println!(" Section 1: Generate 3D skyrmion-like texture");
    println!("======================================================");

    let spins = build_skyrmion_texture();

    println!("Grid: {}×{}×{} = {} points", NX, NY, NZ, N_TOTAL);
    println!(
        "Texture: Bloch-type skyrmion, R={} lattice units",
        SKYRMION_R as usize
    );
    println!("Lattice constant: {:.2e} m ({:.0} nm)", A_LAT, A_LAT * 1e9);

    // Print a quick summary of the texture at the centre and edge.
    let centre_idx = (NZ / 2) * NX * NY + (NY / 2) * NX + NX / 2;
    let edge_idx = (NZ / 2) * NX * NY;
    let s_c = &spins[centre_idx];
    let s_e = &spins[edge_idx];
    println!(
        "  Spin at centre ({},{},{}): ({:.3}, {:.3}, {:.3})",
        NX / 2,
        NY / 2,
        NZ / 2,
        s_c.x,
        s_c.y,
        s_c.z
    );
    println!(
        "  Spin at edge   (0,0,{}):  ({:.3}, {:.3}, {:.3})",
        NZ / 2,
        s_e.x,
        s_e.y,
        s_e.z
    );

    // ─────────────────────────────────────────────────────────────────────────
    // Section 2: Export to VTK (VTU unstructured grid)
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n======================================================");
    println!(" Section 2: VTK (VTU) export");
    println!("======================================================");

    // VtkWriter writes to `{base}_{step:05}.vtu` in the current directory.
    // We direct it to a temp-dir prefix.
    let vtk_base = std::env::temp_dir()
        .join("spintronics_export_vtk")
        .to_string_lossy()
        .into_owned();
    let mut vtk_writer = VtkWriter::new(&vtk_base);
    vtk_writer.write_snapshot(&spins, (NX, NY, NZ))?;

    let vtk_path = PathBuf::from(format!("{}_00000.vtu", vtk_base));
    let vtk_size = file_size_bytes(&vtk_path);
    println!("VTU file: {}", vtk_path.display());
    println!(
        "  size: {} bytes ({:.1} KiB)",
        vtk_size,
        vtk_size as f64 / 1024.0
    );

    // ─────────────────────────────────────────────────────────────────────────
    // Section 3: Export to VTI (VTK ImageData, requires `vti` feature)
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n======================================================");
    println!(" Section 3: VTI export (VTK ImageData, base64 binary)");
    println!("======================================================");

    #[cfg(feature = "vti")]
    {
        let vti_path = temp_path("spintronics_export.vti");
        // VtiWriter::with_uniform_spacing: uniform cubic lattice, spacing = A_LAT.
        let vti_writer = VtiWriter::with_uniform_spacing(NX, NY, NZ, A_LAT);
        vti_writer.write_vector_field(&vti_path, "magnetization", &spins)?;

        let vti_size = file_size_bytes(&vti_path);
        println!("VTI file: {}", vti_path.display());
        println!("  grid spacing: {:.2e} m", A_LAT);
        println!(
            "  size: {} bytes ({:.1} KiB)",
            vti_size,
            vti_size as f64 / 1024.0
        );
        println!("  format: XML + base64-encoded f32 binary (ParaView-ready)");
    }
    #[cfg(not(feature = "vti"))]
    {
        println!("(skipped — build with `--features vti` to enable VTI export)");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Section 4: Export to NetCDF3 (requires `netcdf` feature)
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n======================================================");
    println!(" Section 4: NetCDF3 Classic export");
    println!("======================================================");

    #[cfg(feature = "netcdf")]
    {
        let nc_path = temp_path("spintronics_export.nc");

        let mut nc_writer = NetCdfWriter::new();

        // CF-convention global attributes.
        nc_writer.add_global_attribute("Conventions", "CF-1.10");
        nc_writer.add_global_attribute("institution", "COOLJAPAN — spintronics v0.6.0");
        nc_writer.add_global_attribute("title", "Skyrmion-like magnetisation texture");

        // Dimensions: x, y, z points.
        nc_writer.add_dimension("x", NX);
        nc_writer.add_dimension("y", NY);
        nc_writer.add_dimension("z", NZ);
        // Flat "points" dimension for the component variables.
        nc_writer.add_dimension("points", N_TOTAL);

        // Decompose vector field into three scalar variables (CF convention).
        // Flatten in z→y→x order (matching grid loop above).
        let mx: Vec<f64> = spins.iter().map(|s| s.x).collect();
        let my: Vec<f64> = spins.iter().map(|s| s.y).collect();
        let mz: Vec<f64> = spins.iter().map(|s| s.z).collect();

        nc_writer.add_variable(
            "magnetization_x",
            vec!["points"],
            "1",
            "x-component of unit magnetisation",
            NetCdfData::Float64(mx),
        );
        nc_writer.add_variable(
            "magnetization_y",
            vec!["points"],
            "1",
            "y-component of unit magnetisation",
            NetCdfData::Float64(my),
        );
        nc_writer.add_variable(
            "magnetization_z",
            vec!["points"],
            "1",
            "z-component of unit magnetisation",
            NetCdfData::Float64(mz),
        );

        nc_writer.write(&nc_path)?;

        let nc_size = file_size_bytes(&nc_path);
        let var_names = NetCdfReader::list_variables(&nc_path)?;

        println!("NetCDF3 file: {}", nc_path.display());
        println!("  format: NetCDF3 Classic (CDF-1, XDR big-endian)");
        println!("  variables ({} total): {:?}", var_names.len(), var_names);
        println!(
            "  size: {} bytes ({:.1} KiB)",
            nc_size,
            nc_size as f64 / 1024.0
        );
        println!("  readable with: Xarray, NCO, CDO, Ferret, …");

        // Quick sanity check: read back m_z at the centre point.
        let mz_read = NetCdfReader::read_variable_f64(&nc_path, "magnetization_z")?;
        println!(
            "  verification: m_z[{}] (centre) = {:.4}  (written: {:.4})",
            centre_idx, mz_read[centre_idx], spins[centre_idx].z
        );
    }
    #[cfg(not(feature = "netcdf"))]
    {
        println!("(skipped — build with `--features netcdf` to enable NetCDF export)");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Section 5: Export to Zarr v2 (requires `zarr` feature)
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n======================================================");
    println!(" Section 5: Zarr v2 chunk store export");
    println!("======================================================");

    #[cfg(feature = "zarr")]
    {
        let zarr_root = temp_path("spintronics_export_zarr");

        // Create the Zarr store at zarr_root.
        let mut store = ZarrStore::new_store(&zarr_root)?;

        // We store three 1-D arrays of shape [N_TOTAL], chunked into blocks of 64.
        let chunk_size = 64.min(N_TOTAL);
        let shape = vec![N_TOTAL];
        let chunks = vec![chunk_size];

        store.add_array(
            "magnetization_x",
            shape.clone(),
            chunks.clone(),
            ZarrDtype::Float64,
        );
        store.add_array(
            "magnetization_y",
            shape.clone(),
            chunks.clone(),
            ZarrDtype::Float64,
        );
        store.add_array(
            "magnetization_z",
            shape.clone(),
            chunks.clone(),
            ZarrDtype::Float64,
        );

        // Write metadata (.zarray files) for all arrays.
        store.write_metadata()?;

        // Write the data for each component.
        let mx_data: Vec<f64> = spins.iter().map(|s| s.x).collect();
        let my_data: Vec<f64> = spins.iter().map(|s| s.y).collect();
        let mz_data: Vec<f64> = spins.iter().map(|s| s.z).collect();

        store.write_array_data("magnetization_x", &mx_data)?;
        store.write_array_data("magnetization_y", &my_data)?;
        store.write_array_data("magnetization_z", &mz_data)?;

        // Count chunk files created (all files under zarr_root excluding .zarray).
        let chunk_count = count_chunk_files(&zarr_root);
        let zarr_total_bytes: u64 = dir_total_bytes(&zarr_root);

        println!("Zarr store: {}", zarr_root.display());
        println!("  format: Zarr v2 (no compression, C-order, <f8 little-endian)");
        println!("  arrays created: 3 (magnetization_x, _y, _z)");
        println!(
            "  chunk size: {} elements × 8 bytes = {} bytes/chunk",
            chunk_size,
            chunk_size * 8
        );
        println!("  chunk files written: {}", chunk_count);
        println!(
            "  total store size: {} bytes ({:.1} KiB)",
            zarr_total_bytes,
            zarr_total_bytes as f64 / 1024.0
        );
        println!("  readable with: zarr-python, zarr-js, …");

        // Quick sanity check: read back m_z.
        let mz_zarr = read_array_f64(&zarr_root, "magnetization_z")?;
        println!(
            "  verification: m_z[{}] (centre) = {:.4}  (written: {:.4})",
            centre_idx, mz_zarr[centre_idx], spins[centre_idx].z
        );
    }
    #[cfg(not(feature = "zarr"))]
    {
        println!("(skipped — build with `--features zarr` to enable Zarr export)");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Summary
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n======================================================");
    println!(" Summary");
    println!("======================================================");
    println!("All formats written to: {}", std::env::temp_dir().display());
    println!("Format comparison:");
    println!("  VTU  — best for ParaView unstructured grids");
    println!("  VTI  — best for ParaView regular structured grids (binary efficient)");
    println!("  NC   — best for Xarray/Python scientific analysis (CF conventions)");
    println!("  Zarr — best for large chunked cloud-native data stores");

    Ok(())
}

// ──────────────────────────────────────────────────────────────────────────────
// Private helpers
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(feature = "zarr")]
fn count_chunk_files(root: &std::path::Path) -> usize {
    let mut count = 0usize;
    if let Ok(entries) = fs::read_dir(root) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                count += count_chunk_files(&path);
            } else {
                let name = path.file_name().unwrap_or_default().to_string_lossy();
                // Chunk files are named like "0", "1", "0.0", "1.0", etc.
                if name
                    .chars()
                    .next()
                    .map(|c| c.is_ascii_digit())
                    .unwrap_or(false)
                {
                    count += 1;
                }
            }
        }
    }
    count
}

#[cfg(feature = "zarr")]
fn dir_total_bytes(root: &std::path::Path) -> u64 {
    let mut total = 0u64;
    if let Ok(entries) = fs::read_dir(root) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                total += dir_total_bytes(&path);
            } else if let Ok(m) = fs::metadata(&path) {
                total += m.len();
            }
        }
    }
    total
}
