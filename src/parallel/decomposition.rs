//! Domain decomposition for parallel spin dynamics simulations
//!
//! This module provides spatial decomposition of spin arrays into domains
//! that can be evolved independently in parallel. Ghost cells at domain
//! boundaries carry exchange coupling information between neighbors.
//!
//! ## 1D Decomposition
//!
//! A linear spin chain is split into `N` contiguous domains, each with
//! ghost cells of configurable width at left and right boundaries:
//!
//! ```text
//! Domain 0          Domain 1          Domain 2
//! [g|  spins  |g]   [g|  spins  |g]   [g|  spins  |g]
//!        ↕                ↕                ↕
//!    ghost exchange   ghost exchange   ghost exchange
//! ```
//!
//! ## 2D Decomposition
//!
//! A 2D lattice is decomposed into a grid of rectangular sub-domains,
//! each with ghost regions on all four sides.
//!
//! ## References
//! - M. J. Donahue & D. G. Porter, "OOMMF User's Guide", NIST (2002)
//! - A. Vansteenkiste et al., "MuMax3", AIP Advances 4, 107133 (2014)

use rayon::prelude::*;

use crate::constants::GAMMA;
use crate::error::{Error, Result};
use crate::vector3::Vector3;

/// Boundary data for 2D tile ghost cell exchange: (top_row, bottom_row, left_col, right_col).
type TileBoundary = (
    Vec<Vector3<f64>>,
    Vec<Vector3<f64>>,
    Vec<Vector3<f64>>,
    Vec<Vector3<f64>>,
);

// ---------------------------------------------------------------------------
// 1D Domain Decomposition
// ---------------------------------------------------------------------------

/// A single domain in a 1D decomposition, holding interior spins and
/// ghost cells copied from neighboring domains.
#[derive(Debug, Clone)]
pub struct Domain {
    /// Interior spin vectors for this domain
    pub spins: Vec<Vector3<f64>>,
    /// Ghost cells copied from the left neighbor (empty for leftmost domain)
    pub ghost_left: Vec<Vector3<f64>>,
    /// Ghost cells copied from the right neighbor (empty for rightmost domain)
    pub ghost_right: Vec<Vector3<f64>>,
    /// Index of the first interior spin in the global array
    pub global_offset: usize,
}

impl Domain {
    /// Total number of cells visible to this domain (ghosts + interior)
    pub fn total_cells(&self) -> usize {
        self.ghost_left.len() + self.spins.len() + self.ghost_right.len()
    }

    /// Access a spin by local index where index 0 is the first ghost-left cell.
    /// Returns `None` if out of range.
    pub fn get_spin(&self, local_idx: usize) -> Option<&Vector3<f64>> {
        let gl = self.ghost_left.len();
        let ns = self.spins.len();
        if local_idx < gl {
            self.ghost_left.get(local_idx)
        } else if local_idx < gl + ns {
            self.spins.get(local_idx - gl)
        } else {
            self.ghost_right.get(local_idx - gl - ns)
        }
    }
}

/// 1D domain decomposition of a spin array.
///
/// Splits a linear array of spins into `num_domains` contiguous chunks,
/// each augmented with `ghost_width` ghost cells at each boundary.
#[derive(Debug, Clone)]
pub struct DomainDecomposition {
    /// Number of domains
    pub num_domains: usize,
    /// Number of interior spins per domain (last domain may differ)
    pub domain_size: usize,
    /// Width of ghost region on each side of a domain
    pub ghost_width: usize,
    /// The decomposed domains
    pub domains: Vec<Domain>,
}

impl DomainDecomposition {
    /// Create a new 1D domain decomposition from a global spin array.
    ///
    /// # Arguments
    /// * `spins` - The global spin array to decompose
    /// * `num_domains` - Number of domains to create
    /// * `ghost_width` - Number of ghost cells on each boundary
    ///
    /// # Errors
    /// Returns an error if `num_domains` is zero, `ghost_width` is larger than
    /// any domain, or the spin array is empty.
    pub fn new(spins: &[Vector3<f64>], num_domains: usize, ghost_width: usize) -> Result<Self> {
        if num_domains == 0 {
            return Err(Error::InvalidParameter {
                param: "num_domains".to_string(),
                reason: "must be at least 1".to_string(),
            });
        }
        if spins.is_empty() {
            return Err(Error::InvalidParameter {
                param: "spins".to_string(),
                reason: "spin array must not be empty".to_string(),
            });
        }

        let n = spins.len();
        let base_size = n / num_domains;
        let remainder = n % num_domains;

        // Validate ghost width fits inside smallest domain
        let min_domain = base_size;
        if ghost_width > min_domain && num_domains > 1 {
            return Err(Error::InvalidParameter {
                param: "ghost_width".to_string(),
                reason: format!(
                    "ghost_width {} exceeds smallest domain size {}",
                    ghost_width, min_domain
                ),
            });
        }

        let mut domains = Vec::with_capacity(num_domains);
        let mut offset: usize = 0;

        for d in 0..num_domains {
            // Distribute remainder evenly across the first `remainder` domains
            let size = base_size + if d < remainder { 1 } else { 0 };

            let interior = spins[offset..offset + size].to_vec();

            // Ghost left: copy from global array before this domain
            let ghost_left = if d == 0 {
                Vec::new()
            } else {
                let start = offset.saturating_sub(ghost_width);
                spins[start..offset].to_vec()
            };

            // Ghost right: copy from global array after this domain
            let ghost_right = if d == num_domains - 1 {
                Vec::new()
            } else {
                let end = (offset + size + ghost_width).min(n);
                spins[offset + size..end].to_vec()
            };

            domains.push(Domain {
                spins: interior,
                ghost_left,
                ghost_right,
                global_offset: offset,
            });

            offset += size;
        }

        Ok(Self {
            num_domains,
            domain_size: base_size,
            ghost_width,
            domains,
        })
    }

    /// Total number of interior spins across all domains.
    pub fn total_spins(&self) -> usize {
        self.domains.iter().map(|d| d.spins.len()).sum()
    }

    /// Reconstruct the global spin array from domain interiors.
    pub fn gather(&self) -> Vec<Vector3<f64>> {
        let mut global = Vec::with_capacity(self.total_spins());
        for domain in &self.domains {
            global.extend_from_slice(&domain.spins);
        }
        global
    }

    /// Update ghost cells by copying boundary spins from neighboring domains.
    ///
    /// Each domain's ghost_left is filled from the tail of the left neighbor,
    /// and ghost_right from the head of the right neighbor.  The method is
    /// safe to call from a single thread because it reads from a snapshot of
    /// boundary spins and writes into each domain's ghost regions.
    pub fn update_ghost_cells(&mut self) {
        if self.num_domains <= 1 {
            return;
        }

        // Collect boundary data first (immutable snapshot)
        let boundary_data: Vec<(Vec<Vector3<f64>>, Vec<Vector3<f64>>)> = self
            .domains
            .iter()
            .map(|d| {
                let tail: Vec<Vector3<f64>> = d
                    .spins
                    .iter()
                    .rev()
                    .take(self.ghost_width)
                    .copied()
                    .collect::<Vec<_>>()
                    .into_iter()
                    .rev()
                    .collect();
                let head: Vec<Vector3<f64>> =
                    d.spins.iter().take(self.ghost_width).copied().collect();
                (tail, head)
            })
            .collect();

        // Write ghost cells
        for d in 0..self.num_domains {
            // Ghost left from left neighbor's tail
            if d > 0 {
                self.domains[d].ghost_left = boundary_data[d - 1].0.clone();
            }
            // Ghost right from right neighbor's head
            if d < self.num_domains - 1 {
                self.domains[d].ghost_right = boundary_data[d + 1].1.clone();
            }
        }
    }

    /// Perform one parallel LLG integration step across all domains.
    ///
    /// Each domain independently computes `dm/dt` for its interior spins
    /// using ghost cells for exchange coupling at boundaries, then advances
    /// in time with an Euler step.  Ghost cells are synchronized afterwards.
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `alpha` - Gilbert damping constant
    /// * `exchange_stiffness` - Exchange stiffness constant \[J/m\]
    /// * `cell_size` - Lattice spacing \[m\]
    /// * `dt` - Time step \[s\]
    pub fn parallel_llg_step(
        &mut self,
        h_ext: Vector3<f64>,
        alpha: f64,
        exchange_stiffness: f64,
        cell_size: f64,
        dt: f64,
    ) {
        let ghost_width = self.ghost_width;

        // Parallel evolution: each domain computes and applies dm/dt
        self.domains.par_iter_mut().for_each(|domain| {
            let gl = domain.ghost_left.len();
            let ns = domain.spins.len();

            // Build a contiguous view: ghost_left | spins | ghost_right
            let mut all_spins = Vec::with_capacity(domain.total_cells());
            all_spins.extend_from_slice(&domain.ghost_left);
            all_spins.extend_from_slice(&domain.spins);
            all_spins.extend_from_slice(&domain.ghost_right);

            let mut new_spins = Vec::with_capacity(ns);

            for i in 0..ns {
                let idx = gl + i; // index into all_spins
                let m = all_spins[idx];

                // Exchange field from nearest neighbors
                let exchange_field =
                    compute_exchange_field(&all_spins, idx, exchange_stiffness, cell_size);

                let h_eff = h_ext + exchange_field;

                // LLG equation: dm/dt = -gamma/(1+alpha^2) [m x H + alpha m x (m x H)]
                let dm_dt = llg_torque(m, h_eff, alpha);

                // Euler step + renormalize
                let m_new = (m + dm_dt * dt).normalize();
                new_spins.push(m_new);
            }

            domain.spins = new_spins;
            let _ = (ghost_width, gl); // suppress unused warnings in closure
        });

        // Synchronize ghost cells for next step
        self.update_ghost_cells();
    }

    /// Perform one parallel Heun (improved Euler) LLG step across all domains.
    ///
    /// Two-stage predictor-corrector method for improved accuracy compared
    /// to the plain Euler step in [`Self::parallel_llg_step`].
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `alpha` - Gilbert damping constant
    /// * `exchange_stiffness` - Exchange stiffness \[J/m\]
    /// * `cell_size` - Lattice spacing \[m\]
    /// * `dt` - Time step \[s\]
    pub fn parallel_heun_step(
        &mut self,
        h_ext: Vector3<f64>,
        alpha: f64,
        exchange_stiffness: f64,
        cell_size: f64,
        dt: f64,
    ) {
        let ghost_width = self.ghost_width;

        self.domains.par_iter_mut().for_each(|domain| {
            let gl = domain.ghost_left.len();
            let ns = domain.spins.len();

            // Build contiguous view
            let mut all_spins = Vec::with_capacity(domain.total_cells());
            all_spins.extend_from_slice(&domain.ghost_left);
            all_spins.extend_from_slice(&domain.spins);
            all_spins.extend_from_slice(&domain.ghost_right);

            // Stage 1: compute k1
            let mut k1 = Vec::with_capacity(ns);
            for i in 0..ns {
                let idx = gl + i;
                let m = all_spins[idx];
                let exchange_field =
                    compute_exchange_field(&all_spins, idx, exchange_stiffness, cell_size);
                let h_eff = h_ext + exchange_field;
                k1.push(llg_torque(m, h_eff, alpha));
            }

            // Predictor: m* = normalize(m + k1 * dt)
            let mut predicted = Vec::with_capacity(domain.total_cells());
            predicted.extend_from_slice(&domain.ghost_left);
            for (i, k1_val) in k1.iter().enumerate().take(ns) {
                predicted.push((domain.spins[i] + *k1_val * dt).normalize());
            }
            predicted.extend_from_slice(&domain.ghost_right);

            // Stage 2: compute k2 from predicted state
            let mut k2 = Vec::with_capacity(ns);
            for i in 0..ns {
                let idx = gl + i;
                let m = predicted[idx];
                let exchange_field =
                    compute_exchange_field(&predicted, idx, exchange_stiffness, cell_size);
                let h_eff = h_ext + exchange_field;
                k2.push(llg_torque(m, h_eff, alpha));
            }

            // Corrector: m_new = normalize(m + 0.5*(k1+k2)*dt)
            for i in 0..ns {
                let dm_dt = (k1[i] + k2[i]) * 0.5;
                domain.spins[i] = (domain.spins[i] + dm_dt * dt).normalize();
            }

            let _ = (ghost_width, gl);
        });

        self.update_ghost_cells();
    }
}

// ---------------------------------------------------------------------------
// 2D Domain Decomposition
// ---------------------------------------------------------------------------

/// A single tile in a 2D decomposition.
#[derive(Debug, Clone)]
pub struct Tile2D {
    /// Interior spins stored in row-major order, dimensions `rows x cols`
    pub spins: Vec<Vector3<f64>>,
    /// Number of rows in this tile
    pub rows: usize,
    /// Number of columns in this tile
    pub cols: usize,
    /// Ghost cells on the top edge (one row, width = cols)
    pub ghost_top: Vec<Vector3<f64>>,
    /// Ghost cells on the bottom edge
    pub ghost_bottom: Vec<Vector3<f64>>,
    /// Ghost cells on the left edge (one column, height = rows)
    pub ghost_left: Vec<Vector3<f64>>,
    /// Ghost cells on the right edge
    pub ghost_right: Vec<Vector3<f64>>,
    /// Row index of this tile in the grid
    pub grid_row: usize,
    /// Column index of this tile in the grid
    pub grid_col: usize,
}

impl Tile2D {
    /// Access an interior spin at (row, col).
    pub fn get(&self, row: usize, col: usize) -> Option<&Vector3<f64>> {
        if row < self.rows && col < self.cols {
            Some(&self.spins[row * self.cols + col])
        } else {
            None
        }
    }
}

/// 2D grid decomposition of a rectangular spin lattice.
///
/// The global lattice of size `(global_rows, global_cols)` is tiled into
/// `grid_rows x grid_cols` sub-domains, each with single-cell-wide ghost
/// regions on shared edges.
#[derive(Debug, Clone)]
pub struct DomainDecomposition2D {
    /// Number of tile rows
    pub grid_rows: usize,
    /// Number of tile columns
    pub grid_cols: usize,
    /// Global lattice rows
    pub global_rows: usize,
    /// Global lattice columns
    pub global_cols: usize,
    /// Tiles in row-major order
    pub tiles: Vec<Tile2D>,
}

impl DomainDecomposition2D {
    /// Create a 2D decomposition of a row-major spin lattice.
    ///
    /// # Arguments
    /// * `spins` - Row-major spin array of size `global_rows * global_cols`
    /// * `global_rows` - Number of rows in the global lattice
    /// * `global_cols` - Number of columns in the global lattice
    /// * `grid_rows` - Number of tile rows
    /// * `grid_cols` - Number of tile columns
    ///
    /// # Errors
    /// Returns an error on dimension mismatch or zero grid dimensions.
    pub fn new(
        spins: &[Vector3<f64>],
        global_rows: usize,
        global_cols: usize,
        grid_rows: usize,
        grid_cols: usize,
    ) -> Result<Self> {
        if grid_rows == 0 || grid_cols == 0 {
            return Err(Error::InvalidParameter {
                param: "grid dimensions".to_string(),
                reason: "grid_rows and grid_cols must be at least 1".to_string(),
            });
        }
        if spins.len() != global_rows * global_cols {
            return Err(Error::DimensionMismatch {
                expected: format!(
                    "{}x{} = {}",
                    global_rows,
                    global_cols,
                    global_rows * global_cols
                ),
                actual: format!("{}", spins.len()),
            });
        }

        let base_tile_rows = global_rows / grid_rows;
        let rem_rows = global_rows % grid_rows;
        let base_tile_cols = global_cols / grid_cols;
        let rem_cols = global_cols % grid_cols;

        let idx = |r: usize, c: usize| -> usize { r * global_cols + c };

        let mut tiles = Vec::with_capacity(grid_rows * grid_cols);

        let mut row_offset: usize = 0;
        for gr in 0..grid_rows {
            let tile_rows = base_tile_rows + if gr < rem_rows { 1 } else { 0 };
            let mut col_offset: usize = 0;

            for gc in 0..grid_cols {
                let tile_cols = base_tile_cols + if gc < rem_cols { 1 } else { 0 };

                // Copy interior
                let mut interior = Vec::with_capacity(tile_rows * tile_cols);
                for r in row_offset..row_offset + tile_rows {
                    for c in col_offset..col_offset + tile_cols {
                        interior.push(spins[idx(r, c)]);
                    }
                }

                // Ghost top (row above first interior row)
                let ghost_top = if gr == 0 {
                    Vec::new()
                } else {
                    let r = row_offset - 1;
                    (col_offset..col_offset + tile_cols)
                        .map(|c| spins[idx(r, c)])
                        .collect()
                };

                // Ghost bottom
                let ghost_bottom = if gr == grid_rows - 1 {
                    Vec::new()
                } else {
                    let r = row_offset + tile_rows;
                    (col_offset..col_offset + tile_cols)
                        .map(|c| spins[idx(r, c)])
                        .collect()
                };

                // Ghost left (column to the left of first interior column)
                let ghost_left = if gc == 0 {
                    Vec::new()
                } else {
                    let c = col_offset - 1;
                    (row_offset..row_offset + tile_rows)
                        .map(|r| spins[idx(r, c)])
                        .collect()
                };

                // Ghost right
                let ghost_right = if gc == grid_cols - 1 {
                    Vec::new()
                } else {
                    let c = col_offset + tile_cols;
                    (row_offset..row_offset + tile_rows)
                        .map(|r| spins[idx(r, c)])
                        .collect()
                };

                tiles.push(Tile2D {
                    spins: interior,
                    rows: tile_rows,
                    cols: tile_cols,
                    ghost_top,
                    ghost_bottom,
                    ghost_left,
                    ghost_right,
                    grid_row: gr,
                    grid_col: gc,
                });

                col_offset += tile_cols;
            }
            row_offset += tile_rows;
        }

        Ok(Self {
            grid_rows,
            grid_cols,
            global_rows,
            global_cols,
            tiles,
        })
    }

    /// Total number of interior spins across all tiles.
    pub fn total_spins(&self) -> usize {
        self.tiles.iter().map(|t| t.spins.len()).sum()
    }

    /// Reconstruct the global row-major spin array from tile interiors.
    pub fn gather(&self) -> Vec<Vector3<f64>> {
        let mut global = vec![Vector3::zero(); self.global_rows * self.global_cols];

        for tile in &self.tiles {
            // Compute the row/col offset of this tile in the global lattice
            // by summing tile sizes for tiles above and to the left
            let row_offset = self.row_offset(tile.grid_row);
            let col_offset = self.col_offset(tile.grid_col);

            for r in 0..tile.rows {
                for c in 0..tile.cols {
                    let gr = row_offset + r;
                    let gc = col_offset + c;
                    global[gr * self.global_cols + gc] = tile.spins[r * tile.cols + c];
                }
            }
        }

        global
    }

    /// Update ghost cells for all tiles from neighboring tile boundaries.
    pub fn update_ghost_cells(&mut self) {
        // Snapshot boundary data: for each tile, store (top_row, bottom_row, left_col, right_col)
        let boundary: Vec<TileBoundary> = self
            .tiles
            .iter()
            .map(|t| {
                let top_row: Vec<_> = t.spins.iter().take(t.cols).copied().collect();
                let bottom_row: Vec<_> = t
                    .spins
                    .iter()
                    .skip((t.rows - 1) * t.cols)
                    .take(t.cols)
                    .copied()
                    .collect();
                let left_col: Vec<_> = (0..t.rows).map(|r| t.spins[r * t.cols]).collect();
                let right_col: Vec<_> = (0..t.rows)
                    .map(|r| t.spins[r * t.cols + t.cols - 1])
                    .collect();
                (top_row, bottom_row, left_col, right_col)
            })
            .collect();

        let gc = self.grid_cols;

        for tile in self.tiles.iter_mut() {
            let gr = tile.grid_row;
            let gcol = tile.grid_col;

            // Ghost top from tile above's bottom row
            if gr > 0 {
                let above_idx = (gr - 1) * gc + gcol;
                tile.ghost_top = boundary[above_idx].1.clone();
            }
            // Ghost bottom from tile below's top row
            if gr < self.grid_rows - 1 {
                let below_idx = (gr + 1) * gc + gcol;
                tile.ghost_bottom = boundary[below_idx].0.clone();
            }
            // Ghost left from tile to the left's right column
            if gcol > 0 {
                let left_idx = gr * gc + (gcol - 1);
                tile.ghost_left = boundary[left_idx].3.clone();
            }
            // Ghost right from tile to the right's left column
            if gcol < self.grid_cols - 1 {
                let right_idx = gr * gc + (gcol + 1);
                tile.ghost_right = boundary[right_idx].2.clone();
            }
        }
    }

    /// Helper: compute the global row offset for tile grid-row `gr`.
    fn row_offset(&self, gr: usize) -> usize {
        let base = self.global_rows / self.grid_rows;
        let rem = self.global_rows % self.grid_rows;
        let full = gr.min(rem) * (base + 1);
        let rest = gr.saturating_sub(rem) * base;
        full + rest
    }

    /// Helper: compute the global column offset for tile grid-col `gc`.
    fn col_offset(&self, gc: usize) -> usize {
        let base = self.global_cols / self.grid_cols;
        let rem = self.global_cols % self.grid_cols;
        let full = gc.min(rem) * (base + 1);
        let rest = gc.saturating_sub(rem) * base;
        full + rest
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Compute nearest-neighbor exchange field for spin at `idx` in a 1D array.
///
/// H_exchange = (A / a^2) * (S_{i-1} + S_{i+1} - 2*S_i)
fn compute_exchange_field(
    spins: &[Vector3<f64>],
    idx: usize,
    exchange_stiffness: f64,
    cell_size: f64,
) -> Vector3<f64> {
    let n = spins.len();
    if n < 2 {
        return Vector3::zero();
    }

    let prefactor = exchange_stiffness / (cell_size * cell_size);
    let current = spins[idx];

    let left = if idx > 0 { spins[idx - 1] } else { current };
    let right = if idx < n - 1 { spins[idx + 1] } else { current };

    (left + right - current * 2.0) * prefactor
}

/// Compute the LLG torque dm/dt for a single spin.
///
/// dm/dt = -gamma / (1 + alpha^2) * [m x H_eff + alpha * m x (m x H_eff)]
fn llg_torque(m: Vector3<f64>, h_eff: Vector3<f64>, alpha: f64) -> Vector3<f64> {
    let m_cross_h = m.cross(&h_eff);
    let m_cross_m_cross_h = m.cross(&m_cross_h);
    let prefactor = -GAMMA / (1.0 + alpha * alpha);
    (m_cross_h + m_cross_m_cross_h * alpha) * prefactor
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: create a uniform spin array pointing along +z.
    fn uniform_z_spins(n: usize) -> Vec<Vector3<f64>> {
        vec![Vector3::new(0.0, 0.0, 1.0); n]
    }

    #[test]
    fn test_decomposition_total_size() {
        let spins = uniform_z_spins(100);
        let decomp = DomainDecomposition::new(&spins, 4, 1).expect("decomposition should succeed");

        assert_eq!(decomp.total_spins(), 100);
        assert_eq!(decomp.num_domains, 4);

        // Gather should reproduce the original array
        let gathered = decomp.gather();
        assert_eq!(gathered.len(), 100);
        for (a, b) in gathered.iter().zip(spins.iter()) {
            assert!((a.x - b.x).abs() < 1e-15);
            assert!((a.y - b.y).abs() < 1e-15);
            assert!((a.z - b.z).abs() < 1e-15);
        }
    }

    #[test]
    fn test_decomposition_uneven_split() {
        // 103 spins into 4 domains: 26+26+26+25
        let spins = uniform_z_spins(103);
        let decomp = DomainDecomposition::new(&spins, 4, 1).expect("decomposition should succeed");

        assert_eq!(decomp.total_spins(), 103);

        // First 3 domains get 26, last gets 25
        assert_eq!(decomp.domains[0].spins.len(), 26);
        assert_eq!(decomp.domains[1].spins.len(), 26);
        assert_eq!(decomp.domains[2].spins.len(), 26);
        assert_eq!(decomp.domains[3].spins.len(), 25);
    }

    #[test]
    fn test_ghost_cells_width_1() {
        // Create spins with distinct values so we can verify ghost content
        let spins: Vec<Vector3<f64>> = (0..20).map(|i| Vector3::new(i as f64, 0.0, 0.0)).collect();

        let decomp = DomainDecomposition::new(&spins, 4, 1).expect("decomposition should succeed");

        // Domain 0: no ghost_left, ghost_right = first spin of domain 1
        assert!(decomp.domains[0].ghost_left.is_empty());
        assert_eq!(decomp.domains[0].ghost_right.len(), 1);
        // Domain 0 has spins [0..5], ghost_right should be spin 5
        assert!((decomp.domains[0].ghost_right[0].x - 5.0).abs() < 1e-15);

        // Domain 1: ghost_left = last spin of domain 0, ghost_right = first of domain 2
        assert_eq!(decomp.domains[1].ghost_left.len(), 1);
        assert!((decomp.domains[1].ghost_left[0].x - 4.0).abs() < 1e-15);

        // Last domain: ghost_left present, no ghost_right
        let last = decomp.domains.last().expect("should have domains");
        assert!(last.ghost_right.is_empty());
        assert!(!last.ghost_left.is_empty());
    }

    #[test]
    fn test_ghost_cell_update() {
        let spins: Vec<Vector3<f64>> = (0..20).map(|i| Vector3::new(i as f64, 0.0, 0.0)).collect();

        let mut decomp =
            DomainDecomposition::new(&spins, 4, 1).expect("decomposition should succeed");

        // Modify interior of domain 1 to have new values
        for s in &mut decomp.domains[1].spins {
            s.x += 100.0;
        }

        decomp.update_ghost_cells();

        // Domain 0's ghost_right should now reflect domain 1's first interior spin
        let expected_x = 5.0 + 100.0; // original spin[5] + 100
        assert!(
            (decomp.domains[0].ghost_right[0].x - expected_x).abs() < 1e-15,
            "ghost_right[0] = {}, expected {}",
            decomp.domains[0].ghost_right[0].x,
            expected_x,
        );

        // Domain 2's ghost_left should reflect domain 1's last interior spin
        let d1_last_x = 9.0 + 100.0; // original spin[9] + 100
        assert!(
            (decomp.domains[2].ghost_left[0].x - d1_last_x).abs() < 1e-15,
            "ghost_left[0] = {}, expected {}",
            decomp.domains[2].ghost_left[0].x,
            d1_last_x,
        );
    }

    #[test]
    fn test_parallel_llg_step_matches_serial() {
        // Create a small tilted spin chain
        let n = 40;
        let spins: Vec<Vector3<f64>> = (0..n)
            .map(|i| {
                let angle = 0.1 * (i as f64);
                Vector3::new(angle.sin(), 0.0, angle.cos()).normalize()
            })
            .collect();

        let h_ext = Vector3::new(0.0, 0.0, 1.0);
        let alpha = 0.01;
        let a_ex = 1e-11;
        let cell_size = 1e-9;
        let dt = 1e-14;

        // Serial reference: evolve entire chain as one domain
        let mut serial_decomp =
            DomainDecomposition::new(&spins, 1, 0).expect("single domain decomposition");
        serial_decomp.parallel_llg_step(h_ext, alpha, a_ex, cell_size, dt);
        let serial_result = serial_decomp.gather();

        // Parallel: split into 4 domains with ghost width 1
        let mut par_decomp =
            DomainDecomposition::new(&spins, 4, 1).expect("4-domain decomposition");
        par_decomp.parallel_llg_step(h_ext, alpha, a_ex, cell_size, dt);
        let par_result = par_decomp.gather();

        assert_eq!(serial_result.len(), par_result.len());

        // Interior spins (away from domain boundaries) should match closely
        // Boundary spins may differ slightly due to ghost cell approximation
        let mut max_diff = 0.0_f64;
        for (s, p) in serial_result.iter().zip(par_result.iter()) {
            let diff = (*s - *p).magnitude();
            if diff > max_diff {
                max_diff = diff;
            }
        }

        // For a single Euler step the difference should be very small
        assert!(
            max_diff < 1e-6,
            "max difference between serial and parallel: {:.2e}",
            max_diff,
        );
    }

    #[test]
    fn test_2d_decomposition_total_size() {
        let rows = 12;
        let cols = 15;
        let spins: Vec<Vector3<f64>> = (0..rows * cols)
            .map(|_| Vector3::new(0.0, 0.0, 1.0))
            .collect();

        let decomp = DomainDecomposition2D::new(&spins, rows, cols, 3, 3)
            .expect("2D decomposition should succeed");

        assert_eq!(decomp.total_spins(), rows * cols);

        let gathered = decomp.gather();
        assert_eq!(gathered.len(), rows * cols);
    }

    #[test]
    fn test_2d_ghost_update() {
        let rows = 6;
        let cols = 6;
        let spins: Vec<Vector3<f64>> = (0..rows * cols)
            .map(|i| Vector3::new(i as f64, 0.0, 0.0))
            .collect();

        let mut decomp = DomainDecomposition2D::new(&spins, rows, cols, 2, 2)
            .expect("2D decomposition should succeed");

        // Modify tile (0,0) interior and update ghosts
        for s in &mut decomp.tiles[0].spins {
            s.x += 1000.0;
        }

        decomp.update_ghost_cells();

        // Tile (0,1) ghost_left should now reflect tile (0,0) right column
        let tile_01 = &decomp.tiles[1];
        assert!(
            !tile_01.ghost_left.is_empty(),
            "tile (0,1) should have ghost_left"
        );
        // The ghost values should be the modified values (>= 1000)
        for g in &tile_01.ghost_left {
            assert!(
                g.x >= 1000.0,
                "ghost should reflect modified tile, got {}",
                g.x,
            );
        }
    }

    #[test]
    fn test_error_on_zero_domains() {
        let spins = uniform_z_spins(10);
        let result = DomainDecomposition::new(&spins, 0, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_error_on_empty_spins() {
        let spins: Vec<Vector3<f64>> = Vec::new();
        let result = DomainDecomposition::new(&spins, 2, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_heun_step_preserves_normalization() {
        let n = 20;
        let spins: Vec<Vector3<f64>> = (0..n)
            .map(|i| {
                let angle = 0.2 * (i as f64);
                Vector3::new(angle.sin(), 0.0, angle.cos()).normalize()
            })
            .collect();

        let mut decomp =
            DomainDecomposition::new(&spins, 2, 1).expect("decomposition should succeed");

        decomp.parallel_heun_step(Vector3::new(0.0, 0.0, 1.0), 0.01, 1e-11, 1e-9, 1e-14);

        for domain in &decomp.domains {
            for s in &domain.spins {
                let mag = s.magnitude();
                assert!(
                    (mag - 1.0).abs() < 1e-10,
                    "spin magnitude {} deviates from 1.0",
                    mag,
                );
            }
        }
    }
}
