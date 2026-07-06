//! Frustrated lattice geometries for magnetic systems
//!
//! This module provides lattice generation and manipulation for frustrated
//! magnetic systems, including triangular, kagome, and pyrochlore lattices.
//!
//! # Physics Background
//!
//! Geometric frustration occurs when the lattice geometry prevents simultaneous
//! minimization of all pairwise interactions. The canonical example is an
//! antiferromagnet on a triangular lattice: three spins on a triangle cannot
//! all be antiparallel to each other.
//!
//! The frustration parameter f = |θ_CW| / T_N quantifies the degree of
//! frustration. Materials with f >> 1 are strongly frustrated and may host
//! exotic phases such as spin liquids.
//!
//! # References
//!
//! - A.P. Ramirez, "Strongly Geometrically Frustrated Magnets",
//!   Annu. Rev. Mater. Sci. 24, 453-480 (1994)
//! - L. Balents, "Spin liquids in frustrated magnets",
//!   Nature 464, 199-208 (2010)

use crate::error::{Error, Result};
use crate::vector3::Vector3;

/// Type of frustrated lattice geometry
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LatticeType {
    /// Triangular lattice (2D) - each site has 6 nearest neighbors
    /// Lattice vectors: a1 = (1, 0), a2 = (1/2, √3/2)
    Triangular,

    /// Kagome lattice (2D) - corner-sharing triangles, 3 sublattices
    /// Each site has 4 nearest neighbors
    Kagome,

    /// Pyrochlore lattice (3D) - corner-sharing tetrahedra
    /// Each site has 6 nearest neighbors within tetrahedra
    /// Relevant for spin ice materials (Dy2Ti2O7, Ho2Ti2O7)
    Pyrochlore,
}

/// A frustrated magnetic lattice with spin configurations
///
/// Represents a lattice of interacting magnetic moments (spins) on a
/// geometrically frustrated lattice. The Hamiltonian is the antiferromagnetic
/// Heisenberg model: H = J Σ_{<i,j>} S_i · S_j with J > 0.
#[derive(Debug, Clone)]
pub struct FrustratedLattice {
    /// Type of lattice geometry
    pub lattice_type: LatticeType,
    /// Lattice dimensions (Nx, Ny) for 2D, (Nx, Ny) with Nz layers for 3D
    pub size: (usize, usize),
    /// Exchange coupling constant J \[J\] (positive for antiferromagnetic)
    pub coupling_j: f64,
    /// Spin vectors at each lattice site
    pub spins: Vec<Vector3<f64>>,
    /// Neighbor list: neighbors\[i\] contains indices of nearest neighbors of site i
    pub neighbors: Vec<Vec<usize>>,
    /// Site positions in real space
    pub positions: Vec<Vector3<f64>>,
    /// Lattice constant \[m\]
    pub lattice_constant: f64,
}

/// Deterministic xorshift64 PRNG for Monte Carlo simulations
///
/// Implements the xorshift64 algorithm by George Marsaglia (2003).
/// This avoids external RNG dependencies while providing acceptable
/// statistical quality for physics simulations.
#[derive(Debug, Clone)]
pub struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    /// Create a new PRNG with the given seed
    ///
    /// # Errors
    ///
    /// Returns error if seed is zero (xorshift requires nonzero state).
    pub fn new(seed: u64) -> Result<Self> {
        if seed == 0 {
            return Err(Error::InvalidParameter {
                param: "seed".to_string(),
                reason: "xorshift64 seed must be nonzero".to_string(),
            });
        }
        Ok(Self { state: seed })
    }

    /// Generate next u64 value
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    /// Generate a uniform f64 in [0, 1)
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }

    /// Generate a unit vector on the sphere (Marsaglia method)
    pub fn random_unit_vector(&mut self) -> Vector3<f64> {
        // Marsaglia (1972) method for uniform points on a sphere
        loop {
            let u = 2.0 * self.next_f64() - 1.0;
            let v = 2.0 * self.next_f64() - 1.0;
            let s = u * u + v * v;
            if s < 1.0 && s > 1e-30 {
                let factor = 2.0 * (1.0 - s).sqrt();
                return Vector3::new(u * factor, v * factor, 1.0 - 2.0 * s);
            }
        }
    }
}

impl FrustratedLattice {
    /// Create a new triangular lattice
    ///
    /// Generates a 2D triangular lattice with Nx × Ny sites.
    /// Lattice vectors: a1 = a(1, 0, 0), a2 = a(1/2, √3/2, 0)
    /// Each site has 6 nearest neighbors (with periodic boundary conditions).
    ///
    /// # Arguments
    ///
    /// * `nx` - Number of sites along a1 direction
    /// * `ny` - Number of sites along a2 direction
    /// * `coupling_j` - Exchange coupling J > 0 (antiferromagnetic)
    /// * `lattice_constant` - Lattice spacing \[m\]
    ///
    /// # Errors
    ///
    /// Returns error if dimensions are zero or coupling is non-positive.
    pub fn triangular(
        nx: usize,
        ny: usize,
        coupling_j: f64,
        lattice_constant: f64,
    ) -> Result<Self> {
        validate_lattice_params(nx, ny, coupling_j, lattice_constant)?;

        let n_sites = nx * ny;
        let mut positions = Vec::with_capacity(n_sites);
        let mut neighbors = Vec::with_capacity(n_sites);

        let sqrt3_half = 3.0_f64.sqrt() / 2.0;

        // Generate positions
        for iy in 0..ny {
            for ix in 0..nx {
                let x = (ix as f64 + 0.5 * iy as f64) * lattice_constant;
                let y = iy as f64 * sqrt3_half * lattice_constant;
                positions.push(Vector3::new(x, y, 0.0));
            }
        }

        // Generate neighbor lists with periodic boundary conditions
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = iy * nx + ix;
                let mut nbrs = Vec::with_capacity(6);

                // Six neighbors on triangular lattice:
                // (+1,0), (-1,0), (0,+1), (0,-1), (+1,-1), (-1,+1)
                let offsets: [(i64, i64); 6] = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)];

                for (dx, dy) in offsets {
                    let jx = ((ix as i64 + dx).rem_euclid(nx as i64)) as usize;
                    let jy = ((iy as i64 + dy).rem_euclid(ny as i64)) as usize;
                    let j_idx = jy * nx + jx;
                    if j_idx != idx {
                        nbrs.push(j_idx);
                    }
                }

                // Remove duplicates (can happen for small lattices)
                nbrs.sort_unstable();
                nbrs.dedup();
                let _ = idx; // used above
                neighbors.push(nbrs);
            }
        }

        // Initialize spins pointing along z (will be relaxed by MC)
        let spins = vec![Vector3::unit_z(); n_sites];

        Ok(Self {
            lattice_type: LatticeType::Triangular,
            size: (nx, ny),
            coupling_j,
            spins,
            neighbors,
            positions,
            lattice_constant,
        })
    }

    /// Create a new kagome lattice
    ///
    /// The kagome lattice consists of corner-sharing triangles with 3 sublattices.
    /// Each unit cell contains 3 sites. Total sites = 3 × Nx × Ny.
    /// Each site has 4 nearest neighbors.
    ///
    /// # Arguments
    ///
    /// * `nx` - Number of unit cells along a1
    /// * `ny` - Number of unit cells along a2
    /// * `coupling_j` - Exchange coupling J > 0
    /// * `lattice_constant` - Lattice spacing \[m\]
    ///
    /// # Errors
    ///
    /// Returns error if parameters are invalid.
    pub fn kagome(nx: usize, ny: usize, coupling_j: f64, lattice_constant: f64) -> Result<Self> {
        validate_lattice_params(nx, ny, coupling_j, lattice_constant)?;

        let n_cells = nx * ny;
        let n_sites = 3 * n_cells;
        let mut positions = Vec::with_capacity(n_sites);
        let mut neighbors = Vec::with_capacity(n_sites);

        let sqrt3_half = 3.0_f64.sqrt() / 2.0;
        let a = lattice_constant;

        // Kagome sublattice positions within a unit cell:
        // Sublattice 0: (0, 0)
        // Sublattice 1: (a/2, 0)
        // Sublattice 2: (a/4, a√3/4)
        let sub_offsets = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(a / 2.0, 0.0, 0.0),
            Vector3::new(a / 4.0, a * sqrt3_half / 2.0, 0.0),
        ];

        // Generate positions
        for iy in 0..ny {
            for ix in 0..nx {
                let base_x = (ix as f64 + 0.5 * iy as f64) * a;
                let base_y = iy as f64 * sqrt3_half * a;
                for sub in &sub_offsets {
                    positions.push(Vector3::new(base_x + sub.x, base_y + sub.y, 0.0));
                }
            }
        }

        // Helper to get site index from (cell_x, cell_y, sublattice)
        let site_idx = |cx: usize, cy: usize, sub: usize| -> usize { 3 * (cy * nx + cx) + sub };

        // Generate neighbor lists for kagome lattice
        // Each sublattice site connects to 4 neighbors
        for iy in 0..ny {
            for ix in 0..nx {
                let ixp = (ix + 1) % nx;
                let iyp = (iy + 1) % ny;
                let ixm = (ix + nx - 1) % nx;
                let iym = (iy + ny - 1) % ny;

                // Sublattice 0 neighbors: sub1 in same cell, sub2 in same cell,
                //                         sub1 in (-1,0) cell, sub2 in (0,-1) cell
                let nbrs_0 = vec![
                    site_idx(ix, iy, 1),
                    site_idx(ix, iy, 2),
                    site_idx(ixm, iy, 1),
                    site_idx(ix, iym, 2),
                ];
                neighbors.push(nbrs_0);

                // Sublattice 1 neighbors: sub0 in same cell, sub2 in same cell,
                //                         sub0 in (+1,0) cell, sub2 in (+1,-1) cell
                let nbrs_1 = vec![
                    site_idx(ix, iy, 0),
                    site_idx(ix, iy, 2),
                    site_idx(ixp, iy, 0),
                    site_idx(ixp, iym, 2),
                ];
                neighbors.push(nbrs_1);

                // Sublattice 2 neighbors: sub0 in same cell, sub1 in same cell,
                //                         sub0 in (0,+1) cell, sub1 in (-1,+1) cell
                let nbrs_2 = vec![
                    site_idx(ix, iy, 0),
                    site_idx(ix, iy, 1),
                    site_idx(ix, iyp, 0),
                    site_idx(ixm, iyp, 1),
                ];
                neighbors.push(nbrs_2);
            }
        }

        let spins = vec![Vector3::unit_z(); n_sites];

        Ok(Self {
            lattice_type: LatticeType::Kagome,
            size: (nx, ny),
            coupling_j,
            spins,
            neighbors,
            positions,
            lattice_constant,
        })
    }

    /// Create a new pyrochlore lattice
    ///
    /// The pyrochlore lattice consists of corner-sharing tetrahedra in 3D.
    /// Each cubic unit cell contains 4 sites (FCC with basis).
    /// Total sites = 4 × Nx × Ny × Nz where Nz = min(Nx, Ny).
    ///
    /// This is the lattice relevant for spin ice materials like Dy₂Ti₂O₇.
    ///
    /// # Arguments
    ///
    /// * `nx` - Number of unit cells along x
    /// * `ny` - Number of unit cells along y
    /// * `coupling_j` - Exchange coupling J > 0
    /// * `lattice_constant` - Cubic unit cell size \[m\]
    ///
    /// # Errors
    ///
    /// Returns error if parameters are invalid.
    pub fn pyrochlore(
        nx: usize,
        ny: usize,
        coupling_j: f64,
        lattice_constant: f64,
    ) -> Result<Self> {
        validate_lattice_params(nx, ny, coupling_j, lattice_constant)?;

        let nz = nx.min(ny);
        let n_cells = nx * ny * nz;
        let n_sites = 4 * n_cells;
        let a = lattice_constant;

        // Pyrochlore sublattice positions (corners of tetrahedra in FCC)
        // These are the 4 sites within a conventional cubic unit cell
        let sub_offsets = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(a / 4.0, a / 4.0, 0.0),
            Vector3::new(a / 4.0, 0.0, a / 4.0),
            Vector3::new(0.0, a / 4.0, a / 4.0),
        ];

        let mut positions = Vec::with_capacity(n_sites);
        let mut neighbors = Vec::with_capacity(n_sites);

        // Generate positions
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    for sub in &sub_offsets {
                        positions.push(Vector3::new(
                            ix as f64 * a + sub.x,
                            iy as f64 * a + sub.y,
                            iz as f64 * a + sub.z,
                        ));
                    }
                }
            }
        }

        // Helper to get site index from (cx, cy, cz, sublattice)
        let site_idx = |cx: usize, cy: usize, cz: usize, sub: usize| -> usize {
            4 * ((cz * ny + cy) * nx + cx) + sub
        };

        // Generate neighbor lists
        // Each site in a tetrahedron connects to the other 3 sites in the same
        // tetrahedron, plus 3 sites in neighboring tetrahedra (6 total)
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let ixp = (ix + 1) % nx;
                    let iyp = (iy + 1) % ny;
                    let izp = (iz + 1) % nz;
                    let ixm = (ix + nx - 1) % nx;
                    let iym = (iy + ny - 1) % ny;
                    let izm = (iz + nz - 1) % nz;

                    // Sublattice 0: connects to 1,2,3 in same cell + neighbors
                    let mut nbrs_0 = vec![
                        site_idx(ix, iy, iz, 1),
                        site_idx(ix, iy, iz, 2),
                        site_idx(ix, iy, iz, 3),
                        site_idx(ixm, iym, iz, 1),
                        site_idx(ixm, iy, izm, 2),
                        site_idx(ix, iym, izm, 3),
                    ];
                    nbrs_0.sort_unstable();
                    nbrs_0.dedup();
                    neighbors.push(nbrs_0);

                    // Sublattice 1
                    let mut nbrs_1 = vec![
                        site_idx(ix, iy, iz, 0),
                        site_idx(ix, iy, iz, 2),
                        site_idx(ix, iy, iz, 3),
                        site_idx(ixp, iyp, iz, 0),
                        site_idx(ix, iyp, izm, 3),
                        site_idx(ixp, iy, izm, 2),
                    ];
                    nbrs_1.sort_unstable();
                    nbrs_1.dedup();
                    neighbors.push(nbrs_1);

                    // Sublattice 2
                    let mut nbrs_2 = vec![
                        site_idx(ix, iy, iz, 0),
                        site_idx(ix, iy, iz, 1),
                        site_idx(ix, iy, iz, 3),
                        site_idx(ixp, iy, izp, 0),
                        site_idx(ixp, iy, izm, 1),
                        site_idx(ix, iy, izp, 3),
                    ];
                    nbrs_2.sort_unstable();
                    nbrs_2.dedup();
                    neighbors.push(nbrs_2);

                    // Sublattice 3
                    let mut nbrs_3 = vec![
                        site_idx(ix, iy, iz, 0),
                        site_idx(ix, iy, iz, 1),
                        site_idx(ix, iy, iz, 2),
                        site_idx(ix, iyp, izp, 0),
                        site_idx(ix, iyp, izm, 1),
                        site_idx(ix, iy, izp, 2),
                    ];
                    nbrs_3.sort_unstable();
                    nbrs_3.dedup();
                    neighbors.push(nbrs_3);
                }
            }
        }

        // Initialize spins along local [111] directions (Ising axes for spin ice)
        let ising_axes = [
            Vector3::new(1.0, 1.0, 1.0).normalize(),
            Vector3::new(1.0, -1.0, -1.0).normalize(),
            Vector3::new(-1.0, 1.0, -1.0).normalize(),
            Vector3::new(-1.0, -1.0, 1.0).normalize(),
        ];
        let mut spins = Vec::with_capacity(n_sites);
        for i in 0..n_sites {
            spins.push(ising_axes[i % 4]);
        }

        Ok(Self {
            lattice_type: LatticeType::Pyrochlore,
            size: (nx, ny),
            coupling_j,
            spins,
            neighbors,
            positions,
            lattice_constant,
        })
    }

    /// Number of lattice sites
    pub fn num_sites(&self) -> usize {
        self.spins.len()
    }

    /// Calculate the total energy of the spin configuration
    ///
    /// H = J Σ_{<i,j>} S_i · S_j
    ///
    /// Each bond is counted once.
    pub fn total_energy(&self) -> f64 {
        let mut energy = 0.0;
        for (i, nbrs) in self.neighbors.iter().enumerate() {
            for &j in nbrs {
                if j > i {
                    energy += self.coupling_j * self.spins[i].dot(&self.spins[j]);
                }
            }
        }
        energy
    }

    /// Calculate the energy contribution of a single spin
    ///
    /// E_i = J Σ_{j ∈ neighbors(i)} S_i · S_j
    pub fn site_energy(&self, site: usize) -> f64 {
        let mut energy = 0.0;
        if site < self.neighbors.len() {
            for &j in &self.neighbors[site] {
                energy += self.coupling_j * self.spins[site].dot(&self.spins[j]);
            }
        }
        energy
    }

    /// Calculate the average magnetization vector
    pub fn average_magnetization(&self) -> Vector3<f64> {
        let n = self.spins.len() as f64;
        if n < 1.0 {
            return Vector3::zero();
        }
        let mut sum = Vector3::zero();
        for s in &self.spins {
            sum = sum + *s;
        }
        sum * (1.0 / n)
    }

    /// Calculate the sublattice magnetizations for kagome lattice
    ///
    /// Returns magnetization vectors for the 3 sublattices.
    ///
    /// # Errors
    ///
    /// Returns error if lattice is not kagome type.
    pub fn kagome_sublattice_magnetizations(&self) -> Result<[Vector3<f64>; 3]> {
        if self.lattice_type != LatticeType::Kagome {
            return Err(Error::InvalidParameter {
                param: "lattice_type".to_string(),
                reason: "sublattice magnetizations only defined for kagome lattice".to_string(),
            });
        }

        let n_cells = self.spins.len() / 3;
        let mut mags = [Vector3::zero(); 3];

        for (i, spin) in self.spins.iter().enumerate() {
            let sub = i % 3;
            mags[sub] = mags[sub] + *spin;
        }

        let inv_n = 1.0 / n_cells as f64;
        for m in &mut mags {
            *m = *m * inv_n;
        }

        Ok(mags)
    }

    /// Initialize spins to the 120-degree Néel order (ground state of triangular AFM)
    ///
    /// The three sublattice directions are separated by 120 degrees in the xy-plane:
    /// - Sublattice A: (1, 0, 0)
    /// - Sublattice B: (-1/2, √3/2, 0)
    /// - Sublattice C: (-1/2, -√3/2, 0)
    pub fn set_120_degree_order(&mut self) {
        let sqrt3_half = 3.0_f64.sqrt() / 2.0;
        let directions = [
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(-0.5, sqrt3_half, 0.0),
            Vector3::new(-0.5, -sqrt3_half, 0.0),
        ];

        match self.lattice_type {
            LatticeType::Triangular => {
                let nx = self.size.0;
                for (i, spin) in self.spins.iter_mut().enumerate() {
                    let ix = i % nx;
                    let iy = i / nx;
                    // Sublattice assignment: (ix + 2*iy) mod 3 gives the correct
                    // 3-coloring for the triangular lattice with neighbor offsets
                    // (+1,0), (-1,0), (0,+1), (0,-1), (+1,-1), (-1,+1)
                    let sub = (ix + 2 * iy) % 3;
                    *spin = directions[sub];
                }
            },
            LatticeType::Kagome => {
                // For kagome, use the 3 sublattice indices directly
                for (i, spin) in self.spins.iter_mut().enumerate() {
                    *spin = directions[i % 3];
                }
            },
            LatticeType::Pyrochlore => {
                // Not directly applicable, but set to planar 120-degree as approximation
                for (i, spin) in self.spins.iter_mut().enumerate() {
                    *spin = directions[i % 3];
                }
            },
        }
    }

    /// Initialize spins to a noncoplanar "umbrella" 120-degree order
    ///
    /// Generalizes [`Self::set_120_degree_order`] by canting all three
    /// 120-degree sublattice directions away from the xy-plane by a common
    /// polar angle `polar_angle` (measured from +z), while keeping their
    /// azimuthal directions 120 degrees apart:
    ///
    /// - Sublattice A: (sin θ, 0, cos θ)
    /// - Sublattice B: (-sin θ / 2, sin θ · √3 / 2, cos θ)
    /// - Sublattice C: (-sin θ / 2, -sin θ · √3 / 2, cos θ)
    ///
    /// At `polar_angle = π/2` this reduces exactly to the coplanar
    /// [`Self::set_120_degree_order`] ground state (zero scalar spin
    /// chirality on every plaquette, since three coplanar vectors have a
    /// vanishing scalar triple product). At `polar_angle = 0` all three
    /// sublattices collapse onto +z (collinear, also zero chirality). For
    /// `0 < polar_angle < π/2` (or `π/2 < polar_angle < π`) the configuration
    /// is genuinely noncoplanar and every elementary triangular plaquette
    /// carries the same nonzero chirality *magnitude* -- this is the
    /// canonical "umbrella state" of a triangular/kagome antiferromagnet
    /// canted by an applied field along the C3 axis.
    ///
    /// # Lattice-dependent net chirality
    ///
    /// Whether the per-plaquette chirality has a *uniform sign* (giving a
    /// nonzero net/total chirality summed over the whole lattice) or a
    /// *staggered sign* (giving exact cancellation between the two
    /// elementary-triangle orientations, like a Néel order canceling in net
    /// magnetization) depends on the lattice geometry:
    ///
    /// - **Kagome**: the "up" (same-unit-cell) and "down"
    ///   (corner-sharing) triangles enclose the three sublattices in the
    ///   *same* counterclockwise cyclic order, so this umbrella state carries
    ///   a uniform-sign chirality and a genuinely nonzero net topological
    ///   Hall response (see `crate::frustrated::transport`).
    /// - **Triangular**: the "up" and "down" triangles enclose the three
    ///   sublattices in *opposite* counterclockwise cyclic order for the
    ///   standard `(ix + 2*iy) % 3` 3-coloring used here, so their
    ///   chiralities are equal in magnitude but opposite in sign, and the net
    ///   chirality summed over the whole periodic lattice is exactly zero
    ///   even though every individual plaquette is noncoplanar. A generic
    ///   (non-uniform / disordered) noncoplanar spin texture on the same
    ///   triangular lattice does *not* have this cancellation.
    ///
    /// # Arguments
    ///
    /// * `polar_angle` - Common canting angle θ from +z \[rad\]
    pub fn set_umbrella_order(&mut self, polar_angle: f64) {
        let sin_theta = polar_angle.sin();
        let cos_theta = polar_angle.cos();
        let sqrt3_half = 3.0_f64.sqrt() / 2.0;
        let directions = [
            Vector3::new(sin_theta, 0.0, cos_theta),
            Vector3::new(-0.5 * sin_theta, sin_theta * sqrt3_half, cos_theta),
            Vector3::new(-0.5 * sin_theta, -sin_theta * sqrt3_half, cos_theta),
        ];

        match self.lattice_type {
            LatticeType::Triangular => {
                let nx = self.size.0;
                for (i, spin) in self.spins.iter_mut().enumerate() {
                    let ix = i % nx;
                    let iy = i / nx;
                    // Same 3-coloring as set_120_degree_order
                    let sub = (ix + 2 * iy) % 3;
                    *spin = directions[sub];
                }
            },
            LatticeType::Kagome => {
                // Kagome sites are laid out sublattice-fastest (3 per cell)
                for (i, spin) in self.spins.iter_mut().enumerate() {
                    *spin = directions[i % 3];
                }
            },
            LatticeType::Pyrochlore => {
                // Not directly applicable, but set to canted 120-degree as approximation
                for (i, spin) in self.spins.iter_mut().enumerate() {
                    *spin = directions[i % 3];
                }
            },
        }
    }

    /// Run Metropolis Monte Carlo simulation
    ///
    /// Performs single-spin-flip Metropolis algorithm at the given temperature.
    /// Returns the final energy per spin.
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature \[K\]
    /// * `n_sweeps` - Number of full lattice sweeps
    /// * `seed` - PRNG seed for reproducibility
    ///
    /// # Errors
    ///
    /// Returns error if temperature is negative or seed is zero.
    pub fn metropolis_mc(&mut self, temperature: f64, n_sweeps: usize, seed: u64) -> Result<f64> {
        if temperature < 0.0 {
            return Err(Error::InvalidParameter {
                param: "temperature".to_string(),
                reason: "temperature must be non-negative".to_string(),
            });
        }

        let mut rng = Xorshift64::new(seed)?;
        let n_sites = self.num_sites();
        let beta = if temperature > 1e-30 {
            1.0 / (crate::constants::KB * temperature)
        } else {
            f64::INFINITY
        };

        for _sweep in 0..n_sweeps {
            for _step in 0..n_sites {
                // Pick a random site
                let site = (rng.next_u64() as usize) % n_sites;

                // Calculate current energy contribution
                let old_energy = self.site_energy(site);
                let old_spin = self.spins[site];

                // Propose a new random spin direction
                let new_spin = rng.random_unit_vector();
                self.spins[site] = new_spin;
                let new_energy = self.site_energy(site);

                let delta_e = new_energy - old_energy;

                // Metropolis acceptance criterion
                if delta_e > 0.0 && !beta.is_infinite() {
                    let acceptance = (-beta * delta_e).exp();
                    if rng.next_f64() >= acceptance {
                        // Reject: restore old spin
                        self.spins[site] = old_spin;
                    }
                }
                // If delta_e <= 0 or T=0 and delta_e < 0, always accept
                // If T=0 and delta_e > 0, always reject
                if delta_e > 0.0 && beta.is_infinite() {
                    self.spins[site] = old_spin;
                }
            }
        }

        Ok(self.total_energy() / n_sites as f64)
    }
}

/// Calculate the frustration parameter f = |θ_CW| / T_N
///
/// The frustration parameter quantifies the degree of geometric frustration.
/// For conventional magnets f ≈ 1, while for strongly frustrated magnets f >> 1.
///
/// # Arguments
///
/// * `curie_weiss_temp` - Curie-Weiss temperature θ_CW \[K\] (typically negative for AFM)
/// * `neel_temp` - Néel ordering temperature T_N \[K\] (must be positive)
///
/// # Returns
///
/// The frustration parameter f = |θ_CW| / T_N
///
/// # Errors
///
/// Returns error if T_N <= 0.
pub fn frustration_parameter(curie_weiss_temp: f64, neel_temp: f64) -> Result<f64> {
    if neel_temp <= 0.0 {
        return Err(Error::InvalidParameter {
            param: "neel_temp".to_string(),
            reason: "Néel temperature must be positive".to_string(),
        });
    }
    Ok(curie_weiss_temp.abs() / neel_temp)
}

/// Calculate the mean-field Curie-Weiss temperature for a frustrated lattice
///
/// θ_CW = -z J S(S+1) / (3 k_B)
///
/// where z is the coordination number and S is the spin quantum number.
///
/// # Arguments
///
/// * `coordination_number` - Number of nearest neighbors z
/// * `coupling_j` - Exchange coupling J \[J\]
/// * `spin_s` - Spin quantum number S
pub fn curie_weiss_temperature(coordination_number: usize, coupling_j: f64, spin_s: f64) -> f64 {
    let z = coordination_number as f64;
    -z * coupling_j * spin_s * (spin_s + 1.0) / (3.0 * crate::constants::KB)
}

/// Validate lattice construction parameters
fn validate_lattice_params(
    nx: usize,
    ny: usize,
    coupling_j: f64,
    lattice_constant: f64,
) -> Result<()> {
    if nx == 0 || ny == 0 {
        return Err(Error::InvalidParameter {
            param: "size".to_string(),
            reason: "lattice dimensions must be nonzero".to_string(),
        });
    }
    if coupling_j <= 0.0 {
        return Err(Error::InvalidParameter {
            param: "coupling_j".to_string(),
            reason: "antiferromagnetic coupling must be positive (J > 0)".to_string(),
        });
    }
    if lattice_constant <= 0.0 {
        return Err(Error::InvalidParameter {
            param: "lattice_constant".to_string(),
            reason: "lattice constant must be positive".to_string(),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_triangular_lattice_creation() {
        let lat = FrustratedLattice::triangular(6, 6, 1e-21, 3e-10)
            .expect("failed to create triangular lattice");
        assert_eq!(lat.num_sites(), 36);
        assert_eq!(lat.lattice_type, LatticeType::Triangular);
    }

    #[test]
    fn test_triangular_neighbor_count() {
        let lat = FrustratedLattice::triangular(8, 8, 1e-21, 3e-10)
            .expect("failed to create triangular lattice");
        // Each site on a triangular lattice should have 6 neighbors
        for (i, nbrs) in lat.neighbors.iter().enumerate() {
            assert_eq!(
                nbrs.len(),
                6,
                "site {} has {} neighbors, expected 6",
                i,
                nbrs.len()
            );
        }
    }

    #[test]
    fn test_kagome_lattice_creation() {
        let lat =
            FrustratedLattice::kagome(4, 4, 1e-21, 5e-10).expect("failed to create kagome lattice");
        assert_eq!(lat.num_sites(), 48); // 3 * 4 * 4
        assert_eq!(lat.lattice_type, LatticeType::Kagome);
    }

    #[test]
    fn test_kagome_neighbor_count() {
        let lat =
            FrustratedLattice::kagome(6, 6, 1e-21, 5e-10).expect("failed to create kagome lattice");
        // Each site on kagome lattice has 4 neighbors
        for (i, nbrs) in lat.neighbors.iter().enumerate() {
            assert_eq!(
                nbrs.len(),
                4,
                "site {} has {} neighbors, expected 4",
                i,
                nbrs.len()
            );
        }
    }

    #[test]
    fn test_pyrochlore_lattice_creation() {
        let lat = FrustratedLattice::pyrochlore(3, 3, 1e-21, 1e-9)
            .expect("failed to create pyrochlore lattice");
        // nz = min(3,3) = 3, so 4 * 3 * 3 * 3 = 108
        assert_eq!(lat.num_sites(), 108);
        assert_eq!(lat.lattice_type, LatticeType::Pyrochlore);
    }

    #[test]
    fn test_frustration_parameter_known_materials() {
        // ZnCu3(OH)6Cl2 (herbertsmithite): θ_CW ≈ -300 K, T_N < 0.05 K => f > 6000
        // But more commonly cited: θ_CW ~ -300 K, no ordering => f ~ ∞
        // SrCr8Ga4O19: θ_CW ≈ -500 K, T_f ≈ 3.5 K => f ≈ 143
        let f =
            frustration_parameter(-500.0, 3.5).expect("failed to compute frustration parameter");
        assert!((f - 142.857).abs() < 0.1, "f = {}, expected ~143", f);

        // Conventional antiferromagnet: f ~ 1
        let f_conv =
            frustration_parameter(-100.0, 90.0).expect("failed to compute frustration parameter");
        assert!(
            f_conv > 1.0 && f_conv < 2.0,
            "f = {}, expected ~1.1",
            f_conv
        );
    }

    #[test]
    fn test_frustration_parameter_invalid() {
        let result = frustration_parameter(-100.0, 0.0);
        assert!(result.is_err());
        let result = frustration_parameter(-100.0, -10.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_120_degree_order_energy() {
        // The 120-degree state is the ground state of the classical triangular AFM
        // Energy per bond = J cos(120°) = -J/2
        // For triangular lattice: 3 bonds per site (6 neighbors, each counted once)
        // E/N = 3 * J * cos(120°) = -3J/2
        let mut lat =
            FrustratedLattice::triangular(6, 6, 1.0, 1e-10).expect("failed to create lattice");
        lat.set_120_degree_order();

        let energy_per_site = lat.total_energy() / lat.num_sites() as f64;
        // cos(120°) = -0.5, each site has 6 neighbors, each bond counted once = 3 bonds/site
        // E/N = 3 * J * (-0.5) = -1.5
        assert!(
            (energy_per_site - (-1.5)).abs() < 0.01,
            "energy per site = {}, expected -1.5",
            energy_per_site
        );
    }

    #[test]
    fn test_umbrella_order_matches_120_degree_at_right_angle() {
        // theta = pi/2 should reduce exactly to the coplanar 120-degree state
        let mut lat_umbrella =
            FrustratedLattice::triangular(6, 6, 1.0, 1e-10).expect("failed to create lattice");
        lat_umbrella.set_umbrella_order(std::f64::consts::FRAC_PI_2);

        let mut lat_planar =
            FrustratedLattice::triangular(6, 6, 1.0, 1e-10).expect("failed to create lattice");
        lat_planar.set_120_degree_order();

        for (s_umbrella, s_planar) in lat_umbrella.spins.iter().zip(lat_planar.spins.iter()) {
            assert!((s_umbrella.x - s_planar.x).abs() < 1e-12);
            assert!((s_umbrella.y - s_planar.y).abs() < 1e-12);
            assert!((s_umbrella.z - s_planar.z).abs() < 1e-12);
        }
    }

    #[test]
    fn test_umbrella_order_is_unit_normalized_and_canted() {
        let mut lat =
            FrustratedLattice::triangular(6, 6, 1.0, 1e-10).expect("failed to create lattice");
        let theta = std::f64::consts::FRAC_PI_3; // 60 degrees: away from the coplanar case
        lat.set_umbrella_order(theta);

        for spin in &lat.spins {
            assert!(
                (spin.magnitude() - 1.0).abs() < 1e-12,
                "umbrella spin {:?} is not unit length",
                spin
            );
            assert!(
                (spin.z - theta.cos()).abs() < 1e-12,
                "z-component should equal cos(theta) for every sublattice"
            );
        }

        // Sanity: with three distinct sublattices present, the configuration
        // must be genuinely noncoplanar (not all spins identical or planar)
        let distinct_xy: std::collections::HashSet<(i64, i64)> = lat
            .spins
            .iter()
            .map(|s| ((s.x * 1e6).round() as i64, (s.y * 1e6).round() as i64))
            .collect();
        assert!(
            distinct_xy.len() >= 3,
            "expected 3 distinct sublattice directions, got {}",
            distinct_xy.len()
        );
    }

    #[test]
    fn test_mc_energy_decreases() {
        let mut lat =
            FrustratedLattice::triangular(6, 6, 1e-21, 3e-10).expect("failed to create lattice");

        // Start from random initial state
        let mut rng = Xorshift64::new(42).expect("failed to create rng");
        for spin in lat.spins.iter_mut() {
            *spin = rng.random_unit_vector();
        }
        let initial_energy = lat.total_energy();

        // Run MC at low temperature - energy should decrease
        let _final_e_per_site = lat
            .metropolis_mc(1.0, 100, 12345)
            .expect("MC simulation failed");
        let final_energy = lat.total_energy();

        assert!(
            final_energy <= initial_energy + 1e-30,
            "energy increased: {} -> {}",
            initial_energy,
            final_energy
        );
    }

    #[test]
    fn test_lattice_size_correctness() {
        // Triangular: N = Nx * Ny
        let lat =
            FrustratedLattice::triangular(5, 7, 1e-21, 1e-10).expect("failed to create lattice");
        assert_eq!(lat.num_sites(), 35);

        // Kagome: N = 3 * Nx * Ny
        let lat = FrustratedLattice::kagome(4, 5, 1e-21, 1e-10).expect("failed to create lattice");
        assert_eq!(lat.num_sites(), 60);

        // Pyrochlore: N = 4 * Nx * Ny * min(Nx,Ny)
        let lat =
            FrustratedLattice::pyrochlore(3, 4, 1e-21, 1e-10).expect("failed to create lattice");
        assert_eq!(lat.num_sites(), 144); // 4 * 3 * 4 * 3
    }

    #[test]
    fn test_xorshift64_deterministic() {
        let mut rng1 = Xorshift64::new(42).expect("failed to create rng");
        let mut rng2 = Xorshift64::new(42).expect("failed to create rng");
        for _ in 0..100 {
            assert_eq!(rng1.next_u64(), rng2.next_u64());
        }
    }

    #[test]
    fn test_xorshift64_zero_seed_error() {
        let result = Xorshift64::new(0);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_lattice_params() {
        assert!(FrustratedLattice::triangular(0, 5, 1e-21, 1e-10).is_err());
        assert!(FrustratedLattice::triangular(5, 5, -1e-21, 1e-10).is_err());
        assert!(FrustratedLattice::triangular(5, 5, 1e-21, 0.0).is_err());
    }

    #[test]
    fn test_curie_weiss_temperature() {
        // For a triangular lattice with z=6, J=1e-21, S=1/2:
        // θ_CW = -6 * 1e-21 * 0.5 * 1.5 / (3 * kB) = -4.5e-21 / (3 * 1.38e-23) ≈ -108.7 K
        let theta = curie_weiss_temperature(6, 1e-21, 0.5);
        assert!(theta < 0.0, "CW temperature should be negative for AFM");
        // Rough check of magnitude
        let expected = -6.0 * 1e-21 * 0.5 * 1.5 / (3.0 * crate::constants::KB);
        assert!(
            (theta - expected).abs() / expected.abs() < 1e-10,
            "theta = {}, expected {}",
            theta,
            expected
        );
    }
}
