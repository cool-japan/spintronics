//! Geometric frustration -> topological Hall transport
//!
//! This module wires the scalar spin chirality of a [`FrustratedLattice`]
//! configuration into the emergent-field topological Hall machinery of
//! [`crate::effect::topological_hall`], giving a working
//! chirality -> emergent-field -> Hall-response pipeline for frustrated
//! (triangular/kagome/pyrochlore) spin textures.
//!
//! # Physics Background
//!
//! ## Scalar spin chirality
//!
//! For three spins Sᵢ, Sⱼ, Sₖ on an elementary triangular plaquette, the
//! scalar spin chirality is the signed volume of the parallelepiped they
//! span:
//!
//! χ_ijk = Sᵢ · (Sⱼ × Sₖ)
//!
//! χ vanishes identically for any coplanar spin configuration (in particular
//! the coplanar 120-degree Néel order of the triangular/kagome
//! antiferromagnet, since any three coplanar vectors have a zero scalar
//! triple product), and is nonzero precisely when the local spin texture is
//! noncoplanar -- e.g. a canted "umbrella" 120-degree order
//! ([`FrustratedLattice::set_umbrella_order`]), an all-in/all-out pyrochlore
//! configuration, or a skyrmion/hedgehog texture.
//!
//! ## Plaquette definition
//!
//! A "plaquette" is defined purely combinatorially from the lattice's
//! neighbor adjacency, as an elementary triangle: three mutually-adjacent
//! sites `(i, j, k)`. This is lattice-agnostic and recovers the physically
//! expected triangular motifs for all three `FrustratedLattice` geometries:
//!
//! - **Triangular**: both "up" and "down" elementary triangles of the
//!   triangular net (2 per site on a periodic Nx*Ny lattice).
//! - **Kagome**: both the "up" (same-unit-cell) and "down"
//!   (corner-sharing-between-cells) triangles (2 per unit cell).
//! - **Pyrochlore**: the 4 triangular faces of each corner-sharing
//!   tetrahedron (tetrahedra share only a single vertex, so no additional
//!   triangles are introduced by the inter-tetrahedron bonds).
//!
//! Each plaquette is oriented counterclockwise (as seen from +z, using the
//! in-plane projection of its site positions) so that the chiralities of a
//! coherent (non-random) texture combine constructively rather than
//! cancelling due to an arbitrary index-labelling accident. For the
//! triangular and kagome lattices (both strictly 2D, z = 0 everywhere) this
//! orientation is exact and unambiguous; for the fully 3D pyrochlore lattice
//! the in-plane (x, y) projection used for orientation can degenerate for
//! near-vertical faces, so the *global* sign bookkeeping across the whole
//! pyrochlore lattice should be treated as a documented convention rather
//! than an absolute physical truth (the magnitude of each plaquette's
//! chirality is unaffected either way).
//!
//! ## Net chirality depends on the lattice geometry
//!
//! A spatially *uniform* noncoplanar spin pattern (e.g.
//! [`FrustratedLattice::set_umbrella_order`]) can have every individual
//! plaquette noncoplanar while still summing to zero *net* chirality over
//! the whole periodic lattice -- exactly as a Néel order has nonzero local
//! (staggered) moments but zero net magnetization. Whether this happens
//! depends on whether the two elementary-triangle orientations ("up" and
//! "down") enclose the three sublattices in the same or opposite CCW cyclic
//! order:
//!
//! - **Kagome**: same cyclic order for both orientations -> uniform-sign
//!   chirality -> nonzero net response for the uniform umbrella state.
//! - **Triangular**: opposite cyclic order for the two orientations (under
//!   the standard `(ix + 2*iy) % 3` 3-coloring) -> equal-and-opposite
//!   chirality -> exact net cancellation for the uniform umbrella state,
//!   even though every plaquette is individually noncoplanar. A generic
//!   (disordered) noncoplanar texture on the same triangular lattice does
//!   not have this cancellation and gives a nonzero net response.
//!
//! See the module tests for a worked demonstration of both cases.
//!
//! ## Chirality -> emergent field
//!
//! Following the standard Berg-Lüscher lattice discretization of the
//! topological charge, the chirality together with the pairwise spin dot
//! products determines the solid angle Ω the three spins subtend on the
//! unit sphere (see [`berg_luscher_solid_angle`]). This solid angle is the
//! Berry phase conduction electrons accumulate hopping around the plaquette,
//! which acts on them exactly like an Aharonov-Bohm phase from a real
//! magnetic flux threading the plaquette:
//!
//! B_emergent = Φ₀ Ω / (4π A)
//!
//! where A is the plaquette's real-space area (see
//! [`TopologicalHall::emergent_field_from_solid_angle`]). Summing Ω = 4π over
//! the area of a single skyrmion (Q=1) recovers exactly the "one flux
//! quantum per skyrmion" convention already used throughout
//! `topological_hall.rs`.
//!
//! ## Emergent field -> Hall response
//!
//! The area-averaged emergent field over the whole lattice defines an
//! effective signed topological charge density
//! ([`FrustratedTransport::effective_topological_charge_density`]), which
//! drives the topological Hall resistivity via
//! [`TopologicalHall::hall_resistivity_from_charge_density`] -- the same
//! R₀·n·Q law used for discrete skyrmion lattices, here fed by a
//! continuous/discretized chirality field instead.
//!
//! # References
//!
//! - Y. Taguchi et al., "Spin Chirality, Berry Phase, and Anomalous Hall
//!   Effect in a Frustrated Ferromagnet", Science 291, 2573 (2001)
//! - S. Nakatsuji, N. Kiyohara & T. Higo, "Large anomalous Hall effect in a
//!   non-collinear antiferromagnet at room temperature", Nature 527, 212
//!   (2015)
//! - H. Berg & M. Lüscher, "Definition and statistical distributions of a
//!   topological number in the lattice O(3) sigma-model", Nucl. Phys. B 190,
//!   412 (1981)
//! - A. Van Oosterom & J. Strackee, "The Solid Angle of a Plane Triangle",
//!   IEEE Trans. Biomed. Eng. BME-30, 125 (1983)

use std::f64::consts::PI;

use crate::effect::topological_hall::{TopologicalHall, PHI_0};
use crate::error::{Error, Result};
use crate::frustrated::lattice::{FrustratedLattice, LatticeType};
use crate::vector3::Vector3;

/// Below this scalar-chirality magnitude, the solid-angle formula's `atan2`
/// branch is snapped to zero rather than evaluated (see
/// [`berg_luscher_solid_angle`] for why this regularization is necessary).
const SOLID_ANGLE_CHIRALITY_EPSILON: f64 = 1e-12;

/// Minimum real-space plaquette area \[m²\] below which a triangle is
/// considered degenerate (near-colinear positions) and rejected with an
/// error rather than silently dividing by a near-zero area.
const MIN_PLAQUETTE_AREA: f64 = 1e-30;

/// Scalar spin chirality χ = Sᵢ·(Sⱼ×Sₖ) of an ordered spin triple.
///
/// This is the discrete lattice analogue of the continuum topological
/// charge density; it vanishes identically whenever the three spins are
/// coplanar (in particular for any collinear or planar Néel/120-degree
/// order), and is odd under any orientation-reversing transformation applied
/// to every spin (e.g. mirroring a single shared Cartesian component of all
/// three, which negates chi while leaving the pairwise dot products
/// unchanged).
///
/// # References
/// - Y. Taguchi et al., Science 291, 2573 (2001)
pub fn scalar_spin_chirality(s_i: Vector3<f64>, s_j: Vector3<f64>, s_k: Vector3<f64>) -> f64 {
    s_i.dot(&s_j.cross(&s_k))
}

/// Signed solid angle Ω \[sr\] subtended by three (approximately unit-length)
/// spin directions, via the Van Oosterom-Strackee / Berg-Lüscher triangle
/// formula:
///
/// tan(Ω/2) = χ / (1 + Sᵢ·Sⱼ + Sⱼ·Sₖ + Sₖ·Sᵢ)
///
/// evaluated with `atan2` for the correct branch. Summing Ω over every
/// plaquette of a fine discretization of a skyrmion recovers the standard
/// lattice topological charge Q = (1/4π) ΣΩ.
///
/// # Numerical note
///
/// The `atan2` branch has a removable-singularity ambiguity exactly on the
/// measure-zero locus where the three spins are coplanar (`chi == 0`) *and*
/// the denominator `1 + dot_ij + dot_jk + dot_ki` is negative -- this occurs,
/// e.g., exactly for the coplanar 120-degree triangular Néel order, where the
/// denominator equals -1/2. A vanishingly small floating-point roundoff in
/// `chi` can then round `atan2` to either +π or -π, producing a spurious
/// solid angle near ±2π instead of the physically correct zero (a coplanar
/// plaquette carries no Berry curvature and should contribute nothing).
/// This is regularized by snapping the result to zero whenever `chi` itself
/// is negligible; genuinely large chiralities are unaffected by the
/// regularization regardless of the denominator's sign.
///
/// # Arguments
/// * `chi` - Scalar spin chirality (see [`scalar_spin_chirality`])
/// * `dot_ij`, `dot_jk`, `dot_ki` - Pairwise dot products of the three spins
///
/// # References
/// - H. Berg & M. Lüscher, Nucl. Phys. B 190, 412 (1981)
/// - A. Van Oosterom & J. Strackee, IEEE Trans. Biomed. Eng. BME-30, 125 (1983)
pub fn berg_luscher_solid_angle(chi: f64, dot_ij: f64, dot_jk: f64, dot_ki: f64) -> f64 {
    if chi.abs() < SOLID_ANGLE_CHIRALITY_EPSILON {
        return 0.0;
    }
    2.0 * chi.atan2(1.0 + dot_ij + dot_jk + dot_ki)
}

/// One elementary triangular plaquette of a [`FrustratedLattice`], together
/// with its geometric and chirality-transport data.
#[derive(Debug, Clone, Copy)]
pub struct Plaquette {
    /// Site indices `(i, j, k)` of the plaquette's three vertices, oriented
    /// counterclockwise as seen from +z (using the in-plane projection of
    /// the site positions).
    pub sites: (usize, usize, usize),
    /// Real-space area of the triangle \[m²\].
    pub area: f64,
    /// Scalar spin chirality χ = Sᵢ·(Sⱼ×Sₖ) (dimensionless; magnitude at
    /// most 1 for unit-length spins).
    pub chirality: f64,
    /// Signed solid angle Ω \[sr\] subtended by the three spin directions
    /// (see [`berg_luscher_solid_angle`]).
    pub solid_angle: f64,
    /// Emergent magnetic field \[T\] implied by this plaquette's solid angle
    /// over its real-space area (see
    /// [`TopologicalHall::emergent_field_from_solid_angle`]).
    pub emergent_field: f64,
}

/// Per-plaquette and aggregate topological-Hall transport analysis for a
/// [`FrustratedLattice`] spin configuration.
///
/// Wires the scalar spin chirality of a frustrated-lattice configuration
/// into the emergent-field topological Hall machinery of
/// [`crate::effect::topological_hall`]: for each elementary triangular
/// plaquette it computes the scalar spin chirality, converts it to a solid
/// angle and an emergent magnetic field, and aggregates these into an
/// effective topological charge density that can drive
/// [`TopologicalHall::hall_resistivity_from_charge_density`].
///
/// # Example
///
/// The kagome lattice's "up" and "down" elementary triangles reinforce
/// (rather than cancel) under a uniform umbrella spin texture, giving a
/// genuinely nonzero net chirality and Hall response (see the module
/// documentation for why this differs from the triangular lattice):
///
/// ```
/// use spintronics::frustrated::lattice::FrustratedLattice;
/// use spintronics::frustrated::transport::FrustratedTransport;
/// use spintronics::effect::topological_hall::TopologicalHall;
///
/// let mut lattice = FrustratedLattice::kagome(4, 4, 1.0, 3e-10)
///     .expect("failed to create lattice");
/// lattice.set_umbrella_order(std::f64::consts::FRAC_PI_3);
///
/// let transport = FrustratedTransport::from_lattice(&lattice)
///     .expect("failed to build transport analysis");
/// assert!(transport.num_plaquettes() > 0);
/// assert!(transport.total_chirality().abs() > 0.0);
///
/// let mnsi = TopologicalHall::mnsi();
/// let rho_the = transport.hall_resistivity(&mnsi);
/// assert!(rho_the.is_finite());
/// assert!(rho_the.abs() > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct FrustratedTransport {
    /// Per-plaquette chirality/transport data, one entry per elementary
    /// triangle discovered in the lattice's neighbor adjacency.
    pub plaquettes: Vec<Plaquette>,
}

impl FrustratedTransport {
    /// Build the plaquette chirality/transport decomposition of a
    /// [`FrustratedLattice`] configuration.
    ///
    /// # Errors
    /// Returns an error if the lattice's `spins`/`positions`/`neighbors`
    /// arrays have mismatched lengths, if any neighbor index is out of
    /// bounds, or if a discovered plaquette has a degenerate (near-zero)
    /// real-space area (e.g. colinear positions).
    pub fn from_lattice(lattice: &FrustratedLattice) -> Result<Self> {
        validate_lattice_arrays(lattice)?;

        let raw_triangles = enumerate_plaquettes(&lattice.neighbors);
        let lattice_vectors = periodicity_vectors(lattice);
        let mut plaquettes = Vec::with_capacity(raw_triangles.len());

        for (i0, j0, k0) in raw_triangles {
            // Neighbor bonds that wrap across a periodic boundary have a
            // *naive* Cartesian position difference pointing across the
            // entire simulation cell rather than to the nearest physical
            // image; resolve this with the minimum-image convention before
            // using the displacement for orientation or area.
            let e_i_j = minimum_image_displacement(&lattice.positions, &lattice_vectors, i0, j0);
            let e_i_k = minimum_image_displacement(&lattice.positions, &lattice_vectors, i0, k0);
            let (i, j, k, e1, e2) = orient_ccw(i0, j0, k0, e_i_j, e_i_k);

            let area = 0.5 * e1.cross(&e2).magnitude();
            if !area.is_finite() || area <= MIN_PLAQUETTE_AREA {
                return Err(Error::NumericalError {
                    description: format!(
                        "degenerate plaquette ({}, {}, {}) has near-zero area {:.3e} m^2; \
                         positions may be colinear",
                        i, j, k, area
                    ),
                });
            }

            let s_i = lattice.spins[i];
            let s_j = lattice.spins[j];
            let s_k = lattice.spins[k];

            let chirality = scalar_spin_chirality(s_i, s_j, s_k);
            let dot_ij = s_i.dot(&s_j);
            let dot_jk = s_j.dot(&s_k);
            let dot_ki = s_k.dot(&s_i);
            let solid_angle = berg_luscher_solid_angle(chirality, dot_ij, dot_jk, dot_ki);
            let emergent_field =
                TopologicalHall::emergent_field_from_solid_angle(solid_angle, area)?;

            plaquettes.push(Plaquette {
                sites: (i, j, k),
                area,
                chirality,
                solid_angle,
                emergent_field,
            });
        }

        Ok(Self { plaquettes })
    }

    /// Number of elementary triangular plaquettes found.
    pub fn num_plaquettes(&self) -> usize {
        self.plaquettes.len()
    }

    /// Sum of the scalar spin chirality over all plaquettes.
    pub fn total_chirality(&self) -> f64 {
        self.plaquettes.iter().map(|p| p.chirality).sum()
    }

    /// Mean scalar spin chirality per plaquette (0 if there are none).
    pub fn mean_chirality(&self) -> f64 {
        let n = self.plaquettes.len();
        if n == 0 {
            0.0
        } else {
            self.total_chirality() / n as f64
        }
    }

    /// Sum of the signed solid angle over all plaquettes \[sr\].
    pub fn total_solid_angle(&self) -> f64 {
        self.plaquettes.iter().map(|p| p.solid_angle).sum()
    }

    /// Sum of the real-space area over all plaquettes \[m²\].
    pub fn total_area(&self) -> f64 {
        self.plaquettes.iter().map(|p| p.area).sum()
    }

    /// Effective signed topological charge density \[m⁻²\] of the whole
    /// plaquette decomposition, Q_eff = (ΣΩ) / (4π · ΣA).
    ///
    /// This is the discretized-lattice analogue of a skyrmion density
    /// weighted by topological charge, obtained by accumulating the
    /// Berg-Lüscher solid angle over every plaquette instead of assuming a
    /// discrete skyrmion count -- consistent with the fact that a single
    /// full skyrmion (Q=1) integrates to ΣΩ = 4π over its area.
    ///
    /// Returns `0.0` if the decomposition has zero total area (e.g. no
    /// plaquettes).
    pub fn effective_topological_charge_density(&self) -> f64 {
        let area = self.total_area();
        if area > 0.0 {
            self.total_solid_angle() / (4.0 * PI * area)
        } else {
            0.0
        }
    }

    /// Area-weighted mean emergent magnetic field \[T\] over the whole
    /// plaquette decomposition, Φ₀ · Q_eff.
    pub fn mean_emergent_field(&self) -> f64 {
        PHI_0 * self.effective_topological_charge_density()
    }

    /// Topological Hall resistivity \[Ω·cm\] this plaquette decomposition
    /// would drive in a given material, via
    /// [`TopologicalHall::hall_resistivity_from_charge_density`].
    pub fn hall_resistivity(&self, effect: &TopologicalHall) -> f64 {
        effect.hall_resistivity_from_charge_density(self.effective_topological_charge_density())
    }
}

/// Compute the topological Hall resistivity \[Ω·cm\] driven by the scalar
/// spin chirality of a [`FrustratedLattice`] configuration, in a given
/// [`TopologicalHall`] material.
///
/// Convenience wrapper composing [`FrustratedTransport::from_lattice`] and
/// [`FrustratedTransport::hall_resistivity`] for the common case where only
/// the aggregate resistivity (not the full per-plaquette breakdown) is
/// needed.
///
/// # Errors
/// Propagates any error from [`FrustratedTransport::from_lattice`].
pub fn frustration_hall_response(
    lattice: &FrustratedLattice,
    effect: &TopologicalHall,
) -> Result<f64> {
    Ok(FrustratedTransport::from_lattice(lattice)?.hall_resistivity(effect))
}

/// Enumerate elementary triangular plaquettes from a lattice's neighbor
/// adjacency.
///
/// A plaquette is a triple of mutually-adjacent sites `(i, j, k)` with
/// `i < j < k` -- i.e. an elementary triangle of the neighbor graph. This
/// definition is lattice-agnostic; see the module documentation for the
/// physical triangular motifs it recovers for each `FrustratedLattice`
/// geometry. Runs in O(N * d^2) for coordination number `d`.
fn enumerate_plaquettes(neighbors: &[Vec<usize>]) -> Vec<(usize, usize, usize)> {
    let mut triangles = Vec::new();
    for (i, nbrs_i) in neighbors.iter().enumerate() {
        let mut higher: Vec<usize> = nbrs_i.iter().copied().filter(|&x| x > i).collect();
        higher.sort_unstable();
        higher.dedup();

        for (a, &j) in higher.iter().enumerate() {
            for &k in &higher[(a + 1)..] {
                if neighbors[j].contains(&k) {
                    triangles.push((i, j, k));
                }
            }
        }
    }
    triangles
}

/// Reorder a triangle's vertices to be counterclockwise as seen from +z,
/// using the in-plane (x, y) projection of its two edge-displacement vectors
/// `e_i_j = pos(j) - pos(i)` and `e_i_k = pos(k) - pos(i)` (already resolved
/// to the minimum-image convention by the caller, see
/// [`minimum_image_displacement`]). Returns the reordered vertex indices
/// together with the corresponding (possibly swapped) edge-displacement
/// vectors `(e1, e2)` from the *returned* first vertex, so callers never need
/// to re-derive positions after reordering. See the module documentation for
/// why CCW orientation matters and its limitations for fully 3D lattices.
fn orient_ccw(
    i: usize,
    j: usize,
    k: usize,
    e_i_j: Vector3<f64>,
    e_i_k: Vector3<f64>,
) -> (usize, usize, usize, Vector3<f64>, Vector3<f64>) {
    let normal_z = e_i_j.x * e_i_k.y - e_i_j.y * e_i_k.x;
    if normal_z >= 0.0 {
        (i, j, k, e_i_j, e_i_k)
    } else {
        (i, k, j, e_i_k, e_i_j)
    }
}

/// Bravais (periodicity) vectors of a [`FrustratedLattice`]'s periodic
/// supercell, used to resolve neighbor bonds that wrap across a periodic
/// boundary to their minimum-image displacement.
///
/// The `positions` array stores each site's *unwrapped* Cartesian position
/// within a single copy of the simulation cell; for a bond whose neighbor
/// relationship crosses a periodic boundary (e.g. `ix = nx-1` to `ix = 0`),
/// the raw difference of stored positions points across the *entire*
/// simulation cell instead of to the nearest physical image. Adding/
/// subtracting integer combinations of these vectors recovers the correct
/// short displacement (see [`minimum_image_displacement`]).
fn periodicity_vectors(lattice: &FrustratedLattice) -> Vec<Vector3<f64>> {
    let a = lattice.lattice_constant;
    let (nx, ny) = lattice.size;
    let sqrt3_half = 3.0_f64.sqrt() / 2.0;
    match lattice.lattice_type {
        LatticeType::Triangular | LatticeType::Kagome => {
            // Both geometries repeat on the same underlying triangular
            // Bravais lattice: a1 = a*(1,0), a2 = a*(1/2, sqrt(3)/2).
            vec![
                Vector3::new(nx as f64 * a, 0.0, 0.0),
                Vector3::new(ny as f64 * 0.5 * a, ny as f64 * sqrt3_half * a, 0.0),
            ]
        },
        LatticeType::Pyrochlore => {
            // FrustratedLattice::pyrochlore uses nz = min(nx, ny) cubic cells.
            let nz = nx.min(ny);
            vec![
                Vector3::new(nx as f64 * a, 0.0, 0.0),
                Vector3::new(0.0, ny as f64 * a, 0.0),
                Vector3::new(0.0, 0.0, nz as f64 * a),
            ]
        },
    }
}

/// Minimum-image displacement `positions[j] - positions[i]`, correctly
/// handling neighbor bonds that wrap across a periodic boundary by choosing
/// whichever integer combination of `lattice_vectors` (coefficients in
/// `{-1, 0, 1}`) minimizes the displacement's magnitude. See
/// [`periodicity_vectors`] for why the naive difference can otherwise be
/// wrong.
fn minimum_image_displacement(
    positions: &[Vector3<f64>],
    lattice_vectors: &[Vector3<f64>],
    i: usize,
    j: usize,
) -> Vector3<f64> {
    let naive = positions[j] - positions[i];
    let coeffs: [f64; 3] = [-1.0, 0.0, 1.0];

    let mut best = naive;
    let mut best_mag_sq = naive.magnitude_squared();

    let mut consider = |candidate: Vector3<f64>| {
        let mag_sq = candidate.magnitude_squared();
        if mag_sq < best_mag_sq {
            best = candidate;
            best_mag_sq = mag_sq;
        }
    };

    match lattice_vectors {
        [l1, l2] => {
            for &n1 in &coeffs {
                for &n2 in &coeffs {
                    consider(naive + *l1 * n1 + *l2 * n2);
                }
            }
        },
        [l1, l2, l3] => {
            for &n1 in &coeffs {
                for &n2 in &coeffs {
                    for &n3 in &coeffs {
                        consider(naive + *l1 * n1 + *l2 * n2 + *l3 * n3);
                    }
                }
            }
        },
        _ => {},
    }

    best
}

/// Validate that a lattice's parallel arrays are mutually consistent before
/// indexing into them, so that a malformed [`FrustratedLattice`] (e.g.
/// hand-constructed with mismatched field lengths) produces a descriptive
/// error instead of an out-of-bounds panic.
fn validate_lattice_arrays(lattice: &FrustratedLattice) -> Result<()> {
    let n = lattice.spins.len();
    if lattice.positions.len() != n || lattice.neighbors.len() != n {
        return Err(Error::DimensionMismatch {
            expected: format!("spins/positions/neighbors all of length {}", n),
            actual: format!(
                "spins={}, positions={}, neighbors={}",
                n,
                lattice.positions.len(),
                lattice.neighbors.len()
            ),
        });
    }
    for (site, nbrs) in lattice.neighbors.iter().enumerate() {
        for &j in nbrs {
            if j >= n {
                return Err(Error::DimensionMismatch {
                    expected: format!("neighbor index < {}", n),
                    actual: format!("site {} lists neighbor index {}", site, j),
                });
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frustrated::lattice::Xorshift64;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_3};

    #[test]
    fn test_plaquette_count_triangular() {
        let lat = FrustratedLattice::triangular(6, 6, 1.0, 1e-10)
            .expect("failed to create triangular lattice");
        let transport = FrustratedTransport::from_lattice(&lat).expect("transport failed");
        // Periodic triangular lattice: 2 elementary triangles (1 up + 1 down) per site
        assert_eq!(transport.num_plaquettes(), 2 * lat.num_sites());
    }

    #[test]
    fn test_plaquette_count_kagome() {
        let lat =
            FrustratedLattice::kagome(4, 4, 1.0, 1e-10).expect("failed to create kagome lattice");
        let transport = FrustratedTransport::from_lattice(&lat).expect("transport failed");
        // Periodic kagome lattice: 2 elementary triangles (1 up + 1 down) per unit cell
        let n_cells = 4 * 4;
        assert_eq!(transport.num_plaquettes(), 2 * n_cells);
    }

    #[test]
    fn test_coplanar_120_degree_gives_zero_chirality_and_hall_response() {
        let mut lat = FrustratedLattice::triangular(6, 6, 1.0, 1e-10)
            .expect("failed to create triangular lattice");
        lat.set_120_degree_order();

        let transport = FrustratedTransport::from_lattice(&lat).expect("transport failed");
        assert!(transport.num_plaquettes() > 0);

        for p in &transport.plaquettes {
            assert!(
                p.chirality.abs() < 1e-9,
                "coplanar plaquette {:?} has nonzero chirality {}",
                p.sites,
                p.chirality
            );
            assert!(
                p.solid_angle.abs() < 1e-9,
                "coplanar plaquette {:?} has nonzero solid angle {}",
                p.sites,
                p.solid_angle
            );
        }

        assert!(transport.total_chirality().abs() < 1e-9);
        assert!(transport.total_solid_angle().abs() < 1e-9);

        let mnsi = TopologicalHall::mnsi();
        let rho = transport.hall_resistivity(&mnsi);
        assert!(
            rho.abs() < 1e-9,
            "expected ~zero Hall response, got {}",
            rho
        );
    }

    #[test]
    fn test_triangular_umbrella_per_plaquette_nonzero_but_net_cancels() {
        // On the triangular lattice, the "up" and "down" elementary
        // triangles enclose the three sublattices in *opposite* CCW cyclic
        // order under the standard (ix + 2*iy) % 3 coloring used by
        // set_umbrella_order (see that method's documentation), so a
        // uniform umbrella state has equal-and-opposite chirality on the two
        // triangle orientations: every plaquette is individually
        // noncoplanar, but the net/total over the whole periodic lattice
        // cancels exactly -- analogous to a Neel order having nonzero
        // (staggered) local moments but zero net magnetization.
        let mut lat = FrustratedLattice::triangular(6, 6, 1.0, 1e-10)
            .expect("failed to create triangular lattice");
        lat.set_umbrella_order(FRAC_PI_3); // 60 degrees, away from the coplanar case

        let transport = FrustratedTransport::from_lattice(&lat).expect("transport failed");
        assert!(transport.num_plaquettes() > 0);

        // Analytic value for an elementary triangle with sublattices canted
        // by theta = pi/3: |chi| = (3*sqrt(3)/2) * sin(theta)^2 * cos(theta)
        let theta = FRAC_PI_3;
        let expected_magnitude = 1.5 * 3.0_f64.sqrt() * theta.sin().powi(2) * theta.cos();
        assert!(expected_magnitude > 0.5, "sanity check on analytic value");

        let mut n_positive = 0usize;
        let mut n_negative = 0usize;
        for p in &transport.plaquettes {
            assert!(
                (p.chirality.abs() - expected_magnitude).abs() < 1e-6,
                "plaquette {:?} chirality {} does not match analytic magnitude {}",
                p.sites,
                p.chirality,
                expected_magnitude
            );
            if p.chirality > 0.0 {
                n_positive += 1;
            } else {
                n_negative += 1;
            }
        }

        // Every plaquette is individually noncoplanar and splits evenly into
        // the two (equal-magnitude, opposite-sign) triangle orientations...
        assert_eq!(n_positive + n_negative, transport.num_plaquettes());
        assert!(
            n_positive > 0 && n_negative > 0,
            "expected both triangle orientations to be present"
        );
        assert_eq!(
            n_positive, n_negative,
            "up/down triangles should be equinumerous"
        );

        // ...but the net chirality and Hall response cancel exactly.
        assert!(transport.total_chirality().abs() < 1e-9);
        let mnsi = TopologicalHall::mnsi();
        assert!(transport.hall_resistivity(&mnsi).abs() < 1e-9);
    }

    #[test]
    fn test_triangular_random_noncoplanar_gives_nonzero_net_response() {
        // A generic (non-symmetry-fine-tuned) noncoplanar configuration on
        // the same triangular lattice does not have the up/down cancellation
        // of the uniform umbrella state above, and generically gives a
        // nonzero net chirality and Hall response.
        let mut lat = FrustratedLattice::triangular(6, 6, 1.0, 1e-10)
            .expect("failed to create triangular lattice");
        let mut rng = Xorshift64::new(2024).expect("failed to create rng");
        for spin in lat.spins.iter_mut() {
            *spin = rng.random_unit_vector();
        }

        let transport = FrustratedTransport::from_lattice(&lat).expect("transport failed");
        assert!(transport.total_chirality().abs() > 0.5);

        let mnsi = TopologicalHall::mnsi();
        assert!(transport.hall_resistivity(&mnsi).abs() > 0.0);
    }

    #[test]
    fn test_kagome_umbrella_gives_nonzero_uniform_sign_chirality_and_hall_response() {
        // Unlike the triangular lattice, kagome's "up" (same-unit-cell) and
        // "down" (corner-sharing) triangles enclose the three sublattices in
        // the *same* CCW cyclic order, so the uniform umbrella state carries
        // a uniform-sign chirality and a genuinely nonzero net topological
        // Hall response.
        let mut lat =
            FrustratedLattice::kagome(4, 4, 1.0, 1e-10).expect("failed to create kagome lattice");
        let theta = FRAC_PI_3;
        lat.set_umbrella_order(theta);

        let transport = FrustratedTransport::from_lattice(&lat).expect("transport failed");
        assert!(transport.num_plaquettes() > 0);

        let expected_magnitude = 1.5 * 3.0_f64.sqrt() * theta.sin().powi(2) * theta.cos();
        let first_sign = transport.plaquettes[0].chirality.signum();
        for p in &transport.plaquettes {
            assert!(
                (p.chirality.abs() - expected_magnitude).abs() < 1e-6,
                "plaquette {:?} chirality {} does not match analytic magnitude {}",
                p.sites,
                p.chirality,
                expected_magnitude
            );
            assert_eq!(
                p.chirality.signum(),
                first_sign,
                "plaquette {:?} chirality sign is not coherent with the rest of the lattice",
                p.sites
            );
        }

        assert!(transport.total_chirality().abs() > 1.0);
        assert!(transport.total_solid_angle().abs() > 1.0);
        assert!(transport.mean_emergent_field().is_finite());

        let mnsi = TopologicalHall::mnsi();
        let rho = transport.hall_resistivity(&mnsi);
        assert!(rho.is_finite());
        assert!(
            rho.abs() > 0.0,
            "expected nonzero Hall response for kagome umbrella order"
        );
    }

    #[test]
    fn test_sign_reversal_under_global_mirror() {
        // Uses the kagome lattice, whose uniform umbrella state has a
        // genuinely nonzero net Hall response (see
        // test_kagome_umbrella_gives_nonzero_uniform_sign_chirality_and_hall_response;
        // the triangular lattice's uniform umbrella state cancels in the
        // net, so it is not a suitable baseline for a sign-reversal check).
        let mut lat =
            FrustratedLattice::kagome(4, 4, 1.0, 1e-10).expect("failed to create kagome lattice");
        lat.set_umbrella_order(FRAC_PI_3);

        let mnsi = TopologicalHall::mnsi();
        let transport_before = FrustratedTransport::from_lattice(&lat).expect("transport failed");
        let rho_before = transport_before.hall_resistivity(&mnsi);
        assert!(rho_before.abs() > 0.0);

        // Mirror the z-component of every spin: reverses the global chirality
        for spin in lat.spins.iter_mut() {
            spin.z = -spin.z;
        }

        let transport_after = FrustratedTransport::from_lattice(&lat).expect("transport failed");
        let rho_after = transport_after.hall_resistivity(&mnsi);

        assert!(rho_after.abs() > 0.0);
        assert!(
            (rho_before + rho_after).abs() < 1e-9 * rho_before.abs().max(1e-9),
            "mirroring should reverse the sign of the Hall response: before={:.6e}, after={:.6e}",
            rho_before,
            rho_after
        );
    }

    #[test]
    fn test_frustration_hall_response_matches_struct_method() {
        let mut lat = FrustratedLattice::triangular(5, 5, 1.0, 1e-10)
            .expect("failed to create triangular lattice");
        lat.set_umbrella_order(FRAC_PI_3);

        let mnsi = TopologicalHall::mnsi();
        let via_free_fn = frustration_hall_response(&lat, &mnsi).expect("free fn failed");
        let via_struct = FrustratedTransport::from_lattice(&lat)
            .expect("transport failed")
            .hall_resistivity(&mnsi);

        assert!((via_free_fn - via_struct).abs() < 1e-30);
    }

    #[test]
    fn test_umbrella_at_right_angle_matches_coplanar_zero() {
        // theta = pi/2 reduces set_umbrella_order to the coplanar case, so the
        // Hall response must also vanish there (consistency between the two
        // spin-configuration entry points).
        let mut lat = FrustratedLattice::triangular(6, 6, 1.0, 1e-10)
            .expect("failed to create triangular lattice");
        lat.set_umbrella_order(FRAC_PI_2);

        let transport = FrustratedTransport::from_lattice(&lat).expect("transport failed");
        assert!(transport.total_chirality().abs() < 1e-9);

        let mnsi = TopologicalHall::mnsi();
        assert!(transport.hall_resistivity(&mnsi).abs() < 1e-9);
    }

    #[test]
    fn test_chirality_and_solid_angle_are_bounded_for_relaxed_configuration() {
        let mut lat = FrustratedLattice::triangular(6, 6, 1.0, 1e-10)
            .expect("failed to create triangular lattice");
        let mut rng = Xorshift64::new(7).expect("failed to create rng");
        for spin in lat.spins.iter_mut() {
            *spin = rng.random_unit_vector();
        }
        let _ = lat
            .metropolis_mc(0.5, 20, 99)
            .expect("MC simulation failed");

        let transport = FrustratedTransport::from_lattice(&lat).expect("transport failed");
        for p in &transport.plaquettes {
            assert!(
                p.chirality.abs() <= 1.0 + 1e-9,
                "chirality {} out of the physical [-1,1] bound",
                p.chirality
            );
            assert!(
                p.solid_angle.abs() <= 2.0 * PI + 1e-6,
                "solid angle {} out of the physical [-2*pi,2*pi] bound",
                p.solid_angle
            );
            assert!(p.emergent_field.is_finite());
        }
    }

    #[test]
    fn test_from_lattice_rejects_mismatched_arrays() {
        let mut lat = FrustratedLattice::triangular(3, 3, 1.0, 1e-10)
            .expect("failed to create triangular lattice");
        lat.positions.pop(); // desynchronize positions from spins/neighbors

        assert!(FrustratedTransport::from_lattice(&lat).is_err());
    }

    #[test]
    fn test_from_lattice_rejects_out_of_bounds_neighbor_index() {
        let mut lat = FrustratedLattice::triangular(3, 3, 1.0, 1e-10)
            .expect("failed to create triangular lattice");
        let n = lat.num_sites();
        lat.neighbors[0].push(n + 100); // out-of-bounds neighbor index

        assert!(FrustratedTransport::from_lattice(&lat).is_err());
    }
}
