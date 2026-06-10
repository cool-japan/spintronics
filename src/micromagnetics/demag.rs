//! Demagnetization tensor for finite-difference micromagnetics.
//!
//! ## Physical Background
//!
//! The demagnetizing field arises from the volume and surface magnetic charges of a
//! magnetized body. For a finite-difference grid of uniformly magnetized rectangular
//! prisms (cells), the demagnetizing field at cell i due to cell j is:
//!
//! ```text
//! H_demag(i) = −Σ_j  N̂(i − j) · M(j)
//! ```
//!
//! where N̂ is the **demagnetization tensor** and M(j) = Ms · m̂(j) [A/m].
//!
//! Because the grid is uniform, N̂ depends only on the cell-to-cell displacement,
//! and can be precomputed for all integer offsets.
//!
//! ## Demagnetization Tensor Formula
//!
//! The diagonal component N_xx for a rectangular cell of dimensions dx×dy×dz is
//! computed using the Newell (1993) analytic integration formula. The implementation
//! uses cell-face vertices (i*dx, j*dy, k*dz) with i,j,k ∈ {0,1} and alternating
//! signs, scaled by a factor of 16/(4π V_cell) to match the Newell normalization.
//!
//! For the self-cell (offset = 0,0,0), this formula satisfies:
//! - N_xx + N_yy + N_zz = 1 (demagnetization sum rule)
//! - For a cubic cell: N_xx = N_yy = N_zz = 1/3 (by symmetry)
//! - For a thin plate (dz ≪ dx,dy): N_zz ≫ N_xx ≈ N_yy
//!
//! For non-self cells, the same Newell 8-corner formula is applied with corners
//! at ((p+i−0.5)Δx, (q+j−0.5)Δy, (r+k−0.5)Δz), giving exact results for all
//! offsets including nearest neighbours where the point-dipole approximation fails.
//!
//! ## Simplification
//!
//! Only the diagonal components N_xx, N_yy, N_zz are stored. Off-diagonal terms
//! (N_xy, N_xz, N_yz) are set to zero. These vanish by symmetry for cells on the
//! same coordinate axis; for off-axis cells they represent small corrections.
//!
//! # References
//! - A. J. Newell, W. Williams, D. J. Dunlop, J. Geophys. Res. 98, 9551 (1993)
//! - A. Aharoni, J. Appl. Phys. 83, 3432 (1998) — exact self-demag formulas
//! - M. J. Donahue and D. G. Porter, OOMMF User's Guide, NIST IR 6376 (1999)

use crate::error::{invalid_param, Result};
use crate::vector3::Vector3;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// ─── Newell auxiliary function ───────────────────────────────────────────────

/// Newell auxiliary function for the diagonal demag component N_xx.
///
/// Evaluates the closed-form antiderivative from Newell et al. (1993) at a
/// single corner point (x, y, z). The 8-corner alternating sum of this function
/// gives the demag tensor element via the cell-face formula.
///
/// Uses `atan2(y*z, x*r)` (two-argument arctangent) to correctly handle
/// the sign of the quadrant for all corner positions, including those with
/// negative x-coordinates.
#[inline]
fn f_newell(x: f64, y: f64, z: f64) -> f64 {
    let x2 = x * x;
    let y2 = y * y;
    let z2 = z * z;
    let r = (x2 + y2 + z2).sqrt();
    if r < 1e-30 {
        return 0.0;
    }

    // Term 1: leading (2x²−y²−z²)/6 × r
    let t1 = (2.0 * x2 - y2 - z2) / 6.0 * r;

    // Term 2: y(z²−x²)/6 × ln(y+r)
    // Valid since r ≥ |y| ≥ 0 → y+r ≥ 0 always; strictly positive unless r=y=0.
    let t2 = if y + r > 1e-30 {
        y * (z2 - x2) / 6.0 * (y + r).ln()
    } else {
        0.0
    };

    // Term 3: z(y²−x²)/6 × ln(z+r)
    let t3 = if z + r > 1e-30 {
        z * (y2 - x2) / 6.0 * (z + r).ln()
    } else {
        0.0
    };

    // Term 4: −x·y·z/2 × atan2(y·z, x·r)
    // atan2 correctly handles the π-jump when x changes sign, which is essential
    // for the cell-face formula where x can be 0 or dx.
    let t4 = -x * y * z / 2.0 * f64::atan2(y * z, x * r);

    t1 + t2 + t3 + t4
}

// ─── Self-cell demag factors ─────────────────────────────────────────────────

/// Compute the self-demag factor N_xx (x-component) for a single isolated
/// rectangular prism cell of dimensions dx × dy × dz.
///
/// Uses the cell-face 8-corner formula with the Newell f function and a
/// normalization factor of 16/(4π V_cell). This formula gives the exact
/// self-demagnetization factor for any rectangular prism:
/// - Cubic cell (dx=dy=dz): N_xx = N_yy = N_zz = 1/3.
/// - Thin plate (dz ≪ dx,dy): N_zz → 1, N_xx → 0.
fn self_n_xx(dx: f64, dy: f64, dz: f64) -> f64 {
    // Cell-face corners at (i*dx, j*dy, k*dz) with i,j,k ∈ {0,1}
    // sign = (-1)^(i+j+k)
    let mut acc = 0.0_f64;
    for i in 0..=1_usize {
        for j in 0..=1_usize {
            for k in 0..=1_usize {
                let sign = if (i + j + k) % 2 == 0 { 1.0 } else { -1.0 };
                let xi = i as f64 * dx;
                let yj = j as f64 * dy;
                let zk = k as f64 * dz;
                acc += sign * f_newell(xi, yj, zk);
            }
        }
    }
    // Normalization: 16 / (4π dx dy dz)
    // The factor 16 = 2^4 arises from the cell-face vs cell-centre indexing
    // and ensures Tr(N) = 1 for the self-cell. See module documentation.
    16.0 * acc / (4.0 * std::f64::consts::PI * dx * dy * dz)
}

/// Compute the self-demag factor N_yy by permuting x ↔ y in f_newell.
fn self_n_yy(dx: f64, dy: f64, dz: f64) -> f64 {
    let mut acc = 0.0_f64;
    for i in 0..=1_usize {
        for j in 0..=1_usize {
            for k in 0..=1_usize {
                let sign = if (i + j + k) % 2 == 0 { 1.0 } else { -1.0 };
                let xi = i as f64 * dx;
                let yj = j as f64 * dy;
                let zk = k as f64 * dz;
                // N_yy: permute x ↔ y in f_newell
                acc += sign * f_newell(yj, xi, zk);
            }
        }
    }
    16.0 * acc / (4.0 * std::f64::consts::PI * dx * dy * dz)
}

/// Compute the self-demag factor N_zz by permuting x ↔ z in f_newell.
fn self_n_zz(dx: f64, dy: f64, dz: f64) -> f64 {
    let mut acc = 0.0_f64;
    for i in 0..=1_usize {
        for j in 0..=1_usize {
            for k in 0..=1_usize {
                let sign = if (i + j + k) % 2 == 0 { 1.0 } else { -1.0 };
                let xi = i as f64 * dx;
                let yj = j as f64 * dy;
                let zk = k as f64 * dz;
                // N_zz: permute x ↔ z in f_newell
                acc += sign * f_newell(zk, yj, xi);
            }
        }
    }
    16.0 * acc / (4.0 * std::f64::consts::PI * dx * dy * dz)
}

// ─── Non-self-cell demag: Newell 8-corner formula ────────────────────────────

/// Compute N_xx for a non-self cell at integer grid offset (p, q, r) using the
/// Newell (1993) 8-corner analytic formula.
///
/// For a cell displaced by (p, q, r) grid steps from the target cell, the
/// interaction tensor element is computed as the cell-face alternating corner sum
/// of the Newell antiderivative `f_newell`. The corners of the source cell in the
/// target-centred frame span from `|p|·Δx` to `(|p|+1)·Δx` (and similarly for y, z),
/// i.e., they are located at `(|p| + i)·Δx` for i ∈ {0, 1}:
///
/// ```text
/// N_xx(p,q,r) = 16/(4π V) × Σ_{i,j,k ∈ {0,1}} (−1)^(i+j+k) × f((|p|+i)Δx, (|q|+j)Δy, (|r|+k)Δz)
/// ```
///
/// Using absolute values of the offsets ensures N_xx(p,q,r) = N_xx(−p,−q,−r) by
/// construction, which is required by the symmetry of the dipolar interaction.
/// The prefactor `16/(4πV)` matches the normalization used by `self_n_xx`.
///
/// Physical sign conventions:
/// - N_xx(1,0,0) < 0: a neighbour along x reduces the x-demagnetization.
/// - N_xx(0,1,0) > 0: a neighbour along y increases the x-demagnetization.
/// - For a 2×2×2 cubic grid with uniform magnetization all cells give identical H_demag.
#[inline]
fn newell_n_xx(p: i64, q: i64, r: i64, dx: f64, dy: f64, dz: f64) -> f64 {
    let ap = p.unsigned_abs() as i64;
    let aq = q.unsigned_abs() as i64;
    let ar = r.unsigned_abs() as i64;
    let volume = dx * dy * dz;
    let mut acc = 0.0_f64;
    for i in 0..=1_i64 {
        for j in 0..=1_i64 {
            for k in 0..=1_i64 {
                let sign = if (i + j + k) % 2 == 0 {
                    1.0_f64
                } else {
                    -1.0_f64
                };
                let xi = (ap + i) as f64 * dx;
                let yj = (aq + j) as f64 * dy;
                let zk = (ar + k) as f64 * dz;
                acc += sign * f_newell(xi, yj, zk);
            }
        }
    }
    // Normalization factor 16/(4πV) matches the self_n_xx convention (Newell 1993).
    16.0 * acc / (4.0 * std::f64::consts::PI * volume)
}

/// Compute N_yy for a non-self cell using the Newell 8-corner formula.
///
/// N_yy is obtained from the N_xx antiderivative by swapping the roles of x and y
/// in the `f_newell` evaluation. Corner coordinates along each axis use the absolute
/// value of the corresponding offset component, as for `newell_n_xx`.
#[inline]
fn newell_n_yy(p: i64, q: i64, r: i64, dx: f64, dy: f64, dz: f64) -> f64 {
    let ap = p.unsigned_abs() as i64;
    let aq = q.unsigned_abs() as i64;
    let ar = r.unsigned_abs() as i64;
    let volume = dx * dy * dz;
    let mut acc = 0.0_f64;
    for i in 0..=1_i64 {
        for j in 0..=1_i64 {
            for k in 0..=1_i64 {
                let sign = if (i + j + k) % 2 == 0 {
                    1.0_f64
                } else {
                    -1.0_f64
                };
                let xi = (ap + i) as f64 * dx;
                let yj = (aq + j) as f64 * dy;
                let zk = (ar + k) as f64 * dz;
                // Swap x and y arguments to obtain N_yy from the N_xx antiderivative
                acc += sign * f_newell(yj, xi, zk);
            }
        }
    }
    16.0 * acc / (4.0 * std::f64::consts::PI * volume)
}

/// Compute N_zz for a non-self cell using the Newell 8-corner formula.
///
/// N_zz is obtained by permuting z into the leading argument position of `f_newell`.
/// Corner coordinates use absolute values of the offset components, matching the
/// convention of `newell_n_xx` and `newell_n_yy`.
#[inline]
fn newell_n_zz(p: i64, q: i64, r: i64, dx: f64, dy: f64, dz: f64) -> f64 {
    let ap = p.unsigned_abs() as i64;
    let aq = q.unsigned_abs() as i64;
    let ar = r.unsigned_abs() as i64;
    let volume = dx * dy * dz;
    let mut acc = 0.0_f64;
    for i in 0..=1_i64 {
        for j in 0..=1_i64 {
            for k in 0..=1_i64 {
                let sign = if (i + j + k) % 2 == 0 {
                    1.0_f64
                } else {
                    -1.0_f64
                };
                let xi = (ap + i) as f64 * dx;
                let yj = (aq + j) as f64 * dy;
                let zk = (ar + k) as f64 * dz;
                // Swap z into first argument position to obtain N_zz
                acc += sign * f_newell(zk, xi, yj);
            }
        }
    }
    16.0 * acc / (4.0 * std::f64::consts::PI * volume)
}

// ─── NewellTensor ────────────────────────────────────────────────────────────

/// Precomputed demagnetization tensor for a uniform rectangular grid.
///
/// Stores only the diagonal components N_xx, N_yy, N_zz (off-diagonal are zero).
///
/// All tensor elements (self-cell and non-self cells) are computed using the
/// Newell (1993) 8-corner analytic formula, which gives Tr(N_self) = 1 for any
/// rectangular prism and is exact for all cell-to-cell offsets.
///
/// ## Indexing
///
/// Offset (iu, iv, iw) ∈ [−(nx−1), +(nx−1)] × [−(ny−1), +(ny−1)] × [−(nz−1), +(nz−1)]
/// is stored at flat index:
/// ```text
/// iw_idx * (2ny−1) * (2nx−1) + iv_idx * (2nx−1) + iu_idx
/// ```
/// where `iu_idx = iu + (nx−1)` etc.
///
/// ## Sanity check
///
/// For a cubic cell (dx = dy = dz):
/// - N_xx(0,0,0) = N_yy(0,0,0) = N_zz(0,0,0) = 1/3
/// - Tr(N_self) = N_xx + N_yy + N_zz = 1
///
/// This is verified by `trace_check()` which returns the deviation from 1.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NewellTensor {
    /// Cell dimension along x \[m\]
    pub dx: f64,
    /// Cell dimension along y \[m\]
    pub dy: f64,
    /// Cell dimension along z \[m\]
    pub dz: f64,
    /// Number of cells along x
    pub nx: usize,
    /// Number of cells along y
    pub ny: usize,
    /// Number of cells along z
    pub nz: usize,
    /// N_xx values: shape (2nz−1) × (2ny−1) × (2nx−1), row-major in z-y-x order
    n_xx: Vec<f64>,
    /// N_yy values (same shape)
    n_yy: Vec<f64>,
    /// N_zz values (same shape)
    n_zz: Vec<f64>,
}

impl NewellTensor {
    /// Precompute the demagnetization tensor for an nx × ny × nz grid of cells.
    ///
    /// Covers all cell-to-cell offsets
    /// (iu, iv, iw) ∈ [−(nx−1), +(nx−1)] × [−(ny−1), +(ny−1)] × [−(nz−1), +(nz−1)].
    ///
    /// ## Errors
    /// Returns `Err` if any cell dimension is ≤ 0 or any grid dimension is 0.
    pub fn new(dx: f64, dy: f64, dz: f64, nx: usize, ny: usize, nz: usize) -> Result<Self> {
        if dx <= 0.0 {
            return Err(invalid_param("dx", "cell width must be positive"));
        }
        if dy <= 0.0 {
            return Err(invalid_param("dy", "cell height must be positive"));
        }
        if dz <= 0.0 {
            return Err(invalid_param("dz", "cell depth must be positive"));
        }
        if nx == 0 {
            return Err(invalid_param("nx", "grid dimension must be >= 1"));
        }
        if ny == 0 {
            return Err(invalid_param("ny", "grid dimension must be >= 1"));
        }
        if nz == 0 {
            return Err(invalid_param("nz", "grid dimension must be >= 1"));
        }

        let lu = 2 * nx - 1;
        let lv = 2 * ny - 1;
        let lw = 2 * nz - 1;
        let total = lu * lv * lw;

        let mut n_xx = vec![0.0_f64; total];
        let mut n_yy = vec![0.0_f64; total];
        let mut n_zz = vec![0.0_f64; total];

        for iw in 0..lw {
            let w = iw as i64 - (nz as i64 - 1);
            for iv in 0..lv {
                let v = iv as i64 - (ny as i64 - 1);
                for iu in 0..lu {
                    let u = iu as i64 - (nx as i64 - 1);
                    let flat = iw * lv * lu + iv * lu + iu;

                    if u == 0 && v == 0 && w == 0 {
                        // Self-cell: use the exact Newell cell-face formula.
                        // This gives Tr(N) = 1 and N_xx = N_yy = N_zz = 1/3 for cubes.
                        n_xx[flat] = self_n_xx(dx, dy, dz);
                        n_yy[flat] = self_n_yy(dx, dy, dz);
                        n_zz[flat] = self_n_zz(dx, dy, dz);
                    } else {
                        // Non-self cell: use the exact Newell 8-corner formula.
                        // This is valid for all offsets, including nearest neighbours,
                        // where the former point-dipole approximation was inaccurate.
                        n_xx[flat] = newell_n_xx(u, v, w, dx, dy, dz);
                        n_yy[flat] = newell_n_yy(u, v, w, dx, dy, dz);
                        n_zz[flat] = newell_n_zz(u, v, w, dx, dy, dz);
                    }
                }
            }
        }

        Ok(Self {
            dx,
            dy,
            dz,
            nx,
            ny,
            nz,
            n_xx,
            n_yy,
            n_zz,
        })
    }

    /// Convert a signed cell offset (iu, iv, iw) to a flat array index.
    ///
    /// Returns `None` if the offset is outside the precomputed range.
    #[inline]
    fn flat_index(&self, iu: i64, iv: i64, iw: i64) -> Option<usize> {
        let nx = self.nx as i64;
        let ny = self.ny as i64;
        let nz = self.nz as i64;
        if iu < -(nx - 1) || iu > (nx - 1) {
            return None;
        }
        if iv < -(ny - 1) || iv > (ny - 1) {
            return None;
        }
        if iw < -(nz - 1) || iw > (nz - 1) {
            return None;
        }
        let lu = (2 * self.nx - 1) as i64;
        let lv = (2 * self.ny - 1) as i64;
        let u_idx = iu + (nx - 1);
        let v_idx = iv + (ny - 1);
        let w_idx = iw + (nz - 1);
        Some((w_idx * lv * lu + v_idx * lu + u_idx) as usize)
    }

    /// Return the precomputed N_xx at integer cell offset (iu, iv, iw).
    ///
    /// Returns 0.0 for offsets outside the precomputed range.
    pub fn get_n_xx(&self, iu: i64, iv: i64, iw: i64) -> f64 {
        self.flat_index(iu, iv, iw)
            .map(|idx| self.n_xx[idx])
            .unwrap_or(0.0)
    }

    /// Return the precomputed N_yy at integer cell offset (iu, iv, iw).
    pub fn get_n_yy(&self, iu: i64, iv: i64, iw: i64) -> f64 {
        self.flat_index(iu, iv, iw)
            .map(|idx| self.n_yy[idx])
            .unwrap_or(0.0)
    }

    /// Return the precomputed N_zz at integer cell offset (iu, iv, iw).
    pub fn get_n_zz(&self, iu: i64, iv: i64, iw: i64) -> f64 {
        self.flat_index(iu, iv, iw)
            .map(|idx| self.n_zz[idx])
            .unwrap_or(0.0)
    }

    /// Return the self-demagnetization factors (N_xx, N_yy, N_zz) at offset (0,0,0).
    ///
    /// For a cubic cell: all three equal 1/3 and their sum equals 1.
    /// For a thin plate (dz ≪ dx, dy): N_zz ≫ N_xx ≈ N_yy.
    pub fn self_demag_factors(&self) -> (f64, f64, f64) {
        (
            self.get_n_xx(0, 0, 0),
            self.get_n_yy(0, 0, 0),
            self.get_n_zz(0, 0, 0),
        )
    }

    /// Return |N_xx(0,0,0) + N_yy(0,0,0) + N_zz(0,0,0) − 1|.
    ///
    /// This is the key sanity check: the trace of the self-demagnetization
    /// tensor must equal 1 (demagnetization sum rule). Values below 1e-6
    /// indicate a correct implementation.
    pub fn trace_check(&self) -> f64 {
        let (nxx, nyy, nzz) = self.self_demag_factors();
        (nxx + nyy + nzz - 1.0).abs()
    }

    /// Total number of cells in the grid.
    #[inline]
    pub fn n_cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }
}

// ─── DemagField ──────────────────────────────────────────────────────────────

/// Demagnetizing field evaluator using the precomputed NewellTensor.
///
/// Performs a direct O(N²) spatial convolution over all source-target cell pairs:
/// ```text
/// H_demag_x(i) = −Σ_j N_xx(i − j) × M_x(j)
/// H_demag_y(i) = −Σ_j N_yy(i − j) × M_y(j)
/// H_demag_z(i) = −Σ_j N_zz(i − j) × M_z(j)
/// ```
///
/// Off-diagonal tensor components N_xy, N_xz, N_yz are zero in this model.
///
/// For grids with N ≲ 100 cells this is fast enough for test and validation
/// purposes. Production codes should use an FFT-based convolution (O(N log N)).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DemagField {
    /// The precomputed Newell tensor
    pub tensor: NewellTensor,
}

impl DemagField {
    /// Wrap a precomputed [`NewellTensor`] in a field evaluator.
    pub fn new(tensor: NewellTensor) -> Self {
        Self { tensor }
    }

    /// Compute the demagnetizing field H_demag [A/m] for a given physical
    /// magnetization field M [A/m].
    ///
    /// ## Arguments
    /// - `magnetization` — physical magnetization M at each cell [A/m].
    ///   Length must equal `tensor.n_cells()`.
    ///   Cell indexing: `iz * ny * nx + iy * nx + ix`.
    ///
    /// ## Returns
    /// H_demag [A/m] at each cell (same length and indexing as input).
    ///
    /// Returns zero field (without panicking) if the slice length is wrong.
    pub fn compute(&self, magnetization: &[Vector3<f64>]) -> Vec<Vector3<f64>> {
        let nx = self.tensor.nx;
        let ny = self.tensor.ny;
        let nz = self.tensor.nz;
        let n_cells = nx * ny * nz;

        if magnetization.len() != n_cells {
            return vec![Vector3::zero(); n_cells];
        }

        // Kernel array dimensions (see `NewellTensor::flat_index`):
        //   lu = 2·nx − 1, lv = 2·ny − 1, lw = 2·nz − 1.
        // For any target index ix ∈ [0, nx) and source jx ∈ [0, nx) the offset
        // du = ix − jx is guaranteed to lie within [−(nx−1), nx−1] — exactly the
        // range the kernel arrays store — so `flat_index` could never return None
        // inside this convolution. We therefore index `tensor.n_xx/n_yy/n_zz`
        // directly: no Option, no bounds-branch. The (jz, jy, jx) summation order
        // and per-component accumulation are kept identical to the original
        // brute-force loop, so the floating-point result is bit-for-bit unchanged.
        let lu = 2 * nx - 1;
        let lv = 2 * ny - 1;

        // Compute H_demag for a single target cell whose linear index follows the
        // same ordering as the input: tgt = iz·ny·nx + iy·nx + ix.
        let compute_target = |tgt: usize| -> Vector3<f64> {
            let iz = tgt / (ny * nx);
            let rem = tgt % (ny * nx);
            let iy = rem / nx;
            let ix = rem % nx;

            // u_idx = du + (nx−1) = u_base − jx; with min u_base = nx−1 and
            // max jx = nx−1 this is always ≥ 0, so usize arithmetic never wraps.
            let u_base = ix + (nx - 1);

            let mut hx = 0.0_f64;
            let mut hy = 0.0_f64;
            let mut hz = 0.0_f64;

            for jz in 0..nz {
                // w_idx = dw + (nz−1) = (iz − jz) + (nz−1); always in [0, lw).
                let w_idx = iz + (nz - 1) - jz;
                for jy in 0..ny {
                    let v_idx = iy + (ny - 1) - jy;
                    let row_base = (w_idx * lv + v_idx) * lu;
                    let src_row = (jz * ny + jy) * nx;
                    for jx in 0..nx {
                        let kidx = row_base + (u_base - jx);
                        let m_src = magnetization[src_row + jx];

                        // Diagonal demag: H_α = −N_αα × M_α
                        // Off-diagonal N_αβ (α≠β) are zero.
                        hx -= self.tensor.n_xx[kidx] * m_src.x;
                        hy -= self.tensor.n_yy[kidx] * m_src.y;
                        hz -= self.tensor.n_zz[kidx] * m_src.z;
                    }
                }
            }

            Vector3::new(hx, hy, hz)
        };

        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            (0..n_cells).into_par_iter().map(compute_target).collect()
        }
        #[cfg(not(feature = "parallel"))]
        {
            (0..n_cells).map(compute_target).collect()
        }
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// For a cubic cell (dx = dy = dz) the self-demag trace must equal 1.0
    /// and each diagonal factor must equal 1/3.
    /// This is the demagnetization sum rule: Newell (1993) eq. (1).
    #[test]
    fn test_newell_trace_cubic_cell() {
        let d = 5e-9;
        let tensor =
            NewellTensor::new(d, d, d, 4, 4, 4).expect("tensor construction should succeed");
        let trace_err = tensor.trace_check();
        assert!(
            trace_err < 1e-6,
            "Trace error {trace_err:.2e} exceeds 1e-6 (demagnetization sum rule violated)"
        );
        let (nxx, nyy, nzz) = tensor.self_demag_factors();
        let expected = 1.0 / 3.0;
        assert!(
            (nxx - expected).abs() < 1e-6,
            "N_xx(0,0,0) = {nxx:.8} should be ≈ 1/3 = {expected:.8}"
        );
        assert!(
            (nyy - expected).abs() < 1e-6,
            "N_yy(0,0,0) = {nyy:.8} should be ≈ 1/3"
        );
        assert!(
            (nzz - expected).abs() < 1e-6,
            "N_zz(0,0,0) = {nzz:.8} should be ≈ 1/3"
        );
    }

    /// For a thin platelet (dz ≪ dx, dy) the out-of-plane self-demag N_zz should
    /// dominate over the in-plane factors. This reflects the fact that a flat
    /// disc is strongly demagnetized in the perpendicular direction.
    #[test]
    fn test_newell_thin_plate_demag() {
        let tensor = NewellTensor::new(500e-9, 500e-9, 1e-9, 2, 2, 1).expect("tensor construction");
        let (nxx, _nyy, nzz) = tensor.self_demag_factors();
        assert!(
            nzz > nxx * 5.0,
            "Thin plate: expected N_zz ({nzz:.4}) >> N_xx ({nxx:.4})"
        );
    }

    /// NewellTensor constructor should reject non-positive cell dimensions and
    /// zero grid dimensions.
    #[test]
    fn test_newell_invalid_params() {
        assert!(NewellTensor::new(-1.0, 5e-9, 5e-9, 2, 2, 2).is_err());
        assert!(NewellTensor::new(5e-9, 0.0, 5e-9, 2, 2, 2).is_err());
        assert!(NewellTensor::new(5e-9, 5e-9, 5e-9, 0, 2, 2).is_err());
    }

    /// DemagField output for a uniform x-magnetization should have zero H_y and H_z
    /// (by symmetry of the diagonal-only tensor model).
    #[test]
    fn test_demag_uniform_x_symmetry() {
        let d = 5e-9;
        let tensor = NewellTensor::new(d, d, d, 3, 3, 1).expect("tensor");
        let ms = 800e3_f64;
        let n = tensor.n_cells();
        let mag: Vec<Vector3<f64>> = vec![Vector3::new(ms, 0.0, 0.0); n];
        let demag = DemagField::new(tensor);
        let h = demag.compute(&mag);
        for (i, hi) in h.iter().enumerate() {
            assert!(
                hi.y.abs() < ms * 1e-10,
                "Cell {i}: H_y = {:.3e} should be ≈ 0",
                hi.y
            );
            assert!(
                hi.z.abs() < ms * 1e-10,
                "Cell {i}: H_z = {:.3e} should be ≈ 0",
                hi.z
            );
        }
    }

    /// DemagField should return zero for a zero magnetization field.
    #[test]
    fn test_demag_zero_magnetization() {
        let d = 5e-9;
        let tensor = NewellTensor::new(d, d, d, 2, 2, 2).expect("tensor");
        let n = tensor.n_cells();
        let mag = vec![Vector3::zero(); n];
        let demag = DemagField::new(tensor);
        let h = demag.compute(&mag);
        for hi in &h {
            assert!(hi.x.abs() < 1e-30 && hi.y.abs() < 1e-30 && hi.z.abs() < 1e-30);
        }
    }

    /// get_n_xx should return 0.0 for offsets outside the precomputed range.
    #[test]
    fn test_get_n_xx_out_of_range() {
        let d = 5e-9;
        let tensor = NewellTensor::new(d, d, d, 2, 2, 2).expect("tensor");
        let val = tensor.get_n_xx(10, 0, 0);
        assert_eq!(val, 0.0, "Out-of-range offset should return 0.0");
    }

    /// The Newell f function should return 0 at r = 0 and remain finite near zero.
    #[test]
    fn test_f_newell_zero_and_finite() {
        assert_eq!(f_newell(0.0, 0.0, 0.0), 0.0);
        assert!(f_newell(1e-25, 0.0, 0.0).is_finite());
        assert!(f_newell(5e-9, 2.5e-9, 1e-9).is_finite());
    }

    /// Self demag trace should be 1 for an asymmetric (non-cubic) cell too.
    #[test]
    fn test_newell_trace_asymmetric_cell() {
        let tensor = NewellTensor::new(10e-9, 5e-9, 2e-9, 2, 2, 2).expect("tensor");
        let trace_err = tensor.trace_check();
        assert!(
            trace_err < 1e-6,
            "Asymmetric cell trace error {trace_err:.2e} exceeds 1e-6"
        );
    }

    /// The optimized `DemagField::compute` (direct kernel indexing, optional
    /// rayon parallelism) must reproduce the original brute-force O(N²)
    /// convolution to the last bit. A deliberately NON-cubic grid (nx≠ny≠nz,
    /// dx≠dy≠dz) is used so any kernel-layout / stride bug would surface.
    /// The reference loop below mirrors the pre-optimization implementation
    /// exactly: Option-returning getters and jz→jy→jx summation order.
    #[test]
    fn test_compute_matches_reference_bruteforce() {
        let dx = 4e-9;
        let dy = 5e-9;
        let dz = 6e-9;
        let nx = 3;
        let ny = 4;
        let nz = 2;

        let tensor =
            NewellTensor::new(dx, dy, dz, nx, ny, nz).expect("tensor construction should succeed");
        let n = tensor.n_cells();

        // Deterministic, varied, non-trivial magnetization (no randomness): each
        // component is a distinct trig function of the cell index so no two cells
        // share a value and all three axes differ.
        let mut mag: Vec<Vector3<f64>> = Vec::with_capacity(n);
        for idx in 0..n {
            let t = idx as f64;
            mag.push(Vector3::new(
                (0.7 * t + 0.3).sin() * 8.0e5,
                (0.4 * t - 1.1).cos() * 5.0e5,
                ((t + 1.0) * 0.13).sin() * 6.0e5,
            ));
        }

        let h_fast = DemagField::new(tensor.clone()).compute(&mag);

        // Brute-force reference: byte-for-byte the original O(N²) algorithm.
        let mut h_ref = vec![Vector3::zero(); n];
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let tgt = iz * ny * nx + iy * nx + ix;
                    let mut hx = 0.0_f64;
                    let mut hy = 0.0_f64;
                    let mut hz = 0.0_f64;
                    for jz in 0..nz {
                        for jy in 0..ny {
                            for jx in 0..nx {
                                let src = jz * ny * nx + jy * nx + jx;
                                let m_src = mag[src];
                                let du = ix as i64 - jx as i64;
                                let dv = iy as i64 - jy as i64;
                                let dw = iz as i64 - jz as i64;
                                hx -= tensor.get_n_xx(du, dv, dw) * m_src.x;
                                hy -= tensor.get_n_yy(du, dv, dw) * m_src.y;
                                hz -= tensor.get_n_zz(du, dv, dw) * m_src.z;
                            }
                        }
                    }
                    h_ref[tgt] = Vector3::new(hx, hy, hz);
                }
            }
        }

        assert_eq!(h_fast.len(), h_ref.len(), "output length mismatch");

        // Combined absolute + relative tolerance; results should in fact be
        // bit-identical, so this passes with enormous margin.
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-9 + 1e-9 * b.abs();
        for (i, (fast, refv)) in h_fast.iter().zip(h_ref.iter()).enumerate() {
            assert!(
                close(fast.x, refv.x) && close(fast.y, refv.y) && close(fast.z, refv.z),
                "Cell {i}: optimized H = ({:.6e}, {:.6e}, {:.6e}) \
                 differs from brute-force reference ({:.6e}, {:.6e}, {:.6e})",
                fast.x,
                fast.y,
                fast.z,
                refv.x,
                refv.y,
                refv.z
            );
        }
    }
}
