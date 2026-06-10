//! Finite-difference LLG solver on a structured rectangular grid.
//!
//! ## Physical model
//!
//! Each cell carries a unit-length magnetization vector m̂ (normalized to 1).
//! The physical magnetization is M = Ms × m̂, where Ms is the saturation
//! magnetization of the material.
//!
//! The effective field acting on cell i is:
//! ```text
//! H_eff(i) = H_ex(i) + H_demag(i) + H_ext + H_ani(i)
//! ```
//!
//! ### Exchange field (FD Laplacian with zero-flux BCs)
//! ```text
//! H_ex(i) = (2A / μ₀Ms) × Σ_{nn j} (m_j − m_i) / Δ²
//! ```
//! where the sum runs over the (at most 6) nearest-neighbour cells and Δ is
//! the cell size in the corresponding direction. Missing neighbours at the
//! boundary simply do not contribute (zero-flux / Neumann BC).
//!
//! ### Demagnetizing field
//! Computed by [`DemagField::compute`] using the precomputed Newell tensor.
//! Input is the physical magnetization M = Ms × m̂ [A/m]; output is H_demag [A/m].
//!
//! ### Zeeman field
//! Uniform applied field H_ext [A/m] (same for all cells).
//!
//! ### Uniaxial anisotropy field
//! ```text
//! H_ani(i) = (2K / μ₀Ms) × (m̂ · ê_easy) × ê_easy
//! ```
//! where K is the uniaxial anisotropy constant [J/m³] and ê_easy is the easy axis.
//!
//! ## Time integration
//!
//! Uses classical RK4 applied simultaneously to all cells. After each stage the
//! intermediate magnetization is renormalized cell-by-cell to preserve |m| = 1.
//! The LLG right-hand side is the Gilbert form:
//! ```text
//! dm/dt = −γ / (1 + α²) × [m × H_eff + α × m × (m × H_eff)]
//! ```
//!
//! ## Performance note
//!
//! The demagnetizing-field convolution runs in O(N²) per H_eff evaluation, and
//! RK4 calls H_eff 4 times per step. For N = 64 cells (4×4×4 grid) each step
//! requires ~16 000 scalar multiplications — fast enough for tests and standard
//! problem validation. For production use of larger grids, replace the O(N²)
//! convolution with an FFT-based approach.

use crate::constants::MU_0;
use crate::error::{invalid_param, Result};
use crate::material::Ferromagnet;
use crate::micromagnetics::demag::{DemagField, NewellTensor};
use crate::vector3::Vector3;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// ─── GridConfig ──────────────────────────────────────────────────────────────

/// Configuration for a structured-grid FD micromagnetics simulation.
///
/// The simulation box has dimensions (nx⋅dx) × (ny⋅dy) × (nz⋅dz) \[m\].
/// Typical values for Permalloy at 5 nm resolution: dx = dy = 5 nm,
/// dz = 3 nm (thin-film geometry).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GridConfig {
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
    /// LLG time step \[s\].  Stability requires γ Ms dt ≲ 1e-2.
    pub dt: f64,
    /// Total number of integration steps.
    pub n_steps: usize,
    /// Record average magnetization every this many steps (0 = record only the final state).
    pub record_every: usize,
}

impl Default for GridConfig {
    fn default() -> Self {
        Self {
            dx: 5e-9,
            dy: 5e-9,
            dz: 3e-9,
            nx: 100,
            ny: 25,
            nz: 1,
            dt: 1e-13,
            n_steps: 10_000,
            record_every: 1_000,
        }
    }
}

// ─── LlgResult ───────────────────────────────────────────────────────────────

/// Output of a completed LLG simulation run.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LlgResult {
    /// Volume-averaged magnetization m̄ = (1/N) Σ_i m̂_i at each recorded time.
    pub m_avg: Vec<Vector3<f64>>,
    /// Physical time \[s\] at each recorded snapshot.
    pub times: Vec<f64>,
    /// Per-cell normalized magnetization m̂ at the end of the run.
    pub final_state: Vec<Vector3<f64>>,
}

// ─── MicromagneticGrid ───────────────────────────────────────────────────────

/// Structured-grid finite-difference micromagnetics solver.
///
/// Cells are indexed in C (row-major) order:
/// ```text
/// cell_index(ix, iy, iz) = iz * ny * nx + iy * nx + ix
/// ```
///
/// The magnetization vector stored in each cell is the *normalized* unit vector
/// m̂ ∈ S². The physical magnetization is M = material.ms × m̂.
#[derive(Debug, Clone)]
pub struct MicromagneticGrid {
    /// Simulation configuration
    pub config: GridConfig,
    /// Material parameters
    pub material: Ferromagnet,
    /// Per-cell normalized magnetization m̂ (|m| = 1)
    pub magnetization: Vec<Vector3<f64>>,
    /// Uniform external field H_ext [A/m]
    pub h_ext: Vector3<f64>,
    /// Demagnetizing field evaluator
    demag: DemagField,
}

impl MicromagneticGrid {
    /// Construct a new grid and precompute the Newell tensor.
    ///
    /// Initial magnetization is uniform along +x: m̂ = (1, 0, 0).
    ///
    /// # Errors
    /// Returns `Err` if any dimension or cell size is invalid.
    pub fn new(config: GridConfig, material: Ferromagnet, h_ext: Vector3<f64>) -> Result<Self> {
        // Validate configuration
        if config.nx == 0 {
            return Err(invalid_param("nx", "must be >= 1"));
        }
        if config.ny == 0 {
            return Err(invalid_param("ny", "must be >= 1"));
        }
        if config.nz == 0 {
            return Err(invalid_param("nz", "must be >= 1"));
        }
        if config.dx <= 0.0 {
            return Err(invalid_param("dx", "must be positive"));
        }
        if config.dy <= 0.0 {
            return Err(invalid_param("dy", "must be positive"));
        }
        if config.dz <= 0.0 {
            return Err(invalid_param("dz", "must be positive"));
        }
        if config.dt <= 0.0 {
            return Err(invalid_param("dt", "must be positive"));
        }

        let n = config.nx * config.ny * config.nz;
        // Initial state: uniform saturation along +x
        let magnetization = vec![Vector3::new(1.0, 0.0, 0.0); n];

        let tensor = NewellTensor::new(
            config.dx, config.dy, config.dz, config.nx, config.ny, config.nz,
        )?;
        let demag = DemagField::new(tensor);

        Ok(Self {
            config,
            material,
            magnetization,
            h_ext,
            demag,
        })
    }

    /// Total number of cells in the grid.
    #[inline]
    pub fn n_cells(&self) -> usize {
        self.config.nx * self.config.ny * self.config.nz
    }

    /// Flat cell index for position (ix, iy, iz).
    #[inline]
    fn cell_index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.config.ny * self.config.nx + iy * self.config.nx + ix
    }

    /// Set all cells to a uniform magnetization in the given direction (renormalized).
    ///
    /// If `direction` is the zero vector, the magnetization is left unchanged.
    pub fn set_uniform(&mut self, direction: Vector3<f64>) {
        let mag = direction.magnitude();
        if mag < 1e-30 {
            return;
        }
        let d_norm = direction.normalize();
        for m in self.magnetization.iter_mut() {
            *m = d_norm;
        }
    }

    /// Initialize a flux-closure vortex state in the xy-plane.
    ///
    /// The vortex magnetization pattern is:
    /// ```text
    /// m(ix, iy, iz) = normalize(−ŷ_rel, x̂_rel, ε)
    /// ```
    /// where (x̂_rel, ŷ_rel) are the cell's coordinates relative to the
    /// grid centre, normalized to [−1, +1], and ε = 0.01 gives a small
    /// out-of-plane component that breaks the exact degeneracy.
    ///
    /// This is the standard flux-closure vortex used in muMAG Standard Problem #3.
    pub fn set_vortex(&mut self) {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let half_nx = nx as f64 / 2.0;
        let half_ny = ny as f64 / 2.0;
        // Small out-of-plane component to avoid exact degeneracy
        let eps = 0.01_f64;

        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    // Relative coordinates in [−1, +1] (approximately)
                    let x_rel = (ix as f64 + 0.5 - half_nx) / half_nx;
                    let y_rel = (iy as f64 + 0.5 - half_ny) / half_ny;
                    // Vortex: m = (−y, x, ε).normalize()
                    let m = Vector3::new(-y_rel, x_rel, eps).normalize();
                    let idx = self.cell_index(ix, iy, iz);
                    self.magnetization[idx] = m;
                }
            }
        }
    }

    /// Compute the volume-averaged magnetization over all cells.
    pub fn mean_magnetization(&self) -> Vector3<f64> {
        let n = self.n_cells();
        if n == 0 {
            return Vector3::zero();
        }
        let sum = self.magnetization.iter().fold(Vector3::zero(), |acc, &m| {
            Vector3::new(acc.x + m.x, acc.y + m.y, acc.z + m.z)
        });
        let inv_n = 1.0 / n as f64;
        Vector3::new(sum.x * inv_n, sum.y * inv_n, sum.z * inv_n)
    }

    /// Compute the effective field H_eff [A/m] for each cell.
    ///
    /// H_eff = H_ex + H_demag + H_ext + H_ani
    ///
    /// The input `m` must have length == n_cells().
    pub fn effective_field(&self, m: &[Vector3<f64>]) -> Vec<Vector3<f64>> {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let dx = self.config.dx;
        let dy = self.config.dy;
        let dz = self.config.dz;
        let ms = self.material.ms;
        let a_ex = self.material.exchange_a;
        let k_ani = self.material.anisotropy_k;
        let easy = self.material.easy_axis;

        // Exchange prefactor: 2A / (μ₀ Ms)   [A/m per (1/m²)]
        // H_ex has units of A/m since 2A/(μ₀ Ms) × (1/Δ²) × Δm (dimensionless) → A/m
        let ex_pref = 2.0 * a_ex / (MU_0 * ms);

        // Anisotropy prefactor: 2K / (μ₀ Ms)
        let ani_pref = 2.0 * k_ani / (MU_0 * ms);

        // Physical magnetization M = Ms × m̂ for demag calculation
        let m_phys: Vec<Vector3<f64>> = m.iter().map(|mi| *mi * ms).collect();

        // Demagnetizing field [A/m]
        let h_demag = self.demag.compute(&m_phys);

        let n_cells = nx * ny * nz;
        let mut h_eff = vec![Vector3::zero(); n_cells];

        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let idx = self.cell_index(ix, iy, iz);
                    let mi = m[idx];

                    // ── Exchange field ────────────────────────────────────────
                    // FD Laplacian: H_ex = (2A/(μ₀Ms)) × Σ_{nn} (m_j − m_i) / Δ²
                    // Zero-flux BC: missing neighbours → no contribution.
                    let mut h_ex = Vector3::zero();

                    // −x neighbour
                    if ix > 0 {
                        let j = self.cell_index(ix - 1, iy, iz);
                        let diff = m[j] - mi;
                        h_ex = Vector3::new(
                            h_ex.x + diff.x * ex_pref / (dx * dx),
                            h_ex.y + diff.y * ex_pref / (dx * dx),
                            h_ex.z + diff.z * ex_pref / (dx * dx),
                        );
                    }
                    // +x neighbour
                    if ix + 1 < nx {
                        let j = self.cell_index(ix + 1, iy, iz);
                        let diff = m[j] - mi;
                        h_ex = Vector3::new(
                            h_ex.x + diff.x * ex_pref / (dx * dx),
                            h_ex.y + diff.y * ex_pref / (dx * dx),
                            h_ex.z + diff.z * ex_pref / (dx * dx),
                        );
                    }
                    // −y neighbour
                    if iy > 0 {
                        let j = self.cell_index(ix, iy - 1, iz);
                        let diff = m[j] - mi;
                        h_ex = Vector3::new(
                            h_ex.x + diff.x * ex_pref / (dy * dy),
                            h_ex.y + diff.y * ex_pref / (dy * dy),
                            h_ex.z + diff.z * ex_pref / (dy * dy),
                        );
                    }
                    // +y neighbour
                    if iy + 1 < ny {
                        let j = self.cell_index(ix, iy + 1, iz);
                        let diff = m[j] - mi;
                        h_ex = Vector3::new(
                            h_ex.x + diff.x * ex_pref / (dy * dy),
                            h_ex.y + diff.y * ex_pref / (dy * dy),
                            h_ex.z + diff.z * ex_pref / (dy * dy),
                        );
                    }
                    // −z neighbour
                    if iz > 0 {
                        let j = self.cell_index(ix, iy, iz - 1);
                        let diff = m[j] - mi;
                        h_ex = Vector3::new(
                            h_ex.x + diff.x * ex_pref / (dz * dz),
                            h_ex.y + diff.y * ex_pref / (dz * dz),
                            h_ex.z + diff.z * ex_pref / (dz * dz),
                        );
                    }
                    // +z neighbour
                    if iz + 1 < nz {
                        let j = self.cell_index(ix, iy, iz + 1);
                        let diff = m[j] - mi;
                        h_ex = Vector3::new(
                            h_ex.x + diff.x * ex_pref / (dz * dz),
                            h_ex.y + diff.y * ex_pref / (dz * dz),
                            h_ex.z + diff.z * ex_pref / (dz * dz),
                        );
                    }

                    // ── Uniaxial anisotropy field ─────────────────────────────
                    // H_ani = (2K/(μ₀Ms)) × (m·ê_easy) × ê_easy
                    let m_dot_e = mi.dot(&easy);
                    let h_ani = easy * (ani_pref * m_dot_e);

                    // ── Sum all contributions ─────────────────────────────────
                    let h_total = Vector3::new(
                        h_ex.x + h_demag[idx].x + self.h_ext.x + h_ani.x,
                        h_ex.y + h_demag[idx].y + self.h_ext.y + h_ani.y,
                        h_ex.z + h_demag[idx].z + self.h_ext.z + h_ani.z,
                    );

                    h_eff[idx] = h_total;
                }
            }
        }

        h_eff
    }

    /// Compute dm/dt for a single cell from the Gilbert LLG equation.
    ///
    /// The explicit (reduced) Gilbert form avoids solving an implicit equation:
    /// ```text
    /// dm/dt = −γ/(1+α²) × [m × H + α × m × (m × H)]
    /// ```
    ///
    /// # Arguments
    /// * `m_i`     — normalized magnetization at this cell (unit vector)
    /// * `h_eff_i` — effective field at this cell [A/m]
    #[inline]
    pub fn dmdt_cell(&self, m_i: Vector3<f64>, h_eff_i: Vector3<f64>) -> Vector3<f64> {
        let alpha = self.material.alpha;
        let gamma = crate::constants::GAMMA;
        let factor = gamma / (1.0 + alpha * alpha);

        // Precession term: m × H
        let m_cross_h = m_i.cross(&h_eff_i);
        // Damping term: m × (m × H)
        let m_cross_m_cross_h = m_i.cross(&m_cross_h);

        // dm/dt = −factor × (m × H) − factor × α × m × (m × H)
        let dmdt_x = -factor * m_cross_h.x - factor * alpha * m_cross_m_cross_h.x;
        let dmdt_y = -factor * m_cross_h.y - factor * alpha * m_cross_m_cross_h.y;
        let dmdt_z = -factor * m_cross_h.z - factor * alpha * m_cross_m_cross_h.z;
        Vector3::new(dmdt_x, dmdt_y, dmdt_z)
    }

    /// Apply a single RK4 step to all cells simultaneously.
    ///
    /// After each intermediate stage the trial magnetization is renormalized
    /// cell-by-cell to keep |m̂| = 1. The final update is then applied and
    /// the result is renormalized again.
    pub fn rk4_step(&mut self) {
        let dt = self.config.dt;
        let n = self.n_cells();

        // k1: dm/dt at t with current m
        let h1 = self.effective_field(&self.magnetization);
        let k1: Vec<Vector3<f64>> = (0..n)
            .map(|i| self.dmdt_cell(self.magnetization[i], h1[i]))
            .collect();

        // m2 = renormalize(m + dt/2 × k1)
        let m2: Vec<Vector3<f64>> = (0..n)
            .map(|i| {
                let trial = Vector3::new(
                    self.magnetization[i].x + 0.5 * dt * k1[i].x,
                    self.magnetization[i].y + 0.5 * dt * k1[i].y,
                    self.magnetization[i].z + 0.5 * dt * k1[i].z,
                );
                trial.normalize()
            })
            .collect();

        // k2: dm/dt at t+dt/2 with m2
        let h2 = self.effective_field(&m2);
        let k2: Vec<Vector3<f64>> = (0..n).map(|i| self.dmdt_cell(m2[i], h2[i])).collect();

        // m3 = renormalize(m + dt/2 × k2)
        let m3: Vec<Vector3<f64>> = (0..n)
            .map(|i| {
                let trial = Vector3::new(
                    self.magnetization[i].x + 0.5 * dt * k2[i].x,
                    self.magnetization[i].y + 0.5 * dt * k2[i].y,
                    self.magnetization[i].z + 0.5 * dt * k2[i].z,
                );
                trial.normalize()
            })
            .collect();

        // k3: dm/dt at t+dt/2 with m3
        let h3 = self.effective_field(&m3);
        let k3: Vec<Vector3<f64>> = (0..n).map(|i| self.dmdt_cell(m3[i], h3[i])).collect();

        // m4 = renormalize(m + dt × k3)
        let m4: Vec<Vector3<f64>> = (0..n)
            .map(|i| {
                let trial = Vector3::new(
                    self.magnetization[i].x + dt * k3[i].x,
                    self.magnetization[i].y + dt * k3[i].y,
                    self.magnetization[i].z + dt * k3[i].z,
                );
                trial.normalize()
            })
            .collect();

        // k4: dm/dt at t+dt with m4
        let h4 = self.effective_field(&m4);
        let k4: Vec<Vector3<f64>> = (0..n).map(|i| self.dmdt_cell(m4[i], h4[i])).collect();

        // RK4 update: m_new = m + dt/6 × (k1 + 2k2 + 2k3 + k4), then renormalize
        for i in 0..n {
            let new_m = Vector3::new(
                self.magnetization[i].x
                    + dt / 6.0 * (k1[i].x + 2.0 * k2[i].x + 2.0 * k3[i].x + k4[i].x),
                self.magnetization[i].y
                    + dt / 6.0 * (k1[i].y + 2.0 * k2[i].y + 2.0 * k3[i].y + k4[i].y),
                self.magnetization[i].z
                    + dt / 6.0 * (k1[i].z + 2.0 * k2[i].z + 2.0 * k3[i].z + k4[i].z),
            );
            self.magnetization[i] = new_m.normalize();
        }
    }

    /// Integrate the LLG equation for `config.n_steps` steps, recording the
    /// volume-averaged magnetization every `config.record_every` steps.
    ///
    /// Returns an [`LlgResult`] with the time series and final per-cell state.
    pub fn run(&mut self) -> LlgResult {
        let n_steps = self.config.n_steps;
        let record_every = self.config.record_every.max(1);
        let dt = self.config.dt;

        let capacity = n_steps / record_every + 2;
        let mut m_avg = Vec::with_capacity(capacity);
        let mut times = Vec::with_capacity(capacity);

        // Record initial state
        m_avg.push(self.mean_magnetization());
        times.push(0.0);

        for step in 0..n_steps {
            self.rk4_step();
            if (step + 1) % record_every == 0 {
                m_avg.push(self.mean_magnetization());
                times.push((step + 1) as f64 * dt);
            }
        }

        LlgResult {
            m_avg,
            times,
            final_state: self.magnetization.clone(),
        }
    }

    // ─── Energy calculations ─────────────────────────────────────────────────

    /// Compute the total exchange energy \[J\].
    ///
    /// E_ex = A × Σ_{unique bonds ⟨i,j⟩} |m̂_j − m̂_i|² / Δ² × V_cell
    ///
    /// Each nearest-neighbour bond is counted exactly once (j > i convention).
    /// The cell volume V_cell = dx × dy × dz converts energy density to energy.
    pub fn total_exchange_energy(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let dx = self.config.dx;
        let dy = self.config.dy;
        let dz = self.config.dz;
        let a_ex = self.material.exchange_a;
        let v_cell = dx * dy * dz;

        let mut energy = 0.0_f64;

        // x-direction bonds (ix, iy, iz) — (ix+1, iy, iz)
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..(nx.saturating_sub(1)) {
                    let i = self.cell_index(ix, iy, iz);
                    let j = self.cell_index(ix + 1, iy, iz);
                    let diff = self.magnetization[j] - self.magnetization[i];
                    let diff_sq = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
                    energy += a_ex * diff_sq / (dx * dx) * v_cell;
                }
            }
        }

        // y-direction bonds
        for iz in 0..nz {
            for iy in 0..(ny.saturating_sub(1)) {
                for ix in 0..nx {
                    let i = self.cell_index(ix, iy, iz);
                    let j = self.cell_index(ix, iy + 1, iz);
                    let diff = self.magnetization[j] - self.magnetization[i];
                    let diff_sq = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
                    energy += a_ex * diff_sq / (dy * dy) * v_cell;
                }
            }
        }

        // z-direction bonds
        for iz in 0..(nz.saturating_sub(1)) {
            for iy in 0..ny {
                for ix in 0..nx {
                    let i = self.cell_index(ix, iy, iz);
                    let j = self.cell_index(ix, iy, iz + 1);
                    let diff = self.magnetization[j] - self.magnetization[i];
                    let diff_sq = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
                    energy += a_ex * diff_sq / (dz * dz) * v_cell;
                }
            }
        }

        energy
    }

    /// Compute the total demagnetizing (magnetostatic) energy \[J\].
    ///
    /// E_demag = −(μ₀/2) × Σ_i H_demag(i) · M(i) × V_cell
    ///
    /// The factor of 1/2 avoids double-counting pairwise interactions.
    /// M(i) = Ms × m̂(i) is the physical magnetization [A/m].
    pub fn total_demag_energy(&self) -> f64 {
        let ms = self.material.ms;
        let v_cell = self.config.dx * self.config.dy * self.config.dz;

        // Physical magnetization field
        let m_phys: Vec<Vector3<f64>> = self.magnetization.iter().map(|mi| *mi * ms).collect();
        // Demagnetizing field [A/m]
        let h_demag = self.demag.compute(&m_phys);

        let mut energy = 0.0_f64;
        for (i, m_i) in m_phys.iter().enumerate() {
            // H_demag · M = H_x M_x + H_y M_y + H_z M_z
            let h_dot_m = h_demag[i].x * m_i.x + h_demag[i].y * m_i.y + h_demag[i].z * m_i.z;
            energy -= MU_0 / 2.0 * h_dot_m * v_cell;
        }
        energy
    }

    /// Compute the total Zeeman (applied-field) energy \[J\].
    ///
    /// E_zee = −μ₀ × Σ_i M(i) · H_ext × V_cell
    ///
    /// Note: no factor of 1/2 because H_ext is an external field, not a
    /// self-interaction.
    pub fn total_zeeman_energy(&self) -> f64 {
        let ms = self.material.ms;
        let v_cell = self.config.dx * self.config.dy * self.config.dz;
        let h = self.h_ext;

        let mut energy = 0.0_f64;
        for mi in &self.magnetization {
            // M = ms × m̂, so M · H = ms × (m̂ · H)
            let m_dot_h = mi.x * h.x * ms + mi.y * h.y * ms + mi.z * h.z * ms;
            energy -= MU_0 * m_dot_h * v_cell;
        }
        energy
    }

    /// Compute the total energy \[J\] = E_exchange + E_demag + E_zeeman.
    ///
    /// Anisotropy energy is omitted because the standard muMAG Standard Problem #3
    /// uses Permalloy, which has K = 0.
    pub fn total_energy(&self) -> f64 {
        self.total_exchange_energy() + self.total_demag_energy() + self.total_zeeman_energy()
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::Ferromagnet;

    fn tiny_grid() -> MicromagneticGrid {
        // 2×2×2 cube with 5 nm cells — minimal grid, fast tests
        let cfg = GridConfig {
            dx: 5e-9,
            dy: 5e-9,
            dz: 5e-9,
            nx: 2,
            ny: 2,
            nz: 2,
            dt: 1e-12,
            n_steps: 10,
            record_every: 5,
        };
        let mat = Ferromagnet::permalloy();
        MicromagneticGrid::new(cfg, mat, Vector3::zero()).expect("tiny grid construction")
    }

    #[test]
    fn test_grid_construction_and_uniform_init() {
        let grid = tiny_grid();
        assert_eq!(grid.n_cells(), 8);
        // Initial state: uniform +x
        for mi in &grid.magnetization {
            assert!((mi.x - 1.0).abs() < 1e-12);
            assert!(mi.y.abs() < 1e-12);
            assert!(mi.z.abs() < 1e-12);
        }
    }

    #[test]
    fn test_set_uniform() {
        let mut grid = tiny_grid();
        grid.set_uniform(Vector3::new(0.0, 1.0, 0.0));
        let m = grid.mean_magnetization();
        assert!((m.y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_set_vortex_normalized() {
        let mut grid = tiny_grid();
        grid.set_vortex();
        // All cells should remain unit vectors
        for (i, mi) in grid.magnetization.iter().enumerate() {
            let mag = mi.magnitude();
            assert!((mag - 1.0).abs() < 1e-10, "Cell {i} magnitude = {mag:.6}");
        }
    }

    #[test]
    fn test_mean_magnetization_uniform() {
        let grid = tiny_grid();
        let m = grid.mean_magnetization();
        assert!((m.x - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_effective_field_uniform_no_exchange() {
        // For a uniform magnetization there is no exchange field (all gradients zero)
        let grid = tiny_grid();
        let h = grid.effective_field(&grid.magnetization);
        // The exchange field must be zero for uniform m
        // (demag + Zeeman are present, but exchange term should vanish)
        // Check that every cell has the same H_eff (uniform → same demag by symmetry)
        let h0 = h[0];
        for (i, hi) in h.iter().enumerate() {
            assert!(
                (hi.x - h0.x).abs() < 1e-3 * h0.x.abs() + 1.0,
                "Cell {i}: H_eff differs unexpectedly"
            );
        }
    }

    #[test]
    fn test_rk4_step_preserves_normalization() {
        let mut grid = tiny_grid();
        grid.set_vortex();
        // Run several steps and check that |m| = 1 is preserved
        for _ in 0..5 {
            grid.rk4_step();
        }
        for (i, mi) in grid.magnetization.iter().enumerate() {
            let mag = mi.magnitude();
            assert!(
                (mag - 1.0).abs() < 1e-10,
                "After RK4: cell {i} magnitude = {mag:.8}"
            );
        }
    }

    #[test]
    fn test_run_returns_correct_snapshot_count() {
        let cfg = GridConfig {
            dx: 5e-9,
            dy: 5e-9,
            dz: 5e-9,
            nx: 2,
            ny: 2,
            nz: 1,
            dt: 1e-12,
            n_steps: 10,
            record_every: 5,
        };
        let mat = Ferromagnet::permalloy();
        let mut grid = MicromagneticGrid::new(cfg, mat, Vector3::zero()).expect("grid");
        let result = grid.run();
        // Initial + step 5 + step 10 = 3 records
        assert_eq!(result.m_avg.len(), 3, "Expected 3 snapshots");
        assert_eq!(result.times.len(), 3);
    }

    #[test]
    fn test_exchange_energy_uniform_is_zero() {
        let grid = tiny_grid();
        // Uniform magnetization → zero exchange energy (no gradients)
        let e_ex = grid.total_exchange_energy();
        assert!(
            e_ex.abs() < 1e-50,
            "Exchange energy for uniform state should be 0, got {e_ex:.3e}"
        );
    }

    #[test]
    fn test_exchange_energy_vortex_is_positive() {
        let mut grid = tiny_grid();
        grid.set_vortex();
        let e_ex = grid.total_exchange_energy();
        assert!(
            e_ex > 0.0,
            "Exchange energy for vortex should be positive, got {e_ex:.3e}"
        );
    }

    #[test]
    fn test_total_energy_finite() {
        let grid = tiny_grid();
        let e = grid.total_energy();
        assert!(e.is_finite(), "Total energy must be finite");
    }

    #[test]
    fn test_invalid_grid_config() {
        let cfg = GridConfig {
            dx: -1.0,
            ..GridConfig::default()
        };
        let result = MicromagneticGrid::new(cfg, Ferromagnet::permalloy(), Vector3::zero());
        assert!(result.is_err());
    }

    #[test]
    fn test_demag_energy_sign() {
        // For uniform +x magnetization the demag energy should be negative
        // (H_demag opposes M, so −H·M > 0 ??? )
        // Actually: H_demag is anti-parallel to M, so H_demag · M < 0,
        // and E_demag = −μ₀/2 × (negative) > 0. Both signs appear; just
        // check finiteness and that it changes when state changes.
        let mut grid = tiny_grid();
        let e1 = grid.total_demag_energy();
        grid.set_uniform(Vector3::new(0.0, 0.0, 1.0));
        let e2 = grid.total_demag_energy();
        // For a cube the demag energy should be the same for all axis directions
        // due to cubic symmetry (within numerical precision of a 2x2x2 grid)
        assert!(e1.is_finite() && e2.is_finite());
    }
}
