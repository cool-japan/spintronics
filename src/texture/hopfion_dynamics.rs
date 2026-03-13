//! Hopfion Dynamics — LLG-based time evolution of 3D hopfion spin textures
//!
//! This module implements the Landau-Lifshitz-Gilbert (LLG) equation solver for
//! three-dimensional magnetic textures on a regular cubic lattice, with support
//! for:
//!
//! - Exchange interaction (6-point Laplacian stencil, periodic BC)
//! - Bulk DMI (curl-based antisymmetric exchange)
//! - External Zeeman field
//! - Fourth-order Runge-Kutta (RK4) time integration
//! - Periodic tracking of the Hopf invariant and total energy
//!
//! # Physical Model
//!
//! The effective field acting on each spin is
//!
//! ```text
//! H_eff = H_ex + H_DMI + H_ext
//! ```
//!
//! and the LLG equation of motion reads
//!
//! ```text
//! dm/dt = -γ / (1 + α²) · [m × H_eff + α · m × (m × H_eff)]
//! ```
//!
//! After each RK4 step the spin is re-normalised to |m| = 1.

use std::f64::consts::PI;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::error::Error;
use crate::vector3::Vector3;

/// Configuration for hopfion dynamics simulations
///
/// All fields use SI units unless stated otherwise.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HopfionDynamicsConfig {
    /// Exchange stiffness A \[J/m\]
    pub exchange_a: f64,
    /// Bulk DMI constant D \[J/m²\]
    pub dmi_d: f64,
    /// Gilbert damping parameter α (dimensionless)
    pub alpha: f64,
    /// External magnetic field H_ext \[T\]
    pub h_ext: Vector3<f64>,
    /// Saturation magnetization Mₛ \[A/m\]
    pub ms: f64,
    /// Simulation time step Δt \[s\]
    pub dt: f64,
    /// Lattice grid spacing Δx \[m\]
    pub dx: f64,
}

impl HopfionDynamicsConfig {
    /// Create a new dynamics configuration
    ///
    /// # Arguments
    /// * `exchange_a` – Exchange stiffness \[J/m\]
    /// * `dmi_d`      – Bulk DMI constant \[J/m²\]
    /// * `alpha`      – Gilbert damping (dimensionless ≥ 0)
    /// * `h_ext`      – External field vector \[T\]
    /// * `ms`         – Saturation magnetization \[A/m\]
    /// * `dt`         – Time step \[s\]
    /// * `dx`         – Grid spacing \[m\]
    pub fn new(
        exchange_a: f64,
        dmi_d: f64,
        alpha: f64,
        h_ext: Vector3<f64>,
        ms: f64,
        dt: f64,
        dx: f64,
    ) -> Self {
        Self {
            exchange_a,
            dmi_d,
            alpha,
            h_ext,
            ms,
            dt,
            dx,
        }
    }

    /// Validate physical reasonableness of all parameters
    ///
    /// # Errors
    /// Returns `Error::InvalidParameter` if any of `ms`, `dt`, or `dx` are
    /// non-positive or non-finite.
    pub fn validate(&self) -> Result<(), Error> {
        if self.ms <= 0.0 || !self.ms.is_finite() {
            return Err(Error::InvalidParameter {
                param: "ms".to_string(),
                reason: "Saturation magnetization must be positive and finite".to_string(),
            });
        }
        if self.dt <= 0.0 || !self.dt.is_finite() {
            return Err(Error::InvalidParameter {
                param: "dt".to_string(),
                reason: "Time step must be positive and finite".to_string(),
            });
        }
        if self.dx <= 0.0 || !self.dx.is_finite() {
            return Err(Error::InvalidParameter {
                param: "dx".to_string(),
                reason: "Grid spacing must be positive and finite".to_string(),
            });
        }
        if self.alpha < 0.0 || !self.alpha.is_finite() {
            return Err(Error::InvalidParameter {
                param: "alpha".to_string(),
                reason: "Gilbert damping must be non-negative and finite".to_string(),
            });
        }
        Ok(())
    }
}

// ─── Solver ──────────────────────────────────────────────────────────────────

/// LLG solver for hopfion dynamics on a 3D cubic grid
///
/// Integrates the Landau-Lifshitz-Gilbert equation with exchange, bulk DMI, and
/// Zeeman terms using a standard fourth-order Runge-Kutta scheme.
#[derive(Debug)]
pub struct HopfionDynamicsSolver {
    config: HopfionDynamicsConfig,
    /// Gyromagnetic ratio γ \[rad s⁻¹ T⁻¹\]
    gamma: f64,
}

impl HopfionDynamicsSolver {
    /// Gyromagnetic ratio for free electrons (rad s⁻¹ T⁻¹)
    pub const GAMMA_DEFAULT: f64 = 1.760_859_644e11;

    /// Permeability of free space μ₀ \[T·m/A\]
    const MU0: f64 = 4.0 * PI * 1.0e-7;

    /// Build a solver from a validated configuration.
    ///
    /// # Errors
    /// Propagates any validation errors from [`HopfionDynamicsConfig::validate`].
    pub fn new(config: HopfionDynamicsConfig) -> Result<Self, Error> {
        config.validate()?;
        Ok(Self {
            config,
            gamma: Self::GAMMA_DEFAULT,
        })
    }

    // ── Public entry point ────────────────────────────────────────────────

    /// Run the LLG dynamics on a 3-D magnetization grid.
    ///
    /// # Arguments
    /// * `grid`         – `[nx][ny][nz]` grid of unit-length magnetization vectors
    /// * `num_steps`    – Total number of time-integration steps
    /// * `record_every` – Record Hopf invariant and energy every N steps
    ///
    /// # Returns
    /// A [`HopfionDynamicsResult`] containing the history of Hopf invariant,
    /// total energy, and time sampled at each recording point.
    ///
    /// # Errors
    /// Returns `Error::InvalidParameter` if the grid is empty.
    pub fn run(
        &self,
        mut grid: Vec<Vec<Vec<Vector3<f64>>>>,
        num_steps: usize,
        record_every: usize,
    ) -> Result<HopfionDynamicsResult, Error> {
        // Validate grid dimensions
        let nx = grid.len();
        if nx == 0 {
            return Err(Error::InvalidParameter {
                param: "grid".to_string(),
                reason: "Grid must be non-empty (nx > 0)".to_string(),
            });
        }
        let ny = grid[0].len();
        if ny == 0 {
            return Err(Error::InvalidParameter {
                param: "grid".to_string(),
                reason: "Grid must be non-empty (ny > 0)".to_string(),
            });
        }
        let nz = grid[0][0].len();
        if nz == 0 {
            return Err(Error::InvalidParameter {
                param: "grid".to_string(),
                reason: "Grid must be non-empty (nz > 0)".to_string(),
            });
        }

        let safe_record_every = if record_every == 0 { 1 } else { record_every };

        let capacity = num_steps / safe_record_every + 1;
        let mut hopf_history = Vec::with_capacity(capacity);
        let mut energy_history = Vec::with_capacity(capacity);
        let mut time_vec = Vec::with_capacity(capacity);

        // Record initial state
        let q0 = self.compute_hopf_invariant(&grid, nx, ny, nz);
        let e0 = self.total_energy(&grid, nx, ny, nz);
        hopf_history.push(q0);
        energy_history.push(e0);
        time_vec.push(0.0);

        for step in 1..=num_steps {
            // Compute effective field for the entire grid before updating any spin
            let h_eff = self.compute_heff(&grid, nx, ny, nz);

            // Apply RK4 step to every site
            for ix in 0..nx {
                for iy in 0..ny {
                    for iz in 0..nz {
                        let m_new = self.llg_rk4_step(&grid[ix][iy][iz], &h_eff[ix][iy][iz])?;
                        grid[ix][iy][iz] = m_new.normalize();
                    }
                }
            }

            if step % safe_record_every == 0 {
                let q = self.compute_hopf_invariant(&grid, nx, ny, nz);
                let e = self.total_energy(&grid, nx, ny, nz);
                let t = step as f64 * self.config.dt;
                hopf_history.push(q);
                energy_history.push(e);
                time_vec.push(t);
            }
        }

        Ok(HopfionDynamicsResult {
            hopf_invariant_history: hopf_history,
            energy_history,
            time: time_vec,
        })
    }

    // ── Field computations ────────────────────────────────────────────────

    /// Compute the effective field H_eff = H_exchange + H_DMI + H_external
    /// for every lattice site.
    fn compute_heff(
        &self,
        grid: &[Vec<Vec<Vector3<f64>>>],
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Vec<Vec<Vec<Vector3<f64>>>> {
        let mut h_eff = vec![vec![vec![Vector3::zero(); nz]; ny]; nx];

        // ix/iy/iz are required as neighbour indices passed to exchange/DMI field
        // functions — clippy needless_range_loop is a false positive here.
        #[allow(clippy::needless_range_loop)]
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let h_ex = self.exchange_field(grid, ix, iy, iz, nx, ny, nz);
                    let h_dmi = self.dmi_field(grid, ix, iy, iz, nx, ny, nz);
                    h_eff[ix][iy][iz] = h_ex + h_dmi + self.config.h_ext;
                }
            }
        }
        h_eff
    }

    /// Exchange field via a 6-point Laplacian with periodic boundary conditions.
    ///
    /// ```text
    /// ∇²m ≈ (m[i+1] + m[i-1] + m[j+1] + m[j-1] + m[k+1] + m[k-1] − 6·m) / Δx²
    /// H_ex = (2A / (μ₀·Mₛ)) · ∇²m
    /// ```
    fn exchange_field(
        &self,
        grid: &[Vec<Vec<Vector3<f64>>>],
        ix: usize,
        iy: usize,
        iz: usize,
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Vector3<f64> {
        let dx2 = self.config.dx * self.config.dx;
        let prefactor = 2.0 * self.config.exchange_a / (Self::MU0 * self.config.ms * dx2);

        let m_center = grid[ix][iy][iz];

        // Periodic neighbour indices
        let ixp = (ix + 1) % nx;
        let ixm = (ix + nx - 1) % nx;
        let iyp = (iy + 1) % ny;
        let iym = (iy + ny - 1) % ny;
        let izp = (iz + 1) % nz;
        let izm = (iz + nz - 1) % nz;

        let laplacian = grid[ixp][iy][iz]
            + grid[ixm][iy][iz]
            + grid[ix][iyp][iz]
            + grid[ix][iym][iz]
            + grid[ix][iy][izp]
            + grid[ix][iy][izm]
            + m_center * (-6.0);

        laplacian * prefactor
    }

    /// Bulk DMI field using central-difference curl of the magnetization.
    ///
    /// ```text
    /// curl m = [∂mz/∂y − ∂my/∂z,  ∂mx/∂z − ∂mz/∂x,  ∂my/∂x − ∂mx/∂y]
    /// H_DMI  = (2D / (μ₀·Mₛ)) · curl m
    /// ```
    fn dmi_field(
        &self,
        grid: &[Vec<Vec<Vector3<f64>>>],
        ix: usize,
        iy: usize,
        iz: usize,
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Vector3<f64> {
        let two_dx = 2.0 * self.config.dx;
        let prefactor = 2.0 * self.config.dmi_d / (Self::MU0 * self.config.ms * two_dx);

        // Periodic neighbour indices
        let ixp = (ix + 1) % nx;
        let ixm = (ix + nx - 1) % nx;
        let iyp = (iy + 1) % ny;
        let iym = (iy + ny - 1) % ny;
        let izp = (iz + 1) % nz;
        let izm = (iz + nz - 1) % nz;

        let mxp = grid[ixp][iy][iz];
        let mxm = grid[ixm][iy][iz];
        let myp = grid[ix][iyp][iz];
        let mym = grid[ix][iym][iz];
        let mzp = grid[ix][iy][izp];
        let mzm = grid[ix][iy][izm];

        // Central differences: ∂f/∂x ≈ (f[i+1] − f[i-1]) / (2·dx)
        // Using prefactor = 2D / (μ₀·Mₛ·2dx) so multiply each difference directly
        let dmz_dy = myp.z - mym.z;
        let dmy_dz = mzp.y - mzm.y;
        let dmx_dz = mzp.x - mzm.x;
        let dmz_dx = mxp.z - mxm.z;
        let dmy_dx = mxp.y - mxm.y;
        let dmx_dy = myp.x - mym.x;

        Vector3::new(
            (dmz_dy - dmy_dz) * prefactor,
            (dmx_dz - dmz_dx) * prefactor,
            (dmy_dx - dmx_dy) * prefactor,
        )
    }

    // ── Time integration ──────────────────────────────────────────────────

    /// Fourth-order Runge-Kutta step for a single spin via the LLG equation.
    ///
    /// ```text
    /// dm/dt = −γ / (1 + α²) · [m × H_eff + α · m × (m × H_eff)]
    /// ```
    ///
    /// # Errors
    /// Returns `NumericalError` if the resulting vector contains NaN or ±∞.
    fn llg_rk4_step(&self, m: &Vector3<f64>, h_eff: &Vector3<f64>) -> Result<Vector3<f64>, Error> {
        let dt = self.config.dt;

        let k1 = self.llg_dm_dt(m, h_eff);
        let m1 = (*m + k1 * (dt * 0.5)).normalize();

        let k2 = self.llg_dm_dt(&m1, h_eff);
        let m2 = (*m + k2 * (dt * 0.5)).normalize();

        let k3 = self.llg_dm_dt(&m2, h_eff);
        let m3 = (*m + k3 * dt).normalize();

        let k4 = self.llg_dm_dt(&m3, h_eff);

        let m_new = *m + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);

        if !m_new.x.is_finite() || !m_new.y.is_finite() || !m_new.z.is_finite() {
            return Err(Error::NumericalError {
                description: "LLG RK4 step produced non-finite magnetization".to_string(),
            });
        }

        Ok(m_new)
    }

    /// Evaluate the LLG right-hand side dm/dt for given m and H_eff.
    #[inline]
    fn llg_dm_dt(&self, m: &Vector3<f64>, h_eff: &Vector3<f64>) -> Vector3<f64> {
        let alpha = self.config.alpha;
        let coeff = -self.gamma / (1.0 + alpha * alpha);

        let m_cross_h = m.cross(h_eff);
        let m_cross_m_cross_h = m.cross(&m_cross_h);

        (m_cross_h + m_cross_m_cross_h * alpha) * coeff
    }

    // ── Energy ────────────────────────────────────────────────────────────

    /// Compute the total magnetic energy of the grid \[J\].
    ///
    /// ```text
    /// E = −(μ₀ Mₛ / 2) Σ m·H_ex Δx³   (exchange, factor ½ avoids double-counting)
    ///   − μ₀ Mₛ Σ m·H_ext Δx³          (Zeeman)
    ///   + D Σ ε_DMI Δx³                  (DMI via curl m)
    /// ```
    pub fn total_energy(
        &self,
        grid: &[Vec<Vec<Vector3<f64>>>],
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> f64 {
        let dv = self.config.dx.powi(3);
        let mu0_ms = Self::MU0 * self.config.ms;

        let mut e_exchange = 0.0;
        let mut e_zeeman = 0.0;
        let mut e_dmi = 0.0;

        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let m = grid[ix][iy][iz];

                    // Exchange: −(μ₀ Mₛ / 2) m · H_ex
                    let h_ex = self.exchange_field(grid, ix, iy, iz, nx, ny, nz);
                    e_exchange += m.dot(&h_ex);

                    // Zeeman: −μ₀ Mₛ m · H_ext
                    e_zeeman += m.dot(&self.config.h_ext);

                    // DMI density contribution (D/μ₀Mₛ prefactor already in dmi_field)
                    // H_DMI = (2D / μ₀Mₛ) curl m  →  DMI energy density = D · m · curl_m
                    // We recover D · curl_m by multiplying H_DMI back:
                    // H_DMI · μ₀Mₛ / 2 = D · curl_m
                    // E_DMI = Σ D·curl_m·m·dV = Σ H_DMI · μ₀Mₛ/2 · m · dV
                    let h_dmi = self.dmi_field(grid, ix, iy, iz, nx, ny, nz);
                    e_dmi += m.dot(&h_dmi);
                }
            }
        }

        // Assemble with proper prefactors
        let exchange_energy = -0.5 * mu0_ms * e_exchange * dv;
        let zeeman_energy = -mu0_ms * e_zeeman * dv;
        // DMI: H_dmi already carries 2D/(μ₀Mₛ), so m·H_dmi·μ₀Mₛ/2·dV gives D·m·curl_m·dV
        let dmi_energy = 0.5 * mu0_ms * e_dmi * dv;

        exchange_energy + zeeman_energy + dmi_energy
    }

    // ── Hopf invariant ────────────────────────────────────────────────────

    /// Compute the Hopf invariant for a raw grid (Berry-connection method).
    ///
    /// For grids smaller than 4 in any dimension the invariant is ill-defined;
    /// this function returns 0.0 in that case (uniform / trivial topology).
    ///
    /// The formula used is
    ///
    /// ```text
    /// Q_H = (1 / 4π²) ∫ A · (∇ × A) dV
    /// ```
    ///
    /// where `A_i = m · (∂m/∂xⱼ × ∂m/∂xₖ)` for cyclic (i,j,k).
    pub fn compute_hopf_invariant(
        &self,
        grid: &[Vec<Vec<Vector3<f64>>>],
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> f64 {
        if nx < 4 || ny < 4 || nz < 4 {
            return 0.0;
        }

        let dx = self.config.dx;
        let dv = dx * dx * dx;
        let mut total = 0.0;

        // Interior points: skip outermost shell to safely use ± neighbour
        for ix in 1..nx - 1 {
            for iy in 1..ny - 1 {
                for iz in 1..nz - 1 {
                    let berry_a = Self::berry_connection_at_grid(grid, ix, iy, iz, nx, ny, nz, dx);

                    // Compute Berry connection at the six face-centre neighbours
                    let a_xp = Self::berry_connection_at_grid(grid, ix + 1, iy, iz, nx, ny, nz, dx);
                    let a_xm = Self::berry_connection_at_grid(grid, ix - 1, iy, iz, nx, ny, nz, dx);
                    let a_yp = Self::berry_connection_at_grid(grid, ix, iy + 1, iz, nx, ny, nz, dx);
                    let a_ym = Self::berry_connection_at_grid(grid, ix, iy - 1, iz, nx, ny, nz, dx);
                    let a_zp = Self::berry_connection_at_grid(grid, ix, iy, iz + 1, nx, ny, nz, dx);
                    let a_zm = Self::berry_connection_at_grid(grid, ix, iy, iz - 1, nx, ny, nz, dx);

                    let inv_2dx = 0.5 / dx;
                    let curl_a = Vector3::new(
                        (a_yp.z - a_ym.z - a_zp.y + a_zm.y) * inv_2dx,
                        (a_zp.x - a_zm.x - a_xp.z + a_xm.z) * inv_2dx,
                        (a_xp.y - a_xm.y - a_yp.x + a_ym.x) * inv_2dx,
                    );

                    total += berry_a.dot(&curl_a) * dv;
                }
            }
        }

        total / (4.0 * PI * PI)
    }

    /// Compute the Berry connection vector A at a specific grid site.
    ///
    /// Uses central differences for the spatial derivatives of m with clamped
    /// indices (no wrap) so that boundary calls are safe.
    fn berry_connection_at_grid(
        grid: &[Vec<Vec<Vector3<f64>>>],
        ix: usize,
        iy: usize,
        iz: usize,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
    ) -> Vector3<f64> {
        let ix = ix.min(nx - 1);
        let iy = iy.min(ny - 1);
        let iz = iz.min(nz - 1);

        let ix_lo = if ix > 0 { ix - 1 } else { ix };
        let ix_hi = if ix + 1 < nx { ix + 1 } else { ix };
        let iy_lo = if iy > 0 { iy - 1 } else { iy };
        let iy_hi = if iy + 1 < ny { iy + 1 } else { iy };
        let iz_lo = if iz > 0 { iz - 1 } else { iz };
        let iz_hi = if iz + 1 < nz { iz + 1 } else { iz };

        let m = grid[ix][iy][iz];

        let dx_eff = (ix_hi - ix_lo) as f64 * dx;
        let dy_eff = (iy_hi - iy_lo) as f64 * dx;
        let dz_eff = (iz_hi - iz_lo) as f64 * dx;

        let dm_dx = (grid[ix_hi][iy][iz] - grid[ix_lo][iy][iz]) * (1.0 / dx_eff.max(dx));
        let dm_dy = (grid[ix][iy_hi][iz] - grid[ix][iy_lo][iz]) * (1.0 / dy_eff.max(dx));
        let dm_dz = (grid[ix][iy][iz_hi] - grid[ix][iy][iz_lo]) * (1.0 / dz_eff.max(dx));

        Vector3::new(
            m.dot(&dm_dy.cross(&dm_dz)),
            m.dot(&dm_dz.cross(&dm_dx)),
            m.dot(&dm_dx.cross(&dm_dy)),
        )
    }
}

// ─── Result ───────────────────────────────────────────────────────────────────

/// Result produced by [`HopfionDynamicsSolver::run`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HopfionDynamicsResult {
    /// Hopf invariant Q_H sampled at each recording step
    pub hopf_invariant_history: Vec<f64>,
    /// Total magnetic energy \[J\] sampled at each recording step
    pub energy_history: Vec<f64>,
    /// Simulation time \[s\] corresponding to each sample
    pub time: Vec<f64>,
}

// ─── Helper (module-private) ──────────────────────────────────────────────────

/// Build a uniform `[nx][ny][nz]` grid pointing in +z direction.
///
/// Useful in tests and as a trivial initial condition.
pub fn uniform_grid(nx: usize, ny: usize, nz: usize) -> Vec<Vec<Vec<Vector3<f64>>>> {
    vec![vec![vec![Vector3::unit_z(); nz]; ny]; nx]
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helpers ──────────────────────────────────────────────────────────────

    /// Build a default physically reasonable config for tests.
    fn default_config() -> HopfionDynamicsConfig {
        HopfionDynamicsConfig::new(
            1.0e-11,                     // A  [J/m]
            1.0e-3,                      // D  [J/m²]
            0.1,                         // α
            Vector3::new(0.0, 0.0, 0.1), // H_ext [T]
            8.0e5,                       // Mₛ [A/m]
            1.0e-14,                     // dt [s]
            5.0e-10,                     // dx [m]
        )
    }

    fn small_solver() -> HopfionDynamicsSolver {
        HopfionDynamicsSolver::new(default_config()).expect("default config must be valid")
    }

    // Test 1 ─────────────────────────────────────────────────────────────────
    #[test]
    fn test_config_validation_valid() {
        let cfg = default_config();
        assert!(
            cfg.validate().is_ok(),
            "Valid config should pass validation"
        );
    }

    // Test 2 ─────────────────────────────────────────────────────────────────
    #[test]
    fn test_config_validation_invalid_ms() {
        let mut cfg = default_config();
        cfg.ms = -1.0;
        assert!(
            cfg.validate().is_err(),
            "Negative ms should fail validation"
        );

        cfg.ms = 0.0;
        assert!(cfg.validate().is_err(), "Zero ms should fail validation");

        cfg.ms = f64::INFINITY;
        assert!(
            cfg.validate().is_err(),
            "Infinite ms should fail validation"
        );
    }

    // Test 3 ─────────────────────────────────────────────────────────────────
    #[test]
    fn test_uniform_exchange_field_is_zero() {
        // A spatially uniform spin texture has ∇²m = 0 → H_ex = 0
        let solver = small_solver();
        let nx = 6;
        let ny = 6;
        let nz = 6;
        let grid = uniform_grid(nx, ny, nz);

        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let h_ex = solver.exchange_field(&grid, ix, iy, iz, nx, ny, nz);
                    assert!(
                        h_ex.magnitude() < 1e-30,
                        "Exchange field must vanish on uniform grid, got {:?}",
                        h_ex
                    );
                }
            }
        }
    }

    // Test 4 ─────────────────────────────────────────────────────────────────
    #[test]
    fn test_known_gradient_exchange_field() {
        // Place a single spin pointing along +x in a sea of +z spins on a 5×5×5
        // periodic grid. The Laplacian at (2,2,2) should equal
        //   (m_xp + m_xm + m_yp + m_ym + m_zp + m_zm − 6·m_center) / dx²
        // All neighbours are +z, center is +x.
        let mut cfg = default_config();
        cfg.exchange_a = 1.0e-11;
        let solver = HopfionDynamicsSolver::new(cfg.clone()).expect("valid config");

        let nx = 5;
        let ny = 5;
        let nz = 5;
        let mut grid = uniform_grid(nx, ny, nz);
        grid[2][2][2] = Vector3::unit_x();

        let h_ex = solver.exchange_field(&grid, 2, 2, 2, nx, ny, nz);

        // ∇²m_x = (0+0+0+0+0+0 − 6·1) / dx² = −6/dx²
        let expected_mx =
            -6.0 / (cfg.dx * cfg.dx) * 2.0 * cfg.exchange_a / (HopfionDynamicsSolver::MU0 * cfg.ms);
        assert!(
            (h_ex.x - expected_mx).abs() < expected_mx.abs() * 1e-10,
            "H_ex.x = {} expected = {}",
            h_ex.x,
            expected_mx
        );
        // y and z Laplacian components: neighbouring +z gives (∇²m)_z contribution
        // For the z-component: neighbours are +z, center is 0 → Laplacian = 6
        let expected_mz =
            6.0 / (cfg.dx * cfg.dx) * 2.0 * cfg.exchange_a / (HopfionDynamicsSolver::MU0 * cfg.ms);
        assert!(
            (h_ex.z - expected_mz).abs() < expected_mz.abs() * 1e-10,
            "H_ex.z = {} expected = {}",
            h_ex.z,
            expected_mz
        );
    }

    // Test 5 ─────────────────────────────────────────────────────────────────
    #[test]
    fn test_zeeman_only_heff() {
        // When A=0 and D=0, H_eff must equal h_ext for every site
        let h_target = Vector3::new(0.3, -0.1, 0.5);
        let cfg = HopfionDynamicsConfig::new(0.0, 0.0, 0.1, h_target, 8.0e5, 1.0e-14, 5.0e-10);
        let solver = HopfionDynamicsSolver::new(cfg).expect("valid config");

        let nx = 4;
        let ny = 4;
        let nz = 4;
        let grid = uniform_grid(nx, ny, nz);
        let h_eff = solver.compute_heff(&grid, nx, ny, nz);

        for row in &h_eff {
            for col in row {
                for h in col {
                    let diff = (*h - h_target).magnitude();
                    assert!(
                        diff < 1e-20,
                        "H_eff should equal h_ext when A=D=0, diff={}",
                        diff
                    );
                }
            }
        }
    }

    // Test 6 ─────────────────────────────────────────────────────────────────
    #[test]
    fn test_llg_rk4_magnitude_preserved() {
        // After one RK4 step and normalisation the spin must still be a unit vector.
        let solver = small_solver();
        let m = Vector3::new(
            1.0_f64 / 3.0_f64.sqrt(),
            1.0_f64 / 3.0_f64.sqrt(),
            1.0_f64 / 3.0_f64.sqrt(),
        );
        let h = Vector3::new(0.0, 0.0, 0.05);
        let m_new = solver.llg_rk4_step(&m, &h).expect("RK4 step must succeed");
        let m_norm = m_new.normalize();
        assert!(
            (m_norm.magnitude() - 1.0).abs() < 1e-10,
            "Normalised spin magnitude = {}",
            m_norm.magnitude()
        );
    }

    // Test 7 ─────────────────────────────────────────────────────────────────
    #[test]
    fn test_zero_damping_energy_conservation() {
        // With α = 0 the LLG is a pure precession; energy must be conserved.
        let cfg = HopfionDynamicsConfig::new(
            1.0e-11,
            0.0,
            0.0, // no damping, no DMI
            Vector3::new(0.0, 0.0, 0.05),
            8.0e5,
            5.0e-15, // small dt for accuracy
            5.0e-10,
        );
        let solver = HopfionDynamicsSolver::new(cfg).expect("valid config");

        let nx = 4;
        let ny = 4;
        let nz = 4;
        let grid = uniform_grid(nx, ny, nz);

        let result = solver.run(grid, 20, 10).expect("run must succeed");

        let e0 = result.energy_history[0];
        for &e in &result.energy_history {
            // Allow 1% relative deviation due to discretisation
            let rel = (e - e0).abs() / (e0.abs() + 1e-30);
            assert!(
                rel < 0.01,
                "Energy conservation violated: e0={:.6e}, e={:.6e}, rel={:.4e}",
                e0,
                e,
                rel
            );
        }
    }

    // Test 8 ─────────────────────────────────────────────────────────────────
    #[test]
    fn test_uniform_state_hopf_invariant_zero() {
        // A uniform magnetization state is topologically trivial: Q_H = 0
        let solver = small_solver();
        let nx = 6;
        let ny = 6;
        let nz = 6;
        let grid = uniform_grid(nx, ny, nz);
        let q = solver.compute_hopf_invariant(&grid, nx, ny, nz);
        assert!(
            q.abs() < 1e-10,
            "Uniform state should have Q_H = 0, got {}",
            q
        );
    }

    // Test 9 ─────────────────────────────────────────────────────────────────
    #[test]
    fn test_run_history_length_correct() {
        // Ensure that running N steps with record_every=R produces ceil(N/R)+1
        // (initial + recorded) entries in each history.
        let solver = small_solver();
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let grid = uniform_grid(nx, ny, nz);

        let num_steps = 10;
        let record_every = 5;
        let result = solver
            .run(grid, num_steps, record_every)
            .expect("run must succeed");

        // Initial recording + recordings at steps 5 and 10
        let expected = 1 + num_steps / record_every;
        assert_eq!(
            result.hopf_invariant_history.len(),
            expected,
            "Hopf history length mismatch"
        );
        assert_eq!(
            result.energy_history.len(),
            expected,
            "Energy history length mismatch"
        );
        assert_eq!(result.time.len(), expected, "Time history length mismatch");
    }

    // Test 10 ────────────────────────────────────────────────────────────────
    #[test]
    fn test_energy_is_finite() {
        let solver = small_solver();
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let grid = uniform_grid(nx, ny, nz);
        let e = solver.total_energy(&grid, nx, ny, nz);
        assert!(e.is_finite(), "Total energy must be finite, got {}", e);
    }

    // Test 11 ────────────────────────────────────────────────────────────────
    #[test]
    fn test_periodic_bc_index_wrapping() {
        // The exchange field at corner (0,0,0) uses wrap-around neighbours
        // (nx-1,0,0) etc. On a uniform grid the result must be zero.
        let solver = small_solver();
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let grid = uniform_grid(nx, ny, nz);

        // Corner sites
        for &ix in &[0, nx - 1] {
            for &iy in &[0, ny - 1] {
                for &iz in &[0, nz - 1] {
                    let h = solver.exchange_field(&grid, ix, iy, iz, nx, ny, nz);
                    assert!(
                        h.magnitude() < 1e-30,
                        "Exchange field at corner ({},{},{}) must be 0 on uniform grid",
                        ix,
                        iy,
                        iz
                    );
                }
            }
        }
    }

    // Test 12 ────────────────────────────────────────────────────────────────
    #[test]
    fn test_history_lengths_consistent() {
        let solver = small_solver();
        let grid = uniform_grid(4, 4, 4);
        let result = solver.run(grid, 8, 4).expect("run must succeed");

        let n = result.hopf_invariant_history.len();
        assert_eq!(result.energy_history.len(), n);
        assert_eq!(result.time.len(), n);
    }

    // Test 13 ────────────────────────────────────────────────────────────────
    #[test]
    #[allow(clippy::needless_range_loop)]
    fn test_external_field_alignment() {
        // With large damping and a strong z-field, all spins should relax to +z.
        let h_z = Vector3::new(0.0, 0.0, 1.0e3); // very large field
        let cfg = HopfionDynamicsConfig::new(
            0.0, 0.0, 0.99, // near-critical damping
            h_z, 8.0e5, 1.0e-13, // larger dt to speed relaxation
            5.0e-10,
        );
        let solver = HopfionDynamicsSolver::new(cfg).expect("valid config");

        // Start with spins tilted away from z
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let mut grid = vec![vec![vec![Vector3::zero(); nz]; ny]; nx];
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    grid[ix][iy][iz] = Vector3::new(0.5, 0.5, 0.5_f64.sqrt()).normalize();
                }
            }
        }

        let result = solver.run(grid, 500, 500).expect("run must succeed");

        // Last energy should be lower (more negative) than first → relaxation
        let e_init = result.energy_history[0];
        let e_final = *result.energy_history.last().expect("non-empty");
        assert!(
            e_final <= e_init + 1e-30,
            "Energy should decrease (or stay same) during damped relaxation: \
             e_init={:.6e}, e_final={:.6e}",
            e_init,
            e_final
        );
    }

    // Test 14 ────────────────────────────────────────────────────────────────
    #[test]
    fn test_time_step_spacing() {
        // The recorded times must be spaced by record_every * dt
        let cfg = default_config();
        let dt = cfg.dt;
        let solver = HopfionDynamicsSolver::new(cfg).expect("valid config");

        let grid = uniform_grid(4, 4, 4);
        let num_steps = 9;
        let record_every = 3;
        let result = solver
            .run(grid, num_steps, record_every)
            .expect("run must succeed");

        // times should be [0, 3dt, 6dt, 9dt]
        for (i, &t) in result.time.iter().enumerate() {
            let expected = i as f64 * record_every as f64 * dt;
            assert!(
                (t - expected).abs() < expected.abs() * 1e-10 + 1e-30,
                "Time[{}] = {:.6e}, expected {:.6e}",
                i,
                t,
                expected
            );
        }
    }

    // Test 15 ────────────────────────────────────────────────────────────────
    #[test]
    fn test_small_grid_no_panic() {
        // 2×2×2 grid: Hopf invariant is short-circuited (returns 0), but the
        // dynamics loop itself must not panic.
        let solver = small_solver();
        let grid = uniform_grid(2, 2, 2);
        let result = solver.run(grid, 2, 1).expect("2×2×2 run must succeed");

        // Hopf invariant skipped for small grid → should be 0
        for &q in &result.hopf_invariant_history {
            assert_eq!(q, 0.0, "Q_H on 2×2×2 must be 0 (small-grid bypass)");
        }
        // Energy must still be finite
        for &e in &result.energy_history {
            assert!(e.is_finite(), "Energy must be finite");
        }
    }
}
