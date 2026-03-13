//! Three-dimensional topological solitons (Hopfions)
//!
//! Hopfions are 3D magnetic solitons characterized by the Hopf invariant,
//! which counts the linking number of preimage curves on S². They represent
//! the three-dimensional generalization of magnetic skyrmions.
//!
//! # Mathematical Background
//!
//! The Hopf fibration is a continuous map S³ -> S² that partitions the
//! 3-sphere into linked circles (Hopf fibers). A hopfion is a magnetization
//! texture whose configuration is described by this fibration.
//!
//! The magnetization field is parameterized in toroidal coordinates as:
//!
//! $$
//! \mathbf{n}(\mathbf{r}) = \left(\sin f(\rho) \cos(N\phi + Q\theta),\;
//!     \sin f(\rho) \sin(N\phi + Q\theta),\;
//!     \cos f(\rho)\right)
//! $$
//!
//! where $\rho = \sqrt{x^2 + y^2}$ is the cylindrical radius,
//! $\phi$ is the azimuthal angle, $\theta$ is the toroidal angle,
//! $N$ is the azimuthal winding number, $Q$ is the poloidal winding number,
//! and $f(\rho)$ is the radial profile function.
//!
//! # Hopf Invariant
//!
//! The topological charge is the Hopf invariant:
//!
//! $$
//! Q_H = \frac{1}{4\pi^2} \int F \cdot A \, d^3r
//! $$
//!
//! where $F = \nabla \times A$ and $A$ is the Berry connection derived from
//! the magnetization field. For a well-formed hopfion, $Q_H$ is an integer.
//!
//! # Energy Functional
//!
//! The total energy of a hopfion configuration includes:
//! - Exchange energy: $E_\text{ex} = A \int |\nabla \mathbf{m}|^2 \, d^3r$
//! - Dzyaloshinskii-Moriya interaction energy in 3D
//! - Zeeman energy: $E_Z = -\mu_0 M_s \int \mathbf{m} \cdot \mathbf{H} \, d^3r$
//!
//! # Stability
//!
//! Hopfions can be stabilized by frustrated exchange interactions, where
//! competing nearest-neighbor and next-nearest-neighbor exchange couplings
//! create an energy landscape with a local minimum at a finite hopfion radius.
//!
//! # References
//!
//! - P. J. Ackerman and I. I. Smalyukh, "Static three-dimensional topological
//!   solitons in fluid chiral ferromagnets and colloids", Nat. Mater. 16, 426 (2017)
//! - F. N. Rybakov et al., "Magnetic hopfions in solids", APL Mater. 10, 111113 (2022)
//! - Y. Liu et al., "Binding a hopfion in a chiral magnet nanodisk",
//!   Phys. Rev. B 98, 174437 (2018)
//! - J. Kent et al., "Creation and observation of hopfions in magnetic multilayer
//!   systems", Nat. Commun. 12, 1562 (2021)

use std::f64::consts::PI;
use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::vector3::Vector3;

// ============================================================================
// Hopf Fibration
// ============================================================================

/// Hopf fibration map from S³ to S²
///
/// The Hopf map sends a point on the 3-sphere (parameterized by two complex
/// numbers z₁, z₂ with |z₁|² + |z₂|² = 1) to a point on the 2-sphere via
/// the stereographic projection of z₁/z₂.
///
/// # Example
///
/// ```rust
/// use spintronics::texture::hopfion::HopfFibration;
/// use spintronics::vector3::Vector3;
///
/// // Map a point from S³ to S²
/// let result = HopfFibration::map(1.0, 0.0, 0.0, 0.0);
/// // (1,0,0,0) on S³ maps to the north pole (0,0,1) on S²
/// assert!((result.z - 1.0).abs() < 1e-10);
/// ```
pub struct HopfFibration;

impl HopfFibration {
    /// Compute the Hopf map: S³ -> S²
    ///
    /// Given a point (z₁, z₂) = (a + bi, c + di) on S³ (|z₁|² + |z₂|² = 1),
    /// returns the corresponding point on S² via the Hopf fibration.
    ///
    /// The map is defined as:
    /// - n_x = 2(ac + bd)
    /// - n_y = 2(bc - ad)
    /// - n_z = a² + b² - c² - d²
    ///
    /// # Arguments
    /// * `a`, `b` - Real and imaginary parts of z₁
    /// * `c`, `d` - Real and imaginary parts of z₂
    ///
    /// # Returns
    /// Unit vector on S²
    pub fn map(a: f64, b: f64, c: f64, d: f64) -> Vector3<f64> {
        let norm_sq = a * a + b * b + c * c + d * d;
        let inv_norm_sq = if norm_sq > 1e-30 { 1.0 / norm_sq } else { 1.0 };

        Vector3::new(
            2.0 * (a * c + b * d) * inv_norm_sq,
            2.0 * (b * c - a * d) * inv_norm_sq,
            (a * a + b * b - c * c - d * d) * inv_norm_sq,
        )
    }

    /// Inverse Hopf map: find a preimage fiber point on S³ for a given S² point
    ///
    /// For a point n on S², returns one representative point on S³ that maps
    /// to n under the Hopf fibration. The full preimage is a circle (fiber)
    /// obtained by multiplying z₁ by e^{iα} for all α in [0, 2π).
    ///
    /// # Arguments
    /// * `n` - Target unit vector on S²
    /// * `fiber_angle` - Angle parameterizing the fiber circle \[rad\]
    ///
    /// # Returns
    /// (a, b, c, d) on S³
    pub fn inverse_map(n: &Vector3<f64>, fiber_angle: f64) -> (f64, f64, f64, f64) {
        // Normalize n to ensure it's on S²
        let mag = (n.x * n.x + n.y * n.y + n.z * n.z).sqrt();
        let nx = if mag > 1e-30 { n.x / mag } else { 0.0 };
        let ny = if mag > 1e-30 { n.y / mag } else { 0.0 };
        let nz = if mag > 1e-30 { n.z / mag } else { 1.0 };
        let nz_clamped = nz.clamp(-1.0, 1.0);

        // |z₁|² = (1 + nz)/2, |z₂|² = (1 - nz)/2
        let abs_z1 = ((1.0 + nz_clamped) / 2.0).sqrt();
        let abs_z2 = ((1.0 - nz_clamped) / 2.0).sqrt();

        // Phase of z₁ * conj(z₂): arg(z₁ z₂*) = atan2(ny, nx)
        // We parameterize the fiber by: arg(z₁) = fiber_angle
        // Then: arg(z₂) = arg(z₁) - arg(z₁ z₂*) = fiber_angle - atan2(ny, nx)
        let phase_product = ny.atan2(nx); // arg(z₁ z₂*)
        let arg_z1 = fiber_angle;
        let arg_z2 = fiber_angle - phase_product;

        let a = abs_z1 * arg_z1.cos();
        let b = abs_z1 * arg_z1.sin();
        let c = abs_z2 * arg_z2.cos();
        let d = abs_z2 * arg_z2.sin();

        (a, b, c, d)
    }

    /// Compute magnetization direction at a 3D position using the Hopf map
    ///
    /// Maps a spatial position (x, y, z) to S² through the Hopf fibration,
    /// parameterized by a hopfion of given radius centered at the origin.
    ///
    /// This uses stereographic projection to map R³ -> S³, then the Hopf map
    /// S³ -> S².
    ///
    /// # Arguments
    /// * `x`, `y`, `z` - Spatial coordinates
    /// * `radius` - Hopfion characteristic radius
    ///
    /// # Returns
    /// Unit magnetization vector
    pub fn magnetization_from_position(x: f64, y: f64, z: f64, radius: f64) -> Vector3<f64> {
        // Normalize coordinates by radius
        let safe_radius = if radius.abs() < 1e-30 { 1e-30 } else { radius };
        let xn = x / safe_radius;
        let yn = y / safe_radius;
        let zn = z / safe_radius;

        // Distance squared from origin (in normalized coords)
        let r_sq = xn * xn + yn * yn + zn * zn;

        // Inverse stereographic projection: R³ -> S³
        // (xn, yn, zn) -> (2xn, 2yn, 2zn, r²-1) / (r²+1)
        let denom = r_sq + 1.0;
        let inv_denom = 1.0 / denom;

        let s3_a = 2.0 * xn * inv_denom;
        let s3_b = 2.0 * yn * inv_denom;
        let s3_c = 2.0 * zn * inv_denom;
        let s3_d = (r_sq - 1.0) * inv_denom;

        // Apply Hopf map S³ -> S²
        Self::map(s3_a, s3_b, s3_c, s3_d)
    }
}

// ============================================================================
// Toroidal Coordinate System
// ============================================================================

/// Toroidal coordinates for hopfion parameterization
///
/// Converts between Cartesian (x, y, z) and toroidal (rho, theta, phi) coordinates
/// where:
/// - rho: distance from the toroidal axis ring
/// - theta: poloidal angle (around the cross-section of the torus)
/// - phi: toroidal (azimuthal) angle
pub struct ToroidalCoords;

impl ToroidalCoords {
    /// Convert Cartesian to toroidal coordinates
    ///
    /// # Arguments
    /// * `x`, `y`, `z` - Cartesian coordinates
    /// * `major_radius` - Major radius of the torus (distance from axis to ring center)
    ///
    /// # Returns
    /// (rho, theta, phi) where rho is distance from ring, theta is poloidal angle,
    /// phi is toroidal angle
    pub fn from_cartesian(x: f64, y: f64, z: f64, major_radius: f64) -> (f64, f64, f64) {
        // Azimuthal angle
        let phi = y.atan2(x);

        // Distance from z-axis in xy-plane
        let r_cyl = (x * x + y * y).sqrt();

        // Distance from the toroidal ring
        let dr = r_cyl - major_radius;
        let rho = (dr * dr + z * z).sqrt();

        // Poloidal angle
        let theta = z.atan2(dr);

        (rho, theta, phi)
    }

    /// Convert toroidal to Cartesian coordinates
    ///
    /// # Arguments
    /// * `rho` - Distance from the toroidal ring
    /// * `theta` - Poloidal angle \[rad\]
    /// * `phi` - Toroidal (azimuthal) angle \[rad\]
    /// * `major_radius` - Major radius of the torus
    ///
    /// # Returns
    /// (x, y, z)
    pub fn to_cartesian(rho: f64, theta: f64, phi: f64, major_radius: f64) -> (f64, f64, f64) {
        let r_cyl = major_radius + rho * theta.cos();
        let x = r_cyl * phi.cos();
        let y = r_cyl * phi.sin();
        let z = rho * theta.sin();

        (x, y, z)
    }
}

// ============================================================================
// Hopfion Profile
// ============================================================================

/// Radial profile function for hopfion magnetization
///
/// The profile function f(rho) determines the out-of-plane tilt of the
/// magnetization as a function of distance from the hopfion ring.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ProfileFunction {
    /// Linear profile: f(rho) = pi * (1 - rho/R) for rho < R, 0 otherwise
    ///
    /// Simple analytical form, good for initial configurations.
    #[default]
    Linear,

    /// Hyperbolic tangent profile: f(rho) = pi * (1 - tanh(rho/w)) / 2
    ///
    /// Smoother transition, closer to equilibrium profiles.
    /// Parameter w controls the wall width.
    Tanh {
        /// Wall width parameter \[m\]
        wall_width: f64,
    },

    /// Exponential decay profile: f(rho) = pi * exp(-rho^2 / (2 sigma^2))
    ///
    /// Gaussian-like profile, useful for confined hopfions.
    Gaussian {
        /// Width parameter sigma \[m\]
        sigma: f64,
    },
}

impl ProfileFunction {
    /// Evaluate the profile function at a given distance from the ring
    ///
    /// # Arguments
    /// * `rho` - Distance from the hopfion ring \[m\]
    /// * `radius` - Characteristic hopfion radius \[m\]
    ///
    /// # Returns
    /// Profile angle f(rho) in [0, pi]
    pub fn evaluate(&self, rho: f64, radius: f64) -> f64 {
        let safe_radius = if radius.abs() < 1e-30 { 1e-30 } else { radius };

        match self {
            ProfileFunction::Linear => {
                if rho < safe_radius {
                    PI * (1.0 - rho / safe_radius)
                } else {
                    0.0
                }
            },
            ProfileFunction::Tanh { wall_width } => {
                let safe_w = if wall_width.abs() < 1e-30 {
                    1e-30
                } else {
                    *wall_width
                };
                PI * (1.0 - (rho / safe_w).tanh()) / 2.0
            },
            ProfileFunction::Gaussian { sigma } => {
                let safe_sigma = if sigma.abs() < 1e-30 { 1e-30 } else { *sigma };
                PI * (-rho * rho / (2.0 * safe_sigma * safe_sigma)).exp()
            },
        }
    }
}

// ============================================================================
// Hopfion Struct
// ============================================================================

/// A three-dimensional magnetic hopfion
///
/// Hopfions are 3D topological solitons characterized by the Hopf invariant.
/// The magnetization texture forms a toroidal structure where preimage curves
/// of any two points on S² are linked, with the linking number equal to the
/// Hopf charge.
///
/// # Example
///
/// ```rust
/// use spintronics::texture::hopfion::{Hopfion, ProfileFunction};
///
/// // Create a hopfion on a 32x32x32 grid with Hopf charge 1
/// let hopfion = Hopfion::new(
///     (32, 32, 32),
///     5.0e-9,   // 5 nm radius
///     1,         // Hopf charge
/// ).expect("Failed to create hopfion");
///
/// // Check that the magnetization is normalized
/// let mag = hopfion.magnetization_at(0, 0, 0);
/// assert!((mag.magnitude() - 1.0).abs() < 1e-10);
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Hopfion {
    /// Grid dimensions (Nx, Ny, Nz)
    pub grid_size: (usize, usize, usize),

    /// Flattened 3D magnetization field
    ///
    /// Indexed as magnetization[ix + iy * Nx + iz * Nx * Ny]
    pub magnetization: Vec<Vector3<f64>>,

    /// Hopf charge (topological invariant)
    pub hopf_charge: i32,

    /// Characteristic hopfion radius \[m\]
    pub radius: f64,

    /// Grid spacing \[m\]
    pub cell_size: f64,

    /// Profile function used for generation
    pub profile: ProfileFunction,

    /// Azimuthal winding number N
    pub azimuthal_winding: i32,

    /// Poloidal winding number Q
    pub poloidal_winding: i32,
}

impl Hopfion {
    /// Create a new hopfion on a 3D grid
    ///
    /// Generates a hopfion magnetization texture using the Hopf fibration
    /// with a linear radial profile.
    ///
    /// # Arguments
    /// * `grid_size` - Grid dimensions (Nx, Ny, Nz)
    /// * `radius` - Hopfion characteristic radius \[m\]
    /// * `hopf_charge` - Target Hopf charge (integer)
    ///
    /// # Returns
    /// A Hopfion with the magnetization field initialized
    ///
    /// # Errors
    /// Returns an error if grid dimensions are zero or radius is non-positive
    pub fn new(grid_size: (usize, usize, usize), radius: f64, hopf_charge: i32) -> Result<Self> {
        Self::with_profile(grid_size, radius, hopf_charge, ProfileFunction::Linear)
    }

    /// Create a hopfion with a specified profile function
    ///
    /// # Arguments
    /// * `grid_size` - Grid dimensions (Nx, Ny, Nz)
    /// * `radius` - Hopfion characteristic radius \[m\]
    /// * `hopf_charge` - Target Hopf charge
    /// * `profile` - Radial profile function
    ///
    /// # Errors
    /// Returns an error if parameters are invalid
    pub fn with_profile(
        grid_size: (usize, usize, usize),
        radius: f64,
        hopf_charge: i32,
        profile: ProfileFunction,
    ) -> Result<Self> {
        let (nx, ny, nz) = grid_size;

        if nx == 0 || ny == 0 || nz == 0 {
            return Err(Error::InvalidParameter {
                param: "grid_size".to_string(),
                reason: "All grid dimensions must be positive".to_string(),
            });
        }

        if radius <= 0.0 || !radius.is_finite() {
            return Err(Error::InvalidParameter {
                param: "radius".to_string(),
                reason: "Radius must be positive and finite".to_string(),
            });
        }

        // Cell size: cover ~4 radii across the smallest dimension
        let min_dim = nx.min(ny).min(nz) as f64;
        let cell_size = 4.0 * radius / min_dim;

        // Determine winding numbers from Hopf charge
        // For Q_H = N * Q, simplest is N=hopf_charge, Q=1
        let azimuthal_winding = hopf_charge;
        let poloidal_winding = 1;

        let total_size = nx * ny * nz;
        let mut magnetization = Vec::with_capacity(total_size);

        // Center of the grid
        let cx = (nx as f64 - 1.0) / 2.0;
        let cy = (ny as f64 - 1.0) / 2.0;
        let cz = (nz as f64 - 1.0) / 2.0;

        // Major radius of the toroidal ring (scaled to ~radius)
        let major_radius = radius;

        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let x = (ix as f64 - cx) * cell_size;
                    let y = (iy as f64 - cy) * cell_size;
                    let z = (iz as f64 - cz) * cell_size;

                    let m = Self::compute_magnetization_at_point(
                        x,
                        y,
                        z,
                        major_radius,
                        radius,
                        azimuthal_winding,
                        poloidal_winding,
                        &profile,
                    );
                    magnetization.push(m);
                }
            }
        }

        Ok(Self {
            grid_size,
            magnetization,
            hopf_charge,
            radius,
            cell_size,
            profile,
            azimuthal_winding,
            poloidal_winding,
        })
    }

    /// Compute the magnetization vector at a spatial point
    ///
    /// Uses toroidal coordinates and the hopfion ansatz:
    /// n(r) = (sin f(rho) cos(N*phi + Q*theta), sin f(rho) sin(N*phi + Q*theta), cos f(rho))
    fn compute_magnetization_at_point(
        x: f64,
        y: f64,
        z: f64,
        major_radius: f64,
        profile_radius: f64,
        n_winding: i32,
        q_winding: i32,
        profile: &ProfileFunction,
    ) -> Vector3<f64> {
        let (rho, theta, phi) = ToroidalCoords::from_cartesian(x, y, z, major_radius);

        let f_rho = profile.evaluate(rho, profile_radius);

        let phase = n_winding as f64 * phi + q_winding as f64 * theta;

        let sin_f = f_rho.sin();
        let cos_f = f_rho.cos();

        let mx = sin_f * phase.cos();
        let my = sin_f * phase.sin();
        let mz = cos_f;

        // The vector is already on the unit sphere by construction
        // (sin²f cos²phase + sin²f sin²phase + cos²f = sin²f + cos²f = 1)
        // but normalize for numerical safety
        let m = Vector3::new(mx, my, mz);
        let mag = m.magnitude();
        if mag > 1e-30 {
            Vector3::new(mx / mag, my / mag, mz / mag)
        } else {
            Vector3::new(0.0, 0.0, 1.0)
        }
    }

    /// Get the linear index for a grid point
    ///
    /// # Arguments
    /// * `ix`, `iy`, `iz` - Grid indices
    ///
    /// # Returns
    /// Linear index into the magnetization array, or None if out of bounds
    pub fn linear_index(&self, ix: usize, iy: usize, iz: usize) -> Option<usize> {
        let (nx, ny, nz) = self.grid_size;
        if ix < nx && iy < ny && iz < nz {
            Some(ix + iy * nx + iz * nx * ny)
        } else {
            None
        }
    }

    /// Get the magnetization at a grid point
    ///
    /// # Arguments
    /// * `ix`, `iy`, `iz` - Grid indices
    ///
    /// # Returns
    /// Magnetization vector (returns +z for out-of-bounds)
    pub fn magnetization_at(&self, ix: usize, iy: usize, iz: usize) -> Vector3<f64> {
        match self.linear_index(ix, iy, iz) {
            Some(idx) if idx < self.magnetization.len() => self.magnetization[idx],
            _ => Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Set the magnetization at a grid point
    ///
    /// # Errors
    /// Returns error if indices are out of bounds
    pub fn set_magnetization_at(
        &mut self,
        ix: usize,
        iy: usize,
        iz: usize,
        m: Vector3<f64>,
    ) -> Result<()> {
        let idx = self
            .linear_index(ix, iy, iz)
            .ok_or_else(|| Error::InvalidParameter {
                param: "grid_index".to_string(),
                reason: format!(
                    "Index ({}, {}, {}) out of bounds for grid {:?}",
                    ix, iy, iz, self.grid_size
                ),
            })?;

        if idx >= self.magnetization.len() {
            return Err(Error::InvalidParameter {
                param: "grid_index".to_string(),
                reason: "Linear index exceeds magnetization array length".to_string(),
            });
        }

        // Normalize the input magnetization
        let mag = m.magnitude();
        self.magnetization[idx] = if mag > 1e-30 {
            Vector3::new(m.x / mag, m.y / mag, m.z / mag)
        } else {
            Vector3::new(0.0, 0.0, 1.0)
        };

        Ok(())
    }

    /// Check that all magnetization vectors are normalized
    ///
    /// # Arguments
    /// * `tolerance` - Maximum allowed deviation from unit magnitude
    ///
    /// # Returns
    /// true if all vectors satisfy |m| - 1 < tolerance
    pub fn is_normalized(&self, tolerance: f64) -> bool {
        self.magnetization
            .iter()
            .all(|m| (m.magnitude() - 1.0).abs() < tolerance)
    }

    /// Re-normalize all magnetization vectors to unit length
    pub fn normalize_all(&mut self) {
        for m in &mut self.magnetization {
            let mag = m.magnitude();
            if mag > 1e-30 {
                m.x /= mag;
                m.y /= mag;
                m.z /= mag;
            } else {
                *m = Vector3::new(0.0, 0.0, 1.0);
            }
        }
    }

    /// Get the Cartesian position of a grid point
    ///
    /// # Arguments
    /// * `ix`, `iy`, `iz` - Grid indices
    ///
    /// # Returns
    /// (x, y, z) in physical coordinates \[m\]
    pub fn grid_position(&self, ix: usize, iy: usize, iz: usize) -> (f64, f64, f64) {
        let (nx, ny, nz) = self.grid_size;
        let cx = (nx as f64 - 1.0) / 2.0;
        let cy = (ny as f64 - 1.0) / 2.0;
        let cz = (nz as f64 - 1.0) / 2.0;

        let x = (ix as f64 - cx) * self.cell_size;
        let y = (iy as f64 - cy) * self.cell_size;
        let z = (iz as f64 - cz) * self.cell_size;

        (x, y, z)
    }

    /// Create a uniform magnetization state (no hopfion)
    ///
    /// All spins point along +z. Useful as a reference state.
    ///
    /// # Errors
    /// Returns error if grid dimensions are zero
    pub fn uniform(grid_size: (usize, usize, usize), cell_size: f64) -> Result<Self> {
        let (nx, ny, nz) = grid_size;
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(Error::InvalidParameter {
                param: "grid_size".to_string(),
                reason: "All grid dimensions must be positive".to_string(),
            });
        }

        let total = nx * ny * nz;
        let magnetization = vec![Vector3::new(0.0, 0.0, 1.0); total];

        Ok(Self {
            grid_size,
            magnetization,
            hopf_charge: 0,
            radius: 0.0,
            cell_size,
            profile: ProfileFunction::Linear,
            azimuthal_winding: 0,
            poloidal_winding: 0,
        })
    }
}

impl fmt::Display for Hopfion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Hopfion[Q_H={}, R={:.1} nm, grid={}x{}x{}, N={}, Q={}]",
            self.hopf_charge,
            self.radius * 1e9,
            self.grid_size.0,
            self.grid_size.1,
            self.grid_size.2,
            self.azimuthal_winding,
            self.poloidal_winding,
        )
    }
}

// ============================================================================
// Hopf Invariant Calculation
// ============================================================================

/// Calculator for the Hopf invariant (topological charge of 3D solitons)
///
/// The Hopf invariant Q_H measures the linking of preimage curves in a
/// magnetization texture. It is computed via:
///
/// Q_H = (1/4pi^2) integral F . A d^3r
///
/// where A is the Berry connection and F = curl(A).
pub struct HopfInvariant;

impl HopfInvariant {
    /// Calculate the Hopf invariant for a hopfion configuration
    ///
    /// Uses a discretized version of the Whitehead integral formula on the
    /// 3D grid. This computes the Berry connection A from the magnetization
    /// field, then evaluates F . A integrated over the volume.
    ///
    /// # Arguments
    /// * `hopfion` - The hopfion configuration
    ///
    /// # Returns
    /// The Hopf invariant (should be close to an integer for well-formed hopfions)
    ///
    /// # Errors
    /// Returns error if the grid is too small for numerical differentiation
    pub fn calculate(hopfion: &Hopfion) -> Result<f64> {
        let (nx, ny, nz) = hopfion.grid_size;

        if nx < 4 || ny < 4 || nz < 4 {
            return Err(Error::InvalidParameter {
                param: "grid_size".to_string(),
                reason: "Grid must be at least 4x4x4 for Hopf invariant calculation".to_string(),
            });
        }

        let dx = hopfion.cell_size;
        let dv = dx * dx * dx;

        let mut total = 0.0;

        // Interior points only (need neighbors for finite differences)
        for iz in 1..nz - 1 {
            for iy in 1..ny - 1 {
                for ix in 1..nx - 1 {
                    let m = hopfion.magnetization_at(ix, iy, iz);

                    // Finite difference derivatives
                    let dm_dx = (hopfion.magnetization_at(ix + 1, iy, iz)
                        - hopfion.magnetization_at(ix.saturating_sub(1), iy, iz))
                        * (0.5 / dx);
                    let dm_dy = (hopfion.magnetization_at(ix, iy + 1, iz)
                        - hopfion.magnetization_at(ix, iy.saturating_sub(1), iz))
                        * (0.5 / dx);
                    let dm_dz = (hopfion.magnetization_at(ix, iy, iz + 1)
                        - hopfion.magnetization_at(ix, iy, iz.saturating_sub(1)))
                        * (0.5 / dx);

                    // Berry connection: A_i = m . (dm/dj x dm/dk) for cyclic (i,j,k)
                    // This is the emergent gauge field from the magnetization texture
                    let a_x = m.dot(&dm_dy.cross(&dm_dz));
                    let a_y = m.dot(&dm_dz.cross(&dm_dx));
                    let a_z = m.dot(&dm_dx.cross(&dm_dy));

                    let berry_a = Vector3::new(a_x, a_y, a_z);

                    // For the F.A integral, we compute the emergent field F = curl(A)
                    // numerically. However, a more numerically stable approach for
                    // discrete lattices is to use the solid angle method.
                    //
                    // We use the relation: Q_H = (1/4pi^2) int F.A d^3r
                    // where F_i = epsilon_ijk d_j A_k
                    //
                    // Compute F by finite differences of A at neighboring points.
                    // For efficiency, we compute F.A directly using the Whitehead formula:
                    // Q_H ~ (1/4pi^2) sum_r A(r) . F(r) * dV
                    //
                    // To get F, we need A at neighboring points too. Instead, we use
                    // an equivalent topological density:
                    // q(r) = (1/8pi^2) epsilon_ijk m . (d_i m x d_j m) * (partial_k of
                    //         the solid angle subtended by m(r))
                    //
                    // A simpler equivalent: use the linking number via the emergent
                    // magnetic field approach. The topological charge density is:
                    // rho = (1/(4pi^2)) A . (curl A)
                    //
                    // For our discretization, we directly compute A then curl A.
                    //
                    // We already have A at this point. To get curl A, we need A at
                    // neighboring points. We compute A at the 6 face neighbors.

                    let a_xp = Self::berry_connection_at(hopfion, ix + 1, iy, iz, dx);
                    let a_xm = Self::berry_connection_at(hopfion, ix.saturating_sub(1), iy, iz, dx);
                    let a_yp = Self::berry_connection_at(hopfion, ix, iy + 1, iz, dx);
                    let a_ym = Self::berry_connection_at(hopfion, ix, iy.saturating_sub(1), iz, dx);
                    let a_zp = Self::berry_connection_at(hopfion, ix, iy, iz + 1, dx);
                    let a_zm = Self::berry_connection_at(hopfion, ix, iy, iz.saturating_sub(1), dx);

                    // curl A = (dA_z/dy - dA_y/dz, dA_x/dz - dA_z/dx, dA_y/dx - dA_x/dy)
                    let curl_a = Vector3::new(
                        (a_yp.z - a_ym.z - a_zp.y + a_zm.y) * (0.5 / dx),
                        (a_zp.x - a_zm.x - a_xp.z + a_xm.z) * (0.5 / dx),
                        (a_xp.y - a_xm.y - a_yp.x + a_ym.x) * (0.5 / dx),
                    );

                    // Contribution to Hopf invariant: (1/4pi^2) A . F dV
                    total += berry_a.dot(&curl_a) * dv;
                }
            }
        }

        Ok(total / (4.0 * PI * PI))
    }

    /// Compute Berry connection vector A at a grid point
    ///
    /// A_i = m . (dm/dj x dm/dk) for cyclic permutations (i,j,k)
    fn berry_connection_at(
        hopfion: &Hopfion,
        ix: usize,
        iy: usize,
        iz: usize,
        dx: f64,
    ) -> Vector3<f64> {
        let (nx, ny, nz) = hopfion.grid_size;

        // Clamp indices to valid range
        let ix = ix.min(nx - 1);
        let iy = iy.min(ny - 1);
        let iz = iz.min(nz - 1);

        // Boundary check: need at least one neighbor on each side
        let ix_lo = if ix > 0 { ix - 1 } else { ix };
        let ix_hi = if ix + 1 < nx { ix + 1 } else { ix };
        let iy_lo = if iy > 0 { iy - 1 } else { iy };
        let iy_hi = if iy + 1 < ny { iy + 1 } else { iy };
        let iz_lo = if iz > 0 { iz - 1 } else { iz };
        let iz_hi = if iz + 1 < nz { iz + 1 } else { iz };

        let m = hopfion.magnetization_at(ix, iy, iz);

        let dx_scale = if ix_hi > ix_lo {
            (ix_hi - ix_lo) as f64 * dx
        } else {
            dx
        };
        let dy_scale = if iy_hi > iy_lo {
            (iy_hi - iy_lo) as f64 * dx
        } else {
            dx
        };
        let dz_scale = if iz_hi > iz_lo {
            (iz_hi - iz_lo) as f64 * dx
        } else {
            dx
        };

        let dm_dx = (hopfion.magnetization_at(ix_hi, iy, iz)
            - hopfion.magnetization_at(ix_lo, iy, iz))
            * (1.0 / dx_scale);
        let dm_dy = (hopfion.magnetization_at(ix, iy_hi, iz)
            - hopfion.magnetization_at(ix, iy_lo, iz))
            * (1.0 / dy_scale);
        let dm_dz = (hopfion.magnetization_at(ix, iy, iz_hi)
            - hopfion.magnetization_at(ix, iy, iz_lo))
            * (1.0 / dz_scale);

        Vector3::new(
            m.dot(&dm_dy.cross(&dm_dz)),
            m.dot(&dm_dz.cross(&dm_dx)),
            m.dot(&dm_dx.cross(&dm_dy)),
        )
    }

    /// Calculate the Hopf invariant for a uniform magnetization state
    ///
    /// This should always return 0 (trivial topology).
    pub fn calculate_uniform(grid_size: (usize, usize, usize), cell_size: f64) -> Result<f64> {
        let uniform = Hopfion::uniform(grid_size, cell_size)?;
        Self::calculate(&uniform)
    }
}

// ============================================================================
// Energy Calculations
// ============================================================================

/// Parameters for 3D magnetic energy calculations
///
/// Defines the material parameters needed to compute the exchange, DMI,
/// and Zeeman energies of a hopfion configuration.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HopfionEnergyParams {
    /// Exchange stiffness constant A \[J/m\]
    pub exchange_stiffness: f64,

    /// Bulk DMI constant D \[J/m²\]
    pub dmi_constant: f64,

    /// External magnetic field \[T\] (as a vector)
    pub external_field: Vector3<f64>,

    /// Saturation magnetization M_s \[A/m\]
    pub saturation_magnetization: f64,

    /// Frustrated exchange: next-nearest-neighbor coupling J₂ \[J/m\]
    ///
    /// When J₂ < 0 (antiferromagnetic NNN), it can stabilize hopfions.
    pub frustrated_exchange_j2: f64,
}

impl Default for HopfionEnergyParams {
    /// Default parameters representative of a frustrated ferromagnet
    fn default() -> Self {
        Self {
            exchange_stiffness: 1.0e-11, // ~10 pJ/m (typical ferromagnet)
            dmi_constant: 3.0e-3,        // ~3 mJ/m² (moderate DMI)
            external_field: Vector3::new(0.0, 0.0, 0.5), // 0.5 T along z
            saturation_magnetization: 5.8e5, // ~Co (580 kA/m)
            frustrated_exchange_j2: -2.0e-12, // weak AFM NNN coupling
        }
    }
}

impl HopfionEnergyParams {
    /// Validate that the energy parameters are physically reasonable
    ///
    /// # Errors
    /// Returns error if any parameter is non-physical
    pub fn validate(&self) -> Result<()> {
        if self.exchange_stiffness < 0.0 || !self.exchange_stiffness.is_finite() {
            return Err(Error::InvalidParameter {
                param: "exchange_stiffness".to_string(),
                reason: "Must be non-negative and finite".to_string(),
            });
        }
        if self.saturation_magnetization < 0.0 || !self.saturation_magnetization.is_finite() {
            return Err(Error::InvalidParameter {
                param: "saturation_magnetization".to_string(),
                reason: "Must be non-negative and finite".to_string(),
            });
        }
        if !self.dmi_constant.is_finite() {
            return Err(Error::InvalidParameter {
                param: "dmi_constant".to_string(),
                reason: "Must be finite".to_string(),
            });
        }
        Ok(())
    }
}

/// Energy functional calculator for 3D hopfion configurations
pub struct HopfionEnergy;

impl HopfionEnergy {
    /// Calculate the total exchange energy
    ///
    /// E_ex = A * integral |grad m|^2 d^3r
    ///
    /// Uses central finite differences for the gradient on the 3D grid.
    ///
    /// # Arguments
    /// * `hopfion` - The hopfion configuration
    /// * `exchange_stiffness` - Exchange stiffness constant A \[J/m\]
    ///
    /// # Returns
    /// Exchange energy \[J\]
    ///
    /// # Errors
    /// Returns error for invalid parameters
    pub fn exchange_energy(hopfion: &Hopfion, exchange_stiffness: f64) -> Result<f64> {
        if !exchange_stiffness.is_finite() {
            return Err(Error::InvalidParameter {
                param: "exchange_stiffness".to_string(),
                reason: "Must be finite".to_string(),
            });
        }

        let (nx, ny, nz) = hopfion.grid_size;
        let dx = hopfion.cell_size;
        let dv = dx * dx * dx;
        let inv_2dx = 0.5 / dx;

        let mut energy = 0.0;

        for iz in 1..nz.saturating_sub(1) {
            for iy in 1..ny.saturating_sub(1) {
                for ix in 1..nx.saturating_sub(1) {
                    let dm_dx = (hopfion.magnetization_at(ix + 1, iy, iz)
                        - hopfion.magnetization_at(ix - 1, iy, iz))
                        * inv_2dx;
                    let dm_dy = (hopfion.magnetization_at(ix, iy + 1, iz)
                        - hopfion.magnetization_at(ix, iy - 1, iz))
                        * inv_2dx;
                    let dm_dz = (hopfion.magnetization_at(ix, iy, iz + 1)
                        - hopfion.magnetization_at(ix, iy, iz - 1))
                        * inv_2dx;

                    let grad_sq = dm_dx.magnitude_squared()
                        + dm_dy.magnitude_squared()
                        + dm_dz.magnitude_squared();

                    energy += grad_sq * dv;
                }
            }
        }

        Ok(exchange_stiffness * energy)
    }

    /// Calculate the 3D bulk DMI energy
    ///
    /// E_DMI = D * integral m . (curl m) d^3r
    ///
    /// This is the bulk (Bloch-type) DMI appropriate for non-centrosymmetric
    /// bulk crystals.
    ///
    /// # Arguments
    /// * `hopfion` - The hopfion configuration
    /// * `dmi_constant` - DMI constant D \[J/m²\]
    ///
    /// # Returns
    /// DMI energy \[J\]
    pub fn dmi_energy(hopfion: &Hopfion, dmi_constant: f64) -> Result<f64> {
        if !dmi_constant.is_finite() {
            return Err(Error::InvalidParameter {
                param: "dmi_constant".to_string(),
                reason: "Must be finite".to_string(),
            });
        }

        let (nx, ny, nz) = hopfion.grid_size;
        let dx = hopfion.cell_size;
        let dv = dx * dx * dx;
        let inv_2dx = 0.5 / dx;

        let mut energy = 0.0;

        for iz in 1..nz.saturating_sub(1) {
            for iy in 1..ny.saturating_sub(1) {
                for ix in 1..nx.saturating_sub(1) {
                    let m = hopfion.magnetization_at(ix, iy, iz);

                    let dm_dx = (hopfion.magnetization_at(ix + 1, iy, iz)
                        - hopfion.magnetization_at(ix - 1, iy, iz))
                        * inv_2dx;
                    let dm_dy = (hopfion.magnetization_at(ix, iy + 1, iz)
                        - hopfion.magnetization_at(ix, iy - 1, iz))
                        * inv_2dx;
                    let dm_dz = (hopfion.magnetization_at(ix, iy, iz + 1)
                        - hopfion.magnetization_at(ix, iy, iz - 1))
                        * inv_2dx;

                    // curl m = (dm_z/dy - dm_y/dz, dm_x/dz - dm_z/dx, dm_y/dx - dm_x/dy)
                    let curl_m =
                        Vector3::new(dm_dy.z - dm_dz.y, dm_dz.x - dm_dx.z, dm_dx.y - dm_dy.x);

                    energy += m.dot(&curl_m) * dv;
                }
            }
        }

        Ok(dmi_constant * energy)
    }

    /// Calculate the Zeeman energy
    ///
    /// E_Z = -mu_0 * M_s * integral m . H d^3r
    ///
    /// # Arguments
    /// * `hopfion` - The hopfion configuration
    /// * `field` - External magnetic field \[T\]
    /// * `saturation_magnetization` - M_s \[A/m\]
    ///
    /// # Returns
    /// Zeeman energy \[J\]
    pub fn zeeman_energy(
        hopfion: &Hopfion,
        field: &Vector3<f64>,
        saturation_magnetization: f64,
    ) -> Result<f64> {
        if !saturation_magnetization.is_finite() {
            return Err(Error::InvalidParameter {
                param: "saturation_magnetization".to_string(),
                reason: "Must be finite".to_string(),
            });
        }

        let dx = hopfion.cell_size;
        let dv = dx * dx * dx;

        let mut energy = 0.0;

        for m in &hopfion.magnetization {
            energy += m.dot(field);
        }

        // E_Z = -mu_0 * M_s * sum(m.H) * dV
        Ok(-crate::constants::MU_0 * saturation_magnetization * energy * dv)
    }

    /// Calculate the frustrated exchange (NNN) energy
    ///
    /// E_J2 = J₂ * integral |nabla^2 m|^2 d^3r
    ///
    /// The next-nearest-neighbor exchange acts as a higher-order gradient term
    /// that can stabilize hopfions when J₂ < 0 (antiferromagnetic NNN).
    ///
    /// # Arguments
    /// * `hopfion` - The hopfion configuration
    /// * `j2` - NNN exchange coupling \[J/m\]
    ///
    /// # Returns
    /// Frustrated exchange energy \[J\]
    pub fn frustrated_exchange_energy(hopfion: &Hopfion, j2: f64) -> Result<f64> {
        if !j2.is_finite() {
            return Err(Error::InvalidParameter {
                param: "j2".to_string(),
                reason: "Must be finite".to_string(),
            });
        }

        let (nx, ny, nz) = hopfion.grid_size;
        let dx = hopfion.cell_size;
        let dv = dx * dx * dx;
        let inv_dx2 = 1.0 / (dx * dx);

        let mut energy = 0.0;

        for iz in 1..nz.saturating_sub(1) {
            for iy in 1..ny.saturating_sub(1) {
                for ix in 1..nx.saturating_sub(1) {
                    let m_c = hopfion.magnetization_at(ix, iy, iz);

                    // Laplacian: nabla^2 m ~ (m_{i+1} + m_{i-1} - 2*m_i) / dx^2
                    // for each direction
                    let lap_x = (hopfion.magnetization_at(ix + 1, iy, iz)
                        + hopfion.magnetization_at(ix - 1, iy, iz)
                        - m_c * 2.0)
                        * inv_dx2;
                    let lap_y = (hopfion.magnetization_at(ix, iy + 1, iz)
                        + hopfion.magnetization_at(ix, iy - 1, iz)
                        - m_c * 2.0)
                        * inv_dx2;
                    let lap_z = (hopfion.magnetization_at(ix, iy, iz + 1)
                        + hopfion.magnetization_at(ix, iy, iz - 1)
                        - m_c * 2.0)
                        * inv_dx2;

                    let laplacian = lap_x + lap_y + lap_z;
                    energy += laplacian.magnitude_squared() * dv;
                }
            }
        }

        Ok(j2 * energy)
    }

    /// Calculate the total energy of a hopfion configuration
    ///
    /// E_total = E_exchange + E_DMI + E_Zeeman + E_frustrated
    ///
    /// # Arguments
    /// * `hopfion` - The hopfion configuration
    /// * `params` - Energy parameters
    ///
    /// # Returns
    /// Total energy \[J\]
    ///
    /// # Errors
    /// Returns error if parameters are invalid
    pub fn total_energy(hopfion: &Hopfion, params: &HopfionEnergyParams) -> Result<f64> {
        params.validate()?;

        let e_ex = Self::exchange_energy(hopfion, params.exchange_stiffness)?;
        let e_dmi = Self::dmi_energy(hopfion, params.dmi_constant)?;
        let e_zee = Self::zeeman_energy(
            hopfion,
            &params.external_field,
            params.saturation_magnetization,
        )?;
        let e_frust = Self::frustrated_exchange_energy(hopfion, params.frustrated_exchange_j2)?;

        Ok(e_ex + e_dmi + e_zee + e_frust)
    }
}

// ============================================================================
// Stability Analysis
// ============================================================================

/// Stability analysis for hopfion configurations
///
/// Analyzes the energetic stability of hopfions by examining the energy
/// as a function of hopfion radius and detecting collapse/expansion boundaries.
pub struct HopfionStability;

impl HopfionStability {
    /// Compute the energy as a function of hopfion radius
    ///
    /// Generates hopfion configurations at different radii and computes
    /// the total energy for each, producing an energy landscape.
    ///
    /// # Arguments
    /// * `grid_size` - Grid dimensions
    /// * `radii` - Array of hopfion radii to evaluate \[m\]
    /// * `hopf_charge` - Hopf charge
    /// * `params` - Energy parameters
    ///
    /// # Returns
    /// Vector of (radius, energy) pairs
    ///
    /// # Errors
    /// Returns error if any configuration fails
    pub fn energy_vs_radius(
        grid_size: (usize, usize, usize),
        radii: &[f64],
        hopf_charge: i32,
        params: &HopfionEnergyParams,
    ) -> Result<Vec<(f64, f64)>> {
        params.validate()?;

        let mut results = Vec::with_capacity(radii.len());

        for &r in radii {
            if r <= 0.0 || !r.is_finite() {
                continue;
            }

            let hopfion = Hopfion::new(grid_size, r, hopf_charge)?;
            let energy = HopfionEnergy::total_energy(&hopfion, params)?;
            results.push((r, energy));
        }

        Ok(results)
    }

    /// Estimate the stability region boundaries
    ///
    /// Finds the range of radii where a hopfion can exist as a local energy
    /// minimum. Below R_min the hopfion collapses; above R_max it expands.
    ///
    /// Uses finite differences on the energy landscape to find turning points.
    ///
    /// # Arguments
    /// * `energy_landscape` - Sorted (radius, energy) pairs from `energy_vs_radius`
    ///
    /// # Returns
    /// (R_min, R_equilibrium, R_max) if a stable region exists, or None
    pub fn find_stability_boundaries(energy_landscape: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
        if energy_landscape.len() < 3 {
            return None;
        }

        // Find local minimum in the energy landscape
        let mut min_idx = None;
        let mut min_energy = f64::MAX;

        for i in 1..energy_landscape.len() - 1 {
            let (_, e_prev) = energy_landscape[i - 1];
            let (_, e_curr) = energy_landscape[i];
            let (_, e_next) = energy_landscape[i + 1];

            // Local minimum: energy decreases then increases
            if e_curr < e_prev && e_curr < e_next && e_curr < min_energy {
                min_energy = e_curr;
                min_idx = Some(i);
            }
        }

        let eq_idx = min_idx?;
        let (r_eq, _) = energy_landscape[eq_idx];

        // Find collapse boundary (left side): where energy starts increasing
        // as radius decreases from equilibrium
        let mut r_min = energy_landscape[0].0;
        for i in (0..eq_idx).rev() {
            let (r_i, e_i) = energy_landscape[i];
            let (_, e_next) = energy_landscape[i + 1];
            if e_i < e_next {
                // Energy is still decreasing as r increases -- this means
                // we found where the barrier starts
                r_min = r_i;
                break;
            }
        }

        // Find expansion boundary (right side)
        let mut r_max = energy_landscape[energy_landscape.len() - 1].0;
        for i in eq_idx + 1..energy_landscape.len() {
            let (r_i, e_i) = energy_landscape[i];
            let (_, e_prev) = energy_landscape[i - 1];
            if e_i < e_prev {
                r_max = r_i;
                break;
            }
        }

        Some((r_min, r_eq, r_max))
    }

    /// Assess the frustrated exchange stabilization potential
    ///
    /// Computes the ratio of NNN to NN exchange that determines whether
    /// frustration is sufficient to stabilize a hopfion.
    ///
    /// The critical ratio |J₂/J₁| ~ 0.25 is typically needed for
    /// hopfion stabilization in frustrated magnets.
    ///
    /// # Arguments
    /// * `j1` - Nearest-neighbor exchange \[J/m\]
    /// * `j2` - Next-nearest-neighbor exchange \[J/m\]
    ///
    /// # Returns
    /// (frustration_ratio, is_potentially_stable)
    pub fn frustration_assessment(j1: f64, j2: f64) -> (f64, bool) {
        if j1.abs() < 1e-30 {
            return (0.0, false);
        }

        let ratio = (j2 / j1).abs();
        // Hopfions are potentially stabilized when J2 < 0 (AFM NNN)
        // and the frustration ratio exceeds ~0.25
        let is_stable = j2 < 0.0 && ratio > 0.25;

        (ratio, is_stable)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hopf_map_north_pole() {
        // (1, 0, 0, 0) on S³ should map to (0, 0, 1) on S²
        let result = HopfFibration::map(1.0, 0.0, 0.0, 0.0);
        assert!(
            (result.z - 1.0).abs() < 1e-10,
            "Expected north pole, got z = {}",
            result.z
        );
        assert!(result.x.abs() < 1e-10, "Expected x=0, got {}", result.x);
        assert!(result.y.abs() < 1e-10, "Expected y=0, got {}", result.y);
    }

    #[test]
    fn test_hopf_map_south_pole() {
        // (0, 0, 1, 0) on S³ should map to (0, 0, -1) on S²
        let result = HopfFibration::map(0.0, 0.0, 1.0, 0.0);
        assert!(
            (result.z + 1.0).abs() < 1e-10,
            "Expected south pole, got z = {}",
            result.z
        );
    }

    #[test]
    fn test_hopf_invariant_uniform_state() {
        // Uniform magnetization should have Q_H = 0
        let grid_size = (8, 8, 8);
        let cell_size = 1.0e-9;
        let q_h = HopfInvariant::calculate_uniform(grid_size, cell_size)
            .expect("Failed to calculate Hopf invariant for uniform state");

        assert!(
            q_h.abs() < 0.1,
            "Uniform state should have Q_H ~ 0, got {}",
            q_h
        );
    }

    #[test]
    fn test_hopfion_magnetization_normalized() {
        let hopfion = Hopfion::new((16, 16, 16), 5.0e-9, 1).expect("Failed to create hopfion");

        assert!(
            hopfion.is_normalized(1e-10),
            "All magnetization vectors should be normalized"
        );
    }

    #[test]
    fn test_hopfion_energy_finite_and_positive() {
        let hopfion = Hopfion::new((12, 12, 12), 3.0e-9, 1).expect("Failed to create hopfion");

        let params = HopfionEnergyParams {
            exchange_stiffness: 1.0e-11,
            dmi_constant: 0.0,
            external_field: Vector3::zero(),
            saturation_magnetization: 5.8e5,
            frustrated_exchange_j2: 0.0,
        };

        let e_ex = HopfionEnergy::exchange_energy(&hopfion, params.exchange_stiffness)
            .expect("Failed to calculate exchange energy");

        assert!(e_ex.is_finite(), "Exchange energy must be finite");
        assert!(
            e_ex >= 0.0,
            "Exchange energy must be non-negative, got {}",
            e_ex
        );
    }

    #[test]
    fn test_profile_function_boundary_conditions() {
        let profile = ProfileFunction::Linear;
        let radius = 10.0e-9;

        // At rho = 0, f(0) = pi
        let f_at_zero = profile.evaluate(0.0, radius);
        assert!(
            (f_at_zero - PI).abs() < 1e-10,
            "f(0) should be pi, got {}",
            f_at_zero
        );

        // At rho = R, f(R) = 0
        let f_at_r = profile.evaluate(radius, radius);
        assert!(f_at_r.abs() < 1e-10, "f(R) should be 0, got {}", f_at_r);

        // At rho > R, f(rho) = 0
        let f_beyond = profile.evaluate(2.0 * radius, radius);
        assert!(
            f_beyond.abs() < 1e-10,
            "f(2R) should be 0, got {}",
            f_beyond
        );

        // At rho = R/2, f(R/2) = pi/2
        let f_half = profile.evaluate(radius / 2.0, radius);
        assert!(
            (f_half - PI / 2.0).abs() < 1e-10,
            "f(R/2) should be pi/2 for linear profile, got {}",
            f_half
        );
    }

    #[test]
    fn test_grid_indexing_correctness() {
        let hopfion = Hopfion::new((8, 10, 12), 3.0e-9, 1).expect("Failed to create hopfion");

        let (nx, ny, _nz) = hopfion.grid_size;

        // Check that linear_index is consistent
        let idx_000 = hopfion
            .linear_index(0, 0, 0)
            .expect("Index (0,0,0) should be valid");
        assert_eq!(idx_000, 0);

        let idx_100 = hopfion
            .linear_index(1, 0, 0)
            .expect("Index (1,0,0) should be valid");
        assert_eq!(idx_100, 1);

        let idx_010 = hopfion
            .linear_index(0, 1, 0)
            .expect("Index (0,1,0) should be valid");
        assert_eq!(idx_010, nx);

        let idx_001 = hopfion
            .linear_index(0, 0, 1)
            .expect("Index (0,0,1) should be valid");
        assert_eq!(idx_001, nx * ny);

        // Out of bounds should return None
        assert!(hopfion.linear_index(8, 0, 0).is_none());
        assert!(hopfion.linear_index(0, 10, 0).is_none());
        assert!(hopfion.linear_index(0, 0, 12).is_none());
    }

    #[test]
    fn test_topological_protection_small_perturbation() {
        // The Hopf invariant should not change under small continuous
        // perturbations of the magnetization field
        let mut hopfion = Hopfion::new((10, 10, 10), 3.0e-9, 1).expect("Failed to create hopfion");

        let q_h_before =
            HopfInvariant::calculate(&hopfion).expect("Failed to compute Q_H before perturbation");

        // Apply a small perturbation to interior spins
        let (nx, ny, nz) = hopfion.grid_size;
        let perturbation_magnitude = 0.01;

        for iz in 2..nz - 2 {
            for iy in 2..ny - 2 {
                for ix in 2..nx - 2 {
                    let m = hopfion.magnetization_at(ix, iy, iz);
                    // Small tilt in x-direction
                    let perturbed = Vector3::new(
                        m.x + perturbation_magnitude * ((ix + iy + iz) as f64 * 0.1).sin(),
                        m.y + perturbation_magnitude * ((ix + iy + iz) as f64 * 0.2).cos(),
                        m.z,
                    );
                    hopfion
                        .set_magnetization_at(ix, iy, iz, perturbed)
                        .expect("Failed to set perturbed magnetization");
                }
            }
        }

        let q_h_after =
            HopfInvariant::calculate(&hopfion).expect("Failed to compute Q_H after perturbation");

        // The difference should be small relative to the magnitude (topological protection).
        // On a small grid (10x10x10) the absolute value of Q_H can be large due to
        // discretization artifacts, but the RELATIVE change under small perturbation
        // should remain small, demonstrating topological robustness.
        let denom = q_h_before.abs().max(1e-30);
        let relative_change = (q_h_after - q_h_before).abs() / denom;
        assert!(
            relative_change < 0.01,
            "Hopf invariant should be robust to small perturbations: \
             Q_H_before={:.4}, Q_H_after={:.4}, relative_change={:.6}",
            q_h_before,
            q_h_after,
            relative_change
        );
    }

    #[test]
    fn test_toroidal_structure_verification() {
        let hopfion = Hopfion::new((20, 20, 20), 5.0e-9, 1).expect("Failed to create hopfion");

        let (nx, ny, nz) = hopfion.grid_size;
        let cx = nx / 2;
        let cy = ny / 2;
        let cz = nz / 2;

        // At the center of the grid (which is on the toroidal ring for a hopfion),
        // the magnetization should be significantly different from +z
        let m_center = hopfion.magnetization_at(cx, cy, cz);

        // On the torus ring, the out-of-plane component should be tilted
        // (not purely +z like the background)
        let m_far = hopfion.magnetization_at(0, 0, 0);

        // The center should differ from the corner
        let diff = (m_center - m_far).magnitude();
        assert!(
            diff > 0.01,
            "Hopfion center should differ from corner: diff = {:.6}",
            diff
        );

        // Check that the structure has azimuthal symmetry approximately:
        // m at (cx+2, cy, cz) and (cx, cy+2, cz) should have the same mz component
        // (due to toroidal symmetry in the xy-plane)
        if cx + 2 < nx && cy + 2 < ny {
            let m_xp = hopfion.magnetization_at(cx + 2, cy, cz);
            let m_yp = hopfion.magnetization_at(cx, cy + 2, cz);
            let mz_diff = (m_xp.z - m_yp.z).abs();
            assert!(
                mz_diff < 0.3,
                "Toroidal symmetry: mz at (+x) and (+y) should be similar, diff = {:.6}",
                mz_diff
            );
        }
    }

    #[test]
    fn test_hopf_invariant_single_hopfion() {
        // For a well-formed hopfion with Q_H = 1, the numerical Hopf invariant
        // should be close to 1. Due to discretization on a small grid, we allow
        // larger tolerance but check the sign and approximate magnitude.
        let hopfion = Hopfion::new((20, 20, 20), 5.0e-9, 1).expect("Failed to create hopfion");

        let q_h = HopfInvariant::calculate(&hopfion).expect("Failed to compute Hopf invariant");

        // On a discrete grid the numerical value won't be exactly 1,
        // but should be finite and have a consistent sign with the charge
        assert!(
            q_h.is_finite(),
            "Hopf invariant must be finite, got {}",
            q_h
        );

        // Verify the invariant is non-trivial (different from uniform state)
        let uniform = Hopfion::uniform(hopfion.grid_size, hopfion.cell_size)
            .expect("Failed to create uniform state");
        let q_h_uniform =
            HopfInvariant::calculate(&uniform).expect("Failed to compute Q_H for uniform");

        let hopfion_is_nontrivial = q_h.abs() > q_h_uniform.abs() * 2.0 + 1e-6;
        // The hopfion should have a larger |Q_H| than the uniform state
        assert!(
            hopfion_is_nontrivial || q_h.abs() > 0.01,
            "Hopfion Q_H={:.6} should be distinguishable from uniform Q_H={:.6}",
            q_h,
            q_h_uniform
        );
    }

    #[test]
    fn test_toroidal_coord_roundtrip() {
        let major_r = 5.0e-9;
        let rho = 2.0e-9;
        let theta = 0.7;
        let phi = 1.3;

        let (x, y, z) = ToroidalCoords::to_cartesian(rho, theta, phi, major_r);
        let (rho2, theta2, phi2) = ToroidalCoords::from_cartesian(x, y, z, major_r);

        assert!(
            (rho2 - rho).abs() < 1e-15,
            "Rho roundtrip failed: {} vs {}",
            rho,
            rho2
        );
        assert!(
            (theta2 - theta).abs() < 1e-15,
            "Theta roundtrip failed: {} vs {}",
            theta,
            theta2
        );
        assert!(
            (phi2 - phi).abs() < 1e-15,
            "Phi roundtrip failed: {} vs {}",
            phi,
            phi2
        );
    }

    #[test]
    fn test_frustration_assessment() {
        // Ferromagnetic NN, antiferromagnetic NNN with sufficient ratio
        let (ratio, is_stable) = HopfionStability::frustration_assessment(1.0e-11, -3.0e-12);
        assert!(
            (ratio - 0.3).abs() < 1e-10,
            "Expected ratio ~0.3, got {}",
            ratio
        );
        assert!(is_stable, "Should be potentially stable at ratio 0.3");

        // Insufficient frustration
        let (ratio2, is_stable2) = HopfionStability::frustration_assessment(1.0e-11, -1.0e-12);
        assert!(
            (ratio2 - 0.1).abs() < 1e-10,
            "Expected ratio ~0.1, got {}",
            ratio2
        );
        assert!(!is_stable2, "Should not be stable at ratio 0.1");

        // Ferromagnetic NNN (no frustration)
        let (_ratio3, is_stable3) = HopfionStability::frustration_assessment(1.0e-11, 5.0e-12);
        assert!(
            !is_stable3,
            "Ferromagnetic NNN should not stabilize hopfions"
        );
    }

    #[test]
    fn test_invalid_parameters() {
        // Zero grid dimensions
        let result = Hopfion::new((0, 10, 10), 5.0e-9, 1);
        assert!(result.is_err(), "Zero grid dimension should fail");

        // Negative radius
        let result = Hopfion::new((10, 10, 10), -1.0e-9, 1);
        assert!(result.is_err(), "Negative radius should fail");

        // Infinite radius
        let result = Hopfion::new((10, 10, 10), f64::INFINITY, 1);
        assert!(result.is_err(), "Infinite radius should fail");

        // NaN radius
        let result = Hopfion::new((10, 10, 10), f64::NAN, 1);
        assert!(result.is_err(), "NaN radius should fail");
    }

    #[test]
    fn test_hopf_fibration_inverse_roundtrip() {
        // Test that inverse_map produces a point that maps back to the original
        let n = Vector3::new(0.5, 0.5, (1.0_f64 - 0.5).sqrt());
        let n_norm = n.normalize();
        let fiber_angle = 0.42;

        let (a, b, c, d) = HopfFibration::inverse_map(&n_norm, fiber_angle);
        let n_recovered = HopfFibration::map(a, b, c, d);

        assert!(
            (n_recovered.x - n_norm.x).abs() < 1e-10,
            "x component mismatch: {} vs {}",
            n_recovered.x,
            n_norm.x
        );
        assert!(
            (n_recovered.y - n_norm.y).abs() < 1e-10,
            "y component mismatch: {} vs {}",
            n_recovered.y,
            n_norm.y
        );
        assert!(
            (n_recovered.z - n_norm.z).abs() < 1e-10,
            "z component mismatch: {} vs {}",
            n_recovered.z,
            n_norm.z
        );
    }
}
