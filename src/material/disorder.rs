//! Random anisotropy and disorder models for magnetic materials
//!
//! This module implements disorder models relevant to nanocrystalline and
//! polycrystalline magnetic materials:
//!
//! - **Random Anisotropy Model (RAM)**: Herzer model for effective anisotropy
//!   reduction in nanocrystalline soft magnets, where K_eff = K₁/√N.
//! - **Grain structure generation**: Simplified Voronoi-based grain assignment
//!   with random anisotropy axes per grain.
//! - **Random field disorder**: Gaussian-distributed random fields modeling
//!   dilute impurities and random substitution.
//! - **Surface roughness**: Height profile with RMS roughness and correlation
//!   length, affecting effective surface anisotropy.
//! - **Graded interfaces**: Smooth Linear/Exponential/error-function
//!   transition profiles for a spatially-varying scalar material parameter
//!   (Ms, exchange stiffness A, or anisotropy K) across a compositionally
//!   graded or interdiffused interface region, with optional roughness-like
//!   jitter of the transition position.
//!
//! # References
//!
//! - G. Herzer, IEEE Trans. Magn. 26, 1397 (1990)
//! - R. Alben, J.J. Becker, M.C. Chi, J. Appl. Phys. 49, 1653 (1978)

use crate::error::{Error, Result};
use crate::vector3::Vector3;

// ============================================================================
// Deterministic PRNG (xorshift64)
// ============================================================================

/// A simple deterministic xorshift64 pseudo-random number generator.
///
/// This avoids any dependency on external random number crates while providing
/// reproducible sequences suitable for disorder configuration generation.
#[derive(Debug, Clone)]
pub struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    /// Create a new PRNG with the given seed.
    ///
    /// The seed must be non-zero; if zero is provided it is replaced with a
    /// default non-zero value.
    pub fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 {
                0x5EED_DEAD_BEEF_CAFE
            } else {
                seed
            },
        }
    }

    /// Generate the next raw u64 value.
    pub fn next_u64(&mut self) -> u64 {
        let mut s = self.state;
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        self.state = s;
        s
    }

    /// Generate a uniform f64 in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }

    /// Generate a standard-normal f64 using the Box-Muller transform.
    ///
    /// Returns a pair of independent normal samples.
    pub fn next_normal_pair(&mut self) -> (f64, f64) {
        // Ensure u1 > 0 to avoid ln(0)
        let mut u1 = self.next_f64();
        while u1 <= f64::EPSILON {
            u1 = self.next_f64();
        }
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        (r * theta.cos(), r * theta.sin())
    }

    /// Generate a single standard-normal f64.
    pub fn next_normal(&mut self) -> f64 {
        self.next_normal_pair().0
    }

    /// Generate a uniform random unit vector on the sphere (isotropic).
    ///
    /// Uses the Marsaglia method: pick (x1, x2) uniform in the unit disk,
    /// then map to the sphere.
    pub fn next_unit_vector(&mut self) -> Vector3<f64> {
        loop {
            let x1 = 2.0 * self.next_f64() - 1.0;
            let x2 = 2.0 * self.next_f64() - 1.0;
            let s = x1 * x1 + x2 * x2;
            if s < 1.0 && s > f64::EPSILON {
                let factor = (1.0 - s).sqrt();
                return Vector3::new(2.0 * x1 * factor, 2.0 * x2 * factor, 1.0 - 2.0 * s);
            }
        }
    }
}

// ============================================================================
// Disorder types
// ============================================================================

/// Type of disorder present in the magnetic system.
#[derive(Debug, Clone)]
pub enum DisorderType {
    /// Random anisotropy: each grain has a randomly oriented easy axis.
    RandomAnisotropy,
    /// Random field: spatially varying random magnetic field.
    RandomField,
    /// Site disorder with a given fraction of vacancies.
    SiteDisorder {
        /// Fraction of sites that are vacant (0 to 1).
        vacancy_fraction: f64,
    },
}

/// Configuration for generating a disordered system.
#[derive(Debug, Clone)]
pub struct DisorderConfig {
    /// Type of disorder to apply.
    pub disorder_type: DisorderType,
    /// Disorder strength (units depend on type: J/m³ for anisotropy, T for field).
    pub strength: f64,
    /// Spatial correlation length of the disorder \[m\].
    pub correlation_length: f64,
    /// Seed for the deterministic PRNG.
    pub seed: u64,
}

// ============================================================================
// Random Anisotropy Model
// ============================================================================

/// Random anisotropy model (RAM) following the Herzer scaling theory.
///
/// In nanocrystalline materials the effective anisotropy is averaged over
/// the exchange-correlated volume, leading to:
///
///   K_eff = K₁ / √N
///
/// where N = (L_ex / D)³ is the number of grains within the exchange length
/// and L_ex = √(A / K₁) is the ferromagnetic exchange length.
///
/// # Example
///
/// ```
/// use spintronics::material::disorder::RandomAnisotropyModel;
///
/// let ram = RandomAnisotropyModel::new(
///     100,     // 100 grains
///     10.0e-9, // 10 nm grain size
///     1.0e4,   // K₁ = 10 kJ/m³
///     1.0e-11, // A = 10 pJ/m
///     42,      // seed
/// ).expect("valid parameters");
///
/// assert!(ram.effective_anisotropy < ram.bulk_anisotropy);
/// ```
#[derive(Debug, Clone)]
pub struct RandomAnisotropyModel {
    /// Number of grains in the model.
    pub num_grains: usize,
    /// Average grain diameter \[m\].
    pub grain_size: f64,
    /// Bulk magnetocrystalline anisotropy constant K₁ \[J/m³\].
    pub bulk_anisotropy: f64,
    /// Exchange stiffness constant A \[J/m\].
    pub exchange_stiffness: f64,
    /// Random easy axis direction for each grain (unit vectors).
    pub grain_axes: Vec<Vector3<f64>>,
    /// Effective anisotropy K_eff from Herzer averaging \[J/m³\].
    pub effective_anisotropy: f64,
}

impl RandomAnisotropyModel {
    /// Create a new random anisotropy model.
    ///
    /// # Arguments
    ///
    /// * `num_grains` - Number of grains (must be > 0)
    /// * `grain_size` - Average grain diameter D \[m\] (must be > 0)
    /// * `bulk_anisotropy` - K₁ \[J/m³\] (must be > 0)
    /// * `exchange_stiffness` - A \[J/m\] (must be > 0)
    /// * `seed` - PRNG seed for reproducibility
    ///
    /// # Errors
    ///
    /// Returns `Error::InvalidParameter` if any physical parameter is non-positive.
    pub fn new(
        num_grains: usize,
        grain_size: f64,
        bulk_anisotropy: f64,
        exchange_stiffness: f64,
        seed: u64,
    ) -> Result<Self> {
        if num_grains == 0 {
            return Err(Error::InvalidParameter {
                param: "num_grains".to_string(),
                reason: "must be greater than zero".to_string(),
            });
        }
        if grain_size <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "grain_size".to_string(),
                reason: "must be positive".to_string(),
            });
        }
        if bulk_anisotropy <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "bulk_anisotropy".to_string(),
                reason: "must be positive".to_string(),
            });
        }
        if exchange_stiffness <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "exchange_stiffness".to_string(),
                reason: "must be positive".to_string(),
            });
        }

        let mut rng = Xorshift64::new(seed);
        let grain_axes: Vec<Vector3<f64>> =
            (0..num_grains).map(|_| rng.next_unit_vector()).collect();

        let l_ex = Self::exchange_length_static(exchange_stiffness, bulk_anisotropy);
        let n_grains_in_volume = (l_ex / grain_size).powi(3);
        // Herzer: K_eff = K₁ / √N  (N = number of grains in exchange volume)
        let effective_anisotropy = if n_grains_in_volume > 1.0 {
            bulk_anisotropy / n_grains_in_volume.sqrt()
        } else {
            // If the grain is larger than exchange length, no averaging occurs
            bulk_anisotropy
        };

        Ok(Self {
            num_grains,
            grain_size,
            bulk_anisotropy,
            exchange_stiffness,
            grain_axes,
            effective_anisotropy,
        })
    }

    /// Compute the ferromagnetic exchange length L_ex = √(A / K₁) \[m\].
    pub fn exchange_length(&self) -> f64 {
        Self::exchange_length_static(self.exchange_stiffness, self.bulk_anisotropy)
    }

    /// Static helper for exchange length calculation.
    fn exchange_length_static(exchange_stiffness: f64, bulk_anisotropy: f64) -> f64 {
        (exchange_stiffness / bulk_anisotropy).sqrt()
    }

    /// Number of grains within one exchange volume: N = (L_ex / D)³.
    pub fn grains_in_exchange_volume(&self) -> f64 {
        let l_ex = self.exchange_length();
        (l_ex / self.grain_size).powi(3)
    }

    /// Herzer scaling exponent check: K_eff should scale as K₁⁴ D⁶ / A³.
    ///
    /// Returns the ratio K_eff / (K₁⁴ D⁶ / A³) which should be of order unity
    /// in the regime where D << L_ex.
    pub fn herzer_scaling_ratio(&self) -> f64 {
        let numerator = self.effective_anisotropy;
        let denominator = self.bulk_anisotropy.powi(4) * self.grain_size.powi(6)
            / self.exchange_stiffness.powi(3);
        if denominator.abs() < f64::EPSILON {
            return 0.0;
        }
        numerator / denominator
    }
}

// ============================================================================
// Grain structure
// ============================================================================

/// A single grain in a polycrystalline material.
#[derive(Debug, Clone)]
pub struct Grain {
    /// Position of the grain center in normalized coordinates [0, 1)³.
    pub center: Vector3<f64>,
    /// Local easy axis direction (unit vector).
    pub easy_axis: Vector3<f64>,
    /// Grain diameter \[m\].
    pub diameter: f64,
}

/// A simplified grain structure based on nearest-neighbor assignment to
/// randomly placed grain centers (Voronoi-like tessellation).
#[derive(Debug, Clone)]
pub struct GrainStructure {
    /// Individual grains with their properties.
    pub grains: Vec<Grain>,
    /// System size in each dimension \[m\].
    pub system_size: Vector3<f64>,
}

impl GrainStructure {
    /// Generate a grain structure with `num_grains` randomly placed grain centers.
    ///
    /// Each grain receives a random easy axis and a diameter drawn from a uniform
    /// distribution around `mean_diameter` with spread `diameter_spread`.
    ///
    /// # Arguments
    ///
    /// * `num_grains` - Number of grains to generate
    /// * `system_size` - Physical dimensions of the simulation box \[m\]
    /// * `mean_diameter` - Mean grain diameter \[m\]
    /// * `diameter_spread` - Half-width of the uniform diameter distribution \[m\]
    /// * `seed` - PRNG seed
    ///
    /// # Errors
    ///
    /// Returns an error if parameters are non-positive.
    pub fn generate(
        num_grains: usize,
        system_size: Vector3<f64>,
        mean_diameter: f64,
        diameter_spread: f64,
        seed: u64,
    ) -> Result<Self> {
        if num_grains == 0 {
            return Err(Error::InvalidParameter {
                param: "num_grains".to_string(),
                reason: "must be greater than zero".to_string(),
            });
        }
        if mean_diameter <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "mean_diameter".to_string(),
                reason: "must be positive".to_string(),
            });
        }
        if diameter_spread < 0.0 {
            return Err(Error::InvalidParameter {
                param: "diameter_spread".to_string(),
                reason: "must be non-negative".to_string(),
            });
        }
        if system_size.x <= 0.0 || system_size.y <= 0.0 || system_size.z <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "system_size".to_string(),
                reason: "all dimensions must be positive".to_string(),
            });
        }

        let mut rng = Xorshift64::new(seed);
        let mut grains = Vec::with_capacity(num_grains);

        for _ in 0..num_grains {
            let center = Vector3::new(rng.next_f64(), rng.next_f64(), rng.next_f64());
            let easy_axis = rng.next_unit_vector();
            // Uniform distribution: [mean - spread, mean + spread]
            let diameter = if diameter_spread > f64::EPSILON {
                let u = rng.next_f64();
                let d = mean_diameter + diameter_spread * (2.0 * u - 1.0);
                if d > 0.0 {
                    d
                } else {
                    mean_diameter * 0.1
                }
            } else {
                mean_diameter
            };

            grains.push(Grain {
                center,
                easy_axis,
                diameter,
            });
        }

        Ok(Self {
            grains,
            system_size,
        })
    }

    /// Assign a spatial point (in normalized coordinates [0,1)³) to the nearest grain.
    ///
    /// Returns the grain index or `None` if the structure has no grains.
    pub fn nearest_grain(&self, point: Vector3<f64>) -> Option<usize> {
        if self.grains.is_empty() {
            return None;
        }
        let mut best_idx = 0;
        let mut best_dist_sq = f64::MAX;
        for (i, grain) in self.grains.iter().enumerate() {
            let dx = point.x - grain.center.x;
            let dy = point.y - grain.center.y;
            let dz = point.z - grain.center.z;
            let dist_sq = dx * dx + dy * dy + dz * dz;
            if dist_sq < best_dist_sq {
                best_dist_sq = dist_sq;
                best_idx = i;
            }
        }
        Some(best_idx)
    }
}

// ============================================================================
// Random field disorder
// ============================================================================

/// Random field configuration at discrete lattice sites.
///
/// Each site gets a random field vector drawn from a Gaussian distribution
/// with standard deviation `field_strength` (in Tesla) and optional spatial
/// correlation length.
#[derive(Debug, Clone)]
pub struct RandomFieldDisorder {
    /// Random field vectors at each site \[T\].
    pub fields: Vec<Vector3<f64>>,
    /// RMS field strength \[T\].
    pub field_strength: f64,
    /// Spatial correlation length \[m\] (0 = uncorrelated).
    pub correlation_length: f64,
}

impl RandomFieldDisorder {
    /// Generate random field disorder for `num_sites` lattice points.
    ///
    /// # Arguments
    ///
    /// * `num_sites` - Number of lattice sites
    /// * `field_strength` - RMS amplitude of the random field \[T\]
    /// * `correlation_length` - Spatial correlation length \[m\] (currently
    ///   uncorrelated; correlation is recorded for metadata)
    /// * `seed` - PRNG seed
    ///
    /// # Errors
    ///
    /// Returns an error if `num_sites` is zero or `field_strength` is negative.
    pub fn generate(
        num_sites: usize,
        field_strength: f64,
        correlation_length: f64,
        seed: u64,
    ) -> Result<Self> {
        if num_sites == 0 {
            return Err(Error::InvalidParameter {
                param: "num_sites".to_string(),
                reason: "must be greater than zero".to_string(),
            });
        }
        if field_strength < 0.0 {
            return Err(Error::InvalidParameter {
                param: "field_strength".to_string(),
                reason: "must be non-negative".to_string(),
            });
        }

        let mut rng = Xorshift64::new(seed);
        let mut fields = Vec::with_capacity(num_sites);

        for _ in 0..num_sites {
            let hx = field_strength * rng.next_normal();
            let hy = field_strength * rng.next_normal();
            let hz = field_strength * rng.next_normal();
            fields.push(Vector3::new(hx, hy, hz));
        }

        Ok(Self {
            fields,
            field_strength,
            correlation_length,
        })
    }

    /// Mean magnitude of the random fields.
    pub fn mean_magnitude(&self) -> f64 {
        if self.fields.is_empty() {
            return 0.0;
        }
        let sum: f64 = self
            .fields
            .iter()
            .map(|f| (f.x * f.x + f.y * f.y + f.z * f.z).sqrt())
            .sum();
        sum / self.fields.len() as f64
    }
}

// ============================================================================
// Surface roughness
// ============================================================================

/// Surface roughness profile h(x, y) on a 2D grid.
///
/// The roughness is characterized by an RMS roughness σ and a correlation
/// length ξ. The height profile is generated as Gaussian noise and then
/// optionally smoothed to impose the correlation length.
#[derive(Debug, Clone)]
pub struct SurfaceRoughness {
    /// Height values on a 2D grid \[m\].
    pub height_map: Vec<Vec<f64>>,
    /// Number of grid points in x.
    pub nx: usize,
    /// Number of grid points in y.
    pub ny: usize,
    /// RMS roughness σ \[m\].
    pub rms_roughness: f64,
    /// Correlation length ξ \[m\].
    pub correlation_length: f64,
    /// Grid spacing \[m\].
    pub grid_spacing: f64,
}

impl SurfaceRoughness {
    /// Generate a surface roughness profile on an nx × ny grid.
    ///
    /// # Arguments
    ///
    /// * `nx`, `ny` - Grid dimensions (both must be > 0)
    /// * `grid_spacing` - Distance between grid points \[m\]
    /// * `rms_roughness` - RMS height σ \[m\]
    /// * `correlation_length` - Spatial correlation ξ \[m\]
    /// * `seed` - PRNG seed
    ///
    /// # Errors
    ///
    /// Returns an error if grid dimensions are zero or physical parameters are negative.
    pub fn generate(
        nx: usize,
        ny: usize,
        grid_spacing: f64,
        rms_roughness: f64,
        correlation_length: f64,
        seed: u64,
    ) -> Result<Self> {
        if nx == 0 || ny == 0 {
            return Err(Error::InvalidParameter {
                param: "grid dimensions".to_string(),
                reason: "nx and ny must be greater than zero".to_string(),
            });
        }
        if grid_spacing <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "grid_spacing".to_string(),
                reason: "must be positive".to_string(),
            });
        }
        if rms_roughness < 0.0 {
            return Err(Error::InvalidParameter {
                param: "rms_roughness".to_string(),
                reason: "must be non-negative".to_string(),
            });
        }
        if correlation_length < 0.0 {
            return Err(Error::InvalidParameter {
                param: "correlation_length".to_string(),
                reason: "must be non-negative".to_string(),
            });
        }

        let mut rng = Xorshift64::new(seed);

        // Generate uncorrelated Gaussian noise
        let mut height_map: Vec<Vec<f64>> = Vec::with_capacity(nx);
        for _ in 0..nx {
            let mut row = Vec::with_capacity(ny);
            for _ in 0..ny {
                row.push(rms_roughness * rng.next_normal());
            }
            height_map.push(row);
        }

        // Apply simple Gaussian smoothing if correlation length > grid spacing
        if correlation_length > grid_spacing {
            let kernel_radius = (correlation_length / grid_spacing).ceil() as usize;
            let sigma_grid = correlation_length / grid_spacing;
            let mut smoothed = vec![vec![0.0; ny]; nx];

            for ix in 0..nx {
                for iy in 0..ny {
                    let mut weight_sum = 0.0;
                    let mut val_sum = 0.0;

                    let jx_lo = ix.saturating_sub(kernel_radius);
                    let jx_hi = (ix + kernel_radius + 1).min(nx);
                    let jy_lo = iy.saturating_sub(kernel_radius);
                    let jy_hi = (iy + kernel_radius + 1).min(ny);

                    for (jx, height_row) in height_map.iter().enumerate().take(jx_hi).skip(jx_lo) {
                        for (jy, height_val) in
                            height_row.iter().enumerate().take(jy_hi).skip(jy_lo)
                        {
                            let dx = jx as f64 - ix as f64;
                            let dy = jy as f64 - iy as f64;
                            let r2 = dx * dx + dy * dy;
                            let w = (-r2 / (2.0 * sigma_grid * sigma_grid)).exp();
                            weight_sum += w;
                            val_sum += w * height_val;
                        }
                    }

                    smoothed[ix][iy] = if weight_sum > f64::EPSILON {
                        val_sum / weight_sum
                    } else {
                        height_map[ix][iy]
                    };
                }
            }

            // Rescale to restore target RMS
            let current_rms = Self::compute_rms(&smoothed, nx, ny);
            if current_rms > f64::EPSILON {
                let scale = rms_roughness / current_rms;
                for row in &mut smoothed {
                    for h in row.iter_mut() {
                        *h *= scale;
                    }
                }
            }
            height_map = smoothed;
        }

        Ok(Self {
            height_map,
            nx,
            ny,
            rms_roughness,
            correlation_length,
            grid_spacing,
        })
    }

    /// Compute the RMS of a 2D height map.
    fn compute_rms(map: &[Vec<f64>], nx: usize, ny: usize) -> f64 {
        let n = (nx * ny) as f64;
        if n < 1.0 {
            return 0.0;
        }
        let sum_sq: f64 = map.iter().flat_map(|row| row.iter()).map(|h| h * h).sum();
        (sum_sq / n).sqrt()
    }

    /// Get the actual RMS roughness of the generated height map.
    pub fn actual_rms(&self) -> f64 {
        Self::compute_rms(&self.height_map, self.nx, self.ny)
    }

    /// Estimate the effective surface anisotropy modification due to roughness.
    ///
    /// The roughness-induced effective anisotropy scales as:
    ///   ΔK_s ∝ K_s · (σ / ξ)²
    ///
    /// where K_s is the intrinsic surface anisotropy.
    ///
    /// # Arguments
    ///
    /// * `intrinsic_surface_anisotropy` - K_s \[J/m²\]
    ///
    /// # Returns
    ///
    /// The anisotropy modification ΔK_s \[J/m²\].
    pub fn anisotropy_modification(&self, intrinsic_surface_anisotropy: f64) -> f64 {
        if self.correlation_length < f64::EPSILON {
            return 0.0;
        }
        let ratio = self.rms_roughness / self.correlation_length;
        intrinsic_surface_anisotropy * ratio * ratio
    }
}

// ============================================================================
// Graded (inhomogeneous) interfaces
// ============================================================================

/// Abramowitz-Stegun 7.1.26 polynomial approximation of the Gauss error
/// function `erf(x)`.
///
/// The maximum absolute error over all real `x` is `≤ 1.5×10⁻⁷`, which is
/// more than sufficient for graded-interface material profiles (the
/// physical width/shape of a real interdiffused interface is itself
/// typically known only to a few percent). Implemented locally rather than
/// reusing the `autodiff` module's private `erf_approx` so that graded
/// interfaces remain available under the crate's default feature set.
fn erf_approx(x: f64) -> f64 {
    const A1: f64 = 0.254_829_592;
    const A2: f64 = -0.284_496_736;
    const A3: f64 = 1.421_413_741;
    const A4: f64 = -1.453_152_027;
    const A5: f64 = 1.061_405_429;
    const P: f64 = 0.327_591_1;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let ax = x.abs();
    let t = 1.0 / (1.0 + P * ax);
    let poly = (((A5 * t + A4) * t + A3) * t + A2) * t + A1;
    let y = 1.0 - poly * t * (-(ax * ax)).exp();
    sign * y
}

/// Spatial grading law describing how a scalar material parameter (e.g.
/// saturation magnetisation `Ms`, exchange stiffness `A`, or anisotropy `K`)
/// transitions between two endpoint values across an interface region.
///
/// Writing `d = position − center` and `v̄ = (value_left + value_right) / 2`
/// for the mean of the two endpoints, every law below satisfies the exact
/// point-antisymmetry identity
///
/// ```text
/// value(center + d) + value(center − d) = value_left + value_right
/// ```
///
/// for *every* `d`, not merely in some asymptotic limit. This is what makes
/// the barycenter of any position window that is symmetric about `center`
/// equal exactly to the arithmetic mean of the two endpoints, independent of
/// the transition `width` or the sampling density used to estimate it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum GradingLaw {
    /// Piecewise-linear ramp with *compact support*: exactly equal to
    /// `value_left` / `value_right` once `|position − center| ≥ width / 2`,
    /// and linearly interpolated in between. Here `width` is the full
    /// transition span (end to end).
    Linear,
    /// Two-sided exponential relaxation from each endpoint toward the mean
    /// value, with decay length `width`; continuously differentiable (C¹)
    /// at `center`. The endpoint values are reached only asymptotically:
    /// the residual at `|d| = n · width` is a factor `exp(−n)` of
    /// `(value_right − value_left) / 2`.
    Exponential,
    /// Smooth error-function (erf) sigmoid of characteristic width `width`
    /// — the Boltzmann-Matano solution of the 1D diffusion equation for a
    /// step initial condition, and the standard textbook model for
    /// interdiffusion-broadened magnetic interfaces. Endpoint values are
    /// reached only asymptotically, with residual `erfc(n)` at
    /// `|d| = n · width`.
    ErrorFunction,
}

/// A graded (inhomogeneous) interface: a scalar material parameter that
/// varies smoothly across an interface transition region instead of
/// jumping discontinuously between two bulk values, with optional
/// site-to-site roughness broadening of the transition position.
///
/// The nominal (unperturbed) profile is centered at `center` with
/// characteristic transition `width`, saturating to `value_left` for
/// `position ≪ center` and to `value_right` for `position ≫ center`. The
/// model is generic over *any* scalar material parameter — `Ms`, exchange
/// stiffness `A`, anisotropy `K`, etc. — graded across a compositionally
/// interdiffused or roughness-broadened interface.
///
/// # Example
///
/// ```
/// use spintronics::material::disorder::{GradedInterface, GradingLaw};
///
/// // Ms graded across a 2 nm interdiffusion layer between two ferromagnets.
/// let profile = GradedInterface::generate(
///     GradingLaw::ErrorFunction,
///     0.0,      // interface center at the origin
///     2.0e-9,   // 2 nm characteristic width
///     4.0e5,    // Ms = 4e5 A/m on the left
///     8.0e5,    // Ms = 8e5 A/m on the right
///     1,        // one lateral site (no roughness realization needed)
///     0.0,      // no roughness jitter
///     42,
/// ).expect("valid parameters");
///
/// assert!((profile.value_at(0.0) - 6.0e5).abs() < 1.0);
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GradedInterface {
    /// Grading law used to interpolate between the two endpoint values.
    pub law: GradingLaw,
    /// Nominal position of the interface transition midpoint (same length
    /// unit as the `position` argument of [`GradedInterface::value_at`]).
    pub center: f64,
    /// Characteristic transition width (see [`GradingLaw`] for the exact
    /// meaning per law). Always finite and strictly positive.
    pub width: f64,
    /// Material-parameter value approached for `position ≪ center`.
    pub value_left: f64,
    /// Material-parameter value approached for `position ≫ center`.
    pub value_right: f64,
    /// Per-site random offsets of the local interface `center`, modeling
    /// roughness-induced lateral fluctuation of the graded transition.
    /// Always has length `num_sites` (as passed to `generate`); every
    /// entry is exactly `0.0` when `jitter_rms == 0.0`.
    pub center_offsets: Vec<f64>,
    /// RMS amplitude of the center-position jitter (`0.0` = perfectly
    /// flat, laterally uniform interface).
    pub jitter_rms: f64,
    /// PRNG seed used to generate `center_offsets`.
    pub seed: u64,
}

impl GradedInterface {
    /// Generate a graded interface profile.
    ///
    /// # Arguments
    ///
    /// * `law` - Grading law (see [`GradingLaw`])
    /// * `center` - Nominal transition midpoint position (must be finite)
    /// * `width` - Characteristic transition width (must be finite and > 0)
    /// * `value_left` - Endpoint value for `position ≪ center` (must be finite)
    /// * `value_right` - Endpoint value for `position ≫ center` (must be finite)
    /// * `num_sites` - Number of lateral sites for interfacial-roughness
    ///   realizations (must be > 0; pass `1` if no per-site variation is
    ///   needed)
    /// * `jitter_rms` - RMS roughness-induced jitter of the local interface
    ///   center (must be finite and ≥ 0; `0.0` disables jitter)
    /// * `seed` - PRNG seed for reproducibility
    ///
    /// # Errors
    ///
    /// Returns `Error::InvalidParameter` if `width` is not finite and
    /// strictly positive, if `center`, `value_left`, `value_right`, or
    /// `jitter_rms` is not finite, if `jitter_rms` is negative, or if
    /// `num_sites == 0`.
    pub fn generate(
        law: GradingLaw,
        center: f64,
        width: f64,
        value_left: f64,
        value_right: f64,
        num_sites: usize,
        jitter_rms: f64,
        seed: u64,
    ) -> Result<Self> {
        if !center.is_finite() {
            return Err(Error::InvalidParameter {
                param: "center".to_string(),
                reason: "must be finite".to_string(),
            });
        }
        if !(width.is_finite() && width > 0.0) {
            return Err(Error::InvalidParameter {
                param: "width".to_string(),
                reason: "must be finite and strictly positive".to_string(),
            });
        }
        if !value_left.is_finite() {
            return Err(Error::InvalidParameter {
                param: "value_left".to_string(),
                reason: "must be finite".to_string(),
            });
        }
        if !value_right.is_finite() {
            return Err(Error::InvalidParameter {
                param: "value_right".to_string(),
                reason: "must be finite".to_string(),
            });
        }
        if num_sites == 0 {
            return Err(Error::InvalidParameter {
                param: "num_sites".to_string(),
                reason: "must be greater than zero".to_string(),
            });
        }
        if !(jitter_rms.is_finite() && jitter_rms >= 0.0) {
            return Err(Error::InvalidParameter {
                param: "jitter_rms".to_string(),
                reason: "must be finite and non-negative".to_string(),
            });
        }

        let mut rng = Xorshift64::new(seed);
        let center_offsets: Vec<f64> = (0..num_sites)
            .map(|_| jitter_rms * rng.next_normal())
            .collect();

        Ok(Self {
            law,
            center,
            width,
            value_left,
            value_right,
            center_offsets,
            jitter_rms,
            seed,
        })
    }

    /// Sample the graded material-parameter value at `position`, using the
    /// nominal (unperturbed) interface `center`.
    ///
    /// Implements the grading law selected at construction time; see
    /// [`GradingLaw`] for the exact functional form of each law.
    pub fn value_at(&self, position: f64) -> f64 {
        self.value_at_center(position, self.center)
    }

    /// Sample the graded profile at `position` for lateral site
    /// `site_index`, using the (possibly roughness-jittered) local
    /// interface center at that site.
    ///
    /// An out-of-range `site_index` (`>= self.center_offsets.len()`) is
    /// treated as an unperturbed (zero-jitter) sample, gracefully falling
    /// back to the nominal `center`.
    pub fn value_at_site(&self, position: f64, site_index: usize) -> f64 {
        let local_center =
            self.center + self.center_offsets.get(site_index).copied().unwrap_or(0.0);
        self.value_at_center(position, local_center)
    }

    /// Mean of the two endpoint values: the exact barycenter of any
    /// position window that is symmetric about `center`, for every law.
    pub fn barycenter(&self) -> f64 {
        0.5 * (self.value_left + self.value_right)
    }

    /// Core grading-law evaluation about an arbitrary `center` position
    /// (shared by [`Self::value_at`] and [`Self::value_at_site`]).
    fn value_at_center(&self, position: f64, center: f64) -> f64 {
        let d = position - center;
        let mean = 0.5 * (self.value_left + self.value_right);
        let half_delta = 0.5 * (self.value_right - self.value_left);

        match self.law {
            GradingLaw::Linear => {
                let t = (d / (0.5 * self.width)).clamp(-1.0, 1.0);
                mean + half_delta * t
            },
            GradingLaw::Exponential => {
                if d >= 0.0 {
                    self.value_right - half_delta * (-d / self.width).exp()
                } else {
                    self.value_left + half_delta * (d / self.width).exp()
                }
            },
            GradingLaw::ErrorFunction => mean + half_delta * erf_approx(d / self.width),
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xorshift64_deterministic() {
        let mut rng1 = Xorshift64::new(12345);
        let mut rng2 = Xorshift64::new(12345);
        for _ in 0..100 {
            assert_eq!(rng1.next_u64(), rng2.next_u64());
        }
    }

    #[test]
    fn test_exchange_length_calculation() {
        // For A = 1e-11 J/m, K = 1e4 J/m³:
        // L_ex = sqrt(A/K) = sqrt(1e-11 / 1e4) = sqrt(1e-15) = ~31.6 nm
        let ram =
            RandomAnisotropyModel::new(10, 10.0e-9, 1.0e4, 1.0e-11, 42).expect("valid parameters");
        let l_ex = ram.exchange_length();
        let expected = (1.0e-11_f64 / 1.0e4_f64).sqrt();
        assert!(
            (l_ex - expected).abs() / expected < 1e-10,
            "Exchange length mismatch: got {}, expected {}",
            l_ex,
            expected
        );
    }

    #[test]
    fn test_herzer_scaling_k_eff_less_than_bulk() {
        // With small grains (D << L_ex), K_eff should be much less than K₁
        let ram =
            RandomAnisotropyModel::new(1000, 5.0e-9, 1.0e4, 1.0e-11, 99).expect("valid parameters");

        // L_ex ~ 31.6 nm >> D = 5 nm, so N >> 1
        assert!(
            ram.effective_anisotropy < ram.bulk_anisotropy,
            "K_eff ({}) should be less than K₁ ({}) for small grains",
            ram.effective_anisotropy,
            ram.bulk_anisotropy
        );
    }

    #[test]
    fn test_herzer_scaling_proportionality() {
        // K_eff = K₁ / sqrt(N) where N = (L_ex/D)³
        // For two different grain sizes, check the scaling:
        // K_eff ∝ 1/sqrt((L_ex/D)³) = D^(3/2) / L_ex^(3/2)
        // So K_eff(D2)/K_eff(D1) = (D2/D1)^(3/2) when D << L_ex
        let d1 = 5.0e-9;
        let d2 = 10.0e-9;
        let k1_bulk = 1.0e4;
        let a = 1.0e-11;

        let ram1 = RandomAnisotropyModel::new(100, d1, k1_bulk, a, 1).expect("valid parameters");
        let ram2 = RandomAnisotropyModel::new(100, d2, k1_bulk, a, 1).expect("valid parameters");

        let ratio = ram2.effective_anisotropy / ram1.effective_anisotropy;
        let expected_ratio = (d2 / d1).powf(1.5); // D^(3/2) scaling

        assert!(
            (ratio - expected_ratio).abs() / expected_ratio < 0.01,
            "Herzer scaling ratio: got {}, expected {}",
            ratio,
            expected_ratio
        );
    }

    #[test]
    fn test_random_anisotropy_axes_are_unit_vectors() {
        let ram =
            RandomAnisotropyModel::new(200, 10.0e-9, 5.0e4, 1.0e-11, 7).expect("valid parameters");

        for (i, axis) in ram.grain_axes.iter().enumerate() {
            let norm = (axis.x * axis.x + axis.y * axis.y + axis.z * axis.z).sqrt();
            assert!(
                (norm - 1.0).abs() < 1e-10,
                "Grain axis {} has norm {}, expected 1.0",
                i,
                norm
            );
        }
    }

    #[test]
    fn test_grain_structure_generation() {
        let gs = GrainStructure::generate(
            50,
            Vector3::new(100.0e-9, 100.0e-9, 20.0e-9),
            10.0e-9,
            2.0e-9,
            42,
        )
        .expect("valid parameters");

        assert_eq!(gs.grains.len(), 50);

        // All grain centers should be in [0, 1)
        for grain in &gs.grains {
            assert!(grain.center.x >= 0.0 && grain.center.x < 1.0);
            assert!(grain.center.y >= 0.0 && grain.center.y < 1.0);
            assert!(grain.center.z >= 0.0 && grain.center.z < 1.0);
            assert!(grain.diameter > 0.0);
        }
    }

    #[test]
    fn test_nearest_grain_assignment() {
        let gs = GrainStructure::generate(
            10,
            Vector3::new(100.0e-9, 100.0e-9, 100.0e-9),
            20.0e-9,
            0.0,
            123,
        )
        .expect("valid parameters");

        // A point exactly at a grain center should map to that grain
        let center = gs.grains[0].center;
        let idx = gs.nearest_grain(center).expect("should find nearest grain");
        assert_eq!(idx, 0);
    }

    #[test]
    fn test_random_field_disorder() {
        let rfd = RandomFieldDisorder::generate(500, 0.01, 0.0, 55).expect("valid parameters");

        assert_eq!(rfd.fields.len(), 500);

        // The mean magnitude of a 3D Gaussian with σ=0.01 should be
        // roughly σ·√(8/π) ≈ 0.016 (chi distribution with k=3)
        let mean_mag = rfd.mean_magnitude();
        let expected_order = 0.01 * (8.0_f64 / std::f64::consts::PI).sqrt();
        assert!(
            (mean_mag - expected_order).abs() / expected_order < 0.3,
            "Mean random field magnitude {} deviates too far from expected {}",
            mean_mag,
            expected_order
        );
    }

    #[test]
    fn test_surface_roughness_rms() {
        let sr =
            SurfaceRoughness::generate(64, 64, 1.0e-9, 0.5e-9, 0.0, 77).expect("valid parameters");

        let actual = sr.actual_rms();
        let target = 0.5e-9;
        // With no correlation smoothing, the actual RMS should be close to target
        assert!(
            (actual - target).abs() / target < 0.15,
            "RMS roughness {} deviates too far from target {}",
            actual,
            target
        );
    }

    #[test]
    fn test_surface_anisotropy_modification() {
        let sr = SurfaceRoughness::generate(32, 32, 1.0e-9, 0.3e-9, 5.0e-9, 88)
            .expect("valid parameters");

        let k_s = 1.0e-3; // J/m²
        let delta_k = sr.anisotropy_modification(k_s);

        // ΔK_s = K_s · (σ/ξ)² = 1e-3 · (0.3e-9 / 5e-9)² = 1e-3 · 0.0036 = 3.6e-6
        let expected = k_s * (0.3e-9 / 5.0e-9) * (0.3e-9 / 5.0e-9);
        assert!(
            (delta_k - expected).abs() / expected < 1e-10,
            "Anisotropy modification {} != expected {}",
            delta_k,
            expected
        );
    }

    #[test]
    fn test_invalid_parameters() {
        assert!(RandomAnisotropyModel::new(0, 10.0e-9, 1.0e4, 1.0e-11, 1).is_err());
        assert!(RandomAnisotropyModel::new(10, -1.0, 1.0e4, 1.0e-11, 1).is_err());
        assert!(RandomAnisotropyModel::new(10, 10.0e-9, -1.0, 1.0e-11, 1).is_err());
        assert!(RandomAnisotropyModel::new(10, 10.0e-9, 1.0e4, -1.0, 1).is_err());
        assert!(RandomFieldDisorder::generate(0, 0.01, 0.0, 1).is_err());
        assert!(SurfaceRoughness::generate(0, 10, 1.0e-9, 0.5e-9, 0.0, 1).is_err());
    }

    // ------------------------------------------------------------------
    // Graded interfaces
    // ------------------------------------------------------------------

    #[test]
    fn test_graded_interface_linear_monotonic_and_endpoints() {
        let profile =
            GradedInterface::generate(GradingLaw::Linear, 0.0, 4.0e-9, 1.0e5, 3.0e5, 1, 0.0, 1)
                .expect("valid parameters");

        // Exact endpoint values once |position - center| >= width / 2
        assert!((profile.value_at(-10.0e-9) - 1.0e5).abs() < 1e-6);
        assert!((profile.value_at(10.0e-9) - 3.0e5).abs() < 1e-6);
        // Exact midpoint value
        assert!((profile.value_at(0.0) - 2.0e5).abs() < 1e-6);

        let mut prev = profile.value_at(-10.0e-9);
        for i in 1..=400 {
            let x = -10.0e-9 + i as f64 * (20.0e-9 / 400.0);
            let v = profile.value_at(x);
            assert!(
                v >= prev - 1e-12,
                "Linear profile not monotonic at x = {}: {} < {}",
                x,
                v,
                prev
            );
            prev = v;
        }
    }

    #[test]
    fn test_graded_interface_exponential_monotonic_and_endpoints() {
        let width = 2.0e-9;
        let profile =
            GradedInterface::generate(GradingLaw::Exponential, 0.0, width, 1.0e5, 3.0e5, 1, 0.0, 2)
                .expect("valid parameters");

        let far = 10.0 * width;
        let v_left = profile.value_at(-far);
        let v_right = profile.value_at(far);
        assert!(
            (v_left - 1.0e5).abs() / 1.0e5 < 1e-3,
            "Exponential left endpoint {} != 1e5",
            v_left
        );
        assert!(
            (v_right - 3.0e5).abs() / 3.0e5 < 1e-3,
            "Exponential right endpoint {} != 3e5",
            v_right
        );

        let mut prev = profile.value_at(-far);
        for i in 1..=400 {
            let x = -far + i as f64 * (2.0 * far / 400.0);
            let v = profile.value_at(x);
            assert!(
                v >= prev - 1e-9,
                "Exponential profile not monotonic at x = {}: {} < {}",
                x,
                v,
                prev
            );
            prev = v;
        }
    }

    #[test]
    fn test_graded_interface_error_function_monotonic_and_endpoints() {
        let width = 1.5e-9;
        let profile = GradedInterface::generate(
            GradingLaw::ErrorFunction,
            0.0,
            width,
            2.0e5,
            6.0e5,
            1,
            0.0,
            3,
        )
        .expect("valid parameters");

        let far = 10.0 * width;
        let v_left = profile.value_at(-far);
        let v_right = profile.value_at(far);
        assert!((v_left - 2.0e5).abs() < 1.0, "erf left endpoint {}", v_left);
        assert!(
            (v_right - 6.0e5).abs() < 1.0,
            "erf right endpoint {}",
            v_right
        );

        let mut prev = profile.value_at(-far);
        for i in 1..=400 {
            let x = -far + i as f64 * (2.0 * far / 400.0);
            let v = profile.value_at(x);
            assert!(
                v >= prev - 1e-9,
                "ErrorFunction profile not monotonic at x = {}: {} < {}",
                x,
                v,
                prev
            );
            prev = v;
        }
    }

    #[test]
    fn test_graded_interface_barycenter_symmetric_profile() {
        // For every law, value(center+d) + value(center-d) == value_left +
        // value_right exactly (to floating-point precision), independent of
        // d and width -- this is what makes the mean over any
        // center-symmetric sampling window equal the barycenter exactly.
        for law in [
            GradingLaw::Linear,
            GradingLaw::Exponential,
            GradingLaw::ErrorFunction,
        ] {
            let profile = GradedInterface::generate(law, 5.0e-9, 3.0e-9, 1.0e5, 7.0e5, 1, 0.0, 7)
                .expect("valid parameters");
            let expected_sum = profile.value_left + profile.value_right;
            assert!(
                (profile.barycenter() - 4.0e5).abs() < 1e-6,
                "barycenter mismatch for {:?}",
                law
            );

            for &d in &[0.0, 0.5e-9, 1.0e-9, 3.0e-9, 6.0e-9, 15.0e-9, 30.0e-9] {
                let v_plus = profile.value_at(profile.center + d);
                let v_minus = profile.value_at(profile.center - d);
                let sum = v_plus + v_minus;
                assert!(
                    (sum - expected_sum).abs() / expected_sum.abs() < 1e-6,
                    "{:?}: value(c+{}) + value(c-{}) = {}, expected {}",
                    law,
                    d,
                    d,
                    sum,
                    expected_sum
                );
            }

            // Discretised-mean check over a domain symmetric about the center.
            let n = 201usize; // odd => the sampling grid includes the center itself
            let half_span = 20.0 * profile.width;
            let mut sum = 0.0;
            for i in 0..n {
                let frac = i as f64 / (n - 1) as f64; // 0..1
                let x = profile.center - half_span + frac * 2.0 * half_span;
                sum += profile.value_at(x);
            }
            let mean = sum / n as f64;
            assert!(
                (mean - profile.barycenter()).abs() / profile.barycenter().abs() < 1e-3,
                "{:?}: discretised mean {} != barycenter {}",
                law,
                mean,
                profile.barycenter()
            );
        }
    }

    #[test]
    fn test_graded_interface_invalid_parameters() {
        // width must be finite and strictly positive
        assert!(
            GradedInterface::generate(GradingLaw::Linear, 0.0, 0.0, 0.0, 1.0, 1, 0.0, 1).is_err()
        );
        assert!(
            GradedInterface::generate(GradingLaw::Linear, 0.0, -1.0e-9, 0.0, 1.0, 1, 0.0, 1)
                .is_err()
        );
        assert!(
            GradedInterface::generate(GradingLaw::Linear, 0.0, f64::NAN, 0.0, 1.0, 1, 0.0, 1)
                .is_err()
        );
        // endpoints must be finite
        assert!(GradedInterface::generate(
            GradingLaw::Linear,
            0.0,
            1.0e-9,
            f64::NAN,
            1.0,
            1,
            0.0,
            1
        )
        .is_err());
        assert!(GradedInterface::generate(
            GradingLaw::Linear,
            0.0,
            1.0e-9,
            0.0,
            f64::INFINITY,
            1,
            0.0,
            1
        )
        .is_err());
        // center must be finite
        assert!(GradedInterface::generate(
            GradingLaw::Linear,
            f64::NAN,
            1.0e-9,
            0.0,
            1.0,
            1,
            0.0,
            1
        )
        .is_err());
        // num_sites must be > 0
        assert!(
            GradedInterface::generate(GradingLaw::Linear, 0.0, 1.0e-9, 0.0, 1.0, 0, 0.0, 1)
                .is_err()
        );
        // jitter_rms must be finite and non-negative
        assert!(
            GradedInterface::generate(GradingLaw::Linear, 0.0, 1.0e-9, 0.0, 1.0, 1, -1.0, 1)
                .is_err()
        );
        assert!(GradedInterface::generate(
            GradingLaw::Linear,
            0.0,
            1.0e-9,
            0.0,
            1.0,
            1,
            f64::NAN,
            1
        )
        .is_err());
        // sanity: a valid call succeeds
        assert!(
            GradedInterface::generate(GradingLaw::Linear, 0.0, 1.0e-9, 0.0, 1.0, 1, 0.0, 1).is_ok()
        );
    }

    #[test]
    fn test_graded_interface_roughness_jitter_deterministic_and_scaled() {
        let n_sites = 200usize;
        let jitter_rms = 0.5e-9;
        let p1 = GradedInterface::generate(
            GradingLaw::ErrorFunction,
            0.0,
            2.0e-9,
            1.0e5,
            3.0e5,
            n_sites,
            jitter_rms,
            123,
        )
        .expect("valid parameters");
        let p2 = GradedInterface::generate(
            GradingLaw::ErrorFunction,
            0.0,
            2.0e-9,
            1.0e5,
            3.0e5,
            n_sites,
            jitter_rms,
            123,
        )
        .expect("valid parameters");

        // Same seed => bit-identical jitter realization
        assert_eq!(p1.center_offsets.len(), n_sites);
        for (a, b) in p1.center_offsets.iter().zip(p2.center_offsets.iter()) {
            assert_eq!(a, b);
        }

        // RMS of the realized offsets should be close to the requested jitter_rms
        let mean_sq: f64 = p1.center_offsets.iter().map(|o| o * o).sum::<f64>() / n_sites as f64;
        let actual_rms = mean_sq.sqrt();
        assert!(
            (actual_rms - jitter_rms).abs() / jitter_rms < 0.3,
            "actual jitter RMS {} deviates too far from requested {}",
            actual_rms,
            jitter_rms
        );

        // An out-of-range site index gracefully falls back to the nominal center
        let fallback = p1.value_at_site(0.0, n_sites + 1000);
        let nominal = p1.value_at(0.0);
        assert!((fallback - nominal).abs() < 1e-9);
    }
}
