//! Random Anisotropy Magnet (RAM) — Imry-Ma disorder physics.
//!
//! Implements per-site uniaxial anisotropy with random easy-axis directions:
//!
//! $$H_A = -\sum_i K_i (\hat{n}_i \cdot \mathbf{S}_i)^2$$
//!
//! where K_i is the local anisotropy strength and n̂_i is a random unit vector.
//!
//! # Physical Background
//!
//! In the Imry-Ma (1975) theory, random anisotropy destroys long-range ferromagnetic
//! order in dimensions d < 4. The correlation length scales as:
//!
//! $$\xi_M = \xi_a \left(\frac{A}{K_{\rm rms}}\right)^2$$
//!
//! where A is the exchange stiffness and ξ_a is the grain (correlation) length.
//!
//! The Harris criterion states disorder is relevant when the sample size L > ξ_a.
//!
//! # References
//! - Y. Imry, S. K. Ma, Phys. Rev. Lett. 35, 1399 (1975)
//! - A. B. Harris, J. Phys. C 7, 1671 (1974)
//! - E. M. Chudnovsky, W. M. Saslow, R. A. Serota, Phys. Rev. B 33, 251 (1986)
//! - D. J. Sellmyer, S. Nafis, J. Appl. Phys. 57, 3584 (1985)

use crate::error::{self, Result};
use crate::material::disorder::Xorshift64;
use crate::vector3::Vector3;

/// Vacuum permeability μ₀ [T·m/A].
const MU_0: f64 = 1.256_637_062_12e-6;

// ============================================================================
// Distribution type
// ============================================================================

/// Distribution type for random anisotropy strengths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum RandomAnisotropyDistribution {
    /// Gaussian distribution: K_i ~ N(k_mean, k_std²)
    Gaussian,
    /// Uniform distribution: K_i ~ U(k_mean - k_std√3, k_mean + k_std√3)
    /// This gives variance = k_std² (matching Gaussian variance).
    Uniform,
    /// Fixed easy axes (user-supplied), only strengths are random.
    /// Uses Gaussian sampling for the strength magnitudes.
    FixedAxes,
}

// ============================================================================
// RandomAnisotropy struct
// ============================================================================

/// Random anisotropy model for disordered magnetic materials.
///
/// Generates random easy-axis directions and anisotropy strengths for each site
/// in a spin system, following the Imry-Ma random anisotropy model.
///
/// # Example
///
/// ```rust
/// use spintronics::material::random_anisotropy::{RandomAnisotropy, RandomAnisotropyDistribution};
///
/// let mut model = RandomAnisotropy::new(
///     1e3,        // k_mean [J/m³]
///     1e2,        // k_std  [J/m³]
///     5e-9,       // correlation_length [m]
///     RandomAnisotropyDistribution::Gaussian,
///     42,
/// ).expect("valid parameters");
///
/// let n_sites = 100;
/// let axes      = model.generate_axes(n_sites);
/// let strengths = model.generate_strengths(n_sites);
/// ```
#[derive(Debug, Clone)]
pub struct RandomAnisotropy {
    /// Mean anisotropy strength [J/m³].
    pub k_mean: f64,
    /// Standard deviation of anisotropy strength [J/m³].
    pub k_std: f64,
    /// Grain (correlation) length ξ_a [m].
    pub correlation_length: f64,
    /// Distribution of anisotropy strengths.
    pub distribution: RandomAnisotropyDistribution,
    /// Random seed for reproducibility.
    pub seed: u64,
    /// Internal RNG state.
    rng: Xorshift64,
}

impl RandomAnisotropy {
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Create a new random anisotropy model.
    ///
    /// # Arguments
    /// * `k_mean` — mean anisotropy strength [J/m³]; negative values are
    ///   allowed (unusual anisotropy) but zero is fine too.
    /// * `k_std` — standard deviation of strength [J/m³]; must be ≥ 0.
    /// * `correlation_length` — grain size ξ_a [m]; must be > 0.
    /// * `distribution` — how strengths are distributed across sites.
    /// * `seed` — RNG seed for reproducibility.
    ///
    /// # Errors
    /// Returns an error if `k_std < 0` or `correlation_length <= 0`.
    pub fn new(
        k_mean: f64,
        k_std: f64,
        correlation_length: f64,
        distribution: RandomAnisotropyDistribution,
        seed: u64,
    ) -> Result<Self> {
        if k_std < 0.0 {
            return Err(error::invalid_param(
                "k_std",
                "anisotropy strength standard deviation must be non-negative",
            ));
        }
        if correlation_length <= 0.0 {
            return Err(error::invalid_param(
                "correlation_length",
                "grain correlation length must be strictly positive",
            ));
        }

        Ok(Self {
            k_mean,
            k_std,
            correlation_length,
            distribution,
            seed,
            rng: Xorshift64::new(seed),
        })
    }

    /// Default model for a nanocrystalline soft magnet (e.g., permalloy-like).
    ///
    /// Parameters: k_mean = 1 kJ/m³, k_std = 100 J/m³,
    /// correlation_length = 5 nm, Gaussian distribution, seed = 42.
    pub fn nanocrystalline() -> Self {
        // Safety: parameters are valid by construction.
        Self {
            k_mean: 1e3,
            k_std: 1e2,
            correlation_length: 5e-9,
            distribution: RandomAnisotropyDistribution::Gaussian,
            seed: 42,
            rng: Xorshift64::new(42),
        }
    }

    // -----------------------------------------------------------------------
    // Generation
    // -----------------------------------------------------------------------

    /// Generate `n` random easy-axis unit vectors with isotropic distribution.
    ///
    /// Uses the Marsaglia sphere-point picking method (via [`Xorshift64::next_unit_vector`])
    /// so that all orientations on the unit sphere are equally likely.
    pub fn generate_axes(&mut self, n: usize) -> Vec<Vector3<f64>> {
        (0..n).map(|_| self.rng.next_unit_vector()).collect()
    }

    /// Generate `n` random anisotropy strength values.
    ///
    /// The sampling strategy depends on the chosen distribution:
    ///
    /// * **Gaussian** — K_i = k_mean + k_std · z  where z ~ N(0, 1).
    /// * **Uniform**  — K_i = k_mean + k_std · √3 · (2u − 1)  where u ~ U[0,1).
    ///   This parameterisation preserves the variance k_std² so that Gaussian
    ///   and Uniform runs are comparable in their second moment.
    /// * **FixedAxes** — identical to Gaussian (easy axes are supplied separately
    ///   by the caller via `generate_axes` or a user-provided list).
    pub fn generate_strengths(&mut self, n: usize) -> Vec<f64> {
        let sqrt3 = 3.0_f64.sqrt();
        match self.distribution {
            RandomAnisotropyDistribution::Gaussian | RandomAnisotropyDistribution::FixedAxes => {
                let mut out = Vec::with_capacity(n);
                let mut remaining = n;
                while remaining > 0 {
                    let (z0, z1) = self.rng.next_normal_pair();
                    out.push(self.k_mean + self.k_std * z0);
                    remaining -= 1;
                    if remaining > 0 {
                        out.push(self.k_mean + self.k_std * z1);
                        remaining -= 1;
                    }
                }
                out
            },
            RandomAnisotropyDistribution::Uniform => (0..n)
                .map(|_| {
                    let u = self.rng.next_f64();
                    // (2u − 1) ∈ (−1, 1), so K_i ∈ (k_mean − k_std√3, k_mean + k_std√3)
                    self.k_mean + self.k_std * sqrt3 * (2.0 * u - 1.0)
                })
                .collect(),
        }
    }

    // -----------------------------------------------------------------------
    // Correlation & length scales
    // -----------------------------------------------------------------------

    /// Spatial autocorrelation function C(r) = exp(−r / ξ_a).
    ///
    /// Simple exponential decay model for random anisotropy correlations,
    /// consistent with an Ornstein-Uhlenbeck disorder field.
    pub fn correlation_function(&self, r: f64) -> f64 {
        (-r / self.correlation_length).exp()
    }

    /// Imry-Ma magnetic correlation length:
    ///
    /// ξ_M = ξ_a · (A / K_rms)²
    ///
    /// where K_rms = √(k_mean² + k_std²) is the RMS anisotropy and A is the
    /// exchange stiffness.  This length scale characterises the size of
    /// correlated spin clusters in the random-anisotropy phase.
    ///
    /// # Arguments
    /// * `exchange_a` — exchange stiffness A [J/m], typically ~1×10⁻¹¹ J/m
    ///   for soft ferromagnets.
    pub fn imry_ma_correlation_length(&self, exchange_a: f64) -> f64 {
        let k_rms = self.k_rms();
        if k_rms < f64::EPSILON {
            // No disorder → infinite magnetic correlation length (pure ferromagnet).
            f64::INFINITY
        } else {
            let ratio = exchange_a / k_rms;
            self.correlation_length * ratio * ratio
        }
    }

    /// Harris criterion: disorder is relevant when sample_size > ξ_a.
    ///
    /// Returns `true` when the sample is large enough that random-anisotropy
    /// fluctuations significantly renormalise the effective anisotropy.
    pub fn harris_relevant(&self, sample_size: f64) -> bool {
        sample_size > self.correlation_length
    }

    // -----------------------------------------------------------------------
    // Energetics
    // -----------------------------------------------------------------------

    /// Total uniaxial random-anisotropy energy for a spin configuration:
    ///
    /// E_A = −Σ_i K_i (n̂_i · m_i)²
    ///
    /// where the sum runs over all sites.  Slices must all have the same length.
    ///
    /// # Arguments
    /// * `magnetization` — unit magnetisation vector at each site.
    /// * `axes`          — random easy-axis direction at each site (unit vectors).
    /// * `strengths`     — anisotropy strength K_i at each site [J/m³].
    pub fn anisotropy_energy(
        &self,
        magnetization: &[Vector3<f64>],
        axes: &[Vector3<f64>],
        strengths: &[f64],
    ) -> f64 {
        let n = magnetization.len().min(axes.len()).min(strengths.len());
        let mut energy = 0.0_f64;
        for i in 0..n {
            let cos_theta = axes[i].dot(&magnetization[i]);
            energy -= strengths[i] * cos_theta * cos_theta;
        }
        energy
    }

    /// Effective anisotropy field at each site for use in the LLG equation:
    ///
    /// **H**_eff,i = (2 K_i / (μ₀ Ms)) (n̂_i · m_i) n̂_i
    ///
    /// This is the negative functional derivative of the anisotropy energy
    /// density with respect to the magnetisation, divided by μ₀ Ms:
    ///
    /// H_eff = −(1 / μ₀ Ms) ∂E/∂m = +(2K / μ₀ Ms)(n̂·m) n̂
    ///
    /// The result is in [A/m].
    ///
    /// # Arguments
    /// * `magnetization` — unit magnetisation vector at each site.
    /// * `axes`          — easy-axis unit vectors (from [`generate_axes`]).
    /// * `strengths`     — anisotropy strengths K_i [J/m³] (from [`generate_strengths`]).
    /// * `ms`            — saturation magnetisation [A/m]; must be > 0.
    ///
    /// # Errors
    /// Returns an error when `ms ≤ 0`.
    ///
    /// [`generate_axes`]: RandomAnisotropy::generate_axes
    pub fn effective_field(
        &self,
        magnetization: &[Vector3<f64>],
        axes: &[Vector3<f64>],
        strengths: &[f64],
        ms: f64,
    ) -> Result<Vec<Vector3<f64>>> {
        if ms <= 0.0 {
            return Err(error::invalid_param(
                "ms",
                "saturation magnetisation must be strictly positive",
            ));
        }

        let prefactor_denom = MU_0 * ms;
        let n = magnetization.len().min(axes.len()).min(strengths.len());
        let mut fields = Vec::with_capacity(n);

        for i in 0..n {
            let cos_theta = axes[i].dot(&magnetization[i]);
            // H_eff,i = (2 K_i / (μ₀ Ms)) (n̂_i · m_i) n̂_i
            let scalar = 2.0 * strengths[i] * cos_theta / prefactor_denom;
            fields.push(Vector3::new(
                scalar * axes[i].x,
                scalar * axes[i].y,
                scalar * axes[i].z,
            ));
        }

        Ok(fields)
    }

    // -----------------------------------------------------------------------
    // Derived quantities
    // -----------------------------------------------------------------------

    /// RMS anisotropy strength: K_rms = √(k_mean² + k_std²).
    ///
    /// This is the second moment of the distribution about zero and enters
    /// the Imry-Ma correlation-length formula.
    pub fn k_rms(&self) -> f64 {
        (self.k_mean * self.k_mean + self.k_std * self.k_std).sqrt()
    }

    /// Reset the internal RNG to the original seed.
    ///
    /// Useful for reproducible experiments where the same disorder
    /// configuration must be regenerated without constructing a new model.
    pub fn reset_rng(&mut self) {
        self.rng = Xorshift64::new(self.seed);
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build a default model with Gaussian distribution and seed 7.
    fn make_gaussian(seed: u64) -> RandomAnisotropy {
        RandomAnisotropy::new(
            1e3,  // k_mean  [J/m³]
            1e2,  // k_std   [J/m³]
            5e-9, // xi_a    [m]
            RandomAnisotropyDistribution::Gaussian,
            seed,
        )
        .expect("parameters are valid")
    }

    // -------------------------------------------------------------------------
    // 1. Isotropy of axis distribution
    // -------------------------------------------------------------------------
    #[test]
    fn test_isotropic_axis_average() {
        let mut model = make_gaussian(1);
        let n = 10_000_usize;
        let axes = model.generate_axes(n);

        let sum_x: f64 = axes.iter().map(|v| v.x).sum();
        let sum_y: f64 = axes.iter().map(|v| v.y).sum();
        let sum_z: f64 = axes.iter().map(|v| v.z).sum();

        let mean_x = sum_x / n as f64;
        let mean_y = sum_y / n as f64;
        let mean_z = sum_z / n as f64;

        let net_mag = (mean_x * mean_x + mean_y * mean_y + mean_z * mean_z).sqrt();

        // For a truly isotropic distribution the net magnetisation should vanish.
        // With N = 10 000 and unit vectors the CLT gives std ~ 1/√N ≈ 0.01,
        // so |net| < 0.05 is well within 5σ.
        assert!(
            net_mag < 0.05,
            "axis distribution is not isotropic: |mean| = {net_mag}"
        );
    }

    // -------------------------------------------------------------------------
    // 2. Gaussian strength mean convergence
    // -------------------------------------------------------------------------
    #[test]
    fn test_gaussian_strengths_mean() {
        let mut model = make_gaussian(2);
        let n = 10_000_usize;
        let k_std = model.k_std;
        let k_mean = model.k_mean;

        let strengths = model.generate_strengths(n);
        let sample_mean: f64 = strengths.iter().sum::<f64>() / n as f64;

        // CLT: std of sample mean = k_std / sqrt(N)
        let tolerance = 3.0 * k_std / (n as f64).sqrt();
        let deviation = (sample_mean - k_mean).abs();

        assert!(
            deviation < tolerance,
            "sample mean {sample_mean} deviates from k_mean {k_mean} by {deviation} > 3σ = {tolerance}"
        );
    }

    // -------------------------------------------------------------------------
    // 3. Correlation function is strictly monotone decreasing for r > 0
    // -------------------------------------------------------------------------
    #[test]
    fn test_correlation_decay_monotone() {
        let model = make_gaussian(3);
        // Test a range of positive radii.
        let radii: [f64; 6] = [1e-10, 1e-9, 5e-9, 1e-8, 5e-8, 1e-7];
        for &r in &radii {
            let c_r = model.correlation_function(r);
            let c_2r = model.correlation_function(2.0 * r);
            assert!(
                c_r > c_2r,
                "C({r}) = {c_r} should be > C({d}) = {c_2r}",
                d = 2.0 * r
            );
        }
    }

    // -------------------------------------------------------------------------
    // 4. Correlation function at zero equals 1
    // -------------------------------------------------------------------------
    #[test]
    fn test_correlation_at_zero_is_one() {
        let model = make_gaussian(4);
        let c0 = model.correlation_function(0.0);
        assert!((c0 - 1.0).abs() < 1e-15, "C(0) should be 1.0, got {c0}");
    }

    // -------------------------------------------------------------------------
    // 5. Imry-Ma correlation length is positive for A > 0
    // -------------------------------------------------------------------------
    #[test]
    fn test_imry_ma_length_positive() {
        let model = make_gaussian(5);
        // Try several exchange stiffness values.
        let a_values: [f64; 4] = [1e-13, 1e-12, 1e-11, 1e-10];
        for &exchange_a in &a_values {
            let xi_m = model.imry_ma_correlation_length(exchange_a);
            assert!(
                xi_m > 0.0,
                "Imry-Ma length should be positive for A = {exchange_a}, got {xi_m}"
            );
        }
    }

    // -------------------------------------------------------------------------
    // 6. Harris criterion boundary behaviour
    // -------------------------------------------------------------------------
    #[test]
    fn test_harris_relevance_boundary() {
        let xi_a = 5e-9_f64;
        let model =
            RandomAnisotropy::new(1e3, 1e2, xi_a, RandomAnisotropyDistribution::Gaussian, 6)
                .expect("valid params");

        // Slightly above ξ_a → disorder is relevant.
        let above = xi_a + 1e-12;
        assert!(
            model.harris_relevant(above),
            "disorder should be Harris-relevant for L = {above} > xi_a = {xi_a}"
        );

        // Slightly below ξ_a → disorder is irrelevant.
        let below = xi_a - 1e-12;
        assert!(
            !model.harris_relevant(below),
            "disorder should not be Harris-relevant for L = {below} < xi_a = {xi_a}"
        );
    }

    // -------------------------------------------------------------------------
    // 7. Anisotropy energy is minimised when m_i = n̂_i for all i
    // -------------------------------------------------------------------------
    #[test]
    fn test_anisotropy_energy_aligned() {
        let mut model = make_gaussian(7);
        let n = 50_usize;
        let axes = model.generate_axes(n);
        let strengths = model.generate_strengths(n);

        // When m_i = n̂_i we have (n̂_i · m_i)² = 1, so E = −Σ K_i.
        let e_aligned = model.anisotropy_energy(&axes, &axes, &strengths);
        let e_expected: f64 = -strengths.iter().sum::<f64>();

        let tol = 1e-10 * e_expected.abs().max(1.0);
        assert!(
            (e_aligned - e_expected).abs() < tol,
            "aligned energy {e_aligned} should equal −ΣK_i = {e_expected}"
        );

        // A perpendicular magnetisation should give E = 0 (higher energy).
        // Construct m_perp ⊥ n̂_i by rotating each axis 90° around an
        // arbitrary orthogonal direction.
        let perp_mag: Vec<Vector3<f64>> = axes
            .iter()
            .map(|n_hat| {
                // Find any vector not parallel to n_hat.
                let candidate = if n_hat.x.abs() < 0.9 {
                    Vector3::unit_x()
                } else {
                    Vector3::unit_z()
                };
                // Gram-Schmidt: remove component along n_hat.
                let proj_scale = n_hat.dot(&candidate);
                let perp = Vector3::new(
                    candidate.x - proj_scale * n_hat.x,
                    candidate.y - proj_scale * n_hat.y,
                    candidate.z - proj_scale * n_hat.z,
                );
                perp.normalize()
            })
            .collect();

        let e_perp = model.anisotropy_energy(&perp_mag, &axes, &strengths);
        // E_perp = −Σ K_i · 0 = 0
        let tol_perp = 1e-10 * (strengths.iter().map(|k| k.abs()).sum::<f64>() + 1.0);
        assert!(
            e_perp.abs() < tol_perp,
            "perpendicular energy should be ~0, got {e_perp}"
        );

        // Aligned configuration must have lower (more negative) energy.
        assert!(
            e_aligned <= e_perp,
            "aligned energy {e_aligned} should be ≤ perpendicular energy {e_perp}"
        );
    }

    // -------------------------------------------------------------------------
    // 8. Effective field is non-zero for partially aligned m and n̂
    // -------------------------------------------------------------------------
    #[test]
    fn test_effective_field_nonzero() {
        let mut model = make_gaussian(8);
        let n = 10_usize;
        let axes = model.generate_axes(n);
        let strengths = model.generate_strengths(n);

        // Use axes themselves as the magnetisation so (n̂·m) = 1 everywhere.
        let ms = 8.0e5_f64; // A/m  (typical permalloy Ms)
        let fields = model
            .effective_field(&axes, &axes, &strengths, ms)
            .expect("ms is positive");

        assert_eq!(fields.len(), n);

        // Every field should have a non-zero magnitude when K_i ≠ 0 and m ∥ n̂.
        // With k_mean = 1 kJ/m³ the field magnitude is ~2K/(μ₀ Ms) ≈ 2 kA/m.
        for (i, h) in fields.iter().enumerate() {
            let mag = (h.x * h.x + h.y * h.y + h.z * h.z).sqrt();
            assert!(
                mag > 0.0,
                "effective field at site {i} should be non-zero, got magnitude {mag}"
            );
        }

        // Validate that ms ≤ 0 returns an error.
        let err = model.effective_field(&axes, &axes, &strengths, 0.0);
        assert!(err.is_err(), "effective_field with ms=0 should return Err");

        let err_neg = model.effective_field(&axes, &axes, &strengths, -1.0);
        assert!(
            err_neg.is_err(),
            "effective_field with ms<0 should return Err"
        );
    }
}
