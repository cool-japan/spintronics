//! Point defects, pinning sites, and grain boundary models
//!
//! This module models structural imperfections in magnetic materials that
//! critically influence domain wall dynamics and coercivity:
//!
//! - **Point defects**: Vacancies and substitutional impurities that locally
//!   modify exchange coupling and anisotropy.
//! - **Pinning sites**: Localized potentials that trap domain walls, requiring
//!   a threshold field or current to depin.
//! - **Grain boundaries**: Reduced exchange coupling at inter-grain interfaces
//!   affecting domain wall propagation.
//! - **Concentration scaling**: Coercivity dependence on defect density
//!   following power-law scaling H_c ~ n^α.
//!
//! # References
//!
//! - A. Aharoni, "Introduction to the Theory of Ferromagnetism", Oxford (2000)
//! - H. Kronmüller, M. Fähnle, "Micromagnetism and the Microstructure of
//!   Ferromagnetic Solids", Cambridge (2003)

use crate::constants::MU_0;
use crate::error::{Error, Result};

// ============================================================================
// Defect types
// ============================================================================

/// Type of point defect at a lattice site.
#[derive(Debug, Clone)]
pub enum DefectType {
    /// Vacancy: the magnetic moment at this site is absent.
    Vacancy,
    /// Substitutional impurity with modified local properties.
    Substitutional {
        /// Local anisotropy constant \[J/m³\] (replaces bulk value).
        anisotropy: f64,
        /// Exchange coupling modification factor (1.0 = unchanged).
        exchange_mod: f64,
    },
    /// A localized pinning center for domain walls.
    Pinning {
        /// Pinning potential depth V₀ \[J/m²\].
        strength: f64,
        /// Spatial width of the pinning potential w \[m\].
        width: f64,
    },
}

/// A defect located at a specific lattice site.
#[derive(Debug, Clone)]
pub struct DefectSite {
    /// Index of the lattice site where the defect resides.
    pub position: usize,
    /// The type and parameters of the defect.
    pub defect_type: DefectType,
}

impl DefectSite {
    /// Create a vacancy defect at the given site index.
    pub fn vacancy(position: usize) -> Self {
        Self {
            position,
            defect_type: DefectType::Vacancy,
        }
    }

    /// Create a substitutional defect at the given site index.
    ///
    /// # Arguments
    ///
    /// * `position` - Lattice site index
    /// * `anisotropy` - Local anisotropy \[J/m³\]
    /// * `exchange_mod` - Exchange coupling modification factor
    pub fn substitutional(position: usize, anisotropy: f64, exchange_mod: f64) -> Self {
        Self {
            position,
            defect_type: DefectType::Substitutional {
                anisotropy,
                exchange_mod,
            },
        }
    }

    /// Create a pinning defect at the given site index.
    ///
    /// # Arguments
    ///
    /// * `position` - Lattice site index
    /// * `strength` - Pinning potential V₀ \[J/m²\]
    /// * `width` - Spatial width w \[m\]
    pub fn pinning(position: usize, strength: f64, width: f64) -> Self {
        Self {
            position,
            defect_type: DefectType::Pinning { strength, width },
        }
    }

    /// Check whether this defect is a vacancy.
    pub fn is_vacancy(&self) -> bool {
        matches!(self.defect_type, DefectType::Vacancy)
    }

    /// Return the effective local anisotropy for this site.
    ///
    /// - Vacancy → 0 (no magnetic moment, no anisotropy contribution)
    /// - Substitutional → the substitutional anisotropy value
    /// - Pinning → bulk value is unmodified, returns `bulk_anisotropy`
    pub fn effective_anisotropy(&self, bulk_anisotropy: f64) -> f64 {
        match &self.defect_type {
            DefectType::Vacancy => 0.0,
            DefectType::Substitutional { anisotropy, .. } => *anisotropy,
            DefectType::Pinning { .. } => bulk_anisotropy,
        }
    }

    /// Return the effective exchange coupling modification factor.
    ///
    /// - Vacancy → 0.0 (no exchange through vacancy)
    /// - Substitutional → the exchange_mod factor
    /// - Pinning → 1.0 (unmodified exchange)
    pub fn exchange_factor(&self) -> f64 {
        match &self.defect_type {
            DefectType::Vacancy => 0.0,
            DefectType::Substitutional { exchange_mod, .. } => *exchange_mod,
            DefectType::Pinning { .. } => 1.0,
        }
    }
}

// ============================================================================
// Pinning potential
// ============================================================================

/// Evaluate the Gaussian pinning potential at displacement x from the center.
///
/// V(x) = V₀ · exp(-(x - x₀)² / (2 w²))
///
/// # Arguments
///
/// * `x` - Position along the domain wall propagation direction \[m\]
/// * `x0` - Center of the pinning site \[m\]
/// * `v0` - Pinning potential depth V₀ \[J/m²\]
/// * `w` - Width of the pinning potential \[m\]
///
/// # Returns
///
/// The potential energy density at position x.
pub fn pinning_potential(x: f64, x0: f64, v0: f64, w: f64) -> f64 {
    if w.abs() < f64::EPSILON {
        // Delta-function limit: only nonzero exactly at x0
        return if (x - x0).abs() < f64::EPSILON {
            v0
        } else {
            0.0
        };
    }
    let dx = x - x0;
    v0 * (-dx * dx / (2.0 * w * w)).exp()
}

/// Compute the force from a Gaussian pinning potential (negative gradient).
///
/// F(x) = -dV/dx = V₀ · (x - x₀) / w² · exp(-(x-x₀)²/(2w²))
///
/// # Arguments
///
/// * `x` - Position \[m\]
/// * `x0` - Pinning center \[m\]
/// * `v0` - Pinning strength V₀ \[J/m²\]
/// * `w` - Pinning width \[m\]
pub fn pinning_force(x: f64, x0: f64, v0: f64, w: f64) -> f64 {
    if w.abs() < f64::EPSILON {
        return 0.0;
    }
    let dx = x - x0;
    v0 * dx / (w * w) * (-dx * dx / (2.0 * w * w)).exp()
}

// ============================================================================
// Depinning threshold
// ============================================================================

/// Parameters for computing the depinning field of a domain wall.
#[derive(Debug, Clone)]
pub struct DepinningParams {
    /// Pinning potential depth V₀ \[J/m²\].
    pub pinning_strength: f64,
    /// Saturation magnetization M_s \[A/m\].
    pub saturation_magnetization: f64,
    /// Domain wall width δ_w \[m\].
    pub domain_wall_width: f64,
}

impl DepinningParams {
    /// Create new depinning parameters.
    ///
    /// # Errors
    ///
    /// Returns an error if any parameter is non-positive.
    pub fn new(
        pinning_strength: f64,
        saturation_magnetization: f64,
        domain_wall_width: f64,
    ) -> Result<Self> {
        if pinning_strength <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "pinning_strength".to_string(),
                reason: "must be positive".to_string(),
            });
        }
        if saturation_magnetization <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "saturation_magnetization".to_string(),
                reason: "must be positive".to_string(),
            });
        }
        if domain_wall_width <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "domain_wall_width".to_string(),
                reason: "must be positive".to_string(),
            });
        }
        Ok(Self {
            pinning_strength,
            saturation_magnetization,
            domain_wall_width,
        })
    }

    /// Compute the depinning field H_dep \[A/m\].
    ///
    /// H_dep = V₀ / (2 μ₀ M_s δ_w)
    ///
    /// This is the minimum external field required to move a domain wall
    /// past the pinning center at zero temperature.
    pub fn depinning_field(&self) -> f64 {
        self.pinning_strength
            / (2.0 * MU_0 * self.saturation_magnetization * self.domain_wall_width)
    }

    /// Depinning field with thermal activation (Arrhenius-Néel model).
    ///
    /// At finite temperature the effective depinning field is reduced:
    ///
    ///   H_dep(T) = H_dep(0) · [1 - (k_B T / E_b)^(1/μ)]
    ///
    /// where E_b is the energy barrier and μ ≈ 1 for a simple Gaussian
    /// potential. The barrier is E_b = V₀ · δ_w.
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature \[K\] (must be >= 0)
    ///
    /// # Returns
    ///
    /// The thermally-reduced depinning field \[A/m\]. Returns 0 if thermal
    /// energy exceeds the barrier.
    pub fn depinning_field_thermal(&self, temperature: f64) -> Result<f64> {
        if temperature < 0.0 {
            return Err(Error::InvalidParameter {
                param: "temperature".to_string(),
                reason: "must be non-negative".to_string(),
            });
        }
        let h_dep_0 = self.depinning_field();
        if temperature < f64::EPSILON {
            return Ok(h_dep_0);
        }

        let kb = crate::constants::KB;
        let energy_barrier = self.pinning_strength * self.domain_wall_width;
        let thermal_ratio = kb * temperature / energy_barrier;

        if thermal_ratio >= 1.0 {
            // Thermal energy exceeds barrier — wall is thermally depinned
            Ok(0.0)
        } else {
            Ok(h_dep_0 * (1.0 - thermal_ratio))
        }
    }
}

// ============================================================================
// Grain boundaries
// ============================================================================

/// A grain boundary between two lattice sites with reduced exchange coupling.
///
/// The exchange coupling across the boundary is J_gb = η · J_bulk where
/// η ∈ [0, 1] is the coupling reduction factor.
#[derive(Debug, Clone)]
pub struct GrainBoundary {
    /// The pair of neighboring site indices connected by this boundary.
    pub sites: (usize, usize),
    /// Exchange coupling reduction factor η ∈ [0, 1].
    /// η = 1.0 means no reduction (perfect boundary).
    /// η = 0.0 means complete decoupling.
    pub coupling_reduction: f64,
}

impl GrainBoundary {
    /// Create a new grain boundary.
    ///
    /// # Arguments
    ///
    /// * `site_a` - First site index
    /// * `site_b` - Second site index
    /// * `eta` - Coupling reduction factor η ∈ [0, 1]
    ///
    /// # Errors
    ///
    /// Returns an error if η is outside [0, 1].
    pub fn new(site_a: usize, site_b: usize, eta: f64) -> Result<Self> {
        if !(0.0..=1.0).contains(&eta) {
            return Err(Error::InvalidParameter {
                param: "eta".to_string(),
                reason: "coupling reduction must be in [0, 1]".to_string(),
            });
        }
        Ok(Self {
            sites: (site_a, site_b),
            coupling_reduction: eta,
        })
    }

    /// Compute the effective exchange coupling across this boundary.
    ///
    /// J_gb = η · J_bulk
    ///
    /// # Arguments
    ///
    /// * `bulk_exchange` - Bulk exchange stiffness J_bulk \[J/m\]
    pub fn effective_exchange(&self, bulk_exchange: f64) -> f64 {
        self.coupling_reduction * bulk_exchange
    }

    /// Check if the boundary is completely decoupled.
    pub fn is_decoupled(&self) -> bool {
        self.coupling_reduction < f64::EPSILON
    }
}

// ============================================================================
// Defect collection and concentration effects
// ============================================================================

/// A collection of defects in a magnetic system, providing aggregate analysis.
#[derive(Debug, Clone)]
pub struct DefectCollection {
    /// All defect sites in the system.
    pub defects: Vec<DefectSite>,
    /// Total number of lattice sites in the system.
    pub total_sites: usize,
    /// Grain boundaries.
    pub grain_boundaries: Vec<GrainBoundary>,
}

impl DefectCollection {
    /// Create a new defect collection.
    ///
    /// # Arguments
    ///
    /// * `total_sites` - Total number of sites in the lattice
    ///
    /// # Errors
    ///
    /// Returns an error if `total_sites` is zero.
    pub fn new(total_sites: usize) -> Result<Self> {
        if total_sites == 0 {
            return Err(Error::InvalidParameter {
                param: "total_sites".to_string(),
                reason: "must be greater than zero".to_string(),
            });
        }
        Ok(Self {
            defects: Vec::new(),
            total_sites,
            grain_boundaries: Vec::new(),
        })
    }

    /// Add a defect to the collection.
    ///
    /// # Errors
    ///
    /// Returns an error if the defect position is out of range.
    pub fn add_defect(&mut self, defect: DefectSite) -> Result<()> {
        if defect.position >= self.total_sites {
            return Err(Error::InvalidParameter {
                param: "defect.position".to_string(),
                reason: format!(
                    "position {} exceeds total sites {}",
                    defect.position, self.total_sites
                ),
            });
        }
        self.defects.push(defect);
        Ok(())
    }

    /// Add a grain boundary to the collection.
    pub fn add_grain_boundary(&mut self, gb: GrainBoundary) -> Result<()> {
        if gb.sites.0 >= self.total_sites || gb.sites.1 >= self.total_sites {
            return Err(Error::InvalidParameter {
                param: "grain_boundary.sites".to_string(),
                reason: format!(
                    "site indices ({}, {}) exceed total sites {}",
                    gb.sites.0, gb.sites.1, self.total_sites
                ),
            });
        }
        self.grain_boundaries.push(gb);
        Ok(())
    }

    /// Defect concentration (number of defects / total sites).
    pub fn concentration(&self) -> f64 {
        self.defects.len() as f64 / self.total_sites as f64
    }

    /// Vacancy concentration (number of vacancies / total sites).
    pub fn vacancy_concentration(&self) -> f64 {
        let n_vac = self.defects.iter().filter(|d| d.is_vacancy()).count();
        n_vac as f64 / self.total_sites as f64
    }

    /// Number of pinning defects.
    pub fn num_pinning_sites(&self) -> usize {
        self.defects
            .iter()
            .filter(|d| matches!(d.defect_type, DefectType::Pinning { .. }))
            .count()
    }

    /// Estimate coercivity scaling with defect density.
    ///
    /// A simplified power-law model: H_c = H_c0 · n^α
    ///
    /// where n is the defect concentration and α is the scaling exponent
    /// (typically α ≈ 0.5 for dilute pinning centers).
    ///
    /// # Arguments
    ///
    /// * `h_c0` - Prefactor for coercivity \[A/m\]
    /// * `alpha` - Scaling exponent (dimensionless)
    ///
    /// # Returns
    ///
    /// Estimated coercivity \[A/m\].
    pub fn coercivity_estimate(&self, h_c0: f64, alpha: f64) -> f64 {
        let n = self.concentration();
        if n < f64::EPSILON {
            return 0.0;
        }
        h_c0 * n.powf(alpha)
    }

    /// Compute average exchange modification across all defect sites.
    ///
    /// For non-defect sites the exchange factor is 1.0; for defect sites it
    /// is given by `exchange_factor()`. This returns the site-averaged value.
    pub fn average_exchange_factor(&self) -> f64 {
        if self.total_sites == 0 {
            return 1.0;
        }
        let defect_sum: f64 = self.defects.iter().map(|d| d.exchange_factor()).sum();
        let non_defect_count = self.total_sites - self.defects.len();
        (defect_sum + non_defect_count as f64) / self.total_sites as f64
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pinning_potential_maximum_at_center() {
        let v0 = 1.0e-3; // J/m²
        let w = 5.0e-9; // 5 nm
        let x0 = 50.0e-9; // center at 50 nm

        let v_center = pinning_potential(x0, x0, v0, w);
        let v_away = pinning_potential(x0 + 10.0e-9, x0, v0, w);
        let v_far = pinning_potential(x0 + 50.0e-9, x0, v0, w);

        assert!(
            (v_center - v0).abs() < 1e-15,
            "Potential at center should equal V₀"
        );
        assert!(
            v_center > v_away,
            "Potential should decrease away from center"
        );
        assert!(
            v_away > v_far,
            "Potential should decrease further from center"
        );
    }

    #[test]
    fn test_depinning_field_is_positive() {
        let params =
            DepinningParams::new(1.0e-3, 8.0e5, 20.0e-9).expect("valid depinning parameters");

        let h_dep = params.depinning_field();
        assert!(
            h_dep > 0.0,
            "Depinning field should be positive, got {}",
            h_dep
        );

        // Check the formula: H_dep = V₀ / (2 μ₀ M_s δ_w)
        let expected = 1.0e-3 / (2.0 * MU_0 * 8.0e5 * 20.0e-9);
        assert!(
            (h_dep - expected).abs() / expected < 1e-10,
            "Depinning field mismatch: got {}, expected {}",
            h_dep,
            expected
        );
    }

    #[test]
    fn test_depinning_field_thermal_reduction() {
        let params =
            DepinningParams::new(1.0e-3, 8.0e5, 20.0e-9).expect("valid depinning parameters");

        let h_0 = params.depinning_field();
        let h_300k = params
            .depinning_field_thermal(300.0)
            .expect("valid temperature");

        assert!(
            h_300k < h_0,
            "Thermal depinning field {} should be less than zero-T value {}",
            h_300k,
            h_0
        );
        assert!(h_300k >= 0.0, "Depinning field should be non-negative");
    }

    #[test]
    fn test_vacancy_removes_local_contribution() {
        let defect = DefectSite::vacancy(42);

        assert!(defect.is_vacancy());
        assert!(
            defect.effective_anisotropy(5.0e4).abs() < f64::EPSILON,
            "Vacancy should have zero anisotropy"
        );
        assert!(
            defect.exchange_factor().abs() < f64::EPSILON,
            "Vacancy should have zero exchange"
        );
    }

    #[test]
    fn test_substitutional_defect() {
        let defect = DefectSite::substitutional(10, 2.0e4, 0.5);

        assert!(!defect.is_vacancy());
        assert!(
            (defect.effective_anisotropy(5.0e4) - 2.0e4).abs() < f64::EPSILON,
            "Substitutional should use its own anisotropy"
        );
        assert!(
            (defect.exchange_factor() - 0.5).abs() < f64::EPSILON,
            "Substitutional exchange mod should be 0.5"
        );
    }

    #[test]
    fn test_grain_boundary_reduces_exchange() {
        let gb = GrainBoundary::new(0, 1, 0.3).expect("valid grain boundary");
        let bulk_j = 1.0e-11; // J/m

        let j_gb = gb.effective_exchange(bulk_j);
        assert!(
            j_gb < bulk_j,
            "Grain boundary exchange {} should be less than bulk {}",
            j_gb,
            bulk_j
        );
        assert!(
            (j_gb - 0.3 * bulk_j).abs() / bulk_j < 1e-15,
            "J_gb should be η·J_bulk"
        );

        // Fully decoupled boundary
        let gb_decoupled = GrainBoundary::new(0, 1, 0.0).expect("valid grain boundary");
        assert!(gb_decoupled.is_decoupled());
        assert!(gb_decoupled.effective_exchange(bulk_j).abs() < f64::EPSILON);
    }

    #[test]
    fn test_grain_boundary_invalid_eta() {
        assert!(GrainBoundary::new(0, 1, 1.5).is_err());
        assert!(GrainBoundary::new(0, 1, -0.1).is_err());
    }

    #[test]
    fn test_defect_concentration_scaling() {
        // Build a collection and check coercivity scaling
        let n_sites = 10000;
        let mut collection = DefectCollection::new(n_sites).expect("valid collection");

        // Add 100 pinning defects (1% concentration)
        for i in 0..100 {
            collection
                .add_defect(DefectSite::pinning(i, 1.0e-3, 5.0e-9))
                .expect("valid defect");
        }

        let conc = collection.concentration();
        assert!(
            (conc - 0.01).abs() < 1e-10,
            "Concentration should be 0.01, got {}",
            conc
        );

        // H_c = H_c0 · n^0.5
        let h_c0 = 1.0e6; // prefactor [A/m]
        let alpha = 0.5;
        let h_c = collection.coercivity_estimate(h_c0, alpha);
        let expected = h_c0 * (0.01_f64).powf(0.5);
        assert!(
            (h_c - expected).abs() / expected < 1e-10,
            "Coercivity scaling: got {}, expected {}",
            h_c,
            expected
        );

        // Higher defect density should increase coercivity
        let mut collection2 = DefectCollection::new(n_sites).expect("valid collection");
        for i in 0..400 {
            collection2
                .add_defect(DefectSite::pinning(i, 1.0e-3, 5.0e-9))
                .expect("valid defect");
        }
        let h_c2 = collection2.coercivity_estimate(h_c0, alpha);
        assert!(
            h_c2 > h_c,
            "More defects should increase coercivity: {} vs {}",
            h_c2,
            h_c
        );
    }

    #[test]
    fn test_defect_collection_average_exchange() {
        let mut collection = DefectCollection::new(100).expect("valid collection");

        // Add 10 vacancies (exchange factor = 0)
        for i in 0..10 {
            collection
                .add_defect(DefectSite::vacancy(i))
                .expect("valid defect");
        }

        let avg = collection.average_exchange_factor();
        // 10 sites with factor 0, 90 sites with factor 1 → average = 90/100 = 0.9
        assert!(
            (avg - 0.9).abs() < 1e-10,
            "Average exchange factor should be 0.9, got {}",
            avg
        );
    }

    #[test]
    fn test_pinning_force_direction() {
        let v0 = 1.0e-3;
        let w = 5.0e-9;
        let x0 = 0.0;

        // Force should be positive for x > x0 (restoring toward center)
        // F = V₀ (x-x₀)/w² exp(...)
        let f_right = pinning_force(2.0e-9, x0, v0, w);
        let f_left = pinning_force(-2.0e-9, x0, v0, w);

        assert!(f_right > 0.0, "Force should push right for x > x₀");
        assert!(f_left < 0.0, "Force should push left for x < x₀");

        // Force at center should be zero (by symmetry)
        let f_center = pinning_force(x0, x0, v0, w);
        assert!(
            f_center.abs() < 1e-30,
            "Force at center should be zero, got {}",
            f_center
        );
    }

    #[test]
    fn test_depinning_invalid_params() {
        assert!(DepinningParams::new(0.0, 8.0e5, 20.0e-9).is_err());
        assert!(DepinningParams::new(1.0e-3, 0.0, 20.0e-9).is_err());
        assert!(DepinningParams::new(1.0e-3, 8.0e5, 0.0).is_err());
        assert!(DepinningParams::new(1.0e-3, 8.0e5, -1.0).is_err());
    }
}
