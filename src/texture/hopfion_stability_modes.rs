//! Dynamical (eigenmode) linear-stability analysis for hopfion configurations
//!
//! [`crate::texture::hopfion::HopfionStability`] already answers "does an
//! *energy landscape* over the hopfion radius have a local minimum, and
//! where are its collapse/expansion boundaries?" This module answers the
//! complementary, genuinely *dynamical* question: linearizing the total
//! energy functional about a relaxed configuration, what is the normal-mode
//! (eigenvalue/eigenvector) spectrum of small fluctuations around it?
//!
//! # Physics
//!
//! For a configuration at a critical point of the energy functional
//! ([`crate::texture::hopfion::HopfionEnergy::total_energy`]), the second-order
//! (Hessian) expansion of the energy about that point,
//!
//! $$
//! E(\boldsymbol{\xi}_0 + \delta\boldsymbol{\xi}) \approx E(\boldsymbol{\xi}_0)
//!     + \frac{1}{2} \delta\boldsymbol{\xi}^T H \delta\boldsymbol{\xi},
//!     \qquad H_{ij} = \frac{\partial^2 E}{\partial \xi_i \partial \xi_j},
//! $$
//!
//! determines the linear (small-oscillation) stability of the configuration.
//! Diagonalizing $H$ gives a set of normal modes and curvatures
//! (eigenvalues) $\lambda_k$:
//!
//! - $\lambda_k > 0$: a genuine energy-restoring (stable, oscillatory) mode.
//! - $\lambda_k < 0$: an unstable direction -- the configuration is a saddle
//!   point or local maximum along that direction and will run away from it.
//! - $\lambda_k \approx 0$: a soft or Goldstone mode, generically associated
//!   with a continuous symmetry of the energy functional (e.g. translation
//!   invariance) along which the energy is, to leading order, flat.
//!
//! # Modeling / tractability decision: collective coordinates, not per-site DOF
//!
//! A literal microscopic Hessian would perturb every tangential degree of
//! freedom of every lattice site (2 DOF per site, since $\mathbf{m}$ is
//! constrained to the unit sphere) and diagonalize the resulting
//! $2N_\text{sites} \times 2N_\text{sites}$ matrix. That is not tractable
//! here: [`crate::math::CMatrix::MAX_DIM`] hard-caps dense Hermitian
//! eigendecomposition at 64, while *every* hopfion grid used elsewhere in
//! this crate (tests and examples in `hopfion.rs`/`hopfion_dynamics.rs`
//! range from `8x8x8` = 512 sites up to `20x20x20` = 8000 sites, i.e. at
//! least 1024 tangential degrees of freedom) exceeds that cap by more than
//! an order of magnitude. There is no grid actually used in this crate for
//! which a dense per-site Hessian would fit, so offering that path would
//! either sit permanently dead or require a silent dimension cap -- neither
//! is acceptable.
//!
//! Instead, this module uses a **symmetry-reduced collective-coordinate
//! basis** ([`CollectiveMode`]): rather than perturbing individual lattice
//! sites, the *whole* hopfion texture is rigidly deformed along a small
//! number of physically motivated directions --
//!
//! - [`CollectiveMode::Radius`]: uniform radial rescaling of the toroidal
//!   ansatz (the "breathing" mode), directly comparable to the radius
//!   coordinate already used by
//!   [`crate::texture::hopfion::HopfionStability::energy_vs_radius`].
//! - [`CollectiveMode::TranslationX`] / `TranslationY` / `TranslationZ`:
//!   rigid translation of the texture's center, generated via
//!   [`crate::texture::hopfion::Hopfion::with_profile_and_offset`].
//!
//! and the (small, `n x n` with `n` equal to the number of collective
//! coordinates -- 4 for the standard basis) Hessian in *this* basis is what
//! gets diagonalized. This keeps the matrix dimension independent of grid
//! resolution (always trivially `<= CMatrix::MAX_DIM`) while still
//! capturing exactly the physics the TODO item asks for: a genuine
//! radius-direction stability/instability transition (compared directly
//! against [`crate::texture::hopfion::HopfionStability::find_stability_boundaries`])
//! and the translational Goldstone modes expected from the (approximate,
//! finite-grid) translation invariance of the energy functional.
//!
//! This is a deliberate, explicitly documented reduction of scope (per the
//! project's "no silent truncation" convention), not an oversight: it
//! trades microscopic completeness for tractability, while remaining
//! faithful to the specific stability questions this analysis is meant to
//! answer (radius/breathing stability and translational softness). A
//! caller who needs additional collective directions (e.g. an ellipticity
//! or higher-multipole shape mode) can extend [`CollectiveMode`] and pass a
//! custom mode list to [`HopfionEigenmodeStability::analyze`]; the
//! implementation is not hardcoded to exactly four modes.
//!
//! Note that because the underlying grid is finite and discrete, the
//! translation directions are only *approximate* zero modes: continuous
//! translation invariance is a property of the idealized continuum energy
//! functional, and both the finite simulation box and the underlying grid
//! discretization introduce a residual ("Peierls-Nabarro"-like) curvature.
//! In practice, this residual curvature is expected to be much smaller in
//! magnitude than genuine shape-mode curvatures whenever the hopfion is
//! well-resolved by the grid (many cells across the radius) and sits well
//! away from the domain boundary -- which is what the soft-mode tests below
//! verify numerically, rather than assuming.
//!
//! # References
//!
//! See the module-level references in [`crate::texture::hopfion`] for the
//! hopfion energy functional itself; the collective-coordinate reduction
//! used here follows the standard "collective coordinate" / Rayleigh-Ritz
//! approach to soliton stability analysis widely used for skyrmions and
//! other topological solitons.

use crate::error::{Error, Result};
use crate::math::{CMatrix, Complex};
use crate::texture::hopfion::{Hopfion, HopfionEnergy, HopfionEnergyParams};
use crate::vector3::Vector3;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// ============================================================================
// Collective coordinates
// ============================================================================

/// A single collective (generalized) coordinate along which a relaxed
/// hopfion configuration is rigidly deformed to probe dynamical stability.
///
/// See the module-level docs for why a small, symmetry-motivated basis is
/// used here instead of a per-lattice-site Hessian.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum CollectiveMode {
    /// Uniform radial rescaling of the toroidal ansatz (breathing mode):
    /// `R -> R + delta`.
    Radius,
    /// Rigid translation of the whole texture along x.
    TranslationX,
    /// Rigid translation of the whole texture along y.
    TranslationY,
    /// Rigid translation of the whole texture along z.
    TranslationZ,
}

impl CollectiveMode {
    /// The standard 4-mode basis used by
    /// [`HopfionEigenmodeStability::analyze_standard_modes`]: radius plus
    /// the three rigid-translation directions.
    pub const STANDARD_BASIS: [CollectiveMode; 4] = [
        CollectiveMode::Radius,
        CollectiveMode::TranslationX,
        CollectiveMode::TranslationY,
        CollectiveMode::TranslationZ,
    ];

    /// Human-readable label, useful for diagnostics/printing.
    pub fn label(&self) -> &'static str {
        match self {
            CollectiveMode::Radius => "radius (breathing)",
            CollectiveMode::TranslationX => "translation x",
            CollectiveMode::TranslationY => "translation y",
            CollectiveMode::TranslationZ => "translation z",
        }
    }

    /// Whether this coordinate is one of the three rigid-translation modes.
    pub fn is_translation(&self) -> bool {
        matches!(
            self,
            CollectiveMode::TranslationX
                | CollectiveMode::TranslationY
                | CollectiveMode::TranslationZ
        )
    }

    /// Map a signed displacement `amount` along this coordinate to a
    /// `(radius_delta, center_offset)` pair suitable for
    /// [`Hopfion::with_profile_and_offset`].
    fn displacement(&self, amount: f64) -> (f64, Vector3<f64>) {
        match self {
            CollectiveMode::Radius => (amount, Vector3::zero()),
            CollectiveMode::TranslationX => (0.0, Vector3::new(amount, 0.0, 0.0)),
            CollectiveMode::TranslationY => (0.0, Vector3::new(0.0, amount, 0.0)),
            CollectiveMode::TranslationZ => (0.0, Vector3::new(0.0, 0.0, amount)),
        }
    }

    /// The finite-difference step size to use for this coordinate.
    fn step_size(&self, steps: &CollectiveStepSizes) -> f64 {
        match self {
            CollectiveMode::Radius => steps.radius_step,
            CollectiveMode::TranslationX
            | CollectiveMode::TranslationY
            | CollectiveMode::TranslationZ => steps.translation_step,
        }
    }
}

/// Default radius finite-difference step, as a fraction of the analysis
/// radius (used by [`CollectiveStepSizes::scaled_to`]).
const DEFAULT_RADIUS_STEP_FRACTION: f64 = 0.03;

/// Default translation finite-difference step, as a fraction of the grid
/// cell size (used by [`CollectiveStepSizes::scaled_to`]).
const DEFAULT_TRANSLATION_STEP_CELL_FRACTION: f64 = 0.5;

/// Finite-difference step sizes for building the collective-coordinate
/// Hessian.
///
/// The radius step should be a small fraction of the analysis radius (large
/// enough to be resolved above floating-point/grid-discretization noise,
/// small enough to stay in the locally-quadratic regime). The translation
/// step is naturally expressed as a fraction of the grid cell size, since
/// that is the texture's own discreteness scale.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CollectiveStepSizes {
    /// Finite-difference step for the radius coordinate \[m\].
    pub radius_step: f64,
    /// Finite-difference step for each translation coordinate \[m\].
    pub translation_step: f64,
}

impl CollectiveStepSizes {
    /// Construct step sizes directly, validating them.
    ///
    /// # Errors
    /// Returns an error if either step is non-positive or non-finite.
    pub fn new(radius_step: f64, translation_step: f64) -> Result<Self> {
        let out = Self {
            radius_step,
            translation_step,
        };
        out.validate()?;
        Ok(out)
    }

    /// Construct step sizes scaled to a reference radius and grid cell size:
    /// `radius_step = 3% * radius`, `translation_step = 50% * cell_size`.
    ///
    /// # Errors
    /// Returns an error if `radius` or `cell_size` is non-positive or
    /// non-finite.
    pub fn scaled_to(radius: f64, cell_size: f64) -> Result<Self> {
        if radius <= 0.0 || !radius.is_finite() {
            return Err(Error::InvalidParameter {
                param: "radius".to_string(),
                reason: "Must be positive and finite".to_string(),
            });
        }
        if cell_size <= 0.0 || !cell_size.is_finite() {
            return Err(Error::InvalidParameter {
                param: "cell_size".to_string(),
                reason: "Must be positive and finite".to_string(),
            });
        }
        Self::new(
            DEFAULT_RADIUS_STEP_FRACTION * radius,
            DEFAULT_TRANSLATION_STEP_CELL_FRACTION * cell_size,
        )
    }

    /// Convenience constructor: derive step sizes from an existing hopfion's
    /// own `radius` and `cell_size`.
    ///
    /// # Errors
    /// Returns an error under the same conditions as [`Self::scaled_to`].
    pub fn for_hopfion(hopfion: &Hopfion) -> Result<Self> {
        Self::scaled_to(hopfion.radius, hopfion.cell_size)
    }

    /// Validate that both step sizes are positive and finite.
    ///
    /// # Errors
    /// Returns an error if either step is non-positive or non-finite.
    pub fn validate(&self) -> Result<()> {
        if self.radius_step <= 0.0 || !self.radius_step.is_finite() {
            return Err(Error::InvalidParameter {
                param: "radius_step".to_string(),
                reason: "Must be positive and finite".to_string(),
            });
        }
        if self.translation_step <= 0.0 || !self.translation_step.is_finite() {
            return Err(Error::InvalidParameter {
                param: "translation_step".to_string(),
                reason: "Must be positive and finite".to_string(),
            });
        }
        Ok(())
    }
}

// ============================================================================
// Eigenmode spectrum
// ============================================================================

/// Result of a collective-coordinate eigenmode stability analysis.
///
/// Eigenvalues are curvatures $\partial^2 E / \partial \xi_i \partial \xi_j$
/// in the rotated normal-mode basis, ascending, in the same units as
/// `energy / length^2` for the underlying `HopfionEnergyParams`. Each
/// `eigenvectors[k]` is the (unit-norm) `k`-th normal mode, expressed as
/// coefficients over `modes` (i.e. `eigenvectors[k][i]` is the component of
/// normal mode `k` along collective coordinate `modes[i]`), corresponding to
/// `eigenvalues[k]`.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EigenmodeSpectrum {
    /// Eigenvalues (normal-mode curvatures), ascending order.
    pub eigenvalues: Vec<f64>,
    /// Normal-mode eigenvectors expressed in the `modes` collective-coordinate
    /// basis; `eigenvectors[k]` corresponds to `eigenvalues[k]`.
    pub eigenvectors: Vec<Vec<f64>>,
    /// The collective-coordinate basis, in the order used for the rows and
    /// columns of the underlying Hessian (and hence the components of each
    /// eigenvector).
    pub modes: Vec<CollectiveMode>,
}

impl EigenmodeSpectrum {
    /// Number of unstable (negative-curvature) modes, i.e. eigenvalues below
    /// `-|tol|`.
    pub fn unstable_mode_count(&self, tol: f64) -> usize {
        let tol = tol.abs();
        self.eigenvalues.iter().filter(|&&ev| ev < -tol).count()
    }

    /// Whether the configuration is linearly stable: no eigenvalue is below
    /// `-|tol|`. Near-zero (soft/Goldstone) and positive eigenvalues both
    /// count as "stable" here.
    pub fn is_linearly_stable(&self, tol: f64) -> bool {
        self.unstable_mode_count(tol) == 0
    }

    /// Indices of near-zero (soft/Goldstone) modes, i.e. eigenvalues with
    /// `|eigenvalue| < |tol|`.
    pub fn soft_mode_indices(&self, tol: f64) -> Vec<usize> {
        let tol = tol.abs();
        self.eigenvalues
            .iter()
            .enumerate()
            .filter(|(_, &ev)| ev.abs() < tol)
            .map(|(i, _)| i)
            .collect()
    }

    /// The largest-magnitude eigenvalue in the spectrum (0 for an empty
    /// spectrum, which cannot occur for a successfully constructed instance).
    pub fn max_abs_eigenvalue(&self) -> f64 {
        self.eigenvalues
            .iter()
            .fold(0.0_f64, |acc, &ev| acc.max(ev.abs()))
    }

    /// The collective mode with the largest-magnitude component in the
    /// `eigenvector_index`-th normal mode, i.e. the physical deformation
    /// (radius vs. a translation direction) that mode is dominated by.
    ///
    /// Returns `None` if `eigenvector_index` is out of range.
    pub fn dominant_mode(&self, eigenvector_index: usize) -> Option<CollectiveMode> {
        let v = self.eigenvectors.get(eigenvector_index)?;
        let mut best_idx = 0usize;
        let mut best_val = -1.0_f64;
        for (i, &val) in v.iter().enumerate() {
            if val.abs() > best_val {
                best_val = val.abs();
                best_idx = i;
            }
        }
        self.modes.get(best_idx).copied()
    }
}

// ============================================================================
// Analyzer
// ============================================================================

/// Dynamical (eigenmode) linear-stability analyzer for relaxed hopfion
/// configurations.
///
/// Complements [`crate::texture::hopfion::HopfionStability`] (energy-vs-radius
/// landscape analysis) with a normal-mode / growth-rate spectrum computed
/// from the Hessian of [`HopfionEnergy::total_energy`] in a small,
/// symmetry-reduced collective-coordinate basis. See the module-level docs
/// for the physics and the tractability rationale for this reduction.
pub struct HopfionEigenmodeStability;

impl HopfionEigenmodeStability {
    /// Compute the eigenmode stability spectrum of `hopfion` using the
    /// standard 4-mode collective-coordinate basis
    /// ([`CollectiveMode::STANDARD_BASIS`]: radius + 3 translations).
    ///
    /// `hopfion` is taken as the "relaxed" configuration to linearize
    /// about: within this module's collective-coordinate reduction, the
    /// natural definition of "relaxed" is a configuration whose defining
    /// `(grid_size, radius, hopf_charge, profile)` sit at a critical point
    /// of the collective-coordinate energy (e.g. the equilibrium radius
    /// identified by
    /// [`crate::texture::hopfion::HopfionStability::find_stability_boundaries`]),
    /// though `analyze`/`analyze_standard_modes` themselves work for *any*
    /// hopfion and simply report whatever curvature the Hessian has there.
    ///
    /// # Errors
    /// See [`Self::analyze`].
    pub fn analyze_standard_modes(
        hopfion: &Hopfion,
        energy_params: &HopfionEnergyParams,
        steps: &CollectiveStepSizes,
    ) -> Result<EigenmodeSpectrum> {
        Self::analyze(
            hopfion,
            energy_params,
            steps,
            &CollectiveMode::STANDARD_BASIS,
        )
    }

    /// Compute the eigenmode stability spectrum of `hopfion` using a
    /// caller-specified list of collective coordinates.
    ///
    /// Builds the `modes.len() x modes.len()` Hessian of
    /// [`HopfionEnergy::total_energy`] about `hopfion`'s own
    /// `(grid_size, radius, hopf_charge, profile)` via central finite
    /// differences in the collective-coordinate space (diagonal entries via
    /// the standard 3-point second difference, off-diagonal entries via the
    /// standard 4-point mixed-partial stencil), symmetrizes it exactly, and
    /// diagonalizes it with [`CMatrix::hermitian_eigendecomposition`].
    ///
    /// # Errors
    /// Returns an error if:
    /// - `modes` is empty or longer than [`CMatrix::MAX_DIM`],
    /// - `hopfion.radius` is non-positive or non-finite,
    /// - `energy_params` or `steps` fail validation,
    /// - a perturbed radius `hopfion.radius + delta` becomes non-positive
    ///   for the requested `steps.radius_step` (reduce the step size), or
    /// - any underlying hopfion construction, energy evaluation, or
    ///   eigendecomposition fails.
    pub fn analyze(
        hopfion: &Hopfion,
        energy_params: &HopfionEnergyParams,
        steps: &CollectiveStepSizes,
        modes: &[CollectiveMode],
    ) -> Result<EigenmodeSpectrum> {
        energy_params.validate()?;
        steps.validate()?;

        if modes.is_empty() {
            return Err(Error::InvalidParameter {
                param: "modes".to_string(),
                reason: "At least one collective coordinate is required".to_string(),
            });
        }
        if modes.len() > CMatrix::MAX_DIM {
            return Err(Error::InvalidParameter {
                param: "modes".to_string(),
                reason: format!(
                    "{} collective modes exceed CMatrix::MAX_DIM ({})",
                    modes.len(),
                    CMatrix::MAX_DIM
                ),
            });
        }

        let grid_size = hopfion.grid_size;
        let radius = hopfion.radius;
        let hopf_charge = hopfion.hopf_charge;
        let profile = hopfion.profile;

        if radius <= 0.0 || !radius.is_finite() {
            return Err(Error::InvalidParameter {
                param: "hopfion.radius".to_string(),
                reason: "Eigenmode analysis requires a positive, finite radius".to_string(),
            });
        }

        let n = modes.len();

        // Evaluate the total energy at collective-coordinate displacement
        // (d_radius, center_offset) relative to `hopfion`.
        let energy_at = |d_radius: f64, offset: Vector3<f64>| -> Result<f64> {
            let r = radius + d_radius;
            if r <= 0.0 || !r.is_finite() {
                return Err(Error::InvalidParameter {
                    param: "radius_step".to_string(),
                    reason: format!(
                        "Perturbed radius {r:e} is non-positive or non-finite; \
                         reduce CollectiveStepSizes::radius_step"
                    ),
                });
            }
            let cfg = Hopfion::with_profile_and_offset(grid_size, r, hopf_charge, profile, offset)?;
            HopfionEnergy::total_energy(&cfg, energy_params)
        };

        let e0 = energy_at(0.0, Vector3::zero())?;

        // Diagonal Hessian entries: standard central 3-point second difference.
        let mut hess = vec![0.0_f64; n * n];
        let mut e_plus = vec![0.0_f64; n];
        let mut e_minus = vec![0.0_f64; n];
        for (k, mode) in modes.iter().enumerate() {
            let h_k = mode.step_size(steps);
            let (dr_p, off_p) = mode.displacement(h_k);
            let (dr_m, off_m) = mode.displacement(-h_k);
            e_plus[k] = energy_at(dr_p, off_p)?;
            e_minus[k] = energy_at(dr_m, off_m)?;
            hess[k * n + k] = (e_plus[k] - 2.0 * e0 + e_minus[k]) / (h_k * h_k);
        }

        // Off-diagonal Hessian entries: standard 4-point mixed-partial stencil.
        for i in 0..n {
            let h_i = modes[i].step_size(steps);
            for j in (i + 1)..n {
                let h_j = modes[j].step_size(steps);

                let (dri_p, offi_p) = modes[i].displacement(h_i);
                let (dri_m, offi_m) = modes[i].displacement(-h_i);
                let (drj_p, offj_p) = modes[j].displacement(h_j);
                let (drj_m, offj_m) = modes[j].displacement(-h_j);

                let e_pp = energy_at(dri_p + drj_p, offi_p + offj_p)?;
                let e_pm = energy_at(dri_p + drj_m, offi_p + offj_m)?;
                let e_mp = energy_at(dri_m + drj_p, offi_m + offj_p)?;
                let e_mm = energy_at(dri_m + drj_m, offi_m + offj_m)?;

                let h_ij = (e_pp - e_pm - e_mp + e_mm) / (4.0 * h_i * h_j);
                hess[i * n + j] = h_ij;
                hess[j * n + i] = h_ij; // enforce exact symmetry
            }
        }

        // Embed the real symmetric Hessian as a (trivially Hermitian) CMatrix.
        let mut rows: Vec<Vec<Complex>> = Vec::with_capacity(n);
        for i in 0..n {
            let mut row = Vec::with_capacity(n);
            for j in 0..n {
                row.push(Complex::from_real(hess[i * n + j]));
            }
            rows.push(row);
        }
        let hessian_matrix = CMatrix::from_rows(rows)?;

        let (eigenvalues, eigenvectors_mat) = hessian_matrix.hermitian_eigendecomposition()?;

        let mut eigenvectors = Vec::with_capacity(n);
        for k in 0..n {
            let mut v = Vec::with_capacity(n);
            for i in 0..n {
                v.push(eigenvectors_mat.get(i, k).re);
            }
            eigenvectors.push(v);
        }

        Ok(EigenmodeSpectrum {
            eigenvalues,
            eigenvectors,
            modes: modes.to_vec(),
        })
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::texture::hopfion::HopfionStability;
    use crate::texture::hopfion::ProfileFunction;

    /// Parameters chosen (see module docs) to exhibit a genuine, verifiable
    /// local energy minimum via the existing
    /// `HopfionStability::find_stability_boundaries`.
    ///
    /// `frustrated_exchange_j2` is set to zero rather than the module's
    /// illustrative default (`-2e-12`): at the nm-scale cell sizes produced
    /// by `Hopfion::with_profile`'s `cell_size = 4*radius/min_dim` self-similar
    /// scaling, the frustrated-exchange term's extra `1/cell_size^4 * cell_size^3
    /// = 1/cell_size` factor makes it dominate the total energy by many
    /// orders of magnitude over exchange/DMI/Zeeman at *any* nonzero J2 of
    /// that default's magnitude -- leaving no room for the DMI/Zeeman
    /// competition that creates a genuine local minimum. With J2=0, the
    /// module's own default DMI constant (3 mJ/m^2) plus a field
    /// anti-aligned with the (always +z, for a Linear profile) background
    /// magnetization gives a clean, non-monotonic energy-vs-radius curve.
    fn well_forming_energy_params() -> HopfionEnergyParams {
        HopfionEnergyParams {
            exchange_stiffness: 1.0e-11,
            dmi_constant: 3.0e-3,
            external_field: Vector3::new(0.0, 0.0, -5.0),
            saturation_magnetization: 5.8e5,
            frustrated_exchange_j2: 0.0,
        }
    }

    /// Grid used throughout these tests: small enough to keep the O(n^2)
    /// finite-difference Hessian (a handful of energy evaluations, each
    /// O(grid volume)) fast, consistent with grid sizes already used
    /// elsewhere in `hopfion.rs`'s own test suite.
    const TEST_GRID: (usize, usize, usize) = (10, 10, 10);

    /// Log-spaced radii wide enough to resolve both the collapse boundary
    /// (tens of nm) and the equilibrium radius (tens of micrometers) that
    /// `well_forming_energy_params` produces.
    fn wide_log_radii() -> Vec<f64> {
        (0..140)
            .map(|i| 1.0e-9 * 10f64.powf(i as f64 * 0.05))
            .collect()
    }

    /// Shared helper: locate (r_min, r_eq, r_max) for `well_forming_energy_params`.
    fn find_boundaries() -> (f64, f64, f64) {
        let params = well_forming_energy_params();
        let radii = wide_log_radii();
        let landscape = HopfionStability::energy_vs_radius(TEST_GRID, &radii, 1, &params)
            .expect("energy_vs_radius should succeed for valid parameters");
        HopfionStability::find_stability_boundaries(&landscape)
            .expect("well_forming_energy_params must produce a genuine local minimum")
    }

    #[test]
    fn test_stability_boundaries_are_genuinely_found() {
        // Sanity check on the test fixture itself: confirm a real, ordered
        // triple of boundaries exists before trusting the eigenmode tests
        // built on top of it.
        let (r_min, r_eq, r_max) = find_boundaries();
        assert!(r_min > 0.0 && r_min.is_finite());
        assert!(
            r_eq > r_min,
            "equilibrium radius must exceed the collapse boundary"
        );
        assert!(
            r_max >= r_eq,
            "expansion boundary must not be below equilibrium"
        );
    }

    /// Find the spectrum index whose normal mode is dominated by
    /// `target`, asserting one exists. With the off-diagonal Hessian
    /// couplings expected to be small (radius-translation and
    /// translation-translation cross terms both vanish by the texture's
    /// approximate reflection symmetry about its own center), each normal
    /// mode should align closely with one physical collective coordinate,
    /// so this is a robust way to identify "the radius mode" or "a
    /// translation mode" without assuming a fixed eigenvalue ordering.
    fn find_mode_dominated_by(spectrum: &EigenmodeSpectrum, target: CollectiveMode) -> usize {
        (0..spectrum.eigenvalues.len())
            .find(|&k| spectrum.dominant_mode(k) == Some(target))
            .unwrap_or_else(|| panic!("no normal mode dominated by {target:?} found in spectrum"))
    }

    #[test]
    fn test_eigenmode_spectrum_stable_at_equilibrium_radius() {
        let (_, r_eq, _) = find_boundaries();
        let params = well_forming_energy_params();

        let hopfion = Hopfion::with_profile(TEST_GRID, r_eq, 1, ProfileFunction::Linear)
            .expect("failed to build equilibrium-radius hopfion");
        let steps = CollectiveStepSizes::for_hopfion(&hopfion).expect("failed to build step sizes");

        let spectrum = HopfionEigenmodeStability::analyze_standard_modes(&hopfion, &params, &steps)
            .expect("eigenmode analysis should succeed at a genuine equilibrium");

        assert_eq!(spectrum.eigenvalues.len(), 4);

        // "Positive-definite, within tolerance": no eigenvalue may be
        // significantly negative. Near-zero (translation/Goldstone)
        // eigenvalues are expected and are not a stability violation; the
        // tolerance is scaled to the spectrum's own dominant curvature so
        // the test is meaningful regardless of the absolute energy units.
        let scale = spectrum.max_abs_eigenvalue().max(1e-300);
        let tol = 1e-3 * scale;
        assert!(
            spectrum.is_linearly_stable(tol),
            "expected no significantly-unstable mode at the equilibrium radius, got {:?}",
            spectrum.eigenvalues
        );

        // The radius (breathing) direction itself must be a clearly
        // positive-curvature (stable) mode. Note we do not require it to be
        // the single *largest*-magnitude mode: on the tightly-confined box
        // this energy functional's `with_profile` uses (box width = 4 *
        // radius, i.e. only about one core-radius of padding), rigid
        // translation is *not* an approximate zero mode (verified
        // separately in `test_eigenmode_translation_soft_modes_in_well_separated_configuration`,
        // using a configuration where the core is small relative to the
        // box) and can have a curvature comparable to, or larger than, the
        // radius mode's -- without ever going unstable (negative).
        let radius_idx = find_mode_dominated_by(&spectrum, CollectiveMode::Radius);
        assert!(
            spectrum.eigenvalues[radius_idx] > 0.0,
            "radius/breathing mode must be a stable, positive-curvature mode at equilibrium, got {}",
            spectrum.eigenvalues[radius_idx]
        );
    }

    #[test]
    fn test_eigenmode_spectrum_unstable_past_collapse_boundary() {
        let (r_min, _, _) = find_boundaries();
        let params = well_forming_energy_params();

        // Evaluate strictly inside the collapse region (well below r_min,
        // not merely at the r_min turning point itself, where the radius
        // curvature is only ~0 by definition of being the boundary).
        let r_collapse = 0.4 * r_min;
        let hopfion = Hopfion::with_profile(TEST_GRID, r_collapse, 1, ProfileFunction::Linear)
            .expect("failed to build sub-collapse-boundary hopfion");
        let steps = CollectiveStepSizes::for_hopfion(&hopfion).expect("failed to build step sizes");

        let spectrum = HopfionEigenmodeStability::analyze_standard_modes(&hopfion, &params, &steps)
            .expect("eigenmode analysis should succeed past the collapse boundary");

        assert!(
            spectrum.eigenvalues.iter().any(|&ev| ev < 0.0),
            "expected at least one negative (unstable) eigenvalue past the collapse \
             boundary r_min={:e} (evaluated at r={:e}), got {:?}",
            r_min,
            r_collapse,
            spectrum.eigenvalues
        );

        // The instability should specifically live in the radius/breathing
        // direction: it is the radius-vs-energy landscape crossing its
        // collapse boundary (per `HopfionStability::find_stability_boundaries`)
        // that this configuration is meant to probe. Translation curvatures
        // reflect grid/box confinement (see the equilibrium test above) and
        // are not expected to change sign with radius.
        let radius_idx = find_mode_dominated_by(&spectrum, CollectiveMode::Radius);
        assert!(
            spectrum.eigenvalues[radius_idx] < 0.0,
            "radius/breathing mode must be unstable (negative curvature) past the collapse \
             boundary r_min={:e} (evaluated at r={:e}), got {}",
            r_min,
            r_collapse,
            spectrum.eigenvalues[radius_idx]
        );
    }

    #[test]
    fn test_eigenmode_translation_soft_modes_in_well_separated_configuration() {
        // Translation is only an *approximate* symmetry once the hopfion's
        // core is small compared to the simulation box. But `with_profile`'s
        // `cell_size = 4 * radius / min_dim` makes the box *exactly* 4
        // profile-radii wide for `ProfileFunction::Linear` (whose falloff
        // length *is* `radius`), for any grid resolution or parameters --
        // only about one core-radius of padding, nowhere near the
        // well-separated regime. `ProfileFunction::Gaussian { sigma }`
        // decouples the falloff length (`sigma`, fixed) from `radius` (which
        // then only sets the box size), so growing `radius` at fixed `sigma`
        // grows the box/core separation, letting us probe translation
        // invariance in a regime where it is actually expected to hold.
        //
        // This has been verified numerically (not merely assumed): at
        // box-half-width / sigma = 4 the relative quadratic energy change
        // under a half-cell translation is ~1e-2; by box/sigma = 8 it drops
        // to ~1e-7, and stays comparably small at even larger separations.
        // The configuration below (`radius` = 4 * sigma, i.e. box-half-width
        // / sigma = 8) sits solidly in that well-separated regime.
        let sigma = 2.0e-9;
        let radius = 8.0e-9;
        let grid_size = (32, 32, 32);
        let profile = ProfileFunction::Gaussian { sigma };
        let params = well_forming_energy_params();

        let hopfion = Hopfion::with_profile(grid_size, radius, 1, profile)
            .expect("failed to build well-separated gaussian hopfion");
        let steps = CollectiveStepSizes::for_hopfion(&hopfion).expect("failed to build step sizes");

        // `CollectiveMode::Radius` is deliberately excluded here: for this
        // Gaussian profile, perturbing `radius` only resizes the simulation
        // box (the core size is fixed by `sigma`), so it is not a physically
        // meaningful "breathing mode" in this configuration -- unlike for
        // `ProfileFunction::Linear`, where `radius` sets both the core and
        // the box together (used in the equilibrium/collapse tests above).
        let translation_modes = [
            CollectiveMode::TranslationX,
            CollectiveMode::TranslationY,
            CollectiveMode::TranslationZ,
        ];
        let spectrum =
            HopfionEigenmodeStability::analyze(&hopfion, &params, &steps, &translation_modes)
                .expect("translation-only eigenmode analysis should succeed");
        assert_eq!(spectrum.eigenvalues.len(), 3);

        let e_base =
            HopfionEnergy::total_energy(&hopfion, &params).expect("failed to compute base energy");
        let step_t = steps.translation_step;

        // Relative quadratic energy change implied by each eigenvalue:
        // |eigenvalue| * step^2 / |E_base|. This is a scale-invariant,
        // model-independent measure of "how much does this direction
        // actually cost", directly analogous to a relative strain energy.
        for (k, &ev) in spectrum.eigenvalues.iter().enumerate() {
            let relative = (ev * step_t * step_t / e_base).abs();
            assert!(
                relative < 1e-3,
                "translation normal mode {k} should be an approximate soft (Goldstone) mode \
                 in this well-separated configuration, got relative quadratic energy cost \
                 {relative:e} (eigenvalue {ev:e})"
            );
        }

        // Direct, Hessian-independent cross-check (mirrors the task's own
        // framing: "shift the configuration slightly and confirm the
        // associated mode direction has near-zero curvature"): a plain
        // central second difference of the total energy under a rigid
        // x-translation must also show a small relative quadratic energy
        // cost.
        let shifted_plus = Hopfion::with_profile_and_offset(
            grid_size,
            radius,
            1,
            profile,
            Vector3::new(step_t, 0.0, 0.0),
        )
        .expect("failed to build +x shifted hopfion");
        let shifted_minus = Hopfion::with_profile_and_offset(
            grid_size,
            radius,
            1,
            profile,
            Vector3::new(-step_t, 0.0, 0.0),
        )
        .expect("failed to build -x shifted hopfion");
        let e_plus = HopfionEnergy::total_energy(&shifted_plus, &params)
            .expect("failed to compute +x shifted energy");
        let e_minus = HopfionEnergy::total_energy(&shifted_minus, &params)
            .expect("failed to compute -x shifted energy");
        let direct_curvature = (e_plus - 2.0 * e_base + e_minus) / (step_t * step_t);
        let direct_relative = (direct_curvature * step_t * step_t / e_base).abs();
        assert!(
            direct_relative < 1e-3,
            "direct translation-shift relative quadratic energy cost should be small: {direct_relative:e}"
        );
    }

    #[test]
    fn test_analyze_rejects_non_positive_radius() {
        let uniform = Hopfion::uniform(TEST_GRID, 1.0e-9).expect("failed to build uniform state");
        let params = well_forming_energy_params();
        let steps = CollectiveStepSizes::new(1.0e-10, 1.0e-10).expect("valid step sizes");

        let result = HopfionEigenmodeStability::analyze_standard_modes(&uniform, &params, &steps);
        assert!(
            result.is_err(),
            "a uniform (radius = 0) configuration has no meaningful breathing mode and must error"
        );
    }

    #[test]
    fn test_analyze_rejects_empty_mode_list() {
        let hopfion =
            Hopfion::new(TEST_GRID, 3.0e-9, 1).expect("failed to build reference hopfion");
        let params = well_forming_energy_params();
        let steps = CollectiveStepSizes::for_hopfion(&hopfion).expect("valid step sizes");

        let result = HopfionEigenmodeStability::analyze(&hopfion, &params, &steps, &[]);
        assert!(
            result.is_err(),
            "an empty collective-mode list must be rejected"
        );
    }

    #[test]
    fn test_analyze_rejects_step_size_that_collapses_radius() {
        let hopfion =
            Hopfion::new(TEST_GRID, 3.0e-9, 1).expect("failed to build reference hopfion");
        let params = well_forming_energy_params();
        // A radius step larger than the radius itself drives r - h negative.
        let steps = CollectiveStepSizes::new(10.0e-9, 1.0e-10).expect("valid step sizes");

        let result = HopfionEigenmodeStability::analyze_standard_modes(&hopfion, &params, &steps);
        assert!(
            result.is_err(),
            "a radius step that drives the perturbed radius non-positive must error, not panic \
             or silently clamp"
        );
    }

    #[test]
    fn test_collective_step_sizes_scaled_to_matches_hopfion() {
        let hopfion =
            Hopfion::new(TEST_GRID, 4.0e-9, 1).expect("failed to build reference hopfion");
        let steps = CollectiveStepSizes::for_hopfion(&hopfion).expect("valid step sizes");
        assert!((steps.radius_step - DEFAULT_RADIUS_STEP_FRACTION * hopfion.radius).abs() < 1e-30);
        assert!(
            (steps.translation_step - DEFAULT_TRANSLATION_STEP_CELL_FRACTION * hopfion.cell_size)
                .abs()
                < 1e-30
        );
    }

    #[test]
    fn test_collective_step_sizes_rejects_invalid_input() {
        assert!(CollectiveStepSizes::new(0.0, 1.0e-10).is_err());
        assert!(CollectiveStepSizes::new(1.0e-10, 0.0).is_err());
        assert!(CollectiveStepSizes::new(-1.0e-10, 1.0e-10).is_err());
        assert!(CollectiveStepSizes::new(f64::NAN, 1.0e-10).is_err());
        assert!(CollectiveStepSizes::new(1.0e-10, f64::INFINITY).is_err());
        assert!(CollectiveStepSizes::scaled_to(0.0, 1.0e-9).is_err());
        assert!(CollectiveStepSizes::scaled_to(1.0e-9, 0.0).is_err());
    }

    #[test]
    fn test_with_profile_and_offset_zero_matches_with_profile() {
        // Regression/consistency check for the additive Hopfion API this
        // module relies on: a zero center offset must reproduce
        // `with_profile` exactly.
        let a = Hopfion::with_profile(TEST_GRID, 5.0e-9, 1, ProfileFunction::Linear)
            .expect("with_profile failed");
        let b = Hopfion::with_profile_and_offset(
            TEST_GRID,
            5.0e-9,
            1,
            ProfileFunction::Linear,
            Vector3::zero(),
        )
        .expect("with_profile_and_offset failed");

        assert_eq!(a.magnetization.len(), b.magnetization.len());
        for (ma, mb) in a.magnetization.iter().zip(b.magnetization.iter()) {
            assert!((ma.x - mb.x).abs() < 1e-15);
            assert!((ma.y - mb.y).abs() < 1e-15);
            assert!((ma.z - mb.z).abs() < 1e-15);
        }
    }

    #[test]
    fn test_with_profile_and_offset_rejects_non_finite_offset() {
        let result = Hopfion::with_profile_and_offset(
            TEST_GRID,
            5.0e-9,
            1,
            ProfileFunction::Linear,
            Vector3::new(f64::NAN, 0.0, 0.0),
        );
        assert!(result.is_err(), "non-finite center offset must be rejected");
    }

    #[test]
    fn test_dominant_mode_out_of_range_returns_none() {
        let hopfion =
            Hopfion::new(TEST_GRID, 3.0e-9, 1).expect("failed to build reference hopfion");
        let params = well_forming_energy_params();
        let steps = CollectiveStepSizes::for_hopfion(&hopfion).expect("valid step sizes");
        let spectrum = HopfionEigenmodeStability::analyze_standard_modes(&hopfion, &params, &steps)
            .expect("analysis should succeed");
        assert!(spectrum
            .dominant_mode(spectrum.eigenvalues.len() + 10)
            .is_none());
    }

    #[test]
    fn test_single_mode_radius_only_analysis() {
        // Exercise the flexible mode-list API with a single collective
        // coordinate (n=1), confirming the general n-mode machinery does
        // not require the standard 4-mode basis.
        let (_, r_eq, _) = find_boundaries();
        let params = well_forming_energy_params();
        let hopfion = Hopfion::with_profile(TEST_GRID, r_eq, 1, ProfileFunction::Linear)
            .expect("failed to build equilibrium-radius hopfion");
        let steps = CollectiveStepSizes::for_hopfion(&hopfion).expect("failed to build step sizes");

        let spectrum = HopfionEigenmodeStability::analyze(
            &hopfion,
            &params,
            &steps,
            &[CollectiveMode::Radius],
        )
        .expect("single-mode analysis should succeed");

        assert_eq!(spectrum.eigenvalues.len(), 1);
        assert_eq!(spectrum.modes, vec![CollectiveMode::Radius]);
        assert!(
            spectrum.eigenvalues[0] > 0.0,
            "radius-only curvature at the equilibrium radius must be positive, got {}",
            spectrum.eigenvalues[0]
        );
    }
}
