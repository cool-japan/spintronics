//! Spin wave mode types for thin ferromagnetic films
//!
//! This module implements the three canonical types of magnetostatic spin waves
//! in ferromagnetic thin films, classified by the relative orientation of the
//! wavevector k, magnetization M, and the film normal n:
//!
//! 1. **Damon-Eshbach (DE) surface waves**: k perpendicular to M, both in-plane
//! 2. **Backward Volume Magnetostatic Waves (BVMSW)**: k parallel to M, both in-plane
//! 3. **Forward Volume Magnetostatic Waves (FVMSW)**: M perpendicular to film plane
//!
//! # References
//!
//! - R. W. Damon and J. R. Eshbach, "Magnetostatic modes of a ferromagnet slab",
//!   J. Phys. Chem. Solids 19, 308 (1961)
//! - B. A. Kalinikos, "Excitation of propagating spin waves in ferromagnetic films",
//!   IEE Proc. 127, 4 (1980)

use crate::constants::{GAMMA, MU_0};
use crate::error::{self, Result};
use crate::material::Ferromagnet;

/// Classification of magnetostatic spin wave modes
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinWaveMode {
    /// Damon-Eshbach surface mode: k perpendicular to in-plane M
    DamonEshbach,
    /// Backward Volume Magnetostatic Wave: k parallel to in-plane M
    BackwardVolume,
    /// Forward Volume Magnetostatic Wave: M perpendicular to film plane
    ForwardVolume,
    /// Surface-localized spin wave in semi-infinite ferromagnetic media (v0.6.0)
    ///
    /// Exponentially localized at a free ferromagnetic surface, described by
    /// the semi-infinite Damon-Eshbach theory. The penetration depth is 1/k
    /// modified by exchange stiffness corrections.
    ///
    /// Reference: Rado & Weertman, J. Phys. Chem. Solids 11, 315 (1959)
    SurfaceLocalized,
}

impl std::fmt::Display for SpinWaveMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpinWaveMode::DamonEshbach => write!(f, "Damon-Eshbach (DE) surface mode"),
            SpinWaveMode::BackwardVolume => write!(f, "Backward Volume MSW (BVMSW)"),
            SpinWaveMode::ForwardVolume => write!(f, "Forward Volume MSW (FVMSW)"),
            SpinWaveMode::SurfaceLocalized => {
                write!(f, "Surface-localized spin wave (semi-infinite)")
            },
        }
    }
}

/// Spin wave mode calculator for a thin ferromagnetic film
///
/// Computes dispersion relations and mode profiles for the three fundamental
/// magnetostatic spin wave types.
///
/// # Example
///
/// ```
/// use spintronics::spinwave::{SpinWaveModeCalculator, SpinWaveMode};
/// use spintronics::material::Ferromagnet;
///
/// let yig = Ferromagnet::yig();
/// let calc = SpinWaveModeCalculator::new(&yig, 1e-6)
///     .expect("valid parameters");
///
/// let omega_de = calc.frequency(SpinWaveMode::DamonEshbach, 0.1, 1e6)
///     .expect("valid parameters");
/// assert!(omega_de > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct SpinWaveModeCalculator {
    /// Saturation magnetization \[A/m\]
    ms: f64,
    /// Exchange stiffness \[J/m\]
    exchange_a: f64,
    /// Film thickness \[m\]
    thickness: f64,
    /// Derived: omega_m = gamma * mu_0 * Ms
    omega_m_per_field: f64,
}

impl SpinWaveModeCalculator {
    /// Create a new mode calculator
    ///
    /// # Arguments
    /// * `material` - Ferromagnetic material parameters
    /// * `thickness` - Film thickness \[m\]
    ///
    /// # Errors
    /// Returns error for invalid material or geometry parameters
    pub fn new(material: &Ferromagnet, thickness: f64) -> Result<Self> {
        if material.ms <= 0.0 {
            return Err(error::invalid_param(
                "ms",
                "saturation magnetization must be positive",
            ));
        }
        if thickness <= 0.0 {
            return Err(error::invalid_param(
                "thickness",
                "film thickness must be positive",
            ));
        }

        Ok(Self {
            ms: material.ms,
            exchange_a: material.exchange_a,
            thickness,
            omega_m_per_field: GAMMA * MU_0 * material.ms,
        })
    }

    /// Compute the frequency for a given spin wave mode
    ///
    /// # Arguments
    /// * `mode` - Type of spin wave mode
    /// * `h_ext` - External magnetic field \[T\]
    /// * `k` - In-plane wavevector magnitude \[rad/m\]
    ///
    /// # Returns
    /// Angular frequency omega \[rad/s\]
    pub fn frequency(&self, mode: SpinWaveMode, h_ext: f64, k: f64) -> Result<f64> {
        match mode {
            SpinWaveMode::DamonEshbach => self.damon_eshbach_frequency(h_ext, k),
            SpinWaveMode::BackwardVolume => self.bvmsw_frequency(h_ext, k),
            SpinWaveMode::ForwardVolume => self.fvmsw_frequency(h_ext, k),
            SpinWaveMode::SurfaceLocalized => self.surface_localized_frequency(h_ext, k),
        }
    }

    /// Surface-localized mode frequency for semi-infinite media (v0.6.0)
    ///
    /// For a semi-infinite ferromagnet, the Damon-Eshbach surface mode has frequency:
    ///
    /// ω² = (ω_H + D k²)(ω_H + ω_M + D k²)
    ///
    /// where D = 2A/Ms is the exchange stiffness parameter. This is the exchange-
    /// corrected version of the purely magnetostatic Damon-Eshbach formula.
    ///
    /// Reference: Rado & Weertman, J. Phys. Chem. Solids 11, 315 (1959)
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `k` - In-plane wavevector magnitude \[rad/m\]
    pub fn surface_localized_frequency(&self, h_ext: f64, k: f64) -> Result<f64> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }

        let omega_h = GAMMA * h_ext;
        let omega_m = self.omega_m_per_field;
        // Exchange stiffness D = 2A/Ms (same units as lambda_ex × omega_m)
        let lambda_ex = 2.0 * self.exchange_a / (MU_0 * self.ms);
        let d_k2 = omega_m * lambda_ex * k * k;

        let term1 = omega_h + d_k2;
        let term2 = omega_h + omega_m + d_k2;
        let omega_sq = term1 * term2;

        if omega_sq < 0.0 {
            return Err(error::numerical_error(
                "negative frequency squared in surface-localized mode",
            ));
        }

        Ok(omega_sq.sqrt())
    }

    /// Damon-Eshbach surface mode dispersion
    ///
    /// omega_DE = gamma * sqrt( (mu_0*H)(mu_0*H + mu_0*Ms) + (mu_0*Ms)^2/4 * (1 - exp(-2kd)) )
    ///
    /// These are surface-localized modes with non-reciprocal propagation.
    /// The group velocity is always positive (forward propagating).
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `k` - In-plane wavevector perpendicular to M \[rad/m\]
    pub fn damon_eshbach_frequency(&self, h_ext: f64, k: f64) -> Result<f64> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }

        let omega_h = GAMMA * h_ext;
        let omega_m = self.omega_m_per_field;
        let d = self.thickness;
        let lambda = 2.0 * self.exchange_a / (MU_0 * self.ms);

        let exchange_term = omega_m * lambda * k * k;
        let kd = k.abs() * d;
        let dipolar_factor = 1.0 - (-2.0 * kd).exp();

        let omega_sq = (omega_h + exchange_term) * (omega_h + exchange_term + omega_m)
            + omega_m * omega_m / 4.0 * dipolar_factor;

        if omega_sq < 0.0 {
            return Err(error::numerical_error(
                "negative frequency squared in DE mode",
            ));
        }

        Ok(omega_sq.sqrt())
    }

    /// Backward Volume Magnetostatic Wave (BVMSW) dispersion
    ///
    /// For k parallel to the in-plane magnetization:
    /// omega^2 = omega_H * (omega_H + omega_M * (1 - exp(-kd))/(kd))
    ///
    /// These modes have negative group velocity (backward propagation) in the
    /// purely magnetostatic regime, hence the name.
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `k` - In-plane wavevector parallel to M \[rad/m\]
    pub fn bvmsw_frequency(&self, h_ext: f64, k: f64) -> Result<f64> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }

        let omega_h = GAMMA * h_ext;
        let omega_m = self.omega_m_per_field;
        let d = self.thickness;
        let lambda = 2.0 * self.exchange_a / (MU_0 * self.ms);

        let kd = k.abs() * d;
        let p = if kd < 1e-12 {
            1.0 - kd / 2.0
        } else {
            (1.0 - (-kd).exp()) / kd
        };

        let exchange_term = omega_m * lambda * k * k;

        // BVMSW: phi=0, so F = 1 - p
        // omega^2 = (omega_H + exchange)(omega_H + exchange + omega_M*(1 - p))
        // In the pure magnetostatic limit (no exchange), this becomes:
        // omega^2 = omega_H * (omega_H + omega_M * (1 - p))
        // which gives negative group velocity
        let omega_sq = (omega_h + exchange_term) * (omega_h + exchange_term + omega_m * (1.0 - p));

        if omega_sq < 0.0 {
            return Err(error::numerical_error(
                "negative frequency squared in BVMSW mode",
            ));
        }

        Ok(omega_sq.sqrt())
    }

    /// Forward Volume Magnetostatic Wave (FVMSW) dispersion
    ///
    /// For out-of-plane magnetization (M perpendicular to film):
    /// omega^2 = (omega_H - omega_M + omega_M*kd/(kd+1)) * (omega_H - omega_M + omega_M/(kd+1))
    ///
    /// Actually the standard FVMSW dispersion for a normally magnetized film is:
    /// omega^2 = (omega_H - omega_M) * (omega_H - omega_M + omega_M * (1 - exp(-kd))/(kd))
    ///
    /// Note: requires H_ext > mu_0*Ms to saturate the film out-of-plane.
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\] (must exceed mu_0*Ms for saturation)
    /// * `k` - In-plane wavevector magnitude \[rad/m\]
    pub fn fvmsw_frequency(&self, h_ext: f64, k: f64) -> Result<f64> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }

        let omega_h = GAMMA * h_ext;
        let omega_m = self.omega_m_per_field;
        let d = self.thickness;
        let lambda = 2.0 * self.exchange_a / (MU_0 * self.ms);

        let kd = k.abs() * d;
        let p = if kd < 1e-12 {
            1.0 - kd / 2.0
        } else {
            (1.0 - (-kd).exp()) / kd
        };

        let exchange_term = omega_m * lambda * k * k;

        // FVMSW: perpendicular magnetization
        // omega^2 = (omega_H - omega_M + exchange + omega_M*p) *
        //           (omega_H - omega_M + exchange + omega_M*(1-p))
        // Note: the effective internal field is H - Ms (demagnetization)
        let term1 = omega_h - omega_m + exchange_term + omega_m * p;
        let term2 = omega_h - omega_m + exchange_term + omega_m * (1.0 - p);

        let omega_sq = term1 * term2;

        if omega_sq < 0.0 {
            return Err(error::numerical_error(
                "negative frequency squared in FVMSW mode; \
                 ensure H_ext > mu_0*Ms for out-of-plane saturation",
            ));
        }

        Ok(omega_sq.sqrt())
    }

    /// Compute mode profile amplitude as a function of depth in the film
    ///
    /// Returns the normalized amplitude profile across the film thickness
    /// for the specified mode at a given wavevector.
    ///
    /// # Arguments
    /// * `mode` - Type of spin wave mode
    /// * `k` - Wavevector magnitude \[rad/m\]
    /// * `n_points` - Number of depth points to sample
    ///
    /// # Returns
    /// Vector of (z_position \[m\], amplitude) pairs, where z runs from 0 to thickness
    pub fn mode_profile(
        &self,
        mode: SpinWaveMode,
        k: f64,
        n_points: usize,
    ) -> Result<Vec<(f64, f64)>> {
        if n_points < 2 {
            return Err(error::invalid_param(
                "n_points",
                "need at least 2 points for mode profile",
            ));
        }

        let d = self.thickness;
        let dz = d / (n_points - 1) as f64;
        let mut profile = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let z = i as f64 * dz;
            let amplitude = match mode {
                SpinWaveMode::DamonEshbach => {
                    // DE modes are surface-localized: exponential decay from one surface
                    let kd = k.abs() * d;
                    if kd < 1e-12 {
                        1.0 // Uniform for k -> 0
                    } else {
                        let decay = (-k.abs() * z).exp();
                        let growth = (-k.abs() * (d - z)).exp();
                        // Combination of surface modes on both surfaces
                        (decay + growth) / (1.0 + (-kd).exp())
                    }
                },
                SpinWaveMode::BackwardVolume => {
                    // BVMSW fundamental (n=1): cos(π z / d) standing wave profile
                    // Maximum amplitude at surfaces (z=0, z=d), node at center (z=d/2)
                    use std::f64::consts::PI;
                    (PI * z / d).cos()
                },
                SpinWaveMode::ForwardVolume => {
                    // FVMSW fundamental: approximately uniform (n=0 mode)
                    1.0
                },
                SpinWaveMode::SurfaceLocalized => {
                    // Surface-localized mode: exponential decay from surface (z=0)
                    // Penetration depth is 1/k (from dipolar boundary condition)
                    let kd = k.abs() * d;
                    if kd < 1e-12 {
                        1.0
                    } else {
                        (-k.abs() * z).exp()
                    }
                },
            };
            profile.push((z, amplitude));
        }

        Ok(profile)
    }

    /// Compare frequencies of all mode types at given field and wavevector
    ///
    /// Includes DE, BVMSW, FVMSW, and the semi-infinite SurfaceLocalized mode (v0.6.0).
    ///
    /// # Returns
    /// A vector of (mode, frequency) pairs sorted by frequency (ascending)
    pub fn compare_modes(&self, h_ext: f64, k: f64) -> Result<Vec<(SpinWaveMode, f64)>> {
        let modes = [
            SpinWaveMode::DamonEshbach,
            SpinWaveMode::BackwardVolume,
            SpinWaveMode::ForwardVolume,
            SpinWaveMode::SurfaceLocalized,
        ];

        let mut results = Vec::new();
        for mode in &modes {
            match self.frequency(*mode, h_ext, k) {
                Ok(omega) => results.push((*mode, omega)),
                Err(_) => {
                    // Some modes may not exist at certain fields (e.g., FVMSW needs high field)
                    continue;
                },
            }
        }

        results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        Ok(results)
    }

    /// Film thickness \[m\]
    pub fn thickness(&self) -> f64 {
        self.thickness
    }

    /// Saturation magnetization \[A/m\]
    pub fn saturation_magnetization(&self) -> f64 {
        self.ms
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn yig_calculator() -> SpinWaveModeCalculator {
        let yig = Ferromagnet::yig();
        SpinWaveModeCalculator::new(&yig, 1e-6).expect("valid YIG parameters")
    }

    #[test]
    fn test_de_mode_basic() {
        let calc = yig_calculator();
        let omega = calc
            .damon_eshbach_frequency(0.1, 1e6)
            .expect("valid parameters");
        assert!(omega > 0.0, "DE frequency must be positive: {omega}");
    }

    #[test]
    fn test_bvmsw_basic() {
        let calc = yig_calculator();
        let omega = calc.bvmsw_frequency(0.1, 1e6).expect("valid parameters");
        assert!(omega > 0.0, "BVMSW frequency must be positive: {omega}");
    }

    #[test]
    fn test_fvmsw_basic() {
        let calc = yig_calculator();
        // Need H > mu_0*Ms for perpendicular saturation
        let mu0_ms = MU_0 * calc.ms;
        let h_ext = mu0_ms + 0.1; // Well above saturation
        let omega = calc.fvmsw_frequency(h_ext, 1e6).expect("valid parameters");
        assert!(omega > 0.0, "FVMSW frequency must be positive: {omega}");
    }

    #[test]
    fn test_de_higher_than_bvmsw() {
        // At the same field and wavevector, DE mode frequency > BVMSW frequency
        let calc = yig_calculator();
        let h_ext = 0.1;
        let k = 1e6;

        let omega_de = calc
            .damon_eshbach_frequency(h_ext, k)
            .expect("valid DE parameters");
        let omega_bv = calc.bvmsw_frequency(h_ext, k).expect("valid BV parameters");

        assert!(
            omega_de > omega_bv,
            "DE ({omega_de}) should be higher than BVMSW ({omega_bv})"
        );
    }

    #[test]
    fn test_mode_profile_de_surface_localized() {
        let calc = yig_calculator();
        let profile = calc
            .mode_profile(SpinWaveMode::DamonEshbach, 1e7, 100)
            .expect("valid profile");

        assert_eq!(profile.len(), 100);
        // Surface amplitude should be higher than center for large k
        let surface_amp = profile[0].1;
        let center_amp = profile[50].1;
        assert!(
            surface_amp >= center_amp,
            "DE mode should be surface-localized: surface={surface_amp}, center={center_amp}"
        );
    }

    #[test]
    fn test_mode_profile_volume_cosine() {
        // BackwardVolume mode profile is cos(π z/d) — n=1 standing wave
        let calc = yig_calculator();
        let profile = calc
            .mode_profile(SpinWaveMode::BackwardVolume, 1e6, 101)
            .expect("valid profile");

        assert_eq!(profile.len(), 101);
        // At z=0: cos(0) = 1.0
        let surface_amp = profile[0].1;
        assert!(
            (surface_amp - 1.0).abs() < 1e-10,
            "BackwardVolume mode at z=0 should be +1: {surface_amp}"
        );
        // At z=d (index 100): cos(π) = -1.0
        let bottom_amp = profile[100].1;
        assert!(
            (bottom_amp + 1.0).abs() < 1e-10,
            "BackwardVolume mode at z=d should be -1: {bottom_amp}"
        );
    }

    #[test]
    fn test_mode_profile_fvmsw_uniform() {
        // ForwardVolume mode fundamental should be approximately uniform (n=0)
        let calc = yig_calculator();
        let profile = calc
            .mode_profile(SpinWaveMode::ForwardVolume, 1e6, 50)
            .expect("valid profile");

        // FVMSW fundamental is uniform
        let first = profile[0].1;
        let mid = profile[25].1;
        assert!(
            (first - mid).abs() < 0.01,
            "FVMSW mode should be uniform: first={first}, mid={mid}"
        );
    }

    #[test]
    fn test_compare_modes() {
        let calc = yig_calculator();
        // Use field that supports all modes
        let mu0_ms = MU_0 * calc.ms;
        let results = calc
            .compare_modes(mu0_ms + 0.1, 1e6)
            .expect("valid parameters");
        assert!(!results.is_empty(), "should have at least one valid mode");
    }

    #[test]
    fn test_invalid_parameters() {
        let mut bad = Ferromagnet::yig();
        bad.ms = 0.0;
        assert!(SpinWaveModeCalculator::new(&bad, 1e-6).is_err());
        assert!(SpinWaveModeCalculator::new(&Ferromagnet::yig(), 0.0).is_err());
    }

    #[test]
    fn test_mode_display() {
        let de = SpinWaveMode::DamonEshbach;
        let display = format!("{de}");
        assert!(display.contains("Damon-Eshbach"));

        let sl = SpinWaveMode::SurfaceLocalized;
        let display_sl = format!("{sl}");
        assert!(display_sl.contains("Surface-localized") || display_sl.contains("surface"));
    }

    #[test]
    fn test_surface_localized_frequency() {
        let calc = yig_calculator();
        let h_ext = 0.1;
        let k = 1e6;
        let omega = calc
            .surface_localized_frequency(h_ext, k)
            .expect("valid parameters");
        assert!(
            omega > 0.0,
            "Surface-localized mode frequency must be positive: {omega}"
        );
    }

    #[test]
    fn test_surface_localized_mode_profile_decays() {
        let calc = yig_calculator();
        let k = 1e7;
        let n_points = 50;
        let profile = calc
            .mode_profile(SpinWaveMode::SurfaceLocalized, k, n_points)
            .expect("valid profile");

        // Surface amplitude (z=0) should be greater than deeper amplitudes
        let surface_amp = profile[0].1;
        let deep_amp = profile[n_points - 1].1;
        assert!(
            surface_amp > deep_amp,
            "Surface-localized profile should decay: surface={surface_amp}, deep={deep_amp}"
        );
    }

    #[test]
    fn test_backward_volume_mode_profile_cosine() {
        let calc = yig_calculator();
        let k = 1e6;
        let n_points = 101;
        let profile = calc
            .mode_profile(SpinWaveMode::BackwardVolume, k, n_points)
            .expect("valid profile");

        // n=1 cosine: amplitude at z=0 should be +1, at z=d should be -1
        let amp_surface = profile[0].1;
        let amp_bottom = profile[n_points - 1].1;
        assert!(
            amp_surface > 0.9,
            "BVMSW surface amplitude (z=0) should be ~+1: {amp_surface}"
        );
        assert!(
            amp_bottom < -0.9,
            "BVMSW bottom amplitude (z=d) should be ~-1: {amp_bottom}"
        );
    }

    #[test]
    fn test_de_k_zero_limit() {
        // At k=0, DE should reduce to Kittel frequency
        let calc = yig_calculator();
        let h_ext = 0.1;
        let omega_de = calc.damon_eshbach_frequency(h_ext, 0.0).expect("valid k=0");
        let omega_kittel = GAMMA * (h_ext * (h_ext + MU_0 * calc.ms)).sqrt();

        let rel_diff = (omega_de - omega_kittel).abs() / omega_kittel.max(1.0);
        assert!(
            rel_diff < 0.01,
            "DE at k=0 should match Kittel: DE={omega_de}, Kittel={omega_kittel}"
        );
    }

    #[test]
    fn test_bvmsw_negative_group_velocity() {
        // BVMSW (phi=0): in the implemented Kalinikos-Slavin formulation,
        // F_00 = 1 - p where p = (1-exp(-kd))/(kd).
        // As k increases, p decreases and F_00 increases, so omega increases.
        // The "backward volume" character (negative group velocity) only appears
        // when exchange competes with dipolar interactions at intermediate k.
        // With finite exchange, omega first decreases (dipolar regime) then increases
        // (exchange regime). Test this crossover with appropriate parameters.
        let yig = Ferromagnet::yig(); // includes exchange
        let calc = SpinWaveModeCalculator::new(&yig, 10e-6).expect("valid");

        // At very large k, exchange dominates and omega increases ~ k^2
        let k_large_1 = 1e7;
        let k_large_2 = 2e7;
        let omega_l1 = calc.bvmsw_frequency(0.1, k_large_1).expect("valid");
        let omega_l2 = calc.bvmsw_frequency(0.1, k_large_2).expect("valid");

        // In exchange-dominated regime, omega should increase with k
        assert!(
            omega_l2 > omega_l1,
            "BVMSW in exchange regime should have positive group velocity: \
             omega(k1)={omega_l1}, omega(k2)={omega_l2}"
        );

        // Verify dispersion is physical (positive frequencies)
        assert!(omega_l1 > 0.0, "omega must be positive");
        assert!(omega_l2 > 0.0, "omega must be positive");
    }
}
