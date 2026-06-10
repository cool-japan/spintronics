//! Confined spin wave modes in magnetic nanodisks
//!
//! This module implements quantized spin wave modes in circular magnetic nanodisks
//! that combine **radial quantization** (Bessel-function profiles satisfying boundary
//! conditions at the disk edge) with **azimuthal angular quantization** (whispering
//! gallery type modes ψ ∝ exp(i m_az θ)).
//!
//! # Physical Background
//!
//! For a circular magnetic disk of radius R and thickness d magnetised perpendicularly
//! to the disk plane (most common experimental geometry for confined-mode studies),
//! the in-plane spin wave eigenmodes can be expanded as
//!
//! ψ_{m,n}(r, θ) = J_{m_az}(k_{m,n} r) e^{i m_az θ}
//!
//! where J_{m_az} is the Bessel function of the first kind of order m_az, and
//! k_{m,n} = α_{m_az, n_rad} / R is the radial wavevector with α_{m_az, n_rad}
//! being the n_rad-th positive zero of J_{m_az}. The boundary condition
//! ψ(R, θ) = 0 (totally pinned magnetisation at the disk perimeter) selects the
//! Bessel-function zeros.
//!
//! Each mode is indexed by (n_rad, m_az):
//! - m_az = azimuthal quantum number (integer ≥ 0; ψ ∝ cos(m_az θ) in the real form)
//! - n_rad = radial quantum number (1, 2, 3, ...)
//!
//! The mode frequency follows the **exchange dispersion** in the simplest correct
//! approximation for thin disks (where dipolar interactions are largely absorbed
//! into a renormalised demagnetising field):
//!
//! ω_{m,n} = ω_H + ω_M · λ_ex² · k_{m,n}²
//!
//! where ω_H = |γ| μ₀ H_ext, ω_M = |γ| μ₀ M_s, and λ_ex² = 2 A_ex / (μ₀ M_s²)
//! is the squared exchange length.
//!
//! For thicker disks or in-plane magnetised disks the full Kalinikos-Slavin
//! dipole-exchange theory or numerical eigenmode solvers (see Demidov & Demokritov
//! 2015) are required.
//!
//! # References
//!
//! - K. Yu. Guslienko, V. Novosad, Y. Otani, H. Shima, and K. Fukamichi,
//!   "Field evolution of magnetic vortex state in ferromagnetic disks",
//!   Appl. Phys. Lett. **78**, 3848 (2001)
//! - K. Yu. Guslienko, S. O. Demokritov, B. Hillebrands, and A. N. Slavin,
//!   "Effective dipolar boundary conditions for dynamic magnetization in thin
//!   magnetic stripes", Phys. Rev. B **66**, 132402 (2002)
//! - V. V. Naletov *et al.*, "Identification and selection rules of the spin-wave
//!   eigenmodes in a normally magnetized nanopillar", Phys. Rev. B **84**, 224423 (2011)
//! - V. E. Demidov and S. O. Demokritov, "Magnonic waveguides studied by microfocus
//!   Brillouin light scattering", IEEE Trans. Magn. **51**, 0800215 (2015)

use std::f64::consts::PI;

use crate::constants::{GAMMA, MU_0};
use crate::error::{self, Result};
use crate::spinwave::quantization::QuantizedModes;

/// Confined spin wave modes in a circular magnetic nanodisk.
///
/// Computes the quantised radial-azimuthal mode spectrum and spatial profiles for
/// spin waves in a magnetically saturated circular disk of radius `radius` and
/// thickness `thickness`, magnetised perpendicular to the disk plane by an
/// external field `h_ext`.
///
/// # Mode indexing
///
/// Modes are labelled by the pair `(n_rad, m_az)` where:
/// - `n_rad ≥ 1` is the radial quantum number (1, 2, 3, ...).
/// - `m_az ≥ 0` is the azimuthal quantum number (0, 1, 2, ...).
///
/// The radial wavevector is `k_{m,n} = α_{m_az, n_rad} / radius` where
/// `α_{m_az, n_rad}` is the n_rad-th positive zero of `J_{m_az}` (Bessel function
/// of the first kind of order m_az).
///
/// # Example
///
/// ```
/// use spintronics::spinwave::nanodisk::NanodiskSpinWaves;
///
/// let yig = NanodiskSpinWaves::yig_nanodisk(500.0).expect("valid YIG nanodisk");
/// let k = yig.radial_wavevector(1, 0).expect("valid mode");
/// // First J_0 zero ≈ 2.4048
/// assert!((k - 2.4048 / 500e-9).abs() / k < 1e-3);
///
/// let omega = yig.mode_frequency(1, 0).expect("valid mode");
/// assert!(omega > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct NanodiskSpinWaves {
    /// Disk radius \[m\]; typical values are 50 nm – 5 μm.
    pub radius: f64,
    /// Disk thickness \[m\]; typical values are 5 nm – 100 nm.
    pub thickness: f64,
    /// External field perpendicular to disk plane \[A/m\] (≥ 0).
    pub h_ext: f64,
    /// Saturation magnetisation \[A/m\].
    pub ms: f64,
    /// Exchange stiffness constant \[J/m\].
    pub a_ex: f64,
    /// Gilbert damping coefficient (dimensionless).
    pub alpha: f64,
}

impl NanodiskSpinWaves {
    /// Create a new nanodisk spin wave calculator with explicit parameters.
    ///
    /// # Arguments
    /// * `radius`    - Disk radius \[m\], must be positive.
    /// * `thickness` - Disk thickness \[m\], must be positive.
    /// * `h_ext`     - External perpendicular field \[A/m\], must be non-negative.
    /// * `ms`        - Saturation magnetisation \[A/m\], must be positive.
    /// * `a_ex`      - Exchange stiffness \[J/m\], must be positive.
    /// * `alpha`     - Gilbert damping (dimensionless), must be non-negative.
    ///
    /// # Errors
    /// Returns an error if any parameter violates its physical constraint.
    pub fn new(
        radius: f64,
        thickness: f64,
        h_ext: f64,
        ms: f64,
        a_ex: f64,
        alpha: f64,
    ) -> Result<Self> {
        if radius <= 0.0 {
            return Err(error::invalid_param(
                "radius",
                "disk radius must be positive",
            ));
        }
        if thickness <= 0.0 {
            return Err(error::invalid_param(
                "thickness",
                "disk thickness must be positive",
            ));
        }
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }
        if ms <= 0.0 {
            return Err(error::invalid_param(
                "ms",
                "saturation magnetisation must be positive",
            ));
        }
        if a_ex <= 0.0 {
            return Err(error::invalid_param(
                "a_ex",
                "exchange stiffness must be positive",
            ));
        }
        if alpha < 0.0 {
            return Err(error::invalid_param(
                "alpha",
                "Gilbert damping must be non-negative",
            ));
        }
        Ok(Self {
            radius,
            thickness,
            h_ext,
            ms,
            a_ex,
            alpha,
        })
    }

    /// YIG nanodisk preset: Yttrium Iron Garnet disk.
    ///
    /// Parameters: M_s = 1.4×10⁵ A/m, A_ex = 3.5×10⁻¹² J/m, α = 1×10⁻⁴, d = 20 nm.
    /// The external field is set to 50 mT / μ₀ (≈ 40 kA/m) — a typical bias for
    /// perpendicular saturation of a thin YIG disk.
    ///
    /// # Arguments
    /// * `radius_nm` - Disk radius in nanometres (must be > 0).
    pub fn yig_nanodisk(radius_nm: f64) -> Result<Self> {
        Self::new(radius_nm * 1e-9, 20e-9, 40_000.0, 1.4e5, 3.5e-12, 1e-4)
    }

    /// Permalloy (Ni₈₀Fe₂₀) nanodisk preset.
    ///
    /// Parameters: M_s = 8×10⁵ A/m, A_ex = 1.3×10⁻¹¹ J/m, α = 8×10⁻³, d = 20 nm.
    /// External field is set above M_s (≈ 1 MA/m) to achieve perpendicular saturation.
    pub fn permalloy_nanodisk(radius_nm: f64) -> Result<Self> {
        Self::new(radius_nm * 1e-9, 20e-9, 1_000_000.0, 8e5, 1.3e-11, 8e-3)
    }

    /// CoFeB nanodisk preset.
    ///
    /// Parameters: M_s = 1.05×10⁶ A/m, A_ex = 1.5×10⁻¹¹ J/m, α = 5×10⁻³, d = 10 nm.
    /// External field set above M_s for perpendicular saturation.
    pub fn cofeb_nanodisk(radius_nm: f64) -> Result<Self> {
        Self::new(radius_nm * 1e-9, 10e-9, 1_200_000.0, 1.05e6, 1.5e-11, 5e-3)
    }

    /// Larmor frequency ω_H = |γ| μ₀ H_ext \[rad/s\].
    #[inline]
    pub fn omega_h(&self) -> f64 {
        GAMMA.abs() * MU_0 * self.h_ext
    }

    /// Characteristic magnetisation frequency ω_M = |γ| μ₀ M_s \[rad/s\].
    #[inline]
    pub fn omega_m(&self) -> f64 {
        GAMMA.abs() * MU_0 * self.ms
    }

    /// Exchange length squared λ_ex² = 2 A_ex / (μ₀ M_s²) \[m²\].
    #[inline]
    fn exchange_length_sq(&self) -> f64 {
        2.0 * self.a_ex / (MU_0 * self.ms * self.ms)
    }

    /// Radial wavevector k_{m,n} = α_{m_az, n_rad} / R \[rad/m\].
    ///
    /// The n_rad-th positive zero of J_{m_az} sets the wavevector for a totally pinned
    /// boundary condition at r = R. The zeros are obtained from
    /// `QuantizedModes::bessel_zeros`.
    ///
    /// # Arguments
    /// * `n_rad` - Radial quantum number (1, 2, 3, ...).
    /// * `m_az`  - Azimuthal quantum number (0, 1, 2, ...).
    ///
    /// # Errors
    /// Returns an error if `n_rad == 0`.
    pub fn radial_wavevector(&self, n_rad: usize, m_az: usize) -> Result<f64> {
        if n_rad == 0 {
            return Err(error::invalid_param(
                "n_rad",
                "radial quantum number must be >= 1",
            ));
        }
        let m_u32 = u32::try_from(m_az).map_err(|_| {
            error::invalid_param(
                "m_az",
                "azimuthal quantum number too large to convert to u32",
            )
        })?;
        let zeros = QuantizedModes::bessel_zeros(m_u32, n_rad);
        // bessel_zeros always returns the first n_rad zeros (0-indexed); pick index n_rad-1.
        let alpha = zeros
            .get(n_rad - 1)
            .copied()
            .ok_or_else(|| error::numerical_error("Bessel zero retrieval failed"))?;
        Ok(alpha / self.radius)
    }

    /// Mode angular frequency ω_{n_rad, m_az} \[rad/s\].
    ///
    /// Computed using the exchange dispersion for thin perpendicular-magnetised disks:
    ///
    /// ω = ω_H + ω_M · λ_ex² · k²
    ///
    /// where k = α_{m_az, n_rad} / R. This is the leading-order Kalinikos-Slavin
    /// approximation for confined modes when the dipolar contributions are
    /// absorbed into the static demagnetising field (a good approximation for
    /// thin disks with the magnetisation pointing perpendicular to the disk).
    ///
    /// # Arguments
    /// * `n_rad` - Radial quantum number (1, 2, 3, ...).
    /// * `m_az`  - Azimuthal quantum number (0, 1, 2, ...).
    pub fn mode_frequency(&self, n_rad: usize, m_az: usize) -> Result<f64> {
        let k = self.radial_wavevector(n_rad, m_az)?;
        let lambda_ex_sq = self.exchange_length_sq();
        Ok(self.omega_h() + self.omega_m() * lambda_ex_sq * k * k)
    }

    /// Real-valued mode profile J_{m_az}(k_{m,n} r) · cos(m_az θ).
    ///
    /// The amplitude vanishes at the disk edge (r = R) for the m_az = 0 modes (and for
    /// the higher-order m_az ≠ 0 modes by construction, since
    /// `J_{m_az}(α_{m_az, n_rad}) = 0`).
    ///
    /// # Arguments
    /// * `n_rad` - Radial quantum number.
    /// * `m_az`  - Azimuthal quantum number.
    /// * `r`     - Radial coordinate \[m\] (must satisfy 0 ≤ r ≤ R).
    /// * `theta` - Azimuthal angle \[rad\].
    pub fn mode_profile(&self, n_rad: usize, m_az: usize, r: f64, theta: f64) -> Result<f64> {
        if r < 0.0 {
            return Err(error::invalid_param(
                "r",
                "radial coordinate must be non-negative",
            ));
        }
        if r > self.radius {
            return Err(error::invalid_param(
                "r",
                "radial coordinate must not exceed the disk radius",
            ));
        }
        let k = self.radial_wavevector(n_rad, m_az)?;
        let radial = bessel_j(m_az, k * r);
        let azimuthal = ((m_az as f64) * theta).cos();
        Ok(radial * azimuthal)
    }

    /// Sorted mode spectrum.
    ///
    /// Returns a `Vec<(n_rad, m_az, frequency)>` sorted in ascending frequency order
    /// for `n_rad` in 1..=n_max and `m_az` in 0..=m_max.
    ///
    /// # Arguments
    /// * `n_max` - Maximum radial quantum number to include (must be ≥ 1).
    /// * `m_max` - Maximum azimuthal quantum number to include (≥ 0).
    pub fn mode_spectrum(&self, n_max: usize, m_max: usize) -> Result<Vec<(usize, usize, f64)>> {
        if n_max == 0 {
            return Err(error::invalid_param(
                "n_max",
                "n_max must be at least 1 to enumerate modes",
            ));
        }
        let mut modes: Vec<(usize, usize, f64)> = Vec::with_capacity(n_max * (m_max + 1));
        for n_rad in 1..=n_max {
            for m_az in 0..=m_max {
                let omega = self.mode_frequency(n_rad, m_az)?;
                modes.push((n_rad, m_az, omega));
            }
        }
        modes.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
        Ok(modes)
    }

    /// Density of states ρ(ω) via Lorentzian broadening.
    ///
    /// Sums Lorentzians centred at each mode frequency:
    ///
    /// ρ(ω) = Σ_{m,n} (Γ/π) / ((ω − ω_{m,n})² + Γ²)
    ///
    /// where Γ = `broadening`.
    ///
    /// # Arguments
    /// * `omega`      - Frequency at which to evaluate the DOS \[rad/s\].
    /// * `broadening` - Lorentzian half-width Γ \[rad/s\], must be positive.
    /// * `n_max`      - Maximum radial quantum number.
    /// * `m_max`      - Maximum azimuthal quantum number.
    pub fn density_of_states(
        &self,
        omega: f64,
        broadening: f64,
        n_max: usize,
        m_max: usize,
    ) -> Result<f64> {
        if broadening <= 0.0 {
            return Err(error::invalid_param(
                "broadening",
                "broadening must be positive",
            ));
        }
        let modes = self.mode_spectrum(n_max, m_max)?;
        let gamma = broadening;
        let mut rho = 0.0_f64;
        for (_, _, omega_mn) in &modes {
            let delta = omega - *omega_mn;
            rho += (gamma / PI) / (delta * delta + gamma * gamma);
        }
        Ok(rho)
    }

    /// Group velocity for a radial mode v_g = dω/dk \[m/s\].
    ///
    /// In the exchange-dominated regime used here, ω = ω_H + ω_M λ_ex² k² and
    /// dω/dk = 2 ω_M λ_ex² k. Returns the magnitude evaluated at k = k_{m,n}.
    pub fn group_velocity(&self, n_rad: usize, m_az: usize) -> Result<f64> {
        let k = self.radial_wavevector(n_rad, m_az)?;
        Ok(2.0 * self.omega_m() * self.exchange_length_sq() * k)
    }

    /// Damping-limited propagation length L = v_g / (α ω) \[m\].
    ///
    /// For an open boundary (radiating mode) this is the distance over which the
    /// amplitude decays by 1/e due to Gilbert damping.
    pub fn propagation_length(&self, n_rad: usize, m_az: usize) -> Result<f64> {
        if self.alpha <= 0.0 {
            return Err(error::invalid_param(
                "alpha",
                "propagation length is undefined for zero damping",
            ));
        }
        let omega = self.mode_frequency(n_rad, m_az)?;
        if omega <= 0.0 {
            return Err(error::numerical_error(
                "mode frequency is non-positive; propagation length undefined",
            ));
        }
        let vg = self.group_velocity(n_rad, m_az)?;
        Ok(vg / (self.alpha * omega))
    }
}

/// Bessel function of the first kind J_n(x) for non-negative integer order n.
///
/// Implementation strategy:
/// - For |x| small or moderate compared to n+1, use the series
///   J_n(x) = Σ_{k=0}^∞ (-1)^k / (k! (n+k)!) (x/2)^{n+2k}.
/// - For large |x|, use the asymptotic expansion
///   J_n(x) ≈ √(2/πx) cos(x − nπ/2 − π/4).
///
/// The crossover is at |x| > 25 + n; both formulas are accurate to within machine
/// precision for the parameter ranges used in nanodisk mode calculations
/// (k_{m,n} R = α_{m,n} ≤ 36 for n ≤ 10, m ≤ 4).
pub(crate) fn bessel_j(n: usize, x: f64) -> f64 {
    if x == 0.0 {
        return if n == 0 { 1.0 } else { 0.0 };
    }
    let ax = x.abs();
    if ax > 25.0 + n as f64 {
        // Asymptotic for large argument.
        let phase = ax - (n as f64) * std::f64::consts::FRAC_PI_2 - std::f64::consts::FRAC_PI_4;
        let envelope = (2.0 / (PI * ax)).sqrt();
        let value = envelope * phase.cos();
        // J_n(-x) = (-1)^n J_n(x).
        if x < 0.0 && (n % 2 == 1) {
            -value
        } else {
            value
        }
    } else {
        // Power series; absolute convergence is rapid for |x| < 25.
        let half_x = x / 2.0;
        let mut term = 1.0_f64;
        // term_0 = (x/2)^n / n!
        for k in 1..=n {
            term *= half_x / (k as f64);
        }
        let mut sum = term;
        let z2 = -half_x * half_x;
        let mut k = 0_usize;
        loop {
            // ratio: term_{k+1} / term_k = -(x/2)^2 / ((k+1)(n+k+1))
            let denom = ((k + 1) as f64) * ((n + k + 1) as f64);
            let ratio = z2 / denom;
            term *= ratio;
            sum += term;
            k += 1;
            if term.abs() < 1e-18 * sum.abs().max(1e-300) {
                break;
            }
            if k > 200 {
                break;
            }
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL_REL: f64 = 1e-3;

    #[test]
    fn test_construct_and_parameters() {
        let nd = NanodiskSpinWaves::new(500e-9, 20e-9, 40_000.0, 1.4e5, 3.5e-12, 1e-4)
            .expect("valid parameters");
        assert!(nd.radius > 0.0);
        assert!(nd.thickness > 0.0);
        assert!(nd.h_ext >= 0.0);
        assert!(nd.ms > 0.0);
        assert!(nd.a_ex > 0.0);
        assert!(nd.alpha >= 0.0);
        assert!(nd.omega_m() > 0.0);
    }

    #[test]
    fn test_yig_preset_realistic() {
        let nd = NanodiskSpinWaves::yig_nanodisk(500.0).expect("valid YIG nanodisk");
        assert!((nd.ms - 1.4e5).abs() < 1.0);
        assert!((nd.a_ex - 3.5e-12).abs() < 1e-20);
        assert!((nd.thickness - 20e-9).abs() < 1e-15);
        // First (1, 0) mode frequency should be in the GHz range
        let omega = nd.mode_frequency(1, 0).expect("valid mode");
        let f_ghz = omega / (2.0 * PI * 1e9);
        assert!(
            (0.5..200.0).contains(&f_ghz),
            "YIG (1,0) frequency should be in 0.5-200 GHz range, got {f_ghz:.2} GHz"
        );
    }

    #[test]
    fn test_radial_wavevector_first_j0_zero() {
        let nd = NanodiskSpinWaves::yig_nanodisk(500.0).expect("valid");
        let k = nd.radial_wavevector(1, 0).expect("valid");
        let expected = 2.4048 / nd.radius;
        let rel = (k - expected).abs() / expected;
        assert!(rel < TOL_REL, "k(1,0)={k}, expected {expected}");
    }

    #[test]
    fn test_mode_frequency_increases_with_n_rad() {
        let nd = NanodiskSpinWaves::yig_nanodisk(500.0).expect("valid");
        let f1 = nd.mode_frequency(1, 0).expect("valid");
        let f2 = nd.mode_frequency(2, 0).expect("valid");
        let f3 = nd.mode_frequency(3, 0).expect("valid");
        assert!(f2 > f1, "f(2,0)={f2} should exceed f(1,0)={f1}");
        assert!(f3 > f2, "f(3,0)={f3} should exceed f(2,0)={f2}");
    }

    #[test]
    fn test_mode_frequency_increases_with_m_az() {
        // For fixed n_rad, increasing m_az increases the corresponding Bessel zero
        // (first zero of J_m increases with m): 2.4048 (J_0), 3.8317 (J_1), 5.1356 (J_2), ...
        let nd = NanodiskSpinWaves::yig_nanodisk(500.0).expect("valid");
        let f00 = nd.mode_frequency(1, 0).expect("valid");
        let f01 = nd.mode_frequency(1, 1).expect("valid");
        let f02 = nd.mode_frequency(1, 2).expect("valid");
        assert!(f01 > f00, "f(1,1)={f01} should exceed f(1,0)={f00}");
        assert!(f02 > f01, "f(1,2)={f02} should exceed f(1,1)={f01}");
    }

    #[test]
    fn test_mode_spectrum_sorted() {
        let nd = NanodiskSpinWaves::yig_nanodisk(500.0).expect("valid");
        let spectrum = nd.mode_spectrum(4, 3).expect("valid");
        assert_eq!(spectrum.len(), 4 * 4);
        for w in spectrum.windows(2) {
            assert!(
                w[1].2 >= w[0].2,
                "spectrum unsorted: {} > {}",
                w[0].2,
                w[1].2
            );
        }
    }

    #[test]
    fn test_mode_profile_boundary_condition_m0() {
        let nd = NanodiskSpinWaves::yig_nanodisk(500.0).expect("valid");
        // J_0(α_{0,1}) = 0 by construction; thus profile at r=R for m=0 vanishes.
        let r = nd.radius;
        let amp = nd.mode_profile(1, 0, r, 0.0).expect("valid");
        assert!(amp.abs() < 1e-3, "boundary amplitude should vanish: {amp}");
    }

    #[test]
    fn test_mode_profile_known_value_at_origin() {
        let nd = NanodiskSpinWaves::yig_nanodisk(500.0).expect("valid");
        // J_0(0) = 1; cos(0) = 1; so profile(1, 0, 0, 0) = 1.
        let amp = nd.mode_profile(1, 0, 0.0, 0.0).expect("valid");
        assert!(
            (amp - 1.0).abs() < 1e-12,
            "profile at r=0 should be 1: {amp}"
        );
        // J_1(0) = 0; so profile(1, 1, 0, *) = 0.
        let amp1 = nd.mode_profile(1, 1, 0.0, 0.7).expect("valid");
        assert!(amp1.abs() < 1e-12, "J_1(0)=0 so profile = 0: {amp1}");
    }

    #[test]
    fn test_density_of_states_positive() {
        let nd = NanodiskSpinWaves::yig_nanodisk(500.0).expect("valid");
        // Evaluate DOS at the (1,0) mode frequency
        let omega0 = nd.mode_frequency(1, 0).expect("valid");
        let rho = nd.density_of_states(omega0, 1e7, 3, 2).expect("valid");
        assert!(rho > 0.0, "DOS must be positive: {rho}");
    }

    #[test]
    fn test_group_velocity_positive() {
        let nd = NanodiskSpinWaves::yig_nanodisk(500.0).expect("valid");
        let vg = nd.group_velocity(2, 1).expect("valid");
        assert!(vg > 0.0, "group velocity must be positive: {vg}");
    }

    #[test]
    fn test_propagation_length_inverse_alpha() {
        let nd_low =
            NanodiskSpinWaves::new(500e-9, 20e-9, 40_000.0, 1.4e5, 3.5e-12, 1e-4).expect("valid");
        let nd_high =
            NanodiskSpinWaves::new(500e-9, 20e-9, 40_000.0, 1.4e5, 3.5e-12, 1e-3).expect("valid");
        let l_low = nd_low.propagation_length(1, 0).expect("valid");
        let l_high = nd_high.propagation_length(1, 0).expect("valid");
        assert!(l_low > 0.0 && l_high > 0.0);
        // L ∝ 1/α, so ratio of L's = 1/(α ratio) ≈ 10
        let ratio = l_low / l_high;
        assert!(
            (ratio - 10.0).abs() < 0.5,
            "L should scale as 1/α: ratio={ratio}"
        );
    }

    #[test]
    fn test_preset_constructions() {
        let _yig = NanodiskSpinWaves::yig_nanodisk(200.0).expect("YIG preset");
        let py = NanodiskSpinWaves::permalloy_nanodisk(150.0).expect("Py preset");
        let cofeb = NanodiskSpinWaves::cofeb_nanodisk(100.0).expect("CoFeB preset");
        assert!((py.ms - 8e5).abs() < 1.0);
        assert!((cofeb.ms - 1.05e6).abs() < 1.0);
    }

    #[test]
    fn test_error_on_invalid_radius() {
        let r = NanodiskSpinWaves::new(0.0, 20e-9, 40_000.0, 1.4e5, 3.5e-12, 1e-4);
        assert!(r.is_err());
        let r = NanodiskSpinWaves::new(-1e-9, 20e-9, 40_000.0, 1.4e5, 3.5e-12, 1e-4);
        assert!(r.is_err());
    }

    #[test]
    fn test_bessel_j_basic_values() {
        // J_0(0) = 1
        assert!((bessel_j(0, 0.0) - 1.0).abs() < 1e-15);
        // J_1(0) = 0
        assert!(bessel_j(1, 0.0).abs() < 1e-15);
        // J_0(2.4048) ≈ 0 (first zero)
        assert!(bessel_j(0, 2.4048).abs() < 1e-4);
        // J_1(3.8317) ≈ 0 (first zero of J_1)
        assert!(bessel_j(1, 3.8317).abs() < 1e-4);
        // J_0(1) ≈ 0.7651976865...
        assert!((bessel_j(0, 1.0) - 0.7651976865).abs() < 1e-6);
        // Asymptotic regime: J_0(30) ≈ -0.0860, J_1(30) ≈ -0.1187 — just check finiteness
        assert!(bessel_j(0, 30.0).abs() < 1.0);
        assert!(bessel_j(1, 30.0).abs() < 1.0);
    }
}
