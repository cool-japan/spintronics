//! Spin wave dispersion relations
//!
//! This module implements various spin wave dispersion relations used in magnonics
//! and spin wave physics. It covers the full hierarchy from the simple Kittel formula
//! for uniform precession (k=0) through exchange-dominated spin waves to the complete
//! Kalinikos-Slavin dipole-exchange theory for thin ferromagnetic films.
//!
//! # Physical Background
//!
//! Spin waves (magnons) are collective excitations of the magnetic order in a
//! ferromagnet. Their dispersion relation omega(k) determines the frequency as a
//! function of wavevector, and depends on exchange, dipolar, and Zeeman energies.
//!
//! # References
//!
//! - C. Kittel, "On the Theory of Ferromagnetic Resonance Absorption",
//!   Phys. Rev. 73, 155 (1948)
//! - B. A. Kalinikos and A. N. Slavin, "Theory of dipole-exchange spin wave
//!   spectrum for ferromagnetic films with mixed exchange boundary conditions",
//!   J. Phys. C 19, 7013 (1986)

use crate::constants::{GAMMA, MU_0};
use crate::error::{self, Result};
use crate::material::Ferromagnet;

/// Spin wave dispersion calculator for ferromagnetic thin films
///
/// Encapsulates material parameters and film geometry needed to compute
/// spin wave frequencies across the full wavevector range, from dipolar-dominated
/// (long wavelength) to exchange-dominated (short wavelength) regimes.
///
/// # Example
///
/// ```
/// use spintronics::spinwave::SpinWaveDispersion;
/// use spintronics::material::Ferromagnet;
///
/// let yig = Ferromagnet::yig();
/// let disp = SpinWaveDispersion::new(&yig, 1e-6, 0.0)
///     .expect("valid parameters");
///
/// let omega = disp.kittel_frequency(0.1)
///     .expect("valid field");
/// assert!(omega > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct SpinWaveDispersion {
    /// Saturation magnetization \[A/m\]
    pub ms: f64,
    /// Exchange stiffness constant \[J/m\]
    pub exchange_stiffness: f64,
    /// Film thickness \[m\]
    pub thickness: f64,
    /// Gilbert damping parameter (dimensionless)
    pub alpha: f64,
    /// Exchange length D = 2A / (mu_0 * Ms) [m^2]
    exchange_length_sq: f64,
}

impl SpinWaveDispersion {
    /// Create a new spin wave dispersion calculator from material and geometry
    ///
    /// # Arguments
    /// * `material` - Ferromagnetic material parameters
    /// * `thickness` - Film thickness in meters
    /// * `_gamma_override` - Reserved for future use (set to 0.0 to use default GAMMA)
    ///
    /// # Errors
    /// Returns error if material parameters are non-physical (negative or zero Ms, etc.)
    pub fn new(material: &Ferromagnet, thickness: f64, _gamma_override: f64) -> Result<Self> {
        if material.ms <= 0.0 {
            return Err(error::invalid_param(
                "ms",
                "saturation magnetization must be positive",
            ));
        }
        if material.exchange_a < 0.0 {
            return Err(error::invalid_param(
                "exchange_a",
                "exchange stiffness must be non-negative",
            ));
        }
        if thickness <= 0.0 {
            return Err(error::invalid_param(
                "thickness",
                "film thickness must be positive",
            ));
        }
        if material.alpha < 0.0 {
            return Err(error::invalid_param(
                "alpha",
                "Gilbert damping must be non-negative",
            ));
        }

        let exchange_length_sq = 2.0 * material.exchange_a / (MU_0 * material.ms);

        Ok(Self {
            ms: material.ms,
            exchange_stiffness: material.exchange_a,
            thickness,
            alpha: material.alpha,
            exchange_length_sq,
        })
    }

    /// Kittel formula for uniform ferromagnetic resonance (k=0)
    ///
    /// omega = gamma * mu_0 * sqrt(H * (H + Ms))
    ///
    /// This is the fundamental resonance frequency for a thin film magnetized
    /// in-plane by an external field H (in Tesla, interpreted as mu_0*H_ext).
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field magnitude \[T\] (i.e., mu_0 * H in SI)
    ///
    /// # Returns
    /// Angular frequency omega \[rad/s\]
    ///
    /// # Errors
    /// Returns error if the field is negative
    pub fn kittel_frequency(&self, h_ext: f64) -> Result<f64> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }
        // Convert to internal field: H = h_ext / mu_0 in A/m,
        // then omega = gamma * sqrt( mu_0*H * (mu_0*H + mu_0*Ms) )
        // = gamma * sqrt( h_ext * (h_ext + mu_0*Ms) )
        let mu0_ms = MU_0 * self.ms;
        let omega = GAMMA * (h_ext * (h_ext + mu0_ms)).sqrt();
        Ok(omega)
    }

    /// Exchange spin wave dispersion
    ///
    /// omega(k) = gamma * (mu_0 * H_ext + D * k^2)
    ///
    /// where D = 2A / (mu_0 * Ms) is the exchange length squared.
    /// Valid in the short-wavelength (large k) regime where exchange dominates
    /// over dipolar interactions.
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `k` - Wavevector magnitude \[rad/m\]
    ///
    /// # Returns
    /// Angular frequency omega \[rad/s\]
    ///
    /// # Errors
    /// Returns error if field is negative
    pub fn exchange_dispersion(&self, h_ext: f64, k: f64) -> Result<f64> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }
        let omega = GAMMA * (h_ext + MU_0 * self.exchange_length_sq * k * k);
        Ok(omega)
    }

    /// Kalinikos-Slavin dipole-exchange dispersion for thin ferromagnetic films
    ///
    /// This is the full dispersion relation accounting for both dipolar and exchange
    /// interactions in a thin film geometry. The result depends on the angle phi
    /// between the wavevector k and the in-plane magnetization direction.
    ///
    /// omega^2 = (omega_H + omega_M * lambda * k^2) *
    ///           (omega_H + omega_M * lambda * k^2 + omega_M * F_nn)
    ///
    /// where F_nn is the dipolar matrix element for the n=0 (uniform) mode:
    ///   F_00 = 1 - (1 - exp(-k*d))/(k*d) + sin^2(phi) * (1 - exp(-k*d))/(k*d)
    /// which simplifies for specific angles to known limiting cases.
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `k` - In-plane wavevector magnitude \[rad/m\]
    /// * `phi` - Angle between k and in-plane magnetization \[rad\]
    ///
    /// # Returns
    /// Angular frequency omega \[rad/s\]
    ///
    /// # Errors
    /// Returns error for invalid parameters or if the argument under the square root is negative
    ///
    /// # References
    /// B. A. Kalinikos and A. N. Slavin, J. Phys. C 19, 7013 (1986)
    pub fn kalinikos_slavin(&self, h_ext: f64, k: f64, phi: f64) -> Result<f64> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }

        let omega_h = GAMMA * h_ext;
        let omega_m = GAMMA * MU_0 * self.ms;
        let lambda = self.exchange_length_sq;
        let d = self.thickness;
        let kd = k.abs() * d;

        // Dipolar matrix element F_nn for n=0 mode
        let p = if kd < 1e-12 {
            // Taylor expansion: (1 - exp(-kd)) / kd -> 1 - kd/2 + ...
            1.0 - kd / 2.0
        } else {
            (1.0 - (-kd).exp()) / kd
        };

        let sin2_phi = phi.sin().powi(2);
        let f_nn = 1.0 - p + sin2_phi * p;

        let term1 = omega_h + omega_m * lambda * k * k;
        let term2 = omega_h + omega_m * lambda * k * k + omega_m * f_nn;

        let omega_sq = term1 * term2;
        if omega_sq < 0.0 {
            return Err(error::numerical_error(
                "negative argument under square root in Kalinikos-Slavin dispersion",
            ));
        }

        Ok(omega_sq.sqrt())
    }

    /// Dipolar (magnetostatic) spin wave dispersion in the long-wavelength limit
    ///
    /// In the limit where exchange is negligible (k -> 0, but k*d finite),
    /// the dispersion is purely determined by dipolar interactions:
    ///
    /// omega^2 = omega_H * (omega_H + omega_M) + omega_M^2/4 * (1 - exp(-2*k*d))
    ///
    /// This expression interpolates between the Damon-Eshbach surface mode
    /// and the volume mode depending on k*d.
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `k` - Wavevector magnitude \[rad/m\]
    ///
    /// # Returns
    /// Angular frequency omega \[rad/s\]
    pub fn dipolar_dispersion(&self, h_ext: f64, k: f64) -> Result<f64> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }

        let omega_h = GAMMA * h_ext;
        let omega_m = GAMMA * MU_0 * self.ms;
        let kd = k.abs() * self.thickness;

        let omega_sq =
            omega_h * (omega_h + omega_m) + omega_m * omega_m / 4.0 * (1.0 - (-2.0 * kd).exp());

        if omega_sq < 0.0 {
            return Err(error::numerical_error(
                "negative argument in dipolar dispersion",
            ));
        }

        Ok(omega_sq.sqrt())
    }

    /// Group velocity v_g = d(omega)/dk computed via numerical central difference
    ///
    /// Uses a small finite difference step to approximate the derivative of the
    /// dispersion relation at the given wavevector.
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `k` - Wavevector \[rad/m\]
    /// * `phi` - Angle between k and magnetization \[rad\]
    ///
    /// # Returns
    /// Group velocity \[m/s\]
    pub fn group_velocity(&self, h_ext: f64, k: f64, phi: f64) -> Result<f64> {
        // Use central difference with adaptive step size
        let dk = if k.abs() > 1e-6 {
            k.abs() * 1e-6
        } else {
            1.0 // 1 rad/m for very small k
        };

        let omega_plus = self.kalinikos_slavin(h_ext, k + dk, phi)?;
        let omega_minus = self.kalinikos_slavin(h_ext, (k - dk).max(0.0), phi)?;

        let effective_dk = if k - dk < 0.0 {
            k + dk // forward difference if k-dk would be negative
        } else {
            2.0 * dk
        };

        Ok((omega_plus - omega_minus) / effective_dk)
    }

    /// Spin wave lifetime from Gilbert damping
    ///
    /// tau = 1 / (alpha * omega)
    ///
    /// The Gilbert damping causes exponential decay of spin wave amplitude
    /// with this characteristic time scale.
    ///
    /// # Arguments
    /// * `omega` - Angular frequency of the spin wave \[rad/s\]
    ///
    /// # Returns
    /// Lifetime tau \[s\]
    ///
    /// # Errors
    /// Returns error if omega is non-positive or alpha is zero
    pub fn spin_wave_lifetime(&self, omega: f64) -> Result<f64> {
        if omega <= 0.0 {
            return Err(error::invalid_param(
                "omega",
                "angular frequency must be positive",
            ));
        }
        if self.alpha <= 0.0 {
            return Err(error::invalid_param(
                "alpha",
                "Gilbert damping must be positive for finite lifetime",
            ));
        }

        Ok(1.0 / (self.alpha * omega))
    }

    /// Spin wave propagation length
    ///
    /// l_prop = v_g * tau = v_g / (alpha * omega)
    ///
    /// The characteristic distance a spin wave can travel before its amplitude
    /// decays by a factor of 1/e.
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `k` - Wavevector \[rad/m\]
    /// * `phi` - Angle between k and magnetization \[rad\]
    ///
    /// # Returns
    /// Propagation length \[m\]
    pub fn propagation_length(&self, h_ext: f64, k: f64, phi: f64) -> Result<f64> {
        let vg = self.group_velocity(h_ext, k, phi)?;
        let omega = self.kalinikos_slavin(h_ext, k, phi)?;
        let tau = self.spin_wave_lifetime(omega)?;
        Ok(vg.abs() * tau)
    }

    /// Exchange length squared D = 2A / (mu_0 * Ms) [m^2]
    ///
    /// This is the characteristic length scale at which exchange and dipolar
    /// energies are comparable.
    pub fn exchange_length_squared(&self) -> f64 {
        self.exchange_length_sq
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    fn yig_dispersion() -> SpinWaveDispersion {
        let yig = Ferromagnet::yig();
        SpinWaveDispersion::new(&yig, 1e-6, 0.0).expect("valid YIG parameters")
    }

    fn permalloy_dispersion() -> SpinWaveDispersion {
        let py = Ferromagnet::permalloy();
        SpinWaveDispersion::new(&py, 50e-9, 0.0).expect("valid Permalloy parameters")
    }

    #[test]
    fn test_kittel_yig() {
        // YIG: Ms = 1.4e5 A/m, mu_0*Ms ~ 0.176 T
        // At H_ext = 0.1 T:
        // omega = gamma * sqrt(0.1 * (0.1 + 0.176)) ~ gamma * sqrt(0.0276)
        let disp = yig_dispersion();
        let omega = disp.kittel_frequency(0.1).expect("valid field");

        // Expected: ~2.92e10 rad/s (about 4.65 GHz)
        assert!(omega > 2.0e10, "omega = {omega}");
        assert!(omega < 4.0e10, "omega = {omega}");
    }

    #[test]
    fn test_kittel_permalloy() {
        // Permalloy: Ms = 8.0e5 A/m, mu_0*Ms ~ 1.005 T
        // At H_ext = 0.05 T:
        // omega = gamma * sqrt(0.05 * (0.05 + 1.005)) ~ gamma * sqrt(0.05275)
        let disp = permalloy_dispersion();
        let omega = disp.kittel_frequency(0.05).expect("valid field");

        // Expected: ~4.04e10 rad/s (about 6.4 GHz)
        assert!(omega > 3.0e10, "omega = {omega}");
        assert!(omega < 5.0e10, "omega = {omega}");
    }

    #[test]
    fn test_kittel_zero_field() {
        let disp = yig_dispersion();
        let omega = disp.kittel_frequency(0.0).expect("zero field allowed");
        assert_eq!(omega, 0.0);
    }

    #[test]
    fn test_exchange_quadratic_dispersion() {
        let disp = yig_dispersion();
        let h_ext = 0.1;
        let k1 = 1e7;
        let k2 = 2e7;

        let omega1 = disp.exchange_dispersion(h_ext, k1).expect("valid k1");
        let omega2 = disp.exchange_dispersion(h_ext, k2).expect("valid k2");

        // omega(k) - omega(0) should scale as k^2
        let omega0 = disp.exchange_dispersion(h_ext, 0.0).expect("valid k=0");
        let delta1 = omega1 - omega0;
        let delta2 = omega2 - omega0;

        // delta2/delta1 should be close to (k2/k1)^2 = 4
        let ratio = delta2 / delta1;
        assert!(
            (ratio - 4.0).abs() < 0.01,
            "ratio should be ~4.0, got {ratio}"
        );
    }

    #[test]
    fn test_kalinikos_slavin_phi_0() {
        // phi=0 corresponds to backward volume geometry (k parallel to M)
        let disp = yig_dispersion();
        let omega = disp
            .kalinikos_slavin(0.1, 1e6, 0.0)
            .expect("valid parameters");
        assert!(omega > 0.0, "omega must be positive");
    }

    #[test]
    fn test_kalinikos_slavin_phi_pi2() {
        // phi=pi/2 corresponds to Damon-Eshbach geometry (k perp to M)
        let disp = yig_dispersion();
        let omega = disp
            .kalinikos_slavin(0.1, 1e6, PI / 2.0)
            .expect("valid parameters");
        assert!(omega > 0.0, "omega must be positive");
    }

    #[test]
    fn test_kalinikos_slavin_reduces_to_kittel_at_k0() {
        let disp = yig_dispersion();
        let h_ext = 0.1;

        let omega_kittel = disp.kittel_frequency(h_ext).expect("valid field");
        // At k=0 with phi=pi/2 (Damon-Eshbach geometry), F_nn -> 1,
        // so KS gives sqrt(omega_H * (omega_H + omega_M)) = Kittel frequency
        // Use very small k
        let omega_ks = disp
            .kalinikos_slavin(h_ext, 1e-3, std::f64::consts::FRAC_PI_2)
            .expect("valid small k");

        let rel_diff = (omega_ks - omega_kittel).abs() / omega_kittel.max(1.0);
        assert!(
            rel_diff < 0.01,
            "KS at k~0 should match Kittel: omega_ks={omega_ks}, omega_kittel={omega_kittel}"
        );
    }

    #[test]
    fn test_dipolar_dispersion() {
        let disp = yig_dispersion();
        let omega = disp.dipolar_dispersion(0.1, 1e5).expect("valid parameters");
        assert!(omega > 0.0);
    }

    #[test]
    fn test_group_velocity_positive_de() {
        // Damon-Eshbach modes have positive group velocity
        let disp = yig_dispersion();
        let vg = disp
            .group_velocity(0.1, 1e6, PI / 2.0)
            .expect("valid parameters");
        assert!(
            vg > 0.0,
            "DE mode should have positive group velocity: vg={vg}"
        );
    }

    #[test]
    fn test_spin_wave_lifetime() {
        let disp = yig_dispersion();
        let omega = 2.0 * PI * 5e9; // 5 GHz
        let tau = disp.spin_wave_lifetime(omega).expect("valid omega");

        // YIG alpha = 0.0001, tau = 1/(0.0001 * 2*pi*5e9) ~ 318 ns
        assert!(tau > 100e-9, "YIG lifetime should be > 100 ns: tau={tau}");
        assert!(tau < 1e-6, "YIG lifetime should be < 1 us: tau={tau}");
    }

    #[test]
    fn test_lifetime_scaling_with_damping() {
        let yig = Ferromagnet::yig();
        let omega = 2.0 * PI * 5e9;

        let disp1 = SpinWaveDispersion::new(&yig, 1e-6, 0.0).expect("valid");
        let tau1 = disp1.spin_wave_lifetime(omega).expect("valid omega");

        let py = Ferromagnet::permalloy();
        let disp2 = SpinWaveDispersion::new(&py, 1e-6, 0.0).expect("valid");
        let tau2 = disp2.spin_wave_lifetime(omega).expect("valid omega");

        // tau scales as 1/alpha, so tau_yig/tau_py ~ alpha_py/alpha_yig = 100
        let ratio = tau1 / tau2;
        assert!(
            (ratio - 100.0).abs() < 1.0,
            "lifetime ratio should be ~100: got {ratio}"
        );
    }

    #[test]
    fn test_invalid_negative_field() {
        let disp = yig_dispersion();
        assert!(disp.kittel_frequency(-0.1).is_err());
        assert!(disp.exchange_dispersion(-0.1, 1e6).is_err());
        assert!(disp.kalinikos_slavin(-0.1, 1e6, 0.0).is_err());
    }

    #[test]
    fn test_invalid_parameters() {
        let mut bad_material = Ferromagnet::yig();
        bad_material.ms = -1.0;
        assert!(SpinWaveDispersion::new(&bad_material, 1e-6, 0.0).is_err());

        assert!(SpinWaveDispersion::new(&Ferromagnet::yig(), -1.0, 0.0).is_err());
        assert!(SpinWaveDispersion::new(&Ferromagnet::yig(), 0.0, 0.0).is_err());
    }

    #[test]
    fn test_propagation_length() {
        let disp = yig_dispersion();
        let l = disp
            .propagation_length(0.1, 1e6, PI / 2.0)
            .expect("valid parameters");
        // YIG propagation length can be mm-scale
        assert!(l > 1e-6, "propagation length should be > 1 um: l={l}");
    }
}
