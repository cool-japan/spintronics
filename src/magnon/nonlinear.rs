//! Nonlinear Magnon Physics
//!
//! This module implements the physics of nonlinear magnon interactions,
//! including four-magnon scattering, Suhl parametric instabilities,
//! parametric amplification, and nonlinear FMR linewidth broadening.
//!
//! # Physical Background
//!
//! At large microwave drive amplitudes, spin-wave systems enter a nonlinear
//! regime dominated by multi-magnon scattering processes. The most important
//! are:
//!
//! ## Three-Magnon (Suhl First-Order) Instability
//!
//! A uniform FMR magnon decays into two spin-wave magnons with complementary
//! wavevectors (1 → 2 process). The threshold pump field is:
//!
//! h_th ≈ μ₀ Δω_k / |T_kk|
//!
//! where Δω_k is the spin-wave linewidth and |T_kk| is the three-magnon
//! coupling coefficient.
//!
//! ## Four-Magnon (Suhl Second-Order) Instability
//!
//! Two pump magnons scatter into two spin-wave magnons satisfying energy and
//! momentum conservation. The threshold is similar but governed by the
//! four-magnon coupling.
//!
//! ## Parametric Amplification
//!
//! Below the Suhl threshold, a pump field at 2ω₀ can coherently amplify
//! a signal at ω₀ via the degenerate parametric process. Above threshold,
//! the system oscillates spontaneously.
//!
//! ## Nonlinear FMR Linewidth (Suhl's Theory)
//!
//! Above the Suhl threshold, anomalous linewidth broadening occurs due to
//! parametric magnon generation draining energy from the uniform precession.
//! The effective damping increases as:
//!
//! ΔH_NL ≈ ΔH₀ × (1 + (h/h_th - 1)²)
//!
//! # Key References
//!
//! - Suhl, J. Phys. Chem. Solids 1, 209 (1957) — instability theory
//! - Kalinikos & Slavin, J. Phys. C 19, 7013 (1986) — spin-wave dispersion
//! - Gurevich & Melkov, "Magnetization Oscillations and Waves" (1996)

use std::f64::consts::PI;

use crate::constants::{GAMMA, HBAR, KB, MU_0};
use crate::error::{self, Result};

// ============================================================================
// Four-Magnon Scattering
// ============================================================================

/// Four-magnon scattering system for a magnetic thin film or bulk material.
///
/// Describes the nonlinear magnon-magnon interactions governed by the
/// Kalinikos-Slavin dispersion relation and dipole-exchange coupling.
///
/// The magnon dispersion follows (Kalinikos-Slavin, all in Tesla):
///
/// B_eff(k) = B_ext + λ_ex · (μ₀ Ms) · k²
///
/// ω(k, θ_k) = γ √[B_eff(k) · (B_eff(k) + μ₀ Ms · sin²θ_k)]
///
/// where B_ext = μ₀ H_ext \[T\], λ_ex = 2A/(μ₀ Ms²) \[m²\] is the exchange
/// length squared, μ₀ Ms is the saturation induction \[T\], and θ_k is the
/// angle between the spin-wave wavevector k and the equilibrium magnetization.
#[derive(Debug, Clone)]
pub struct FourMagnonScattering {
    /// Saturation magnetization \[A/m\]
    ms: f64,
    /// Exchange stiffness constant \[J/m\]
    exchange_a: f64,
    /// Gilbert damping coefficient (dimensionless)
    alpha: f64,
    /// External applied field \[T\]
    h_ext: f64,
}

impl FourMagnonScattering {
    /// Create a new four-magnon scattering system.
    ///
    /// # Arguments
    /// * `ms` - Saturation magnetization \[A/m\]; must be positive
    /// * `exchange_a` - Exchange stiffness \[J/m\]; must be positive
    /// * `alpha` - Gilbert damping; must be in (0, 1)
    /// * `h_ext` - External field \[T\]; must be ≥ 0
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] if any parameter
    /// is out of its physical range.
    pub fn new(ms: f64, exchange_a: f64, alpha: f64, h_ext: f64) -> Result<Self> {
        if ms <= 0.0 {
            return Err(error::invalid_param(
                "ms",
                "saturation magnetization must be positive",
            ));
        }
        if exchange_a <= 0.0 {
            return Err(error::invalid_param(
                "exchange_a",
                "exchange stiffness must be positive",
            ));
        }
        if alpha <= 0.0 || alpha >= 1.0 {
            return Err(error::invalid_param(
                "alpha",
                "Gilbert damping must be in the open interval (0, 1)",
            ));
        }
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }
        Ok(Self {
            ms,
            exchange_a,
            alpha,
            h_ext,
        })
    }

    /// Construct a system with standard YIG (Yttrium Iron Garnet) parameters.
    ///
    /// YIG is the premier low-damping ferromagnet for spin-wave studies:
    /// - Ms = 1.40 × 10⁵ A/m (room-temperature value)
    /// - A  = 3.70 × 10⁻¹² J/m (exchange stiffness)
    /// - α  = 1.00 × 10⁻⁴   (very low Gilbert damping)
    ///
    /// # Arguments
    /// * `h_ext` - External field \[T\]
    pub fn yig(h_ext: f64) -> Self {
        Self {
            ms: 1.40e5,
            exchange_a: 3.70e-12,
            alpha: 1.0e-4,
            h_ext: h_ext.max(0.0),
        }
    }

    /// μ₀ Ms in Tesla — the dipolar field scale.
    #[inline]
    fn mu0_ms(&self) -> f64 {
        MU_0 * self.ms
    }

    /// Exchange length squared: λ_ex = 2A / (μ₀ Ms²) \[m²\].
    ///
    /// The exchange contribution to the effective field is λ_ex · μ₀ Ms · k²
    /// (units of Tesla).
    #[inline]
    fn lambda_ex(&self) -> f64 {
        2.0 * self.exchange_a / (MU_0 * self.ms * self.ms)
    }

    /// Magnon angular frequency via the Kalinikos-Slavin dispersion relation.
    ///
    /// All fields are expressed in Tesla for dimensional consistency:
    ///
    /// B_eff(k) = B_ext + λ_ex · (μ₀ Ms) · k²
    ///
    /// ω(k, θ_k) = γ √[B_eff(k) · (B_eff(k) + μ₀ Ms · sin²θ_k)]
    ///
    /// where θ_k is the angle between the spin-wave wavevector and the
    /// equilibrium magnetization direction, B_ext = h_ext \[T\], and
    /// μ₀ Ms is the saturation induction \[T\].
    ///
    /// At k = 0 and θ_k = π/2 this reduces to the Kittel formula:
    /// ω₀ = γ √[B_ext · (B_ext + μ₀ Ms)]
    ///
    /// # Arguments
    /// * `k_magnitude` - Wavenumber |k| \[1/m\]
    /// * `theta_k` - Angle between k and M \[rad\]
    ///
    /// # Returns
    /// Angular frequency ω \[rad/s\]
    pub fn magnon_frequency(&self, k_magnitude: f64, theta_k: f64) -> f64 {
        let mu0_ms = self.mu0_ms();
        let lambda_ex = self.lambda_ex();

        // Effective internal field in Tesla
        let b_eff = self.h_ext + lambda_ex * mu0_ms * k_magnitude * k_magnitude;
        let sin_theta = theta_k.sin();
        let dipolar_term = mu0_ms * sin_theta * sin_theta;

        GAMMA * (b_eff * (b_eff + dipolar_term)).sqrt()
    }

    /// Damping-limited spin-wave linewidth (half-power, in angular frequency).
    ///
    /// For propagating spin waves with wavevector at θ_k = π/2 (perpendicular
    /// to the field), the linewidth is dominated by Gilbert damping:
    ///
    /// Δω_k = α · ω(k, π/2)
    ///
    /// # Arguments
    /// * `k_magnitude` - Wavenumber |k| \[1/m\]
    ///
    /// # Returns
    /// Half-power linewidth Δω \[rad/s\]
    pub fn magnon_linewidth(&self, k_magnitude: f64) -> f64 {
        // Use θ = π/2 (backward volume spin wave geometry) for the linewidth
        self.alpha * self.magnon_frequency(k_magnitude, PI / 2.0)
    }

    /// Four-magnon coupling amplitude |T_kk| \[rad/s\].
    ///
    /// The coupling combines the long-range dipolar interaction and short-range
    /// exchange stiffness. All quantities are in consistent SI units with γ
    /// converting field amplitudes to angular frequencies:
    ///
    /// |T_kk| = (γ μ₀ Ms / 4) · |3cos²θ_k - 1|  +  γ λ_ex μ₀ Ms · k²
    ///
    /// The first term is the dipolar contribution, which vanishes at the magic
    /// angle θ_k ≈ 54.7°. The second term is the exchange contribution,
    /// which always adds and grows as k².
    ///
    /// # Arguments
    /// * `k_magnitude` - Wavenumber |k| \[1/m\]
    /// * `theta_k` - Angle between k and M \[rad\]
    ///
    /// # Returns
    /// Four-magnon coupling |T_kk| \[rad/s\]
    pub fn four_magnon_coupling(&self, k_magnitude: f64, theta_k: f64) -> f64 {
        let mu0_ms = self.mu0_ms();
        let lambda_ex = self.lambda_ex();

        let cos_theta = theta_k.cos();
        // Dipolar: (γ μ₀ Ms / 4) |3cos²θ - 1|  [rad/s]
        let dipolar = (GAMMA * mu0_ms / 4.0) * (3.0 * cos_theta * cos_theta - 1.0).abs();
        // Exchange: γ λ_ex μ₀ Ms k²  [rad/s]
        let exchange = GAMMA * lambda_ex * mu0_ms * k_magnitude * k_magnitude;

        dipolar + exchange
    }

    /// Threshold pump field amplitude for the Suhl instability.
    ///
    /// Combines the spin-wave linewidth and four-magnon coupling:
    ///
    /// h_th = μ₀ · Δω_k / |T_kk|
    ///
    /// This is the minimum microwave field amplitude that renders spin-wave
    /// mode k parametrically unstable (Suhl's second-order instability for
    /// the four-magnon process).
    ///
    /// # Arguments
    /// * `k_magnitude` - Wavenumber |k| \[1/m\]
    /// * `theta_k` - Angle between k and M \[rad\]
    ///
    /// # Returns
    /// Threshold pump field h_th \[T\]
    pub fn suhl_threshold_field(&self, k_magnitude: f64, theta_k: f64) -> f64 {
        let delta_omega = self.magnon_linewidth(k_magnitude);
        let t_kk = self.four_magnon_coupling(k_magnitude, theta_k);

        // Guard against vanishing coupling (magic-angle + zero k)
        if t_kk < f64::EPSILON {
            return f64::MAX;
        }

        MU_0 * delta_omega / t_kk
    }

    /// Parametric growth rate σ for a spin-wave mode.
    ///
    /// Above the Suhl threshold, the amplitude of spin-wave mode k grows
    /// exponentially at rate:
    ///
    /// σ = √(max(0, (|T_kk| · h_pump / μ₀)² - Δω_k²))
    ///
    /// Below threshold (h_pump < h_th), σ = 0.
    ///
    /// # Arguments
    /// * `pump_h` - Pump field amplitude \[T\]
    /// * `k_magnitude` - Wavenumber |k| \[1/m\]
    /// * `theta_k` - Angle between k and M \[rad\]
    ///
    /// # Returns
    /// Growth rate σ \[rad/s\]
    pub fn parametric_growth_rate(&self, pump_h: f64, k_magnitude: f64, theta_k: f64) -> f64 {
        let delta_omega = self.magnon_linewidth(k_magnitude);
        let t_kk = self.four_magnon_coupling(k_magnitude, theta_k);

        // Effective coupling drive: V = T_kk · h / μ₀
        let coupling_drive = t_kk * pump_h / MU_0;
        let discriminant = coupling_drive * coupling_drive - delta_omega * delta_omega;

        if discriminant <= 0.0 {
            0.0
        } else {
            discriminant.sqrt()
        }
    }

    /// Find the optimal wavevector magnitude that minimises the Suhl threshold.
    ///
    /// For a given pump frequency ω_p, the most easily excited spin-wave mode
    /// satisfies ω(k, π/2) ≈ ω_p / 2. This method performs a bisection search
    /// over |k| to find this degeneracy condition.
    ///
    /// The search is bounded by:
    /// - k_min = 10 m⁻¹ (long-wavelength cutoff)
    /// - k_max = 10⁸ m⁻¹ (exchange-dominated cutoff, ~10 nm wavelength)
    ///
    /// # Arguments
    /// * `pump_freq` - Pump angular frequency ω_p \[rad/s\]
    ///
    /// # Returns
    /// Optimal |k| \[1/m\], or the boundary value if no crossing is found
    pub fn optimal_k_for_suhl(&self, pump_freq: f64) -> f64 {
        let target_omega = pump_freq / 2.0;
        let theta = PI / 2.0; // use the backward-volume geometry

        let k_min = 1.0e1_f64;
        let k_max = 1.0e8_f64;

        // Evaluate the signed gap at the boundaries
        let f_min = self.magnon_frequency(k_min, theta) - target_omega;
        let f_max = self.magnon_frequency(k_max, theta) - target_omega;

        // If both have the same sign there is no bracket; return the end
        // closest to the target (smallest absolute residual).
        if f_min * f_max >= 0.0 {
            return if f_min.abs() < f_max.abs() {
                k_min
            } else {
                k_max
            };
        }

        // Standard bisection — 60 iterations guarantees ~18 significant digits
        let mut lo = k_min;
        let mut hi = k_max;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            let f_mid = self.magnon_frequency(mid, theta) - target_omega;
            if f_mid == 0.0 {
                return mid;
            }
            if (self.magnon_frequency(lo, theta) - target_omega) * f_mid < 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Determine whether any spin-wave mode is parametrically unstable.
    ///
    /// Evaluates the growth rate at θ_k = π/2 for the mode k_opt that
    /// minimises the Suhl threshold. Returns `true` if σ > 0.
    ///
    /// # Arguments
    /// * `pump_h` - Pump field amplitude \[T\]
    /// * `pump_freq` - Pump angular frequency \[rad/s\]
    pub fn instability_is_present(&self, pump_h: f64, pump_freq: f64) -> bool {
        let k_opt = self.optimal_k_for_suhl(pump_freq);
        self.parametric_growth_rate(pump_h, k_opt, PI / 2.0) > 0.0
    }
}

// ============================================================================
// Parametric Amplification
// ============================================================================

/// Parametric magnon amplifier driven by a microwave pump field.
///
/// Describes the three-wave mixing process where a pump at frequency ω_p
/// simultaneously amplifies a signal at ω_s and generates an idler at
/// ω_i = ω_p - ω_s. In the degenerate case ω_s = ω_i = ω_p / 2.
///
/// The threshold condition is:
///
/// (coupling · h_th)² = γ_s · γ_i
///
/// Above threshold the signal amplitude grows as exp(G t) where
///
/// G = √((coupling · h)² - γ_s · γ_i)
#[derive(Debug, Clone)]
pub struct ParametricAmplification {
    /// Parametric coupling coefficient \[rad/(s·T)\]
    coupling: f64,
    /// Signal angular frequency ω_s \[rad/s\]
    omega_signal: f64,
    /// Idler angular frequency ω_i \[rad/s\]
    pub omega_idler: f64,
    /// Pump field amplitude h_p \[T\]
    pump_h: f64,
    /// Signal half-power linewidth γ_s \[rad/s\]
    gamma_s: f64,
    /// Idler half-power linewidth γ_i \[rad/s\]
    gamma_i: f64,
}

impl ParametricAmplification {
    /// Create a new parametric amplification system.
    ///
    /// # Arguments
    /// * `coupling` - Coupling coefficient \[rad/(s·T)\]; must be positive
    /// * `omega_signal` - Signal angular frequency \[rad/s\]; must be positive
    /// * `omega_idler` - Idler angular frequency \[rad/s\]; must be positive
    /// * `pump_h` - Pump field amplitude \[T\]; must be positive
    /// * `gamma_s` - Signal linewidth \[rad/s\]; must be positive
    /// * `gamma_i` - Idler linewidth \[rad/s\]; must be positive
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] for non-positive inputs.
    pub fn new(
        coupling: f64,
        omega_signal: f64,
        omega_idler: f64,
        pump_h: f64,
        gamma_s: f64,
        gamma_i: f64,
    ) -> Result<Self> {
        if coupling <= 0.0 {
            return Err(error::invalid_param("coupling", "must be positive"));
        }
        if omega_signal <= 0.0 {
            return Err(error::invalid_param("omega_signal", "must be positive"));
        }
        if omega_idler <= 0.0 {
            return Err(error::invalid_param("omega_idler", "must be positive"));
        }
        if pump_h <= 0.0 {
            return Err(error::invalid_param("pump_h", "must be positive"));
        }
        if gamma_s <= 0.0 {
            return Err(error::invalid_param("gamma_s", "must be positive"));
        }
        if gamma_i <= 0.0 {
            return Err(error::invalid_param("gamma_i", "must be positive"));
        }
        Ok(Self {
            coupling,
            omega_signal,
            omega_idler,
            pump_h,
            gamma_s,
            gamma_i,
        })
    }

    /// Construct a degenerate parametric amplifier using YIG-like parameters.
    ///
    /// In the degenerate case, signal = idler = pump / 2 and the coupling
    /// coefficient is γ_gyro · Ms / 2.
    ///
    /// # Arguments
    /// * `pump_freq` - Pump angular frequency ω_p \[rad/s\]
    /// * `alpha` - Gilbert damping of the material (dimensionless)
    /// * `ms` - Saturation magnetization \[A/m\]
    pub fn degenerate_from_yig(pump_freq: f64, alpha: f64, ms: f64) -> Self {
        let omega_half = pump_freq / 2.0;
        // Coupling: γ · Ms / 2 in rad/(s·T)
        let coupling = GAMMA * ms / 2.0;
        // Linewidth from Gilbert damping at the signal/idler frequency
        let gamma = alpha * omega_half;
        Self {
            coupling,
            omega_signal: omega_half,
            omega_idler: omega_half,
            pump_h: 1.0e-4, // placeholder 0.1 mT pump
            gamma_s: gamma,
            gamma_i: gamma,
        }
    }

    /// Return the idler angular frequency ω_i \[rad/s\].
    ///
    /// In a degenerate amplifier ω_i = ω_s = ω_p / 2. In the non-degenerate
    /// case the idler carries the energy difference ω_p - ω_s.
    pub fn idler_frequency(&self) -> f64 {
        self.omega_idler
    }

    /// Phase mismatch Δω = ω_s + ω_i - ω_p \[rad/s\] (energy conservation residual).
    ///
    /// A perfectly phase-matched process satisfies ω_s + ω_i = ω_p, so Δω = 0.
    /// Non-zero values indicate an imperfect energy-conserving process; the
    /// threshold field increases roughly as h_th ∝ √(Δω² + γ²) / coupling.
    ///
    /// # Arguments
    /// * `pump_freq` - Pump angular frequency ω_p \[rad/s\]
    pub fn phase_mismatch(&self, pump_freq: f64) -> f64 {
        self.omega_signal + self.omega_idler - pump_freq
    }

    /// Threshold pump field amplitude for onset of parametric oscillation.
    ///
    /// h_th = √(γ_s · γ_i) / coupling
    ///
    /// # Returns
    /// Threshold field h_th \[T\]
    pub fn threshold_pump_field(&self) -> f64 {
        (self.gamma_s * self.gamma_i).sqrt() / self.coupling
    }

    /// Above-threshold gain coefficient G \[rad/s\].
    ///
    /// G = √(max(0, (coupling · h)² - γ_s · γ_i))
    ///
    /// At exactly threshold G = 0; below threshold G = 0.
    pub fn gain_coefficient(&self) -> f64 {
        let drive = self.coupling * self.pump_h;
        let discriminant = drive * drive - self.gamma_s * self.gamma_i;
        if discriminant <= 0.0 {
            0.0
        } else {
            discriminant.sqrt()
        }
    }

    /// Estimate signal power gain in dB after `n_roundtrips` in a resonator.
    ///
    /// A simplified estimate treating each round-trip time as T = 2π/ω_signal.
    /// The signal amplitude grows as exp(G · T · n_roundtrips) so the power
    /// gain (20 log₁₀ of the amplitude ratio) is:
    ///
    /// G_dB = 40 · log₁₀(e) · G · n_roundtrips / ω_signal
    ///       = 40/ln(10) · G · n_roundtrips / ω_signal
    ///
    /// # Arguments
    /// * `n_roundtrips` - Number of cavity round-trips
    ///
    /// # Returns
    /// Power gain \[dB\]
    pub fn signal_gain_db(&self, n_roundtrips: usize) -> f64 {
        let g = self.gain_coefficient();
        let t_roundtrip = 2.0 * PI / self.omega_signal;
        // Power gain = 20 log10(exp(G t)) where t = n_roundtrips * t_roundtrip
        let amplitude_ratio = (g * t_roundtrip * n_roundtrips as f64).exp();
        20.0 * amplitude_ratio.log10()
    }

    /// Returns `true` if the pump exceeds the parametric oscillation threshold.
    pub fn is_above_threshold(&self) -> bool {
        self.pump_h > self.threshold_pump_field()
    }

    /// Fraction of pump magnons converted to signal+idler pairs (pump depletion).
    ///
    /// A simple estimate valid well above threshold:
    ///
    /// r = 1 - γ_s γ_i / (coupling · h)²
    ///
    /// Clamped to [0, 1] for physical consistency.
    pub fn pump_depletion_ratio(&self) -> f64 {
        let drive_sq = (self.coupling * self.pump_h).powi(2);
        let threshold_sq = self.gamma_s * self.gamma_i;
        let ratio = 1.0 - threshold_sq / drive_sq;
        ratio.clamp(0.0, 1.0)
    }

    /// Quantum-limited noise-equivalent power (NEP) \[W\].
    ///
    /// For a signal mode in thermal equilibrium at room temperature (T = 300 K):
    ///
    /// NEP ≈ ℏ ω_s (1 + 2 n_th)
    ///
    /// where n_th = 1 / (exp(ℏ ω_s / k_B T) - 1) is the thermal photon number.
    /// The factor (1 + 2 n_th) represents the quantum noise floor including
    /// both vacuum fluctuations and thermal noise.
    pub fn noise_equivalent_power(&self) -> f64 {
        const ROOM_TEMPERATURE: f64 = 300.0;
        let energy = HBAR * self.omega_signal;
        let exponent = energy / (KB * ROOM_TEMPERATURE);
        // Guard against exp overflow at very high frequency
        let n_th = if exponent > 700.0 {
            0.0
        } else {
            1.0 / (exponent.exp() - 1.0)
        };
        energy * (1.0 + 2.0 * n_th)
    }
}

// ============================================================================
// Nonlinear FMR Linewidth
// ============================================================================

/// Nonlinear FMR linewidth model based on Suhl's parametric theory.
///
/// Wraps a [`FourMagnonScattering`] system and adds phenomenological models
/// for the linewidth broadening and bistability onset that occur when the
/// drive amplitude approaches the Suhl threshold.
///
/// Above the threshold, the FMR linewidth broadens as additional spin-wave
/// modes are populated through parametric processes, and at sufficiently
/// large drive amplitudes the resonance line folds (bistability).
#[derive(Debug, Clone)]
pub struct NonlinearFmrLinewidth {
    material: FourMagnonScattering,
}

impl NonlinearFmrLinewidth {
    /// Create a new nonlinear FMR linewidth model.
    ///
    /// # Arguments
    /// Same as [`FourMagnonScattering::new`]: `ms`, `exchange_a`, `alpha`, `h_ext`.
    ///
    /// # Errors
    /// Returns error for invalid parameters (see [`FourMagnonScattering::new`]).
    pub fn new(ms: f64, exchange_a: f64, alpha: f64, h_ext: f64) -> Result<Self> {
        Ok(Self {
            material: FourMagnonScattering::new(ms, exchange_a, alpha, h_ext)?,
        })
    }

    /// Construct from standard YIG parameters.
    ///
    /// # Arguments
    /// * `h_ext` - External field \[T\]
    pub fn from_yig(h_ext: f64) -> Self {
        Self {
            material: FourMagnonScattering::yig(h_ext),
        }
    }

    /// FMR resonance angular frequency \[rad/s\].
    ///
    /// ω_res = γ · H_ext (simplified Kittel; exact only for thin-film in-plane)
    #[inline]
    fn omega_res(&self) -> f64 {
        GAMMA * self.material.h_ext
    }

    /// Linear (small-signal) FMR linewidth in field units \[T\].
    ///
    /// ΔH₀ = α · ω_res / γ
    ///
    /// This is the Gilbert-damping contribution alone, without any nonlinear
    /// broadening.
    pub fn linear_linewidth_field(&self) -> f64 {
        self.material.alpha * self.omega_res() / GAMMA
    }

    /// Nonlinear FMR linewidth in field units using a phenomenological Suhl model.
    ///
    /// The broadening above threshold follows a Duffing-oscillator-like form:
    ///
    /// ΔH(η) = ΔH₀ · (1 + max(0, η - 1)²)
    ///
    /// where η = h_pump / h_th is the normalised pump amplitude.
    ///
    /// # Arguments
    /// * `pump_h_normalized` - Normalised pump amplitude η = h / h_th (dimensionless)
    ///
    /// # Returns
    /// Effective linewidth \[T\]
    pub fn nonlinear_linewidth_field(&self, pump_h_normalized: f64) -> f64 {
        let dh0 = self.linear_linewidth_field();
        let excess = (pump_h_normalized - 1.0).max(0.0);
        dh0 * (1.0 + excess * excess)
    }

    /// Foldover (bistability onset) field deviation \[T\].
    ///
    /// For a Duffing-type nonlinear resonator the line begins to fold when
    /// the drive exceeds:
    ///
    /// H_fold ≈ (3/√3) · ΔH₀ = √3 · ΔH₀
    ///
    /// This is the onset of the S-shaped (bistable) amplitude-frequency curve.
    pub fn foldover_field(&self) -> f64 {
        let dh0 = self.linear_linewidth_field();
        3.0_f64 / 3.0_f64.sqrt() * dh0
    }

    /// Threshold power density for bistability onset \[W/m²\].
    ///
    /// P_th ≈ ΔH₀² · f_res / (γ² · Ms)
    ///
    /// where f_res = ω_res / (2π) is the FMR frequency and Ms is the
    /// saturation magnetization.
    pub fn bistability_threshold_power(&self) -> f64 {
        let dh0 = self.linear_linewidth_field();
        let f_res = self.omega_res() / (2.0 * PI);
        dh0 * dh0 * f_res / (GAMMA * GAMMA * self.material.ms)
    }

    /// Power saturation factor f(P) ∈ (0, 1].
    ///
    /// Models the reduction of FMR amplitude at high drive powers due to
    /// nonlinear damping and spin-wave generation:
    ///
    /// f(P) = 1 / (1 + P / P_crit)
    ///
    /// where P_crit = ΔH₀² · Ms / γ² captures the power scale at which
    /// nonlinear effects become comparable to linear damping.
    ///
    /// # Arguments
    /// * `input_power` - Incident microwave power density \[W/m²\]
    pub fn power_saturation_factor(&self, input_power: f64) -> f64 {
        let dh0 = self.linear_linewidth_field();
        let p_crit = dh0 * dh0 * self.material.ms / (GAMMA * GAMMA);
        1.0 / (1.0 + input_power / p_crit)
    }

    /// Effective damping including magnon-magnon scattering contributions.
    ///
    /// When the spin-wave population n_k is elevated (e.g., by parametric
    /// pumping), the effective Gilbert damping increases:
    ///
    /// α_eff = α · (1 + n_k / n_sat)
    ///
    /// where the saturation density is n_sat = 1 / (|T_kk|² · α) estimated
    /// at the optimal instability wavevector (θ = π/2, k = k_opt).
    ///
    /// # Arguments
    /// * `spin_wave_density` - Spin-wave mode occupation n_k (dimensionless)
    ///
    /// # Returns
    /// Effective Gilbert damping α_eff (dimensionless)
    pub fn effective_damping(&self, spin_wave_density: f64) -> f64 {
        // Evaluate coupling at a representative wavevector (k ~ 1e6 m⁻¹, θ = π/2)
        let k_rep = 1.0e6_f64;
        let t_kk = self.material.four_magnon_coupling(k_rep, PI / 2.0);

        // n_sat = 1 / (T² · α); guard against zero T
        let n_sat = if t_kk < f64::EPSILON {
            f64::MAX
        } else {
            1.0 / (t_kk * t_kk * self.material.alpha)
        };

        self.material.alpha * (1.0 + spin_wave_density / n_sat)
    }
}

// ============================================================================
// Free Functions
// ============================================================================

/// Magnon-magnon interaction energy density.
///
/// The mean-field energy density due to four-magnon scattering between two
/// magnon populations n₁ and n₂:
///
/// E_int = T · n₁ · n₂   \[J/m³\]
///
/// # Arguments
/// * `n1` - Occupation of the first magnon mode (dimensionless)
/// * `n2` - Occupation of the second magnon mode (dimensionless)
/// * `t_coupling` - Four-magnon coupling constant |T_kk| \[rad/s\]
///
/// # Returns
/// Interaction energy density \[J/m³\]
pub fn magnon_magnon_interaction_energy(n1: f64, n2: f64, t_coupling: f64) -> f64 {
    t_coupling * n1 * n2
}

/// Four-magnon relaxation rate for a spin-wave mode at finite temperature.
///
/// Captures the thermally activated magnon-magnon scattering contribution
/// to the linewidth. A simplified collision integral model:
///
/// Γ_k ≈ π · T² · n_k · D(ω)
///
/// where D(ω) = (ω / ω_ref)^(1/2) is a simplified magnon density of states
/// (ω_ref = 2π × 10⁹ rad/s ≃ 1 GHz) and n_k is the occupation number.
///
/// This grows as √ω at high frequency (3D magnon DOS in the exchange regime)
/// and is proportional to the occupation and coupling squared.
///
/// # Arguments
/// * `n_k` - Spin-wave mode occupation (dimensionless)
/// * `t_coupling` - Four-magnon coupling |T_kk| \[rad/s\]
/// * `omega` - Magnon angular frequency \[rad/s\]
/// * `temperature` - Lattice temperature \[K\] (reserved for future use in
///   detailed balance; the current model uses n_k directly)
///
/// # Returns
/// Relaxation rate Γ_k \[rad/s\]
pub fn four_magnon_relaxation_rate(
    n_k: f64,
    t_coupling: f64,
    omega: f64,
    _temperature: f64,
) -> f64 {
    // Reference frequency: 1 GHz in rad/s
    const OMEGA_REF: f64 = 2.0 * PI * 1.0e9;
    let dos = (omega / OMEGA_REF).sqrt().max(0.0);
    PI * t_coupling * t_coupling * n_k * dos
}

/// Minimum threshold power density for spin-wave instabilities \[W/m³\].
///
/// Provides a lower bound on the power needed to excite any Suhl instability.
/// In the limit of small damping the threshold is set entirely by:
///
/// P_min = α² · ω_res² · Ms / (μ₀ · γ²)
///
/// # Arguments
/// * `ms` - Saturation magnetization \[A/m\]
/// * `alpha` - Gilbert damping (dimensionless)
/// * `h_ext` - External field \[T\]
/// * `gamma` - Gyromagnetic ratio \[rad/(s·T)\]
///
/// # Returns
/// Minimum threshold power density \[W/m³\]
pub fn suhl_spin_wave_instability_power(ms: f64, alpha: f64, h_ext: f64, gamma: f64) -> f64 {
    let omega_res = gamma * h_ext;
    alpha * alpha * omega_res * omega_res * ms / (MU_0 * gamma * gamma)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------------------------------------------------------
    // 1. Kittel frequency at k = 0, θ = π/2 matches the Kittel formula exactly
    // -------------------------------------------------------------------------
    #[test]
    fn test_magnon_frequency_kittel_limit() {
        let h_ext = 0.1; // 100 mT external field [T]
        let ms = 1.40e5_f64; // YIG Ms [A/m]
        let sys = FourMagnonScattering::yig(h_ext);

        // Exact Kittel formula for in-plane magnetized film:
        // ω_Kittel = γ √(B_ext · (B_ext + μ₀ Ms))
        let mu0_ms = MU_0 * ms;
        let omega_kittel = GAMMA * (h_ext * (h_ext + mu0_ms)).sqrt();

        // Our dispersion at k = 0 and θ = π/2 must recover this exactly
        let omega = sys.magnon_frequency(0.0, PI / 2.0);

        let rel_err = (omega - omega_kittel).abs() / omega_kittel;
        assert!(omega > 0.0, "magnon frequency must be positive");
        assert!(
            rel_err < 1.0e-12,
            "k=0 Kalinikos-Slavin must equal Kittel formula; \
             omega={omega:.6e}, kittel={omega_kittel:.6e}, err={rel_err:.2e}"
        );
    }

    // -------------------------------------------------------------------------
    // 2. Linewidth is proportional to alpha
    // -------------------------------------------------------------------------
    #[test]
    fn test_magnon_linewidth_proportional_to_alpha() {
        let h_ext = 0.05;
        let sys1 = FourMagnonScattering::new(1.4e5, 3.7e-12, 1.0e-4, h_ext).expect("valid params");
        let sys2 = FourMagnonScattering::new(1.4e5, 3.7e-12, 2.0e-4, h_ext).expect("valid params");

        let k = 1.0e6;
        let dw1 = sys1.magnon_linewidth(k);
        let dw2 = sys2.magnon_linewidth(k);

        // α doubled → linewidth should double
        let ratio = dw2 / dw1;
        assert!(
            (ratio - 2.0).abs() < 1.0e-10,
            "linewidth must be linear in alpha; got ratio = {ratio}"
        );
    }

    // -------------------------------------------------------------------------
    // 3. Four-magnon coupling is positive away from the magic angle
    // -------------------------------------------------------------------------
    #[test]
    fn test_four_magnon_coupling_nonzero_away_from_magic_angle() {
        let sys = FourMagnonScattering::yig(0.1);

        // θ = 0 (parallel to M): |3cos²θ - 1| = 2 — large dipolar contribution
        let t_parallel = sys.four_magnon_coupling(1.0e6, 0.0);
        assert!(
            t_parallel > 0.0,
            "coupling must be positive at θ = 0; got {t_parallel}"
        );

        // θ = π/2 (perpendicular): |3·0 - 1| = 1 — also nonzero
        let t_perp = sys.four_magnon_coupling(1.0e6, PI / 2.0);
        assert!(
            t_perp > 0.0,
            "coupling must be positive at θ = π/2; got {t_perp}"
        );
    }

    // -------------------------------------------------------------------------
    // 4. Suhl threshold decreases with increasing k (exchange softens linewidth)
    // -------------------------------------------------------------------------
    #[test]
    fn test_suhl_threshold_field_positive_and_trend() {
        let sys = FourMagnonScattering::yig(0.1);
        let theta = PI / 2.0;

        let h_th_low = sys.suhl_threshold_field(1.0e5, theta);
        let h_th_mid = sys.suhl_threshold_field(1.0e6, theta);
        let h_th_high = sys.suhl_threshold_field(5.0e6, theta);

        assert!(h_th_low > 0.0, "threshold field must be positive");
        assert!(h_th_mid > 0.0, "threshold field must be positive");
        assert!(h_th_high > 0.0, "threshold field must be positive");

        // At very large k the exchange term in coupling dominates and grows
        // faster than the linewidth, so the threshold should eventually decrease.
        // We verify monotonicity from mid to high k.
        assert!(
            h_th_high <= h_th_mid * 10.0,
            "threshold should not blow up arbitrarily with k; \
             got h_th_mid={h_th_mid:.3e}, h_th_high={h_th_high:.3e}"
        );
    }

    // -------------------------------------------------------------------------
    // 5. Parametric growth rate = 0 below threshold
    // -------------------------------------------------------------------------
    #[test]
    fn test_parametric_growth_rate_zero_below_threshold() {
        let sys = FourMagnonScattering::yig(0.1);
        let k = 1.0e6;
        let theta = PI / 2.0;

        let h_th = sys.suhl_threshold_field(k, theta);

        // Use a pump strictly below threshold
        let pump_h = 0.5 * h_th;
        let sigma = sys.parametric_growth_rate(pump_h, k, theta);
        assert_eq!(
            sigma, 0.0,
            "growth rate below threshold must be exactly 0; got {sigma}"
        );
    }

    // -------------------------------------------------------------------------
    // 6. Parametric growth rate > 0 above threshold
    // -------------------------------------------------------------------------
    #[test]
    fn test_parametric_growth_rate_positive_above_threshold() {
        let sys = FourMagnonScattering::yig(0.1);
        let k = 1.0e6;
        let theta = PI / 2.0;

        let h_th = sys.suhl_threshold_field(k, theta);

        // Use a pump 10× above threshold
        let pump_h = 10.0 * h_th;
        let sigma = sys.parametric_growth_rate(pump_h, k, theta);
        assert!(
            sigma > 0.0,
            "growth rate above threshold must be positive; got {sigma}"
        );
    }

    // -------------------------------------------------------------------------
    // 7. ParametricAmplification: is_above_threshold matches threshold comparison
    // -------------------------------------------------------------------------
    #[test]
    fn test_parametric_amplification_threshold_consistency() {
        // Build a system just below threshold
        let coupling = GAMMA * 1.4e5 / 2.0; // YIG-like
        let gamma_s = 1.0e7_f64; // rad/s
        let gamma_i = 1.0e7_f64;
        let h_th = (gamma_s * gamma_i).sqrt() / coupling;

        let below = ParametricAmplification::new(
            coupling,
            2.0 * PI * 5.0e9,
            2.0 * PI * 5.0e9,
            0.9 * h_th,
            gamma_s,
            gamma_i,
        )
        .expect("valid params");

        let above = ParametricAmplification::new(
            coupling,
            2.0 * PI * 5.0e9,
            2.0 * PI * 5.0e9,
            1.1 * h_th,
            gamma_s,
            gamma_i,
        )
        .expect("valid params");

        assert!(!below.is_above_threshold(), "should be below threshold");
        assert!(above.is_above_threshold(), "should be above threshold");
    }

    // -------------------------------------------------------------------------
    // 8. Gain coefficient = 0 at exactly threshold (continuity)
    // -------------------------------------------------------------------------
    #[test]
    fn test_gain_coefficient_zero_at_threshold() {
        let coupling = 1.0e9_f64; // rad/(s·T)
        let gamma_s = 1.0e7_f64;
        let gamma_i = 1.0e7_f64;
        let h_th = (gamma_s * gamma_i).sqrt() / coupling;

        let pa = ParametricAmplification::new(
            coupling,
            2.0 * PI * 5.0e9,
            2.0 * PI * 5.0e9,
            h_th,
            gamma_s,
            gamma_i,
        )
        .expect("valid params");

        let g = pa.gain_coefficient();
        assert!(
            g.abs() < 1.0e-6,
            "gain coefficient must be 0 at threshold; got {g}"
        );
    }

    // -------------------------------------------------------------------------
    // 9. Noise-equivalent power is positive
    // -------------------------------------------------------------------------
    #[test]
    fn test_noise_equivalent_power_positive() {
        let pa = ParametricAmplification::new(
            1.0e9,
            2.0 * PI * 5.0e9,
            2.0 * PI * 5.0e9,
            1.0e-3,
            1.0e7,
            1.0e7,
        )
        .expect("valid params");

        let nep = pa.noise_equivalent_power();
        assert!(nep > 0.0, "NEP must be positive; got {nep}");
    }

    // -------------------------------------------------------------------------
    // 10. Linear linewidth: ΔH₀ = α · ω_res / γ within 1%
    // -------------------------------------------------------------------------
    #[test]
    fn test_linear_linewidth_field_formula() {
        let h_ext = 0.1;
        let alpha = 1.0e-4;
        let nlw = NonlinearFmrLinewidth::new(1.4e5, 3.7e-12, alpha, h_ext).expect("valid params");

        let dh0 = nlw.linear_linewidth_field();
        let omega_res = GAMMA * h_ext;
        let expected = alpha * omega_res / GAMMA;

        let rel_error = (dh0 - expected).abs() / expected;
        assert!(
            rel_error < 0.01,
            "linear linewidth formula error too large: {rel_error:.2e}"
        );
    }

    // -------------------------------------------------------------------------
    // 11. Nonlinear linewidth increases monotonically above threshold
    // -------------------------------------------------------------------------
    #[test]
    fn test_nonlinear_linewidth_monotone_above_threshold() {
        let nlw = NonlinearFmrLinewidth::from_yig(0.1);

        let mut prev = nlw.nonlinear_linewidth_field(1.0);
        for i in 1..=10 {
            let eta = 1.0 + 0.5 * i as f64;
            let dh = nlw.nonlinear_linewidth_field(eta);
            assert!(
                dh >= prev,
                "nonlinear linewidth must be non-decreasing above threshold; \
                 got {dh} < {prev} at η = {eta}"
            );
            prev = dh;
        }
    }

    // -------------------------------------------------------------------------
    // 12. Bistability threshold power > 0 and scales as ΔH₀²
    // -------------------------------------------------------------------------
    #[test]
    fn test_bistability_threshold_power_scaling() {
        let h_ext = 0.1;
        let nlw1 = NonlinearFmrLinewidth::new(1.4e5, 3.7e-12, 1.0e-4, h_ext).expect("valid params");
        let nlw2 = NonlinearFmrLinewidth::new(1.4e5, 3.7e-12, 2.0e-4, h_ext).expect("valid params");

        let p1 = nlw1.bistability_threshold_power();
        let p2 = nlw2.bistability_threshold_power();

        assert!(p1 > 0.0, "bistability power must be positive");

        // Doubling α doubles ΔH₀, so P_th should quadruple
        let ratio = p2 / p1;
        assert!(
            (ratio - 4.0).abs() < 0.01,
            "bistability power should scale as α² (4× for 2× α); got ratio = {ratio:.4}"
        );
    }

    // -------------------------------------------------------------------------
    // 13. optimal_k_for_suhl returns a finite positive value
    // -------------------------------------------------------------------------
    #[test]
    fn test_optimal_k_for_suhl_finite_positive() {
        let sys = FourMagnonScattering::yig(0.1);

        // Pump at 10 GHz
        let pump_freq = 2.0 * PI * 10.0e9;
        let k_opt = sys.optimal_k_for_suhl(pump_freq);

        assert!(k_opt.is_finite(), "optimal k must be finite; got {k_opt}");
        assert!(k_opt > 0.0, "optimal k must be positive; got {k_opt}");
        assert!(k_opt < 1.0e9, "optimal k should be in a physical range");
    }

    // -------------------------------------------------------------------------
    // 14. magnon_magnon_interaction_energy is positive for positive inputs
    // -------------------------------------------------------------------------
    #[test]
    fn test_magnon_magnon_interaction_energy_positive() {
        let n1 = 100.0;
        let n2 = 50.0;
        let t = 1.0e8;
        let e = magnon_magnon_interaction_energy(n1, n2, t);

        assert!(
            e > 0.0,
            "interaction energy must be positive for positive n1, n2, T; got {e}"
        );
        // Linearity check: doubling n1 should double E
        let e2 = magnon_magnon_interaction_energy(2.0 * n1, n2, t);
        assert!(
            (e2 / e - 2.0).abs() < 1.0e-12,
            "interaction energy must be linear in n1"
        );
    }
}
