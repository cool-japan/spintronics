//! Optomagnonics — Light-Magnon Interaction in Ferromagnetic Insulators
//!
//! Optomagnonics couples optical photons to magnons via the magneto-optical Kerr
//! and Faraday effects (spin-orbit interaction). The key microscopic mechanism is
//! inelastic Brillouin light scattering (BLS): an incident photon is scattered into
//! a shifted photon and a magnon (Stokes) or absorbs a magnon (anti-Stokes).
//!
//! This module provides:
//! - [`BrillouinScattering`] — single-mode BLS scattering rate and cross-section.
//! - [`OptomagnonicCoupling`] — single-magnon optomagnonic coupling and Kerr shift.
//! - [`MicrowaveToOptical`] — three-mode microwave-magnon-optical transducer efficiency.
//! - [`MagnonicFrequencyComb`] — parametric magnon comb: equally-spaced magnon lines.
//!
//! # References
//!
//! - A. Osada et al., Phys. Rev. Lett. **116**, 223601 (2016) — optomagnonic cavity with YIG
//! - X. Zhang et al., Sci. Adv. **2**, e1501286 (2016) — cavity optomagnonics
//! - R. Hisatomi et al., Phys. Rev. B **93**, 174427 (2016) — bidirectional conversion
//! - J. A. Haigh et al., Phys. Rev. A **92**, 063845 (2015) — magnonic Kerr nonlinearity
//! - N. Hashemi et al., Phys. Rev. Research (2020) — efficiency formula for transducers

use std::f64::consts::PI;

use crate::error::{self, Result};

// ---------------------------------------------------------------------------
// BrillouinScattering
// ---------------------------------------------------------------------------

/// Brillouin light scattering: inelastic interaction between optical photons and magnons.
///
/// The optomagnonic coupling Hamiltonian (single mode) is:
///
///   H_om = ℏ g_om (a†_opt a_opt)(b + b†)
///
/// where a_opt is the optical photon operator and b the magnon operator.
/// The scattering rate (from Fermi's golden rule) for Stokes/anti-Stokes is:
///
///   Γ = g_om² n_magnon / (κ² + Δ²)
///
/// and the Brillouin differential cross-section peaks at 90° scattering angle.
///
/// # Examples
///
/// ```rust
/// use spintronics::cavity::BrillouinScattering;
///
/// let bls = BrillouinScattering::yig_optical();
/// let rate = bls.scattering_rate(0.0, 1.0);
/// assert!(rate > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct BrillouinScattering {
    /// Magnon frequency \[rad/s\].
    pub magnon_frequency: f64,
    /// Optical photon frequency \[rad/s\].
    pub optical_frequency: f64,
    /// Optomagnonic coupling constant g_om \[rad/s\].
    pub coupling: f64,
}

impl BrillouinScattering {
    /// Create a new `BrillouinScattering` instance.
    ///
    /// # Errors
    ///
    /// Returns `Err` if any frequency or coupling is negative.
    pub fn new(magnon_frequency: f64, optical_frequency: f64, coupling: f64) -> Result<Self> {
        if magnon_frequency < 0.0 {
            return Err(error::invalid_param(
                "magnon_frequency",
                "magnon frequency must be non-negative",
            ));
        }
        if optical_frequency < 0.0 {
            return Err(error::invalid_param(
                "optical_frequency",
                "optical photon frequency must be non-negative",
            ));
        }
        if coupling < 0.0 {
            return Err(error::invalid_param(
                "coupling",
                "optomagnonic coupling must be non-negative",
            ));
        }
        Ok(Self {
            magnon_frequency,
            optical_frequency,
            coupling,
        })
    }

    /// YIG optical BLS preset.
    ///
    /// - ω_m = 2π × 10 GHz  (uniform Kittel mode)
    /// - ω_opt = 2π × 193 THz  (1550 nm telecom band)
    /// - g_om = 2π × 1 MHz   (typical YIG sphere value)
    pub fn yig_optical() -> Self {
        Self {
            magnon_frequency: 2.0 * PI * 10.0e9,
            optical_frequency: 2.0 * PI * 193.0e12,
            coupling: 2.0 * PI * 1.0e6,
        }
    }

    /// Stokes/anti-Stokes scattering rate.
    ///
    ///   Γ = g_om² · n_magnon / (κ_eff² + detuning²)
    ///
    /// where κ_eff is implicitly contained in the linewidth denominator.
    /// Here we use a simplified single-pole susceptibility:
    ///
    ///   Γ = g_om² · n_magnon / (1 + (detuning / g_om)²)
    ///
    /// which gives units of \[rad/s\] · n_magnon (dimensionless magnon occupation).
    ///
    /// The full expression requires knowledge of the optical cavity linewidth κ_o;
    /// the formula above assumes κ_o = g_om for the normalisation.  For actual
    /// experiments, replace the denominator with `(κ_o² + detuning²)`.
    ///
    /// # Arguments
    ///
    /// - `detuning`: laser-cavity detuning Δ = ω_laser − ω_opt \[rad/s\].
    /// - `magnon_density`: average magnon occupation ⟨b†b⟩ (dimensionless).
    pub fn scattering_rate(&self, detuning: f64, magnon_density: f64) -> f64 {
        let g = self.coupling;
        let denom = g.powi(2) + detuning.powi(2);
        g.powi(2) * magnon_density / denom
    }

    /// Output photon frequency for Stokes (ω_out = ω_opt − ω_m) or anti-Stokes
    /// (ω_out = ω_opt + ω_m) scattering.
    ///
    /// # Arguments
    ///
    /// - `stokes`: `true` for Stokes (magnon emission), `false` for anti-Stokes (magnon absorption).
    pub fn frequency_shift(&self, stokes: bool) -> f64 {
        if stokes {
            self.optical_frequency - self.magnon_frequency
        } else {
            self.optical_frequency + self.magnon_frequency
        }
    }

    /// Relative Brillouin differential cross-section as a function of scattering angle.
    ///
    /// Models the angular dependence of inelastic BLS from a spherical YIG sample:
    ///
    ///   σ(θ) ∝ 1 / (1 + (2θ/π − 1)²)
    ///
    /// This Lorentzian-on-angle form peaks at θ = π/2 (90°, backscattering geometry
    /// commonly used in BLS spectroscopy) and returns values in [0, 1].
    ///
    /// # Arguments
    ///
    /// - `scattering_angle_rad`: scattering angle θ ∈ [0, π] \[rad\].
    pub fn cross_section(&self, scattering_angle_rad: f64) -> f64 {
        // Normalised Lorentzian centred at θ = π/2.
        let x = 2.0 * scattering_angle_rad / PI - 1.0;
        1.0 / (1.0 + x.powi(2))
    }
}

// ---------------------------------------------------------------------------
// OptomagnonicCoupling
// ---------------------------------------------------------------------------

/// Effective optomagnonic coupling strength and Kerr nonlinearity.
///
/// The single-magnon optomagnonic coupling g_om quantifies the frequency shift of
/// an optical mode per magnon.  It is enhanced by the square root of the intra-cavity
/// photon number n_opt via the effective coupling g_eff = g_om·√n_opt.
///
/// The magnonic Kerr nonlinearity shifts the magnon frequency by Δω = K·n_photons,
/// where K is the Kerr coefficient.
///
/// # Examples
///
/// ```rust
/// use spintronics::cavity::OptomagnonicCoupling;
///
/// let oc = OptomagnonicCoupling::yig();
/// let g_eff = oc.effective_coupling(1000.0);
/// assert!(g_eff > oc.g_om);
/// ```
#[derive(Debug, Clone)]
pub struct OptomagnonicCoupling {
    /// Single-magnon optomagnonic coupling g_om \[rad/s\].
    pub g_om: f64,
    /// Magnonic Kerr coefficient K [rad/s per intra-cavity photon].
    pub kerr_strength: f64,
}

impl OptomagnonicCoupling {
    /// Create a new `OptomagnonicCoupling`.
    ///
    /// # Errors
    ///
    /// Returns `Err` if `g_om` or `kerr_strength` is negative.
    pub fn new(g_om: f64, kerr_strength: f64) -> Result<Self> {
        if g_om < 0.0 {
            return Err(error::invalid_param(
                "g_om",
                "optomagnonic coupling must be non-negative",
            ));
        }
        if kerr_strength < 0.0 {
            return Err(error::invalid_param(
                "kerr_strength",
                "Kerr coefficient must be non-negative",
            ));
        }
        Ok(Self {
            g_om,
            kerr_strength,
        })
    }

    /// YIG optomagnonic preset.
    ///
    /// - g_om = 2π × 1 MHz   (Haigh et al. 2015 value for 1 mm YIG sphere)
    /// - K    = 2π × 100 Hz  (magnonic Kerr — small but measurable)
    pub fn yig() -> Self {
        Self {
            g_om: 2.0 * PI * 1.0e6,
            kerr_strength: 2.0 * PI * 100.0,
        }
    }

    /// Effective coupling enhanced by intra-cavity optical field.
    ///
    ///   g_eff = g_om · √(n_opt)
    ///
    /// # Arguments
    ///
    /// - `n_optical_photons`: mean intra-cavity photon number n_opt ≥ 0.
    pub fn effective_coupling(&self, n_optical_photons: f64) -> f64 {
        self.g_om * n_optical_photons.max(0.0).sqrt()
    }

    /// Kerr frequency shift per intra-cavity photon.
    ///
    ///   Δω = K · n_photons
    pub fn kerr_shift(&self, n_photons: f64) -> f64 {
        self.kerr_strength * n_photons
    }
}

// ---------------------------------------------------------------------------
// MicrowaveToOptical
// ---------------------------------------------------------------------------

/// Three-mode microwave ↔ magnon ↔ optical transducer.
///
/// Mediates coherent conversion between microwave and optical photons via an
/// intermediate magnon mode.  The conversion efficiency follows the
/// Hashemi-Mahmoodian formula (generalised impedance-matching condition):
///
///   η = 4 C_me C_mo / (1 + C_me + C_mo)²
///
/// where the cooperativities are:
///
///   C_me = g_me² / (κ_m · γ_m)   (microwave-magnon)
///   C_mo = g_mo² / (κ_o · γ_m)   (optical-magnon)
///
/// Maximum efficiency η = 1 is reached when C_me = C_mo = C → ∞; at unit
/// cooperativity C_me = C_mo = 1 one obtains η = 4/9 ≈ 0.444.
///
/// # Examples
///
/// ```rust
/// use spintronics::cavity::MicrowaveToOptical;
///
/// let t = MicrowaveToOptical::yig_telecom();
/// let eta = t.conversion_efficiency();
/// assert!(eta >= 0.0 && eta <= 1.0);
/// ```
#[derive(Debug, Clone)]
pub struct MicrowaveToOptical {
    /// Microwave photon frequency \[rad/s\].
    pub microwave_freq: f64,
    /// Magnon mode frequency \[rad/s\].
    pub magnon_freq: f64,
    /// Optical photon frequency \[rad/s\].
    pub optical_freq: f64,
    /// Microwave-magnon coupling g_me \[rad/s\].
    pub g_me: f64,
    /// Optical-magnon coupling g_mo \[rad/s\].
    pub g_mo: f64,
    /// Microwave cavity linewidth κ_m \[rad/s\].
    pub kappa_m: f64,
    /// Optical cavity linewidth κ_o \[rad/s\].
    pub kappa_o: f64,
    /// Magnon damping rate γ_m \[rad/s\].
    pub gamma_m: f64,
}

impl MicrowaveToOptical {
    /// Create a new `MicrowaveToOptical` transducer.
    ///
    /// # Errors
    ///
    /// Returns `Err` if any frequency, coupling, or decay rate is negative.
    pub fn new(
        microwave_freq: f64,
        magnon_freq: f64,
        optical_freq: f64,
        g_me: f64,
        g_mo: f64,
        kappa_m: f64,
        kappa_o: f64,
        gamma_m: f64,
    ) -> Result<Self> {
        let params = [
            ("microwave_freq", microwave_freq),
            ("magnon_freq", magnon_freq),
            ("optical_freq", optical_freq),
            ("g_me", g_me),
            ("g_mo", g_mo),
            ("kappa_m", kappa_m),
            ("kappa_o", kappa_o),
            ("gamma_m", gamma_m),
        ];
        for (name, val) in &params {
            if *val < 0.0 {
                return Err(error::invalid_param(
                    name,
                    "transducer parameter must be non-negative",
                ));
            }
        }
        Ok(Self {
            microwave_freq,
            magnon_freq,
            optical_freq,
            g_me,
            g_mo,
            kappa_m,
            kappa_o,
            gamma_m,
        })
    }

    /// YIG-based microwave-to-telecom transducer (representative lab parameters).
    ///
    /// - f_MW   = 10 GHz (microwave cavity)
    /// - f_m    ≈ 10 GHz (Kittel mode at 80 kA/m bias)
    /// - f_opt  = 193 THz (telecom C-band)
    /// - g_me   = 2π × 100 MHz (strong MW-magnon coupling)
    /// - g_mo   = 2π × 1 MHz   (weak optical coupling)
    /// - κ_m   = 2π × 1 MHz
    /// - κ_o   = 2π × 10 MHz  (higher loss in optical cavity)
    /// - γ_m   ≈ 2π × 1 MHz   (YIG magnon linewidth at room temperature)
    pub fn yig_telecom() -> Self {
        Self {
            microwave_freq: 2.0 * PI * 10.0e9,
            magnon_freq: 2.0 * PI * 10.0e9,
            optical_freq: 2.0 * PI * 193.0e12,
            g_me: 2.0 * PI * 100.0e6,
            g_mo: 2.0 * PI * 1.0e6,
            kappa_m: 2.0 * PI * 1.0e6,
            kappa_o: 2.0 * PI * 10.0e6,
            gamma_m: 2.0 * PI * 1.0e6,
        }
    }

    /// Microwave-magnon cooperativity C_me = g_me² / (κ_m · γ_m).
    fn cooperativity_me(&self) -> f64 {
        self.g_me.powi(2) / (self.kappa_m * self.gamma_m)
    }

    /// Optical-magnon cooperativity C_mo = g_mo² / (κ_o · γ_m).
    fn cooperativity_mo(&self) -> f64 {
        self.g_mo.powi(2) / (self.kappa_o * self.gamma_m)
    }

    /// Quantum transduction efficiency η ∈ [0, 1].
    ///
    ///   η = 4 C_me C_mo / (1 + C_me + C_mo)²
    pub fn conversion_efficiency(&self) -> f64 {
        let c_me = self.cooperativity_me();
        let c_mo = self.cooperativity_mo();
        let denom = (1.0 + c_me + c_mo).powi(2);
        if denom == 0.0 {
            return 0.0;
        }
        (4.0 * c_me * c_mo / denom).min(1.0)
    }

    /// Effective transduction bandwidth γ_eff = γ_m · (1 + C_me + C_mo).
    pub fn bandwidth(&self) -> f64 {
        let c_me = self.cooperativity_me();
        let c_mo = self.cooperativity_mo();
        self.gamma_m * (1.0 + c_me + c_mo)
    }

    /// `true` when C_me ≈ C_mo (impedance-matched condition, within 10% relative tolerance).
    pub fn is_impedance_matched(&self) -> bool {
        let c_me = self.cooperativity_me();
        let c_mo = self.cooperativity_mo();
        if c_me + c_mo == 0.0 {
            return true;
        }
        ((c_me - c_mo) / (c_me + c_mo)).abs() < 0.1
    }
}

// ---------------------------------------------------------------------------
// MagnonicFrequencyComb
// ---------------------------------------------------------------------------

/// Magnonic frequency comb: an equally-spaced array of magnon spectral lines.
///
/// Analogous to optical frequency combs, a magnonic comb arises from parametric
/// driving of a nonlinear magnon system.  The comb consists of N lines separated
/// by the repetition frequency f_rep, centred on the carrier f_carrier, with
/// Lorentzian linewidths broadened by decoherence.
///
/// Phase noise follows the standard Leeson model:
///
///   L(f) = decoherence² / (2π f)²   [dBc/Hz, linearised]
///
/// which decreases as 1/f² away from the carrier.
///
/// # Examples
///
/// ```rust
/// use spintronics::cavity::MagnonicFrequencyComb;
///
/// let comb = MagnonicFrequencyComb::yig_microwave(0.1); // 100 MHz repetition rate
/// let spec = comb.comb_spectrum();
/// assert_eq!(spec.len(), comb.n_lines);
/// ```
#[derive(Debug, Clone)]
pub struct MagnonicFrequencyComb {
    /// Repetition rate (line spacing) f_rep \[Hz\].
    pub f_repetition: f64,
    /// Number of comb lines.
    pub n_lines: usize,
    /// Linewidth broadening / decoherence parameter \[Hz\].
    pub decoherence: f64,
    /// Carrier (centre) frequency \[Hz\].
    pub carrier: f64,
}

impl MagnonicFrequencyComb {
    /// Create a new `MagnonicFrequencyComb`.
    ///
    /// # Errors
    ///
    /// Returns `Err` if `f_repetition` or `carrier` is negative, `n_lines` is zero,
    /// or `decoherence` is negative.
    pub fn new(f_repetition: f64, n_lines: usize, decoherence: f64, carrier: f64) -> Result<Self> {
        if f_repetition < 0.0 {
            return Err(error::invalid_param(
                "f_repetition",
                "repetition rate must be non-negative",
            ));
        }
        if n_lines == 0 {
            return Err(error::invalid_param(
                "n_lines",
                "comb must have at least one line",
            ));
        }
        if decoherence < 0.0 {
            return Err(error::invalid_param(
                "decoherence",
                "decoherence must be non-negative",
            ));
        }
        if carrier < 0.0 {
            return Err(error::invalid_param(
                "carrier",
                "carrier frequency must be non-negative",
            ));
        }
        Ok(Self {
            f_repetition,
            n_lines,
            decoherence,
            carrier,
        })
    }

    /// Typical YIG microwave comb.
    ///
    /// - Carrier: 10 GHz (Kittel mode)
    /// - n_lines: 21 lines  (±10 sidebands)
    /// - decoherence: 1 MHz (YIG linewidth)
    /// - f_rep: `f_rep_ghz` GHz
    pub fn yig_microwave(f_rep_ghz: f64) -> Self {
        Self {
            f_repetition: f_rep_ghz * 1.0e9,
            n_lines: 21,
            decoherence: 1.0e6,
            carrier: 10.0e9,
        }
    }

    /// Comb spectral lines \[Hz\].
    ///
    /// Returns `n_lines` frequencies:
    ///
    ///   f_n = carrier + n · f_repetition
    ///
    /// for n ∈ {-⌊n_lines/2⌋, …, +⌊n_lines/2⌋} (centred at carrier).
    /// If `n_lines` is even the range is symmetric around carrier ± f_rep/2.
    pub fn comb_spectrum(&self) -> Vec<f64> {
        let half = (self.n_lines / 2) as i64;
        let start = -half;
        // For even n_lines we produce n_lines points starting from -(n_lines/2).
        (start..start + self.n_lines as i64)
            .map(|n| self.carrier + n as f64 * self.f_repetition)
            .collect()
    }

    /// Phase noise spectral density at offset frequency `offset_freq_hz` \[Hz\].
    ///
    ///   L(f) = decoherence² / (2π f)²   (Lorentzian / Leeson, −20 dB/decade)
    ///
    /// Returns linear power spectral density (not dBc/Hz).
    ///
    /// # Arguments
    ///
    /// - `offset_freq_hz`: offset from carrier \[Hz\]; must be > 0.
    pub fn phase_noise(&self, offset_freq_hz: f64) -> f64 {
        if offset_freq_hz <= 0.0 {
            return f64::INFINITY;
        }
        let omega_offset = 2.0 * PI * offset_freq_hz;
        self.decoherence.powi(2) / omega_offset.powi(2)
    }

    /// Coherence time τ_c = 1 / (π · decoherence) \[s\].
    ///
    /// Estimated from the Fourier-transform relation between a Lorentzian lineshape
    /// of FWHM decoherence and its time-domain exponential decay.
    pub fn coherence_time(&self) -> f64 {
        if self.decoherence == 0.0 {
            return f64::INFINITY;
        }
        1.0 / (PI * self.decoherence)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // 1. Scattering rate is positive for positive magnon density.
    #[test]
    fn test_scattering_rate_positive() {
        let bls = BrillouinScattering::yig_optical();
        let rate = bls.scattering_rate(0.0, 5.0);
        assert!(rate > 0.0, "scattering rate must be positive, got {rate}");
    }

    // 2. Stokes output frequency is lower than optical carrier.
    #[test]
    fn test_stokes_lower_frequency() {
        let bls = BrillouinScattering::yig_optical();
        let f_stokes = bls.frequency_shift(true);
        assert!(
            f_stokes < bls.optical_frequency,
            "Stokes frequency {f_stokes} should be below optical {freq}",
            freq = bls.optical_frequency
        );
    }

    // 3. Anti-Stokes output frequency is higher than optical carrier.
    #[test]
    fn test_antistokes_higher_frequency() {
        let bls = BrillouinScattering::yig_optical();
        let f_as = bls.frequency_shift(false);
        assert!(
            f_as > bls.optical_frequency,
            "anti-Stokes frequency {f_as} should be above optical {freq}",
            freq = bls.optical_frequency
        );
    }

    // 4. Cross-section peaks (maximum) at 90°.
    #[test]
    fn test_cross_section_at_pi_over_2() {
        let bls = BrillouinScattering::yig_optical();
        let at_90 = bls.cross_section(PI / 2.0);
        let at_30 = bls.cross_section(PI / 6.0);
        let at_150 = bls.cross_section(5.0 * PI / 6.0);
        // 90° should give the Lorentzian peak value = 1.0
        assert!(
            approx(at_90, 1.0, 1e-12),
            "cross section at 90° should be 1.0, got {at_90}"
        );
        assert!(at_90 > at_30, "σ(90°) should exceed σ(30°)");
        assert!(at_90 > at_150, "σ(90°) should exceed σ(150°)");
    }

    // 5. Conversion efficiency ∈ [0, 1].
    #[test]
    fn test_conversion_efficiency_in_0_1() {
        let t = MicrowaveToOptical::yig_telecom();
        let eta = t.conversion_efficiency();
        assert!(
            (0.0..=1.0).contains(&eta),
            "efficiency {eta} must be in [0, 1]"
        );
    }

    // 6. At unit cooperativity (C_me = C_mo = 1): η = 4/9 ≈ 0.4444.
    #[test]
    fn test_conversion_efficiency_unit_cooperativity() {
        // η = 4*1*1/(1+1+1)^2 = 4/9
        // Choose parameters such that C_me = g_me^2/(kappa_m*gamma_m) = 1
        // and C_mo = g_mo^2/(kappa_o*gamma_m) = 1.
        let gamma_m = 1.0e6_f64;
        let kappa_m = 1.0e6_f64;
        let kappa_o = 1.0e6_f64;
        // g_me^2 = kappa_m * gamma_m  → g_me = sqrt(1e12) = 1e6
        let g_me = (kappa_m * gamma_m).sqrt();
        let g_mo = (kappa_o * gamma_m).sqrt();
        let t = MicrowaveToOptical::new(
            1.0e10,
            1.0e10,
            2.0 * PI * 193.0e12,
            g_me,
            g_mo,
            kappa_m,
            kappa_o,
            gamma_m,
        )
        .unwrap();
        let eta = t.conversion_efficiency();
        let expected = 4.0 / 9.0;
        assert!(
            approx(eta, expected, 1e-10),
            "η at unit cooperativity should be 4/9={expected:.6}, got {eta:.6}"
        );
    }

    // 7. Bandwidth increases with cooperativity.
    #[test]
    fn test_bandwidth_increases_with_c() {
        let t_low = MicrowaveToOptical::new(
            1.0e10,
            1.0e10,
            2.0 * PI * 193.0e12,
            1.0e5,
            1.0e5,
            1.0e6,
            1.0e6,
            1.0e6,
        )
        .unwrap();
        let t_high = MicrowaveToOptical::new(
            1.0e10,
            1.0e10,
            2.0 * PI * 193.0e12,
            1.0e7,
            1.0e7,
            1.0e6,
            1.0e6,
            1.0e6,
        )
        .unwrap();
        assert!(
            t_high.bandwidth() > t_low.bandwidth(),
            "higher cooperativity should yield larger bandwidth"
        );
    }

    // 8. comb_spectrum has exactly n_lines entries.
    #[test]
    fn test_comb_spectrum_length() {
        let comb = MagnonicFrequencyComb::new(100.0e6, 15, 1.0e6, 10.0e9).unwrap();
        let spec = comb.comb_spectrum();
        assert_eq!(
            spec.len(),
            15,
            "comb spectrum should have 15 lines, got {}",
            spec.len()
        );
    }

    // 9. Consecutive comb lines are separated by f_repetition.
    #[test]
    fn test_comb_spacing_equals_f_rep() {
        let f_rep = 500.0e6_f64;
        let comb = MagnonicFrequencyComb::new(f_rep, 11, 1.0e6, 10.0e9).unwrap();
        let spec = comb.comb_spectrum();
        for i in 1..spec.len() {
            let spacing = spec[i] - spec[i - 1];
            assert!(
                approx(spacing, f_rep, 1.0),
                "spacing at index {i} is {spacing}, expected {f_rep}"
            );
        }
    }

    // 10. Phase noise decreases as offset frequency increases (1/f² scaling).
    #[test]
    fn test_phase_noise_decreases_with_offset() {
        let comb = MagnonicFrequencyComb::yig_microwave(0.1);
        let l_low = comb.phase_noise(1.0e3);
        let l_mid = comb.phase_noise(1.0e6);
        let l_high = comb.phase_noise(1.0e9);
        assert!(
            l_low > l_mid && l_mid > l_high,
            "phase noise should decrease with offset: L(1kHz)={l_low:.3e}, L(1MHz)={l_mid:.3e}, L(1GHz)={l_high:.3e}"
        );
    }
}
