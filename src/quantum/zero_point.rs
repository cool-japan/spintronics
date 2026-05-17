//! Zero-point fluctuations of quantised spin-wave modes.
//!
//! In a confined magnetic nanostructure the spin-wave spectrum is discrete.
//! Each mode n carries a *zero-point* (vacuum) energy `E_n = ℏω_n / 2` and a
//! corresponding amplitude of quantum fluctuations
//!
//! ```text
//! u_n = √(ℏ / (2 m_eff ω_n))   [m]
//! ```
//!
//! analogous to a quantum harmonic oscillator with effective mass `m_eff`.
//!
//! At finite temperature the free energy is
//!
//! ```text
//! F = Σ_n ℏω_n [1/2 + n_BE(ω_n, T)]
//! n_BE = 1 / [exp(ℏω_n / kT) - 1]
//! ```
//!
//! The derivative of F with respect to the confinement length gives a
//! Casimir-like magnon pressure (quantum confinement force).
//!
//! # Physical Context
//!
//! Zero-point spin fluctuations are responsible for the *quantum depletion* of
//! the ordered moment (Anderson's sublattice magnetisation reduction), and for
//! fluctuation-driven effects such as order-by-disorder in frustrated magnets.
//!
//! # References
//!
//! - P. W. Anderson, *Phys. Rev.* **86**, 694 (1952) — quantum corrections in AFMs.
//! - T. Holstein, H. Primakoff, *Phys. Rev.* **58**, 1098 (1940).

use crate::constants::{HBAR, KB};
use crate::error::{self, Result};
use crate::spinwave::QuantizedModes;

// ---------------------------------------------------------------------------
// Main struct
// ---------------------------------------------------------------------------

/// Zero-point fluctuation amplitudes and vacuum energetics for a set of
/// discrete spin-wave modes in a confined nanostructure.
///
/// Each mode is characterised by its frequency `ω_n` \[rad/s\].  The effective
/// mass `m_eff` \[kg\] sets the oscillator energy scale through the harmonic
/// oscillator analogy; in practice it is often taken as the magnon effective
/// mass or the volume-weighted spin-inertia.
pub struct ZeroPointFluctuations {
    /// Mode index and angular frequency: (mode_index, ω_n \[rad/s\]).
    pub modes: Vec<(usize, f64)>,
    /// Effective mass parameter \[kg\].
    pub m_eff: f64,
    /// Spin quantum number S \[dimensionless\].
    pub spin: f64,
    /// Effective mode volume V \[m³\].
    pub volume: f64,
}

impl ZeroPointFluctuations {
    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    /// Create a zero-point fluctuation calculator from explicit mode list.
    ///
    /// # Arguments
    ///
    /// * `modes`  — pairs (mode_index, ω_n \[rad/s\]) with all ω_n > 0.
    /// * `m_eff`  — effective mass \[kg\], must be positive.
    /// * `spin`   — spin quantum number S \[dimensionless\], must be positive.
    /// * `volume` — effective mode volume \[m³\], must be positive.
    ///
    /// # Errors
    ///
    /// Returns [`crate::error::Error::InvalidParameter`] if any parameter
    /// violates its positivity constraint or any frequency is ≤ 0.
    pub fn new(modes: Vec<(usize, f64)>, m_eff: f64, spin: f64, volume: f64) -> Result<Self> {
        if m_eff <= 0.0 {
            return Err(error::invalid_param(
                "m_eff",
                "effective mass must be positive",
            ));
        }
        if volume <= 0.0 {
            return Err(error::invalid_param(
                "volume",
                "mode volume must be positive",
            ));
        }
        if spin <= 0.0 {
            return Err(error::invalid_param(
                "spin",
                "spin quantum number must be positive",
            ));
        }
        for (idx, omega) in &modes {
            if *omega <= 0.0 {
                return Err(error::invalid_param(
                    "omega",
                    &format!("mode {idx}: frequency must be positive, got {omega}"),
                ));
            }
        }
        Ok(Self {
            modes,
            m_eff,
            spin,
            volume,
        })
    }

    /// Construct from an existing [`QuantizedModes`] stripe geometry calculator.
    ///
    /// Calls [`QuantizedModes::stripe_mode_frequencies`] to obtain `(mode_n, ω_n)`
    /// pairs and stores them for subsequent fluctuation calculations.
    ///
    /// # Arguments
    ///
    /// * `quantized_modes` — pre-built [`QuantizedModes`] for the nanostructure.
    /// * `h_ext`  — external magnetic field \[T\] for the dispersion.
    /// * `n_modes` — number of modes to extract (1, 2, ..., n_modes).
    /// * `spin`   — spin quantum number S \[dimensionless\].
    /// * `m_eff`  — effective mass \[kg\].
    /// * `volume` — effective mode volume \[m³\].
    ///
    /// # Errors
    ///
    /// Returns an error if the underlying frequency computation fails or if
    /// the parameters are invalid.
    pub fn from_quantized_modes(
        quantized_modes: &QuantizedModes,
        h_ext: f64,
        n_modes: usize,
        spin: f64,
        m_eff: f64,
        volume: f64,
    ) -> Result<Self> {
        // stripe_mode_frequencies returns Vec<(usize, f64)> = (mode_number, omega)
        let freq_pairs = quantized_modes.stripe_mode_frequencies(h_ext, n_modes)?;
        // Build (mode_index, omega) list — use the returned mode_number as index
        let modes: Vec<(usize, f64)> = freq_pairs.into_iter().collect();
        Self::new(modes, m_eff, spin, volume)
    }

    // -----------------------------------------------------------------------
    // Mode access
    // -----------------------------------------------------------------------

    /// Number of discrete modes stored.
    pub fn n_modes(&self) -> usize {
        self.modes.len()
    }

    // -----------------------------------------------------------------------
    // Physical observables
    // -----------------------------------------------------------------------

    /// Zero-point amplitude of mode `mode_idx`: `u_n = √(ℏ / (2 m_eff ω_n))` \[m\].
    ///
    /// This is the root-mean-square displacement amplitude of a quantum
    /// harmonic oscillator with effective mass `m_eff` at frequency `ω_n`
    /// in its ground state.
    ///
    /// # Arguments
    ///
    /// * `mode_idx` — index into `self.modes` (0-based position in the stored list).
    ///
    /// # Returns
    ///
    /// `√(ℏ / (2 m_eff ω_n))` \[m\].
    ///
    /// # Errors
    ///
    /// Returns an error if `mode_idx ≥ n_modes()`.
    pub fn zero_point_amplitude(&self, mode_idx: usize) -> Result<f64> {
        if mode_idx >= self.modes.len() {
            return Err(error::invalid_param(
                "mode_idx",
                &format!(
                    "index {mode_idx} out of range for {} modes",
                    self.modes.len()
                ),
            ));
        }
        let (_n, omega) = self.modes[mode_idx];
        Ok((HBAR / (2.0 * self.m_eff * omega)).sqrt())
    }

    /// Ground-state (vacuum) energy: `E_0 = (1/2) Σ_n ℏ ω_n` \[J\].
    ///
    /// The sum of zero-point energies over all stored modes — equivalent to
    /// the Casimir energy of the confined magnon gas at zero temperature.
    ///
    /// # Returns
    ///
    /// `(ℏ/2) Σ_n ω_n` \[J\].
    pub fn ground_state_energy(&self) -> f64 {
        let omega_sum: f64 = self.modes.iter().map(|&(_, omega)| omega).sum();
        0.5 * HBAR * omega_sum
    }

    /// Vacuum fluctuation density of the x-component of spin.
    ///
    /// Returns `Σ_n u_n² / V` where `u_n = √(ℏ/(2 m_eff ω_n))` is the
    /// zero-point amplitude of mode n and V = `self.volume`.  This has
    /// dimensions of m⁻¹ (fluctuation amplitude squared per volume), acting as
    /// a dimensionless fluctuation density when multiplied by a length scale.
    ///
    /// # Returns
    ///
    /// `Σ_n u_n² / volume` \[m⁻¹\].
    ///
    /// # Errors
    ///
    /// Returns an error if any mode amplitude computation fails.
    pub fn vacuum_fluctuation_sx(&self) -> Result<f64> {
        let mut sum = 0.0;
        for idx in 0..self.modes.len() {
            let un = self.zero_point_amplitude(idx)?;
            sum += un * un;
        }
        Ok(sum / self.volume)
    }

    /// Magnon free energy at temperature T \[J\].
    ///
    /// ```text
    /// F(T) = Σ_n ℏω_n [1/2 + n_BE(ω_n, T)]
    /// n_BE = 1 / [exp(ℏω_n / (k_B T)) − 1]
    /// ```
    ///
    /// At `T → 0`, `n_BE → 0` and `F → E_0 = (ℏ/2) Σ_n ω_n`.
    ///
    /// # Arguments
    ///
    /// * `temperature` — temperature T \[K\], must be ≥ 0.
    ///
    /// # Returns
    ///
    /// Free energy \[J\].
    ///
    /// # Errors
    ///
    /// Returns an error if `temperature < 0`.
    pub fn casimir_free_energy(&self, temperature: f64) -> Result<f64> {
        if temperature < 0.0 {
            return Err(error::invalid_param(
                "temperature",
                "temperature must be non-negative",
            ));
        }
        let mut f = 0.0;
        for &(_, omega) in &self.modes {
            if temperature < 1e-30 {
                // T → 0: Bose-Einstein = 0, only zero-point contribution
                f += 0.5 * HBAR * omega;
            } else {
                let x = HBAR * omega / (KB * temperature);
                // Numerically stable n_BE for x > 0
                let n_be = if x > 300.0 {
                    // Avoid overflow: n_BE ≈ exp(-x) → 0
                    (-x).exp()
                } else {
                    1.0 / (x.exp() - 1.0)
                };
                f += HBAR * omega * (0.5 + n_be);
            }
        }
        Ok(f)
    }

    /// Numerical Casimir-like magnon pressure from confinement length variation.
    ///
    /// In a stripe geometry with standing-wave quantisation, stretching the
    /// length by `dl` shifts each mode frequency as
    ///
    /// ```text
    /// ω_n(L + dL) ≈ ω_n(L) · L / (L + dL)
    /// ```
    ///
    /// (because `k_n = nπ/L`, so `ω ∝ k ∝ 1/L` in the exchange-dominated
    /// regime).  The finite-difference derivative gives a confinement force:
    ///
    /// ```text
    /// dF/dL ≈ [F(ω·L/(L+dL)) − F(ω)] / dL
    /// ```
    ///
    /// Note: this is evaluated at `T = 0` (ground-state Casimir pressure).
    ///
    /// # Arguments
    ///
    /// * `length` — current stripe length L \[m\], must be positive.
    /// * `dl`     — length increment \[m\], must be > 0.
    ///
    /// # Returns
    ///
    /// `dF/dL` \[J/m = N\], i.e. the confinement force (negative = attractive).
    ///
    /// # Errors
    ///
    /// Returns an error if `length ≤ 0` or `dl ≤ 0`.
    pub fn casimir_pressure_change(&self, length: f64, dl: f64) -> Result<f64> {
        if length <= 0.0 {
            return Err(error::invalid_param(
                "length",
                "confinement length must be positive",
            ));
        }
        if dl <= 0.0 {
            return Err(error::invalid_param(
                "dl",
                "length increment dl must be positive",
            ));
        }

        // Original ground state energy at L
        let f0 = self.ground_state_energy();

        // Shifted modes: ω_n → ω_n * L/(L+dl)
        let scale = length / (length + dl);
        let shifted_modes: Vec<(usize, f64)> = self
            .modes
            .iter()
            .map(|&(n, omega)| (n, omega * scale))
            .collect();

        // Ground-state energy at L + dl (T=0)
        let f1: f64 = shifted_modes
            .iter()
            .map(|&(_, omega)| 0.5 * HBAR * omega)
            .sum();

        // Finite-difference derivative dF/dL
        Ok((f1 - f0) / dl)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const ME_EFF: f64 = 9.109e-31; // electron mass as proxy [kg]
    const VOL: f64 = 1e-21; // 1 μm³ [m³]

    fn make_zpf(n: usize) -> ZeroPointFluctuations {
        // Construct n modes with frequencies ω_1, 2ω_1, ..., n*ω_1
        let omega0 = 1e11; // 100 GHz base frequency
        let modes: Vec<(usize, f64)> = (1..=n).map(|i| (i, i as f64 * omega0)).collect();
        ZeroPointFluctuations::new(modes, ME_EFF, 1.0, VOL).expect("valid params")
    }

    // ── amplitude positive ────────────────────────────────────────────────────

    #[test]
    fn test_amplitude_positive() {
        let zpf = make_zpf(3);
        for i in 0..3 {
            let u = zpf.zero_point_amplitude(i).expect("valid index");
            assert!(u > 0.0, "amplitude at mode {i} should be positive, got {u}");
        }
    }

    // ── amplitude scales as 1/sqrt(omega) ─────────────────────────────────────

    #[test]
    fn test_amplitude_scales_as_inverse_sqrt_omega() {
        let zpf = make_zpf(2);
        let u0 = zpf.zero_point_amplitude(0).expect("valid");
        let u1 = zpf.zero_point_amplitude(1).expect("valid");
        // u0 is at omega, u1 at 2*omega: ratio should be sqrt(2)
        let ratio = u0 / u1;
        let expected = 2.0_f64.sqrt();
        assert!(
            (ratio - expected).abs() < 1e-10,
            "u_0/u_1 should be sqrt(2), got {ratio}"
        );
    }

    // ── ground_state_energy positive ──────────────────────────────────────────

    #[test]
    fn test_ground_state_energy_positive() {
        let zpf = make_zpf(5);
        let e0 = zpf.ground_state_energy();
        assert!(e0 > 0.0, "ground state energy must be positive, got {e0}");
    }

    // ── vacuum_fluctuation positive ───────────────────────────────────────────

    #[test]
    fn test_vacuum_fluctuation_positive() {
        let zpf = make_zpf(4);
        let sx = zpf.vacuum_fluctuation_sx().expect("valid");
        assert!(sx > 0.0, "vacuum fluctuation Sx must be positive, got {sx}");
    }

    // ── Casimir free energy at T≈0 equals ground state energy ─────────────────

    #[test]
    fn test_casimir_at_zero_t_equals_ground_state() {
        let zpf = make_zpf(3);
        let e0 = zpf.ground_state_energy();
        // At T → 0: use T = 1e-50 K (effectively 0)
        let f_low = zpf.casimir_free_energy(1e-50).expect("valid");
        let rel_err = (f_low - e0).abs() / e0;
        assert!(
            rel_err < 1e-6,
            "Casimir F(T→0) should equal E_0, rel err={rel_err}"
        );
    }

    // ── Casimir free energy classical limit ───────────────────────────────────

    #[test]
    fn test_casimir_high_t_classical_limit() {
        // At high T: n_BE ≈ kT/(ℏω) → F ≈ Σ ℏω * kT/(ℏω) = N*kT
        let zpf = make_zpf(3);
        let t = 1e6; // extreme high T
        let f_high = zpf.casimir_free_energy(t).expect("valid");
        let n_modes = zpf.n_modes() as f64;
        let classical = n_modes * KB * t;
        // At such high T the thermal term dominates; ratio should approach 1
        let ratio = f_high / classical;
        assert!(
            (ratio - 1.0).abs() < 0.05,
            "high-T limit: F/(NkT) = {ratio}, expected ≈1"
        );
    }

    // ── from_quantized_modes round-trip ──────────────────────────────────────

    #[test]
    fn test_from_quantized_modes_round_trip_n_modes() {
        use crate::material::Ferromagnet;
        use crate::spinwave::{NanostructureGeometry, QuantizedModes};

        let py = Ferromagnet::permalloy();
        let geo = NanostructureGeometry::Stripe {
            width: 1e-6,
            length: 10e-6,
        };
        let qm = QuantizedModes::new(&py, geo, 50e-9).expect("valid");
        let n_req = 5;
        let zpf = ZeroPointFluctuations::from_quantized_modes(&qm, 0.05, n_req, 1.0, ME_EFF, VOL)
            .expect("valid zpf");
        assert_eq!(zpf.n_modes(), n_req, "round-trip must preserve mode count");
    }

    // ── parameter validation ──────────────────────────────────────────────────

    #[test]
    fn test_negative_m_eff_rejected() {
        let modes = vec![(1, 1e11)];
        assert!(ZeroPointFluctuations::new(modes, -1.0, 1.0, VOL).is_err());
    }

    #[test]
    fn test_negative_volume_rejected() {
        let modes = vec![(1, 1e11)];
        assert!(ZeroPointFluctuations::new(modes, ME_EFF, 1.0, -1.0).is_err());
    }

    // ── out-of-bounds mode index ──────────────────────────────────────────────

    #[test]
    fn test_invalid_mode_idx() {
        let zpf = make_zpf(2);
        assert!(zpf.zero_point_amplitude(5).is_err());
    }

    // ── casimir free energy increases with T ──────────────────────────────────

    #[test]
    fn test_casimir_free_energy_increases_with_t() {
        let zpf = make_zpf(3);
        let f_low = zpf.casimir_free_energy(1.0).expect("valid");
        let f_mid = zpf.casimir_free_energy(100.0).expect("valid");
        let f_high = zpf.casimir_free_energy(10000.0).expect("valid");
        assert!(f_mid > f_low, "F should increase with T");
        assert!(f_high > f_mid, "F should increase with T");
    }
}
