//! Tavis-Cummings / Dicke Model — Collective Spin-Cavity Coupling
//!
//! Models N spin-1/2 emitters (magnons, qubits, etc.) collectively coupled to a
//! single cavity mode. In the Holstein-Primakoff / low-excitation limit the collective
//! spin behaves as a bosonic mode b with coupling g_N = g·√N, leading to a
//! pair of coupled harmonic oscillators amenable to analytic treatment.
//!
//! The mean-field equations of motion (rotating-frame, dissipative) are:
//!
//!   da/dt = -iω_c a - (κ/2) a - i g_N b + f(t)
//!   db/dt = -iω_m b - (γ/2) b - i g_N a
//!
//! where a is the cavity photon amplitude and b is the collective magnon amplitude.
//!
//! # References
//!
//! - M. Tavis, F. W. Cummings, Phys. Rev. **170**, 379 (1968)
//! - R. H. Dicke, Phys. Rev. **93**, 99 (1954)
//! - K. Hepp, E. H. Lieb, Ann. Phys. **76**, 360 (1973) — superradiant transition
//! - Y. Zhang et al., Phys. Rev. Lett. **113**, 156401 (2014) — magnon-polariton strong coupling

use std::f64::consts::PI;

use crate::constants::{GAMMA, MU_0};
use crate::error::{self, Result};
use crate::math::Complex;

/// Tavis-Cummings model: N spin-1/2 emitters coupled to a single cavity mode.
///
/// In the Holstein-Primakoff / low-excitation limit the N spins are treated as
/// a single collective bosonic magnon mode b with enhanced coupling g_N = g·√N.
///
/// The system is characterised by:
/// - A photon mode a at frequency ω_c with decay rate κ.
/// - A collective magnon mode b at frequency ω_m with damping rate γ.
/// - Collective coupling g_N = g·√N between the two modes.
///
/// # Examples
///
/// ```rust
/// use spintronics::cavity::TavisCummings;
///
/// let tc = TavisCummings::yig_ensemble(100);
/// let (omega_lo, omega_hi) = tc.polariton_frequencies();
/// assert!(omega_hi > omega_lo);
/// ```
#[derive(Debug, Clone)]
pub struct TavisCummings {
    /// Number of spins in the ensemble.
    pub n_spins: usize,
    /// Cavity angular frequency ω_c [rad/s].
    pub omega_cavity: f64,
    /// Magnon (spin) angular frequency ω_m [rad/s].
    pub omega_magnon: f64,
    /// Single-spin coupling constant g [rad/s].
    pub coupling_g: f64,
    /// Cavity field decay rate κ [rad/s].
    pub kappa: f64,
    /// Spin damping rate γ [rad/s].
    pub gamma_damp: f64,
    /// Internal: cavity field complex amplitude a(t).
    cavity_amp: Complex,
    /// Internal: collective magnon complex amplitude b(t).
    magnon_amp: Complex,
}

impl TavisCummings {
    /// Create a new Tavis-Cummings system.
    ///
    /// # Errors
    ///
    /// Returns `Err` if `n_spins` is zero, or if any frequency / rate is negative.
    pub fn new(
        n_spins: usize,
        omega_cavity: f64,
        omega_magnon: f64,
        coupling_g: f64,
        kappa: f64,
        gamma_damp: f64,
    ) -> Result<Self> {
        if n_spins == 0 {
            return Err(error::invalid_param(
                "n_spins",
                "ensemble must contain at least one spin",
            ));
        }
        if omega_cavity < 0.0 {
            return Err(error::invalid_param(
                "omega_cavity",
                "cavity frequency must be non-negative",
            ));
        }
        if omega_magnon < 0.0 {
            return Err(error::invalid_param(
                "omega_magnon",
                "magnon frequency must be non-negative",
            ));
        }
        if coupling_g < 0.0 {
            return Err(error::invalid_param(
                "coupling_g",
                "coupling constant must be non-negative",
            ));
        }
        if kappa < 0.0 {
            return Err(error::invalid_param(
                "kappa",
                "cavity decay rate must be non-negative",
            ));
        }
        if gamma_damp < 0.0 {
            return Err(error::invalid_param(
                "gamma_damp",
                "spin damping rate must be non-negative",
            ));
        }
        Ok(Self {
            n_spins,
            omega_cavity,
            omega_magnon,
            coupling_g,
            kappa,
            gamma_damp,
            cavity_amp: Complex::ZERO,
            magnon_amp: Complex::ZERO,
        })
    }

    /// YIG ensemble preset matching the `HybridSystem::yig_cavity()` parameters.
    ///
    /// Uses a **fixed per-spin coupling** g = 2π × 100 MHz (same for all N), so that
    /// the collective coupling g_N = g·√N scales naturally with the ensemble size:
    ///
    /// - ω_c = 2π × 10 GHz
    /// - ω_m = GAMMA × μ₀ × H_bias  (H_bias = 80 kA/m)
    /// - g (per-spin) = 2π × 100 MHz  (single-spin coupling as in `HybridSystem::yig_cavity()`)
    /// - κ = 2π × 1 MHz
    /// - γ = ω_m × α  (α = 0.0001, ultra-low YIG damping)
    pub fn yig_ensemble(n_spins: usize) -> Self {
        let n = n_spins.max(1);
        let omega_cavity = 2.0 * PI * 10.0e9;
        let h_bias = 8.0e4_f64; // 80 kA/m
        let omega_magnon = GAMMA * MU_0 * h_bias;
        // Fixed per-spin coupling g, matching HybridSystem::yig_cavity().
        // The collective coupling then grows as g_N = g * sqrt(N).
        let coupling_g = 2.0 * PI * 100.0e6;
        let kappa = 2.0 * PI * 1.0e6;
        let alpha = 0.0001_f64;
        let gamma_damp = alpha * omega_magnon;
        Self {
            n_spins: n,
            omega_cavity,
            omega_magnon,
            coupling_g,
            kappa,
            gamma_damp,
            cavity_amp: Complex::ZERO,
            magnon_amp: Complex::ZERO,
        }
    }

    /// Collective coupling g_N = g·√N [rad/s].
    ///
    /// This is the effective magnon-photon coupling after bosonic enhancement by
    /// the square root of the spin ensemble size.
    pub fn collective_coupling(&self) -> f64 {
        self.coupling_g * (self.n_spins as f64).sqrt()
    }

    /// Polariton frequencies from 2×2 matrix diagonalisation.
    ///
    /// The non-Hermitian matrix eigenfrequencies of the lossless dressed system are:
    ///
    ///   ω_± = (ω_c + ω_m)/2 ± √[(Δ/2)² + g_N²]
    ///
    /// where Δ = ω_c − ω_m is the detuning.
    ///
    /// Returns `(omega_lower, omega_upper)`.
    pub fn polariton_frequencies(&self) -> (f64, f64) {
        let g_n = self.collective_coupling();
        let delta = self.omega_cavity - self.omega_magnon;
        let discriminant = ((delta / 2.0).powi(2) + g_n.powi(2)).sqrt();
        let center = (self.omega_cavity + self.omega_magnon) / 2.0;
        (center - discriminant, center + discriminant)
    }

    /// Superradiant phase-transition critical coupling (Hepp-Lieb-Wang threshold).
    ///
    /// g_c = √(ω_c · ω_m) / 2
    ///
    /// The Dicke superradiant phase transition occurs when g_N > g_c.
    pub fn superradiant_threshold(&self) -> f64 {
        (self.omega_cavity * self.omega_magnon).sqrt() / 2.0
    }

    /// Cooperativity C = 4 g_N² / (κ · γ).
    ///
    /// The strong-coupling regime corresponds to C >> 1.
    pub fn cooperativity(&self) -> f64 {
        let g_n = self.collective_coupling();
        4.0 * g_n.powi(2) / (self.kappa * self.gamma_damp)
    }

    /// Frequency detuning Δ = ω_c − ω_m [rad/s].
    pub fn detuning(&self) -> f64 {
        self.omega_cavity - self.omega_magnon
    }

    /// Current cavity field amplitude a(t) (complex).
    pub fn cavity_amplitude(&self) -> Complex {
        self.cavity_amp
    }

    /// Current collective magnon amplitude b(t) (complex).
    pub fn magnon_amplitude(&self) -> Complex {
        self.magnon_amp
    }

    /// Advance the mean-field equations of motion by one Euler step of size `dt` [s].
    ///
    /// Equations (lossless rotating-frame with dissipation):
    ///
    ///   da/dt = −iω_c a − (κ/2) a − i g_N b + drive
    ///   db/dt = −iω_m b − (γ/2) b − i g_N a
    ///
    /// Uses the explicit Euler method; for accurate trajectories choose dt << 1/ω_c.
    pub fn mean_field_evolve(&mut self, dt: f64, drive: Complex) {
        let g_n = self.collective_coupling();
        let a = self.cavity_amp;
        let b = self.magnon_amp;

        // da/dt = -i ω_c a  - (κ/2) a  - i g_N b  + drive
        let da = a
            .mul_i()
            .scale(-self.omega_cavity)
            .sub(&a.scale(self.kappa / 2.0))
            .sub(&b.mul_i().scale(g_n))
            .add(&drive);

        // db/dt = -i ω_m b  - (γ/2) b  - i g_N a
        let db = b
            .mul_i()
            .scale(-self.omega_magnon)
            .sub(&b.scale(self.gamma_damp / 2.0))
            .sub(&a.mul_i().scale(g_n));

        self.cavity_amp = a.add(&da.scale(dt));
        self.magnon_amp = b.add(&db.scale(dt));
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Relative tolerance helper.
    fn rel_close(a: f64, b: f64, tol: f64) -> bool {
        if b == 0.0 {
            return a.abs() < tol;
        }
        ((a - b) / b).abs() < tol
    }

    // 1. Collective coupling scales as √N.
    #[test]
    fn test_collective_coupling_scales_sqrt_n() {
        let n: usize = 100;
        let tc_n = TavisCummings::yig_ensemble(n);
        let tc_4n = TavisCummings::yig_ensemble(4 * n);
        let ratio = tc_4n.collective_coupling() / tc_n.collective_coupling();
        // g_N ∝ √N  →  g_N(4N)/g_N(N) = √(4N)/√N = 2
        assert!((ratio - 2.0).abs() < 1e-10, "ratio = {ratio}, expected 2.0");
    }

    // 2. At resonance (Δ=0) the polariton splitting equals 2·g_N.
    #[test]
    fn test_polariton_anticrossing_at_resonance() {
        let omega = 2.0 * PI * 10.0e9;
        let g = 2.0 * PI * 100.0e6;
        let kappa = 2.0 * PI * 1.0e6;
        let gamma = 2.0 * PI * 1.0e6;
        let tc = TavisCummings::new(1, omega, omega, g, kappa, gamma).unwrap();
        let (lo, hi) = tc.polariton_frequencies();
        let g_n = tc.collective_coupling();
        assert!(
            rel_close(hi - lo, 2.0 * g_n, 1e-12),
            "splitting = {}, 2g_N = {}",
            hi - lo,
            2.0 * g_n
        );
    }

    // 3. Superradiant threshold is strictly positive for physical parameters.
    #[test]
    fn test_superradiant_threshold_positive() {
        let tc = TavisCummings::yig_ensemble(50);
        let g_c = tc.superradiant_threshold();
        assert!(g_c > 0.0, "threshold must be positive, got {g_c}");
        // g_c = sqrt(ω_c · ω_m) / 2 — verify formula directly
        let expected = (tc.omega_cavity * tc.omega_magnon).sqrt() / 2.0;
        assert!(
            rel_close(g_c, expected, 1e-14),
            "threshold mismatch: {} vs {}",
            g_c,
            expected
        );
    }

    // 4. Cooperativity matches analytic formula 4 g_N² / (κ γ).
    #[test]
    fn test_cooperativity_formula() {
        let omega_c = 2.0 * PI * 10.0e9;
        let omega_m = 2.0 * PI * 9.0e9;
        let g = 2.0 * PI * 50.0e6;
        let kappa = 2.0 * PI * 2.0e6;
        let gamma = 2.0 * PI * 1.0e6;
        let n: usize = 16;
        let tc = TavisCummings::new(n, omega_c, omega_m, g, kappa, gamma).unwrap();
        let g_n = g * (n as f64).sqrt();
        let expected_c = 4.0 * g_n.powi(2) / (kappa * gamma);
        assert!(
            rel_close(tc.cooperativity(), expected_c, 1e-12),
            "C = {}, expected = {}",
            tc.cooperativity(),
            expected_c
        );
    }

    // 5. After 1000 zero-drive steps the amplitudes remain finite.
    #[test]
    fn test_mean_field_evolve_finite() {
        let mut tc = TavisCummings::yig_ensemble(100);
        // Seed a non-zero initial cavity amplitude.
        tc.cavity_amp = Complex::new(1.0e-3, 0.0);
        let dt = 1.0e-13_f64; // 0.1 ps
        for _ in 0..1000 {
            tc.mean_field_evolve(dt, Complex::ZERO);
        }
        assert!(
            tc.cavity_amplitude().is_finite(),
            "cavity amplitude became non-finite"
        );
        assert!(
            tc.magnon_amplitude().is_finite(),
            "magnon amplitude became non-finite"
        );
    }

    // 6. n_spins=1 cooperativity is the same order of magnitude as HybridSystem.
    #[test]
    fn test_yig_ensemble_single_spin_matches_hybrid() {
        let tc = TavisCummings::yig_ensemble(1);
        let c = tc.cooperativity();
        // The HybridSystem::yig_cavity() has C >> 10.
        // With n=1 and the per-spin g chosen so that g_N = 2π·100 MHz,
        // C should be large (>> 1 in strong-coupling regime).
        assert!(
            c > 1.0,
            "expected cooperativity > 1 for single-spin YIG, got {c}"
        );
    }

    // 7. Constructing with n_spins=0 must return an error.
    #[test]
    fn test_n_spins_zero_errors() {
        let result = TavisCummings::new(0, 1.0e10, 1.0e10, 1.0e8, 1.0e6, 1.0e5);
        assert!(result.is_err(), "expected Err for n_spins=0");
    }

    // 8. Large detuning: polariton frequencies approach ω_c and ω_m.
    #[test]
    fn test_large_detuning_polaritons() {
        // Δ = 100 · g_N  →  ω_± ≈ ω_c, ω_m
        let g = 2.0 * PI * 10.0e6;
        let omega_c = 2.0 * PI * 10.0e9;
        let omega_m = 2.0 * PI * 8.0e9; // Δ = 2π · 2 GHz >> g
        let kappa = 2.0 * PI * 1.0e6;
        let gamma = 2.0 * PI * 1.0e5;
        let tc = TavisCummings::new(1, omega_c, omega_m, g, kappa, gamma).unwrap();
        let (lo, hi) = tc.polariton_frequencies();
        // Lower polariton should be close to the smaller bare frequency (ω_m),
        // upper should be close to ω_c.  Tolerance: 5% of the gap.
        let gap = omega_c - omega_m;
        assert!(
            (lo - omega_m).abs() < 0.05 * gap,
            "lo polariton {lo:.3e} not near omega_m {omega_m:.3e}"
        );
        assert!(
            (hi - omega_c).abs() < 0.05 * gap,
            "hi polariton {hi:.3e} not near omega_c {omega_c:.3e}"
        );
    }

    // 9. Detuning sign: positive when ω_c > ω_m.
    #[test]
    fn test_detuning_sign() {
        let omega_c = 2.0 * PI * 10.0e9;
        let omega_m = 2.0 * PI * 8.0e9;
        let tc = TavisCummings::new(1, omega_c, omega_m, 1.0e8, 1.0e6, 1.0e5).unwrap();
        assert!(
            tc.detuning() > 0.0,
            "detuning should be positive when omega_c > omega_m"
        );
    }

    // 10. A non-zero drive increases the cavity amplitude from zero.
    #[test]
    fn test_mean_field_drive_increases_amplitude() {
        let mut tc = TavisCummings::yig_ensemble(100);
        let drive = Complex::new(1.0e6, 0.0); // strong drive
        let dt = 1.0e-13_f64;
        let initial_norm = tc.cavity_amplitude().norm();
        for _ in 0..200 {
            tc.mean_field_evolve(dt, drive);
        }
        let final_norm = tc.cavity_amplitude().norm();
        assert!(
            final_norm > initial_norm,
            "drive should increase cavity amplitude: initial={initial_norm}, final={final_norm}"
        );
    }
}
