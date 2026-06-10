//! Exchange Bias Effect at FM/AFM Interfaces
//!
//! Exchange bias arises when a ferromagnet (FM) is exchange-coupled to an antiferromagnet (AFM)
//! across a sharp interface. After field-cooling through the AFM Néel temperature, the FM
//! hysteresis loop shifts along the applied-field axis — the hallmark of exchange bias.
//!
//! ## Physical Origin
//!
//! The interfacial exchange coupling pins some AFM spins antiparallel to the FM. On reversing
//! the applied field the FM must overcome not only its own anisotropy but also this interface
//! exchange energy, producing an asymmetric switching and a net loop shift H_EB.
//!
//! ## Meiklejohn–Bean Model (rigid AFM limit)
//!
//! Interface exchange energy density: J_EB [J/m²] (positive → prefers the configuration that
//! gives a left shift under our sign convention).
//!
//! Loop shift:   H_EB = -J_EB / (μ₀ M_FM t_FM)    [A/m]
//!
//! Enhanced coercive field contribution from bias pinning:
//!   ΔH_c = |J_EB| / (μ₀ M_FM t_FM)               [A/m]
//!
//! ## Key References
//!
//! - W. H. Meiklejohn, C. P. Bean, "New Magnetic Anisotropy", Phys. Rev. 102, 1413 (1956).
//! - J. Nogués, I. K. Schuller, "Exchange bias",
//!   J. Magn. Magn. Mater. **192**, 203–232 (1999).
//! - A. E. Berkowitz, K. Takano, "Exchange anisotropy — a review",
//!   J. Magn. Magn. Mater. **200**, 552–570 (1999).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::constants::MU_0;

// =============================================================================
// Material-specific interface exchange energy density constants
// =============================================================================

/// Interface exchange energy density for IrMn(7 nm)/Co bilayer [J/m²].
///
/// Representative value for room-temperature Co(~5 nm)/IrMn(7 nm) systems; the
/// exact value depends strongly on IrMn thickness, texture, and deposition conditions.
///
/// # References
/// Nogués & Schuller (1999), Table 1.
pub const J_EB_IRMN_CO: f64 = 0.12e-3; // 0.12 mJ/m²

/// Interface exchange energy density for CoO/Co bilayer [J/m²].
///
/// Classic system studied by Meiklejohn and Bean; strong AFM superexchange in CoO gives
/// a larger bias than IrMn at low temperature.
///
/// # References
/// Meiklejohn & Bean (1956); Nogués & Schuller (1999).
pub const J_EB_COO_CO: f64 = 0.20e-3; // 0.20 mJ/m²

/// Interface exchange energy density for FeMn/NiFe bilayer [J/m²].
///
/// The classic antiferromagnet FeMn (face-centered-cubic, equiatomic) is frequently
/// used with Permalloy (Ni₈₀Fe₂₀) for spin-valve spin-valve read heads.
///
/// # References
/// Nogués & Schuller (1999), Table 1.
pub const J_EB_FEMN_NIFE: f64 = 0.15e-3; // 0.15 mJ/m²

/// Interface exchange energy density for BiFeO₃/CoFe bilayer [J/m²].
///
/// Multiferroic BiFeO₃ as AFM; large room-temperature bias owing to the strong
/// G-type AFM order and electrically tunable anisotropy.
///
/// # References
/// Nogués & Schuller (1999); You et al., Phys. Rev. Lett. 2013.
pub const J_EB_BFCO_COFE: f64 = 0.40e-3; // 0.40 mJ/m²

// =============================================================================
// ExchangeBias struct
// =============================================================================

/// Exchange bias system: a ferromagnet / antiferromagnet bilayer after field cooling.
///
/// Implements the Meiklejohn–Bean rigid-AFM model extended with:
/// - Stoner–Wohlfarth-based hysteresis loop computation,
/// - Hoffmann training-effect formula,
/// - Power-law temperature dependence.
///
/// All fields are in SI units.
///
/// # References
/// - J. Nogués, I. K. Schuller, J. Magn. Magn. Mater. **192**, 203 (1999).
/// - A. Hoffmann, Phys. Rev. Lett. **93**, 097203 (2004).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ExchangeBias {
    /// Interface exchange energy density [J/m²].
    ///
    /// Positive value: the interface coupling stabilises the configuration that
    /// produces a *left* shift of the hysteresis loop (H_EB < 0 in standard convention).
    pub j_eb: f64,

    /// FM saturation magnetisation M_s [A/m].
    pub m_fm: f64,

    /// FM layer thickness t_FM \[m\].
    pub t_fm: f64,

    /// FM uniaxial anisotropy constant K_FM [J/m³].
    pub k_fm: f64,

    /// Blocking temperature T_B \[K\].
    ///
    /// Below T_B the loop shift is non-zero. T_B < T_N of the AFM in practice.
    pub t_b: f64,

    /// AFM Néel temperature T_N \[K\].
    pub t_n: f64,
}

impl ExchangeBias {
    /// Construct an [`ExchangeBias`] with all parameters specified explicitly.
    pub fn new(j_eb: f64, m_fm: f64, t_fm: f64, k_fm: f64, t_b: f64, t_n: f64) -> Self {
        Self {
            j_eb,
            m_fm,
            t_fm,
            k_fm,
            t_b,
            t_n,
        }
    }

    /// IrMn / Co bilayer preset.
    ///
    /// AFM layer assumed thick enough to fully set the bias (`t_afm ≳ 6 nm`).
    /// `t_co_nm` specifies the Co ferromagnet thickness in nanometres.
    ///
    /// Material parameters:
    /// - M_s(Co) ≈ 1.4 × 10⁶ A/m
    /// - K_FM(Co) ≈ 4.5 × 10⁵ J/m³  (bulk; reduced slightly at interface)
    /// - T_B(IrMn) ≈ 500 K, T_N(IrMn) ≈ 690 K
    ///
    /// # Arguments
    /// * `t_co_nm` - Co thickness in nanometres (typical 3–10 nm).
    pub fn irmn_co(t_co_nm: f64) -> Self {
        Self {
            j_eb: J_EB_IRMN_CO,
            m_fm: 1.4e6, // A/m — Co saturation magnetisation
            t_fm: t_co_nm * 1e-9,
            k_fm: 4.5e5, // J/m³ — Co uniaxial anisotropy
            t_b: 500.0,  // K   — IrMn blocking temp (Nogués 1999)
            t_n: 690.0,  // K   — IrMn Néel temp
        }
    }

    /// CoO / Co bilayer preset (Meiklejohn–Bean archetype).
    ///
    /// Uses a canonical 5 nm Co FM layer, representative of the original Meiklejohn & Bean
    /// particle-in-oxide geometry approximated as a continuous film.
    pub fn coo_co() -> Self {
        Self {
            j_eb: J_EB_COO_CO,
            m_fm: 1.4e6, // A/m — Co
            t_fm: 5e-9,  // 5 nm Co
            k_fm: 4.5e5, // J/m³ — Co
            t_b: 200.0,  // K   — CoO blocking temp
            t_n: 291.0,  // K   — CoO Néel temp (bulk)
        }
    }

    /// FeMn / NiFe (Permalloy) bilayer preset.
    ///
    /// Representative of 1990s spin-valve read heads. FeMn(10 nm)/NiFe(8 nm) stack.
    pub fn femn_nife() -> Self {
        Self {
            j_eb: J_EB_FEMN_NIFE,
            m_fm: 8.0e5, // A/m — NiFe saturation magnetisation
            t_fm: 8e-9,  // 8 nm NiFe
            k_fm: 2.5e3, // J/m³ — NiFe (very soft)
            t_b: 430.0,  // K   — FeMn blocking temp
            t_n: 490.0,  // K   — FeMn Néel temp
        }
    }

    // =========================================================================
    // Core physics: Meiklejohn–Bean model
    // =========================================================================

    /// Loop shift field H_EB [A/m].
    ///
    /// H_EB = -J_EB / (μ₀ M_FM t_FM)
    ///
    /// Negative return value means the loop shifts to *lower* (more negative) applied fields,
    /// which is the experimentally observed direction for a positive J_EB (interface prefers the
    /// FM magnetisation antiparallel to the field-cooled direction after reversal).
    pub fn loop_shift_field(&self) -> f64 {
        -self.j_eb / (MU_0 * self.m_fm * self.t_fm)
    }

    /// Additional coercive field contribution from the exchange-bias pinning [A/m].
    ///
    /// In the Meiklejohn–Bean model the switching asymmetry widens the loop by
    /// ΔH_c = |J_EB| / (μ₀ M_FM t_FM) on each side.
    pub fn coercive_enhancement(&self) -> f64 {
        self.j_eb.abs() / (MU_0 * self.m_fm * self.t_fm)
    }

    /// Anisotropy field H_K = 2 K_FM / (μ₀ M_FM) [A/m].
    ///
    /// Governs the intrinsic Stoner–Wohlfarth switching field in the absence of exchange bias.
    fn anisotropy_field(&self) -> f64 {
        2.0 * self.k_fm / (MU_0 * self.m_fm)
    }

    // =========================================================================
    // Training effect — Hoffmann (2004) model
    // =========================================================================

    /// Loop shift after `n` hysteresis cycles (training effect) [A/m].
    ///
    /// Follows the Hoffmann square-root decay:
    ///
    /// H_EB(n) = H_EB_∞ + (H_EB_1 - H_EB_∞) / √n
    ///
    /// The asymptotic value H_EB_∞ is taken as 80 % of the first-cycle value,
    /// consistent with the typical ~20 % total training decay measured in IrMn/Co.
    ///
    /// # Arguments
    /// * `n` - Cycle number (1-based; n = 1 returns the field-cooled loop shift).
    ///
    /// # Panics
    /// Panics if `n == 0` (cycle numbers start at 1).
    ///
    /// # References
    /// A. Hoffmann, Phys. Rev. Lett. **93**, 097203 (2004).
    pub fn training_field(&self, n: u32) -> f64 {
        assert!(n >= 1, "cycle number must be ≥ 1");
        let h_eb_1 = self.loop_shift_field();
        // Asymptotic (saturated) bias: 80 % of the first-cycle value (retained fraction).
        let h_eb_inf = 0.80 * h_eb_1;
        h_eb_inf + (h_eb_1 - h_eb_inf) / (n as f64).sqrt()
    }

    // =========================================================================
    // Temperature dependence
    // =========================================================================

    /// Temperature-dependent loop shift H_EB(T) [A/m].
    ///
    /// Power-law scaling above 0 K and below the blocking temperature:
    ///
    /// H_EB(T) = H_EB(0) × [1 - (T / T_B)]^(3/2)   for T < T_B
    /// H_EB(T) = 0                                     for T ≥ T_B
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]. Must be non-negative.
    pub fn loop_shift_at_temperature(&self, temperature: f64) -> f64 {
        if temperature >= self.t_b {
            return 0.0;
        }
        let t_reduced = temperature / self.t_b;
        let h_eb_0 = self.loop_shift_field();
        h_eb_0 * (1.0 - t_reduced).powf(1.5)
    }

    // =========================================================================
    // Hysteresis loop — Stoner–Wohlfarth + exchange bias
    // =========================================================================

    /// Compute the full hysteresis loop as (H [A/m], m/M_s) pairs.
    ///
    /// Uses the rigid-AFM Stoner–Wohlfarth model for the easy-axis geometry.
    /// The rectangular-loop approximation is exact for fields applied exactly
    /// along the uniaxial easy axis.
    ///
    /// - Descending branch (H: +H_max → -H_max): m = +1 until the FM switches at
    ///   H_sw^- = -(H_K + |H_EB|) for H_EB < 0, then m = -1.
    /// - Ascending branch  (H: -H_max → +H_max): m = -1 until the FM switches at
    ///   H_sw^+ = H_K - |H_EB| for H_EB < 0 (left shift reduces the positive
    ///   switching field).
    ///
    /// Returns 2 × `n_steps` points, first the descending branch then the ascending.
    ///
    /// # Arguments
    /// * `h_max`    - Maximum applied field magnitude [A/m]. Should exceed H_sw.
    /// * `n_steps`  - Number of field steps per half-sweep (≥ 2).
    pub fn compute_hysteresis_loop(&self, h_max: f64, n_steps: usize) -> Vec<(f64, f64)> {
        assert!(n_steps >= 2, "n_steps must be ≥ 2");

        let h_k = self.anisotropy_field();
        // Bias field (negative for positive J_EB → left shift).
        let h_eb_field = self.loop_shift_field();

        // Switching fields along the easy axis (Stoner–Wohlfarth + EB).
        // Descending (positive → negative sweep): switch at H_sw_desc.
        let h_sw_desc = -(h_k + h_eb_field.abs());
        // Ascending (negative → positive sweep): switch at H_sw_asc.
        // Left shift means the FM needs less field to re-saturate in +direction.
        let h_sw_asc = h_k - h_eb_field.abs();

        let mut loop_data: Vec<(f64, f64)> = Vec::with_capacity(2 * n_steps);

        // ── Descending branch: H from +h_max down to -h_max ──────────────────
        let mut m_desc = 1.0_f64; // start saturated positive
        for i in 0..n_steps {
            let h = h_max - (2.0 * h_max / (n_steps - 1) as f64) * i as f64;
            if h < h_sw_desc && m_desc > 0.0 {
                m_desc = -1.0; // irreversible switch to negative saturation
            }
            loop_data.push((h, m_desc));
        }

        // ── Ascending branch: H from -h_max up to +h_max ─────────────────────
        let mut m_asc = -1.0_f64; // start saturated negative
        for i in 0..n_steps {
            let h = -h_max + (2.0 * h_max / (n_steps - 1) as f64) * i as f64;
            if h > h_sw_asc && m_asc < 0.0 {
                m_asc = 1.0; // irreversible switch to positive saturation
            }
            loop_data.push((h, m_asc));
        }

        loop_data
    }

    /// Extract loop shift H_EB and coercive field H_c from a computed hysteresis loop.
    ///
    /// Scans the descending and ascending branches separately for the field at which the
    /// normalised magnetisation crosses zero (the coercive fields H_c^- and H_c^+).
    ///
    /// H_c  = (|H_c^+| + |H_c^-|) / 2   [average half-width]
    /// H_EB = (H_c^+ + H_c^-)    / 2     [centre of the loop — signed shift]
    ///
    /// # Arguments
    /// * `loop_data` - Slice of (H, m/M_s) pairs as returned by [`Self::compute_hysteresis_loop`].
    ///
    /// # Returns
    /// `(H_EB [A/m], H_c [A/m])`
    pub fn extract_loop_parameters(&self, loop_data: &[(f64, f64)]) -> (f64, f64) {
        let n = loop_data.len();
        if n < 4 {
            return (0.0, 0.0);
        }

        let half = n / 2;
        let descending = &loop_data[..half];
        let ascending = &loop_data[half..];

        let h_c_desc = Self::zero_crossing_field(descending);
        let h_c_asc = Self::zero_crossing_field(ascending);

        let h_eb = 0.5 * (h_c_asc + h_c_desc);
        let h_c = 0.5 * (h_c_asc.abs() + h_c_desc.abs());
        (h_eb, h_c)
    }

    /// Linear interpolation to find the H at which m/M_s crosses zero in a branch.
    fn zero_crossing_field(branch: &[(f64, f64)]) -> f64 {
        for window in branch.windows(2) {
            let (h0, m0) = window[0];
            let (h1, m1) = window[1];
            // Sign change: one point positive, one negative (or exactly zero).
            if m0 * m1 <= 0.0 && (m1 - m0).abs() > 1e-30 {
                // Linear interpolation.
                let t = -m0 / (m1 - m0);
                return h0 + t * (h1 - h0);
            }
        }
        // If no crossing found, return the midpoint H as a fallback.
        0.5 * (branch.first().map(|p| p.0).unwrap_or(0.0)
            + branch.last().map(|p| p.0).unwrap_or(0.0))
    }
}

// =============================================================================
// LoopShiftResult — analytical summary
// =============================================================================

/// Summary of exchange-bias loop parameters derived from the Meiklejohn–Bean model.
///
/// # References
/// - J. Nogués, I. K. Schuller, J. Magn. Magn. Mater. **192**, 203 (1999).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LoopShiftResult {
    /// Loop shift H_EB [A/m] (negative = left shift in standard convention).
    pub h_eb: f64,

    /// Coercive field H_c [A/m] (average of |H_c^+| and |H_c^-|, always positive).
    pub h_c: f64,

    /// Dimensionless ratio H_EB / H_c.
    ///
    /// Characterises how strongly the bias shifts the loop relative to its width.
    /// Values near ±1 indicate a strongly asymmetric loop (one switching field is
    /// near zero); values near 0 indicate the bias is negligible compared to the
    /// anisotropy.
    pub relative_loop_shift: f64,
}

impl LoopShiftResult {
    /// Compute all loop parameters analytically from an [`ExchangeBias`] system.
    ///
    /// Uses the Stoner–Wohlfarth easy-axis expressions:
    ///   H_EB = -J_EB / (μ₀ M_FM t_FM)
    ///   H_c  = 2 K_FM / (μ₀ M_FM)      (intrinsic anisotropy field = coercive field in SW model)
    pub fn from_exchange_bias(eb: &ExchangeBias) -> Self {
        let h_eb = eb.loop_shift_field();
        let h_c = eb.anisotropy_field();
        let relative_loop_shift = if h_c.abs() > 1e-30 { h_eb / h_c } else { 0.0 };
        Self {
            h_eb,
            h_c,
            relative_loop_shift,
        }
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    // ── Unit: loop_shift_field ────────────────────────────────────────────────

    #[test]
    fn test_loop_shift_field_sign() {
        // Positive J_EB must give a negative H_EB (left shift).
        let eb = ExchangeBias::irmn_co(5.0);
        let h_eb = eb.loop_shift_field();
        assert!(
            h_eb < 0.0,
            "H_EB should be negative for positive J_EB, got {h_eb}"
        );
    }

    #[test]
    fn test_loop_shift_field_magnitude() {
        // H_EB = -J_EB / (μ₀ M_s t_FM)
        // Use J_EB = 0.12e-3, M_s = 1.4e6, t_FM = 5e-9.
        let eb = ExchangeBias::irmn_co(5.0);
        let expected = -J_EB_IRMN_CO / (MU_0 * 1.4e6 * 5e-9);
        let got = eb.loop_shift_field();
        assert!(
            (got - expected).abs() < TOL * expected.abs(),
            "H_EB magnitude mismatch: got {got}, expected {expected}"
        );
    }

    #[test]
    fn test_loop_shift_scales_inversely_with_thickness() {
        // Doubling FM thickness halves H_EB.
        let eb5 = ExchangeBias::irmn_co(5.0);
        let eb10 = ExchangeBias::irmn_co(10.0);
        let ratio = eb5.loop_shift_field() / eb10.loop_shift_field();
        assert!(
            (ratio - 2.0).abs() < 1e-6,
            "H_EB should scale as 1/t_FM; ratio = {ratio}"
        );
    }

    // ── Unit: coercive_enhancement ───────────────────────────────────────────

    #[test]
    fn test_coercive_enhancement_positive() {
        let eb = ExchangeBias::coo_co();
        assert!(eb.coercive_enhancement() > 0.0);
    }

    #[test]
    fn test_coercive_enhancement_equals_abs_h_eb() {
        // ΔH_c = |J_EB| / (μ₀ M_s t_FM) = |H_EB|.
        let eb = ExchangeBias::femn_nife();
        let delta_hc = eb.coercive_enhancement();
        let h_eb_abs = eb.loop_shift_field().abs();
        assert!(
            (delta_hc - h_eb_abs).abs() < TOL * h_eb_abs,
            "coercive enhancement should equal |H_EB|: Δhc={delta_hc}, |H_EB|={h_eb_abs}"
        );
    }

    // ── Unit: training_field ─────────────────────────────────────────────────

    #[test]
    fn test_training_first_cycle() {
        // n = 1 must return exactly H_EB_1 (= loop_shift_field).
        let eb = ExchangeBias::irmn_co(5.0);
        let h_eb_1 = eb.loop_shift_field();
        let h_train_1 = eb.training_field(1);
        assert!(
            (h_train_1 - h_eb_1).abs() < TOL * h_eb_1.abs(),
            "training_field(1) should equal loop_shift_field(); got {h_train_1} vs {h_eb_1}"
        );
    }

    #[test]
    fn test_training_monotone_decay() {
        // |H_EB(n)| should decrease monotonically with n.
        let eb = ExchangeBias::irmn_co(5.0);
        let vals: Vec<f64> = (1u32..=10).map(|n| eb.training_field(n).abs()).collect();
        for w in vals.windows(2) {
            assert!(
                w[1] <= w[0] + 1e-30,
                "training effect should be monotonically decreasing: {vals:?}"
            );
        }
    }

    #[test]
    fn test_training_approaches_asymptote() {
        // training_field(large n) should approach 80 % of H_EB_1.
        let eb = ExchangeBias::irmn_co(5.0);
        let h_eb_1 = eb.loop_shift_field();
        let h_eb_inf = 0.80 * h_eb_1;
        let h_train_1000 = eb.training_field(1000);
        assert!(
            (h_train_1000 - h_eb_inf).abs() < 0.01 * h_eb_inf.abs(),
            "training_field(1000) should be ≈ 0.80×H_EB_1; got {h_train_1000} vs {h_eb_inf}"
        );
    }

    // ── Unit: loop_shift_at_temperature ─────────────────────────────────────

    #[test]
    fn test_temperature_zero() {
        // At T = 0 K the shift must equal the full H_EB_0.
        let eb = ExchangeBias::irmn_co(5.0);
        let h_0 = eb.loop_shift_field();
        let h_at_0 = eb.loop_shift_at_temperature(0.0);
        assert!(
            (h_at_0 - h_0).abs() < TOL * h_0.abs(),
            "T=0 should give full H_EB(0); got {h_at_0} vs {h_0}"
        );
    }

    #[test]
    fn test_temperature_above_blocking() {
        // At T ≥ T_B the bias vanishes.
        let eb = ExchangeBias::irmn_co(5.0);
        assert_eq!(eb.loop_shift_at_temperature(eb.t_b), 0.0);
        assert_eq!(eb.loop_shift_at_temperature(eb.t_b + 100.0), 0.0);
    }

    #[test]
    fn test_temperature_monotone_decrease() {
        // |H_EB(T)| should decrease monotonically as T increases toward T_B.
        let eb = ExchangeBias::irmn_co(5.0);
        let temps: Vec<f64> = (0..10).map(|i| i as f64 * eb.t_b / 10.0).collect();
        let h_vals: Vec<f64> = temps
            .iter()
            .map(|&t| eb.loop_shift_at_temperature(t).abs())
            .collect();
        for w in h_vals.windows(2) {
            assert!(
                w[1] <= w[0] + 1e-30,
                "|H_EB(T)| should decrease with T: {h_vals:?}"
            );
        }
    }

    #[test]
    fn test_temperature_power_law_exponent() {
        // Verify the 3/2 power-law numerically at T = T_B / 2.
        let eb = ExchangeBias::irmn_co(5.0);
        let h_0 = eb.loop_shift_field();
        let h_half = eb.loop_shift_at_temperature(eb.t_b / 2.0);
        // expected = h_0 × (1 - 0.5)^(3/2) = h_0 × 0.5^1.5 ≈ h_0 × 0.35355.
        let expected = h_0 * (0.5_f64).powf(1.5);
        assert!(
            (h_half - expected).abs() < 1e-10 * expected.abs(),
            "T=T_B/2 power law: got {h_half}, expected {expected}"
        );
    }

    // ── Unit: hysteresis loop ─────────────────────────────────────────────────

    #[test]
    fn test_hysteresis_loop_length() {
        let eb = ExchangeBias::irmn_co(5.0);
        let h_max = 3.0e6_f64;
        let n_steps = 200_usize;
        let loop_data = eb.compute_hysteresis_loop(h_max, n_steps);
        assert_eq!(loop_data.len(), 2 * n_steps);
    }

    #[test]
    fn test_hysteresis_loop_saturation_ends() {
        // At the extreme positive field the magnetisation should be +1.
        let eb = ExchangeBias::irmn_co(5.0);
        let h_max = 5.0e6_f64;
        let loop_data = eb.compute_hysteresis_loop(h_max, 500);
        // First point of descending branch: H = +h_max → m = +1.
        assert!(
            (loop_data[0].1 - 1.0).abs() < 1e-9,
            "First point m should be +1; got {}",
            loop_data[0].1
        );
        // First point of ascending branch: H = -h_max → m = -1.
        let asc_start = loop_data[500].1;
        assert!(
            (asc_start + 1.0).abs() < 1e-9,
            "Ascending start m should be -1; got {asc_start}"
        );
    }

    #[test]
    fn test_hysteresis_loop_shift_direction() {
        // The extracted H_EB from the loop should be negative (left shift).
        let eb = ExchangeBias::irmn_co(5.0);
        let h_max = 5.0e6_f64;
        let loop_data = eb.compute_hysteresis_loop(h_max, 1000);
        let (h_eb_extracted, _h_c) = eb.extract_loop_parameters(&loop_data);
        assert!(
            h_eb_extracted < 0.0,
            "Extracted H_EB should be negative (left shift); got {h_eb_extracted}"
        );
    }

    #[test]
    fn test_hysteresis_extract_consistency_with_analytical() {
        // The extracted H_EB should agree with the analytical loop_shift_field()
        // to within 15 %. The discretisation error is bounded by ½ × step ≈ step/2,
        // and since H_EB ≪ H_K, step/2 is an 8–10 % error at n_steps=2000.
        // Using 15 % gives comfortable margin over the worst-case half-step error.
        let eb = ExchangeBias::irmn_co(5.0);
        let h_max = 5.0e6_f64;
        let n_steps = 2000_usize;
        let loop_data = eb.compute_hysteresis_loop(h_max, n_steps);
        let (h_eb_extracted, _h_c) = eb.extract_loop_parameters(&loop_data);
        let h_eb_analytical = eb.loop_shift_field();
        let rel_err = (h_eb_extracted - h_eb_analytical).abs() / h_eb_analytical.abs();
        assert!(
            rel_err < 0.15,
            "Extracted H_EB should match analytical to 15% (discretisation bound); rel_err = {rel_err:.4}"
        );
    }

    // ── Unit: LoopShiftResult ────────────────────────────────────────────────

    #[test]
    fn test_loop_shift_result_from_exchange_bias() {
        let eb = ExchangeBias::irmn_co(5.0);
        let result = LoopShiftResult::from_exchange_bias(&eb);
        assert!(result.h_eb < 0.0, "H_EB should be negative");
        assert!(result.h_c > 0.0, "H_c should be positive");
        assert!(
            result.relative_loop_shift < 0.0,
            "relative loop shift should be negative"
        );
        // |relative shift| should be << 1 for typical parameters (bias << anisotropy field in IrMn/Co).
        assert!(
            result.relative_loop_shift.abs() < 1.0,
            "relative loop shift magnitude should be < 1; got {}",
            result.relative_loop_shift
        );
    }

    #[test]
    fn test_loop_shift_result_relative_shift_formula() {
        let eb = ExchangeBias::irmn_co(5.0);
        let result = LoopShiftResult::from_exchange_bias(&eb);
        let expected_ratio = eb.loop_shift_field() / eb.anisotropy_field();
        assert!(
            (result.relative_loop_shift - expected_ratio).abs() < 1e-12,
            "relative_loop_shift mismatch: got {}, expected {expected_ratio}",
            result.relative_loop_shift
        );
    }

    // ── Material presets sanity ───────────────────────────────────────────────

    #[test]
    fn test_presets_physical_range() {
        // All presets should give H_EB in the range 1 kA/m – 10 MA/m for typical thicknesses.
        for eb in [
            ExchangeBias::coo_co(),
            ExchangeBias::femn_nife(),
            ExchangeBias::irmn_co(5.0),
        ] {
            let h_eb_abs = eb.loop_shift_field().abs();
            assert!(
                h_eb_abs > 1e3 && h_eb_abs < 1e7,
                "H_EB = {h_eb_abs} A/m is outside plausible range [1 kA/m, 10 MA/m]"
            );
        }
    }

    #[test]
    fn test_neel_above_blocking_for_all_presets() {
        for eb in [
            ExchangeBias::coo_co(),
            ExchangeBias::femn_nife(),
            ExchangeBias::irmn_co(5.0),
        ] {
            assert!(
                eb.t_n > eb.t_b,
                "T_N must exceed T_B: T_N={}, T_B={}",
                eb.t_n,
                eb.t_b
            );
        }
    }
}
