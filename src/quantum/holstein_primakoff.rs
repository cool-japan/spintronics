//! Holstein-Primakoff transformation for spin-S operators.
//!
//! The Holstein-Primakoff (HP) transformation maps spin operators to bosonic
//! creation and annihilation operators in a systematic 1/S expansion.  For a
//! spin-S site with quantisation axis along z the exact identities are:
//!
//! ```text
//! S^z  = S - a†a
//! S^+  = √(2S - a†a) · a  ≈ √(2S) a  (linear order)
//! S^-  = a† · √(2S - a†a) ≈ √(2S) a†  (linear order)
//! ```
//!
//! where `a†`, `a` are canonical bosonic operators satisfying `[a, a†] = 1`.
//!
//! The linear (spin-wave) approximation is valid when the average boson occupation
//! `<n> = <a†a>` is much smaller than 2S.  The next correction enters at order 1/S
//! and shifts the effective excitation energy.
//!
//! # References
//!
//! T. Holstein and H. Primakoff, *Phys. Rev.* **58**, 1098 (1940).

use crate::error::{self, Result};

// ---------------------------------------------------------------------------
// Expansion order
// ---------------------------------------------------------------------------

/// Order of the 1/S expansion in the Holstein-Primakoff transformation.
///
/// `Linear` keeps only the leading-order term `S^± ≈ √(2S) a^(†)`.  `Quadratic`
/// includes the O(1/S) correction that accounts for the finite bosonic space.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum HpOrder {
    /// Leading-order (spin-wave) approximation: `S^± ≈ √(2S) a^(†)`.
    Linear,
    /// Next-to-leading-order: includes `−n(n−1)/(4S)` correction.
    Quadratic,
}

// ---------------------------------------------------------------------------
// Main struct
// ---------------------------------------------------------------------------

/// Holstein-Primakoff transformation parameters for a spin-S lattice site.
///
/// Stores the spin quantum number S and the desired expansion order.  All
/// coefficient methods are analytic in S, so no matrix operations are needed.
///
/// # Validity
///
/// The HP expansion is an asymptotic series in powers of 1/(2S).  Results are
/// physically meaningful only when the fractional boson occupation `n/(2S)` is
/// small (typically < 0.1 for quantitative accuracy).
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HolsteinPrimakoff {
    /// Spin quantum number S (positive half-integer: 1/2, 1, 3/2, …).
    pub spin: f64,
    /// Order of the 1/S expansion used for S^± coefficients.
    pub expansion_order: HpOrder,
}

impl HolsteinPrimakoff {
    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    /// Create a Holstein-Primakoff transform for spin `spin`.
    ///
    /// # Arguments
    ///
    /// * `spin` — spin quantum number S \[dimensionless\].  Must be a positive
    ///   multiple of 1/2 within tolerance 10⁻⁹.
    ///
    /// # Errors
    ///
    /// Returns [`crate::error::Error::InvalidParameter`] if `spin ≤ 0` or if
    /// `spin` is not a half-integer.
    pub fn new(spin: f64) -> Result<Self> {
        if spin <= 0.0 {
            return Err(error::invalid_param(
                "spin",
                "spin quantum number must be positive",
            ));
        }
        // Check half-integer: 2S must be an integer within 1e-9
        let twice_spin = 2.0 * spin;
        let rounded = twice_spin.round();
        if (twice_spin - rounded).abs() > 1e-9 {
            return Err(error::invalid_param(
                "spin",
                "spin must be a positive multiple of 1/2 (half-integer or integer)",
            ));
        }
        Ok(Self {
            spin,
            expansion_order: HpOrder::Linear,
        })
    }

    /// Set the expansion order, returning `self` for builder-style chaining.
    ///
    /// # Arguments
    ///
    /// * `order` — [`HpOrder::Linear`] or [`HpOrder::Quadratic`].
    pub fn with_order(self, order: HpOrder) -> Self {
        Self {
            expansion_order: order,
            ..self
        }
    }

    /// Ferromagnetic default: spin-1, linear expansion (spin-wave approximation).
    ///
    /// Suitable for simple ferromagnets such as Fe, Ni in which S ≈ 1 is a
    /// common effective value.
    pub fn ferromagnet_default() -> Self {
        Self {
            spin: 1.0,
            expansion_order: HpOrder::Linear,
        }
    }

    /// Antiferromagnetic default: spin-1/2, quadratic expansion.
    ///
    /// Spin-1/2 is the natural choice for localized-moment antiferromagnets
    /// (Cu²⁺, Mn²⁺ in CuGeO₃, MnO, etc.).  The quadratic order includes the
    /// leading quantum correction beyond linear spin-wave theory.
    pub fn antiferromagnet_default() -> Self {
        Self {
            spin: 0.5,
            expansion_order: HpOrder::Quadratic,
        }
    }

    // -----------------------------------------------------------------------
    // Physical observables
    // -----------------------------------------------------------------------

    /// Expectation value of S^z for boson occupation `n_boson`.
    ///
    /// Exact in the HP representation: `⟨S^z⟩ = S − n`.
    ///
    /// # Arguments
    ///
    /// * `n_boson` — boson occupation number n \[dimensionless, ≥ 0\].
    ///
    /// # Returns
    ///
    /// `S − n_boson` \[dimensionless\].
    #[inline]
    pub fn transform_sz(&self, n_boson: f64) -> f64 {
        self.spin - n_boson
    }

    /// Leading-order coefficient of S^+ in the HP expansion: √(2S).
    ///
    /// In the linear approximation `S^+ ≈ c · a` where `c = √(2S)`.
    /// This coefficient appears as the magnon amplitude in spin-wave theory.
    ///
    /// # Returns
    ///
    /// `√(2S)` \[dimensionless\].
    #[inline]
    pub fn transform_splus_coeff(&self) -> f64 {
        (2.0 * self.spin).sqrt()
    }

    /// Leading-order coefficient of S^− in the HP expansion: √(2S).
    ///
    /// By Hermitian conjugation, `S^− ≈ c · a†` with the same coefficient
    /// `c = √(2S)` as for `S^+`.
    ///
    /// # Returns
    ///
    /// `√(2S)` \[dimensionless\].
    #[inline]
    pub fn transform_sminus_coeff(&self) -> f64 {
        (2.0 * self.spin).sqrt()
    }

    /// Quadratic correction to the S^+ S^− product from the HP expansion.
    ///
    /// At next-to-leading order the exact HP factor expands as
    ///
    /// ```text
    /// √(2S − n) ≈ √(2S) [1 − n/(4S) − n(n−1)/(16S²) − …]
    /// ```
    ///
    /// The leading correction to `S^+S^-` (i.e., to the magnon number operator)
    /// is `−n(n−1)/(4S)`.  This is relevant in dense-magnon regimes.
    ///
    /// # Arguments
    ///
    /// * `n_boson` — boson occupation n \[dimensionless\].
    ///
    /// # Returns
    ///
    /// Correction `−n(n−1)/(4S)` \[dimensionless\].
    pub fn quadratic_correction(&self, n_boson: f64) -> f64 {
        -n_boson * (n_boson - 1.0) / (4.0 * self.spin)
    }

    /// Check the commutator identity \[S^+, S^-\] = 2 S^z to leading HP order.
    ///
    /// At linear order `[S^+, S^-] = 2S − 2a†a = 2S^z`.  The fractional
    /// deviation `|([S^+,S^-] − 2S^z)/(2S)|` is computed at two boson
    /// occupation values `n1` (for S^z in the commutator) and `n2` (for
    /// the explicit S^z reference).
    ///
    /// For very low occupation (`n1, n2 ≪ 2S`) the commutator is satisfied
    /// to within the linear approximation.  At high occupation (n ≈ 2S) the
    /// higher-order HP terms become important and the linear approximation fails.
    ///
    /// # Arguments
    ///
    /// * `n1` — boson occupation entering the commutator evaluation.
    /// * `n2` — boson occupation entering the explicit S^z evaluation.
    /// * `tol` — tolerance for the fractional deviation.
    ///
    /// # Returns
    ///
    /// `true` if `|deviation| < tol`.
    pub fn commutator_check(&self, n1: f64, n2: f64, tol: f64) -> bool {
        // [S^+, S^-] in linear HP = 2(S - n1)
        let commutator_val = 2.0 * (self.spin - n1);
        // 2*Sz at occupation n2
        let two_sz = 2.0 * self.transform_sz(n2);
        let deviation = (commutator_val - two_sz).abs() / (2.0 * self.spin);
        deviation < tol
    }

    /// Maximum meaningful boson occupation: `2S`.
    ///
    /// The HP boson space is restricted to occupation numbers `n = 0, 1, …, 2S`.
    /// States with `n > 2S` are unphysical — the HP Hilbert space has dimension
    /// `2S + 1`.
    ///
    /// # Returns
    ///
    /// `2S` \[dimensionless\].
    pub fn maximum_n_boson(&self) -> f64 {
        2.0 * self.spin
    }

    /// Occupation threshold above which the linear HP approximation degrades.
    ///
    /// The linear spin-wave expansion is reliable for `n/(2S) ≪ 1`.  A
    /// practical warning threshold is 10 % of S (i.e., 5 % of 2S):
    /// `n_warn = 0.1 * S`.
    ///
    /// # Returns
    ///
    /// `0.1 * S` \[dimensionless\].
    pub fn validity_warning_threshold(&self) -> f64 {
        0.1 * self.spin
    }

    /// Fractional occupation `n / (2S)`.
    ///
    /// Should satisfy `validity_fraction ≪ 1` for the linear approximation to
    /// be quantitatively valid.  A value above 0.1 signals a strongly quantum
    /// regime.
    ///
    /// # Arguments
    ///
    /// * `n_boson` — boson occupation n \[dimensionless\].
    ///
    /// # Returns
    ///
    /// `n / (2S)` \[dimensionless\].
    pub fn validity_fraction(&self, n_boson: f64) -> f64 {
        n_boson / (2.0 * self.spin)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── spin validation ──────────────────────────────────────────────────────

    #[test]
    fn test_negative_spin_rejected() {
        assert!(HolsteinPrimakoff::new(-0.5).is_err());
    }

    #[test]
    fn test_zero_spin_rejected() {
        assert!(HolsteinPrimakoff::new(0.0).is_err());
    }

    #[test]
    fn test_non_half_integer_rejected() {
        // 0.3 is not a multiple of 0.5
        assert!(HolsteinPrimakoff::new(0.3).is_err());
        // 1.7 is not a multiple of 0.5
        assert!(HolsteinPrimakoff::new(1.7).is_err());
    }

    #[test]
    fn test_valid_half_integer_spins() {
        assert!(HolsteinPrimakoff::new(0.5).is_ok());
        assert!(HolsteinPrimakoff::new(1.0).is_ok());
        assert!(HolsteinPrimakoff::new(1.5).is_ok());
        assert!(HolsteinPrimakoff::new(2.0).is_ok());
        assert!(HolsteinPrimakoff::new(5.0).is_ok());
    }

    // ── transform_sz ─────────────────────────────────────────────────────────

    #[test]
    fn test_transform_sz_at_n_zero_returns_spin() {
        let hp = HolsteinPrimakoff::new(1.5).expect("valid");
        let sz = hp.transform_sz(0.0);
        assert!((sz - 1.5).abs() < 1e-12, "Sz at n=0 should equal S={}", 1.5);
    }

    #[test]
    fn test_transform_sz_at_n_equals_spin_returns_zero() {
        // For spin S, at n = S we get Sz = S - S = 0
        let hp = HolsteinPrimakoff::new(2.0).expect("valid");
        let sz = hp.transform_sz(2.0);
        assert!(sz.abs() < 1e-12, "Sz at n=S should be 0, got {}", sz);
    }

    // ── splus / sminus coefficients ───────────────────────────────────────────

    #[test]
    fn test_splus_sminus_coefficients_equal_sqrt_2s() {
        for &spin in &[0.5_f64, 1.0, 1.5, 2.0, 5.0] {
            let hp = HolsteinPrimakoff::new(spin).expect("valid");
            let expected = (2.0 * spin).sqrt();
            let cp = hp.transform_splus_coeff();
            let cm = hp.transform_sminus_coeff();
            assert!(
                (cp - expected).abs() < 1e-12,
                "S^+ coeff mismatch for S={}: got {cp}, expected {expected}",
                spin
            );
            assert!(
                (cm - expected).abs() < 1e-12,
                "S^- coeff mismatch for S={}: got {cm}, expected {expected}",
                spin
            );
        }
    }

    // ── quadratic correction ──────────────────────────────────────────────────

    #[test]
    fn test_quadratic_correction_negative_for_n_greater_than_one() {
        let hp = HolsteinPrimakoff::new(2.0).expect("valid");
        // n=2: correction = -2*(2-1)/(4*2) = -0.25
        let corr = hp.quadratic_correction(2.0);
        assert!(
            corr < 0.0,
            "quadratic correction at n=2 should be negative, got {}",
            corr
        );
        assert!((corr - (-0.25)).abs() < 1e-12);
    }

    #[test]
    fn test_quadratic_correction_zero_at_n_zero() {
        let hp = HolsteinPrimakoff::new(1.0).expect("valid");
        assert!(hp.quadratic_correction(0.0).abs() < 1e-12);
    }

    #[test]
    fn test_quadratic_correction_zero_at_n_one() {
        let hp = HolsteinPrimakoff::new(1.0).expect("valid");
        // -1*(1-1)/(4S) = 0
        assert!(hp.quadratic_correction(1.0).abs() < 1e-12);
    }

    // ── commutator check ──────────────────────────────────────────────────────

    #[test]
    fn test_commutator_check_passes_at_low_density() {
        let hp = HolsteinPrimakoff::new(5.0).expect("valid");
        // At low occupation (n=0), commutator is exact
        assert!(hp.commutator_check(0.0, 0.0, 1e-10));
    }

    #[test]
    fn test_commutator_check_fails_at_high_density() {
        let hp = HolsteinPrimakoff::new(1.0).expect("valid");
        // At n1 >> n2 the linear approximation breaks down
        assert!(!hp.commutator_check(1.8, 0.0, 1e-3));
    }

    // ── preset constructors ───────────────────────────────────────────────────

    #[test]
    fn test_ferromagnet_default_spin_one() {
        let hp = HolsteinPrimakoff::ferromagnet_default();
        assert!((hp.spin - 1.0).abs() < 1e-12);
        assert_eq!(hp.expansion_order, HpOrder::Linear);
    }

    #[test]
    fn test_antiferromagnet_default_spin_half() {
        let hp = HolsteinPrimakoff::antiferromagnet_default();
        assert!((hp.spin - 0.5).abs() < 1e-12);
        assert_eq!(hp.expansion_order, HpOrder::Quadratic);
    }

    // ── maximum occupation and validity ──────────────────────────────────────

    #[test]
    fn test_maximum_n_boson_equals_2s() {
        let hp = HolsteinPrimakoff::new(3.0).expect("valid");
        assert!((hp.maximum_n_boson() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_validity_fraction() {
        let hp = HolsteinPrimakoff::new(2.0).expect("valid");
        // n=1 for S=2 → fraction = 1/4 = 0.25
        let frac = hp.validity_fraction(1.0);
        assert!(
            (frac - 0.25).abs() < 1e-12,
            "validity fraction mismatch: {frac}"
        );
    }
}
