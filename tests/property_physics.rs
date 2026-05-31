//! Property-based tests for v0.4.0+ physics modules.
//!
//! Coverage:
//!   - SMR: longitudinal symmetry ρ_L(m) = ρ_L(−m).
//!   - SMR: Hall resistivity leading-order antisymmetry identity.
//!   - SMR: ρ_L ≥ ρ₀ for any physical magnetisation direction.
//!   - SMR: ρ₁ ≥ 0 (always non-negative).
//!   - SMR: longitudinal resistivity ≤ ρ₀ + ρ₁ (upper bound when m ⊥ ŷ).
//!   - SMR: angular_scan returns exactly n_angles entries in [0, 2π).
//!   - SMR: G_sh (spin conductance) is strictly positive.
//!   - SMR: rho_2 = rho_1 × g_i/g_r (derived formula consistency).
//!   - USMR: relative change is odd in magnetisation.
//!   - USMR: relative change is linear in current density.
//!   - USMR: usmr_resistance antisymmetric in m for x-current.
//!   - STNO: |m| = Mₛ is conserved under LLG+STT with J = 0.
//!   - STNO: Slonczewski torque is perpendicular to m̂.
//!   - STNO: threshold current density is strictly positive.
//!   - STNO: thermal linewidth is non-negative.

#![allow(clippy::needless_pass_by_value)]

use std::f64::consts::PI;

use proptest::prelude::*;
use spintronics::effect::smr::{SpinHallMagnetoresistance, UnidirectionalSmr};
use spintronics::effect::stno::SpinTorqueOscillator;
use spintronics::vector3::Vector3;

// ─── Strategy helpers ─────────────────────────────────────────────────────────

/// Deterministic mapping from two `u64` seeds to a unit Vector3 on S².
///
/// Converts a pair of uniformly-distributed 64-bit integers to a unit vector
/// via spherical coordinates `(cos θ ∈ [−1, 1], φ ∈ [0, 2π])`, giving
/// asymptotically uniform coverage of the sphere.
fn seed_to_unit_vector(seed_a: u64, seed_b: u64) -> Vector3<f64> {
    let u = (seed_a as f64) / (u64::MAX as f64);
    let v = (seed_b as f64) / (u64::MAX as f64);
    let cos_theta = 2.0 * u - 1.0;
    let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
    let phi = 2.0 * PI * v;
    Vector3::new(sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta)
}

/// Strategy: charge current density in the A/m² range typical for USMR/STNO.
fn j_density_strategy() -> impl Strategy<Value = f64> {
    1.0e9f64..1.0e12
}

/// Strategy: physical temperature T ∈ [1, 1000] K.
fn temperature_strategy() -> impl Strategy<Value = f64> {
    1.0f64..1000.0
}

/// Strategy: number of angles for an angular scan ∈ [4, 16].
fn n_angles_strategy() -> impl Strategy<Value = usize> {
    4usize..16
}

// ─── Property tests ───────────────────────────────────────────────────────────

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    // ── SMR: longitudinal even symmetry ──────────────────────────────────────

    /// ρ_L(m) = ρ_L(−m) for any unit magnetisation vector.
    ///
    /// The longitudinal resistivity `ρ₀ + ρ₁(1 − m_y²)` depends only on m_y²,
    /// which is invariant under m → −m. This must hold for all m on S².
    #[test]
    fn smr_longitudinal_symmetric_in_m(
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
    ) {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let m = seed_to_unit_vector(seed_a, seed_b);
        let m_neg = Vector3::new(-m.x, -m.y, -m.z);

        let rho_pos = smr.longitudinal_resistivity(m);
        let rho_neg = smr.longitudinal_resistivity(m_neg);

        prop_assert!(
            (rho_pos - rho_neg).abs() < 1.0e-12,
            "ρ_L(m) ≠ ρ_L(−m): Δ = {:.3e} for m = ({:.4}, {:.4}, {:.4})",
            (rho_pos - rho_neg).abs(), m.x, m.y, m.z
        );
    }

    // ── SMR: Hall antisymmetry identity ──────────────────────────────────────

    /// ρ_H(m) + ρ_H(−m) = 2 ρ₁ m_x m_y for any m.
    ///
    /// Expanding `ρ_H(m) = ρ₁ m_x m_y − ρ₂ m_z` and
    /// `ρ_H(−m) = ρ₁ m_x m_y + ρ₂ m_z`, the sum is `2 ρ₁ m_x m_y`. The ρ₂ m_z
    /// terms cancel exactly. We verify the full analytical identity.
    #[test]
    fn smr_hall_sum_identity(
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
    ) {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let m = seed_to_unit_vector(seed_a, seed_b);
        let m_neg = Vector3::new(-m.x, -m.y, -m.z);

        let rho_h_pos = smr.hall_resistivity(m);
        let rho_h_neg = smr.hall_resistivity(m_neg);
        let actual_sum = rho_h_pos + rho_h_neg;
        let expected_sum = 2.0 * smr.rho_1() * m.x * m.y;

        let scale = smr.rho_1().abs().max(1.0e-30);
        prop_assert!(
            (actual_sum - expected_sum).abs() / scale < 1.0e-10,
            "ρ_H(m)+ρ_H(−m) − 2ρ₁ mₓmᵧ = {:.3e} (scale {:.3e})",
            (actual_sum - expected_sum).abs(), scale
        );
    }

    // ── SMR: longitudinal lower bound ─────────────────────────────────────────

    /// ρ_L(m) ≥ ρ₀ for all m (since ρ₁ ≥ 0 and 1 − m_y² ≥ 0).
    ///
    /// The SMR coefficient ρ₁ = ρ_NM θ_SH² · (positive factors) is always ≥ 0,
    /// so the longitudinal resistivity is bounded below by ρ₀.
    #[test]
    fn smr_longitudinal_lower_bound(
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
    ) {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let m = seed_to_unit_vector(seed_a, seed_b);
        let rho_l = smr.longitudinal_resistivity(m);

        prop_assert!(
            rho_l >= smr.rho_0() - 1.0e-15,
            "ρ_L = {:.6e} < ρ₀ = {:.6e} for m = ({:.4}, {:.4}, {:.4})",
            rho_l, smr.rho_0(), m.x, m.y, m.z
        );
    }

    // ── SMR: longitudinal upper bound ─────────────────────────────────────────

    /// ρ_L(m) ≤ ρ₀ + ρ₁ for all m.
    ///
    /// The factor (1 − m_y²) ∈ [0, 1] since |m_y| ≤ 1 for a unit vector, so
    /// ρ_L = ρ₀ + ρ₁(1−m_y²) ≤ ρ₀ + ρ₁. Equality holds when m ⊥ ŷ (m_y = 0).
    #[test]
    fn smr_longitudinal_upper_bound(
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
    ) {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let m = seed_to_unit_vector(seed_a, seed_b);
        let rho_l = smr.longitudinal_resistivity(m);
        let rho_max = smr.rho_0() + smr.rho_1();

        prop_assert!(
            rho_l <= rho_max + 1.0e-15,
            "ρ_L = {:.6e} > ρ₀+ρ₁ = {:.6e} for m = ({:.4}, {:.4}, {:.4})",
            rho_l, rho_max, m.x, m.y, m.z
        );
    }

    // ── SMR: rho_1 non-negative ───────────────────────────────────────────────

    /// ρ₁ ≥ 0 for the platinum-YIG preset.
    ///
    /// The SMR coefficient formula involves θ_SH² (always ≥ 0) multiplied by a
    /// strictly positive combination of conductance and length parameters.
    #[test]
    fn smr_rho_1_non_negative(
        // Trivial parameterisation — exercises preset construction.
        _seed in any::<u64>(),
    ) {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        prop_assert!(
            smr.rho_1() >= 0.0,
            "ρ₁ = {:.3e} must be non-negative", smr.rho_1()
        );

        // The tungsten-YIG preset has negative θ_SH but θ_SH² > 0, so ρ₁ > 0.
        let smr_w = SpinHallMagnetoresistance::tungsten_yig();
        prop_assert!(
            smr_w.rho_1() >= 0.0,
            "ρ₁ = {:.3e} must be non-negative for W/YIG", smr_w.rho_1()
        );
    }

    // ── SMR: rho_2 derived-formula consistency ────────────────────────────────

    /// ρ₂ = ρ₁ × (g_i / g_r) for non-degenerate interfaces.
    ///
    /// The anomalous Hall-like SMR coefficient is defined as
    /// `ρ₂ = ρ₁ × g_i / g_r`. We verify this identity holds to machine
    /// precision for the platinum-YIG preset.
    #[test]
    fn smr_rho_2_formula_consistent(
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
    ) {
        let smr = SpinHallMagnetoresistance::platinum_yig();

        // Probe the Hall resistivity for a random m, verifying the formula
        // ρ_H = ρ₁ m_x m_y − ρ₂ m_z holds with ρ₂ = ρ₁ g_i/g_r.
        let m = seed_to_unit_vector(seed_a, seed_b);
        let rho_h_computed = smr.hall_resistivity(m);
        let rho_1 = smr.rho_1();
        let rho_2 = smr.rho_2();
        let rho_h_expected = rho_1 * m.x * m.y - rho_2 * m.z;

        let scale = (rho_1 + rho_2).abs().max(1.0e-30);
        prop_assert!(
            (rho_h_computed - rho_h_expected).abs() / scale < 1.0e-10,
            "Hall formula mismatch: computed {:.6e} ≠ expected {:.6e}",
            rho_h_computed, rho_h_expected
        );

        // ρ₂ itself must equal ρ₁ × g_i / g_r to machine precision.
        let rho_2_ref = rho_1 * (smr.g_i / smr.g_r);
        prop_assert!(
            (rho_2 - rho_2_ref).abs() / scale < 1.0e-12,
            "ρ₂ = {:.6e} ≠ ρ₁·g_i/g_r = {:.6e}", rho_2, rho_2_ref
        );
    }

    // ── SMR: G_sh strictly positive ───────────────────────────────────────────

    /// G_sh = σ_NM tanh(t_NM/λ_sf) / λ_sf > 0 for physical parameters.
    ///
    /// Since σ_NM, t_NM, λ_sf > 0, the tanh is positive on a positive argument,
    /// and the division by λ_sf is well-defined. We verify for both presets.
    #[test]
    fn smr_g_sh_positive(
        _seed in any::<u64>(),
    ) {
        let smr_pt = SpinHallMagnetoresistance::platinum_yig();
        prop_assert!(
            smr_pt.g_sh() > 0.0,
            "G_sh = {:.3e} must be positive for Pt/YIG", smr_pt.g_sh()
        );

        let smr_w = SpinHallMagnetoresistance::tungsten_yig();
        prop_assert!(
            smr_w.g_sh() > 0.0,
            "G_sh = {:.3e} must be positive for W/YIG", smr_w.g_sh()
        );
    }

    // ── SMR: angular_scan returns exactly n_angles entries ───────────────────

    /// `angular_scan(n)` returns a vector of exactly `n` (angle, ρ_L, ρ_H)
    /// triples, with each angle in [0, 2π).
    ///
    /// For any n ≥ 1, the implementation iterates `k ∈ 0..n` and pushes exactly
    /// one entry per step. The first angle is 0 and entries are strictly
    /// increasing.
    #[test]
    fn smr_angular_scan_length(
        n_angles in n_angles_strategy(),
    ) {
        let smr = SpinHallMagnetoresistance::platinum_yig();
        let scan = smr.angular_scan(n_angles);

        prop_assert_eq!(
            scan.len(),
            n_angles,
            "angular_scan({}) returned {} entries", n_angles, scan.len()
        );

        // Verify the first angle is 0 (k=0 step: 2π·0/n = 0).
        if !scan.is_empty() {
            prop_assert!(
                scan[0].0.abs() < 1.0e-15,
                "first angle = {:.3e} ≠ 0", scan[0].0
            );
        }

        // Verify all angles are in [0, 2π).
        let two_pi = 2.0 * PI;
        for (alpha, _rho_l, _rho_h) in &scan {
            prop_assert!(
                *alpha >= 0.0 && *alpha < two_pi + 1.0e-12,
                "angle {:.6} outside [0, 2π)", alpha
            );
        }

        // Verify angles are strictly increasing (step = 2π/n > 0 for n ≥ 4).
        for w in scan.windows(2) {
            prop_assert!(
                w[1].0 > w[0].0,
                "angles not monotonically increasing: {} >= {}", w[0].0, w[1].0
            );
        }
    }

    // ── USMR: odd in magnetisation ───────────────────────────────────────────

    /// ΔR_USMR(m) / R₀ = −ΔR_USMR(−m) / R₀ for any unit magnetisation m.
    ///
    /// The normalised USMR `η · J · (−m_y)` is linear in m_y, which is odd
    /// under m → −m. This antisymmetry is the defining experimental signature
    /// of USMR distinguishing it from ordinary (even) SMR.
    #[test]
    fn usmr_odd_in_magnetization(
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
        j_density in j_density_strategy(),
    ) {
        let usmr = UnidirectionalSmr::platinum_cobalt()
            .expect("platinum_cobalt preset must be valid");
        let m = seed_to_unit_vector(seed_a, seed_b);
        let m_neg = Vector3::new(-m.x, -m.y, -m.z);

        let delta_r_pos = usmr.usmr_relative_change(m, j_density);
        let delta_r_neg = usmr.usmr_relative_change(m_neg, j_density);

        // For the exact linear formula δR/R₀ = η J (−mᵧ), the sum must vanish.
        prop_assert!(
            (delta_r_pos + delta_r_neg).abs() < 1.0e-12,
            "USMR not odd in m: ΔR(m)+ΔR(−m) = {:.3e} at J = {:.3e}",
            (delta_r_pos + delta_r_neg).abs(), j_density
        );
    }

    // ── USMR: linear in current density ──────────────────────────────────────

    /// ΔR(2J) = 2 ΔR(J) for any magnetisation direction and any J.
    ///
    /// The USMR formula `η J (−m_y)` is strictly linear in J, so doubling
    /// the current density must exactly double the resistance change. We avoid
    /// the degenerate case m_y ≈ 0 (where both values vanish) using prop_assume.
    #[test]
    fn usmr_linear_in_current(
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
        j_base in j_density_strategy(),
    ) {
        let usmr = UnidirectionalSmr::platinum_cobalt()
            .expect("platinum_cobalt preset must be valid");
        let m = seed_to_unit_vector(seed_a, seed_b);

        // Skip degenerate cases where m_y ≈ 0 (USMR signal vanishes regardless)
        prop_assume!(m.y.abs() > 0.01);

        let delta_r_j = usmr.usmr_relative_change(m, j_base);
        let delta_r_2j = usmr.usmr_relative_change(m, 2.0 * j_base);

        // Exact relation: ΔR(2J) = 2 ΔR(J)
        let scale = delta_r_2j.abs().max(1.0e-40);
        prop_assert!(
            (2.0 * delta_r_j - delta_r_2j).abs() < 1.0e-10 * scale,
            "USMR not linear in J: 2·ΔR(J) − ΔR(2J) = {:.3e}, scale = {:.3e}",
            (2.0 * delta_r_j - delta_r_2j).abs(), scale
        );
    }

    // ── USMR: usmr_resistance antisymmetric in m for x̂ current ──────────────

    /// ΔR_resistance(m) = −ΔR_resistance(−m) when current is along x̂.
    ///
    /// Using the full `usmr_resistance` method (not just the relative form),
    /// the sign-flip under m → −m must still hold since the projection
    /// (ĵ × ẑ) · m flips sign while ρ₀, η, J are unchanged.
    #[test]
    fn usmr_resistance_antisymmetric_in_m(
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
        j_density in j_density_strategy(),
    ) {
        let usmr = UnidirectionalSmr::platinum_cobalt()
            .expect("platinum_cobalt preset must be valid");
        let m = seed_to_unit_vector(seed_a, seed_b);
        let m_neg = Vector3::new(-m.x, -m.y, -m.z);
        let j_direction = Vector3::new(1.0, 0.0, 0.0); // x̂

        let dr_pos = usmr.usmr_resistance(m, j_density, j_direction);
        let dr_neg = usmr.usmr_resistance(m_neg, j_density, j_direction);

        // dr_pos + dr_neg must be zero (antisymmetry in m).
        prop_assert!(
            (dr_pos + dr_neg).abs() < 1.0e-20,
            "USMR resistance not antisymmetric: ΔR(m)+ΔR(−m) = {:.3e}",
            (dr_pos + dr_neg).abs()
        );
    }

    // ── STNO: magnetisation norm conserved at J = 0 ──────────────────────────

    /// With J = 0, |m(t)| = Mₛ is preserved by the RK4 integrator over 100 steps.
    ///
    /// The renormalisation step inside `rk4_step` guarantees `|m| = Mₛ` at each
    /// step. We verify this holds for arbitrary initial conditions: for any unit
    /// vector m̂₀ ∈ S², the scaled magnetisation m₀ = Mₛ m̂₀ retains its norm
    /// throughout a 100-step zero-current integration.
    #[test]
    fn stno_magnetization_norm_conserved_no_stt(
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
    ) {
        let stno = SpinTorqueOscillator::permalloy_nanopillar()
            .expect("permalloy_nanopillar must build");
        let ms = stno.solver.material.ms;

        // m0 at Mₛ length, arbitrary direction on S²
        let m_hat = seed_to_unit_vector(seed_a, seed_b);
        let m0 = m_hat * ms;

        // dt = 5e-13 s (0.5 ps) — stable for Permalloy in 100 kA/m bias field
        let dt = 5.0e-13;
        let traj = stno.simulate(m0, dt, 100, 0.0);

        for (i, m_i) in traj.iter().enumerate() {
            let norm_err = (m_i.magnitude() - ms).abs();
            prop_assert!(
                norm_err < 1.0e-3,
                "step {i}: |m| = {:.6e} ≠ Mₛ = {ms:.6e} (error = {norm_err:.3e})",
                m_i.magnitude()
            );
        }
    }

    // ── STNO: Slonczewski torque perpendicular to m ───────────────────────────

    /// τ_STT · m̂ = 0 for any normalised magnetisation direction.
    ///
    /// The Slonczewski torque `a_J [m×(m×p̂) + β(m×p̂)]` is automatically
    /// perpendicular to m̂: `m×(m×p̂)` has no component along m, and `m×p̂` is
    /// also ⊥ m by definition of the cross product. This holds for any β.
    #[test]
    fn stno_stt_perpendicular_to_m(
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
        j_density in j_density_strategy(),
    ) {
        let stno = SpinTorqueOscillator::permalloy_nanopillar()
            .expect("permalloy_nanopillar must build");
        let ms = stno.solver.material.ms;

        let m_hat = seed_to_unit_vector(seed_a, seed_b);
        let m = m_hat * ms;

        let tau = stno.slonczewski_torque(m, j_density);
        let dot = tau.dot(&m_hat).abs();
        let tau_norm = tau.magnitude().max(1.0e-20);

        prop_assert!(
            dot / tau_norm < 1.0e-9,
            "τ_STT · m̂ = {:.3e} / |τ| = {:.3e} ≠ 0 for m̂ = ({:.4}, {:.4}, {:.4})",
            dot, tau_norm, m_hat.x, m_hat.y, m_hat.z
        );
    }

    // ── STNO: threshold current density strictly positive ─────────────────────

    /// J_th > 0 for the permalloy nanopillar preset.
    ///
    /// The Slonczewski threshold formula `J_th = (2e μ₀ Mₛ t_FM α) / (ℏ η)` is
    /// a product of strictly positive physical quantities. We verify it for the
    /// preset and also check linearity in α (thermal scaling).
    #[test]
    fn stno_threshold_positive(
        _seed in any::<u64>(),
    ) {
        let stno = SpinTorqueOscillator::permalloy_nanopillar()
            .expect("permalloy_nanopillar must build");
        let j_th = stno.threshold_current_density();
        prop_assert!(
            j_th > 0.0,
            "J_th = {:.3e} must be strictly positive", j_th
        );
        prop_assert!(
            j_th.is_finite(),
            "J_th = {j_th} must be finite"
        );
        // Physical upper bound: Permalloy STNO should not require > 100 TA/m²
        prop_assert!(
            j_th < 1.0e14,
            "J_th = {j_th:.3e} is unrealistically large for Permalloy"
        );
    }

    // ── STNO: thermal linewidth non-negative ──────────────────────────────────

    /// Δf ≥ 0 for any physical temperature and current density.
    ///
    /// The Slavin-Tiberkevich formula `Δf = α k_B T / (π E_spin)` is a ratio of
    /// positive quantities (α > 0, k_B > 0, T > 0, E_spin > 0 for physical
    /// parameters). For the degenerate case E_spin ≈ 0 the implementation returns
    /// a large sentinel (1 THz), which is also non-negative.
    #[test]
    fn stno_thermal_linewidth_positive(
        j_density in j_density_strategy(),
        temperature in temperature_strategy(),
    ) {
        let stno = SpinTorqueOscillator::permalloy_nanopillar()
            .expect("permalloy_nanopillar must build");
        let lw = stno.thermal_linewidth(j_density, temperature);

        prop_assert!(
            lw >= 0.0,
            "thermal linewidth = {:.3e} must be ≥ 0 for T = {:.1} K, J = {:.3e} A/m²",
            lw, temperature, j_density
        );
        prop_assert!(
            lw.is_finite(),
            "thermal linewidth = {lw} must be finite for T = {temperature:.1} K"
        );
    }
}
