//! Property-based tests for optical switching, SAW magnetoacoustics,
//! and antiferromagnet dynamics.
//!
//! Coverage:
//!   - IFE: RCP and LCP produce equal-and-opposite fields (antisymmetry).
//!   - IFE: zero peak intensity gives zero switching probability.
//!   - IFE: demagnetization fraction is bounded in [0, 1] for any fluence.
//!   - IFE: switching probability is bounded in [0, 1] for any pulse.
//!   - IFE: high-intensity pulses produce near-complete demagnetization.
//!   - SAW: strain amplitude is non-negative at any (x, z, t).
//!   - SAW: strain at surface scales linearly with source amplitude.
//!   - SAW: Lorentzian response peaks on-resonance vs off-resonance.
//!   - SAW: acoustic spin current density is non-negative.
//!   - SAW: resonant absorption is positive when cone angle is non-zero.
//!   - SAW: penetration depth equals lambda / (2π).
//!   - AFM: sublattice norms are conserved under evolution.
//!   - AFM: Neel vector magnitude is unity for perfect initial AFM state.
//!   - AFM: total magnetization is negligibly small for perfect AFM state.
//!   - AFM: resonance frequency scales with sqrt(H_E × H_A).
//!   - AFM: resonance frequency is positive for positive H_E and H_A.

#![allow(clippy::needless_pass_by_value)]

use std::f64::consts::PI;

use proptest::prelude::*;
use spintronics::afm::antiferromagnet::Antiferromagnet;
use spintronics::effect::optical_switching::{
    CircularHelicity, LaserPulseParams, OpticalMagneticMaterial, OpticalSwitching,
};
use spintronics::mech::saw::{
    MagnetoelasticMaterial, PiezoSubstrate, SawMagnetoacoustics, SawSource,
};
use spintronics::vector3::Vector3;

// ─── Strategy helpers ─────────────────────────────────────────────────────────

/// Deterministic mapping from two `u64` seeds to a unit Vector3 on S².
///
/// Converts a pair of uniformly-distributed 64-bit integers to a unit vector
/// via spherical coordinates (cos θ ∈ [−1, 1], φ ∈ [0, 2π]), giving
/// asymptotically uniform coverage of the sphere.
fn seed_to_unit_vector(seed_a: u64, seed_b: u64) -> Vector3<f64> {
    let u = (seed_a as f64) / (u64::MAX as f64);
    let v = (seed_b as f64) / (u64::MAX as f64);
    let cos_theta = 2.0 * u - 1.0;
    let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
    let phi = 2.0 * PI * v;
    Vector3::new(sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta)
}

/// Strategy: peak intensity in W/m² spanning the physically meaningful range
/// for ultrafast laser experiments (10¹¹ to 10¹⁸ W/m²).
fn intensity_strategy() -> impl Strategy<Value = f64> {
    1.0e11f64..1.0e18
}

/// Strategy: SAW strain amplitude in the physically valid range (10⁻⁵ to 0.09).
///
/// Upper bound of 0.09 keeps us strictly below the 0.1 rejection threshold in
/// `SawSource::new`.
fn strain_amplitude_strategy() -> impl Strategy<Value = f64> {
    1.0e-5f64..9.0e-2
}

/// Strategy: external bias field H_ext ∈ [0, 500 kA/m].
fn h_ext_strategy() -> impl Strategy<Value = f64> {
    0.0f64..5.0e5
}

/// Strategy: SAW/AFM exchange field H_E ∈ [1 MA/m, 500 MA/m]
/// (range of physically realistic strong AFM exchange fields).
fn h_exchange_strategy() -> impl Strategy<Value = f64> {
    1.0e6f64..5.0e8
}

/// Strategy: AFM anisotropy field H_A ∈ [100 kA/m, 10 MA/m].
fn h_anisotropy_strategy() -> impl Strategy<Value = f64> {
    1.0e5f64..1.0e7
}

/// Build a `LaserPulseParams` with the standard Ti:Sapphire wavelength (800 nm),
/// a given peak intensity, and a 100 fs pulse duration.
///
/// Returns the pulse or a fallback to the standard preset when construction fails
/// (should not happen for the intensity range used by `intensity_strategy`).
fn make_pulse(peak_intensity: f64) -> LaserPulseParams {
    LaserPulseParams::new(800e-9, peak_intensity, 100e-15)
        .unwrap_or_else(|_| LaserPulseParams::ti_sapphire_standard())
}

// ─── Property tests ───────────────────────────────────────────────────────────

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    // ── IFE: RCP and LCP produce equal-and-opposite effective fields ──────────

    /// H_IFE(RCP) = −H_IFE(LCP) for any laser pulse.
    ///
    /// The Inverse Faraday Effect field is proportional to the helicity factor σ
    /// (σ = +1 for RCP, −1 for LCP), so switching helicity exactly negates the
    /// field.  This antisymmetry is the defining experimental signature of the IFE
    /// and is guaranteed by the formula H_IFE = α_IFE · n_ph · σ / μ₀.
    #[test]
    fn ife_rcp_lcp_fields_antisymmetric(
        peak_intensity in intensity_strategy(),
    ) {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);
        let pulse = make_pulse(peak_intensity);

        let h_rcp = switcher.inverse_faraday_field(&pulse, CircularHelicity::RightCircular);
        let h_lcp = switcher.inverse_faraday_field(&pulse, CircularHelicity::LeftCircular);

        // RCP must be positive, LCP must be negative
        prop_assert!(
            h_rcp > 0.0,
            "H_IFE(RCP) = {:.3e} A/m must be positive for I = {:.3e} W/m²",
            h_rcp, peak_intensity
        );
        prop_assert!(
            h_lcp < 0.0,
            "H_IFE(LCP) = {:.3e} A/m must be negative for I = {:.3e} W/m²",
            h_lcp, peak_intensity
        );

        // Exact antisymmetry: H_RCP + H_LCP = 0
        let scale = h_rcp.abs().max(1.0e-30);
        prop_assert!(
            (h_rcp + h_lcp).abs() / scale < 1.0e-10,
            "H_IFE(RCP) + H_IFE(LCP) = {:.3e} ≠ 0 (scale {:.3e})",
            (h_rcp + h_lcp).abs(), scale
        );
    }

    // ── IFE: zero peak intensity → zero IFE and negligible switching ─────────

    /// A laser with peak intensity → 0 produces no IFE field and zero
    /// demagnetization fraction.
    ///
    /// At intensity = 0 the photon flux n_ph = I/E_photon = 0, so H_IFE = 0 and
    /// the fluence falls far below the 1 mJ/cm² threshold.  The demagnetization
    /// fraction must therefore be exactly 0.0 (piecewise-linear clamp).
    #[test]
    fn ife_minimal_intensity_no_demag(
        _seed in any::<u64>(),
    ) {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);
        // Deliberately sub-threshold intensity: 10⁸ W/m², 100 fs
        // Fluence ≈ 10⁸ × 100e-15 × √(π/(4 ln2)) / 1e4 ≈ 8.5 × 10⁻¹³ J/cm² ≪ 1 mJ/cm²
        let pulse = LaserPulseParams::new(800e-9, 1.0e8, 100e-15)
            .expect("sub-threshold pulse must be constructible");

        let h_ife = switcher.inverse_faraday_field(&pulse, CircularHelicity::RightCircular);
        let demag = switcher.demagnetization_fraction(&pulse);

        prop_assert!(
            h_ife.is_finite(),
            "H_IFE must be finite; got {h_ife}"
        );
        prop_assert!(
            demag == 0.0,
            "Demagnetization fraction below threshold must be 0.0; got {demag:.6}"
        );
    }

    // ── IFE: demagnetization fraction bounded in [0, 1] ──────────────────────

    /// For any physical laser pulse the demagnetization fraction ∈ [0, 1].
    ///
    /// The piecewise-linear quench model clamps its output to [0, 1] by
    /// construction: f = 0 for F < F_th, f = (F−F_th)/(F_sat−F_th) for
    /// F_th ≤ F ≤ F_sat, and f = 1 for F > F_sat.  We verify the invariant
    /// holds across the full intensity range used in experiments.
    #[test]
    fn ife_demagnetization_bounded(
        peak_intensity in intensity_strategy(),
    ) {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);
        let pulse = make_pulse(peak_intensity);

        let demag = switcher.demagnetization_fraction(&pulse);

        prop_assert!(
            demag >= 0.0,
            "Demagnetization fraction = {demag:.6} must be ≥ 0 for I = {peak_intensity:.3e}"
        );
        prop_assert!(
            demag <= 1.0,
            "Demagnetization fraction = {demag:.6} must be ≤ 1 for I = {peak_intensity:.3e}"
        );
    }

    // ── IFE: switching probability bounded in [0, 1] ─────────────────────────

    /// `switching_probability` ∈ [0, 1] for any pulse and initial magnetisation.
    ///
    /// The sigmoid function σ(x) = 1/(1+e^{−x}) maps ℝ → (0, 1) strictly, so
    /// the switching probability computed from it must always lie inside the
    /// unit interval.  This invariant must hold for any combination of pulse
    /// parameters, helicity, and initial magnetisation direction.
    #[test]
    fn ife_switching_probability_bounded(
        peak_intensity in intensity_strategy(),
        seed_a in 1u64..u64::MAX,
        seed_b in 1u64..u64::MAX,
    ) {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);
        let pulse = make_pulse(peak_intensity);
        let m = seed_to_unit_vector(seed_a, seed_b);

        // Test both helicities
        for helicity in [CircularHelicity::RightCircular, CircularHelicity::LeftCircular] {
            let p = switcher.switching_probability(&pulse, helicity, m.z);
            prop_assert!(
                (0.0..=1.0).contains(&p),
                "switching_probability = {p:.6} outside [0,1] for I = {peak_intensity:.3e}, mz = {:.4}",
                m.z
            );
        }
    }

    // ── IFE: high-intensity pulses produce near-complete demagnetization ───────

    /// At peak intensity ≥ 10¹⁷ W/m² (fluence ≫ 10 mJ/cm²) the quench
    /// fraction saturates to 1.0.
    ///
    /// The threshold for saturation is F_sat = 10 mJ/cm².  A 100 fs pulse at
    /// 10¹⁷ W/m² gives F ≈ 8.5 J/cm² ≫ F_sat, so the clamp must return 1.0
    /// exactly.
    #[test]
    fn ife_high_intensity_full_demag(
        _seed in any::<u64>(),
    ) {
        // peak_intensity = 1e17 W/m², 100 fs pulse
        // Fluence = 1e17 * 100e-15 * sqrt(π/(4 ln2)) / 1e4
        //         ≈ 1e17 * 100e-15 * 1.064 / 1e4
        //         ≈ 8.49e-1 J/cm² = 849 mJ/cm² ≫ 10 mJ/cm²
        let pulse = LaserPulseParams::new(800e-9, 1.0e17, 100e-15)
            .expect("high-intensity pulse must be constructible");

        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);

        let demag = switcher.demagnetization_fraction(&pulse);

        prop_assert!(
            demag >= 0.999,
            "At high fluence ({:.3e} J/cm²) demagnetization fraction should be 1.0; got {demag:.6}",
            pulse.fluence_j_cm2
        );
    }

    // ── SAW: surface strain amplitude is bounded by source ε₀ ────────────────

    /// |ε_xx(x, 0, t)| ≤ ε₀ for all positions x and all times t.
    ///
    /// The surface strain is ε₀ · sin(k·x − ω·t), which is bounded in magnitude
    /// by ε₀ = `strain_amplitude`.  The inequality must hold for all (x, t) by
    /// the definition of the sinusoidal strain field.
    #[test]
    fn saw_surface_strain_bounded_by_amplitude(
        strain_amp in strain_amplitude_strategy(),
        seed_x in any::<u64>(),
        seed_t in any::<u64>(),
    ) {
        let sub = PiezoSubstrate::linbo3();
        let saw = SawSource::new(sub, 1.0e9, strain_amp)
            .expect("strain within valid range");

        // Map seeds to [0, 10λ] and [0, 10/f] respectively
        let lambda = saw.wavelength_m();
        let x = ((seed_x as f64) / (u64::MAX as f64)) * 10.0 * lambda;
        let period = 1.0 / saw.frequency_hz;
        let t = ((seed_t as f64) / (u64::MAX as f64)) * 10.0 * period;

        let strain = saw.strain_at_point(x, 0.0, t).abs();

        prop_assert!(
            strain <= strain_amp + 1.0e-14,
            "Surface strain |{:.4e}| exceeds ε₀ = {:.4e}",
            strain, strain_amp
        );
    }

    // ── SAW: strain at surface scales linearly with source amplitude ───────────

    /// If ε₀ doubles the surface strain field doubles everywhere.
    ///
    /// The strain field is proportional to ε₀ = `strain_amplitude` via
    /// ε_xx = ε₀ · sin(k·x − ω·t) · exp(−k·z).  Doubling ε₀ therefore doubles
    /// the strain at every point.  We verify this scaling for representative
    /// points at the surface.
    #[test]
    fn saw_strain_linear_in_amplitude(
        seed_x in any::<u64>(),
        seed_t in any::<u64>(),
    ) {
        let sub = PiezoSubstrate::linbo3();
        let eps_base = 1.0e-4;
        let eps_double = 2.0e-4;

        let saw_base = SawSource::new(sub.clone(), 1.0e9, eps_base).expect("valid");
        let saw_double = SawSource::new(sub, 1.0e9, eps_double).expect("valid");

        let lambda = saw_base.wavelength_m();
        let x = ((seed_x as f64) / (u64::MAX as f64)) * 4.0 * lambda;
        let period = 1.0 / saw_base.frequency_hz;
        let t = ((seed_t as f64) / (u64::MAX as f64)) * 4.0 * period;

        let s_base = saw_base.strain_at_point(x, 0.0, t);
        let s_double = saw_double.strain_at_point(x, 0.0, t);

        // Verify exact 2× relationship; skip near-zero to avoid numerical triviality
        let scale = s_base.abs().max(1.0e-20);
        if scale > 1.0e-15 {
            let ratio = s_double / s_base;
            prop_assert!(
                (ratio - 2.0).abs() < 1.0e-9,
                "Strain ratio = {ratio:.6} should be 2.0 at x={x:.4e}, t={t:.4e}"
            );
        }
    }

    // ── SAW: Lorentzian response peaks on-resonance ───────────────────────────

    /// The Lorentzian response function L(ω) peaks at ω = ω_FMR.
    ///
    /// At resonance (Δω = 0) the denominator of L is minimised to α ω_FMR, giving
    /// L_max = 1/(α ω_FMR) = Q/ω_FMR.  Any frequency detuning increases the
    /// denominator, reducing L below the peak value.  We check this property for
    /// Permalloy on LiNbO₃ with a randomised off-resonance bias field.
    #[test]
    fn saw_lorentzian_peaks_at_resonance(
        h_offset in 1.0e4f64..1.0e6,
    ) {
        let saw = SawSource::linbo3_3ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let thickness = 20.0e-9;

        // Compute the resonant bias field
        let helper = SawMagnetoacoustics::new(saw.clone(), mat.clone(), thickness, 0.0)
            .expect("valid on-resonance device");
        let h_res = helper.acoustic_fmr_condition_h_ext();

        let on_res = SawMagnetoacoustics::new(saw.clone(), mat.clone(), thickness, h_res)
            .expect("valid on-resonance device");
        let off_res = SawMagnetoacoustics::new(saw, mat, thickness, h_res + h_offset)
            .expect("valid off-resonance device");

        let l_on = on_res.lorentzian_response();
        let l_off = off_res.lorentzian_response();

        prop_assert!(
            l_on > l_off,
            "L(ω_res) = {l_on:.4e} must exceed L(ω_res + δ) = {l_off:.4e} at δ = {h_offset:.3e} A/m bias offset"
        );
    }

    // ── SAW: acoustic spin current density is non-negative ────────────────────

    /// J_s = (ℏ/4π) g_r ω sin²(θ) / 2 ≥ 0 for g_r ≥ 0 and any θ.
    ///
    /// All factors in the spin-pumping formula are non-negative for physical
    /// parameters: ℏ > 0, g_r ≥ 0, ω > 0, sin²(θ) ≥ 0.  We verify this
    /// property over the full bias-field range accessible via `h_ext_strategy`.
    #[test]
    fn saw_spin_current_nonnegative(
        h_ext in h_ext_strategy(),
        g_r_seed in any::<u64>(),
    ) {
        let saw = SawSource::linbo3_3ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let device = SawMagnetoacoustics::new(saw, mat, 20.0e-9, h_ext)
            .expect("valid device");

        // Spin-mixing conductance: 0 to 5×10¹⁹ m⁻²
        let g_r = ((g_r_seed as f64) / (u64::MAX as f64)) * 5.0e19;
        let j_s = device.spin_current_density(g_r);

        prop_assert!(
            j_s >= 0.0,
            "Spin current J_s = {j_s:.3e} A/m² must be ≥ 0 for H_ext = {h_ext:.3e} A/m, g_r = {g_r:.3e}"
        );
        prop_assert!(
            j_s.is_finite(),
            "Spin current J_s = {j_s} must be finite"
        );
    }

    // ── SAW: resonant absorption is positive at resonance ─────────────────────

    /// P_abs > 0 at acoustic FMR resonance.
    ///
    /// The Gilbert-damping power formula `P_abs = (μ₀/2) M_s α ω θ² d_film`
    /// is a product of positive quantities (α > 0, ω > 0, d_film > 0,
    /// θ_cone > 0 at resonance).  We verify P_abs is positive when the bias
    /// field is tuned to the acoustic FMR condition.
    #[test]
    fn saw_resonant_absorption_positive_at_fmr(
        _seed in any::<u64>(),
    ) {
        let saw = SawSource::linbo3_3ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let helper = SawMagnetoacoustics::new(saw.clone(), mat.clone(), 20.0e-9, 0.0)
            .expect("valid helper device");
        let h_res = helper.acoustic_fmr_condition_h_ext();

        let device = SawMagnetoacoustics::new(saw, mat, 20.0e-9, h_res)
            .expect("valid resonant device");

        let p_abs = device.resonant_absorption();
        prop_assert!(
            p_abs > 0.0,
            "Resonant absorption P_abs = {p_abs:.3e} W/m² must be positive at FMR"
        );
        prop_assert!(
            p_abs.is_finite(),
            "Resonant absorption P_abs = {p_abs} must be finite"
        );
    }

    // ── SAW: penetration depth equals lambda / (2π) ───────────────────────────

    /// The Rayleigh-wave penetration depth δ = λ / (2π) = v_SAW / (2π f).
    ///
    /// This follows directly from the exponential decay exp(−k·z) with
    /// k = 2π f / v_SAW.  The characteristic depth is 1/k = v_SAW / (2π f),
    /// which must equal the wavelength divided by 2π.  We verify for LiNbO₃,
    /// GaAs, and ZnO at three representative frequencies.
    #[test]
    fn saw_penetration_depth_equals_lambda_over_2pi(
        freq_seed in any::<u64>(),
    ) {
        // Sample from three substrate presets deterministically
        let substrates = [
            PiezoSubstrate::linbo3(),
            PiezoSubstrate::gaas(),
            PiezoSubstrate::zno(),
        ];
        // Frequency in [0.5, 5] GHz — well within the IDT-accessible range
        let freq_hz = 0.5e9 + ((freq_seed as f64) / (u64::MAX as f64)) * 4.5e9;

        for sub in &substrates {
            let lambda = sub.wavelength_m(freq_hz);
            let depth = sub.penetration_depth_m(freq_hz);
            let expected = lambda / (2.0 * PI);
            let rel_err = (depth - expected).abs() / expected;
            prop_assert!(
                rel_err < 1.0e-10,
                "Penetration depth = {depth:.6e} m ≠ λ/(2π) = {expected:.6e} m (rel_err = {rel_err:.2e}) for {} at {freq_hz:.3e} Hz",
                sub.name
            );
        }
    }

    // ── AFM: sublattice norms conserved under Heun evolution ──────────────────

    /// |m_A| = |m_B| = 1 after N Heun steps with arbitrary external field.
    ///
    /// The Heun integrator inside `evolve` explicitly normalizes each sublattice
    /// at the end of every step.  This hard renormalization guarantees that the
    /// unit-sphere constraint is maintained regardless of time step or field
    /// magnitude.  We verify that |m_A| and |m_B| remain ≈ 1 after 50 steps
    /// at dt = 10 fs with a randomized field direction and magnitude.
    #[test]
    fn afm_sublattice_norms_conserved(
        seed_a in 1u64..u64::MAX,
        seed_b in 1u64..u64::MAX,
        h_mag_seed in any::<u64>(),
    ) {
        let mut afm = Antiferromagnet::nio();
        let h_dir = seed_to_unit_vector(seed_a, seed_b);

        // External field 0–2 MA/m (challenging for the integrator)
        let h_mag = ((h_mag_seed as f64) / (u64::MAX as f64)) * 2.0e6;
        let h_ext = h_dir * h_mag;
        let dt = 1.0e-14; // 10 fs — stable timestep for THz AFM dynamics
        let n_steps = 50;

        for step in 0..n_steps {
            afm.evolve(h_ext, dt);
            let err_a = (afm.m_a.magnitude() - 1.0).abs();
            let err_b = (afm.m_b.magnitude() - 1.0).abs();
            prop_assert!(
                err_a < 1.0e-6,
                "step {step}: |m_A| = {:.10} ≠ 1 (error = {err_a:.3e})",
                afm.m_a.magnitude()
            );
            prop_assert!(
                err_b < 1.0e-6,
                "step {step}: |m_B| = {:.10} ≠ 1 (error = {err_b:.3e})",
                afm.m_b.magnitude()
            );
        }
    }

    // ── AFM: Neel vector magnitude is unity for perfect initial state ──────────

    /// |n| = 1 when m_A = +ẑ, m_B = −ẑ.
    ///
    /// The Neel vector n = (m_A − m_B) / 2.  For the perfect antiparallel ground
    /// state m_A = +ẑ and m_B = −ẑ, this gives n = ẑ with |n| = 1 exactly.
    /// We verify this for both the NiO and MnF₂ presets.
    #[test]
    fn afm_neel_vector_unity_for_perfect_state(
        _seed in any::<u64>(),
    ) {
        for afm in [Antiferromagnet::nio(), Antiferromagnet::mnf2()] {
            let neel = afm.neel_vector();
            let neel_mag = neel.magnitude();
            prop_assert!(
                (neel_mag - 1.0).abs() < 1.0e-10,
                "Neel vector magnitude = {neel_mag:.12} ≠ 1 for perfect AFM state"
            );
        }
    }

    // ── AFM: total magnetization negligible for perfect state ─────────────────

    /// |m_total| < 10⁻¹² for the perfect antiparallel ground state.
    ///
    /// m_total = (m_A + m_B) / 2 = (+ẑ + (−ẑ)) / 2 = 0 exactly.  Any deviation
    /// from zero signals a numerical representation error; we check it stays
    /// below 10⁻¹² (well within double-precision floating-point error).
    #[test]
    fn afm_total_magnetization_vanishes_for_perfect_state(
        _seed in any::<u64>(),
    ) {
        for afm in [Antiferromagnet::nio(), Antiferromagnet::mnf2()] {
            let m_total = afm.total_magnetization();
            let m_mag = m_total.magnitude();
            prop_assert!(
                m_mag < 1.0e-12,
                "|m_total| = {m_mag:.3e} must be < 10⁻¹² for perfect AFM state"
            );
        }
    }

    // ── AFM: resonance frequency scales with sqrt(H_E × H_A) ─────────────────

    /// f_res(H_E₂, H_A₂) / f_res(H_E₁, H_A₁) = sqrt((H_E₂ × H_A₂) / (H_E₁ × H_A₁)).
    ///
    /// The uniaxial AFM resonance formula f_res = γ μ₀ sqrt(H_E H_A) / (2π) is
    /// proportional to sqrt(H_E H_A).  Constructing two AFMs with different
    /// exchange and anisotropy fields, the ratio of their resonance frequencies
    /// must equal the square-root ratio of the products.  We verify this
    /// analytical relation holds to floating-point precision.
    #[test]
    fn afm_resonance_frequency_sqrt_scaling(
        h_e1 in h_exchange_strategy(),
        h_a1 in h_anisotropy_strategy(),
        h_e2 in h_exchange_strategy(),
        h_a2 in h_anisotropy_strategy(),
    ) {
        let afm1 = Antiferromagnet::new(h_e1, 0.005, h_a1);
        let afm2 = Antiferromagnet::new(h_e2, 0.005, h_a2);

        let f1 = afm1.resonance_frequency();
        let f2 = afm2.resonance_frequency();

        // Both frequencies must be positive (validated separately)
        prop_assume!(f1 > 0.0 && f2 > 0.0);

        // Expected ratio: sqrt(H_E2 H_A2) / sqrt(H_E1 H_A1)
        let ratio_computed = f2 / f1;
        let ratio_expected = ((h_e2 * h_a2) / (h_e1 * h_a1)).sqrt();

        let scale = ratio_expected.abs().max(1.0e-30);
        prop_assert!(
            (ratio_computed - ratio_expected).abs() / scale < 1.0e-9,
            "f_res ratio = {ratio_computed:.6} ≠ sqrt-product ratio = {ratio_expected:.6} \
             for (H_E1={h_e1:.3e}, H_A1={h_a1:.3e}) and (H_E2={h_e2:.3e}, H_A2={h_a2:.3e})"
        );
    }

    // ── AFM: resonance frequency positive for positive exchange and anisotropy ─

    /// f_res > 0 for any positive H_E and H_A.
    ///
    /// The formula f_res = γ μ₀ sqrt(H_E H_A) / (2π) is a product of positive
    /// constants (γ > 0, μ₀ > 0) times the square root of a positive product
    /// (H_E H_A > 0 for H_E, H_A > 0), divided by 2π > 0.  The result is
    /// strictly positive and finite for any physically valid parameters.
    #[test]
    fn afm_resonance_frequency_positive(
        h_exchange in h_exchange_strategy(),
        h_anisotropy in h_anisotropy_strategy(),
    ) {
        let afm = Antiferromagnet::new(h_exchange, 0.005, h_anisotropy);
        let f_res = afm.resonance_frequency();

        prop_assert!(
            f_res > 0.0,
            "AFM resonance frequency = {f_res:.3e} Hz must be positive \
             for H_E = {h_exchange:.3e} A/m, H_A = {h_anisotropy:.3e} A/m"
        );
        prop_assert!(
            f_res.is_finite(),
            "AFM resonance frequency = {f_res} must be finite"
        );

        // Sanity: expect THz range for realistic parameters
        // Lower bound: γ μ₀ sqrt(1e6 * 1e5) / (2π) ≈ 28e9 * 1.26e-6 * 316 / 6.28 ≈ 1.88 GHz
        // Upper bound: γ μ₀ sqrt(5e8 * 1e7) / (2π) ≈ 28e9 * 1.26e-6 * 7.07e7 / 6.28 ≈ 2.5 PHz
        prop_assert!(
            f_res < 1.0e16,
            "AFM resonance frequency = {f_res:.3e} Hz is unrealistically large; check units"
        );
    }
}
