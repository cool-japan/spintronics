//! Property-based tests for the altermagnet electronic band model
//! (`AltermagnetBandModel`, `BlochHamiltonian`, `KuboBerry`).
//!
//! Coverage:
//!   - `H_sigma(k)` (via the `BlochHamiltonian` trait) is Hermitian for random
//!     k-points and random material parameters, any symmetry, with SOC on.
//!   - Berry curvature oddness under the combined spin + momentum + Neel-order
//!     inversion: `Omega_sigma(k; M) = -Omega_{-sigma}(-k; -M)`, for any symmetry.
//!   - Band-sum-to-zero: `Omega_Lower(k) + Omega_Upper(k) = 0` for a fixed spin,
//!     any symmetry.
//!   - Zero net spin polarization of the occupied Fermi sea for d-wave
//!     altermagnets without spin-orbit coupling, for any material parameters,
//!     Fermi energy, or Neel angle.
//!   - `KuboBerry` (finite-difference Kubo formula) agrees with the closed-form
//!     `AltermagnetBandModel::berry_curvature` at random k-points when SOC is on
//!     (nonzero curvature expected), and both vanish together when SOC is off.
//!   - Generalized spin-group relation `E_up(k) = E_down(R k)` (`R` = the
//!     `pi / harmonic_order` sublattice-swap rotation) for any of the d-wave,
//!     g-wave, or i-wave symmetries, without SOC.

#![allow(clippy::needless_pass_by_value)]
// proptest depends on rusty-fork -> wait-timeout, which has no wasm32 backend
// (no process-fork model on wasm32-unknown-unknown); this suite is native-only.
#![cfg(not(target_arch = "wasm32"))]

use std::f64::consts::TAU;

use proptest::prelude::*;
use spintronics::altermagnet::{
    AltermagnetBandModel, AltermagneticSymmetry, AltermagnetSpinHamiltonian, Band,
    BlochHamiltonian, KuboBerry, Spin,
};

// ---------------------------------------------------------------------------
// Strategy helpers
// ---------------------------------------------------------------------------

/// Strategy: one of the three altermagnetic symmetries (d/g/i-wave).
fn symmetry_strategy() -> impl Strategy<Value = AltermagneticSymmetry> {
    prop_oneof![
        Just(AltermagneticSymmetry::DWave),
        Just(AltermagneticSymmetry::GWave),
        Just(AltermagneticSymmetry::IWave),
    ]
}

/// Strategy: momentum component in a moderate range around the zone centre.
fn k_strategy() -> impl Strategy<Value = f64> {
    -2.5f64..2.5
}

/// Strategy: kinetic hopping scale `t_kin > 0`.
fn t_kin_strategy() -> impl Strategy<Value = f64> {
    0.3f64..2.5
}

/// Strategy: sublattice hybridization amplitude `delta_hyb >= 0`.
fn delta_hyb_strategy() -> impl Strategy<Value = f64> {
    0.0f64..1.0
}

/// Strategy: Neel exchange splitting `M >= 0`.
fn exchange_m_strategy() -> impl Strategy<Value = f64> {
    0.0f64..1.5
}

/// Strategy: altermagnetic form-factor amplitude `t_am >= 0`.
fn t_am_strategy() -> impl Strategy<Value = f64> {
    0.0f64..1.5
}

/// Strategy: nonzero spin-orbit-coupling amplitude (for SOC-on tests).
fn lambda_soc_nonzero_strategy() -> impl Strategy<Value = f64> {
    0.1f64..0.8
}

/// Strategy: Fermi energy spanning both bands' typical range.
fn fermi_energy_strategy() -> impl Strategy<Value = f64> {
    -2.0f64..3.0
}

/// Strategy: Neel-vector in-plane angle, full period.
fn neel_angle_strategy() -> impl Strategy<Value = f64> {
    0.0f64..TAU
}

/// Build a validated `AltermagnetBandModel` with unit lattice constant.
///
/// All parameters produced by the strategies above satisfy
/// `AltermagnetBandModel::new`'s validation (`t_kin > 0`, `delta_hyb >= 0`,
/// `exchange_m >= 0`, `t_am >= 0`, `a_lattice > 0`, finite `fermi_energy`), so
/// construction cannot fail.
#[allow(clippy::too_many_arguments)]
fn build_model(
    symmetry: AltermagneticSymmetry,
    t_kin: f64,
    delta_hyb: f64,
    exchange_m: f64,
    t_am: f64,
    lambda_soc: f64,
    fermi_energy: f64,
    neel_angle: f64,
) -> AltermagnetBandModel {
    AltermagnetBandModel::new(
        symmetry,
        t_kin,
        delta_hyb,
        exchange_m,
        t_am,
        lambda_soc,
        fermi_energy,
        neel_angle,
        1.0,
    )
    .expect("strategies only produce parameters within AltermagnetBandModel::new's valid range")
}

// ---------------------------------------------------------------------------
// Property tests: 32 cases per property, matching the crate's other
// tests/property_*.rs suites.
// ---------------------------------------------------------------------------

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    /// `H_sigma(k)` (accessed through the `BlochHamiltonian` trait) is
    /// Hermitian at any random momentum, for any symmetry, spin, and SOC
    /// strength.
    #[test]
    fn hamiltonian_is_hermitian(
        symmetry in symmetry_strategy(),
        t_kin in t_kin_strategy(),
        delta_hyb in delta_hyb_strategy(),
        exchange_m in exchange_m_strategy(),
        t_am in t_am_strategy(),
        lambda_soc in lambda_soc_nonzero_strategy(),
        fermi_energy in fermi_energy_strategy(),
        neel_angle in neel_angle_strategy(),
        kx in k_strategy(),
        ky in k_strategy(),
    ) {
        let model = build_model(
            symmetry, t_kin, delta_hyb, exchange_m, t_am, lambda_soc, fermi_energy, neel_angle,
        );
        for spin in [Spin::Up, Spin::Down] {
            let h = AltermagnetSpinHamiltonian::new(&model, spin);
            let matrix = h.hamiltonian_at(kx, ky).expect("2x2 construction cannot fail");
            for i in 0..matrix.n() {
                for j in 0..matrix.n() {
                    let hij = matrix.get(i, j);
                    let hji_conj = matrix.get(j, i).conj();
                    prop_assert!(
                        (hij.re - hji_conj.re).abs() < 1e-9 && (hij.im - hji_conj.im).abs() < 1e-9,
                        "H[{i}][{j}]={hij:?} != conj(H[{j}][{i}])={hji_conj:?} at k=({kx},{ky}), spin={spin:?}"
                    );
                }
            }
        }
    }

    /// Berry curvature oddness under the combined spin + momentum + Neel-order
    /// inversion: `Omega_sigma(k; M) = -Omega_{-sigma}(-k; -M)`.
    ///
    /// See `src/altermagnet/band_model.rs`'s
    /// `berry_curvature_odd_under_combined_spin_k_neel_flip` unit test for the
    /// derivation from the `Theta = i*tau_y*K` operator identity documented on
    /// `AltermagnetBandModel::spin_d_vector`. A literal same-spin
    /// `Omega_sigma(k) = -Omega_sigma(-k)` (holding `M` fixed) does *not* hold
    /// in general once `delta_hyb` and SOC are both nonzero -- this is the
    /// corrected, generally-true form.
    #[test]
    fn berry_curvature_odd_under_combined_flip(
        symmetry in symmetry_strategy(),
        t_kin in t_kin_strategy(),
        delta_hyb in delta_hyb_strategy(),
        exchange_m in exchange_m_strategy(),
        t_am in t_am_strategy(),
        lambda_soc in lambda_soc_nonzero_strategy(),
        neel_angle in neel_angle_strategy(),
        kx in k_strategy(),
        ky in k_strategy(),
    ) {
        let model = build_model(
            symmetry, t_kin, delta_hyb, exchange_m, t_am, lambda_soc, 0.5, neel_angle,
        );
        let mut reversed = model.clone();
        reversed.exchange_m = -model.exchange_m;

        for band in [Band::Lower, Band::Upper] {
            let omega_up = model.berry_curvature(kx, ky, Spin::Up, band);
            let omega_down_flipped = reversed.berry_curvature(-kx, -ky, Spin::Down, band);
            let scale = 1.0 + omega_up.abs().max(omega_down_flipped.abs());
            prop_assert!(
                (omega_up + omega_down_flipped).abs() < 1e-6 * scale,
                "Omega_Up(k;M) + Omega_Down(-k;-M) = {} (expected 0)",
                omega_up + omega_down_flipped
            );
        }
    }

    /// Band-sum-to-zero: for a fixed spin, `Omega_Lower(k) + Omega_Upper(k) = 0`.
    #[test]
    fn berry_curvature_band_sum_to_zero(
        symmetry in symmetry_strategy(),
        t_kin in t_kin_strategy(),
        delta_hyb in delta_hyb_strategy(),
        exchange_m in exchange_m_strategy(),
        t_am in t_am_strategy(),
        lambda_soc in lambda_soc_nonzero_strategy(),
        neel_angle in neel_angle_strategy(),
        kx in k_strategy(),
        ky in k_strategy(),
    ) {
        let model = build_model(
            symmetry, t_kin, delta_hyb, exchange_m, t_am, lambda_soc, 0.5, neel_angle,
        );
        for spin in [Spin::Up, Spin::Down] {
            let lower = model.berry_curvature(kx, ky, spin, Band::Lower);
            let upper = model.berry_curvature(kx, ky, spin, Band::Upper);
            let scale = 1.0 + lower.abs().max(upper.abs());
            prop_assert!(
                (lower + upper).abs() < 1e-9 * scale,
                "Omega_Lower + Omega_Upper = {} (expected 0)",
                lower + upper
            );
        }
    }

    /// Zero net spin polarization of the Fermi sea for d-wave altermagnets
    /// without SOC, for any material parameters, Fermi energy, or Neel angle.
    ///
    /// (The g-wave/i-wave presets are only *approximately* compensated at
    /// specific Fermi energies -- see
    /// `band_model::tests::net_spin_polarization_vanishes_for_presets` and
    /// `AltermagnetBandModel::sublattice_swap_rotation`'s docs for why only
    /// `n = 2` has an exact, parameter-independent grid symmetry.)
    #[test]
    fn net_spin_polarization_vanishes_dwave_no_soc(
        t_kin in t_kin_strategy(),
        delta_hyb in delta_hyb_strategy(),
        exchange_m in exchange_m_strategy(),
        t_am in t_am_strategy(),
        fermi_energy in fermi_energy_strategy(),
        neel_angle in neel_angle_strategy(),
    ) {
        let model = build_model(
            AltermagneticSymmetry::DWave, t_kin, delta_hyb, exchange_m, t_am, 0.0, fermi_energy, neel_angle,
        );
        let polarization = model.net_spin_polarization(31, 2.5).expect("valid grid parameters");
        prop_assert!(
            polarization.abs() < 1e-9,
            "net_spin_polarization = {polarization} (expected ~0 for a d-wave altermagnet without SOC)"
        );
    }

    /// `KuboBerry` agrees with the closed-form `berry_curvature` at random
    /// k-points when SOC is on (nonzero curvature expected).
    #[test]
    fn kubo_matches_closed_form_with_soc(
        lambda_soc in lambda_soc_nonzero_strategy(),
        kx in -1.2f64..1.2,
        ky in -1.2f64..1.2,
    ) {
        let model = build_model(
            AltermagneticSymmetry::DWave, 1.0, 0.3, 0.8, 0.6, lambda_soc, 0.5, 0.0,
        );
        for spin in [Spin::Up, Spin::Down] {
            let h = AltermagnetSpinHamiltonian::new(&model, spin);
            let kubo = KuboBerry::new(&h);
            for (band_idx, band) in [Band::Lower, Band::Upper].into_iter().enumerate() {
                let closed = model.berry_curvature(kx, ky, spin, band);
                let numeric = kubo.curvature_at(kx, ky, band_idx).expect("valid band index");
                let scale = closed.abs().max(numeric.abs()).max(1e-6);
                prop_assert!(
                    (closed - numeric).abs() / scale < 5e-3,
                    "closed={closed}, kubo={numeric} at k=({kx},{ky}), spin={spin:?}, band={band:?}"
                );
            }
        }
    }

    /// Both the closed form and `KuboBerry` vanish together when SOC is off.
    #[test]
    fn kubo_and_closed_form_vanish_together_without_soc(
        kx in -1.5f64..1.5,
        ky in -1.5f64..1.5,
        t_am in t_am_strategy(),
        exchange_m in exchange_m_strategy(),
    ) {
        let model = build_model(
            AltermagneticSymmetry::DWave, 1.0, 0.3, exchange_m, t_am, 0.0, 0.5, 0.0,
        );
        let h = AltermagnetSpinHamiltonian::new(&model, Spin::Up);
        let kubo = KuboBerry::new(&h);
        let closed = model.berry_curvature(kx, ky, Spin::Up, Band::Lower);
        let numeric = kubo.curvature_at(kx, ky, 0).expect("valid band index");
        prop_assert!(closed.abs() < 1e-9, "closed form should vanish without SOC: {closed}");
        prop_assert!(numeric.abs() < 1e-5, "Kubo cross-check should vanish without SOC: {numeric}");
    }

    /// Generalized spin-group relation `E_up(k) = E_down(R k)` for any harmonic
    /// order (d/g/i-wave), without SOC, where `R` is the `pi/n` sublattice-swap
    /// rotation (`n` = `AltermagneticSymmetry::harmonic_order`).
    #[test]
    fn spin_group_relation_generalized_rotation(
        symmetry in symmetry_strategy(),
        t_kin in t_kin_strategy(),
        delta_hyb in delta_hyb_strategy(),
        exchange_m in exchange_m_strategy(),
        t_am in t_am_strategy(),
        fermi_energy in fermi_energy_strategy(),
        neel_angle in neel_angle_strategy(),
        kx in k_strategy(),
        ky in k_strategy(),
    ) {
        let model = build_model(
            symmetry, t_kin, delta_hyb, exchange_m, t_am, 0.0, fermi_energy, neel_angle,
        );
        let (rkx, rky) = model.sublattice_swap_rotation(kx, ky);
        let (up, _) = model.spin_bands(kx, ky);
        let (_, down_rotated) = model.spin_bands(rkx, rky);
        prop_assert!(
            (up.lower - down_rotated.lower).abs() < 1e-8,
            "E_up,lower(k)={} != E_down,lower(Rk)={}", up.lower, down_rotated.lower
        );
        prop_assert!(
            (up.upper - down_rotated.upper).abs() < 1e-8,
            "E_up,upper(k)={} != E_down,upper(Rk)={}", up.upper, down_rotated.upper
        );
    }
}
