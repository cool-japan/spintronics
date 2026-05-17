//! Differentiable physics functions for spintronic model fitting.
//!
//! Every function here is built entirely from [`Var`] arithmetic, so gradients
//! w.r.t. any free parameter are obtained automatically by calling
//! [`Tape::backward`] after the forward pass.
//!
//! # Included functions
//!
//! | Function | Physics |
//! |----------|---------|
//! | [`kittel_frequency_diff`] | FMR resonance frequency (Kittel, 1948) |
//! | [`zeeman_energy_diff`]    | Zeeman coupling energy |
//! | [`exchange_energy_diff`]  | Heisenberg nearest-neighbour exchange |
//! | [`dmi_energy_diff`]       | Interfacial Dzyaloshinskii–Moriya interaction |
//! | [`anisotropy_energy_diff`]| Uniaxial magnetic anisotropy |
//! | [`llg_torque_norm_diff`]  | LLG precessional torque magnitude |
//!
//! # References
//!
//! - C. Kittel, "On the Theory of Ferromagnetic Resonance Absorption",
//!   *Phys. Rev.* **73**, 155 (1948)
//! - I. E. Dzyaloshinsky, "A Thermodynamic Theory of 'Weak' Ferromagnetism …",
//!   *J. Phys. Chem. Solids* **4**, 241 (1958)
//! - T. Moriya, "Anisotropic Superexchange Interaction and Weak Ferromagnetism",
//!   *Phys. Rev.* **120**, 91 (1960)
//! - W. F. Brown, *Micromagnetics*, Wiley 1963

use crate::autodiff::tape::{Tape, Var};
use crate::constants::{GAMMA, MU_0};

// ─── Kittel FMR frequency ─────────────────────────────────────────────────────

/// Differentiable Kittel ferromagnetic resonance frequency.
///
/// The Kittel formula for an in-plane magnetized thin film:
///
/// ω = |γ| · μ₀ · √[ H_ext · (H_ext + M_s) ]
///
/// where
/// - `h_ext`: external field amplitude [A/m]
/// - `ms`: saturation magnetisation [A/m]
/// - `GAMMA`: gyromagnetic ratio |γ| [rad/(s·T)]
/// - `MU_0`: vacuum permeability [H/m]
///
/// Returns angular frequency ω [rad/s].
///
/// # Differentiability
/// Gradients are available w.r.t. both `ms` and `h_ext` after `tape.backward`.
pub fn kittel_frequency_diff<'t>(_tape: &'t Tape, ms: Var<'t>, h_ext: Var<'t>) -> Var<'t> {
    let gamma = GAMMA.abs();
    let mu0 = MU_0;
    // ω = γ·μ₀·√(h·(h + ms))
    let h_plus_ms = h_ext + ms;
    let product = h_ext * h_plus_ms;
    let sq = product.sqrt();
    sq * (gamma * mu0)
}

// ─── Zeeman energy ────────────────────────────────────────────────────────────

/// Differentiable Zeeman energy density.
///
/// E_Zeeman = −μ₀ · ms · (m · H)
///           = −μ₀ · ms · (mx·Hx + my·Hy + mz·Hz)
///
/// The field components `(hx, hy, hz)` are treated as constants (f64).
/// The magnetisation direction `(mx, my, mz)` and the saturation magnetisation
/// `ms` are differentiable [`Var`] inputs.
///
/// Returns energy density [J/m³].
pub fn zeeman_energy_diff<'t>(
    _tape: &'t Tape,
    mx: Var<'t>,
    my: Var<'t>,
    mz: Var<'t>,
    hx: f64,
    hy: f64,
    hz: f64,
    ms: Var<'t>,
) -> Var<'t> {
    let mu0 = MU_0;
    // m·H = mx*hx + my*hy + mz*hz
    let mdoth = mx * hx + my * hy + mz * hz;
    // E = -μ₀ · ms · (m·H)
    let neg_mu0 = -mu0;
    ms * mdoth * neg_mu0
}

// ─── Exchange energy ──────────────────────────────────────────────────────────

/// Differentiable Heisenberg exchange energy for a nearest-neighbour pair.
///
/// The bilinear Heisenberg Hamiltonian for a single bond:
///
/// E_ex = −(2 A_ex / a²) · (m₁ · m₂)
///
/// where
/// - `a_ex`: exchange stiffness constant [J/m], differentiable
/// - `lattice_const`: lattice parameter `a` [m], constant
/// - `m₁ = (m1x, m1y, m1z)`, `m₂ = (m2x, m2y, m2z)`: unit-vector spin
///   directions, differentiable
///
/// Returns energy [J] (assuming unit volume / bond).
pub fn exchange_energy_diff<'t>(
    _tape: &'t Tape,
    m1x: Var<'t>,
    m1y: Var<'t>,
    m1z: Var<'t>,
    m2x: Var<'t>,
    m2y: Var<'t>,
    m2z: Var<'t>,
    a_ex: Var<'t>,
    lattice_const: f64,
) -> Var<'t> {
    // m₁ · m₂
    let dot = m1x * m2x + m1y * m2y + m1z * m2z;
    // −2 A_ex / a²  (scalar part)
    let coeff = -2.0 / (lattice_const * lattice_const);
    a_ex * dot * coeff
}

// ─── DMI energy ───────────────────────────────────────────────────────────────

/// Differentiable interfacial Dzyaloshinskii–Moriya interaction energy.
///
/// Interfacial DMI with a DM vector along ẑ:
///
/// E_DMI = D · (m₁ × m₂)_z = D · (m1x·m2y − m1y·m2x)
///
/// The DMI constant `d` [J/m²] is a differentiable variable.
/// The spin directions `m₁`, `m₂` are also differentiable.
///
/// Returns energy density [J/m²] (per interface area).
pub fn dmi_energy_diff<'t>(
    _tape: &'t Tape,
    m1x: Var<'t>,
    m1y: Var<'t>,
    _m1z: Var<'t>,
    m2x: Var<'t>,
    m2y: Var<'t>,
    _m2z: Var<'t>,
    d: Var<'t>,
) -> Var<'t> {
    // (m₁ × m₂)_z = m1x·m2y − m1y·m2x
    let cross_z = m1x * m2y - m1y * m2x;
    d * cross_z
}

// ─── Uniaxial anisotropy energy ───────────────────────────────────────────────

/// Differentiable uniaxial magnetocrystalline anisotropy energy.
///
/// E_aniso = −Ku · mz²
///
/// where `ku` is the uniaxial anisotropy constant [J/m³] (differentiable)
/// and `mz` is the z-component of the unit magnetisation vector (differentiable).
///
/// Positive `ku` favours magnetisation along ẑ.
///
/// Returns energy density [J/m³].
pub fn anisotropy_energy_diff<'t>(_tape: &'t Tape, mz: Var<'t>, ku: Var<'t>) -> Var<'t> {
    // E = -Ku * mz²
    let mz_sq = mz * mz;
    ku * mz_sq * (-1.0)
}

// ─── LLG torque magnitude ─────────────────────────────────────────────────────

/// Differentiable LLG precessional torque magnitude.
///
/// Returns |dm/dt|_prec ≈ |γ| · |m × H_eff|, where
///
/// |m × H_eff|² = (my·Hz − mz·Hy)² + (mz·Hx − mx·Hz)² + (mx·Hy − my·Hx)²
///
/// The effective field components `(heff_x, heff_y, heff_z)` are constants;
/// the magnetisation components `(mx, my, mz)` and damping `alpha` are
/// differentiable (though damping only enters prefactors in the full LLG —
/// here it is used to return the full |dm/dt| including the damping torque
/// term).  The `gamma` parameter is the gyromagnetic ratio magnitude.
///
/// Returns |dm/dt| in [rad/s] (precessional part only).
pub fn llg_torque_norm_diff<'t>(
    _tape: &'t Tape,
    mx: Var<'t>,
    my: Var<'t>,
    mz: Var<'t>,
    heff_x: f64,
    heff_y: f64,
    heff_z: f64,
    _alpha: Var<'t>,
    gamma: f64,
) -> Var<'t> {
    // (m × H_eff) components
    let cross_x = my * heff_z - mz * heff_y;
    let cross_y = mz * heff_x - mx * heff_z;
    let cross_z = mx * heff_y - my * heff_x;
    // |m × H_eff|²
    let norm_sq = cross_x * cross_x + cross_y * cross_y + cross_z * cross_z;
    // |γ| · |m × H_eff|
    norm_sq.sqrt() * gamma
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::autodiff::tape::finite_diff_grad;
    use crate::constants::{GAMMA, MU_0};

    // Tolerance for physics comparisons (double precision)
    const TOL: f64 = 1e-9;

    // 1. kittel_frequency_diff: value matches non-AD Kittel formula
    #[test]
    fn test_kittel_value_matches_formula() {
        let ms_val = 1.4e5_f64; // A/m (e.g. YIG-like)
        let h_val = 1.592e5_f64; // A/m (~200 mT)

        let tape = Tape::new();
        let ms = Var::leaf(&tape, ms_val);
        let h = Var::leaf(&tape, h_val);
        let omega = kittel_frequency_diff(&tape, ms, h);

        // Reference: ω = |γ|·μ₀·√(H·(H+Ms))
        let expected = GAMMA.abs() * MU_0 * (h_val * (h_val + ms_val)).sqrt();
        assert!(
            (omega.value() - expected).abs() < TOL * expected.abs().max(1.0),
            "Kittel value mismatch: got {}, expected {}",
            omega.value(),
            expected
        );
    }

    // 2. dω/dms at known point (finite-difference check)
    #[test]
    fn test_kittel_grad_wrt_ms() {
        let ms_val = 1.4e5_f64;
        let h_val = 1.592e5_f64;

        // AD gradient
        let tape = Tape::new();
        let ms = Var::leaf(&tape, ms_val);
        let h = Var::leaf(&tape, h_val);
        let omega = kittel_frequency_diff(&tape, ms, h);
        tape.backward(omega);
        let ad_grad = ms.grad();

        // Finite difference w.r.t. ms
        let fd_grad = finite_diff_grad(
            |ms_i| {
                let t2 = Tape::new();
                let ms2 = Var::leaf(&t2, ms_i);
                let h2 = Var::leaf(&t2, h_val);
                kittel_frequency_diff(&t2, ms2, h2).value()
            },
            ms_val,
            ms_val * 1e-5,
        );

        let rel_err = (ad_grad - fd_grad).abs() / fd_grad.abs().max(1.0);
        assert!(
            rel_err < 1e-5,
            "dω/dms AD ({}) vs FD ({}), rel_err={}",
            ad_grad,
            fd_grad,
            rel_err
        );
    }

    // 3. dω/dh_ext at known point (finite-difference check)
    #[test]
    fn test_kittel_grad_wrt_hext() {
        let ms_val = 1.4e5_f64;
        let h_val = 1.592e5_f64;

        let tape = Tape::new();
        let ms = Var::leaf(&tape, ms_val);
        let h = Var::leaf(&tape, h_val);
        let omega = kittel_frequency_diff(&tape, ms, h);
        tape.backward(omega);
        let ad_grad = h.grad();

        let fd_grad = finite_diff_grad(
            |hi| {
                let t2 = Tape::new();
                let ms2 = Var::leaf(&t2, ms_val);
                let h2 = Var::leaf(&t2, hi);
                kittel_frequency_diff(&t2, ms2, h2).value()
            },
            h_val,
            h_val * 1e-5,
        );

        let rel_err = (ad_grad - fd_grad).abs() / fd_grad.abs().max(1.0);
        assert!(
            rel_err < 1e-5,
            "dω/dh AD ({}) vs FD ({}), rel_err={}",
            ad_grad,
            fd_grad,
            rel_err
        );
    }

    // 4. zeeman_energy_diff: value matches analytical
    #[test]
    fn test_zeeman_value() {
        let ms_val = 1.4e5_f64;
        let (mx_val, my_val, mz_val) = (0.6, 0.8, 0.0);
        let (hx, hy, hz) = (1e5_f64, 0.0, 0.0);

        let tape = Tape::new();
        let mx = Var::leaf(&tape, mx_val);
        let my = Var::leaf(&tape, my_val);
        let mz = Var::leaf(&tape, mz_val);
        let ms = Var::leaf(&tape, ms_val);
        let e = zeeman_energy_diff(&tape, mx, my, mz, hx, hy, hz, ms);

        // E = -μ₀·ms·(m·H)
        let mdoth = mx_val * hx + my_val * hy + mz_val * hz;
        let expected = -MU_0 * ms_val * mdoth;
        assert!(
            (e.value() - expected).abs() < TOL * expected.abs().max(1.0),
            "Zeeman energy: got {}, expected {}",
            e.value(),
            expected
        );
    }

    // 5. exchange_energy_diff: orthogonal spins → dot product = 0 → E = 0
    #[test]
    fn test_exchange_orthogonal_spins() {
        let tape = Tape::new();
        let m1x = Var::leaf(&tape, 1.0);
        let m1y = Var::leaf(&tape, 0.0);
        let m1z = Var::leaf(&tape, 0.0);
        let m2x = Var::leaf(&tape, 0.0);
        let m2y = Var::leaf(&tape, 1.0);
        let m2z = Var::leaf(&tape, 0.0);
        let a_ex = Var::leaf(&tape, 3.5e-12_f64); // typical exchange stiffness
        let lattice_const = 3.5e-10_f64; // 3.5 Å

        let e = exchange_energy_diff(&tape, m1x, m1y, m1z, m2x, m2y, m2z, a_ex, lattice_const);
        // m₁ · m₂ = 0 → E = 0
        assert!(
            e.value().abs() < 1e-30,
            "exchange energy of orthogonal spins should be 0"
        );
    }

    // 6. exchange_energy_diff: parallel spins → E = -2·A/a²
    #[test]
    fn test_exchange_parallel_spins() {
        let a_ex_val = 3.5e-12_f64;
        let a = 3.5e-10_f64;

        let tape = Tape::new();
        let m1x = Var::leaf(&tape, 1.0);
        let m1y = Var::leaf(&tape, 0.0);
        let m1z = Var::leaf(&tape, 0.0);
        let m2x = Var::leaf(&tape, 1.0);
        let m2y = Var::leaf(&tape, 0.0);
        let m2z = Var::leaf(&tape, 0.0);
        let a_ex = Var::leaf(&tape, a_ex_val);

        let e = exchange_energy_diff(&tape, m1x, m1y, m1z, m2x, m2y, m2z, a_ex, a);
        let expected = -2.0 * a_ex_val / (a * a); // m₁·m₂ = 1
        assert!(
            (e.value() - expected).abs() < 1e-3 * expected.abs().max(1.0),
            "parallel exchange: got {}, expected {}",
            e.value(),
            expected
        );
    }

    // 7. dmi_energy_diff: value and gradient
    #[test]
    fn test_dmi_value_and_grad() {
        // m₁ = (1, 0, 0), m₂ = (0, 1, 0) → (m₁ × m₂)_z = 1
        let d_val = 1.5e-3_f64;

        let tape = Tape::new();
        let m1x = Var::leaf(&tape, 1.0);
        let m1y = Var::leaf(&tape, 0.0);
        let m1z = Var::leaf(&tape, 0.0);
        let m2x = Var::leaf(&tape, 0.0);
        let m2y = Var::leaf(&tape, 1.0);
        let m2z = Var::leaf(&tape, 0.0);
        let d = Var::leaf(&tape, d_val);

        let e = dmi_energy_diff(&tape, m1x, m1y, m1z, m2x, m2y, m2z, d);
        tape.backward(e);

        let expected_val = d_val * 1.0; // D · 1
        assert!(
            (e.value() - expected_val).abs() < TOL * expected_val.abs().max(1.0),
            "DMI value: got {}, expected {}",
            e.value(),
            expected_val
        );
        // ∂E/∂D = (m₁ × m₂)_z = 1.0
        assert!(
            (d.grad() - 1.0).abs() < 1e-12,
            "DMI ∂E/∂D should be 1.0, got {}",
            d.grad()
        );
    }

    // 8. anisotropy_energy_diff: value for mz=1 should give -Ku
    #[test]
    fn test_anisotropy_value() {
        let ku_val = 5e4_f64; // J/m³ (typical PMA value)

        let tape = Tape::new();
        let mz = Var::leaf(&tape, 1.0);
        let ku = Var::leaf(&tape, ku_val);
        let e = anisotropy_energy_diff(&tape, mz, ku);

        let expected = -ku_val * 1.0; // -Ku·mz²
        assert!(
            (e.value() - expected).abs() < TOL * expected.abs().max(1.0),
            "anisotropy: got {}, expected {}",
            e.value(),
            expected
        );
    }

    // 9. anisotropy gradient check: ∂E/∂mz = -2·Ku·mz
    #[test]
    fn test_anisotropy_grad_mz() {
        let ku_val = 5e4_f64;
        let mz_val = 0.7_f64;

        let tape = Tape::new();
        let mz = Var::leaf(&tape, mz_val);
        let ku = Var::leaf(&tape, ku_val);
        let e = anisotropy_energy_diff(&tape, mz, ku);
        tape.backward(e);

        let expected_grad = -2.0 * ku_val * mz_val; // ∂(-Ku·mz²)/∂mz
        assert!(
            (mz.grad() - expected_grad).abs() < 1.0,
            "anisotropy ∂E/∂mz: got {}, expected {}",
            mz.grad(),
            expected_grad
        );
    }

    // 10. llg_torque_norm_diff: m along x, H along y → |m×H| = |H|
    #[test]
    fn test_llg_torque_norm_value() {
        let h_mag = 1e5_f64; // A/m

        let tape = Tape::new();
        let mx = Var::leaf(&tape, 1.0);
        let my = Var::leaf(&tape, 0.0);
        let mz = Var::leaf(&tape, 0.0);
        let alpha = Var::leaf(&tape, 0.01);

        let torque = llg_torque_norm_diff(&tape, mx, my, mz, 0.0, h_mag, 0.0, alpha, GAMMA.abs());
        tape.backward(torque);

        // |m × H_eff| where m=(1,0,0), H=(0,H_mag,0) → cross = (0,0,H_mag), norm = H_mag
        let expected = GAMMA.abs() * h_mag;
        let rel_err = (torque.value() - expected).abs() / expected;
        assert!(
            rel_err < 1e-10,
            "LLG torque: got {}, expected {}",
            torque.value(),
            expected
        );
    }

    // 11. Total energy gradient via chain rule
    #[test]
    fn test_total_energy_gradient() {
        // E_total = E_zeeman + E_anisotropy, differentiate w.r.t. mz
        let ku_val = 4.0e4_f64;
        let ms_val = 8.0e5_f64;
        let hz_val = 1.0e5_f64;
        let mz_val = 0.5_f64;

        let tape = Tape::new();
        let mx = Var::leaf(&tape, 0.0);
        let my = Var::leaf(&tape, 0.0);
        let mz = Var::leaf(&tape, mz_val);
        let ms = Var::leaf(&tape, ms_val);
        let ku = Var::leaf(&tape, ku_val);

        let e_z = zeeman_energy_diff(&tape, mx, my, mz, 0.0, 0.0, hz_val, ms);
        let e_a = anisotropy_energy_diff(&tape, mz, ku);
        // Total loss = e_z + e_a  (we'll differentiate the sum)
        let total = e_z + e_a;
        tape.backward(total);

        // ∂E_zeeman/∂mz = -μ₀·ms·Hz
        // ∂E_aniso/∂mz  = -2·Ku·mz
        let expected = -MU_0 * ms_val * hz_val + (-2.0 * ku_val * mz_val);
        let rel_err = (mz.grad() - expected).abs() / expected.abs().max(1.0);
        assert!(
            rel_err < 1e-8,
            "total ∂E/∂mz: got {}, expected {}",
            mz.grad(),
            expected
        );
    }
}
