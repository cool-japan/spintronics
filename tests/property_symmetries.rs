//! Property-based tests for SO(3) symmetries and physical invariances.
//!
//! Coverage:
//!   - Rodrigues rotation matrices are orthogonal: R R^T = I.
//!   - Rotation determinant equals +1 (proper rotation).
//!   - Dot product is SO(3)-invariant: (R a) . (R b) = a . b.
//!   - Magnitude is SO(3)-invariant: |R a| = |a|.
//!   - Cross product is SO(3)-equivariant: R (a x b) = (R a) x (R b).
//!   - Heisenberg exchange energy E_ex = -J (m1 . m2) is SO(3)-invariant.
//!   - Uniaxial anisotropy E_a = -K (m . n)^2 is invariant when m and n co-rotate.
//!   - Zeeman energy E_Z = -mu_0 ms (m . H) is invariant when m and H co-rotate.
//!   - calc_dm_dt is SO(3)-equivariant: R . dm/dt(m, h) = dm/dt(R m, R h).
//!   - Dot product is parity-even: (-m1) . (-m2) = m1 . m2.
//!   - Time-reversal: dm/dt(m, h) = -dm/dt(-m, -h) at alpha = 0.
//!   - Linearity in field: dm/dt(m, h1 + h2) = dm/dt(m, h1) + dm/dt(m, h2).
//!   - Rotation composition: (R1 R2) v = R1 (R2 v).
//!   - Rotation inverse: R^T R v = v.
//!   - Rotation trace bound: tr(R) in [-1, 3].

#![allow(clippy::needless_pass_by_value)]
#![allow(clippy::needless_range_loop)]

use std::f64::consts::PI;

use proptest::prelude::*;
use spintronics::constants::{GAMMA, MU_0};
use spintronics::dynamics::llg::{anisotropy_energy, calc_dm_dt, zeeman_energy};
use spintronics::vector3::Vector3;

// ---------------------------------------------------------------------------
// Strategy helpers
// ---------------------------------------------------------------------------

/// Deterministic mapping from two `u64` seeds to a unit Vector3 on `S^2`.
fn seed_to_unit_vector(seed_a: u64, seed_b: u64) -> Vector3<f64> {
    let u = (seed_a as f64) / (u64::MAX as f64);
    let v = (seed_b as f64) / (u64::MAX as f64);
    let cos_theta = 2.0 * u - 1.0;
    let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
    let phi = 2.0 * PI * v;
    Vector3::new(sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta)
}

/// Strategy: unit vector uniformly distributed on `S^2`.
fn unit_vector_strategy() -> impl Strategy<Value = Vector3<f64>> {
    (any::<u64>(), any::<u64>()).prop_map(|(a, b)| seed_to_unit_vector(a, b))
}

/// Rodrigues' rotation formula: build a 3x3 rotation matrix from an axis-angle.
fn rodrigues_rotation(axis: Vector3<f64>, angle: f64) -> [[f64; 3]; 3] {
    let c = angle.cos();
    let s = angle.sin();
    let one_m_c = 1.0 - c;
    let ux = axis.x;
    let uy = axis.y;
    let uz = axis.z;
    [
        [
            c + ux * ux * one_m_c,
            ux * uy * one_m_c - uz * s,
            ux * uz * one_m_c + uy * s,
        ],
        [
            uy * ux * one_m_c + uz * s,
            c + uy * uy * one_m_c,
            uy * uz * one_m_c - ux * s,
        ],
        [
            uz * ux * one_m_c - uy * s,
            uz * uy * one_m_c + ux * s,
            c + uz * uz * one_m_c,
        ],
    ]
}

/// Apply a 3x3 matrix to a Vector3.
fn apply_rotation(r: &[[f64; 3]; 3], v: Vector3<f64>) -> Vector3<f64> {
    Vector3::new(
        r[0][0] * v.x + r[0][1] * v.y + r[0][2] * v.z,
        r[1][0] * v.x + r[1][1] * v.y + r[1][2] * v.z,
        r[2][0] * v.x + r[2][1] * v.y + r[2][2] * v.z,
    )
}

/// Compute the transpose of a 3x3 matrix.
fn transpose(r: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [r[0][0], r[1][0], r[2][0]],
        [r[0][1], r[1][1], r[2][1]],
        [r[0][2], r[1][2], r[2][2]],
    ]
}

/// Multiply two 3x3 matrices: `a b`.
fn matmul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut sum = 0.0;
            for k in 0..3 {
                sum += a[i][k] * b[k][j];
            }
            out[i][j] = sum;
        }
    }
    out
}

/// Determinant of a 3x3 matrix via cofactor expansion along the first row.
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Strategy: a uniformly distributed rotation matrix in SO(3).
///
/// Sampled by uniform axis on `S^2` plus angle in `[0, 2 pi]`. Not the Haar
/// measure on SO(3), but sufficient to explore the manifold for property tests.
fn rotation_matrix_strategy() -> impl Strategy<Value = [[f64; 3]; 3]> {
    (any::<u64>(), any::<u64>(), any::<u64>()).prop_map(|(a, b, c)| {
        let axis = seed_to_unit_vector(a, b);
        let angle = 2.0 * PI * ((c as f64) / (u64::MAX as f64));
        rodrigues_rotation(axis, angle)
    })
}

// ---------------------------------------------------------------------------
// Proptest configuration: 32 cases per property keeps CI fast.
// ---------------------------------------------------------------------------

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    /// Rodrigues rotation matrix is orthogonal: R R^T = I.
    #[test]
    fn rotation_matrix_orthogonal(r in rotation_matrix_strategy()) {
        let rt = transpose(&r);
        let prod = matmul(&r, &rt);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                prop_assert!((prod[i][j] - expected).abs() < 1.0e-12,
                    "R R^T[{}][{}] = {}, expected {}", i, j, prod[i][j], expected);
            }
        }
    }

    /// Determinant of a proper rotation matrix equals +1.
    #[test]
    fn rotation_determinant_unity(r in rotation_matrix_strategy()) {
        let d = det3(&r);
        prop_assert!((d - 1.0).abs() < 1.0e-12,
            "det(R) = {}, expected 1.0", d);
    }

    /// Dot product is invariant under simultaneous rotation:
    /// (R a) . (R b) = a . b.
    #[test]
    fn dot_product_so3_invariant(
        a in unit_vector_strategy(),
        b in unit_vector_strategy(),
        r in rotation_matrix_strategy(),
    ) {
        let lhs = apply_rotation(&r, a).dot(&apply_rotation(&r, b));
        let rhs = a.dot(&b);
        prop_assert!((lhs - rhs).abs() < 1.0e-12);
    }

    /// Magnitude is invariant under rotation: |R a| = |a|.
    #[test]
    fn magnitude_so3_invariant(a in unit_vector_strategy(), r in rotation_matrix_strategy()) {
        let rotated = apply_rotation(&r, a);
        prop_assert!((rotated.magnitude() - a.magnitude()).abs() < 1.0e-12);
    }

    /// Cross product is equivariant under proper rotation:
    /// R (a x b) = (R a) x (R b).
    #[test]
    fn cross_product_so3_equivariant(
        a in unit_vector_strategy(),
        b in unit_vector_strategy(),
        r in rotation_matrix_strategy(),
    ) {
        let lhs = apply_rotation(&r, a.cross(&b));
        let rhs = apply_rotation(&r, a).cross(&apply_rotation(&r, b));
        prop_assert!((lhs - rhs).magnitude() < 1.0e-11);
    }

    /// Heisenberg exchange energy E_ex = -J (m1 . m2) is SO(3)-invariant.
    #[test]
    fn exchange_energy_so3_invariant(
        m1 in unit_vector_strategy(),
        m2 in unit_vector_strategy(),
        r in rotation_matrix_strategy(),
        seed_j in any::<u64>(),
    ) {
        // Exchange constant J spans [1e-22, 1e-20] J (typical magnitudes for
        // a localised Heisenberg interaction).
        let j_couple = 1.0e-22 + ((seed_j as f64) / (u64::MAX as f64)) * 9.9e-21;
        let e_before = -j_couple * m1.dot(&m2);
        let m1_rot = apply_rotation(&r, m1);
        let m2_rot = apply_rotation(&r, m2);
        let e_after = -j_couple * m1_rot.dot(&m2_rot);
        prop_assert!((e_before - e_after).abs() / j_couple.abs() < 1.0e-10);
    }

    /// Uniaxial anisotropy E_a = -K (m . n)^2 is invariant when both m and n
    /// rotate together.
    #[test]
    fn uniaxial_anisotropy_invariant(
        m in unit_vector_strategy(),
        n in unit_vector_strategy(),
        r in rotation_matrix_strategy(),
        seed_k in any::<u64>(),
    ) {
        // Anisotropy constant K spans [1e4, 1e6] J/m^3 (typical hard magnets).
        let k_anis = 1.0e4 + ((seed_k as f64) / (u64::MAX as f64)) * 9.9e5;
        let e_before = anisotropy_energy(m, n, k_anis);
        let m_rot = apply_rotation(&r, m);
        let n_rot = apply_rotation(&r, n);
        let e_after = anisotropy_energy(m_rot, n_rot, k_anis);
        prop_assert!((e_before - e_after).abs() / k_anis < 1.0e-10);
    }

    /// Zeeman energy E_Z = -mu_0 ms (m . H) is invariant when m and H co-rotate.
    #[test]
    fn zeeman_energy_so3_invariant(
        m in unit_vector_strategy(),
        h_dir in unit_vector_strategy(),
        r in rotation_matrix_strategy(),
        seed_h in any::<u64>(),
        seed_ms in any::<u64>(),
    ) {
        let h_mag = 0.01 + ((seed_h as f64) / (u64::MAX as f64)) * 10.0;
        let ms = 1.0e4 + ((seed_ms as f64) / (u64::MAX as f64)) * 1.99e6;
        let h = h_dir * h_mag;
        let e_before = zeeman_energy(m, h, ms);
        let m_rot = apply_rotation(&r, m);
        let h_rot = apply_rotation(&r, h);
        let e_after = zeeman_energy(m_rot, h_rot, ms);
        let scale = (MU_0 * ms * h_mag).abs().max(1.0e-20);
        prop_assert!((e_before - e_after).abs() / scale < 1.0e-10);
    }

    /// LLG right-hand side is SO(3)-equivariant:
    /// R . dm/dt(m, h, gamma, alpha) = dm/dt(R m, R h, gamma, alpha).
    #[test]
    fn llg_rhs_so3_equivariant(
        m in unit_vector_strategy(),
        h_dir in unit_vector_strategy(),
        r in rotation_matrix_strategy(),
        alpha in 0.0f64..0.5,
        seed_h in any::<u64>(),
    ) {
        let h_mag = 0.01 + ((seed_h as f64) / (u64::MAX as f64)) * 10.0;
        let h = h_dir * h_mag;
        let dm_orig = calc_dm_dt(m, h, GAMMA, alpha);
        let dm_orig_rot = apply_rotation(&r, dm_orig);
        let dm_rot = calc_dm_dt(apply_rotation(&r, m), apply_rotation(&r, h), GAMMA, alpha);
        // dm/dt scales like GAMMA * |H|, so normalise tolerance accordingly.
        let scale = (GAMMA * h_mag).abs().max(1.0);
        prop_assert!((dm_orig_rot - dm_rot).magnitude() / scale < 1.0e-9,
            "LLG not SO(3)-equivariant: residual {} (scale {})",
            (dm_orig_rot - dm_rot).magnitude(), scale);
    }

    /// Parity: dot product is even under simultaneous reflection
    /// (m1 -> -m1, m2 -> -m2).
    #[test]
    fn dot_product_parity_even(
        a in unit_vector_strategy(),
        b in unit_vector_strategy(),
    ) {
        let original = a.dot(&b);
        let parity = (a * -1.0).dot(&(b * -1.0));
        prop_assert!((original - parity).abs() < 1.0e-12);
    }

    /// Time-reversal at alpha = 0: `dm/dt(m, h) = dm/dt(-m, -h)`.
    ///
    /// Under time reversal both `m -> -m` and `h -> -h` (axial vectors flip
    /// sign). The Larmor equation `-gamma (m x h)` transforms as
    /// `-gamma ((-m) x (-h)) = -gamma (m x h)`, so the RHS is invariant.
    /// Combined with `d(-m)/d(-t) = dm/dt`, this gives the equality above.
    /// (At nonzero alpha the damping term `m x (m x h)` picks up an overall
    /// sign change, so the dissipative LLG is not time-reversal symmetric.)
    #[test]
    fn llg_time_reversal_zero_damping(
        m in unit_vector_strategy(),
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
    ) {
        let h = seed_to_unit_vector(seed_a, seed_b);
        let dm = calc_dm_dt(m, h, GAMMA, 0.0);
        let dm_reversed = calc_dm_dt(m * -1.0, h * -1.0, GAMMA, 0.0);
        let diff = dm - dm_reversed;
        let scale = dm.magnitude().max(1.0);
        prop_assert!(diff.magnitude() / scale < 1.0e-10,
            "Time-reversal symmetry broken at alpha=0: |dm - dm_rev|/scale = {}",
            diff.magnitude() / scale);
    }

    /// Linearity in field at alpha = 0:
    /// dm/dt(m, h1 + h2) = dm/dt(m, h1) + dm/dt(m, h2).
    ///
    /// (For alpha > 0, the LLG RHS is non-linear in h due to the m x (m x h)
    /// term, which is linear in h, so this property actually holds at all alpha.)
    #[test]
    fn llg_linearity_in_field(
        m in unit_vector_strategy(),
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
        seed_c in any::<u64>(),
        seed_d in any::<u64>(),
        alpha in 0.0f64..0.3,
    ) {
        let h1 = seed_to_unit_vector(seed_a, seed_b);
        let h2 = seed_to_unit_vector(seed_c, seed_d);
        let dm_sum = calc_dm_dt(m, h1 + h2, GAMMA, alpha);
        let dm1 = calc_dm_dt(m, h1, GAMMA, alpha);
        let dm2 = calc_dm_dt(m, h2, GAMMA, alpha);
        let diff = dm_sum - dm1 - dm2;
        let scale = GAMMA.abs().max(1.0);
        prop_assert!(diff.magnitude() / scale < 1.0e-9,
            "LLG not linear in H: residual {}", diff.magnitude());
    }

    /// Rotation composition: (R1 R2) v = R1 (R2 v).
    #[test]
    fn rotation_composition(
        r1 in rotation_matrix_strategy(),
        r2 in rotation_matrix_strategy(),
        v in unit_vector_strategy(),
    ) {
        let r12 = matmul(&r1, &r2);
        let lhs = apply_rotation(&r12, v);
        let rhs = apply_rotation(&r1, apply_rotation(&r2, v));
        prop_assert!((lhs - rhs).magnitude() < 1.0e-12);
    }

    /// R^T R v = v: orthogonal matrices satisfy R^{-1} = R^T.
    #[test]
    fn rotation_inverse(
        r in rotation_matrix_strategy(),
        v in unit_vector_strategy(),
    ) {
        let rt = transpose(&r);
        let v_rot = apply_rotation(&r, v);
        let v_back = apply_rotation(&rt, v_rot);
        prop_assert!((v_back - v).magnitude() < 1.0e-12);
    }

    /// Trace of a rotation matrix lies in `[-1, 3]` and equals `1 + 2 cos(theta)`,
    /// where `theta` is the rotation angle.
    #[test]
    fn rotation_trace_bounded(r in rotation_matrix_strategy()) {
        let tr = r[0][0] + r[1][1] + r[2][2];
        let bounds = (-1.0 - 1.0e-12)..=(3.0 + 1.0e-12);
        prop_assert!(bounds.contains(&tr), "tr(R) = {} outside [-1, 3]", tr);
    }
}
