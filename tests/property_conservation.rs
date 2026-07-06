//! Property-based tests for conservation laws under LLG / Heisenberg dynamics.
//!
//! Uses `proptest` to generate random valid parameters (within physical bounds)
//! and verify that fundamental invariants hold across the parameter space:
//!
//!   - |m| = 1 conserved by LLG under RK4 (within numerical tolerance)
//!   - |m| = 1 conserved by Heun predictor-corrector
//!   - |m| = 1 conserved by Euler (post-normalization)
//!   - Energy conservation by RK4 at alpha = 0 (Larmor precession)
//!   - LLG damping aligns m with H_eff in the long-time limit
//!   - Larmor precession reverses sign under H -> -H
//!   - Vector3 algebraic identities (cross anti-commutativity, triple product cyclic)
//!   - Lagrange identity |a x b|^2 + (a . b)^2 = |a|^2 |b|^2
//!   - Normalization idempotency: v.normalize().magnitude() == 1
//!   - At alpha = 0, dm/dt perpendicular to both m and h_eff
//!
//! Each `proptest!` test runs `cases` random trials. We configure 32 cases per
//! property for fast CI; proptest's shrinking still produces minimal failing
//! examples on regression.

#![allow(clippy::needless_pass_by_value)]
// proptest depends on rusty-fork -> wait-timeout, which has no wasm32 backend
// (no process-fork model on wasm32-unknown-unknown); this suite is native-only.
#![cfg(not(target_arch = "wasm32"))]

use std::f64::consts::PI;

use proptest::prelude::*;
use spintronics::constants::{GAMMA, MU_0};
use spintronics::dynamics::llg::{calc_dm_dt, zeeman_energy, LlgSolver};
use spintronics::vector3::Vector3;

// ---------------------------------------------------------------------------
// Strategy helpers
// ---------------------------------------------------------------------------

/// Marsaglia-style deterministic conversion from `(u64, u64)` to a unit Vector3.
///
/// Maps a pair of uniformly-distributed integers to a unit vector via spherical
/// coordinates `(theta, phi)`. The transformation `cos(theta) ~ U(-1, 1)` plus
/// `phi ~ U(0, 2*pi)` yields uniform sampling on `S^2`.
fn seed_to_unit_vector(seed_a: u64, seed_b: u64) -> Vector3<f64> {
    let u = (seed_a as f64) / (u64::MAX as f64);
    let v = (seed_b as f64) / (u64::MAX as f64);
    let cos_theta = 2.0 * u - 1.0;
    let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
    let phi = 2.0 * PI * v;
    Vector3::new(sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta)
}

/// Strategy: uniformly distributed unit Vector3 on S^2.
fn unit_vector_strategy() -> impl Strategy<Value = Vector3<f64>> {
    (any::<u64>(), any::<u64>()).prop_map(|(a, b)| seed_to_unit_vector(a, b))
}

/// Strategy: an arbitrary (non-unit) Vector3 with components in `[-1, 1]`,
/// guaranteed to have magnitude >= 1e-3 (avoids the degenerate zero vector).
fn nonzero_vector_strategy() -> impl Strategy<Value = Vector3<f64>> {
    (any::<u64>(), any::<u64>(), any::<u64>()).prop_map(|(a, b, c)| {
        let x = (a as f64) / (u64::MAX as f64) * 2.0 - 1.0;
        let y = (b as f64) / (u64::MAX as f64) * 2.0 - 1.0;
        let z = (c as f64) / (u64::MAX as f64) * 2.0 - 1.0;
        let v = Vector3::new(x, y, z);
        // Guarantee non-degeneracy by adding a tiny offset along z.
        if v.magnitude() < 1.0e-3 {
            Vector3::new(v.x, v.y, v.z + 1.0)
        } else {
            v
        }
    })
}

/// Strategy: small valid time step (seconds).
fn dt_strategy() -> impl Strategy<Value = f64> {
    1.0e-15f64..1.0e-13
}

/// Strategy: physically reasonable Gilbert damping.
fn alpha_strategy() -> impl Strategy<Value = f64> {
    0.001f64..0.5
}

/// Strategy: applied field magnitude (Tesla), single component.
fn h_field_strategy() -> impl Strategy<Value = f64> {
    1.0e-3f64..1.0e1
}

/// Strategy: saturation magnetization (A/m).
fn ms_strategy() -> impl Strategy<Value = f64> {
    1.0e4f64..2.0e6
}

/// Strategy: small number of integration steps for tests (keeps runtime modest).
fn n_steps_strategy() -> impl Strategy<Value = usize> {
    5usize..30
}

// ---------------------------------------------------------------------------
// Proptest configuration: 32 cases per property keeps CI fast while still
// exploring the parameter space meaningfully.
// ---------------------------------------------------------------------------

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    /// |m| = 1 conserved by RK4 over short integration windows.
    ///
    /// `step_rk4` renormalises after each step, so this checks the integrator
    /// returns a unit vector even for arbitrary fields / damping / time steps.
    #[test]
    fn norm_conserved_rk4(
        m0 in unit_vector_strategy(),
        alpha in alpha_strategy(),
        h_z in h_field_strategy(),
        dt in dt_strategy(),
        n_steps in n_steps_strategy(),
    ) {
        let mut m = m0;
        let h = Vector3::new(0.0, 0.0, h_z);
        let solver = LlgSolver::new(alpha, dt);
        for _ in 0..n_steps {
            m = solver.step_rk4(m, |_| h);
        }
        let mag = m.magnitude();
        prop_assert!((mag - 1.0).abs() < 1.0e-10,
            "|m| drifted to {} after {} RK4 steps", mag, n_steps);
    }

    /// |m| = 1 conserved by Heun (predictor-corrector) integrator.
    #[test]
    fn norm_conserved_heun(
        m0 in unit_vector_strategy(),
        alpha in alpha_strategy(),
        h_x in h_field_strategy(),
        h_z in h_field_strategy(),
        dt in dt_strategy(),
        n_steps in n_steps_strategy(),
    ) {
        let mut m = m0;
        let h = Vector3::new(h_x, 0.0, h_z);
        let solver = LlgSolver::new(alpha, dt);
        for _ in 0..n_steps {
            m = solver.step_heun(m, |_| h);
        }
        prop_assert!((m.magnitude() - 1.0).abs() < 1.0e-10);
    }

    /// |m| = 1 conserved by Euler (post-normalization).
    #[test]
    fn norm_conserved_euler(
        m0 in unit_vector_strategy(),
        alpha in alpha_strategy(),
        h_z in h_field_strategy(),
        dt in dt_strategy(),
        n_steps in n_steps_strategy(),
    ) {
        let mut m = m0;
        let h = Vector3::new(0.0, 0.0, h_z);
        let solver = LlgSolver::new(alpha, dt);
        for _ in 0..n_steps {
            m = solver.step_euler(m, h);
        }
        prop_assert!((m.magnitude() - 1.0).abs() < 1.0e-10);
    }

    /// Zeeman energy E = -mu_0 ms (m . h) at alpha = 0 should be conserved by
    /// RK4 over short windows (Larmor precession is energy-conserving).
    #[test]
    fn zeeman_energy_conserved_zero_damping(
        m0 in unit_vector_strategy(),
        h_z in h_field_strategy(),
        dt in dt_strategy(),
        ms in ms_strategy(),
        n_steps in n_steps_strategy(),
    ) {
        let h = Vector3::new(0.0, 0.0, h_z);
        let solver = LlgSolver::new(0.0, dt);
        let mut m = m0;
        let e0 = zeeman_energy(m, h, ms);
        for _ in 0..n_steps {
            m = solver.step_rk4(m, |_| h);
        }
        let e1 = zeeman_energy(m, h, ms);
        // Relative tolerance scales with field magnitude; symplectic-like
        // behaviour of normalised RK4 should give |dE/E| << 1e-3 for n <= 30.
        let scale = (MU_0 * ms * h_z).abs().max(1.0e-20);
        prop_assert!((e1 - e0).abs() / scale < 1.0e-2,
            "Energy drift {} over {} steps (scale {})", e1 - e0, n_steps, scale);
    }

    /// With damping alpha > 0, the magnetisation must move toward (m . h_eff)
    /// growing (or staying constant): the projection along the field cannot
    /// decrease over a damping-only test starting from m perpendicular to h.
    #[test]
    fn damping_aligns_m_with_field(
        alpha in 0.05f64..0.3,
        h_z in 1.0f64..5.0,
        dt in 1.0e-13f64..5.0e-13,
    ) {
        let h = Vector3::new(0.0, 0.0, h_z);
        // Start with m perpendicular to h, so initial projection is 0.
        let m0 = Vector3::new(1.0, 0.0, 0.0);
        let solver = LlgSolver::new(alpha, dt);
        let mut m = m0;
        // A few thousand steps of integration to overcome precession period.
        for _ in 0..2000 {
            m = solver.step_rk4(m, |_| h);
        }
        let proj_final = m.dot(&h) / h.magnitude();
        // After many steps, the z-projection must have grown (proj initially 0,
        // damping drives m toward +z since gamma > 0 in this convention).
        prop_assert!(proj_final > 0.0,
            "Damping did not align m with h: m.z/|h| = {}", proj_final);
    }

    /// Larmor precession reverses sign under H -> -H.
    ///
    /// dm/dt (m, h, 0, 0) = -dm/dt (m, -h, 0, 0).
    /// Relative tolerance scaled by `|dm/dt|`.
    #[test]
    fn larmor_reverses_under_field_flip(
        m in unit_vector_strategy(),
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
    ) {
        let h = seed_to_unit_vector(seed_a, seed_b);
        let dm_forward = calc_dm_dt(m, h, GAMMA, 0.0);
        let dm_reversed = calc_dm_dt(m, h * -1.0, GAMMA, 0.0);
        let sum = dm_forward + dm_reversed;
        let scale = dm_forward.magnitude().max(1.0);
        prop_assert!(sum.magnitude() / scale < 1.0e-10,
            "Larmor not antisymmetric in H: |sum|/scale = {}",
            sum.magnitude() / scale);
    }

    /// At alpha = 0, dm/dt is perpendicular to m (norm-preserving).
    ///
    /// Tolerance is relative to `|dm/dt|`, which itself scales like
    /// `GAMMA * |h|` (here `|h| = 1` Tesla, so the natural scale is GAMMA).
    #[test]
    fn dm_dt_perpendicular_to_m_zero_damping(
        m in unit_vector_strategy(),
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
    ) {
        let h = seed_to_unit_vector(seed_a, seed_b);
        let dm = calc_dm_dt(m, h, GAMMA, 0.0);
        let dm_mag = dm.magnitude().max(1.0);
        prop_assert!(dm.dot(&m).abs() / dm_mag < 1.0e-10,
            "dm/dt . m / |dm/dt| = {} not negligible at alpha=0",
            dm.dot(&m) / dm_mag);
    }

    /// At alpha = 0, dm/dt is also perpendicular to h_eff (pure precession,
    /// no work done by the field).
    ///
    /// Uses relative tolerance scaled by `|dm/dt| * |h|`.
    #[test]
    fn dm_dt_perpendicular_to_h_zero_damping(
        m in unit_vector_strategy(),
        seed_a in any::<u64>(),
        seed_b in any::<u64>(),
    ) {
        let h = seed_to_unit_vector(seed_a, seed_b);
        let dm = calc_dm_dt(m, h, GAMMA, 0.0);
        let scale = (dm.magnitude() * h.magnitude()).max(1.0);
        prop_assert!(dm.dot(&h).abs() / scale < 1.0e-10,
            "dm/dt . h / scale = {} not negligible at alpha=0",
            dm.dot(&h) / scale);
    }

    /// Cross product is anti-commutative: a x b = -(b x a).
    #[test]
    fn cross_product_anticommutative(
        a in unit_vector_strategy(),
        b in unit_vector_strategy(),
    ) {
        let lhs = a.cross(&b);
        let rhs = b.cross(&a) * -1.0;
        prop_assert!((lhs.x - rhs.x).abs() < 1.0e-12);
        prop_assert!((lhs.y - rhs.y).abs() < 1.0e-12);
        prop_assert!((lhs.z - rhs.z).abs() < 1.0e-12);
    }

    /// Triple product invariant under cyclic permutation:
    /// a . (b x c) = b . (c x a) = c . (a x b).
    #[test]
    fn triple_product_cyclic(
        a in unit_vector_strategy(),
        b in unit_vector_strategy(),
        c in unit_vector_strategy(),
    ) {
        let lhs = a.cross(&b).dot(&c);
        let mid = b.cross(&c).dot(&a);
        let rhs = c.cross(&a).dot(&b);
        prop_assert!((lhs - mid).abs() < 1.0e-12);
        prop_assert!((mid - rhs).abs() < 1.0e-12);
    }

    /// Lagrange's identity: |a x b|^2 + (a . b)^2 = |a|^2 |b|^2.
    #[test]
    fn lagrange_identity(
        a in unit_vector_strategy(),
        b in unit_vector_strategy(),
    ) {
        let cross_sq = a.cross(&b).magnitude_squared();
        let dot = a.dot(&b);
        let lhs = cross_sq + dot * dot;
        let rhs = a.magnitude_squared() * b.magnitude_squared();
        prop_assert!((lhs - rhs).abs() < 1.0e-12);
    }

    /// Normalize is idempotent: v.normalize().magnitude() ~ 1.
    #[test]
    fn normalize_idempotent(v in nonzero_vector_strategy()) {
        let n = v.normalize();
        prop_assert!((n.magnitude() - 1.0).abs() < 1.0e-12);
        // Applying normalize again should not change the vector.
        let nn = n.normalize();
        prop_assert!((nn.magnitude() - 1.0).abs() < 1.0e-12);
        prop_assert!((nn - n).magnitude() < 1.0e-12);
    }
}
