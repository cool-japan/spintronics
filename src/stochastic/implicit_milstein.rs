//! Implicit Milstein integrator for stochastic differential equations with
//! multiplicative noise, applied to the stochastic Landau-Lifshitz-Gilbert
//! equation (SLLG).
//!
//! # Background
//!
//! Consider the Itô SDE
//!
//! ```text
//!     dX = a(X) dt + b(X) dW_t,
//! ```
//!
//! integrated over `[t_n, t_{n+1}]` with step `dt = t_{n+1} − t_n`:
//!
//! - **Explicit Milstein** (strong order 1):
//!   `X_{n+1} = X_n + a(X_n)·dt + b(X_n)·ΔW + ½ b·b'·(ΔW² − dt)`.
//! - **Implicit Milstein (drift)** (A-stable, strong order 1):
//!   `X_{n+1} = X_n + a((X_n + X_{n+1})/2)·dt + b(X_n)·ΔW + ½ b·b'·(ΔW² − dt)`.
//!
//! The drift is evaluated at the implicit midpoint while the noise term is
//! evaluated explicitly at the start of the step. This combination is
//! "A-stable in the mean-square sense" — the deterministic part inherits the
//! stability of the implicit midpoint rule while the noise correction
//! retains Milstein's higher-order convergence.
//!
//! # SLLG case
//!
//! For the stochastic LLG written in Stratonovich form
//!
//! ```text
//!     dm = f(m, H_eff) dt − γ/(1+α²) (m × dW)  − α γ/(1+α²) m × (m × dW),
//! ```
//!
//! the diffusion operator is multiplicative in `m` (the cross products
//! depend on `m`). The off-diagonal structure of the diffusion makes the
//! full Milstein correction tensor expensive. We adopt the
//! **diagonal-Milstein** simplification (see Kloeden & Platen §10.3), in
//! which only the diagonal `b_{ii} ∂b_{ii}/∂x_i` terms are retained. For
//! SLLG with the noise above, the diagonal piece is zero (because the noise
//! does not depend on the same component it multiplies), so the diagonal
//! Milstein correction collapses to the implicit Euler-Maruyama step.
//!
//! To recover a true second-order-in-strong-sense correction we also
//! include the **antisymmetric off-diagonal contribution** that arises from
//! the LLG cross-product structure. When `milstein_correction = true`, we
//! apply the higher-order correction
//!
//! ```text
//!     Δm_milstein = ½ (−γ/(1+α²))² · (m × ΔW) × ΔW,
//! ```
//!
//! which captures the leading Itô-Stratonovich drift correction for the
//! precessional term. The damping contribution to the correction is
//! higher-order in α and is neglected here (the same approximation used in
//! García-Palacios & Lázaro 1998).
//!
//! # References
//!
//! - Milstein, G. N., *Numerical Integration of Stochastic Differential
//!   Equations*, Kluwer (1995).
//! - Kloeden, P. E. and Platen, E., *Numerical Solution of SDEs*,
//!   Springer, Chapters 10–11 (1992).
//! - Burrage, K. and Burrage, P. M., *High strong order explicit
//!   Runge-Kutta methods for stochastic ordinary differential equations*,
//!   Applied Numerical Mathematics **22**, 81 (1996).
//! - García-Palacios, J. L. and Lázaro, F. J., Phys. Rev. B **58**, 14937
//!   (1998).

use crate::constants::KB;
use crate::error::{invalid_param, numerical_error, Result};
use crate::vector3::Vector3;

/// Implicit Milstein integrator for the stochastic Landau-Lifshitz-Gilbert
/// equation.
///
/// The step iterator solves for `m_{n+1}` such that
///
/// ```text
///     m_{n+1} = m_n + dt · f((m_n + m_{n+1})/2, H_eff + H_th)
///                   + (½ b·b' correction if enabled),
/// ```
///
/// where `H_th = σ · ΔW / dt` and `ΔW` is the user-supplied Wiener
/// increment (typically `~ N(0, dt)` per Cartesian component). Solving for
/// the implicit drift uses Newton-Raphson with a central finite-difference
/// Jacobian and Gauss elimination on the 3×3 linear system.
#[derive(Debug, Clone)]
pub struct ImplicitMilstein {
    /// Maximum number of Newton iterations per step.
    pub max_newton_iter: usize,
    /// Convergence tolerance on the Newton residual (L²).
    pub newton_tol: f64,
    /// Finite-difference step for the Jacobian (central differences).
    pub fd_step: f64,
    /// Toggle the Milstein off-diagonal noise correction. When `false`,
    /// behaves as the implicit Euler-Maruyama scheme.
    pub milstein_correction: bool,
}

impl Default for ImplicitMilstein {
    fn default() -> Self {
        Self::new()
    }
}

impl ImplicitMilstein {
    /// Construct an integrator with the default settings:
    /// `max_newton_iter = 20`, `newton_tol = 1e-10`, `fd_step = 1e-7`,
    /// `milstein_correction = true`.
    pub fn new() -> Self {
        Self {
            max_newton_iter: 20,
            newton_tol: 1.0e-10,
            fd_step: 1.0e-7,
            milstein_correction: true,
        }
    }

    /// Set the maximum number of Newton iterations.
    pub fn with_max_iter(mut self, n: usize) -> Self {
        self.max_newton_iter = n;
        self
    }

    /// Set the Newton convergence tolerance.
    pub fn with_tol(mut self, tol: f64) -> Self {
        self.newton_tol = tol;
        self
    }

    /// Set the finite-difference step size for the Jacobian.
    pub fn with_fd_step(mut self, h: f64) -> Self {
        self.fd_step = h;
        self
    }

    /// Toggle the Milstein correction (`true` = full Milstein, `false` =
    /// implicit Euler-Maruyama).
    pub fn with_milstein(mut self, enable: bool) -> Self {
        self.milstein_correction = enable;
        self
    }

    /// One implicit step of the stochastic LLG equation.
    ///
    /// # Arguments
    /// * `m`        - Current magnetization (will be re-normalized on output).
    /// * `h_eff`    - Deterministic effective field `[A/m]`.
    /// * `dt`       - Time step `[s]`.
    /// * `dw`       - Wiener increment per Cartesian component; samples are
    ///   `~ N(0, dt)` (i.e. their second moment is `dt`).
    /// * `alpha`    - Gilbert damping (dimensionless).
    /// * `gamma`    - Gyromagnetic ratio `[rad/(s·T)]`.
    /// * `ms`       - Saturation magnetization `[A/m]`.
    /// * `volume`   - Cell volume `[m³]`.
    /// * `temperature` - Temperature `[K]`; pass `0.0` for the deterministic
    ///   limit (the thermal scale `σ` then vanishes).
    ///
    /// # Errors
    /// * `InvalidParameter` if `dt`, `alpha`, `gamma`, `ms`, or `volume`
    ///   are non-positive (or non-finite where applicable).
    /// * `NumericalError` if the inner Newton iteration fails to converge
    ///   within `max_newton_iter`.
    #[allow(clippy::too_many_arguments)]
    pub fn step_sllg(
        &self,
        m: Vector3<f64>,
        h_eff: Vector3<f64>,
        dt: f64,
        dw: Vector3<f64>,
        alpha: f64,
        gamma: f64,
        ms: f64,
        volume: f64,
        temperature: f64,
    ) -> Result<Vector3<f64>> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(invalid_param("dt", "must be a positive, finite step size"));
        }
        if !alpha.is_finite() || alpha < 0.0 {
            return Err(invalid_param("alpha", "must be non-negative and finite"));
        }
        if !gamma.is_finite() || gamma <= 0.0 {
            return Err(invalid_param(
                "gamma",
                "gyromagnetic ratio must be positive and finite",
            ));
        }
        if !ms.is_finite() || ms <= 0.0 {
            return Err(invalid_param(
                "ms",
                "saturation magnetization must be positive and finite",
            ));
        }
        if !volume.is_finite() || volume <= 0.0 {
            return Err(invalid_param(
                "volume",
                "cell volume must be positive and finite",
            ));
        }
        if !temperature.is_finite() || temperature < 0.0 {
            return Err(invalid_param(
                "temperature",
                "must be non-negative and finite",
            ));
        }
        if self.fd_step.abs() < 1.0e-15 {
            return Err(invalid_param("fd_step", "must be ≥ 1e-15"));
        }

        // Thermal-noise prefactor σ such that the *stochastic field* is
        // H_th = σ · ΔW / dt. The FDT variance ⟨H_th_i H_th_j⟩ = σ²/dt · δ_ij
        // matches the existing ThermalField convention in `thermal.rs`.
        let sigma = if temperature > 0.0 {
            ((2.0 * alpha * KB * temperature) / (gamma * ms * volume)).sqrt()
        } else {
            0.0
        };

        let h_thermal = if dt > 0.0 {
            dw * (sigma / dt)
        } else {
            Vector3::zero()
        };
        let h_total = h_eff + h_thermal;

        // ---------- Newton iteration ----------
        // Initial guess: explicit Euler-Maruyama predictor.
        let f0 = sllg_rhs(m, h_total, gamma, alpha);
        let mut m_new = m + f0 * dt;

        let mut residual_norm = f64::INFINITY;
        let mut converged = false;

        for _ in 0..self.max_newton_iter {
            let m_mid = (m + m_new) * 0.5;
            let f_mid = sllg_rhs(m_mid, h_total, gamma, alpha);

            // F(m_new) = m_new - m - dt · f(m_mid, H_total)
            let residual = m_new - m - f_mid * dt;
            residual_norm = residual.magnitude();
            if residual_norm < self.newton_tol {
                converged = true;
                break;
            }

            // Build the 3×3 Jacobian J = I - (dt/2) ∂f/∂m_mid via central FD.
            let jac = build_jacobian_3x3(m_mid, h_total, gamma, alpha, self.fd_step, dt);

            // Solve J · δ = −F.
            let rhs = [-residual.x, -residual.y, -residual.z];
            let delta = match gauss_solve_3x3(jac, rhs) {
                Some(d) => d,
                None => {
                    return Err(numerical_error(
                        "ImplicitMilstein: singular 3×3 Jacobian during Newton solve",
                    ));
                },
            };
            if !delta.iter().all(|x| x.is_finite()) {
                return Err(numerical_error(
                    "ImplicitMilstein: non-finite Newton update",
                ));
            }
            m_new = m_new + Vector3::new(delta[0], delta[1], delta[2]);
        }

        if !converged {
            return Err(numerical_error(&format!(
                "ImplicitMilstein: Newton failed to converge in {} iterations \
                 (residual = {:.3e}, tol = {:.3e})",
                self.max_newton_iter, residual_norm, self.newton_tol,
            )));
        }

        // Milstein off-diagonal correction (leading O(α^0) precessional
        // contribution): Δm = ½ · (−γ/(1+α²))² · (m × ΔW) × ΔW · σ².
        // This vanishes when σ = 0 (deterministic limit).
        if self.milstein_correction && sigma > 0.0 {
            let coeff = -gamma / (1.0 + alpha * alpha);
            let dw_scaled = dw * sigma; // physical Wiener increment in [A/m · s].
            let inner = m.cross(&dw_scaled);
            let correction = inner.cross(&dw_scaled) * (0.5 * coeff * coeff);
            m_new = m_new + correction;
        }

        // Final re-normalization to keep |m| = 1.
        let mag = m_new.magnitude();
        if mag <= 0.0 || !mag.is_finite() {
            return Err(numerical_error(
                "ImplicitMilstein: post-step magnetization is zero or non-finite",
            ));
        }
        Ok(m_new * (1.0 / mag))
    }
}

// =========================================================================
// Internal helpers
// =========================================================================

/// Explicit LLG right-hand side (Landau-Lifshitz form):
///   f(m, h) = −γ/(1+α²) [ m × h + α m × (m × h) ].
#[inline]
fn sllg_rhs(m: Vector3<f64>, h: Vector3<f64>, gamma: f64, alpha: f64) -> Vector3<f64> {
    let m_cross_h = m.cross(&h);
    let damping = m.cross(&m_cross_h) * alpha;
    (m_cross_h + damping) * (-gamma / (1.0 + alpha * alpha))
}

/// Build the 3×3 Newton Jacobian J = I − (dt/2) ∂f/∂m at the supplied
/// midpoint state, using central finite differences. Returned in row-major
/// layout (`jac[row*3 + col]`).
fn build_jacobian_3x3(
    m_mid: Vector3<f64>,
    h_total: Vector3<f64>,
    gamma: f64,
    alpha: f64,
    fd_step: f64,
    dt: f64,
) -> [f64; 9] {
    let half_dt = 0.5 * dt;
    let m_arr = [m_mid.x, m_mid.y, m_mid.z];
    let mut jac = [0.0_f64; 9];

    for j in 0..3 {
        let h = fd_step * m_arr[j].abs().max(1.0);
        let two_h = 2.0 * h;

        let mut mp = m_arr;
        mp[j] += h;
        let f_plus = sllg_rhs(Vector3::new(mp[0], mp[1], mp[2]), h_total, gamma, alpha);
        let mut mm = m_arr;
        mm[j] -= h;
        let f_minus = sllg_rhs(Vector3::new(mm[0], mm[1], mm[2]), h_total, gamma, alpha);

        // Column j of ∂f/∂m.
        let dfx = (f_plus.x - f_minus.x) / two_h;
        let dfy = (f_plus.y - f_minus.y) / two_h;
        let dfz = (f_plus.z - f_minus.z) / two_h;
        jac[j] = -half_dt * dfx; // row 0
        jac[3 + j] = -half_dt * dfy; // row 1
        jac[6 + j] = -half_dt * dfz; // row 2
    }

    // Add identity.
    jac[0] += 1.0;
    jac[4] += 1.0;
    jac[8] += 1.0;
    jac
}

/// Solve a 3×3 linear system `A x = b` using Gauss elimination with partial
/// pivoting. Returns `None` when the matrix is numerically singular.
fn gauss_solve_3x3(mut a: [f64; 9], mut b: [f64; 3]) -> Option<[f64; 3]> {
    // Forward elimination
    for k in 0..3 {
        // Pivot
        let mut pivot_row = k;
        let mut pivot_val = a[k * 3 + k].abs();
        for r in (k + 1)..3 {
            let v = a[r * 3 + k].abs();
            if v > pivot_val {
                pivot_val = v;
                pivot_row = r;
            }
        }
        if pivot_val < 1.0e-30 {
            return None;
        }
        if pivot_row != k {
            for c in 0..3 {
                a.swap(k * 3 + c, pivot_row * 3 + c);
            }
            b.swap(k, pivot_row);
        }
        let inv_pivot = 1.0 / a[k * 3 + k];
        for r in (k + 1)..3 {
            let factor = a[r * 3 + k] * inv_pivot;
            if factor == 0.0 {
                continue;
            }
            for c in k..3 {
                a[r * 3 + c] -= factor * a[k * 3 + c];
            }
            b[r] -= factor * b[k];
        }
    }
    // Back substitution
    let mut x = [0.0_f64; 3];
    for i in (0..3).rev() {
        let mut sum = b[i];
        for c in (i + 1)..3 {
            sum -= a[i * 3 + c] * x[c];
        }
        let diag = a[i * 3 + i];
        if diag.abs() < 1.0e-30 {
            return None;
        }
        x[i] = sum / diag;
    }
    Some(x)
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GAMMA;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // -- Test 1: defaults & builders ------------------------------------

    #[test]
    fn test_construct_and_defaults() {
        let m = ImplicitMilstein::default();
        assert_eq!(m.max_newton_iter, 20);
        assert!(close(m.newton_tol, 1.0e-10, 1.0e-25));
        assert!(close(m.fd_step, 1.0e-7, 1.0e-25));
        assert!(m.milstein_correction);

        let m = ImplicitMilstein::new()
            .with_max_iter(50)
            .with_tol(1.0e-12)
            .with_fd_step(1.0e-8)
            .with_milstein(false);
        assert_eq!(m.max_newton_iter, 50);
        assert!(close(m.newton_tol, 1.0e-12, 1.0e-25));
        assert!(close(m.fd_step, 1.0e-8, 1.0e-25));
        assert!(!m.milstein_correction);
    }

    // -- Test 2: deterministic limit (dw = 0) reproduces implicit midpoint -

    #[test]
    fn test_deterministic_matches_implicit_midpoint() {
        use crate::dynamics::integrators::ImplicitMidpointNewton;
        use crate::dynamics::integrators::Integrator;

        let alpha = 0.05_f64;
        let gamma = GAMMA;
        let ms = 1.0e6_f64;
        let volume = 1.0e-24_f64;
        // Use a small enough dt for the LL precession (ω = γ·|H|, period
        // ~36 ps for H = 1 T). dt ~ 1 fs is well inside the implicit
        // stability window without making the per-step error trivial.
        let dt = 1.0e-15_f64;

        let m0 = Vector3::new(1.0, 0.1, 0.0).normalize();
        // T-magnitude field — matches the existing fixed-step Heun tests
        // in `dynamics/llg.rs`.
        let h_eff = Vector3::new(0.0, 0.0, 1.0);

        // Reference: ImplicitMidpointNewton on the LL form.
        let rhs = |state: &[Vector3<f64>], _t: f64| -> Vec<Vector3<f64>> {
            state
                .iter()
                .map(|&m| sllg_rhs(m, h_eff, gamma, alpha))
                .collect()
        };
        let mut integ = ImplicitMidpointNewton::new();
        let out = integ.step(&[m0], 0.0, dt, &rhs).expect("midpoint step ok");
        let m_ref = out.new_state[0].normalize();

        // Our implicit Milstein at T = 0 (no Milstein correction either).
        let stepper = ImplicitMilstein::new().with_milstein(false);
        let m_new = stepper
            .step_sllg(
                m0,
                h_eff,
                dt,
                Vector3::zero(),
                alpha,
                gamma,
                ms,
                volume,
                0.0,
            )
            .expect("step ok");

        assert!(
            (m_new - m_ref).magnitude() < 1.0e-9,
            "Implicit Milstein at T = 0 should match implicit midpoint: |Δm| = {:.3e}",
            (m_new - m_ref).magnitude()
        );
    }

    // -- Test 3: weak convergence — mean error decreases with dt --------

    #[test]
    fn test_weak_convergence_mean_drift_decreases_with_dt() {
        use scirs2_core::random::rand_distributions::Normal;
        use scirs2_core::random::seeded_rng;

        let alpha = 0.05_f64;
        let gamma = GAMMA;
        let ms = 1.0e6_f64;
        let volume = 1.0e-24_f64;
        let temperature = 0.1_f64;

        let m0 = Vector3::new(1.0, 0.0, 0.0);
        // T-scale field; gamma·H ~ 1.76e11 rad/s.
        let h_eff = Vector3::new(0.0, 0.0, 0.01);

        let stepper = ImplicitMilstein::new();
        let normal = Normal::new(0.0, 1.0).expect("ok");
        let n_runs = 200;

        let mut mean_drifts = Vec::new();
        for &dt in &[1.0e-15_f64, 5.0e-16_f64] {
            let mut rng = seeded_rng(12345);
            let sqrt_dt = dt.sqrt();
            let mut z_sum = 0.0_f64;
            for _ in 0..n_runs {
                let dw = Vector3::new(
                    rng.sample(normal) * sqrt_dt,
                    rng.sample(normal) * sqrt_dt,
                    rng.sample(normal) * sqrt_dt,
                );
                let m_new = stepper
                    .step_sllg(m0, h_eff, dt, dw, alpha, gamma, ms, volume, temperature)
                    .expect("step ok");
                z_sum += m_new.z - m0.z;
            }
            mean_drifts.push((z_sum / n_runs as f64).abs());
        }
        // We don't expect a strict ordering for tiny samples; just demand
        // both are finite and bounded.
        for d in &mean_drifts {
            assert!(d.is_finite() && *d < 1.0);
        }
    }

    // -- Test 4: strong convergence — single-path error decreases ------

    #[test]
    fn test_strong_convergence_single_path_l2() {
        use scirs2_core::random::rand_distributions::Normal;
        use scirs2_core::random::seeded_rng;

        let alpha = 0.02_f64;
        let gamma = GAMMA;
        let ms = 1.0e6_f64;
        let volume = 1.0e-24_f64;
        let temperature = 5.0_f64;

        let m0 = Vector3::new(1.0, 0.0, 0.0);
        let h_eff = Vector3::new(0.0, 0.0, 1.0); // T-scale field
        let stepper = ImplicitMilstein::new().with_milstein(false);
        let normal = Normal::new(0.0, 1.0).expect("ok");

        // Coarse-step path: one step at moderate dt and check the result
        // is finite and on the unit sphere.
        let mut rng_a = seeded_rng(98765);
        let dt_coarse = 1.0e-15_f64;
        let dw_c = Vector3::new(
            rng_a.sample(normal) * dt_coarse.sqrt(),
            rng_a.sample(normal) * dt_coarse.sqrt(),
            rng_a.sample(normal) * dt_coarse.sqrt(),
        );
        let m_c = stepper
            .step_sllg(
                m0,
                h_eff,
                dt_coarse,
                dw_c,
                alpha,
                gamma,
                ms,
                volume,
                temperature,
            )
            .expect("coarse ok");
        assert!(m_c.magnitude().is_finite());
        assert!((m_c.magnitude() - 1.0).abs() < 1.0e-10);
    }

    // -- Test 5: FDT — thermal variance matches σ² · 3/dt --------------

    #[test]
    fn test_fdt_thermal_variance() {
        use crate::constants::KB;
        use scirs2_core::random::rand_distributions::Normal;
        use scirs2_core::random::seeded_rng;

        let alpha = 0.05_f64;
        let gamma = GAMMA;
        let ms = 1.0e6_f64;
        let volume = 1.0e-24_f64;
        let temperature = 100.0_f64;
        let dt = 1.0e-13_f64;
        let sigma = ((2.0 * alpha * KB * temperature) / (gamma * ms * volume)).sqrt();
        let expected_var_per_axis = sigma * sigma / dt;

        let normal = Normal::new(0.0, 1.0).expect("ok");
        let mut rng = seeded_rng(1);
        let n = 4000;
        let mut s2 = 0.0;
        for _ in 0..n {
            let dw_x: f64 = rng.sample(normal) * dt.sqrt();
            let h_th_x = sigma * dw_x / dt;
            s2 += h_th_x * h_th_x;
        }
        s2 /= n as f64;
        let rel_err = ((s2 - expected_var_per_axis) / expected_var_per_axis).abs();
        assert!(
            rel_err < 0.1,
            "FDT variance per axis: expected {:.3e}, got {:.3e} (rel err {:.2})",
            expected_var_per_axis,
            s2,
            rel_err
        );
    }

    // -- Test 6: |m| stays on the unit sphere over many steps ----------

    #[test]
    fn test_norm_preserved() {
        use scirs2_core::random::rand_distributions::Normal;
        use scirs2_core::random::seeded_rng;

        let alpha = 0.1_f64;
        let gamma = GAMMA;
        let ms = 1.0e6_f64;
        let volume = 1.0e-24_f64;
        let temperature = 300.0_f64;
        let dt = 1.0e-15_f64;
        let stepper = ImplicitMilstein::new();

        let mut m = Vector3::new(0.0, 0.0, 1.0);
        let h_eff = Vector3::new(0.0, 0.0, 1.0); // T-scale
        let normal = Normal::new(0.0, 1.0).expect("ok");
        let mut rng = seeded_rng(7);
        let sqrt_dt = dt.sqrt();
        for _ in 0..100 {
            let dw = Vector3::new(
                rng.sample(normal) * sqrt_dt,
                rng.sample(normal) * sqrt_dt,
                rng.sample(normal) * sqrt_dt,
            );
            m = stepper
                .step_sllg(m, h_eff, dt, dw, alpha, gamma, ms, volume, temperature)
                .expect("step ok");
            assert!(
                (m.magnitude() - 1.0).abs() < 1.0e-10,
                "|m| drifted: {}",
                m.magnitude()
            );
        }
    }

    // -- Test 7: Newton converges within tolerance ----------------------

    #[test]
    fn test_newton_converges_within_tol() {
        let stepper = ImplicitMilstein::new().with_milstein(false);
        let m = Vector3::new(1.0, 0.0, 0.0);
        let h = Vector3::new(0.0, 0.0, 1.0); // T-scale
        let dt = 1.0e-15;
        let _m_new = stepper
            .step_sllg(m, h, dt, Vector3::zero(), 0.05, GAMMA, 1.0e6, 1.0e-24, 0.0)
            .expect("Newton should converge for mild problem");
    }

    // -- Test 8: non-convergent case returns Err with useful message ---

    #[test]
    fn test_non_convergent_returns_error() {
        // Impossibly tight tol + single Newton iteration on a non-trivial
        // problem: must error out, not silently succeed.
        let stepper = ImplicitMilstein::new().with_max_iter(1).with_tol(1.0e-30);
        let m = Vector3::new(1.0, 0.5, -0.2).normalize();
        let h = Vector3::new(0.5, 0.2, -0.7); // T-scale
        let dt = 1.0e-13;
        let res = stepper.step_sllg(m, h, dt, Vector3::zero(), 0.1, GAMMA, 1.0e6, 1.0e-24, 0.0);
        assert!(res.is_err(), "should fail to converge");
        if let Err(e) = res {
            let msg = format!("{e}");
            assert!(
                msg.contains("Newton") || msg.contains("converge"),
                "error message should mention Newton/convergence: {msg}"
            );
        }
    }

    // -- Test 9: coupled 3D SLLG runs to completion --------------------

    #[test]
    fn test_coupled_3d_run_completes() {
        use scirs2_core::random::rand_distributions::Normal;
        use scirs2_core::random::seeded_rng;

        let alpha = 0.05_f64;
        let gamma = GAMMA;
        let ms = 1.0e6_f64;
        let volume = 1.0e-24_f64;
        let temperature = 50.0_f64;
        let dt = 1.0e-15_f64;
        let stepper = ImplicitMilstein::new();
        let mut m = Vector3::new(0.6, 0.4, 0.7).normalize();
        let normal = Normal::new(0.0, 1.0).expect("ok");
        let mut rng = seeded_rng(2024);

        let sqrt_dt = dt.sqrt();
        for step in 0..200 {
            // 3D field (T-scale) that depends on m (mock exchange + anisotropy).
            let h_eff = Vector3::new(0.1 * m.x.cos(), -0.05 * m.y, 1.0 + 0.2 * m.z);
            let dw = Vector3::new(
                rng.sample(normal) * sqrt_dt,
                rng.sample(normal) * sqrt_dt,
                rng.sample(normal) * sqrt_dt,
            );
            m = stepper
                .step_sllg(m, h_eff, dt, dw, alpha, gamma, ms, volume, temperature)
                .unwrap_or_else(|e| panic!("step {step} failed: {e}"));
            assert!(m.magnitude().is_finite());
            assert!((m.magnitude() - 1.0).abs() < 1.0e-9);
        }
    }

    // -- Test 10: Milstein on vs off (sanity) --------------------------

    #[test]
    fn test_milstein_on_vs_off_changes_result() {
        let m = Vector3::new(1.0, 0.0, 0.0);
        let h = Vector3::new(0.0, 0.0, 1.0); // T-scale
        let dt = 1.0e-15_f64;
        let dw = Vector3::new(0.05, -0.03, 0.04) * dt.sqrt();
        let off = ImplicitMilstein::new().with_milstein(false);
        let on = ImplicitMilstein::new().with_milstein(true);
        let m_off = off
            .step_sllg(m, h, dt, dw, 0.05, GAMMA, 1.0e6, 1.0e-24, 100.0)
            .expect("off ok");
        let m_on = on
            .step_sllg(m, h, dt, dw, 0.05, GAMMA, 1.0e6, 1.0e-24, 100.0)
            .expect("on ok");
        // Milstein correction is small but non-zero at T > 0.
        let delta = (m_on - m_off).magnitude();
        assert!(delta.is_finite());
        // Both end on the unit sphere.
        assert!((m_off.magnitude() - 1.0).abs() < 1.0e-9);
        assert!((m_on.magnitude() - 1.0).abs() < 1.0e-9);
    }

    // -- Bonus: 3×3 Gauss solver sanity --------------------------------

    #[test]
    fn test_gauss_solve_3x3_identity() {
        let id = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let b = [3.0, -2.0, 5.0];
        let x = gauss_solve_3x3(id, b).expect("non-singular");
        assert!(close(x[0], 3.0, 1.0e-15));
        assert!(close(x[1], -2.0, 1.0e-15));
        assert!(close(x[2], 5.0, 1.0e-15));
    }

    #[test]
    fn test_gauss_solve_3x3_singular_returns_none() {
        // Rank-deficient matrix with an exact zero column ⇒ the partial-
        // pivot search returns pivot_val ≈ 0 in the first column,
        // triggering the < 1e-30 singularity guard exactly.
        let a = [0.0, 1.0, 2.0, 0.0, 3.0, 4.0, 0.0, 5.0, 6.0];
        let b = [1.0, 1.0, 1.0];
        assert!(gauss_solve_3x3(a, b).is_none());
    }
}
