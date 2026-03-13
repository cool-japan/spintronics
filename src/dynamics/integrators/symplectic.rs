//! Structure-preserving symplectic integrators: Velocity Verlet, Yoshida 4th-order, Forest-Ruth.

use super::rhs_fn::{check_nan, Integrator, IntegratorOutput, RhsFn};
use crate::error::{Error, Result};
use crate::vector3::Vector3;

// =========================================================================
// 3a. Velocity Verlet (2nd-order symplectic)
// =========================================================================

/// Velocity Verlet integrator (2nd-order symplectic / Stormer-Verlet).
///
/// For a system written in the form  dq/dt = p,  dp/dt = F(q, t),
/// this integrator preserves the symplectic structure and provides excellent
/// long-time energy conservation. The state is treated as pairs:
/// even indices are positions (q), odd indices are momenta (p).
///
/// In magnetization dynamics the "position" and "momentum" interpretation
/// maps to the two components of the precession equation when cast in
/// Hamiltonian form.
pub struct VelocityVerlet {
    _private: (),
}

impl VelocityVerlet {
    /// Create a new Velocity Verlet integrator.
    pub fn new() -> Self {
        Self { _private: () }
    }
}

impl Default for VelocityVerlet {
    fn default() -> Self {
        Self::new()
    }
}

impl Integrator for VelocityVerlet {
    fn step(
        &mut self,
        state: &[Vector3<f64>],
        t: f64,
        dt: f64,
        f: &RhsFn<'_>,
    ) -> Result<IntegratorOutput> {
        if state.len() % 2 != 0 {
            return Err(Error::DimensionMismatch {
                expected: "even number of state components (q, p pairs)".to_string(),
                actual: format!("{}", state.len()),
            });
        }

        let half = state.len() / 2;
        let a_old = f(state, t);

        // Half-step velocity: v(t + dt/2) = v(t) + a(t) * dt/2
        // Full-step position: q(t + dt)   = q(t) + v(t + dt/2) * dt
        let mut new_state = state.to_vec();

        for i in 0..half {
            // velocity half-step
            new_state[half + i] = state[half + i] + a_old[half + i] * (dt * 0.5);
            // position full-step
            new_state[i] = state[i] + new_state[half + i] * dt;
        }

        // Evaluate acceleration at new position
        let a_new = f(&new_state, t + dt);

        // Complete velocity step
        for i in 0..half {
            new_state[half + i] = new_state[half + i] + a_new[half + i] * (dt * 0.5);
        }

        check_nan(&new_state)?;

        Ok(IntegratorOutput {
            new_state,
            error_estimate: None,
            suggested_dt: None,
        })
    }
}

// =========================================================================
// 3b. Yoshida 4th-order symplectic integrator
// =========================================================================

/// Yoshida 4th-order symplectic integrator.
///
/// A composition method that applies the leapfrog (Stormer-Verlet) integrator
/// with specially chosen sub-steps to achieve 4th-order accuracy while
/// preserving the symplectic structure.
///
/// Reference: H. Yoshida, "Construction of higher order symplectic
/// integrators", Phys. Lett. A 150, 262 (1990).
///
/// State layout: first half = positions (q), second half = momenta (p).
pub struct Yoshida4 {
    _private: (),
}

impl Yoshida4 {
    /// Create a new Yoshida 4th-order integrator.
    pub fn new() -> Self {
        Self { _private: () }
    }
}

impl Default for Yoshida4 {
    fn default() -> Self {
        Self::new()
    }
}

// Yoshida 4th-order coefficients
const YOSHIDA_W1: f64 = 1.351_207_191_959_657_6;
const YOSHIDA_W0: f64 = -1.702_414_383_919_315;

impl Integrator for Yoshida4 {
    fn step(
        &mut self,
        state: &[Vector3<f64>],
        t: f64,
        dt: f64,
        f: &RhsFn<'_>,
    ) -> Result<IntegratorOutput> {
        if state.len() % 2 != 0 {
            return Err(Error::DimensionMismatch {
                expected: "even number of state components (q, p pairs)".to_string(),
                actual: format!("{}", state.len()),
            });
        }

        // Yoshida 4th-order symmetric composition of leapfrog:
        // Pattern: drift-kick-drift-kick-drift-kick-drift (PQPQPQP)
        // Position drift coefficients: c1=c4=w1/2, c2=c3=(w0+w1)/2
        // Momentum kick coefficients: d1=d3=w1, d2=w0
        let c1 = YOSHIDA_W1 / 2.0;
        let c2 = (YOSHIDA_W0 + YOSHIDA_W1) / 2.0;
        let c3 = c2;
        let c4 = c1;
        let d1 = YOSHIDA_W1;
        let d2 = YOSHIDA_W0;
        let d3 = YOSHIDA_W1;

        let c_coeffs = [c1, c2, c3, c4];
        let d_coeffs = [d1, d2, d3];

        let half = state.len() / 2;
        let mut cur = state.to_vec();
        let mut time = t;

        // First drift
        for i in 0..half {
            cur[i] = cur[i] + cur[half + i] * (c_coeffs[0] * dt);
        }
        time += c_coeffs[0] * dt;

        for step_idx in 0..3 {
            // Kick: p += d_i * dt * F(q, t)
            let accel = f(&cur, time);
            for i in 0..half {
                cur[half + i] = cur[half + i] + accel[half + i] * (d_coeffs[step_idx] * dt);
            }

            // Drift: q += c_{i+1} * dt * p
            for i in 0..half {
                cur[i] = cur[i] + cur[half + i] * (c_coeffs[step_idx + 1] * dt);
            }
            time += c_coeffs[step_idx + 1] * dt;
        }

        check_nan(&cur)?;

        Ok(IntegratorOutput {
            new_state: cur,
            error_estimate: None,
            suggested_dt: None,
        })
    }
}

// =========================================================================
// 3c. Forest-Ruth 4th-order symplectic integrator
// =========================================================================

/// Forest-Ruth 4th-order symplectic integrator.
///
/// An alternative 4th-order composition method with slightly different
/// stability properties than Yoshida. Uses 4 stages.
///
/// Reference: E. Forest and R. D. Ruth, "Fourth-order symplectic
/// integration", Physica D 43, 105 (1990).
///
/// State layout: first half = positions (q), second half = momenta (p).
pub struct ForestRuth {
    _private: (),
}

impl ForestRuth {
    /// Create a new Forest-Ruth 4th-order integrator.
    pub fn new() -> Self {
        Self { _private: () }
    }
}

impl Default for ForestRuth {
    fn default() -> Self {
        Self::new()
    }
}

// Forest-Ruth coefficients
const FR_THETA: f64 = 1.351_207_191_959_657_6; // 2^(1/3) / (2 - 2^(1/3))

impl Integrator for ForestRuth {
    fn step(
        &mut self,
        state: &[Vector3<f64>],
        t: f64,
        dt: f64,
        f: &RhsFn<'_>,
    ) -> Result<IntegratorOutput> {
        if state.len() % 2 != 0 {
            return Err(Error::DimensionMismatch {
                expected: "even number of state components (q, p pairs)".to_string(),
                actual: format!("{}", state.len()),
            });
        }

        let half = state.len() / 2;
        let theta = FR_THETA;
        let mut cur = state.to_vec();

        // Step 1: x += theta * dt/2 * v
        for i in 0..half {
            cur[i] = cur[i] + cur[half + i] * (theta * dt * 0.5);
        }

        // Step 2: v += theta * dt * F
        let a1 = f(&cur, t + theta * dt * 0.5);
        for i in 0..half {
            cur[half + i] = cur[half + i] + a1[half + i] * (theta * dt);
        }

        // Step 3: x += (1 - theta) * dt/2 * v
        for i in 0..half {
            cur[i] = cur[i] + cur[half + i] * ((1.0 - theta) * dt * 0.5);
        }

        // Step 4: v += (1 - 2*theta) * dt * F
        let a2 = f(&cur, t + dt * 0.5);
        for i in 0..half {
            cur[half + i] = cur[half + i] + a2[half + i] * ((1.0 - 2.0 * theta) * dt);
        }

        // Step 5: x += (1 - theta) * dt/2 * v
        for i in 0..half {
            cur[i] = cur[i] + cur[half + i] * ((1.0 - theta) * dt * 0.5);
        }

        // Step 6: v += theta * dt * F
        let a3 = f(&cur, t + (1.0 - theta * 0.5) * dt);
        for i in 0..half {
            cur[half + i] = cur[half + i] + a3[half + i] * (theta * dt);
        }

        // Step 7: x += theta * dt/2 * v
        for i in 0..half {
            cur[i] = cur[i] + cur[half + i] * (theta * dt * 0.5);
        }

        check_nan(&cur)?;

        Ok(IntegratorOutput {
            new_state: cur,
            error_estimate: None,
            suggested_dt: None,
        })
    }
}
