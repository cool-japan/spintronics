//! Dormand-Prince embedded Runge-Kutta integrators: RK5(4) and RK8(7).

use super::rhs_fn::{
    check_nan, max_error_norm, vec_add_scaled, vec_combine, Integrator, IntegratorOutput, RhsFn,
};
use crate::error::Result;
use crate::vector3::Vector3;

// =========================================================================
// 1. Dormand-Prince RK5(4) — the classic "ode45"
// =========================================================================

/// Dormand-Prince 5(4) embedded Runge-Kutta integrator.
///
/// This is the method underlying MATLAB's `ode45`. It uses 7 stages with the
/// First-Same-As-Last (FSAL) property: the last stage evaluation of one step
/// becomes the first stage of the next step, effectively giving 6 new function
/// evaluations per step.
///
/// The embedded 4th-order solution provides a local error estimate that drives
/// adaptive step size selection.
///
/// # Butcher Tableau
///
/// Uses the standard Dormand-Prince coefficients (1980).
pub struct DormandPrince45 {
    /// Last stage evaluation from previous step (FSAL reuse).
    fsal_k7: Option<Vec<Vector3<f64>>>,
}

impl DormandPrince45 {
    /// Create a new Dormand-Prince 5(4) integrator.
    pub fn new() -> Self {
        Self { fsal_k7: None }
    }
}

impl Default for DormandPrince45 {
    fn default() -> Self {
        Self::new()
    }
}

// Butcher tableau nodes
const DP45_C2: f64 = 1.0 / 5.0;
const DP45_C3: f64 = 3.0 / 10.0;
const DP45_C4: f64 = 4.0 / 5.0;
const DP45_C5: f64 = 8.0 / 9.0;
// c6 = 1, c7 = 1

// Butcher tableau coefficients (a_ij)
const DP45_A21: f64 = 1.0 / 5.0;
const DP45_A31: f64 = 3.0 / 40.0;
const DP45_A32: f64 = 9.0 / 40.0;
const DP45_A41: f64 = 44.0 / 45.0;
const DP45_A42: f64 = -56.0 / 15.0;
const DP45_A43: f64 = 32.0 / 9.0;
const DP45_A51: f64 = 19372.0 / 6561.0;
const DP45_A52: f64 = -25360.0 / 2187.0;
const DP45_A53: f64 = 64448.0 / 6561.0;
const DP45_A54: f64 = -212.0 / 729.0;
const DP45_A61: f64 = 9017.0 / 3168.0;
const DP45_A62: f64 = -355.0 / 33.0;
const DP45_A63: f64 = 46732.0 / 5247.0;
const DP45_A64: f64 = 49.0 / 176.0;
const DP45_A65: f64 = -5103.0 / 18656.0;
// Row 7 of the Butcher tableau (identical to B weights due to FSAL property)
const _DP45_A71: f64 = 35.0 / 384.0;
// a72 = 0
const _DP45_A73: f64 = 500.0 / 1113.0;
const _DP45_A74: f64 = 125.0 / 192.0;
const _DP45_A75: f64 = -2187.0 / 6784.0;
const _DP45_A76: f64 = 11.0 / 84.0;

// 5th-order weights (same as a7i row due to FSAL)
const DP45_B1: f64 = 35.0 / 384.0;
// b2 = 0
const DP45_B3: f64 = 500.0 / 1113.0;
const DP45_B4: f64 = 125.0 / 192.0;
const DP45_B5: f64 = -2187.0 / 6784.0;
const DP45_B6: f64 = 11.0 / 84.0;
// b7 = 0

// 4th-order weights (for error estimation)
const DP45_E1: f64 = 5179.0 / 57600.0;
// e2 = 0
const DP45_E3: f64 = 7571.0 / 16695.0;
const DP45_E4: f64 = 393.0 / 640.0;
const DP45_E5: f64 = -92097.0 / 339200.0;
const DP45_E6: f64 = 187.0 / 2100.0;
const DP45_E7: f64 = 1.0 / 40.0;

impl Integrator for DormandPrince45 {
    fn step(
        &mut self,
        state: &[Vector3<f64>],
        t: f64,
        dt: f64,
        f: &RhsFn<'_>,
    ) -> Result<IntegratorOutput> {
        let n = state.len();

        // Stage 1: reuse FSAL if available
        let k1 = match self.fsal_k7.take() {
            Some(k) if k.len() == n => k,
            _ => f(state, t),
        };

        // Stage 2
        let y2 = vec_add_scaled(state, &k1, dt * DP45_A21);
        let k2 = f(&y2, t + dt * DP45_C2);

        // Stage 3
        let y3 = vec_combine(state, &[&k1, &k2], &[dt * DP45_A31, dt * DP45_A32]);
        let k3 = f(&y3, t + dt * DP45_C3);

        // Stage 4
        let y4 = vec_combine(
            state,
            &[&k1, &k2, &k3],
            &[dt * DP45_A41, dt * DP45_A42, dt * DP45_A43],
        );
        let k4 = f(&y4, t + dt * DP45_C4);

        // Stage 5
        let y5 = vec_combine(
            state,
            &[&k1, &k2, &k3, &k4],
            &[dt * DP45_A51, dt * DP45_A52, dt * DP45_A53, dt * DP45_A54],
        );
        let k5 = f(&y5, t + dt * DP45_C5);

        // Stage 6
        let y6 = vec_combine(
            state,
            &[&k1, &k2, &k3, &k4, &k5],
            &[
                dt * DP45_A61,
                dt * DP45_A62,
                dt * DP45_A63,
                dt * DP45_A64,
                dt * DP45_A65,
            ],
        );
        let k6 = f(&y6, t + dt);

        // 5th-order solution (stage 7 position)
        let y_5th = vec_combine(
            state,
            &[&k1, &k3, &k4, &k5, &k6],
            &[
                dt * DP45_B1,
                dt * DP45_B3,
                dt * DP45_B4,
                dt * DP45_B5,
                dt * DP45_B6,
            ],
        );

        // Stage 7 (FSAL — evaluate at the 5th-order solution)
        let k7 = f(&y_5th, t + dt);

        // 4th-order solution for error estimation
        let y_4th = vec_combine(
            state,
            &[&k1, &k3, &k4, &k5, &k6, &k7],
            &[
                dt * DP45_E1,
                dt * DP45_E3,
                dt * DP45_E4,
                dt * DP45_E5,
                dt * DP45_E6,
                dt * DP45_E7,
            ],
        );

        let error = max_error_norm(&y_5th, &y_4th);

        check_nan(&y_5th)?;

        // Store k7 for FSAL reuse
        self.fsal_k7 = Some(k7);

        // PI controller for step size (Gustafsson)
        let safety = 0.9;
        let p = 5.0; // order of the higher method
        let suggested_dt = if error > 0.0 {
            // Standard adaptive formula: dt_new = safety * dt * (tol / err)^(1/p)
            // We return a ratio-based suggestion; the AdaptiveIntegrator uses it.
            // Here we just embed a reasonable suggestion assuming tolerance ~= error
            // for maximum reuse. The adaptive wrapper overrides this.
            dt * safety * (1.0 / (error + 1e-30)).powf(1.0 / (p + 1.0)).min(5.0)
        } else {
            dt * 2.0 // error is zero, double the step
        };

        Ok(IntegratorOutput {
            new_state: y_5th,
            error_estimate: Some(error),
            suggested_dt: Some(suggested_dt),
        })
    }
}

// =========================================================================
// 2. Dormand-Prince RK8(7) — 13 stages, high precision
// =========================================================================

/// Dormand-Prince 8(7) embedded Runge-Kutta integrator.
///
/// A 13-stage method delivering 8th-order accuracy with a 7th-order embedded
/// solution for error estimation. This is the method of choice when very high
/// precision is required for long-time dynamics (e.g., computing Lyapunov
/// exponents or adiabatic invariants).
///
/// The coefficients are from:
/// Dormand & Prince, "A family of embedded Runge-Kutta formulae",
/// J. Comput. Appl. Math. 6, 19-26 (1980) — extended 8(7) pair.
pub struct DormandPrince87 {
    _private: (),
}

impl DormandPrince87 {
    /// Create a new Dormand-Prince 8(7) integrator.
    pub fn new() -> Self {
        Self { _private: () }
    }
}

impl Default for DormandPrince87 {
    fn default() -> Self {
        Self::new()
    }
}

/// Butcher tableau for the Dormand-Prince 8(7) method.
///
/// The nodes (c values) for the 13 stages.
const DP87_C: [f64; 13] = [
    0.0,
    1.0 / 18.0,
    1.0 / 12.0,
    1.0 / 8.0,
    5.0 / 16.0,
    3.0 / 8.0,
    59.0 / 400.0,
    93.0 / 200.0,
    5490023248.0 / 9719169821.0,
    13.0 / 20.0,
    1201146811.0 / 1299019798.0,
    1.0,
    1.0,
];

/// 8th-order weights.
const DP87_B8: [f64; 13] = [
    14005451.0 / 335480064.0,
    0.0,
    0.0,
    0.0,
    0.0,
    -59238493.0 / 1068277825.0,
    181606767.0 / 758867731.0,
    561292985.0 / 797845732.0,
    -1041891430.0 / 1371343529.0,
    760417239.0 / 1151165299.0,
    118820643.0 / 751138087.0,
    -528747749.0 / 2220607170.0,
    1.0 / 4.0,
];

/// 7th-order weights (for error estimation).
const DP87_B7: [f64; 13] = [
    13451932.0 / 455176623.0,
    0.0,
    0.0,
    0.0,
    0.0,
    -808719846.0 / 976000145.0,
    1757004468.0 / 5645159321.0,
    656045339.0 / 265891186.0,
    -3867574721.0 / 1518517206.0,
    465885868.0 / 322736535.0,
    53011238.0 / 667516719.0,
    2.0 / 45.0,
    0.0,
];

impl Integrator for DormandPrince87 {
    fn step(
        &mut self,
        state: &[Vector3<f64>],
        t: f64,
        dt: f64,
        f: &RhsFn<'_>,
    ) -> Result<IntegratorOutput> {
        // Full Butcher tableau a-matrix (only non-zero entries used per stage).
        // We store them inline to avoid large static arrays while keeping the
        // code readable.

        let k1 = f(state, t);

        // Stage 2
        let y2 = vec_add_scaled(state, &k1, dt / 18.0);
        let k2 = f(&y2, t + dt * DP87_C[1]);

        // Stage 3
        let y3 = vec_combine(state, &[&k1, &k2], &[dt / 48.0, dt / 16.0]);
        let k3 = f(&y3, t + dt * DP87_C[2]);

        // Stage 4
        let y4 = vec_combine(state, &[&k1, &k3], &[dt / 32.0, dt * 3.0 / 32.0]);
        let k4 = f(&y4, t + dt * DP87_C[3]);

        // Stage 5
        let y5 = vec_combine(
            state,
            &[&k1, &k3, &k4],
            &[dt * 5.0 / 16.0, dt * (-75.0 / 64.0), dt * 75.0 / 64.0],
        );
        let k5 = f(&y5, t + dt * DP87_C[4]);

        // Stage 6
        let y6 = vec_combine(
            state,
            &[&k1, &k4, &k5],
            &[dt * 3.0 / 80.0, dt * 3.0 / 16.0, dt * 3.0 / 20.0],
        );
        let k6 = f(&y6, t + dt * DP87_C[5]);

        // Stage 7
        let y7 = vec_combine(
            state,
            &[&k1, &k4, &k5, &k6],
            &[
                dt * 29443841.0 / 614563906.0,
                dt * 77736538.0 / 692538347.0,
                dt * (-28693883.0 / 1125000000.0),
                dt * 23124283.0 / 1800000000.0,
            ],
        );
        let k7 = f(&y7, t + dt * DP87_C[6]);

        // Stage 8
        let y8 = vec_combine(
            state,
            &[&k1, &k4, &k5, &k6, &k7],
            &[
                dt * 16016141.0 / 946692911.0,
                dt * 61564180.0 / 158732637.0,
                dt * 22789713.0 / 633445777.0,
                dt * 545815736.0 / 2771057229.0,
                dt * (-180193667.0 / 1043307555.0),
            ],
        );
        let k8 = f(&y8, t + dt * DP87_C[7]);

        // Stage 9
        let y9 = vec_combine(
            state,
            &[&k1, &k4, &k5, &k6, &k7, &k8],
            &[
                dt * 39632708.0 / 573591083.0,
                dt * (-433636366.0 / 683701615.0),
                dt * (-421739975.0 / 2616292301.0),
                dt * 100302831.0 / 723423059.0,
                dt * 790204164.0 / 839813087.0,
                dt * 800635310.0 / 3783071287.0,
            ],
        );
        let k9 = f(&y9, t + dt * DP87_C[8]);

        // Stage 10
        let y10 = vec_combine(
            state,
            &[&k1, &k4, &k5, &k6, &k7, &k8, &k9],
            &[
                dt * 246121993.0 / 1340847787.0,
                dt * (-37695042795.0 / 15268766246.0),
                dt * (-309121744.0 / 1061227803.0),
                dt * (-12992083.0 / 490766935.0),
                dt * 6005943493.0 / 2108947869.0,
                dt * 393006217.0 / 1396673457.0,
                dt * 123872331.0 / 1001029789.0,
            ],
        );
        let k10 = f(&y10, t + dt * DP87_C[9]);

        // Stage 11
        let y11 = vec_combine(
            state,
            &[&k1, &k4, &k5, &k6, &k7, &k8, &k9, &k10],
            &[
                dt * (-1028468189.0 / 846180014.0),
                dt * 8478235783.0 / 508512852.0,
                dt * 1311729495.0 / 1432422823.0,
                dt * (-10304129995.0 / 1701304382.0),
                dt * (-48777925059.0 / 3047939560.0),
                dt * 15336726248.0 / 1032824649.0,
                dt * (-45442868181.0 / 3398467696.0),
                dt * 3065993473.0 / 597172653.0,
            ],
        );
        let k11 = f(&y11, t + dt * DP87_C[10]);

        // Stage 12
        let y12 = vec_combine(
            state,
            &[&k1, &k4, &k5, &k6, &k7, &k8, &k9, &k10, &k11],
            &[
                dt * 185892177.0 / 718116043.0,
                dt * (-3185094517.0 / 667107341.0),
                dt * (-477755414.0 / 1098053517.0),
                dt * (-703635378.0 / 230739211.0),
                dt * 5731566787.0 / 1027545527.0,
                dt * 5232866602.0 / 850066563.0,
                dt * (-4093664535.0 / 808688257.0),
                dt * 3962137247.0 / 1805957418.0,
                dt * 65686358.0 / 487910083.0,
            ],
        );
        let k12 = f(&y12, t + dt * DP87_C[11]);

        // Stage 13
        let y13 = vec_combine(
            state,
            &[&k1, &k4, &k5, &k6, &k7, &k8, &k9, &k10, &k11, &k12],
            &[
                dt * 403863854.0 / 491063109.0,
                dt * (-5068492393.0 / 434740067.0),
                dt * (-411421997.0 / 543043805.0),
                dt * 652783627.0 / 914296604.0,
                dt * 11173962825.0 / 925320556.0,
                dt * (-13158990841.0 / 6184727034.0),
                dt * 3936647629.0 / 1978049680.0,
                dt * (-160528059.0 / 685178525.0),
                dt * 248638103.0 / 1413531060.0,
                dt * 0.0,
            ],
        );
        let k13 = f(&y13, t + dt * DP87_C[12]);

        // Collect all stages
        let ks: [&[Vector3<f64>]; 13] = [
            &k1, &k2, &k3, &k4, &k5, &k6, &k7, &k8, &k9, &k10, &k11, &k12, &k13,
        ];

        // 8th-order solution
        let n = state.len();
        let mut y8th = state.to_vec();
        for i in 0..n {
            let mut sum = Vector3::zero();
            for (stage, &bw) in ks.iter().zip(DP87_B8.iter()) {
                sum = sum + stage[i] * bw;
            }
            y8th[i] = y8th[i] + sum * dt;
        }

        // 7th-order solution
        let mut y7th = state.to_vec();
        for i in 0..n {
            let mut sum = Vector3::zero();
            for (stage, &bw) in ks.iter().zip(DP87_B7.iter()) {
                sum = sum + stage[i] * bw;
            }
            y7th[i] = y7th[i] + sum * dt;
        }

        let error = max_error_norm(&y8th, &y7th);

        check_nan(&y8th)?;

        let safety = 0.9;
        let p = 8.0;
        let suggested_dt = if error > 0.0 {
            dt * safety * (1.0 / (error + 1e-30)).powf(1.0 / (p + 1.0)).min(5.0)
        } else {
            dt * 2.0
        };

        Ok(IntegratorOutput {
            new_state: y8th,
            error_estimate: Some(error),
            suggested_dt: Some(suggested_dt),
        })
    }
}
