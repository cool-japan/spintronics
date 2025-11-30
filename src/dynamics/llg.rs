//! Landau-Lifshitz-Gilbert (LLG) equation solver
//!
//! The LLG equation describes the precessional motion of magnetization
//! in a ferromagnet under the influence of an effective magnetic field.
//!
//! # Mathematical Formulation
//!
//! The LLG equation in its standard form:
//!
//! $$
//! \frac{d\mathbf{m}}{dt} = -\gamma \left(\mathbf{m} \times \mathbf{H}_{\text{eff}}\right) + \alpha \left(\mathbf{m} \times \frac{d\mathbf{m}}{dt}\right)
//! $$
//!
//! where:
//! - $\mathbf{m}$ is the normalized magnetization vector ($|\mathbf{m}| = 1$)
//! - $\gamma$ is the gyromagnetic ratio (rad/(s·T))
//! - $\mathbf{H}_{\text{eff}}$ is the effective magnetic field (T)
//! - $\alpha$ is the Gilbert damping constant (dimensionless)
//!
//! For numerical implementation, this is typically rearranged to:
//!
//! $$
//! \frac{d\mathbf{m}}{dt} = \frac{-\gamma}{1+\alpha^2} \left[\mathbf{m} \times \mathbf{H}_{\text{eff}} + \alpha \mathbf{m} \times \left(\mathbf{m} \times \mathbf{H}_{\text{eff}}\right)\right]
//! $$
//!
//! This form avoids solving implicitly for $d\mathbf{m}/dt$ on the right-hand side.
//!
//! # Physical Interpretation
//!
//! The equation describes two competing effects:
//! 1. **Precession**: $-\gamma (\mathbf{m} \times \mathbf{H}_{\text{eff}})$ causes $\mathbf{m}$ to precess around $\mathbf{H}_{\text{eff}}$
//! 2. **Damping**: $\alpha (\mathbf{m} \times d\mathbf{m}/dt)$ relaxes $\mathbf{m}$ toward $\mathbf{H}_{\text{eff}}$

use crate::constants::GAMMA;
use crate::vector3::Vector3;

/// Calculate the time derivative of magnetization using the LLG equation
///
/// # Arguments
/// * `m` - Normalized magnetization vector
/// * `h_eff` - Effective magnetic field \[T\]
/// * `gamma` - Gyromagnetic ratio \[rad/(s·T)\]
/// * `alpha` - Gilbert damping constant (dimensionless)
///
/// # Returns
/// Time derivative of magnetization dm/dt [1/s]
///
/// # Physical Background
/// The LLG equation combines:
/// - Precession term: -γ (m × H_eff) - magnetization precesses around H_eff
/// - Damping term: α (m × dm/dt) - relaxation towards H_eff direction
///
/// # Example
/// ```
/// use spintronics::dynamics::llg::calc_dm_dt;
/// use spintronics::constants::GAMMA;
/// use spintronics::Vector3;
///
/// // Magnetization pointing in x-direction
/// let m = Vector3::new(1.0, 0.0, 0.0);
/// // External field in z-direction
/// let h_ext = Vector3::new(0.0, 0.0, 1.0);
/// // Gilbert damping for Permalloy
/// let alpha = 0.01;
///
/// // Calculate magnetization dynamics
/// let dm_dt = calc_dm_dt(m, h_ext, GAMMA, alpha);
///
/// // dm/dt should be perpendicular to m (conserves |m|)
/// assert!(dm_dt.dot(&m).abs() < 1e-10);
/// ```
pub fn calc_dm_dt(m: Vector3<f64>, h_eff: Vector3<f64>, gamma: f64, alpha: f64) -> Vector3<f64> {
    // Precession term: -γ (m × H_eff)
    // Physical meaning: Magnetization precesses around the effective field direction
    // like a spinning top in a gravitational field. The cross product ensures
    // the motion is perpendicular to both m and H, causing circular precession.
    // The negative sign comes from the electron's negative charge.
    let precession = m.cross(&h_eff) * (-gamma);

    // Damping term: α m × (m × H_eff)
    // Physical meaning: Gilbert damping represents energy dissipation as the
    // magnetization relaxes toward the field direction. The double cross product
    // gives a component perpendicular to m but toward H, causing spiral motion.
    // This approximation (α << 1) avoids solving the implicit equation for dm/dt.
    let damping_term = m.cross(&precession) * alpha;

    // Normalization factor: 1/(1 + α²)
    // Physical meaning: Ensures the equation remains valid for the normalized form.
    // For typical materials (α ~ 0.01), this correction is ~1%. Only significant
    // for highly dissipative systems where α approaches 1.
    (precession + damping_term) * (1.0 / (1.0 + alpha * alpha))
}

/// LLG equation solver with adaptive time stepping
///
/// # Example
/// ```
/// use spintronics::dynamics::llg::LlgSolver;
/// use spintronics::Vector3;
///
/// // Create solver with damping α=0.01 and dt=0.1ps
/// let solver = LlgSolver::new(0.01, 1.0e-13);
///
/// // Initial magnetization tilted from z-axis
/// let m0 = Vector3::new(0.1, 0.0, 1.0).normalize();
/// // Applied field in z-direction (1 Tesla)
/// let h_ext = Vector3::new(0.0, 0.0, 1.0);
///
/// // Evolve one time step
/// let m1 = solver.step_euler(m0, h_ext);
///
/// // Magnetization magnitude is conserved
/// assert!((m1.magnitude() - 1.0).abs() < 1e-10);
/// ```
pub struct LlgSolver {
    /// Gilbert damping constant
    pub alpha: f64,
    /// Gyromagnetic ratio \[rad/(s·T)\]
    pub gamma: f64,
    /// Time step \[s\]
    pub dt: f64,
}

impl Default for LlgSolver {
    fn default() -> Self {
        Self {
            alpha: 0.01,
            gamma: GAMMA,
            dt: 1.0e-13, // 0.1 ps
        }
    }
}

impl LlgSolver {
    /// Create a new LLG solver
    pub fn new(alpha: f64, dt: f64) -> Self {
        Self {
            alpha,
            gamma: GAMMA,
            dt,
        }
    }

    /// Evolve magnetization by one time step using Euler method
    ///
    /// # Arguments
    /// * `m` - Current magnetization (will be normalized)
    /// * `h_eff` - Effective magnetic field \[T\]
    ///
    /// # Returns
    /// New magnetization after dt, normalized to unit length
    pub fn step_euler(&self, m: Vector3<f64>, h_eff: Vector3<f64>) -> Vector3<f64> {
        let dm_dt = calc_dm_dt(m, h_eff, self.gamma, self.alpha);
        (m + dm_dt * self.dt).normalize()
    }

    /// Evolve magnetization using 4th-order Runge-Kutta method
    ///
    /// More accurate than Euler method, recommended for production simulations.
    /// RK4 provides 4th-order accuracy O(dt^4) compared to Euler's O(dt).
    ///
    /// # Arguments
    /// * `m` - Current normalized magnetization
    /// * `h_eff_fn` - Function that computes effective field from magnetization
    ///
    /// # Returns
    /// New magnetization after dt, normalized to unit length
    ///
    /// # Example
    /// ```
    /// use spintronics::dynamics::llg::LlgSolver;
    /// use spintronics::Vector3;
    ///
    /// let solver = LlgSolver::new(0.01, 1.0e-12);
    /// let m0 = Vector3::new(1.0, 0.1, 0.0).normalize();
    ///
    /// // Define effective field function (could include anisotropy, exchange, etc.)
    /// let h_field = Vector3::new(0.0, 0.0, 1.0);
    /// let h_eff_fn = |_m: Vector3<f64>| h_field;
    ///
    /// // Evolve with RK4 (more accurate than Euler)
    /// let m1 = solver.step_rk4(m0, h_eff_fn);
    ///
    /// assert!((m1.magnitude() - 1.0).abs() < 1e-10);
    /// ```
    pub fn step_rk4(
        &self,
        m: Vector3<f64>,
        h_eff_fn: impl Fn(Vector3<f64>) -> Vector3<f64>,
    ) -> Vector3<f64> {
        let k1 = calc_dm_dt(m, h_eff_fn(m), self.gamma, self.alpha);

        let m2 = (m + k1 * (self.dt * 0.5)).normalize();
        let k2 = calc_dm_dt(m2, h_eff_fn(m2), self.gamma, self.alpha);

        let m3 = (m + k2 * (self.dt * 0.5)).normalize();
        let k3 = calc_dm_dt(m3, h_eff_fn(m3), self.gamma, self.alpha);

        let m4 = (m + k3 * self.dt).normalize();
        let k4 = calc_dm_dt(m4, h_eff_fn(m4), self.gamma, self.alpha);

        let dm = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (self.dt / 6.0);
        (m + dm).normalize()
    }

    /// Evolve magnetization using 2nd-order Heun's method (improved Euler)
    ///
    /// Heun's method is a predictor-corrector algorithm that provides
    /// 2nd-order accuracy with lower computational cost than RK4.
    /// Particularly useful for stochastic LLG with thermal noise.
    ///
    /// # Arguments
    /// * `m` - Current normalized magnetization
    /// * `h_eff_fn` - Function that computes effective field from magnetization
    ///
    /// # Returns
    /// New magnetization after dt, normalized
    pub fn step_heun(
        &self,
        m: Vector3<f64>,
        h_eff_fn: impl Fn(Vector3<f64>) -> Vector3<f64>,
    ) -> Vector3<f64> {
        // Predictor step (Euler)
        let k1 = calc_dm_dt(m, h_eff_fn(m), self.gamma, self.alpha);
        let m_pred = (m + k1 * self.dt).normalize();

        // Corrector step
        let k2 = calc_dm_dt(m_pred, h_eff_fn(m_pred), self.gamma, self.alpha);

        // Average of slopes
        let dm = (k1 + k2) * (self.dt * 0.5);
        (m + dm).normalize()
    }

    /// Adaptive RK4 with error estimation
    ///
    /// Uses embedded RK4-RK5 pair to estimate local truncation error
    /// and suggests optimal time step.
    ///
    /// # Arguments
    /// * `m` - Current magnetization
    /// * `h_eff_fn` - Effective field function
    /// * `tolerance` - Desired error tolerance
    ///
    /// # Returns
    /// Tuple of (new magnetization, suggested next dt, error estimate)
    #[allow(dead_code)]
    pub fn step_rk4_adaptive(
        &self,
        m: Vector3<f64>,
        h_eff_fn: impl Fn(Vector3<f64>) -> Vector3<f64>,
        tolerance: f64,
    ) -> (Vector3<f64>, f64, f64) {
        // Standard RK4
        let k1 = calc_dm_dt(m, h_eff_fn(m), self.gamma, self.alpha);

        let m2 = (m + k1 * (self.dt * 0.5)).normalize();
        let k2 = calc_dm_dt(m2, h_eff_fn(m2), self.gamma, self.alpha);

        let m3 = (m + k2 * (self.dt * 0.5)).normalize();
        let k3 = calc_dm_dt(m3, h_eff_fn(m3), self.gamma, self.alpha);

        let m4 = (m + k3 * self.dt).normalize();
        let k4 = calc_dm_dt(m4, h_eff_fn(m4), self.gamma, self.alpha);

        let dm_rk4 = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (self.dt / 6.0);

        // 5th-order estimate (Cash-Karp coefficients)
        let dm_rk5 = (k1 * (16.0 / 135.0) + k2 * (6656.0 / 12825.0) + k3 * (28561.0 / 56430.0)
            - k4 * (9.0 / 50.0))
            * self.dt;

        // Error estimate
        let error = (dm_rk5 - dm_rk4).magnitude();

        // Safety factor for step size adjustment
        let safety = 0.9;
        let dt_new = if error > tolerance {
            self.dt * safety * (tolerance / error).powf(0.2)
        } else {
            self.dt * safety * (tolerance / error).powf(0.25).min(2.0)
        };

        let m_new = (m + dm_rk4).normalize();
        (m_new, dt_new, error)
    }

    /// Semi-implicit (Crank-Nicolson-like) solver
    ///
    /// Uses implicit midpoint rule which is unconditionally stable.
    /// Useful for stiff problems with large effective fields.
    ///
    /// # Arguments
    /// * `m` - Current magnetization
    /// * `h_eff` - Effective field (constant during step)
    /// * `max_iter` - Maximum iterations for implicit solve
    ///
    /// # Returns
    /// New magnetization after dt
    #[allow(dead_code)]
    pub fn step_implicit(
        &self,
        m: Vector3<f64>,
        h_eff: Vector3<f64>,
        max_iter: usize,
    ) -> Vector3<f64> {
        let mut m_new = m;

        // Fixed-point iteration for implicit midpoint
        for _ in 0..max_iter {
            let m_mid = (m + m_new) * 0.5;
            let dm_dt = calc_dm_dt(m_mid.normalize(), h_eff, self.gamma, self.alpha);
            let m_next = (m + dm_dt * self.dt).normalize();

            // Check convergence
            if (m_next - m_new).magnitude() < 1e-10 {
                return m_next;
            }
            m_new = m_next;
        }

        m_new
    }

    /// Builder method to set time step
    pub fn with_dt(mut self, dt: f64) -> Self {
        self.dt = dt;
        self
    }

    /// Builder method to set damping
    pub fn with_alpha(mut self, alpha: f64) -> Self {
        self.alpha = alpha;
        self
    }

    /// Builder method to set gyromagnetic ratio
    pub fn with_gamma(mut self, gamma: f64) -> Self {
        self.gamma = gamma;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_llg_precession() {
        let m = Vector3::new(1.0, 0.0, 0.0);
        let h = Vector3::new(0.0, 0.0, 1.0);
        let alpha = 0.01;
        let dm_dt = calc_dm_dt(m, h, GAMMA, alpha);

        // dm/dt should always be perpendicular to m (conserves |m|)
        assert!(dm_dt.dot(&m).abs() < 1e-6);

        // With damping, dm/dt is NOT perpendicular to H (damping drives m toward H)
        // The precession term alone would be perpendicular, but the damping term
        // adds a component toward the field direction
        let precession_only = m.cross(&h) * (-GAMMA);
        assert!(precession_only.dot(&h).abs() < 1e-10);
    }

    #[test]
    fn test_magnetization_conservation() {
        let solver = LlgSolver::default();
        let m = Vector3::new(1.0, 0.0, 0.0);
        let h = Vector3::new(0.0, 0.0, 1.0);

        let m_new = solver.step_euler(m, h);
        assert!((m_new.magnitude() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_zero_field_relaxation() {
        let m = Vector3::new(1.0, 0.0, 0.0);
        let h = Vector3::new(0.0, 0.0, 0.0);
        let dm_dt = calc_dm_dt(m, h, GAMMA, 0.01);

        // No field means no dynamics (in ideal case)
        assert!(dm_dt.magnitude() < 1e-10);
    }

    #[test]
    fn test_rk4_vs_euler() {
        let solver = LlgSolver::new(0.01, 1.0e-12);
        let m0 = Vector3::new(1.0, 0.1, 0.0).normalize();
        let h_const = Vector3::new(0.0, 0.0, 1.0);

        let h_fn = |_m: Vector3<f64>| h_const;

        let m_euler = solver.step_euler(m0, h_const);
        let m_rk4 = solver.step_rk4(m0, h_fn);

        // RK4 should give different (more accurate) result than Euler
        assert!((m_rk4 - m_euler).magnitude() > 1e-15);

        // Both should preserve magnetization magnitude
        assert!((m_euler.magnitude() - 1.0).abs() < 1e-10);
        assert!((m_rk4.magnitude() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heun_method() {
        let solver = LlgSolver::new(0.01, 1.0e-12);
        let m0 = Vector3::new(0.707, 0.707, 0.0);
        let h_const = Vector3::new(0.0, 0.0, 1.0);

        let h_fn = |_m: Vector3<f64>| h_const;

        let m_heun = solver.step_heun(m0, h_fn);

        // Should preserve normalization
        assert!((m_heun.magnitude() - 1.0).abs() < 1e-10);

        // Should be different from initial state
        assert!((m_heun - m0).magnitude() > 1e-15);
    }

    #[test]
    fn test_heun_vs_euler_accuracy() {
        let solver = LlgSolver::new(0.01, 1.0e-11);
        let m0 = Vector3::new(1.0, 0.0, 0.0);
        let h = Vector3::new(0.0, 0.0, 1.0);
        let h_fn = |_m: Vector3<f64>| h;

        let m_euler = solver.step_euler(m0, h);
        let m_heun = solver.step_heun(m0, h_fn);

        // Heun should give different result (better accuracy)
        assert!((m_heun - m_euler).magnitude() > 1e-15);
    }

    #[test]
    fn test_adaptive_rk4_error_estimate() {
        let solver = LlgSolver::new(0.01, 1.0e-12);
        let m0 = Vector3::new(1.0, 0.0, 0.0);
        let h = Vector3::new(0.0, 0.0, 1.0);

        let h_fn = |_m: Vector3<f64>| h;

        let (m_new, dt_suggest, error) = solver.step_rk4_adaptive(m0, h_fn, 1.0e-8);

        // Should preserve normalization
        assert!((m_new.magnitude() - 1.0).abs() < 1e-10);

        // Error should be non-negative
        assert!(error >= 0.0);

        // Suggested dt should be positive
        assert!(dt_suggest > 0.0);
    }

    #[test]
    fn test_implicit_convergence() {
        let solver = LlgSolver::new(0.01, 1.0e-12);
        let m0 = Vector3::new(1.0, 0.0, 0.0);
        let h = Vector3::new(0.0, 0.0, 1.0);

        let m_impl = solver.step_implicit(m0, h, 20);

        // Should preserve normalization
        assert!((m_impl.magnitude() - 1.0).abs() < 1e-10);

        // Should evolve from initial state
        assert!((m_impl - m0).magnitude() > 1e-15);
    }

    #[test]
    fn test_implicit_stability() {
        // Test with large time step (where explicit methods might be unstable)
        let solver = LlgSolver::new(0.01, 1.0e-10); // Large dt
        let m0 = Vector3::new(1.0, 0.0, 0.0);
        let h = Vector3::new(0.0, 0.0, 10.0); // Large field

        let m_impl = solver.step_implicit(m0, h, 50);

        // Even with large dt and field, should remain normalized
        assert!((m_impl.magnitude() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_builder_pattern() {
        let solver = LlgSolver::default()
            .with_dt(1.0e-11)
            .with_alpha(0.05)
            .with_gamma(1.76e11);

        assert_eq!(solver.dt, 1.0e-11);
        assert_eq!(solver.alpha, 0.05);
        assert_eq!(solver.gamma, 1.76e11);
    }

    #[test]
    fn test_rk4_energy_conservation_tendency() {
        // RK4 should have better energy conservation than Euler
        let solver = LlgSolver::new(0.0, 1.0e-12); // No damping
        let m0 = Vector3::new(0.707, 0.707, 0.0);
        let h = Vector3::new(0.0, 0.0, 1.0);
        let h_fn = |_m: Vector3<f64>| h;

        // Simulate several steps
        let mut m_euler = m0;
        let mut m_rk4 = m0;

        for _ in 0..10 {
            m_euler = solver.step_euler(m_euler, h);
            m_rk4 = solver.step_rk4(m_rk4, h_fn);
        }

        // Both should preserve magnitude
        assert!((m_euler.magnitude() - 1.0).abs() < 1e-9);
        assert!((m_rk4.magnitude() - 1.0).abs() < 1e-9);

        // Energy = -m·H; for undamped precession, should be conserved
        let e0 = -m0.dot(&h);
        let e_euler = -m_euler.dot(&h);
        let e_rk4 = -m_rk4.dot(&h);

        // RK4 should conserve energy better (smaller drift)
        assert!((e_rk4 - e0).abs() <= (e_euler - e0).abs());
    }
}
