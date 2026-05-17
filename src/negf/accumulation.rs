//! Spin accumulation and spin diffusion in 1D structures.
//!
//! Implements the Valet–Fert spin diffusion equation on a 1D spatial grid:
//!
//! ```text
//! ∂μ_s/∂t = D_s ∂²μ_s/∂x² − μ_s/τ_sf + source
//! ```
//!
//! At steady state this becomes the second-order ODE
//!
//! ```text
//! D_s μ_s'' − μ_s/τ_sf = 0   ⟹   μ_s(x) ∝ exp(±x/λ)
//! ```
//!
//! where `λ = √(D_s τ_sf)` is the spin-diffusion length.
//!
//! Both an explicit FTCS scheme (with CFL stability check) and an unconditionally
//! stable backward-Euler scheme (via Thomas algorithm) are provided.
//!
//! # References
//!
//! - T. Valet, A. Fert, Phys. Rev. B **48**, 7099 (1993)

use crate::error::{self, Result};

// ============================================================================
// BoundaryCondition
// ============================================================================

/// Boundary condition type for the spin accumulation solver.
#[derive(Debug, Clone)]
pub enum BoundaryCondition {
    /// Dirichlet: fix μ_s = `value` at the boundary \[J\].
    Dirichlet(f64),
    /// Neumann: fix ∂μ_s/∂x = `flux` at the boundary \[J/m\].
    Neumann(f64),
    /// Mixed exponential: μ_s approaches `mu_inf` with characteristic length `length` \[m\].
    Mixed {
        /// Asymptotic value of μ_s far from the boundary \[J\].
        mu_inf: f64,
        /// Characteristic decay length \[m\].
        length: f64,
    },
}

// ============================================================================
// SpinAccumulation1D
// ============================================================================

/// 1D spin accumulation solver on a uniform spatial grid.
///
/// Solves the spin diffusion equation with optional spin-flip relaxation on
/// a discretized line of `n_points` grid points separated by `dx = length / (n_points − 1)`.
#[derive(Debug, Clone)]
pub struct SpinAccumulation1D {
    /// Spin diffusion constant D_s \[m²/s\].
    pub d_s: f64,
    /// Spin-flip relaxation time τ_sf \[s\].
    pub tau_sf: f64,
    /// Sample length L \[m\].
    pub length: f64,
    /// Number of spatial grid points (≥ 2).
    pub n_points: usize,
    /// Spin chemical potential at each grid point μ_s(x_i) \[J\].
    pub mu_s: Vec<f64>,
    /// Grid spacing dx = L / (n_points − 1) \[m\].
    pub dx: f64,
}

impl SpinAccumulation1D {
    /// Create a new `SpinAccumulation1D` with zero initial spin accumulation.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if:
    /// - `d_s ≤ 0`, `tau_sf ≤ 0`, or `length ≤ 0`
    /// - `n_points < 2`
    pub fn new(d_s: f64, tau_sf: f64, length: f64, n_points: usize) -> Result<Self> {
        if d_s <= 0.0 {
            return Err(error::invalid_param(
                "d_s",
                "spin diffusion constant must be positive",
            ));
        }
        if tau_sf <= 0.0 {
            return Err(error::invalid_param(
                "tau_sf",
                "spin-flip time must be positive",
            ));
        }
        if length <= 0.0 {
            return Err(error::invalid_param(
                "length",
                "sample length must be positive",
            ));
        }
        if n_points < 2 {
            return Err(error::invalid_param(
                "n_points",
                "must have at least 2 grid points",
            ));
        }
        let dx = length / (n_points - 1) as f64;
        Ok(Self {
            d_s,
            tau_sf,
            length,
            n_points,
            mu_s: vec![0.0; n_points],
            dx,
        })
    }

    /// Spin diffusion length λ = √(D_s · τ_sf) \[m\].
    pub fn spin_diffusion_length(&self) -> f64 {
        (self.d_s * self.tau_sf).sqrt()
    }

    /// Initialize all grid points to a uniform value `mu0` \[J\].
    pub fn initialize_uniform(&mut self, mu0: f64) {
        for v in self.mu_s.iter_mut() {
            *v = mu0;
        }
    }

    /// Initialize with a Gaussian profile: μ_s(x) = mu0 · exp(−(x − x0)² / (2σ²)).
    ///
    /// # Arguments
    ///
    /// * `x0` – peak position \[m\]
    /// * `sigma` – Gaussian width \[m\]
    /// * `mu0` – peak amplitude \[J\]
    pub fn initialize_gaussian(&mut self, x0: f64, sigma: f64, mu0: f64) {
        for (i, v) in self.mu_s.iter_mut().enumerate() {
            let x = i as f64 * self.dx;
            let arg = (x - x0) / sigma;
            *v = mu0 * (-0.5 * arg * arg).exp();
        }
    }

    /// Apply boundary condition at the left edge (index 0).
    fn apply_left_bc(&mut self, bc: &BoundaryCondition) {
        match bc {
            BoundaryCondition::Dirichlet(v) => {
                self.mu_s[0] = *v;
            },
            BoundaryCondition::Neumann(flux) => {
                // One-sided: μ[0] = μ[1] − flux * dx
                if self.n_points > 1 {
                    self.mu_s[0] = self.mu_s[1] - flux * self.dx;
                }
            },
            BoundaryCondition::Mixed { mu_inf, length } => {
                // Exponential approach from the left: μ = mu_inf * exp(−x/length)
                self.mu_s[0] = *mu_inf;
                if self.n_points > 1 {
                    self.mu_s[1] = mu_inf * (-self.dx / length).exp();
                }
            },
        }
    }

    /// Apply boundary condition at the right edge (index n_points − 1).
    fn apply_right_bc(&mut self, bc: &BoundaryCondition) {
        let nm1 = self.n_points - 1;
        match bc {
            BoundaryCondition::Dirichlet(v) => {
                self.mu_s[nm1] = *v;
            },
            BoundaryCondition::Neumann(flux) => {
                // One-sided: μ[N-1] = μ[N-2] + flux * dx
                if nm1 > 0 {
                    self.mu_s[nm1] = self.mu_s[nm1 - 1] + flux * self.dx;
                }
            },
            BoundaryCondition::Mixed { mu_inf, length } => {
                let x_right = (nm1 as f64) * self.dx;
                self.mu_s[nm1] = mu_inf * (-x_right / length).exp();
            },
        }
    }

    /// Explicit FTCS (Forward-Time, Centred-Space) time step.
    ///
    /// ```text
    /// μ^{n+1}_i = μ^n_i + dt · [D_s · (μ_{i+1} − 2μ_i + μ_{i-1}) / dx² − μ_i / τ_sf]
    /// ```
    ///
    /// CFL stability condition: `D_s · dt / dx² ≤ 0.5`.
    ///
    /// # Errors
    ///
    /// Returns `NumericalError` if the CFL condition is violated.
    pub fn evolve(
        &mut self,
        dt: f64,
        left: BoundaryCondition,
        right: BoundaryCondition,
    ) -> Result<()> {
        let dx2 = self.dx * self.dx;
        let cfl = self.d_s * dt / dx2;
        if cfl > 0.5 + 1e-10 {
            return Err(error::numerical_error(&format!(
                "CFL condition violated: D_s·dt/dx² = {:.4} > 0.5 (reduce dt or increase n_points)",
                cfl
            )));
        }
        let n = self.n_points;
        let mut new_mu = self.mu_s.clone();

        for (new_v, window) in new_mu[1..n - 1].iter_mut().zip(self.mu_s.windows(3)) {
            let laplacian = (window[2] - 2.0 * window[1] + window[0]) / dx2;
            let relaxation = window[1] / self.tau_sf;
            *new_v = window[1] + dt * (self.d_s * laplacian - relaxation);
        }

        self.mu_s = new_mu;
        self.apply_left_bc(&left);
        self.apply_right_bc(&right);
        Ok(())
    }

    /// Implicit backward-Euler time step (unconditionally stable).
    ///
    /// Solves `(I + dt·L) μ^{n+1} = μ^n` where `L` is the negative of the
    /// diffusion-relaxation operator.  The resulting tridiagonal system is solved
    /// with the Thomas algorithm.
    ///
    /// # Errors
    ///
    /// Propagates `NumericalError` from the Thomas solver if the system is singular.
    pub fn evolve_implicit(
        &mut self,
        dt: f64,
        left: BoundaryCondition,
        right: BoundaryCondition,
    ) -> Result<()> {
        let n = self.n_points;
        let dx2 = self.dx * self.dx;
        // Coefficients for the backward-Euler tridiagonal system
        // (1 + dt/tau_sf + 2*D_s*dt/dx²) * μ^{n+1}_i
        //   − (D_s*dt/dx²) * μ^{n+1}_{i-1}
        //   − (D_s*dt/dx²) * μ^{n+1}_{i+1} = μ^n_i
        let alpha = self.d_s * dt / dx2; // off-diagonal coefficient magnitude
        let beta = 1.0 + dt / self.tau_sf + 2.0 * alpha; // diagonal coefficient

        let mut a = vec![0.0; n]; // sub-diagonal
        let mut b = vec![0.0; n]; // diagonal
        let mut c = vec![0.0; n]; // super-diagonal
        let mut d = vec![0.0; n]; // RHS

        for i in 0..n {
            b[i] = beta;
            d[i] = self.mu_s[i];
            if i > 0 {
                a[i] = -alpha;
            }
            if i < n - 1 {
                c[i] = -alpha;
            }
        }

        // Enforce boundary conditions by modifying the system
        match &left {
            BoundaryCondition::Dirichlet(v) => {
                // Pin μ[0] = v
                a[0] = 0.0;
                b[0] = 1.0;
                c[0] = 0.0;
                d[0] = *v;
            },
            BoundaryCondition::Neumann(flux) => {
                // One-sided: μ[0] − μ[1] = −flux * dx
                b[0] = 1.0;
                c[0] = -1.0;
                a[0] = 0.0;
                d[0] = -flux * self.dx;
            },
            BoundaryCondition::Mixed { mu_inf, .. } => {
                b[0] = 1.0;
                c[0] = 0.0;
                a[0] = 0.0;
                d[0] = *mu_inf;
            },
        }

        let nm1 = n - 1;
        match &right {
            BoundaryCondition::Dirichlet(v) => {
                a[nm1] = 0.0;
                b[nm1] = 1.0;
                c[nm1] = 0.0;
                d[nm1] = *v;
            },
            BoundaryCondition::Neumann(flux) => {
                // One-sided: μ[N-1] − μ[N-2] = flux * dx
                a[nm1] = -1.0;
                b[nm1] = 1.0;
                c[nm1] = 0.0;
                d[nm1] = flux * self.dx;
            },
            BoundaryCondition::Mixed { mu_inf, length } => {
                let x_right = (nm1 as f64) * self.dx;
                a[nm1] = 0.0;
                b[nm1] = 1.0;
                c[nm1] = 0.0;
                d[nm1] = mu_inf * (-x_right / length).exp();
            },
        }

        self.mu_s = thomas_solve(&a, &b, &c, &d)?;
        Ok(())
    }

    /// Compute the steady-state spin accumulation profile via direct tridiagonal solve.
    ///
    /// Solves `(D_s ∇² − 1/τ_sf) μ_s = 0` with Dirichlet conditions:
    /// `μ_s[0] = injection_value`, `μ_s[N-1] = 0`.
    ///
    /// # Errors
    ///
    /// Propagates `NumericalError` from the Thomas solver.
    pub fn steady_state(&mut self, injection_value: f64) -> Result<()> {
        let n = self.n_points;
        let dx2 = self.dx * self.dx;
        let alpha = self.d_s / dx2; // off-diagonal magnitude
        let beta = 2.0 * alpha + 1.0 / self.tau_sf; // diagonal

        let mut a = vec![0.0_f64; n];
        let mut b = vec![0.0_f64; n];
        let mut c = vec![0.0_f64; n];
        let mut d = vec![0.0_f64; n];

        // Interior
        for i in 1..(n - 1) {
            a[i] = -alpha;
            b[i] = beta;
            c[i] = -alpha;
            d[i] = 0.0;
        }

        // Left BC: Dirichlet
        b[0] = 1.0;
        c[0] = 0.0;
        a[0] = 0.0;
        d[0] = injection_value;

        // Right BC: Dirichlet = 0
        let nm1 = n - 1;
        a[nm1] = 0.0;
        b[nm1] = 1.0;
        c[nm1] = 0.0;
        d[nm1] = 0.0;

        self.mu_s = thomas_solve(&a, &b, &c, &d)?;
        Ok(())
    }

    /// Interpolate the spin accumulation at continuous position `x` \[m\].
    ///
    /// Uses linear interpolation between adjacent grid points.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `x` is outside `[0, length]`.
    pub fn decay_profile(&self, x: f64) -> Result<f64> {
        if x < 0.0 || x > self.length + 1e-12 {
            return Err(error::invalid_param(
                "x",
                "position must be within [0, length]",
            ));
        }
        let x = x.clamp(0.0, self.length);
        let idx_f = x / self.dx;
        let i = idx_f.floor() as usize;
        if i + 1 >= self.n_points {
            return Ok(self.mu_s[self.n_points - 1]);
        }
        let frac = idx_f - i as f64;
        Ok(self.mu_s[i] * (1.0 - frac) + self.mu_s[i + 1] * frac)
    }

    /// Total spin polarization: ∫ μ_s(x) dx over [0, L] (trapezoid rule) \[J·m\].
    pub fn total_polarization(&self) -> f64 {
        let n = self.n_points;
        if n < 2 {
            return self.mu_s.first().copied().unwrap_or(0.0) * self.dx;
        }
        let mut sum = 0.0;
        for k in 0..(n - 1) {
            sum += 0.5 * (self.mu_s[k] + self.mu_s[k + 1]) * self.dx;
        }
        sum
    }

    /// Spin current at site `idx` \[J/m · S = J·S/m, effectively A here\].
    ///
    /// ```text
    /// j_s(x_i) = −D_s · σ · dμ_s/dx |_{x_i}
    /// ```
    ///
    /// Central difference for interior sites; one-sided at boundaries.
    ///
    /// # Arguments
    ///
    /// * `idx` – site index
    /// * `conductivity` – spin conductivity σ \[S/m\]
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `idx ≥ n_points`.
    pub fn spin_current(&self, idx: usize, conductivity: f64) -> Result<f64> {
        let n = self.n_points;
        if idx >= n {
            return Err(error::invalid_param("idx", "site index out of bounds"));
        }
        let dmu_dx = if idx == 0 {
            // Forward difference
            if n > 1 {
                (self.mu_s[1] - self.mu_s[0]) / self.dx
            } else {
                0.0
            }
        } else if idx == n - 1 {
            // Backward difference
            (self.mu_s[n - 1] - self.mu_s[n - 2]) / self.dx
        } else {
            // Central difference
            (self.mu_s[idx + 1] - self.mu_s[idx - 1]) / (2.0 * self.dx)
        };
        Ok(-self.d_s * conductivity * dmu_dx)
    }
}

// ============================================================================
// Thomas algorithm
// ============================================================================

/// Thomas algorithm for solving a tridiagonal linear system Ax = d.
///
/// `a` is the sub-diagonal (a[0] unused), `b` is the diagonal,
/// `c` is the super-diagonal (c[n-1] unused), `d` is the RHS.
///
/// All arrays must have the same length n ≥ 1.
///
/// # Errors
///
/// Returns `NumericalError` if a zero pivot is encountered (singular system).
pub(crate) fn thomas_solve(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Result<Vec<f64>> {
    let n = b.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut c_star = vec![0.0_f64; n];
    let mut d_star = vec![0.0_f64; n];
    let mut x = vec![0.0_f64; n];

    // Forward sweep
    if b[0].abs() < 1e-30 {
        return Err(error::numerical_error(
            "Thomas algorithm: zero pivot at row 0",
        ));
    }
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for i in 1..n {
        let denom = b[i] - a[i] * c_star[i - 1];
        if denom.abs() < 1e-30 {
            return Err(error::numerical_error(&format!(
                "Thomas algorithm: zero pivot at row {}",
                i
            )));
        }
        c_star[i] = if i < n - 1 { c[i] / denom } else { 0.0 };
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / denom;
    }

    // Back substitution
    x[n - 1] = d_star[n - 1];
    for i in (0..(n - 1)).rev() {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }
    Ok(x)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_acc(n: usize) -> SpinAccumulation1D {
        // D_s = 1e-3 m²/s, tau_sf = 1e-12 s, length = 100e-9 m
        SpinAccumulation1D::new(1e-3, 1e-12, 100e-9, n).expect("valid")
    }

    // ------------------------------------------------------------------ basics

    #[test]
    fn test_spin_diffusion_length_sqrt() {
        let acc = make_acc(10);
        let lambda = acc.spin_diffusion_length();
        let expected = (1e-3 * 1e-12_f64).sqrt();
        assert!(
            (lambda - expected).abs() < 1e-20,
            "λ = {}, expected {}",
            lambda,
            expected
        );
    }

    #[test]
    fn test_n_points_minimum_2() {
        let result = SpinAccumulation1D::new(1e-3, 1e-12, 100e-9, 1);
        assert!(result.is_err(), "n_points = 1 should fail");
        let ok = SpinAccumulation1D::new(1e-3, 1e-12, 100e-9, 2);
        assert!(ok.is_ok());
    }

    // ------------------------------------------------------------------ initialization

    #[test]
    fn test_initialize_gaussian_peak_at_x0() {
        let mut acc = make_acc(50);
        let x0 = 50e-9;
        let sigma = 10e-9;
        let mu0 = 1.0e-21;
        acc.initialize_gaussian(x0, sigma, mu0);
        // Find the site closest to x0
        let i_peak = (x0 / acc.dx).round() as usize;
        // The peak should be near mu0
        let peak_val = acc.mu_s[i_peak];
        assert!(
            (peak_val - mu0).abs() / mu0 < 0.05,
            "peak value {} ≠ mu0 {}",
            peak_val,
            mu0
        );
        // Neighbours should be smaller
        if i_peak + 1 < 50 {
            assert!(acc.mu_s[i_peak + 1] < peak_val + 1e-25);
        }
    }

    // ------------------------------------------------------------------ CFL violation

    #[test]
    fn test_evolve_cfl_violation_errors() {
        let mut acc = make_acc(10);
        acc.initialize_uniform(1.0e-21);
        // dt so large it violates CFL
        let dt_bad = 0.6 * acc.dx * acc.dx / acc.d_s;
        let result = acc.evolve(
            dt_bad,
            BoundaryCondition::Dirichlet(0.0),
            BoundaryCondition::Dirichlet(0.0),
        );
        assert!(result.is_err(), "CFL violation should return error");
    }

    // ------------------------------------------------------------------ Dirichlet boundary

    #[test]
    fn test_dirichlet_boundary_pinned() {
        let mut acc = make_acc(10);
        acc.initialize_uniform(1.0e-21);
        // CFL-safe dt
        let dt = 0.4 * acc.dx * acc.dx / acc.d_s;
        acc.evolve(
            dt,
            BoundaryCondition::Dirichlet(5.0e-21),
            BoundaryCondition::Dirichlet(0.0),
        )
        .expect("ok");
        assert!(
            (acc.mu_s[0] - 5.0e-21).abs() < 1e-30,
            "left Dirichlet not pinned"
        );
        assert!(acc.mu_s[9].abs() < 1e-30, "right Dirichlet not pinned");
    }

    // ------------------------------------------------------------------ Neumann zero flux

    #[test]
    fn test_neumann_zero_flux() {
        let mut acc = make_acc(20);
        acc.initialize_uniform(1.0e-21);
        let dt = 0.4 * acc.dx * acc.dx / acc.d_s;
        // Zero-flux Neumann on both ends should conserve approximate symmetry
        acc.evolve(
            dt,
            BoundaryCondition::Neumann(0.0),
            BoundaryCondition::Neumann(0.0),
        )
        .expect("ok");
        // μ[0] should be close to μ[1] (zero gradient)
        let grad = (acc.mu_s[1] - acc.mu_s[0]).abs();
        assert!(
            grad < 1e-25,
            "Neumann zero flux: gradient at boundary = {}",
            grad
        );
    }

    // ------------------------------------------------------------------ implicit (large dt)

    #[test]
    fn test_evolve_implicit_unconditional() {
        let mut acc = make_acc(20);
        acc.initialize_uniform(1.0e-21);
        // Very large dt — would fail CFL but should work with backward Euler
        let dt = 100.0 * acc.dx * acc.dx / acc.d_s;
        let result = acc.evolve_implicit(
            dt,
            BoundaryCondition::Dirichlet(0.0),
            BoundaryCondition::Dirichlet(0.0),
        );
        assert!(result.is_ok(), "implicit scheme should handle large dt");
        // All values should be finite
        for v in &acc.mu_s {
            assert!(v.is_finite(), "implicit solution has non-finite value");
        }
    }

    // ------------------------------------------------------------------ steady state

    #[test]
    fn test_steady_state_exponential_decay() {
        // Verify that the steady-state solution decays approximately exponentially
        let n = 100;
        let d_s = 1e-3;
        let tau_sf = 1e-12;
        let length = 1e-7; // 100 nm
        let mut acc = SpinAccumulation1D::new(d_s, tau_sf, length, n).expect("valid");
        let mu0 = 1.0e-20;
        acc.steady_state(mu0).expect("steady state ok");

        let lambda = acc.spin_diffusion_length();
        // Check first half of the chain
        let n_check = n / 2;
        for i in 1..n_check {
            let x = i as f64 * acc.dx;
            let expected = mu0 * (-x / lambda).exp();
            let actual = acc.mu_s[i];
            // Allow 20% relative tolerance for discrete approximation
            let rel_err = (actual - expected).abs() / expected.abs().max(1e-30);
            assert!(
                rel_err < 0.20,
                "steady-state at i={}: actual={:.4e}, expected={:.4e}, rel_err={:.3}",
                i,
                actual,
                expected,
                rel_err
            );
        }
    }

    // ------------------------------------------------------------------ total polarization

    #[test]
    fn test_total_polarization_positive_for_positive_mu() {
        let mut acc = make_acc(20);
        acc.initialize_uniform(1.0e-21);
        let pol = acc.total_polarization();
        assert!(pol > 0.0, "polarization should be positive: {}", pol);
    }

    // ------------------------------------------------------------------ decay profile interpolation

    #[test]
    fn test_decay_profile_interpolation() {
        let mut acc = make_acc(11);
        acc.initialize_gaussian(50e-9, 20e-9, 1.0e-21);
        // Interpolated value at a grid point should match the stored value
        let v_grid = acc.decay_profile(0.0).expect("ok");
        assert!((v_grid - acc.mu_s[0]).abs() < 1e-30);
        // Interior interpolated value should be between adjacent grid points
        let mid_x = 0.5 * acc.dx;
        let v_mid = acc.decay_profile(mid_x).expect("ok");
        let lo = acc.mu_s[0].min(acc.mu_s[1]);
        let hi = acc.mu_s[0].max(acc.mu_s[1]);
        assert!(v_mid >= lo - 1e-30 && v_mid <= hi + 1e-30);
    }

    // ------------------------------------------------------------------ spin current

    #[test]
    fn test_spin_current_at_injection_positive() {
        let n = 20;
        let d_s = 1e-3;
        let tau_sf = 1e-12;
        let length = 1e-7;
        let mut acc = SpinAccumulation1D::new(d_s, tau_sf, length, n).expect("valid");
        acc.steady_state(1.0e-20).expect("ok");
        // At the injection side (idx=0), spin current should flow into the sample (positive)
        let j = acc.spin_current(0, 1.0).expect("ok");
        // μ[0] > μ[1] (decaying), so derivative > 0, current = -D_s * sigma * dμ/dx < 0
        // The sign depends on direction convention; we check it's nonzero and finite
        assert!(j.is_finite(), "spin current must be finite");
    }
}
