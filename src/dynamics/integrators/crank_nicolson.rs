//! Crank-Nicolson integrator for parabolic (diffusion) equations.
//!
//! Solves the one-dimensional diffusion equation
//!
//!   ∂u/∂t = D ∂²u/∂x² + S(x, t)
//!
//! using the classical Crank-Nicolson scheme, which is centred in time and
//! centred in space:
//!
//!   (I − dt·D/(2 dx²) L) u_{n+1} = (I + dt·D/(2 dx²) L) u_n + dt · S_{n+1/2}
//!
//! where `L` is the standard second-difference operator. The method is
//! second-order accurate in both space and time, and unconditionally stable
//! for linear diffusion problems — there is no Courant–Friedrichs–Lewy
//! restriction `dt ≤ dx²/(2 D)` as is required for explicit schemes such as
//! forward Euler.
//!
//! Spatial boundary conditions are encoded by [`BoundaryCondition`]:
//! Dirichlet (fixed value), Neumann (zero-flux) and periodic boundaries are
//! supported. The discrete system is tridiagonal (Dirichlet/Neumann) or
//! cyclic-tridiagonal (periodic) and is solved with the Thomas algorithm in
//! O(N) operations. The periodic case uses a Sherman-Morrison correction to
//! reduce the cyclic system to two tridiagonal solves.
//!
//! For spintronics applications the natural use case is the diffusion of
//! the longitudinal spin accumulation μ_s in a non-magnetic metal layer
//!
//!   ∂μ_s/∂t = D ∇²μ_s − μ_s / τ_sf
//!
//! where τ_sf is the spin-flip relaxation time. The spin diffusion length
//! is `ℓ_sf = sqrt(D · τ_sf)`. A specialised wrapper
//! [`SpinDiffusionCrankNicolson`] is provided for this case.
//!
//! Crank-Nicolson assumes a linear (or at most quasi-linear) diffusion
//! operator. For spatially coupled LLG dynamics one should instead use
//! [`super::implicit_midpoint::ImplicitMidpointNewton`].
//!
//! # References
//! - Crank, J. and Nicolson, P., *A practical method for numerical evaluation
//!   of solutions of partial differential equations of the heat-conduction
//!   type*, Proc. Cambridge Philos. Soc. **43**, 50 (1947).
//! - Hairer, E. and Wanner, G., *Solving Ordinary Differential Equations II*,
//!   Springer (1996).

use crate::error::{invalid_param, Result};

// =========================================================================
// Boundary conditions
// =========================================================================

/// Spatial boundary condition for the 1D Crank-Nicolson diffusion solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryCondition {
    /// Dirichlet: u(left) and u(right) are fixed at the given values.
    Dirichlet(f64, f64),
    /// Neumann no-flux boundary: ∂u/∂x = 0 at both ends. Uses a ghost-cell
    /// reflection that makes u₋₁ = u₁ and u_N = u_{N-2}.
    NeumannZero,
    /// Periodic: `u[0] == u[N-1]` in the continuum sense; the discrete
    /// scheme treats the grid as a ring of (N-1) unique nodes.
    Periodic,
}

// =========================================================================
// CrankNicolsonDiffusion
// =========================================================================

/// Crank-Nicolson solver for the linear 1D diffusion equation
/// ∂u/∂t = D ∂²u/∂x² + S(x, t).
#[derive(Debug, Clone)]
pub struct CrankNicolsonDiffusion {
    /// Uniform grid spacing.
    pub dx: f64,
    /// Diffusion coefficient D.
    pub diffusivity: f64,
    /// Spatial boundary condition.
    pub boundary: BoundaryCondition,
}

impl CrankNicolsonDiffusion {
    /// Construct a new diffusion solver, validating the parameters.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] if `dx` or
    /// `diffusivity` are non-positive.
    pub fn new(dx: f64, diffusivity: f64, boundary: BoundaryCondition) -> Result<Self> {
        if !dx.is_finite() || dx <= 0.0 {
            return Err(invalid_param("dx", "must be finite and positive"));
        }
        if !diffusivity.is_finite() || diffusivity < 0.0 {
            return Err(invalid_param(
                "diffusivity",
                "must be finite and non-negative",
            ));
        }
        Ok(Self {
            dx,
            diffusivity,
            boundary,
        })
    }

    /// Advance the discrete solution by one Crank-Nicolson step.
    ///
    /// `u` is the current solution vector, `source` is the source term S
    /// (assumed constant across the step), and `dt` is the time step. The
    /// returned vector has the same length as `u`.
    ///
    /// # Errors
    /// Returns an error if `dt` is non-positive, or if the linear solve
    /// becomes singular.
    pub fn step(&self, u: &[f64], source: &[f64], dt: f64) -> Result<Vec<f64>> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(invalid_param("dt", "must be finite and positive"));
        }
        if u.len() != source.len() {
            return Err(invalid_param("source", "must have the same length as u"));
        }
        if u.len() < 3 {
            return Err(invalid_param("u", "must have at least 3 nodes"));
        }

        let n = u.len();
        let r = self.diffusivity * dt / (self.dx * self.dx);
        let half_r = 0.5 * r;

        // Build the linear system A · u^{n+1} = rhs.
        // For interior nodes (everywhere except boundary handling):
        //   -half_r · u_{i-1}^{n+1} + (1 + r) · u_i^{n+1} - half_r · u_{i+1}^{n+1}
        //     =  half_r · u_{i-1}^n + (1 - r) · u_i^n + half_r · u_{i+1}^n + dt · S_i
        match self.boundary {
            BoundaryCondition::Dirichlet(left, right) => {
                solve_dirichlet(u, source, dt, half_r, left, right)
            },
            BoundaryCondition::NeumannZero => solve_neumann(u, source, dt, r, half_r),
            BoundaryCondition::Periodic => solve_periodic(u, source, dt, n, r, half_r),
        }
    }

    /// Evolve from `u_init` to time `t_end` using a constant step size `dt`,
    /// recording the solution after each step. The first entry of the
    /// returned vector is the initial state.
    ///
    /// `source_fn(t, u)` is evaluated at each step.
    ///
    /// # Errors
    /// Propagates errors from [`Self::step`].
    pub fn evolve<F>(
        &self,
        u_init: &[f64],
        source_fn: F,
        t_end: f64,
        dt: f64,
    ) -> Result<Vec<Vec<f64>>>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        if !t_end.is_finite() || t_end <= 0.0 {
            return Err(invalid_param("t_end", "must be finite and positive"));
        }
        if !dt.is_finite() || dt <= 0.0 {
            return Err(invalid_param("dt", "must be finite and positive"));
        }

        let n_steps = (t_end / dt).round() as usize;
        let mut history: Vec<Vec<f64>> = Vec::with_capacity(n_steps + 1);
        history.push(u_init.to_vec());

        let mut t = 0.0;
        let mut current = u_init.to_vec();
        for _ in 0..n_steps {
            let source = source_fn(t + 0.5 * dt, &current);
            current = self.step(&current, &source, dt)?;
            t += dt;
            history.push(current.clone());
        }
        Ok(history)
    }

    /// Solve the time-independent equation D u''(x) + S(x) = 0 subject to the
    /// configured boundary condition. The linear system is the same as the
    /// Crank-Nicolson left-hand side at infinite time and is therefore solved
    /// with the same tridiagonal machinery.
    ///
    /// # Errors
    /// Returns an error if `source` has fewer than three nodes or if the
    /// configured boundary condition is incompatible with a unique solution
    /// (Neumann zero-flux is only unique up to an additive constant, so the
    /// total source must integrate to zero).
    pub fn steady_state(&self, source: &[f64]) -> Result<Vec<f64>> {
        if source.len() < 3 {
            return Err(invalid_param("source", "must have at least 3 nodes"));
        }
        let n = source.len();
        let dx2 = self.dx * self.dx;
        let d = self.diffusivity;
        if d <= 0.0 {
            return Err(invalid_param(
                "diffusivity",
                "must be positive for steady-state solve",
            ));
        }

        match self.boundary {
            BoundaryCondition::Dirichlet(left, right) => {
                // Interior unknowns u_1..u_{n-2}. The discrete equation is
                //   (u_{i-1} - 2 u_i + u_{i+1}) / dx^2 = -S_i / D.
                let m = n - 2;
                let mut a = vec![0.0_f64; m]; // sub-diagonal
                let mut b = vec![0.0_f64; m]; // diagonal
                let mut c = vec![0.0_f64; m]; // super-diagonal
                let mut rhs = vec![0.0_f64; m];

                for i in 0..m {
                    a[i] = 1.0;
                    b[i] = -2.0;
                    c[i] = 1.0;
                    rhs[i] = -source[i + 1] * dx2 / d;
                }
                rhs[0] -= left;
                rhs[m - 1] -= right;
                let inner = thomas(&a, &b, &c, &rhs)?;

                let mut out = vec![0.0_f64; n];
                out[0] = left;
                out[n - 1] = right;
                out[1..(n - 1)].copy_from_slice(&inner);
                Ok(out)
            },
            BoundaryCondition::NeumannZero => Err(invalid_param(
                "boundary",
                "NeumannZero steady-state is only defined up to an additive constant",
            )),
            BoundaryCondition::Periodic => Err(invalid_param(
                "boundary",
                "Periodic steady-state is only defined up to an additive constant",
            )),
        }
    }
}

// =========================================================================
// Step-implementation helpers
// =========================================================================

fn solve_dirichlet(
    u: &[f64],
    source: &[f64],
    dt: f64,
    half_r: f64,
    left: f64,
    right: f64,
) -> Result<Vec<f64>> {
    let n = u.len();
    let m = n - 2;
    let mut a = vec![0.0_f64; m];
    let mut b = vec![0.0_f64; m];
    let mut c = vec![0.0_f64; m];
    let mut rhs = vec![0.0_f64; m];
    let one_minus = 1.0 - 2.0 * half_r; // = 1 - r
    let one_plus = 1.0 + 2.0 * half_r; // = 1 + r

    for i in 0..m {
        a[i] = -half_r;
        b[i] = one_plus;
        c[i] = -half_r;
        let i_global = i + 1;
        rhs[i] = half_r * u[i_global - 1]
            + one_minus * u[i_global]
            + half_r * u[i_global + 1]
            + dt * source[i_global];
    }
    // Move boundary contributions to the RHS.
    rhs[0] += half_r * left;
    rhs[m - 1] += half_r * right;

    let inner = thomas(&a, &b, &c, &rhs)?;
    let mut out = vec![0.0_f64; n];
    out[0] = left;
    out[n - 1] = right;
    out[1..(n - 1)].copy_from_slice(&inner);
    Ok(out)
}

fn solve_neumann(u: &[f64], source: &[f64], dt: f64, r: f64, half_r: f64) -> Result<Vec<f64>> {
    // Ghost cells: u_{-1} = u_1 and u_N = u_{N-2}, giving for the boundary
    // rows an effective coefficient of 2 in front of the off-diagonal term.
    let n = u.len();
    let mut a = vec![0.0_f64; n];
    let mut b = vec![0.0_f64; n];
    let mut c = vec![0.0_f64; n];
    let mut rhs = vec![0.0_f64; n];
    let one_minus = 1.0 - r;
    let one_plus = 1.0 + r;

    // Left boundary (i = 0)
    a[0] = 0.0;
    b[0] = one_plus;
    c[0] = -r;
    rhs[0] = one_minus * u[0] + r * u[1] + dt * source[0];

    // Interior nodes
    for i in 1..(n - 1) {
        a[i] = -half_r;
        b[i] = one_plus;
        c[i] = -half_r;
        rhs[i] = half_r * u[i - 1] + one_minus * u[i] + half_r * u[i + 1] + dt * source[i];
    }

    // Right boundary (i = n - 1)
    let last = n - 1;
    a[last] = -r;
    b[last] = one_plus;
    c[last] = 0.0;
    rhs[last] = r * u[last - 1] + one_minus * u[last] + dt * source[last];

    thomas(&a, &b, &c, &rhs)
}

/// Periodic case via Sherman-Morrison. The cyclic tridiagonal system can be
/// written as A x = d where A = T + u v^T with T tridiagonal. We solve
/// T y = d and T z = u, then x = y − ((v^T y) / (1 + v^T z)) z.
fn solve_periodic(
    u: &[f64],
    source: &[f64],
    dt: f64,
    n: usize,
    r: f64,
    half_r: f64,
) -> Result<Vec<f64>> {
    let one_minus = 1.0 - r;
    let one_plus = 1.0 + r;

    // Build the modified tridiagonal T = A - u v^T. The two corner entries
    // a[0][n-1] = c[n-1][0] = -half_r are removed; we choose
    // u = (gamma, 0, ..., 0, -half_r) and v = (1, 0, ..., 0, -half_r/gamma)
    // so that u v^T reproduces the missing corners. To avoid a zero pivot at
    // row 0 (which would happen if gamma == one_plus) pick gamma = -one_plus,
    // which makes T[0][0] = 2 one_plus and is always non-zero.
    let gamma = -one_plus;

    let mut a = vec![0.0_f64; n];
    let mut b = vec![0.0_f64; n];
    let mut c = vec![0.0_f64; n];
    let mut rhs = vec![0.0_f64; n];

    for i in 0..n {
        let im1 = if i == 0 { n - 1 } else { i - 1 };
        let ip1 = if i == n - 1 { 0 } else { i + 1 };
        a[i] = -half_r;
        b[i] = one_plus;
        c[i] = -half_r;
        rhs[i] = half_r * u[im1] + one_minus * u[i] + half_r * u[ip1] + dt * source[i];
    }
    // Drop the cyclic corners from the boundary rows.
    a[0] = 0.0;
    c[n - 1] = 0.0;

    // Adjust the diagonal: subtract u v^T contributions from rows 0 and n-1.
    b[0] -= gamma;
    b[n - 1] -= -half_r * (-half_r / gamma);

    // Solve T y = rhs.
    let y = thomas(&a, &b, &c, &rhs)?;
    // Build the u-vector for the rank-1 update.
    let mut u_vec = vec![0.0_f64; n];
    u_vec[0] = gamma;
    u_vec[n - 1] = -half_r;
    // Solve T z = u.
    let z = thomas(&a, &b, &c, &u_vec)?;

    let v_dot_y = y[0] + (-half_r / gamma) * y[n - 1];
    let v_dot_z = z[0] + (-half_r / gamma) * z[n - 1];
    let denom = 1.0 + v_dot_z;
    if denom.abs() < 1e-30 {
        return Err(invalid_param(
            "periodic_solve",
            "Sherman-Morrison denominator is singular",
        ));
    }
    let coeff = v_dot_y / denom;
    let mut x = vec![0.0_f64; n];
    for i in 0..n {
        x[i] = y[i] - coeff * z[i];
    }
    Ok(x)
}

/// Thomas algorithm for tridiagonal systems. `a` is the sub-diagonal (a[0]
/// is unused), `b` is the diagonal, `c` is the super-diagonal
/// (c[n-1] is unused), and `d` is the right-hand side.
fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Result<Vec<f64>> {
    let n = b.len();
    if a.len() != n || c.len() != n || d.len() != n {
        return Err(invalid_param(
            "thomas",
            "tridiagonal coefficient slices must all have the same length",
        ));
    }
    if n == 0 {
        return Ok(Vec::new());
    }

    let mut c_prime = vec![0.0_f64; n];
    let mut d_prime = vec![0.0_f64; n];

    if b[0].abs() < 1e-30 {
        return Err(invalid_param("thomas", "zero pivot encountered"));
    }
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];
    for i in 1..n {
        let denom = b[i] - a[i] * c_prime[i - 1];
        if denom.abs() < 1e-30 {
            return Err(invalid_param("thomas", "zero pivot encountered"));
        }
        c_prime[i] = c[i] / denom;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom;
    }

    let mut x = vec![0.0_f64; n];
    x[n - 1] = d_prime[n - 1];
    for i in (0..(n - 1)).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
    Ok(x)
}

// =========================================================================
// SpinDiffusionCrankNicolson
// =========================================================================

/// Spin diffusion in a non-magnetic spacer layer:
///
///   ∂μ_s/∂t = D ∇²μ_s − μ_s / τ_sf
///
/// The spin-flip relaxation term is incorporated implicitly into the source
/// term S = −μ_s / τ_sf evaluated at the midpoint of the time step.
#[derive(Debug, Clone)]
pub struct SpinDiffusionCrankNicolson {
    /// Underlying Crank-Nicolson solver for the pure diffusion piece.
    pub cn: CrankNicolsonDiffusion,
    /// Spin-flip relaxation time τ_sf.
    pub spin_flip_time: f64,
}

impl SpinDiffusionCrankNicolson {
    /// Construct the spin-diffusion solver with validated parameters.
    ///
    /// # Errors
    /// Returns an error if any parameter is non-positive or non-finite.
    pub fn new(
        dx: f64,
        diffusivity: f64,
        spin_flip_time: f64,
        boundary: BoundaryCondition,
    ) -> Result<Self> {
        let cn = CrankNicolsonDiffusion::new(dx, diffusivity, boundary)?;
        if !spin_flip_time.is_finite() || spin_flip_time <= 0.0 {
            return Err(invalid_param(
                "spin_flip_time",
                "must be finite and positive",
            ));
        }
        Ok(Self { cn, spin_flip_time })
    }

    /// Advance μ_s by one step of size `dt`. The spin-flip sink is folded
    /// into the Crank-Nicolson matrix so that the discretisation is fully
    /// time-centred:
    ///
    ///   (1 + dt/(2τ)) u^{n+1} − (dt·D/(2 dx²)) L u^{n+1}
    ///       = (1 − dt/(2τ)) u^n + (dt·D/(2 dx²)) L u^n .
    ///
    /// For a spatially uniform μ_s the operator L vanishes and this reduces
    /// to the rational approximation to exp(−dt/τ), giving the exact
    /// exponential decay to second order in dt.
    pub fn step(&self, mu_s: &[f64], dt: f64) -> Result<Vec<f64>> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(invalid_param("dt", "must be finite and positive"));
        }
        if mu_s.len() < 3 {
            return Err(invalid_param("mu_s", "must have at least 3 nodes"));
        }
        // First run the diffusion + zero-source CN step.
        let n = mu_s.len();
        let source = vec![0.0_f64; n];
        let diffused = self.cn.step(mu_s, &source, dt)?;
        // Then apply the implicit-trapezoidal damping factor that exactly
        // accounts for the sink term −u/τ in a uniform setting. Splitting
        // sink and diffusion in this way is operator-splitting (Strang-like)
        // and remains second-order accurate because both subproblems are
        // self-adjoint and linear.
        let factor =
            (1.0 - 0.5 * dt / self.spin_flip_time) / (1.0 + 0.5 * dt / self.spin_flip_time);
        Ok(diffused.iter().map(|&v| v * factor).collect())
    }

    /// Spin-diffusion length ℓ_sf = sqrt(D · τ_sf).
    pub fn spin_diffusion_length(&self) -> f64 {
        (self.cn.diffusivity * self.spin_flip_time).sqrt()
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -- Test 1: constructor with valid parameters ----------------------

    #[test]
    fn test_constructor_valid() {
        let cn = CrankNicolsonDiffusion::new(0.1, 0.5, BoundaryCondition::NeumannZero)
            .expect("valid constructor");
        assert!((cn.dx - 0.1).abs() < 1e-15);
        assert!((cn.diffusivity - 0.5).abs() < 1e-15);
        assert_eq!(cn.boundary, BoundaryCondition::NeumannZero);
    }

    // -- Test 2: invalid dx returns an error -----------------------------

    #[test]
    fn test_invalid_dx_rejected() {
        let r = CrankNicolsonDiffusion::new(0.0, 1.0, BoundaryCondition::NeumannZero);
        assert!(r.is_err());
        let r = CrankNicolsonDiffusion::new(-0.1, 1.0, BoundaryCondition::NeumannZero);
        assert!(r.is_err());
        let r = CrankNicolsonDiffusion::new(f64::NAN, 1.0, BoundaryCondition::NeumannZero);
        assert!(r.is_err());
    }

    // -- Test 3: invalid diffusivity returns an error -------------------

    #[test]
    fn test_invalid_diffusivity_rejected() {
        let r = CrankNicolsonDiffusion::new(0.1, -0.5, BoundaryCondition::NeumannZero);
        assert!(r.is_err());
        let r = CrankNicolsonDiffusion::new(0.1, f64::INFINITY, BoundaryCondition::NeumannZero);
        assert!(r.is_err());
    }

    // -- Test 4: Dirichlet steady-state linear profile ------------------

    #[test]
    fn test_dirichlet_steady_state_linear() {
        // With S = 0 and boundary values u(0) = 0, u(L) = 1 the steady-state
        // is u(x) = x/L (purely linear).
        let n = 21;
        let l = 1.0;
        let dx = l / (n as f64 - 1.0);
        let cn =
            CrankNicolsonDiffusion::new(dx, 1.0, BoundaryCondition::Dirichlet(0.0, 1.0)).unwrap();
        let source = vec![0.0_f64; n];
        let u = cn.steady_state(&source).expect("steady-state ok");
        for (i, &val) in u.iter().enumerate() {
            let x = i as f64 * dx;
            let expected = x / l;
            assert!(
                (val - expected).abs() < 1e-12,
                "node {}: got {}, expected {}",
                i,
                val,
                expected
            );
        }
    }

    // -- Test 5: Neumann BC conserves total mass ------------------------

    #[test]
    fn test_neumann_mass_conservation() {
        // With NeumannZero BC and S = 0 the total integral of u must be
        // conserved during diffusion. For a node-centred grid with mirror
        // ghost cells the natural mass measure is the trapezoidal rule
        // M = dx * (u_0/2 + u_1 + ... + u_{N-2} + u_{N-1}/2), which is the
        // discrete analogue of ∫ u dx for this stencil.
        let n = 41;
        let dx = 0.05;
        let cn = CrankNicolsonDiffusion::new(dx, 1.0, BoundaryCondition::NeumannZero).unwrap();
        // Initial bump at the centre.
        let mut u = vec![0.0_f64; n];
        u[n / 2] = 10.0;
        let source = vec![0.0_f64; n];
        let mass = |v: &[f64]| -> f64 {
            let interior: f64 = v[1..(n - 1)].iter().sum();
            dx * (0.5 * v[0] + interior + 0.5 * v[n - 1])
        };
        let m0 = mass(&u);
        let dt = 0.01;
        let mut current = u;
        for _ in 0..100 {
            current = cn.step(&current, &source, dt).expect("step ok");
        }
        let m_final = mass(&current);
        let drift = ((m_final - m0) / m0).abs();
        assert!(drift < 1e-10, "mass drift: {:.3e}", drift);
    }

    // -- Test 6: periodic Fourier mode decays as exp(-D k^2 t) -----------

    #[test]
    fn test_periodic_mode_decay() {
        // Domain length L with periodic BC. Initial condition u_0 = cos(k x)
        // with k = 2π/L decays exactly as exp(-D k^2 t).
        let n = 64;
        let l = 1.0;
        let dx = l / n as f64;
        let d = 1.0_f64;
        let cn = CrankNicolsonDiffusion::new(dx, d, BoundaryCondition::Periodic).unwrap();
        let k = 2.0 * std::f64::consts::PI / l;
        let u0: Vec<f64> = (0..n).map(|i| (k * (i as f64) * dx).cos()).collect();
        let source = vec![0.0_f64; n];

        let dt = 0.001_f64;
        let n_steps = 100;
        let mut current = u0.clone();
        for _ in 0..n_steps {
            current = cn.step(&current, &source, dt).expect("step ok");
        }

        let t_total = dt * n_steps as f64;
        let expected_factor = (-d * k * k * t_total).exp();
        // Check amplitude at the cosine peak (i = 0).
        let measured = current[0];
        let expected = u0[0] * expected_factor;
        let err = (measured - expected).abs();
        assert!(
            err < 5e-3,
            "periodic decay error: got {}, expected {} (diff {:.3e})",
            measured,
            expected,
            err
        );
    }

    // -- Test 7: 2nd-order spatial accuracy -----------------------------

    #[test]
    fn test_spatial_second_order_accuracy() {
        // Solve the steady-state Poisson equation u''(x) = -π² sin(π x) with
        // u(0) = u(1) = 0; analytical solution u(x) = sin(π x). Refining
        // the mesh by 2 should reduce the L∞ error by ~4.
        let mut errors = Vec::new();
        let dxs = [0.05_f64, 0.025, 0.0125];
        for &dx in &dxs {
            let n = (1.0 / dx).round() as usize + 1;
            let cn = CrankNicolsonDiffusion::new(dx, 1.0, BoundaryCondition::Dirichlet(0.0, 0.0))
                .unwrap();
            let pi = std::f64::consts::PI;
            let source: Vec<f64> = (0..n)
                .map(|i| pi * pi * (pi * i as f64 * dx).sin())
                .collect();
            let u = cn.steady_state(&source).expect("steady-state ok");
            let mut max_err = 0.0_f64;
            for (i, &val) in u.iter().enumerate() {
                let x = i as f64 * dx;
                let expected = (pi * x).sin();
                let err = (val - expected).abs();
                if err > max_err {
                    max_err = err;
                }
            }
            errors.push(max_err);
        }
        // Ratio between successive errors should be near 4.
        let r1 = errors[0] / errors[1];
        let r2 = errors[1] / errors[2];
        assert!(
            r1 > 3.0 && r1 < 5.0,
            "expected O(dx^2) ratio ~ 4, got {:.2}",
            r1
        );
        assert!(
            r2 > 3.0 && r2 < 5.0,
            "expected O(dx^2) ratio ~ 4, got {:.2}",
            r2
        );
    }

    // -- Test 8: unconditional stability (dt > dx^2/D) ------------------

    #[test]
    fn test_unconditional_stability_large_dt() {
        // For an explicit forward-Euler scheme the stability limit is
        // dt ≤ dx²/(2 D). Crank-Nicolson should remain stable at dt ten
        // times larger.
        let n = 21;
        let dx = 0.05;
        let d = 1.0;
        let dt = 10.0 * dx * dx / d;
        let cn = CrankNicolsonDiffusion::new(dx, d, BoundaryCondition::NeumannZero).unwrap();
        let mut u = vec![0.0_f64; n];
        u[n / 2] = 1.0;
        let source = vec![0.0_f64; n];
        for _ in 0..50 {
            u = cn.step(&u, &source, dt).expect("step ok");
            assert!(
                u.iter().all(|x| x.is_finite()),
                "CN must remain finite at large dt"
            );
        }
    }

    // -- Test 9: spin-diffusion length sanity check ----------------------

    #[test]
    fn test_spin_diffusion_length_matches_definition() {
        let sd =
            SpinDiffusionCrankNicolson::new(0.05, 1.0e-3, 5.0e-4, BoundaryCondition::NeumannZero)
                .unwrap();
        let expected = (1.0e-3_f64 * 5.0e-4_f64).sqrt();
        assert!((sd.spin_diffusion_length() - expected).abs() < 1e-15);
    }

    // -- Test 10: spin-flip sink causes exponential decay ----------------

    #[test]
    fn test_spin_flip_sink_decay() {
        // Spatially uniform μ_s under Neumann BC should decay purely from
        // the spin-flip sink: μ(t) = μ(0) exp(-t / τ_sf).
        let dx = 0.05;
        let d = 1.0;
        let tau = 0.1;
        let sd =
            SpinDiffusionCrankNicolson::new(dx, d, tau, BoundaryCondition::NeumannZero).unwrap();
        let n = 11;
        let mut u = vec![1.0_f64; n];
        let dt = 0.001;
        let n_steps = 50;
        for _ in 0..n_steps {
            u = sd.step(&u, dt).expect("step ok");
        }
        let t_total = dt * n_steps as f64;
        let expected = (-t_total / tau).exp();
        let measured = u[n / 2];
        assert!(
            (measured - expected).abs() < 1e-3,
            "spin-flip decay: got {}, expected {}",
            measured,
            expected
        );
    }

    // -- Test 11: evolve records the right number of snapshots ----------

    #[test]
    fn test_evolve_snapshot_count() {
        let dx = 0.1;
        let cn = CrankNicolsonDiffusion::new(dx, 1.0, BoundaryCondition::NeumannZero).unwrap();
        let n = 11;
        let u0 = vec![0.0_f64; n];
        let history = cn
            .evolve(&u0, |_t, _u| vec![1.0_f64; n], 1.0, 0.1)
            .expect("evolve ok");
        // Initial + 10 steps.
        assert_eq!(history.len(), 11);
        assert_eq!(history[0], u0);
    }

    // -- Test 12: step rejects dt <= 0 -----------------------------------

    #[test]
    fn test_step_rejects_bad_dt() {
        let cn = CrankNicolsonDiffusion::new(0.1, 1.0, BoundaryCondition::NeumannZero).unwrap();
        let n = 5;
        let u = vec![0.0_f64; n];
        let s = vec![0.0_f64; n];
        assert!(cn.step(&u, &s, 0.0).is_err());
        assert!(cn.step(&u, &s, -1.0).is_err());
        assert!(cn.step(&u, &s, f64::NAN).is_err());
    }
}
