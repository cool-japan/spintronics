//! Linear system solvers for FEM
//!
//! Provides iterative solvers for sparse linear systems arising from
//! finite element discretization. Uses scirs2-linalg for advanced solvers.

use crate::fem::assembly::SparseMatrix;

/// Solver type selection
#[derive(Debug, Clone, Copy)]
pub enum SolverType {
    /// Conjugate Gradient (for symmetric positive definite matrices)
    CG,
    /// BiConjugate Gradient Stabilized (for general non-symmetric matrices)
    BiCGSTAB,
    /// Simple iterative solver (Jacobi method)
    Jacobi,
    /// Successive Over-Relaxation (faster convergence than Jacobi)
    SOR,
}

/// Preconditioner type for iterative solvers
#[derive(Debug, Clone, Copy)]
pub enum Preconditioner {
    /// No preconditioning
    None,
    /// Jacobi (diagonal) preconditioner
    Jacobi,
    /// Symmetric Successive Over-Relaxation preconditioner
    SSOR,
}

/// Solver parameters
#[derive(Debug, Clone)]
pub struct SolverParams {
    /// Maximum number of iterations
    pub max_iter: usize,
    /// Convergence tolerance
    pub tolerance: f64,
    /// Relaxation parameter for SOR (typically 1.0 < omega < 2.0)
    pub omega: f64,
    /// Preconditioner to use
    pub preconditioner: Preconditioner,
}

impl Default for SolverParams {
    fn default() -> Self {
        Self {
            max_iter: 1000,
            tolerance: 1e-6,
            omega: 1.5,
            preconditioner: Preconditioner::None,
        }
    }
}

/// Apply preconditioner M^{-1} r
fn apply_preconditioner(
    a: &SparseMatrix,
    r: &[f64],
    precond: Preconditioner,
    omega: f64,
) -> Vec<f64> {
    match precond {
        Preconditioner::None => r.to_vec(),
        Preconditioner::Jacobi => apply_jacobi_preconditioner(a, r),
        Preconditioner::SSOR => apply_ssor_preconditioner(a, r, omega),
    }
}

/// Apply Jacobi (diagonal) preconditioner: z = D^{-1} r
fn apply_jacobi_preconditioner(a: &SparseMatrix, r: &[f64]) -> Vec<f64> {
    let n = a.nrows;
    let mut z = vec![0.0; n];

    // Extract diagonal
    let mut diag = vec![1.0; n];
    for ((&i, &j), &val) in a.rows.iter().zip(&a.cols).zip(&a.vals) {
        if i == j {
            diag[i] = val.max(1e-14); // Avoid division by zero
        }
    }

    // z = D^{-1} r
    for i in 0..n {
        z[i] = r[i] / diag[i];
    }

    z
}

/// Apply SSOR preconditioner
fn apply_ssor_preconditioner(a: &SparseMatrix, r: &[f64], omega: f64) -> Vec<f64> {
    let n = a.nrows;
    let mut z = vec![0.0; n];

    // Extract diagonal and lower triangular part
    let mut diag = vec![1.0; n];
    for ((&i, &j), &val) in a.rows.iter().zip(&a.cols).zip(&a.vals) {
        if i == j {
            diag[i] = val.max(1e-14);
        }
    }

    // Forward sweep: solve (D/omega + L) z_temp = r
    let mut z_temp = vec![0.0; n];
    for i in 0..n {
        let mut sum = 0.0;
        for ((&row, &col), &val) in a.rows.iter().zip(&a.cols).zip(&a.vals) {
            if row == i && col < i {
                sum += val * z_temp[col];
            }
        }
        z_temp[i] = omega * (r[i] - sum) / diag[i];
    }

    // Backward sweep: solve (D/omega + U) z = D z_temp
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for ((&row, &col), &val) in a.rows.iter().zip(&a.cols).zip(&a.vals) {
            if row == i && col > i {
                sum += val * z[col];
            }
        }
        z[i] = z_temp[i] - omega * sum / diag[i];
    }

    z
}

/// Solve linear system Ax = b using specified iterative solver
///
/// # Arguments
/// * `a` - System matrix (sparse)
/// * `b` - Right-hand side vector
/// * `solver_type` - Type of iterative solver to use
///
/// # Returns
/// Solution vector x
pub fn solve_linear_system(
    a: &SparseMatrix,
    b: &[f64],
    solver_type: SolverType,
) -> Result<Vec<f64>, String> {
    solve_linear_system_with_params(a, b, solver_type, SolverParams::default())
}

/// Solve linear system Ax = b with custom solver parameters
///
/// # Arguments
/// * `a` - System matrix (sparse)
/// * `b` - Right-hand side vector
/// * `solver_type` - Type of iterative solver to use
/// * `params` - Solver parameters (tolerance, max iterations, etc.)
///
/// # Returns
/// Solution vector x
pub fn solve_linear_system_with_params(
    a: &SparseMatrix,
    b: &[f64],
    solver_type: SolverType,
    params: SolverParams,
) -> Result<Vec<f64>, String> {
    if b.len() != a.nrows {
        return Err(format!(
            "Dimension mismatch: A is {}x{}, b has length {}",
            a.nrows,
            a.ncols,
            b.len()
        ));
    }

    match solver_type {
        SolverType::CG => solve_cg(a, b, &params),
        SolverType::BiCGSTAB => solve_bicgstab(a, b, &params),
        SolverType::Jacobi => solve_jacobi(a, b, &params),
        SolverType::SOR => solve_sor(a, b, &params),
    }
}

/// Preconditioned Conjugate Gradient solver for symmetric positive definite systems
fn solve_cg(a: &SparseMatrix, b: &[f64], params: &SolverParams) -> Result<Vec<f64>, String> {
    let n = a.nrows;
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();

    // Apply preconditioner to initial residual
    let mut z = apply_preconditioner(a, &r, params.preconditioner, params.omega);
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);

    for _iter in 0..params.max_iter {
        // Check convergence
        let r_norm = dot(&r, &r).sqrt();
        if r_norm < params.tolerance {
            break;
        }

        let ap = a.matvec(&p);
        let alpha = rz_old / dot(&p, &ap);

        // x = x + alpha * p
        for i in 0..n {
            x[i] += alpha * p[i];
        }

        // r = r - alpha * Ap
        for i in 0..n {
            r[i] -= alpha * ap[i];
        }

        // Apply preconditioner: z = M^{-1} r
        z = apply_preconditioner(a, &r, params.preconditioner, params.omega);

        let rz_new = dot(&r, &z);
        if rz_new.abs() < params.tolerance * params.tolerance {
            break;
        }

        let beta = rz_new / rz_old;

        // p = z + beta * p
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }

        rz_old = rz_new;
    }

    Ok(x)
}

/// BiConjugate Gradient Stabilized solver for general systems
fn solve_bicgstab(a: &SparseMatrix, b: &[f64], params: &SolverParams) -> Result<Vec<f64>, String> {
    let n = a.nrows;
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let r0 = r.clone();
    let mut rho = 1.0;
    let mut alpha = 1.0;
    let mut omega = 1.0;
    let mut v = vec![0.0; n];
    let mut p = vec![0.0; n];

    for _iter in 0..params.max_iter {
        let rho_new = dot(&r0, &r);

        if rho_new.abs() < params.tolerance {
            break;
        }

        let beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;

        // p = r + beta * (p - omega * v)
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        v = a.matvec(&p);
        alpha = rho / dot(&r0, &v);

        // s = r - alpha * v
        let mut s = vec![0.0; n];
        for i in 0..n {
            s[i] = r[i] - alpha * v[i];
        }

        // Check for early convergence
        if dot(&s, &s).sqrt() < params.tolerance {
            for i in 0..n {
                x[i] += alpha * p[i];
            }
            break;
        }

        let t = a.matvec(&s);
        omega = dot(&t, &s) / dot(&t, &t);

        // x = x + alpha * p + omega * s
        for i in 0..n {
            x[i] += alpha * p[i] + omega * s[i];
        }

        // r = s - omega * t
        for i in 0..n {
            r[i] = s[i] - omega * t[i];
        }

        if dot(&r, &r).sqrt() < params.tolerance {
            break;
        }
    }

    Ok(x)
}

/// Jacobi iterative solver
fn solve_jacobi(a: &SparseMatrix, b: &[f64], params: &SolverParams) -> Result<Vec<f64>, String> {
    let n = a.nrows;
    let mut x = vec![0.0; n];

    // Extract diagonal
    let mut diag = vec![1.0; n];
    for ((&i, &j), &val) in a.rows.iter().zip(&a.cols).zip(&a.vals) {
        if i == j {
            diag[i] = val.max(1e-10); // Avoid division by zero
        }
    }

    for _iter in 0..params.max_iter {
        let ax = a.matvec(&x);
        let mut residual_norm = 0.0;

        for i in 0..n {
            let r = b[i] - ax[i];
            x[i] += r / diag[i];
            residual_norm += r * r;
        }

        if residual_norm.sqrt() < params.tolerance {
            break;
        }
    }

    Ok(x)
}

/// Successive Over-Relaxation (SOR) solver
fn solve_sor(a: &SparseMatrix, b: &[f64], params: &SolverParams) -> Result<Vec<f64>, String> {
    let n = a.nrows;
    let mut x = vec![0.0; n];

    // Extract diagonal
    let mut diag = vec![1.0; n];
    for ((&i, &j), &val) in a.rows.iter().zip(&a.cols).zip(&a.vals) {
        if i == j {
            diag[i] = val.max(1e-10);
        }
    }

    for _iter in 0..params.max_iter {
        let mut residual_norm = 0.0;

        for i in 0..n {
            let ax_i = a.matvec_row(&x, i);
            let r = b[i] - ax_i;
            x[i] += params.omega * r / diag[i];
            residual_norm += r * r;
        }

        if residual_norm.sqrt() < params.tolerance {
            break;
        }
    }

    Ok(x)
}

/// Dot product of two vectors
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_simple() {
        // Diagonal system
        let mut a = SparseMatrix::new(3, 3);
        a.add_entry(0, 0, 2.0);
        a.add_entry(1, 1, 2.0);
        a.add_entry(2, 2, 2.0);

        let b = vec![2.0, 4.0, 6.0];
        let x = solve_linear_system(&a, &b, SolverType::Jacobi).unwrap();

        // Solution should be [1.0, 2.0, 3.0]
        assert!((x[0] - 1.0).abs() < 0.01);
        assert!((x[1] - 2.0).abs() < 0.01);
        assert!((x[2] - 3.0).abs() < 0.01);
    }

    #[test]
    fn test_solve_cg() {
        // Symmetric positive definite system
        let mut a = SparseMatrix::new(3, 3);
        a.add_entry(0, 0, 4.0);
        a.add_entry(0, 1, 1.0);
        a.add_entry(1, 0, 1.0);
        a.add_entry(1, 1, 3.0);
        a.add_entry(1, 2, 1.0);
        a.add_entry(2, 1, 1.0);
        a.add_entry(2, 2, 2.0);

        let b = vec![5.0, 5.0, 3.0];
        let x = solve_linear_system(&a, &b, SolverType::CG).unwrap();

        // Verify Ax = b
        let ax = a.matvec(&x);
        for i in 0..3 {
            assert!((ax[i] - b[i]).abs() < 0.01, "CG solver residual too large");
        }
    }

    #[test]
    fn test_solve_bicgstab() {
        // Non-symmetric system
        let mut a = SparseMatrix::new(3, 3);
        a.add_entry(0, 0, 3.0);
        a.add_entry(0, 1, 1.0);
        a.add_entry(1, 0, 2.0);
        a.add_entry(1, 1, 4.0);
        a.add_entry(1, 2, 1.0);
        a.add_entry(2, 1, 1.0);
        a.add_entry(2, 2, 2.0);

        let b = vec![4.0, 7.0, 3.0];
        let x = solve_linear_system(&a, &b, SolverType::BiCGSTAB).unwrap();

        // Verify Ax = b
        let ax = a.matvec(&x);
        for i in 0..3 {
            assert!(
                (ax[i] - b[i]).abs() < 0.1,
                "BiCGSTAB solver residual too large"
            );
        }
    }

    #[test]
    fn test_solve_sor() {
        // Test SOR solver
        let mut a = SparseMatrix::new(3, 3);
        a.add_entry(0, 0, 4.0);
        a.add_entry(0, 1, -1.0);
        a.add_entry(1, 0, -1.0);
        a.add_entry(1, 1, 4.0);
        a.add_entry(1, 2, -1.0);
        a.add_entry(2, 1, -1.0);
        a.add_entry(2, 2, 4.0);

        let b = vec![3.0, 2.0, 3.0];
        let x = solve_linear_system(&a, &b, SolverType::SOR).unwrap();

        // Verify Ax ≈ b
        let ax = a.matvec(&x);
        for i in 0..3 {
            assert!((ax[i] - b[i]).abs() < 0.1, "SOR solver residual too large");
        }
    }

    #[test]
    fn test_solver_params() {
        // Test custom solver parameters
        let params = SolverParams {
            max_iter: 100,
            tolerance: 1e-8,
            omega: 1.8,
            preconditioner: Preconditioner::None,
        };

        let mut a = SparseMatrix::new(2, 2);
        a.add_entry(0, 0, 2.0);
        a.add_entry(1, 1, 2.0);

        let b = vec![4.0, 6.0];
        let x = solve_linear_system_with_params(&a, &b, SolverType::Jacobi, params).unwrap();

        assert!((x[0] - 2.0).abs() < 1e-6);
        assert!((x[1] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_preconditioned_cg() {
        // Test CG with Jacobi preconditioner
        let mut a = SparseMatrix::new(3, 3);
        a.add_entry(0, 0, 4.0);
        a.add_entry(0, 1, 1.0);
        a.add_entry(1, 0, 1.0);
        a.add_entry(1, 1, 3.0);
        a.add_entry(1, 2, 1.0);
        a.add_entry(2, 1, 1.0);
        a.add_entry(2, 2, 2.0);

        let b = vec![5.0, 5.0, 3.0];

        let params_jacobi = SolverParams {
            max_iter: 100,
            tolerance: 1e-6,
            omega: 1.0,
            preconditioner: Preconditioner::Jacobi,
        };

        let x = solve_linear_system_with_params(&a, &b, SolverType::CG, params_jacobi).unwrap();

        // Verify Ax = b
        let ax = a.matvec(&x);
        for i in 0..3 {
            assert!(
                (ax[i] - b[i]).abs() < 0.01,
                "Preconditioned CG residual too large"
            );
        }
    }

    #[test]
    fn test_ssor_preconditioner() {
        // Test CG with SSOR preconditioner
        let mut a = SparseMatrix::new(3, 3);
        a.add_entry(0, 0, 4.0);
        a.add_entry(0, 1, -1.0);
        a.add_entry(1, 0, -1.0);
        a.add_entry(1, 1, 4.0);
        a.add_entry(1, 2, -1.0);
        a.add_entry(2, 1, -1.0);
        a.add_entry(2, 2, 4.0);

        let b = vec![3.0, 2.0, 3.0];

        let params_ssor = SolverParams {
            max_iter: 100,
            tolerance: 1e-6,
            omega: 1.2,
            preconditioner: Preconditioner::SSOR,
        };

        let x = solve_linear_system_with_params(&a, &b, SolverType::CG, params_ssor).unwrap();

        // Verify Ax ≈ b
        let ax = a.matvec(&x);
        for i in 0..3 {
            assert!(
                (ax[i] - b[i]).abs() < 0.01,
                "SSOR preconditioned CG residual too large"
            );
        }
    }
}
