//! Simplified Navier-Stokes solver for liquid metal flows
//!
//! Implements incompressible Navier-Stokes equations for spin-hydrodynamics

use crate::vector3::Vector3;

/// 2D fluid field
#[derive(Debug, Clone)]
pub struct FluidField {
    /// Velocity field \[m/s\]
    pub velocity: Vec<Vec<Vector3<f64>>>,

    /// Pressure field \[Pa\]
    pub pressure: Vec<Vec<f64>>,

    /// Grid dimension (x direction)
    pub nx: usize,
    /// Grid dimension (y direction)
    pub ny: usize,

    /// Grid spacing in x \[m\]
    pub dx: f64,
    /// Grid spacing in y \[m\]
    pub dy: f64,
}

impl FluidField {
    /// Create a new fluid field
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64) -> Self {
        Self {
            velocity: vec![vec![Vector3::new(0.0, 0.0, 0.0); ny]; nx],
            pressure: vec![vec![0.0; ny]; nx],
            nx,
            ny,
            dx,
            dy,
        }
    }

    /// Set initial velocity field
    pub fn set_velocity(&mut self, i: usize, j: usize, v: Vector3<f64>) {
        if i < self.nx && j < self.ny {
            self.velocity[i][j] = v;
        }
    }
}

/// Simplified Navier-Stokes solver
///
/// Solves incompressible N-S equations:
/// ∂v/∂t + (v·∇)v = -∇p/ρ + ν∇²v + f
/// ∇·v = 0 (incompressibility)
#[derive(Debug, Clone)]
pub struct NavierStokes {
    /// Kinematic viscosity [m²/s]
    /// Mercury: 1.15e-7, Galinstan: 3.4e-7
    pub viscosity: f64,

    /// Density [kg/m³]
    /// Mercury: 13534, Galinstan: 6440
    pub density: f64,

    /// Time step \[s\]
    pub dt: f64,
}

impl NavierStokes {
    /// Create solver for mercury
    pub fn mercury(dt: f64) -> Self {
        Self {
            viscosity: 1.15e-7,
            density: 13534.0,
            dt,
        }
    }

    /// Create solver for galinstan (liquid metal alloy)
    pub fn galinstan(dt: f64) -> Self {
        Self {
            viscosity: 3.4e-7,
            density: 6440.0,
            dt,
        }
    }

    /// Evolve velocity field by one time step (simplified Euler)
    ///
    /// This is a highly simplified solver for demonstration.
    /// Production code should use proper CFD methods (projection method, etc.)
    #[allow(clippy::needless_range_loop)]
    pub fn step(&self, field: &mut FluidField, external_force: Vector3<f64>) {
        let mut new_velocity = field.velocity.clone();

        for i in 1..field.nx - 1 {
            for j in 1..field.ny - 1 {
                let v = field.velocity[i][j];

                // Advection: -(v·∇)v (simplified)
                let advection = self.calculate_advection(field, i, j);

                // Diffusion: ν∇²v
                let diffusion = self.calculate_diffusion(field, i, j);

                // Pressure gradient: -∇p/ρ (simplified, assuming pressure equilibrium)
                let pressure_grad = Vector3::new(0.0, 0.0, 0.0);

                // Total acceleration
                let dv_dt = advection * (-1.0)
                    + diffusion * self.viscosity
                    + pressure_grad
                    + external_force;

                new_velocity[i][j] = v + dv_dt * self.dt;
            }
        }

        field.velocity = new_velocity;
    }

    /// Calculate advection term (v·∇)v
    #[allow(dead_code)]
    fn calculate_advection(&self, field: &FluidField, i: usize, j: usize) -> Vector3<f64> {
        let v = field.velocity[i][j];

        // Upwind scheme (very simple)
        let dv_dx = if v.x > 0.0 {
            (v - field.velocity[i - 1][j]) * (1.0 / field.dx)
        } else {
            (field.velocity[i + 1][j] - v) * (1.0 / field.dx)
        };

        let dv_dy = if v.y > 0.0 {
            (v - field.velocity[i][j - 1]) * (1.0 / field.dy)
        } else {
            (field.velocity[i][j + 1] - v) * (1.0 / field.dy)
        };

        Vector3::new(
            v.x * dv_dx.x + v.y * dv_dy.x,
            v.x * dv_dx.y + v.y * dv_dy.y,
            0.0,
        )
    }

    /// Calculate diffusion term ∇²v
    fn calculate_diffusion(&self, field: &FluidField, i: usize, j: usize) -> Vector3<f64> {
        let v = field.velocity[i][j];
        let v_xp = field.velocity[i + 1][j];
        let v_xm = field.velocity[i - 1][j];
        let v_yp = field.velocity[i][j + 1];
        let v_ym = field.velocity[i][j - 1];

        // Laplacian: ∇²v = ∂²v/∂x² + ∂²v/∂y²
        let d2v_dx2 = (v_xp + v_xm - v * 2.0) * (1.0 / field.dx.powi(2));
        let d2v_dy2 = (v_yp + v_ym - v * 2.0) * (1.0 / field.dy.powi(2));

        d2v_dx2 + d2v_dy2
    }

    /// Calculate Reynolds number
    ///
    /// Re = V L / ν
    ///
    /// where V is characteristic velocity and L is characteristic length
    pub fn reynolds_number(&self, velocity: f64, length: f64) -> f64 {
        velocity * length / self.viscosity
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fluid_field_creation() {
        let field = FluidField::new(10, 10, 1.0e-3, 1.0e-3);
        assert_eq!(field.nx, 10);
        assert_eq!(field.ny, 10);
    }

    #[test]
    fn test_navier_stokes_mercury() {
        let ns = NavierStokes::mercury(1.0e-6);
        assert!((ns.density - 13534.0).abs() < 1.0);
    }

    #[test]
    fn test_reynolds_number() {
        let ns = NavierStokes::galinstan(1.0e-6);
        let re = ns.reynolds_number(0.1, 1.0e-3); // 0.1 m/s, 1 mm

        assert!(re > 0.0);
        // For these values, should be laminar (Re < 2000)
        assert!(re < 2000.0);
    }

    #[test]
    fn test_diffusion_calculation() {
        let ns = NavierStokes::mercury(1.0e-6);
        let mut field = FluidField::new(5, 5, 1.0e-3, 1.0e-3);

        // Set a velocity bump in the center
        field.set_velocity(2, 2, Vector3::new(1.0, 0.0, 0.0));

        let diffusion = ns.calculate_diffusion(&field, 2, 2);

        // Diffusion should smooth out the bump (negative)
        assert!(diffusion.x < 0.0);
    }
}
