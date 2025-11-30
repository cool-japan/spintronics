//! Vorticity calculation for fluid fields
//!
//! Vorticity ω = ∇ × v measures local rotation in a fluid.

use crate::vector3::Vector3;

/// Vorticity calculator for discretized velocity fields
pub struct VorticityCalculator {
    /// Grid spacing in x direction \[m\]
    pub dx: f64,

    /// Grid spacing in y direction \[m\]
    pub dy: f64,

    /// Grid spacing in z direction \[m\]
    pub dz: f64,
}

impl VorticityCalculator {
    /// Create a new vorticity calculator
    pub fn new(dx: f64, dy: f64, dz: f64) -> Self {
        Self { dx, dy, dz }
    }

    /// Create for uniform grid
    pub fn uniform(grid_spacing: f64) -> Self {
        Self::new(grid_spacing, grid_spacing, grid_spacing)
    }

    /// Calculate vorticity at a point in a 3D velocity field
    ///
    /// ω = ∇ × v = (∂v_z/∂y - ∂v_y/∂z, ∂v_x/∂z - ∂v_z/∂x, ∂v_y/∂x - ∂v_x/∂y)
    ///
    /// # Arguments
    /// * `velocity_field` - 3D array of velocity vectors
    /// * `i`, `j`, `k` - Grid indices
    ///
    /// # Returns
    /// Vorticity vector [1/s]
    pub fn calculate_3d(
        &self,
        velocity_field: &[Vec<Vec<Vector3<f64>>>],
        i: usize,
        j: usize,
        k: usize,
    ) -> Vector3<f64> {
        let nx = velocity_field.len();
        if nx == 0 || i == 0 || i >= nx - 1 {
            return Vector3::new(0.0, 0.0, 0.0);
        }

        let ny = velocity_field[0].len();
        if ny == 0 || j == 0 || j >= ny - 1 {
            return Vector3::new(0.0, 0.0, 0.0);
        }

        let nz = velocity_field[0][0].len();
        if nz == 0 || k == 0 || k >= nz - 1 {
            return Vector3::new(0.0, 0.0, 0.0);
        }

        // Central differences for derivatives
        let _v = velocity_field[i][j][k];

        // ∂v_z/∂y
        let dvz_dy =
            (velocity_field[i][j + 1][k].z - velocity_field[i][j - 1][k].z) / (2.0 * self.dy);
        // ∂v_y/∂z
        let dvy_dz =
            (velocity_field[i][j][k + 1].y - velocity_field[i][j][k - 1].y) / (2.0 * self.dz);

        // ∂v_x/∂z
        let dvx_dz =
            (velocity_field[i][j][k + 1].x - velocity_field[i][j][k - 1].x) / (2.0 * self.dz);
        // ∂v_z/∂x
        let dvz_dx =
            (velocity_field[i + 1][j][k].z - velocity_field[i - 1][j][k].z) / (2.0 * self.dx);

        // ∂v_y/∂x
        let dvy_dx =
            (velocity_field[i + 1][j][k].y - velocity_field[i - 1][j][k].y) / (2.0 * self.dx);
        // ∂v_x/∂y
        let dvx_dy =
            (velocity_field[i][j + 1][k].x - velocity_field[i][j - 1][k].x) / (2.0 * self.dy);

        Vector3::new(dvz_dy - dvy_dz, dvx_dz - dvz_dx, dvy_dx - dvx_dy)
    }

    /// Calculate vorticity for a 2D flow (z-component only)
    ///
    /// ω_z = ∂v_y/∂x - ∂v_x/∂y
    ///
    /// This is the most common case for planar flows
    pub fn calculate_2d(&self, velocity_field: &[Vec<Vector3<f64>>], i: usize, j: usize) -> f64 {
        let nx = velocity_field.len();
        if nx == 0 || i == 0 || i >= nx - 1 {
            return 0.0;
        }

        let ny = velocity_field[0].len();
        if ny == 0 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // ∂v_y/∂x
        let dvy_dx = (velocity_field[i + 1][j].y - velocity_field[i - 1][j].y) / (2.0 * self.dx);

        // ∂v_x/∂y
        let dvx_dy = (velocity_field[i][j + 1].x - velocity_field[i][j - 1].x) / (2.0 * self.dy);

        dvy_dx - dvx_dy
    }

    /// Calculate circulation around a closed loop
    ///
    /// Γ = ∮ v · dl
    ///
    /// By Stokes' theorem: Γ = ∫∫ (∇ × v) · dA = ∫∫ ω · dA
    pub fn circulation(&self, velocity_path: &[(Vector3<f64>, Vector3<f64>)]) -> f64 {
        let mut circulation = 0.0;

        for (v, dl) in velocity_path {
            circulation += v.dot(dl);
        }

        circulation
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculator_creation() {
        let calc = VorticityCalculator::uniform(1.0e-3);
        assert_eq!(calc.dx, 1.0e-3);
        assert_eq!(calc.dy, 1.0e-3);
        assert_eq!(calc.dz, 1.0e-3);
    }

    #[test]
    fn test_2d_vorticity_rigid_rotation() {
        let calc = VorticityCalculator::uniform(0.1);

        // Create a rigidly rotating velocity field: v = ω × r
        // For rotation about z-axis: v_x = -ω y, v_y = ω x
        let omega = 2.0; // rad/s
        let mut field = vec![vec![Vector3::new(0.0, 0.0, 0.0); 5]; 5];

        #[allow(clippy::needless_range_loop)]
        for i in 0..5 {
            for j in 0..5 {
                let x = (i as f64 - 2.0) * 0.1;
                let y = (j as f64 - 2.0) * 0.1;
                field[i][j] = Vector3::new(-omega * y, omega * x, 0.0);
            }
        }

        let vort = calc.calculate_2d(&field, 2, 2);

        // Rigid rotation should give vorticity = 2ω
        assert!((vort - 2.0 * omega).abs() < 0.5);
    }

    #[test]
    fn test_zero_vorticity_uniform_flow() {
        let calc = VorticityCalculator::uniform(0.1);

        // Uniform flow has zero vorticity
        let field = vec![vec![Vector3::new(1.0, 0.0, 0.0); 5]; 5];

        let vort = calc.calculate_2d(&field, 2, 2);
        assert!(vort.abs() < 1e-10);
    }
}
