//! Topological charge calculations
//!
//! Implements calculation of topological invariants like the skyrmion number.

use crate::vector3::Vector3;
use std::f64::consts::PI;

/// Topological charge calculator
pub struct TopologicalCharge;

impl TopologicalCharge {
    /// Calculate topological charge density at a point
    ///
    /// n(r) = (1/4π) m · (∂m/∂x × ∂m/∂y)
    ///
    /// # Arguments
    /// * `m` - Magnetization at the point
    /// * `dm_dx` - Spatial derivative in x direction
    /// * `dm_dy` - Spatial derivative in y direction
    ///
    /// # Returns
    /// Topological charge density [1/m²]
    pub fn density(m: Vector3<f64>, dm_dx: Vector3<f64>, dm_dy: Vector3<f64>) -> f64 {
        let cross = dm_dx.cross(&dm_dy);
        m.dot(&cross) / (4.0 * PI)
    }
}

/// Calculate skyrmion number for a 2D magnetization field
///
/// Numerically integrates the topological charge density over a region.
///
/// # Arguments
/// * `magnetization` - 2D array of magnetization vectors
/// * `dx` - Grid spacing in x direction \[m\]
/// * `dy` - Grid spacing in y direction \[m\]
///
/// # Returns
/// Skyrmion number (should be close to an integer)
pub fn calculate_skyrmion_number(magnetization: &[Vec<Vector3<f64>>], dx: f64, dy: f64) -> f64 {
    let nx = magnetization.len();
    if nx == 0 {
        return 0.0;
    }
    let ny = magnetization[0].len();
    if ny == 0 {
        return 0.0;
    }

    let mut total_charge = 0.0;

    // Use central differences for interior points
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let m = magnetization[i][j];

            // Calculate derivatives using central differences
            let dm_dx = (magnetization[i + 1][j] - magnetization[i - 1][j]) * (0.5 / dx);
            let dm_dy = (magnetization[i][j + 1] - magnetization[i][j - 1]) * (0.5 / dy);

            let charge_density = TopologicalCharge::density(m, dm_dx, dm_dy);

            total_charge += charge_density * dx * dy;
        }
    }

    total_charge
}

/// Calculate winding number along a closed path
///
/// Useful for characterizing domain walls and vortices.
pub fn winding_number(path: &[Vector3<f64>]) -> f64 {
    if path.len() < 3 {
        return 0.0;
    }

    let mut total_angle = 0.0;

    for i in 0..path.len() {
        let m1 = path[i];
        let m2 = path[(i + 1) % path.len()];

        // Calculate angle between consecutive magnetizations (in-plane)
        let angle1 = m1.y.atan2(m1.x);
        let angle2 = m2.y.atan2(m2.x);

        let mut delta = angle2 - angle1;

        // Unwrap angle
        if delta > PI {
            delta -= 2.0 * PI;
        } else if delta < -PI {
            delta += 2.0 * PI;
        }

        total_angle += delta;
    }

    total_angle / (2.0 * PI)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_charge_density_uniform() {
        let m = Vector3::new(0.0, 0.0, 1.0);
        let dm_dx = Vector3::new(0.0, 0.0, 0.0);
        let dm_dy = Vector3::new(0.0, 0.0, 0.0);

        let density = TopologicalCharge::density(m, dm_dx, dm_dy);
        assert!(density.abs() < 1e-20);
    }

    #[test]
    fn test_skyrmion_number_uniform() {
        // Create uniform magnetization
        let mag = vec![vec![Vector3::new(0.0, 0.0, 1.0); 10]; 10];

        let q = calculate_skyrmion_number(&mag, 1.0e-9, 1.0e-9);
        assert!(q.abs() < 0.1); // Should be ~0 for uniform
    }

    #[test]
    fn test_winding_number_vortex() {
        // Create a simple vortex path (circle in xy plane)
        let n_points = 8;
        let mut path = Vec::new();

        for i in 0..n_points {
            let angle = 2.0 * PI * (i as f64) / (n_points as f64);
            path.push(Vector3::new(angle.cos(), angle.sin(), 0.0));
        }

        let w = winding_number(&path);
        assert!((w - 1.0).abs() < 0.2); // Should be ~1 for vortex
    }
}
