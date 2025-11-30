//! Finite element definitions and shape functions
//!
//! Provides linear triangular and tetrahedral elements
//! with Lagrange basis functions.

use crate::vector3::Vector3;

/// Linear triangular element (P1)
///
/// Uses 3-node linear basis functions:
/// N₁(ξ,η) = 1 - ξ - η
/// N₂(ξ,η) = ξ
/// N₃(ξ,η) = η
#[derive(Debug, Clone)]
pub struct TriangularElement {
    /// Node coordinates
    pub nodes: [Vector3<f64>; 3],
}

impl TriangularElement {
    /// Create new triangular element
    pub fn new(nodes: [Vector3<f64>; 3]) -> Self {
        Self { nodes }
    }

    /// Evaluate shape functions at point (ξ, η)
    pub fn shape_functions(&self, xi: f64, eta: f64) -> [f64; 3] {
        [1.0 - xi - eta, xi, eta]
    }

    /// Shape function gradients in physical coordinates
    ///
    /// Returns ∇N in (x,y) coordinates
    pub fn shape_gradients(&self) -> [[f64; 2]; 3] {
        let p0 = self.nodes[0];
        let p1 = self.nodes[1];
        let p2 = self.nodes[2];

        // Jacobian matrix
        let dx_dxi = p1.x - p0.x;
        let dy_dxi = p1.y - p0.y;
        let dx_deta = p2.x - p0.x;
        let dy_deta = p2.y - p0.y;

        let det_j = dx_dxi * dy_deta - dx_deta * dy_dxi;

        // Inverse Jacobian
        let inv_j11 = dy_deta / det_j;
        let inv_j12 = -dx_deta / det_j;
        let inv_j21 = -dy_dxi / det_j;
        let inv_j22 = dx_dxi / det_j;

        // Gradients in reference coordinates
        let dn_dxi = [-1.0, 1.0, 0.0];
        let dn_deta = [-1.0, 0.0, 1.0];

        // Transform to physical coordinates
        let mut gradients = [[0.0; 2]; 3];
        for i in 0..3 {
            gradients[i][0] = inv_j11 * dn_dxi[i] + inv_j12 * dn_deta[i];
            gradients[i][1] = inv_j21 * dn_dxi[i] + inv_j22 * dn_deta[i];
        }

        gradients
    }

    /// Element area
    pub fn area(&self) -> f64 {
        let v1 = self.nodes[1] - self.nodes[0];
        let v2 = self.nodes[2] - self.nodes[0];
        0.5 * v1.cross(&v2).magnitude()
    }
}

/// Linear tetrahedral element (P1)
///
/// Uses 4-node linear basis functions
#[derive(Debug, Clone)]
pub struct TetrahedralElement {
    /// Node coordinates
    pub nodes: [Vector3<f64>; 4],
}

impl TetrahedralElement {
    /// Create new tetrahedral element
    pub fn new(nodes: [Vector3<f64>; 4]) -> Self {
        Self { nodes }
    }

    /// Evaluate shape functions at point (ξ, η, ζ)
    pub fn shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> [f64; 4] {
        [1.0 - xi - eta - zeta, xi, eta, zeta]
    }

    /// Shape function gradients in physical coordinates
    pub fn shape_gradients(&self) -> [[f64; 3]; 4] {
        let p0 = self.nodes[0];
        let p1 = self.nodes[1];
        let p2 = self.nodes[2];
        let p3 = self.nodes[3];

        // Jacobian matrix
        let j = [
            [p1.x - p0.x, p2.x - p0.x, p3.x - p0.x],
            [p1.y - p0.y, p2.y - p0.y, p3.y - p0.y],
            [p1.z - p0.z, p2.z - p0.z, p3.z - p0.z],
        ];

        // Determinant
        let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
            - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
            + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);

        // Inverse Jacobian (simplified for 3x3)
        let inv_j = [
            [
                (j[1][1] * j[2][2] - j[1][2] * j[2][1]) / det,
                (j[0][2] * j[2][1] - j[0][1] * j[2][2]) / det,
                (j[0][1] * j[1][2] - j[0][2] * j[1][1]) / det,
            ],
            [
                (j[1][2] * j[2][0] - j[1][0] * j[2][2]) / det,
                (j[0][0] * j[2][2] - j[0][2] * j[2][0]) / det,
                (j[0][2] * j[1][0] - j[0][0] * j[1][2]) / det,
            ],
            [
                (j[1][0] * j[2][1] - j[1][1] * j[2][0]) / det,
                (j[0][1] * j[2][0] - j[0][0] * j[2][1]) / det,
                (j[0][0] * j[1][1] - j[0][1] * j[1][0]) / det,
            ],
        ];

        // Gradients in reference coordinates
        let dn_dxi = [-1.0, 1.0, 0.0, 0.0];
        let dn_deta = [-1.0, 0.0, 1.0, 0.0];
        let dn_dzeta = [-1.0, 0.0, 0.0, 1.0];

        // Transform to physical coordinates
        let mut gradients = [[0.0; 3]; 4];
        #[allow(clippy::needless_range_loop)]
        for i in 0..4 {
            for j_idx in 0..3 {
                gradients[i][j_idx] = inv_j[j_idx][0] * dn_dxi[i]
                    + inv_j[j_idx][1] * dn_deta[i]
                    + inv_j[j_idx][2] * dn_dzeta[i];
            }
        }

        gradients
    }

    /// Element volume
    pub fn volume(&self) -> f64 {
        let v1 = self.nodes[1] - self.nodes[0];
        let v2 = self.nodes[2] - self.nodes[0];
        let v3 = self.nodes[3] - self.nodes[0];

        // Volume = 1/6 * |triple product|
        (1.0 / 6.0) * v1.cross(&v2).dot(&v3).abs()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_triangular_element() {
        let nodes = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let elem = TriangularElement::new(nodes);

        // Area should be 0.5
        assert!((elem.area() - 0.5).abs() < 1e-10);

        // Shape functions at center
        let n = elem.shape_functions(1.0 / 3.0, 1.0 / 3.0);
        for &val in &n {
            assert!((val - 1.0 / 3.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_tetrahedral_element() {
        let nodes = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let elem = TetrahedralElement::new(nodes);

        // Volume should be 1/6
        assert!((elem.volume() - 1.0 / 6.0).abs() < 1e-10);

        // Shape functions at center
        let n = elem.shape_functions(0.25, 0.25, 0.25);
        for &val in &n {
            assert!((val - 0.25).abs() < 1e-10);
        }
    }
}
