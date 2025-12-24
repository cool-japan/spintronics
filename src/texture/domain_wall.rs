//! Magnetic domain walls
//!
//! Domain walls are interfaces between regions of different magnetization.
//! They are characterized by their width, energy, and type (Bloch vs Néel).

use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::vector3::Vector3;

/// Domain wall type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum WallType {
    /// Bloch wall (magnetization rotates in-plane perpendicular to wall normal)
    Bloch,
    /// Néel wall (magnetization rotates in-plane parallel to wall normal)
    Neel,
}

/// Magnetic domain wall
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DomainWall {
    /// Position of wall center \[m\]
    pub center: f64,

    /// Wall width parameter \[m\]
    pub width: f64,

    /// Wall type
    pub wall_type: WallType,

    /// Propagation direction normal to wall
    pub normal: Vector3<f64>,
}

impl Default for DomainWall {
    /// Default domain wall: Bloch type with 10 nm width at origin
    fn default() -> Self {
        Self::new(0.0, 10.0e-9, WallType::Bloch)
    }
}

impl DomainWall {
    /// Create a new domain wall
    pub fn new(center: f64, width: f64, wall_type: WallType) -> Self {
        Self {
            center,
            width,
            wall_type,
            normal: Vector3::new(1.0, 0.0, 0.0),
        }
    }

    /// Calculate domain wall width from material parameters
    ///
    /// δ = √(A/K)
    ///
    /// where A is exchange stiffness and K is anisotropy
    pub fn calculate_width(exchange_a: f64, anisotropy_k: f64) -> f64 {
        (exchange_a / anisotropy_k).sqrt()
    }

    /// Get magnetization profile through the wall
    ///
    /// # Arguments
    /// * `x` - Position along wall normal \[m\]
    ///
    /// # Returns
    /// Magnetization vector (normalized)
    pub fn magnetization_at(&self, x: f64) -> Vector3<f64> {
        // Distance from wall center
        let xi = (x - self.center) / self.width;

        // Out-of-plane angle (for head-to-head walls)
        // θ(x) = 2 arctan(exp(x/δ))
        let theta = 2.0 * xi.exp().atan();

        match self.wall_type {
            WallType::Bloch => {
                // Bloch wall: rotation in yz plane
                Vector3::new(theta.cos(), theta.sin(), 0.0)
            },
            WallType::Neel => {
                // Néel wall: rotation in xz plane
                Vector3::new(theta.cos(), 0.0, theta.sin())
            },
        }
    }

    /// Calculate domain wall energy per unit area \[J/m²\]
    ///
    /// σ = 4√(AK)
    pub fn energy_density(exchange_a: f64, anisotropy_k: f64) -> f64 {
        4.0 * (exchange_a * anisotropy_k).sqrt()
    }

    /// Calculate velocity under spin-transfer torque
    ///
    /// v = (μ_B P / e M_s δ) × j
    ///
    /// where P is spin polarization, j is current density
    #[allow(dead_code)]
    pub fn velocity_stt(
        &self,
        current_density: f64,
        spin_polarization: f64,
        saturation_mag: f64,
    ) -> f64 {
        use crate::constants::{E_CHARGE, MU_B};

        let prefactor = MU_B * spin_polarization / (E_CHARGE * saturation_mag * self.width);
        prefactor * current_density
    }

    /// Builder method to set wall center position
    pub fn with_center(mut self, center: f64) -> Self {
        self.center = center;
        self
    }

    /// Builder method to set wall width
    pub fn with_width(mut self, width: f64) -> Self {
        self.width = width;
        self
    }

    /// Builder method to set wall type
    pub fn with_type(mut self, wall_type: WallType) -> Self {
        self.wall_type = wall_type;
        self
    }

    /// Builder method to set normal direction
    pub fn with_normal(mut self, normal: Vector3<f64>) -> Self {
        self.normal = normal.normalize();
        self
    }
}

impl fmt::Display for WallType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WallType::Bloch => write!(f, "Bloch"),
            WallType::Neel => write!(f, "Néel"),
        }
    }
}

impl fmt::Display for DomainWall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "DomainWall[{}]: δ={:.1} nm, center={:.1} nm",
            self.wall_type,
            self.width * 1e9,
            self.center * 1e9
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wall_width_calculation() {
        let a_ex = 1.0e-11; // J/m
        let k_u = 1.0e5; // J/m³

        let width = DomainWall::calculate_width(a_ex, k_u);
        assert!(width > 0.0);
        assert!(width < 1.0e-6); // Should be nm scale
    }

    #[test]
    fn test_magnetization_profile() {
        let wall = DomainWall::new(0.0, 10.0e-9, WallType::Bloch);

        // Far left: should point in one direction
        let m_left = wall.magnetization_at(-100.0e-9);
        assert!(m_left.x.abs() > 0.9);

        // Far right: should point in opposite direction
        let m_right = wall.magnetization_at(100.0e-9);
        assert!((m_right.x - (-1.0)).abs() < 0.1 || (m_right.x - 1.0).abs() < 0.1);

        // At center: should be rotated
        let m_center = wall.magnetization_at(0.0);
        assert!(m_center.y.abs() > 0.5 || m_center.z.abs() > 0.5);
    }

    #[test]
    fn test_wall_energy() {
        let a_ex = 1.0e-11;
        let k_u = 1.0e5;

        let energy = DomainWall::energy_density(a_ex, k_u);
        assert!(energy > 0.0);
        assert!(energy < 1.0); // Should be mJ/m² scale
    }

    #[test]
    fn test_default_domain_wall() {
        let dw = DomainWall::default();

        // Default should be Bloch type
        assert_eq!(dw.wall_type, WallType::Bloch);
        // Default width should be 10 nm
        assert_eq!(dw.width, 10.0e-9);
        // Default center should be at origin
        assert_eq!(dw.center, 0.0);
        // Default normal should be along x
        assert_eq!(dw.normal.x, 1.0);
        assert_eq!(dw.normal.y, 0.0);
        assert_eq!(dw.normal.z, 0.0);
    }

    #[test]
    fn test_wall_types() {
        let bloch = DomainWall::new(0.0, 10.0e-9, WallType::Bloch);
        let neel = DomainWall::new(0.0, 10.0e-9, WallType::Neel);

        let m_bloch = bloch.magnetization_at(0.0);
        let m_neel = neel.magnetization_at(0.0);

        // Bloch wall has y-component, Néel has z-component
        assert!(m_bloch.y.abs() > m_neel.y.abs());
        assert!(m_neel.z.abs() > m_bloch.z.abs());
    }
}
