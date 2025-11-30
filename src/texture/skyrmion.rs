//! Magnetic skyrmions
//!
//! Skyrmions are topologically protected magnetic textures characterized by
//! a topological charge (skyrmion number).
//!
//! # Mathematical Formulation
//!
//! The topological charge (skyrmion number) is defined as:
//!
//! $$
//! Q = \frac{1}{4\pi} \int \mathbf{m} \cdot \left(\frac{\partial \mathbf{m}}{\partial x} \times \frac{\partial \mathbf{m}}{\partial y}\right) dx\,dy
//! $$
//!
//! where $\mathbf{m}(\mathbf{r})$ is the normalized magnetization field.
//!
//! For a single skyrmion, $Q = \pm 1$ depending on chirality.
//!
//! # Skyrmion Profile
//!
//! The magnetization distribution is typically parameterized as:
//!
//! $$
//! \mathbf{m}(r, \phi) = \left(\sin\Theta(r)\cos\Phi(\phi), \sin\Theta(r)\sin\Phi(\phi), \cos\Theta(r)\right)
//! $$
//!
//! where:
//! - $r$ is the radial distance from center
//! - $\phi$ is the azimuthal angle
//! - $\Theta(r)$ is the out-of-plane angle: $\Theta(r) = \pi \tanh\left(\frac{r-R}{\lambda}\right)$
//! - $\Phi(\phi)$ is the in-plane angle
//!
//! For Néel-type skyrmions: $\Phi(\phi) = \phi$ (radial)
//!
//! For Bloch-type skyrmions: $\Phi(\phi) = \phi + \pi/2$ (circular)
//!
//! # Stability
//!
//! Skyrmions are stabilized by the Dzyaloshinskii-Moriya interaction (DMI):
//!
//! $$
//! E_{\text{DMI}} = D \int \mathbf{m} \cdot (\nabla \times \mathbf{m}) \, dV
//! $$
//!
//! where $D$ is the DMI constant. The sign of $D$ determines the chirality.
//!
//! # References
//!
//! - T. H. R. Skyrme, "A unified field theory of mesons and baryons",
//!   Nucl. Phys. 31, 556 (1962)
//! - A. N. Bogdanov and D. A. Yablonskii, "Thermodynamically stable vortices
//!   in magnetically ordered crystals", Sov. Phys. JETP 68, 101 (1989)
//! - N. Nagaosa and Y. Tokura, "Topological properties and dynamics of
//!   magnetic skyrmions", Nat. Nanotechnol. 8, 899 (2013)

use crate::vector3::Vector3;
use std::f64::consts::PI;

/// Skyrmion helicity types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Helicity {
    /// Néel-type (radial)
    Neel,
    /// Bloch-type (circular)
    Bloch,
}

/// Chirality (handedness) of skyrmion
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Chirality {
    /// Clockwise rotation
    Clockwise,
    /// Counter-clockwise rotation
    CounterClockwise,
}

/// Single magnetic skyrmion
///
/// # Example
/// ```
/// use spintronics::texture::skyrmion::{Skyrmion, Helicity, Chirality};
///
/// // Create a Néel-type skyrmion at origin with 50 nm radius
/// let skyrmion = Skyrmion::new(
///     (0.0, 0.0),           // center position
///     50.0e-9,              // 50 nm radius
///     Helicity::Neel,       // Néel-type (radial)
///     Chirality::CounterClockwise,
/// );
///
/// // Calculate magnetization at skyrmion center
/// let m_center = skyrmion.magnetization_at(0.0, 0.0, 10.0e-9);
///
/// // At center, m_z should be negative (pointing down)
/// assert!(m_center.z < 0.0);
///
/// // Skyrmion number should be -1
/// assert_eq!(skyrmion.topological_charge, -1);
/// ```
#[derive(Debug, Clone)]
pub struct Skyrmion {
    /// Center position \[m\]
    pub center: (f64, f64),

    /// Skyrmion radius \[m\]
    pub radius: f64,

    /// Helicity type
    pub helicity: Helicity,

    /// Chirality
    pub chirality: Chirality,

    /// Topological charge (skyrmion number)
    /// Typically ±1 for skyrmions
    pub topological_charge: i32,
}

impl Skyrmion {
    /// Create a new skyrmion
    pub fn new(center: (f64, f64), radius: f64, helicity: Helicity, chirality: Chirality) -> Self {
        let topological_charge = match chirality {
            Chirality::CounterClockwise => -1,
            Chirality::Clockwise => 1,
        };

        Self {
            center,
            radius,
            helicity,
            chirality,
            topological_charge,
        }
    }

    /// Calculate magnetization at a given position
    ///
    /// Uses the standard skyrmion profile with a domain wall width parameter.
    ///
    /// # Arguments
    /// * `x`, `y` - Position \[m\]
    /// * `wall_width` - Domain wall width parameter \[m\]
    ///
    /// # Returns
    /// Magnetization vector (normalized)
    pub fn magnetization_at(&self, x: f64, y: f64, wall_width: f64) -> Vector3<f64> {
        // Distance from skyrmion center
        let dx = x - self.center.0;
        let dy = y - self.center.1;
        let r = (dx * dx + dy * dy).sqrt();

        // Polar angle: φ = atan2(y, x)
        // Physical meaning: Azimuthal position around the skyrmion center
        let phi = dy.atan2(dx);

        // Out-of-plane component: m_z(r)
        // Physical meaning: Skyrmions have m_z pointing DOWN at center (m_z = -1)
        // and UP at infinity (m_z = +1), with a smooth transition at radius R.
        // The tanh profile gives a domain wall: m_z(r) = tanh((r-R)/λ)
        // where λ is the wall width. This is energetically favorable due to
        // competition between exchange (prefers smooth) and DMI (prefers twist).
        let mz = if wall_width > 0.0 {
            ((r - self.radius) / wall_width).tanh()
        } else if r < self.radius {
            -1.0 // Sharp transition: down inside
        } else {
            1.0 // Sharp transition: up outside
        };

        // In-plane component: |m_xy| = √(1 - m_z²)
        // Physical meaning: Magnetization is always normalized |m| = 1, so when
        // m_z varies from -1 to +1, the in-plane component traces out a cone.
        // At the domain wall (m_z = 0), in-plane component is maximum (|m_xy| = 1).
        let m_inplane = (1.0 - mz * mz).sqrt().max(0.0);

        // In-plane angle: θ(φ)
        // Physical meaning: Determines how in-plane magnetization rotates around center
        //   - Néel type: θ = φ (radial, like wheel spokes pointing outward)
        //   - Bloch type: θ = φ + π/2 (circular, like vortex swirling)
        // Néel skyrmions are stabilized by interfacial DMI (e.g., Pt/Co/AlOx)
        // Bloch skyrmions appear in bulk materials with bulk DMI (e.g., MnSi, FeGe)
        let theta = match self.helicity {
            Helicity::Neel => phi,             // Radial (DMI at interface)
            Helicity::Bloch => phi + PI / 2.0, // Circular (bulk DMI)
        };

        // Chirality sign
        // Physical meaning: Determines rotation direction (handedness)
        //   - Counter-clockwise: left-handed rotation (Q = -1)
        //   - Clockwise: right-handed rotation (Q = +1)
        // Set by sign of DMI constant: D > 0 or D < 0
        let sign = match self.chirality {
            Chirality::CounterClockwise => 1.0,
            Chirality::Clockwise => -1.0,
        };

        // Final magnetization: m(r,φ) = (m_x, m_y, m_z)
        // where m_x = m_xy cos(θ), m_y = m_xy sin(θ)
        Vector3::new(
            m_inplane * (sign * theta).cos(),
            m_inplane * (sign * theta).sin(),
            mz,
        )
    }

    /// Calculate the topological charge density at a point
    ///
    /// n(r) = (1/4π) m · (∂m/∂x × ∂m/∂y)
    pub fn topological_charge_density(&self, x: f64, y: f64, wall_width: f64, dx: f64) -> f64 {
        // Calculate magnetization and derivatives numerically
        let m = self.magnetization_at(x, y, wall_width);
        let m_xp = self.magnetization_at(x + dx, y, wall_width);
        let m_yp = self.magnetization_at(x, y + dx, wall_width);

        let dm_dx = (m_xp - m) * (1.0 / dx);
        let dm_dy = (m_yp - m) * (1.0 / dx);

        let cross = dm_dx.cross(&dm_dy);
        m.dot(&cross) / (4.0 * PI)
    }
}

/// Skyrmion lattice (array of skyrmions)
#[derive(Debug, Clone)]
pub struct SkymionLattice {
    /// Skyrmions in the lattice
    pub skyrmions: Vec<Skyrmion>,

    /// Lattice constant \[m\]
    pub lattice_constant: f64,

    /// Lattice type
    pub lattice_type: LatticeType,
}

/// Lattice arrangement type
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LatticeType {
    /// Square lattice
    Square,
    /// Hexagonal (triangular) lattice
    Hexagonal,
}

impl SkymionLattice {
    /// Create a square lattice of skyrmions
    pub fn square(
        nx: usize,
        ny: usize,
        lattice_constant: f64,
        skyrmion_radius: f64,
        helicity: Helicity,
        chirality: Chirality,
    ) -> Self {
        let mut skyrmions = Vec::new();

        for i in 0..nx {
            for j in 0..ny {
                let x = i as f64 * lattice_constant;
                let y = j as f64 * lattice_constant;
                skyrmions.push(Skyrmion::new((x, y), skyrmion_radius, helicity, chirality));
            }
        }

        Self {
            skyrmions,
            lattice_constant,
            lattice_type: LatticeType::Square,
        }
    }

    /// Create a hexagonal lattice of skyrmions
    pub fn hexagonal(
        nx: usize,
        ny: usize,
        lattice_constant: f64,
        skyrmion_radius: f64,
        helicity: Helicity,
        chirality: Chirality,
    ) -> Self {
        let mut skyrmions = Vec::new();
        let dy = lattice_constant * (3.0_f64).sqrt() / 2.0;

        for i in 0..nx {
            for j in 0..ny {
                let x = i as f64 * lattice_constant + (j % 2) as f64 * lattice_constant / 2.0;
                let y = j as f64 * dy;
                skyrmions.push(Skyrmion::new((x, y), skyrmion_radius, helicity, chirality));
            }
        }

        Self {
            skyrmions,
            lattice_constant,
            lattice_type: LatticeType::Hexagonal,
        }
    }

    /// Calculate total topological charge
    pub fn total_topological_charge(&self) -> i32 {
        self.skyrmions.iter().map(|s| s.topological_charge).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_skyrmion_creation() {
        let sk = Skyrmion::new(
            (0.0, 0.0),
            10.0e-9,
            Helicity::Neel,
            Chirality::CounterClockwise,
        );

        assert_eq!(sk.topological_charge, -1);
    }

    #[test]
    fn test_magnetization_at_center() {
        let sk = Skyrmion::new(
            (0.0, 0.0),
            10.0e-9,
            Helicity::Neel,
            Chirality::CounterClockwise,
        );

        let m = sk.magnetization_at(0.0, 0.0, 2.0e-9);

        // At center, should point down
        assert!(m.z < 0.0);

        // Should be normalized
        assert!((m.magnitude() - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_magnetization_far_from_center() {
        let sk = Skyrmion::new(
            (0.0, 0.0),
            10.0e-9,
            Helicity::Neel,
            Chirality::CounterClockwise,
        );

        let m = sk.magnetization_at(100.0e-9, 0.0, 2.0e-9);

        // Far away, should point up
        assert!(m.z > 0.5);
    }

    #[test]
    fn test_square_lattice() {
        let lattice = SkymionLattice::square(
            3,
            3,
            50.0e-9,
            10.0e-9,
            Helicity::Neel,
            Chirality::CounterClockwise,
        );

        assert_eq!(lattice.skyrmions.len(), 9);
        assert_eq!(lattice.total_topological_charge(), -9);
    }

    #[test]
    fn test_hexagonal_lattice() {
        let lattice = SkymionLattice::hexagonal(
            4,
            4,
            50.0e-9,
            10.0e-9,
            Helicity::Bloch,
            Chirality::Clockwise,
        );

        assert_eq!(lattice.skyrmions.len(), 16);
        assert_eq!(lattice.lattice_type, LatticeType::Hexagonal);
    }
}
