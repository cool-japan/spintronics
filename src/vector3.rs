//! Simple 3D vector implementation for spintronics
//!
//! This is a temporary internal implementation until scirs2-linalg is available.

use std::fmt;
use std::ops::{Add, Mul, Sub};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A 3D vector with generic type T
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Vector3<T> {
    /// X component
    pub x: T,
    /// Y component
    pub y: T,
    /// Z component
    pub z: T,
}

impl<T> Vector3<T> {
    /// Create a new 3D vector
    pub const fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }
}

impl Vector3<f64> {
    /// Create a zero vector
    pub const fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Create a unit vector along the x-axis
    pub const fn unit_x() -> Self {
        Self {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Create a unit vector along the y-axis
    pub const fn unit_y() -> Self {
        Self {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }
    }

    /// Create a unit vector along the z-axis
    pub const fn unit_z() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        }
    }

    /// Calculate the dot product with another vector
    #[inline]
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Calculate the cross product with another vector
    #[inline]
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Calculate the magnitude (length) of the vector
    #[inline]
    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Calculate the squared magnitude (avoids sqrt for performance)
    #[inline]
    pub fn magnitude_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Return a normalized (unit) vector in the same direction
    #[inline]
    pub fn normalize(&self) -> Self {
        let mag = self.magnitude();
        if mag > 0.0 {
            Self {
                x: self.x / mag,
                y: self.y / mag,
                z: self.z / mag,
            }
        } else {
            *self
        }
    }

    /// Check if the vector is normalized (unit length)
    ///
    /// Uses a tolerance of 1e-10 for floating point comparison
    #[inline]
    pub fn is_normalized(&self) -> bool {
        (self.magnitude_squared() - 1.0).abs() < 1e-10
    }

    /// Calculate angle between two vectors in radians
    ///
    /// Returns angle in range \[0, π\]
    #[inline]
    pub fn angle_between(&self, other: &Self) -> f64 {
        let dot = self.dot(other);
        let mags = self.magnitude() * other.magnitude();
        if mags > 0.0 {
            (dot / mags).clamp(-1.0, 1.0).acos()
        } else {
            0.0
        }
    }

    /// Project this vector onto another vector
    ///
    /// Returns the component of self in the direction of other
    #[inline]
    pub fn project(&self, other: &Self) -> Self {
        let other_mag_sq = other.magnitude_squared();
        if other_mag_sq > 0.0 {
            *other * (self.dot(other) / other_mag_sq)
        } else {
            Self::zero()
        }
    }
}

impl Add for Vector3<f64> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vector3<f64> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<f64> for Vector3<f64> {
    type Output = Self;

    fn mul(self, scalar: f64) -> Self {
        Self {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
        }
    }
}

impl fmt::Display for Vector3<f64> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.6}, {:.6}, {:.6})", self.x, self.y, self.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cross_product() {
        let v1 = Vector3::new(1.0, 0.0, 0.0);
        let v2 = Vector3::new(0.0, 1.0, 0.0);
        let result = v1.cross(&v2);
        assert!((result.x - 0.0).abs() < 1e-10);
        assert!((result.y - 0.0).abs() < 1e-10);
        assert!((result.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize() {
        let v = Vector3::new(3.0, 4.0, 0.0);
        let normalized = v.normalize();
        assert!((normalized.magnitude() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_zero() {
        let v = Vector3::zero();
        assert_eq!(v.x, 0.0);
        assert_eq!(v.y, 0.0);
        assert_eq!(v.z, 0.0);
    }

    #[test]
    fn test_unit_vectors() {
        let ux = Vector3::unit_x();
        let uy = Vector3::unit_y();
        let uz = Vector3::unit_z();

        assert_eq!(ux.x, 1.0);
        assert_eq!(ux.y, 0.0);
        assert_eq!(ux.z, 0.0);

        assert_eq!(uy.x, 0.0);
        assert_eq!(uy.y, 1.0);
        assert_eq!(uy.z, 0.0);

        assert_eq!(uz.x, 0.0);
        assert_eq!(uz.y, 0.0);
        assert_eq!(uz.z, 1.0);

        // All should be normalized
        assert!(ux.is_normalized());
        assert!(uy.is_normalized());
        assert!(uz.is_normalized());
    }

    #[test]
    fn test_magnitude_squared() {
        let v = Vector3::new(3.0, 4.0, 0.0);
        assert_eq!(v.magnitude_squared(), 25.0);
        assert!((v.magnitude() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_is_normalized() {
        let v1 = Vector3::new(1.0, 0.0, 0.0);
        assert!(v1.is_normalized());

        let v2 = Vector3::new(3.0, 4.0, 0.0);
        assert!(!v2.is_normalized());

        let v3 = v2.normalize();
        assert!(v3.is_normalized());
    }

    #[test]
    fn test_angle_between() {
        let v1 = Vector3::new(1.0, 0.0, 0.0);
        let v2 = Vector3::new(0.0, 1.0, 0.0);

        // Should be 90 degrees (π/2 radians)
        let angle = v1.angle_between(&v2);
        assert!((angle - std::f64::consts::FRAC_PI_2).abs() < 1e-10);

        // Parallel vectors
        let v3 = Vector3::new(2.0, 0.0, 0.0);
        let angle2 = v1.angle_between(&v3);
        assert!(angle2.abs() < 1e-10);

        // Anti-parallel vectors
        let v4 = Vector3::new(-1.0, 0.0, 0.0);
        let angle3 = v1.angle_between(&v4);
        assert!((angle3 - std::f64::consts::PI).abs() < 1e-10);
    }

    #[test]
    fn test_project() {
        let v1 = Vector3::new(3.0, 4.0, 0.0);
        let v2 = Vector3::new(1.0, 0.0, 0.0);

        // Project v1 onto v2 (should give (3, 0, 0))
        let proj = v1.project(&v2);
        assert!((proj.x - 3.0).abs() < 1e-10);
        assert!(proj.y.abs() < 1e-10);
        assert!(proj.z.abs() < 1e-10);

        // Project onto diagonal
        let v3 = Vector3::new(1.0, 1.0, 0.0);
        let proj2 = v1.project(&v3);
        // Projection should be (3.5, 3.5, 0)
        assert!((proj2.x - 3.5).abs() < 1e-10);
        assert!((proj2.y - 3.5).abs() < 1e-10);
        assert!(proj2.z.abs() < 1e-10);
    }
}
