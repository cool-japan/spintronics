//! Spin interface properties
//!
//! This module implements the spin-mixing conductance model that is central
//! to Prof. Saitoh's research on spin pumping phenomena.

use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::vector3::Vector3;

/// Spin interface between ferromagnet and normal metal
///
/// The spin-mixing conductance g_r is a key parameter introduced in:
/// Saitoh et al., "Conversion of spin current into charge current at room
/// temperature: Inverse spin-Hall effect", Appl. Phys. Lett. 88, 182509 (2006)
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpinInterface {
    /// Real part of spin-mixing conductance [Ω⁻¹ m⁻²]
    ///
    /// This parameter determines the efficiency of spin current transmission
    /// across the ferromagnet/normal-metal interface.
    /// Typical values: 1e18 to 1e20 for YIG/Pt, Py/Pt interfaces
    pub g_r: f64,

    /// Imaginary part of spin-mixing conductance [Ω⁻¹ m⁻²]
    ///
    /// Related to interface-induced magnetic anisotropy
    pub g_i: f64,

    /// Interface normal vector (normalized)
    ///
    /// Points from ferromagnet to normal metal
    pub normal: Vector3<f64>,

    /// Interface area \[m²\]
    pub area: f64,
}

impl Default for SpinInterface {
    fn default() -> Self {
        Self {
            g_r: 1.0e19,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12, // 1 μm²
        }
    }
}

impl SpinInterface {
    /// Create a YIG/Pt interface
    ///
    /// Parameters based on experimental values from Saitoh group's research
    pub fn yig_pt() -> Self {
        Self {
            g_r: 1.0e19,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12,
        }
    }

    /// Create a Permalloy/Pt interface
    pub fn py_pt() -> Self {
        Self {
            g_r: 5.0e19,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12,
        }
    }

    /// Create a CoFeB/Pt interface
    ///
    /// Common in perpendicular magnetic anisotropy (PMA) systems
    pub fn cofeb_pt() -> Self {
        Self {
            g_r: 4.0e19,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12,
        }
    }

    /// Create a Co/Pt interface
    ///
    /// Strong PMA and spin-orbit coupling
    pub fn co_pt() -> Self {
        Self {
            g_r: 6.0e19,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12,
        }
    }

    /// Create a Fe/Pt interface
    ///
    /// Large spin polarization
    pub fn fe_pt() -> Self {
        Self {
            g_r: 5.5e19,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12,
        }
    }

    /// Create a YIG/Ta interface
    ///
    /// Tantalum with negative spin Hall angle
    pub fn yig_ta() -> Self {
        Self {
            g_r: 9.0e18,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12,
        }
    }

    /// Create a Permalloy/Ta interface
    ///
    /// Commonly used for negative spin Hall studies
    pub fn py_ta() -> Self {
        Self {
            g_r: 4.5e19,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12,
        }
    }

    /// Create a CoFeB/Ta interface
    ///
    /// Perpendicular anisotropy with Ta underlayer
    pub fn cofeb_ta() -> Self {
        Self {
            g_r: 3.8e19,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12,
        }
    }

    /// Create a CoFeB/W interface
    ///
    /// Tungsten with giant negative spin Hall angle
    pub fn cofeb_w() -> Self {
        Self {
            g_r: 4.2e19,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12,
        }
    }

    /// Create a Py/W interface
    ///
    /// Strong spin-orbit coupling
    pub fn py_w() -> Self {
        Self {
            g_r: 5.2e19,
            g_i: 0.0,
            normal: Vector3::new(0.0, 1.0, 0.0),
            area: 1.0e-12,
        }
    }

    /// Builder method to set spin-mixing conductance
    pub fn with_g_r(mut self, g_r: f64) -> Self {
        self.g_r = g_r;
        self
    }

    /// Builder method to set imaginary part of spin-mixing conductance
    pub fn with_g_i(mut self, g_i: f64) -> Self {
        self.g_i = g_i;
        self
    }

    /// Builder method to set interface normal vector
    pub fn with_normal(mut self, normal: Vector3<f64>) -> Self {
        self.normal = normal.normalize();
        self
    }

    /// Builder method to set interface area
    pub fn with_area(mut self, area: f64) -> Self {
        self.area = area;
        self
    }
}

impl fmt::Display for SpinInterface {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SpinInterface(g_r={:.2e} Ω⁻¹m⁻², g_i={:.2e} Ω⁻¹m⁻², area={:.2e} m²)",
            self.g_r, self.g_i, self.area
        )
    }
}

impl super::traits::InterfaceMaterial for SpinInterface {
    fn spin_mixing_conductance(&self) -> f64 {
        self.g_r
    }

    fn interface_area(&self) -> f64 {
        self.area
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_interface() {
        let interface = SpinInterface::default();
        assert!(interface.g_r > 0.0);
        assert!((interface.normal.magnitude() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_yig_pt_interface() {
        let interface = SpinInterface::yig_pt();
        assert!((interface.g_r - 1.0e19).abs() < 1.0);
    }

    #[test]
    fn test_py_pt_interface() {
        let interface = SpinInterface::py_pt();
        assert!((interface.g_r - 5.0e19).abs() < 1.0);
    }

    #[test]
    fn test_cofeb_pt_interface() {
        let interface = SpinInterface::cofeb_pt();
        assert!((interface.g_r - 4.0e19).abs() < 1.0);
    }

    #[test]
    fn test_co_pt_interface() {
        let interface = SpinInterface::co_pt();
        assert!((interface.g_r - 6.0e19).abs() < 1.0);
    }

    #[test]
    fn test_fe_pt_interface() {
        let interface = SpinInterface::fe_pt();
        assert!((interface.g_r - 5.5e19).abs() < 1.0);
    }

    #[test]
    fn test_builder_methods() {
        let custom = SpinInterface::yig_pt()
            .with_g_r(2.0e19)
            .with_g_i(1.0e18)
            .with_area(5.0e-13);

        assert_eq!(custom.g_r, 2.0e19);
        assert_eq!(custom.g_i, 1.0e18);
        assert_eq!(custom.area, 5.0e-13);
    }

    #[test]
    fn test_with_normal() {
        let custom = SpinInterface::default().with_normal(Vector3::new(1.0, 1.0, 1.0));

        // Normal should be normalized
        assert!((custom.normal.magnitude() - 1.0).abs() < 1e-10);
    }
}
