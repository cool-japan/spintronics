//! Trait hierarchy for magnetic materials
//!
//! This module defines common traits that magnetic materials can implement,
//! enabling polymorphic behavior and consistent interfaces across different
//! material types.

use crate::vector3::Vector3;

/// Core trait for materials with magnetic properties
///
/// Implemented by ferromagnets, antiferromagnets, and 2D magnetic materials.
pub trait MagneticMaterial {
    /// Returns the saturation magnetization [A/m]
    fn saturation_magnetization(&self) -> f64;

    /// Returns the Gilbert damping constant (dimensionless)
    fn damping(&self) -> f64;

    /// Returns the exchange stiffness [J/m]
    fn exchange_stiffness(&self) -> f64;

    /// Returns the anisotropy constant [J/m³]
    fn anisotropy(&self) -> f64;

    /// Returns the easy axis direction (normalized)
    fn easy_axis(&self) -> Vector3<f64>;

    /// Check if the material is a soft magnet (low anisotropy)
    fn is_soft(&self) -> bool {
        self.anisotropy().abs() < 1e3
    }

    /// Check if the material is a hard magnet (high anisotropy)
    fn is_hard(&self) -> bool {
        self.anisotropy().abs() > 1e5
    }
}

/// Trait for materials that can perform spin-charge conversion
///
/// Implemented by heavy metals, topological insulators, and Weyl semimetals.
pub trait SpinChargeConverter {
    /// Returns the spin Hall angle (dimensionless)
    ///
    /// θ_SH = σ_SH / σ_charge
    fn spin_hall_angle(&self) -> f64;

    /// Returns the spin Hall conductivity [S/m]
    fn spin_hall_conductivity(&self) -> f64;

    /// Returns the conversion efficiency metric
    fn conversion_efficiency(&self) -> f64 {
        self.spin_hall_angle().abs()
    }
}

/// Trait for materials with temperature-dependent properties
pub trait TemperatureDependent {
    /// Returns the Curie (or Néel) temperature \[K\]
    fn critical_temperature(&self) -> f64;

    /// Calculate magnetization at given temperature [A/m]
    fn magnetization_at_temperature(&self, temperature: f64) -> f64;

    /// Check if the material is magnetically ordered at given temperature
    fn is_ordered_at(&self, temperature: f64) -> bool {
        temperature < self.critical_temperature()
    }

    /// Calculate the reduced temperature T/T_c
    fn reduced_temperature(&self, temperature: f64) -> f64 {
        temperature / self.critical_temperature()
    }
}

/// Trait for topological materials with protected surface states
pub trait TopologicalMaterial {
    /// Returns the bulk bandgap \[eV\]
    fn bulk_gap(&self) -> f64;

    /// Returns the Fermi velocity of surface/edge states [m/s]
    fn surface_fermi_velocity(&self) -> f64;

    /// Check if the material has well-defined surface states
    fn has_surface_states(&self) -> bool {
        self.bulk_gap() > 0.0 && self.surface_fermi_velocity() > 0.0
    }
}

/// Trait for materials in FM/NM bilayer structures
pub trait InterfaceMaterial {
    /// Returns the spin-mixing conductance \[Ω⁻¹ m⁻²\]
    fn spin_mixing_conductance(&self) -> f64;

    /// Returns the interface area \[m²\]
    fn interface_area(&self) -> f64;

    /// Calculate the total spin conductance [Ω⁻¹]
    fn total_spin_conductance(&self) -> f64 {
        self.spin_mixing_conductance() * self.interface_area()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MockFerromagnet {
        ms: f64,
        alpha: f64,
        k: f64,
        a: f64,
    }

    impl MagneticMaterial for MockFerromagnet {
        fn saturation_magnetization(&self) -> f64 {
            self.ms
        }

        fn damping(&self) -> f64 {
            self.alpha
        }

        fn exchange_stiffness(&self) -> f64 {
            self.a
        }

        fn anisotropy(&self) -> f64 {
            self.k
        }

        fn easy_axis(&self) -> Vector3<f64> {
            Vector3::new(0.0, 0.0, 1.0)
        }
    }

    #[test]
    fn test_magnetic_material_trait() {
        let soft = MockFerromagnet {
            ms: 8e5,
            alpha: 0.01,
            k: 100.0,
            a: 1e-11,
        };
        assert!(soft.is_soft());
        assert!(!soft.is_hard());

        let hard = MockFerromagnet {
            ms: 1e6,
            alpha: 0.01,
            k: 5e5,
            a: 2e-11,
        };
        assert!(!hard.is_soft());
        assert!(hard.is_hard());
    }

    struct MockTemperatureMaterial {
        tc: f64,
        ms0: f64,
    }

    impl TemperatureDependent for MockTemperatureMaterial {
        fn critical_temperature(&self) -> f64 {
            self.tc
        }

        fn magnetization_at_temperature(&self, temperature: f64) -> f64 {
            if temperature >= self.tc {
                0.0
            } else {
                self.ms0 * (1.0 - (temperature / self.tc).powf(0.5))
            }
        }
    }

    #[test]
    fn test_temperature_dependent_trait() {
        let material = MockTemperatureMaterial {
            tc: 500.0,
            ms0: 1e6,
        };

        assert!(material.is_ordered_at(300.0));
        assert!(!material.is_ordered_at(600.0));
        assert!((material.reduced_temperature(250.0) - 0.5).abs() < 1e-10);
    }
}
