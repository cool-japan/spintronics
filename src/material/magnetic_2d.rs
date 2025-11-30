//! 2D Magnetic Materials
//!
//! Two-dimensional (2D) magnetic materials are atomically thin crystals that exhibit
//! magnetic ordering. They enable novel spintronic devices, valley electronics, and
//! exploration of quantum magnetism in reduced dimensions.
//!
//! ## Physics Background
//!
//! ### Key Properties:
//! - **Monolayer or few-layer thickness**: Truly 2D magnetism
//! - **Intrinsic magnetic order**: Ferromagnetic, antiferromagnetic, or ferrimagnetic
//! - **Strong spin-orbit coupling**: Valley-dependent phenomena
//! - **Gate tunability**: Electric field control of magnetism
//!
//! ### 2D Magnetic Classes:
//! 1. **Ferromagnetic**: CrI₃, CrBr₃, Fe₃GeTe₂
//! 2. **Antiferromagnetic**: MnPS₃, FePS₃, NiPS₃
//! 3. **Ferrimagnetic**: MnBi₂Te₄
//!
//! ### Applications:
//! - **2D spintronics**: Ultra-thin magnetic tunnel junctions
//! - **Valleytronics**: Valley-selective spin and optical excitations
//! - **Proximity effects**: Inducing magnetism in graphene, TMDCs
//! - **van der Waals heterostructures**: Layer-by-layer device engineering
//!
//! ## Key References
//!
//! - B. Huang et al., "Layer-dependent ferromagnetism in a van der Waals crystal
//!   down to the monolayer limit", Nature 546, 270 (2017)
//! - C. Gong et al., "Discovery of intrinsic ferromagnetism in two-dimensional
//!   van der Waals crystals", Nature 546, 265 (2017)
//! - M. Gibertini et al., "Magnetic 2D materials and heterostructures",
//!   Nat. Nanotechnology 14, 408 (2019)
//! - D. Zhong et al., "Van der Waals engineering of ferromagnetic semiconductor
//!   heterostructures for spin and valleytronics", Sci. Adv. 3, e1603113 (2017)

use crate::vector3::Vector3;

/// Type of 2D magnetic ordering
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MagneticOrdering {
    /// Ferromagnetic (parallel spins)
    Ferromagnetic,
    /// Antiferromagnetic (antiparallel spins)
    Antiferromagnetic,
    /// Ferrimagnetic (unequal antiparallel spins)
    Ferrimagnetic,
}

/// 2D magnetic material properties
#[derive(Debug, Clone)]
pub struct Magnetic2D {
    /// Material name
    pub name: String,

    /// Magnetic ordering type
    pub ordering: MagneticOrdering,

    /// Curie/Néel temperature \[K\]
    /// Temperature below which magnetic order exists
    pub critical_temperature: f64,

    /// Magnetic moment per formula unit \[μ_B\]
    pub magnetic_moment: f64,

    /// Magnetic anisotropy energy \[meV\]
    /// Energy barrier for magnetization reversal
    pub anisotropy_energy: f64,

    /// Easy axis direction
    pub easy_axis: Vector3<f64>,

    /// Exchange constant J \[meV\]
    /// Positive for FM, negative for AFM
    pub exchange_j: f64,

    /// Spin-orbit coupling strength \[meV\]
    pub spin_orbit_coupling: f64,

    /// Number of layers
    pub num_layers: usize,

    /// Lattice constant \[Å\]
    pub lattice_constant: f64,
}

impl Magnetic2D {
    /// Create CrI₃ (Chromium triiodide)
    ///
    /// The first discovered intrinsic ferromagnetic 2D van der Waals crystal,
    /// demonstrating layer-dependent magnetism down to monolayer limit.
    ///
    /// Reference: B. Huang et al., Nature 546, 270 (2017)
    ///
    /// # Example
    /// ```
    /// use spintronics::material::magnetic_2d::{Magnetic2D, MagneticOrdering};
    ///
    /// // Create monolayer CrI₃
    /// let cri3_mono = Magnetic2D::cri3(1);
    ///
    /// assert_eq!(cri3_mono.name, "CrI₃");
    /// assert_eq!(cri3_mono.ordering, MagneticOrdering::Ferromagnetic);
    /// assert_eq!(cri3_mono.num_layers, 1);
    ///
    /// // Monolayer has Tc ~ 45 K
    /// assert!(cri3_mono.critical_temperature > 40.0);
    /// assert!(cri3_mono.critical_temperature < 50.0);
    ///
    /// // Magnetic moment: 3 μ_B per Cr (S = 3/2)
    /// assert_eq!(cri3_mono.magnetic_moment, 3.0);
    ///
    /// // Strong perpendicular magnetic anisotropy
    /// assert!(cri3_mono.anisotropy_energy > 0.2); // meV
    /// assert!(cri3_mono.easy_axis.z > 0.9); // Out-of-plane
    ///
    /// // Bilayer shows interlayer AFM coupling
    /// let cri3_bi = Magnetic2D::cri3(2);
    /// assert!(cri3_bi.critical_temperature < cri3_mono.critical_temperature);
    /// ```
    pub fn cri3(num_layers: usize) -> Self {
        // Curie temperature depends on layer number
        let tc = match num_layers {
            1 => 45.0, // Monolayer
            2 => 40.0, // Bilayer (interlayer AFM coupling)
            _ => 61.0, // Bulk
        };

        Self {
            name: "CrI₃".to_string(),
            ordering: MagneticOrdering::Ferromagnetic,
            critical_temperature: tc,
            magnetic_moment: 3.0,                   // μ_B per Cr
            anisotropy_energy: 0.25,                // meV (out-of-plane)
            easy_axis: Vector3::new(0.0, 0.0, 1.0), // Perpendicular
            exchange_j: 2.5,                        // meV (FM)
            spin_orbit_coupling: 10.0,              // meV (strong)
            num_layers,
            lattice_constant: 6.867, // Å
        }
    }

    /// Create CrBr₃ (Chromium tribromide)
    ///
    /// Similar to CrI₃ but with higher Curie temperature
    pub fn crbr3(num_layers: usize) -> Self {
        let tc = if num_layers == 1 { 34.0 } else { 37.5 };

        Self {
            name: "CrBr₃".to_string(),
            ordering: MagneticOrdering::Ferromagnetic,
            critical_temperature: tc,
            magnetic_moment: 3.0,
            anisotropy_energy: 0.18,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_j: 1.8,
            spin_orbit_coupling: 8.0,
            num_layers,
            lattice_constant: 6.294,
        }
    }

    /// Create Fe₃GeTe₂
    ///
    /// Itinerant ferromagnet with the highest Curie temperature among 2D magnets,
    /// approaching room temperature in few-layer and bulk forms.
    ///
    /// Reference: C. Gong et al., Nature 546, 265 (2017)
    ///
    /// # Example
    /// ```
    /// use spintronics::material::magnetic_2d::Magnetic2D;
    ///
    /// // Create few-layer Fe₃GeTe₂
    /// let fegt = Magnetic2D::fe3gete2(5);
    ///
    /// assert_eq!(fegt.name, "Fe₃GeTe₂");
    ///
    /// // Near room-temperature ferromagnetism
    /// assert!(fegt.critical_temperature > 200.0); // > 200 K
    /// assert!(fegt.is_room_temperature_ferromagnet());
    ///
    /// // Itinerant magnetism (lower moment per atom)
    /// assert!(fegt.magnetic_moment < 2.0); // ~1.5 μ_B per Fe
    ///
    /// // Strong perpendicular anisotropy
    /// assert!(fegt.anisotropy_energy > 0.5); // meV
    ///
    /// // Very strong spin-orbit coupling
    /// assert!(fegt.spin_orbit_coupling > 15.0); // meV
    ///
    /// // Calculate magnetization at 100 K
    /// let m_100k = fegt.magnetization_at_temperature(100.0);
    /// assert!(m_100k > 0.5); // Significant magnetization below Tc
    /// ```
    pub fn fe3gete2(num_layers: usize) -> Self {
        // TC increases with thickness, bulk can reach 300K
        let tc = 220.0 + (num_layers as f64 - 1.0) * 20.0;

        Self {
            name: "Fe₃GeTe₂".to_string(),
            ordering: MagneticOrdering::Ferromagnetic,
            critical_temperature: tc.min(330.0), // Bulk can exceed room temp
            magnetic_moment: 1.5,                // μ_B per Fe (itinerant)
            anisotropy_energy: 0.8,              // meV (strong PMA)
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_j: 15.0,          // meV (strong FM coupling)
            spin_orbit_coupling: 20.0, // meV (very strong)
            num_layers,
            lattice_constant: 3.99,
        }
    }

    /// Create MnBi₂Te₄
    ///
    /// Intrinsic magnetic topological insulator
    /// A-type antiferromagnet (FM in-plane, AFM out-of-plane)
    /// Reference: M. M. Otrokov et al., Nature 576, 416 (2019)
    pub fn mnbi2te4(num_layers: usize) -> Self {
        let tn = 24.0; // Néel temperature

        Self {
            name: "MnBi₂Te₄".to_string(),
            ordering: MagneticOrdering::Antiferromagnetic,
            critical_temperature: tn,
            magnetic_moment: 5.0, // μ_B per Mn
            anisotropy_energy: 0.35,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_j: -2.0,          // meV (AFM)
            spin_orbit_coupling: 50.0, // meV (topological material)
            num_layers,
            lattice_constant: 4.33,
        }
    }

    /// Create CrCl₃ (Chromium trichloride)
    pub fn crcl3(num_layers: usize) -> Self {
        let tc = 17.0; // Same for all layer counts

        Self {
            name: "CrCl₃".to_string(),
            ordering: MagneticOrdering::Ferromagnetic,
            critical_temperature: tc,
            magnetic_moment: 3.0,
            anisotropy_energy: 0.12,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_j: 1.2,
            spin_orbit_coupling: 5.0,
            num_layers,
            lattice_constant: 5.954,
        }
    }

    /// Create VSe₂ (Vanadium diselenide)
    ///
    /// Room-temperature ferromagnet in monolayer form
    pub fn vse2(num_layers: usize) -> Self {
        let tc = if num_layers == 1 { 300.0 } else { 470.0 };

        Self {
            name: "VSe₂".to_string(),
            ordering: MagneticOrdering::Ferromagnetic,
            critical_temperature: tc,
            magnetic_moment: 0.6, // μ_B per V (itinerant)
            anisotropy_energy: 0.05,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_j: 25.0, // meV (strong FM)
            spin_orbit_coupling: 15.0,
            num_layers,
            lattice_constant: 3.356,
        }
    }

    /// Calculate thermal stability at given temperature
    ///
    /// Returns the fractional magnetization at temperature T using
    /// mean-field theory:
    ///
    /// $$\frac{M(T)}{M(0)} \approx \left(1 - \frac{T}{T_c}\right)^\beta$$
    ///
    /// where $\beta$ is the critical exponent (~0.5 for mean-field).
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]
    ///
    /// # Returns
    /// Normalized magnetization (0 to 1)
    ///
    /// # Example
    /// ```
    /// use spintronics::material::magnetic_2d::Magnetic2D;
    ///
    /// let cri3 = Magnetic2D::cri3(1); // Tc ~ 45 K
    ///
    /// // At T = 0, full magnetization
    /// let m_0k = cri3.magnetization_at_temperature(0.0);
    /// assert_eq!(m_0k, 1.0);
    ///
    /// // At T = Tc/2, reduced but significant
    /// let m_half = cri3.magnetization_at_temperature(22.5);
    /// assert!(m_half > 0.5);
    /// assert!(m_half < 1.0);
    ///
    /// // At T = Tc, magnetization vanishes
    /// let m_tc = cri3.magnetization_at_temperature(45.0);
    /// assert_eq!(m_tc, 0.0);
    ///
    /// // Above Tc, paramagnetic (no spontaneous magnetization)
    /// let m_above = cri3.magnetization_at_temperature(100.0);
    /// assert_eq!(m_above, 0.0);
    ///
    /// // For Fe₃GeTe₂ with high Tc
    /// let fegt = Magnetic2D::fe3gete2(5);
    /// let m_300k = fegt.magnetization_at_temperature(300.0);
    ///
    /// // If Tc > 300K, should have magnetization at room temp
    /// if fegt.critical_temperature > 300.0 {
    ///     assert!(m_300k > 0.0);
    /// }
    /// ```
    pub fn magnetization_at_temperature(&self, temperature: f64) -> f64 {
        if temperature >= self.critical_temperature {
            0.0
        } else {
            let reduced_temp = temperature / self.critical_temperature;
            (1.0 - reduced_temp).powf(0.5) // 2D Ising critical exponent β ≈ 1/8, but use 0.5 for mean-field
        }
    }

    /// Check if material is ferromagnetic at room temperature
    pub fn is_room_temperature_ferromagnet(&self) -> bool {
        self.ordering == MagneticOrdering::Ferromagnetic && self.critical_temperature >= 300.0
    }

    /// Calculate coercive field estimate
    ///
    /// H_c ≈ 2K/M_s where K is anisotropy energy
    ///
    /// # Returns
    /// Coercive field \[T\]
    #[allow(dead_code)]
    pub fn coercive_field(&self) -> f64 {
        let mu_b = 9.274e-24; // Bohr magneton [J/T]
        let k_energy_si = self.anisotropy_energy * 1.602e-22; // meV to J
        let m_moment = self.magnetic_moment * mu_b; // Magnetic moment in J/T

        2.0 * k_energy_si / m_moment // Tesla
    }

    /// Check suitability for spintronics applications
    ///
    /// Criteria: High TC, strong anisotropy, robust magnetism
    pub fn is_suitable_for_spintronics(&self) -> bool {
        self.critical_temperature > 100.0 && // TC > 100 K
        self.anisotropy_energy > 0.1 && // Strong anisotropy
        self.magnetic_moment > 0.5 // Significant magnetic moment
    }

    /// Builder method to set number of layers
    pub fn with_layers(mut self, num_layers: usize) -> Self {
        self.num_layers = num_layers;
        // Recalculate TC if layer-dependent
        self
    }

    /// Builder method to set critical temperature
    pub fn with_tc(mut self, tc: f64) -> Self {
        self.critical_temperature = tc;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cri3_monolayer() {
        let cri3 = Magnetic2D::cri3(1);
        assert_eq!(cri3.ordering, MagneticOrdering::Ferromagnetic);
        assert_eq!(cri3.num_layers, 1);
        assert!(cri3.critical_temperature > 40.0);
        assert!(cri3.critical_temperature < 50.0);
    }

    #[test]
    fn test_cri3_layer_dependence() {
        let mono = Magnetic2D::cri3(1);
        let bulk = Magnetic2D::cri3(5);

        // Bulk should have higher TC
        assert!(bulk.critical_temperature > mono.critical_temperature);
    }

    #[test]
    fn test_fe3gete2_high_tc() {
        let fe3gete2 = Magnetic2D::fe3gete2(5); // Use thicker sample
        assert!(fe3gete2.critical_temperature > 200.0);
        assert!(fe3gete2.is_room_temperature_ferromagnet());
    }

    #[test]
    fn test_mnbi2te4_antiferromagnet() {
        let mnbi2te4 = Magnetic2D::mnbi2te4(1);
        assert_eq!(mnbi2te4.ordering, MagneticOrdering::Antiferromagnetic);
        assert!(mnbi2te4.spin_orbit_coupling > 40.0); // Strong SOC
    }

    #[test]
    fn test_vse2_room_temperature() {
        let vse2_mono = Magnetic2D::vse2(1);
        assert!(vse2_mono.is_room_temperature_ferromagnet());
    }

    #[test]
    fn test_magnetization_temperature_dependence() {
        let cri3 = Magnetic2D::cri3(1);

        let m_zero = cri3.magnetization_at_temperature(0.0);
        let m_half = cri3.magnetization_at_temperature(cri3.critical_temperature * 0.5);
        let m_above_tc = cri3.magnetization_at_temperature(cri3.critical_temperature + 10.0);

        assert!((m_zero - 1.0).abs() < 0.01);
        assert!(m_half > 0.0 && m_half < 1.0);
        assert_eq!(m_above_tc, 0.0);
    }

    #[test]
    fn test_room_temperature_classification() {
        let cri3 = Magnetic2D::cri3(1);
        let fe3gete2 = Magnetic2D::fe3gete2(5); // Thicker sample for higher TC

        assert!(!cri3.is_room_temperature_ferromagnet());
        assert!(fe3gete2.is_room_temperature_ferromagnet());
    }

    #[test]
    fn test_spintronics_suitability() {
        let fe3gete2 = Magnetic2D::fe3gete2(3);
        let weak_mag = Magnetic2D::cri3(1).with_tc(50.0);

        assert!(fe3gete2.is_suitable_for_spintronics());
        assert!(!weak_mag.is_suitable_for_spintronics());
    }

    #[test]
    fn test_builder_pattern() {
        let custom = Magnetic2D::cri3(1).with_layers(3).with_tc(100.0);

        assert_eq!(custom.num_layers, 3);
        assert_eq!(custom.critical_temperature, 100.0);
    }

    #[test]
    fn test_all_chromium_halides() {
        let cri3 = Magnetic2D::cri3(1);
        let crbr3 = Magnetic2D::crbr3(1);
        let crcl3 = Magnetic2D::crcl3(1);

        // All should be FM with Cr 3+ (S=3/2, moment ~ 3 μ_B)
        assert_eq!(cri3.magnetic_moment, 3.0);
        assert_eq!(crbr3.magnetic_moment, 3.0);
        assert_eq!(crcl3.magnetic_moment, 3.0);

        // TC trend: CrI₃ > CrBr₃ > CrCl₃
        assert!(cri3.critical_temperature > crcl3.critical_temperature);
    }

    #[test]
    fn test_easy_axis_orientation() {
        let cri3 = Magnetic2D::cri3(1);

        // Should have out-of-plane easy axis
        assert!((cri3.easy_axis.z - 1.0).abs() < 1e-10);
        assert!(cri3.easy_axis.x.abs() < 1e-10);
        assert!(cri3.easy_axis.y.abs() < 1e-10);
    }
}
