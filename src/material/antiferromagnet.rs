//! Antiferromagnetic Materials
//!
//! Antiferromagnets (AFMs) are materials with antiparallel spin ordering resulting
//! in zero net magnetization. Despite lacking macroscopic magnetization, AFMs exhibit
//! rich spintronic phenomena and are promising for ultrafast, THz-frequency devices.
//!
//! ## Physics Background
//!
//! ### Antiferromagnetic Order:
//! - Sublattice magnetizations: **M_A** = -**M_B**
//! - Néel vector: **L** = **M_A** - **M_B** (AFM order parameter)
//! - No net magnetization: **M** = **M_A** + **M_B** = 0
//!
//! ### Key Advantages:
//! 1. **THz dynamics**: Natural resonance frequencies ~0.1-10 THz
//! 2. **No stray fields**: Ideal for high-density integration
//! 3. **Robust to magnetic fields**: Stable in large fields
//! 4. **Ultrafast switching**: Inertia-free dynamics
//!
//! ### Spintronic Phenomena:
//! - **Anisotropic Magnetoresistance (AMR)**: Angle-dependent resistance
//! - **Spin Hall Magnetoresistance**: Interface effect with heavy metals
//! - **Néel Spin-Orbit Torque**: Current-induced switching
//! - **THz Emission**: From spin dynamics
//!
//! ## Key References
//!
//! - T. Jungwirth et al., "Antiferromagnetic spintronics",
//!   Nat. Nanotechnology 11, 231 (2016)
//! - V. Baltz et al., "Antiferromagnetic spintronics",
//!   Rev. Mod. Phys. 90, 015005 (2018)
//! - P. Wadley et al., "Electrical switching of an antiferromagnet",
//!   Science 351, 587 (2016) - CuMnAs
//! - S. Y. Bodnar et al., "Writing and reading antiferromagnetic Mn₂Au by
//!   Néel spin-orbit torques", Nat. Commun. 9, 348 (2018)

use crate::vector3::Vector3;

/// Antiferromagnetic crystal structure
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AfmStructure {
    /// Type-I: Simple antiferromagnetic (e.g., MnO, NiO)
    TypeI,
    /// Type-II: With canted moments (weak ferromagnetism)
    TypeII,
    /// Collinear: Spins along single axis
    Collinear,
    /// Non-collinear: Complex spin structure
    NonCollinear,
}

/// Antiferromagnetic material properties
#[derive(Debug, Clone)]
pub struct Antiferromagnet {
    /// Material name
    pub name: String,

    /// AFM structure type
    pub structure: AfmStructure,

    /// Néel temperature \[K\]
    /// Temperature below which AFM order exists
    pub neel_temperature: f64,

    /// Sublattice magnetization \[A/m\]
    pub sublattice_magnetization: f64,

    /// Exchange field \[T\]
    /// Effective field from exchange interaction
    pub exchange_field: f64,

    /// Anisotropy field \[T\]
    /// Determines easy axis
    pub anisotropy_field: f64,

    /// Easy axis direction
    pub easy_axis: Vector3<f64>,

    /// AFM resonance frequency \[THz\]
    pub resonance_frequency: f64,

    /// Electrical conductivity \[S/m\]
    pub electrical_conductivity: f64,

    /// Spin Hall angle (for interfacial effects)
    pub spin_hall_angle: f64,
}

impl Antiferromagnet {
    /// Create NiO (Nickel Oxide)
    ///
    /// Prototypical insulating antiferromagnet
    /// Reference: S. A. Wolf et al., Science 294, 1488 (2001)
    pub fn nio() -> Self {
        Self {
            name: "NiO".to_string(),
            structure: AfmStructure::TypeII,
            neel_temperature: 523.0,                            // K
            sublattice_magnetization: 1.77e5,                   // A/m
            exchange_field: 660.0,                              // T (huge!)
            anisotropy_field: 0.02,                             // T
            easy_axis: Vector3::new(1.0, 1.0, 1.0).normalize(), // [111]
            resonance_frequency: 1.0,                           // THz
            electrical_conductivity: 1.0e-10,                   // Insulator
            spin_hall_angle: 0.0,                               // Insulator
        }
    }

    /// Create Mn₂Au (Metallic antiferromagnet)
    ///
    /// Room-temperature electrical switching
    /// Reference: S. Y. Bodnar et al., Nat. Commun. 9, 348 (2018)
    pub fn mn2au() -> Self {
        Self {
            name: "Mn₂Au".to_string(),
            structure: AfmStructure::Collinear,
            neel_temperature: 1500.0, // K (well above RT)
            sublattice_magnetization: 2.0e6,
            exchange_field: 800.0,
            anisotropy_field: 1.0,
            easy_axis: Vector3::new(0.0, 0.0, 1.0), // c-axis
            resonance_frequency: 0.9,
            electrical_conductivity: 1.0e6, // Metallic
            spin_hall_angle: 0.02,          // Weak SOC
        }
    }

    /// Create CuMnAs
    ///
    /// First demonstrated electrical AFM switching
    /// Reference: P. Wadley et al., Science 351, 587 (2016)
    pub fn cumnas() -> Self {
        Self {
            name: "CuMnAs".to_string(),
            structure: AfmStructure::Collinear,
            neel_temperature: 480.0,
            sublattice_magnetization: 1.5e6,
            exchange_field: 700.0,
            anisotropy_field: 0.5,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            resonance_frequency: 0.8,
            electrical_conductivity: 5.0e5,
            spin_hall_angle: 0.015,
        }
    }

    /// Create IrMn₃ (Iridium Manganese)
    ///
    /// Common pinning layer in spin valves
    pub fn irmn3() -> Self {
        Self {
            name: "IrMn₃".to_string(),
            structure: AfmStructure::TypeI,
            neel_temperature: 960.0,
            sublattice_magnetization: 1.2e6,
            exchange_field: 500.0,
            anisotropy_field: 10.0, // Strong for pinning
            easy_axis: Vector3::new(1.0, 0.0, 0.0),
            resonance_frequency: 0.5,
            electrical_conductivity: 2.0e6,
            spin_hall_angle: 0.05, // Strong SOC from Ir
        }
    }

    /// Create α-Fe₂O₃ (Hematite)
    ///
    /// Weak ferromagnetism, common mineral
    pub fn fe2o3() -> Self {
        Self {
            name: "α-Fe₂O₃".to_string(),
            structure: AfmStructure::TypeII, // Canted
            neel_temperature: 948.0,         // Morin transition at 260 K
            sublattice_magnetization: 2.5e6,
            exchange_field: 600.0,
            anisotropy_field: 0.05,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            resonance_frequency: 0.6,
            electrical_conductivity: 1.0e-6, // Semiconductor
            spin_hall_angle: 0.0,
        }
    }

    /// Create MnF₂ (Manganese Fluoride)
    ///
    /// Model Ising antiferromagnet
    pub fn mnf2() -> Self {
        Self {
            name: "MnF₂".to_string(),
            structure: AfmStructure::TypeI,
            neel_temperature: 67.0, // Low TN
            sublattice_magnetization: 5.0e5,
            exchange_field: 55.0,
            anisotropy_field: 8.8, // Strong uniaxial
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            resonance_frequency: 0.25,
            electrical_conductivity: 1.0e-12, // Insulator
            spin_hall_angle: 0.0,
        }
    }

    /// Calculate AFM resonance frequency
    ///
    /// ω_R = γ √(2H_E H_A)
    ///
    /// # Returns
    /// Resonance frequency \[THz\]
    pub fn calculate_resonance_frequency(&self) -> f64 {
        let gamma = 2.8e10; // Hz/T (gyromagnetic ratio)
        let omega = gamma * (2.0 * self.exchange_field * self.anisotropy_field).sqrt();
        omega * 1e-12 // Convert to THz
    }

    /// Calculate exchange stiffness constant
    ///
    /// A = M_s H_E a / 2
    ///
    /// # Returns
    /// Exchange stiffness \[J/m\]
    pub fn exchange_stiffness(&self) -> f64 {
        let mu_0 = 4.0 * std::f64::consts::PI * 1e-7;
        let a = 3.0e-10; // Typical lattice constant \[m\]
        let h_e_si = self.exchange_field * mu_0;

        self.sublattice_magnetization * h_e_si * a / 2.0
    }

    /// Check if suitable for room-temperature applications
    pub fn is_room_temperature_stable(&self) -> bool {
        self.neel_temperature > 300.0
    }

    /// Check if suitable for electrical switching
    pub fn is_electrically_switchable(&self) -> bool {
        // Metallic AFMs with sufficient SOC
        self.electrical_conductivity > 1.0e4 && self.spin_hall_angle > 0.01
    }

    /// Check if suitable for THz spintronics
    pub fn is_thz_active(&self) -> bool {
        self.resonance_frequency > 0.1 && // Above 100 GHz
        self.neel_temperature > 77.0 // Stable at liquid N2 temp
    }

    /// Calculate Néel spin-orbit torque efficiency
    ///
    /// For switching applications
    ///
    /// # Returns
    /// Effective torque efficiency
    pub fn neel_sot_efficiency(&self) -> f64 {
        // Proportional to spin Hall angle and structure
        let structure_factor = match self.structure {
            AfmStructure::Collinear => 1.0,
            AfmStructure::NonCollinear => 0.8,
            _ => 0.5,
        };

        self.spin_hall_angle * structure_factor
    }

    /// Builder method to set Néel temperature
    pub fn with_neel_temperature(mut self, tn: f64) -> Self {
        self.neel_temperature = tn;
        self
    }

    /// Builder method to set easy axis
    pub fn with_easy_axis(mut self, axis: Vector3<f64>) -> Self {
        self.easy_axis = axis.normalize();
        self
    }

    /// Builder method to set sublattice magnetization
    pub fn with_sublattice_magnetization(mut self, ms: f64) -> Self {
        self.sublattice_magnetization = ms;
        self
    }

    /// Builder method to set exchange field
    pub fn with_exchange_field(mut self, h_ex: f64) -> Self {
        self.exchange_field = h_ex;
        self
    }

    /// Builder method to set anisotropy field
    pub fn with_anisotropy_field(mut self, h_k: f64) -> Self {
        self.anisotropy_field = h_k;
        self
    }

    /// Builder method to set resonance frequency
    pub fn with_resonance_frequency(mut self, freq: f64) -> Self {
        self.resonance_frequency = freq;
        self
    }

    /// Builder method to set spin Hall angle
    pub fn with_spin_hall_angle(mut self, theta_sh: f64) -> Self {
        self.spin_hall_angle = theta_sh;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nio() {
        let nio = Antiferromagnet::nio();
        assert_eq!(nio.structure, AfmStructure::TypeII);
        assert!(nio.neel_temperature > 500.0);
        assert!(nio.electrical_conductivity < 1.0e-8); // Insulator
    }

    #[test]
    fn test_mn2au_metallic() {
        let mn2au = Antiferromagnet::mn2au();
        assert!(mn2au.is_room_temperature_stable());
        assert!(mn2au.is_electrically_switchable());
        assert!(mn2au.electrical_conductivity > 1.0e5);
    }

    #[test]
    fn test_cumnas_switching() {
        let cumnas = Antiferromagnet::cumnas();
        assert!(cumnas.is_electrically_switchable());
        assert!(cumnas.is_room_temperature_stable());
    }

    #[test]
    fn test_irmn3_pinning() {
        let irmn3 = Antiferromagnet::irmn3();
        // Strong anisotropy for pinning
        assert!(irmn3.anisotropy_field > 5.0);
        // High Néel temperature
        assert!(irmn3.neel_temperature > 900.0);
    }

    #[test]
    fn test_mnf2_low_tn() {
        let mnf2 = Antiferromagnet::mnf2();
        assert!(!mnf2.is_room_temperature_stable());
        assert!(mnf2.neel_temperature < 100.0);
        // Strong anisotropy
        assert!(mnf2.anisotropy_field > 5.0);
    }

    #[test]
    fn test_resonance_frequency() {
        let nio = Antiferromagnet::nio();
        let f_r = nio.calculate_resonance_frequency();

        // Should be in THz range
        assert!(f_r > 0.1);
        assert!(f_r < 10.0);
    }

    #[test]
    fn test_exchange_stiffness() {
        let nio = Antiferromagnet::nio();
        let a_ex = nio.exchange_stiffness();

        // Should be positive and reasonable
        assert!(a_ex > 0.0);
        assert!(a_ex.is_finite());
    }

    #[test]
    fn test_thz_activity() {
        let nio = Antiferromagnet::nio();
        let mnf2 = Antiferromagnet::mnf2();

        assert!(nio.is_thz_active());
        assert!(!mnf2.is_thz_active()); // Too low TN for RT operation
    }

    #[test]
    fn test_electrical_switchability() {
        let mn2au = Antiferromagnet::mn2au();
        let nio = Antiferromagnet::nio();

        assert!(mn2au.is_electrically_switchable());
        assert!(!nio.is_electrically_switchable()); // Insulator
    }

    #[test]
    fn test_neel_sot_efficiency() {
        let mn2au = Antiferromagnet::mn2au();
        let nio = Antiferromagnet::nio();

        let eff_mn2au = mn2au.neel_sot_efficiency();
        let eff_nio = nio.neel_sot_efficiency();

        // Metallic AFM should have higher efficiency
        assert!(eff_mn2au > eff_nio);
        assert!(eff_mn2au > 0.01);
    }

    #[test]
    fn test_room_temperature_stability() {
        let mn2au = Antiferromagnet::mn2au();
        let cumnas = Antiferromagnet::cumnas();
        let mnf2 = Antiferromagnet::mnf2();

        assert!(mn2au.is_room_temperature_stable());
        assert!(cumnas.is_room_temperature_stable());
        assert!(!mnf2.is_room_temperature_stable());
    }

    #[test]
    fn test_builder_pattern() {
        let custom = Antiferromagnet::nio()
            .with_neel_temperature(600.0)
            .with_easy_axis(Vector3::new(1.0, 0.0, 0.0));

        assert_eq!(custom.neel_temperature, 600.0);
        assert!((custom.easy_axis.x - 1.0).abs() < 1e-10);
    }
}
