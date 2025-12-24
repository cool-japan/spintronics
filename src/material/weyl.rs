//! Weyl Semimetal Materials
//!
//! Weyl semimetals are topological quantum materials where conduction and valence
//! bands touch at discrete points (Weyl nodes) in momentum space. These nodes act
//! as monopoles of Berry curvature, leading to exotic phenomena like the chiral
//! anomaly and giant anomalous Hall effect.
//!
//! ## Physics Background
//!
//! ### Weyl Nodes:
//! - Point degeneracies with linear dispersion E(k) = ℏv_F |k|
//! - Characterized by chirality (monopole charge): ±1
//! - Protected by topology (cannot be gapped without breaking symmetry)
//!
//! ### Key Phenomena:
//! 1. **Chiral Anomaly**: Parallel E and B fields cause charge pumping between nodes
//! 2. **Fermi Arc Surface States**: Open arcs connecting Weyl nodes of opposite chirality
//! 3. **Anomalous Hall Effect**: Large intrinsic Hall conductivity from Berry curvature
//! 4. **Negative Magnetoresistance**: From chiral anomaly
//!
//! ### Material Classes:
//! - **Type-I**: Weyl nodes with point-like Fermi surface
//! - **Type-II**: Tilted Weyl cones with finite Fermi surface
//! - **Magnetic**: Time-reversal broken, minimum 2 Weyl nodes
//! - **Non-magnetic**: Inversion broken, minimum 4 Weyl nodes
//!
//! ## Key References
//!
//! - S.-Y. Xu et al., "Discovery of a Weyl fermion semimetal and topological Fermi arcs",
//!   Science 349, 613 (2015) - TaAs
//! - B. Q. Lv et al., "Experimental Discovery of Weyl Semimetal TaAs",
//!   Phys. Rev. X 5, 031013 (2015)
//! - N. P. Armitage et al., "Weyl and Dirac semimetals in three-dimensional solids",
//!   Rev. Mod. Phys. 90, 015001 (2018)
//! - E. Liu et al., "Giant anomalous Hall effect in a ferromagnetic kagome-lattice semimetal",
//!   Nat. Phys. 14, 1125 (2018) - Co₃Sn₂S₂

use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::vector3::Vector3;

/// Weyl semimetal classification
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum WeylType {
    /// Type-I: Weyl nodes with point-like Fermi surface
    TypeI,
    /// Type-II: Tilted Weyl cones with finite Fermi surface
    TypeII,
}

/// Magnetic state
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum MagneticState {
    /// Non-magnetic (inversion symmetry broken)
    Nonmagnetic,
    /// Ferromagnetic (time-reversal symmetry broken)
    Ferromagnetic,
}

/// Weyl semimetal material properties
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct WeylSemimetal {
    /// Material name
    pub name: String,

    /// Weyl type (I or II)
    pub weyl_type: WeylType,

    /// Magnetic state
    pub magnetic_state: MagneticState,

    /// Number of Weyl nodes (pairs)
    pub num_weyl_nodes: usize,

    /// Fermi velocity at Weyl nodes \[m/s\]
    pub fermi_velocity: f64,

    /// Energy separation between Weyl nodes \[eV\]
    pub node_separation_energy: f64,

    /// Momentum separation between nodes [Å⁻¹]
    pub node_separation_momentum: f64,

    /// Anomalous Hall conductivity \[S/m\]
    pub anomalous_hall_conductivity: f64,

    /// Electrical conductivity \[S/m\]
    pub electrical_conductivity: f64,

    /// Carrier density [m⁻³]
    pub carrier_density: f64,
}

impl Default for WeylSemimetal {
    /// Default to TaAs parameters
    fn default() -> Self {
        Self::taas()
    }
}

impl WeylSemimetal {
    /// Create TaAs (Tantalum Arsenide)
    ///
    /// First experimentally confirmed Weyl semimetal
    /// Reference: S.-Y. Xu et al., Science 349, 613 (2015)
    pub fn taas() -> Self {
        Self {
            name: "TaAs".to_string(),
            weyl_type: WeylType::TypeI,
            magnetic_state: MagneticState::Nonmagnetic,
            num_weyl_nodes: 12,                 // 24 Weyl points total
            fermi_velocity: 6.0e5,              // m/s
            node_separation_energy: 0.02,       // eV
            node_separation_momentum: 0.03,     // Å⁻¹
            anomalous_hall_conductivity: 1.0e4, // S/m
            electrical_conductivity: 5.0e5,     // S/m
            carrier_density: 1.0e25,            // m⁻³
        }
    }

    /// Create Co₃Sn₂S₂ (Magnetic Weyl semimetal)
    ///
    /// Giant anomalous Hall effect in kagome lattice
    /// Reference: E. Liu et al., Nat. Phys. 14, 1125 (2018)
    pub fn co3sn2s2() -> Self {
        Self {
            name: "Co₃Sn₂S₂".to_string(),
            weyl_type: WeylType::TypeI,
            magnetic_state: MagneticState::Ferromagnetic,
            num_weyl_nodes: 2, // Minimum for magnetic WSM
            fermi_velocity: 4.0e5,
            node_separation_energy: 0.06, // eV
            node_separation_momentum: 0.1,
            anomalous_hall_conductivity: 1.5e6, // Giant AHE!
            electrical_conductivity: 2.0e6,
            carrier_density: 5.0e26,
        }
    }

    /// Create WTe₂ (Type-II Weyl semimetal)
    ///
    /// Tilted Weyl cones, large magnetoresistance
    /// Reference: A. A. Soluyanov et al., Nature 527, 495 (2015)
    pub fn wte2() -> Self {
        Self {
            name: "WTe₂".to_string(),
            weyl_type: WeylType::TypeII,
            magnetic_state: MagneticState::Nonmagnetic,
            num_weyl_nodes: 4, // 8 Weyl points
            fermi_velocity: 3.0e5,
            node_separation_energy: 0.05,
            node_separation_momentum: 0.02,
            anomalous_hall_conductivity: 5.0e4,
            electrical_conductivity: 1.0e6,
            carrier_density: 2.0e25,
        }
    }

    /// Create NbAs (Niobium Arsenide)
    ///
    /// Similar to TaAs, clean Weyl semimetal
    pub fn nbas() -> Self {
        Self {
            name: "NbAs".to_string(),
            weyl_type: WeylType::TypeI,
            magnetic_state: MagneticState::Nonmagnetic,
            num_weyl_nodes: 12,
            fermi_velocity: 5.5e5,
            node_separation_energy: 0.018,
            node_separation_momentum: 0.028,
            anomalous_hall_conductivity: 8.0e3,
            electrical_conductivity: 4.5e5,
            carrier_density: 8.0e24,
        }
    }

    /// Create TaP (Tantalum Phosphide)
    pub fn tap() -> Self {
        Self {
            name: "TaP".to_string(),
            weyl_type: WeylType::TypeI,
            magnetic_state: MagneticState::Nonmagnetic,
            num_weyl_nodes: 12,
            fermi_velocity: 6.5e5,
            node_separation_energy: 0.025,
            node_separation_momentum: 0.032,
            anomalous_hall_conductivity: 1.2e4,
            electrical_conductivity: 6.0e5,
            carrier_density: 1.2e25,
        }
    }

    /// Calculate anomalous Hall angle
    ///
    /// θ_AH = σ_AH / σ_xx
    ///
    /// # Returns
    /// Anomalous Hall angle (dimensionless)
    pub fn anomalous_hall_angle(&self) -> f64 {
        self.anomalous_hall_conductivity / self.electrical_conductivity
    }

    /// Calculate chiral anomaly magnetoconductivity
    ///
    /// Δσ = C · (E · B)
    ///
    /// # Arguments
    /// * `e_field` - Electric field [V/m]
    /// * `b_field` - Magnetic field \[T\]
    ///
    /// # Returns
    /// Additional conductivity from chiral anomaly \[S/m\]
    pub fn chiral_anomaly_conductivity(&self, e_field: Vector3<f64>, b_field: Vector3<f64>) -> f64 {
        // Chiral anomaly coefficient depends on node separation
        let e_dot_b = e_field.dot(&b_field);

        // Empirical scaling
        let coeff = self.node_separation_momentum * 1e10 * 1e-6; // Convert units

        coeff * e_dot_b.abs()
    }

    /// Calculate negative magnetoresistance from chiral anomaly
    ///
    /// MR = [ρ(B) - ρ(0)] / ρ(0)
    ///
    /// # Arguments
    /// * `b_field_magnitude` - Magnetic field strength \[T\]
    ///
    /// # Returns
    /// Magnetoresistance (negative for chiral anomaly)
    pub fn magnetoresistance(&self, b_field_magnitude: f64) -> f64 {
        // Negative MR from chiral anomaly
        // MR ∝ -B²
        let mr_coeff = 0.01 * self.node_separation_momentum; // Empirical
        -mr_coeff * b_field_magnitude * b_field_magnitude
    }

    /// Estimate Berry curvature magnitude near Weyl node
    ///
    /// Ω ~ k / |k|³ near node
    ///
    /// # Arguments
    /// * `k_distance` - Distance from Weyl node in k-space [Å⁻¹]
    ///
    /// # Returns
    /// Berry curvature magnitude \[Å²\]
    pub fn berry_curvature(&self, k_distance: f64) -> f64 {
        if k_distance < 1e-6 {
            // Singular at node
            1.0e6
        } else {
            // Ω ~ 1/k²
            1.0 / (k_distance * k_distance)
        }
    }

    /// Check if suitable for spintronics applications
    pub fn is_suitable_for_spintronics(&self) -> bool {
        // Strong anomalous Hall effect
        self.anomalous_hall_angle() > 0.01 ||
        // Magnetic Weyl semimetals are particularly useful
        self.magnetic_state == MagneticState::Ferromagnetic
    }

    /// Builder method to set number of Weyl nodes
    pub fn with_num_nodes(mut self, num: usize) -> Self {
        self.num_weyl_nodes = num;
        self
    }

    /// Builder method to set Fermi velocity
    pub fn with_fermi_velocity(mut self, vf: f64) -> Self {
        self.fermi_velocity = vf;
        self
    }
}

impl fmt::Display for WeylType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WeylType::TypeI => write!(f, "Type-I"),
            WeylType::TypeII => write!(f, "Type-II"),
        }
    }
}

impl fmt::Display for MagneticState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MagneticState::Nonmagnetic => write!(f, "Non-magnetic"),
            MagneticState::Ferromagnetic => write!(f, "Ferromagnetic"),
        }
    }
}

impl fmt::Display for WeylSemimetal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} [{}, {}]: {} Weyl nodes, v_F={:.2e} m/s, σ_AH={:.2e} S/m",
            self.name,
            self.weyl_type,
            self.magnetic_state,
            self.num_weyl_nodes * 2,
            self.fermi_velocity,
            self.anomalous_hall_conductivity
        )
    }
}

impl super::traits::SpinChargeConverter for WeylSemimetal {
    fn spin_hall_angle(&self) -> f64 {
        self.anomalous_hall_angle()
    }

    fn spin_hall_conductivity(&self) -> f64 {
        self.anomalous_hall_conductivity
    }
}

impl super::traits::TopologicalMaterial for WeylSemimetal {
    fn bulk_gap(&self) -> f64 {
        // Weyl semimetals are gapless at Weyl nodes
        // Return the energy separation between nodes as a characteristic scale
        self.node_separation_energy
    }

    fn surface_fermi_velocity(&self) -> f64 {
        self.fermi_velocity
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_taas() {
        let taas = WeylSemimetal::taas();
        assert_eq!(taas.name, "TaAs");
        assert_eq!(taas.weyl_type, WeylType::TypeI);
        assert_eq!(taas.magnetic_state, MagneticState::Nonmagnetic);
        assert_eq!(taas.num_weyl_nodes, 12);
    }

    #[test]
    fn test_co3sn2s2_magnetic() {
        let co3sn2s2 = WeylSemimetal::co3sn2s2();
        assert_eq!(co3sn2s2.magnetic_state, MagneticState::Ferromagnetic);
        // Magnetic WSM has minimum 2 nodes
        assert!(co3sn2s2.num_weyl_nodes >= 2);
        // Giant AHE
        assert!(co3sn2s2.anomalous_hall_conductivity > 1.0e6);
    }

    #[test]
    fn test_wte2_type_ii() {
        let wte2 = WeylSemimetal::wte2();
        assert_eq!(wte2.weyl_type, WeylType::TypeII);
        assert_eq!(wte2.magnetic_state, MagneticState::Nonmagnetic);
    }

    #[test]
    fn test_anomalous_hall_angle() {
        let taas = WeylSemimetal::taas();
        let theta_ah = taas.anomalous_hall_angle();

        // Should be small but measurable
        assert!(theta_ah > 0.0);
        assert!(theta_ah < 0.1);
    }

    #[test]
    fn test_co3sn2s2_giant_ahe() {
        let co3sn2s2 = WeylSemimetal::co3sn2s2();
        let theta_ah = co3sn2s2.anomalous_hall_angle();

        // Giant AHE should have large angle
        assert!(theta_ah > 0.1);
    }

    #[test]
    fn test_chiral_anomaly() {
        let taas = WeylSemimetal::taas();
        let e_field = Vector3::new(1000.0, 0.0, 0.0); // V/m
        let b_field = Vector3::new(1.0, 0.0, 0.0); // T (parallel)

        let delta_sigma = taas.chiral_anomaly_conductivity(e_field, b_field);

        // Should enhance conductivity
        assert!(delta_sigma > 0.0);
    }

    #[test]
    fn test_chiral_anomaly_perpendicular() {
        let taas = WeylSemimetal::taas();
        let e_field = Vector3::new(1000.0, 0.0, 0.0);
        let b_field = Vector3::new(0.0, 1.0, 0.0); // Perpendicular

        let delta_sigma = taas.chiral_anomaly_conductivity(e_field, b_field);

        // E ⊥ B → no chiral anomaly
        assert!(delta_sigma.abs() < 1e-10);
    }

    #[test]
    fn test_negative_magnetoresistance() {
        let taas = WeylSemimetal::taas();
        let b_field = 5.0; // T

        let mr = taas.magnetoresistance(b_field);

        // Should be negative (chiral anomaly effect)
        assert!(mr < 0.0);
    }

    #[test]
    fn test_magnetoresistance_quadratic() {
        let taas = WeylSemimetal::taas();
        let b1 = 1.0;
        let b2 = 2.0;

        let mr1 = taas.magnetoresistance(b1);
        let mr2 = taas.magnetoresistance(b2);

        // MR ∝ B²
        assert!((mr2 / mr1 - 4.0).abs() < 0.01);
    }

    #[test]
    fn test_berry_curvature_divergence() {
        let taas = WeylSemimetal::taas();

        let omega_far = taas.berry_curvature(0.1); // Far from node
        let omega_near = taas.berry_curvature(0.01); // Near node

        // Berry curvature increases near Weyl node
        assert!(omega_near > omega_far);
    }

    #[test]
    fn test_spintronics_suitability() {
        let co3sn2s2 = WeylSemimetal::co3sn2s2(); // Magnetic with giant AHE
        let taas = WeylSemimetal::taas(); // Nonmagnetic with small AHE

        assert!(co3sn2s2.is_suitable_for_spintronics());
        // TaAs has small AHE but still useful
        assert!(taas.anomalous_hall_conductivity > 0.0);
    }

    #[test]
    fn test_builder_pattern() {
        let custom = WeylSemimetal::taas()
            .with_num_nodes(24)
            .with_fermi_velocity(1.0e6);

        assert_eq!(custom.num_weyl_nodes, 24);
        assert_eq!(custom.fermi_velocity, 1.0e6);
    }

    #[test]
    fn test_material_family() {
        let taas = WeylSemimetal::taas();
        let nbas = WeylSemimetal::nbas();
        let tap = WeylSemimetal::tap();

        // All should be Type-I nonmagnetic
        assert_eq!(taas.weyl_type, WeylType::TypeI);
        assert_eq!(nbas.weyl_type, WeylType::TypeI);
        assert_eq!(tap.weyl_type, WeylType::TypeI);

        // All should have same number of nodes
        assert_eq!(taas.num_weyl_nodes, nbas.num_weyl_nodes);
        assert_eq!(taas.num_weyl_nodes, tap.num_weyl_nodes);
    }
}
