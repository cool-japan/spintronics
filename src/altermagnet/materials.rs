//! Altermagnet material database
//!
//! Altermagnets are a newly classified magnetic phase (2022-2024) distinct from
//! both ferromagnets and antiferromagnets. Key properties include:
//!
//! - **Zero net magnetization** (like antiferromagnets)
//! - **Non-relativistic spin splitting** in k-space (unlike AFMs)
//! - **d-wave, g-wave, or i-wave symmetry** of the spin splitting
//! - **No spin-orbit coupling required** for spin-split bands
//!
//! ## Materials
//!
//! | Material | Symmetry | T_N (K) | Spin splitting (eV) |
//! |----------|----------|---------|---------------------|
//! | RuO2     | d-wave   | ~300    | ~1.4                |
//! | CrSb     | g-wave   | ~700    | ~1.2                |
//! | MnTe     | d-wave   | ~307    | ~0.9                |
//! | Fe2O3    | d-wave   | ~960    | ~0.6                |
//!
//! ## References
//!
//! - L. Smejkal et al., "Emerging Research Landscape of Altermagnetism",
//!   Phys. Rev. X 12, 040501 (2022)
//! - L. Smejkal et al., "Beyond Conventional Ferromagnetism and Antiferromagnetism",
//!   Phys. Rev. X 12, 031042 (2022)

use std::f64::consts::PI;
use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::constants::KB;
use crate::error::{Error, Result};

/// Symmetry classification of the altermagnetic spin splitting in k-space.
///
/// The symmetry determines the angular dependence of the spin splitting:
/// - d-wave: 2-fold rotational symmetry with nodes at 45 degrees
/// - g-wave: 4-fold rotational symmetry with nodes at 22.5 degrees
/// - i-wave: 6-fold rotational symmetry with nodes at 15 degrees
///
/// The spin splitting Δ(φ) follows:
/// - d-wave: Δ₀ cos(2φ)
/// - g-wave: Δ₀ cos(4φ)
/// - i-wave: Δ₀ cos(6φ)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum AltermagneticSymmetry {
    /// d-wave symmetry: 2-fold rotational symmetry of splitting.
    /// Nodes at 45, 135, 225, 315 degrees in the Brillouin zone.
    /// Most common symmetry type (e.g., RuO2, MnTe, Fe2O3).
    DWave,
    /// g-wave symmetry: 4-fold rotational symmetry of splitting.
    /// Nodes at 22.5, 67.5, 112.5, ... degrees.
    /// Example: CrSb.
    GWave,
    /// i-wave symmetry: 6-fold rotational symmetry of splitting.
    /// Nodes at 15, 45, 75, ... degrees.
    /// Theoretically predicted, fewer experimental examples.
    IWave,
}

impl AltermagneticSymmetry {
    /// Return the harmonic order of the spin splitting.
    ///
    /// d-wave -> 2, g-wave -> 4, i-wave -> 6
    pub fn harmonic_order(&self) -> u32 {
        match self {
            AltermagneticSymmetry::DWave => 2,
            AltermagneticSymmetry::GWave => 4,
            AltermagneticSymmetry::IWave => 6,
        }
    }

    /// Return the angular positions (in radians) where the spin splitting
    /// has nodes (zero crossings) within [0, 2π).
    ///
    /// For a cos(n*φ) dependence, nodes occur at φ = π/(2n) + k*π/n.
    pub fn node_angles(&self) -> Vec<f64> {
        let n = self.harmonic_order() as f64;
        let num_nodes = (2.0 * n) as usize;
        (0..num_nodes)
            .map(|k| PI / (2.0 * n) + (k as f64) * PI / n)
            .collect()
    }

    /// Compute the angular factor of spin splitting at angle φ (radians).
    ///
    /// Returns cos(n*φ) where n is the harmonic order.
    /// The full splitting is Δ(φ) = Δ₀ * angular_factor(φ).
    pub fn angular_factor(&self, phi: f64) -> f64 {
        let n = self.harmonic_order() as f64;
        (n * phi).cos()
    }

    /// Number of nodes in the full 2π angular range.
    pub fn num_nodes(&self) -> usize {
        2 * self.harmonic_order() as usize
    }
}

impl fmt::Display for AltermagneticSymmetry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AltermagneticSymmetry::DWave => write!(f, "d-wave"),
            AltermagneticSymmetry::GWave => write!(f, "g-wave"),
            AltermagneticSymmetry::IWave => write!(f, "i-wave"),
        }
    }
}

/// Altermagnet material with physical parameters.
///
/// Represents an altermagnetic material with zero net magnetization but
/// non-relativistic spin splitting in k-space. The spin splitting has
/// characteristic angular symmetry (d-wave, g-wave, or i-wave).
///
/// # Example
///
/// ```rust
/// use spintronics::altermagnet::{Altermagnet, AltermagneticSymmetry};
///
/// let ruo2 = Altermagnet::ruo2();
/// assert_eq!(ruo2.symmetry, AltermagneticSymmetry::DWave);
/// assert!(ruo2.spin_splitting > 1.0); // Large spin splitting in eV
///
/// // Verify zero net magnetization
/// assert!((ruo2.net_magnetization()).abs() < 1e-10);
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Altermagnet {
    /// Material name
    pub name: &'static str,
    /// Neel temperature \[K\]
    pub neel_temperature: f64,
    /// Maximum spin splitting energy \[eV\]
    pub spin_splitting: f64,
    /// Symmetry type of the spin splitting
    pub symmetry: AltermagneticSymmetry,
    /// Crystal structure description
    pub crystal_structure: &'static str,
    /// Sublattice magnetization magnitude \[A/m\]
    ///
    /// Each sublattice has this magnetization magnitude, but they point
    /// in opposite directions yielding zero net magnetization.
    pub sublattice_magnetization: f64,
    /// Lattice constant \[m\]
    pub lattice_constant: f64,
}

impl Altermagnet {
    /// Create RuO2 (ruthenium dioxide) - a prototypical d-wave altermagnet.
    ///
    /// Rutile crystal structure. The altermagnetic classification of RuO2
    /// has been debated, with the Neel temperature reported around 300 K
    /// in some studies, while others question the magnetic ordering.
    ///
    /// Reference: L. Smejkal et al., Phys. Rev. X 12, 040501 (2022)
    pub fn ruo2() -> Self {
        Self {
            name: "RuO2",
            neel_temperature: 300.0,
            spin_splitting: 1.4,
            symmetry: AltermagneticSymmetry::DWave,
            crystal_structure: "Rutile (P4_2/mnm)",
            sublattice_magnetization: 1.0e5,
            lattice_constant: 4.49e-10,
        }
    }

    /// Create CrSb (chromium antimonide) - a g-wave altermagnet.
    ///
    /// NiAs-type crystal structure with high Neel temperature and large
    /// spin splitting. The g-wave symmetry is experimentally confirmed
    /// via ARPES measurements.
    ///
    /// Reference: J. Krempasky et al., Nature 626, 517 (2024)
    pub fn crsb() -> Self {
        Self {
            name: "CrSb",
            neel_temperature: 700.0,
            spin_splitting: 1.2,
            symmetry: AltermagneticSymmetry::GWave,
            crystal_structure: "NiAs-type (P6_3/mmc)",
            sublattice_magnetization: 2.4e5,
            lattice_constant: 4.12e-10,
        }
    }

    /// Create MnTe (manganese telluride) - a d-wave altermagnet.
    ///
    /// NiAs-type crystal structure. First material where altermagnetic
    /// spin splitting was directly observed via ARPES.
    ///
    /// Reference: J. Krempasky et al., Nature 626, 517 (2024)
    pub fn mnte() -> Self {
        Self {
            name: "MnTe",
            neel_temperature: 307.0,
            spin_splitting: 0.9,
            symmetry: AltermagneticSymmetry::DWave,
            crystal_structure: "NiAs-type (P6_3/mmc)",
            sublattice_magnetization: 3.5e5,
            lattice_constant: 4.14e-10,
        }
    }

    /// Create Fe2O3 (hematite) - a d-wave altermagnet candidate.
    ///
    /// Corundum crystal structure with very high Neel temperature.
    /// Below the Morin transition (263 K), the spins are aligned along
    /// the c-axis; above it, they cant slightly (weak ferromagnetism).
    /// The altermagnetic properties are most relevant in the collinear
    /// antiferromagnetic phase below the Morin transition.
    ///
    /// Reference: L. Smejkal et al., Phys. Rev. X 12, 031042 (2022)
    pub fn fe2o3() -> Self {
        Self {
            name: "Fe2O3",
            neel_temperature: 960.0,
            spin_splitting: 0.6,
            symmetry: AltermagneticSymmetry::DWave,
            crystal_structure: "Corundum (R-3c)",
            sublattice_magnetization: 8.0e5,
            lattice_constant: 5.04e-10,
        }
    }

    /// Create a custom altermagnet material.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Neel temperature is not positive
    /// - Spin splitting is negative
    /// - Sublattice magnetization is not positive
    /// - Lattice constant is not positive
    pub fn custom(
        name: &'static str,
        neel_temperature: f64,
        spin_splitting: f64,
        symmetry: AltermagneticSymmetry,
        crystal_structure: &'static str,
        sublattice_magnetization: f64,
        lattice_constant: f64,
    ) -> Result<Self> {
        if neel_temperature <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "neel_temperature".to_string(),
                reason: "Neel temperature must be positive".to_string(),
            });
        }
        if spin_splitting < 0.0 {
            return Err(Error::InvalidParameter {
                param: "spin_splitting".to_string(),
                reason: "Spin splitting must be non-negative".to_string(),
            });
        }
        if sublattice_magnetization <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "sublattice_magnetization".to_string(),
                reason: "Sublattice magnetization must be positive".to_string(),
            });
        }
        if lattice_constant <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "lattice_constant".to_string(),
                reason: "Lattice constant must be positive".to_string(),
            });
        }
        Ok(Self {
            name,
            neel_temperature,
            spin_splitting,
            symmetry,
            crystal_structure,
            sublattice_magnetization,
            lattice_constant,
        })
    }

    /// Net magnetization of the altermagnet.
    ///
    /// By definition, altermagnets have zero net magnetization because
    /// the two sublattices have equal and opposite magnetic moments.
    /// This is a fundamental property distinguishing them from ferromagnets.
    pub fn net_magnetization(&self) -> f64 {
        0.0
    }

    /// Compute the spin splitting at a given angle φ in the Brillouin zone \[eV\].
    ///
    /// The angular dependence follows the symmetry of the altermagnet:
    /// Δ(φ) = Δ₀ * cos(n*φ), where n is the harmonic order.
    ///
    /// # Arguments
    /// * `phi` - Angle in the Brillouin zone \[radians\]
    pub fn splitting_at_angle(&self, phi: f64) -> f64 {
        self.spin_splitting * self.symmetry.angular_factor(phi)
    }

    /// Compute the temperature-dependent spin splitting \[eV\].
    ///
    /// The spin splitting follows the sublattice order parameter, which
    /// vanishes at the Neel temperature following a mean-field power law:
    ///
    /// Δ(T) = Δ₀ * (1 - T/T_N)^β
    ///
    /// where β ≈ 0.35 is the critical exponent for a 3D Heisenberg
    /// antiferromagnet (between mean-field 0.5 and Ising 0.33).
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]
    ///
    /// # Errors
    ///
    /// Returns an error if temperature is negative.
    pub fn splitting_at_temperature(&self, temperature: f64) -> Result<f64> {
        if temperature < 0.0 {
            return Err(Error::InvalidParameter {
                param: "temperature".to_string(),
                reason: "Temperature must be non-negative".to_string(),
            });
        }
        if temperature >= self.neel_temperature {
            return Ok(0.0);
        }
        // 3D Heisenberg critical exponent
        let beta = 0.35;
        let reduced_t = temperature / self.neel_temperature;
        Ok(self.spin_splitting * (1.0 - reduced_t).powf(beta))
    }

    /// Compute the angle-resolved and temperature-dependent spin splitting \[eV\].
    ///
    /// Combines angular symmetry and thermal suppression:
    /// Δ(φ, T) = Δ₀ * cos(n*φ) * (1 - T/T_N)^β
    ///
    /// # Arguments
    /// * `phi` - Angle in the Brillouin zone \[radians\]
    /// * `temperature` - Temperature \[K\]
    ///
    /// # Errors
    ///
    /// Returns an error if temperature is negative.
    pub fn splitting_at_angle_and_temperature(&self, phi: f64, temperature: f64) -> Result<f64> {
        let thermal_splitting = self.splitting_at_temperature(temperature)?;
        Ok(thermal_splitting * self.symmetry.angular_factor(phi))
    }

    /// Estimate the sublattice order parameter at a given temperature.
    ///
    /// The order parameter L(T) = M_sublattice(T) / M_sublattice(0) follows
    /// a power law near T_N: L(T) ∝ (1 - T/T_N)^β.
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]
    ///
    /// # Errors
    ///
    /// Returns an error if temperature is negative.
    pub fn order_parameter(&self, temperature: f64) -> Result<f64> {
        if temperature < 0.0 {
            return Err(Error::InvalidParameter {
                param: "temperature".to_string(),
                reason: "Temperature must be non-negative".to_string(),
            });
        }
        if temperature >= self.neel_temperature {
            return Ok(0.0);
        }
        let beta = 0.35;
        let reduced_t = temperature / self.neel_temperature;
        Ok((1.0 - reduced_t).powf(beta))
    }

    /// Estimate the magnetic susceptibility in the paramagnetic phase \[1/T\].
    ///
    /// Above T_N, the susceptibility follows a Curie-Weiss law:
    /// χ = C / (T + T_N)
    ///
    /// where C is the Curie constant estimated from the sublattice
    /// magnetization. Note the +T_N (antiferromagnetic Weiss constant).
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\], must be above T_N
    ///
    /// # Errors
    ///
    /// Returns an error if temperature is below T_N or negative.
    pub fn paramagnetic_susceptibility(&self, temperature: f64) -> Result<f64> {
        if temperature < 0.0 {
            return Err(Error::InvalidParameter {
                param: "temperature".to_string(),
                reason: "Temperature must be non-negative".to_string(),
            });
        }
        if temperature < self.neel_temperature {
            return Err(Error::InvalidParameter {
                param: "temperature".to_string(),
                reason: "Paramagnetic susceptibility is only valid above T_N".to_string(),
            });
        }
        // Estimate Curie constant from sublattice magnetization
        // C ≈ μ₀ M_sub² / (3 k_B) (simplified mean-field)
        let mu_0 = crate::constants::MU_0;
        let curie_constant =
            mu_0 * self.sublattice_magnetization * self.sublattice_magnetization / (3.0 * KB);
        Ok(curie_constant / (temperature + self.neel_temperature))
    }

    /// Validate the material parameters.
    ///
    /// Checks that all parameters are physically reasonable.
    ///
    /// # Errors
    ///
    /// Returns an error if any parameter is out of physical range.
    pub fn validate(&self) -> Result<()> {
        if self.neel_temperature <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "neel_temperature".to_string(),
                reason: "Neel temperature must be positive".to_string(),
            });
        }
        if self.spin_splitting < 0.0 {
            return Err(Error::InvalidParameter {
                param: "spin_splitting".to_string(),
                reason: "Spin splitting must be non-negative".to_string(),
            });
        }
        if self.sublattice_magnetization <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "sublattice_magnetization".to_string(),
                reason: "Sublattice magnetization must be positive".to_string(),
            });
        }
        if self.lattice_constant <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "lattice_constant".to_string(),
                reason: "Lattice constant must be positive".to_string(),
            });
        }
        if self.neel_temperature > 2000.0 {
            return Err(Error::InvalidParameter {
                param: "neel_temperature".to_string(),
                reason: "Neel temperature unreasonably high (>2000 K)".to_string(),
            });
        }
        if self.spin_splitting > 10.0 {
            return Err(Error::InvalidParameter {
                param: "spin_splitting".to_string(),
                reason: "Spin splitting unreasonably large (>10 eV)".to_string(),
            });
        }
        Ok(())
    }
}

impl fmt::Display for Altermagnet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Altermagnet({}, {}, T_N={:.0} K, Δ={:.2} eV, {})",
            self.name,
            self.symmetry,
            self.neel_temperature,
            self.spin_splitting,
            self.crystal_structure
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ruo2_properties() {
        let ruo2 = Altermagnet::ruo2();
        assert_eq!(ruo2.name, "RuO2");
        assert_eq!(ruo2.symmetry, AltermagneticSymmetry::DWave);
        assert!((ruo2.neel_temperature - 300.0).abs() < f64::EPSILON);
        assert!((ruo2.spin_splitting - 1.4).abs() < f64::EPSILON);
        ruo2.validate().expect("RuO2 parameters should be valid");
    }

    #[test]
    fn test_crsb_properties() {
        let crsb = Altermagnet::crsb();
        assert_eq!(crsb.name, "CrSb");
        assert_eq!(crsb.symmetry, AltermagneticSymmetry::GWave);
        assert!((crsb.neel_temperature - 700.0).abs() < f64::EPSILON);
        assert!((crsb.spin_splitting - 1.2).abs() < f64::EPSILON);
        crsb.validate().expect("CrSb parameters should be valid");
    }

    #[test]
    fn test_mnte_properties() {
        let mnte = Altermagnet::mnte();
        assert_eq!(mnte.name, "MnTe");
        assert_eq!(mnte.symmetry, AltermagneticSymmetry::DWave);
        assert!((mnte.neel_temperature - 307.0).abs() < f64::EPSILON);
        assert!((mnte.spin_splitting - 0.9).abs() < f64::EPSILON);
        mnte.validate().expect("MnTe parameters should be valid");
    }

    #[test]
    fn test_fe2o3_properties() {
        let fe2o3 = Altermagnet::fe2o3();
        assert_eq!(fe2o3.name, "Fe2O3");
        assert_eq!(fe2o3.symmetry, AltermagneticSymmetry::DWave);
        assert!((fe2o3.neel_temperature - 960.0).abs() < f64::EPSILON);
        assert!((fe2o3.spin_splitting - 0.6).abs() < f64::EPSILON);
        fe2o3.validate().expect("Fe2O3 parameters should be valid");
    }

    #[test]
    fn test_zero_net_magnetization() {
        // Fundamental property: all altermagnets have zero net magnetization
        let materials = [
            Altermagnet::ruo2(),
            Altermagnet::crsb(),
            Altermagnet::mnte(),
            Altermagnet::fe2o3(),
        ];
        for mat in &materials {
            assert!(
                mat.net_magnetization().abs() < 1e-10,
                "{} should have zero net magnetization",
                mat.name
            );
        }
    }

    #[test]
    fn test_dwave_nodes_at_45_degrees() {
        let ruo2 = Altermagnet::ruo2();
        // d-wave: nodes at 45, 135, 225, 315 degrees
        let node_angles = [PI / 4.0, 3.0 * PI / 4.0, 5.0 * PI / 4.0, 7.0 * PI / 4.0];
        for &angle in &node_angles {
            let splitting = ruo2.splitting_at_angle(angle);
            assert!(
                splitting.abs() < 1e-10,
                "d-wave splitting should be zero at {:.1} degrees, got {}",
                angle.to_degrees(),
                splitting
            );
        }
        // Maxima at 0, 90, 180, 270 degrees
        let max_angle = 0.0;
        let splitting_max = ruo2.splitting_at_angle(max_angle);
        assert!(
            (splitting_max.abs() - ruo2.spin_splitting).abs() < 1e-10,
            "d-wave splitting should be maximal at 0 degrees"
        );
    }

    #[test]
    fn test_gwave_nodes_at_22_5_degrees() {
        let crsb = Altermagnet::crsb();
        // g-wave: first node at 22.5 degrees = PI/8
        let node_angle = PI / 8.0;
        let splitting = crsb.splitting_at_angle(node_angle);
        assert!(
            splitting.abs() < 1e-10,
            "g-wave splitting should be zero at 22.5 degrees, got {}",
            splitting
        );
        // Check node count
        assert_eq!(crsb.symmetry.num_nodes(), 8);
    }

    #[test]
    fn test_splitting_vanishes_at_neel_temperature() {
        let materials = [
            Altermagnet::ruo2(),
            Altermagnet::crsb(),
            Altermagnet::mnte(),
            Altermagnet::fe2o3(),
        ];
        for mat in &materials {
            let splitting = mat
                .splitting_at_temperature(mat.neel_temperature)
                .expect("Should compute splitting at T_N");
            assert!(
                splitting.abs() < 1e-10,
                "{}: splitting should vanish at T_N={} K, got {}",
                mat.name,
                mat.neel_temperature,
                splitting
            );
            // Also vanishes above T_N
            let splitting_above = mat
                .splitting_at_temperature(mat.neel_temperature + 100.0)
                .expect("Should compute splitting above T_N");
            assert!(
                splitting_above.abs() < 1e-10,
                "{}: splitting should vanish above T_N",
                mat.name
            );
        }
    }

    #[test]
    fn test_temperature_dependence_monotonic() {
        let ruo2 = Altermagnet::ruo2();
        let mut prev_splitting = ruo2
            .splitting_at_temperature(0.0)
            .expect("Should compute splitting at T=0");
        // Splitting should decrease monotonically with temperature
        for t_step in 1..=100 {
            let t = (t_step as f64) * ruo2.neel_temperature / 100.0;
            let splitting = ruo2
                .splitting_at_temperature(t)
                .expect("Should compute splitting");
            assert!(
                splitting <= prev_splitting + 1e-15,
                "Splitting should decrease monotonically: at T={} K got {}, prev was {}",
                t,
                splitting,
                prev_splitting
            );
            prev_splitting = splitting;
        }
    }

    #[test]
    fn test_order_parameter_range() {
        let mnte = Altermagnet::mnte();
        let op_zero = mnte.order_parameter(0.0).expect("Should compute at T=0");
        assert!(
            (op_zero - 1.0).abs() < 1e-10,
            "Order parameter should be 1.0 at T=0"
        );
        let op_tn = mnte
            .order_parameter(mnte.neel_temperature)
            .expect("Should compute at T_N");
        assert!(op_tn.abs() < 1e-10, "Order parameter should be 0.0 at T_N");
    }

    #[test]
    fn test_symmetry_node_angles() {
        let d_nodes = AltermagneticSymmetry::DWave.node_angles();
        assert_eq!(d_nodes.len(), 4);
        assert!((d_nodes[0] - PI / 4.0).abs() < 1e-10);

        let g_nodes = AltermagneticSymmetry::GWave.node_angles();
        assert_eq!(g_nodes.len(), 8);
        assert!((g_nodes[0] - PI / 8.0).abs() < 1e-10);

        let i_nodes = AltermagneticSymmetry::IWave.node_angles();
        assert_eq!(i_nodes.len(), 12);
        assert!((i_nodes[0] - PI / 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_custom_material_validation() {
        // Valid custom material
        let mat = Altermagnet::custom(
            "TestMat",
            500.0,
            0.8,
            AltermagneticSymmetry::DWave,
            "Test structure",
            2.0e5,
            4.0e-10,
        );
        assert!(mat.is_ok());

        // Invalid: negative Neel temperature
        let invalid = Altermagnet::custom(
            "Bad",
            -100.0,
            0.8,
            AltermagneticSymmetry::DWave,
            "Test",
            2.0e5,
            4.0e-10,
        );
        assert!(invalid.is_err());

        // Invalid: negative spin splitting
        let invalid = Altermagnet::custom(
            "Bad",
            300.0,
            -0.5,
            AltermagneticSymmetry::DWave,
            "Test",
            2.0e5,
            4.0e-10,
        );
        assert!(invalid.is_err());
    }

    #[test]
    fn test_negative_temperature_error() {
        let ruo2 = Altermagnet::ruo2();
        assert!(ruo2.splitting_at_temperature(-10.0).is_err());
        assert!(ruo2.order_parameter(-10.0).is_err());
        assert!(ruo2.paramagnetic_susceptibility(-10.0).is_err());
    }

    #[test]
    fn test_paramagnetic_susceptibility_above_tn() {
        let ruo2 = Altermagnet::ruo2();
        // Should work above T_N
        let chi = ruo2
            .paramagnetic_susceptibility(600.0)
            .expect("Should compute susceptibility above T_N");
        assert!(chi > 0.0, "Susceptibility must be positive");

        // Should fail below T_N
        assert!(ruo2.paramagnetic_susceptibility(200.0).is_err());
    }

    #[test]
    fn test_symmetry_comparison_between_types() {
        // Higher symmetry order means more nodes and faster oscillation
        let d = AltermagneticSymmetry::DWave;
        let g = AltermagneticSymmetry::GWave;
        let i = AltermagneticSymmetry::IWave;

        assert!(d.harmonic_order() < g.harmonic_order());
        assert!(g.harmonic_order() < i.harmonic_order());
        assert!(d.num_nodes() < g.num_nodes());
        assert!(g.num_nodes() < i.num_nodes());
    }

    #[test]
    fn test_display_formatting() {
        let ruo2 = Altermagnet::ruo2();
        let display = format!("{}", ruo2);
        assert!(display.contains("RuO2"));
        assert!(display.contains("d-wave"));
        assert!(display.contains("300"));

        let d = AltermagneticSymmetry::DWave;
        assert_eq!(format!("{}", d), "d-wave");
        let g = AltermagneticSymmetry::GWave;
        assert_eq!(format!("{}", g), "g-wave");
    }

    #[test]
    fn test_combined_angle_and_temperature_splitting() {
        let mnte = Altermagnet::mnte();
        // At T=0, angle=0: should give full splitting
        let full = mnte
            .splitting_at_angle_and_temperature(0.0, 0.0)
            .expect("Should compute at T=0, phi=0");
        assert!(
            (full - mnte.spin_splitting).abs() < 1e-10,
            "At T=0, phi=0 should give full splitting"
        );
        // At T=0, node angle: should give zero
        let at_node = mnte
            .splitting_at_angle_and_temperature(PI / 4.0, 0.0)
            .expect("Should compute at node");
        assert!(
            at_node.abs() < 1e-10,
            "At node angle should give zero regardless of temperature"
        );
        // At T=T_N, any angle: should give zero
        let at_tn = mnte
            .splitting_at_angle_and_temperature(0.0, mnte.neel_temperature)
            .expect("Should compute at T_N");
        assert!(at_tn.abs() < 1e-10, "At T_N should give zero");
    }
}
