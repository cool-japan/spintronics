//! Magnetic Multilayer Structures
//!
//! This module implements magnetic multilayer systems including:
//! - Giant Magnetoresistance (GMR) spin valves
//! - Tunnel Magnetoresistance (TMR) magnetic tunnel junctions
//! - Synthetic Antiferromagnets (SAF)
//! - Exchange-biased structures
//!
//! ## Physics Background
//!
//! ### Giant Magnetoresistance (GMR):
//! Discovered by Fert and Grünberg (2007 Nobel Prize), GMR occurs in
//! ferromagnet/non-magnet/ferromagnet trilayers:
//!
//! ΔR/R = (R_AP - R_P) / R_P
//!
//! where R_P and R_AP are resistances for parallel and antiparallel configurations.
//!
//! ### Tunnel Magnetoresistance (TMR):
//! Similar to GMR but with insulating barrier:
//!
//! TMR = 2P₁P₂ / (1 - P₁P₂)
//!
//! where P₁, P₂ are spin polarizations of the electrodes.
//!
//! ### Synthetic Antiferromagnet (SAF):
//! Two ferromagnetic layers coupled antiferromagnetically through
//! RKKY interaction via thin spacer layer.
//!
//! ## Key References
//!
//! - M. N. Baibich et al., "Giant Magnetoresistance of (001)Fe/(001)Cr Magnetic
//!   Superlattices", Phys. Rev. Lett. 61, 2472 (1988)
//! - G. Binasch et al., "Enhanced magnetoresistance in layered magnetic
//!   structures with antiferromagnetic interlayer exchange", Phys. Rev. B 39, 4828 (1989)
//! - J. S. Moodera et al., "Large Magnetoresistance at Room Temperature in
//!   Ferromagnetic Thin Film Tunnel Junctions", Phys. Rev. Lett. 74, 3273 (1995)
//! - S. S. P. Parkin et al., "Giant tunnelling magnetoresistance at room
//!   temperature with MgO (100) tunnel barriers", Nat. Mater. 3, 862 (2004)

use crate::material::Ferromagnet;

/// Type of magnetic multilayer structure
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MultilayerType {
    /// GMR spin valve (FM/NM/FM)
    GmrSpinValve,
    /// TMR magnetic tunnel junction (FM/I/FM)
    TmrJunction,
    /// Synthetic antiferromagnet (FM/NM/FM with RKKY coupling)
    SyntheticAntiferromagnet,
    /// Exchange-biased structure (FM/AFM)
    ExchangeBiased,
}

/// Spacer layer material
#[derive(Debug, Clone)]
pub struct SpacerLayer {
    /// Material name (e.g., "Cu", "Cr", "MgO", "Al₂O₃")
    pub material: String,

    /// Thickness \[nm\]
    pub thickness: f64,

    /// Resistivity [Ω·m]
    pub resistivity: f64,

    /// Is insulating (for TMR)
    pub is_insulator: bool,

    /// Spin diffusion length \[nm\] (for metals)
    pub spin_diffusion_length: f64,
}

impl SpacerLayer {
    /// Create Cu (Copper) spacer for GMR
    pub fn copper(thickness: f64) -> Self {
        Self {
            material: "Cu".to_string(),
            thickness,
            resistivity: 1.7e-8, // Ω·m
            is_insulator: false,
            spin_diffusion_length: 500.0, // nm (very long)
        }
    }

    /// Create Cr (Chromium) spacer for SAF
    pub fn chromium(thickness: f64) -> Self {
        Self {
            material: "Cr".to_string(),
            thickness,
            resistivity: 1.3e-7,
            is_insulator: false,
            spin_diffusion_length: 5.0, // nm (short)
        }
    }

    /// Create MgO barrier for TMR
    pub fn mgo(thickness: f64) -> Self {
        Self {
            material: "MgO".to_string(),
            thickness,
            resistivity: 1.0e14, // Very high (insulator)
            is_insulator: true,
            spin_diffusion_length: 0.0, // N/A for insulator
        }
    }

    /// Create Al₂O₃ barrier for TMR
    pub fn al2o3(thickness: f64) -> Self {
        Self {
            material: "Al₂O₃".to_string(),
            thickness,
            resistivity: 1.0e13,
            is_insulator: true,
            spin_diffusion_length: 0.0,
        }
    }
}

/// Magnetic multilayer structure
#[derive(Debug, Clone)]
pub struct MagneticMultilayer {
    /// Structure type
    pub structure_type: MultilayerType,

    /// Bottom ferromagnetic layer (fixed/reference)
    pub bottom_layer: Ferromagnet,

    /// Bottom layer thickness \[nm\]
    pub bottom_thickness: f64,

    /// Spacer layer
    pub spacer: SpacerLayer,

    /// Top ferromagnetic layer (free)
    pub top_layer: Ferromagnet,

    /// Top layer thickness \[nm\]
    pub top_thickness: f64,

    /// RKKY coupling strength \[J/m²\] (for SAF)
    pub rkky_coupling: f64,

    /// Exchange bias field \[T\] (for exchange-biased)
    pub exchange_bias: f64,
}

impl MagneticMultilayer {
    /// Create GMR spin valve: Co/Cu/Co
    ///
    /// Reference: M. N. Baibich et al., PRL 61, 2472 (1988)
    pub fn gmr_co_cu_co() -> Self {
        Self {
            structure_type: MultilayerType::GmrSpinValve,
            bottom_layer: Ferromagnet::cobalt(),
            bottom_thickness: 10.0, // nm
            spacer: SpacerLayer::copper(2.0),
            top_layer: Ferromagnet::cobalt(),
            top_thickness: 5.0,
            rkky_coupling: 0.0,
            exchange_bias: 0.0,
        }
    }

    /// Create GMR spin valve: Py/Cu/Py (Permalloy-based)
    pub fn gmr_py_cu_py() -> Self {
        Self {
            structure_type: MultilayerType::GmrSpinValve,
            bottom_layer: Ferromagnet::permalloy(),
            bottom_thickness: 8.0,
            spacer: SpacerLayer::copper(2.5),
            top_layer: Ferromagnet::permalloy(),
            top_thickness: 4.0,
            rkky_coupling: 0.0,
            exchange_bias: 0.0,
        }
    }

    /// Create TMR junction: CoFeB/MgO/CoFeB
    ///
    /// High TMR at room temperature
    /// Reference: S. S. P. Parkin et al., Nat. Mater. 3, 862 (2004)
    pub fn tmr_cofeb_mgo_cofeb() -> Self {
        Self {
            structure_type: MultilayerType::TmrJunction,
            bottom_layer: Ferromagnet::cofeb(),
            bottom_thickness: 3.0,
            spacer: SpacerLayer::mgo(1.2), // ~1.2 nm optimal
            top_layer: Ferromagnet::cofeb(),
            top_thickness: 2.0,
            rkky_coupling: 0.0,
            exchange_bias: 0.0,
        }
    }

    /// Create SAF structure: CoFe/Ru/CoFe
    ///
    /// Strong antiferromagnetic coupling through RKKY
    pub fn saf_cofe_ru_cofe() -> Self {
        let mut cofe = Ferromagnet::cofe();
        cofe.ms = 1.4e6; // Adjusted for SAF

        Self {
            structure_type: MultilayerType::SyntheticAntiferromagnet,
            bottom_layer: cofe.clone(),
            bottom_thickness: 4.0,
            spacer: SpacerLayer {
                material: "Ru".to_string(),
                thickness: 0.8, // First RKKY peak
                resistivity: 7.1e-8,
                is_insulator: false,
                spin_diffusion_length: 10.0,
            },
            top_layer: cofe,
            top_thickness: 4.0,
            rkky_coupling: -1.0e-3, // J/m² (antiferromagnetic)
            exchange_bias: 0.0,
        }
    }

    /// Calculate GMR ratio
    ///
    /// GMR = (R_AP - R_P) / R_P
    ///
    /// # Arguments
    /// * `angle` - Angle between magnetizations \[rad\]
    ///
    /// # Returns
    /// Magnetoresistance ratio
    pub fn gmr_ratio(&self, angle: f64) -> f64 {
        if self.structure_type != MultilayerType::GmrSpinValve {
            return 0.0;
        }

        // Simplified Valet-Fert model
        let gmr_max = 0.15; // 15% typical for Co/Cu/Co
        gmr_max * (1.0 - angle.cos()) / 2.0
    }

    /// Calculate TMR ratio
    ///
    /// TMR = 2P₁P₂ / (1 - P₁P₂) × (1 - cos θ) / 2
    ///
    /// # Arguments
    /// * `angle` - Angle between magnetizations \[rad\]
    ///
    /// # Returns
    /// TMR ratio
    pub fn tmr_ratio(&self, angle: f64) -> f64 {
        if self.structure_type != MultilayerType::TmrJunction {
            return 0.0;
        }

        // Spin polarizations
        let p1 = 0.7; // CoFeB with MgO
        let p2 = 0.7;

        let tmr_max = 2.0 * p1 * p2 / (1.0 - p1 * p2);
        tmr_max * (1.0 - angle.cos()) / 2.0
    }

    /// Calculate RKKY coupling field
    ///
    /// H_RKKY = J_RKKY / (μ₀ M_s t)
    ///
    /// # Returns
    /// RKKY coupling field \[T\]
    pub fn rkky_field(&self) -> f64 {
        if self.rkky_coupling == 0.0 {
            return 0.0;
        }

        let mu_0 = 4.0 * std::f64::consts::PI * 1e-7;
        let t = self.top_thickness * 1e-9; // nm to m
        let ms = self.top_layer.ms;

        self.rkky_coupling / (mu_0 * ms * t)
    }

    /// Calculate total resistance
    ///
    /// # Arguments
    /// * `angle` - Angle between magnetizations \[rad\]
    /// * `area` - Junction area \[m²\]
    ///
    /// # Returns
    /// Resistance \[Ω\]
    pub fn resistance(&self, angle: f64, area: f64) -> f64 {
        let r0 = 1.0e-12; // Base resistance-area product [Ω·m²]

        let mr = match self.structure_type {
            MultilayerType::GmrSpinValve => self.gmr_ratio(angle),
            MultilayerType::TmrJunction => self.tmr_ratio(angle),
            _ => 0.0,
        };

        r0 * (1.0 + mr) / area
    }

    /// Check if SAF is properly antiferromagnetically coupled
    pub fn is_properly_coupled_saf(&self) -> bool {
        self.structure_type == MultilayerType::SyntheticAntiferromagnet &&
        self.rkky_coupling < 0.0 && // Antiferromagnetic
        self.spacer.thickness < 2.0 // Reasonable spacer thickness
    }

    /// Builder: set RKKY coupling
    pub fn with_rkky_coupling(mut self, j_rkky: f64) -> Self {
        self.rkky_coupling = j_rkky;
        self
    }

    /// Builder: set spacer thickness
    pub fn with_spacer_thickness(mut self, thickness: f64) -> Self {
        self.spacer.thickness = thickness;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gmr_co_cu_co() {
        let gmr = MagneticMultilayer::gmr_co_cu_co();
        assert_eq!(gmr.structure_type, MultilayerType::GmrSpinValve);
        assert!(!gmr.spacer.is_insulator);
        assert_eq!(gmr.spacer.material, "Cu");
    }

    #[test]
    fn test_tmr_cofeb_mgo() {
        let tmr = MagneticMultilayer::tmr_cofeb_mgo_cofeb();
        assert_eq!(tmr.structure_type, MultilayerType::TmrJunction);
        assert!(tmr.spacer.is_insulator);
        assert_eq!(tmr.spacer.material, "MgO");
    }

    #[test]
    fn test_saf_structure() {
        let saf = MagneticMultilayer::saf_cofe_ru_cofe();
        assert_eq!(saf.structure_type, MultilayerType::SyntheticAntiferromagnet);
        assert!(saf.rkky_coupling < 0.0); // Antiferromagnetic
        assert!(saf.is_properly_coupled_saf());
    }

    #[test]
    fn test_gmr_ratio_parallel() {
        let gmr = MagneticMultilayer::gmr_co_cu_co();
        let ratio_p = gmr.gmr_ratio(0.0); // Parallel
        assert!(ratio_p.abs() < 0.01); // Nearly zero
    }

    #[test]
    fn test_gmr_ratio_antiparallel() {
        let gmr = MagneticMultilayer::gmr_co_cu_co();
        let ratio_ap = gmr.gmr_ratio(std::f64::consts::PI); // Antiparallel
        assert!(ratio_ap > 0.1); // Should be ~15%
    }

    #[test]
    fn test_tmr_ratio() {
        let tmr = MagneticMultilayer::tmr_cofeb_mgo_cofeb();
        let ratio_ap = tmr.tmr_ratio(std::f64::consts::PI);
        // TMR should be much larger than GMR
        assert!(ratio_ap > 1.0); // > 100%
    }

    #[test]
    fn test_rkky_field() {
        let saf = MagneticMultilayer::saf_cofe_ru_cofe();
        let h_rkky = saf.rkky_field();
        assert!(h_rkky < 0.0); // Negative for AFM coupling
        assert!(h_rkky.abs() > 0.01); // Significant field
    }

    #[test]
    fn test_resistance_calculation() {
        let gmr = MagneticMultilayer::gmr_co_cu_co();
        let area = 1.0e-14; // 100 nm × 100 nm

        let r_p = gmr.resistance(0.0, area);
        let r_ap = gmr.resistance(std::f64::consts::PI, area);

        assert!(r_ap > r_p); // Antiparallel has higher resistance
    }

    #[test]
    fn test_builder_pattern() {
        let saf = MagneticMultilayer::saf_cofe_ru_cofe()
            .with_rkky_coupling(-2.0e-3)
            .with_spacer_thickness(0.9);

        assert_eq!(saf.rkky_coupling, -2.0e-3);
        assert_eq!(saf.spacer.thickness, 0.9);
    }
}
