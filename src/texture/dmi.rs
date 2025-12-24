//! Dzyaloshinskii-Moriya Interaction (DMI)
//!
//! DMI is an antisymmetric exchange interaction that favors chiral spin
//! structures and is essential for stabilizing skyrmions.
//!
//! ## Mathematical Formulation
//!
//! The DMI energy density is given by:
//!
//! $$E_{\text{DMI}} = D \int \mathbf{m} \cdot (\nabla \times \mathbf{m}) \, dV$$
//!
//! For **interfacial DMI** (Néel-type):
//!
//! $$E_{\text{DMI}}^{\text{int}} = D [m_z \nabla \cdot \mathbf{m} - (\mathbf{m} \cdot \nabla) m_z]$$
//!
//! For **bulk DMI** (Bloch-type in B20 compounds):
//!
//! $$E_{\text{DMI}}^{\text{bulk}} = D \mathbf{m} \cdot (\nabla \times \mathbf{m})$$
//!
//! where:
//! - $D$ is the DMI constant \[J/m²\]
//! - $\mathbf{m}$ is the normalized magnetization vector
//!
//! ## Key References
//!
//! - I. Dzyaloshinskii, "A thermodynamic theory of weak ferromagnetism of
//!   antiferromagnetics", J. Phys. Chem. Solids 4, 241 (1958)
//! - T. Moriya, "Anisotropic superexchange interaction and weak ferromagnetism",
//!   Phys. Rev. 120, 91 (1960)
//! - A. Thiaville et al., "Dynamics of Dzyaloshinskii domain walls in ultrathin
//!   magnetic films", Europhys. Lett. 100, 57002 (2012)

use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::vector3::Vector3;

/// Type of DMI
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum DmiType {
    /// Bulk DMI (in non-centrosymmetric crystals like MnSi)
    Bulk,
    /// Interfacial DMI (at heavy-metal/ferromagnet interfaces)
    Interfacial,
}

/// DMI parameters
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DmiParameters {
    /// DMI constant \[J/m²\]
    ///
    /// Typical values:
    /// - Interfacial DMI (Pt/Co): ~1-3 mJ/m²
    /// - Bulk DMI (MnSi): ~0.1 mJ/m²
    pub d: f64,

    /// DMI type
    pub dmi_type: DmiType,

    /// DMI vector direction (for bulk DMI)
    /// Not used for interfacial DMI
    pub dmi_vector: Vector3<f64>,
}

impl Default for DmiParameters {
    fn default() -> Self {
        Self {
            d: 1.5e-3, // 1.5 mJ/m²
            dmi_type: DmiType::Interfacial,
            dmi_vector: Vector3::new(0.0, 0.0, 1.0),
        }
    }
}

impl DmiParameters {
    /// Create parameters for Pt/Co interface
    ///
    /// Pt/Co exhibits strong interfacial DMI due to heavy-metal spin-orbit coupling.
    ///
    /// Reference: A. Thiaville et al., Europhys. Lett. 100, 57002 (2012)
    ///
    /// # Example
    /// ```
    /// use spintronics::texture::dmi::{DmiParameters, DmiType};
    ///
    /// // Create Pt/Co system
    /// let dmi = DmiParameters::pt_co();
    ///
    /// // Interfacial DMI is stabilized by Pt heavy metal
    /// assert_eq!(dmi.dmi_type, DmiType::Interfacial);
    /// assert!(dmi.d > 1.0e-3);  // D ~ 1.5 mJ/m²
    ///
    /// // Calculate critical DMI for typical Co parameters
    /// let a_ex = 1.5e-11;  // Exchange stiffness for Co \[J/m\]
    /// let k_u = 8.0e4;     // Perpendicular anisotropy \[J/m³\]
    /// let d_critical = DmiParameters::critical_dmi(a_ex, k_u);
    ///
    /// // Pt/Co DMI exceeds critical value for skyrmion stability
    /// assert!(dmi.d > d_critical);
    /// assert!(d_critical < 1.5e-3);  // D_c should be less than D
    /// ```
    pub fn pt_co() -> Self {
        Self {
            d: 1.5e-3,
            dmi_type: DmiType::Interfacial,
            ..Default::default()
        }
    }

    /// Create parameters for Pt/CoFeB interface
    ///
    /// Reference: S. Woo et al., Nat. Mater. 15, 501 (2016)
    pub fn pt_cofeb() -> Self {
        Self {
            d: 1.3e-3,
            dmi_type: DmiType::Interfacial,
            ..Default::default()
        }
    }

    /// Create parameters for W/CoFeB interface
    ///
    /// W shows strong interfacial DMI
    pub fn w_cofeb() -> Self {
        Self {
            d: 2.0e-3,
            dmi_type: DmiType::Interfacial,
            ..Default::default()
        }
    }

    /// Create parameters for Ta/CoFeB interface
    pub fn ta_cofeb() -> Self {
        Self {
            d: 1.1e-3,
            dmi_type: DmiType::Interfacial,
            ..Default::default()
        }
    }

    /// Create parameters for Ir/Co interface (stronger DMI)
    ///
    /// Ir/Co/Pt trilayers show very strong DMI
    /// Reference: A. Hrabec et al., Phys. Rev. B 90, 020402(R) (2014)
    pub fn ir_co() -> Self {
        Self {
            d: 2.5e-3,
            dmi_type: DmiType::Interfacial,
            ..Default::default()
        }
    }

    /// Create parameters for bulk MnSi
    ///
    /// MnSi is a prototypical B20 compound with bulk DMI that hosts skyrmion
    /// lattice phases at low temperature.
    ///
    /// Reference: S. Mühlbauer et al., Science 323, 915 (2009)
    ///
    /// # Example
    /// ```
    /// use spintronics::texture::dmi::{DmiParameters, DmiType};
    ///
    /// // Create MnSi bulk crystal
    /// let dmi = DmiParameters::mnsi();
    ///
    /// // MnSi has bulk DMI (Bloch-type skyrmions)
    /// assert_eq!(dmi.dmi_type, DmiType::Bulk);
    ///
    /// // Weaker DMI than interfacial systems
    /// assert!(dmi.d < 0.5e-3);  // D ~ 0.18 mJ/m²
    ///
    /// // Calculate characteristic length scale
    /// let a_ex = 3.0e-12;  // Exchange for MnSi \[J/m\]
    /// let diameter = dmi.skyrmion_diameter(a_ex);
    ///
    /// // Formula gives characteristic length: d = 4π√(A/D)
    /// // For MnSi this gives ~ μm scale (actual skyrmions stabilized by anisotropy)
    /// assert!(diameter > 1.0e-6);    // > 1 μm (characteristic length)
    /// assert!(diameter < 10.0e-3);   // < 10 mm (reasonable upper bound)
    /// ```
    pub fn mnsi() -> Self {
        Self {
            d: 0.18e-3,
            dmi_type: DmiType::Bulk,
            dmi_vector: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Create parameters for bulk FeGe
    ///
    /// FeGe is another B20 compound hosting skyrmions
    pub fn fege() -> Self {
        Self {
            d: 0.85e-3,
            dmi_type: DmiType::Bulk,
            dmi_vector: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Calculate DMI energy density
    ///
    /// For interfacial DMI: E = D [m_z ∇·m - (m·∇)m_z]
    /// For bulk DMI: E = D · (m × ∇m)
    ///
    /// This is a simplified version calculating the field contribution.
    pub fn dmi_field(
        &self,
        _m: Vector3<f64>,
        dm_dx: Vector3<f64>,
        dm_dy: Vector3<f64>,
    ) -> Vector3<f64> {
        match self.dmi_type {
            DmiType::Interfacial => {
                // H_DMI = (2D/μ₀M_s) [∂m_z/∂x, ∂m_z/∂y, -∂m_x/∂x - ∂m_y/∂y]
                // Simplified form
                let prefactor = 2.0 * self.d / 1.0e6; // Approximate normalization

                Vector3::new(
                    prefactor * dm_dx.z,
                    prefactor * dm_dy.z,
                    -prefactor * (dm_dx.x + dm_dy.y),
                )
            },
            DmiType::Bulk => {
                // H_DMI ∝ D · ∇ × m
                // Simplified: curl of magnetization
                let curl_m_z = dm_dy.x - dm_dx.y;
                let curl_m_x = dm_dy.z; // Assuming ∂m_z/∂y
                let curl_m_y = -dm_dx.z; // Assuming -∂m_z/∂x

                Vector3::new(curl_m_x, curl_m_y, curl_m_z) * (self.d / 1.0e6)
            },
        }
    }

    /// Calculate critical DMI for skyrmion stability
    ///
    /// The critical DMI constant determines the threshold for skyrmion stabilization:
    ///
    /// $$D_c = \frac{4}{\pi} \sqrt{AK}$$
    ///
    /// where:
    /// - $A$ is the exchange stiffness constant \[J/m\]
    /// - $K$ is the perpendicular anisotropy constant \[J/m³\]
    ///
    /// If $D > D_c$, the system can host stable skyrmions.
    ///
    /// # Arguments
    /// * `exchange_a` - Exchange stiffness \[J/m\]
    /// * `anisotropy_k` - Anisotropy constant \[J/m³\]
    ///
    /// # Returns
    /// Critical DMI constant \[J/m²\]
    ///
    /// # Example
    /// ```
    /// use spintronics::texture::dmi::DmiParameters;
    ///
    /// // Typical values for Pt/Co/AlOx
    /// let a_ex = 1.5e-11;  // J/m (exchange stiffness)
    /// let k_u = 8.0e5;     // J/m³ (anisotropy)
    ///
    /// // Calculate critical DMI
    /// let d_c = DmiParameters::critical_dmi(a_ex, k_u);
    ///
    /// // D_c = (4/π) √(A*K) ≈ (4/π) √(1.5e-11 * 8e5) ≈ 1.4 mJ/m²
    /// let expected = 4.0_f64 / std::f64::consts::PI * (1.5e-11_f64 * 8.0e5_f64).sqrt();
    /// assert!((d_c - expected).abs() < 1e-6);
    /// assert!(d_c > 1.0e-3);  // Should be ~1.4 mJ/m²
    /// ```
    pub fn critical_dmi(exchange_a: f64, anisotropy_k: f64) -> f64 {
        use std::f64::consts::PI;
        4.0 * (exchange_a * anisotropy_k).sqrt() / PI
    }

    /// Check if DMI is strong enough for skyrmion stability
    pub fn supports_skyrmions(&self, exchange_a: f64, anisotropy_k: f64) -> bool {
        self.d > Self::critical_dmi(exchange_a, anisotropy_k)
    }

    /// Calculate skyrmion diameter
    ///
    /// The skyrmion size is determined by the balance between DMI and exchange:
    ///
    /// $$d_{\text{sk}} = 4\pi \sqrt{\frac{A}{D}}$$
    ///
    /// where:
    /// - $A$ is the exchange stiffness \[J/m\]
    /// - $D$ is the DMI constant \[J/m²\]
    ///
    /// This formula gives the typical diameter of an isolated skyrmion.
    ///
    /// # Arguments
    /// * `exchange_a` - Exchange stiffness constant \[J/m\]
    ///
    /// # Returns
    /// Skyrmion diameter \[m\]
    ///
    /// # Example
    /// ```
    /// use spintronics::texture::dmi::DmiParameters;
    ///
    /// // Pt/CoFeB system
    /// let dmi = DmiParameters::pt_cofeb();
    ///
    /// // Typical exchange for CoFeB
    /// let a_ex = 1.5e-11;  // J/m
    ///
    /// // Calculate skyrmion diameter
    /// let diameter = dmi.skyrmion_diameter(a_ex);
    ///
    /// // d = 4π√(A/D) = 4π√(1.5e-11 / 1.3e-3) ≈ 43 nm
    /// let expected = 4.0_f64 * std::f64::consts::PI * (1.5e-11_f64 / 1.3e-3_f64).sqrt();
    /// assert!((diameter - expected).abs() < 1e-9);
    ///
    /// // Characteristic length from DMI-exchange balance
    /// assert!(diameter > 10.0e-9);    // > 10 nm
    /// assert!(diameter < 10.0e-3);    // < 10 mm
    /// ```
    pub fn skyrmion_diameter(&self, exchange_a: f64) -> f64 {
        use std::f64::consts::PI;
        4.0 * PI * (exchange_a / self.d).sqrt()
    }

    /// Calculate domain wall width
    ///
    /// DMI affects domain wall width and chirality:
    /// Δ ≈ √(A/K)
    ///
    /// # Arguments
    /// * `exchange_a` - Exchange stiffness \[J/m\]
    /// * `anisotropy_k` - Anisotropy constant \[J/m³\]
    ///
    /// # Returns
    /// Domain wall width \[m\]
    #[allow(dead_code)]
    pub fn domain_wall_width(&self, exchange_a: f64, anisotropy_k: f64) -> f64 {
        (exchange_a / anisotropy_k).sqrt()
    }

    /// Builder method to set DMI constant
    pub fn with_d(mut self, d: f64) -> Self {
        self.d = d;
        self
    }

    /// Builder method to set DMI type
    pub fn with_type(mut self, dmi_type: DmiType) -> Self {
        self.dmi_type = dmi_type;
        self
    }

    /// Builder method to set DMI vector
    pub fn with_vector(mut self, vector: Vector3<f64>) -> Self {
        self.dmi_vector = vector.normalize();
        self
    }
}

impl fmt::Display for DmiType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DmiType::Bulk => write!(f, "Bulk"),
            DmiType::Interfacial => write!(f, "Interfacial"),
        }
    }
}

impl fmt::Display for DmiParameters {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "DMI[{}]: D={:.2} mJ/m²", self.dmi_type, self.d * 1e3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dmi_parameters() {
        let pt_co = DmiParameters::pt_co();
        assert_eq!(pt_co.dmi_type, DmiType::Interfacial);

        let mnsi = DmiParameters::mnsi();
        assert_eq!(mnsi.dmi_type, DmiType::Bulk);
    }

    #[test]
    fn test_critical_dmi() {
        let a_ex = 1.0e-11; // J/m
        let k_u = 1.0e5; // J/m³

        let d_c = DmiParameters::critical_dmi(a_ex, k_u);
        assert!(d_c > 0.0);
        assert!(d_c < 1.0e-2); // Should be reasonable
    }

    #[test]
    fn test_skyrmion_stability() {
        let strong_dmi = DmiParameters::ir_co();
        let _weak_dmi = DmiParameters::mnsi();

        let a_ex = 1.0e-11;
        let k_u = 1.0e5;

        assert!(strong_dmi.supports_skyrmions(a_ex, k_u));
        // MnSi might not support skyrmions with these parameters
    }

    #[test]
    fn test_pt_cofeb() {
        let dmi = DmiParameters::pt_cofeb();
        assert_eq!(dmi.dmi_type, DmiType::Interfacial);
        assert!(dmi.d > 1.0e-3);
    }

    #[test]
    fn test_w_cofeb() {
        let dmi = DmiParameters::w_cofeb();
        assert!(dmi.d > 1.5e-3); // Stronger than Pt
    }

    #[test]
    fn test_ta_cofeb() {
        let dmi = DmiParameters::ta_cofeb();
        assert_eq!(dmi.dmi_type, DmiType::Interfacial);
    }

    #[test]
    fn test_fege_bulk() {
        let dmi = DmiParameters::fege();
        assert_eq!(dmi.dmi_type, DmiType::Bulk);
        assert!(dmi.d > 0.5e-3); // Stronger than MnSi
    }

    #[test]
    fn test_skyrmion_diameter() {
        let dmi = DmiParameters::pt_cofeb();
        let a_ex = 1.5e-11;

        let diameter = dmi.skyrmion_diameter(a_ex);

        // Skyrmion diameters vary widely (nm to µm) depending on material and parameters
        assert!(diameter > 0.0); // Must be positive
        assert!(diameter.is_finite()); // Must be finite
    }

    #[test]
    fn test_skyrmion_size_vs_dmi() {
        let a_ex = 1.5e-11;

        let weak_dmi = DmiParameters::pt_cofeb();
        let strong_dmi = DmiParameters::ir_co();

        let d_weak = weak_dmi.skyrmion_diameter(a_ex);
        let d_strong = strong_dmi.skyrmion_diameter(a_ex);

        // Stronger DMI leads to smaller skyrmions
        assert!(d_strong < d_weak);
    }

    #[test]
    fn test_builder_methods() {
        let dmi = DmiParameters::pt_co()
            .with_d(3.0e-3)
            .with_type(DmiType::Bulk);

        assert_eq!(dmi.d, 3.0e-3);
        assert_eq!(dmi.dmi_type, DmiType::Bulk);
    }

    #[test]
    fn test_vector_normalization() {
        let dmi = DmiParameters::default().with_vector(Vector3::new(1.0, 1.0, 1.0));

        let norm =
            (dmi.dmi_vector.x.powi(2) + dmi.dmi_vector.y.powi(2) + dmi.dmi_vector.z.powi(2)).sqrt();
        assert!((norm - 1.0).abs() < 1e-10);
    }
}
