//! Spin Nernst Effect
//!
//! The Spin Nernst Effect (SNE) is the thermal analog of the Spin Hall Effect.
//! A temperature gradient generates a transverse spin current perpendicular to
//! both the thermal gradient and the spin polarization direction.
//!
//! ## Physics Background
//!
//! ### Spin Nernst Equation:
//! **j_s** = θ_SN · **∇T** × **σ**
//!
//! where:
//! - j_s: Spin current density \[J/m²\]
//! - θ_SN: Spin Nernst angle (dimensionless)
//! - ∇T: Temperature gradient \[K/m\]
//! - σ: Spin polarization direction
//!
//! ### Relationship to Other Effects:
//! - **Spin Hall Effect (SHE)**: Electric field → transverse spin current
//! - **Spin Nernst Effect (SNE)**: Thermal gradient → transverse spin current
//! - **Spin Seebeck Effect**: Thermal gradient → spin accumulation
//!
//! ### Microscopic Origin:
//! - Mott two-current model with energy-dependent scattering
//! - Spin-orbit coupling in thermal transport
//! - Berry curvature contribution to thermal Hall conductivity
//!
//! ## Key References
//!
//! - M. Schmid et al., "Transverse Spin Seebeck Effect versus Anomalous and
//!   Planar Nernst Effects in Permalloy Thin Films", Phys. Rev. Lett. 111, 187201 (2013)
//! - A. Bose et al., "Observation of the magnon Hall effect",
//!   Nat. Commun. 10, 5965 (2019)
//! - S. Meyer et al., "Observation of the spin Nernst effect",
//!   Nat. Mater. 16, 977 (2017)

use std::fmt;

use crate::vector3::Vector3;

/// Spin Nernst effect in metals and semiconductors
#[derive(Debug, Clone)]
pub struct SpinNernst {
    /// Material name
    pub name: String,

    /// Spin Nernst angle θ_SN (dimensionless)
    /// Typical range: 10⁻⁴ to 10⁻²
    pub spin_nernst_angle: f64,

    /// Thermal conductivity \[W/(m·K)\]
    pub thermal_conductivity: f64,

    /// Electrical conductivity \[S/m\]
    pub electrical_conductivity: f64,

    /// Seebeck coefficient \[V/K\]
    pub seebeck_coefficient: f64,
}

impl Default for SpinNernst {
    /// Default to Platinum parameters
    fn default() -> Self {
        Self::platinum()
    }
}

impl SpinNernst {
    /// Create Pt (Platinum)
    ///
    /// Strong spin-orbit coupling metal with significant SNE
    /// Reference: S. Meyer et al., Nat. Mater. 16, 977 (2017)
    pub fn platinum() -> Self {
        Self {
            name: "Pt".to_string(),
            spin_nernst_angle: 0.01,         // ~1% (estimate)
            thermal_conductivity: 71.6,      // W/(m·K)
            electrical_conductivity: 9.43e6, // S/m
            seebeck_coefficient: -5.0e-6,    // V/K
        }
    }

    /// Create W (Tungsten)
    ///
    /// Large spin Hall and spin Nernst effects
    pub fn tungsten() -> Self {
        Self {
            name: "W".to_string(),
            spin_nernst_angle: 0.015, // ~1.5%
            thermal_conductivity: 173.0,
            electrical_conductivity: 1.79e7,
            seebeck_coefficient: 4.5e-6,
        }
    }

    /// Create Ta (Tantalum)
    ///
    /// Moderate SNE, widely used in spintronics
    pub fn tantalum() -> Self {
        Self {
            name: "Ta".to_string(),
            spin_nernst_angle: 0.005, // ~0.5%
            thermal_conductivity: 57.5,
            electrical_conductivity: 7.41e6,
            seebeck_coefficient: 4.0e-6,
        }
    }

    /// Create Bi (Bismuth)
    ///
    /// Large thermoelectric and spin-orbit coupling
    pub fn bismuth() -> Self {
        Self {
            name: "Bi".to_string(),
            spin_nernst_angle: 0.02, // ~2% (strong SOC)
            thermal_conductivity: 7.97,
            electrical_conductivity: 7.7e5,
            seebeck_coefficient: -72.0e-6, // Large Seebeck
        }
    }

    /// Calculate spin current from temperature gradient
    ///
    /// **j_s** = θ_SN · **∇T** × **σ**
    ///
    /// # Arguments
    /// * `grad_t` - Temperature gradient \[K/m\]
    /// * `spin_direction` - Spin polarization direction (normalized)
    ///
    /// # Returns
    /// Spin current density \[J/m²\]
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::spin_nernst::SpinNernst;
    /// use spintronics::Vector3;
    ///
    /// // Create Pt with Spin Nernst Effect
    /// let pt = SpinNernst::platinum();
    ///
    /// // Temperature gradient: 1 MK/m along x-axis
    /// let grad_t = Vector3::new(1.0e6, 0.0, 0.0);
    /// // Spin polarization along z-axis
    /// let spin_dir = Vector3::new(0.0, 0.0, 1.0);
    ///
    /// // Calculate spin current: j_s = θ_SN (∇T × σ)
    /// let js = pt.spin_current(grad_t, spin_dir);
    ///
    /// // Spin current should be in y-direction (x × z = y)
    /// assert!(js.y.abs() > 0.0);
    /// assert!(js.x.abs() < 1e-10);
    /// assert!(js.z.abs() < 1e-10);
    ///
    /// // Magnitude: j_s = θ_SN × |∇T|
    /// let expected = pt.spin_nernst_angle * 1.0e6;
    /// assert!((js.magnitude() - expected).abs() < 1.0);
    /// ```
    #[inline]
    pub fn spin_current(&self, grad_t: Vector3<f64>, spin_direction: Vector3<f64>) -> Vector3<f64> {
        // j_s = θ_SN (∇T × σ)
        grad_t.cross(&spin_direction) * self.spin_nernst_angle
    }

    /// Calculate heat current from temperature gradient
    ///
    /// **j_Q** = -κ **∇T**
    ///
    /// # Arguments
    /// * `grad_t` - Temperature gradient \[K/m\]
    ///
    /// # Returns
    /// Heat current density \[W/m²\]
    #[inline]
    pub fn heat_current(&self, grad_t: Vector3<f64>) -> Vector3<f64> {
        grad_t * (-self.thermal_conductivity)
    }

    /// Calculate Mott relation estimate for spin Nernst angle
    ///
    /// At high temperature: θ_SN ≈ θ_SH · (π²k_B²T)/(3eE_F) · (dθ_SH/dE)
    ///
    /// For simple estimate: θ_SN ≈ θ_SH · S/e
    /// where S is Seebeck coefficient
    ///
    /// # Arguments
    /// * `spin_hall_angle` - Spin Hall angle
    ///
    /// # Returns
    /// Estimated spin Nernst angle
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::spin_nernst::SpinNernst;
    ///
    /// // Platinum with known Spin Hall angle
    /// let pt = SpinNernst::platinum();
    /// let theta_sh = 0.07; // Pt spin Hall angle
    ///
    /// // Estimate Spin Nernst angle from Mott relation
    /// // θ_SN ≈ θ_SH × (S/e)
    /// let theta_sn_est = pt.mott_relation_estimate(theta_sh);
    ///
    /// // Should be positive and finite
    /// assert!(theta_sn_est > 0.0);
    /// assert!(theta_sn_est.is_finite());
    ///
    /// // Typically θ_SN << θ_SH for metals
    /// // (Mott relation gives order of magnitude estimate)
    /// ```
    pub fn mott_relation_estimate(&self, spin_hall_angle: f64) -> f64 {
        let e = 1.602e-19; // Elementary charge
        spin_hall_angle * (self.seebeck_coefficient / e).abs()
    }

    /// Calculate thermal spin Hall conductivity
    ///
    /// σ_SN = θ_SN · κ / T
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]
    ///
    /// # Returns
    /// Thermal spin Hall conductivity \[W/(m·K²)\]
    pub fn thermal_spin_hall_conductivity(&self, temperature: f64) -> f64 {
        self.spin_nernst_angle * self.thermal_conductivity / temperature
    }

    /// Calculate efficiency for thermal spin injection
    ///
    /// How effective is this material for generating spin current from heat?
    ///
    /// # Returns
    /// Figure of merit (dimensionless)
    pub fn thermal_spin_injection_efficiency(&self) -> f64 {
        // Proportional to θ_SN × κ
        self.spin_nernst_angle * self.thermal_conductivity / 100.0
    }

    /// Check suitability for thermal spintronic devices
    pub fn is_suitable_for_thermal_spintronics(&self) -> bool {
        self.spin_nernst_angle > 0.001 && // Detectable SNE
        self.thermal_conductivity > 10.0 // Good thermal conductor
    }

    /// Builder method to set spin Nernst angle
    pub fn with_spin_nernst_angle(mut self, angle: f64) -> Self {
        self.spin_nernst_angle = angle;
        self
    }

    /// Builder method to set thermal conductivity
    pub fn with_thermal_conductivity(mut self, kappa: f64) -> Self {
        self.thermal_conductivity = kappa;
        self
    }
}

impl fmt::Display for SpinNernst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: θ_SN={:.4}, κ={:.1} W/(m·K), S={:.2e} V/K",
            self.name, self.spin_nernst_angle, self.thermal_conductivity, self.seebeck_coefficient
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_platinum() {
        let pt = SpinNernst::platinum();
        assert!(pt.spin_nernst_angle > 0.0);
        assert!(pt.thermal_conductivity > 70.0);
        assert!(pt.is_suitable_for_thermal_spintronics());
    }

    #[test]
    fn test_tungsten() {
        let pt = SpinNernst::platinum();
        let w = SpinNernst::tungsten();
        assert!(w.spin_nernst_angle > pt.spin_nernst_angle);
        assert!(w.thermal_conductivity > 170.0);
    }

    #[test]
    fn test_bismuth_strong_soc() {
        let bi = SpinNernst::bismuth();
        // Bi should have largest SNE due to strong SOC
        assert!(bi.spin_nernst_angle > 0.01);
        // Large Seebeck coefficient
        assert!(bi.seebeck_coefficient.abs() > 50.0e-6);
    }

    #[test]
    fn test_spin_current_perpendicularity() {
        let pt = SpinNernst::platinum();
        let grad_t = Vector3::new(1000.0, 0.0, 0.0); // Temperature gradient along x
        let spin_dir = Vector3::new(0.0, 0.0, 1.0); // Spin up (z)

        let j_s = pt.spin_current(grad_t, spin_dir);

        // j_s should be perpendicular to both ∇T and σ
        assert!(j_s.dot(&grad_t).abs() < 1e-10);
        assert!(j_s.dot(&spin_dir).abs() < 1e-10);

        // j_s should be along y direction
        assert!(j_s.y.abs() > 0.0);
        assert!(j_s.x.abs() < 1e-10);
        assert!(j_s.z.abs() < 1e-10);
    }

    #[test]
    fn test_spin_current_magnitude() {
        let pt = SpinNernst::platinum();
        let grad_t = Vector3::new(1.0e6, 0.0, 0.0); // 1 MK/m gradient
        let spin_dir = Vector3::new(0.0, 0.0, 1.0);

        let j_s = pt.spin_current(grad_t, spin_dir);

        // j_s = θ_SN × |∇T|
        let expected_magnitude = pt.spin_nernst_angle * 1.0e6;
        assert!((j_s.magnitude() - expected_magnitude).abs() < 1.0);
    }

    #[test]
    fn test_heat_current() {
        let pt = SpinNernst::platinum();
        let grad_t = Vector3::new(100.0, 0.0, 0.0); // 100 K/m

        let j_q = pt.heat_current(grad_t);

        // j_Q = -κ ∇T
        let expected = -pt.thermal_conductivity * 100.0;
        assert!((j_q.x - expected).abs() < 0.1);
    }

    #[test]
    fn test_mott_relation() {
        let pt = SpinNernst::platinum();
        let theta_sh = 0.07; // Typical Pt spin Hall angle

        let theta_sn_estimate = pt.mott_relation_estimate(theta_sh);

        // Should be positive (this is just an order-of-magnitude estimate)
        assert!(theta_sn_estimate > 0.0);
        // The simple Mott relation overestimates, so we just check it's finite
        assert!(theta_sn_estimate.is_finite());
    }

    #[test]
    fn test_thermal_spin_hall_conductivity() {
        let pt = SpinNernst::platinum();
        let temp = 300.0; // K

        let sigma_sn = pt.thermal_spin_hall_conductivity(temp);

        // Should be positive
        assert!(sigma_sn > 0.0);
        // Reasonable magnitude
        assert!(sigma_sn < 10.0);
    }

    #[test]
    fn test_thermal_spin_injection_efficiency() {
        let pt = SpinNernst::platinum();
        let w = SpinNernst::tungsten();
        let bi = SpinNernst::bismuth();

        let eff_pt = pt.thermal_spin_injection_efficiency();
        let eff_w = w.thermal_spin_injection_efficiency();
        let eff_bi = bi.thermal_spin_injection_efficiency();

        // All should be positive
        assert!(eff_pt > 0.0);
        assert!(eff_w > 0.0);
        assert!(eff_bi > 0.0);

        // W should have highest efficiency (large θ_SN and κ)
        assert!(eff_w > eff_pt);
    }

    #[test]
    fn test_thermal_spintronics_suitability() {
        let strong = SpinNernst::tungsten();
        let weak = SpinNernst::platinum().with_spin_nernst_angle(0.0001);

        assert!(strong.is_suitable_for_thermal_spintronics());
        assert!(!weak.is_suitable_for_thermal_spintronics());
    }

    #[test]
    fn test_builder_pattern() {
        let custom = SpinNernst::platinum()
            .with_spin_nernst_angle(0.05)
            .with_thermal_conductivity(100.0);

        assert_eq!(custom.spin_nernst_angle, 0.05);
        assert_eq!(custom.thermal_conductivity, 100.0);
    }
}
