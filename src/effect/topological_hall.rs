//! Topological Hall Effect
//!
//! The Topological Hall Effect (THE) is an emergent electromagnetic phenomenon
//! arising from non-coplanar spin textures in magnetic materials. Conduction
//! electrons moving through these textures accumulate a Berry phase, leading to
//! a transverse Hall voltage even in the absence of an external magnetic field.
//!
//! ## Physics Background
//!
//! ### Berry Phase Origin:
//! Electrons traversing a skyrmion texture experience an effective magnetic field:
//!
//! **B_eff** = (Φ₀/2) · **m** · (∂_x **m** × ∂_y **m**)
//!
//! where Φ₀ = h/e is the flux quantum and **m** is the unit magnetization.
//!
//! ### Topological Hall Resistivity:
//! ρ_THE = R₀ · n_sk · Q
//!
//! where:
//! - R₀: Material-dependent coefficient [Ω·cm]
//! - n_sk: Skyrmion density [m⁻²]
//! - Q: Topological charge (winding number, typically ±1)
//!
//! ### Relationship to Skyrmions:
//! - Each skyrmion contributes one quantum of emergent flux
//! - THE is proportional to skyrmion density
//! - Sign of THE indicates skyrmion chirality
//!
//! ## Key References
//!
//! - A. Neubauer et al., "Topological Hall Effect in the A Phase of MnSi",
//!   Phys. Rev. Lett. 102, 186602 (2009)
//! - N. Kanazawa et al., "Large Topological Hall Effect in a Short-Period
//!   Helimagnet MnGe", Phys. Rev. Lett. 106, 156603 (2011)
//! - S. Huang et al., "Topological Hall effect in thin films of the
//!   Heusler compound Co₂MnGa", Phys. Rev. B 95, 075133 (2017)

use std::f64::consts::PI;

/// Topological Hall effect in skyrmion-hosting materials
#[derive(Debug, Clone)]
pub struct TopologicalHall {
    /// Material name
    pub name: String,

    /// Normal Hall coefficient R₀ [Ω·cm/T]
    pub hall_coefficient: f64,

    /// Anomalous Hall coefficient R_s [Ω·cm]
    pub anomalous_hall_coeff: f64,

    /// Resistivity [Ω·cm]
    pub resistivity: f64,

    /// Saturation magnetization \[A/m\]
    pub magnetization: f64,

    /// Typical skyrmion size \[nm\]
    pub skyrmion_diameter: f64,
}

impl TopologicalHall {
    /// Create MnSi (Manganese Silicide)
    ///
    /// First material where THE was clearly identified in the skyrmion A-phase.
    ///
    /// Reference: A. Neubauer et al., PRL 102, 186602 (2009)
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::topological_hall::TopologicalHall;
    ///
    /// // Create MnSi with skyrmion lattice phase
    /// let mnsi = TopologicalHall::mnsi();
    ///
    /// assert_eq!(mnsi.name, "MnSi");
    /// assert!(mnsi.skyrmion_diameter > 10.0); // ~18 nm skyrmions
    ///
    /// // Typical skyrmion density in MnSi: ~10^14 m^-2
    /// let n_sk = 1.0e14; // Skyrmion density [m^-2]
    /// let q = -1.0;      // Topological charge
    ///
    /// // Calculate topological Hall resistivity
    /// let rho_the = mnsi.topological_hall_resistivity(n_sk, q);
    ///
    /// // Should be measurable (~nΩ·cm range)
    /// assert!(rho_the.abs() > 1e-10);
    /// assert!(mnsi.is_the_measurable(n_sk));
    /// ```
    pub fn mnsi() -> Self {
        Self {
            name: "MnSi".to_string(),
            hall_coefficient: 1.5e-10,    // Ω·cm/T
            anomalous_hall_coeff: 2.0e-8, // Ω·cm
            resistivity: 1.0e-5,          // Ω·cm
            magnetization: 1.8e5,         // A/m
            skyrmion_diameter: 18.0,      // nm
        }
    }

    /// Create MnGe (Manganese Germanium)
    ///
    /// Large THE in short-period helimagnet
    /// Reference: N. Kanazawa et al., PRL 106, 156603 (2011)
    pub fn mnge() -> Self {
        Self {
            name: "MnGe".to_string(),
            hall_coefficient: 2.0e-10,
            anomalous_hall_coeff: 3.0e-8,
            resistivity: 2.0e-5,
            magnetization: 1.5e5,
            skyrmion_diameter: 3.0, // nm (very small!)
        }
    }

    /// Create FeGe (Iron Germanium)
    ///
    /// B20 chiral magnet with stable skyrmion phase
    pub fn fege() -> Self {
        Self {
            name: "FeGe".to_string(),
            hall_coefficient: 1.8e-10,
            anomalous_hall_coeff: 2.5e-8,
            resistivity: 1.2e-5,
            magnetization: 3.84e5,   // A/m
            skyrmion_diameter: 70.0, // nm
        }
    }

    /// Create Co₂MnGa Heusler compound
    ///
    /// THE in thin films
    /// Reference: S. Huang et al., PRB 95, 075133 (2017)
    pub fn co2mnga() -> Self {
        Self {
            name: "Co₂MnGa".to_string(),
            hall_coefficient: 1.0e-10,
            anomalous_hall_coeff: 1.5e-8,
            resistivity: 5.0e-6,
            magnetization: 6.0e5,
            skyrmion_diameter: 50.0,
        }
    }

    /// Calculate topological Hall resistivity
    ///
    /// The topological Hall resistivity arises from the Berry phase that
    /// conduction electrons acquire when traversing skyrmion textures:
    ///
    /// $$\rho_{\text{THE}} = \alpha \cdot n_{\text{sk}} \cdot Q$$
    ///
    /// where:
    /// - $\alpha$ is a material-dependent coefficient [Ω·cm]
    /// - $n_{\text{sk}}$ is the skyrmion density [m⁻²]
    /// - $Q$ is the topological charge (winding number, typically ±1)
    ///
    /// # Arguments
    /// * `skyrmion_density` - Skyrmion density [m⁻²]
    /// * `topological_charge` - Q (typically ±1)
    ///
    /// # Returns
    /// Topological Hall resistivity [Ω·cm]
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::topological_hall::TopologicalHall;
    ///
    /// let mnsi = TopologicalHall::mnsi();
    ///
    /// // Skyrmion lattice density from Neubauer 2009
    /// // Lattice constant a ≈ 18 nm → density ≈ 1/(a²)
    /// let lattice_constant = 18e-9; // m
    /// let n_sk = 1.0 / (lattice_constant * lattice_constant); // ~3×10^15 m^-2
    ///
    /// // Each skyrmion has Q = -1 (negative chirality)
    /// let q = -1.0;
    ///
    /// // Calculate THE resistivity
    /// let rho_the = mnsi.topological_hall_resistivity(n_sk, q);
    ///
    /// // THE resistivity should be non-zero and finite
    /// assert!(rho_the.abs() > 0.0);
    /// assert!(rho_the.is_finite());
    ///
    /// println!("MnSi THE resistivity: {:.2e} Ω·cm", rho_the);
    /// ```
    pub fn topological_hall_resistivity(
        &self,
        skyrmion_density: f64,
        topological_charge: f64,
    ) -> f64 {
        // Empirical: ρ_THE ∝ R₀ × n_sk × Q
        // With proper unit conversion: m^-2 to cm^-2
        let alpha = self.hall_coefficient * 1.0e-4; // Ω·cm/T × cm²/m²
        alpha * skyrmion_density * topological_charge
    }

    /// Calculate topological Hall angle
    ///
    /// tan(θ_THE) = ρ_THE / ρ_xx
    ///
    /// # Arguments
    /// * `skyrmion_density` - Skyrmion density [m⁻²]
    /// * `topological_charge` - Q (typically ±1)
    ///
    /// # Returns
    /// Topological Hall angle \[rad\]
    pub fn topological_hall_angle(&self, skyrmion_density: f64, topological_charge: f64) -> f64 {
        let rho_the = self.topological_hall_resistivity(skyrmion_density, topological_charge);
        (rho_the / self.resistivity).atan()
    }

    /// Calculate emergent magnetic field from skyrmion texture
    ///
    /// Each skyrmion creates an emergent magnetic flux quantum:
    ///
    /// $$B_{\text{eff}} = \frac{\Phi_0 \cdot |Q|}{\pi r^2}$$
    ///
    /// where:
    /// - $\Phi_0 = h/e$ is the magnetic flux quantum (4.136×10⁻¹⁵ V·s)
    /// - $Q$ is the topological charge
    /// - $r$ is the skyrmion radius
    ///
    /// This effective field is what electrons "see" when moving through
    /// the skyrmion texture due to Berry phase effects.
    ///
    /// # Arguments
    /// * `topological_charge` - Skyrmion winding number
    ///
    /// # Returns
    /// Effective magnetic field \[T\]
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::topological_hall::TopologicalHall;
    ///
    /// // MnSi with 18 nm diameter skyrmions
    /// let mnsi = TopologicalHall::mnsi();
    ///
    /// // Calculate emergent field from single skyrmion
    /// let b_eff = mnsi.emergent_magnetic_field(-1.0); // Q = -1
    ///
    /// // For 18 nm skyrmion: B_eff ~ 10-100 mT range
    /// assert!(b_eff > 0.0);
    /// assert!(b_eff < 1.0); // Should be sub-Tesla
    ///
    /// println!("Emergent B-field: {:.1} mT", b_eff * 1e3);
    ///
    /// // Smaller skyrmions → stronger emergent field
    /// let mnge = TopologicalHall::mnge();
    /// let b_eff_mnge = mnge.emergent_magnetic_field(-1.0);
    ///
    /// // MnGe has 3 nm skyrmions → much stronger field
    /// assert!(b_eff_mnge > b_eff);
    /// ```
    pub fn emergent_magnetic_field(&self, topological_charge: f64) -> f64 {
        let phi_0 = 4.136e-15; // h/e in V·s
        let radius = self.skyrmion_diameter * 0.5 * 1e-9; // m
        let area = PI * radius * radius;

        (phi_0 * topological_charge.abs() / area) * 1e-4 // Convert to Tesla
    }

    /// Estimate skyrmion density from Hall resistivity measurement
    ///
    /// # Arguments
    /// * `measured_rho_the` - Measured topological Hall resistivity [Ω·cm]
    /// * `assumed_charge` - Assumed topological charge (±1)
    ///
    /// # Returns
    /// Estimated skyrmion density [m⁻²]
    pub fn estimate_skyrmion_density(&self, measured_rho_the: f64, assumed_charge: f64) -> f64 {
        let alpha = self.hall_coefficient * 1.0e-4;
        measured_rho_the / (alpha * assumed_charge)
    }

    /// Calculate Hall voltage from skyrmion lattice
    ///
    /// The Hall voltage in a skyrmion-hosting sample:
    ///
    /// $$V_H = \frac{\rho_{\text{THE}}}{t} \cdot I$$
    ///
    /// where $t$ is the sample thickness and $I$ is the applied current.
    ///
    /// # Arguments
    /// * `current` - Applied current \[A\]
    /// * `thickness` - Sample thickness \[nm\]
    /// * `skyrmion_density` - Skyrmion density [m⁻²]
    /// * `topological_charge` - Q
    ///
    /// # Returns
    /// Hall voltage \[V\]
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::topological_hall::TopologicalHall;
    ///
    /// let fege = TopologicalHall::fege();
    ///
    /// // Experimental conditions for FeGe thin film
    /// let current = 1.0e-6;      // 1 μA applied current
    /// let thickness = 100.0;     // 100 nm film
    /// let n_sk = 5.0e14;         // Skyrmion density [m^-2]
    /// let q = -1.0;              // Topological charge
    ///
    /// // Calculate THE voltage
    /// let v_hall = fege.hall_voltage(current, thickness, n_sk, q);
    ///
    /// // Topological Hall voltage should be non-zero and finite
    /// assert!(v_hall.abs() > 0.0);
    /// assert!(v_hall.is_finite());
    ///
    /// println!("THE voltage: {:.3e} V", v_hall);
    /// ```
    pub fn hall_voltage(
        &self,
        current: f64,
        thickness: f64,
        skyrmion_density: f64,
        topological_charge: f64,
    ) -> f64 {
        let rho_the = self.topological_hall_resistivity(skyrmion_density, topological_charge);
        let thickness_cm = thickness * 1e-7; // nm to cm
        (rho_the / thickness_cm) * current
    }

    /// Check if THE is measurable
    pub fn is_the_measurable(&self, skyrmion_density: f64) -> bool {
        let rho_the = self.topological_hall_resistivity(skyrmion_density, 1.0);
        // THE should be at least 1% of background resistivity
        rho_the / self.resistivity > 0.01
    }

    /// Builder method to set skyrmion diameter
    pub fn with_skyrmion_diameter(mut self, diameter: f64) -> Self {
        self.skyrmion_diameter = diameter;
        self
    }

    /// Builder method to set Hall coefficient
    pub fn with_hall_coefficient(mut self, r0: f64) -> Self {
        self.hall_coefficient = r0;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mnsi() {
        let mnsi = TopologicalHall::mnsi();
        assert_eq!(mnsi.name, "MnSi");
        assert!(mnsi.skyrmion_diameter > 15.0);
        assert!(mnsi.skyrmion_diameter < 20.0);
    }

    #[test]
    fn test_mnge_small_skyrmions() {
        let mnge = TopologicalHall::mnge();
        // MnGe has the smallest skyrmions
        assert!(mnge.skyrmion_diameter < 5.0);
    }

    #[test]
    fn test_topological_hall_resistivity() {
        let mnsi = TopologicalHall::mnsi();
        let n_sk = 1.0e14; // skyrmions/m²
        let q = 1.0; // Topological charge

        let rho_the = mnsi.topological_hall_resistivity(n_sk, q);

        // Should be positive and reasonable (THE is a small effect)
        assert!(rho_the > 0.0);
        assert!(rho_the.is_finite());
    }

    #[test]
    fn test_topological_hall_angle() {
        let mnsi = TopologicalHall::mnsi();
        let n_sk = 5.0e14; // high density

        let theta = mnsi.topological_hall_angle(n_sk, 1.0);

        // Should be positive (typically a small angle)
        assert!(theta > 0.0);
        assert!(theta.is_finite());
    }

    #[test]
    fn test_emergent_magnetic_field() {
        let mnsi = TopologicalHall::mnsi();
        let b_eff = mnsi.emergent_magnetic_field(1.0);

        // Should be significant (mT range)
        assert!(b_eff > 0.0);
        assert!(b_eff.is_finite());
    }

    #[test]
    fn test_estimate_skyrmion_density() {
        let mnsi = TopologicalHall::mnsi();
        let n_sk_actual = 1.0e14; // m⁻²

        // Calculate THE for this density
        let rho_the = mnsi.topological_hall_resistivity(n_sk_actual, 1.0);

        // Estimate density from resistivity
        let n_sk_estimated = mnsi.estimate_skyrmion_density(rho_the, 1.0);

        // Should recover the original density
        assert!((n_sk_estimated - n_sk_actual).abs() / n_sk_actual < 0.01);
    }

    #[test]
    fn test_hall_voltage() {
        let mnsi = TopologicalHall::mnsi();
        let current = 1.0e-6; // 1 μA
        let thickness = 100.0; // nm
        let n_sk = 1.0e14;

        let v_h = mnsi.hall_voltage(current, thickness, n_sk, 1.0);

        // Should be measurable (nV to μV range)
        assert!(v_h > 0.0);
        assert!(v_h.is_finite());
    }

    #[test]
    fn test_chirality_sign() {
        let mnsi = TopologicalHall::mnsi();
        let n_sk = 1.0e14;

        let rho_plus = mnsi.topological_hall_resistivity(n_sk, 1.0);
        let rho_minus = mnsi.topological_hall_resistivity(n_sk, -1.0);

        // Opposite chirality should give opposite sign
        assert!((rho_plus + rho_minus).abs() < 1e-15);
    }

    #[test]
    fn test_the_measurability() {
        let mnsi = TopologicalHall::mnsi();

        let high_density = 1.0e15; // Lots of skyrmions
        let low_density = 1.0e6; // Very few

        assert!(mnsi.is_the_measurable(high_density));
        assert!(!mnsi.is_the_measurable(low_density));
    }

    #[test]
    fn test_density_dependence() {
        let mnsi = TopologicalHall::mnsi();
        let n_sk1 = 1.0e14;
        let n_sk2 = 2.0e14;

        let rho1 = mnsi.topological_hall_resistivity(n_sk1, 1.0);
        let rho2 = mnsi.topological_hall_resistivity(n_sk2, 1.0);

        // Should scale linearly with density
        assert!((rho2 / rho1 - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_builder_pattern() {
        let custom = TopologicalHall::mnsi()
            .with_skyrmion_diameter(100.0)
            .with_hall_coefficient(5.0e-10);

        assert_eq!(custom.skyrmion_diameter, 100.0);
        assert_eq!(custom.hall_coefficient, 5.0e-10);
    }
}
