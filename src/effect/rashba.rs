//! Rashba Effect in 2D Electron Gas
//!
//! The Rashba effect is a momentum-dependent spin splitting of electronic bands
//! in systems lacking inversion symmetry, particularly at interfaces and surfaces.
//! It's a key mechanism for electrical control of spin in spintronics devices.
//!
//! ## Physics Background
//!
//! ### Rashba Hamiltonian:
//! H_R = α_R (σ × k) · ẑ
//!
//! where:
//! - α_R: Rashba parameter [eV·Å]
//! - σ: Pauli spin matrices
//! - k: in-plane momentum
//! - ẑ: surface normal
//!
//! ### Spin-Momentum Locking:
//! The Rashba effect locks the electron spin perpendicular to its momentum,
//! enabling conversion between charge and spin currents.
//!
//! ### Applications:
//! - **Spin-orbit torques**: Current-induced magnetization switching
//! - **Spin Hall effect**: Charge current → transverse spin current
//! - **Edelstein effect**: Charge accumulation → spin polarization
//! - **Datta-Das spin FET**: Gate-tunable spin precession
//!
//! ## Key References
//!
//! - E. I. Rashba, "Properties of semiconductors with an extremum loop",
//!   Sov. Phys. Solid State 2, 1109 (1960)
//! - A. Manchon et al., "New perspectives for Rashba spin-orbit coupling",
//!   Nat. Mater. 14, 871 (2015)
//! - J. C. R. Sánchez et al., "Spin-to-charge conversion using Rashba coupling",
//!   Nat. Commun. 4, 2944 (2013)

use crate::constants::HBAR;
use crate::vector3::Vector3;

/// Rashba spin-orbit coupling interface/material
#[derive(Debug, Clone)]
pub struct RashbaSystem {
    /// Material/interface name
    pub name: String,

    /// Rashba parameter α_R [eV·Å]
    /// Typical values: 0.1-10 eV·Å
    pub alpha_r: f64,

    /// Fermi energy \[eV\]
    pub fermi_energy: f64,

    /// Effective mass (in units of electron mass m_e)
    pub effective_mass: f64,

    /// Carrier density [m⁻²]
    pub carrier_density: f64,

    /// Surface normal direction
    pub surface_normal: Vector3<f64>,
}

impl RashbaSystem {
    /// Create Ag/Bi interface
    ///
    /// Giant Rashba splitting at metal-Bi interface
    /// Reference: C. R. Ast et al., Phys. Rev. Lett. 98, 186807 (2007)
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::rashba::RashbaSystem;
    ///
    /// // Create Ag/Bi interface with giant Rashba coupling
    /// let system = RashbaSystem::ag_bi();
    ///
    /// // Ag/Bi has one of the largest known Rashba parameters
    /// assert!(system.alpha_r > 3.0);  // eV·Å
    /// assert_eq!(system.name, "Ag/Bi");
    ///
    /// // Check it's suitable for spintronics
    /// assert!(system.is_suitable_for_spintronics());
    /// ```
    pub fn ag_bi() -> Self {
        Self {
            name: "Ag/Bi".to_string(),
            alpha_r: 3.05, // eV·Å (giant Rashba)
            fermi_energy: 0.1,
            effective_mass: 0.15,
            carrier_density: 5.0e17, // m⁻²
            surface_normal: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Create Au(111) surface
    ///
    /// Strong surface state with Rashba splitting
    /// Reference: S. LaShell et al., Phys. Rev. Lett. 77, 3419 (1996)
    pub fn au_111() -> Self {
        Self {
            name: "Au(111)".to_string(),
            alpha_r: 0.33, // eV·Å
            fermi_energy: 0.4,
            effective_mass: 0.25,
            carrier_density: 1.0e18,
            surface_normal: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Create BiTeI crystal
    ///
    /// Bulk Rashba semiconductor
    /// Reference: K. Ishizaka et al., Nat. Mater. 10, 521 (2011)
    pub fn bitei() -> Self {
        Self {
            name: "BiTeI".to_string(),
            alpha_r: 3.8, // eV·Å (largest bulk Rashba)
            fermi_energy: 0.05,
            effective_mass: 0.15,
            carrier_density: 3.0e17,
            surface_normal: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Create GaAs/AlGaAs 2DEG
    ///
    /// Tunable Rashba via gate voltage
    /// Reference: J. Nitta et al., Phys. Rev. Lett. 78, 1335 (1997)
    pub fn gaas_2deg(gate_voltage: f64) -> Self {
        // α_R ∝ electric field ∝ gate voltage
        let alpha_r = 0.07 * (1.0 + gate_voltage / 1.0); // eV·Å

        Self {
            name: format!("GaAs/AlGaAs (V_g = {:.2} V)", gate_voltage),
            alpha_r,
            fermi_energy: 0.02,
            effective_mass: 0.067,
            carrier_density: 1.0e16,
            surface_normal: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Create LAO/STO interface
    ///
    /// Oxide 2DEG with strong Rashba coupling
    /// Reference: A. D. Caviglia et al., Phys. Rev. Lett. 104, 126803 (2010)
    pub fn lao_sto() -> Self {
        Self {
            name: "LaAlO₃/SrTiO₃".to_string(),
            alpha_r: 0.1, // eV·Å (tunable)
            fermi_energy: 0.03,
            effective_mass: 1.0, // Heavy 3d electrons
            carrier_density: 5.0e16,
            surface_normal: Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Calculate spin splitting energy at given momentum
    ///
    /// The Rashba Hamiltonian splits spin-degenerate bands:
    ///
    /// $$\Delta E_R = 2\alpha_R |k|$$
    ///
    /// where:
    /// - $\alpha_R$ is the Rashba parameter [eV·Å]
    /// - $k$ is the in-plane momentum [Å⁻¹]
    ///
    /// # Arguments
    /// * `k_momentum` - In-plane momentum vector [Å⁻¹]
    ///
    /// # Returns
    /// Spin splitting energy \[eV\]
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::rashba::RashbaSystem;
    /// use spintronics::Vector3;
    ///
    /// // Au(111) surface state
    /// let system = RashbaSystem::au_111();
    ///
    /// // Momentum: k = 0.1 Å⁻¹ in x-direction
    /// let k = Vector3::new(0.1, 0.0, 0.0);
    ///
    /// // Calculate spin splitting
    /// let splitting = system.spin_splitting(k);
    ///
    /// // ΔE = 2 * α_R * |k| = 2 * 0.33 * 0.1 = 0.066 eV
    /// let expected = 2.0 * 0.33 * 0.1;
    /// assert!((splitting - expected).abs() < 1e-10);
    /// ```
    pub fn spin_splitting(&self, k_momentum: Vector3<f64>) -> f64 {
        // Project momentum onto plane perpendicular to surface normal
        let k_inplane = k_momentum - self.surface_normal * k_momentum.dot(&self.surface_normal);
        2.0 * self.alpha_r * k_inplane.magnitude()
    }

    /// Calculate spin texture (momentum-dependent spin direction)
    ///
    /// The Rashba effect creates spin-momentum locking:
    ///
    /// $$\mathbf{s} \propto \hat{z} \times \mathbf{k}$$
    ///
    /// For momentum **k**, the spin is perpendicular to both the surface
    /// normal and the momentum direction, forming a hedgehog texture in
    /// momentum space.
    ///
    /// # Arguments
    /// * `k_momentum` - In-plane momentum [Å⁻¹]
    ///
    /// # Returns
    /// Spin polarization direction (normalized)
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::rashba::RashbaSystem;
    /// use spintronics::Vector3;
    ///
    /// let system = RashbaSystem::bitei();
    ///
    /// // Electron moving in +x direction
    /// let k_x = Vector3::new(0.1, 0.0, 0.0);
    /// let spin_x = system.spin_texture(k_x);
    ///
    /// // Spin should point in +y direction (ẑ × x̂ = ŷ)
    /// assert!(spin_x.y > 0.9);
    /// assert!(spin_x.x.abs() < 1e-10);
    ///
    /// // Electron moving in +y direction
    /// let k_y = Vector3::new(0.0, 0.1, 0.0);
    /// let spin_y = system.spin_texture(k_y);
    ///
    /// // Spin should point in -x direction (ẑ × ŷ = -x̂)
    /// assert!(spin_y.x < -0.9);
    /// assert!(spin_y.y.abs() < 1e-10);
    /// ```
    pub fn spin_texture(&self, k_momentum: Vector3<f64>) -> Vector3<f64> {
        let k_inplane = k_momentum - self.surface_normal * k_momentum.dot(&self.surface_normal);
        self.surface_normal.cross(&k_inplane).normalize()
    }

    /// Calculate Edelstein effect (charge-to-spin conversion)
    ///
    /// The Edelstein (inverse Rashba-Edelstein) effect converts charge
    /// current into non-equilibrium spin polarization:
    ///
    /// $$\mathbf{s} = \lambda_E \cdot \mathbf{j}_c$$
    ///
    /// where:
    /// - $\lambda_E = \alpha_R / (2E_F)$ is the Edelstein length
    /// - $\mathbf{j}_c$ is the charge current density
    ///
    /// The spin direction is perpendicular to the current flow.
    ///
    /// # Arguments
    /// * `current_density` - Charge current density \[A/m²\]
    /// * `current_direction` - Direction of current flow
    ///
    /// # Returns
    /// Spin density [J·s/m³]
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::rashba::RashbaSystem;
    /// use spintronics::Vector3;
    ///
    /// // LAO/STO interface
    /// let system = RashbaSystem::lao_sto();
    ///
    /// // Apply current in x-direction: 1 MA/m²
    /// let j_c = 1.0e6;  // A/m²
    /// let current_dir = Vector3::new(1.0, 0.0, 0.0);
    ///
    /// // Calculate induced spin density
    /// let spin_density = system.edelstein_spin_density(j_c, current_dir);
    ///
    /// // Spin should be perpendicular to current (in y-direction)
    /// assert!(spin_density.y.abs() > 0.0);
    /// assert!(spin_density.x.abs() < 1e-20);
    ///
    /// // Spin density magnitude proportional to current
    /// assert!(spin_density.magnitude() > 0.0);
    /// ```
    pub fn edelstein_spin_density(
        &self,
        current_density: f64,
        current_direction: Vector3<f64>,
    ) -> Vector3<f64> {
        // Edelstein length
        let lambda_e = self.alpha_r / (2.0 * self.fermi_energy); // Å
        let lambda_e_si = lambda_e * 1e-10; // m

        // Spin direction: perpendicular to current
        let spin_direction = self.surface_normal.cross(&current_direction.normalize());

        spin_direction * (lambda_e_si * current_density)
    }

    /// Calculate inverse Edelstein effect (spin-to-charge conversion)
    ///
    /// Charge current: j_c = (2e/ℏ) · α_R · (n × s)
    ///
    /// # Arguments
    /// * `spin_density` - Non-equilibrium spin density [J·s/m³]
    ///
    /// # Returns
    /// Charge current density \[A/m²\]
    pub fn inverse_edelstein_current(&self, spin_density: Vector3<f64>) -> Vector3<f64> {
        let e = 1.602e-19; // Elementary charge \[C\]
        let alpha_r_si = self.alpha_r * 1.602e-19 * 1e-10; // Convert eV·Å to J·m

        // j_c ∝ n × s
        let current_direction = self.surface_normal.cross(&spin_density);

        current_direction * (2.0 * e / HBAR * alpha_r_si)
    }

    /// Calculate Fermi momentum
    ///
    /// k_F = (2π n_s)^(1/2) for 2D system
    ///
    /// # Returns
    /// Fermi momentum [Å⁻¹]
    pub fn fermi_momentum(&self) -> f64 {
        let k_f_si = (2.0 * std::f64::consts::PI * self.carrier_density).sqrt(); // m⁻¹
        k_f_si * 1e-10 // Å⁻¹
    }

    /// Calculate spin precession length (Datta-Das)
    ///
    /// The spin precession length determines the distance over which
    /// an electron's spin rotates by 2π due to Rashba spin-orbit coupling:
    ///
    /// $$L_{SO} = \frac{\pi \hbar^2}{m^* \alpha_R}$$
    ///
    /// This is the fundamental length scale for the Datta-Das spin FET.
    ///
    /// # Returns
    /// Spin precession length \[nm\]
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::rashba::RashbaSystem;
    ///
    /// // GaAs 2DEG with gate voltage
    /// let system = RashbaSystem::gaas_2deg(1.0);  // 1V gate
    ///
    /// // Calculate spin precession length
    /// let l_so = system.spin_precession_length();
    ///
    /// // For GaAs, typically 100-1000 nm range
    /// assert!(l_so > 10.0);   // nm
    /// assert!(l_so < 10000.0); // nm
    /// ```
    pub fn spin_precession_length(&self) -> f64 {
        let m_e = 9.109e-31; // kg
        let m_eff = self.effective_mass * m_e;
        let alpha_r_si = self.alpha_r * 1.602e-19 * 1e-10; // J·m

        let l_so = std::f64::consts::PI * HBAR * HBAR / (m_eff * alpha_r_si);
        l_so * 1e9 // nm
    }

    /// Check if suitable for spin-orbitronics
    pub fn is_suitable_for_spintronics(&self) -> bool {
        self.alpha_r > 0.1 && // Strong Rashba coupling
        self.fermi_energy > 0.01 // Metallic regime
    }

    /// Builder method to set Rashba parameter
    pub fn with_alpha_r(mut self, alpha_r: f64) -> Self {
        self.alpha_r = alpha_r;
        self
    }

    /// Builder method to set carrier density
    pub fn with_carrier_density(mut self, density: f64) -> Self {
        self.carrier_density = density;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ag_bi() {
        let rashba = RashbaSystem::ag_bi();
        assert!(rashba.alpha_r > 3.0);
        assert!(rashba.is_suitable_for_spintronics());
    }

    #[test]
    fn test_au_111() {
        let rashba = RashbaSystem::au_111();
        assert!(rashba.alpha_r > 0.3);
        assert!(rashba.alpha_r < 0.4);
    }

    #[test]
    fn test_bitei() {
        let bitei = RashbaSystem::bitei();
        assert!(bitei.alpha_r > 3.5); // Largest bulk Rashba
        assert!(bitei.is_suitable_for_spintronics());
    }

    #[test]
    fn test_gaas_gate_tuning() {
        let gaas_0v = RashbaSystem::gaas_2deg(0.0);
        let gaas_1v = RashbaSystem::gaas_2deg(1.0);

        // Gate voltage should increase Rashba parameter
        assert!(gaas_1v.alpha_r > gaas_0v.alpha_r);
    }

    #[test]
    fn test_spin_splitting() {
        let rashba = RashbaSystem::ag_bi();
        let k = Vector3::new(0.1, 0.0, 0.0); // k_x = 0.1 Å⁻¹

        let splitting = rashba.spin_splitting(k);

        // ΔE = 2 α_R |k|
        assert!((splitting - 2.0 * rashba.alpha_r * 0.1).abs() < 1e-6);
    }

    #[test]
    fn test_spin_texture_perpendicularity() {
        let rashba = RashbaSystem::au_111();
        let k = Vector3::new(1.0, 0.0, 0.0); // k along x

        let spin = rashba.spin_texture(k);

        // Spin should be perpendicular to k
        assert!(spin.dot(&k).abs() < 1e-10);
        // Spin should be normalized
        assert!((spin.magnitude() - 1.0).abs() < 1e-10);
        // For k_x, spin should be along ±y
        assert!(spin.y.abs() > 0.9);
    }

    #[test]
    fn test_edelstein_effect() {
        let rashba = RashbaSystem::ag_bi();
        let j_c = 1.0e10; // A/m²
        let j_dir = Vector3::new(1.0, 0.0, 0.0);

        let spin_density = rashba.edelstein_spin_density(j_c, j_dir);

        // Should be non-zero
        assert!(spin_density.magnitude() > 0.0);
        // Should be perpendicular to current
        assert!(spin_density.dot(&j_dir).abs() < 1e-10);
    }

    #[test]
    fn test_inverse_edelstein() {
        let rashba = RashbaSystem::bitei();
        let spin_dens = Vector3::new(0.0, 1.0e10, 0.0);

        let current = rashba.inverse_edelstein_current(spin_dens);

        // Should produce current
        assert!(current.magnitude() > 0.0);
        // Current perpendicular to spin and normal
        assert!(current.dot(&spin_dens).abs() < 1e-5);
    }

    #[test]
    fn test_fermi_momentum() {
        let rashba = RashbaSystem::au_111();
        let k_f = rashba.fermi_momentum();

        // Should be positive
        assert!(k_f > 0.0);
        // Reasonable range for 2DEG
        assert!(k_f > 0.01);
        assert!(k_f < 1.0);
    }

    #[test]
    fn test_spin_precession_length() {
        let rashba = RashbaSystem::gaas_2deg(1.0);
        let l_so = rashba.spin_precession_length();

        // Should be in nanometer range
        assert!(l_so > 1.0); // nm
        assert!(l_so < 1000.0);
    }

    #[test]
    fn test_spintronics_suitability() {
        let strong = RashbaSystem::ag_bi();
        let weak = RashbaSystem::au_111().with_alpha_r(0.01);

        assert!(strong.is_suitable_for_spintronics());
        assert!(!weak.is_suitable_for_spintronics());
    }

    #[test]
    fn test_builder_pattern() {
        let custom = RashbaSystem::ag_bi()
            .with_alpha_r(5.0)
            .with_carrier_density(1.0e18);

        assert_eq!(custom.alpha_r, 5.0);
        assert_eq!(custom.carrier_density, 1.0e18);
    }
}
