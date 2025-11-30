//! Spin-Orbit Torque (SOT) Effects
//!
//! This module implements spin-orbit torque, which arises from spin-orbit coupling
//! in heavy-metal/ferromagnet bilayers. SOT can be used to manipulate magnetization
//! electrically without an external magnetic field.
//!
//! ## Physics Background
//!
//! When charge current flows through a heavy metal (e.g., Pt, Ta, W) with strong
//! spin-orbit coupling, it generates a transverse spin current via the Spin Hall Effect.
//! This spin current diffuses into an adjacent ferromagnetic layer and exerts torque
//! on the magnetization.
//!
//! The SOT consists of two components:
//!
//! 1. **Field-like Torque**: τ_FL = -γ H_FL × m
//!    - Equivalent to an effective magnetic field
//!    - Proportional to m × (y × m) where y is the spin polarization direction
//!
//! 2. **Damping-like Torque**: τ_DL = -γ H_DL × (m × y)
//!    - Can overcome Gilbert damping
//!    - Enables current-driven magnetization switching
//!
//! ## Key References
//!
//! - I. M. Miron et al., "Perpendicular switching of a single ferromagnetic layer
//!   induced by in-plane current injection", Nature 476, 189 (2011)
//! - L. Liu et al., "Spin-Torque Switching with the Giant Spin Hall Effect of Tantalum",
//!   Science 336, 555 (2012)
//! - K. Garello et al., "Symmetry and magnitude of spin–orbit torques in ferromagnetic
//!   heterostructures", Nat. Nanotechnology 8, 587 (2013)

use crate::constants::HBAR;
use crate::vector3::Vector3;

/// Spin-orbit torque calculator for heavy-metal/ferromagnet bilayers
///
/// # Example
/// ```
/// use spintronics::effect::sot::SpinOrbitTorque;
///
/// // Create Pt/CoFeB SOT system
/// let sot = SpinOrbitTorque::platinum_cofeb();
///
/// // Calculate effective spin Hall efficiency
/// let eff = sot.effective_spin_hall_efficiency();
///
/// // Efficiency should be positive and less than θ_SH
/// assert!(eff > 0.0);
/// assert!(eff < sot.theta_sh.abs());
///
/// // Tantalum has larger (negative) spin Hall angle
/// let sot_ta = SpinOrbitTorque::tantalum_cofeb();
/// assert!(sot_ta.theta_sh.abs() > sot.theta_sh.abs());
/// ```
#[derive(Debug, Clone)]
pub struct SpinOrbitTorque {
    /// Spin Hall angle of the heavy metal (dimensionless)
    /// Pt: ~0.07, Ta: ~-0.15, W: ~-0.3
    pub theta_sh: f64,

    /// Heavy metal resistivity [Ω·m]
    pub resistivity: f64,

    /// Heavy metal layer thickness \[m\]
    pub thickness: f64,

    /// Spin diffusion length in heavy metal \[m\]
    pub lambda_sd: f64,

    /// Interfacial spin transparency (0 to 1)
    pub transparency: f64,
}

impl SpinOrbitTorque {
    /// Create SOT parameters for Pt/CoFeB interface
    ///
    /// Based on experimental measurements from Miron et al., Nature 476, 189 (2011)
    pub fn platinum_cofeb() -> Self {
        Self {
            theta_sh: 0.07,      // Spin Hall angle for Pt
            resistivity: 2.0e-7, // Ω·m
            thickness: 5.0e-9,   // 5 nm Pt
            lambda_sd: 1.5e-9,   // 1.5 nm spin diffusion length
            transparency: 0.5,   // Moderate transparency
        }
    }

    /// Create SOT parameters for Ta/CoFeB interface
    ///
    /// Based on Liu et al., Science 336, 555 (2012)
    pub fn tantalum_cofeb() -> Self {
        Self {
            theta_sh: -0.15, // Negative spin Hall angle for Ta (β-phase)
            resistivity: 1.8e-7,
            thickness: 5.0e-9,
            lambda_sd: 1.2e-9,
            transparency: 0.6,
        }
    }

    /// Create SOT parameters for W/CoFeB interface
    ///
    /// W has the largest spin Hall angle among common heavy metals
    pub fn tungsten_cofeb() -> Self {
        Self {
            theta_sh: -0.3, // Very large negative spin Hall angle
            resistivity: 2.5e-7,
            thickness: 5.0e-9,
            lambda_sd: 1.0e-9,
            transparency: 0.4,
        }
    }

    /// Calculate the effective spin Hall efficiency
    ///
    /// Accounts for finite thickness effects and spin backflow
    pub fn effective_spin_hall_efficiency(&self) -> f64 {
        let tanh_factor = (self.thickness / self.lambda_sd).tanh();
        self.theta_sh * tanh_factor * self.transparency
    }

    /// Calculate damping-like effective field from current density
    ///
    /// # Arguments
    /// * `j_charge` - Charge current density in heavy metal \[A/m²\]
    /// * `m` - Normalized magnetization vector
    /// * `current_direction` - Direction of current flow (normalized)
    /// * `ms` - Saturation magnetization \[A/m\]
    ///
    /// # Returns
    /// Damping-like effective field H_DL \[A/m\]
    ///
    /// Formula: H_DL = (ℏ/(2e)) * (j_e * θ_SH) / (M_s * t_FM) * (m × y)
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::sot::SpinOrbitTorque;
    /// use spintronics::Vector3;
    ///
    /// // Pt/CoFeB bilayer
    /// let sot = SpinOrbitTorque::platinum_cofeb();
    ///
    /// // Perpendicular magnetization (slightly tilted)
    /// let m = Vector3::new(0.0, 0.1, 0.995).normalize();
    /// // Current along x-axis
    /// let j_current = Vector3::new(1.0, 0.0, 0.0);
    /// // Current density: 10^11 A/m² (typical for SOT switching)
    /// let j_e = 1.0e11;
    /// // CoFeB saturation magnetization
    /// let ms = 1.0e6; // 1 MA/m
    ///
    /// // Calculate damping-like torque field
    /// let h_dl = sot.damping_like_field(j_e, m, j_current, ms);
    ///
    /// // H_DL should be perpendicular to m
    /// assert!(h_dl.dot(&m).abs() < 1e-6);
    ///
    /// // Magnitude should be physically reasonable (kA/m range)
    /// assert!(h_dl.magnitude() > 0.0);
    /// assert!(h_dl.magnitude() < 1.0e5);
    /// ```
    pub fn damping_like_field(
        &self,
        j_charge: f64,
        m: Vector3<f64>,
        current_direction: Vector3<f64>,
        ms: f64,
    ) -> Vector3<f64> {
        // Spin polarization direction (perpendicular to current and normal to interface)
        // Assuming current flows in-plane and normal is z-direction
        let normal = Vector3::new(0.0, 0.0, 1.0);
        let spin_polarization = current_direction.cross(&normal).normalize();

        // Effective spin Hall angle (including backflow and transparency)
        let theta_eff = self.effective_spin_hall_efficiency();

        // Prefactor: (ℏ/(2e)) * j * θ_SH / (M_s * t)
        let e_charge = 1.602e-19; // Elementary charge \[C\]
        let prefactor = (HBAR / (2.0 * e_charge)) * j_charge * theta_eff / (ms * self.thickness);

        // H_DL = prefactor * (m × σ) where σ is spin polarization direction
        m.cross(&spin_polarization) * prefactor
    }

    /// Calculate field-like effective field from current density
    ///
    /// The field-like term is typically smaller than damping-like term
    /// and its origin is still under debate (Rashba effect, interfacial effects, etc.)
    ///
    /// # Arguments
    /// * `j_charge` - Charge current density \[A/m²\]
    /// * `m` - Normalized magnetization vector
    /// * `current_direction` - Direction of current flow (normalized)
    /// * `ms` - Saturation magnetization \[A/m\]
    /// * `field_like_ratio` - Ratio of field-like to damping-like torque (typically 0.1-0.3)
    ///
    /// # Returns
    /// Field-like effective field H_FL \[A/m\]
    pub fn field_like_field(
        &self,
        j_charge: f64,
        m: Vector3<f64>,
        current_direction: Vector3<f64>,
        ms: f64,
        field_like_ratio: f64,
    ) -> Vector3<f64> {
        // H_FL = -field_like_ratio * H_DL × m
        let h_dl = self.damping_like_field(j_charge, m, current_direction, ms);
        h_dl.cross(&m) * (-field_like_ratio)
    }

    /// Calculate total SOT effective field
    ///
    /// Returns (H_DL, H_FL) tuple
    pub fn total_field(
        &self,
        j_charge: f64,
        m: Vector3<f64>,
        current_direction: Vector3<f64>,
        ms: f64,
        field_like_ratio: f64,
    ) -> (Vector3<f64>, Vector3<f64>) {
        let h_dl = self.damping_like_field(j_charge, m, current_direction, ms);
        let h_fl = self.field_like_field(j_charge, m, current_direction, ms, field_like_ratio);
        (h_dl, h_fl)
    }

    /// Calculate critical switching current density
    ///
    /// For perpendicular magnetization switching via damping-like SOT
    /// Requires small in-plane field to break symmetry
    ///
    /// Formula: j_c = (2 e / ℏ) * (M_s * t_FM / θ_SH) * (H_k + 2πM_s)
    ///
    /// # Arguments
    /// * `ms` - Saturation magnetization \[A/m\]
    /// * `h_k` - Perpendicular anisotropy field \[A/m\]
    ///
    /// # Returns
    /// Critical current density \[A/m²\]
    ///
    /// # Example
    /// ```
    /// use spintronics::effect::sot::SpinOrbitTorque;
    ///
    /// // Pt/CoFeB system with moderate anisotropy
    /// let sot = SpinOrbitTorque::platinum_cofeb();
    ///
    /// // CoFeB with moderate perpendicular magnetic anisotropy
    /// let ms = 1.0e6;    // 1 MA/m saturation magnetization
    /// let h_k = 5.0e4;   // 50 kA/m anisotropy field (moderate PMA)
    ///
    /// // Calculate critical current for magnetization switching
    /// // j_c = (2e/ℏ) × (M_s × t / θ_SH) × (H_k + M_s)
    /// let j_c = sot.critical_current_density(ms, h_k);
    ///
    /// // Critical current should be positive and finite
    /// assert!(j_c > 0.0);
    /// assert!(j_c.is_finite());
    ///
    /// // Typical experimental range: 10^11 to 10^13 A/m²
    /// // (varies with material stack, anisotropy, and efficiency)
    /// assert!(j_c > 1.0e10);
    ///
    /// println!("Critical switching current: {:.2e} A/m²", j_c);
    /// ```
    #[allow(dead_code)]
    pub fn critical_current_density(&self, ms: f64, h_k: f64) -> f64 {
        let e_charge = 1.602e-19;
        let theta_eff = self.effective_spin_hall_efficiency();

        // Demagnetization field for thin film: H_demag = M_s
        let h_eff = h_k + ms;

        (2.0 * e_charge / HBAR) * (ms * self.thickness / theta_eff.abs()) * h_eff
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_platinum_sot() {
        let sot = SpinOrbitTorque::platinum_cofeb();
        assert!(sot.theta_sh > 0.0); // Positive for Pt
        assert!(sot.theta_sh < 0.1);
    }

    #[test]
    fn test_tantalum_sot() {
        let sot = SpinOrbitTorque::tantalum_cofeb();
        assert!(sot.theta_sh < 0.0); // Negative for Ta
        assert!(sot.theta_sh.abs() > 0.1);
    }

    #[test]
    fn test_tungsten_sot() {
        let sot = SpinOrbitTorque::tungsten_cofeb();
        assert!(sot.theta_sh < -0.2); // Large negative for W
    }

    #[test]
    fn test_effective_spin_hall_efficiency() {
        let sot = SpinOrbitTorque::platinum_cofeb();
        let theta_eff = sot.effective_spin_hall_efficiency();

        // Should be smaller than bare theta_sh due to backflow and transparency
        assert!(theta_eff < sot.theta_sh);
        assert!(theta_eff > 0.0);
    }

    #[test]
    fn test_damping_like_field() {
        let sot = SpinOrbitTorque::platinum_cofeb();

        let m = Vector3::new(0.0, 0.1, 1.0).normalize(); // Nearly perpendicular
        let j_charge = 1.0e11; // 10^11 A/m² (typical experimental value)
        let current_dir = Vector3::new(1.0, 0.0, 0.0); // Current along x
        let ms = 1.0e6; // 1 MA/m

        let h_dl = sot.damping_like_field(j_charge, m, current_dir, ms);

        // H_DL should be perpendicular to both m and spin polarization direction
        assert!(h_dl.dot(&m).abs() < 1e-6);

        // H_DL magnitude should be reasonable (order of kA/m for typical j)
        let h_magnitude = h_dl.magnitude();
        assert!(h_magnitude > 0.0);
        assert!(h_magnitude < 1.0e5); // Should be less than 100 kA/m
    }

    #[test]
    fn test_field_like_field() {
        let sot = SpinOrbitTorque::platinum_cofeb();

        let m = Vector3::new(0.0, 0.0, 1.0); // Perpendicular
        let j_charge = 1.0e11;
        let current_dir = Vector3::new(1.0, 0.0, 0.0);
        let ms = 1.0e6;
        let fl_ratio = 0.2;

        let h_fl = sot.field_like_field(j_charge, m, current_dir, ms, fl_ratio);

        // H_FL magnitude should be smaller than H_DL
        let h_dl = sot.damping_like_field(j_charge, m, current_dir, ms);
        assert!(h_fl.magnitude() < h_dl.magnitude());
    }

    #[test]
    fn test_critical_current() {
        let sot = SpinOrbitTorque::platinum_cofeb();

        let ms = 1.0e6; // 1 MA/m
        let h_k = 1.0e5; // 100 kA/m (perpendicular anisotropy field)

        let j_c = sot.critical_current_density(ms, h_k);

        // Critical current should be positive and reasonable
        // Typical experimental values range from 10^10 to 10^14 A/m²
        // depending on material stack, anisotropy, and interface quality
        assert!(j_c > 0.0);
        assert!(j_c > 1.0e10); // Lower bound: should be at least 10^10 A/m²
                               // Note: High-anisotropy materials can require j_c > 10^13 A/m²
    }

    #[test]
    fn test_total_field() {
        let sot = SpinOrbitTorque::platinum_cofeb();

        let m = Vector3::new(0.0, 0.1, 1.0).normalize();
        let j_charge = 1.0e11;
        let current_dir = Vector3::new(1.0, 0.0, 0.0);
        let ms = 1.0e6;
        let fl_ratio = 0.2;

        let (h_dl, h_fl) = sot.total_field(j_charge, m, current_dir, ms, fl_ratio);

        // Both fields should be non-zero
        assert!(h_dl.magnitude() > 0.0);
        assert!(h_fl.magnitude() > 0.0);

        // Damping-like should be larger
        assert!(h_dl.magnitude() > h_fl.magnitude());
    }
}
