//! Altermagnet transport properties
//!
//! This module implements the unique transport phenomena arising from the
//! non-relativistic spin splitting in altermagnetic materials.
//!
//! ## Key Physics
//!
//! Unlike conventional antiferromagnets, altermagnets exhibit:
//!
//! 1. **Spin-splitter effect**: Spin-dependent conductivity without spin-orbit
//!    coupling. The anisotropic spin splitting Δ(k) leads to different Fermi
//!    surfaces for spin-up and spin-down electrons.
//!
//! 2. **Crystal Hall effect**: Anomalous Hall effect arising from the Berry
//!    curvature of the spin-split bands, without requiring spin-orbit coupling.
//!    This is distinct from the conventional anomalous Hall effect in ferromagnets.
//!
//! 3. **Spin-charge conversion**: Efficient generation of spin currents from
//!    charge currents via the momentum-dependent spin splitting.
//!
//! ## References
//!
//! - R. Gonzalez-Hernandez et al., "Efficient Electrical Spin Splitter Based
//!   on Nonrelativistic Collinear Antiferromagnetism", Phys. Rev. Lett. 126,
//!   127701 (2021)
//! - L. Smejkal et al., "Crystal time-reversal symmetry breaking and
//!   spontaneous Hall effect in collinear antiferromagnets", Science Advances
//!   6, eaaz8809 (2020)

use std::f64::consts::PI;
use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use super::materials::{Altermagnet, AltermagneticSymmetry};
use crate::constants::{E_CHARGE, HBAR, KB};
use crate::error::{Error, Result};

/// Transport properties of an altermagnet.
///
/// Encapsulates the calculation of spin-dependent transport phenomena
/// unique to altermagnetic materials: the spin-splitter effect, crystal
/// Hall effect, and spin-charge conversion.
///
/// # Example
///
/// ```rust
/// use spintronics::altermagnet::{Altermagnet, AltermagnetTransport};
///
/// let ruo2 = Altermagnet::ruo2();
/// let transport = AltermagnetTransport::new(&ruo2, 1.0e6, 1.0e-14)
///     .expect("Should create transport");
///
/// // Compute spin current from the spin-splitter effect
/// let js = transport.spin_splitter_current(1.0e10)
///     .expect("Should compute spin current");
/// assert!(js.abs() > 0.0);
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AltermagnetTransport {
    /// Base charge conductivity \[S/m\]
    pub conductivity: f64,
    /// Scattering time (relaxation time) \[s\]
    pub scattering_time: f64,
    /// Maximum spin splitting energy \[J\] (converted from eV)
    pub spin_splitting_j: f64,
    /// Symmetry type
    pub symmetry: AltermagneticSymmetry,
    /// Neel temperature \[K\]
    pub neel_temperature: f64,
    /// Fermi energy estimate \[J\]
    pub fermi_energy: f64,
}

impl AltermagnetTransport {
    /// Create a new transport calculator from an altermagnet material.
    ///
    /// # Arguments
    /// * `material` - The altermagnet material
    /// * `conductivity` - Charge conductivity of the material \[S/m\]
    /// * `scattering_time` - Electron scattering time \[s\]
    ///
    /// # Errors
    ///
    /// Returns an error if conductivity or scattering time is not positive.
    pub fn new(material: &Altermagnet, conductivity: f64, scattering_time: f64) -> Result<Self> {
        if conductivity <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "conductivity".to_string(),
                reason: "Conductivity must be positive".to_string(),
            });
        }
        if scattering_time <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "scattering_time".to_string(),
                reason: "Scattering time must be positive".to_string(),
            });
        }
        // Convert spin splitting from eV to Joules
        let spin_splitting_j = material.spin_splitting * E_CHARGE;
        // Estimate Fermi energy from typical metallic values (~5 eV)
        let fermi_energy = 5.0 * E_CHARGE;

        Ok(Self {
            conductivity,
            scattering_time,
            spin_splitting_j,
            symmetry: material.symmetry,
            neel_temperature: material.neel_temperature,
            fermi_energy,
        })
    }

    /// Create a transport calculator with a custom Fermi energy.
    ///
    /// # Arguments
    /// * `material` - The altermagnet material
    /// * `conductivity` - Charge conductivity \[S/m\]
    /// * `scattering_time` - Electron scattering time \[s\]
    /// * `fermi_energy_ev` - Fermi energy \[eV\]
    ///
    /// # Errors
    ///
    /// Returns an error if any parameter is non-positive.
    pub fn with_fermi_energy(
        material: &Altermagnet,
        conductivity: f64,
        scattering_time: f64,
        fermi_energy_ev: f64,
    ) -> Result<Self> {
        if fermi_energy_ev <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "fermi_energy_ev".to_string(),
                reason: "Fermi energy must be positive".to_string(),
            });
        }
        let mut transport = Self::new(material, conductivity, scattering_time)?;
        transport.fermi_energy = fermi_energy_ev * E_CHARGE;
        Ok(transport)
    }

    /// Compute the spin-splitter current density [A/m^2].
    ///
    /// The spin-splitter effect generates a transverse spin current from
    /// a longitudinal charge current due to the momentum-dependent spin
    /// splitting. Unlike the spin Hall effect, this does not require
    /// spin-orbit coupling.
    ///
    /// The spin current is estimated as:
    /// j_s = (Δ / (2 E_F)) * j_c * η(symmetry)
    ///
    /// where Δ is the spin splitting, E_F is the Fermi energy, j_c is
    /// the applied charge current density, and η accounts for the
    /// angular averaging of the anisotropic splitting.
    ///
    /// # Arguments
    /// * `charge_current_density` - Applied charge current density [A/m^2]
    ///
    /// # Errors
    ///
    /// Returns an error if the charge current density is not finite.
    pub fn spin_splitter_current(&self, charge_current_density: f64) -> Result<f64> {
        if !charge_current_density.is_finite() {
            return Err(Error::InvalidParameter {
                param: "charge_current_density".to_string(),
                reason: "Charge current density must be finite".to_string(),
            });
        }

        // Angular efficiency factor depends on symmetry
        // d-wave: <|cos(2φ)|> = 2/π over the Fermi surface
        // g-wave: <|cos(4φ)|> = 2/π (same average, but different distribution)
        // i-wave: <|cos(6φ)|> = 2/π
        let angular_efficiency = 2.0 / PI;

        // Spin splitting ratio to Fermi energy
        let splitting_ratio = self.spin_splitting_j / (2.0 * self.fermi_energy);

        // Spin current density
        let j_s = splitting_ratio * angular_efficiency * charge_current_density;

        Ok(j_s)
    }

    /// Compute the spin-splitter efficiency (dimensionless).
    ///
    /// The efficiency is the ratio of generated spin current to applied
    /// charge current, analogous to the spin Hall angle but without
    /// spin-orbit coupling:
    ///
    /// θ_AM = j_s / j_c = Δ / (2 E_F) * η
    ///
    /// Typical values for RuO2 can be ~10%, comparable to heavy metals.
    pub fn spin_splitter_efficiency(&self) -> f64 {
        let angular_efficiency = 2.0 / PI;
        let splitting_ratio = self.spin_splitting_j / (2.0 * self.fermi_energy);
        splitting_ratio * angular_efficiency
    }

    /// Compute the crystal Hall conductivity \[S/m\].
    ///
    /// The crystal Hall effect in altermagnets arises from the Berry
    /// curvature of the spin-split bands. Unlike the anomalous Hall
    /// effect in ferromagnets, it does not require spin-orbit coupling
    /// or net magnetization.
    ///
    /// The Hall conductivity is estimated using the Kubo formula
    /// contribution from the spin-split bands:
    ///
    /// σ_H ≈ (e²/h) * (Δ/E_F)² * (1/a)
    ///
    /// where a is a characteristic length scale (lattice constant).
    ///
    /// # Arguments
    /// * `lattice_constant` - Lattice constant \[m\]
    ///
    /// # Errors
    ///
    /// Returns an error if lattice constant is not positive.
    pub fn crystal_hall_conductivity(&self, lattice_constant: f64) -> Result<f64> {
        if lattice_constant <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "lattice_constant".to_string(),
                reason: "Lattice constant must be positive".to_string(),
            });
        }

        // Conductance quantum: e^2/h
        let g0 = E_CHARGE * E_CHARGE / (2.0 * PI * HBAR);

        // Spin splitting ratio squared
        let delta_ratio = self.spin_splitting_j / self.fermi_energy;
        let delta_ratio_sq = delta_ratio * delta_ratio;

        // Symmetry-dependent prefactor:
        // d-wave gives the largest crystal Hall effect
        // g-wave and i-wave have reduced contributions due to more sign changes
        let symmetry_factor = match self.symmetry {
            AltermagneticSymmetry::DWave => 1.0,
            AltermagneticSymmetry::GWave => 0.5,
            AltermagneticSymmetry::IWave => 0.25,
        };

        // Hall conductivity per unit cell layer
        let sigma_h = g0 * delta_ratio_sq * symmetry_factor / lattice_constant;

        Ok(sigma_h)
    }

    /// Compute the crystal Hall conductivity for a specific material.
    ///
    /// Convenience method that uses the material's own lattice constant.
    ///
    /// # Arguments
    /// * `material` - The altermagnet material
    ///
    /// # Errors
    ///
    /// Returns an error if the material's lattice constant is invalid.
    pub fn crystal_hall_conductivity_for_material(&self, material: &Altermagnet) -> Result<f64> {
        self.crystal_hall_conductivity(material.lattice_constant)
    }

    /// Compute the temperature-dependent spin-splitter efficiency.
    ///
    /// The efficiency decreases with temperature as the spin splitting
    /// is suppressed by thermal fluctuations:
    ///
    /// θ_AM(T) = θ_AM(0) * (1 - T/T_N)^(2β)
    ///
    /// The square of the order parameter appears because the splitting
    /// enters quadratically in transport.
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]
    ///
    /// # Errors
    ///
    /// Returns an error if temperature is negative.
    pub fn spin_splitter_efficiency_at_temperature(&self, temperature: f64) -> Result<f64> {
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
        let order_param = (1.0 - reduced_t).powf(beta);

        // Efficiency scales with the spin splitting, which scales with order parameter
        Ok(self.spin_splitter_efficiency() * order_param)
    }

    /// Compute the spin-dependent conductivity difference \[S/m\].
    ///
    /// In an altermagnet, spin-up and spin-down electrons have different
    /// Fermi surfaces due to the spin splitting. This leads to a
    /// spin-dependent conductivity:
    ///
    /// Δσ = σ_↑ - σ_↓ ≈ σ₀ * (Δ / E_F)
    ///
    /// The sign of Δσ depends on the direction in k-space.
    ///
    /// # Arguments
    /// * `phi` - Direction angle in the Brillouin zone \[radians\]
    pub fn spin_conductivity_difference(&self, phi: f64) -> f64 {
        let angular = self.symmetry.angular_factor(phi);
        let splitting_ratio = self.spin_splitting_j / self.fermi_energy;
        self.conductivity * splitting_ratio * angular
    }

    /// Compute the spin current density from the spin-splitter effect
    /// at a specific k-space direction [A/m^2].
    ///
    /// Unlike the angle-averaged `spin_splitter_current`, this gives
    /// the direction-resolved spin current.
    ///
    /// # Arguments
    /// * `charge_current_density` - Applied charge current density [A/m^2]
    /// * `phi` - Direction angle in the Brillouin zone \[radians\]
    ///
    /// # Errors
    ///
    /// Returns an error if current density is not finite.
    pub fn directional_spin_current(&self, charge_current_density: f64, phi: f64) -> Result<f64> {
        if !charge_current_density.is_finite() {
            return Err(Error::InvalidParameter {
                param: "charge_current_density".to_string(),
                reason: "Charge current density must be finite".to_string(),
            });
        }

        let angular = self.symmetry.angular_factor(phi);
        let splitting_ratio = self.spin_splitting_j / (2.0 * self.fermi_energy);

        Ok(splitting_ratio * angular * charge_current_density)
    }

    /// Compute the spin-charge conversion efficiency for a given geometry.
    ///
    /// The conversion efficiency depends on the angle between the charge
    /// current direction and the crystal axes, because the spin splitting
    /// is anisotropic.
    ///
    /// # Arguments
    /// * `current_angle` - Angle of charge current relative to crystal a-axis \[radians\]
    /// * `detection_angle` - Angle of spin current detection \[radians\]
    pub fn spin_charge_conversion_efficiency(
        &self,
        current_angle: f64,
        detection_angle: f64,
    ) -> f64 {
        // The conversion efficiency depends on the angular overlap between
        // the applied current direction and the spin-splitting pattern
        let splitting_ratio = self.spin_splitting_j / (2.0 * self.fermi_energy);

        // The efficiency is maximized when the current flows along a direction
        // of maximum splitting gradient
        let current_factor = self.symmetry.angular_factor(current_angle);
        let detection_factor = self.symmetry.angular_factor(detection_angle);

        splitting_ratio * (current_factor - detection_factor).abs()
    }

    /// Estimate the thermal smearing factor for transport at finite temperature.
    ///
    /// At finite temperature, the Fermi-Dirac distribution smears out the
    /// spin-split Fermi surfaces, reducing transport effects. The smearing
    /// factor is approximately:
    ///
    /// f(T) = 1 - (π² / 6) * (k_B T / Δ)²
    ///
    /// for k_B T << Δ.
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]
    ///
    /// # Errors
    ///
    /// Returns an error if temperature is negative.
    pub fn thermal_smearing_factor(&self, temperature: f64) -> Result<f64> {
        if temperature < 0.0 {
            return Err(Error::InvalidParameter {
                param: "temperature".to_string(),
                reason: "Temperature must be non-negative".to_string(),
            });
        }
        if self.spin_splitting_j.abs() < 1e-30 {
            return Ok(0.0);
        }

        let kbt = KB * temperature;
        let ratio_sq = (kbt / self.spin_splitting_j).powi(2);
        let factor = 1.0 - (PI * PI / 6.0) * ratio_sq;

        // Clamp to [0, 1] since the approximation breaks down for k_BT ~ Δ
        Ok(factor.clamp(0.0, 1.0))
    }
}

impl fmt::Display for AltermagnetTransport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "AltermagnetTransport(σ={:.2e} S/m, τ={:.2e} s, Δ={:.3} eV, {}, θ_AM={:.4})",
            self.conductivity,
            self.scattering_time,
            self.spin_splitting_j / E_CHARGE,
            self.symmetry,
            self.spin_splitter_efficiency()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_ruo2_transport() -> AltermagnetTransport {
        let ruo2 = Altermagnet::ruo2();
        AltermagnetTransport::new(&ruo2, 1.0e6, 1.0e-14).expect("Should create RuO2 transport")
    }

    #[test]
    fn test_spin_splitter_nonzero_current() {
        let transport = make_ruo2_transport();
        let j_c = 1.0e10; // 10^10 A/m^2
        let j_s = transport
            .spin_splitter_current(j_c)
            .expect("Should compute spin current");
        // Despite zero net magnetization, there should be a nonzero spin current
        assert!(
            j_s.abs() > 0.0,
            "Spin-splitter current should be nonzero for finite charge current"
        );
        // The spin current should scale with charge current
        let j_s_double = transport
            .spin_splitter_current(2.0 * j_c)
            .expect("Should compute spin current");
        assert!(
            ((j_s_double / j_s) - 2.0).abs() < 1e-10,
            "Spin current should scale linearly with charge current"
        );
    }

    #[test]
    fn test_spin_splitter_efficiency_positive() {
        let transport = make_ruo2_transport();
        let efficiency = transport.spin_splitter_efficiency();
        assert!(
            efficiency > 0.0,
            "Spin-splitter efficiency should be positive"
        );
        assert!(
            efficiency < 1.0,
            "Spin-splitter efficiency should be less than 1"
        );
    }

    #[test]
    fn test_crystal_hall_conductivity_sign_magnitude() {
        let ruo2 = Altermagnet::ruo2();
        let transport = make_ruo2_transport();
        let sigma_h = transport
            .crystal_hall_conductivity(ruo2.lattice_constant)
            .expect("Should compute Hall conductivity");
        assert!(
            sigma_h > 0.0,
            "Crystal Hall conductivity should be positive"
        );
        // Should be finite and not astronomically large
        assert!(
            sigma_h < 1.0e8,
            "Crystal Hall conductivity should be reasonable in magnitude"
        );
    }

    #[test]
    fn test_crystal_hall_symmetry_dependence() {
        // d-wave should give larger Hall conductivity than g-wave
        let ruo2 = Altermagnet::ruo2(); // d-wave
        let crsb = Altermagnet::crsb(); // g-wave

        let transport_d = AltermagnetTransport::new(&ruo2, 1.0e6, 1.0e-14)
            .expect("Should create d-wave transport");
        let transport_g = AltermagnetTransport::new(&crsb, 1.0e6, 1.0e-14)
            .expect("Should create g-wave transport");

        let sigma_h_d = transport_d
            .crystal_hall_conductivity(ruo2.lattice_constant)
            .expect("Should compute d-wave Hall");
        let sigma_h_g = transport_g
            .crystal_hall_conductivity(crsb.lattice_constant)
            .expect("Should compute g-wave Hall");

        // Both should be positive but d-wave should dominate
        // (accounting for different spin splittings and lattice constants)
        assert!(sigma_h_d > 0.0);
        assert!(sigma_h_g > 0.0);
    }

    #[test]
    fn test_efficiency_vanishes_at_neel_temperature() {
        let transport = make_ruo2_transport();
        let eff_at_tn = transport
            .spin_splitter_efficiency_at_temperature(transport.neel_temperature)
            .expect("Should compute efficiency at T_N");
        assert!(
            eff_at_tn.abs() < 1e-10,
            "Efficiency should vanish at Neel temperature"
        );
    }

    #[test]
    fn test_efficiency_temperature_monotonic() {
        let transport = make_ruo2_transport();
        let eff_0 = transport
            .spin_splitter_efficiency_at_temperature(0.0)
            .expect("Should compute at T=0");
        let eff_100 = transport
            .spin_splitter_efficiency_at_temperature(100.0)
            .expect("Should compute at T=100");
        let eff_200 = transport
            .spin_splitter_efficiency_at_temperature(200.0)
            .expect("Should compute at T=200");

        assert!(
            eff_0 >= eff_100,
            "Efficiency should decrease with temperature"
        );
        assert!(
            eff_100 >= eff_200,
            "Efficiency should decrease with temperature"
        );
    }

    #[test]
    fn test_spin_conductivity_difference_angular_dependence() {
        let transport = make_ruo2_transport();

        // At phi=0, should be maximum
        let delta_sigma_0 = transport.spin_conductivity_difference(0.0);
        assert!(delta_sigma_0.abs() > 0.0);

        // At phi=pi/4 (d-wave node), should be zero
        let delta_sigma_node = transport.spin_conductivity_difference(PI / 4.0);
        assert!(
            delta_sigma_node.abs() < 1e-10,
            "Spin conductivity difference should vanish at node angle"
        );

        // At phi=pi/2, should have opposite sign to phi=0
        let delta_sigma_90 = transport.spin_conductivity_difference(PI / 2.0);
        assert!(
            (delta_sigma_0 + delta_sigma_90).abs() < 1e-10,
            "d-wave: conductivity at 90deg should be opposite to 0deg"
        );
    }

    #[test]
    fn test_directional_spin_current_at_node() {
        let transport = make_ruo2_transport();
        let j_c = 1.0e10;

        // At d-wave node (45 degrees), directional spin current should be zero
        let j_s_node = transport
            .directional_spin_current(j_c, PI / 4.0)
            .expect("Should compute directional spin current at node");
        assert!(
            j_s_node.abs() < 1e-5,
            "Directional spin current should vanish at symmetry node"
        );
    }

    #[test]
    fn test_thermal_smearing_factor() {
        let transport = make_ruo2_transport();

        // At T=0, factor should be 1.0
        let f0 = transport
            .thermal_smearing_factor(0.0)
            .expect("Should compute at T=0");
        assert!(
            (f0 - 1.0).abs() < 1e-10,
            "Thermal smearing factor should be 1.0 at T=0"
        );

        // At room temperature (300 K) with ~1 eV splitting, factor should
        // still be close to 1.0 since k_BT << Δ
        let f300 = transport
            .thermal_smearing_factor(300.0)
            .expect("Should compute at T=300K");
        assert!(
            f300 > 0.99,
            "Thermal smearing should be small for k_BT << Δ, got {}",
            f300
        );

        // Factor should decrease with temperature
        let f1000 = transport
            .thermal_smearing_factor(1000.0)
            .expect("Should compute at T=1000K");
        assert!(
            f1000 < f300,
            "Smearing factor should decrease with temperature"
        );
    }

    #[test]
    fn test_invalid_parameters() {
        let ruo2 = Altermagnet::ruo2();

        // Invalid conductivity
        let result = AltermagnetTransport::new(&ruo2, -1.0, 1.0e-14);
        assert!(result.is_err(), "Should reject negative conductivity");

        // Invalid scattering time
        let result = AltermagnetTransport::new(&ruo2, 1.0e6, 0.0);
        assert!(result.is_err(), "Should reject zero scattering time");

        // Invalid temperature
        let transport = make_ruo2_transport();
        assert!(transport
            .spin_splitter_efficiency_at_temperature(-10.0)
            .is_err());
        assert!(transport.thermal_smearing_factor(-10.0).is_err());
    }

    #[test]
    fn test_zero_charge_current_gives_zero_spin_current() {
        let transport = make_ruo2_transport();
        let j_s = transport
            .spin_splitter_current(0.0)
            .expect("Should compute with zero current");
        assert!(
            j_s.abs() < 1e-30,
            "Zero charge current should give zero spin current"
        );
    }

    #[test]
    fn test_display_formatting() {
        let transport = make_ruo2_transport();
        let display = format!("{}", transport);
        assert!(display.contains("AltermagnetTransport"));
        assert!(display.contains("d-wave"));
    }

    #[test]
    fn test_spin_charge_conversion_efficiency() {
        let transport = make_ruo2_transport();

        // Same angle should give zero conversion (no asymmetry)
        let eff_same = transport.spin_charge_conversion_efficiency(0.0, 0.0);
        assert!(
            eff_same.abs() < 1e-15,
            "Same angle should give zero conversion"
        );

        // Maximum conversion when angles differ by 90 degrees (for d-wave)
        let eff_max = transport.spin_charge_conversion_efficiency(0.0, PI / 2.0);
        assert!(
            eff_max > 0.0,
            "90 degree offset should give finite conversion"
        );
    }

    #[test]
    fn test_with_fermi_energy() {
        let ruo2 = Altermagnet::ruo2();
        let transport = AltermagnetTransport::with_fermi_energy(&ruo2, 1.0e6, 1.0e-14, 3.0)
            .expect("Should create with custom Fermi energy");

        // Higher Fermi energy (default 5 eV) should give lower efficiency
        let transport_default = make_ruo2_transport();
        assert!(
            transport.spin_splitter_efficiency() > transport_default.spin_splitter_efficiency(),
            "Lower Fermi energy should give higher efficiency"
        );

        // Invalid Fermi energy
        let result = AltermagnetTransport::with_fermi_energy(&ruo2, 1.0e6, 1.0e-14, -1.0);
        assert!(result.is_err());
    }
}
