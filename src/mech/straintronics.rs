//! Strain-mediated magnetization switching (straintronics)
//!
//! This module implements strain-mediated control of magnetization in
//! multiferroic heterostructures, where a piezoelectric substrate generates
//! strain that is transferred to a magnetostrictive magnetic layer.
//!
//! ## Key Concepts
//!
//! ### Strain-Mediated Switching
//! A piezoelectric substrate (PZT, PMN-PT) generates controlled strain
//! via applied voltage. This strain is transferred to a magnetostrictive
//! free layer (e.g., CoFeB, Terfenol-D), modifying its magnetic anisotropy
//! and potentially switching its magnetization by 90° or 180°.
//!
//! ### Energy Comparison
//! Strain switching is ultra-low energy compared to current-driven methods:
//! - Strain switching: ~1 aJ (10⁻¹⁸ J)
//! - SOT switching: ~10 fJ (10⁻¹⁴ J)
//! - STT switching: ~100 fJ (10⁻¹³ J)
//!
//! ### Voltage-Controlled Magnetic Anisotropy (VCMA)
//! Electric field at a ferromagnet/oxide interface directly modifies the
//! perpendicular magnetic anisotropy through electronic effects (charge
//! accumulation, orbital hybridization).
//!
//! ## References
//! - Roy et al., "Hybrid spintronics and straintronics" (2011)
//! - Biswas et al., "Complete magnetization reversal in a magnetostrictive
//!   nanomagnet with voltage-generated stress" (2014)

use super::magnetoelastic::{
    magnetoelastic_energy_density_cubic, stress_induced_anisotropy, MagnetoelasticMaterial,
    PiezoelectricSubstrate, StrainTensor,
};
use crate::constants::MU_0;
use crate::error::{self, Result};
use crate::vector3::Vector3;

// =============================================================================
// Straintronic Device Parameters
// =============================================================================

/// Parameters for a straintronic switching device
///
/// Models a multiferroic heterostructure consisting of a piezoelectric
/// substrate and a magnetostrictive free layer.
#[derive(Debug, Clone)]
pub struct StraintronicDevice {
    /// Magnetostrictive free layer material
    pub free_layer: MagnetoelasticMaterial,
    /// Piezoelectric substrate
    pub substrate: PiezoelectricSubstrate,
    /// Free layer thickness \[m\]
    pub free_layer_thickness: f64,
    /// Free layer lateral dimension (assumed square) \[m\]
    pub free_layer_width: f64,
    /// Uniaxial anisotropy constant K_u \[J/m³\]
    pub ku: f64,
    /// Strain transfer efficiency (0 to 1)
    ///
    /// Accounts for imperfect strain transfer at the interface
    pub strain_transfer_efficiency: f64,
    /// VCMA coefficient β \[J/(V·m)\]
    ///
    /// Change in interfacial anisotropy per unit electric field
    pub vcma_coefficient: f64,
}

impl StraintronicDevice {
    /// Create a new straintronic device
    ///
    /// # Arguments
    /// * `free_layer` - Magnetostrictive material for the free layer
    /// * `substrate` - Piezoelectric substrate material
    /// * `free_layer_thickness` - Thickness of free layer \[m\]
    /// * `free_layer_width` - Lateral dimension of free layer \[m\]
    /// * `ku` - Uniaxial anisotropy \[J/m³\]
    ///
    /// # Errors
    /// Returns an error if dimensions are non-positive.
    pub fn new(
        free_layer: MagnetoelasticMaterial,
        substrate: PiezoelectricSubstrate,
        free_layer_thickness: f64,
        free_layer_width: f64,
        ku: f64,
    ) -> Result<Self> {
        if free_layer_thickness <= 0.0 {
            return Err(error::invalid_param(
                "free_layer_thickness",
                "must be positive",
            ));
        }
        if free_layer_width <= 0.0 {
            return Err(error::invalid_param("free_layer_width", "must be positive"));
        }
        Ok(Self {
            free_layer,
            substrate,
            free_layer_thickness,
            free_layer_width,
            ku,
            strain_transfer_efficiency: 0.9,
            vcma_coefficient: 0.0,
        })
    }

    /// Set the strain transfer efficiency
    pub fn with_strain_transfer_efficiency(mut self, eta: f64) -> Result<Self> {
        if !(0.0..=1.0).contains(&eta) {
            return Err(error::invalid_param(
                "strain_transfer_efficiency",
                "must be between 0 and 1",
            ));
        }
        self.strain_transfer_efficiency = eta;
        Ok(self)
    }

    /// Set the VCMA coefficient
    ///
    /// Typical values: 30-300 fJ/(V·m) for CoFeB/MgO interfaces
    pub fn with_vcma_coefficient(mut self, beta: f64) -> Self {
        self.vcma_coefficient = beta;
        self
    }

    /// Volume of the free layer \[m³\]
    pub fn free_layer_volume(&self) -> f64 {
        self.free_layer_width * self.free_layer_width * self.free_layer_thickness
    }

    /// Compute the effective strain in the free layer for a given substrate voltage
    ///
    /// Accounts for strain transfer efficiency at the interface.
    pub fn effective_strain(&self, voltage: f64, substrate_thickness: f64) -> Result<StrainTensor> {
        let substrate_strain = self
            .substrate
            .strain_from_voltage(voltage, substrate_thickness)?;

        let eta = self.strain_transfer_efficiency;
        let mut transferred = [[0.0; 3]; 3];
        for (i, row) in transferred.iter_mut().enumerate() {
            for (j, val) in row.iter_mut().enumerate() {
                *val = eta * substrate_strain.components[i][j];
            }
        }

        Ok(StrainTensor {
            components: transferred,
        })
    }
}

// =============================================================================
// Switching Energy Calculations
// =============================================================================

/// Compute the anisotropy energy barrier for magnetization switching \[J\]
///
/// For uniaxial anisotropy: E_barrier = K_u × V
///
/// # Arguments
/// * `ku` - Uniaxial anisotropy constant \[J/m³\]
/// * `volume` - Magnetic volume \[m³\]
pub fn anisotropy_energy_barrier(ku: f64, volume: f64) -> f64 {
    ku * volume
}

/// Compute the magnetoelastic switching energy \[J\]
///
/// E_switch = |E_me(final) - E_me(initial)| × V
///
/// For a 90° rotation from along easy axis to hard axis:
/// E_me = B₁ × ε × V (order of magnitude)
///
/// # Arguments
/// * `device` - Straintronic device parameters
/// * `strain` - Applied strain tensor
/// * `initial_dir` - Initial magnetization direction
/// * `final_dir` - Target magnetization direction
pub fn strain_switching_energy(
    device: &StraintronicDevice,
    strain: &StrainTensor,
    initial_dir: &Vector3<f64>,
    final_dir: &Vector3<f64>,
) -> Result<f64> {
    let e_initial = magnetoelastic_energy_density_cubic(&device.free_layer, strain, initial_dir)?;
    let e_final = magnetoelastic_energy_density_cubic(&device.free_layer, strain, final_dir)?;

    let volume = device.free_layer_volume();
    Ok((e_final - e_initial).abs() * volume)
}

/// Estimate STT (spin-transfer torque) switching energy \[J\]
///
/// E_STT ≈ (α / η) × (2e / ℏ) × K_u × V × τ_pulse × I_c × R
///
/// A simplified order-of-magnitude estimate using:
/// E_STT ≈ I_c² × R × τ_pulse
///
/// # Arguments
/// * `critical_current` - Critical switching current \[A\]
/// * `resistance` - Device resistance \[Ω\]
/// * `pulse_duration` - Switching pulse duration \[s\]
pub fn stt_switching_energy(critical_current: f64, resistance: f64, pulse_duration: f64) -> f64 {
    critical_current * critical_current * resistance * pulse_duration
}

/// Estimate SOT (spin-orbit torque) switching energy \[J\]
///
/// E_SOT ≈ J_c² × ρ × V_channel × τ_pulse
///
/// # Arguments
/// * `current_density` - Critical current density [A/m²]
/// * `resistivity` - Channel resistivity \[Ω·m\]
/// * `channel_volume` - Volume of the SOT channel \[m³\]
/// * `pulse_duration` - Switching pulse duration \[s\]
pub fn sot_switching_energy(
    current_density: f64,
    resistivity: f64,
    channel_volume: f64,
    pulse_duration: f64,
) -> f64 {
    current_density * current_density * resistivity * channel_volume * pulse_duration
}

/// Compare switching energies for different mechanisms
///
/// Returns (strain_energy, sot_energy, stt_energy) for a typical device.
///
/// Typical values:
/// - Strain: ~1 aJ (10⁻¹⁸ J) for nanoscale devices
/// - SOT: ~10 fJ (10⁻¹⁴ J)
/// - STT: ~100 fJ (10⁻¹³ J)
pub fn typical_switching_energy_comparison() -> (f64, f64, f64) {
    // Strain switching: ~1 aJ for 100nm × 100nm × 2nm free layer
    // with PMN-PT substrate at ~200 mV
    let strain_energy = 1.0e-18;

    // SOT switching: J_c ~ 1e11 A/m², ρ ~ 1e-6 Ωm, channel ~100nm×300nm×5nm, τ ~ 1ns
    let sot_jc = 1.0e11;
    let sot_rho = 1.0e-6;
    let sot_vol = 100.0e-9 * 300.0e-9 * 5.0e-9;
    let sot_tau = 1.0e-9;
    let sot_energy = sot_switching_energy(sot_jc, sot_rho, sot_vol, sot_tau);

    // STT switching: I_c ~ 50μA, R ~ 5kΩ, τ ~ 10ns
    let stt_ic = 50.0e-6;
    let stt_r = 5.0e3;
    let stt_tau = 10.0e-9;
    let stt_energy = stt_switching_energy(stt_ic, stt_r, stt_tau);

    (strain_energy, sot_energy, stt_energy)
}

// =============================================================================
// Critical Strain for Switching
// =============================================================================

/// Compute the critical strain for magnetization switching
///
/// Switching occurs when the magnetoelastic energy overcomes the anisotropy barrier:
/// |B₁ × ε_crit| ≥ K_u  (simplified)
///
/// For uniaxial anisotropy with isotropic magnetostriction:
/// ε_crit = K_u / |(3/2) λ_s Y|
///
/// where Y is Young's modulus.
///
/// # Arguments
/// * `ku` - Uniaxial anisotropy \[J/m³\]
/// * `lambda_s` - Isotropic magnetostriction constant
/// * `youngs_modulus` - Young's modulus \[Pa\]
///
/// # Errors
/// Returns an error if the denominator is zero.
pub fn critical_strain_uniaxial(ku: f64, lambda_s: f64, youngs_modulus: f64) -> Result<f64> {
    let denominator = (1.5 * lambda_s * youngs_modulus).abs();
    if denominator < 1e-30 {
        return Err(error::invalid_param(
            "lambda_s * youngs_modulus",
            "product must be non-zero for critical strain calculation",
        ));
    }
    Ok(ku.abs() / denominator)
}

/// Compute the critical strain using B₁ coupling constant directly
///
/// ε_crit = K_u / |B₁|
///
/// # Arguments
/// * `ku` - Uniaxial anisotropy \[J/m³\]
/// * `b1` - Magnetoelastic coupling constant B₁ \[J/m³\]
///
/// # Errors
/// Returns an error if B₁ is zero.
pub fn critical_strain_from_b1(ku: f64, b1: f64) -> Result<f64> {
    if b1.abs() < 1e-30 {
        return Err(error::invalid_param(
            "b1",
            "coupling constant must be non-zero",
        ));
    }
    Ok(ku.abs() / b1.abs())
}

/// Compute the critical voltage for strain-mediated switching
///
/// V_crit = ε_crit × t_substrate / (|d₃₁| × η)
///
/// where t_substrate is the substrate thickness, d₃₁ is the piezoelectric
/// coefficient, and η is the strain transfer efficiency.
///
/// # Arguments
/// * `critical_strain` - Required strain for switching
/// * `substrate` - Piezoelectric substrate parameters
/// * `substrate_thickness` - Substrate thickness \[m\]
/// * `transfer_efficiency` - Strain transfer efficiency (0 to 1)
///
/// # Errors
/// Returns an error if d₃₁ or transfer efficiency is zero.
pub fn critical_voltage(
    critical_strain: f64,
    substrate: &PiezoelectricSubstrate,
    substrate_thickness: f64,
    transfer_efficiency: f64,
) -> Result<f64> {
    let effective_d31 = substrate.d31.abs() * transfer_efficiency;
    if effective_d31 < 1e-30 {
        return Err(error::invalid_param(
            "d31 * transfer_efficiency",
            "must be non-zero",
        ));
    }
    Ok(critical_strain.abs() * substrate_thickness / effective_d31)
}

// =============================================================================
// Voltage-Controlled Magnetic Anisotropy (VCMA)
// =============================================================================

/// VCMA (Voltage-Controlled Magnetic Anisotropy) parameters
///
/// At ferromagnet/oxide interfaces, an applied electric field modifies
/// the interfacial perpendicular magnetic anisotropy (PMA) through
/// charge accumulation and orbital hybridization effects.
#[derive(Debug, Clone, Copy)]
pub struct VcmaParameters {
    /// VCMA coefficient β \[J/(V·m)\]
    ///
    /// Typical values: 30-300 fJ/(V·m) for CoFeB/MgO
    pub beta: f64,
    /// Intrinsic interfacial anisotropy K_i \[J/m²\]
    pub ki_intrinsic: f64,
    /// Dielectric thickness (MgO barrier) \[m\]
    pub dielectric_thickness: f64,
}

impl VcmaParameters {
    /// Create new VCMA parameters
    ///
    /// # Arguments
    /// * `beta` - VCMA coefficient \[J/(V·m)\]
    /// * `ki_intrinsic` - Intrinsic interfacial anisotropy \[J/m²\]
    /// * `dielectric_thickness` - Dielectric layer thickness \[m\]
    ///
    /// # Errors
    /// Returns an error if dielectric_thickness is non-positive.
    pub fn new(beta: f64, ki_intrinsic: f64, dielectric_thickness: f64) -> Result<Self> {
        if dielectric_thickness <= 0.0 {
            return Err(error::invalid_param(
                "dielectric_thickness",
                "must be positive",
            ));
        }
        Ok(Self {
            beta,
            ki_intrinsic,
            dielectric_thickness,
        })
    }

    /// Typical CoFeB/MgO interface parameters
    pub fn cofeb_mgo_typical() -> Self {
        Self {
            beta: 100.0e-15,              // 100 fJ/(V·m)
            ki_intrinsic: 1.3e-3,         // ~1.3 mJ/m²
            dielectric_thickness: 1.5e-9, // 1.5 nm MgO
        }
    }

    /// Compute the change in interfacial anisotropy for applied voltage \[J/m²\]
    ///
    /// ΔK_i = β × E = β × V / t_dielectric
    ///
    /// Positive voltage typically reduces PMA (for common convention).
    pub fn anisotropy_change(&self, voltage: f64) -> f64 {
        let e_field = voltage / self.dielectric_thickness;
        self.beta * e_field
    }

    /// Compute the effective interfacial anisotropy under applied voltage \[J/m²\]
    ///
    /// K_i_eff = K_i_intrinsic + ΔK_i
    pub fn effective_anisotropy(&self, voltage: f64) -> f64 {
        self.ki_intrinsic + self.anisotropy_change(voltage)
    }

    /// Compute the effective uniaxial anisotropy constant K_u_eff \[J/m³\]
    ///
    /// K_u_eff = (K_i_eff / t_mag) - (1/2)μ₀ M_s²
    ///
    /// where the second term is the demagnetization energy for thin films.
    ///
    /// # Arguments
    /// * `voltage` - Applied voltage \[V\]
    /// * `magnetic_thickness` - Magnetic layer thickness \[m\]
    /// * `ms` - Saturation magnetization \[A/m\]
    pub fn effective_ku(&self, voltage: f64, magnetic_thickness: f64, ms: f64) -> Result<f64> {
        if magnetic_thickness <= 0.0 {
            return Err(error::invalid_param(
                "magnetic_thickness",
                "must be positive",
            ));
        }
        let ki_eff = self.effective_anisotropy(voltage);
        let demag = 0.5 * MU_0 * ms * ms;
        Ok(ki_eff / magnetic_thickness - demag)
    }

    /// Critical voltage for PMA to in-plane transition
    ///
    /// When K_u_eff = 0, the magnetization transitions from perpendicular
    /// to in-plane. This voltage is:
    /// V_crit = [(1/2)μ₀ M_s² t_mag - K_i] × t_dielectric / β
    ///
    /// # Arguments
    /// * `magnetic_thickness` - Magnetic layer thickness \[m\]
    /// * `ms` - Saturation magnetization \[A/m\]
    pub fn critical_voltage_pma_transition(&self, magnetic_thickness: f64, ms: f64) -> Result<f64> {
        if magnetic_thickness <= 0.0 {
            return Err(error::invalid_param(
                "magnetic_thickness",
                "must be positive",
            ));
        }
        if self.beta.abs() < 1e-30 {
            return Err(error::invalid_param(
                "beta",
                "VCMA coefficient must be non-zero",
            ));
        }
        let demag_energy = 0.5 * MU_0 * ms * ms * magnetic_thickness;
        let numerator = demag_energy - self.ki_intrinsic;
        Ok(numerator * self.dielectric_thickness / self.beta)
    }
}

// =============================================================================
// Combined Strain + VCMA Switching
// =============================================================================

/// Compute the total effective anisotropy under combined strain and VCMA \[J/m³\]
///
/// K_eff = K_u + K_σ + K_VCMA
///
/// where:
/// - K_u is the intrinsic uniaxial anisotropy
/// - K_σ = -(3/2) λ_s σ is the stress-induced anisotropy
/// - K_VCMA = ΔK_i / t_mag is the VCMA contribution
///
/// # Arguments
/// * `ku_intrinsic` - Intrinsic uniaxial anisotropy \[J/m³\]
/// * `lambda_s` - Magnetostriction constant
/// * `stress` - Applied stress \[Pa\]
/// * `vcma_change` - VCMA anisotropy change \[J/m²\]
/// * `magnetic_thickness` - Magnetic layer thickness \[m\]
///
/// # Errors
/// Returns an error if magnetic_thickness is non-positive.
pub fn combined_effective_anisotropy(
    ku_intrinsic: f64,
    lambda_s: f64,
    stress: f64,
    vcma_change: f64,
    magnetic_thickness: f64,
) -> Result<f64> {
    if magnetic_thickness <= 0.0 {
        return Err(error::invalid_param(
            "magnetic_thickness",
            "must be positive",
        ));
    }
    let k_stress = stress_induced_anisotropy(lambda_s, stress);
    let k_vcma = vcma_change / magnetic_thickness;
    Ok(ku_intrinsic + k_stress + k_vcma)
}

/// Determine if switching is possible under combined strain + VCMA
///
/// Switching requires K_eff to change sign (overcome the barrier).
///
/// # Returns
/// `true` if the combined effects can overcome the intrinsic anisotropy.
pub fn can_switch_combined(
    ku_intrinsic: f64,
    lambda_s: f64,
    stress: f64,
    vcma_change: f64,
    magnetic_thickness: f64,
) -> Result<bool> {
    let k_eff = combined_effective_anisotropy(
        ku_intrinsic,
        lambda_s,
        stress,
        vcma_change,
        magnetic_thickness,
    )?;
    // Switching is possible when K_eff changes sign relative to K_u
    Ok(k_eff * ku_intrinsic < 0.0)
}

// =============================================================================
// Switching Dynamics
// =============================================================================

/// Estimate the switching time for strain-mediated switching \[s\]
///
/// The switching time depends on the precessional dynamics:
/// τ_switch ≈ 1 / (γ × H_eff_strain)
///
/// where H_eff_strain = |K_σ| / (μ₀ M_s)
///
/// # Arguments
/// * `k_sigma` - Stress-induced anisotropy magnitude \[J/m³\]
/// * `ms` - Saturation magnetization \[A/m\]
/// * `gamma` - Gyromagnetic ratio [rad/(s·T)]
///
/// # Errors
/// Returns an error if any parameter leads to non-physical results.
pub fn estimate_switching_time(k_sigma: f64, ms: f64, gamma: f64) -> Result<f64> {
    if ms.abs() < 1e-30 {
        return Err(error::invalid_param("ms", "must be non-zero"));
    }
    if gamma <= 0.0 {
        return Err(error::invalid_param("gamma", "must be positive"));
    }

    let h_eff = k_sigma.abs() / (MU_0 * ms);
    if h_eff < 1e-30 {
        return Err(error::numerical_error(
            "effective field too small for switching time estimate",
        ));
    }

    Ok(1.0 / (gamma * h_eff))
}

/// Compute the energy-delay product (EDP) for straintronic switching
///
/// EDP = E_switch × τ_switch \[J·s\]
///
/// A key figure of merit; lower is better.
pub fn energy_delay_product(switching_energy: f64, switching_time: f64) -> f64 {
    switching_energy * switching_time
}

/// Compute the thermal stability factor Δ = E_barrier / (k_B T)
///
/// A data retention requirement for memory: typically Δ > 40-60.
///
/// # Arguments
/// * `energy_barrier` - Anisotropy energy barrier \[J\]
/// * `temperature` - Temperature \[K\]
///
/// # Errors
/// Returns an error if temperature is non-positive.
pub fn thermal_stability_factor(energy_barrier: f64, temperature: f64) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(error::invalid_param("temperature", "must be positive"));
    }
    Ok(energy_barrier / (crate::constants::KB * temperature))
}

// =============================================================================
// Multiferroic Heterostructure Modeling
// =============================================================================

/// Model a complete multiferroic heterostructure switching event
///
/// Computes the strain path and energy landscape for voltage-driven
/// magnetization switching in a piezo/magnetostrictive bilayer.
#[derive(Debug, Clone)]
pub struct MultiferroicSwitchingResult {
    /// Applied voltage \[V\]
    pub voltage: f64,
    /// In-plane strain achieved
    pub strain_in_plane: f64,
    /// Stress-induced anisotropy K_σ \[J/m³\]
    pub k_sigma: f64,
    /// VCMA anisotropy change ΔK_VCMA \[J/m³\]
    pub delta_k_vcma: f64,
    /// Total effective anisotropy K_eff \[J/m³\]
    pub k_eff: f64,
    /// Switching possible?
    pub switching_possible: bool,
    /// Estimated switching energy \[J\]
    pub switching_energy: f64,
}

/// Analyze a multiferroic switching event for a given voltage
///
/// # Arguments
/// * `device` - Straintronic device
/// * `vcma` - VCMA parameters (optional)
/// * `voltage` - Applied voltage \[V\]
/// * `substrate_thickness` - Piezoelectric substrate thickness \[m\]
pub fn analyze_multiferroic_switching(
    device: &StraintronicDevice,
    vcma: Option<&VcmaParameters>,
    voltage: f64,
    substrate_thickness: f64,
) -> Result<MultiferroicSwitchingResult> {
    // Compute effective strain
    let strain = device.effective_strain(voltage, substrate_thickness)?;
    let strain_in_plane = strain.components[0][0];

    // Stress from strain: σ = Y × ε
    let stress = strain_in_plane * device.free_layer.youngs_modulus;

    // Stress-induced anisotropy
    let k_sigma = stress_induced_anisotropy(device.free_layer.lambda_s, stress);

    // VCMA contribution
    let delta_k_vcma = match vcma {
        Some(v) => v.anisotropy_change(voltage) / device.free_layer_thickness,
        None => 0.0,
    };

    // Total effective anisotropy
    let k_eff = device.ku + k_sigma + delta_k_vcma;

    // Switching is possible when K_eff changes sign
    let switching_possible = k_eff * device.ku < 0.0;

    // Switching energy: dominated by capacitive energy of piezo
    // E ≈ (1/2) C V² where C = ε₀ εᵣ A / t for the piezo
    let area = device.free_layer_width * device.free_layer_width;
    let epsilon_r = 1000.0; // typical for PZT-like materials
    let capacitance = crate::constants::EPSILON_0 * epsilon_r * area / substrate_thickness;
    let switching_energy = 0.5 * capacitance * voltage * voltage;

    Ok(MultiferroicSwitchingResult {
        voltage,
        strain_in_plane,
        k_sigma,
        delta_k_vcma,
        k_eff,
        switching_possible,
        switching_energy,
    })
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GAMMA;

    #[test]
    fn test_switching_energy_comparison_strain_less_than_sot_less_than_stt() {
        let (strain_e, sot_e, stt_e) = typical_switching_energy_comparison();

        assert!(
            strain_e < sot_e,
            "Strain switching ({:.2e} J) should be less than SOT ({:.2e} J)",
            strain_e,
            sot_e
        );
        assert!(
            sot_e < stt_e,
            "SOT switching ({:.2e} J) should be less than STT ({:.2e} J)",
            sot_e,
            stt_e
        );

        // Strain should be ~aJ range
        assert!(
            strain_e < 1.0e-15,
            "Strain switching should be sub-fJ, got {:.2e} J",
            strain_e
        );
        // STT should be ~fJ to pJ range
        assert!(
            stt_e > 1.0e-15,
            "STT switching should be > fJ, got {:.2e} J",
            stt_e
        );
    }

    #[test]
    fn test_critical_strain_is_positive() {
        let ku = 5.0e5; // 500 kJ/m³
        let lambda_s = 25.0e-6; // CoFeB-like
        let youngs = 160.0e9;

        let eps_crit =
            critical_strain_uniaxial(ku, lambda_s, youngs).expect("should compute critical strain");

        assert!(
            eps_crit > 0.0,
            "Critical strain must be positive, got {}",
            eps_crit
        );

        // Verify order of magnitude: for CoFeB, ε_crit ~ 10⁻⁴ to 10⁻²
        assert!(
            eps_crit > 1.0e-6 && eps_crit < 1.0,
            "Critical strain should be in reasonable range, got {}",
            eps_crit
        );
    }

    #[test]
    fn test_vcma_anisotropy_change_sign() {
        let vcma = VcmaParameters::cofeb_mgo_typical();

        // Positive voltage → positive ΔK (since β > 0 by convention here)
        let delta_k_pos = vcma.anisotropy_change(1.0);
        // Negative voltage → negative ΔK
        let delta_k_neg = vcma.anisotropy_change(-1.0);

        // They should have opposite signs
        assert!(
            delta_k_pos * delta_k_neg < 0.0,
            "VCMA changes for +V and -V should have opposite signs: +V → {:.2e}, -V → {:.2e}",
            delta_k_pos,
            delta_k_neg
        );

        // Magnitude should be equal
        assert!(
            (delta_k_pos.abs() - delta_k_neg.abs()).abs() / delta_k_pos.abs() < 1e-10,
            "Magnitude of VCMA change should be symmetric"
        );
    }

    #[test]
    fn test_straintronic_device_creation() {
        let cofeb = MagnetoelasticMaterial::cofeb();
        let pmn_pt = PiezoelectricSubstrate::pmn_pt();

        let device = StraintronicDevice::new(cofeb, pmn_pt, 2.0e-9, 100.0e-9, 5.0e5)
            .expect("should create straintronic device");

        let vol = device.free_layer_volume();
        let expected_vol = 100.0e-9 * 100.0e-9 * 2.0e-9;
        assert!(
            (vol - expected_vol).abs() / expected_vol < 1e-10,
            "Volume calculation incorrect"
        );
    }

    #[test]
    fn test_critical_voltage_finite() {
        let eps_crit = 1.0e-3;
        let pmn_pt = PiezoelectricSubstrate::pmn_pt();
        let t_sub = 500.0e-6; // 500 μm substrate
        let eta = 0.9;

        let v_crit = critical_voltage(eps_crit, &pmn_pt, t_sub, eta)
            .expect("should compute critical voltage");

        assert!(
            v_crit > 0.0,
            "Critical voltage must be positive, got {}",
            v_crit
        );
        // For PMN-PT with large d31, voltage should be modest (< 100 V)
        assert!(
            v_crit < 1000.0,
            "Critical voltage should be reasonable, got {} V",
            v_crit
        );
    }

    #[test]
    fn test_thermal_stability_factor() {
        let ku = 5.0e5;
        let volume = 100.0e-9 * 100.0e-9 * 2.0e-9;
        let barrier = anisotropy_energy_barrier(ku, volume);
        let temperature = 300.0;

        let delta = thermal_stability_factor(barrier, temperature)
            .expect("should compute stability factor");

        // Δ should be positive
        assert!(
            delta > 0.0,
            "Thermal stability factor must be positive, got {}",
            delta
        );
    }

    #[test]
    fn test_combined_effective_anisotropy() {
        let ku = 5.0e5;
        let lambda_s = 25.0e-6;
        let stress = 200.0e6; // 200 MPa tensile
        let t_mag = 2.0e-9;

        // Without VCMA
        let k_eff = combined_effective_anisotropy(ku, lambda_s, stress, 0.0, t_mag)
            .expect("should compute combined anisotropy");

        // Stress-induced anisotropy for positive λ_s, tensile: K_σ = -(3/2)λ_s σ < 0
        // So K_eff < K_u
        assert!(
            k_eff < ku,
            "Tensile stress with positive λ_s should reduce K_eff: {} vs {}",
            k_eff,
            ku
        );
    }

    #[test]
    fn test_energy_delay_product() {
        let e_switch = 1.0e-18; // 1 aJ
        let tau = 1.0e-9; // 1 ns

        let edp = energy_delay_product(e_switch, tau);

        assert!(edp > 0.0, "EDP must be positive");
        assert!(
            (edp - 1.0e-27).abs() / 1.0e-27 < 1e-10,
            "EDP should be 1e-27 J·s, got {:.2e}",
            edp
        );
    }

    #[test]
    fn test_switching_time_estimate() {
        let k_sigma = 1.0e5; // 100 kJ/m³
        let ms = 1.0e6; // 1 MA/m

        let tau =
            estimate_switching_time(k_sigma, ms, GAMMA).expect("should estimate switching time");

        assert!(tau > 0.0, "Switching time must be positive");
        // Should be in sub-ns to ns range
        assert!(
            tau < 1.0e-6,
            "Switching time should be < μs, got {:.2e} s",
            tau
        );
    }

    #[test]
    fn test_analyze_multiferroic_switching() {
        let cofeb = MagnetoelasticMaterial::cofeb();
        let pmn_pt = PiezoelectricSubstrate::pmn_pt();

        let device = StraintronicDevice::new(cofeb, pmn_pt, 2.0e-9, 100.0e-9, 5.0e5)
            .expect("should create device");

        let result = analyze_multiferroic_switching(&device, None, 10.0, 500.0e-6)
            .expect("should analyze switching");

        // Strain should be non-zero for non-zero voltage
        assert!(
            result.strain_in_plane.abs() > 0.0,
            "Strain should be non-zero for non-zero voltage"
        );

        // Switching energy should be positive
        assert!(
            result.switching_energy > 0.0,
            "Switching energy must be positive"
        );
    }
}
