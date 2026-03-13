//! Spin Density Waves (SDW)
//!
//! This module implements the physics of spin density waves, which arise from
//! Fermi surface nesting in itinerant antiferromagnets. Chromium (Cr) is the
//! canonical SDW material.
//!
//! # Physical Background
//!
//! In metals with nested Fermi surfaces, an instability can drive the
//! formation of a static spin modulation:
//!
//! M(r) = M₀ cos(Q·r + φ)
//!
//! where Q is the nesting vector connecting parallel portions of the Fermi
//! surface. This is analogous to the BCS instability in superconductors but
//! in the spin channel.
//!
//! # Key Features
//!
//! - SDW order parameter with amplitude, nesting vector, and phase
//! - Mean-field gap equation (BCS-like self-consistency)
//! - Chromium-specific SDW properties (T_N = 311 K, spin-flip at 123 K)
//! - Transport anomaly at the Neel temperature
//! - SDW condensation and elastic energies
//!
//! # References
//!
//! - Fawcett, Rev. Mod. Phys. 60, 209 (1988) -- comprehensive review of Cr SDW
//! - Overhauser, Phys. Rev. 128, 1437 (1962) -- original SDW theory

use std::f64::consts::PI;

use crate::constants::KB;
use crate::error::{self, Result};
use crate::vector3::Vector3;

// ============================================================================
// SDW Order Parameter
// ============================================================================

/// Spin density wave state.
///
/// Represents a sinusoidal modulation of the spin density:
///
/// M(r) = M₀ cos(Q·r + φ)
///
/// The SDW is characterized by its amplitude, wavevector (nesting vector),
/// phase, energy gap, and transition temperature.
#[derive(Debug, Clone)]
pub struct SpinDensityWave {
    /// SDW amplitude M₀ \[A/m\]
    pub amplitude: f64,

    /// Nesting vector Q \[1/m\]
    ///
    /// The wavevector of the spin modulation, determined by Fermi surface
    /// geometry. For Cr: Q ≈ (2π/a)(1-δ, 0, 0) with δ ≈ 0.05.
    pub nesting_vector: Vector3<f64>,

    /// SDW phase φ \[rad\]
    pub phase: f64,

    /// Energy gap Δ \[eV\]
    ///
    /// The gap that opens at the Fermi surface due to the SDW formation.
    /// Related to the amplitude by Δ = U·M₀.
    pub gap: f64,

    /// Neel temperature T_N \[K\]
    ///
    /// Temperature above which the SDW order is destroyed.
    pub neel_temperature: f64,
}

impl SpinDensityWave {
    /// Create a new spin density wave state.
    ///
    /// # Arguments
    /// * `amplitude` - SDW amplitude M₀ \[A/m\]
    /// * `nesting_vector` - Nesting vector Q \[1/m\]
    /// * `phase` - SDW phase φ \[rad\]
    /// * `gap` - Energy gap Δ \[eV\]
    /// * `neel_temperature` - Neel temperature T_N \[K\]
    ///
    /// # Errors
    /// Returns error for non-physical parameters.
    pub fn new(
        amplitude: f64,
        nesting_vector: Vector3<f64>,
        phase: f64,
        gap: f64,
        neel_temperature: f64,
    ) -> Result<Self> {
        if amplitude < 0.0 {
            return Err(error::invalid_param("amplitude", "must be non-negative"));
        }
        if gap < 0.0 {
            return Err(error::invalid_param("gap", "must be non-negative"));
        }
        if neel_temperature <= 0.0 {
            return Err(error::invalid_param("neel_temperature", "must be positive"));
        }

        Ok(Self {
            amplitude,
            nesting_vector,
            phase,
            gap,
            neel_temperature,
        })
    }

    /// Create a Chromium SDW with standard parameters.
    ///
    /// Cr is the canonical incommensurate SDW material:
    /// - T_N = 311 K
    /// - Q ≈ (2π/a)(0.95, 0, 0) with a = 2.88 Å
    /// - Δ ≈ 0.12 eV at T = 0
    /// - M₀ ≈ 0.62 μ_B per atom
    pub fn chromium() -> Result<Self> {
        let lattice_constant = 2.88e-10; // Cr lattice constant [m]
        let delta = 0.05; // incommensurability parameter

        // Nesting vector Q = (2π/a)(1 - δ) along x
        let q_magnitude = 2.0 * PI / lattice_constant * (1.0 - delta);
        let nesting_vector = Vector3::new(q_magnitude, 0.0, 0.0);

        // SDW amplitude: ~0.62 μ_B per atom → convert to A/m
        // M₀ ≈ 0.62 μ_B × n_atoms, where n_atoms = 2/a³ for BCC Cr
        // For bulk: M₀ ≈ 5.0e4 A/m (approximate)
        let amplitude = 5.0e4; // A/m

        Self::new(
            amplitude,
            nesting_vector,
            0.0,   // phase = 0
            0.12,  // gap in eV
            311.0, // T_N = 311 K
        )
    }

    /// Evaluate the spin density at position r.
    ///
    /// M(r) = M₀ cos(Q·r + φ)
    ///
    /// Returns the local magnetization magnitude (can be negative for
    /// antiferromagnetic regions).
    pub fn magnetization_at(&self, r: &Vector3<f64>) -> f64 {
        let q_dot_r = self.nesting_vector.dot(r);
        self.amplitude * (q_dot_r + self.phase).cos()
    }

    /// Calculate the SDW wavelength λ = 2π / |Q|.
    pub fn wavelength(&self) -> f64 {
        let q_mag = self.nesting_vector.dot(&self.nesting_vector).sqrt();
        if q_mag > 0.0 {
            2.0 * PI / q_mag
        } else {
            f64::INFINITY
        }
    }

    /// Check if the SDW is commensurate with the lattice.
    ///
    /// An SDW is commensurate if Q = (2π/a)(p/q, 0, 0) where p, q are
    /// small integers. We check if Q·a/(2π) is close to a rational number.
    ///
    /// # Arguments
    /// * `lattice_constant` - Lattice parameter a \[m\]
    /// * `tolerance` - How close to rational to consider commensurate
    pub fn is_commensurate(&self, lattice_constant: f64, tolerance: f64) -> bool {
        let q_mag = self.nesting_vector.dot(&self.nesting_vector).sqrt();
        let q_ratio = q_mag * lattice_constant / (2.0 * PI);

        // Check if close to a simple rational p/q with q ≤ 4
        for denom in 1..=4 {
            let numer = (q_ratio * denom as f64).round();
            let rational = numer / denom as f64;
            if (q_ratio - rational).abs() < tolerance {
                return true;
            }
        }
        false
    }
}

// ============================================================================
// SDW Gap Equation
// ============================================================================

/// Mean-field SDW gap equation solver.
///
/// The SDW gap follows a BCS-like self-consistency equation:
///
/// 1 = U · N(0) · ∫₀^{ω_D} dε tanh(√(ε² + Δ²)/(2k_BT)) / √(ε² + Δ²)
///
/// where U is the electron-electron interaction, N(0) is the density of
/// states at the Fermi level, and ω_D is an energy cutoff.
#[derive(Debug, Clone)]
pub struct SdwGapSolver {
    /// Electron-electron interaction strength U·N(0) (dimensionless coupling)
    pub coupling: f64,

    /// Energy cutoff (Debye-like) \[eV\]
    pub cutoff_energy: f64,

    /// Neel temperature \[K\]
    pub neel_temperature: f64,

    /// Number of integration points for numerical quadrature
    pub n_points: usize,
}

impl SdwGapSolver {
    /// Create a new SDW gap solver.
    ///
    /// # Arguments
    /// * `coupling` - Dimensionless coupling U·N(0)
    /// * `cutoff_energy` - Energy cutoff \[eV\]
    /// * `neel_temperature` - Neel temperature \[K\]
    ///
    /// # Errors
    /// Returns error for non-physical parameters.
    pub fn new(coupling: f64, cutoff_energy: f64, neel_temperature: f64) -> Result<Self> {
        if coupling <= 0.0 {
            return Err(error::invalid_param("coupling", "must be positive"));
        }
        if cutoff_energy <= 0.0 {
            return Err(error::invalid_param("cutoff_energy", "must be positive"));
        }
        if neel_temperature <= 0.0 {
            return Err(error::invalid_param("neel_temperature", "must be positive"));
        }

        Ok(Self {
            coupling,
            cutoff_energy,
            neel_temperature,
            n_points: 200,
        })
    }

    /// Estimate the zero-temperature gap from BCS-like weak-coupling formula.
    ///
    /// Δ(0) ≈ 2 ω_D exp(-1/λ)
    ///
    /// where λ = U·N(0) is the dimensionless coupling.
    pub fn zero_temperature_gap(&self) -> f64 {
        2.0 * self.cutoff_energy * (-1.0 / self.coupling).exp()
    }

    /// Calculate the SDW gap at a given temperature using BCS-like temperature
    /// dependence.
    ///
    /// Near T_N, the gap vanishes as:
    ///   Δ(T) ≈ Δ(0) · √(1 - (T/T_N))   (mean-field)
    ///
    /// More precisely, we use a smooth interpolation that captures:
    /// - Δ(0) = full gap at T = 0
    /// - Δ(T_N) = 0 at the Neel temperature
    /// - BCS-like temperature dependence
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]
    ///
    /// # Returns
    /// The gap Δ(T) in eV.
    pub fn gap_at_temperature(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return self.zero_temperature_gap();
        }
        if temperature >= self.neel_temperature {
            return 0.0;
        }

        let t_ratio = temperature / self.neel_temperature;
        let delta_0 = self.zero_temperature_gap();

        // BCS-like interpolation: Δ(T) = Δ(0) · tanh(α √(T_N/T - 1))
        // where α ≈ 1.74 (from BCS theory, 2Δ(0)/(k_B T_c) ≈ 3.53)
        let alpha = 1.74;
        let arg = alpha * (1.0 / t_ratio - 1.0).sqrt();

        delta_0 * arg.tanh()
    }

    /// Self-consistent gap equation evaluation.
    ///
    /// Computes the right-hand side of:
    ///
    /// 1/λ = ∫₀^{ω_D} dε tanh(√(ε² + Δ²)/(2k_BT)) / √(ε² + Δ²)
    ///
    /// # Arguments
    /// * `gap` - Trial gap value \[eV\]
    /// * `temperature` - Temperature \[K\]
    ///
    /// # Returns
    /// The integrated kernel value. Self-consistency requires this equals 1/λ.
    pub fn gap_equation_rhs(&self, gap: f64, temperature: f64) -> f64 {
        let n = self.n_points;
        let de = self.cutoff_energy / n as f64;
        let mut integral = 0.0;

        // Convert temperature to eV for consistency
        let k_b_t_ev = KB * temperature / (1.602_176_634e-19); // K → eV

        for i in 0..n {
            // Use midpoint rule, avoiding ε = 0 singularity when gap = 0
            let epsilon = (i as f64 + 0.5) * de;
            let e_k = (epsilon * epsilon + gap * gap).sqrt();

            let tanh_arg = if k_b_t_ev > 1e-30 {
                e_k / (2.0 * k_b_t_ev)
            } else {
                // T → 0 limit: tanh → 1
                f64::MAX
            };

            let tanh_val = if tanh_arg > 30.0 {
                1.0
            } else {
                tanh_arg.tanh()
            };

            integral += tanh_val / e_k;
        }

        integral * de
    }

    /// Solve the gap equation self-consistently at a given temperature.
    ///
    /// Uses bisection to find Δ such that the gap equation is satisfied.
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]
    /// * `max_iterations` - Maximum number of bisection iterations
    ///
    /// # Returns
    /// Self-consistent gap Δ(T) in eV.
    ///
    /// # Errors
    /// Returns error if convergence fails.
    pub fn solve_gap(&self, temperature: f64, max_iterations: usize) -> Result<f64> {
        if temperature >= self.neel_temperature {
            return Ok(0.0);
        }

        let target = 1.0 / self.coupling;
        let delta_0 = self.zero_temperature_gap();

        // Bisection: find Δ where gap_equation_rhs(Δ, T) = 1/λ
        let mut lo = 0.0_f64;
        let mut hi = delta_0 * 1.5; // generous upper bound

        for _ in 0..max_iterations {
            let mid = (lo + hi) / 2.0;
            let rhs = self.gap_equation_rhs(mid, temperature);

            if (rhs - target).abs() < 1e-8 {
                return Ok(mid);
            }

            // The RHS is a decreasing function of Δ
            if rhs > target {
                lo = mid;
            } else {
                hi = mid;
            }

            if (hi - lo) < 1e-14 * delta_0 {
                return Ok((lo + hi) / 2.0);
            }
        }

        Ok((lo + hi) / 2.0)
    }
}

// ============================================================================
// Chromium SDW Properties
// ============================================================================

/// Chromium-specific SDW parameters and properties.
///
/// Cr is the archetypal itinerant antiferromagnet with an incommensurate
/// SDW below T_N = 311 K and a spin-flip transition at T_SF = 123 K.
pub struct ChromiumSdw;

impl ChromiumSdw {
    /// Neel temperature of Cr \[K\]
    pub const NEEL_TEMPERATURE: f64 = 311.0;

    /// Spin-flip transition temperature \[K\]
    ///
    /// Below T_SF, the SDW polarization rotates from transverse to longitudinal
    pub const SPIN_FLIP_TEMPERATURE: f64 = 123.0;

    /// Lattice constant of BCC Cr \[m\]
    pub const LATTICE_CONSTANT: f64 = 2.88e-10;

    /// Incommensurability parameter δ
    ///
    /// Q = (2π/a)(1 - δ, 0, 0)
    pub const INCOMMENSURABILITY: f64 = 0.05;

    /// SDW gap at T = 0 \[eV\]
    pub const ZERO_TEMP_GAP: f64 = 0.12;

    /// Calculate the Cr nesting vector.
    ///
    /// Q = (2π/a)(1 - δ) x̂
    pub fn nesting_vector() -> Vector3<f64> {
        let q = 2.0 * PI / Self::LATTICE_CONSTANT * (1.0 - Self::INCOMMENSURABILITY);
        Vector3::new(q, 0.0, 0.0)
    }

    /// Calculate the nesting vector in units of 2π/a.
    ///
    /// Returns Q / (2π/a) which should be approximately (0.95, 0, 0).
    pub fn nesting_vector_reduced() -> Vector3<f64> {
        let factor = 1.0 - Self::INCOMMENSURABILITY;
        Vector3::new(factor, 0.0, 0.0)
    }

    /// Determine SDW polarization at a given temperature.
    ///
    /// - T < T_SF: longitudinal (M ∥ Q) -- spin-flip phase
    /// - T_SF < T < T_N: transverse (M ⊥ Q)
    /// - T > T_N: paramagnetic (no SDW)
    pub fn polarization_type(temperature: f64) -> SdwPolarization {
        if temperature > Self::NEEL_TEMPERATURE {
            SdwPolarization::None
        } else if temperature > Self::SPIN_FLIP_TEMPERATURE {
            SdwPolarization::Transverse
        } else {
            SdwPolarization::Longitudinal
        }
    }
}

/// SDW polarization relative to the nesting vector Q.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SdwPolarization {
    /// No SDW (paramagnetic state above T_N)
    None,
    /// Transverse: M ⊥ Q (between T_SF and T_N in Cr)
    Transverse,
    /// Longitudinal: M ∥ Q (below T_SF in Cr)
    Longitudinal,
}

// ============================================================================
// Transport in SDW State
// ============================================================================

/// Calculate the resistivity anomaly at the SDW transition.
///
/// The SDW gap removes states from the Fermi surface, reducing the
/// density of carriers and increasing resistivity. The anomalous
/// contribution to resistivity is approximately:
///
/// Δρ/ρ_0 ≈ -f · (Δ(T)/(k_B T_N))²
///
/// where f is the fraction of the Fermi surface that is gapped.
///
/// # Arguments
/// * `gap_ev` - SDW gap at the given temperature \[eV\]
/// * `neel_temperature` - Neel temperature \[K\]
/// * `fermi_surface_fraction` - Fraction of Fermi surface affected by nesting (0 to 1)
///
/// # Returns
/// Relative resistivity change Δρ/ρ₀ (negative means conductivity decreases,
/// i.e. resistivity increases due to reduced scattering phase space, or
/// in the SDW state the gapping leads to anomalous behavior).
///
/// # Errors
/// Returns error for non-physical parameters.
pub fn resistivity_anomaly(
    gap_ev: f64,
    neel_temperature: f64,
    fermi_surface_fraction: f64,
) -> Result<f64> {
    if neel_temperature <= 0.0 {
        return Err(error::invalid_param("neel_temperature", "must be positive"));
    }
    if !(0.0..=1.0).contains(&fermi_surface_fraction) {
        return Err(error::invalid_param(
            "fermi_surface_fraction",
            "must be between 0 and 1",
        ));
    }

    // k_B T_N in eV
    let kb_tn_ev = KB * neel_temperature / 1.602_176_634e-19;

    if kb_tn_ev < 1e-30 {
        return Ok(0.0);
    }

    // The gap opening removes carriers → resistivity increase in some channels,
    // but the net effect is a kink/anomaly. Using a simplified model:
    let ratio = gap_ev / kb_tn_ev;
    let delta_rho = -fermi_surface_fraction * ratio * ratio;

    Ok(delta_rho)
}

// ============================================================================
// SDW Energy
// ============================================================================

/// Calculate the SDW condensation energy density.
///
/// The condensation energy is the energy gained by forming the SDW state:
///
/// E_cond = -½ N(0) Δ²
///
/// where N(0) is the density of states at the Fermi level and Δ is the gap.
/// This is always negative (energetically favorable).
///
/// # Arguments
/// * `density_of_states` - N(0), density of states at Fermi level [states/(eV·m³)]
/// * `gap_ev` - SDW gap Δ \[eV\]
///
/// # Returns
/// Condensation energy density \[J/m³\] (negative = favorable).
///
/// # Errors
/// Returns error for non-physical parameters.
pub fn condensation_energy(density_of_states: f64, gap_ev: f64) -> Result<f64> {
    if density_of_states <= 0.0 {
        return Err(error::invalid_param(
            "density_of_states",
            "must be positive",
        ));
    }
    if gap_ev < 0.0 {
        return Err(error::invalid_param("gap_ev", "must be non-negative"));
    }

    let ev_to_j = 1.602_176_634e-19;
    // E_cond = -½ N(0) Δ², convert Δ from eV to J, N(0) from states/(eV·m³) to states/(J·m³)
    // N(0) [states/(eV·m³)] × 1/eV_to_J = N(0) [states/(J·m³)]
    // Δ² [eV²] × eV_to_J² = Δ² [J²]
    // Result = -½ × N(0)/eV_to_J × Δ²×eV_to_J² = -½ × N(0) × Δ² × eV_to_J
    let e_cond = -0.5 * density_of_states * gap_ev * gap_ev * ev_to_j;

    Ok(e_cond)
}

/// Calculate the elastic energy associated with lattice distortion in the SDW state.
///
/// The SDW can couple to the lattice via magnetoelastic effects, producing
/// a strain wave at wavevector 2Q (for a sinusoidal SDW). The elastic
/// energy density is:
///
/// E_elastic = ½ C ε²
///
/// where C is the elastic modulus and ε is the strain amplitude.
///
/// # Arguments
/// * `elastic_modulus` - Elastic modulus C \[Pa\]
/// * `strain_amplitude` - Lattice strain ε (dimensionless)
///
/// # Returns
/// Elastic energy density \[J/m³\] (positive).
///
/// # Errors
/// Returns error for non-physical parameters.
pub fn elastic_energy(elastic_modulus: f64, strain_amplitude: f64) -> Result<f64> {
    if elastic_modulus <= 0.0 {
        return Err(error::invalid_param("elastic_modulus", "must be positive"));
    }

    Ok(0.5 * elastic_modulus * strain_amplitude * strain_amplitude)
}

/// Calculate the total SDW energy (condensation + elastic) per unit volume.
///
/// # Arguments
/// * `density_of_states` - N(0) [states/(eV·m³)]
/// * `gap_ev` - SDW gap \[eV\]
/// * `elastic_modulus` - Elastic modulus \[Pa\]
/// * `strain_amplitude` - Strain amplitude (dimensionless)
///
/// # Returns
/// Total energy density \[J/m³\].
pub fn total_sdw_energy(
    density_of_states: f64,
    gap_ev: f64,
    elastic_modulus: f64,
    strain_amplitude: f64,
) -> Result<f64> {
    let e_cond = condensation_energy(density_of_states, gap_ev)?;
    let e_elastic = elastic_energy(elastic_modulus, strain_amplitude)?;
    Ok(e_cond + e_elastic)
}

// ============================================================================
// SDW Susceptibility
// ============================================================================

/// Calculate the static spin susceptibility enhancement near the SDW instability.
///
/// The Lindhard susceptibility diverges at the nesting vector Q as T → T_N:
///
/// χ(Q, T) = χ₀ / (1 - U·χ₀)
///
/// where χ₀ = N(0) · ln(1.13 ω_D / (k_B T)) is the bare susceptibility
/// (Lindhard function at perfect nesting) and U is the interaction.
///
/// Near T_N this diverges as ~1/(T - T_N), signaling the SDW instability.
///
/// # Arguments
/// * `density_of_states` - N(0) [states/(eV·m³)]
/// * `interaction_u` - Electron-electron interaction U \[eV\]
/// * `cutoff_energy` - Energy cutoff ω_D \[eV\]
/// * `temperature` - Temperature \[K\]
///
/// # Returns
/// Enhanced susceptibility [states/(eV·m³)].
///
/// # Errors
/// Returns error at or very near the divergence.
pub fn enhanced_susceptibility(
    density_of_states: f64,
    interaction_u: f64,
    cutoff_energy: f64,
    temperature: f64,
) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(error::invalid_param("temperature", "must be positive"));
    }
    if density_of_states <= 0.0 {
        return Err(error::invalid_param(
            "density_of_states",
            "must be positive",
        ));
    }

    let kb_t_ev = KB * temperature / 1.602_176_634e-19;
    if kb_t_ev < 1e-30 {
        return Err(error::numerical_error(
            "temperature too low for susceptibility calculation",
        ));
    }

    let log_arg = 1.13 * cutoff_energy / kb_t_ev;
    if log_arg <= 0.0 {
        return Err(error::numerical_error(
            "invalid argument for logarithm in susceptibility",
        ));
    }

    let chi_0 = density_of_states * log_arg.ln();
    let denominator = 1.0 - interaction_u * chi_0;

    if denominator.abs() < 1e-10 {
        return Err(error::numerical_error(
            "susceptibility diverges at SDW transition",
        ));
    }

    Ok(chi_0 / denominator)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sdw_gap_vanishes_at_neel_temperature() {
        let solver = SdwGapSolver::new(0.5, 1.0, 311.0).expect("SdwGapSolver::new should succeed");

        let gap_at_tn = solver.gap_at_temperature(311.0);
        assert!(
            gap_at_tn.abs() < 1e-15,
            "SDW gap must vanish at T_N, got {}",
            gap_at_tn
        );

        // Also check above T_N
        let gap_above = solver.gap_at_temperature(400.0);
        assert!(
            gap_above.abs() < 1e-15,
            "SDW gap must be zero above T_N, got {}",
            gap_above
        );
    }

    #[test]
    fn test_cr_nesting_vector_approximately_095() {
        let q_reduced = ChromiumSdw::nesting_vector_reduced();

        assert!(
            (q_reduced.x - 0.95).abs() < 0.001,
            "Cr nesting vector x-component should be ~0.95 in units of 2π/a, got {}",
            q_reduced.x
        );
        assert!(
            q_reduced.y.abs() < 1e-15,
            "Cr nesting vector y-component should be 0"
        );
        assert!(
            q_reduced.z.abs() < 1e-15,
            "Cr nesting vector z-component should be 0"
        );
    }

    #[test]
    fn test_sdw_condensation_energy_is_negative() {
        // Condensation energy should be negative (favorable)
        let n0 = 1e28; // states/(eV·m³), typical for transition metals
        let gap = 0.12; // eV

        let e_cond = condensation_energy(n0, gap).expect("condensation_energy should succeed");
        assert!(
            e_cond < 0.0,
            "SDW condensation energy must be negative (favorable), got {}",
            e_cond
        );
    }

    #[test]
    fn test_sdw_magnetization_oscillates() {
        let sdw = SpinDensityWave::chromium().expect("Chromium SDW creation should succeed");

        let wavelength = sdw.wavelength();
        assert!(wavelength > 0.0, "SDW wavelength must be positive");

        // Check that magnetization oscillates: sample at 0, λ/4, λ/2
        let r0 = Vector3::new(0.0, 0.0, 0.0);
        let r_quarter = Vector3::new(wavelength / 4.0, 0.0, 0.0);
        let r_half = Vector3::new(wavelength / 2.0, 0.0, 0.0);

        let m0 = sdw.magnetization_at(&r0);
        let m_quarter = sdw.magnetization_at(&r_quarter);
        let m_half = sdw.magnetization_at(&r_half);

        // At r=0 (phase=0): M = M₀
        assert!(
            (m0 - sdw.amplitude).abs() < sdw.amplitude * 1e-10,
            "M(0) should equal M₀"
        );

        // At r=λ/4: M ≈ 0 (cos(π/2) = 0)
        assert!(
            m_quarter.abs() < sdw.amplitude * 1e-10,
            "M(λ/4) should be approximately 0, got {}",
            m_quarter
        );

        // At r=λ/2: M = -M₀
        assert!(
            (m_half + sdw.amplitude).abs() < sdw.amplitude * 1e-10,
            "M(λ/2) should equal -M₀"
        );
    }

    #[test]
    fn test_chromium_sdw_parameters() {
        let sdw = SpinDensityWave::chromium().expect("Chromium SDW creation should succeed");

        assert!(
            (sdw.neel_temperature - 311.0).abs() < 1e-10,
            "Cr T_N should be 311 K"
        );
        assert!((sdw.gap - 0.12).abs() < 1e-10, "Cr gap should be 0.12 eV");
        assert!(sdw.amplitude > 0.0, "Cr SDW amplitude should be positive");
    }

    #[test]
    fn test_chromium_polarization_phases() {
        // T < T_SF = 123 K → longitudinal
        assert_eq!(
            ChromiumSdw::polarization_type(100.0),
            SdwPolarization::Longitudinal,
            "Below T_SF should be longitudinal"
        );

        // T_SF < T < T_N → transverse
        assert_eq!(
            ChromiumSdw::polarization_type(200.0),
            SdwPolarization::Transverse,
            "Between T_SF and T_N should be transverse"
        );

        // T > T_N → none
        assert_eq!(
            ChromiumSdw::polarization_type(400.0),
            SdwPolarization::None,
            "Above T_N should be paramagnetic (no SDW)"
        );
    }

    #[test]
    fn test_resistivity_anomaly_sign() {
        let delta_rho =
            resistivity_anomaly(0.12, 311.0, 0.3).expect("resistivity_anomaly should succeed");

        // The anomaly should be negative in our model
        assert!(
            delta_rho < 0.0,
            "Resistivity anomaly should be negative (gap reduces scattering channels)"
        );
    }

    #[test]
    fn test_resistivity_anomaly_zero_gap() {
        let delta_rho =
            resistivity_anomaly(0.0, 311.0, 0.3).expect("resistivity_anomaly should succeed");

        assert!(
            delta_rho.abs() < 1e-15,
            "Resistivity anomaly should be zero when gap is zero"
        );
    }

    #[test]
    fn test_elastic_energy_positive() {
        let e_el = elastic_energy(3.5e11, 1e-5).expect("elastic_energy should succeed");
        assert!(e_el > 0.0, "Elastic energy should be positive");
    }

    #[test]
    fn test_total_sdw_energy() {
        let n0 = 1e28;
        let gap = 0.12;
        let c_modulus = 3.5e11;
        let strain = 1e-6; // very small strain

        let e_total =
            total_sdw_energy(n0, gap, c_modulus, strain).expect("total_sdw_energy should succeed");

        // Condensation energy dominates for small strain → total should be negative
        let e_cond = condensation_energy(n0, gap).expect("condensation_energy should succeed");
        let e_el = elastic_energy(c_modulus, strain).expect("elastic_energy should succeed");

        assert!(
            (e_total - (e_cond + e_el)).abs() < 1e-30,
            "Total energy should be sum of condensation and elastic"
        );

        // For small strain, condensation should dominate
        assert!(
            e_total < 0.0,
            "Total SDW energy should be negative for small strain"
        );
    }

    #[test]
    fn test_sdw_gap_temperature_dependence() {
        let solver = SdwGapSolver::new(0.5, 1.0, 311.0).expect("SdwGapSolver::new should succeed");

        let gap_0 = solver.gap_at_temperature(0.0);
        let gap_150 = solver.gap_at_temperature(150.0);
        let gap_300 = solver.gap_at_temperature(300.0);

        // Gap should decrease with temperature
        assert!(
            gap_0 > gap_150,
            "Gap at T=0 should be larger than at T=150K"
        );
        assert!(
            gap_150 > gap_300,
            "Gap at T=150K should be larger than at T=300K"
        );
        assert!(gap_0 > 0.0, "Zero-temperature gap should be positive");
    }

    #[test]
    fn test_sdw_incommensurability() {
        let sdw = SpinDensityWave::chromium().expect("Chromium SDW creation should succeed");

        // Cr SDW is incommensurate
        let is_comm = sdw.is_commensurate(ChromiumSdw::LATTICE_CONSTANT, 0.01);
        assert!(!is_comm, "Cr SDW should be incommensurate (δ = 0.05)");

        // A commensurate SDW with Q = 2π/a should register as commensurate
        let comm_sdw = SpinDensityWave::new(
            1e4,
            Vector3::new(2.0 * PI / 3e-10, 0.0, 0.0),
            0.0,
            0.1,
            300.0,
        )
        .expect("commensurate SDW creation should succeed");
        let is_comm2 = comm_sdw.is_commensurate(3e-10, 0.01);
        assert!(is_comm2, "SDW with Q = 2π/a should be commensurate");
    }

    #[test]
    fn test_invalid_sdw_parameters() {
        // Negative amplitude
        assert!(SpinDensityWave::new(-1.0, Vector3::new(1.0, 0.0, 0.0), 0.0, 0.1, 300.0).is_err());

        // Negative gap
        assert!(SpinDensityWave::new(1.0, Vector3::new(1.0, 0.0, 0.0), 0.0, -0.1, 300.0).is_err());

        // Negative T_N
        assert!(SpinDensityWave::new(1.0, Vector3::new(1.0, 0.0, 0.0), 0.0, 0.1, -1.0).is_err());

        // Invalid gap solver parameters
        assert!(SdwGapSolver::new(0.0, 1.0, 311.0).is_err());
        assert!(SdwGapSolver::new(0.5, 0.0, 311.0).is_err());
    }
}
