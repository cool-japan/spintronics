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
//! - Time-dependent Ginzburg-Landau (TDGL) relaxational dynamics of the
//!   order-parameter amplitude toward the self-consistent gap
//!
//! # References
//!
//! - Fawcett, Rev. Mod. Phys. 60, 209 (1988) -- comprehensive review of Cr SDW
//! - Overhauser, Phys. Rev. 128, 1437 (1962) -- original SDW theory
//! - Hohenberg & Halperin, Rev. Mod. Phys. 49, 435 (1977) -- Model A
//!   (non-conserved, purely relaxational) order-parameter dynamics

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
// TDGL Relaxational Dynamics
// ============================================================================

/// Time-dependent Ginzburg-Landau (TDGL) relaxational dynamics for the SDW
/// order-parameter amplitude.
///
/// # Physical background
///
/// The machinery above ([`SdwGapSolver`], [`condensation_energy`],
/// [`elastic_energy`]) describes the *equilibrium* SDW state at a given
/// temperature but says nothing about how the order parameter reaches that
/// equilibrium dynamically -- e.g. relaxing after a quench, or approaching
/// equilibrium from an arbitrary initial condition. This is supplied by the
/// standard "Model A" (non-conserved, purely dissipative) time-dependent
/// Ginzburg-Landau equation of motion for the gap amplitude Δ (the modulus
/// of the complex SDW order parameter Δe^{iφ}, in the same sense that
/// [`SpinDensityWave::gap`] plays this role):
///
/// dΔ/dt = -Γ (∂F/∂Δ)
///
/// where Γ > 0 (`relaxation_rate`) is a phenomenological kinetic coefficient
/// and F(Δ, T) is a Landau free energy functional whose *minimum*
/// reproduces the self-consistent BCS-like gap already computed by
/// [`SdwGapSolver::gap_at_temperature`].
///
/// # Free energy construction
///
/// F is the standard quartic (φ⁴) Landau functional
///
/// F(Δ, T) = (b/4) (Δ² - Δ_eq(T)²)²
///
/// where Δ_eq(T) = `gap_solver.gap_at_temperature(T)` is the pre-existing
/// self-consistent equilibrium amplitude and b = `quartic_stiffness` > 0 is a
/// fixed Landau quartic coefficient. Expanding the square gives exactly the
/// textbook a(T)Δ²/2 + bΔ⁴/4 Landau expansion with quadratic coefficient
/// a(T) = -b·Δ_eq(T)², i.e. the standard way a phenomenological Landau
/// theory is calibrated against a microscopic (here BCS mean-field) result:
/// a(T) < 0 for T < T_N, so Δ = 0 is an unstable local *maximum* (the
/// paramagnetic state is unstable to SDW formation, as physically
/// required), while a(T) = 0 identically once T ≥ T_N (Δ_eq = 0), leaving
/// the pure positive quartic well (b/4)Δ⁴ whose unique minimum is at Δ = 0
/// (no SDW above the Neel temperature).
///
/// The quartic term is implemented via the module's existing
/// [`elastic_energy`] harmonic-energy formula applied to the auxiliary
/// "generalized strain" x = Δ² - Δ_eq(T)² (a standard Landau-theory device:
/// since F must be even in Δ by the Δ → -Δ symmetry of the order parameter,
/// any smooth F is naturally a harmonic function of Δ² near its minimum):
///
/// F(Δ, T) = `elastic_energy`(b/2, Δ² - Δ_eq(T)²)
///
/// The physical *condensation* energy gained at a given Δ (reusing
/// [`condensation_energy`] directly) is tracked separately, via
/// [`SdwRelaxationDynamics::condensation_energy_at`], as a diagnostic. It is
/// **not** summed into F above: doing so would shift F's minimum away from
/// Δ_eq(T) by a T-independent offset (adding a raw -½N(0)Δ² term moves the
/// stationary point to Δ² = Δ_eq(T)² + N(0)/(2b) rather than Δ_eq(T)²
/// itself), which would silently detune the relaxation target away from the
/// self-consistent gap. Keeping it as a separate report avoids that.
///
/// # References
///
/// - Standard Ginzburg-Landau theory of a continuous order-parameter
///   transition, e.g. Chaikin & Lubensky, "Principles of Condensed Matter
///   Physics", ch. 4.
/// - Hohenberg & Halperin, Rev. Mod. Phys. 49, 435 (1977) -- Model A
///   dynamics.
#[derive(Debug, Clone)]
pub struct SdwRelaxationDynamics {
    /// Self-consistent BCS-like gap solver providing the equilibrium target
    /// Δ_eq(T) = `gap_solver.gap_at_temperature(T)`.
    pub gap_solver: SdwGapSolver,

    /// Density of states at the Fermi level N(0) \[states/(eV·m³)\], used
    /// only for the [`condensation_energy_at`](Self::condensation_energy_at)
    /// diagnostic. Does not affect the relaxation dynamics itself.
    pub density_of_states: f64,

    /// Landau quartic stiffness b \[eV⁻²\] controlling the curvature of the
    /// free-energy well around its minimum, and hence (together with
    /// `relaxation_rate`) the relaxation timescale.
    pub quartic_stiffness: f64,

    /// TDGL relaxation rate Γ \[s⁻¹\].
    pub relaxation_rate: f64,
}

impl SdwRelaxationDynamics {
    /// Create new TDGL relaxational dynamics for an SDW order parameter.
    ///
    /// # Arguments
    /// * `gap_solver` - self-consistent BCS-like gap solver (supplies the
    ///   equilibrium target Δ_eq(T))
    /// * `density_of_states` - N(0) \[states/(eV·m³)\], used only for the
    ///   condensation-energy diagnostic
    /// * `quartic_stiffness` - Landau quartic coefficient b \[eV⁻²\], must be
    ///   positive
    /// * `relaxation_rate` - TDGL relaxation rate Γ \[s⁻¹\], must be positive
    ///
    /// # Errors
    /// Returns error for non-physical (non-positive) parameters.
    pub fn new(
        gap_solver: SdwGapSolver,
        density_of_states: f64,
        quartic_stiffness: f64,
        relaxation_rate: f64,
    ) -> Result<Self> {
        if density_of_states <= 0.0 {
            return Err(error::invalid_param(
                "density_of_states",
                "must be positive",
            ));
        }
        if quartic_stiffness <= 0.0 {
            return Err(error::invalid_param(
                "quartic_stiffness",
                "must be positive",
            ));
        }
        if relaxation_rate <= 0.0 {
            return Err(error::invalid_param("relaxation_rate", "must be positive"));
        }

        Ok(Self {
            gap_solver,
            density_of_states,
            quartic_stiffness,
            relaxation_rate,
        })
    }

    /// Self-consistent equilibrium amplitude Δ_eq(T) \[eV\], i.e. the fixed
    /// point that the relaxation dynamics converges to.
    ///
    /// Thin convenience wrapper around
    /// `self.gap_solver.gap_at_temperature(temperature)`.
    pub fn equilibrium_amplitude(&self, temperature: f64) -> f64 {
        self.gap_solver.gap_at_temperature(temperature)
    }

    /// Landau free energy density F(Δ, T) = (b/4)(Δ² - Δ_eq(T)²)² \[eV²\]
    /// (reduced units; see the struct-level documentation).
    ///
    /// This is the functional whose gradient drives the relaxational
    /// dynamics in [`relax_to_equilibrium`](Self::relax_to_equilibrium). It
    /// has a single global minimum F = 0 at Δ = Δ_eq(T) and, for T < T_N, an
    /// unstable local maximum at Δ = 0.
    ///
    /// # Errors
    /// Propagates errors from the underlying [`elastic_energy`] call (only
    /// possible if `quartic_stiffness` were non-positive, which is already
    /// excluded by [`SdwRelaxationDynamics::new`]).
    pub fn free_energy_density(&self, amplitude: f64, temperature: f64) -> Result<f64> {
        let eq = self.equilibrium_amplitude(temperature);
        let strain = amplitude * amplitude - eq * eq;
        elastic_energy(0.5 * self.quartic_stiffness, strain)
    }

    /// Gradient ∂F/∂Δ = b·Δ·(Δ² - Δ_eq(T)²) \[eV\] of the Landau free energy
    /// with respect to the amplitude Δ.
    pub fn free_energy_gradient(&self, amplitude: f64, temperature: f64) -> f64 {
        let eq = self.equilibrium_amplitude(temperature);
        self.quartic_stiffness * amplitude * (amplitude * amplitude - eq * eq)
    }

    /// Condensation energy density gained at amplitude Δ \[J/m³\], reusing
    /// the module's existing BCS condensation-energy formula
    /// ([`condensation_energy`]).
    ///
    /// This is a physical diagnostic tracked alongside, but independently
    /// of, the reduced-unit Landau free energy used to drive the relaxation
    /// dynamics -- see the struct-level documentation for why it is not
    /// folded directly into [`free_energy_density`](Self::free_energy_density).
    ///
    /// # Errors
    /// Propagates errors from [`condensation_energy`] (non-physical inputs).
    pub fn condensation_energy_at(&self, amplitude: f64) -> Result<f64> {
        condensation_energy(self.density_of_states, amplitude)
    }

    /// Right-hand side of the TDGL equation, dΔ/dt = f(Δ) = -Γ·∂F/∂Δ.
    fn relaxation_rhs(&self, amplitude: f64, temperature: f64) -> f64 {
        -self.relaxation_rate * self.free_energy_gradient(amplitude, temperature)
    }

    /// Integrate the TDGL relaxation dΔ/dt = -Γ·∂F/∂Δ starting from
    /// `initial_amplitude` at fixed `temperature`, using classical 4th-order
    /// Runge-Kutta (RK4).
    ///
    /// Returns the amplitude trajectory \[eV\] with
    /// `trajectory[0] == initial_amplitude` and `trajectory[n_steps]` the
    /// amplitude after `n_steps` steps of size `dt`. For stable `dt`, the
    /// trajectory converges to `self.equilibrium_amplitude(temperature)` as
    /// `n_steps` grows -- this is the self-consistent gap already computed
    /// by [`SdwGapSolver::gap_at_temperature`].
    ///
    /// # Stability
    /// This is an explicit integrator: `dt` must be small compared to the
    /// natural relaxation timescale near equilibrium,
    /// τ ≈ 1 / (2 Γ b Δ_eq(T)²), or the iteration becomes numerically
    /// unstable (RK4's linear stability region is larger than explicit
    /// Euler's, but still finite). Note that τ *diverges* as T → T_N
    /// (Δ_eq → 0): this is the physically correct "critical slowing down"
    /// of order-parameter relaxation near a continuous phase transition, and
    /// callers relaxing close to T_N should budget more steps accordingly.
    ///
    /// # Arguments
    /// * `initial_amplitude` - starting amplitude Δ(0) \[eV\], must be
    ///   non-negative
    /// * `temperature` - fixed bath temperature \[K\] for the relaxation
    /// * `dt` - integration time step \[s\], must be positive
    /// * `n_steps` - number of RK4 steps to take
    ///
    /// # Errors
    /// Returns an error if `initial_amplitude` is negative or `dt` is not
    /// positive.
    pub fn relax_to_equilibrium(
        &self,
        initial_amplitude: f64,
        temperature: f64,
        dt: f64,
        n_steps: usize,
    ) -> Result<Vec<f64>> {
        if initial_amplitude < 0.0 {
            return Err(error::invalid_param(
                "initial_amplitude",
                "must be non-negative",
            ));
        }
        if dt <= 0.0 {
            return Err(error::invalid_param("dt", "must be positive"));
        }

        let mut trajectory = Vec::with_capacity(n_steps + 1);
        trajectory.push(initial_amplitude);

        let mut amplitude = initial_amplitude;
        for _ in 0..n_steps {
            let k1 = self.relaxation_rhs(amplitude, temperature);
            let k2 = self.relaxation_rhs(amplitude + 0.5 * dt * k1, temperature);
            let k3 = self.relaxation_rhs(amplitude + 0.5 * dt * k2, temperature);
            let k4 = self.relaxation_rhs(amplitude + dt * k3, temperature);

            amplitude += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            trajectory.push(amplitude);
        }

        Ok(trajectory)
    }
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

    // ========================================================================
    // TDGL relaxational dynamics tests
    // ========================================================================

    /// Standard test-fixture `SdwRelaxationDynamics`, reusing the same
    /// gap-solver parameters (coupling = 0.5, cutoff = 1 eV, T_N = 311 K) as
    /// the static-SDW tests above.
    fn test_relaxation_dynamics() -> SdwRelaxationDynamics {
        let gap_solver =
            SdwGapSolver::new(0.5, 1.0, 311.0).expect("SdwGapSolver::new should succeed");
        SdwRelaxationDynamics::new(gap_solver, 1e28, 1.0, 1.0)
            .expect("SdwRelaxationDynamics::new should succeed")
    }

    #[test]
    fn test_sdw_relaxation_converges_to_self_consistent_gap() {
        let dynamics = test_relaxation_dynamics();

        let temperature = 150.0;
        let target = dynamics.equilibrium_amplitude(temperature);
        assert!(
            target > 0.0,
            "equilibrium amplitude should be positive well below T_N"
        );

        // Start well displaced below equilibrium (a small seed amplitude)
        // and integrate for many decay times.
        let initial = 0.05 * target;
        let trajectory = dynamics
            .relax_to_equilibrium(initial, temperature, 0.05, 50_000)
            .expect("relax_to_equilibrium should succeed");

        let relaxed = *trajectory
            .last()
            .expect("trajectory must contain at least the initial amplitude");
        let rel_err = (relaxed - target).abs() / target;
        assert!(
            rel_err < 1e-6,
            "relaxed amplitude {relaxed} should converge to the self-consistent gap {target} \
             (relative error {rel_err})"
        );

        // Sanity: relaxing down from above equilibrium converges to the same point.
        let trajectory_above = dynamics
            .relax_to_equilibrium(1.8 * target, temperature, 0.05, 50_000)
            .expect("relax_to_equilibrium should succeed");
        let relaxed_above = *trajectory_above
            .last()
            .expect("trajectory must contain at least the initial amplitude");
        let rel_err_above = (relaxed_above - target).abs() / target;
        assert!(
            rel_err_above < 1e-6,
            "relaxing from above equilibrium should also converge to {target}, got {relaxed_above}"
        );
    }

    #[test]
    fn test_sdw_relaxation_free_energy_monotonically_decreases() {
        let dynamics = test_relaxation_dynamics();
        let temperature = 200.0;
        let target = dynamics.equilibrium_amplitude(temperature);
        assert!(
            target > 0.0,
            "equilibrium amplitude should be positive at T=200K"
        );

        // Pure relaxation (no driving), started both below and above the
        // equilibrium amplitude.
        for &initial in &[0.1 * target, 1.7 * target] {
            let trajectory = dynamics
                .relax_to_equilibrium(initial, temperature, 0.01, 2_000)
                .expect("relax_to_equilibrium should succeed");

            let free_energies: Vec<f64> = trajectory
                .iter()
                .map(|&amplitude| {
                    dynamics
                        .free_energy_density(amplitude, temperature)
                        .expect("free_energy_density should succeed")
                })
                .collect();

            for window in free_energies.windows(2) {
                let (prev, next) = (window[0], window[1]);
                assert!(
                    next <= prev + 1e-15,
                    "pure relaxation must not increase the free energy: {prev} -> {next} \
                     (initial amplitude {initial})"
                );
            }

            // The free energy should also have genuinely decreased overall
            // (not merely stayed flat), confirming the dynamics actually moved.
            let first = *free_energies.first().expect("non-empty trajectory");
            let last = *free_energies.last().expect("non-empty trajectory");
            assert!(
                last < first,
                "free energy should strictly decrease from {first} to {last} over the relaxation"
            );
        }
    }

    #[test]
    fn test_sdw_relaxation_amplitude_vanishes_near_neel_temperature() {
        let dynamics = test_relaxation_dynamics();

        // The self-consistent target itself should continuously vanish
        // approaching T_N (already exercised for `gap_at_temperature`
        // directly in `test_sdw_gap_temperature_dependence`; re-verify the
        // trend holds through the new `equilibrium_amplitude` wrapper too).
        let t_well_below = dynamics.equilibrium_amplitude(100.0);
        let t_mid = dynamics.equilibrium_amplitude(250.0);
        let t_just_below = dynamics.equilibrium_amplitude(305.0);
        let t_closer = dynamics.equilibrium_amplitude(309.0);

        assert!(
            t_well_below > t_mid,
            "target should shrink monotonically approaching T_N"
        );
        assert!(
            t_mid > t_just_below,
            "target should shrink monotonically approaching T_N"
        );
        assert!(
            t_just_below > t_closer,
            "target should shrink monotonically approaching T_N"
        );
        assert!(
            t_closer < 0.2 * t_well_below,
            "target amplitude should be much smaller just below T_N than well below T_N"
        );

        // Now actually relax the dynamics (not just read off the target) at
        // a temperature well below T_N and one just below T_N, confirming
        // the *relaxed* amplitude reproduces the same vanishing trend.
        // Critical slowing down (tau ~ 1/Delta_eq^2) means the near-T_N case
        // needs many more steps to reach comparable relative convergence.
        let seed = 0.05 * t_well_below;

        let relaxed_well_below = *dynamics
            .relax_to_equilibrium(seed, 100.0, 0.05, 50_000)
            .expect("relax_to_equilibrium should succeed")
            .last()
            .expect("non-empty trajectory");

        let relaxed_just_below = *dynamics
            .relax_to_equilibrium(seed, 305.0, 0.05, 2_000_000)
            .expect("relax_to_equilibrium should succeed")
            .last()
            .expect("non-empty trajectory");

        assert!(
            relaxed_just_below < relaxed_well_below,
            "relaxed amplitude near T_N ({relaxed_just_below}) should be much smaller than \
             well below T_N ({relaxed_well_below})"
        );
        assert!(
            (relaxed_just_below - t_just_below).abs() / t_just_below < 1e-2,
            "relaxed amplitude near T_N should converge close to its (small) self-consistent \
             target: relaxed {relaxed_just_below}, target {t_just_below}"
        );

        // At and above T_N the target -- and hence the amplitude that
        // relaxation converges to -- is exactly zero.
        let at_tn = dynamics.equilibrium_amplitude(311.0);
        let above_tn = dynamics.equilibrium_amplitude(400.0);
        assert!(at_tn.abs() < 1e-15, "target must vanish exactly at T_N");
        assert!(above_tn.abs() < 1e-15, "target must vanish above T_N");
    }

    #[test]
    fn test_sdw_relaxation_dynamics_invalid_parameters() {
        let gap_solver =
            SdwGapSolver::new(0.5, 1.0, 311.0).expect("SdwGapSolver::new should succeed");

        // Non-positive density of states, quartic stiffness, relaxation rate.
        assert!(SdwRelaxationDynamics::new(gap_solver.clone(), 0.0, 1.0, 1.0).is_err());
        assert!(SdwRelaxationDynamics::new(gap_solver.clone(), 1e28, 0.0, 1.0).is_err());
        assert!(SdwRelaxationDynamics::new(gap_solver.clone(), 1e28, 1.0, 0.0).is_err());

        let dynamics = SdwRelaxationDynamics::new(gap_solver, 1e28, 1.0, 1.0)
            .expect("SdwRelaxationDynamics::new should succeed");

        // Negative initial amplitude and non-positive time step must be rejected.
        assert!(dynamics
            .relax_to_equilibrium(-1.0, 150.0, 0.05, 10)
            .is_err());
        assert!(dynamics.relax_to_equilibrium(0.1, 150.0, 0.0, 10).is_err());
        assert!(dynamics
            .relax_to_equilibrium(0.1, 150.0, -0.05, 10)
            .is_err());
    }
}
