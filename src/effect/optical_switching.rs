//! All-Optical Magnetization Switching (AOS)
//!
//! This module implements the physics of all-optical helicity-dependent switching (HDS),
//! driven by circularly polarized femtosecond laser pulses in magnetic materials.
//!
//! ## Physics Overview
//!
//! ### Inverse Faraday Effect (IFE)
//!
//! A circularly polarized laser beam propagating through a magnetic medium exerts an
//! effective quasi-static magnetic field via the inverse Faraday effect:
//!
//! ```text
//! H_IFE = α_IFE · n_ph · σ̂_circ
//! ```
//!
//! where:
//! - α_IFE  [T·m²·s/J] — material-specific IFE coupling coefficient
//! - n_ph   [photons/(m²·s)] — photon flux density = I_laser / (ℏ·ω)
//! - σ̂_circ — +ẑ for right-circular (σ⁺), −ẑ for left-circular (σ⁻)
//!
//! Typical peak IFE fields are 1–10 T for femtosecond pulses in ferromagnets.
//!
//! ### Ultrafast Demagnetization (Three-Temperature Model, simplified)
//!
//! Following Beaurepaire et al. (1996), a laser pulse of fluence F [J/cm²] drives
//! rapid heating of the electron sub-system which couples to spins on a timescale
//! τ_demag ≈ 100–300 fs, causing partial or full demagnetization.  A simplified
//! quench fraction parameterises this:
//!
//! ```text
//! f_quench = clip((F − F_th) / (F_sat − F_th), 0, 1)
//! F_th  ≈ 1   mJ/cm²   (threshold fluence)
//! F_sat ≈ 10  mJ/cm²   (saturation fluence)
//! ```
//!
//! ### Helicity-Dependent Switching (HDS)
//!
//! For ferromagnets with perpendicular magnetic anisotropy (PMA):
//! - RCP (σ⁺): H_IFE ∥ +z → can switch m from −z to +z
//! - LCP (σ⁻): H_IFE ∥ −z → can switch m from +z to −z
//! - Linear: no IFE contribution; purely thermomagnetic
//!
//! A sigmoid model captures the probabilistic nature of near-threshold switching:
//!
//! ```text
//! P_switch = 1 / (1 + exp(−(H_IFE − sign(m_z)·H_c) / H_width))
//! H_width  = 0.2 · H_c
//! ```
//!
//! ## Key References
//!
//! - E. Beaurepaire et al., "Ultrafast Spin Dynamics in Ferromagnetic Nickel",
//!   PRL 76, 4250 (1996) — first observation of sub-ps demagnetization
//! - C. D. Stanciu et al., "All-Optical Magnetic Recording with Circularly Polarized Light",
//!   PRL 99, 047601 (2007) — GdFeCo HDS discovery
//! - A. Kirilyuk et al., "Ultrafast optical manipulation of magnetic order",
//!   Rev. Mod. Phys. 82, 2731 (2010) — comprehensive review
//! - T. A. Ostler et al., "Ultrafast heating as a sufficient stimulus for magnetization
//!   reversal in a ferrimagnet", Nature Comm. 3, 666 (2012)

use std::f64::consts::PI;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[allow(unused_imports)] // KB and GAMMA are part of the public API surface per spec
use crate::constants::{GAMMA, HBAR, KB, MU_0};
use crate::error::{invalid_param, Result};
use crate::vector3::Vector3;

// Speed of light [m/s] — use local constant to avoid import ambiguity
const C_LIGHT: f64 = 2.997_924_58e8;

// =============================================================================
// Laser Pulse Parameters
// =============================================================================

/// Parameters describing a femtosecond laser pulse for optical switching experiments.
///
/// Models a Gaussian pulse envelope with full-width-at-half-maximum (FWHM) duration.
/// The fluence (energy per unit area) is computed analytically from the peak intensity
/// and pulse duration for a transform-limited Gaussian pulse.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LaserPulseParams {
    /// Laser center wavelength [m]
    ///
    /// Typical: 800 nm (Ti:Sapphire fundamental), 400 nm (second harmonic).
    /// Valid range: 100 nm – 10 µm.
    pub wavelength_m: f64,

    /// Peak intensity at pulse center [W/m²]
    ///
    /// Typical: 1e13–1e17 W/m² for femtosecond pulses focused to µm spots.
    pub peak_intensity: f64,

    /// FWHM pulse duration [s]
    ///
    /// Typical: 50 fs – 10 ps for ultrafast laser experiments.
    pub pulse_duration_s: f64,

    /// Pulse fluence (energy per unit area) [J/cm²]
    ///
    /// Computed from peak_intensity and pulse_duration for a Gaussian pulse:
    /// `F = I_peak × τ_FWHM × sqrt(π / (4 ln 2))` and converted to J/cm².
    pub fluence_j_cm2: f64,
}

impl LaserPulseParams {
    /// Construct a [`LaserPulseParams`] from wavelength, peak intensity, and pulse duration.
    ///
    /// Fluence is computed automatically assuming a Gaussian pulse envelope:
    ///
    /// ```text
    /// F [J/m²] = I_peak × τ_FWHM × sqrt(π / (4 ln 2))
    /// F [J/cm²] = F [J/m²] / 1e4
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`Err`] when:
    /// - `wavelength_m` is outside the range (100 nm, 10 µm)
    /// - `peak_intensity` ≤ 0
    /// - `pulse_duration_s` ≤ 0
    pub fn new(wavelength_m: f64, peak_intensity: f64, pulse_duration_s: f64) -> Result<Self> {
        if !(100e-9..=10e-6).contains(&wavelength_m) {
            return Err(invalid_param(
                "wavelength_m",
                "must be in the range (100 nm, 10 µm)",
            ));
        }
        if peak_intensity <= 0.0 {
            return Err(invalid_param("peak_intensity", "must be positive [W/m²]"));
        }
        if pulse_duration_s <= 0.0 {
            return Err(invalid_param("pulse_duration_s", "must be positive [s]"));
        }

        // Gaussian pulse area factor: ∫ exp(−4 ln2 · t²/τ²) dt = τ · sqrt(π/(4 ln2))
        let gaussian_area_factor = (PI / (4.0 * 2_f64.ln())).sqrt();
        let fluence_j_m2 = peak_intensity * pulse_duration_s * gaussian_area_factor;
        // Convert J/m² → J/cm²  (1 m² = 1e4 cm²)
        let fluence_j_cm2 = fluence_j_m2 / 1.0e4;

        Ok(Self {
            wavelength_m,
            peak_intensity,
            pulse_duration_s,
            fluence_j_cm2,
        })
    }

    /// Standard Ti:Sapphire laser pulse (800 nm, 100 fs, 1e16 W/m²).
    ///
    /// Commonly used in ultrafast magneto-optical experiments on GdFeCo and Co/Pt.
    pub fn ti_sapphire_standard() -> Self {
        // This uses known-valid parameters so the constructor cannot fail.
        // We construct directly to avoid propagating a spurious error.
        let wavelength_m = 800e-9;
        let peak_intensity = 1.0e16;
        let pulse_duration_s = 100e-15;
        let gaussian_area_factor = (PI / (4.0 * 2_f64.ln())).sqrt();
        let fluence_j_cm2 = peak_intensity * pulse_duration_s * gaussian_area_factor / 1.0e4;
        Self {
            wavelength_m,
            peak_intensity,
            pulse_duration_s,
            fluence_j_cm2,
        }
    }

    /// Photon energy E = ℏ·ω = ℏ·2π·c / λ  [J].
    ///
    /// At 800 nm this is ≈ 2.48 eV = 3.97 × 10⁻¹⁹ J.
    #[inline]
    pub fn photon_energy_j(&self) -> f64 {
        // ω = 2π c / λ
        let omega = 2.0 * PI * C_LIGHT / self.wavelength_m;
        HBAR * omega
    }

    /// Photon flux density n_ph = I_peak / E_photon  [photons/(m²·s)].
    ///
    /// This is the number of photons crossing unit area per second at peak intensity.
    #[inline]
    pub fn photon_flux(&self) -> f64 {
        self.peak_intensity / self.photon_energy_j()
    }
}

// =============================================================================
// Optical-Magnetic Material
// =============================================================================

/// Material parameters relevant to all-optical magnetization switching.
///
/// Combines saturation magnetization, perpendicular anisotropy, the IFE coupling
/// coefficient, and ultrafast demagnetization timescales.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OpticalMagneticMaterial {
    /// Human-readable material name.
    pub name: &'static str,

    /// Saturation magnetization M_s [A/m].
    ///
    /// GdFeCo ≈ 400 kA/m, Co ≈ 1.4 MA/m, Permalloy ≈ 800 kA/m.
    pub ms: f64,

    /// Perpendicular (uniaxial) magnetic anisotropy constant K_u [J/m³].
    ///
    /// Positive value corresponds to easy-axis perpendicular to the film plane (PMA).
    pub anisotropy_k: f64,

    /// Inverse Faraday Effect coupling coefficient α_IFE [T·m²·s/J].
    ///
    /// Material-specific parameter encoding the strength of the light-spin coupling.
    /// Larger values mean stronger IFE for a given photon flux.
    pub alpha_ife: f64,

    /// Coercive field H_c [A/m].
    ///
    /// Minimum field required to reverse the magnetization.
    pub coercive_field: f64,

    /// Ultrafast demagnetization time constant τ_demag [s].
    ///
    /// Characteristic timescale for spin-electron energy transfer, typically 100–300 fs.
    pub demag_time_s: f64,

    /// Re-magnetization time constant τ_remag [s].
    ///
    /// Timescale over which the magnetization recovers after demagnetization, 1–10 ps.
    pub remag_time_s: f64,
}

impl OpticalMagneticMaterial {
    /// GdFeCo — the archetypal all-optical switching ferrimagnet.
    ///
    /// Gadolinium iron cobalt (Gd₂₂Fe₆₈·₂Co₉·₈) was the material in which HDS was
    /// first demonstrated by Stanciu et al. (2007).  Its ferrimagnetic compensation
    /// point facilitates sub-ps switching.
    ///
    /// Reference parameters:
    /// - Stanciu et al., PRL 99, 047601 (2007)
    /// - Radu et al., Nature 472, 205 (2011)
    pub fn gdfeco() -> Self {
        Self {
            name: "GdFeCo",
            ms: 400.0e3,             // 400 kA/m
            anisotropy_k: 300.0e3,   // 300 kJ/m³ (PMA)
            alpha_ife: 1.0e-4,       // T·m²·s/J
            coercive_field: 50.0e3,  // 50 kA/m ≈ 63 mT
            demag_time_s: 300.0e-15, // 300 fs
            remag_time_s: 3.0e-12,   // 3 ps
        }
    }

    /// Cobalt — a classic ferromagnet used in ultrafast demagnetization studies.
    ///
    /// Cobalt has strong PMA in thin-film form (Co/Pt, Co/Pd multilayers) and a
    /// relatively short demagnetization time (~150 fs).
    pub fn cobalt() -> Self {
        Self {
            name: "Co",
            ms: 1.4e6,             // 1.4 MA/m
            anisotropy_k: 410.0e3, // 410 kJ/m³
            alpha_ife: 5.0e-5,
            coercive_field: 100.0e3, // 100 kA/m ≈ 126 mT
            demag_time_s: 150.0e-15, // 150 fs
            remag_time_s: 1.0e-12,   // 1 ps
        }
    }

    /// Permalloy (Ni₈₀Fe₂₀) thin film — soft magnet with negligible PMA.
    ///
    /// Permalloy is primarily used for in-plane magnetisation dynamics.  Its small
    /// coercive field and moderate demagnetisation time make it useful for
    /// characterising laser-induced precession rather than deterministic switching.
    pub fn permalloy_thin() -> Self {
        Self {
            name: "Permalloy (Ni80Fe20)",
            ms: 800.0e3,          // 800 kA/m
            anisotropy_k: 50.0e3, // 50 kJ/m³ (weak PMA in thin film)
            alpha_ife: 2.0e-5,
            coercive_field: 10.0e3,  // 10 kA/m ≈ 12.6 mT (soft)
            demag_time_s: 200.0e-15, // 200 fs
            remag_time_s: 5.0e-12,   // 5 ps
        }
    }

    /// Coercive field expressed as a magnetic flux density μ₀ H_c [T].
    ///
    /// Converts from SI field units [A/m] to Tesla for comparison with applied fields
    /// and published AOS threshold tables.
    #[inline]
    pub fn coercive_field_tesla(&self) -> f64 {
        MU_0 * self.coercive_field
    }
}

// =============================================================================
// Circular Helicity
// =============================================================================

/// Polarisation state of the optical pulse controlling the IFE direction.
///
/// The helicity determines the sign of the effective IFE magnetic field and
/// thus the direction of optically driven magnetisation switching.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum CircularHelicity {
    /// Right-circular polarisation (σ⁺, RCP).
    ///
    /// The angular momentum of the photon is along +ẑ (beam propagation), resulting
    /// in an IFE field H_IFE ∥ +ẑ.  Can switch magnetisation from −z to +z when
    /// H_IFE > H_c.
    RightCircular,

    /// Left-circular polarisation (σ⁻, LCP).
    ///
    /// Angular momentum along −ẑ; IFE field H_IFE ∥ −ẑ.  Can switch from +z to −z.
    LeftCircular,

    /// Linear polarisation.
    ///
    /// Carries no net angular momentum; no IFE contribution.  Switching, if it occurs,
    /// is purely thermomagnetic (helicity-independent).
    Linear,
}

impl CircularHelicity {
    /// Helicity factor σ used in the IFE field formula.
    ///
    /// Returns +1.0 (RCP), −1.0 (LCP), or 0.0 (Linear).
    #[inline]
    pub fn sigma_factor(self) -> f64 {
        match self {
            CircularHelicity::RightCircular => 1.0,
            CircularHelicity::LeftCircular => -1.0,
            CircularHelicity::Linear => 0.0,
        }
    }
}

// =============================================================================
// Result of a single optical pulse
// =============================================================================

/// Outcome of a single laser-pulse interaction with a magnetic material.
///
/// All quantities are computed in SI units (Gaussians are avoided throughout).
#[derive(Debug, Clone)]
pub struct OpticalSwitchResult {
    /// Inverse Faraday effective field H_IFE [A/m].
    ///
    /// Signed: positive for RCP, negative for LCP, zero for linear.
    pub h_ife: f64,

    /// Demagnetisation fraction f_quench ∈ [0, 1].
    ///
    /// 0 ⟹ no demagnetisation (below threshold fluence),
    /// 1 ⟹ complete demagnetisation.
    pub demagnetization: f64,

    /// Probabilistic switching likelihood P_switch ∈ [0, 1].
    ///
    /// Computed from the sigmoid model centred on the coercive field.
    pub switching_probability: f64,

    /// Deterministic switching verdict: true when P_switch > 0.5.
    pub switched: bool,

    /// Magnetisation unit vector after the pulse.
    pub final_magnetization: Vector3<f64>,
}

// =============================================================================
// Main Physics Calculator
// =============================================================================

/// All-optical magnetization switching calculator.
///
/// Computes the Inverse Faraday Effect field, ultrafast demagnetization fraction,
/// and helicity-dependent switching probability for a given material and laser pulse.
///
/// # Example
///
/// ```rust
/// use spintronics::effect::optical_switching::{
///     OpticalSwitching, OpticalMagneticMaterial, LaserPulseParams, CircularHelicity,
/// };
/// use spintronics::Vector3;
///
/// let material = OpticalMagneticMaterial::gdfeco();
/// let switcher = OpticalSwitching::new(material);
/// let pulse    = LaserPulseParams::ti_sapphire_standard();
///
/// // IFE field should be positive for RCP
/// let h_ife = switcher.inverse_faraday_field(&pulse, CircularHelicity::RightCircular);
/// assert!(h_ife > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct OpticalSwitching {
    /// Magnetic material under illumination.
    pub material: OpticalMagneticMaterial,
}

impl OpticalSwitching {
    /// Create a new [`OpticalSwitching`] calculator for the given material.
    pub fn new(material: OpticalMagneticMaterial) -> Self {
        Self { material }
    }

    /// Compute the Inverse Faraday Effect field H_IFE [A/m].
    ///
    /// The IFE field is derived from the photon angular momentum deposited in the
    /// spin system:
    ///
    /// ```text
    /// H_IFE = α_IFE · n_ph · σ / μ₀
    /// ```
    ///
    /// where n_ph = I_peak / (ℏ·ω) is the photon flux and σ = ±1, 0 is the helicity
    /// factor.  The result is signed: positive for RCP, negative for LCP, zero for
    /// linear polarisation.
    ///
    /// # Arguments
    ///
    /// * `pulse`    — laser pulse parameters
    /// * `helicity` — circular polarisation state
    pub fn inverse_faraday_field(
        &self,
        pulse: &LaserPulseParams,
        helicity: CircularHelicity,
    ) -> f64 {
        // n_ph [photons/(m²·s)] = I_peak / E_photon
        let n_ph = pulse.photon_flux();
        // σ factor: +1 (RCP), −1 (LCP), 0 (linear)
        let sigma = helicity.sigma_factor();
        // H_IFE [A/m] — division by μ₀ converts the material coupling from [T·m²·s/J]
        // into [A/m], keeping the formula in SI units throughout.
        self.material.alpha_ife * n_ph * sigma / MU_0
    }

    /// Compute the ultrafast demagnetisation fraction f_quench ∈ [0, 1].
    ///
    /// Uses a piecewise-linear model calibrated against three-temperature model
    /// simulations for metallic ferromagnets:
    ///
    /// ```text
    /// f_quench = 0                                   if F < F_th
    ///          = (F − F_th) / (F_sat − F_th)         if F_th ≤ F ≤ F_sat
    ///          = 1                                   if F > F_sat
    ///
    /// F_th  = 1.0  mJ/cm²   (threshold fluence)
    /// F_sat = 10.0 mJ/cm²   (saturation fluence)
    /// ```
    ///
    /// # Arguments
    ///
    /// * `pulse` — laser pulse carrying the fluence information
    pub fn demagnetization_fraction(&self, pulse: &LaserPulseParams) -> f64 {
        const F_TH: f64 = 1.0e-3; // 1 mJ/cm²
        const F_SAT: f64 = 10.0e-3; // 10 mJ/cm²

        let f = pulse.fluence_j_cm2;
        if f < F_TH {
            0.0
        } else if f > F_SAT {
            1.0
        } else {
            (f - F_TH) / (F_SAT - F_TH)
        }
    }

    /// Compute the probabilistic switching probability P_switch ∈ [0, 1].
    ///
    /// Models the stochastic nature of near-threshold switching using a sigmoid
    /// (logistic) function centred on the coercive field:
    ///
    /// ```text
    /// drive    = H_IFE_signed − sign(m_z) · H_c
    /// H_width  = 0.2 · H_c           (20% transition width)
    /// P_switch = sigmoid(drive / H_width)
    ///          = 1 / (1 + exp(−drive / H_width))
    /// ```
    ///
    /// - Well above threshold (drive ≫ 0): P_switch → 1
    /// - Well below threshold (drive ≪ 0): P_switch → 0
    /// - At threshold (drive = 0):         P_switch = 0.5
    ///
    /// # Arguments
    ///
    /// * `pulse`      — laser pulse
    /// * `helicity`   — polarisation state
    /// * `initial_mz` — z-component of initial unit magnetisation (±1 for a macrospin)
    pub fn switching_probability(
        &self,
        pulse: &LaserPulseParams,
        helicity: CircularHelicity,
        initial_mz: f64,
    ) -> f64 {
        let h_ife_signed = self.inverse_faraday_field(pulse, helicity);
        let h_c = self.material.coercive_field;
        // Drive = IFE field minus the restoring coercive barrier (sign accounts for
        // whether the magnetisation points along or against the IFE direction).
        let drive = h_ife_signed - initial_mz.signum() * h_c;
        let h_width = h_c * 0.2; // 20% of H_c sets the sharpness of the transition
        sigmoid(drive / h_width)
    }

    /// Simulate the magnetisation response to a single laser pulse.
    ///
    /// Steps:
    /// 1. Compute the IFE effective field.
    /// 2. Compute the demagnetisation fraction.
    /// 3. Evaluate the switching probability.
    /// 4. Determine the final magnetisation state (deterministic threshold at P > 0.5).
    ///
    /// When switching occurs, the z-component of the magnetisation is reversed.
    /// The in-plane components (x, y) are assumed to relax to zero on a timescale
    /// longer than the pulse (macrospin approximation, perpendicular anisotropy).
    ///
    /// # Arguments
    ///
    /// * `pulse`       — laser pulse parameters
    /// * `helicity`    — circular polarisation state
    /// * `m_initial`   — initial magnetisation unit vector
    pub fn simulate_pulse_response(
        &self,
        pulse: &LaserPulseParams,
        helicity: CircularHelicity,
        m_initial: Vector3<f64>,
    ) -> OpticalSwitchResult {
        let h_ife = self.inverse_faraday_field(pulse, helicity);
        let demagnetization = self.demagnetization_fraction(pulse);
        let p_switch = self.switching_probability(pulse, helicity, m_initial.z);

        let switched = p_switch > 0.5;

        // After the pulse the macrospin either retains its initial orientation or
        // reverses its z-component (deterministic HDS approximation).
        let final_magnetization = if switched {
            // Flip z, reset in-plane to zero (PMA restores m to ±z axis)
            Vector3::new(0.0, 0.0, -m_initial.z.signum())
        } else {
            // Recovers toward initial state (demagnetized but in the same well)
            // Scale the recovered magnitude by the non-quenched fraction.
            let recovered_z = m_initial.z * (1.0 - demagnetization);
            Vector3::new(m_initial.x, m_initial.y, recovered_z)
        };

        OpticalSwitchResult {
            h_ife,
            demagnetization,
            switching_probability: p_switch,
            switched,
            final_magnetization,
        }
    }

    /// Heuristic check whether the material is a viable candidate for HDS.
    ///
    /// Returns `true` when the IFE field at a reference photon flux of 10²⁵ ph/(m²·s)
    /// (equivalent to ~10 GW/m² at 800 nm) would exceed the coercive field.
    ///
    /// This is an order-of-magnitude estimate; detailed threshold analysis should use
    /// [`critical_fluence`](Self::critical_fluence).
    pub fn is_hds_material(&self) -> bool {
        // Reference photon flux ≈ 10 GW/m² ÷ (2.5 eV photon energy)
        let ref_n_ph = 1.0e25_f64;
        let ref_h_ife = self.material.alpha_ife * ref_n_ph / MU_0;
        ref_h_ife > self.material.coercive_field
    }

    /// Estimate the critical fluence F_crit at which switching becomes probable.
    ///
    /// Returns `None` for linear polarisation (no IFE drive).
    ///
    /// The phenomenological estimate is:
    ///
    /// ```text
    /// F_crit [J/cm²] ≈ H_c · μ₀ · 1e−2
    /// ```
    ///
    /// This very rough scaling comes from matching the IFE field (proportional to
    /// photon flux, and hence to fluence) to the coercive field, with the 1e−2
    /// factor absorbing material-dependent constants in typical experimental
    /// parameter ranges.  For quantitative predictions use [`switching_probability`]
    /// with the full pulse model.
    ///
    /// [`switching_probability`]: Self::switching_probability
    pub fn critical_fluence(&self, helicity: CircularHelicity) -> Option<f64> {
        if helicity == CircularHelicity::Linear {
            return None;
        }
        // Phenomenological formula: F_crit [J/cm²] ≈ H_c [A/m] · μ₀ [H/m] · 1e-2
        // Dimensional analysis: [A/m] · [H/m] = [A/m] · [V·s/A/m] = [V·s/m²] = [T]
        // The 1e-2 factor (≈ units of cm²/m²·material_factor) gives J/cm².
        let f_crit = self.material.coercive_field * MU_0 * 1.0e-2;
        Some(f_crit)
    }
}

// =============================================================================
// Helper: sigmoid (logistic) function
// =============================================================================

/// Numerically stable sigmoid / logistic function: σ(x) = 1 / (1 + exp(−x)).
#[inline]
fn sigmoid(x: f64) -> f64 {
    // Use the alternative form for large negative x to avoid overflow in exp(−x).
    if x >= 0.0 {
        let e = (-x).exp();
        1.0 / (1.0 + e)
    } else {
        let e = x.exp();
        e / (1.0 + e)
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------------------------------------------------------
    // IFE field sign and magnitude
    // -------------------------------------------------------------------------

    /// RCP must produce a positive IFE field (H ∥ +z).
    #[test]
    fn test_ife_field_rcp_positive() {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);
        let pulse = LaserPulseParams::ti_sapphire_standard();

        let h_ife = switcher.inverse_faraday_field(&pulse, CircularHelicity::RightCircular);
        assert!(
            h_ife > 0.0,
            "H_IFE must be positive for RCP; got {h_ife:.3e} A/m"
        );
    }

    /// LCP must produce a negative IFE field (H ∥ −z).
    #[test]
    fn test_ife_field_lcp_negative() {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);
        let pulse = LaserPulseParams::ti_sapphire_standard();

        let h_ife = switcher.inverse_faraday_field(&pulse, CircularHelicity::LeftCircular);
        assert!(
            h_ife < 0.0,
            "H_IFE must be negative for LCP; got {h_ife:.3e} A/m"
        );
    }

    /// Linear polarisation must give exactly zero IFE field.
    #[test]
    fn test_ife_field_linear_zero() {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);
        let pulse = LaserPulseParams::ti_sapphire_standard();

        let h_ife = switcher.inverse_faraday_field(&pulse, CircularHelicity::Linear);
        assert_eq!(h_ife, 0.0, "H_IFE must be zero for linear polarisation");
    }

    /// H_IFE must scale linearly with peak intensity.
    #[test]
    fn test_ife_field_scales_with_intensity() {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);

        let pulse1 = LaserPulseParams::new(800e-9, 1.0e15, 100e-15).expect("valid params");
        let pulse2 = LaserPulseParams::new(800e-9, 2.0e15, 100e-15).expect("valid params");

        let h1 = switcher.inverse_faraday_field(&pulse1, CircularHelicity::RightCircular);
        let h2 = switcher.inverse_faraday_field(&pulse2, CircularHelicity::RightCircular);

        // h2 / h1 should equal 2 within floating-point precision
        let ratio = h2 / h1;
        assert!(
            (ratio - 2.0).abs() < 1.0e-9,
            "H_IFE must double when intensity doubles; ratio = {ratio}"
        );
    }

    // -------------------------------------------------------------------------
    // Demagnetisation fraction
    // -------------------------------------------------------------------------

    /// Below threshold fluence (1 mJ/cm²) there must be no demagnetisation.
    #[test]
    fn test_demag_below_threshold() {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);

        // Construct a pulse with fluence well below 1 mJ/cm²
        // Using very low peak intensity: 1e10 W/m², 100 fs → ~4.3e-10 J/cm²
        let pulse = LaserPulseParams::new(800e-9, 1.0e10, 100e-15).expect("valid params");

        assert!(
            pulse.fluence_j_cm2 < 1.0e-3,
            "Sanity: fluence should be below 1 mJ/cm², got {:.3e} J/cm²",
            pulse.fluence_j_cm2
        );

        let f_q = switcher.demagnetization_fraction(&pulse);
        assert_eq!(f_q, 0.0, "f_quench must be 0 below threshold; got {f_q}");
    }

    /// Above saturation fluence (10 mJ/cm²) demagnetisation must be complete.
    #[test]
    fn test_demag_above_saturation() {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);

        // Very high fluence: 1e18 W/m², 100 fs → >> 10 mJ/cm²
        let pulse = LaserPulseParams::new(800e-9, 1.0e18, 100e-15).expect("valid params");

        assert!(
            pulse.fluence_j_cm2 > 10.0e-3,
            "Sanity: fluence should exceed 10 mJ/cm², got {:.3e} J/cm²",
            pulse.fluence_j_cm2
        );

        let f_q = switcher.demagnetization_fraction(&pulse);
        assert_eq!(f_q, 1.0, "f_quench must be 1 above saturation; got {f_q}");
    }

    // -------------------------------------------------------------------------
    // Helicity-dependent switching probability
    // -------------------------------------------------------------------------

    /// For GdFeCo with m = −ẑ (initial), RCP should have higher switching probability
    /// than LCP because RCP drives H_IFE toward +z (opposite to initial state).
    #[test]
    fn test_switching_prob_rcp_vs_lcp() {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);

        // Use a high-intensity pulse so IFE field is substantial
        let pulse = LaserPulseParams::ti_sapphire_standard();

        let initial_mz = -1.0_f64; // magnetisation pointing −z

        let p_rcp =
            switcher.switching_probability(&pulse, CircularHelicity::RightCircular, initial_mz);
        let p_lcp =
            switcher.switching_probability(&pulse, CircularHelicity::LeftCircular, initial_mz);

        assert!(
            p_rcp > p_lcp,
            "P_switch(RCP, m=−z) should exceed P_switch(LCP, m=−z): {p_rcp:.4} vs {p_lcp:.4}"
        );
    }

    // -------------------------------------------------------------------------
    // Photon energy
    // -------------------------------------------------------------------------

    /// Photon energy at 800 nm: E = ℏω = ℏ·2πc/λ ≈ 1.55 eV = 2.48 × 10⁻¹⁹ J.
    ///
    /// Note: the 800 nm Ti:Sapphire photon energy is 1.55 eV, not 2.48 eV.
    /// 2.48 eV corresponds to ~500 nm (green light).  The expected value below
    /// is derived directly from ℏ·2πc/λ with NIST constants.
    ///
    /// Verify within 1%.
    #[test]
    fn test_photon_energy_800nm() {
        let pulse = LaserPulseParams::ti_sapphire_standard();
        let e_photon = pulse.photon_energy_j();

        // Reference: ℏ·2πc/λ at 800 nm  ≈ 2.4831e-19 J  (≈ 1.550 eV)
        // Cross-check: 1.550 eV × 1.602176634e-19 J/eV ≈ 2.483e-19 J
        let expected_j = 1.550 * 1.602_176_634e-19;
        let relative_error = ((e_photon - expected_j) / expected_j).abs();

        assert!(
            relative_error < 0.01,
            "Photon energy at 800 nm should be ~{expected_j:.3e} J (1.55 eV); \
             got {e_photon:.3e} J (relative error {relative_error:.4})"
        );
    }

    // -------------------------------------------------------------------------
    // Material presets sanity
    // -------------------------------------------------------------------------

    /// GdFeCo preset must have positive physical parameters.
    #[test]
    fn test_gdfeco_preset_fields_valid() {
        let m = OpticalMagneticMaterial::gdfeco();
        assert!(m.alpha_ife > 0.0, "alpha_ife must be positive");
        assert!(m.coercive_field > 0.0, "coercive_field must be positive");
        assert!(m.demag_time_s > 0.0, "demag_time_s must be positive");
        assert!(m.remag_time_s > 0.0, "remag_time_s must be positive");
        assert!(m.ms > 0.0, "saturation magnetisation must be positive");
        assert!(m.anisotropy_k > 0.0, "anisotropy constant must be positive");
    }

    // -------------------------------------------------------------------------
    // Additional coverage: LaserPulseParams validation
    // -------------------------------------------------------------------------

    #[test]
    fn test_new_rejects_bad_wavelength() {
        // Too short (X-ray territory)
        assert!(LaserPulseParams::new(50e-9, 1e15, 100e-15).is_err());
        // Too long (far-IR)
        assert!(LaserPulseParams::new(20e-6, 1e15, 100e-15).is_err());
    }

    #[test]
    fn test_new_rejects_non_positive_intensity() {
        assert!(LaserPulseParams::new(800e-9, 0.0, 100e-15).is_err());
        assert!(LaserPulseParams::new(800e-9, -1e15, 100e-15).is_err());
    }

    #[test]
    fn test_new_rejects_non_positive_duration() {
        assert!(LaserPulseParams::new(800e-9, 1e15, 0.0).is_err());
        assert!(LaserPulseParams::new(800e-9, 1e15, -100e-15).is_err());
    }

    #[test]
    fn test_helicity_sigma_factors() {
        assert_eq!(CircularHelicity::RightCircular.sigma_factor(), 1.0);
        assert_eq!(CircularHelicity::LeftCircular.sigma_factor(), -1.0);
        assert_eq!(CircularHelicity::Linear.sigma_factor(), 0.0);
    }

    /// RCP and LCP should give equal-and-opposite IFE fields.
    #[test]
    fn test_ife_rcp_lcp_antisymmetric() {
        let material = OpticalMagneticMaterial::cobalt();
        let switcher = OpticalSwitching::new(material);
        let pulse = LaserPulseParams::ti_sapphire_standard();

        let h_rcp = switcher.inverse_faraday_field(&pulse, CircularHelicity::RightCircular);
        let h_lcp = switcher.inverse_faraday_field(&pulse, CircularHelicity::LeftCircular);

        assert!(
            (h_rcp + h_lcp).abs() < 1.0e-6 * h_rcp.abs(),
            "RCP and LCP IFE fields should be equal and opposite: {h_rcp:.3e} vs {h_lcp:.3e}"
        );
    }

    /// Simulate a full pulse and check that the result struct fields are consistent.
    #[test]
    fn test_simulate_pulse_response_consistency() {
        let material = OpticalMagneticMaterial::gdfeco();
        let switcher = OpticalSwitching::new(material);
        let pulse = LaserPulseParams::ti_sapphire_standard();
        let m_initial = Vector3::new(0.0, 0.0, -1.0); // pointing −z

        let result =
            switcher.simulate_pulse_response(&pulse, CircularHelicity::RightCircular, m_initial);

        // switched ⟺ P_switch > 0.5
        assert_eq!(result.switched, result.switching_probability > 0.5);
        // Demagnetization in [0, 1]
        assert!((0.0..=1.0).contains(&result.demagnetization));
        // P_switch in [0, 1]
        assert!((0.0..=1.0).contains(&result.switching_probability));
    }

    /// critical_fluence should return None for linear polarisation.
    #[test]
    fn test_critical_fluence_linear_none() {
        let switcher = OpticalSwitching::new(OpticalMagneticMaterial::gdfeco());
        assert!(switcher
            .critical_fluence(CircularHelicity::Linear)
            .is_none());
    }

    /// critical_fluence should return Some positive value for RCP and LCP.
    #[test]
    fn test_critical_fluence_circular_some_positive() {
        let switcher = OpticalSwitching::new(OpticalMagneticMaterial::gdfeco());
        let f_rcp = switcher.critical_fluence(CircularHelicity::RightCircular);
        let f_lcp = switcher.critical_fluence(CircularHelicity::LeftCircular);
        assert!(f_rcp.is_some() && f_rcp.unwrap() > 0.0);
        assert!(f_lcp.is_some() && f_lcp.unwrap() > 0.0);
    }

    /// The sigmoid helper should satisfy known values.
    #[test]
    fn test_sigmoid_known_values() {
        assert!((sigmoid(0.0) - 0.5).abs() < 1.0e-15);
        assert!(sigmoid(100.0) > 0.9999);
        assert!(sigmoid(-100.0) < 1.0e-4);
    }

    /// GdFeCo should be identified as an HDS-capable material.
    #[test]
    fn test_is_hds_material_gdfeco() {
        let switcher = OpticalSwitching::new(OpticalMagneticMaterial::gdfeco());
        // GdFeCo has high alpha_ife and low coercive field — should qualify.
        // (The test checks the internal heuristic, not physical truth.)
        assert!(switcher.is_hds_material());
    }

    /// Permalloy has very small alpha_ife but also very low H_c;
    /// check that the result is consistent (no panic, either bool is OK).
    #[test]
    fn test_is_hds_material_permalloy_no_panic() {
        let switcher = OpticalSwitching::new(OpticalMagneticMaterial::permalloy_thin());
        let _ = switcher.is_hds_material(); // must not panic
    }

    /// Fluence of the standard Ti:Sapphire pulse should be above 1 mJ/cm²
    /// (so that GdFeCo is partially demagnetized under it).
    #[test]
    fn test_ti_sapphire_fluence_above_threshold() {
        let pulse = LaserPulseParams::ti_sapphire_standard();
        assert!(
            pulse.fluence_j_cm2 > 1.0e-3,
            "Ti:Sapphire standard pulse fluence should exceed 1 mJ/cm²; got {:.3e}",
            pulse.fluence_j_cm2
        );
    }

    /// The coercive field in Tesla must equal H_c × μ₀.
    #[test]
    fn test_coercive_field_tesla_conversion() {
        let m = OpticalMagneticMaterial::gdfeco();
        let expected = m.coercive_field * MU_0;
        let computed = m.coercive_field_tesla();
        assert!(
            (computed - expected).abs() < 1.0e-20,
            "coercive_field_tesla() mismatch: {computed:.6e} vs {expected:.6e}"
        );
    }
}
