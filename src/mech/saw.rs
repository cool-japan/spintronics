//! Surface Acoustic Wave (SAW) Magnetoacoustics
//!
//! This module implements the physics of Surface Acoustic Wave (SAW) driven
//! magnetization dynamics in ferromagnetic thin films grown on piezoelectric
//! substrates — a sub-field known as **magneto-acoustics** or **acoustic
//! spintronics**.
//!
//! ## Physical Background
//!
//! ### Rayleigh SAWs on Piezoelectric Substrates
//! Rayleigh surface acoustic waves are elastic waves confined within roughly
//! one wavelength of the surface of a piezoelectric crystal (LiNbO₃, ZnO,
//! GaAs, …). An inter-digital transducer (IDT) converts an RF voltage into
//! a SAW at frequency f_SAW with:
//!
//! - Phase velocity  v_SAW  [m/s] (material-dependent; ~3000–4000 m/s for LiNbO₃)
//! - Wavelength      λ_SAW = v_SAW / f_SAW
//! - Penetration     δ ≈ λ_SAW / (2π) — exponential depth decay
//! - Strain amplitude ε₀ (dimensionless; ~10⁻⁵–10⁻³ for IDT-generated SAWs)
//!
//! For propagation along x̂ the dominant strain component is:
//!   ε_xx(x,z,t) = ε₀ · sin(k·x − ω·t) · exp(−k·z)
//! with k = 2π·f / v_SAW, ω = 2π·f.  The Poisson coupling gives
//!   ε_zz ≈ −ν · ε_xx.
//!
//! ### Magnetoelastic Coupling
//! When a ferromagnetic film is grown on the piezoelectric substrate the
//! time-periodic strain acts on the magnetization through the magnetoelastic
//! energy density:
//!   E_me = B₁ (ε_xx m_x² + ε_yy m_y² + ε_zz m_z²)
//!          + B₂ (ε_xy m_x m_y + …)
//!
//! where B₁, B₂ [J/m³] are the magnetoelastic coupling constants.
//! Differentiating with respect to **m** gives an effective driving field:
//!   H_me = 2 B₁ ε_xx / (μ₀ M_s)  (simplified, along dominant ε_xx axis)
//!
//! ### SAW-Driven FMR (Acoustic FMR)
//! When the SAW frequency equals the FMR frequency the magnetization enters
//! resonant precession described by a Lorentzian amplitude:
//!   m_⊥(ω) ∝ 1 / sqrt[(1 − (ω/ω_FMR)²)² + α²]
//!
//! The resonance half-linewidth is set by the Gilbert damping α, giving a
//! quality factor Q = 1/α.
//!
//! ### Acoustic Spin Pumping
//! Precessing magnetization in an FM/NM bilayer pumps a spin current:
//!   J_s ∝ (ℏ/4π) g_r ω_SAW sin²(θ_cone)
//! This enables generation of pure spin currents without microwave cavities.
//!
//! ## Reference Values
//! | Material  | Parameter | Value       |
//! |-----------|-----------|-------------|
//! | LiNbO₃   | v_SAW     | 3985 m/s    |
//! | LiNbO₃   | K²        | 5.5 %       |
//! | GaAs      | v_SAW     | 2900 m/s    |
//! | ZnO       | v_SAW     | 2700 m/s    |
//! | Nickel    | B₁        | −8.5 MJ/m³  |
//! | Permalloy | B₁        | +6.7 MJ/m³  |
//! | Cobalt    | B₁        | −8.0 MJ/m³  |
//!
//! ## References
//! - Weiler et al., PRL **106**, 117601 (2011) — acoustic FMR
//! - Thevenard et al., PRB **93**, 134430 (2016) — SAW-driven precession
//! - Gowtham et al., PRB **94**, 014436 (2016) — acoustic spin pumping
//! - Xu et al., PRB **102**, 075432 (2020) — magneto-acoustic coupling review

use std::f64::consts::PI;

use crate::constants::{GAMMA, HBAR, KB, MU_0};
use crate::error::{invalid_param, Result};
use crate::vector3::Vector3;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// Keep KB in scope — it is exported so that downstream code can scale thermal
// fluctuations without re-importing constants, but we also need it here for
// tests.
#[allow(dead_code)]
const _KB_USED: f64 = KB;

// =============================================================================
// PiezoSubstrate
// =============================================================================

/// A piezoelectric substrate that supports Rayleigh surface acoustic waves.
///
/// Different crystal cuts and orientations give different SAW velocities and
/// electromechanical coupling strengths.  The three most common substrates in
/// magneto-acoustic experiments are LiNbO₃, GaAs, and ZnO.
///
/// # Physical Parameters
///
/// - `saw_velocity` — phase velocity of the Rayleigh SAW [m/s]
/// - `electromechanical_k2` — dimensionless coupling coefficient K²;
///   governs the transduction efficiency of IDT devices.
/// - `poisson_ratio` — ν; connects ε_zz to ε_xx through ε_zz ≈ −ν ε_xx.
/// - `density` — mass density ρ [kg/m³] (needed for acoustic impedance).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PiezoSubstrate {
    /// Human-readable material name.
    ///
    /// Skipped during serialization because `&'static str` is not `Deserialize`-
    /// compatible with the borrowed `'de` lifetime.
    #[cfg_attr(feature = "serde", serde(skip))]
    pub name: &'static str,

    /// Rayleigh SAW phase velocity v_SAW [m/s].
    ///
    /// Material-dependent; typical range 2700–4000 m/s for common substrates.
    pub saw_velocity: f64,

    /// Electromechanical coupling coefficient K² (dimensionless).
    ///
    /// Determines IDT transduction efficiency.  LiNbO₃ is preferred for
    /// high-coupling (K² ≈ 5.5 %) applications.
    pub electromechanical_k2: f64,

    /// Poisson's ratio ν (dimensionless).
    ///
    /// Relates out-of-plane strain to in-plane strain:
    ///   ε_zz ≈ −ν · ε_xx
    pub poisson_ratio: f64,

    /// Mass density ρ [kg/m³].
    ///
    /// Used for acoustic impedance Z = ρ · v_SAW.
    pub density: f64,
}

impl PiezoSubstrate {
    /// Lithium Niobate (LiNbO₃) — 128°-rotated Y-cut, X-propagation.
    ///
    /// The workhorse substrate for magneto-acoustic experiments owing to its
    /// large electromechanical coupling (K² ≈ 5.5 %) and well-characterised
    /// surface acoustic properties.
    ///
    /// Parameters:
    /// - v_SAW = 3985 m/s
    /// - K²    = 0.055
    /// - ν     = 0.25
    /// - ρ     = 4640 kg/m³
    pub fn linbo3() -> Self {
        Self {
            name: "LiNbO3",
            saw_velocity: 3985.0,
            electromechanical_k2: 0.055,
            poisson_ratio: 0.25,
            density: 4640.0,
        }
    }

    /// Gallium Arsenide (GaAs) — (001) cut.
    ///
    /// Lower coupling than LiNbO₃ but monolithic integration with III-V
    /// semiconductor quantum wells is possible.
    ///
    /// Parameters:
    /// - v_SAW = 2900 m/s
    /// - K²    = 0.007
    /// - ν     = 0.31
    /// - ρ     = 5320 kg/m³
    pub fn gaas() -> Self {
        Self {
            name: "GaAs",
            saw_velocity: 2900.0,
            electromechanical_k2: 0.007,
            poisson_ratio: 0.31,
            density: 5320.0,
        }
    }

    /// Zinc Oxide (ZnO) — c-axis perpendicular to surface.
    ///
    /// Grown as a thin-film transducer on arbitrary substrates; useful for
    /// heterogeneous integration.
    ///
    /// Parameters:
    /// - v_SAW = 2700 m/s
    /// - K²    = 0.012
    /// - ν     = 0.36
    /// - ρ     = 5600 kg/m³
    pub fn zno() -> Self {
        Self {
            name: "ZnO",
            saw_velocity: 2700.0,
            electromechanical_k2: 0.012,
            poisson_ratio: 0.36,
            density: 5600.0,
        }
    }

    /// SAW wavelength at a given excitation frequency [m].
    ///
    /// λ_SAW = v_SAW / f
    ///
    /// # Arguments
    /// * `frequency_hz` — SAW frequency f [Hz].
    #[inline]
    pub fn wavelength_m(&self, frequency_hz: f64) -> f64 {
        self.saw_velocity / frequency_hz
    }

    /// SAW wavenumber at a given excitation frequency [rad/m].
    ///
    /// k = 2π f / v_SAW = 2π / λ_SAW
    ///
    /// # Arguments
    /// * `frequency_hz` — SAW frequency f [Hz].
    #[inline]
    pub fn wavenumber(&self, frequency_hz: f64) -> f64 {
        2.0 * PI * frequency_hz / self.saw_velocity
    }

    /// Rayleigh-wave penetration depth [m].
    ///
    /// The elastic displacement decays as exp(−k·z), so the characteristic
    /// penetration depth is:
    ///   δ = 1/k = v_SAW / (2π f) = λ_SAW / (2π)
    ///
    /// # Arguments
    /// * `frequency_hz` — SAW frequency f [Hz].
    #[inline]
    pub fn penetration_depth_m(&self, frequency_hz: f64) -> f64 {
        // δ = v_SAW / (2π f)
        self.saw_velocity / (2.0 * PI * frequency_hz)
    }

    /// Acoustic impedance Z = ρ · v_SAW [Pa·s/m = N·s/m³].
    #[inline]
    pub fn acoustic_impedance(&self) -> f64 {
        self.density * self.saw_velocity
    }
}

// =============================================================================
// MagnetoelasticMaterial  (SAW-focused — distinct from mech::magnetoelastic)
// =============================================================================

/// Ferromagnetic thin-film material characterised by its magneto-acoustic
/// coupling constants.
///
/// This struct is *distinct* from [`super::magnetoelastic::MagnetoelasticMaterial`]
/// which focuses on magnetostriction and static strain.  Here we carry the full
/// set of parameters needed for dynamic SAW-driven FMR and spin pumping:
/// saturation magnetisation, Gilbert damping, uniaxial anisotropy energy, and
/// both magnetoelastic coupling constants B₁ and B₂.
///
/// # Magnetoelastic coupling constants
/// B₁ and B₂ enter the magnetoelastic energy density:
///   E_me = B₁(ε_xx m_x² + ε_yy m_y² + ε_zz m_z²)
///          + B₂(ε_xy m_x m_y + ε_yz m_y m_z + ε_xz m_x m_z)
///
/// They are related to the magnetostriction constants λ_100, λ_111 and the
/// elastic stiffness constants of the crystal.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MagnetoelasticMaterial {
    /// Human-readable material name.
    ///
    /// Skipped during serialization because `&'static str` is not `Deserialize`-
    /// compatible with the borrowed `'de` lifetime.
    #[cfg_attr(feature = "serde", serde(skip))]
    pub name: &'static str,

    /// Saturation magnetisation M_s [A/m].
    pub ms: f64,

    /// Gilbert damping parameter α (dimensionless, 0 < α ≪ 1 typically).
    ///
    /// Sets the FMR linewidth and the quality factor Q = 1/α.
    pub alpha: f64,

    /// Uniaxial anisotropy constant K [J/m³].
    ///
    /// Positive K → easy axis along the symmetry axis.
    pub anisotropy_k: f64,

    /// First magnetoelastic coupling constant B₁ [J/m³].
    ///
    /// Couples normal strains (ε_xx, ε_yy, ε_zz) to the squared direction
    /// cosines of the magnetisation.
    pub b1: f64,

    /// Second magnetoelastic coupling constant B₂ [J/m³].
    ///
    /// Couples shear strains (ε_xy, ε_yz, ε_xz) to products of direction
    /// cosines.
    pub b2: f64,
}

impl MagnetoelasticMaterial {
    /// Nickel (Ni) — face-centred cubic ferromagnet.
    ///
    /// Negative B₁ means compressive strain along the magnetisation direction
    /// lowers the energy.  Widely used in early magneto-acoustic experiments.
    ///
    /// Parameters (literature values):
    /// - M_s = 490 kA/m
    /// - α   = 0.045
    /// - K   = 5 kJ/m³  (small uniaxial contribution)
    /// - B₁  = −8.5 MJ/m³
    /// - B₂  = −9.2 MJ/m³
    pub fn nickel() -> Self {
        Self {
            name: "Nickel",
            ms: 490.0e3,
            alpha: 0.045,
            anisotropy_k: 5.0e3,
            b1: -8.5e6,
            b2: -9.2e6,
        }
    }

    /// Permalloy (Py = Ni₈₀Fe₂₀) — soft ferromagnetic alloy.
    ///
    /// Nearly zero anisotropy and low damping make Py the standard spin-pump
    /// and FMR test material.  Its positive B₁ is opposite in sign to Ni,
    /// making it complementary for acoustic phase tuning.
    ///
    /// Parameters:
    /// - M_s = 800 kA/m
    /// - α   = 0.010
    /// - K   = 200 J/m³  (weak shape/interface contribution)
    /// - B₁  = +6.7 MJ/m³
    /// - B₂  = +9.0 MJ/m³
    pub fn permalloy() -> Self {
        Self {
            name: "Permalloy",
            ms: 800.0e3,
            alpha: 0.01,
            anisotropy_k: 200.0,
            b1: 6.7e6,
            b2: 9.0e6,
        }
    }

    /// Cobalt (Co) — hexagonally close-packed ferromagnet.
    ///
    /// Large uniaxial anisotropy, moderate damping.  Used in perpendicular
    /// magnetic anisotropy devices and acoustic switching experiments.
    ///
    /// Parameters:
    /// - M_s = 1.4 MA/m
    /// - α   = 0.010
    /// - K   = 410 kJ/m³
    /// - B₁  = −8.0 MJ/m³
    /// - B₂  = −9.0 MJ/m³
    pub fn cobalt() -> Self {
        Self {
            name: "Cobalt",
            ms: 1.4e6,
            alpha: 0.01,
            anisotropy_k: 410.0e3,
            b1: -8.0e6,
            b2: -9.0e6,
        }
    }

    /// Effective anisotropy field H_k [A/m].
    ///
    /// H_k = 2 K / (μ₀ M_s)
    ///
    /// This is the field an out-of-plane (or hard-axis) bias must overcome to
    /// rotate the magnetisation away from the easy axis.
    pub fn anisotropy_field(&self) -> f64 {
        2.0 * self.anisotropy_k / (MU_0 * self.ms)
    }

    /// FMR frequency [Hz] for out-of-plane magnetisation saturation.
    ///
    /// Kittel formula (simplified — uniform mode, out-of-plane saturation):
    ///   f_FMR = γ / (2π) · μ₀ · (H_ext + H_k)
    ///
    /// Here H_ext and H_k are in A/m; μ₀ converts to Tesla before multiplying
    /// by γ [rad/(s·T)] to obtain angular frequency.
    ///
    /// For the full in-plane Kittel formula see [`SawSpinWaveExcitation`].
    ///
    /// # Arguments
    /// * `h_ext` — External applied field [A/m].
    pub fn fmr_frequency_hz(&self, h_ext: f64) -> f64 {
        GAMMA / (2.0 * PI) * MU_0 * (h_ext + self.anisotropy_field())
    }

    /// Quality factor Q = 1 / α.
    ///
    /// The Q-factor sets the FMR linewidth Δω = α ω_FMR and the peak
    /// precession amplitude at resonance.
    ///
    /// # Arguments
    /// * `_h_ext` — External applied field (not used in the simple α model;
    ///   included for API consistency with frequency-dependent damping models).
    pub fn quality_factor(&self, _h_ext: f64) -> f64 {
        1.0 / self.alpha
    }
}

// =============================================================================
// SawSource
// =============================================================================

/// A source of Rayleigh surface acoustic waves on a piezoelectric substrate.
///
/// This struct bundles the substrate properties with the IDT operating
/// parameters: drive frequency, strain amplitude, and the propagation
/// direction unit vector.
///
/// The instantaneous in-plane strain at position (x, z) and time t is:
///   ε_xx(x,z,t) = ε₀ · sin(k·x − ω·t) · exp(−k·z)
///
/// where k = 2π f / v_SAW and ω = 2π f.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SawSource {
    /// Piezoelectric substrate on which the SAW propagates.
    pub substrate: PiezoSubstrate,

    /// SAW drive frequency f [Hz].
    pub frequency_hz: f64,

    /// Peak strain amplitude ε₀ (dimensionless).
    ///
    /// Typical range: 10⁻⁵ to 10⁻³ for IDT-generated SAWs.
    /// Values ≥ 0.1 are unphysical (fracture / nonlinear regime).
    pub strain_amplitude: f64,

    /// Unit vector along the SAW propagation direction (lab frame).
    ///
    /// Defaults to x̂ = (1, 0, 0).
    pub propagation_dir: Vector3<f64>,
}

impl SawSource {
    /// Create a new `SawSource` with explicit propagation direction x̂.
    ///
    /// # Validation
    /// - `frequency_hz` must be positive.
    /// - `strain_amplitude` must be positive and < 0.1 (physical bound).
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] if any validation fails.
    pub fn new(
        substrate: PiezoSubstrate,
        frequency_hz: f64,
        strain_amplitude: f64,
    ) -> Result<Self> {
        if frequency_hz <= 0.0 {
            return Err(invalid_param(
                "frequency_hz",
                "SAW frequency must be positive",
            ));
        }
        if strain_amplitude <= 0.0 {
            return Err(invalid_param(
                "strain_amplitude",
                "strain amplitude must be positive",
            ));
        }
        if strain_amplitude >= 0.1 {
            return Err(invalid_param(
                "strain_amplitude",
                "strain amplitude ≥ 0.1 is unphysical (nonlinear / fracture regime)",
            ));
        }
        Ok(Self {
            substrate,
            frequency_hz,
            strain_amplitude,
            propagation_dir: Vector3::unit_x(),
        })
    }

    /// Pre-built: LiNbO₃ at 1 GHz with ε₀ = 10⁻⁴.
    ///
    /// Typical SAW-FMR experiment at 1 GHz — well below Ni and Co FMR but
    /// near Py FMR at low bias fields.
    pub fn linbo3_1ghz() -> Self {
        Self {
            substrate: PiezoSubstrate::linbo3(),
            frequency_hz: 1.0e9,
            strain_amplitude: 1.0e-4,
            propagation_dir: Vector3::unit_x(),
        }
    }

    /// Pre-built: LiNbO₃ at 3 GHz with ε₀ = 10⁻⁴.
    ///
    /// Closer to the zero-field FMR of Permalloy (~3 GHz); used to study
    /// resonant acoustic spin pumping.
    pub fn linbo3_3ghz() -> Self {
        Self {
            substrate: PiezoSubstrate::linbo3(),
            frequency_hz: 3.0e9,
            strain_amplitude: 1.0e-4,
            propagation_dir: Vector3::unit_x(),
        }
    }

    /// SAW wavelength [m] at the current drive frequency.
    ///
    /// λ = v_SAW / f
    #[inline]
    pub fn wavelength_m(&self) -> f64 {
        self.substrate.wavelength_m(self.frequency_hz)
    }

    /// Angular frequency ω = 2π f [rad/s].
    #[inline]
    pub fn angular_frequency(&self) -> f64 {
        2.0 * PI * self.frequency_hz
    }

    /// Instantaneous ε_xx strain at position (x, z) and time t.
    ///
    /// ε_xx(x,z,t) = ε₀ · sin(k·x − ω·t) · exp(−k·z)
    ///
    /// - z = 0 → surface; z < 0 is below surface (exponential decay).
    /// - z > 0 is above the surface → strain is super-exponential (unphysical
    ///   for Rayleigh wave; result is still returned but callers should restrict
    ///   z ≥ 0).
    ///
    /// # Arguments
    /// * `x` — Position along the propagation direction [m].
    /// * `z` — Depth below the surface [m]; 0 = surface, increasing into bulk.
    /// * `t` — Time [s].
    pub fn strain_at_point(&self, x: f64, z: f64, t: f64) -> f64 {
        let k = self.substrate.wavenumber(self.frequency_hz);
        let omega = self.angular_frequency();
        // Exponential depth decay characteristic of Rayleigh waves
        self.strain_amplitude * (k * x - omega * t).sin() * (-k * z).exp()
    }
}

// =============================================================================
// SawMagnetoacoustics
// =============================================================================

/// Main calculator for SAW–ferromagnet coupling.
///
/// Combines a [`SawSource`] and a [`MagnetoelasticMaterial`] to evaluate:
/// - The magnetoelastic driving field amplitude H_me
/// - The FMR frequency and detuning from the SAW frequency
/// - The Lorentzian precession response function
/// - The cone angle of SAW-driven magnetisation precession
/// - The acoustic spin current density (acoustic spin pumping)
/// - The resonant SAW power absorption per unit area
///
/// # Example
///
/// ```rust
/// use spintronics::mech::saw::{SawMagnetoacoustics, SawSource, MagnetoelasticMaterial};
///
/// let saw = SawSource::linbo3_3ghz();
/// let mat = MagnetoelasticMaterial::permalloy();
/// let h_res = {
///     let tmp = SawMagnetoacoustics::new(saw.clone(), mat.clone(), 20e-9, 0.0)
///         .expect("valid parameters");
///     tmp.acoustic_fmr_condition_h_ext()
/// };
/// let device = SawMagnetoacoustics::new(saw, mat, 20e-9, h_res)
///     .expect("valid parameters");
/// let theta = device.precession_cone_angle_rad();
/// assert!(theta > 0.0);
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SawMagnetoacoustics {
    /// SAW source (substrate + drive parameters).
    pub saw: SawSource,

    /// Ferromagnetic thin film material.
    pub material: MagnetoelasticMaterial,

    /// Ferromagnetic film thickness d [m].
    pub film_thickness_m: f64,

    /// External bias field H_ext [A/m].
    ///
    /// Used to tune the FMR frequency to match the SAW frequency.
    pub h_ext: f64,
}

impl SawMagnetoacoustics {
    /// Construct a new `SawMagnetoacoustics` coupling calculator.
    ///
    /// # Validation
    /// - `film_thickness_m` must be strictly positive.
    /// - `h_ext` must be non-negative (demagnetising-field effects ignored).
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] on validation failure.
    pub fn new(
        saw: SawSource,
        material: MagnetoelasticMaterial,
        film_thickness_m: f64,
        h_ext: f64,
    ) -> Result<Self> {
        if film_thickness_m <= 0.0 {
            return Err(invalid_param(
                "film_thickness_m",
                "ferromagnetic film thickness must be positive",
            ));
        }
        if h_ext < 0.0 {
            return Err(invalid_param(
                "h_ext",
                "external bias field must be non-negative",
            ));
        }
        Ok(Self {
            saw,
            material,
            film_thickness_m,
            h_ext,
        })
    }

    /// Effective magnetoelastic driving field amplitude H_me [A/m].
    ///
    /// Derived from the leading-order magnetoelastic energy density under
    /// strain ε_xx = ε₀:
    ///
    ///   H_me = |2 B₁ ε₀| / (μ₀ M_s)
    ///
    /// This field drives magnetisation precession in the same way a microwave
    /// field h_rf does in conventional FMR.
    pub fn magnetoelastic_field_amplitude(&self) -> f64 {
        (2.0 * self.material.b1.abs() * self.saw.strain_amplitude) / (MU_0 * self.material.ms)
    }

    /// FMR frequency f_FMR [Hz] at the current bias field.
    ///
    /// Uses the simplified Kittel formula (out-of-plane saturation / uniform mode):
    ///   f_FMR = γ / (2π) · (H_ext + H_k)
    pub fn fmr_frequency(&self) -> f64 {
        self.material.fmr_frequency_hz(self.h_ext)
    }

    /// Frequency detuning Δf = f_SAW − f_FMR [Hz].
    ///
    /// Δf = 0 → exact acoustic FMR resonance.
    /// |Δf| ≫ α f_FMR → far off-resonance, exponentially small response.
    pub fn resonance_detuning(&self) -> f64 {
        self.saw.frequency_hz - self.fmr_frequency()
    }

    /// Lorentzian response function L(ω) [1/(rad/s)].
    ///
    /// The driven magnetisation amplitude is proportional to:
    ///   L(ω) = 1 / (ω_FMR · sqrt[(ω/ω_FMR − 1)² + α²])
    ///
    /// At exact resonance (ω = ω_FMR) this reduces to L = 1/(α ω_FMR) = Q/ω_FMR,
    /// giving the peak precession amplitude.
    ///
    /// # Physical note
    /// The square-root form is correct for the *amplitude* Lorentzian (not the
    /// power Lorentzian which has a squared denominator).  This is consistent
    /// with how a driven harmonic oscillator responds at first order.
    pub fn lorentzian_response(&self) -> f64 {
        let omega_fmr = 2.0 * PI * self.fmr_frequency();
        let omega_saw = self.saw.angular_frequency();
        let alpha = self.material.alpha;
        let detuning_sq = (omega_saw / omega_fmr - 1.0).powi(2) + alpha.powi(2);
        1.0 / (detuning_sq.sqrt() * omega_fmr)
    }

    /// Magnetisation precession cone angle θ_cone [rad].
    ///
    /// The driven precession amplitude is:
    ///   θ_cone = H_me · L(ω)
    ///
    /// where H_me is the magnetoelastic field amplitude [A/m] and L(ω) is the
    /// Lorentzian response [1/(rad/s)].  The product has units of
    /// [A/m · s/rad] = [A·s/(m·rad)] which, via γ = rad/(s·T) and
    /// μ₀ M_s, simplifies to a dimensionless angle (radians).
    ///
    /// The result is clamped to [0, π/4] to maintain physical validity
    /// (small-angle precession assumption underlying the Lorentzian model).
    pub fn precession_cone_angle_rad(&self) -> f64 {
        let h_me = self.magnetoelastic_field_amplitude();
        let theta_raw = h_me * self.lorentzian_response();
        // Clamp to physically valid small-angle regime [0, π/4]
        theta_raw.clamp(0.0, PI / 4.0)
    }

    /// Acoustic spin current density J_s [A/m²].
    ///
    /// A precessing magnetisation in an FM/NM bilayer pumps a DC spin current
    /// via the spin-pumping mechanism (Tserkovnyak et al., 2002):
    ///
    ///   J_s = (ℏ / 4π) · g_r · ω_SAW · sin²(θ_cone) / 2
    ///
    /// where g_r [m⁻²] is the real part of the spin-mixing conductance per
    /// unit area at the FM/NM interface.
    ///
    /// # Arguments
    /// * `g_r` — Spin-mixing conductance (real part) [m⁻²].
    ///   Typical values: ~10¹⁸–10¹⁹ m⁻² for Py/Pt.
    pub fn spin_current_density(&self, g_r: f64) -> f64 {
        let theta = self.precession_cone_angle_rad();
        (HBAR / (4.0 * PI)) * g_r * self.saw.angular_frequency() * theta.sin().powi(2) / 2.0
    }

    /// Resonant SAW power absorption per unit area P_abs [W/m²].
    ///
    /// The FM film absorbs SAW power when the magnetisation precesses under
    /// the alternating magnetoelastic torque.  The Gilbert-damping power
    /// dissipation rate for a precessing thin film is:
    ///
    ///   P_abs = (μ₀ / 2) · M_s · α · ω_FMR · θ_cone² · d_film
    ///
    /// where d_film [m] is the film thickness (converting volume density to
    /// areal density).
    ///
    /// Strictly valid for θ_cone ≪ 1 (small-angle, linear response).
    pub fn resonant_absorption(&self) -> f64 {
        let theta = self.precession_cone_angle_rad();
        let omega_fmr = 2.0 * PI * self.fmr_frequency();
        0.5 * MU_0
            * self.material.ms
            * self.material.alpha
            * omega_fmr
            * theta.powi(2)
            * self.film_thickness_m
    }

    /// Bias field H_ext [A/m] required to tune FMR onto the SAW frequency.
    ///
    /// Inverts the simplified Kittel relation f_FMR = γ μ₀ (H_ext + H_k)/(2π):
    ///   H_ext = ω_SAW / (γ μ₀) − H_k
    ///
    /// The factor μ₀ is needed because γ is in rad/(s·T) while H is in A/m.
    ///
    /// Returns 0 if the formula would give a negative result (occurs when
    /// f_SAW < f_FMR(H=0), i.e. the zero-field FMR already exceeds the SAW
    /// frequency — no positive bias can bring them into resonance in this
    /// simplified model).
    pub fn acoustic_fmr_condition_h_ext(&self) -> f64 {
        let h_ext_raw =
            self.saw.angular_frequency() / (GAMMA * MU_0) - self.material.anisotropy_field();
        // External field cannot be negative in this model
        h_ext_raw.max(0.0)
    }
}

// =============================================================================
// SawSpinWaveExcitation
// =============================================================================

/// SAW excitation of spin waves (magnons) in a magnetic film.
///
/// A SAW with wavenumber k_SAW imprints a periodic magnetoelastic potential
/// with the *same* wavenumber on the magnetic film.  This can resonantly
/// excite a spin wave with k_SW = k_SAW when the dispersion relations
/// intersect — a phenomenon exploited for efficient acoustic-to-magnonic
/// transduction.
///
/// The excited spin-wave frequency follows the Kittel dispersion for the
/// uniform mode (k = 0 limit) or the dipole-exchange dispersion for finite k.
/// Here we use the out-of-plane simplified dispersion for pedagogical clarity:
///   ω_SW(k) ≈ γ (H_ext + H_k)   [independent of k for uniform mode]
///
/// Phase-velocity matching Δv = v_SAW − v_SW(k_SAW) determines whether the
/// coupling is resonant (Δv → 0) or dispersive.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SawSpinWaveExcitation {
    /// SAW source.
    pub saw: SawSource,

    /// Magnetic film material.
    pub material: MagnetoelasticMaterial,
}

impl SawSpinWaveExcitation {
    /// Create a new `SawSpinWaveExcitation` coupling.
    pub fn new(saw: SawSource, material: MagnetoelasticMaterial) -> Self {
        Self { saw, material }
    }

    /// Excited spin-wave wavenumber [rad/m].
    ///
    /// Momentum conservation forces:
    ///   k_SW = k_SAW = 2π f / v_SAW
    ///
    /// This is the key advantage of SAW transduction: the spin-wave wavevector
    /// is set precisely by the IDT geometry (through the SAW wavelength), not
    /// by an external antenna.
    pub fn excited_spinwave_k(&self) -> f64 {
        self.saw.substrate.wavenumber(self.saw.frequency_hz)
    }

    /// Excited spin-wave frequency f_SW [Hz] at the given bias field.
    ///
    /// In the simplified out-of-plane Kittel picture the spin-wave frequency
    /// for the uniform mode is:
    ///   ω_SW = γ μ₀ (H_ext + H_k)
    ///
    /// where H_ext and H_k are in A/m and μ₀ converts to Tesla.
    /// This is independent of k in the uniform-mode limit; the full
    /// dipole-exchange dispersion would add a k-dependent exchange term.
    ///
    /// # Arguments
    /// * `h_ext` — External bias field [A/m].
    ///
    /// # Returns
    /// Spin-wave frequency f_SW [Hz].
    pub fn excited_spinwave_frequency(&self, h_ext: f64) -> f64 {
        let omega_sw = GAMMA * MU_0 * (h_ext + self.material.anisotropy_field());
        omega_sw / (2.0 * PI)
    }

    /// Phase-velocity mismatch Δv [m/s] between SAW and spin wave.
    ///
    /// A large mismatch means the SAW and spin wave are not in resonance; the
    /// interaction is then dispersive and the coupling is weak.  Near
    /// Δv → 0 the magneto-acoustic interaction length diverges.
    ///
    /// Computed as:
    ///   Δv = |f_SAW − f_SW| · λ_SAW
    ///
    /// which has dimensions [Hz · m] = [m/s].
    ///
    /// # Arguments
    /// * `h_ext` — External bias field [A/m].
    pub fn phase_velocity_mismatch(&self, h_ext: f64) -> f64 {
        let f_saw = self.saw.frequency_hz;
        let f_sw = self.excited_spinwave_frequency(h_ext);
        (f_saw - f_sw).abs() * self.saw.wavelength_m()
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------------------------------------------------------
    // PiezoSubstrate
    // -------------------------------------------------------------------------

    #[test]
    fn test_saw_wavelength_linbo3_1ghz() {
        let sub = PiezoSubstrate::linbo3();
        let lambda = sub.wavelength_m(1.0e9);
        // Expected: 3985 m/s / 1e9 Hz = 3.985 μm
        let expected = 3.985e-6;
        let rel_err = (lambda - expected).abs() / expected;
        assert!(
            rel_err < 0.01,
            "LiNbO₃ wavelength at 1 GHz should be ~3.985 μm, got {:.4e} m (rel_err={:.4})",
            lambda,
            rel_err
        );
    }

    #[test]
    fn test_saw_wavenumber_linbo3() {
        let sub = PiezoSubstrate::linbo3();
        let k = sub.wavenumber(1.0e9);
        // k = 2π / λ = 2π / 3.985e-6
        let expected = 2.0 * PI / 3.985e-6;
        let rel_err = (k - expected).abs() / expected;
        assert!(
            rel_err < 1.0e-6,
            "Wavenumber mismatch: got {:.6e} expected {:.6e}",
            k,
            expected
        );
    }

    #[test]
    fn test_penetration_depth_is_lambda_over_2pi() {
        let sub = PiezoSubstrate::linbo3();
        let f = 2.0e9;
        let depth = sub.penetration_depth_m(f);
        let lambda = sub.wavelength_m(f);
        let expected = lambda / (2.0 * PI);
        let rel_err = (depth - expected).abs() / expected;
        assert!(
            rel_err < 1.0e-10,
            "Penetration depth should equal λ/(2π), rel_err = {:.2e}",
            rel_err
        );
    }

    #[test]
    fn test_gaas_slower_than_linbo3() {
        let linbo3 = PiezoSubstrate::linbo3();
        let gaas = PiezoSubstrate::gaas();
        assert!(
            gaas.saw_velocity < linbo3.saw_velocity,
            "GaAs SAW velocity ({}) should be less than LiNbO₃ ({})",
            gaas.saw_velocity,
            linbo3.saw_velocity
        );
    }

    // -------------------------------------------------------------------------
    // MagnetoelasticMaterial
    // -------------------------------------------------------------------------

    #[test]
    fn test_magnetoelastic_field_nickel() {
        // H_me = |2 B₁ ε₀| / (μ₀ M_s)
        // Ni: B₁ = -8.5e6, M_s = 490e3, ε₀ = 1e-4
        // H_me = 2 * 8.5e6 * 1e-4 / (4π×10⁻⁷ * 490e3)
        //      = 1700 / (4π×10⁻⁷ * 490e3)
        //      ≈ 1700 / 0.6158 ≈ 2760 A/m  (of order 10³ A/m, < 10000 A/m)
        let saw = SawSource::linbo3_1ghz();
        let mat = MagnetoelasticMaterial::nickel();
        let dev = SawMagnetoacoustics::new(saw, mat, 20.0e-9, 0.0).expect("valid parameters");
        let h_me = dev.magnetoelastic_field_amplitude();
        // Physical: should be a positive field in the A/m range, < 1 MA/m
        assert!(h_me > 0.0, "H_me must be positive, got {}", h_me);
        assert!(
            h_me < 1.0e6,
            "H_me should be < 1 MA/m for ε₀=1e-4, got {:.3e} A/m",
            h_me
        );
    }

    #[test]
    fn test_fmr_frequency_permalloy_zero_field() {
        // At H_ext = 0 the FMR is driven by the anisotropy field alone.
        // H_k(Py) = 2*200/(4π×10⁻⁷ * 800e3) ≈ 398 A/m
        // f_FMR = γ/(2π) * 398 ≈ 11.2 MHz  — very low; real Py also has a
        // demagnetising contribution not included in this simplified model.
        let mat = MagnetoelasticMaterial::permalloy();
        let f = mat.fmr_frequency_hz(0.0);
        // The simplified model gives a low frequency; we just check it is > 0
        // and < 100 GHz (any value ≫ 100 GHz would indicate a unit error).
        assert!(
            f > 0.0,
            "FMR frequency must be positive at zero field, got {}",
            f
        );
        assert!(
            f < 100.0e9,
            "FMR frequency should be < 100 GHz, got {:.3e} Hz",
            f
        );
    }

    #[test]
    fn test_quality_factor_equals_one_over_alpha() {
        let mat = MagnetoelasticMaterial::nickel();
        let q = mat.quality_factor(0.0);
        let expected = 1.0 / mat.alpha;
        assert!(
            (q - expected).abs() < 1.0e-10,
            "Q factor should be 1/α = {}, got {}",
            expected,
            q
        );
    }

    // -------------------------------------------------------------------------
    // SawSource
    // -------------------------------------------------------------------------

    #[test]
    fn test_saw_source_new_rejects_negative_frequency() {
        let result = SawSource::new(PiezoSubstrate::linbo3(), -1.0e9, 1.0e-4);
        assert!(result.is_err(), "Negative frequency should be rejected");
    }

    #[test]
    fn test_saw_source_new_rejects_large_strain() {
        let result = SawSource::new(PiezoSubstrate::linbo3(), 1.0e9, 0.5);
        assert!(result.is_err(), "Strain ≥ 0.1 should be rejected");
    }

    #[test]
    fn test_strain_at_surface_matches_amplitude() {
        // At z=0, x=λ/4, t=0: strain = ε₀ * sin(k * λ/4) * 1 = ε₀ * sin(π/2) = ε₀
        let saw = SawSource::linbo3_1ghz();
        let lambda = saw.wavelength_m();
        let x_quarter = lambda / 4.0;
        let strain = saw.strain_at_point(x_quarter, 0.0, 0.0);
        let expected = saw.strain_amplitude;
        let rel_err = (strain - expected).abs() / expected;
        assert!(
            rel_err < 1.0e-10,
            "Surface strain at λ/4 should equal ε₀ = {:.4e}, got {:.4e} (rel_err={:.2e})",
            expected,
            strain,
            rel_err
        );
    }

    #[test]
    fn test_strain_decays_with_depth() {
        let saw = SawSource::linbo3_1ghz();
        let lambda = saw.wavelength_m();
        let x_quarter = lambda / 4.0;
        let s_surface = saw.strain_at_point(x_quarter, 0.0, 0.0);
        let depth = saw.substrate.penetration_depth_m(saw.frequency_hz);
        let s_depth = saw.strain_at_point(x_quarter, depth, 0.0);
        // At depth δ = 1/k the strain should decay by factor exp(-1) ≈ 0.368
        let ratio = s_depth / s_surface;
        let expected_ratio = (-1.0_f64).exp();
        assert!(
            (ratio - expected_ratio).abs() < 1.0e-6,
            "Strain ratio at depth δ should be e⁻¹ = {:.4}, got {:.4}",
            expected_ratio,
            ratio
        );
    }

    // -------------------------------------------------------------------------
    // SawMagnetoacoustics
    // -------------------------------------------------------------------------

    #[test]
    fn test_resonance_at_zero_detuning() {
        // Build a device on-resonance: choose H_ext so that f_FMR = f_SAW.
        let saw = SawSource::linbo3_3ghz();
        let mat = MagnetoelasticMaterial::permalloy();

        // First compute the resonant H_ext
        let helper = SawMagnetoacoustics::new(saw.clone(), mat.clone(), 20.0e-9, 0.0)
            .expect("valid parameters");
        let h_res = helper.acoustic_fmr_condition_h_ext();

        let on_res = SawMagnetoacoustics::new(saw.clone(), mat.clone(), 20.0e-9, h_res)
            .expect("valid parameters");
        // Off-resonance: shift H by 10× H_k
        let h_off = h_res + 100.0 * mat.anisotropy_field() + 1.0e5;
        let off_res = SawMagnetoacoustics::new(saw, mat, 20.0e-9, h_off).expect("valid parameters");

        let theta_on = on_res.precession_cone_angle_rad();
        let theta_off = off_res.precession_cone_angle_rad();
        assert!(
            theta_on > theta_off,
            "On-resonance cone angle ({:.4e} rad) should exceed off-resonance ({:.4e} rad)",
            theta_on,
            theta_off
        );
    }

    #[test]
    fn test_precession_angle_grows_with_strain() {
        let sub = PiezoSubstrate::linbo3();
        let mat = MagnetoelasticMaterial::permalloy();

        let saw_small = SawSource::new(sub.clone(), 1.0e9, 1.0e-5).expect("valid");
        let saw_large = SawSource::new(sub, 1.0e9, 1.0e-3).expect("valid");

        let dev_small = SawMagnetoacoustics::new(saw_small, mat.clone(), 20.0e-9, 0.0)
            .expect("valid parameters");
        let dev_large =
            SawMagnetoacoustics::new(saw_large, mat, 20.0e-9, 0.0).expect("valid parameters");

        let theta_small = dev_small.precession_cone_angle_rad();
        let theta_large = dev_large.precession_cone_angle_rad();

        assert!(
            theta_large > theta_small,
            "Larger strain (1e-3) should give larger cone angle ({:.4e}) than small strain (1e-5) ({:.4e})",
            theta_large,
            theta_small
        );
    }

    #[test]
    fn test_spin_current_positive_at_resonance() {
        let saw = SawSource::linbo3_3ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let helper =
            SawMagnetoacoustics::new(saw.clone(), mat.clone(), 20.0e-9, 0.0).expect("valid");
        let h_res = helper.acoustic_fmr_condition_h_ext();
        let device = SawMagnetoacoustics::new(saw, mat, 20.0e-9, h_res).expect("valid parameters");

        // Typical Py/Pt spin-mixing conductance g_r ≈ 4.18 × 10¹⁸ m⁻²
        let g_r = 4.18e18_f64;
        let j_s = device.spin_current_density(g_r);

        assert!(
            j_s > 0.0,
            "Spin current density must be positive at resonance, got {:.3e} A/m²",
            j_s
        );
    }

    #[test]
    fn test_spin_current_zero_for_zero_g_r() {
        let saw = SawSource::linbo3_1ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let device = SawMagnetoacoustics::new(saw, mat, 20.0e-9, 0.0).expect("valid");
        let j_s = device.spin_current_density(0.0);
        assert!(
            j_s.abs() < 1.0e-60,
            "Spin current with g_r=0 should be ~0, got {:.3e}",
            j_s
        );
    }

    #[test]
    fn test_acoustic_fmr_field_positive_for_linbo3_1ghz() {
        // At 1 GHz the SAW frequency is low; for Py with α = 0.01 we might need
        // a small positive H_ext to hit resonance.
        let saw = SawSource::linbo3_1ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let device = SawMagnetoacoustics::new(saw, mat, 20.0e-9, 0.0).expect("valid");
        let h_cond = device.acoustic_fmr_condition_h_ext();
        // Must be non-negative (clamped)
        assert!(
            h_cond >= 0.0,
            "acoustic_fmr_condition_h_ext must be ≥ 0, got {}",
            h_cond
        );
    }

    #[test]
    fn test_resonant_absorption_positive() {
        let saw = SawSource::linbo3_3ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let helper =
            SawMagnetoacoustics::new(saw.clone(), mat.clone(), 20.0e-9, 0.0).expect("valid");
        let h_res = helper.acoustic_fmr_condition_h_ext();
        let device = SawMagnetoacoustics::new(saw, mat, 20.0e-9, h_res).expect("valid");
        let p_abs = device.resonant_absorption();
        assert!(
            p_abs > 0.0,
            "Resonant absorption must be positive, got {:.3e} W/m²",
            p_abs
        );
    }

    #[test]
    fn test_film_thickness_validation() {
        let saw = SawSource::linbo3_1ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let result = SawMagnetoacoustics::new(saw, mat, -5.0e-9, 0.0);
        assert!(
            result.is_err(),
            "Negative film thickness should be rejected"
        );
    }

    #[test]
    fn test_h_ext_validation() {
        let saw = SawSource::linbo3_1ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let result = SawMagnetoacoustics::new(saw, mat, 20.0e-9, -1.0);
        assert!(result.is_err(), "Negative H_ext should be rejected");
    }

    // -------------------------------------------------------------------------
    // SawSpinWaveExcitation
    // -------------------------------------------------------------------------

    #[test]
    fn test_excited_spinwave_k_equals_saw_k() {
        let saw = SawSource::linbo3_1ghz();
        let k_saw = saw.substrate.wavenumber(saw.frequency_hz);
        let mat = MagnetoelasticMaterial::permalloy();
        let exc = SawSpinWaveExcitation::new(saw, mat);
        let k_sw = exc.excited_spinwave_k();
        assert!(
            (k_sw - k_saw).abs() < 1.0e-6,
            "Spin-wave k ({:.4e}) must equal SAW k ({:.4e})",
            k_sw,
            k_saw
        );
    }

    #[test]
    fn test_spinwave_frequency_increases_with_field() {
        let saw = SawSource::linbo3_1ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let exc = SawSpinWaveExcitation::new(saw, mat);
        let f_low = exc.excited_spinwave_frequency(0.0);
        let f_high = exc.excited_spinwave_frequency(1.0e5);
        assert!(
            f_high > f_low,
            "Spin-wave frequency should increase with bias field: f(0)={:.3e} f(1e5)={:.3e}",
            f_low,
            f_high
        );
    }

    #[test]
    fn test_phase_velocity_mismatch_zero_at_resonance() {
        // When f_SAW = f_SW, the mismatch should be (nearly) zero.
        let saw = SawSource::linbo3_1ghz();
        let mat = MagnetoelasticMaterial::permalloy();
        let exc = SawSpinWaveExcitation::new(saw.clone(), mat.clone());

        // Find H_ext such that f_SW = f_SAW
        // f_SW = GAMMA * MU_0 * (H_ext + H_k) / (2π) = f_SAW
        // H_ext = f_SAW * 2π / (GAMMA * MU_0) - H_k
        let h_k = mat.anisotropy_field();
        let h_res = saw.frequency_hz * 2.0 * PI / (GAMMA * MU_0) - h_k;
        let h_res_clamped = h_res.max(0.0);

        let mismatch = exc.phase_velocity_mismatch(h_res_clamped);
        // At resonance the mismatch should be small (bounded by floating-point)
        assert!(
            mismatch < 1.0e3,
            "Phase velocity mismatch at resonance should be near zero, got {:.3e} m/s",
            mismatch
        );
    }

    #[test]
    fn test_lorentzian_response_peaks_at_resonance() {
        let saw = SawSource::linbo3_3ghz();
        let mat = MagnetoelasticMaterial::permalloy();

        let helper =
            SawMagnetoacoustics::new(saw.clone(), mat.clone(), 20.0e-9, 0.0).expect("valid");
        let h_res = helper.acoustic_fmr_condition_h_ext();
        let on_res =
            SawMagnetoacoustics::new(saw.clone(), mat.clone(), 20.0e-9, h_res).expect("valid");
        // Far detuned: add 10× M_s worth of field
        let h_far = h_res + 10.0 * mat.ms;
        let off_res = SawMagnetoacoustics::new(saw, mat, 20.0e-9, h_far).expect("valid");

        let l_on = on_res.lorentzian_response();
        let l_off = off_res.lorentzian_response();
        assert!(
            l_on > l_off,
            "Lorentzian response on-resonance ({:.3e}) should exceed off-resonance ({:.3e})",
            l_on,
            l_off
        );
    }

    #[test]
    fn test_cobalt_higher_fmr_than_permalloy() {
        // Co has larger anisotropy → higher H_k → higher f_FMR at H_ext = 0
        let py = MagnetoelasticMaterial::permalloy();
        let co = MagnetoelasticMaterial::cobalt();
        let f_py = py.fmr_frequency_hz(0.0);
        let f_co = co.fmr_frequency_hz(0.0);
        assert!(
            f_co > f_py,
            "Co FMR ({:.3e} Hz) should exceed Py FMR ({:.3e} Hz) at zero field",
            f_co,
            f_py
        );
    }
}
