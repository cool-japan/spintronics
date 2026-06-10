//! AC Spin Pumping and Spin Battery
//!
//! Implements the theory of spin current generation from a precessing ferromagnet
//! driven at or near ferromagnetic resonance (FMR), including:
//!
//! - **DC spin pumping** (Tserkovnyak 2002): rectified spin current from steady
//!   precession, `J_s^DC ∝ ⟨m × dm/dt⟩`.
//! - **2ω AC component**: oscillatory spin current at twice the drive frequency,
//!   equal in magnitude to the DC component at resonance.
//! - **Backflow correction** (Tserkovnyak RMP 2005): finite NM thickness reduces
//!   the effective spin mixing conductance via spin-diffusion back-action.
//! - **Spin accumulation / spin battery**: DC spin current builds a non-equilibrium
//!   spin chemical potential in the NM.
//! - **Gilbert damping enhancement**: spin pumping adds a radiative damping channel
//!   proportional to G_r.
//! - **ISHE readout**: the ISHE converts the spin battery output into a measurable
//!   open-circuit voltage.
//!
//! # Key References
//!
//! - Y. Tserkovnyak et al., "Enhanced Gilbert Damping in Thin Ferromagnetic Films",
//!   *Phys. Rev. Lett.* **88**, 117601 (2002).
//! - Y. Tserkovnyak et al., "Nonlocal magnetization dynamics in ferromagnetic
//!   heterostructures", *Rev. Mod. Phys.* **77**, 1375 (2005).
//! - E. Saitoh et al., "Conversion of spin current into charge current at room
//!   temperature: Inverse spin-Hall effect", *Appl. Phys. Lett.* **88**, 182509 (2006).
//! - O. Mosendz et al., "Quantifying Spin Hall Angles from Spin Pumping",
//!   *Phys. Rev. Lett.* **104**, 046601 (2010).

use std::f64::consts::PI;

use crate::constants::{GAMMA, HBAR, KB, MU_0};
use crate::effect::ishe::InverseSpinHall;
use crate::error::{invalid_param, Result};
use crate::material::{Ferromagnet, SpinInterface};
use crate::transport::diffusion::SpinDiffusion;
use crate::vector3::Vector3;

// ---------------------------------------------------------------------------
// AcSpinPumping
// ---------------------------------------------------------------------------

/// AC/DC spin pumping at a ferromagnet / normal-metal interface.
///
/// Captures the full Tserkovnyak picture: the precessing FM acts as a *spin
/// battery* that injects a spin current into the adjacent NM. The effective
/// pumping is reduced by spin-diffusion backflow from the finite-thickness NM
/// slab, and the extracted spin current depends on whether the system is driven
/// at FMR (where DC and 2ω components are equal) or off-resonance.
///
/// # Physical picture
///
/// At FMR the cone-angle θ is set by the rf-drive amplitude and the linewidth
/// `Δω = α ω_FMR`. The pumped spin-current density is
///
/// ```text
/// J_s^DC = (ℏ / 8π) · G_r_eff · ω_FMR · sin²θ        [A m⁻²]
/// J_s^{2ω} = J_s^DC                                   (at resonance)
/// ```
///
/// The backflow-corrected mixing conductance is
///
/// ```text
/// G_r_eff = G_r / (1 + G_r / G_NM)
/// G_NM    = σ_NM · tanh(t_NM / λ_sf) / λ_sf           [Ω⁻¹ m⁻²]
/// ```
#[derive(Debug, Clone)]
pub struct AcSpinPumping {
    /// FM/NM interface parameters (spin-mixing conductance, area).
    pub interface: SpinInterface,
    /// Spin-diffusion properties of the NM layer (λ_sf, σ).
    pub diffusion: SpinDiffusion,
    /// NM layer thickness \[m\].
    pub t_nm: f64,
    /// FM layer thickness \[m\].
    pub t_fm: f64,
}

impl AcSpinPumping {
    /// Construct a new `AcSpinPumping` with explicit parameters.
    ///
    /// # Errors
    /// Returns `Err` if `t_nm` or `t_fm` are non-positive.
    pub fn new(
        interface: SpinInterface,
        diffusion: SpinDiffusion,
        t_nm: f64,
        t_fm: f64,
    ) -> Result<Self> {
        if t_nm <= 0.0 {
            return Err(invalid_param("t_nm", "NM thickness must be positive"));
        }
        if t_fm <= 0.0 {
            return Err(invalid_param("t_fm", "FM thickness must be positive"));
        }
        Ok(Self {
            interface,
            diffusion,
            t_nm,
            t_fm,
        })
    }

    /// YIG (5 μm) / Pt (7 nm) preset — the canonical spin-pumping system.
    ///
    /// Parameters are representative of Saitoh-group experiments.
    pub fn yig_pt() -> Self {
        Self {
            interface: SpinInterface::yig_pt(),
            diffusion: SpinDiffusion::platinum(),
            t_nm: 7.0e-9, // 7 nm Pt
            t_fm: 5.0e-6, // 5 μm YIG
        }
    }

    // -----------------------------------------------------------------------
    // Backflow correction
    // -----------------------------------------------------------------------

    /// Spin conductance of the NM slab \[Ω⁻¹ m⁻²\].
    ///
    /// Accounts for the finite spin-diffusion length in the NM:
    ///
    /// ```text
    /// G_NM = σ_NM · tanh(t_NM / λ_sf) / λ_sf
    /// ```
    ///
    /// For a slab much thicker than λ_sf, tanh → 1 and G_NM → σ / λ_sf.
    /// For a very thin slab, G_NM → σ · t_NM / λ_sf².
    pub fn g_nm(&self) -> f64 {
        let ratio = self.t_nm / self.diffusion.lambda_sf;
        self.diffusion.sigma * ratio.tanh() / self.diffusion.lambda_sf
    }

    /// Backflow-corrected real part of spin-mixing conductance \[Ω⁻¹ m⁻²\].
    ///
    /// The NM spin conductance `G_NM` acts in parallel with the interface
    /// conductance `G_r`, reducing the net spin injection:
    ///
    /// ```text
    /// G_r_eff = G_r / (1 + G_r / G_NM)
    /// ```
    pub fn effective_g_r(&self) -> f64 {
        let g_r = self.interface.g_r;
        let g_nm = self.g_nm();
        g_r / (1.0 + g_r / g_nm)
    }

    // -----------------------------------------------------------------------
    // Spin current
    // -----------------------------------------------------------------------

    /// Instantaneous DC spin-current density pumped across the interface.
    ///
    /// Following Tserkovnyak (2002), the pumped spin current per unit area is
    ///
    /// ```text
    /// J_s = (ℏ / 4π) · G_r_eff · (m × dm/dt)
    /// ```
    ///
    /// The returned vector is oriented along **m × dm/dt** (the instantaneous
    /// spin-polarisation direction) with magnitude `(ℏ/4π) · G_r_eff · |m × dm/dt|`.
    ///
    /// # Arguments
    /// * `m`    — unit magnetisation vector (dimensionless).
    /// * `dmdt` — time derivative of the magnetisation \[s⁻¹\] (already in units
    ///   where |m| = 1, so |dm/dt| ≈ ω sin θ at precession angle θ).
    ///
    /// # Returns
    /// Spin-current density vector \[ℏ · A m⁻² / (2e)\] — the conventional
    /// "spin-current" unit used in Tserkovnyak's formulation. Numerically this
    /// is `≈ J_charge / (2e/ℏ)`, i.e. the factor `ℏ/(4π)` carries the units.
    pub fn dc_spin_current(&self, m: Vector3<f64>, dmdt: Vector3<f64>) -> Vector3<f64> {
        let cross = m.cross(&dmdt);
        cross * (HBAR / (4.0 * PI) * self.effective_g_r())
    }

    /// Magnitude of the DC spin current at ferromagnetic resonance \[A m⁻²\].
    ///
    /// Under uniform circular precession at cone angle θ:
    ///
    /// ```text
    /// |J_s^DC| = (ℏ / 8π) · G_r_eff · ω_FMR · sin²θ
    /// ```
    ///
    /// The factor of `1/8` (vs `1/4` in the general formula) arises from
    /// time-averaging the rectified signal `⟨|m × dm/dt|⟩ = ω · sin²θ / 2`.
    ///
    /// # Arguments
    /// * `omega_fmr`  — FMR angular frequency \[rad s⁻¹\].
    /// * `cone_angle` — precession cone half-angle θ \[rad\].
    pub fn dc_magnitude_at_fmr(&self, omega_fmr: f64, cone_angle: f64) -> f64 {
        let sin_sq = cone_angle.sin().powi(2);
        (HBAR / (8.0 * PI)) * self.effective_g_r() * omega_fmr * sin_sq
    }

    /// Magnitude of the 2ω AC spin-current component \[A m⁻²\].
    ///
    /// At resonance, the oscillatory component at twice the drive frequency
    /// has the same magnitude as the rectified (DC) component:
    ///
    /// ```text
    /// |J_s^{2ω}| = |J_s^DC|
    /// ```
    ///
    /// This equivalence holds exactly for circular precession and is an
    /// important experimental signature that distinguishes spin-pumping from
    /// other spin-current sources.
    ///
    /// # Arguments
    /// * `omega_fmr`  — FMR angular frequency \[rad s⁻¹\].
    /// * `cone_angle` — precession cone half-angle θ \[rad\].
    pub fn ac_component_2omega(&self, omega_fmr: f64, cone_angle: f64) -> f64 {
        // Identical to DC magnitude at resonance — they are degenerate for
        // circular precession.
        self.dc_magnitude_at_fmr(omega_fmr, cone_angle)
    }

    // -----------------------------------------------------------------------
    // Damping enhancement
    // -----------------------------------------------------------------------

    /// Additional Gilbert damping constant from spin pumping (dimensionless).
    ///
    /// The radiation of spin angular momentum into the NM opens a new
    /// dissipation channel, effectively increasing the measured FMR linewidth:
    ///
    /// ```text
    /// Δα = (γ · ℏ · G_r_eff) / (4π · μ₀ · M_s · t_FM)
    /// ```
    ///
    /// Here `γ = GAMMA` is the electron gyromagnetic ratio, `M_s` is the
    /// saturation magnetisation of the FM, and `t_FM` is its thickness. The
    /// backflow-corrected `G_r_eff` is used.
    ///
    /// # Arguments
    /// * `material` — ferromagnet whose damping is being enhanced.
    pub fn damping_enhancement(&self, material: &Ferromagnet) -> f64 {
        let g_eff = self.effective_g_r();
        (GAMMA * HBAR * g_eff) / (4.0 * PI * MU_0 * material.ms * self.t_fm)
    }

    // -----------------------------------------------------------------------
    // Spin accumulation (spin battery)
    // -----------------------------------------------------------------------

    /// Spin accumulation (spin chemical potential difference) in the NM \[J\].
    ///
    /// The DC spin current injected into the NM builds up a non-equilibrium
    /// spin population. In the linear-response (spin-battery) regime the
    /// spin accumulation is:
    ///
    /// ```text
    /// μ_s = J_s^DC · R_sf
    /// ```
    ///
    /// where the spin-resistance of the NM slab is
    ///
    /// ```text
    /// R_sf = λ_sf / (σ_NM · A)        [Ω]
    /// ```
    ///
    /// Multiplying current density (A m⁻²) by `A` (m²) gives current (A),
    /// then multiplying by `R_sf` (Ω) gives a voltage (V = J/C).  In spin
    /// physics the natural unit is the spin chemical potential in Joules,
    /// obtained by treating the charge-equivalent voltage as an energy
    /// `μ_s \[J\] = J_s [A/m²] · A \[m²\] · R_sf [Ω] / A [m²]`
    ///         `= J_s [A/m²] · λ_sf / σ_NM`.
    ///
    /// The returned value carries SI units of \[J\] (divide by `E_CHARGE` for
    /// electron-volt equivalent).
    ///
    /// # Arguments
    /// * `omega_fmr`  — FMR angular frequency \[rad s⁻¹\].
    /// * `cone_angle` — precession cone half-angle θ \[rad\].
    pub fn spin_accumulation_dc(&self, omega_fmr: f64, cone_angle: f64) -> f64 {
        let j_s = self.dc_magnitude_at_fmr(omega_fmr, cone_angle);
        // R_sf per unit area: λ_sf / σ  [Ω · m² / m² = Ω ... but we work in
        // volumetric units so the area cancels and we get J/(A/m²) = J·m²/A]
        // μ_s [J] = J_s [A/m²] × (λ_sf / σ)   where λ_sf/σ has dimensions of
        // (m) / (S/m) = m / (A/(V·m)) = V·m² / A = Ω·m² → times A/m² = V = J/C
        // → treating as energy per unit charge, this is already in volts;
        //   we return in joules by noting 1 J/C × 1 C/charge = energy.
        //   In practice the literature quotes μ_s in μV or μeV; we return SI J.
        j_s * self.diffusion.lambda_sf / self.diffusion.sigma
    }
}

// ---------------------------------------------------------------------------
// SpinBattery
// ---------------------------------------------------------------------------

/// Spin battery: an `AcSpinPumping` source read out by an ISHE detector.
///
/// The pumped spin current flows through the NM layer, builds a spin
/// accumulation (spin voltage), and is converted to an open-circuit charge
/// voltage by the inverse spin-Hall effect in the same NM strip.
///
/// # Geometry
///
/// ```text
///  ┌──────── FM (YIG) ────────┐
///  │  precessing magnetization │
///  └──────── interface ────────┘
///  ┌──────── NM (Pt) strip ───┐  ← J_s flows ⊥ interface (z)
///  │  ┌────────────────────┐  │  ← σ̂ in plane (e.g. x)
///  │  │  V_ISHE measured   │  │  ← E_ISHE along y
///  │  └────────────────────┘  │
///  └──────────────────────────┘
///         ←── sample_length ──→
///               ↑ sample_width
/// ```
///
/// The ISHE voltage is measured between the two ends of the NM strip along
/// the direction perpendicular to both J_s and σ̂, giving
///
/// ```text
/// V_ISHE = θ_SH · ρ_NM · J_s^DC · L
/// ```
#[derive(Debug, Clone)]
pub struct SpinBattery {
    /// The spin-pumping source (FM/NM interface + NM diffusion + geometry).
    pub pumping: AcSpinPumping,
    /// The ISHE detector (spin Hall angle, resistivity of the NM).
    pub detector: InverseSpinHall,
    /// Width of the NM strip (perpendicular to J_s, along σ̂) \[m\].
    pub sample_width: f64,
    /// Length of the NM strip along the ISHE voltage direction \[m\].
    pub sample_length: f64,
}

impl SpinBattery {
    /// Construct a `SpinBattery` with explicit parameters.
    ///
    /// # Errors
    /// Returns `Err` if `width` or `length` are non-positive.
    pub fn new(
        pumping: AcSpinPumping,
        detector: InverseSpinHall,
        width: f64,
        length: f64,
    ) -> Result<Self> {
        if width <= 0.0 {
            return Err(invalid_param("width", "sample width must be positive"));
        }
        if length <= 0.0 {
            return Err(invalid_param("length", "sample length must be positive"));
        }
        Ok(Self {
            pumping,
            detector,
            sample_width: width,
            sample_length: length,
        })
    }

    /// Canonical YIG/Pt spin battery preset.
    ///
    /// Dimensions representative of Saitoh-group thin-film devices:
    /// - NM strip width  = 3 mm
    /// - NM strip length = 10 mm
    pub fn yig_pt_battery() -> Self {
        Self {
            pumping: AcSpinPumping::yig_pt(),
            detector: InverseSpinHall::platinum(),
            sample_width: 3.0e-3,
            sample_length: 10.0e-3,
        }
    }

    // -----------------------------------------------------------------------
    // ISHE voltages
    // -----------------------------------------------------------------------

    /// Open-circuit ISHE voltage at FMR \[V\].
    ///
    /// Uses the time-averaged DC spin current at the given cone angle and
    /// resonance frequency, then applies the ISHE conversion:
    ///
    /// ```text
    /// V_ISHE = θ_SH · ρ_NM · |J_s^DC| · L
    /// ```
    ///
    /// where `L = sample_length` is the distance between the voltage contacts
    /// along the ISHE electric field direction.
    ///
    /// # Arguments
    /// * `omega_fmr`  — FMR angular frequency \[rad s⁻¹\].
    /// * `cone_angle` — precession cone half-angle θ \[rad\].
    pub fn ishe_voltage(&self, omega_fmr: f64, cone_angle: f64) -> f64 {
        let j_s = self.pumping.dc_magnitude_at_fmr(omega_fmr, cone_angle);
        self.detector.theta_sh * self.detector.rho * j_s * self.sample_length
    }

    /// Open-circuit ISHE voltage from the instantaneous magnetisation state.
    ///
    /// Computes the spin current directly from `m` and `dm/dt` via
    /// [`AcSpinPumping::dc_spin_current`], takes its magnitude, and converts
    /// via ISHE.  This form is useful when integrating the full LLG trajectory.
    ///
    /// # Arguments
    /// * `m`    — unit magnetisation vector.
    /// * `dmdt` — time derivative of `m` \[s⁻¹\].
    pub fn open_circuit_voltage(&self, m: Vector3<f64>, dmdt: Vector3<f64>) -> f64 {
        let j_s_vec = self.pumping.dc_spin_current(m, dmdt);
        let j_s_mag = j_s_vec.magnitude();
        self.detector.theta_sh * self.detector.rho * j_s_mag * self.sample_length
    }

    // -----------------------------------------------------------------------
    // Signal quality
    // -----------------------------------------------------------------------

    /// Estimate the signal-to-noise ratio at the ISHE output.
    ///
    /// Compares the FMR spin-pumping signal `V_ISHE` against the Johnson–
    /// Nyquist thermal noise voltage of the NM strip at temperature `T`.
    ///
    /// The thermal noise in a bandwidth `B = 1 Hz` is
    ///
    /// ```text
    /// V_noise = sqrt(4 · k_B · T · R_NM · B)
    /// ```
    ///
    /// where the NM strip resistance is estimated from the detector resistivity,
    /// strip length, and an assumed strip thickness of 10 nm:
    ///
    /// ```text
    /// R_NM = ρ_NM · L / (W · d_NM)       d_NM = 10 nm
    /// ```
    ///
    /// # Arguments
    /// * `omega_fmr`   — FMR angular frequency \[rad s⁻¹\].
    /// * `cone_angle`  — precession cone half-angle θ \[rad\].
    /// * `temperature` — measurement temperature \[K\].
    ///
    /// # Returns
    /// Dimensionless SNR = V_signal / V_noise (larger is better).
    pub fn signal_to_noise_estimate(
        &self,
        omega_fmr: f64,
        cone_angle: f64,
        temperature: f64,
    ) -> f64 {
        let v_signal = self.ishe_voltage(omega_fmr, cone_angle);

        // Assumed NM strip thickness for resistance estimate
        let d_nm = 10.0e-9_f64; // 10 nm
        let r_nm = self.detector.rho * self.sample_length / (self.sample_width * d_nm);

        // Johnson–Nyquist noise voltage in 1 Hz bandwidth
        let bandwidth = 1.0_f64; // Hz
        let v_noise = (4.0 * KB * temperature * r_nm * bandwidth).sqrt();

        // Guard against division by zero (e.g. T = 0 K)
        if v_noise == 0.0 {
            return f64::INFINITY;
        }

        v_signal.abs() / v_noise
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn default_pumping() -> AcSpinPumping {
        AcSpinPumping::yig_pt()
    }

    fn default_battery() -> SpinBattery {
        SpinBattery::yig_pt_battery()
    }

    /// DC spin current must be collinear with m × dm/dt and therefore
    /// orthogonal to both m and dm/dt individually (for |m| = 1 and
    /// dm/dt ⊥ m, as in uniform precession).
    #[test]
    fn test_dc_current_direction() {
        let pumping = default_pumping();

        // Uniform circular precession in the x-z plane at angle PI/6
        let theta = PI / 6.0;
        let m = Vector3::new(theta.sin(), 0.0, theta.cos());
        // dm/dt is tangent to the precession circle — for ccw precession about z:
        // dm/dt = ω (ẑ × m) = ω (-m_y ẑ + m_z ŷ - ...) but for m in xz-plane:
        // dm/dt ∝ ŷ × m rotated → simplest: dm/dt = ω (m × ẑ) (LLG precession sense)
        let omega = 2.0 * PI * 10.0e9; // 10 GHz
        let dmdt = m.cross(&Vector3::new(0.0, 0.0, 1.0)) * omega; // ω (m × ẑ)

        let j_s = pumping.dc_spin_current(m, dmdt);

        // J_s is along m × dm/dt — orthogonal to both m and dm/dt
        let dot_m = j_s.dot(&m).abs();
        let dot_dmdt = j_s.dot(&dmdt).abs();

        // Allow small floating-point residuals relative to the vector magnitudes
        let j_mag = j_s.magnitude();
        let m_mag = m.magnitude();
        let d_mag = dmdt.magnitude();

        assert!(
            dot_m < 1.0e-6 * j_mag * m_mag,
            "J_s not orthogonal to m: dot = {dot_m:.3e}"
        );
        assert!(
            dot_dmdt < 1.0e-6 * j_mag * d_mag,
            "J_s not orthogonal to dm/dt: dot = {dot_dmdt:.3e}"
        );
    }

    /// Backflow always reduces the effective mixing conductance.
    #[test]
    fn test_backflow_reduces_g_r() {
        let pumping = default_pumping();
        let g_r_bare = pumping.interface.g_r;
        let g_r_eff = pumping.effective_g_r();

        assert!(
            g_r_eff <= g_r_bare,
            "G_r_eff ({g_r_eff:.3e}) should be ≤ G_r ({g_r_bare:.3e})"
        );
        assert!(g_r_eff > 0.0, "G_r_eff must be positive");
    }

    /// Spin-pumping damping enhancement is always positive for physical parameters.
    #[test]
    fn test_damping_enhancement_positive() {
        let pumping = default_pumping();
        let yig = Ferromagnet::yig();
        let delta_alpha = pumping.damping_enhancement(&yig);

        assert!(
            delta_alpha > 0.0,
            "Δα must be positive, got {delta_alpha:.3e}"
        );
    }

    /// ISHE voltage is positive (non-zero) at FMR for a standard cone angle.
    #[test]
    fn test_ishe_voltage_positive_at_fmr() {
        let battery = default_battery();
        let omega_fmr = 10.0e9 * 2.0 * PI; // 10 GHz
        let cone_angle = PI / 4.0;

        let v = battery.ishe_voltage(omega_fmr, cone_angle);

        assert!(
            v > 0.0,
            "V_ISHE must be positive for Pt (positive θ_SH), got {v:.3e} V"
        );
        assert!(v.is_finite(), "V_ISHE must be finite, got {v}");
    }

    /// G_NM must be positive and finite for physical parameters.
    #[test]
    fn test_g_nm_formula() {
        let pumping = default_pumping();
        let g_nm = pumping.g_nm();

        assert!(g_nm > 0.0, "G_NM must be positive, got {g_nm:.3e}");
        assert!(g_nm.is_finite(), "G_NM must be finite, got {g_nm}");

        // Sanity: G_NM should be in a physically reasonable range
        // For Pt: σ ≈ 9.43e6 S/m, λ_sf ≈ 10 nm, t_NM = 7 nm
        // tanh(0.7) ≈ 0.604, G_NM ≈ 9.43e6 * 0.604 / 10e-9 ≈ 5.7e14 m⁻²
        assert!(
            g_nm > 1.0e12 && g_nm < 1.0e18,
            "G_NM = {g_nm:.3e} outside plausible range [1e12, 1e18]"
        );
    }

    /// At FMR, the 2ω AC component equals the DC component exactly.
    #[test]
    fn test_ac_equals_dc_at_fmr() {
        let pumping = default_pumping();
        let omega_fmr = 2.0 * PI * 10.0e9;
        let cone_angle = PI / 8.0;

        let j_dc = pumping.dc_magnitude_at_fmr(omega_fmr, cone_angle);
        let j_2w = pumping.ac_component_2omega(omega_fmr, cone_angle);

        assert!(
            (j_dc - j_2w).abs() < 1.0e-30,
            "AC 2ω component ({j_2w:.3e}) must equal DC ({j_dc:.3e}) at resonance"
        );
    }

    /// `new()` rejects non-positive thicknesses.
    #[test]
    fn test_new_validates_thicknesses() {
        let iface = SpinInterface::yig_pt();
        let diff = SpinDiffusion::platinum();

        assert!(
            AcSpinPumping::new(iface.clone(), diff.clone(), -1.0e-9, 5.0e-6).is_err(),
            "Should reject t_nm <= 0"
        );
        assert!(
            AcSpinPumping::new(iface.clone(), diff.clone(), 7.0e-9, 0.0).is_err(),
            "Should reject t_fm <= 0"
        );
        assert!(
            AcSpinPumping::new(iface, diff, 7.0e-9, 5.0e-6).is_ok(),
            "Should accept valid thicknesses"
        );
    }

    /// SpinBattery::new() rejects non-positive dimensions.
    #[test]
    fn test_battery_new_validates_dimensions() {
        let pumping = AcSpinPumping::yig_pt();
        let detector = InverseSpinHall::platinum();

        assert!(
            SpinBattery::new(pumping.clone(), detector.clone(), -3.0e-3, 10.0e-3).is_err(),
            "Should reject non-positive width"
        );
        assert!(
            SpinBattery::new(pumping.clone(), detector.clone(), 3.0e-3, 0.0).is_err(),
            "Should reject non-positive length"
        );
        assert!(
            SpinBattery::new(pumping, detector, 3.0e-3, 10.0e-3).is_ok(),
            "Should accept valid dimensions"
        );
    }

    /// Spin accumulation is positive and finite for physical parameters.
    #[test]
    fn test_spin_accumulation_positive() {
        let pumping = default_pumping();
        let omega_fmr = 2.0 * PI * 10.0e9;
        let cone_angle = PI / 4.0;

        let mu_s = pumping.spin_accumulation_dc(omega_fmr, cone_angle);

        assert!(
            mu_s > 0.0,
            "Spin accumulation must be positive, got {mu_s:.3e}"
        );
        assert!(mu_s.is_finite(), "Spin accumulation must be finite");
    }

    /// SNR is positive and finite for room temperature.
    #[test]
    fn test_snr_positive_at_room_temperature() {
        let battery = default_battery();
        let omega_fmr = 2.0 * PI * 10.0e9;
        let cone_angle = PI / 4.0;
        let temperature = 300.0;

        let snr = battery.signal_to_noise_estimate(omega_fmr, cone_angle, temperature);

        assert!(snr > 0.0, "SNR must be positive, got {snr:.3e}");
        assert!(snr.is_finite(), "SNR must be finite");
    }

    /// The DC spin current magnitude matches the expected analytical formula.
    #[test]
    fn test_dc_magnitude_formula() {
        let pumping = default_pumping();
        let omega = 2.0 * PI * 10.0e9;
        let theta = PI / 4.0;

        let j_dc = pumping.dc_magnitude_at_fmr(omega, theta);
        let g_eff = pumping.effective_g_r();
        let expected = (HBAR / (8.0 * PI)) * g_eff * omega * theta.sin().powi(2);

        assert!(
            (j_dc - expected).abs() < 1.0e-40,
            "DC magnitude {j_dc:.6e} ≠ expected {expected:.6e}"
        );
    }
}
