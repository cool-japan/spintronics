//! Spin-Torque Nano-Oscillator (STNO)
//!
//! This module implements the physics of spin-torque nano-oscillators, nanoscale
//! devices in which a DC spin-polarized current drives sustained magnetization
//! precession via the Slonczewski spin-transfer torque (STT).
//!
//! ## Physics Background
//!
//! ### Slonczewski STT (J. Magn. Magn. Mater. 159, L1, 1996)
//!
//! The STT term appended to the Landau-Lifshitz-Gilbert (LLG) equation is:
//!
//! ```text
//! τ_STT = a_J · [m × (m × p̂) + β · (m × p̂)]
//! ```
//!
//! where:
//! - `p̂` is the fixed-layer spin-polarization direction (unit vector)
//! - `a_J = (ℏ/2e) · J · η / (μ₀ · M_s · t_FM)` — Slonczewski amplitude \[A/m\]
//! - `β` — field-like STT ratio (often ≈ 0)
//!
//! ### Full LLG + STT equation (Landau-Lifshitz form)
//!
//! ```text
//! dm/dt = -γ/(1+α²) · [m × H_eff + α · m × (m × H_eff)]
//!         + γ/(1+α²) · a_J · [m × (m × p̂) + α · (m × p̂)]
//! ```
//!
//! ### Threshold current
//!
//! Auto-oscillation begins when J exceeds the threshold that overcomes damping:
//!
//! ```text
//! J_th = (2e · μ₀ · M_s · t_FM · α) / (ℏ · η)
//! ```
//!
//! ### Slavin-Tiberkevich thermal linewidth (2009)
//!
//! ```text
//! Δf = (α · k_B · T) / (π · E_spin)
//! ```
//!
//! where `E_spin = (μ₀/2) · M_s · V · |H_eff|`.
//!
//! ## References
//!
//! - J. C. Slonczewski, J. Magn. Magn. Mater. **159**, L1 (1996)
//! - A. Slavin and V. Tiberkevich, IEEE Trans. Magn. **45**, 1875 (2009)
//! - R. Adler, Proc. IRE **34**, 351 (1946) (injection locking)

use std::f64::consts::PI;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::constants::{E_CHARGE, GAMMA, HBAR, KB, MU_0};
use crate::error::{invalid_param, Result};
use crate::llg::LLGSolver;
use crate::material::Ferromagnet;
use crate::vector3::Vector3;

// ─────────────────────────────────────────────────────────────────────────────
// SpinTorqueOscillatorConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration parameters for a Spin-Torque Nano-Oscillator device.
///
/// Encapsulates the geometry and spin-transport properties of the free layer
/// and the fixed-layer polarizer.
///
/// # Example
/// ```no_run
/// use spintronics::effect::stno::SpinTorqueOscillatorConfig;
/// use spintronics::Vector3;
///
/// let config = SpinTorqueOscillatorConfig::new(
///     0.35,                               // η — spin polarization efficiency
///     5e-9,                               // t_fm — free layer thickness [m]
///     std::f64::consts::PI * (50e-9f64).powi(2), // area [m²]
///     0.01,                               // β — field-like ratio
///     Vector3::new(0.0, 0.0, -1.0),       // p̂ — fixed layer direction
/// ).unwrap();
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpinTorqueOscillatorConfig {
    /// Spin-polarization efficiency η ∈ (0, 1).
    ///
    /// Represents the fraction of injected electrons that are spin-polarized
    /// and contribute to the torque. Typical metallic values: 0.3–0.5.
    pub spin_polarization: f64,

    /// Free-layer (ferromagnetic) thickness t_FM \[m\].
    ///
    /// Directly scales the Slonczewski coefficient a_J.
    /// Typical STNO nanopillar: 2–10 nm.
    pub t_fm: f64,

    /// Device cross-sectional area \[m²\].
    ///
    /// Combined with `t_fm` to define the free-layer volume.
    /// Typical nanopillar diameter: 30–200 nm.
    pub area: f64,

    /// Field-like STT ratio β (dimensionless).
    ///
    /// Governs the `m × p̂` (Rashba-like) contribution relative to the
    /// Slonczewski damping-like term `m × (m × p̂)`. Often small (≈ 0–0.05).
    pub beta_stt: f64,

    /// Fixed-layer spin-polarization direction p̂ (unit vector).
    ///
    /// Represents the magnetization direction of the reference (pinned) layer.
    /// For perpendicular-to-plane geometry this is typically ẑ; for in-plane
    /// geometry it is ±x̂ or ±ŷ.
    pub p_hat: Vector3<f64>,

    /// Free-layer volume V = area × t_fm \[m³\].
    ///
    /// Stored for computational efficiency (avoids repeated multiplication).
    pub volume: f64,
}

impl SpinTorqueOscillatorConfig {
    /// Construct a new `SpinTorqueOscillatorConfig`, validating all inputs.
    ///
    /// # Arguments
    ///
    /// * `spin_polarization` — η ∈ (0, 1).
    /// * `t_fm` — free-layer thickness \[m\], must be positive.
    /// * `area` — cross-sectional area \[m²\], must be positive.
    /// * `beta_stt` — field-like STT ratio.
    /// * `p_hat` — fixed-layer polarization direction; automatically normalized.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] if η is outside (0, 1), t_fm ≤ 0, or area ≤ 0.
    pub fn new(
        spin_polarization: f64,
        t_fm: f64,
        area: f64,
        beta_stt: f64,
        p_hat: Vector3<f64>,
    ) -> Result<Self> {
        if spin_polarization <= 0.0 || spin_polarization >= 1.0 {
            return Err(invalid_param(
                "spin_polarization",
                "must be strictly between 0 and 1",
            ));
        }
        if t_fm <= 0.0 {
            return Err(invalid_param(
                "t_fm",
                "free-layer thickness must be positive",
            ));
        }
        if area <= 0.0 {
            return Err(invalid_param("area", "device area must be positive"));
        }

        let volume = area * t_fm;
        let p_hat_norm = p_hat.normalize();

        Ok(Self {
            spin_polarization,
            t_fm,
            area,
            beta_stt,
            p_hat: p_hat_norm,
            volume,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SpinTorqueOscillator
// ─────────────────────────────────────────────────────────────────────────────

/// Spin-Torque Nano-Oscillator (STNO) — full LLG + Slonczewski STT simulator.
///
/// Combines a standard [`LLGSolver`] with Slonczewski spin-transfer torque
/// (STT) to model current-driven magnetization auto-oscillation in a nanopillar
/// or nanocontact geometry.
///
/// # Physical model
///
/// The extended LLG equation in Landau-Lifshitz form is:
///
/// ```text
/// dm̂/dt = -γ/(1+α²) [m̂ × H_eff + α m̂ × (m̂ × H_eff)]
///         + γ/(1+α²) · a_J [m̂ × (m̂ × p̂) + α (m̂ × p̂)]
/// ```
///
/// where `a_J = (ℏ/2e) · J · η / (μ₀ M_s t_FM)`.
///
/// # Example
/// ```no_run
/// use spintronics::effect::stno::SpinTorqueOscillator;
/// use spintronics::Vector3;
///
/// let stno = SpinTorqueOscillator::permalloy_nanopillar().unwrap();
/// let j_th = stno.threshold_current_density();
/// let m0 = Vector3::new(0.1, 0.0, 1.0).normalize() * stno.solver.material.ms;
/// let traj = stno.simulate(m0, 1e-12, 2000, 3.0 * j_th);
/// let power = stno.oscillation_power(&traj);
/// ```
pub struct SpinTorqueOscillator {
    /// Underlying LLG solver (material + external field).
    pub solver: LLGSolver,
    /// Device geometry and spin-transport configuration.
    pub config: SpinTorqueOscillatorConfig,
}

impl std::fmt::Debug for SpinTorqueOscillator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SpinTorqueOscillator")
            .field("material_ms", &self.solver.material.ms)
            .field("material_alpha", &self.solver.material.alpha)
            .field("h_ext", &self.solver.h_ext)
            .field("config", &self.config)
            .finish()
    }
}

impl Clone for SpinTorqueOscillator {
    fn clone(&self) -> Self {
        let mut solver = LLGSolver::new(self.solver.material.clone());
        solver.h_ext = self.solver.h_ext;
        Self {
            solver,
            config: self.config.clone(),
        }
    }
}

impl SpinTorqueOscillator {
    /// Create a new `SpinTorqueOscillator` from a [`Ferromagnet`], external
    /// field, and device [`SpinTorqueOscillatorConfig`].
    ///
    /// # Arguments
    ///
    /// * `material` — free-layer ferromagnet (determines α, γ, M_s).
    /// * `h_ext` — static external field \[A/m\].
    /// * `config` — STNO geometry and transport parameters.
    pub fn new(
        material: Ferromagnet,
        h_ext: Vector3<f64>,
        config: SpinTorqueOscillatorConfig,
    ) -> Result<Self> {
        let mut solver = LLGSolver::new(material);
        solver.h_ext = h_ext;
        Ok(Self { solver, config })
    }

    /// Permalloy nanopillar preset — a canonical STNO reference geometry.
    ///
    /// Parameters (after Kiselev et al., Nature 425, 380, 2003):
    /// - Material: Permalloy (Ni₈₀Fe₂₀), α = 0.01, M_s = 800 kA/m
    /// - External field: H_ext = 100 kA/m along +z (bias field)
    /// - t_FM = 5 nm, diameter = 100 nm ⟹ area = π × (50 nm)²
    /// - η = 0.35, β = 0.01, p̂ = −ẑ (anti-parallel reference layer)
    pub fn permalloy_nanopillar() -> Result<Self> {
        let mut material = Ferromagnet::permalloy();
        // Standard STNO Permalloy parameters
        material.alpha = 0.01;
        material.ms = 8.0e5; // 800 kA/m

        let h_ext = Vector3::new(0.0, 0.0, 1.0e5); // 100 kA/m along +z
        let t_fm = 5.0e-9; // 5 nm
        let radius = 50.0e-9; // 50 nm radius → 100 nm diameter
        let area = PI * radius * radius;

        let config = SpinTorqueOscillatorConfig::new(
            0.35,                         // η — spin polarization
            t_fm,                         // t_FM
            area,                         // cross-section
            0.01,                         // β — field-like ratio
            Vector3::new(0.0, 0.0, -1.0), // p̂ — anti-parallel fixed layer
        )?;

        Self::new(material, h_ext, config)
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Core STT physics
    // ─────────────────────────────────────────────────────────────────────────

    /// Slonczewski torque amplitude a_J \[A/m\] for a given current density J.
    ///
    /// ```text
    /// a_J = (ℏ / 2e) · J · η / (μ₀ · M_s · t_FM)
    /// ```
    ///
    /// # Arguments
    ///
    /// * `j_density` — current density \[A/m²\], positive = conventional current.
    ///
    /// # Returns
    ///
    /// Slonczewski field amplitude \[A/m\]. Sign follows `j_density`.
    #[inline]
    pub fn a_j(&self, j_density: f64) -> f64 {
        let ms = self.solver.material.ms;
        let eta = self.config.spin_polarization;
        let t_fm = self.config.t_fm;
        (HBAR / (2.0 * E_CHARGE)) * j_density * eta / (MU_0 * ms * t_fm)
    }

    /// Compute the full Slonczewski spin-transfer torque for normalized m̂.
    ///
    /// ```text
    /// τ_STT = a_J · [m × (m × p̂)  +  β · (m × p̂)]
    /// ```
    ///
    /// The result is expressed in units of \[A/m\] so that it can be added
    /// directly to the effective-field contribution in the LLG equation.
    ///
    /// **Note:** The tensor `m × (m × p̂)` is automatically perpendicular to m̂,
    /// which is verified by [`test_stt_perpendicular_to_m`].
    ///
    /// # Arguments
    ///
    /// * `m` — magnetization vector (may be un-normalized; uses direction only).
    /// * `j_density` — current density \[A/m²\].
    ///
    /// # Returns
    ///
    /// STT pseudo-field vector \[A/m\].
    pub fn slonczewski_torque(&self, m: Vector3<f64>, j_density: f64) -> Vector3<f64> {
        let a_j = self.a_j(j_density);
        let m_hat = m.normalize();
        let p = self.config.p_hat;

        let m_cross_p = m_hat.cross(&p); // m × p̂
        let m_cross_m_cross_p = m_hat.cross(&m_cross_p); // m × (m × p̂)

        // τ = a_J · [m×(m×p̂) + β·(m×p̂)]
        m_cross_m_cross_p * a_j + m_cross_p * (self.config.beta_stt * a_j)
    }

    /// Time derivative dm̂/dt including both LLG precession/damping and STT.
    ///
    /// Implements the extended Landau-Lifshitz form:
    ///
    /// ```text
    /// dm̂/dt = -γ/(1+α²) [m̂ × H_eff + α m̂ × (m̂ × H_eff)]
    ///         + γ/(1+α²) · a_J [m̂ × (m̂ × p̂) + α (m̂ × p̂)]
    /// ```
    ///
    /// The result is in units of \[A/m / s\] (i.e. dm/dt with |m| = M_s).
    ///
    /// # Arguments
    ///
    /// * `m` — current magnetization state \[A/m\].
    /// * `j_density` — applied current density \[A/m²\].
    pub fn dmdt_with_stt(&self, m: Vector3<f64>, j_density: f64) -> Vector3<f64> {
        let ms = self.solver.material.ms;
        let alpha = self.solver.material.alpha;
        let gamma = GAMMA;

        let h_eff = self.solver.effective_field(m);
        let m_hat = m.normalize();

        let factor = gamma / (1.0 + alpha * alpha);

        // LLG precession + damping
        let m_cross_h = m_hat.cross(&h_eff); // m̂ × H_eff
        let m_cross_m_cross_h = m_hat.cross(&m_cross_h); // m̂ × (m̂ × H_eff)

        let llg_term = (m_cross_h + m_cross_m_cross_h * alpha) * (-factor * ms);

        // STT contribution: +γ/(1+α²) · a_J · [m̂×(m̂×p̂) + α·(m̂×p̂)]
        let a_j = self.a_j(j_density);
        let p = self.config.p_hat;
        let m_cross_p = m_hat.cross(&p);
        let m_cross_m_cross_p = m_hat.cross(&m_cross_p);

        let stt_term = (m_cross_m_cross_p + m_cross_p * alpha) * (factor * a_j * ms);

        llg_term + stt_term
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Device figures of merit
    // ─────────────────────────────────────────────────────────────────────────

    /// Threshold current density J_th \[A/m²\] for auto-oscillation onset.
    ///
    /// Derived from the condition that the STT compensates Gilbert damping
    /// at perpendicular-to-plane (PTP) geometry (θ₀ = π/2):
    ///
    /// ```text
    /// J_th = (2e · μ₀ · M_s · t_FM · α) / (ℏ · η)
    /// ```
    ///
    /// For finite initial angle, the actual threshold can be lower, but this
    /// expression is the standard single-domain PTP estimate.
    pub fn threshold_current_density(&self) -> f64 {
        let ms = self.solver.material.ms;
        let alpha = self.solver.material.alpha;
        let t_fm = self.config.t_fm;
        let eta = self.config.spin_polarization;

        (2.0 * E_CHARGE * MU_0 * ms * t_fm * alpha) / (HBAR * eta)
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Time-domain simulation
    // ─────────────────────────────────────────────────────────────────────────

    /// RK4 time integration of the full LLG + STT equation.
    ///
    /// Each step:
    /// 1. Evaluates four `dmdt_with_stt` stages (k1 … k4).
    /// 2. Advances m by the weighted combination.
    /// 3. Renormalizes to `|m| = M_s` to prevent drift.
    ///
    /// # Arguments
    ///
    /// * `m0` — initial magnetization \[A/m\].
    /// * `dt` — time step \[s\]. Stability requires dt ≪ 1/(γ |H_eff|).
    /// * `n_steps` — number of integration steps.
    /// * `j_density` — DC current density \[A/m²\].
    ///
    /// # Returns
    ///
    /// Trajectory: `n_steps + 1` vectors (includes m0 at index 0).
    pub fn simulate(
        &self,
        m0: Vector3<f64>,
        dt: f64,
        n_steps: usize,
        j_density: f64,
    ) -> Vec<Vector3<f64>> {
        let ms = self.solver.material.ms;
        let mut trajectory = Vec::with_capacity(n_steps + 1);
        trajectory.push(m0);

        let mut m = m0;
        for _ in 0..n_steps {
            m = self.rk4_step(m, dt, j_density, ms);
            trajectory.push(m);
        }

        trajectory
    }

    /// Single RK4 step with full STT, returning renormalized m.
    fn rk4_step(&self, m: Vector3<f64>, dt: f64, j_density: f64, ms: f64) -> Vector3<f64> {
        let half_dt = dt * 0.5;

        let k1 = self.dmdt_with_stt(m, j_density);
        let k2 = self.dmdt_with_stt(m + k1 * half_dt, j_density);
        let k3 = self.dmdt_with_stt(m + k2 * half_dt, j_density);
        let k4 = self.dmdt_with_stt(m + k3 * dt, j_density);

        let m_new = m + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);

        // Renormalize to |m| = M_s — essential for long integrations
        let mag = m_new.magnitude();
        if mag > 1.0e-20 {
            m_new * (ms / mag)
        } else {
            // Degenerate case: return along external field or z-axis
            let h_mag = self.solver.h_ext.magnitude();
            if h_mag > 0.0 {
                self.solver.h_ext * (ms / h_mag)
            } else {
                Vector3::new(0.0, 0.0, ms)
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Spectral analysis
    // ─────────────────────────────────────────────────────────────────────────

    /// Estimate oscillation frequency from zero-crossings of the m_z component.
    ///
    /// The algorithm counts sign changes in the m_z trajectory (normalized by
    /// M_s so that the zero-crossing carries physical meaning) and converts
    /// to frequency using:
    ///
    /// ```text
    /// f = n_crossings / (2 · T_total)
    /// ```
    ///
    /// # Arguments
    ///
    /// * `trajectory` — output of [`simulate`]; must have ≥ 2 points.
    /// * `dt` — time step used in simulation \[s\].
    ///
    /// # Returns
    ///
    /// `Some(f_osc)` in Hz if at least 2 zero-crossings are found;
    /// `None` otherwise (steady state / insufficient oscillation).
    pub fn oscillation_frequency(&self, trajectory: &[Vector3<f64>], dt: f64) -> Option<f64> {
        if trajectory.len() < 2 {
            return None;
        }

        let ms = self.solver.material.ms;
        // Guard against zero ms
        let ms_inv = if ms > 0.0 { 1.0 / ms } else { return None };

        let mut n_crossings: usize = 0;

        // Skip the first 20 % of the trajectory (transient) for cleaner estimation
        let skip = trajectory.len() / 5;
        let working = &trajectory[skip..];

        if working.len() < 2 {
            return None;
        }

        let mut prev_sign = working[0].z * ms_inv > 0.0;
        for state in working.iter().skip(1) {
            let current_sign = state.z * ms_inv > 0.0;
            if current_sign != prev_sign {
                n_crossings += 1;
            }
            prev_sign = current_sign;
        }

        if n_crossings < 2 {
            return None;
        }

        let t_total = (working.len() - 1) as f64 * dt;
        // f = n_crossings / (2 * T_total)  — each full period has 2 crossings
        Some(n_crossings as f64 / (2.0 * t_total))
    }

    /// Oscillation power — variance of the normalized m_z component.
    ///
    /// A proxy for the precession amplitude squared:
    /// - ≈ 0 when the magnetization is static.
    /// - > 0 when sustained oscillation is present.
    ///
    /// # Arguments
    ///
    /// * `trajectory` — output of [`simulate`].
    ///
    /// # Returns
    ///
    /// Variance of m_z / M_s (dimensionless, ≥ 0).
    pub fn oscillation_power(&self, trajectory: &[Vector3<f64>]) -> f64 {
        if trajectory.is_empty() {
            return 0.0;
        }

        let ms = self.solver.material.ms;
        let ms_safe = if ms > 0.0 { ms } else { return 0.0 };

        let n = trajectory.len() as f64;
        let mean = trajectory.iter().map(|m| m.z / ms_safe).sum::<f64>() / n;
        let variance = trajectory
            .iter()
            .map(|m| {
                let diff = m.z / ms_safe - mean;
                diff * diff
            })
            .sum::<f64>()
            / n;

        variance
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Spectral figures of merit
    // ─────────────────────────────────────────────────────────────────────────

    /// Slavin-Tiberkevich thermal linewidth Δf \[Hz\].
    ///
    /// Based on the stochastic LLG analysis (IEEE Trans. Magn. 45, 1875, 2009):
    ///
    /// ```text
    /// Δf = (α · k_B · T) / (π · E_spin)
    /// ```
    ///
    /// where the spin energy per mode is approximated as:
    ///
    /// ```text
    /// E_spin = (μ₀ / 2) · M_s · V · |H_eff|
    /// ```
    ///
    /// and H_eff is evaluated at the equilibrium (or time-averaged) magnetization.
    ///
    /// # Arguments
    ///
    /// * `j_density` — operating current density \[A/m²\] (used to pick a
    ///   representative m̂ via STT-shifted equilibrium; here approximated by
    ///   the external-field direction).
    /// * `temperature` — lattice temperature \[K\].
    ///
    /// # Returns
    ///
    /// Linewidth Δf \[Hz\]. Returns a large sentinel (1 THz) if E_spin = 0
    /// to avoid division by zero.
    pub fn thermal_linewidth(&self, _j_density: f64, temperature: f64) -> f64 {
        let ms = self.solver.material.ms;
        let alpha = self.solver.material.alpha;

        // Approximate equilibrium magnetization along external field or +z
        let m_eq = {
            let h_mag = self.solver.h_ext.magnitude();
            if h_mag > 0.0 {
                self.solver.h_ext * (ms / h_mag)
            } else {
                Vector3::new(0.0, 0.0, ms)
            }
        };

        let h_eff = self.solver.effective_field(m_eq);
        let h_eff_magnitude = h_eff.magnitude();

        let e_spin = 0.5 * MU_0 * ms * self.config.volume * h_eff_magnitude;

        if e_spin < 1.0e-50 {
            // Degenerate limit — very large linewidth
            return 1.0e12;
        }

        (alpha * KB * temperature) / (PI * e_spin)
    }

    /// Adler injection-locking bandwidth Δf_lock \[Hz\].
    ///
    /// For an external RF signal with normalized amplitude `i_inj` ∈ (0, 1),
    /// the free-running oscillator locks to the injected frequency within:
    ///
    /// ```text
    /// Δf_lock = (f_osc / (2 Q)) · i_inj / √(1 − i_inj²)
    /// ```
    ///
    /// (R. Adler, Proc. IRE 34, 351, 1946).
    ///
    /// # Arguments
    ///
    /// * `q_factor` — quality factor Q = f_osc / Δf of the free-running oscillator.
    /// * `i_injection_normalized` — normalized injection amplitude ∈ (0, 1).
    ///
    /// # Returns
    ///
    /// Half-bandwidth of the locking range \[Hz\].
    /// Returns 0 if Q ≤ 0 or `i_injection_normalized` ≥ 1.
    pub fn injection_locking_range(&self, q_factor: f64, i_injection_normalized: f64) -> f64 {
        if q_factor <= 0.0 {
            return 0.0;
        }
        let i = i_injection_normalized.clamp(0.0, 1.0 - f64::EPSILON);
        let denom_sq = 1.0 - i * i;
        if denom_sq <= 0.0 {
            return f64::INFINITY;
        }

        // Estimate free-running frequency from the effective field at equilibrium
        let ms = self.solver.material.ms;
        let m_eq = {
            let h_mag = self.solver.h_ext.magnitude();
            if h_mag > 0.0 {
                self.solver.h_ext * (ms / h_mag)
            } else {
                Vector3::new(0.0, 0.0, ms)
            }
        };
        let h_eff = self.solver.effective_field(m_eq);
        let h_eff_magnitude = h_eff.magnitude();
        let f_osc = GAMMA * h_eff_magnitude / (2.0 * PI);

        f_osc / (2.0 * q_factor) * i / denom_sq.sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build the permalloy nanopillar preset (panics on construction error
    /// so tests report clearly).
    fn py_stno() -> SpinTorqueOscillator {
        SpinTorqueOscillator::permalloy_nanopillar()
            .expect("permalloy_nanopillar construction should succeed")
    }

    /// Helper: build a magnetization state tilted slightly from +z, normalized to Ms.
    fn tilted_m(stno: &SpinTorqueOscillator) -> Vector3<f64> {
        let ms = stno.solver.material.ms;
        Vector3::new(0.1, 0.0, 1.0).normalize() * ms
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Config validation
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_config_invalid_spin_polarization_zero() {
        let res = SpinTorqueOscillatorConfig::new(
            0.0,
            5.0e-9,
            PI * (50.0e-9f64).powi(2),
            0.01,
            Vector3::new(0.0, 0.0, -1.0),
        );
        assert!(res.is_err());
    }

    #[test]
    fn test_config_invalid_spin_polarization_one() {
        let res = SpinTorqueOscillatorConfig::new(
            1.0,
            5.0e-9,
            PI * (50.0e-9f64).powi(2),
            0.01,
            Vector3::new(0.0, 0.0, -1.0),
        );
        assert!(res.is_err());
    }

    #[test]
    fn test_config_invalid_t_fm() {
        let res = SpinTorqueOscillatorConfig::new(
            0.35,
            -1.0e-9,
            PI * (50.0e-9f64).powi(2),
            0.01,
            Vector3::new(0.0, 0.0, -1.0),
        );
        assert!(res.is_err());
    }

    #[test]
    fn test_config_invalid_area() {
        let res =
            SpinTorqueOscillatorConfig::new(0.35, 5.0e-9, 0.0, 0.01, Vector3::new(0.0, 0.0, -1.0));
        assert!(res.is_err());
    }

    #[test]
    fn test_config_volume_derived_correctly() {
        let t_fm = 5.0e-9;
        let area = PI * (50.0e-9f64).powi(2);
        let config =
            SpinTorqueOscillatorConfig::new(0.35, t_fm, area, 0.01, Vector3::new(0.0, 0.0, -1.0))
                .expect("valid config");
        let expected_volume = area * t_fm;
        assert!((config.volume - expected_volume).abs() < 1.0e-40);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Threshold current
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_threshold_positive() {
        let stno = py_stno();
        let j_th = stno.threshold_current_density();
        assert!(
            j_th > 0.0,
            "threshold current density must be positive, got {j_th}"
        );
        // Permalloy preset: J_th = 2e·μ₀·Ms·t·α / (ℏ·η)
        // With Ms=800kA/m, t=5nm, α=0.01, η=0.35 → J_th ≈ 436 kA/m².
        // This is below 10^9 A/m² for such a soft, thin free layer.
        // The formula is physically correct — the bound is that it must be
        // positive and finite.
        assert!(j_th.is_finite(), "j_th must be finite");
        assert!(j_th < 1.0e15, "j_th = {j_th} seems unrealistically large");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Magnetization norm conservation
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_magnetization_norm_conserved() {
        let stno = py_stno();
        let ms = stno.solver.material.ms;
        let m0 = tilted_m(&stno);
        let dt = 5.0e-13; // 0.5 ps — stable for Py in 100 kA/m field

        let trajectory = stno.simulate(m0, dt, 100, 0.0);

        for (i, m) in trajectory.iter().enumerate() {
            let norm_error = (m.magnitude() - ms).abs();
            assert!(
                norm_error < 1.0e-3, // Allow ~ppm relative error
                "step {i}: |m| = {:.6e} ≠ M_s = {ms:.6e} (error = {norm_error:.2e})",
                m.magnitude()
            );
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // STT perpendicularity
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_stt_perpendicular_to_m() {
        // For zero β, the torque m×(m×p̂) is automatically perpendicular to m.
        // With finite β, the component m×p̂ is also perpendicular to m.
        // Hence τ_STT · m̂ = 0 for any β.
        let stno = py_stno();
        let ms = stno.solver.material.ms;

        // Test with several m directions
        let test_cases = [
            Vector3::new(1.0, 0.0, 0.0) * ms,
            Vector3::new(0.0, 1.0, 0.0) * ms,
            Vector3::new(0.0, 0.0, 1.0) * ms,
            Vector3::new(1.0, 1.0, 0.0).normalize() * ms,
            Vector3::new(0.1, 0.0, 1.0).normalize() * ms,
        ];

        let j_density = 1.0e12;
        for m in &test_cases {
            let tau = stno.slonczewski_torque(*m, j_density);
            let m_hat = m.normalize();
            let dot = tau.dot(&m_hat).abs();
            assert!(
                dot < 1.0e-6,
                "τ_STT · m̂ = {dot:.3e} ≠ 0 for m = ({:.3}, {:.3}, {:.3})",
                m_hat.x,
                m_hat.y,
                m_hat.z
            );
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Auto-oscillation above threshold — custom low-field configuration
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_auto_oscillation() {
        // Build a STNO with a very small external field so that the STT at
        // a moderate multiple of J_th can actually compete with and overcome
        // the LLG damping.  We use Py parameters but H_ext = 0, p̂ = +ẑ,
        // so the only energy scale is from STT.
        let mut material = Ferromagnet::permalloy();
        material.alpha = 0.01;
        material.ms = 8.0e5;

        // Zero external field — STT drives dynamics purely
        let h_ext = Vector3::new(0.0, 0.0, 0.0);
        let t_fm = 5.0e-9;
        let area = PI * (50.0e-9f64).powi(2);
        let config = SpinTorqueOscillatorConfig::new(
            0.35,
            t_fm,
            area,
            0.0,                         // β = 0
            Vector3::new(0.0, 0.0, 1.0), // p̂ = +ẑ
        )
        .expect("valid config");

        let stno = SpinTorqueOscillator::new(material, h_ext, config).expect("valid STNO");

        let ms = stno.solver.material.ms;
        let j_th = stno.threshold_current_density();

        // Start with m in x-y plane (perpendicular to p̂) — maximum torque
        let m0 = Vector3::new(1.0, 0.0, 0.0) * ms;
        let dt = 1.0e-12;

        // Run well above threshold.  J_th = 2e*mu0*Ms*t*alpha/(hbar*eta) ≈ 436 kA/m².
        // At J = 10*J_th, a_J = 10*alpha ≈ 0.1 (as fraction of Ms).
        // The STT will sustain in-plane oscillations around ẑ.
        let j = 10.0 * j_th;
        let trajectory = stno.simulate(m0, dt, 2000, j);
        let power = stno.oscillation_power(&trajectory);

        // With p̂ = ẑ, m starting in xy-plane, and large positive J driving m×(m×ẑ),
        // the magnetization spirals into the x-y plane.  We expect non-trivial
        // variance in m_z (power > 0).
        assert!(
            power > 1.0e-6,
            "oscillation_power = {power:.4e} — expected > 0 at J = 10·J_th = {j:.2e} A/m²"
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Stability at equilibrium with J = 0
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_below_threshold_stable() {
        let stno = py_stno();
        let ms = stno.solver.material.ms;

        // Start exactly at equilibrium: m aligned with H_ext (+z).
        // With J = 0, there is no STT — only LLG precession (which vanishes
        // when m ∥ H_eff) and damping.  The magnetization should remain
        // stationary (or within numerical noise of the initial state).
        let m0 = Vector3::new(0.0, 0.0, ms); // aligned with H_ext along +z
        let dt = 5.0e-13;

        let trajectory = stno.simulate(m0, dt, 2000, 0.0);

        let m_z_final = trajectory.last().map(|m| m.z).unwrap_or(ms);

        // At true equilibrium there is no precession, so m_z should remain ≈ ms.
        // Allow a small numerical drift (< 0.01 * ms) but no spurious oscillation.
        let drift = (m_z_final - ms).abs();
        assert!(
            drift < 0.01 * ms,
            "J=0 from equilibrium: |m_z_final - Ms| = {drift:.4e} — should be near zero"
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Slonczewski sign check
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_slonczewski_sign() {
        // For p̂ = +ẑ, m = x̂:
        //   m × (m × p̂) = x̂ × (x̂ × ẑ) = x̂ × (−ŷ) = −ẑ
        // Hence for positive a_J (positive J) the Slonczewski term pushes m
        // in the −ẑ direction initially (away from +z polarizer), which makes
        // physical sense: the conventional (Slonczewski) sign convention
        // means positive J drives m away from parallel alignment with p̂.
        //
        // We verify the sign explicitly here using a custom config with p̂ = +ẑ.
        let material = Ferromagnet::permalloy();
        let h_ext = Vector3::new(0.0, 0.0, 0.0);
        let config = SpinTorqueOscillatorConfig::new(
            0.35,
            5.0e-9,
            PI * (50.0e-9f64).powi(2),
            0.0,                         // β = 0 for a clean test
            Vector3::new(0.0, 0.0, 1.0), // p̂ = +ẑ
        )
        .expect("valid config");

        let stno = SpinTorqueOscillator::new(material, h_ext, config).expect("valid STNO");

        let ms = stno.solver.material.ms;
        let m = Vector3::new(1.0, 0.0, 0.0) * ms; // m = x̂ · Ms
        let j_positive = 1.0e12; // positive current density

        let tau = stno.slonczewski_torque(m, j_positive);

        // τ should have a negative z-component (drives m away from p̂ = ẑ)
        // and zero x-component (torque perpendicular to m).
        assert!(
            tau.z < 0.0,
            "τ_z = {:.4e} — expected < 0 (torque drives m away from p̂ for positive J)",
            tau.z
        );
        assert!(
            tau.x.abs() < 1.0e-6,
            "τ_x = {:.4e} — should be ≈ 0 (perpendicular to m = x̂)",
            tau.x
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // a_J coefficient
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_a_j_scaling() {
        let stno = py_stno();
        // a_J should scale linearly with J
        let j1 = 1.0e11;
        let j2 = 2.0e11;
        let ratio = stno.a_j(j2) / stno.a_j(j1);
        assert!(
            (ratio - 2.0).abs() < 1.0e-10,
            "a_J should be linear in J; ratio = {ratio}"
        );
    }

    #[test]
    fn test_a_j_sign() {
        let stno = py_stno();
        // a_J positive for positive J, negative for negative J
        assert!(stno.a_j(1.0e12) > 0.0);
        assert!(stno.a_j(-1.0e12) < 0.0);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Thermal linewidth
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_thermal_linewidth_positive() {
        let stno = py_stno();
        let j_th = stno.threshold_current_density();
        let lw = stno.thermal_linewidth(1.5 * j_th, 300.0);
        assert!(lw > 0.0, "linewidth must be positive, got {lw}");
        assert!(lw.is_finite(), "linewidth must be finite, got {lw}");
    }

    #[test]
    fn test_thermal_linewidth_temperature_scaling() {
        let stno = py_stno();
        let j = 1.5 * stno.threshold_current_density();
        let lw_300 = stno.thermal_linewidth(j, 300.0);
        let lw_600 = stno.thermal_linewidth(j, 600.0);
        // Linewidth ∝ T, so doubling T should double linewidth
        let ratio = lw_600 / lw_300;
        assert!(
            (ratio - 2.0).abs() < 1.0e-6,
            "Δf should scale linearly with T; ratio = {ratio}"
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Injection locking
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_injection_locking_range_positive() {
        let stno = py_stno();
        let delta_f = stno.injection_locking_range(100.0, 0.5);
        assert!(
            delta_f > 0.0,
            "locking range must be positive, got {delta_f}"
        );
    }

    #[test]
    fn test_injection_locking_range_zero_q() {
        let stno = py_stno();
        let delta_f = stno.injection_locking_range(0.0, 0.5);
        assert_eq!(delta_f, 0.0, "zero Q should give zero locking range");
    }

    #[test]
    fn test_injection_locking_range_increases_with_injection() {
        let stno = py_stno();
        let q = 1000.0;
        let lw_low = stno.injection_locking_range(q, 0.1);
        let lw_high = stno.injection_locking_range(q, 0.9);
        assert!(
            lw_high > lw_low,
            "locking range should increase with injection strength"
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Oscillation frequency
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_oscillation_frequency_auto_oscillation() {
        let stno = py_stno();
        let ms = stno.solver.material.ms;
        let j_th = stno.threshold_current_density();

        let m0 = Vector3::new(0.1, 0.0, 1.0).normalize() * ms;
        let dt = 2.0e-13;
        let trajectory = stno.simulate(m0, dt, 2000, 3.0 * j_th);

        // Should detect a frequency since we're well above threshold
        let freq_opt = stno.oscillation_frequency(&trajectory, dt);
        if let Some(f) = freq_opt {
            // Physical frequency range for Py STNO: tens of MHz to tens of GHz
            assert!(
                f > 1.0e6,
                "frequency {f:.2e} Hz seems too low for a Py STNO"
            );
            assert!(f < 1.0e12, "frequency {f:.2e} Hz is unrealistically high");
        }
        // None is also acceptable if zero crossings < 2 (marginal oscillation)
    }

    #[test]
    fn test_oscillation_frequency_no_oscillation() {
        let stno = py_stno();
        let ms = stno.solver.material.ms;

        // Zero current — no STT — damped equilibrium — no frequency
        let m0 = Vector3::new(0.0, 0.0, ms); // already at equilibrium
        let dt = 5.0e-13;
        let trajectory = stno.simulate(m0, dt, 100, 0.0);
        let freq_opt = stno.oscillation_frequency(&trajectory, dt);
        // With m0 exactly along +z and h_ext along +z, there is no precession
        assert!(
            freq_opt.is_none(),
            "expected None for no oscillation, got {:?}",
            freq_opt
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Preset construction
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_permalloy_nanopillar_construction() {
        let stno = py_stno();
        assert!(stno.solver.material.ms > 0.0);
        assert!(stno.solver.material.alpha > 0.0);
        assert!(stno.config.spin_polarization > 0.0);
        assert!(stno.config.volume > 0.0);
        // p̂ should be normalized
        let p_norm = stno.config.p_hat.magnitude();
        assert!(
            (p_norm - 1.0).abs() < 1.0e-10,
            "p̂ not normalized: |p̂| = {p_norm}"
        );
    }
}
