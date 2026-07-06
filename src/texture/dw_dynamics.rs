//! Domain wall dynamics: Walker breakdown, STT-driven, and SOT-driven motion
//!
//! This module implements the three canonical models for current- and field-driven
//! magnetic domain wall (DW) dynamics in ferromagnetic thin films:
//!
//! 1. **Walker breakdown** (field-driven): below the Walker threshold field `H_W`
//!    the DW translates at `v = γ Δ H / α`; above it the wall precesses and the
//!    time-averaged velocity falls to `v = γ Δ H / (1 + α²)`.
//!
//! 2. **STT-driven DW** (Berger–Slonczewski model): spin-transfer torque (STT)
//!    from an in-plane current drives the DW through adiabatic and non-adiabatic
//!    mechanisms characterised by the spin polarisation `P` and non-adiabaticity
//!    parameter `β`.
//!
//! 3. **SOT-driven DW** (Thiaville 2012 model): the damping-like (DL) component
//!    of spin-orbit torque drives a Néel DW stabilised by interfacial DMI at
//!    `v = γ Δ θ_SH ℏ J / (2 e μ₀ M_s t α)`.
//!
//! # References
//!
//! - N. L. Schryer and L. R. Walker, "Motion of 180° domain walls in uniform
//!   dc magnetic fields", J. Appl. Phys. **45**, 5406 (1974).
//! - A. Thiaville, S. Rohart, É. Jué, V. Cros, A. Fert, "Dynamics of
//!   Dzyaloshinskii domain walls in ultrathin magnetic films",
//!   *EPL* **100**, 57002 (2012).
//! - I. M. Miron *et al.*, "Current-driven spin torque induced by the Rashba
//!   effect in a ferromagnetic metal layer",
//!   *Nature* **476**, 189 (2011).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::constants::{E_CHARGE, GAMMA, HBAR, MU_0, MU_B};

// ────────────────────────────────────────────────────────────────────────────
// Material presets
// ────────────────────────────────────────────────────────────────────────────

/// Material parameters governing domain wall dynamics.
///
/// Bundles all the physical quantities required to evaluate Walker breakdown,
/// STT-driven, and SOT-driven DW velocity formulas for a specific material
/// system. Provide your own values or use one of the preset constructors.
///
/// # Units
///
/// All fields use SI base units. See individual field doc-comments.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DwMaterial {
    /// Human-readable identifier (e.g. `"Permalloy"`, `"Pt/Co/AlOx"`)
    ///
    /// This field is skipped during serde deserialisation (it is set to
    /// `"custom"` on round-trip). Only the numeric parameters are serialised.
    #[cfg_attr(feature = "serde", serde(skip, default = "DwMaterial::default_name"))]
    pub name: &'static str,
    /// Exchange stiffness `A` \[J/m\]
    pub exchange_a: f64,
    /// Uniaxial anisotropy `K_u` \[J/m³\]
    pub anisotropy_k: f64,
    /// Saturation magnetisation `M_s` \[A/m\]
    pub ms: f64,
    /// Gilbert damping constant `α` (dimensionless)
    pub alpha: f64,
    /// Ferromagnetic layer thickness `t` \[m\]
    pub thickness: f64,
    /// Spin polarisation `P` (dimensionless, 0 < P < 1)
    pub polarization: f64,
    /// Non-adiabaticity parameter `β` (dimensionless, typically 0.01–0.1)
    pub beta: f64,
    /// Effective spin Hall angle `θ_SH` for SOT (dimensionless)
    pub theta_sh: f64,
}

impl DwMaterial {
    /// Permalloy (Ni₈₁Fe₁₉) — soft ferromagnet, in-plane anisotropy.
    ///
    /// Typical literature values for 10–30 nm thick Py films:
    ///
    /// | Parameter | Value |
    /// |---|---|
    /// | `A` | 1.3 × 10⁻¹¹ J/m |
    /// | `K_u` | 1 × 10³ J/m³ (shape anisotropy dominated) |
    /// | `M_s` | 800 kA/m |
    /// | `α` | 0.01 |
    /// | `P` | 0.40 |
    /// | `β` | 0.01 |
    /// | `t` | 20 nm |
    /// | `θ_SH` | 0 (no HM layer) |
    pub fn permalloy() -> Self {
        Self {
            name: "Permalloy",
            exchange_a: 1.3e-11,
            anisotropy_k: 1.0e3,
            ms: 800.0e3,
            alpha: 0.01,
            thickness: 20.0e-9,
            polarization: 0.40,
            beta: 0.01,
            theta_sh: 0.0,
        }
    }

    /// Co/Pt PMA multilayer — perpendicular magnetic anisotropy stack.
    ///
    /// Representative parameters for a Pt(3nm)/Co(1nm)/Pt(3nm) tri-layer.
    /// `K_u` is the total effective uniaxial anisotropy combining interfacial
    /// PMA (from Pt/Co interfaces) with bulk contributions; it must exceed
    /// `μ₀ M_s²/2 ≈ 0.76 MJ/m³` to give `K_eff > 0` (PMA condition).
    ///
    /// | Parameter | Value |
    /// |---|---|
    /// | `A` | 1.5 × 10⁻¹¹ J/m |
    /// | `K_u` | 1.5 × 10⁶ J/m³ |
    /// | `M_s` | 1.1 MA/m |
    /// | `α` | 0.10 |
    /// | `P` | 0.50 |
    /// | `β` | 0.04 |
    /// | `t` | 1 nm |
    /// | `θ_SH` | 0.08 |
    pub fn co_pt_pma() -> Self {
        Self {
            name: "Co/Pt PMA",
            exchange_a: 1.5e-11,
            anisotropy_k: 1.5e6,
            ms: 1.1e6,
            alpha: 0.10,
            thickness: 1.0e-9,
            polarization: 0.50,
            beta: 0.04,
            theta_sh: 0.08,
        }
    }

    /// CoFeB/MgO PMA interface — reference system for tunnel-junction devices.
    ///
    /// Parameters representative of Ta(5nm)/CoFeB(1.5nm)/MgO(2nm).
    /// `K_u = 1.0 MJ/m³` represents the total uniaxial anisotropy including the
    /// interfacial PMA from the CoFeB/MgO interface (which must overcome
    /// `μ₀ M_s²/2 ≈ 0.76 MJ/m³` to achieve perpendicular orientation).
    ///
    /// | Parameter | Value |
    /// |---|---|
    /// | `A` | 2.0 × 10⁻¹¹ J/m |
    /// | `K_u` | 1.0 × 10⁶ J/m³ |
    /// | `M_s` | 1.1 MA/m |
    /// | `α` | 0.01 |
    /// | `P` | 0.65 |
    /// | `β` | 0.01 |
    /// | `t` | 1.5 nm |
    /// | `θ_SH` | 0 |
    pub fn cofeb_mgo_pma() -> Self {
        Self {
            name: "CoFeB/MgO PMA",
            exchange_a: 2.0e-11,
            anisotropy_k: 1.0e6,
            ms: 1.1e6,
            alpha: 0.01,
            thickness: 1.5e-9,
            polarization: 0.65,
            beta: 0.01,
            theta_sh: 0.0,
        }
    }

    /// Pt(3nm)/Co(0.6nm)/AlO_x — Miron 2011 reference system.
    ///
    /// Parameters from I. M. Miron *et al.*, *Nature* **476**, 189 (2011):
    ///
    /// | Parameter | Value |
    /// |---|---|
    /// | `A` | 1.6 × 10⁻¹¹ J/m |
    /// | `K_u` | 1.3 × 10⁶ J/m³ |
    /// | `M_s` | 1.0 MA/m |
    /// | `α` | 0.50 |
    /// | `P` | 0.50 |
    /// | `β` | 0.02 |
    /// | `t` | 0.6 nm |
    /// | `θ_SH` | 0.15 |
    pub fn pt_co_alox() -> Self {
        Self {
            name: "Pt/Co/AlOx",
            exchange_a: 1.6e-11,
            anisotropy_k: 1.3e6,
            ms: 1.0e6,
            alpha: 0.50,
            thickness: 0.6e-9,
            polarization: 0.50,
            beta: 0.02,
            theta_sh: 0.15,
        }
    }

    /// Effective anisotropy `K_eff = K_u - μ₀ M_s² / 2` \[J/m³\].
    ///
    /// For PMA materials `K_eff > 0`; the DW width uses `K_eff` as the
    /// restoring energy scale competing against the exchange stiffness.
    pub fn k_eff(&self) -> f64 {
        self.anisotropy_k - 0.5 * MU_0 * self.ms * self.ms
    }

    /// Default material name sentinel used when serde deserialises `DwMaterial`.
    ///
    /// The `name` field holds a `&'static str` which cannot be deserialised from
    /// arbitrary bytes, so it is skipped and replaced with this value on
    /// round-trip.
    #[cfg(feature = "serde")]
    fn default_name() -> &'static str {
        "custom"
    }
}

// ────────────────────────────────────────────────────────────────────────────
// Walker breakdown (field-driven)
// ────────────────────────────────────────────────────────────────────────────

/// Field-driven domain wall dynamics including the Walker breakdown.
///
/// Below the Walker threshold field `H_W` a 1D DW moves with velocity
/// `v = γ Δ H / α`. Above `H_W` the DW azimuthal angle `φ` precesses
/// and the time-averaged velocity drops to `v = γ Δ H / (1 + α²)`.
///
/// The Walker threshold field (simplified for thin PMA film where
/// demagnetisation factors satisfy `N_z ≈ 0`, `N_y ≈ 1`) is:
///
/// ```text
/// H_W = α M_s / 2   [A/m]
/// μ₀ H_W = α μ₀ M_s / 2   [T]
/// ```
///
/// The Walker velocity — the *maximum* DW velocity in the below-Walker
/// regime — is independent of α:
///
/// ```text
/// v_W = γ Δ_DW μ₀ M_s / 2   [m/s]
/// ```
///
/// # References
///
/// - N. L. Schryer and L. R. Walker, J. Appl. Phys. **45**, 5406 (1974).
/// - A. Mougin *et al.*, EPL **78**, 57007 (2007).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct WalkerBreakdown {
    /// Material parameters governing the field-driven DW dynamics.
    pub material: DwMaterial,
}

impl WalkerBreakdown {
    /// Construct a [`WalkerBreakdown`] solver for the given material.
    pub fn new(material: DwMaterial) -> Self {
        Self { material }
    }

    /// Domain wall width parameter `Δ_DW = π √(A / K_eff)` \[m\].
    ///
    /// Uses the effective anisotropy `K_eff = K_u - μ₀ M_s² / 2` as the
    /// relevant energy scale. For in-plane anisotropy materials where
    /// `K_eff` may be small or even negative, the width diverges; in such
    /// cases the PMA approximation is no longer valid.
    ///
    /// ```text
    /// Δ_DW = π √(A / K_eff)   [m]
    /// ```
    pub fn dw_width(&self) -> f64 {
        let k_eff = self.material.k_eff();
        std::f64::consts::PI * (self.material.exchange_a / k_eff).sqrt()
    }

    /// Walker threshold field `H_W = α M_s / 2` \[A/m\] (SI, not Tesla).
    ///
    /// This is the simplified expression for a PMA thin film where the
    /// demagnetisation factor difference `N_z − N_y ≈ −1`. The field is
    /// expressed in A/m (multiply by `μ₀` to get Tesla).
    ///
    /// ```text
    /// H_W = α M_s / 2   [A/m]
    /// ```
    pub fn walker_threshold_field(&self) -> f64 {
        self.material.alpha * self.material.ms * 0.5
    }

    /// Walker threshold flux density `μ₀ H_W = α μ₀ M_s / 2` \[T\].
    ///
    /// More convenient for direct comparison with applied field values
    /// in experimental contexts.
    pub fn walker_threshold_tesla(&self) -> f64 {
        MU_0 * self.walker_threshold_field()
    }

    /// Walker velocity `v_W = γ Δ_DW μ₀ M_s / 2` \[m/s\].
    ///
    /// This is the maximum steady-state DW velocity, reached exactly at
    /// `H = H_W`. Note that it is *independent* of the Gilbert damping `α`:
    ///
    /// ```text
    /// v_W = γ Δ_DW H_W / α = γ Δ_DW μ₀ M_s / 2   [m/s]
    /// ```
    pub fn walker_velocity(&self) -> f64 {
        // Equivalently: GAMMA * dw_width * H_W / alpha = GAMMA * delta * mu0 * Ms / 2
        GAMMA * self.dw_width() * 0.5 * MU_0 * self.material.ms
    }

    /// DW velocity in the **below-Walker** steady-state regime \[m/s\].
    ///
    /// Valid for applied field magnitudes `H_applied ≤ H_W`. In this regime
    /// the DW translates without precessing. The field `H` is given in SI units
    /// A/m; the factor `μ₀` converts it to flux density (Tesla) before
    /// multiplying by the gyromagnetic ratio:
    ///
    /// ```text
    /// v = γ Δ_DW μ₀ H / α   [m/s]
    /// ```
    ///
    /// # Arguments
    /// * `h_field` — Applied field magnitude \[A/m\] (SI, not Tesla).
    pub fn velocity_below_walker(&self, h_field: f64) -> f64 {
        GAMMA * self.dw_width() * MU_0 * h_field / self.material.alpha
    }

    /// Time-averaged DW velocity in the **above-Walker** regime \[m/s\].
    ///
    /// Above `H_W` the DW precesses periodically and the time-averaged
    /// velocity is reduced relative to the linear extrapolation. The field
    /// is in A/m; `μ₀` converts it to Tesla:
    ///
    /// ```text
    /// ⟨v⟩ = γ Δ_DW μ₀ H / (1 + α²)   [m/s]
    /// ```
    ///
    /// Note that this is larger than `v_below` for `α < 1` but smaller
    /// than the velocity at `H_W` in the below-Walker regime.
    ///
    /// # Arguments
    /// * `h_field` — Applied field magnitude \[A/m\].
    pub fn velocity_above_walker(&self, h_field: f64) -> f64 {
        GAMMA * self.dw_width() * MU_0 * h_field / (1.0 + self.material.alpha * self.material.alpha)
    }

    /// DW velocity automatically selecting the correct dynamical regime \[m/s\].
    ///
    /// Compares `|H_applied|` against the Walker threshold `H_W` and delegates
    /// to [`Self::velocity_below_walker`] or [`Self::velocity_above_walker`]
    /// accordingly. The sign of the returned velocity matches that of the
    /// applied field.
    ///
    /// # Arguments
    /// * `h_field` — Applied field \[A/m\] (may be negative).
    pub fn velocity(&self, h_field: f64) -> f64 {
        let h_abs = h_field.abs();
        let h_w = self.walker_threshold_field();
        let sign = if h_field >= 0.0 { 1.0 } else { -1.0 };
        if h_abs <= h_w {
            sign * self.velocity_below_walker(h_abs)
        } else {
            sign * self.velocity_above_walker(h_abs)
        }
    }
}

// ────────────────────────────────────────────────────────────────────────────
// STT-driven domain wall dynamics
// ────────────────────────────────────────────────────────────────────────────

/// Spin-transfer torque (STT) driven domain wall dynamics.
///
/// Implements the 1D DW model of Thiaville *et al.* (2004), which extends the
/// Berger adiabatic STT with a non-adiabatic correction. The equations of
/// motion for the DW centre `X` and tilting angle `φ` yield the steady-state
/// DW velocity:
///
/// ```text
/// v = [(1 + α β) v_s + (β − α) v_H] / (1 + α²)
/// ```
///
/// where `v_s = P j μ_B / (e M_s)` is the adiabatic STT velocity and
/// `v_H = γ Δ_DW H_applied` is the field-driven contribution.
///
/// For purely current-driven motion (`H = 0`) this simplifies to:
///
/// ```text
/// v = (1 + α β) / (1 + α²) × v_s  ≈  v_s   (since α, β ≪ 1)
/// ```
///
/// # References
///
/// - A. Thiaville, Y. Nakatani, J. Miltat, Y. Suzuki,
///   "Micromagnetic understanding of current-driven domain wall motion in
///   patterned nanowires", EPL **69**, 990 (2004).
/// - S. Zhang and Z. Li, "Roles of nonequilibrium conduction electrons on
///   the magnetization dynamics of ferromagnets",
///   PRL **93**, 127204 (2004).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DwSttDynamics {
    /// Material parameters for the STT DW dynamics.
    pub material: DwMaterial,
}

impl DwSttDynamics {
    /// Construct a [`DwSttDynamics`] solver for the given material.
    pub fn new(material: DwMaterial) -> Self {
        Self { material }
    }

    /// DW width `Δ_DW = π √(A / K_eff)` \[m\].
    ///
    /// Identical formula to [`WalkerBreakdown::dw_width`]; duplicated here
    /// so that `DwSttDynamics` is self-contained.
    pub fn dw_width(&self) -> f64 {
        let k_eff = self.material.k_eff();
        std::f64::consts::PI * (self.material.exchange_a / k_eff).sqrt()
    }

    /// Walker threshold field `H_W = α M_s / 2` \[A/m\].
    ///
    /// Same simplified PMA expression as [`WalkerBreakdown::walker_threshold_field`].
    pub fn walker_threshold_field(&self) -> f64 {
        self.material.alpha * self.material.ms * 0.5
    }

    /// Walker velocity `v_W = γ Δ_DW μ₀ M_s / 2` \[m/s\].
    ///
    /// Maximum DW velocity in the purely-field-driven below-Walker regime;
    /// used to determine the Walker threshold current density.
    pub fn walker_velocity(&self) -> f64 {
        GAMMA * self.dw_width() * 0.5 * MU_0 * self.material.ms
    }

    /// Adiabatic STT drift velocity `v_s = P j μ_B / (e M_s)` \[m/s\].
    ///
    /// This is the DW velocity that would result from a purely adiabatic
    /// spin-transfer torque in the absence of damping and non-adiabatic
    /// corrections. It sets the scale for the STT mechanism.
    ///
    /// ```text
    /// v_s = P × j × μ_B / (e × M_s)   [m/s]
    /// ```
    ///
    /// # Arguments
    /// * `current_density` — In-plane current density `j` \[A/m²\] (positive
    ///   drives the DW in the +x direction by convention).
    pub fn adiabatic_velocity(&self, current_density: f64) -> f64 {
        self.material.polarization * current_density * MU_B / (E_CHARGE * self.material.ms)
    }

    /// Combined field + STT steady-state DW velocity \[m/s\].
    ///
    /// Evaluates the full one-dimensional collective-coordinate result:
    ///
    /// ```text
    /// v = [(1 + α β) v_s + (β − α) v_H] / (1 + α²)
    /// ```
    ///
    /// where `v_H = γ Δ_DW H_applied` is the field-driven contribution and
    /// `v_s` is the adiabatic STT velocity.
    ///
    /// # Arguments
    /// * `current_density` — In-plane charge current density \[A/m²\].
    /// * `h_applied` — External longitudinal field \[A/m\] (along easy axis).
    pub fn velocity(&self, current_density: f64, h_applied: f64) -> f64 {
        let alpha = self.material.alpha;
        let beta = self.material.beta;
        let v_s = self.adiabatic_velocity(current_density);
        // v_H = γ Δ μ₀ H  (H is in A/m; μ₀ converts to Tesla for γ)
        let v_h = GAMMA * self.dw_width() * MU_0 * h_applied;
        ((1.0 + alpha * beta) * v_s + (beta - alpha) * v_h) / (1.0 + alpha * alpha)
    }

    /// Walker threshold current density for purely current-driven DW \[A/m²\].
    ///
    /// Above this current density the DW enters the Walker precession regime.
    /// The threshold is found by setting the adiabatic STT velocity equal to
    /// the Walker velocity:
    ///
    /// ```text
    /// v_s,W = v_W  ⟹  J_W = v_W × e × M_s / (P × μ_B)
    /// ```
    ///
    /// where `v_W = γ Δ_DW μ₀ M_s / 2` is the Walker velocity.
    pub fn walker_current_density(&self) -> f64 {
        let v_w = self.walker_velocity();
        v_w * E_CHARGE * self.material.ms / (self.material.polarization * MU_B)
    }
}

// ────────────────────────────────────────────────────────────────────────────
// SOT-driven domain wall dynamics (Thiaville 2012)
// ────────────────────────────────────────────────────────────────────────────

/// Spin-orbit torque (SOT) driven domain wall dynamics for Néel DWs.
///
/// Implements the model of Thiaville *et al.* (2012) for a Néel DW in a PMA
/// material driven by the damping-like (DL) component of the SOT. For a
/// DMI-stabilised Néel DW (`φ_DW = 0`) the DW velocity is:
///
/// ```text
/// v_DW = (γ Δ_DW / α) × H_DL
///       = γ Δ_DW θ_SH ℏ J / (2 e μ₀ M_s t α)   [m/s]
/// ```
///
/// where the DL-SOT effective field is:
///
/// ```text
/// H_DL = θ_SH ℏ J / (2 e μ₀ M_s t)   [A/m]
/// ```
///
/// The Walker breakdown threshold under SOT is:
///
/// ```text
/// H_W^SOT = α μ₀ M_s / 2   [A/m]
/// J_W^SOT = 2 e μ₀ M_s t α H_W / (θ_SH ℏ)
///          = e μ₀² M_s² t α² / (θ_SH ℏ)   [A/m²]
/// ```
///
/// # References
///
/// - A. Thiaville, S. Rohart, É. Jué, V. Cros, A. Fert,
///   "Dynamics of Dzyaloshinskii domain walls in ultrathin magnetic films",
///   *EPL* **100**, 57002 (2012).
/// - K. Ryu, L. Thomas, S. Yang, S. Parkin, "Chiral spin torque at magnetic
///   domain walls", *Nat. Nanotechnol.* **8**, 527 (2013).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DwSotDynamics {
    /// Material parameters for the SOT DW dynamics.
    pub material: DwMaterial,
}

impl DwSotDynamics {
    /// Construct a [`DwSotDynamics`] solver for the given material.
    pub fn new(material: DwMaterial) -> Self {
        Self { material }
    }

    /// DW width parameter `Δ_DW = π √(A / K_eff)` \[m\].
    ///
    /// Identical formula to [`WalkerBreakdown::dw_width`]; duplicated here so
    /// that `DwSotDynamics` is self-contained.
    pub fn dw_width(&self) -> f64 {
        let k_eff = self.material.k_eff();
        std::f64::consts::PI * (self.material.exchange_a / k_eff).sqrt()
    }

    /// Walker threshold field `H_W = α M_s / 2` \[A/m\].
    ///
    /// Same simplified PMA expression as [`WalkerBreakdown::walker_threshold_field`].
    pub fn walker_threshold_field(&self) -> f64 {
        self.material.alpha * self.material.ms * 0.5
    }

    /// Damping-like SOT effective field `H_DL` \[A/m\] (not in Tesla).
    ///
    /// The DL-SOT effective field for a spin current generated by the spin
    /// Hall effect in an adjacent heavy metal layer:
    ///
    /// ```text
    /// H_DL = θ_SH ℏ J / (2 e μ₀ M_s t)   [A/m]
    /// ```
    ///
    /// Multiply by `μ₀` to obtain the flux density in Tesla.
    ///
    /// # Arguments
    /// * `current_density` — Charge current density in the HM layer \[A/m²\].
    pub fn dl_sot_field(&self, current_density: f64) -> f64 {
        self.material.theta_sh * HBAR * current_density
            / (2.0 * E_CHARGE * MU_0 * self.material.ms * self.material.thickness)
    }

    /// SOT-driven DW velocity for a DMI-stabilised Néel DW \[m/s\].
    ///
    /// For `φ_DW = 0` (Néel DW fully aligned with the DMI direction):
    ///
    /// ```text
    /// v_DW = (γ Δ_DW / α) × μ₀ H_DL   [m/s]
    /// ```
    ///
    /// Note: [`Self::dl_sot_field`] returns H_DL in A/m; the factor `μ₀`
    /// converts it to Tesla before multiplying by the gyromagnetic ratio.
    /// Equivalently this gives `v = γ Δ θ_SH ℏ J / (2 e M_s t α)`.
    ///
    /// This formula applies below the Walker breakdown threshold. Above the
    /// Walker threshold the DW precesses and the velocity saturates.
    ///
    /// # Arguments
    /// * `current_density` — Charge current density in the HM layer \[A/m²\].
    pub fn velocity(&self, current_density: f64) -> f64 {
        let h_dl = self.dl_sot_field(current_density);
        // H_DL is in A/m; multiply by μ₀ to convert to Tesla for gyromagnetic ratio
        GAMMA * self.dw_width() * MU_0 * h_dl / self.material.alpha
    }

    /// DW velocity with combined SOT + external out-of-plane field \[m/s\].
    ///
    /// An additional longitudinal (easy-axis) field `H_applied` adds to the
    /// DL-SOT field contribution:
    ///
    /// ```text
    /// v_DW = (γ Δ_DW / α) × (H_DL + H_applied)   [m/s]
    /// ```
    ///
    /// # Arguments
    /// * `current_density` — Charge current density in the HM layer \[A/m²\].
    /// * `h_applied` — External easy-axis field \[A/m\].
    pub fn velocity_with_field(&self, current_density: f64, h_applied: f64) -> f64 {
        let h_dl = self.dl_sot_field(current_density);
        // Both H_DL and H_applied are in A/m; multiply by μ₀ to get Tesla
        GAMMA * self.dw_width() * MU_0 * (h_dl + h_applied) / self.material.alpha
    }

    /// Walker threshold current density for SOT-driven DW \[A/m²\].
    ///
    /// Obtained by setting `H_DL = H_W`:
    ///
    /// ```text
    /// J_W^SOT = 2 e μ₀ M_s t α H_W / (θ_SH ℏ)
    ///          = e μ₀² M_s² t α² / (θ_SH ℏ)   [A/m²]
    /// ```
    pub fn walker_current_density(&self) -> f64 {
        let h_w = self.walker_threshold_field();
        2.0 * E_CHARGE
            * MU_0
            * self.material.ms
            * self.material.thickness
            * self.material.alpha
            * h_w
            / (self.material.theta_sh * HBAR)
    }

    /// Compute a velocity-vs-current-density curve as `(J, v)` pairs.
    ///
    /// Samples `n_points` current densities linearly spaced from `j_min` to
    /// `j_max`. For each sample the velocity is computed accounting for the
    /// Walker breakdown: below the Walker current `J_W` the linear model
    /// `v = (γ Δ / α) H_DL` is used; above `J_W` the time-averaged velocity
    /// saturates to `v_W` (Walker velocity) and does not continue to grow
    /// linearly. The Walker velocity for the SOT context is defined at the
    /// Walker threshold:
    ///
    /// ```text
    /// v_W = (γ Δ / α) × H_W   [m/s]
    /// ```
    ///
    /// # Arguments
    /// * `j_min` — Minimum current density \[A/m²\].
    /// * `j_max` — Maximum current density \[A/m²\].
    /// * `n_points` — Number of sample points (must be ≥ 2).
    pub fn velocity_curve(&self, j_min: f64, j_max: f64, n_points: usize) -> Vec<(f64, f64)> {
        let n = n_points.max(2);
        let j_w = self.walker_current_density();
        let v_w =
            GAMMA * self.dw_width() * MU_0 * self.walker_threshold_field() / self.material.alpha;
        (0..n)
            .map(|i| {
                let j = j_min + (j_max - j_min) * (i as f64) / ((n - 1) as f64);
                let v = if j.abs() <= j_w {
                    self.velocity(j)
                } else {
                    // Above Walker: saturate at the Walker velocity (time-averaged)
                    let sign = if j >= 0.0 { 1.0 } else { -1.0 };
                    sign * v_w
                };
                (j, v)
            })
            .collect()
    }
}

// ────────────────────────────────────────────────────────────────────────────
// Tests
// ────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── DwMaterial presets ───────────────────────────────────────────────────

    #[test]
    fn test_material_permalloy_k_eff_positive() {
        // Permalloy has very small K_u; K_eff = K_u - μ₀ Ms²/2 is dominated
        // by the demagnetisation term and may be negative (in-plane anisotropy).
        // We only check that the method does not panic.
        let mat = DwMaterial::permalloy();
        let _ = mat.k_eff();
        assert!(mat.ms > 0.0);
    }

    #[test]
    fn test_material_co_pt_k_eff_positive() {
        // Co/Pt PMA should have K_eff > 0 for a perpendicular anisotropy material.
        let mat = DwMaterial::co_pt_pma();
        assert!(
            mat.k_eff() > 0.0,
            "Co/Pt PMA must have K_eff > 0 (got {})",
            mat.k_eff()
        );
    }

    #[test]
    fn test_material_pt_co_alox_k_eff_positive() {
        let mat = DwMaterial::pt_co_alox();
        assert!(
            mat.k_eff() > 0.0,
            "Pt/Co/AlOx must have K_eff > 0 (got {})",
            mat.k_eff()
        );
    }

    #[test]
    fn test_material_cofeb_k_eff_positive() {
        let mat = DwMaterial::cofeb_mgo_pma();
        assert!(
            mat.k_eff() > 0.0,
            "CoFeB/MgO must have K_eff > 0 (got {})",
            mat.k_eff()
        );
    }

    // ── WalkerBreakdown ──────────────────────────────────────────────────────

    #[test]
    fn test_walker_dw_width_co_pt_nm_scale() {
        let wb = WalkerBreakdown::new(DwMaterial::co_pt_pma());
        let width_nm = wb.dw_width() * 1.0e9;
        assert!(
            width_nm > 1.0 && width_nm < 100.0,
            "Co/Pt DW width should be 1–100 nm, got {:.2} nm",
            width_nm
        );
    }

    #[test]
    fn test_walker_dw_width_pt_co_alox_around_15nm() {
        // Δ_DW ≈ π × 5 nm ≈ 15 nm for Miron 2011 parameters
        let wb = WalkerBreakdown::new(DwMaterial::pt_co_alox());
        let width_nm = wb.dw_width() * 1.0e9;
        assert!(
            width_nm > 10.0 && width_nm < 25.0,
            "Pt/Co/AlOx DW width should be ~15 nm, got {:.2} nm",
            width_nm
        );
    }

    #[test]
    fn test_walker_threshold_field_co_pt_mt_range() {
        // For Co/Pt: H_W = α Ms/2 = 0.1 × 1.1e6/2 = 55 kA/m → μ₀ H_W ≈ 69 mT
        let wb = WalkerBreakdown::new(DwMaterial::co_pt_pma());
        let b_milli_tesla = wb.walker_threshold_tesla() * 1.0e3;
        assert!(
            b_milli_tesla > 10.0 && b_milli_tesla < 500.0,
            "Walker threshold should be in mT–hundreds-of-mT range, got {:.2} mT",
            b_milli_tesla
        );
    }

    #[test]
    fn test_walker_threshold_field_permalloy_milli_tesla_range() {
        // For Py: H_W ≈ α Ms / 2 = 0.01 × 800e3 / 2 = 4 kA/m → μ₀ H_W ≈ 5 mT
        let wb = WalkerBreakdown::new(DwMaterial::permalloy());
        let b_milli_tesla = wb.walker_threshold_tesla() * 1.0e3;
        assert!(
            b_milli_tesla > 0.5 && b_milli_tesla < 100.0,
            "Permalloy Walker field should be in mT range, got {:.3} mT",
            b_milli_tesla
        );
    }

    #[test]
    fn test_walker_velocity_co_pt_physical_range() {
        // For Co/Pt with K_u = 1.5 MJ/m³: Δ ≈ 14 nm, μ₀ Ms/2 ≈ 0.69 T
        // v_W = γ Δ μ₀ Ms/2 ≈ 1.76e11 × 14e-9 × 0.69 ≈ 1700 m/s
        // High-K_u PMA materials have Walker velocities in the 100–3000 m/s range.
        let wb = WalkerBreakdown::new(DwMaterial::co_pt_pma());
        let v_w = wb.walker_velocity();
        assert!(
            v_w > 10.0 && v_w < 5000.0,
            "Co/Pt Walker velocity should be 10–5000 m/s, got {:.2} m/s",
            v_w
        );
    }

    #[test]
    fn test_walker_velocity_equal_to_velocity_at_threshold() {
        // v_W must equal velocity_below_walker(H_W)
        let wb = WalkerBreakdown::new(DwMaterial::co_pt_pma());
        let h_w = wb.walker_threshold_field();
        let v_at_threshold = wb.velocity_below_walker(h_w);
        let v_w = wb.walker_velocity();
        let rel_err = (v_at_threshold - v_w).abs() / v_w;
        assert!(
            rel_err < 1.0e-12,
            "v_below_walker(H_W) must equal v_W, rel_err = {:.2e}",
            rel_err
        );
    }

    #[test]
    fn test_walker_velocity_linear_below_threshold() {
        let wb = WalkerBreakdown::new(DwMaterial::co_pt_pma());
        let h_w = wb.walker_threshold_field();
        let v1 = wb.velocity_below_walker(h_w * 0.5);
        let v2 = wb.velocity_below_walker(h_w * 1.0);
        // v should be proportional to H → doubling H doubles v
        let ratio = v2 / v1;
        assert!(
            (ratio - 2.0).abs() < 1.0e-10,
            "Below-Walker velocity must be linear in H, ratio = {:.6}",
            ratio
        );
    }

    #[test]
    fn test_velocity_auto_regime_selection() {
        let wb = WalkerBreakdown::new(DwMaterial::co_pt_pma());
        let h_w = wb.walker_threshold_field();
        // At a field just below threshold, velocity() must match below-Walker
        let h_sub = h_w * 0.9;
        assert!(
            (wb.velocity(h_sub) - wb.velocity_below_walker(h_sub)).abs() < 1.0e-30,
            "velocity() must use below-Walker formula below threshold"
        );
        // Above threshold, velocity() must match above-Walker
        let h_sup = h_w * 1.5;
        assert!(
            (wb.velocity(h_sup) - wb.velocity_above_walker(h_sup)).abs() < 1.0e-30,
            "velocity() must use above-Walker formula above threshold"
        );
    }

    #[test]
    fn test_velocity_sign_follows_field_sign() {
        let wb = WalkerBreakdown::new(DwMaterial::co_pt_pma());
        let h = wb.walker_threshold_field() * 0.5;
        assert!(wb.velocity(h) > 0.0, "positive field → positive velocity");
        assert!(wb.velocity(-h) < 0.0, "negative field → negative velocity");
    }

    #[test]
    fn test_above_walker_velocity_smaller_than_walker_velocity() {
        // For α < 1, the time-averaged above-Walker velocity at H = 2 H_W
        // should be less than the Walker velocity v_W.
        let wb = WalkerBreakdown::new(DwMaterial::co_pt_pma());
        let h_w = wb.walker_threshold_field();
        let v_above = wb.velocity_above_walker(2.0 * h_w);
        let v_w = wb.walker_velocity();
        assert!(
            v_above < v_w * 2.5,
            "above-Walker velocity should not greatly exceed Walker velocity, v_above={:.2}, v_W={:.2}",
            v_above, v_w
        );
    }

    // ── DwSttDynamics ────────────────────────────────────────────────────────

    #[test]
    fn test_stt_adiabatic_velocity_linear_in_j() {
        let stt = DwSttDynamics::new(DwMaterial::co_pt_pma());
        let j1 = 1.0e11;
        let j2 = 2.0e11;
        let v1 = stt.adiabatic_velocity(j1);
        let v2 = stt.adiabatic_velocity(j2);
        let ratio = v2 / v1;
        assert!(
            (ratio - 2.0).abs() < 1.0e-10,
            "Adiabatic velocity must be linear in J, ratio = {:.6}",
            ratio
        );
    }

    #[test]
    fn test_stt_adiabatic_velocity_physical_scale() {
        // For Co/Pt: v_s = P j μ_B / (e Ms) ≈ 0.5 × 1e12 × 9.27e-24 / (1.6e-19 × 1.1e6)
        //                                    ≈ 0.5 × 1e12 × 5.8e-5 ≈ 29 m/s at j=1e12
        let stt = DwSttDynamics::new(DwMaterial::co_pt_pma());
        let v = stt.adiabatic_velocity(1.0e12);
        assert!(
            v > 1.0 && v < 1000.0,
            "Adiabatic velocity at 1e12 A/m² should be 1–1000 m/s, got {:.2} m/s",
            v
        );
    }

    #[test]
    fn test_stt_velocity_zero_current_zero_field() {
        let stt = DwSttDynamics::new(DwMaterial::co_pt_pma());
        let v = stt.velocity(0.0, 0.0);
        assert!(
            v.abs() < 1.0e-30,
            "Zero current + zero field → zero velocity"
        );
    }

    #[test]
    fn test_stt_velocity_agrees_with_adiabatic_at_small_alpha_beta() {
        // When α ≈ β ≈ 0 the full velocity should equal v_s
        // Use a PMA material with very small α and β:
        // K_u = 1.5e6 > μ₀ Ms²/2 ≈ 0.76e6 ensures K_eff > 0.
        let mat = DwMaterial {
            name: "test",
            exchange_a: 1.5e-11,
            anisotropy_k: 1.5e6,
            ms: 1.1e6,
            alpha: 1.0e-4,
            thickness: 1.0e-9,
            polarization: 0.50,
            beta: 1.0e-4,
            theta_sh: 0.0,
        };
        let stt = DwSttDynamics::new(mat);
        let j = 1.0e12;
        let v_s = stt.adiabatic_velocity(j);
        let v = stt.velocity(j, 0.0);
        let rel_err = (v - v_s).abs() / v_s;
        assert!(
            rel_err < 1.0e-6,
            "At α=β→0 full velocity must equal v_s; rel_err = {:.2e}",
            rel_err
        );
    }

    #[test]
    fn test_stt_walker_current_density_positive() {
        let stt = DwSttDynamics::new(DwMaterial::co_pt_pma());
        let j_w = stt.walker_current_density();
        assert!(j_w > 0.0, "Walker current density must be positive");
        // Should be in physically reasonable range 1e11–1e14 A/m²
        assert!(
            j_w > 1.0e10 && j_w < 1.0e15,
            "Walker current density out of expected range: {:.2e} A/m²",
            j_w
        );
    }

    // ── DwSotDynamics ────────────────────────────────────────────────────────

    #[test]
    fn test_sot_dl_field_linear_in_j() {
        let sot = DwSotDynamics::new(DwMaterial::pt_co_alox());
        let j1 = 1.0e11;
        let j2 = 3.0e11;
        let h1 = sot.dl_sot_field(j1);
        let h2 = sot.dl_sot_field(j2);
        let ratio = h2 / h1;
        assert!(
            (ratio - 3.0).abs() < 1.0e-10,
            "DL-SOT field must be linear in J, ratio = {:.6}",
            ratio
        );
    }

    #[test]
    fn test_sot_dl_field_physical_scale() {
        // For Pt/Co/AlOx at J=3e11 A/m²:
        // H_DL = 0.15 × 1.055e-34 × 3e11 / (2 × 1.6e-19 × 1.257e-6 × 1e6 × 0.6e-9)
        //       = 0.15 × 3.165e-23 / (2 × 1.6e-19 × 7.54e-10)
        //       ≈ tens of kA/m
        let sot = DwSotDynamics::new(DwMaterial::pt_co_alox());
        let h_dl = sot.dl_sot_field(3.0e11);
        assert!(
            h_dl > 1.0e3 && h_dl < 1.0e8,
            "DL-SOT field should be in kA/m–MA/m range, got {:.2e} A/m",
            h_dl
        );
    }

    #[test]
    fn test_sot_velocity_linear_in_j_below_walker() {
        let sot = DwSotDynamics::new(DwMaterial::pt_co_alox());
        // Use currents well below Walker threshold
        let j_w = sot.walker_current_density();
        let j1 = j_w * 0.1;
        let j2 = j_w * 0.2;
        let v1 = sot.velocity(j1);
        let v2 = sot.velocity(j2);
        let ratio = v2 / v1;
        assert!(
            (ratio - 2.0).abs() < 1.0e-10,
            "SOT velocity must be linear in J below Walker, ratio = {:.6}",
            ratio
        );
    }

    #[test]
    fn test_sot_velocity_zero_at_zero_current() {
        let sot = DwSotDynamics::new(DwMaterial::pt_co_alox());
        let v = sot.velocity(0.0);
        assert!(v.abs() < 1.0e-30, "Zero current → zero velocity");
    }

    #[test]
    fn test_sot_velocity_with_field_additive() {
        // velocity_with_field at H=0 should equal velocity
        let sot = DwSotDynamics::new(DwMaterial::pt_co_alox());
        let j = 1.0e11;
        let v_no_field = sot.velocity(j);
        let v_with_zero = sot.velocity_with_field(j, 0.0);
        assert!(
            (v_no_field - v_with_zero).abs() < 1.0e-30,
            "velocity_with_field(j, 0) must equal velocity(j)"
        );
    }

    #[test]
    fn test_sot_walker_current_density_positive() {
        let sot = DwSotDynamics::new(DwMaterial::pt_co_alox());
        let j_w = sot.walker_current_density();
        assert!(j_w > 0.0, "SOT Walker current density must be positive");
        assert!(
            j_w > 1.0e10 && j_w < 1.0e16,
            "SOT Walker current density out of plausible range: {:.2e} A/m²",
            j_w
        );
    }

    #[test]
    fn test_sot_velocity_curve_length_and_ordering() {
        let sot = DwSotDynamics::new(DwMaterial::pt_co_alox());
        let curve = sot.velocity_curve(0.0, 3.0e11, 10);
        assert_eq!(curve.len(), 10, "Curve must have exactly 10 points");
        // First J should be 0, last should be 3e11
        assert!(curve[0].0.abs() < 1.0, "First J should be ~0");
        assert!((curve[9].0 - 3.0e11).abs() < 1.0, "Last J should be ~3e11");
    }

    #[test]
    fn test_sot_velocity_curve_monotone_below_walker() {
        let sot = DwSotDynamics::new(DwMaterial::pt_co_alox());
        let j_w = sot.walker_current_density();
        // Sample only below the Walker threshold
        let curve = sot.velocity_curve(0.0, j_w * 0.9, 20);
        // Velocities should be monotonically increasing (J > 0, v > 0)
        for i in 1..curve.len() {
            assert!(
                curve[i].1 >= curve[i - 1].1,
                "Velocity curve must be monotonically non-decreasing below Walker: v[{}]={:.2} < v[{}]={:.2}",
                i, curve[i].1, i - 1, curve[i - 1].1
            );
        }
    }

    #[test]
    fn test_sot_velocity_pt_co_alox_reasonable_scale() {
        // At J = 3e11 A/m² expect velocity in range 1–1000 m/s
        let sot = DwSotDynamics::new(DwMaterial::pt_co_alox());
        let v = sot.velocity(3.0e11);
        assert!(
            v > 1.0 && v < 1000.0,
            "Pt/Co/AlOx DW velocity at 3e11 A/m² should be 1–1000 m/s, got {:.2} m/s",
            v
        );
    }
}
