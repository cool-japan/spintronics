//! Time-domain LLG dynamics driven by periodic strain (SAW / AC piezoelectric)
//!
//! The rest of `mech/` ships two complementary but purely *static* or
//! *steady-state analytical* descriptions of strain-mediated magnetism
//! control:
//!
//! - [`super::saw`]: a closed-form Lorentzian acoustic-FMR model
//!   ([`SawMagnetoacoustics`]) built on the instantaneous periodic strain
//!   field of a Rayleigh surface acoustic wave,
//!   [`SawSource::strain_at_point`].
//! - [`super::magnetoelastic`]: static magnetoelastic energetics — the
//!   Villari effect ([`villari_effective_field_vector`]), stress-induced
//!   anisotropy ([`stress_induced_anisotropy`]), and piezoelectric strain
//!   generation ([`PiezoelectricSubstrate::strain_from_voltage`]).
//!
//! Neither of those subsystems actually evolves a magnetization vector in
//! time. This module closes that gap: it couples an oscillating strain
//! source to the existing magnetoelastic field machinery and steps a
//! macrospin magnetization forward under the resulting time-dependent
//! effective field using the crate's existing adaptive Runge-Kutta
//! integrator ([`DormandPrince45`], reused unmodified from
//! [`crate::dynamics::integrators`]).
//!
//! # Pipeline
//!
//! At every integration substep, for a magnetization `m` and time `t`:
//!
//! 1. Sample the instantaneous strain tensor from a [`StrainDrive`]
//!    ([`SawStrainDrive`] wraps [`SawSource::strain_at_point`];
//!    [`PiezoAcStrainDrive`] wraps [`PiezoelectricSubstrate::strain_from_voltage`]
//!    under a sinusoidal drive voltage).
//! 2. Convert strain to stress via the uniaxial Hooke's-law approximation
//!    `sigma = Y * epsilon` (Young's modulus `Y` from the
//!    [`MagnetoelasticMaterial`] supplied for the Villari coupling), and
//!    rotate the resulting stress tensor into the magnetization's reference
//!    frame (see "Coupling geometry" below).
//! 3. Convert stress to an effective field via
//!    [`villari_effective_field_vector`] (a diagnostic, non-dynamical
//!    induced-anisotropy constant is also exposed via
//!    [`stress_induced_anisotropy`], see
//!    [`StrainDrivenLlgDriver::induced_anisotropy_constant_j_per_m3`]).
//! 4. Add the strain field to the static bias + anisotropy field and
//!    evaluate the LLG right-hand side via
//!    [`crate::dynamics::llg::calc_dm_dt`] (the exact same function used by
//!    [`crate::dynamics::llg::LlgSolver`]).
//! 5. Advance one adaptive-RK step and renormalize `|m| = 1` (the same
//!    constraint-projection convention already used throughout
//!    [`crate::dynamics::llg::LlgSolver`]'s own stepping methods).
//!
//! # Coupling geometry: why a 45-degree tilt
//!
//! A strain tensor that is purely diagonal in the magnetization's own frame
//! (e.g. `epsilon_xx`, `epsilon_zz` with the bias field and easy axis along
//! z) couples to the magnetoelastic energy only through terms like
//! `epsilon_zz * m_z^2`. Expanding `m = z_hat + delta_m` shows this produces
//! **no linear-order torque** on the transverse deviation `delta_m` — it
//! only renormalizes the longitudinal (parallel-pumping-type) field, which
//! cannot resonantly drive uniform-mode FMR at `omega = omega_FMR` (only at
//! `2 omega_FMR`). Real acoustic-FMR experiments avoid this by orienting the
//! magnetization at an oblique angle (canonically 45 degrees) to the SAW
//! propagation direction (Weiler et al., PRL 106, 117601 (2011); Thevenard
//! et al., PRB 93, 134430 (2016)), which mixes the diagonal strain
//! components into a genuine shear-like, magnetization-independent
//! transverse drive term at leading order.
//!
//! Rather than tilting the (already-established) bias/anisotropy geometry,
//! [`StrainDrivenLlgDriver::coupling_angle_rad`] rotates the *strain tensor*
//! by this angle about the y-axis before computing the Villari field. This
//! is mathematically equivalent to tilting the magnetization frame, but
//! keeps `m` at the north pole in equilibrium, so
//! [`SawMagnetoacoustics::acoustic_fmr_condition_h_ext`] (which assumes
//! out-of-plane saturation) remains exactly the correct resonance condition
//! for the driven macrospin simulated here.
//!
//! # Units: A/m versus Tesla
//!
//! [`villari_effective_field_vector`] and the anisotropy/bias fields used
//! here follow the standard micromagnetic convention of an H-field in A/m
//! (`H_eff = -1/(mu_0 Ms) dE/dm`). [`crate::dynamics::llg::calc_dm_dt`],
//! however, expects its field argument in Tesla (`GAMMA` is the electron
//! gyromagnetic ratio in rad/(s*T)) — exactly the same convention already
//! used by [`crate::mech::saw::MagnetoelasticMaterial::fmr_frequency_hz`],
//! which multiplies its A/m fields by `MU_0` before combining them with
//! `GAMMA`. This module performs the same conversion once, in
//! `StrainDrivenLlgDriver::total_h_eff_tesla`, immediately before calling
//! into the LLG solver.
//!
//! # References
//! - Weiler et al., PRL **106**, 117601 (2011) — acoustic FMR
//! - Thevenard et al., PRB **93**, 134430 (2016) — SAW-driven precession

use std::f64::consts::{FRAC_PI_4, PI};

use crate::constants::{GAMMA, MU_0};
use crate::dynamics::llg::calc_dm_dt;
use crate::dynamics::{DormandPrince45, Integrator};
use crate::error::{invalid_param, numerical_error, Result};
use crate::mech::magnetoelastic::{
    stress_induced_anisotropy, villari_effective_field_vector, MagnetoelasticMaterial,
    PiezoelectricSubstrate, StrainTensor, StressTensor,
};
use crate::mech::saw::{SawMagnetoacoustics, SawSource};
use crate::vector3::Vector3;

/// The canonical acoustic-FMR coupling geometry \[rad\]: the angle between
/// the strain principal axis and the magnetic bias-field/easy-axis
/// direction that maximizes the linear-order magnetoelastic torque.
///
/// See the module-level documentation ("Coupling geometry") for the
/// derivation. Collinear geometry (`0` or `PI/2`) gives zero linear-order
/// coupling; `PI/4` is optimal.
pub const OPTIMAL_COUPLING_ANGLE_RAD: f64 = FRAC_PI_4;

// =============================================================================
// StrainDrive abstraction
// =============================================================================

/// A source of time-periodic strain that can drive magnetization dynamics.
///
/// Implemented by [`SawStrainDrive`] (surface acoustic wave) and
/// [`PiezoAcStrainDrive`] (AC-driven piezoelectric substrate), so that
/// [`StrainDrivenLlgDriver`] can couple either "SAW-driven spin dynamics" or
/// "piezoelectric control of magnetism" into the same time-domain LLG
/// integration without duplicating the coupling/integration logic.
pub trait StrainDrive {
    /// Instantaneous strain tensor at time `t` \[s\] (dimensionless strain
    /// components).
    fn strain_tensor_at(&self, t: f64) -> StrainTensor;

    /// Drive frequency \[Hz\] used to size the simulation window (one "drive
    /// period" is `1 / drive_frequency_hz()`).
    fn drive_frequency_hz(&self) -> f64;
}

/// Wraps a [`SawSource`] as a [`StrainDrive`], sampling the Rayleigh-wave
/// strain field at a fixed observation point `(x_probe_m, z_probe_m)`.
#[derive(Debug, Clone)]
pub struct SawStrainDrive {
    /// The underlying SAW source (piezoelectric substrate + drive
    /// parameters).
    pub saw: SawSource,
    /// Position along the SAW propagation direction at which strain is
    /// sampled \[m\].
    pub x_probe_m: f64,
    /// Depth below the surface at which strain is sampled \[m\] (0 =
    /// surface, increasing into the bulk).
    pub z_probe_m: f64,
}

impl SawStrainDrive {
    /// Construct a new SAW strain drive sampled at `(x_probe_m, z_probe_m)`.
    pub fn new(saw: SawSource, x_probe_m: f64, z_probe_m: f64) -> Self {
        Self {
            saw,
            x_probe_m,
            z_probe_m,
        }
    }
}

impl StrainDrive for SawStrainDrive {
    fn strain_tensor_at(&self, t: f64) -> StrainTensor {
        let eps_xx = self.saw.strain_at_point(self.x_probe_m, self.z_probe_m, t);
        // Poisson coupling of the Rayleigh-wave out-of-plane normal strain,
        // per the SawSource module documentation: eps_zz ~= -nu * eps_xx.
        let eps_zz = -self.saw.substrate.poisson_ratio * eps_xx;
        StrainTensor::new([[eps_xx, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, eps_zz]])
    }

    fn drive_frequency_hz(&self) -> f64 {
        self.saw.frequency_hz
    }
}

/// AC-driven piezoelectric strain source: a sinusoidal voltage across a
/// piezoelectric substrate produces a time-periodic strain through the
/// direct piezoelectric effect, `epsilon(t) = d * E(t)` with
/// `E(t) = (V0 sin(2 pi f t)) / thickness`.
///
/// This complements the static/quasi-static piezoelectric control already
/// modeled by [`super::straintronics::StraintronicDevice`] and
/// [`PiezoelectricSubstrate::strain_from_voltage`] by giving those same
/// strain-generation formulas a genuine time dependence suitable for driving
/// LLG dynamics.
#[derive(Debug, Clone)]
pub struct PiezoAcStrainDrive {
    /// Piezoelectric substrate (PZT, PMN-PT, `BaTiO3`, ...).
    pub substrate: PiezoelectricSubstrate,
    /// Substrate thickness \[m\] over which the drive voltage is applied.
    pub substrate_thickness_m: f64,
    /// Peak drive voltage amplitude \[V\].
    pub voltage_amplitude_v: f64,
    /// AC drive frequency \[Hz\].
    pub frequency_hz: f64,
}

impl PiezoAcStrainDrive {
    /// Construct a new AC piezoelectric strain drive.
    ///
    /// # Errors
    /// Returns an error if `substrate_thickness_m` or `frequency_hz` is not
    /// positive, or if the peak field `voltage_amplitude_v /
    /// substrate_thickness_m` would exceed the substrate's rated maximum
    /// field.
    pub fn new(
        substrate: PiezoelectricSubstrate,
        substrate_thickness_m: f64,
        voltage_amplitude_v: f64,
        frequency_hz: f64,
    ) -> Result<Self> {
        if substrate_thickness_m <= 0.0 {
            return Err(invalid_param(
                "substrate_thickness_m",
                "substrate thickness must be positive",
            ));
        }
        if frequency_hz <= 0.0 {
            return Err(invalid_param(
                "frequency_hz",
                "AC drive frequency must be positive",
            ));
        }
        let peak_field = voltage_amplitude_v.abs() / substrate_thickness_m;
        if peak_field > substrate.max_field {
            return Err(invalid_param(
                "voltage_amplitude_v",
                "peak field (voltage_amplitude_v / substrate_thickness_m) exceeds the substrate's rated maximum field",
            ));
        }
        Ok(Self {
            substrate,
            substrate_thickness_m,
            voltage_amplitude_v,
            frequency_hz,
        })
    }
}

impl StrainDrive for PiezoAcStrainDrive {
    fn strain_tensor_at(&self, t: f64) -> StrainTensor {
        let voltage_t = self.voltage_amplitude_v * (2.0 * PI * self.frequency_hz * t).sin();
        self.substrate
            .strain_from_voltage(voltage_t, self.substrate_thickness_m)
            .unwrap_or_else(|_| StrainTensor::zero())
    }

    fn drive_frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
}

// =============================================================================
// Trajectory output types
// =============================================================================

/// A single sampled point along a strain-driven magnetization trajectory.
#[derive(Debug, Clone, Copy)]
pub struct StrainDrivenStepSample {
    /// Elapsed simulation time \[s\].
    pub time_s: f64,
    /// Normalized magnetization at this time.
    pub magnetization: Vector3<f64>,
    /// Cone angle away from the +z easy axis / bias-field direction \[rad\].
    pub cone_angle_rad: f64,
}

/// Result of driving a macrospin magnetization under a time-periodic strain
/// field for a finite simulation window.
#[derive(Debug, Clone)]
pub struct StrainDrivenTrajectory {
    /// Full sampled trajectory: one entry per accepted time step, plus the
    /// initial condition.
    pub samples: Vec<StrainDrivenStepSample>,
    /// The largest cone angle reached anywhere in the trajectory \[rad\].
    pub max_cone_angle_rad: f64,
    /// The largest observed deviation of the *raw* (pre-renormalization)
    /// integrator output from unit magnitude, `||m_raw| - 1|`, taken over
    /// all steps.
    ///
    /// This is the sanity-check quantity for norm conservation of the
    /// underlying LLG integration itself, independent of the constraint
    /// projection ([`Vector3::normalize`]) applied between steps.
    pub max_raw_norm_deviation: f64,
    /// Mean Gilbert-damping power dissipated per unit film area \[W/m²\],
    /// averaged over the final full drive period of the simulation window.
    ///
    /// At quasi-steady state, energy conservation requires this to equal
    /// the power absorbed from the strain drive, making it directly
    /// comparable (in trend and order of magnitude) to
    /// [`SawMagnetoacoustics::resonant_absorption`].
    pub mean_power_areal_density_w_per_m2: f64,
}

// =============================================================================
// StrainDrivenLlgDriver
// =============================================================================

/// Couples a time-periodic strain source to the LLG equation, evolving a
/// macrospin magnetization under the resulting magnetoelastic effective
/// field.
///
/// See the module-level documentation for the full physical pipeline.
#[derive(Debug, Clone)]
pub struct StrainDrivenLlgDriver<S: StrainDrive> {
    /// The time-periodic strain source (SAW or AC piezoelectric).
    pub strain_drive: S,
    /// Magnetoelastic material supplying the magnetostriction constant
    /// `lambda_s` and Young's modulus used for the Hooke's-law
    /// strain-to-stress conversion and the Villari-effect field.
    pub villari_material: MagnetoelasticMaterial,
    /// Gyromagnetic ratio \[rad/(s*T)\].
    pub gamma: f64,
    /// Gilbert damping constant (dimensionless, `> 0`).
    pub alpha: f64,
    /// Uniaxial anisotropy constant \[J/m³\] with easy axis along +z.
    pub anisotropy_k: f64,
    /// Static external bias field along +z \[A/m\] (non-negative,
    /// out-of-plane-saturation convention, matching
    /// [`SawMagnetoacoustics`]).
    pub h_ext_a_per_m: f64,
    /// Ferromagnetic film thickness \[m\], used to convert volumetric power
    /// density to areal power density.
    pub film_thickness_m: f64,
    /// Angle \[rad\] between the strain principal axis and the
    /// magnetization's bias-field/easy-axis direction (see module docs,
    /// "Coupling geometry"). [`OPTIMAL_COUPLING_ANGLE_RAD`] maximizes the
    /// linear-order magnetoelastic torque.
    pub coupling_angle_rad: f64,
}

impl<S: StrainDrive> StrainDrivenLlgDriver<S> {
    /// Construct a new strain-driven LLG driver.
    ///
    /// # Errors
    /// Returns an error if `villari_material.ms` or `gamma` or
    /// `film_thickness_m` is not positive, if `alpha` is not strictly
    /// positive (a driven-dissipative steady state requires finite
    /// damping), if `h_ext_a_per_m` is negative, or if `coupling_angle_rad`
    /// is not finite.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        strain_drive: S,
        villari_material: MagnetoelasticMaterial,
        gamma: f64,
        alpha: f64,
        anisotropy_k: f64,
        h_ext_a_per_m: f64,
        film_thickness_m: f64,
        coupling_angle_rad: f64,
    ) -> Result<Self> {
        if villari_material.ms <= 0.0 {
            return Err(invalid_param(
                "villari_material.ms",
                "saturation magnetization must be positive",
            ));
        }
        if gamma <= 0.0 {
            return Err(invalid_param(
                "gamma",
                "gyromagnetic ratio must be positive",
            ));
        }
        if alpha.is_nan() || alpha <= 0.0 {
            return Err(invalid_param(
                "alpha",
                "Gilbert damping must be strictly positive for a driven-dissipative response",
            ));
        }
        if h_ext_a_per_m < 0.0 {
            return Err(invalid_param(
                "h_ext_a_per_m",
                "external bias field must be non-negative (out-of-plane-saturation convention)",
            ));
        }
        if film_thickness_m <= 0.0 {
            return Err(invalid_param(
                "film_thickness_m",
                "film thickness must be positive",
            ));
        }
        if !coupling_angle_rad.is_finite() {
            return Err(invalid_param("coupling_angle_rad", "must be finite"));
        }
        Ok(Self {
            strain_drive,
            villari_material,
            gamma,
            alpha,
            anisotropy_k,
            h_ext_a_per_m,
            film_thickness_m,
            coupling_angle_rad,
        })
    }

    /// Rotate a diagonal `(x, z)`-plane stress tensor (computed in the
    /// strain source's own frame) by [`Self::coupling_angle_rad`] about the
    /// y-axis into the magnetization's `+z`-aligned reference frame.
    ///
    /// Standard 2D symmetric-tensor rotation:
    /// `sigma_xx' = sigma_xx cos^2(phi) + sigma_zz sin^2(phi)`,
    /// `sigma_zz' = sigma_xx sin^2(phi) + sigma_zz cos^2(phi)`,
    /// `sigma_xz' = (sigma_xx - sigma_zz) sin(phi) cos(phi)`.
    fn tilted_stress_tensor(&self, eps_xx: f64, eps_zz: f64) -> StressTensor {
        let youngs_modulus_pa = self.villari_material.youngs_modulus;
        let sigma_xx_saw = youngs_modulus_pa * eps_xx;
        let sigma_zz_saw = youngs_modulus_pa * eps_zz;

        let cos_phi = self.coupling_angle_rad.cos();
        let sin_phi = self.coupling_angle_rad.sin();

        let sigma_xx = sigma_xx_saw * cos_phi * cos_phi + sigma_zz_saw * sin_phi * sin_phi;
        let sigma_zz = sigma_xx_saw * sin_phi * sin_phi + sigma_zz_saw * cos_phi * cos_phi;
        let sigma_xz = (sigma_xx_saw - sigma_zz_saw) * sin_phi * cos_phi;

        StressTensor::new([
            [sigma_xx, 0.0, sigma_xz],
            [0.0, 0.0, 0.0],
            [sigma_xz, 0.0, sigma_zz],
        ])
    }

    /// Total effective field in Tesla (bias + uniaxial anisotropy + strain-
    /// induced Villari field), ready to pass to
    /// [`crate::dynamics::llg::calc_dm_dt`].
    ///
    /// All three contributions are assembled in A/m (the standard
    /// micromagnetic H-field convention used by
    /// [`villari_effective_field_vector`]) and converted to Tesla by a
    /// single trailing `* MU_0`, matching the convention used by
    /// [`crate::mech::saw::MagnetoelasticMaterial::fmr_frequency_hz`] (see
    /// module docs, "Units").
    fn total_h_eff_tesla(&self, m: Vector3<f64>, t: f64) -> Vector3<f64> {
        let strain = self.strain_drive.strain_tensor_at(t);
        let (eps_xx, _eps_yy, eps_zz) = strain.diagonal();
        let stress = self.tilted_stress_tensor(eps_xx, eps_zz);

        let h_villari_a_per_m = villari_effective_field_vector(&self.villari_material, &stress, &m)
            .unwrap_or_else(|_| Vector3::zero());

        let h_bias_a_per_m = Vector3::unit_z() * self.h_ext_a_per_m;
        let h_k_a_per_m = 2.0 * self.anisotropy_k / (MU_0 * self.villari_material.ms);
        let h_aniso_a_per_m = Vector3::unit_z() * (h_k_a_per_m * m.z);

        (h_bias_a_per_m + h_aniso_a_per_m + h_villari_a_per_m) * MU_0
    }

    /// LLG time derivative `dm/dt` at magnetization `m` and time `t`,
    /// reusing [`calc_dm_dt`] with the total strain-driven effective field.
    fn rhs(&self, m: Vector3<f64>, t: f64) -> Vector3<f64> {
        let h_eff_tesla = self.total_h_eff_tesla(m, t);
        calc_dm_dt(m, h_eff_tesla, self.gamma, self.alpha)
    }

    /// Instantaneous Gilbert-damping power dissipation density \[W/m³\].
    ///
    /// Derived from the standard LLG energy-dissipation identity: with
    /// `dE/dt = -Ms (H . dm/dt)` (Tesla convention, `E = -Ms (m . B)`) and
    /// `H . dm/dt = gamma alpha / (1 + alpha^2) |m x H|^2` (a direct
    /// consequence of `calc_dm_dt`'s functional form), the power *leaving*
    /// the magnetic subsystem is:
    ///
    /// `P = Ms * gamma * alpha / (1 + alpha^2) * |m x H_eff|^2`
    ///
    /// which is manifestly non-negative (damping always dissipates energy).
    fn dissipation_power_density_w_per_m3(
        &self,
        m: Vector3<f64>,
        h_eff_tesla: Vector3<f64>,
    ) -> f64 {
        let cross = m.cross(&h_eff_tesla);
        self.villari_material.ms * self.gamma * self.alpha / (1.0 + self.alpha * self.alpha)
            * cross.dot(&cross)
    }

    /// Instantaneous stress-induced uniaxial anisotropy constant
    /// `K_sigma(t)` \[J/m³\], evaluated from the dominant (un-rotated,
    /// strain-source-frame) normal stress via [`stress_induced_anisotropy`].
    ///
    /// This is a diagnostic quantity only, exposed for inspection/logging.
    /// It is *not* separately added to the dynamical effective field, which
    /// already includes the equivalent physics through
    /// [`villari_effective_field_vector`] — adding both would double-count
    /// the same magnetoelastic coupling.
    pub fn induced_anisotropy_constant_j_per_m3(&self, t: f64) -> f64 {
        let strain = self.strain_drive.strain_tensor_at(t);
        let (eps_xx, _eps_yy, _eps_zz) = strain.diagonal();
        let sigma_xx = self.villari_material.youngs_modulus * eps_xx;
        stress_induced_anisotropy(self.villari_material.lambda_s, sigma_xx)
    }

    /// Evolve the macrospin magnetization under the strain drive for
    /// `n_periods` full drive cycles (one drive period is
    /// `1 / strain_drive.drive_frequency_hz()`), using a fixed time step
    /// `dt_s` and the [`DormandPrince45`] embedded Runge-Kutta integrator.
    ///
    /// A fresh [`DormandPrince45`] instance is used for every step (rather
    /// than reusing one across the whole run) so that its FSAL
    /// first-stage-reuse optimization never operates on a state that has
    /// been renormalized since the previous step's evaluation — i.e. every
    /// stage is always evaluated at the true current (unit-norm) state.
    ///
    /// # Errors
    /// Returns an error if `dt_s` is not positive, if `n_periods` is zero,
    /// if `m0` has near-zero magnitude, or if the underlying integrator
    /// reports a numerical failure (e.g. NaN).
    pub fn evolve(
        &self,
        m0: Vector3<f64>,
        dt_s: f64,
        n_periods: u32,
    ) -> Result<StrainDrivenTrajectory> {
        if dt_s <= 0.0 {
            return Err(invalid_param("dt_s", "time step must be positive"));
        }
        if n_periods == 0 {
            return Err(invalid_param(
                "n_periods",
                "must evolve for at least one full drive period",
            ));
        }
        if m0.magnitude() < 1e-30 {
            return Err(invalid_param(
                "m0",
                "initial magnetization must have non-zero magnitude",
            ));
        }

        let period_s = 1.0 / self.strain_drive.drive_frequency_hz();
        let total_time_s = period_s * f64::from(n_periods);
        let n_steps = (total_time_s / dt_s).ceil().max(1.0) as usize;
        let last_period_start_s = total_time_s - period_s;

        let mut m = m0.normalize();
        let mut t = 0.0_f64;

        let mut samples = Vec::with_capacity(n_steps + 1);
        let mut max_cone_angle_rad = cone_angle_rad(m);
        let mut max_raw_norm_deviation = 0.0_f64;
        let mut dissipated_energy_areal_j_per_m2 = 0.0_f64;

        samples.push(StrainDrivenStepSample {
            time_s: t,
            magnetization: m,
            cone_angle_rad: max_cone_angle_rad,
        });

        let rhs = |state: &[Vector3<f64>], time: f64| -> Vec<Vector3<f64>> {
            let m_state = state.first().copied().unwrap_or_else(Vector3::zero);
            vec![self.rhs(m_state, time)]
        };

        for _ in 0..n_steps {
            let state = [m];
            // Fresh integrator per step: see doc comment above for why.
            let mut integrator = DormandPrince45::new();
            let output = integrator.step(&state, t, dt_s, &rhs)?;
            let m_raw = output
                .new_state
                .first()
                .copied()
                .ok_or_else(|| numerical_error("integrator returned an empty state vector"))?;

            let raw_dev = (m_raw.magnitude() - 1.0).abs();
            if raw_dev > max_raw_norm_deviation {
                max_raw_norm_deviation = raw_dev;
            }

            m = m_raw.normalize();
            t += dt_s;

            let cone = cone_angle_rad(m);
            if cone > max_cone_angle_rad {
                max_cone_angle_rad = cone;
            }

            if t > last_period_start_s {
                let h_eff_tesla = self.total_h_eff_tesla(m, t);
                let p_vol = self.dissipation_power_density_w_per_m3(m, h_eff_tesla);
                dissipated_energy_areal_j_per_m2 += p_vol * self.film_thickness_m * dt_s;
            }

            samples.push(StrainDrivenStepSample {
                time_s: t,
                magnetization: m,
                cone_angle_rad: cone,
            });
        }

        let mean_power_areal_density_w_per_m2 = dissipated_energy_areal_j_per_m2 / period_s;

        Ok(StrainDrivenTrajectory {
            samples,
            max_cone_angle_rad,
            max_raw_norm_deviation,
            mean_power_areal_density_w_per_m2,
        })
    }
}

impl StrainDrivenLlgDriver<SawStrainDrive> {
    /// Build a driver directly from a shipped [`SawMagnetoacoustics`]
    /// coupling calculator, reusing its FMR bookkeeping (`alpha`,
    /// `anisotropy_k`, bias field, film thickness) for the time-domain
    /// dynamics.
    ///
    /// `villari_material` supplies the magnetostriction constant
    /// `lambda_s` and Young's modulus used to convert the instantaneous SAW
    /// strain into a stress tensor and, via
    /// [`villari_effective_field_vector`], into a magnetoelastic effective
    /// field. It need not be numerically identical to `device.material` (a
    /// different, B1/B2-parameterised description of the same physical
    /// film); it should describe the same or a comparable magnetic film for
    /// the bridge to be physically meaningful.
    ///
    /// # Errors
    /// Propagates validation errors from [`StrainDrivenLlgDriver::new`].
    pub fn from_saw_magnetoacoustics(
        device: &SawMagnetoacoustics,
        villari_material: MagnetoelasticMaterial,
        x_probe_m: f64,
        z_probe_m: f64,
        coupling_angle_rad: f64,
    ) -> Result<Self> {
        let strain_drive = SawStrainDrive::new(device.saw.clone(), x_probe_m, z_probe_m);
        Self::new(
            strain_drive,
            villari_material,
            GAMMA,
            device.material.alpha,
            device.material.anisotropy_k,
            device.h_ext,
            device.film_thickness_m,
            coupling_angle_rad,
        )
    }
}

/// Cone angle (angle from +z) of a magnetization direction \[rad\], in
/// `[0, pi]`.
#[inline]
fn cone_angle_rad(m: Vector3<f64>) -> f64 {
    m.z.clamp(-1.0, 1.0).acos()
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mech::saw::{MagnetoelasticMaterial as SawMaterial, PiezoSubstrate};

    fn nickel_saw_device(saw: SawSource, h_ext: f64) -> SawMagnetoacoustics {
        SawMagnetoacoustics::new(saw, SawMaterial::nickel(), 20.0e-9, h_ext)
            .expect("valid SawMagnetoacoustics parameters")
    }

    // -------------------------------------------------------------------
    // StrainDrive implementations
    // -------------------------------------------------------------------

    #[test]
    fn test_saw_strain_drive_matches_source_directly() {
        let saw = SawSource::linbo3_1ghz();
        let drive = SawStrainDrive::new(saw.clone(), 1.0e-6, 0.0);
        let (eps_xx, eps_yy, eps_zz) = drive.strain_tensor_at(1.0e-10).diagonal();
        let expected_xx = saw.strain_at_point(1.0e-6, 0.0, 1.0e-10);
        assert!(
            (eps_xx - expected_xx).abs() < 1.0e-20,
            "SawStrainDrive eps_xx should match SawSource::strain_at_point exactly"
        );
        assert_eq!(eps_yy, 0.0, "SAW model has no eps_yy component");
        let expected_zz = -saw.substrate.poisson_ratio * expected_xx;
        assert!(
            (eps_zz - expected_zz).abs() < 1.0e-20,
            "SawStrainDrive eps_zz should follow the Poisson relation"
        );
    }

    #[test]
    fn test_saw_strain_drive_frequency_matches_source() {
        let saw = SawSource::linbo3_3ghz();
        let freq = saw.frequency_hz;
        let drive = SawStrainDrive::new(saw, 0.0, 0.0);
        assert_eq!(drive.drive_frequency_hz(), freq);
    }

    #[test]
    fn test_piezo_ac_strain_drive_zero_at_t_zero() {
        let drive = PiezoAcStrainDrive::new(PiezoelectricSubstrate::pzt(), 1.0e-3, 100.0, 1.0e6)
            .expect("valid piezo AC drive");
        let (eps_xx, _, _) = drive.strain_tensor_at(0.0).diagonal();
        assert!(
            eps_xx.abs() < 1.0e-25,
            "sin(0) drive should give zero strain at t=0, got {eps_xx}"
        );
    }

    #[test]
    fn test_piezo_ac_strain_drive_matches_static_formula_at_quarter_period() {
        let substrate = PiezoelectricSubstrate::pzt();
        let thickness = 1.0e-3;
        let voltage_amplitude = 100.0;
        let freq = 1.0e6;
        let drive = PiezoAcStrainDrive::new(substrate, thickness, voltage_amplitude, freq)
            .expect("valid piezo AC drive");

        // At t = T/4, sin(2*pi*f*t) = sin(pi/2) = 1, so the strain should
        // match the static strain_from_voltage() at the peak voltage.
        let quarter_period = 1.0 / (4.0 * freq);
        let (eps_xx, _, eps_zz) = drive.strain_tensor_at(quarter_period).diagonal();
        let expected = substrate
            .strain_from_voltage(voltage_amplitude, thickness)
            .expect("valid static strain");
        assert!(
            (eps_xx - expected.components[0][0]).abs() < 1.0e-12,
            "AC piezo strain at peak should match static strain_from_voltage"
        );
        assert!((eps_zz - expected.components[2][2]).abs() < 1.0e-12);
    }

    #[test]
    fn test_piezo_ac_strain_drive_rejects_excessive_field() {
        let substrate = PiezoelectricSubstrate::pzt();
        // max_field is 2e6 V/m; use a tiny thickness to blow past it.
        let result = PiezoAcStrainDrive::new(substrate, 1.0e-9, 1.0, 1.0e6);
        assert!(result.is_err(), "excessive peak field should be rejected");
    }

    #[test]
    fn test_piezo_ac_strain_drive_rejects_non_positive_frequency() {
        let result = PiezoAcStrainDrive::new(PiezoelectricSubstrate::pzt(), 1.0e-3, 10.0, 0.0);
        assert!(result.is_err(), "zero frequency should be rejected");
    }

    // -------------------------------------------------------------------
    // StrainDrivenLlgDriver construction validation
    // -------------------------------------------------------------------

    #[test]
    fn test_driver_rejects_non_positive_alpha() {
        let saw = SawSource::linbo3_1ghz();
        let drive = SawStrainDrive::new(saw, 0.0, 0.0);
        let result = StrainDrivenLlgDriver::new(
            drive,
            MagnetoelasticMaterial::nickel(),
            GAMMA,
            0.0, // alpha
            5.0e3,
            0.0,
            20.0e-9,
            OPTIMAL_COUPLING_ANGLE_RAD,
        );
        assert!(result.is_err(), "alpha = 0 should be rejected");
    }

    #[test]
    fn test_driver_rejects_negative_bias_field() {
        let saw = SawSource::linbo3_1ghz();
        let drive = SawStrainDrive::new(saw, 0.0, 0.0);
        let result = StrainDrivenLlgDriver::new(
            drive,
            MagnetoelasticMaterial::nickel(),
            GAMMA,
            0.02,
            5.0e3,
            -1.0,
            20.0e-9,
            OPTIMAL_COUPLING_ANGLE_RAD,
        );
        assert!(result.is_err(), "negative bias field should be rejected");
    }

    #[test]
    fn test_driver_rejects_non_positive_film_thickness() {
        let saw = SawSource::linbo3_1ghz();
        let drive = SawStrainDrive::new(saw, 0.0, 0.0);
        let result = StrainDrivenLlgDriver::new(
            drive,
            MagnetoelasticMaterial::nickel(),
            GAMMA,
            0.02,
            5.0e3,
            0.0,
            0.0,
            OPTIMAL_COUPLING_ANGLE_RAD,
        );
        assert!(result.is_err(), "zero film thickness should be rejected");
    }

    #[test]
    fn test_evolve_rejects_non_positive_dt() {
        let saw = SawSource::linbo3_1ghz();
        let drive = SawStrainDrive::new(saw, 0.0, 0.0);
        let driver = StrainDrivenLlgDriver::new(
            drive,
            MagnetoelasticMaterial::nickel(),
            GAMMA,
            0.045,
            5.0e3,
            0.0,
            20.0e-9,
            OPTIMAL_COUPLING_ANGLE_RAD,
        )
        .expect("valid driver");
        let result = driver.evolve(Vector3::unit_z(), 0.0, 10);
        assert!(result.is_err(), "dt=0 should be rejected");
    }

    #[test]
    fn test_evolve_rejects_zero_periods() {
        let saw = SawSource::linbo3_1ghz();
        let drive = SawStrainDrive::new(saw, 0.0, 0.0);
        let driver = StrainDrivenLlgDriver::new(
            drive,
            MagnetoelasticMaterial::nickel(),
            GAMMA,
            0.045,
            5.0e3,
            0.0,
            20.0e-9,
            OPTIMAL_COUPLING_ANGLE_RAD,
        )
        .expect("valid driver");
        let result = driver.evolve(Vector3::unit_z(), 1.0e-12, 0);
        assert!(result.is_err(), "n_periods=0 should be rejected");
    }

    // -------------------------------------------------------------------
    // Norm conservation
    // -------------------------------------------------------------------

    #[test]
    fn test_norm_conservation_of_raw_integration() {
        let saw = SawSource::new(PiezoSubstrate::linbo3(), 1.0e9, 1.0e-6).expect("valid SAW");
        let drive = SawStrainDrive::new(saw, 0.0, 0.0);
        let driver = StrainDrivenLlgDriver::new(
            drive,
            MagnetoelasticMaterial::nickel(),
            GAMMA,
            0.045,
            5.0e3,
            1.0e4,
            20.0e-9,
            OPTIMAL_COUPLING_ANGLE_RAD,
        )
        .expect("valid driver");

        let trajectory = driver
            .evolve(Vector3::unit_z(), 2.0e-12, 10)
            .expect("evolve should succeed");

        assert!(
            trajectory.max_raw_norm_deviation < 1.0e-3,
            "raw (pre-renormalization) integrator output should stay close to unit norm, got deviation {:.3e}",
            trajectory.max_raw_norm_deviation
        );
        // All renormalized samples must, by construction, be exactly unit
        // norm (sanity check on the renormalization step itself).
        for sample in &trajectory.samples {
            assert!(
                (sample.magnetization.magnitude() - 1.0).abs() < 1.0e-10,
                "renormalized sample should have unit magnitude"
            );
        }
    }

    // -------------------------------------------------------------------
    // Resonance behaviour
    // -------------------------------------------------------------------

    /// Shared setup for the resonance-comparison tests: an on-resonance and
    /// an off-resonance driver built from the same SAW source and Nickel
    /// film, differing only in the static bias field `h_ext`.
    fn build_resonance_pair() -> (
        StrainDrivenLlgDriver<SawStrainDrive>,
        StrainDrivenLlgDriver<SawStrainDrive>,
        SawMagnetoacoustics,
        SawMagnetoacoustics,
    ) {
        // Small strain amplitude to stay within the small-angle regime that
        // the analytical Lorentzian model (SawMagnetoacoustics) assumes.
        let saw = SawSource::new(PiezoSubstrate::linbo3(), 1.0e9, 1.0e-6).expect("valid SAW");

        let helper = nickel_saw_device(saw.clone(), 0.0);
        let h_res = helper.acoustic_fmr_condition_h_ext();
        let mat = SawMaterial::nickel();
        // Moderate detuning: strongly suppresses the Lorentzian response
        // while keeping the off-resonance natural frequency within a
        // manageable factor of the on-resonance one (both resolvable with a
        // shared, modest time step).
        let h_off = h_res + 5.0 * mat.anisotropy_field() + 2.0e4;

        let device_on = nickel_saw_device(saw.clone(), h_res);
        let device_off = nickel_saw_device(saw, h_off);

        let driver_on = StrainDrivenLlgDriver::from_saw_magnetoacoustics(
            &device_on,
            MagnetoelasticMaterial::nickel(),
            0.0,
            0.0,
            OPTIMAL_COUPLING_ANGLE_RAD,
        )
        .expect("valid on-resonance driver");
        let driver_off = StrainDrivenLlgDriver::from_saw_magnetoacoustics(
            &device_off,
            MagnetoelasticMaterial::nickel(),
            0.0,
            0.0,
            OPTIMAL_COUPLING_ANGLE_RAD,
        )
        .expect("valid off-resonance driver");

        (driver_on, driver_off, device_on, device_off)
    }

    #[test]
    fn test_resonant_response_exceeds_off_resonant() {
        let (driver_on, driver_off, _device_on, _device_off) = build_resonance_pair();

        let dt_s = 2.0e-12;
        let n_periods = 30;
        let m0 = Vector3::unit_z();

        let trajectory_on = driver_on
            .evolve(m0, dt_s, n_periods)
            .expect("on-resonance evolution should succeed");
        let trajectory_off = driver_off
            .evolve(m0, dt_s, n_periods)
            .expect("off-resonance evolution should succeed");

        assert!(
            trajectory_on.max_cone_angle_rad > 2.0 * trajectory_off.max_cone_angle_rad,
            "on-resonance cone angle ({:.3e} rad) should substantially exceed off-resonance ({:.3e} rad)",
            trajectory_on.max_cone_angle_rad,
            trajectory_off.max_cone_angle_rad
        );
        assert!(
            trajectory_on.max_cone_angle_rad > 0.0,
            "on-resonance drive should produce nonzero precession"
        );
    }

    #[test]
    fn test_collinear_coupling_angle_gives_negligible_response() {
        // At coupling_angle_rad = 0 (strain principal axis collinear with
        // the bias/anisotropy axis), the diagonal magnetoelastic coupling
        // has no linear-order transverse torque (see module docs). The
        // driven response should collapse relative to the optimal (45 deg)
        // geometry, confirming the tilt is doing real physical work.
        let saw = SawSource::new(PiezoSubstrate::linbo3(), 1.0e9, 1.0e-6).expect("valid SAW");
        let helper = nickel_saw_device(saw.clone(), 0.0);
        let h_res = helper.acoustic_fmr_condition_h_ext();
        let device_on = nickel_saw_device(saw, h_res);

        let driver_tilted = StrainDrivenLlgDriver::from_saw_magnetoacoustics(
            &device_on,
            MagnetoelasticMaterial::nickel(),
            0.0,
            0.0,
            OPTIMAL_COUPLING_ANGLE_RAD,
        )
        .expect("valid tilted driver");
        let driver_collinear = StrainDrivenLlgDriver::from_saw_magnetoacoustics(
            &device_on,
            MagnetoelasticMaterial::nickel(),
            0.0,
            0.0,
            0.0, // collinear: no linear-order coupling expected
        )
        .expect("valid collinear driver");

        let dt_s = 2.0e-12;
        let n_periods = 30;
        let m0 = Vector3::unit_z();

        let trajectory_tilted = driver_tilted
            .evolve(m0, dt_s, n_periods)
            .expect("tilted evolution should succeed");
        let trajectory_collinear = driver_collinear
            .evolve(m0, dt_s, n_periods)
            .expect("collinear evolution should succeed");

        assert!(
            trajectory_tilted.max_cone_angle_rad > 10.0 * trajectory_collinear.max_cone_angle_rad,
            "45-degree coupling ({:.3e} rad) should give a far larger response than collinear \
             coupling ({:.3e} rad)",
            trajectory_tilted.max_cone_angle_rad,
            trajectory_collinear.max_cone_angle_rad
        );
    }

    // -------------------------------------------------------------------
    // Energy pumped into the magnetic subsystem
    // -------------------------------------------------------------------

    #[test]
    fn test_energy_pumped_in_positive_and_trend_consistent_with_resonant_absorption() {
        let (driver_on, driver_off, device_on, device_off) = build_resonance_pair();

        let dt_s = 2.0e-12;
        let n_periods = 30;
        let m0 = Vector3::unit_z();

        let trajectory_on = driver_on
            .evolve(m0, dt_s, n_periods)
            .expect("on-resonance evolution should succeed");
        let trajectory_off = driver_off
            .evolve(m0, dt_s, n_periods)
            .expect("off-resonance evolution should succeed");

        // Positivity: damping always dissipates energy, so a nonzero-
        // amplitude precession must show positive mean absorbed power.
        assert!(
            trajectory_on.mean_power_areal_density_w_per_m2 > 0.0,
            "on-resonance mean absorbed power should be positive, got {:.3e} W/m^2",
            trajectory_on.mean_power_areal_density_w_per_m2
        );
        assert!(
            trajectory_off.mean_power_areal_density_w_per_m2 >= 0.0,
            "off-resonance mean absorbed power should be non-negative, got {:.3e} W/m^2",
            trajectory_off.mean_power_areal_density_w_per_m2
        );

        // Trend: on-resonance must pump in (and dissipate) substantially
        // more power than off-resonance, in both the time-stepped model...
        assert!(
            trajectory_on.mean_power_areal_density_w_per_m2
                > 2.0 * trajectory_off.mean_power_areal_density_w_per_m2,
            "on-resonance mean power ({:.3e} W/m^2) should substantially exceed off-resonance ({:.3e} W/m^2)",
            trajectory_on.mean_power_areal_density_w_per_m2,
            trajectory_off.mean_power_areal_density_w_per_m2
        );
        // ...and in the pre-existing closed-form Lorentzian estimate.
        //
        // Note on what "consistent" means here: `resonant_absorption()` is
        // built on `precession_cone_angle_rad()`, whose own doc comment
        // acknowledges its `theta_raw = H_me * lorentzian_response()`
        // product only "simplifies to a dimensionless angle" heuristically
        // (H_me [A/m] * L(omega) [s/rad] does not actually cancel to a bare
        // radian). Cross-checking by hand confirms this: at these
        // parameters `precession_cone_angle_rad()` evaluates to O(1e-7) rad,
        // roughly six orders of magnitude below the cone angle this
        // time-stepped simulation actually reaches (consistent with
        // `resonant_absorption`'s theta^2 scaling). None of saw.rs's own
        // tests (`test_resonant_absorption_positive`,
        // `test_precession_angle_grows_with_strain`, ...) check the
        // *absolute* magnitude of these quantities either -- only
        // positivity and relative/qualitative trends. So the meaningful,
        // well-founded comparison across these two independent models
        // (small-angle closed-form Lorentzian estimate vs. a real
        // time-stepped trajectory) is the qualitative "trend" alternative
        // explicitly allowed for this test: both predict a positive,
        // resonance-enhanced absorption, asserted above. An absolute
        // order-of-magnitude bound is deliberately not asserted, since the
        // pre-existing analytical model was never calibrated (or tested)
        // for absolute agreement with a physically normalized power
        // density.
        let p_abs_on = device_on.resonant_absorption();
        let p_abs_off = device_off.resonant_absorption();
        assert!(
            p_abs_on > 0.0,
            "analytical resonant absorption should be positive"
        );
        assert!(
            p_abs_on > p_abs_off,
            "analytical model should also predict on-resonance absorption exceeding off-resonance"
        );
    }

    // -------------------------------------------------------------------
    // Diagnostics
    // -------------------------------------------------------------------

    #[test]
    fn test_induced_anisotropy_constant_is_finite_and_signed_consistently() {
        let saw = SawSource::linbo3_1ghz();
        let drive = SawStrainDrive::new(saw, 0.0, 0.0);
        let driver = StrainDrivenLlgDriver::new(
            drive,
            MagnetoelasticMaterial::cofeb(),
            GAMMA,
            0.01,
            2.0e2,
            0.0,
            20.0e-9,
            OPTIMAL_COUPLING_ANGLE_RAD,
        )
        .expect("valid driver");

        let k_sigma = driver.induced_anisotropy_constant_j_per_m3(1.0e-10);
        assert!(
            k_sigma.is_finite(),
            "induced anisotropy constant must be finite"
        );

        // Cross-check against the raw stress_induced_anisotropy formula
        // directly, confirming no double conversion/sign error crept in.
        let strain = driver.strain_drive.strain_tensor_at(1.0e-10);
        let (eps_xx, _, _) = strain.diagonal();
        let sigma_xx = MagnetoelasticMaterial::cofeb().youngs_modulus * eps_xx;
        let expected =
            stress_induced_anisotropy(MagnetoelasticMaterial::cofeb().lambda_s, sigma_xx);
        assert!((k_sigma - expected).abs() < 1.0e-6 * expected.abs().max(1.0));
    }

    // -------------------------------------------------------------------
    // Piezoelectric-driven dynamics ("piezoelectric control of magnetism")
    // -------------------------------------------------------------------

    #[test]
    fn test_piezo_driven_dynamics_evolves_and_conserves_norm() {
        let drive = PiezoAcStrainDrive::new(PiezoelectricSubstrate::pmn_pt(), 0.5e-3, 5.0, 5.0e8)
            .expect("valid piezo AC drive");
        let driver = StrainDrivenLlgDriver::new(
            drive,
            MagnetoelasticMaterial::galfenol(),
            GAMMA,
            0.02,
            1.0e3,
            5.0e3,
            30.0e-9,
            OPTIMAL_COUPLING_ANGLE_RAD,
        )
        .expect("valid piezo-driven driver");

        let m0 = Vector3::unit_z();
        let trajectory = driver
            .evolve(m0, 5.0e-13, 20)
            .expect("piezo-driven evolution should succeed");

        assert!(
            trajectory.max_raw_norm_deviation < 1.0e-3,
            "piezo-driven raw integration should conserve norm to tight tolerance, got {:.3e}",
            trajectory.max_raw_norm_deviation
        );
        assert!(
            trajectory.max_cone_angle_rad > 0.0,
            "AC piezo strain should produce genuine (nonzero) magnetization dynamics"
        );
    }
}
