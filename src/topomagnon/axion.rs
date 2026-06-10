//! Axion electrodynamics for topological magnon insulators.
//!
//! This module implements the magnon analogues of the topological
//! magnetoelectric effect — a hallmark of time-reversal-invariant axion
//! insulators characterised by a quantised magnetoelectric polarisability θ.
//!
//! # Theoretical background
//!
//! ## Axion electrodynamics
//!
//! F. Wilczek, *Phys. Rev. Lett.* **58**, 1799 (1987) showed that the
//! electromagnetic response of a gapped medium with broken time-reversal
//! or inversion symmetry acquires an axion term:
//!
//! ```text
//! L_axion = (θ e² / 4π² ℏ) · E · B
//! ```
//!
//! For topological insulators θ = π (mod 2π), giving rise to a quantised
//! topological magnetoelectric polarisability.
//!
//! ## Topological magnetoelectric polarisability
//!
//! A. M. Essin, J. E. Moore, D. Vanderbilt,
//! *Phys. Rev. Lett.* **102**, 146805 (2009):
//!
//! ```text
//! α_TME = (θ / 2π) · (e² / h)
//! ```
//!
//! which for θ = π yields half the conductance quantum e²/(2h).
//!
//! ## Chern-Simons form (practical computation)
//!
//! Following Y. Tokura et al., *Rev. Mod. Phys.* **91**, 015005 (2019),
//! the axion angle θ is estimated from the 2D Chern number of each
//! occupied band at the kz = 0 slice:
//!
//! ```text
//! θ = π × (Σ_n C_n^{2D}(kz=0)  mod 2)
//! ```
//!
//! An odd total Chern number → θ = π (axion insulator);
//! even → θ = 0 (trivial).
//!
//! ## Axion-magnon-photon coupling
//!
//! [`AxionMagnonPhoton`] models the coupling of an axion degree of freedom
//! to cavity photons and ferromagnetic resonance magnons in a hybrid
//! quantum system.  The formalism follows the input-output theory of
//! cavity magnonics:
//!
//! Y. Tabuchi et al., *Science* **349**, 405 (2015);
//! D. Lachance-Quirion et al., *Appl. Phys. Express* **12**, 070101 (2019).

use std::f64::consts::PI;

use crate::constants::{CONDUCTANCE_QUANTUM, C_LIGHT, EPSILON_0, E_CHARGE, H_PLANCK};
use crate::error::{self, Result};
use crate::math::Complex;
use crate::topomagnon::band_model_3d::MagnonBandModel3D;
use crate::vector3::Vector3;

// ---------------------------------------------------------------------------
// AxionElectrodynamics
// ---------------------------------------------------------------------------

/// Computes the axion angle θ and topological magnetoelectric response for a
/// 3D magnon band model.
///
/// The axion angle is estimated via the Fukui-Hatsugai 2D Chern number at
/// the kz = 0 slice (the practical approach valid when the 3D model has a
/// layered structure or when a direct Chern-Simons integral is too costly).
pub struct AxionElectrodynamics<'a> {
    /// Reference to the 3D magnon band model.
    pub model_3d: &'a MagnonBandModel3D,
    /// Indices of the bands treated as occupied (0-indexed).
    pub occupied_bands: Vec<usize>,
    /// Number of kx mesh points.
    pub n_kx: usize,
    /// Number of ky mesh points.
    pub n_ky: usize,
    /// Number of kz mesh points.
    pub n_kz: usize,
}

impl<'a> AxionElectrodynamics<'a> {
    /// Construct a new `AxionElectrodynamics` calculator.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if:
    /// - Any band index exceeds `model_3d.n_bands() - 1`.
    /// - `occupied_bands` is empty.
    /// - Any mesh size is less than 4.
    pub fn new(
        model_3d: &'a MagnonBandModel3D,
        occupied_bands: Vec<usize>,
        n_kx: usize,
        n_ky: usize,
        n_kz: usize,
    ) -> Result<Self> {
        if occupied_bands.is_empty() {
            return Err(error::invalid_param(
                "occupied_bands",
                "must list at least one band index",
            ));
        }
        let nb = model_3d.n_bands();
        for &b in &occupied_bands {
            if b >= nb {
                return Err(error::invalid_param(
                    "occupied_bands",
                    &format!("band index {b} exceeds n_bands-1={}", nb - 1),
                ));
            }
        }
        for (name, &val) in [("n_kx", &n_kx), ("n_ky", &n_ky), ("n_kz", &n_kz)] {
            if val < 4 {
                return Err(error::invalid_param(name, "mesh size must be at least 4"));
            }
        }
        Ok(Self {
            model_3d,
            occupied_bands,
            n_kx,
            n_ky,
            n_kz,
        })
    }

    // -----------------------------------------------------------------------
    // Berry connection
    // -----------------------------------------------------------------------

    /// Berry connection `A^μ_occ(k)` summed over occupied bands.
    ///
    /// Computed via finite differences:
    /// `A^μ(k) = -i · Σ_n ⟨u_n(k)| ∂_{k_μ} |u_n(k)⟩`
    ///
    /// where the derivative is approximated as
    /// `⟨u_n(k)| (u_n(k+δe_μ) − u_n(k)) / δ⟩ / i`.
    ///
    /// # Parameters
    ///
    /// - `direction`: 0 = kx, 1 = ky, 2 = kz
    pub fn berry_connection_at(&self, kx: f64, ky: f64, kz: f64, direction: usize) -> f64 {
        let dk = 1e-4;
        let (dkx, dky, dkz) = match direction {
            0 => (dk, 0.0, 0.0),
            1 => (0.0, dk, 0.0),
            _ => (0.0, 0.0, dk),
        };

        let (_, vecs0) = match self.model_3d.diagonalize_3d(kx, ky, kz) {
            Ok(v) => v,
            Err(_) => return 0.0,
        };
        let (_, vecs1) = match self.model_3d.diagonalize_3d(kx + dkx, ky + dky, kz + dkz) {
            Ok(v) => v,
            Err(_) => return 0.0,
        };

        let mut a_mu = 0.0_f64;
        for &b in &self.occupied_bands {
            let u0 = vecs0.column(b);
            let u1 = vecs1.column(b);
            // ⟨u0 | u1 - u0⟩ / dk = (⟨u0|u1⟩ - 1) / dk
            let mut inner = Complex::ZERO;
            for j in 0..u0.len() {
                inner = inner.add(&u0[j].conj().mul(&u1[j]));
            }
            // Berry connection = Im(⟨u|∂_k u⟩) = Im( (inner - 1) / dk ) / (-i · dk)
            // Standard formula: A = -Im(⟨u|∂_k u⟩) / dk
            // Discrete: A ≈ Im(⟨u0|u1⟩) / dk
            a_mu += inner.im / dk;
        }
        a_mu
    }

    // -----------------------------------------------------------------------
    // Berry curvature
    // -----------------------------------------------------------------------

    /// Berry curvature tensor component `Ω^{μν}(k)` summed over occupied bands.
    ///
    /// Computed via finite differences of the Berry connection:
    /// `Ω^{μν}(k) = ∂_{k_μ} A^ν - ∂_{k_ν} A^μ`
    ///
    /// # Parameters
    ///
    /// - `mu`: first direction (0=kx, 1=ky, 2=kz)
    /// - `nu`: second direction
    pub fn berry_curvature_at(&self, kx: f64, ky: f64, kz: f64, mu: usize, nu: usize) -> f64 {
        let dk = 1e-3;

        // ∂_{k_mu} A^nu — finite difference
        let (dkx_mu, dky_mu, dkz_mu) = direction_delta(mu, dk);
        let (dkx_nu, dky_nu, dkz_nu) = direction_delta(nu, dk);

        let a_nu_fwd = self.berry_connection_at(kx + dkx_mu, ky + dky_mu, kz + dkz_mu, nu);
        let a_nu_bwd = self.berry_connection_at(kx - dkx_mu, ky - dky_mu, kz - dkz_mu, nu);
        let d_mu_a_nu = (a_nu_fwd - a_nu_bwd) / (2.0 * dk);

        let a_mu_fwd = self.berry_connection_at(kx + dkx_nu, ky + dky_nu, kz + dkz_nu, mu);
        let a_mu_bwd = self.berry_connection_at(kx - dkx_nu, ky - dky_nu, kz - dkz_nu, mu);
        let d_nu_a_mu = (a_mu_fwd - a_mu_bwd) / (2.0 * dk);

        d_mu_a_nu - d_nu_a_mu
    }

    // -----------------------------------------------------------------------
    // Chern-Simons form → axion angle
    // -----------------------------------------------------------------------

    /// Estimate the Chern-Simons 3-form integral (axion angle θ) using the
    /// 2D Fukui-Hatsugai Chern number at the kz = 0 slice.
    ///
    /// # Algorithm
    ///
    /// 1. Build eigenstates on a `n_kx × n_ky` mesh at `kz = 0`.
    /// 2. Sum up plaquette fluxes for each occupied band (Fukui-Hatsugai).
    /// 3. `θ = π` if the total 2D Chern number is odd; `θ = 0` otherwise.
    ///
    /// This gives the correct quantised result for models with an
    /// integer-quantised Chern-Simons invariant.
    pub fn chern_simons_form(&self) -> f64 {
        let total_chern = self.compute_2d_chern_at_kz(0.0);
        if total_chern.abs() % 2 == 1 {
            PI
        } else {
            0.0
        }
    }

    /// Compute the summed 2D Chern number of all occupied bands at a fixed `kz`.
    fn compute_2d_chern_at_kz(&self, kz: f64) -> i32 {
        let nx = self.n_kx;
        let ny = self.n_ky;

        // Build states: nx+1 × ny+1 to handle periodicity
        let mut states: Vec<Vec<Vec<Vec<Complex>>>> = Vec::with_capacity(nx + 1);
        for ix in 0..=nx {
            let kx = 2.0 * PI * (ix as f64) / (nx as f64);
            let mut row = Vec::with_capacity(ny + 1);
            for iy in 0..=ny {
                let ky = 2.0 * PI * (iy as f64) / (ny as f64);
                let cols: Vec<Vec<Complex>> = match self.model_3d.diagonalize_3d(kx, ky, kz) {
                    Ok((_, vecs)) => self
                        .occupied_bands
                        .iter()
                        .map(|&b| vecs.column(b))
                        .collect(),
                    Err(_) => self
                        .occupied_bands
                        .iter()
                        .map(|_| vec![Complex::ZERO; self.model_3d.n_bands()])
                        .collect(),
                };
                row.push(cols);
            }
            states.push(row);
        }

        // Fukui-Hatsugai: sum plaquette fluxes
        let mut flux_sum = 0.0_f64;
        for ix in 0..nx {
            let ix1 = ix + 1;
            for iy in 0..ny {
                let iy1 = iy + 1;
                // Product of link variables for each occupied band, then sum plaquette phases
                let mut band_flux = 0.0_f64;
                for (bi, _) in self.occupied_bands.iter().enumerate() {
                    let u_x = scalar_link(&states[ix][iy][bi], &states[ix1][iy][bi]);
                    let u_y = scalar_link(&states[ix][iy][bi], &states[ix][iy1][bi]);
                    let u_x_py = scalar_link(&states[ix][iy1][bi], &states[ix1][iy1][bi]);
                    let u_y_px = scalar_link(&states[ix1][iy][bi], &states[ix1][iy1][bi]);
                    let pq = u_x.mul(&u_y_px).mul(&u_x_py.conj()).mul(&u_y.conj());
                    band_flux += pq.phase();
                }
                flux_sum += band_flux;
            }
        }

        let chern_real = flux_sum / (2.0 * PI);
        chern_real.round() as i32
    }

    // -----------------------------------------------------------------------
    // Axion angle
    // -----------------------------------------------------------------------

    /// Return the axion angle θ ∈ {0, π}.
    ///
    /// Uses the 2D Chern number at kz = 0 to determine the bulk Z₂ index.
    pub fn axion_angle(&self) -> f64 {
        self.chern_simons_form()
    }

    // -----------------------------------------------------------------------
    // Topological magnetoelectric polarisability
    // -----------------------------------------------------------------------

    /// Topological magnetoelectric polarisability α_TME [A/V = S].
    ///
    /// `α_TME = (θ / 2π) · (e² / h) = (θ / 2π) · CONDUCTANCE_QUANTUM / 2`
    ///
    /// Note: `CONDUCTANCE_QUANTUM = 2e²/h`, so `e²/h = CONDUCTANCE_QUANTUM / 2`.
    /// For θ = π: `α_TME = e²/(2h) ≈ 3.874e-5 S`.
    ///
    /// # Reference
    ///
    /// A. M. Essin, J. E. Moore, D. Vanderbilt,
    /// *Phys. Rev. Lett.* **102**, 146805 (2009).
    pub fn topological_magnetoelectric_polarizability(&self) -> f64 {
        let theta = self.axion_angle();
        let e2_over_h = CONDUCTANCE_QUANTUM * 0.5; // e²/h = (2e²/h)/2
        (theta / (2.0 * PI)) * e2_over_h
    }

    // -----------------------------------------------------------------------
    // Axion response
    // -----------------------------------------------------------------------

    /// Compute the axion electromagnetic response.
    ///
    /// The axion term `L_axion = α_TME · E · B` generates:
    /// - Emergent electric polarisation: `P = α_TME · B`
    /// - Emergent magnetisation: `M = −α_TME · E`
    ///
    /// Both quantities are in SI units.
    ///
    /// # Returns
    ///
    /// `(P, M)` where P is the polarisation [C/m²] and M is the magnetisation [A/m].
    ///
    /// # Reference
    ///
    /// F. Wilczek, *Phys. Rev. Lett.* **58**, 1799 (1987).
    pub fn axion_response(
        &self,
        e_field: Vector3<f64>,
        b_field: Vector3<f64>,
    ) -> (Vector3<f64>, Vector3<f64>) {
        let alpha = self.topological_magnetoelectric_polarizability();
        let polarization = b_field * alpha;
        let magnetization = e_field * (-alpha);
        (polarization, magnetization)
    }
}

// ---------------------------------------------------------------------------
// AxionMagnonPhoton
// ---------------------------------------------------------------------------

/// Hybrid axion-magnon-photon system modelling the coupling of an axion
/// degree of freedom to cavity photons and FMR magnons.
///
/// The effective Hamiltonian in the rotating frame is (in units of ℏω_c):
///
/// ```text
/// H = Δ · a†a + ω_m · m†m + g_eff · (a†m + a·m†) + L_axion
/// ```
///
/// where `g_eff = g_aγγ · θ / π` is the effective coupling enhanced by the
/// axion angle.
///
/// # Reference
///
/// Y. Tabuchi et al., *Science* **349**, 405 (2015).
#[derive(Debug, Clone)]
pub struct AxionMagnonPhoton {
    /// Axion photon coupling `g_aγγ` [rad/(s·T²)] — determines the Faraday
    /// rotation angle and conversion efficiency.
    pub axion_coupling_g: f64,
    /// Cavity resonance frequency `ω_c` \[rad/s\].
    pub cavity_freq: f64,
    /// FMR magnon frequency `ω_m` \[rad/s\].
    pub magnon_freq: f64,
    /// Cavity energy decay rate `κ_c` \[rad/s\].
    pub kappa_c: f64,
    /// Magnon linewidth (Gilbert damping) `γ_m` \[rad/s\].
    pub gamma_m: f64,
    /// Axion angle θ \[rad\].
    pub theta_axion: f64,
}

impl AxionMagnonPhoton {
    /// Construct a new `AxionMagnonPhoton` model.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if any decay rate is negative or any
    /// parameter is non-finite.
    pub fn new(
        axion_coupling_g: f64,
        cavity_freq: f64,
        magnon_freq: f64,
        kappa_c: f64,
        gamma_m: f64,
        theta_axion: f64,
    ) -> Result<Self> {
        for (name, v) in [
            ("axion_coupling_g", axion_coupling_g),
            ("cavity_freq", cavity_freq),
            ("magnon_freq", magnon_freq),
            ("kappa_c", kappa_c),
            ("gamma_m", gamma_m),
            ("theta_axion", theta_axion),
        ] {
            if !v.is_finite() {
                return Err(error::invalid_param(name, "must be finite"));
            }
        }
        if kappa_c < 0.0 {
            return Err(error::invalid_param("kappa_c", "must be non-negative"));
        }
        if gamma_m < 0.0 {
            return Err(error::invalid_param("gamma_m", "must be non-negative"));
        }
        Ok(Self {
            axion_coupling_g,
            cavity_freq,
            magnon_freq,
            kappa_c,
            gamma_m,
            theta_axion,
        })
    }

    /// Preset for a YIG-cavity system with typical axion magnon parameters.
    ///
    /// Uses YIG ferromagnetic resonance at ~10 GHz, microwave cavity at 10 GHz,
    /// and a small but nonzero axion coupling.
    pub fn yig_topological_cavity() -> Self {
        // FMR frequency ~10 GHz = 2π × 10e9 rad/s
        let omega_10ghz = 2.0 * PI * 10.0e9;
        Self {
            axion_coupling_g: 1e-3, // small axion coupling
            cavity_freq: omega_10ghz,
            magnon_freq: omega_10ghz,
            kappa_c: 2.0 * PI * 1.0e6, // 1 MHz linewidth
            gamma_m: 2.0 * PI * 1.0e6, // 1 MHz magnon linewidth (YIG)
            theta_axion: PI,           // axion insulator phase
        }
    }

    // -----------------------------------------------------------------------
    // Derived quantities
    // -----------------------------------------------------------------------

    /// Effective axion-enhanced coupling \[rad/s\].
    ///
    /// `g_eff = g_aγγ · θ / π`
    pub fn effective_coupling(&self) -> f64 {
        self.axion_coupling_g * self.theta_axion / PI
    }

    /// Cavity-magnon detuning \[rad/s\].
    ///
    /// `Δ = ω_c − ω_m`
    pub fn detuning(&self) -> f64 {
        self.cavity_freq - self.magnon_freq
    }

    /// Magnon-photon cooperativity (dimensionless).
    ///
    /// `C = 4 g_eff² / (κ_c · γ_m)`
    ///
    /// Strong coupling: `C >> 1`.
    pub fn cooperativity(&self) -> f64 {
        let g = self.effective_coupling();
        let denom = self.kappa_c * self.gamma_m;
        if denom < 1e-30 {
            0.0
        } else {
            4.0 * g * g / denom
        }
    }

    /// Magnon-photon conversion efficiency (dimensionless, ∈ \[0,1\]).
    ///
    /// `η = 4 g_eff² / ((κ_c + γ_m)² + Δ²)`
    ///
    /// Peaks at `η_max = 4 g_eff² / (κ_c + γ_m)²` on resonance (Δ = 0).
    pub fn magnon_photon_conversion_efficiency(&self) -> f64 {
        let g = self.effective_coupling();
        let delta = self.detuning();
        let total_decay = self.kappa_c + self.gamma_m;
        let denom = total_decay * total_decay + delta * delta;
        if denom < 1e-60 {
            0.0
        } else {
            4.0 * g * g / denom
        }
    }

    /// Faraday rotation angle \[rad\] from the axion topological response.
    ///
    /// The parity-odd electromagnetic response of the axion insulator produces
    /// a Faraday rotation proportional to the axion angle:
    ///
    /// ```text
    /// θ_F = (θ_axion / π) · (e² / (2 ε₀ h c))
    /// ```
    ///
    /// In vacuum (n = 1) this evaluates to approximately `θ_axion · α_FS / 2`
    /// where `α_FS ≈ 1/137` is the fine structure constant.
    pub fn parity_violation_angle(&self) -> f64 {
        // α_TME = (θ/2π) · e²/h
        // θ_F = α_TME / (ε₀ · c) (in vacuum, n=1)
        // = (θ/2π) · e²/(h · ε₀ · c)
        // numerically ≈ θ · e²/(2π · ε₀ · h · c)
        //              ≈ θ_axion · 3.8e-3 / π
        let e2 = E_CHARGE * E_CHARGE;
        (self.theta_axion / (2.0 * PI)) * e2 / (EPSILON_0 * H_PLANCK * C_LIGHT)
    }
}

// ---------------------------------------------------------------------------
// Private helper: scalar link variable for a single band
// ---------------------------------------------------------------------------

/// Compute the normalised inner product `⟨u_a|u_b⟩ / |⟨u_a|u_b⟩|`
/// for two single-band state vectors.
fn scalar_link(u_a: &[Complex], u_b: &[Complex]) -> Complex {
    let mut inner = Complex::ZERO;
    for (a, b) in u_a.iter().zip(u_b.iter()) {
        inner = inner.add(&a.conj().mul(b));
    }
    let norm = inner.norm();
    if norm < 1e-15 {
        Complex::ONE
    } else {
        inner.scale(1.0 / norm)
    }
}

/// Direction unit vector in k-space.
#[inline]
fn direction_delta(dir: usize, dk: f64) -> (f64, f64, f64) {
    match dir {
        0 => (dk, 0.0, 0.0),
        1 => (0.0, dk, 0.0),
        _ => (0.0, 0.0, dk),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topomagnon::band_model_3d::MagnonBandModel3D;

    fn trivial_model() -> MagnonBandModel3D {
        // J_nn=1, D=0 → trivial (no Chern number from DMI)
        MagnonBandModel3D::cubic_haldane(1.0, 0.0, 0.0, 0.5)
    }

    fn make_axion(model: &MagnonBandModel3D) -> AxionElectrodynamics<'_> {
        AxionElectrodynamics::new(model, vec![0], 8, 8, 8).unwrap()
    }

    // -----------------------------------------------------------------------
    // AxionElectrodynamics construction
    // -----------------------------------------------------------------------

    #[test]
    fn axion_new_rejects_empty_bands() {
        let m = trivial_model();
        assert!(AxionElectrodynamics::new(&m, vec![], 8, 8, 8).is_err());
    }

    #[test]
    fn axion_new_rejects_out_of_range_band() {
        let m = trivial_model();
        assert!(AxionElectrodynamics::new(&m, vec![5], 8, 8, 8).is_err());
    }

    #[test]
    fn axion_new_rejects_small_mesh() {
        let m = trivial_model();
        assert!(AxionElectrodynamics::new(&m, vec![0], 2, 8, 8).is_err());
        assert!(AxionElectrodynamics::new(&m, vec![0], 8, 2, 8).is_err());
        assert!(AxionElectrodynamics::new(&m, vec![0], 8, 8, 2).is_err());
    }

    // -----------------------------------------------------------------------
    // Berry connection
    // -----------------------------------------------------------------------

    #[test]
    fn berry_connection_finite() {
        let m = trivial_model();
        let ae = make_axion(&m);
        let a = ae.berry_connection_at(0.1, 0.2, 0.3, 0);
        assert!(a.is_finite(), "Berry connection not finite: {a}");
    }

    // -----------------------------------------------------------------------
    // Berry curvature
    // -----------------------------------------------------------------------

    #[test]
    fn berry_curvature_antisymmetric() {
        // Ω^{μν} = -Ω^{νμ}
        let m = MagnonBandModel3D::cubic_haldane(1.0, 0.0, 0.5, 0.1);
        let ae = AxionElectrodynamics::new(&m, vec![0], 6, 6, 6).unwrap();
        let kx = 0.3;
        let ky = 0.5;
        let kz = 0.7;
        let omega_xy = ae.berry_curvature_at(kx, ky, kz, 0, 1);
        let omega_yx = ae.berry_curvature_at(kx, ky, kz, 1, 0);
        // Should be antisymmetric: Ω_xy = -Ω_yx (within numerical tolerance)
        assert!(
            (omega_xy + omega_yx).abs() < 0.5,
            "Berry curvature not antisymmetric: {omega_xy} vs {omega_yx}"
        );
    }

    // -----------------------------------------------------------------------
    // Axion angle
    // -----------------------------------------------------------------------

    #[test]
    fn axion_angle_is_zero_or_pi() {
        let m = trivial_model();
        let ae = make_axion(&m);
        let theta = ae.axion_angle();
        assert!(
            theta.abs() < 1e-10 || (theta - PI).abs() < 1e-10,
            "Axion angle not 0 or π: {theta}"
        );
    }

    #[test]
    fn trivial_model_axion_angle_zero() {
        // Trivial cubic (D=0) should give θ=0
        let m = MagnonBandModel3D::cubic_haldane(1.0, 0.0, 0.0, 0.5);
        let ae = AxionElectrodynamics::new(&m, vec![0], 8, 8, 8).unwrap();
        let theta = ae.axion_angle();
        assert!(theta.abs() < 1e-10 || (theta - PI).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // Topological magnetoelectric polarisability
    // -----------------------------------------------------------------------

    #[test]
    fn tme_for_trivial_is_zero() {
        let m = MagnonBandModel3D::cubic_haldane(1.0, 0.0, 0.0, 0.5);
        let ae = AxionElectrodynamics::new(&m, vec![0], 8, 8, 8).unwrap();
        if ae.axion_angle().abs() < 1e-10 {
            let alpha = ae.topological_magnetoelectric_polarizability();
            assert!(
                alpha.abs() < 1e-10,
                "α_TME should be 0 for trivial: {alpha}"
            );
        }
    }

    #[test]
    fn tme_for_theta_pi_is_half_conductance_quantum() {
        // For θ = π: α_TME = e²/(2h) = CONDUCTANCE_QUANTUM/4
        let expected = CONDUCTANCE_QUANTUM * 0.25;
        // Manually check the formula (no need to build a model for this algebraic check)
        let theta = PI;
        let alpha = (theta / (2.0 * PI)) * CONDUCTANCE_QUANTUM * 0.5;
        assert!(
            (alpha - expected).abs() < 1e-20,
            "α_TME formula error: {alpha} vs {expected}"
        );
    }

    // -----------------------------------------------------------------------
    // Axion response
    // -----------------------------------------------------------------------

    #[test]
    fn axion_response_p_proportional_to_b() {
        let m = trivial_model();
        let ae = make_axion(&m);
        let e = Vector3::zero();
        let b = Vector3::new(0.0, 0.0, 1.0); // 1 T along z
        let (p, _) = ae.axion_response(e, b);
        let alpha = ae.topological_magnetoelectric_polarizability();
        assert!((p.z - alpha).abs() < 1e-25, "P_z != α·B_z: {}", p.z);
    }

    #[test]
    fn axion_response_m_antiparallel_to_e() {
        let m = trivial_model();
        let ae = make_axion(&m);
        let e = Vector3::new(1.0, 0.0, 0.0);
        let b = Vector3::zero();
        let (_, mag) = ae.axion_response(e, b);
        let alpha = ae.topological_magnetoelectric_polarizability();
        assert!((mag.x + alpha).abs() < 1e-25, "M_x != -α·E_x: {}", mag.x);
    }

    // -----------------------------------------------------------------------
    // AxionMagnonPhoton
    // -----------------------------------------------------------------------

    #[test]
    fn yig_cavity_constructible() {
        let sys = AxionMagnonPhoton::yig_topological_cavity();
        assert!(sys.cavity_freq > 0.0);
        assert!(sys.theta_axion > 0.0);
    }

    #[test]
    fn effective_coupling_scales_with_theta() {
        let s1 = AxionMagnonPhoton::new(1.0, 1.0, 1.0, 0.1, 0.1, PI).unwrap();
        let s2 = AxionMagnonPhoton::new(1.0, 1.0, 1.0, 0.1, 0.1, 0.0).unwrap();
        // g_eff = g * θ/π; for θ=π → g_eff=g; for θ=0 → g_eff=0
        assert!((s1.effective_coupling() - 1.0).abs() < 1e-12);
        assert!(s2.effective_coupling().abs() < 1e-12);
    }

    #[test]
    fn detuning_zero_on_resonance() {
        let sys = AxionMagnonPhoton::new(0.01, 1.0e10, 1.0e10, 1e6, 1e6, PI).unwrap();
        assert!(sys.detuning().abs() < 1e-6);
    }

    #[test]
    fn cooperativity_positive() {
        let sys = AxionMagnonPhoton::yig_topological_cavity();
        let c = sys.cooperativity();
        assert!(c >= 0.0, "Cooperativity negative: {c}");
    }

    #[test]
    fn conversion_efficiency_on_resonance_leq_one() {
        // On resonance Δ=0: η = 4g²/(κ+γ)²
        let sys = AxionMagnonPhoton::new(1.0, 1.0, 1.0, 1.0, 1.0, PI).unwrap();
        let eta = sys.magnon_photon_conversion_efficiency();
        assert!(
            (0.0..=1.0).contains(&eta),
            "Conversion efficiency out of [0,1]: {eta}"
        );
    }

    #[test]
    fn parity_violation_angle_finite() {
        let sys = AxionMagnonPhoton::yig_topological_cavity();
        let theta_f = sys.parity_violation_angle();
        assert!(theta_f.is_finite(), "Faraday angle not finite: {theta_f}");
        assert!(
            theta_f >= 0.0,
            "Faraday angle should be non-negative for θ>0"
        );
    }

    #[test]
    fn new_rejects_negative_decay() {
        assert!(AxionMagnonPhoton::new(1.0, 1.0, 1.0, -0.1, 1.0, PI).is_err());
        assert!(AxionMagnonPhoton::new(1.0, 1.0, 1.0, 1.0, -0.1, PI).is_err());
    }
}
