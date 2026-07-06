//! Momentum-resolved electronic band model for collinear altermagnets.
//!
//! This module provides [`AltermagnetBandModel`], a continuum k·p tight-binding
//! model that resolves the non-relativistic spin splitting of an altermagnet
//! directly in momentum space. Unlike the scalar [`Altermagnet`] material
//! description, the band model exposes per-spin band energies, closed-form Berry
//! curvature, and zero-temperature Fermi-sea Hall/spin-Hall integrals.
//!
//! # Physical model
//!
//! The Néel order is collinear along `ẑ` and spin-orbit coupling (SOC) is off by
//! default, so spin is a good quantum number and the 4×4 Bloch Hamiltonian in the
//! basis `(A↑, B↑, A↓, B↓)` is block diagonal:
//!
//! ```text
//! H(k) = diag( H_up(k), H_down(k) )
//! H_sigma(k) = eps0(k)*tau0 + d_sigma(k) . tau
//! d_sigma(k) = ( Re(gamma), -Im(gamma), eps_a(k) - sigma*M )
//! ```
//!
//! with sublattice-space Pauli matrices `tau` and `sigma = +1` for spin up,
//! `sigma = -1` for spin down. The altermagnetic form factor
//! `eps_a(k) = t_am*(a|k|)^2 * cos(n*(phi - phi_neel))` carries the crystal
//! symmetry via the harmonic order `n` (2/4/6 for d/g/i-wave). Because the
//! spin-splitting term enters `d_z` with a sign set by `sigma`, spin-up and
//! spin-down bands split at generic `k` while the net magnetization stays zero.
//!
//! # References
//!
//! - L. Šmejkal, J. Sinova, T. Jungwirth, "Emerging Research Landscape of
//!   Altermagnetism", Phys. Rev. X 12, 040501 (2022).
//! - L. Šmejkal et al., "Crystal Hall effect in collinear antiferromagnets",
//!   Sci. Adv. 6, eaaz8809 (2020).

use crate::altermagnet::materials::{Altermagnet, AltermagneticSymmetry};
use crate::error::{Error, Result};
use crate::math::{CMatrix, Complex};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Electron spin projection along the collinear Néel axis (`ẑ`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Spin {
    /// Spin up, `sigma = +1`.
    Up,
    /// Spin down, `sigma = -1`.
    Down,
}

impl Spin {
    /// Signed value `sigma` used throughout the Hamiltonian (`+1` up, `-1` down).
    #[inline]
    pub fn sigma(self) -> f64 {
        match self {
            Spin::Up => 1.0,
            Spin::Down => -1.0,
        }
    }
}

/// Selects one of the two eigen-bands of a per-spin 2×2 sector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Band {
    /// Lower band `E_minus = eps0(k) - |d_sigma(k)|`.
    Lower,
    /// Upper band `E_plus = eps0(k) + |d_sigma(k)|`.
    Upper,
}

/// The two eigen-energies of a single spin sector at a given `k` point \[eV\].
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpinBands {
    /// Lower band energy `eps0(k) - |d_sigma(k)|` \[eV\].
    pub lower: f64,
    /// Upper band energy `eps0(k) + |d_sigma(k)|` \[eV\].
    pub upper: f64,
}

/// Momentum-resolved two-sublattice electronic band model of a collinear
/// altermagnet.
///
/// See the [module documentation](self) for the full Hamiltonian. Construct one
/// from a scalar [`Altermagnet`] preset via [`AltermagnetBandModel::from_altermagnet`]
/// or use the named constructors [`AltermagnetBandModel::mnte`],
/// [`AltermagnetBandModel::crsb`], and [`AltermagnetBandModel::ruo2`].
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AltermagnetBandModel {
    /// Angular symmetry of the altermagnetic spin splitting.
    pub symmetry: AltermagneticSymmetry,
    /// Kinetic hopping scale entering `eps0(k) = t_kin*(a|k|)^2` \[eV\].
    pub t_kin: f64,
    /// Sublattice hybridization amplitude (`gamma`, real by default) \[eV\].
    pub delta_hyb: f64,
    /// Néel exchange splitting `M` \[eV\].
    pub exchange_m: f64,
    /// Altermagnetic form-factor amplitude `t_am` \[eV\].
    pub t_am: f64,
    /// Optional spin-orbit coupling strength (`0` disables SOC) \[eV\].
    pub lambda_soc: f64,
    /// Fermi energy used by the Fermi-sea transport integrals \[eV\].
    pub fermi_energy: f64,
    /// Néel-vector in-plane angle `phi_neel` \[rad\].
    pub neel_angle: f64,
    /// Lattice constant `a` \[m\].
    pub a_lattice: f64,
}

impl AltermagnetBandModel {
    /// Construct a band model from explicit parameters, validating physical ranges.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] if `t_kin <= 0`, `delta_hyb < 0`,
    /// `exchange_m < 0`, `t_am < 0`, `a_lattice <= 0`, or `fermi_energy` is not finite.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        symmetry: AltermagneticSymmetry,
        t_kin: f64,
        delta_hyb: f64,
        exchange_m: f64,
        t_am: f64,
        lambda_soc: f64,
        fermi_energy: f64,
        neel_angle: f64,
        a_lattice: f64,
    ) -> Result<Self> {
        if t_kin <= 0.0 || t_kin.is_nan() {
            return Err(Error::InvalidParameter {
                param: "t_kin".to_string(),
                reason: "kinetic hopping must be strictly positive".to_string(),
            });
        }
        if delta_hyb < 0.0 {
            return Err(Error::InvalidParameter {
                param: "delta_hyb".to_string(),
                reason: "hybridization amplitude must be non-negative".to_string(),
            });
        }
        if exchange_m < 0.0 {
            return Err(Error::InvalidParameter {
                param: "exchange_m".to_string(),
                reason: "exchange splitting must be non-negative".to_string(),
            });
        }
        if t_am < 0.0 {
            return Err(Error::InvalidParameter {
                param: "t_am".to_string(),
                reason: "altermagnetic amplitude must be non-negative".to_string(),
            });
        }
        if a_lattice <= 0.0 || a_lattice.is_nan() {
            return Err(Error::InvalidParameter {
                param: "a_lattice".to_string(),
                reason: "lattice constant must be strictly positive".to_string(),
            });
        }
        if !fermi_energy.is_finite() {
            return Err(Error::InvalidParameter {
                param: "fermi_energy".to_string(),
                reason: "Fermi energy must be finite".to_string(),
            });
        }
        Ok(Self {
            symmetry,
            t_kin,
            delta_hyb,
            exchange_m,
            t_am,
            lambda_soc,
            fermi_energy,
            neel_angle,
            a_lattice,
        })
    }

    /// Build a band model from a scalar [`Altermagnet`] material preset.
    ///
    /// Uses `t_kin = 1.0`, `delta_hyb = 0.3`, `exchange_m = t_am = spin_splitting`,
    /// `lambda_soc = 0`, `fermi_energy = 0.5`, and `neel_angle = 0`.
    pub fn from_altermagnet(m: &Altermagnet) -> Self {
        // The preset parameters always satisfy `new`'s validation, so the model is
        // constructed directly (spin_splitting >= 0 and lattice_constant > 0).
        Self {
            symmetry: m.symmetry,
            t_kin: 1.0,
            delta_hyb: 0.3,
            exchange_m: m.spin_splitting,
            t_am: m.spin_splitting,
            lambda_soc: 0.0,
            fermi_energy: 0.5,
            neel_angle: 0.0,
            a_lattice: m.lattice_constant,
        }
    }

    /// MnTe (manganese telluride), a d-wave altermagnet band model.
    pub fn mnte() -> Self {
        Self::from_altermagnet(&Altermagnet::mnte())
    }

    /// CrSb (chromium antimonide), a g-wave altermagnet band model.
    pub fn crsb() -> Self {
        Self::from_altermagnet(&Altermagnet::crsb())
    }

    /// RuO2 (ruthenium dioxide), a d-wave altermagnet band model.
    pub fn ruo2() -> Self {
        Self::from_altermagnet(&Altermagnet::ruo2())
    }

    /// Spin-independent kinetic term `eps0(k) = t_kin*(a|k|)^2` \[eV\].
    #[inline]
    fn eps0(&self, kx: f64, ky: f64) -> f64 {
        let k_mag = (kx * kx + ky * ky).sqrt();
        let ak = self.a_lattice * k_mag;
        self.t_kin * ak * ak
    }

    /// Altermagnetic form factor
    /// `eps_a(k) = t_am*(a|k|)^2*cos(n*(phi - phi_neel))` \[eV\].
    #[inline]
    fn eps_a(&self, kx: f64, ky: f64) -> f64 {
        let k_mag = (kx * kx + ky * ky).sqrt();
        let ak = self.a_lattice * k_mag;
        let phi = ky.atan2(kx);
        let n = self.symmetry.harmonic_order() as f64;
        self.t_am * ak * ak * (n * (phi - self.neel_angle)).cos()
    }

    /// Sublattice hybridization `gamma(k)` (a real constant by default).
    #[inline]
    fn gamma(&self, _kx: f64, _ky: f64) -> Complex {
        Complex::from_real(self.delta_hyb)
    }

    /// The three-component sublattice-space vector `d_sigma(k)` for the given spin.
    ///
    /// With SOC off (`lambda_soc = 0`) this is exactly the collinear-Néel form
    /// `d = (Re(gamma), -Im(gamma), eps_a(k) - sigma*M)`, so `dy = 0` identically.
    ///
    /// When SOC is on, each `tau` channel picks up a spin-locked structure factor
    /// `dx += sigma*lambda_soc*b_x(k)`, `dy += sigma*lambda_soc*b_y(k)`,
    /// `dz += sigma*lambda_soc*b_z(k)`, with
    /// `b_x(k) = sin(a*kx)` (odd in k), `b_y(k) = cos(a*kx) + cos(a*ky)` (EVEN in
    /// k), and `b_z(k) = sin(a*kx) + sin(a*ky)` (odd in k). This mixed parity —
    /// `x`/`z` odd, `y` even — is required by the time-reversal-to-Néel-reversal
    /// operator identity `Theta H(k; M) Theta^-1 = H(-k; -M)` with
    /// `Theta = i*sigma_y*K`: the `tau_y` channel is purely imaginary, so complex
    /// conjugation `K` forces its multiplying structure factor to be even for the
    /// identity to close consistently, while the real `tau_x`/`tau_z` channels
    /// must be odd. A fully-odd choice for `b_y` would instead force the
    /// Néel-reversal Hall response to be even (not odd) in `M`, incorrectly
    /// making it vanish by symmetry when combined with time-reversal. The
    /// spin-locked chirality of this winding is what makes the crystal Hall
    /// response reverse under Néel reversal (`M -> -M`) while vanishing without
    /// SOC.
    pub fn spin_d_vector(&self, kx: f64, ky: f64, spin: Spin) -> (f64, f64, f64) {
        let sigma = spin.sigma();
        let g = self.gamma(kx, ky);
        let sx = (self.a_lattice * kx).sin();
        let sy = (self.a_lattice * ky).sin();
        let cx = (self.a_lattice * kx).cos();
        let cy = (self.a_lattice * ky).cos();
        let dx = g.re + sigma * self.lambda_soc * sx;
        let dy = -g.im + sigma * self.lambda_soc * (cx + cy);
        let dz = self.eps_a(kx, ky) - sigma * self.exchange_m + sigma * self.lambda_soc * (sx + sy);
        (dx, dy, dz)
    }

    /// Magnitude `|d_sigma(k)|` of the sublattice-space vector.
    #[inline]
    fn d_norm(&self, kx: f64, ky: f64, spin: Spin) -> f64 {
        let (dx, dy, dz) = self.spin_d_vector(kx, ky, spin);
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Explicit 2×2 per-spin Bloch Hamiltonian `H_sigma(k)` as a complex matrix.
    ///
    /// Expands `H_sigma(k) = eps0(k)*tau0 + d_sigma(k).tau` in the Pauli basis
    /// (`tau0` = identity, `tau_x/y/z` the Pauli matrices):
    ///
    /// ```text
    /// H = [ eps0 + dz     dx - i*dy ]
    ///     [ dx + i*dy     eps0 - dz ]
    /// ```
    ///
    /// Built from exactly the same private `eps0` helper and
    /// [`spin_d_vector`](Self::spin_d_vector) used by the closed-form eigenvalue
    /// path ([`band_energy`](Self::band_energy), the private `d_norm` helper,
    /// [`berry_curvature`](Self::berry_curvature)) -- there is exactly one
    /// underlying computation of `(eps0, dx, dy, dz)` in this module, so the
    /// matrix and closed-form representations can never silently diverge from
    /// each other. Because `tau_x, tau_y, tau_z` satisfy the Pauli identity
    /// `(d.tau)^2 = |d|^2 * I`, this matrix's eigenvalues are exactly
    /// `eps0(k) ± |d_sigma(k)|`, matching [`band_energy`](Self::band_energy).
    ///
    /// This is the type used to implement
    /// [`BlochHamiltonian`](crate::altermagnet::bloch_hamiltonian::BlochHamiltonian)
    /// for a single spin sector via
    /// [`AltermagnetSpinHamiltonian`](crate::altermagnet::bloch_hamiltonian::AltermagnetSpinHamiltonian),
    /// and is what
    /// [`KuboBerry`](crate::altermagnet::kubo_berry::KuboBerry) diagonalizes to
    /// numerically cross-check [`berry_curvature`](Self::berry_curvature).
    ///
    /// # Errors
    ///
    /// Propagates the (infallible for `n = 2`) error from
    /// [`CMatrix::from_rows`].
    pub fn hamiltonian_matrix(&self, kx: f64, ky: f64, spin: Spin) -> Result<CMatrix> {
        let e0 = self.eps0(kx, ky);
        let (dx, dy, dz) = self.spin_d_vector(kx, ky, spin);
        CMatrix::from_rows(vec![
            vec![Complex::new(e0 + dz, 0.0), Complex::new(dx, -dy)],
            vec![Complex::new(dx, dy), Complex::new(e0 - dz, 0.0)],
        ])
    }

    /// Rotate momentum `(kx, ky)` by the sublattice-swap angle `pi / n`, where
    /// `n` is [`AltermagneticSymmetry::harmonic_order`].
    ///
    /// This is the correct, general `R` entering the spin-group invariant
    /// `E_up(k) = E_down(R k)` documented in the [module docs](self): rotating
    /// momentum by `pi/n` sends `phi -> phi + pi/n`, so
    /// `eps_a(Rk) = t_am*(a|k|)^2*cos(n*phi + pi - n*phi_neel) = -eps_a(k)` for
    /// **any** harmonic order `n` and **any** `phi_neel` (`cos(x + pi) = -cos(x)`),
    /// while `eps0(Rk) = eps0(k)` because a proper rotation preserves `|k|`.
    /// Combined with `d_sigma = (delta_hyb, 0, eps_a(k) - sigma*M)` at
    /// `lambda_soc = 0`, this gives
    /// `|d_up(k)|^2 = delta_hyb^2 + (eps_a(k)-M)^2 = delta_hyb^2 + (eps_a(Rk)+M)^2 = |d_down(Rk)|^2`,
    /// hence `E_up(k) = E_down(Rk)` exactly.
    ///
    /// For `n = 2` (d-wave) this rotation happens to coincide with reflecting
    /// `kx` and `ky` into each other up to the sign convention used directly by
    /// the `spin_group_relation_dwave` test (both send `phi -> const - phi` or
    /// `phi -> phi + pi/2`, and `cos(2 * (...))` cannot tell the two apart).
    /// For `n = 4` (g-wave, e.g. CrSb) and `n = 6` (i-wave) a literal
    /// `(kx, ky) -> (ky, kx)` coordinate swap does **not** implement this
    /// symmetry (`cos(4*(pi/2 - phi)) = cos(4 phi) = +eps_a(k)/t_am/(ak)^2`,
    /// not `-1` times it) -- the genuine `pi/n` rotation computed here is
    /// required, and is what generalizes the relation to every harmonic order.
    pub fn sublattice_swap_rotation(&self, kx: f64, ky: f64) -> (f64, f64) {
        let alpha = std::f64::consts::PI / f64::from(self.symmetry.harmonic_order());
        let (s, c) = alpha.sin_cos();
        (kx * c - ky * s, kx * s + ky * c)
    }

    /// Energy of a single `(spin, band)` eigenstate at `k` \[eV\].
    ///
    /// `E = eps0(k) ± |d_sigma(k)|` with `+` for [`Band::Upper`] and `-` for
    /// [`Band::Lower`].
    pub fn band_energy(&self, kx: f64, ky: f64, spin: Spin, band: Band) -> f64 {
        let e0 = self.eps0(kx, ky);
        let dn = self.d_norm(kx, ky, spin);
        match band {
            Band::Upper => e0 + dn,
            Band::Lower => e0 - dn,
        }
    }

    /// The `(up, down)` band pair at momentum `k`.
    ///
    /// Each [`SpinBands`] holds the analytic eigenvalues
    /// `eps0(k) ± |d_sigma(k)|` of the corresponding 2×2 spin sector.
    pub fn spin_bands(&self, kx: f64, ky: f64) -> (SpinBands, SpinBands) {
        let e0 = self.eps0(kx, ky);
        let dn_up = self.d_norm(kx, ky, Spin::Up);
        let dn_dn = self.d_norm(kx, ky, Spin::Down);
        let up = SpinBands {
            lower: e0 - dn_up,
            upper: e0 + dn_up,
        };
        let down = SpinBands {
            lower: e0 - dn_dn,
            upper: e0 + dn_dn,
        };
        (up, down)
    }

    /// Spin splitting of the upper band at `k`: `|d_up(k)| - |d_down(k)|` \[eV\].
    ///
    /// This vanishes along the altermagnetic nodal directions returned by
    /// [`AltermagneticSymmetry::node_angles`], where `eps_a(k) = 0`.
    pub fn spin_splitting_at(&self, kx: f64, ky: f64) -> f64 {
        self.d_norm(kx, ky, Spin::Up) - self.d_norm(kx, ky, Spin::Down)
    }

    /// Berry curvature `Omega_sigma^band(k)` of a single `(spin, band)` state.
    ///
    /// Uses the exact closed form for a two-band `d·tau` Hamiltonian,
    /// `Omega = ∓ (1/2) d·(∂_kx d × ∂_ky d) / |d|^3` (`-` for [`Band::Upper`],
    /// `+` for [`Band::Lower`]). The `k`-derivatives of the three-vector `d` are
    /// evaluated by central finite differences with step `h = 1e-6*max(1, |k|)`.
    /// Returns `0` at band-touching points where `|d| = 0`.
    pub fn berry_curvature(&self, kx: f64, ky: f64, spin: Spin, band: Band) -> f64 {
        let k_mag = (kx * kx + ky * ky).sqrt();
        let h = 1e-6 * k_mag.max(1.0);

        let d = self.spin_d_vector(kx, ky, spin);
        let dp_x = self.spin_d_vector(kx + h, ky, spin);
        let dm_x = self.spin_d_vector(kx - h, ky, spin);
        let dp_y = self.spin_d_vector(kx, ky + h, spin);
        let dm_y = self.spin_d_vector(kx, ky - h, spin);

        let inv = 1.0 / (2.0 * h);
        let ddx = (
            (dp_x.0 - dm_x.0) * inv,
            (dp_x.1 - dm_x.1) * inv,
            (dp_x.2 - dm_x.2) * inv,
        );
        let ddy = (
            (dp_y.0 - dm_y.0) * inv,
            (dp_y.1 - dm_y.1) * inv,
            (dp_y.2 - dm_y.2) * inv,
        );

        let cross = (
            ddx.1 * ddy.2 - ddx.2 * ddy.1,
            ddx.2 * ddy.0 - ddx.0 * ddy.2,
            ddx.0 * ddy.1 - ddx.1 * ddy.0,
        );
        let triple = d.0 * cross.0 + d.1 * cross.1 + d.2 * cross.2;
        let dn = (d.0 * d.0 + d.1 * d.1 + d.2 * d.2).sqrt();
        let dn3 = dn * dn * dn;
        if dn3 <= f64::MIN_POSITIVE {
            return 0.0;
        }
        let base = triple / dn3;
        match band {
            Band::Upper => -0.5 * base,
            Band::Lower => 0.5 * base,
        }
    }

    /// Validate the disk-integration grid parameters.
    fn check_grid(n_grid: usize, k_max: f64) -> Result<()> {
        if n_grid == 0 {
            return Err(Error::InvalidParameter {
                param: "n_grid".to_string(),
                reason: "grid resolution must be at least 1".to_string(),
            });
        }
        if k_max <= 0.0 || !k_max.is_finite() {
            return Err(Error::InvalidParameter {
                param: "k_max".to_string(),
                reason: "Brillouin-zone radius must be a finite positive number".to_string(),
            });
        }
        Ok(())
    }

    /// Zero-temperature Fermi-sea Berry sum over a disk of radius `k_max`.
    ///
    /// When `spin_weighted` is false this returns the charge (crystal) Hall
    /// conductivity; when true, each spin contribution is weighted by `sigma`,
    /// giving the spin Hall conductivity.
    fn fermi_sea_berry(&self, n_grid: usize, k_max: f64, spin_weighted: bool) -> Result<f64> {
        Self::check_grid(n_grid, k_max)?;
        let dk = 2.0 * k_max / n_grid as f64;
        let cell_area = dk * dk;
        let prefactor = cell_area / (2.0 * std::f64::consts::PI).powi(2);
        let k_max_sq = k_max * k_max;

        let mut acc = 0.0;
        for i in 0..n_grid {
            let kx = -k_max + (i as f64 + 0.5) * dk;
            for j in 0..n_grid {
                let ky = -k_max + (j as f64 + 0.5) * dk;
                if kx * kx + ky * ky > k_max_sq {
                    continue;
                }
                for spin in [Spin::Up, Spin::Down] {
                    let weight = if spin_weighted { spin.sigma() } else { 1.0 };
                    for band in [Band::Lower, Band::Upper] {
                        let energy = self.band_energy(kx, ky, spin, band);
                        if self.fermi_energy - energy >= 0.0 {
                            acc += weight * self.berry_curvature(kx, ky, spin, band);
                        }
                    }
                }
            }
        }
        Ok(prefactor * acc)
    }

    /// Crystal (charge) Hall conductivity from the occupied Berry curvature.
    ///
    /// Sums Berry curvature over both spins and both bands inside the Fermi disk.
    /// This is identically zero (up to grid resolution) without SOC, and becomes
    /// nonzero and Néel-sign-reversing once `lambda_soc != 0`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] if `n_grid == 0` or `k_max` is not a
    /// finite positive number.
    pub fn crystal_hall_conductivity(&self, n_grid: usize, k_max: f64) -> Result<f64> {
        self.fermi_sea_berry(n_grid, k_max, false)
    }

    /// Spin Hall conductivity: the spin-current Berry sum (`sigma`-weighted).
    ///
    /// Same Fermi-sea integral as [`crystal_hall_conductivity`](Self::crystal_hall_conductivity)
    /// but with the up contribution weighted by `+1` and the down contribution by
    /// `-1`, tracking the spin-current difference rather than the charge sum.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] if `n_grid == 0` or `k_max` is not a
    /// finite positive number.
    pub fn spin_hall_conductivity(&self, n_grid: usize, k_max: f64) -> Result<f64> {
        self.fermi_sea_berry(n_grid, k_max, true)
    }

    /// Net spin polarization of the occupied states inside the Fermi disk.
    ///
    /// Computes `(N_up - N_down)/(N_up + N_down)` over all occupied `(spin, band)`
    /// states. For a compensated altermagnet this is zero within grid resolution.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] if `n_grid == 0` or `k_max` is not a
    /// finite positive number.
    pub fn net_spin_polarization(&self, n_grid: usize, k_max: f64) -> Result<f64> {
        Self::check_grid(n_grid, k_max)?;
        let dk = 2.0 * k_max / n_grid as f64;
        let k_max_sq = k_max * k_max;

        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for i in 0..n_grid {
            let kx = -k_max + (i as f64 + 0.5) * dk;
            for j in 0..n_grid {
                let ky = -k_max + (j as f64 + 0.5) * dk;
                if kx * kx + ky * ky > k_max_sq {
                    continue;
                }
                for spin in [Spin::Up, Spin::Down] {
                    let sigma = spin.sigma();
                    for band in [Band::Lower, Band::Upper] {
                        let energy = self.band_energy(kx, ky, spin, band);
                        if self.fermi_energy - energy >= 0.0 {
                            numerator += sigma;
                            denominator += 1.0;
                        }
                    }
                }
            }
        }
        if denominator == 0.0 {
            return Ok(0.0);
        }
        Ok(numerator / denominator)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (a) Both spin sectors are degenerate when the Néel exchange vanishes.
    #[test]
    fn spin_degenerate_without_exchange() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.0,
            0.7,
            0.0,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let (up, down) = model.spin_bands(0.4, 0.2);
        assert!((up.lower - down.lower).abs() < 1e-12);
        assert!((up.upper - down.upper).abs() < 1e-12);
    }

    /// (b) Ordinary-AFM limit: `t_am = 0` restores spin degeneracy everywhere.
    #[test]
    fn spin_degenerate_afm_limit() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.9,
            0.0,
            0.0,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let (up, down) = model.spin_bands(0.4, 0.2);
        assert!((up.lower - down.lower).abs() < 1e-12);
        assert!((up.upper - down.upper).abs() < 1e-12);
    }

    /// (c) The spin splitting vanishes along the altermagnetic nodal directions.
    #[test]
    fn splitting_vanishes_at_nodes() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.0,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let radius = 0.35;
        for angle in model.symmetry.node_angles() {
            let kx = radius * angle.cos();
            let ky = radius * angle.sin();
            assert!(
                model.spin_splitting_at(kx, ky).abs() < 1e-9,
                "splitting at node angle {angle} should vanish"
            );
        }
        // Off-node the splitting is genuinely nonzero (maximal along phi = 0).
        assert!(model.spin_splitting_at(radius, 0.0).abs() > 1e-3);
    }

    /// (d) Without SOC each spin band is even in `k`: `E(k) == E(-k)`.
    #[test]
    fn bands_even_in_k() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.0,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let (kx, ky) = (0.37, 0.21);
        let (up_p, down_p) = model.spin_bands(kx, ky);
        let (up_m, down_m) = model.spin_bands(-kx, -ky);
        assert!((up_p.lower - up_m.lower).abs() < 1e-12);
        assert!((up_p.upper - up_m.upper).abs() < 1e-12);
        assert!((down_p.lower - down_m.lower).abs() < 1e-12);
        assert!((down_p.upper - down_m.upper).abs() < 1e-12);
    }

    /// (e) d-wave spin-group relation: `E_up(kx, ky) == E_down(ky, kx)`.
    #[test]
    fn spin_group_relation_dwave() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.0,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let (kx, ky) = (0.4, 0.15);
        let (up, _) = model.spin_bands(kx, ky);
        let (_, down_swapped) = model.spin_bands(ky, kx);
        assert!((up.lower - down_swapped.lower).abs() < 1e-12);
        assert!((up.upper - down_swapped.upper).abs() < 1e-12);
    }

    /// (f) The crystal (charge) Hall conductivity is zero without SOC.
    #[test]
    fn crystal_hall_zero_without_soc() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.5,
            0.0,
            0.2,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let hall = model
            .crystal_hall_conductivity(41, 2.5)
            .expect("hall integral");
        assert!(
            hall.abs() < 1e-12,
            "charge Hall must vanish without SOC: {hall}"
        );
    }

    /// (g) With SOC the crystal Hall is nonzero and reverses under Néel reversal.
    #[test]
    fn crystal_hall_flips_under_neel_reversal_with_soc() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.5,
            0.4,
            0.2,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let (n_grid, k_max) = (61, 2.5);
        let hall = model
            .crystal_hall_conductivity(n_grid, k_max)
            .expect("hall integral");
        let mut reversed = model.clone();
        reversed.exchange_m = -model.exchange_m;
        let hall_reversed = reversed
            .crystal_hall_conductivity(n_grid, k_max)
            .expect("hall integral");
        assert!(
            hall.abs() > 1e-9,
            "crystal Hall must be nonzero with SOC: {hall}"
        );
        assert!(
            (hall + hall_reversed).abs() < 1e-6 * (1.0 + hall.abs()),
            "crystal Hall must reverse under Neel reversal: {hall} vs {hall_reversed}"
        );
    }

    /// (h) With SOC on, the d-vector components carry the required mixed k-parity:
    /// `dx`/`dz` pick up an odd-in-k SOC contribution, `dy` picks up an even-in-k
    /// one. Required by the time-reversal-to-Néel-reversal operator identity
    /// documented on `spin_d_vector` — a fully-odd `dy` would be physically wrong.
    #[test]
    fn soc_term_has_correct_mixed_parity() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.4,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let (kx, ky) = (0.4, 0.15);
        let (dx_p, dy_p, dz_p) = model.spin_d_vector(kx, ky, Spin::Up);
        let (dx_m, dy_m, dz_m) = model.spin_d_vector(-kx, -ky, Spin::Up);

        // dx is odd in k: the antisymmetric part (the SOC contribution) is nonzero.
        assert!(
            (dx_p - dx_m).abs() > 1e-6,
            "dx must carry an odd-in-k SOC component: dx(k)={dx_p}, dx(-k)={dx_m}"
        );
        // dy is even in k: the antisymmetric part vanishes, the symmetric part does not.
        assert!(
            (dy_p - dy_m).abs() < 1e-9,
            "dy must be even in k: dy(k)={dy_p}, dy(-k)={dy_m}"
        );
        assert!(
            (dy_p + dy_m).abs() > 1e-6,
            "dy must carry a nonzero even-in-k SOC component: dy(k)={dy_p}, dy(-k)={dy_m}"
        );
        // dz is odd in k: the antisymmetric part (the SOC contribution) is nonzero.
        assert!(
            (dz_p - dz_m).abs() > 1e-6,
            "dz must carry an odd-in-k SOC component: dz(k)={dz_p}, dz(-k)={dz_m}"
        );
    }

    /// (i) SOC breaks the non-relativistic spin-group relation
    /// `E_up(kx, ky) == E_down(ky, kx)` verified at `lambda_soc = 0` in test (e).
    /// This is physically correct, not a bug: that relation is a
    /// non-relativistic-limit symmetry of the collinear (spin-diagonal)
    /// Hamiltonian, and SOC couples spin to the lattice, generically lifting it.
    #[test]
    fn spin_group_relation_broken_by_soc() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.4,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let (kx, ky) = (0.4, 0.15);
        let (up, _) = model.spin_bands(kx, ky);
        let (_, down_swapped) = model.spin_bands(ky, kx);
        assert!(
            (up.lower - down_swapped.lower).abs() > 1e-6,
            "SOC must break the spin-group relation for the lower band: {} vs {}",
            up.lower,
            down_swapped.lower
        );
        assert!(
            (up.upper - down_swapped.upper).abs() > 1e-6,
            "SOC must break the spin-group relation for the upper band: {} vs {}",
            up.upper,
            down_swapped.upper
        );
    }

    /// (j) The Berry curvature vanishes pointwise without SOC and is clearly
    /// nonzero with SOC, isolating the mechanism from the Fermi-sea integral
    /// already checked in test (g).
    #[test]
    fn berry_curvature_pointwise_zero_without_soc_nonzero_with_soc() {
        let base = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.0,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let (kx, ky) = (0.4, 0.15);
        let omega_no_soc = base.berry_curvature(kx, ky, Spin::Up, Band::Lower);
        assert!(
            omega_no_soc.abs() < 1e-9,
            "Berry curvature must vanish pointwise without SOC: {omega_no_soc}"
        );

        let mut with_soc = base.clone();
        with_soc.lambda_soc = 0.4;
        let omega_soc = with_soc.berry_curvature(kx, ky, Spin::Up, Band::Lower);
        assert!(
            omega_soc.abs() > 1e-6,
            "Berry curvature must be clearly nonzero with SOC: {omega_soc}"
        );
    }

    /// (k) `hamiltonian_matrix` is Hermitian at representative `k`-points, both
    /// spins, with SOC on (the general case: `dy != 0`, so off-diagonal terms
    /// are genuinely complex, not just real).
    #[test]
    fn hamiltonian_matrix_is_hermitian() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::GWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.4,
            0.5,
            0.2,
            1.0,
        )
        .expect("parameters are valid");
        for &(kx, ky) in &[(0.4, 0.15), (-0.7, 0.9), (1.1, -0.3), (0.0, 0.0)] {
            for spin in [Spin::Up, Spin::Down] {
                let h = model
                    .hamiltonian_matrix(kx, ky, spin)
                    .expect("2x2 matrix construction cannot fail");
                for i in 0..h.n() {
                    for j in 0..h.n() {
                        let hij = h.get(i, j);
                        let hji_conj = h.get(j, i).conj();
                        assert!(
                            (hij.re - hji_conj.re).abs() < 1e-12
                                && (hij.im - hji_conj.im).abs() < 1e-12,
                            "H[{i}][{j}]={hij:?} != conj(H[{j}][{i}])={hji_conj:?} at k=({kx},{ky}), spin={spin:?}"
                        );
                    }
                }
            }
        }
    }

    /// (l) `hamiltonian_matrix`'s eigenvalues match the closed-form `band_energy`
    /// exactly (to Jacobi-eigensolver tolerance), with and without SOC. This is
    /// the "never silently diverge" guarantee promised by `hamiltonian_matrix`'s
    /// docs: both paths are built from the same `eps0`/`spin_d_vector` call, so
    /// this test would catch any future edit that breaks that invariant.
    #[test]
    fn hamiltonian_matrix_eigenvalues_match_band_energy() {
        for lambda_soc in [0.0, 0.4] {
            let model = AltermagnetBandModel::new(
                AltermagneticSymmetry::GWave,
                1.0,
                0.3,
                0.8,
                0.6,
                lambda_soc,
                0.5,
                0.2,
                1.0,
            )
            .expect("parameters are valid");
            for &(kx, ky) in &[(0.4, 0.15), (-0.7, 0.9), (1.1, -0.3)] {
                for spin in [Spin::Up, Spin::Down] {
                    let h = model.hamiltonian_matrix(kx, ky, spin).expect("valid matrix");
                    let (evals, _) = h
                        .hermitian_eigendecomposition()
                        .expect("2x2 Hermitian eigendecomposition always converges");
                    let expect_lower = model.band_energy(kx, ky, spin, Band::Lower);
                    let expect_upper = model.band_energy(kx, ky, spin, Band::Upper);
                    assert!(
                        (evals[0] - expect_lower).abs() < 1e-9,
                        "lower eigenvalue {} != band_energy {} (soc={lambda_soc}, k=({kx},{ky}), spin={spin:?})",
                        evals[0], expect_lower
                    );
                    assert!(
                        (evals[1] - expect_upper).abs() < 1e-9,
                        "upper eigenvalue {} != band_energy {} (soc={lambda_soc}, k=({kx},{ky}), spin={spin:?})",
                        evals[1], expect_upper
                    );
                }
            }
        }
    }

    /// (m) Band-sum-to-zero: for a fixed spin, the Berry curvature of the two
    /// bands cancels exactly, `Omega_Upper(k) + Omega_Lower(k) = 0`. This holds
    /// by construction (`berry_curvature` computes a single `base` value and
    /// returns `+/- 0.5*base` for Lower/Upper), but is worth a standing
    /// regression test since it is the "Band-sum-to-zero" invariant used
    /// elsewhere (e.g. the Fermi-sea integrals summing over both bands).
    #[test]
    fn berry_curvature_band_sum_to_zero() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::GWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.4,
            0.5,
            0.2,
            1.0,
        )
        .expect("parameters are valid");
        for &(kx, ky) in &[(0.4, 0.15), (-0.7, 0.9), (1.1, -0.3)] {
            for spin in [Spin::Up, Spin::Down] {
                let lower = model.berry_curvature(kx, ky, spin, Band::Lower);
                let upper = model.berry_curvature(kx, ky, spin, Band::Upper);
                assert!(
                    (lower + upper).abs() < 1e-10,
                    "Omega_Lower + Omega_Upper = {} (expected 0) at k=({kx},{ky}), spin={spin:?}",
                    lower + upper
                );
            }
        }
    }

    /// (n) Berry curvature oddness under the combined spin + momentum + Neel-order
    /// inversion: `Omega_sigma(k; M) = -Omega_{-sigma}(-k; -M)`.
    ///
    /// This is the rigorously correct general form of Berry-curvature "oddness"
    /// for this model: a literal same-spin `Omega_sigma(k) = -Omega_sigma(-k)`
    /// (holding `M` fixed) does **not** hold in general once both `delta_hyb`
    /// (a k-independent hybridization offset shared by both spins) and SOC are
    /// nonzero simultaneously -- verified numerically to fail by an amount
    /// comparable to the curvature itself, not floating-point noise. The
    /// combined relation tested here, in contrast, follows directly from the
    /// `Theta = i*tau_y*K` operator identity already documented on
    /// `spin_d_vector` (`Theta H(k;M) Theta^-1 = H(-k;-M)`) once the implicit
    /// spin-relabelling identity `H_sigma(k;M) = H_{-sigma}(k;-M) . conj` is
    /// folded in: `H_{-sigma}(-k;-M)` equals the complex conjugate of
    /// `H_sigma(k;M)` (same `eps0`, `dx`, `dz`; `dy` flips sign), and Berry
    /// curvature flips sign under this antiunitary-like conjugation. Confirmed
    /// numerically to finite-difference precision (~1e-11) for generic
    /// (delta_hyb, exchange_m, t_am, lambda_soc, neel_angle) and multiple k.
    #[test]
    fn berry_curvature_odd_under_combined_spin_k_neel_flip() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.4,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        let mut reversed = model.clone();
        reversed.exchange_m = -model.exchange_m;

        for &(kx, ky) in &[(0.4, 0.15), (0.9, -0.6), (-0.2, 1.3)] {
            for band in [Band::Lower, Band::Upper] {
                let omega_up_k_m = model.berry_curvature(kx, ky, Spin::Up, band);
                let omega_down_mk_negm = reversed.berry_curvature(-kx, -ky, Spin::Down, band);
                assert!(
                    (omega_up_k_m + omega_down_mk_negm).abs() < 1e-8,
                    "Omega_Up(k;M) + Omega_Down(-k;-M) = {} (expected 0) at k=({kx},{ky}), band={band:?}",
                    omega_up_k_m + omega_down_mk_negm
                );

                let omega_down_k_m = model.berry_curvature(kx, ky, Spin::Down, band);
                let omega_up_mk_negm = reversed.berry_curvature(-kx, -ky, Spin::Up, band);
                assert!(
                    (omega_down_k_m + omega_up_mk_negm).abs() < 1e-8,
                    "Omega_Down(k;M) + Omega_Up(-k;-M) = {} (expected 0) at k=({kx},{ky}), band={band:?}",
                    omega_down_k_m + omega_up_mk_negm
                );
            }
        }
    }

    /// (o) `net_spin_polarization` is a real, asserted-in-a-test invariant (not
    /// merely printed by the `altermagnet_band_structure` example): both the
    /// MnTe (d-wave) and CrSb (g-wave) presets are compensated (zero net spin
    /// polarization of the occupied Fermi sea) at their default parameters.
    #[test]
    fn net_spin_polarization_vanishes_for_presets() {
        for model in [AltermagnetBandModel::mnte(), AltermagnetBandModel::crsb()] {
            let k_max = 0.8 / model.a_lattice;
            let polarization = model
                .net_spin_polarization(41, k_max)
                .expect("valid grid parameters");
            assert!(
                polarization.abs() < 1e-9,
                "{} preset should be compensated (zero net spin polarization), got {polarization}",
                model.symmetry
            );
        }
    }

    /// (p) Generalized spin-group relation for g-wave (CrSb, harmonic order 4):
    /// rotating momentum by the sublattice-swap angle `pi/n` swaps the up/down
    /// band energies, `E_up(k) = E_down(R k)`. Unlike the `n = 2` (d-wave) case
    /// in test (e), a literal `(kx, ky) -> (ky, kx)` coordinate swap does *not*
    /// implement this symmetry for `n = 4` (see `sublattice_swap_rotation`'s
    /// docs) -- the genuine `pi/n` rotation is required, and is what generalizes
    /// the relation to every harmonic order.
    #[test]
    fn spin_group_relation_gwave_generalized_rotation() {
        let model = AltermagnetBandModel::new(
            AltermagneticSymmetry::GWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.0,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid");
        assert_eq!(model.symmetry.harmonic_order(), 4);

        let (kx, ky) = (0.4, 0.15);
        let (rkx, rky) = model.sublattice_swap_rotation(kx, ky);
        let (up, _) = model.spin_bands(kx, ky);
        let (_, down_rotated) = model.spin_bands(rkx, rky);
        assert!(
            (up.lower - down_rotated.lower).abs() < 1e-9,
            "E_up,lower(k)={} != E_down,lower(Rk)={}",
            up.lower,
            down_rotated.lower
        );
        assert!(
            (up.upper - down_rotated.upper).abs() < 1e-9,
            "E_up,upper(k)={} != E_down,upper(Rk)={}",
            up.upper,
            down_rotated.upper
        );

        // A literal coordinate swap is NOT the correct rotation for g-wave: it
        // must fail to reproduce the relation (this documents the pitfall the
        // general `sublattice_swap_rotation` avoids).
        let (_, down_swapped) = model.spin_bands(ky, kx);
        assert!(
            (up.lower - down_swapped.lower).abs() > 1e-3,
            "a literal (kx,ky)->(ky,kx) swap should NOT satisfy the g-wave spin-group relation"
        );
    }
}
