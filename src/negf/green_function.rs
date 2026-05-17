//! Green's function and transport machinery for 1D mesoscopic conductors.
//!
//! # Physical Background
//!
//! The retarded Green's function is defined by
//!
//! ```text
//! G^R(E) = [ (E + iη)I − H − Σ_L^R − Σ_R^R ]^{−1}
//! ```
//!
//! The Landauer–Büttiker transmission coefficient is
//!
//! ```text
//! T(E) = Tr[ Γ_L G^R Γ_R G^A ]
//! ```
//!
//! and the current
//!
//! ```text
//! I = (e/h) ∫ T(E) [f_L(E) − f_R(E)] dE
//! ```
//!
//! # References
//!
//! - S. Datta, *Electronic Transport in Mesoscopic Systems* (Cambridge, 1995)
//! - M. P. López Sancho et al., J. Phys. F **15**, 851 (1985)

use crate::constants::{CONDUCTANCE_QUANTUM, E_CHARGE, H_PLANCK, KB};
use crate::error::{self, Result};
use crate::material::Xorshift64;
use crate::math::{CMatrix, Complex};

// ============================================================================
// Hamiltonian1D
// ============================================================================

/// Tight-binding Hamiltonian for a 1D chain of `n_sites` lattice sites.
///
/// The matrix representation is tridiagonal:
/// `H[i][i] = onsite[i]`,  `H[i][i±1] = hopping` (real).
#[derive(Debug, Clone)]
pub struct Hamiltonian1D {
    /// On-site energies \[eV\].
    pub onsite: Vec<f64>,
    /// Nearest-neighbour hopping parameter \[eV\].
    pub hopping: f64,
    /// Number of lattice sites (1 – 64).
    pub n_sites: usize,
}

impl Hamiltonian1D {
    /// Create a Hamiltonian from explicit on-site energies and a uniform hopping.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `n_sites` is 0 or > 64, or if `hopping` is not finite.
    pub fn new(onsite: Vec<f64>, hopping: f64) -> Result<Self> {
        let n = onsite.len();
        if n == 0 {
            return Err(error::invalid_param(
                "onsite",
                "must have at least one site",
            ));
        }
        if n > 64 {
            return Err(error::invalid_param("onsite", "n_sites must be ≤ 64"));
        }
        if !hopping.is_finite() {
            return Err(error::invalid_param("hopping", "must be a finite number"));
        }
        Ok(Self {
            onsite,
            hopping,
            n_sites: n,
        })
    }

    /// Construct a uniform chain: all on-site energies = `e0`, hopping = `t`.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`new`](Self::new).
    pub fn from_uniform(n_sites: usize, e0: f64, t: f64) -> Result<Self> {
        Self::new(vec![e0; n_sites], t)
    }

    /// Add Gaussian disorder to the on-site energies using a seeded `Xorshift64` PRNG.
    ///
    /// Each on-site energy is perturbed by an independent draw from N(0, σ²).
    ///
    /// # Arguments
    ///
    /// * `sigma` – disorder strength \[eV\]
    /// * `seed`  – PRNG seed (deterministic, reproducible)
    pub fn with_disorder(mut self, sigma: f64, seed: u64) -> Self {
        let mut rng = Xorshift64::new(seed);
        for e in self.onsite.iter_mut() {
            *e += sigma * rng.next_normal();
        }
        self
    }

    /// Convert the Hamiltonian to an N×N complex matrix.
    ///
    /// The result is the tridiagonal Hermitian matrix with
    /// `H[i][i] = onsite[i]` and `H[i][i±1] = hopping` (all real).
    pub fn to_matrix(&self) -> CMatrix {
        let n = self.n_sites;
        let mut m = CMatrix::zeros(n);
        for i in 0..n {
            m.set(i, i, Complex::from_real(self.onsite[i]));
            if i + 1 < n {
                let t = Complex::from_real(self.hopping);
                m.set(i, i + 1, t);
                m.set(i + 1, i, t);
            }
        }
        m
    }

    /// Return the on-site energy at site `idx` \[eV\].
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `idx ≥ n_sites`.
    pub fn site_energy(&self, idx: usize) -> Result<f64> {
        if idx >= self.n_sites {
            return Err(error::invalid_param("idx", "site index out of bounds"));
        }
        Ok(self.onsite[idx])
    }

    /// Return the number of lattice sites.
    pub fn n_sites(&self) -> usize {
        self.n_sites
    }
}

// ============================================================================
// LeadSelfEnergy
// ============================================================================

/// Wide-band-limit self-energy for a semi-infinite lead.
///
/// In the wide-band limit the retarded self-energy is purely imaginary and
/// energy-independent:
///
/// ```text
/// Σ^R = −iΓ/2
/// ```
///
/// where Γ is the level-broadening width \[eV\].
#[derive(Debug, Clone)]
pub struct LeadSelfEnergy {
    /// Level-broadening width Γ \[eV\].
    pub gamma: f64,
    /// Fermi level of the lead \[eV\].
    pub e_fermi: f64,
}

impl LeadSelfEnergy {
    /// Create a new lead self-energy.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `gamma ≤ 0`.
    pub fn new(gamma: f64, e_fermi: f64) -> Result<Self> {
        if gamma <= 0.0 {
            return Err(error::invalid_param(
                "gamma",
                "coupling width must be positive",
            ));
        }
        Ok(Self { gamma, e_fermi })
    }

    /// Construct with `e_fermi = 0` (symmetric bias reference).
    pub fn wide_band(gamma: f64) -> Self {
        Self {
            gamma,
            e_fermi: 0.0,
        }
    }

    /// Retarded self-energy Σ^R = −iΓ/2 (wide-band limit, energy-independent).
    ///
    /// Returns a purely imaginary `Complex`.
    pub fn sigma_r(&self, _energy: f64) -> Complex {
        Complex::new(0.0, -self.gamma / 2.0)
    }

    /// Advanced self-energy Σ^A = (Σ^R)* = +iΓ/2.
    pub fn sigma_a(&self, energy: f64) -> Complex {
        self.sigma_r(energy).conj()
    }

    /// Coupling matrix Γ = i(Σ^R − Σ^A) = γ (scalar in the wide-band limit) \[eV\].
    pub fn gamma_matrix(&self) -> f64 {
        self.gamma
    }
}

// ============================================================================
// SanchoRubio
// ============================================================================

/// Iterative surface Green's function for a semi-infinite 1D tight-binding lead.
///
/// Uses the Sancho–Rubio decimation algorithm (López Sancho et al., J. Phys. F 15, 851, 1985)
/// to evaluate the surface Green's function of a semi-infinite 1D chain.
#[derive(Debug, Clone)]
pub struct SanchoRubio {
    /// Inter-cell hopping parameter in the lead \[eV\].
    pub t_lead: f64,
    /// On-site energy in the lead \[eV\].
    pub e_onsite: f64,
    /// Maximum number of decimation iterations.
    pub max_iter: usize,
    /// Convergence tolerance for the transfer matrix norm.
    pub tol: f64,
}

impl SanchoRubio {
    /// Create a `SanchoRubio` solver with default iteration limit (100) and tolerance (1e-8).
    pub fn new(t_lead: f64, e_onsite: f64) -> Self {
        Self {
            t_lead,
            e_onsite,
            max_iter: 100,
            tol: 1e-8,
        }
    }

    /// Compute the surface Green's function at complex energy `E + iη`.
    ///
    /// The scalar Sancho–Rubio iteration for a 1D chain:
    ///
    /// ```text
    /// g_n  = 1 / (z − ε_n)
    /// ε_s_{n+1} = ε_s_n + α_n · g_n · β_n
    /// ε_{n+1}   = ε_n + α_n·g_n·β_n + β_n·g_n·α_n
    /// α_{n+1}   = α_n · g_n · α_n
    /// β_{n+1}   = β_n · g_n · β_n
    /// ```
    ///
    /// Terminates when |α| < `tol`.  Returns `g_s = 1 / (z − ε_s_final)`.
    ///
    /// # Errors
    ///
    /// Returns `NumericalError` if the iteration does not converge.
    pub fn surface_gf(&self, energy: f64, eta: f64) -> Result<Complex> {
        let z = Complex::new(energy, eta); // complex energy
        let mut epsilon_s = Complex::from_real(self.e_onsite);
        let mut epsilon = Complex::from_real(self.e_onsite);
        let mut alpha = Complex::from_real(self.t_lead);
        let mut beta = Complex::from_real(self.t_lead);

        for _iter in 0..self.max_iter {
            // g_n = 1 / (z - epsilon_n)
            let denom = z.sub(&epsilon);
            let g_n = Complex::ONE.div(&denom);

            // ε_s_{n+1} = ε_s_n + α_n * g_n * β_n
            let alpha_g = alpha.mul(&g_n);
            let beta_g = beta.mul(&g_n);
            let delta_s = alpha_g.mul(&beta);
            epsilon_s = epsilon_s.add(&delta_s);

            // ε_{n+1} = ε_n + α_n*g_n*β_n + β_n*g_n*α_n
            let delta_bulk = alpha_g.mul(&beta).add(&beta_g.mul(&alpha));
            epsilon = epsilon.add(&delta_bulk);

            // α_{n+1} = α_n * g_n * α_n
            let alpha_new = alpha_g.mul(&alpha);
            // β_{n+1} = β_n * g_n * β_n
            let beta_new = beta_g.mul(&beta);

            alpha = alpha_new;
            beta = beta_new;

            if alpha.norm() < self.tol {
                // Converged: compute surface GF
                let gs_denom = z.sub(&epsilon_s);
                return Ok(Complex::ONE.div(&gs_denom));
            }
        }

        Err(error::numerical_error(
            "SanchoRubio: surface GF iteration did not converge",
        ))
    }
}

// ============================================================================
// GreenFunction
// ============================================================================

/// Retarded/advanced Green's function for a finite 1D tight-binding chain
/// coupled to two wide-band leads.
///
/// The retarded Green's function is:
///
/// ```text
/// G^R(E) = [ (E + iη)I − H − Σ_L^R·|0><0| − Σ_R^R·|N-1><N-1| ]^{-1}
/// ```
#[derive(Debug, Clone)]
pub struct GreenFunction {
    /// Tight-binding Hamiltonian of the scattering region.
    pub hamiltonian: Hamiltonian1D,
    /// Self-energy of the left lead.
    pub sigma_l: LeadSelfEnergy,
    /// Self-energy of the right lead.
    pub sigma_r: LeadSelfEnergy,
    /// Broadening parameter η \[eV\] (small positive number for numerical stability).
    pub eta: f64,
}

impl GreenFunction {
    /// Create a new `GreenFunction`.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `eta < 0`.
    pub fn new(h: Hamiltonian1D, sl: LeadSelfEnergy, sr: LeadSelfEnergy, eta: f64) -> Result<Self> {
        if eta < 0.0 {
            return Err(error::invalid_param(
                "eta",
                "broadening must be non-negative",
            ));
        }
        Ok(Self {
            hamiltonian: h,
            sigma_l: sl,
            sigma_r: sr,
            eta,
        })
    }

    /// Build the effective Hamiltonian matrix `(E + iη)I − H − Σ_L − Σ_R`.
    fn effective_hamiltonian(&self, energy: f64) -> CMatrix {
        let n = self.hamiltonian.n_sites;
        let h = self.hamiltonian.to_matrix();
        // (E + iη)I
        let e_plus_eta = Complex::new(energy, self.eta);
        let mut eff = CMatrix::zeros(n);
        for i in 0..n {
            for j in 0..n {
                let hij = h.get(i, j);
                let diag_contrib = if i == j { e_plus_eta } else { Complex::ZERO };
                eff.set(i, j, diag_contrib.sub(&hij));
            }
        }
        // Subtract Σ_L at [0][0]
        let sl = self.sigma_l.sigma_r(energy);
        let cur00 = eff.get(0, 0);
        eff.set(0, 0, cur00.sub(&sl));

        // Subtract Σ_R at [N-1][N-1]
        let nm1 = n - 1;
        let sr = self.sigma_r.sigma_r(energy);
        let cur_nn = eff.get(nm1, nm1);
        eff.set(nm1, nm1, cur_nn.sub(&sr));
        eff
    }

    /// Retarded Green's function G^R(E) = [(E+iη)I − H − Σ_L − Σ_R]^{−1}.
    ///
    /// # Errors
    ///
    /// Returns `NumericalError` if the matrix inversion fails (singular matrix).
    pub fn g_retarded(&self, energy: f64) -> Result<CMatrix> {
        let eff = self.effective_hamiltonian(energy);
        eff.inverse()
    }

    /// Advanced Green's function G^A(E) = [G^R(E)]†.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`g_retarded`](Self::g_retarded).
    pub fn g_advanced(&self, energy: f64) -> Result<CMatrix> {
        let gr = self.g_retarded(energy)?;
        Ok(gr.conj_transpose())
    }

    /// Coupling matrix Γ_L (diagonal, non-zero only at [0][0]).
    ///
    /// Γ_L = i(Σ_L^R − Σ_L^A) = γ_L  (wide-band limit) \[eV\].
    pub fn gamma_l_matrix(&self) -> CMatrix {
        let n = self.hamiltonian.n_sites;
        let mut m = CMatrix::zeros(n);
        m.set(0, 0, Complex::from_real(self.sigma_l.gamma_matrix()));
        m
    }

    /// Coupling matrix Γ_R (diagonal, non-zero only at [N-1][N-1]).
    ///
    /// Γ_R = i(Σ_R^R − Σ_R^A) = γ_R \[eV\].
    pub fn gamma_r_matrix(&self) -> CMatrix {
        let n = self.hamiltonian.n_sites;
        let nm1 = n - 1;
        let mut m = CMatrix::zeros(n);
        m.set(nm1, nm1, Complex::from_real(self.sigma_r.gamma_matrix()));
        m
    }

    /// Spectral function A(E) = i[G^R(E) − G^A(E)].
    ///
    /// # Errors
    ///
    /// Propagates errors from [`g_retarded`](Self::g_retarded).
    pub fn spectral_function(&self, energy: f64) -> Result<CMatrix> {
        let gr = self.g_retarded(energy)?;
        let ga = gr.conj_transpose();
        let diff = gr.sub(&ga)?;
        // A = i * (G^R - G^A)
        let n = self.hamiltonian.n_sites;
        let mut a = CMatrix::zeros(n);
        for i in 0..n {
            for j in 0..n {
                // multiply by i: (a+bi)*i = -b + ai
                let v = diff.get(i, j).mul_i();
                a.set(i, j, v);
            }
        }
        Ok(a)
    }

    /// Local density of states at site `site` \[eV^{-1}\].
    ///
    /// LDOS(E, i) = A(E)[i][i] / (2π)
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `site ≥ n_sites`, or propagates `spectral_function` errors.
    pub fn ldos(&self, energy: f64, site: usize) -> Result<f64> {
        if site >= self.hamiltonian.n_sites {
            return Err(error::invalid_param("site", "site index out of bounds"));
        }
        let a = self.spectral_function(energy)?;
        Ok(a.get(site, site).re / (2.0 * std::f64::consts::PI))
    }

    /// Total density of states DOS(E) = Tr[A(E)] / (2π) \[eV^{-1}\].
    ///
    /// # Errors
    ///
    /// Propagates errors from [`spectral_function`](Self::spectral_function).
    pub fn density_of_states(&self, energy: f64) -> Result<f64> {
        let a = self.spectral_function(energy)?;
        Ok(a.trace().re / (2.0 * std::f64::consts::PI))
    }

    /// Lesser Green's function G^<(E) = G^R · Σ^< · G^A.
    ///
    /// In the wide-band limit:
    /// ```text
    /// Σ^< = i·f_L·Γ_L (at [0,0]) + i·f_R·Γ_R (at [N-1,N-1])
    /// ```
    ///
    /// # Errors
    ///
    /// Propagates errors from [`g_retarded`](Self::g_retarded) or matrix arithmetic.
    pub fn g_lesser(&self, energy: f64, mu_l: f64, mu_r: f64, temperature: f64) -> Result<CMatrix> {
        let gr = self.g_retarded(energy)?;
        let ga = gr.conj_transpose();
        let n = self.hamiltonian.n_sites;
        let nm1 = n - 1;

        let f_l = TransportCalculator::fermi_dirac(energy, mu_l, temperature);
        let f_r = TransportCalculator::fermi_dirac(energy, mu_r, temperature);

        // Build Σ^< as sparse matrix (only two non-zero entries)
        let mut sigma_lesser = CMatrix::zeros(n);
        // Σ^<[0][0]     = i * f_L * Γ_L
        let sl_val = Complex::new(0.0, f_l * self.sigma_l.gamma_matrix());
        sigma_lesser.set(0, 0, sl_val);
        // Σ^<[N-1][N-1] = i * f_R * Γ_R
        let sr_val = Complex::new(0.0, f_r * self.sigma_r.gamma_matrix());
        sigma_lesser.set(nm1, nm1, sr_val);

        // G^< = G^R · Σ^< · G^A
        let tmp = gr.matmul(&sigma_lesser)?;
        tmp.matmul(&ga)
    }
}

// ============================================================================
// TransportCalculator
// ============================================================================

/// Landauer–Büttiker transport calculator for a 1D mesoscopic conductor.
///
/// Integrates the transmission function over an energy window to compute
/// current, differential conductance, and zero-bias conductance.
#[derive(Debug, Clone)]
pub struct TransportCalculator {
    /// Green's function object for the scattering region.
    pub gf: GreenFunction,
    /// Number of energy points in the integration grid.
    pub n_energy: usize,
    /// Lower bound of the energy integration window \[eV\].
    pub e_min: f64,
    /// Upper bound of the energy integration window \[eV\].
    pub e_max: f64,
}

impl TransportCalculator {
    /// Create a new `TransportCalculator`.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `e_min ≥ e_max` or `n_energy < 2`.
    pub fn new(gf: GreenFunction, e_min: f64, e_max: f64, n_energy: usize) -> Result<Self> {
        if e_min >= e_max {
            return Err(error::invalid_param("e_min/e_max", "e_min must be < e_max"));
        }
        if n_energy < 2 {
            return Err(error::invalid_param("n_energy", "must be at least 2"));
        }
        Ok(Self {
            gf,
            n_energy,
            e_min,
            e_max,
        })
    }

    /// Transmission coefficient T(E) at a single energy (Landauer formula).
    ///
    /// For the wide-band 1D case:
    /// ```text
    /// T(E) = Γ_L · Γ_R · |G^R[0][N-1]|²
    /// ```
    ///
    /// This is equivalent to Tr[Γ_L G^R Γ_R G^A] for the corner-coupling geometry.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`GreenFunction::g_retarded`].
    pub fn transmission(&self, energy: f64) -> Result<f64> {
        let gr = self.gf.g_retarded(energy)?;
        let n = self.gf.hamiltonian.n_sites;
        let nm1 = n - 1;
        // G^R[0][N-1]: transmission amplitude across the chain
        let g0n = gr.get(0, nm1);
        let gamma_l = self.gf.sigma_l.gamma_matrix();
        let gamma_r = self.gf.sigma_r.gamma_matrix();
        let t = gamma_l * gamma_r * g0n.norm_sq();
        // Clamp to [0, 1] for numerical safety (should be satisfied physically)
        Ok(t.clamp(0.0, 1.0 + 1e-10))
    }

    /// Compute the transmission curve over the energy grid.
    ///
    /// Returns `n_energy` pairs `(E, T(E))` evenly spaced in `[e_min, e_max]`.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`transmission`](Self::transmission).
    pub fn transmission_curve(&self) -> Result<Vec<(f64, f64)>> {
        let de = (self.e_max - self.e_min) / (self.n_energy - 1) as f64;
        let mut curve = Vec::with_capacity(self.n_energy);
        for k in 0..self.n_energy {
            let e = self.e_min + k as f64 * de;
            let t = self.transmission(e)?;
            curve.push((e, t));
        }
        Ok(curve)
    }

    /// Current I(V) via the Landauer formula (trapezoid integration) \[A\].
    ///
    /// ```text
    /// I = (e/h) ∫ T(E) [f_L(E) − f_R(E)] dE
    /// ```
    ///
    /// Symmetric bias: μ_L = +V/2, μ_R = −V/2.
    /// Positive V_bias → left lead at higher chemical potential → positive current.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`transmission`](Self::transmission).
    pub fn current(&self, v_bias: f64, temperature: f64) -> Result<f64> {
        let mu_l = v_bias / 2.0;
        let mu_r = -v_bias / 2.0;
        let de = (self.e_max - self.e_min) / (self.n_energy - 1) as f64;
        let prefactor = E_CHARGE / H_PLANCK;

        let mut values = Vec::with_capacity(self.n_energy);
        for k in 0..self.n_energy {
            let e = self.e_min + k as f64 * de;
            let t = self.transmission(e)?;
            let fl = Self::fermi_dirac(e, mu_l, temperature);
            let fr = Self::fermi_dirac(e, mu_r, temperature);
            values.push(t * (fl - fr));
        }

        // Trapezoid rule
        let mut integral = 0.0;
        for k in 0..(self.n_energy - 1) {
            integral += 0.5 * (values[k] + values[k + 1]) * de;
        }
        Ok(prefactor * integral)
    }

    /// Differential conductance dI/dV at `v_bias` by central difference \[S\].
    ///
    /// Uses step `dV = 1e-4 · max(|v_bias|, 0.001)`.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`current`](Self::current).
    pub fn differential_conductance(&self, v_bias: f64, temperature: f64) -> Result<f64> {
        let dv = 1e-4 * v_bias.abs().max(0.001);
        let i_plus = self.current(v_bias + dv, temperature)?;
        let i_minus = self.current(v_bias - dv, temperature)?;
        Ok((i_plus - i_minus) / (2.0 * dv))
    }

    /// Zero-bias conductance G(V→0) = G₀ · T(E_F) \[S\].
    ///
    /// Uses the Fermi level of the left lead as the reference energy.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`transmission`](Self::transmission).
    pub fn zero_bias_conductance(&self, _temperature: f64) -> Result<f64> {
        let e_f = self.gf.sigma_l.e_fermi;
        let t = self.transmission(e_f)?;
        Ok(CONDUCTANCE_QUANTUM * t)
    }

    /// Fermi–Dirac distribution f(E, μ, T) = 1 / (exp((E−μ)/(k_B T)) + 1).
    ///
    /// Clipped to avoid numerical overflow for large arguments.
    /// At T = 0 (or very small T), uses a step-function approximation.
    pub fn fermi_dirac(energy: f64, mu: f64, temperature: f64) -> f64 {
        if temperature < 1e-10 {
            // T → 0 limit: step function
            return if energy < mu {
                1.0
            } else if energy > mu {
                0.0
            } else {
                0.5
            };
        }
        let x = (energy - mu) / (KB * temperature / E_CHARGE); // dimensionless (energies in eV, KB*T in J → convert)
                                                               // KB [J/K] * T [K] / E_CHARGE [C] = KB*T in eV
                                                               // x = (E[eV] - mu[eV]) / (KB*T [eV])
        let x_clamped = x.clamp(-500.0, 500.0);
        1.0 / (x_clamped.exp() + 1.0)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_gf(n: usize) -> GreenFunction {
        let h = Hamiltonian1D::from_uniform(n, 0.0, 1.0).expect("valid");
        let sl = LeadSelfEnergy::new(0.5, 0.0).expect("valid");
        let sr = LeadSelfEnergy::new(0.5, 0.0).expect("valid");
        GreenFunction::new(h, sl, sr, 1e-3).expect("valid")
    }

    fn make_tc(n: usize) -> TransportCalculator {
        let gf = make_simple_gf(n);
        TransportCalculator::new(gf, -3.0, 3.0, 100).expect("valid")
    }

    // ------------------------------------------------------------------ basics

    #[test]
    fn test_from_uniform_size() {
        let h = Hamiltonian1D::from_uniform(5, 0.0, 1.0).expect("valid");
        assert_eq!(h.n_sites, 5);
        assert_eq!(h.onsite.len(), 5);
        let m = h.to_matrix();
        assert_eq!(m.n(), 5);
    }

    #[test]
    fn test_with_disorder_seeded_reproducible() {
        let h1 = Hamiltonian1D::from_uniform(8, 0.0, 1.0)
            .expect("valid")
            .with_disorder(0.1, 42);
        let h2 = Hamiltonian1D::from_uniform(8, 0.0, 1.0)
            .expect("valid")
            .with_disorder(0.1, 42);
        for i in 0..8 {
            assert!((h1.onsite[i] - h2.onsite[i]).abs() < 1e-15);
        }
        // Different seed gives different result
        let h3 = Hamiltonian1D::from_uniform(8, 0.0, 1.0)
            .expect("valid")
            .with_disorder(0.1, 99);
        let mut differs = false;
        for i in 0..8 {
            if (h1.onsite[i] - h3.onsite[i]).abs() > 1e-12 {
                differs = true;
                break;
            }
        }
        assert!(differs, "different seeds should produce different disorder");
    }

    // --------------------------------------------------------- self-energy

    #[test]
    fn test_sigma_r_purely_imaginary() {
        let se = LeadSelfEnergy::new(0.3, 0.0).expect("valid");
        let sr = se.sigma_r(0.5);
        // Real part should be zero
        assert!(sr.re.abs() < 1e-15);
        // Im part should be -gamma/2
        assert!((sr.im + 0.15).abs() < 1e-15);
    }

    #[test]
    fn test_gamma_matrix_positive() {
        let se = LeadSelfEnergy::new(0.4, 0.0).expect("valid");
        assert!(se.gamma_matrix() > 0.0);
    }

    // --------------------------------------------------------- Sancho–Rubio

    #[test]
    fn test_surface_gf_converges() {
        let sr = SanchoRubio::new(1.0, 0.0);
        // Energy in the band and slightly above it
        let g = sr.surface_gf(1.5, 1e-3).expect("converges");
        assert!(g.is_finite());
    }

    // --------------------------------------------------------- G^R consistency

    #[test]
    fn test_g_retarded_inverse_consistency() {
        // Check that G^R · [(E+iη)I - H - Σ] ≈ I  (Frobenius norm of deviation)
        let n = 4;
        let gf = make_simple_gf(n);
        let energy = 0.5;
        let gr = gf.g_retarded(energy).expect("invertible");
        let eff = gf.effective_hamiltonian(energy);
        let prod = gr.matmul(&eff).expect("matmul ok");
        let eye = CMatrix::eye(n);
        let diff = prod.sub(&eye).expect("sub ok");
        let frob = diff.frobenius_norm();
        assert!(frob < 1e-10, "Frobenius deviation = {}", frob);
    }

    #[test]
    fn test_g_advanced_is_conj_transpose_of_g_retarded() {
        let gf = make_simple_gf(3);
        let gr = gf.g_retarded(0.2).expect("ok");
        let ga = gf.g_advanced(0.2).expect("ok");
        let ga_expected = gr.conj_transpose();
        for i in 0..3 {
            for j in 0..3 {
                let d = ga.get(i, j).sub(&ga_expected.get(i, j));
                assert!(d.norm() < 1e-14, "G^A != (G^R)† at ({}, {})", i, j);
            }
        }
    }

    // --------------------------------------------------------- spectral function

    #[test]
    fn test_spectral_function_is_hermitian() {
        // A = i(G^R - G^A) must satisfy A = A†
        let gf = make_simple_gf(4);
        let a = gf.spectral_function(0.0).expect("ok");
        let n = 4;
        for i in 0..n {
            for j in 0..n {
                let a_ij = a.get(i, j);
                let a_ji_conj = a.get(j, i).conj();
                let d = a_ij.sub(&a_ji_conj);
                assert!(d.norm() < 1e-12, "A not Hermitian at ({}, {})", i, j);
            }
        }
    }

    #[test]
    fn test_dos_positive() {
        let gf = make_simple_gf(5);
        // DOS should be positive inside the band
        let dos = gf.density_of_states(0.0).expect("ok");
        assert!(dos > 0.0, "DOS = {}", dos);
    }

    // --------------------------------------------------------- transmission

    #[test]
    fn test_transmission_in_range_0_to_1() {
        let tc = make_tc(5);
        let t = tc.transmission(0.0).expect("ok");
        assert!((0.0..=1.0 + 1e-8).contains(&t), "T = {}", t);
    }

    #[test]
    fn test_transmission_curve_length() {
        let tc = make_tc(5);
        let curve = tc.transmission_curve().expect("ok");
        assert_eq!(curve.len(), 100);
    }

    // --------------------------------------------------------- transport

    #[test]
    fn test_current_direction_follows_bias() {
        let tc = make_tc(5);
        let i_pos = tc.current(0.5, 300.0).expect("ok");
        let i_neg = tc.current(-0.5, 300.0).expect("ok");
        assert!(
            i_pos > 0.0,
            "positive bias should give positive current: I = {}",
            i_pos
        );
        assert!(
            i_neg < 0.0,
            "negative bias should give negative current: I = {}",
            i_neg
        );
    }

    // --------------------------------------------------------- Fermi–Dirac

    #[test]
    fn test_fermi_dirac_at_mu_is_half() {
        let f = TransportCalculator::fermi_dirac(0.5, 0.5, 300.0);
        assert!((f - 0.5).abs() < 1e-12, "f(μ, μ, T) = {}", f);
    }

    // --------------------------------------------------------- zero-bias conductance

    #[test]
    fn test_zero_bias_conductance_positive() {
        let tc = make_tc(3);
        let g = tc.zero_bias_conductance(300.0).expect("ok");
        assert!(g >= 0.0, "G₀ = {}", g);
    }
}
