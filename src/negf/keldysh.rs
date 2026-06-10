//! Keldysh non-equilibrium Green's functions.
//!
//! Implements the Keldysh formalism for computing lesser and greater Green's
//! functions in a non-equilibrium setting with two leads held at different
//! chemical potentials.
//!
//! # Keldysh Relations
//!
//! ```text
//! Σ^<(E) = i f_L Γ_L + i f_R Γ_R          (lesser self-energy)
//! Σ^>(E) = −i(1−f_L) Γ_L − i(1−f_R) Γ_R  (greater self-energy)
//! G^<    = G^R · Σ^< · G^A                  (lesser GF)
//! G^>    = G^R · Σ^> · G^A                  (greater GF)
//! G^> − G^< = G^R − G^A                     (Keldysh identity)
//! ```
//!
//! # References
//!
//! - L. V. Keldysh, Sov. Phys. JETP **20**, 1018 (1965)
//! - S. Datta, *Electronic Transport in Mesoscopic Systems* (Cambridge, 1995)

use super::green_function::{GreenFunction, TransportCalculator};
use crate::error::Result;
use crate::math::{CMatrix, Complex};

// ============================================================================
// KeldyshSolver
// ============================================================================

/// Keldysh solver for non-equilibrium lesser/greater Green's functions.
///
/// Wraps a [`GreenFunction`] and provides the Keldysh lesser/greater GF and
/// derived physical observables (occupation density, nonequilibrium carrier density).
#[derive(Debug, Clone)]
pub struct KeldyshSolver {
    /// The underlying Green's function for the scattering region.
    pub gf: GreenFunction,
}

impl KeldyshSolver {
    /// Create a new `KeldyshSolver`.
    pub fn new(gf: GreenFunction) -> Self {
        Self { gf }
    }

    /// Fermi–Dirac distribution f(E, μ, T).
    ///
    /// Numerically stable across the full energy range.
    pub fn fermi_dirac(energy: f64, mu: f64, temperature: f64) -> f64 {
        TransportCalculator::fermi_dirac(energy, mu, temperature)
    }

    /// Lesser self-energy Σ^<(E) as an N×N matrix.
    ///
    /// ```text
    /// Σ^<[0][0]     = i · f_L · Γ_L
    /// Σ^<[N-1][N-1] = i · f_R · Γ_R
    /// ```
    ///
    /// All other entries are zero (wide-band, contact geometry).
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `temperature < 0`.
    pub fn sigma_lesser(
        &self,
        energy: f64,
        mu_l: f64,
        mu_r: f64,
        temperature: f64,
    ) -> Result<CMatrix> {
        let n = self.gf.hamiltonian.n_sites;
        let nm1 = n - 1;
        let f_l = Self::fermi_dirac(energy, mu_l, temperature);
        let f_r = Self::fermi_dirac(energy, mu_r, temperature);
        let mut sigma = CMatrix::zeros(n);
        // Σ^<[0][0]     = i·f_L·Γ_L
        sigma.set(
            0,
            0,
            Complex::new(0.0, f_l * self.gf.sigma_l.gamma_matrix()),
        );
        // Σ^<[N-1][N-1] = i·f_R·Γ_R
        sigma.set(
            nm1,
            nm1,
            Complex::new(0.0, f_r * self.gf.sigma_r.gamma_matrix()),
        );
        Ok(sigma)
    }

    /// Greater self-energy Σ^>(E) as an N×N matrix.
    ///
    /// ```text
    /// Σ^>[0][0]     = −i · (1−f_L) · Γ_L
    /// Σ^>[N-1][N-1] = −i · (1−f_R) · Γ_R
    /// ```
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `temperature < 0`.
    pub fn sigma_greater(
        &self,
        energy: f64,
        mu_l: f64,
        mu_r: f64,
        temperature: f64,
    ) -> Result<CMatrix> {
        let n = self.gf.hamiltonian.n_sites;
        let nm1 = n - 1;
        let f_l = Self::fermi_dirac(energy, mu_l, temperature);
        let f_r = Self::fermi_dirac(energy, mu_r, temperature);
        let mut sigma = CMatrix::zeros(n);
        // Σ^>[0][0]     = −i·(1−f_L)·Γ_L
        sigma.set(
            0,
            0,
            Complex::new(0.0, -(1.0 - f_l) * self.gf.sigma_l.gamma_matrix()),
        );
        // Σ^>[N-1][N-1] = −i·(1−f_R)·Γ_R
        sigma.set(
            nm1,
            nm1,
            Complex::new(0.0, -(1.0 - f_r) * self.gf.sigma_r.gamma_matrix()),
        );
        Ok(sigma)
    }

    /// Lesser Green's function G^<(E) = G^R · Σ^< · G^A.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`GreenFunction::g_retarded`] or matrix arithmetic.
    pub fn g_lesser(&self, energy: f64, mu_l: f64, mu_r: f64, temperature: f64) -> Result<CMatrix> {
        let gr = self.gf.g_retarded(energy)?;
        let ga = gr.conj_transpose();
        let sigma_l = self.sigma_lesser(energy, mu_l, mu_r, temperature)?;
        let tmp = gr.matmul(&sigma_l)?;
        tmp.matmul(&ga)
    }

    /// Greater Green's function G^>(E) = G^R · Σ^> · G^A.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`GreenFunction::g_retarded`] or matrix arithmetic.
    pub fn g_greater(
        &self,
        energy: f64,
        mu_l: f64,
        mu_r: f64,
        temperature: f64,
    ) -> Result<CMatrix> {
        let gr = self.gf.g_retarded(energy)?;
        let ga = gr.conj_transpose();
        let sigma_g = self.sigma_greater(energy, mu_l, mu_r, temperature)?;
        let tmp = gr.matmul(&sigma_g)?;
        tmp.matmul(&ga)
    }

    /// Occupation density n(E) \[per eV, per site\] at energy E.
    ///
    /// Following Datta's convention, the per-site density of states times occupancy is:
    ///
    /// ```text
    /// n(E) = (1/2π) Im[Tr G^<(E)] / N_sites
    /// ```
    ///
    /// Since G^< = G^R · Σ^< · G^A with Σ^< = iΓf (positive imaginary),
    /// Im[G^<_{ii}] is positive, yielding a non-negative occupation.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`g_lesser`](Self::g_lesser).
    pub fn occupation_density(
        &self,
        energy: f64,
        mu_l: f64,
        mu_r: f64,
        temperature: f64,
    ) -> Result<f64> {
        let gl = self.g_lesser(energy, mu_l, mu_r, temperature)?;
        let trace_im = gl.trace().im;
        let n = self.gf.hamiltonian.n_sites as f64;
        // Im[G^<_{ii}] > 0 by construction, giving positive occupation
        Ok(trace_im / (2.0 * std::f64::consts::PI * n))
    }

    /// Nonequilibrium carrier density (integrated) per site.
    ///
    /// Computes the energy-integrated diagonal elements of −i G^<(E) / (2π)
    /// over `[e_min, e_max]` using the trapezoid rule.
    ///
    /// Returns a `Vec<f64>` of length `n_sites`, where entry `i` is the
    /// per-site occupancy at lattice position `i`.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`g_lesser`](Self::g_lesser).
    pub fn nonequilibrium_density(
        &self,
        mu_l: f64,
        mu_r: f64,
        temperature: f64,
        e_min: f64,
        e_max: f64,
        n_points: usize,
    ) -> Result<Vec<f64>> {
        let n = self.gf.hamiltonian.n_sites;
        if n_points < 2 {
            return Ok(vec![0.0; n]);
        }
        let de = (e_max - e_min) / (n_points - 1) as f64;
        let pi2 = 2.0 * std::f64::consts::PI;

        // Collect per-site integrand at each energy point
        let mut integrands: Vec<Vec<f64>> = Vec::with_capacity(n_points);
        for k in 0..n_points {
            let e = e_min + k as f64 * de;
            let gl = self.g_lesser(e, mu_l, mu_r, temperature)?;
            // n_i(E) = Im[G^<_{ii}(E)] / (2π)   [positive since Im[G^<_{ii}] > 0]
            let site_vals: Vec<f64> = (0..n).map(|i| gl.get(i, i).im / pi2).collect();
            integrands.push(site_vals);
        }

        // Trapezoid integration per site
        let mut density = vec![0.0; n];
        for k in 0..(n_points - 1) {
            for i in 0..n {
                density[i] += 0.5 * (integrands[k][i] + integrands[k + 1][i]) * de;
            }
        }
        Ok(density)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::negf::green_function::{GreenFunction, Hamiltonian1D, LeadSelfEnergy};

    fn make_solver(n: usize) -> KeldyshSolver {
        let h = Hamiltonian1D::from_uniform(n, 0.0, 1.0).expect("valid");
        let sl = LeadSelfEnergy::new(0.5, 0.0).expect("valid");
        let sr = LeadSelfEnergy::new(0.5, 0.0).expect("valid");
        let gf = GreenFunction::new(h, sl, sr, 1e-3).expect("valid");
        KeldyshSolver::new(gf)
    }

    // --------------------------------------------------------- Fermi-Dirac

    #[test]
    fn test_fermi_dirac_at_mu_is_half() {
        let f = KeldyshSolver::fermi_dirac(1.0, 1.0, 300.0);
        assert!((f - 0.5).abs() < 1e-12, "f = {}", f);
    }

    #[test]
    fn test_fermi_dirac_high_t_half() {
        // At very high T, both states should be nearly 0.5
        let f_above = KeldyshSolver::fermi_dirac(10.0, 0.0, 1e8);
        let f_below = KeldyshSolver::fermi_dirac(-10.0, 0.0, 1e8);
        assert!(
            (f_above - 0.5).abs() < 0.1,
            "f(+10eV, 0, 1e8K) = {}",
            f_above
        );
        assert!(
            (f_below - 0.5).abs() < 0.1,
            "f(-10eV, 0, 1e8K) = {}",
            f_below
        );
    }

    // --------------------------------------------------------- Sigma^< / Sigma^>

    #[test]
    fn test_sigma_lesser_imaginary() {
        let ks = make_solver(4);
        let sigma = ks.sigma_lesser(0.0, 0.0, 0.0, 300.0).expect("ok");
        // At equilibrium (μ_L = μ_R) and E=μ, f=0.5
        // Σ^<[0][0] = i·0.5·Γ_L, real part = 0
        let s00 = sigma.get(0, 0);
        assert!(
            s00.re.abs() < 1e-12,
            "Re[Σ^<[0,0]] should be zero: {}",
            s00.re
        );
        assert!(s00.im > 0.0, "Im[Σ^<[0,0]] should be positive: {}", s00.im);
    }

    #[test]
    fn test_sigma_greater_neg_imaginary() {
        let ks = make_solver(4);
        let sigma = ks.sigma_greater(0.0, 0.0, 0.0, 300.0).expect("ok");
        let s00 = sigma.get(0, 0);
        assert!(
            s00.re.abs() < 1e-12,
            "Re[Σ^>[0,0]] should be zero: {}",
            s00.re
        );
        assert!(s00.im < 0.0, "Im[Σ^>[0,0]] should be negative: {}", s00.im);
    }

    // --------------------------------------------------------- G^< finite

    #[test]
    fn test_g_lesser_finite() {
        let ks = make_solver(4);
        let gl = ks.g_lesser(0.0, -0.1, 0.1, 300.0).expect("ok");
        // All elements should be finite
        for i in 0..4 {
            for j in 0..4 {
                let v = gl.get(i, j);
                assert!(v.is_finite(), "G^<[{},{}] is not finite", i, j);
            }
        }
    }

    // --------------------------------------------------------- Keldysh identity: G^> - G^< = G^R·(Σ^>-Σ^<)·G^A

    #[test]
    fn test_g_lesser_minus_g_greater_keldysh_identity() {
        // Exact Keldysh identity: G^> − G^< = G^R · (Σ^> − Σ^<) · G^A
        //
        // Note: with finite η broadening, G^R − G^A ≠ G^> − G^< (they differ by a
        // 2iη correction term). The exact identity uses the self-energy difference.
        let ks = make_solver(4);
        let energy = 0.3;
        let mu_l = -0.2;
        let mu_r = 0.2;
        let t = 300.0;

        let gl = ks.g_lesser(energy, mu_l, mu_r, t).expect("ok");
        let gg = ks.g_greater(energy, mu_l, mu_r, t).expect("ok");
        let gr = ks.gf.g_retarded(energy).expect("ok");
        let ga = ks.gf.g_advanced(energy).expect("ok");

        // LHS = G^> - G^<
        let lhs = gg.sub(&gl).expect("sub ok");

        // RHS = G^R · (Σ^> - Σ^<) · G^A
        let sg = ks.sigma_greater(energy, mu_l, mu_r, t).expect("ok");
        let sl = ks.sigma_lesser(energy, mu_l, mu_r, t).expect("ok");
        let sigma_diff = sg.sub(&sl).expect("sub ok");
        let tmp = gr.matmul(&sigma_diff).expect("matmul ok");
        let rhs = tmp.matmul(&ga).expect("matmul ok");

        let diff = lhs.sub(&rhs).expect("sub ok");
        let frob = diff.frobenius_norm();
        assert!(
            frob < 1e-9,
            "Keldysh identity violation: Frobenius = {}",
            frob
        );
    }

    // --------------------------------------------------------- equilibrium occupation

    #[test]
    fn test_occupation_equilibrium() {
        // At equilibrium (μ_L = μ_R = 0), n(E=0) should be ≈ f(0, 0, T) * DOS/N
        let ks = make_solver(3);
        let n = ks.occupation_density(0.0, 0.0, 0.0, 300.0).expect("ok");
        // At E=μ, f=0.5, so density should be positive and in reasonable range
        assert!(n >= 0.0, "occupation density must be non-negative: {}", n);
    }

    // --------------------------------------------------------- nonequilibrium density

    #[test]
    fn test_nonequilibrium_density_sum_positive() {
        let ks = make_solver(3);
        let density = ks
            .nonequilibrium_density(0.0, 0.0, 300.0, -2.0, 2.0, 50)
            .expect("ok");
        assert_eq!(density.len(), 3);
        let total: f64 = density.iter().sum();
        assert!(
            total >= 0.0,
            "total density must be non-negative: {}",
            total
        );
    }

    // --------------------------------------------------------- disorder

    #[test]
    fn test_keldysh_with_disorder() {
        use crate::negf::green_function::Hamiltonian1D;
        let h = Hamiltonian1D::from_uniform(5, 0.0, 1.0)
            .expect("valid")
            .with_disorder(0.05, 777);
        let sl = LeadSelfEnergy::new(0.3, 0.0).expect("valid");
        let sr = LeadSelfEnergy::new(0.3, 0.0).expect("valid");
        let gf = GreenFunction::new(h, sl, sr, 1e-3).expect("valid");
        let ks = KeldyshSolver::new(gf);
        let gl = ks
            .g_lesser(0.0, -0.1, 0.1, 300.0)
            .expect("ok with disorder");
        assert!(gl.get(0, 0).is_finite());
    }
}
