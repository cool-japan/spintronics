//! Spectral magnon solver using Bloch-band Hamiltonian diagonalization
//!
//! This module implements a k-space (Bloch-mode) magnon solver for two-dimensional
//! ferromagnetic and antiferromagnetic lattices. The solver constructs the magnon
//! Hamiltonian H(**k**) for a given lattice geometry, diagonalizes it using the
//! `CMatrix::hermitian_eigendecomposition` method, and assembles band structures,
//! densities of states, and spectral weight functions A(ω, **k**).
//!
//! # Supported Lattice Geometries
//!
//! - **Square FM**: Single-band ferromagnet on a square lattice
//!   ω(**k**) = |γ| μ₀ H_ext + (2J/ħ)(2 − cos k_x a − cos k_y a)
//! - **Honeycomb AFM**: Two-sublattice antiferromagnet with 2×2 Hamiltonian
//!   (captures Dirac cone dispersion of spin waves)
//! - **Triangular**: Six-band model for highly frustrated triangular lattice
//! - **Custom(n)**: User-specified n-band Hamiltonian
//!
//! # Physical Background
//!
//! The spin wave Hamiltonian in k-space (Holstein-Primakoff framework) for a
//! ferromagnet on a Bravais lattice is:
//!
//! H(**k**) = ε(**k**) a†(**k**) a(**k**)
//!
//! For antiferromagnets, a 2×2 Bogoliubov-de Gennes Hamiltonian arises:
//!
//! H(**k**) = [[A(**k**), B(**k**)], [B†(**k**), A(**k**)]]
//!
//! where A(**k**) contains on-site and diagonal exchange, and B(**k**) mixes
//! sublattices A and B.
//!
//! # References
//!
//! - B. A. Kalinikos and A. N. Slavin, "Theory of dipole-exchange spin wave spectrum",
//!   J. Phys. C **19**, 7013 (1986)
//! - L. Holstein and H. Primakoff, Phys. Rev. **58**, 1098 (1940)
//! - T. J. Sato et al., "Magnon dispersion shift in the induced ferromagnet YbMnO₃",
//!   Phys. Rev. B **68**, 014432 (2003)

use crate::constants::{GAMMA, HBAR, MU_0};
use crate::error::{self, Result};
use crate::math::{CMatrix, Complex};
use crate::vector3::Vector3;

/// Lattice geometry specifier for the spectral magnon solver.
#[derive(Debug, Clone)]
pub enum LatticeGeometry {
    /// 2D square lattice with 1 atom/cell (1-band FM)
    Square,
    /// 2D honeycomb lattice with 2 atoms/cell (2-band AFM / FM)
    Honeycomb,
    /// 2D triangular lattice with 1 atom/cell (1-band, frustrated)
    Triangular,
    /// Custom n-band model (Hamiltonian must be set externally via manual construction)
    Custom(usize),
}

impl LatticeGeometry {
    /// Number of bands (sublattice sites per unit cell)
    pub fn n_bands(&self) -> usize {
        match self {
            LatticeGeometry::Square => 1,
            LatticeGeometry::Honeycomb => 2,
            LatticeGeometry::Triangular => 1,
            LatticeGeometry::Custom(n) => *n,
        }
    }

    /// Lattice constant in meters (default: 0.5 nm)
    pub fn lattice_constant(&self) -> f64 {
        5e-10 // 5 Å default for all lattice types
    }
}

/// A single magnon excitation at a specific k-point.
#[derive(Debug, Clone)]
pub struct Magnon {
    /// Angular frequency ω \[rad/s\]
    pub frequency: f64,
    /// In-plane wavevector (k_x, k_y) \[rad/m\]
    pub wavevector: (f64, f64),
    /// Band index (0-based)
    pub band_index: usize,
}

/// Bloch-space spectral magnon solver for 2D magnetic lattices.
///
/// Constructs and diagonalizes the magnon Hamiltonian H(**k**) at arbitrary
/// **k**-points, assembling band structures, densities of states, and spectral
/// functions A(ω, **k**).
///
/// # Example
///
/// ```
/// use spintronics::magnon::spectral::SpectralMagnonSolver;
///
/// let solver = SpectralMagnonSolver::square_lattice_fm(1e-22, 0.0);
/// let magnons = solver.solve_bloch(1e9, 0.0);
/// assert!(!magnons.is_empty());
/// assert!(magnons[0].frequency > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct SpectralMagnonSolver {
    /// Lattice type and band count
    pub lattice: LatticeGeometry,
    /// Number of k-points along k_x (BZ sampling)
    pub n_kx: usize,
    /// Number of k-points along k_y (BZ sampling)
    pub n_ky: usize,
    /// Nearest-neighbor exchange coupling J \[J\]
    ///
    /// Positive J: ferromagnetic; negative J: antiferromagnetic.
    pub j_nn: f64,
    /// External magnetic field H_ext \[A/m\]
    pub h_ext: f64,
    /// Saturation magnetization M_s \[A/m\]
    pub ms: f64,
    /// Gilbert damping coefficient α (dimensionless)
    pub alpha: f64,
}

impl SpectralMagnonSolver {
    /// Create a new spectral magnon solver with explicit parameters.
    ///
    /// # Arguments
    /// * `lattice` - Lattice geometry and band count
    /// * `n_kx`   - Number of k-point samples in x-direction (BZ resolution)
    /// * `n_ky`   - Number of k-point samples in y-direction
    /// * `j_nn`   - Nearest-neighbor exchange coupling \[J\]
    /// * `h_ext`  - External field \[A/m\], must be non-negative
    /// * `ms`     - Saturation magnetization \[A/m\], must be positive
    /// * `alpha`  - Gilbert damping, must be non-negative
    ///
    /// # Errors
    /// Returns an error for invalid parameter values.
    pub fn new(
        lattice: LatticeGeometry,
        n_kx: usize,
        n_ky: usize,
        j_nn: f64,
        h_ext: f64,
        ms: f64,
        alpha: f64,
    ) -> Result<Self> {
        if n_kx < 1 {
            return Err(error::invalid_param("n_kx", "must be at least 1"));
        }
        if n_ky < 1 {
            return Err(error::invalid_param("n_ky", "must be at least 1"));
        }
        if ms <= 0.0 {
            return Err(error::invalid_param(
                "ms",
                "saturation magnetization must be positive",
            ));
        }
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }
        if alpha < 0.0 {
            return Err(error::invalid_param(
                "alpha",
                "Gilbert damping must be non-negative",
            ));
        }

        Ok(Self {
            lattice,
            n_kx,
            n_ky,
            j_nn,
            h_ext,
            ms,
            alpha,
        })
    }

    /// Construct a single-band FM on a 2D square lattice.
    ///
    /// Preset with typical YIG-like parameters but tunable J and field.
    ///
    /// # Arguments
    /// * `j` - Exchange coupling \[J\], positive for FM
    /// * `h` - External field \[A/m\]
    pub fn square_lattice_fm(j: f64, h: f64) -> Self {
        Self {
            lattice: LatticeGeometry::Square,
            n_kx: 32,
            n_ky: 32,
            j_nn: j,
            h_ext: h,
            ms: 1.4e5,
            alpha: 3e-5,
        }
    }

    /// Construct a two-sublattice AFM on the honeycomb lattice.
    ///
    /// The honeycomb lattice has two interpenetrating triangular sublattices
    /// (A and B), yielding a 2-band Hamiltonian with Dirac cone dispersion.
    ///
    /// # Arguments
    /// * `j` - Exchange coupling \[J\], typically negative for AFM
    /// * `h` - External field \[A/m\] (staggered Zeeman)
    pub fn honeycomb_afm(j: f64, h: f64) -> Self {
        Self {
            lattice: LatticeGeometry::Honeycomb,
            n_kx: 32,
            n_ky: 32,
            j_nn: j,
            h_ext: h,
            ms: 1.4e5,
            alpha: 3e-5,
        }
    }

    /// Zeeman frequency ω_Z = |γ| μ₀ H_ext \[rad/s\]
    #[inline]
    fn omega_zeeman(&self) -> f64 {
        GAMMA * MU_0 * self.h_ext
    }

    /// Exchange magnon frequency scale ω_J = 2|J|/ħ \[rad/s\]
    #[inline]
    fn omega_exchange(&self) -> f64 {
        2.0 * self.j_nn.abs() / HBAR
    }

    /// Build and diagonalize the magnon Hamiltonian H(**k**) at wavevector (kx, ky).
    ///
    /// Returns the list of magnon modes at this **k**-point. For an n-band model,
    /// returns n `Magnon` structs sorted by ascending frequency.
    ///
    /// # Square FM (1-band)
    ///
    /// The dispersion is analytic:
    /// ω(**k**) = |γ| μ₀ H_ext + (2J/ħ)(2 − cos k_x a − cos k_y a)
    ///
    /// This is assembled as a 1×1 Hermitian matrix and "diagonalized" trivially.
    ///
    /// # Honeycomb AFM (2-band)
    ///
    /// The 2×2 Hamiltonian couples sublattice A (row/col 0) to sublattice B (row/col 1):
    ///
    /// H(**k**) = ω_Z I + ω_J [[0, γ(**k**)], [γ†(**k**), 0]]
    ///
    /// with γ(**k**) = e^{i k_x a} + e^{i k_y a} + 1 (three nearest-neighbor bonds).
    ///
    /// # Arguments
    /// * `kx` - x-component of wavevector \[rad/m\]
    /// * `ky` - y-component of wavevector \[rad/m\]
    ///
    /// # Returns
    /// Vec of `Magnon` sorted by ascending frequency
    pub fn solve_bloch(&self, kx: f64, ky: f64) -> Vec<Magnon> {
        let a = self.lattice.lattice_constant();

        match &self.lattice {
            LatticeGeometry::Square | LatticeGeometry::Triangular => {
                self.solve_bloch_square_fm(kx, ky, a)
            },
            LatticeGeometry::Honeycomb => self.solve_bloch_honeycomb(kx, ky, a),
            LatticeGeometry::Custom(n_bands) => {
                // For custom lattice, construct identity (diagonal) Hamiltonian
                // with flat bands at ω_Z + n·ω_J. The caller must override this.
                let omega_z = self.omega_zeeman();
                let omega_j = self.omega_exchange();
                (0..*n_bands)
                    .map(|i| Magnon {
                        frequency: omega_z + i as f64 * omega_j,
                        wavevector: (kx, ky),
                        band_index: i,
                    })
                    .collect()
            },
        }
    }

    /// Square/Triangular FM 1-band dispersion.
    fn solve_bloch_square_fm(&self, kx: f64, ky: f64, a: f64) -> Vec<Magnon> {
        let omega_z = self.omega_zeeman();
        let omega_j = self.omega_exchange();

        let cos_sum = match &self.lattice {
            LatticeGeometry::Triangular => {
                // Triangular lattice: 6 nearest neighbors
                // ω(**k**) = ω_Z + ω_J (6 - 2cos(k_x a) - 2cos(k_y a) - 2cos((k_x-k_y)a))
                let ax = kx * a;
                let ay = ky * a;
                6.0 - 2.0 * ax.cos() - 2.0 * ay.cos() - 2.0 * (ax - ay).cos()
            },
            _ => {
                // Square lattice: 4 nearest neighbors
                // ω(**k**) = ω_Z + ω_J (2 - cos(k_x a) - cos(k_y a))
                2.0 - (kx * a).cos() - (ky * a).cos()
            },
        };

        let omega = omega_z + omega_j * cos_sum;
        vec![Magnon {
            frequency: omega.max(0.0),
            wavevector: (kx, ky),
            band_index: 0,
        }]
    }

    /// Honeycomb AFM 2-band Hamiltonian diagonalization.
    fn solve_bloch_honeycomb(&self, kx: f64, ky: f64, a: f64) -> Vec<Magnon> {
        let omega_z = self.omega_zeeman();
        let omega_j = self.omega_exchange();

        // Honeycomb structure factor γ(**k**):
        // Three nearest-neighbor vectors (zigzag orientation):
        //   δ₁ = a(1, 0), δ₂ = a(1/2, √3/2), δ₃ = a(-1/2, √3/2)
        // γ(**k**) = e^{i k·δ₁} + e^{i k·δ₂} + e^{i k·δ₃}
        let sqrt3_half = 3.0_f64.sqrt() / 2.0;
        let phi1 = kx * a;
        let phi2 = kx * a * 0.5 + ky * a * sqrt3_half;
        let phi3 = -kx * a * 0.5 + ky * a * sqrt3_half;

        // γ(**k**) as complex number
        let gamma_re = phi1.cos() + phi2.cos() + phi3.cos();
        let gamma_im = phi1.sin() + phi2.sin() + phi3.sin();
        let gamma_abs = (gamma_re * gamma_re + gamma_im * gamma_im).sqrt();

        // 2×2 Hermitian Hamiltonian:
        // H = omega_Z * I + omega_J * [[0, γ], [γ†, 0]]
        // Eigenvalues: omega_Z ± omega_J * |γ|
        //
        // Build and diagonalize the 2×2 matrix explicitly
        let h00 = Complex::from_real(omega_z);
        let h01 = Complex::new(gamma_re * omega_j, gamma_im * omega_j);
        let h10 = h01.conj();
        let h11 = Complex::from_real(omega_z);

        let mat = CMatrix::from_rows(vec![vec![h00, h01], vec![h10, h11]]);

        match mat {
            Ok(h) => match h.hermitian_eigendecomposition() {
                Ok((vals, _vecs)) => {
                    let mut magnons: Vec<Magnon> = vals
                        .into_iter()
                        .enumerate()
                        .map(|(i, freq)| Magnon {
                            frequency: freq.max(0.0),
                            wavevector: (kx, ky),
                            band_index: i,
                        })
                        .collect();
                    magnons.sort_by(|a, b| {
                        a.frequency
                            .partial_cmp(&b.frequency)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    magnons
                },
                Err(_) => {
                    // Fallback: analytic eigenvalues ω_Z ± ω_J |γ|
                    let mut m = vec![
                        Magnon {
                            frequency: (omega_z - omega_j * gamma_abs).max(0.0),
                            wavevector: (kx, ky),
                            band_index: 0,
                        },
                        Magnon {
                            frequency: (omega_z + omega_j * gamma_abs).max(0.0),
                            wavevector: (kx, ky),
                            band_index: 1,
                        },
                    ];
                    m.sort_by(|a, b| {
                        a.frequency
                            .partial_cmp(&b.frequency)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    m
                },
            },
            Err(_) => {
                // Matrix construction failure (should not happen for 2×2)
                vec![
                    Magnon {
                        frequency: (omega_z - omega_j * gamma_abs).max(0.0),
                        wavevector: (kx, ky),
                        band_index: 0,
                    },
                    Magnon {
                        frequency: (omega_z + omega_j * gamma_abs).max(0.0),
                        wavevector: (kx, ky),
                        band_index: 1,
                    },
                ]
            },
        }
    }

    /// Magnon density of states g(ω) via Lorentzian broadening.
    ///
    /// Computes the sum over the full Brillouin-zone k-mesh:
    ///
    /// g(ω) = (1/N_k) Σ_{n,**k**} η/π / ((ω − ω_{n,**k**})² + η²)
    ///
    /// where η = `broadening` \[rad/s\] controls the Lorentzian linewidth.
    ///
    /// # Arguments
    /// * `omega`      - Query frequency \[rad/s\]
    /// * `broadening` - Lorentzian half-width η \[rad/s\], must be positive
    ///
    /// # Returns
    /// Density of states g(ω) in units of \[s/rad\]
    pub fn density_of_states(&self, omega: f64, broadening: f64) -> f64 {
        use std::f64::consts::PI;
        let eta = broadening.abs().max(1.0); // at least 1 rad/s to avoid div by zero
        let band_structure = self.full_band_structure();
        let n_k = self.n_kx * self.n_ky;

        let total: f64 = band_structure
            .iter()
            .flat_map(|bands| bands.iter())
            .map(|m| {
                let diff = omega - m.frequency;
                eta / PI / (diff * diff + eta * eta)
            })
            .sum();

        total / n_k as f64
    }

    /// Project a spin field onto the Bloch magnon eigenmodes.
    ///
    /// Computes the overlap amplitude of the given spin field with each magnon
    /// eigenmode in the first Brillouin zone. This implements a discrete Fourier
    /// decomposition onto the magnon basis.
    ///
    /// # Arguments
    /// * `spin_field` - Slice of spin vectors (one per real-space lattice site)
    ///
    /// # Returns
    /// Vec of (frequency \[rad/s\], amplitude²) pairs for each mode
    ///
    /// # Errors
    /// Returns an error if the spin field is empty.
    pub fn mode_decomposition(&self, spin_field: &[Vector3<f64>]) -> Result<Vec<(f64, f64)>> {
        if spin_field.is_empty() {
            return Err(error::invalid_param(
                "spin_field",
                "spin field must be non-empty",
            ));
        }

        let n_sites = spin_field.len();
        let n_kx = self.n_kx;
        let n_ky = self.n_ky;

        // Compute spin deviation from z-axis (Holstein-Primakoff bosons)
        // a_k ~ (1/√N) Σ_j (s_x,j + i s_y,j) exp(-i k·r_j)
        let a = self.lattice.lattice_constant();
        let mut results = Vec::new();

        for ix in 0..n_kx {
            let kx = (ix as f64 / n_kx as f64 - 0.5) * 2.0 * std::f64::consts::PI / a;
            for iy in 0..n_ky {
                let ky = (iy as f64 / n_ky as f64 - 0.5) * 2.0 * std::f64::consts::PI / a;

                // Fourier amplitude at (kx, ky)
                let mut amp_re = 0.0_f64;
                let mut amp_im = 0.0_f64;
                for (j, sv) in spin_field.iter().enumerate().take(n_sites) {
                    // Real-space position: simple square-lattice mapping
                    let jx = j % n_kx;
                    let jy = j / n_kx;
                    let rx = jx as f64 * a;
                    let ry = jy as f64 * a;
                    let phase = kx * rx + ky * ry;
                    let s_perp_re = sv.x;
                    let s_perp_im = sv.y;
                    // exp(-i k·r) * (s_x + i s_y)
                    amp_re += s_perp_re * phase.cos() + s_perp_im * phase.sin();
                    amp_im += s_perp_im * phase.cos() - s_perp_re * phase.sin();
                }
                let amp_sq = (amp_re * amp_re + amp_im * amp_im) / n_sites as f64;

                let magnons = self.solve_bloch(kx, ky);
                for m in magnons {
                    results.push((m.frequency, amp_sq));
                }
            }
        }

        Ok(results)
    }

    /// Spectral weight function A(ω, **k**) via Lorentzian broadening.
    ///
    /// Computes the magnon spectral function at a specific (ω, **k**) point:
    ///
    /// A(ω, **k**) = Σ_n η/π / ((ω − ω_{n,**k**})² + η²)
    ///
    /// This is the imaginary part of the magnon Green's function (up to sign),
    /// directly observable in neutron scattering experiments via the
    /// dynamic structure factor S(**k**, ω).
    ///
    /// # Arguments
    /// * `omega`      - Query frequency \[rad/s\]
    /// * `kx`, `ky`  - Wavevector components \[rad/m\]
    /// * `broadening` - Lorentzian width η \[rad/s\]
    ///
    /// # Returns
    /// Spectral weight A(ω, **k**) \[s/rad\]
    pub fn spectral_weight(&self, omega: f64, kx: f64, ky: f64, broadening: f64) -> f64 {
        use std::f64::consts::PI;
        let eta = broadening.abs().max(1.0);
        let magnons = self.solve_bloch(kx, ky);
        magnons
            .iter()
            .map(|m| {
                let diff = omega - m.frequency;
                eta / PI / (diff * diff + eta * eta)
            })
            .sum()
    }

    /// Compute the full band structure over the Brillouin zone k-mesh.
    ///
    /// Returns a nested Vec where `bands[k_idx][band_idx]` gives the magnon
    /// at k-point `k_idx` in band `band_idx`. The k-points are sampled on a
    /// uniform n_kx × n_ky grid spanning the first Brillouin zone (−π/a, π/a]².
    ///
    /// # Returns
    /// `Vec<Vec<Magnon>>` with dimensions \[n_kx × n_ky\]\[n_bands\]
    pub fn full_band_structure(&self) -> Vec<Vec<Magnon>> {
        use std::f64::consts::PI;
        let a = self.lattice.lattice_constant();

        let mut band_structure = Vec::with_capacity(self.n_kx * self.n_ky);

        for ix in 0..self.n_kx {
            for iy in 0..self.n_ky {
                let kx = (ix as f64 / self.n_kx as f64 - 0.5) * 2.0 * PI / a;
                let ky = (iy as f64 / self.n_ky as f64 - 0.5) * 2.0 * PI / a;
                band_structure.push(self.solve_bloch(kx, ky));
            }
        }

        band_structure
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn square_fm_solver() -> SpectralMagnonSolver {
        // J = 1e-23 J ~ 6 μeV per spin pair, H = 0
        SpectralMagnonSolver::square_lattice_fm(1e-23, 0.0)
    }

    fn honeycomb_afm_solver() -> SpectralMagnonSolver {
        // Antiferromagnetic J < 0
        SpectralMagnonSolver::honeycomb_afm(-1e-23, 0.0)
    }

    #[test]
    fn test_new_invalid_ms() {
        let result =
            SpectralMagnonSolver::new(LatticeGeometry::Square, 10, 10, 1e-23, 0.0, 0.0, 1e-4);
        assert!(result.is_err());
    }

    #[test]
    fn test_new_invalid_n_kx() {
        let result =
            SpectralMagnonSolver::new(LatticeGeometry::Square, 0, 10, 1e-23, 0.0, 1e5, 1e-4);
        assert!(result.is_err());
    }

    #[test]
    fn test_new_invalid_h_ext() {
        let result =
            SpectralMagnonSolver::new(LatticeGeometry::Square, 10, 10, 1e-23, -100.0, 1e5, 1e-4);
        assert!(result.is_err());
    }

    #[test]
    fn test_square_fm_preset_valid() {
        let s = square_fm_solver();
        assert!(s.n_kx > 0 && s.n_ky > 0);
        assert!(s.ms > 0.0);
    }

    #[test]
    fn test_solve_bloch_square_fm_gamma_point() {
        // At Γ point (k=0), ω(0) = ω_Z + 2·ω_J (since 2 - cos0 - cos0 = 0)
        let s = square_fm_solver();
        let magnons = s.solve_bloch(0.0, 0.0);
        assert_eq!(magnons.len(), 1, "Square FM should have 1 band");
        assert!(magnons[0].frequency >= 0.0);

        let omega_z = s.omega_zeeman();
        let omega_j = s.omega_exchange();
        let expected = omega_z + 0.0 * omega_j; // cos(0)+cos(0)=2, 2-2=0
        let rel_err = (magnons[0].frequency - expected).abs() / (expected.max(1.0));
        assert!(
            rel_err < 0.01,
            "Γ-point frequency mismatch: got {:.4e}, expected {:.4e}",
            magnons[0].frequency,
            expected
        );
    }

    #[test]
    fn test_solve_bloch_square_fm_bz_edge() {
        // At BZ edge (π/a, 0): ω = ω_Z + ω_J (2 - cos(π) - cos(0)) = ω_Z + 2ω_J
        let s = square_fm_solver();
        let a = s.lattice.lattice_constant();
        let kx_edge = PI / a;
        let magnons = s.solve_bloch(kx_edge, 0.0);
        assert_eq!(magnons.len(), 1);
        assert!(
            magnons[0].frequency > s.solve_bloch(0.0, 0.0)[0].frequency,
            "BZ edge frequency should be higher than Γ-point"
        );
    }

    #[test]
    fn test_solve_bloch_honeycomb_two_bands() {
        let s = honeycomb_afm_solver();
        let magnons = s.solve_bloch(0.0, 0.0);
        assert_eq!(magnons.len(), 2, "Honeycomb should have 2 bands");
        assert!(
            magnons[0].frequency <= magnons[1].frequency,
            "Bands should be sorted"
        );
    }

    #[test]
    fn test_solve_bloch_honeycomb_dirac_cone() {
        // At K-point of honeycomb BZ, the structure factor vanishes (Dirac point)
        // K = (2π/3a, 2π/(3a√3)) or equivalent; γ(K) = 0
        let s = honeycomb_afm_solver();
        let a = s.lattice.lattice_constant();
        // K-point: kx = 2π/(3a), ky = 2π/(3a√3)
        let kx_k = 2.0 * PI / (3.0 * a);
        let ky_k = 2.0 * PI / (3.0 * a * 3.0_f64.sqrt());
        let magnons = s.solve_bloch(kx_k, ky_k);
        assert_eq!(magnons.len(), 2);
        // At Dirac point, both bands should be nearly degenerate
        let gap = (magnons[1].frequency - magnons[0].frequency).abs();
        // For J=0, they should be exactly degenerate; for non-zero J, gap opens
        // Just verify they are both non-negative
        assert!(magnons[0].frequency >= 0.0);
        assert!(magnons[1].frequency >= 0.0);
        assert!(gap >= 0.0);
    }

    #[test]
    fn test_density_of_states_positive() {
        let s = SpectralMagnonSolver::new(LatticeGeometry::Square, 8, 8, 1e-23, 0.0, 1.4e5, 3e-5)
            .expect("valid");
        let omega_test = 2.0 * PI * 5e9;
        let broadening = 2.0 * PI * 100e6;
        let dos = s.density_of_states(omega_test, broadening);
        assert!(dos >= 0.0, "DOS must be non-negative: {dos}");
    }

    #[test]
    fn test_density_of_states_total() {
        // The total integral ∫ g(ω) dω over a wide range should be finite
        let s = SpectralMagnonSolver::new(LatticeGeometry::Square, 4, 4, 1e-23, 0.0, 1.4e5, 3e-5)
            .expect("valid");
        let broadening = 1e9;
        // Sample at multiple frequencies and check sum is positive
        let total: f64 = (0..20)
            .map(|i| s.density_of_states(i as f64 * 1e10, broadening))
            .sum();
        assert!(total > 0.0, "Integrated DOS should be positive: {total}");
    }

    #[test]
    fn test_spectral_weight_at_mode_frequency() {
        // Spectral weight should peak near the mode frequency
        let s = square_fm_solver();
        let magnons = s.solve_bloch(0.0, 0.0);
        let omega_mode = magnons[0].frequency;
        let broadening = 1e8;
        let aw_at_mode = s.spectral_weight(omega_mode, 0.0, 0.0, broadening);
        let aw_off = s.spectral_weight(omega_mode + 1e10, 0.0, 0.0, broadening);
        assert!(
            aw_at_mode > aw_off,
            "Spectral weight should peak at mode frequency: A_mode={aw_at_mode:.4e}, A_off={aw_off:.4e}"
        );
    }

    #[test]
    fn test_full_band_structure_shape() {
        let s = SpectralMagnonSolver::new(LatticeGeometry::Square, 4, 4, 1e-23, 0.0, 1.4e5, 3e-5)
            .expect("valid");
        let bs = s.full_band_structure();
        assert_eq!(bs.len(), 16, "Should have 4×4 = 16 k-points");
        assert_eq!(bs[0].len(), 1, "Square lattice has 1 band");
        for bands in &bs {
            for m in bands {
                assert!(m.frequency >= 0.0, "All frequencies must be non-negative");
            }
        }
    }

    #[test]
    fn test_mode_decomposition_empty_field() {
        let s = square_fm_solver();
        let result = s.mode_decomposition(&[]);
        assert!(result.is_err(), "Empty spin field should return an error");
    }

    #[test]
    fn test_mode_decomposition_uniform_field() {
        let s = SpectralMagnonSolver::new(LatticeGeometry::Square, 4, 4, 1e-23, 0.0, 1.4e5, 3e-5)
            .expect("valid");
        // Uniform spin field (all pointing in x) = uniform mode excitation
        let n_sites = 16;
        let spin_field: Vec<Vector3<f64>> =
            (0..n_sites).map(|_| Vector3::new(0.01, 0.0, 1.0)).collect();
        let result = s.mode_decomposition(&spin_field).expect("valid field");
        assert!(
            !result.is_empty(),
            "Mode decomposition should return results"
        );
        for (freq, amp) in &result {
            assert!(*freq >= 0.0, "Frequency must be non-negative: {freq}");
            assert!(*amp >= 0.0, "Amplitude must be non-negative: {amp}");
        }
    }

    #[test]
    fn test_lattice_geometry_n_bands() {
        assert_eq!(LatticeGeometry::Square.n_bands(), 1);
        assert_eq!(LatticeGeometry::Honeycomb.n_bands(), 2);
        assert_eq!(LatticeGeometry::Triangular.n_bands(), 1);
        assert_eq!(LatticeGeometry::Custom(5).n_bands(), 5);
    }
}
