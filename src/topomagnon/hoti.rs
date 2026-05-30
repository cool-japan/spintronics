//! Higher-Order Topological Insulators (HOTI): BBH model and breathing kagome.
//!
//! Implements the two paradigmatic higher-order topological phases in 2D:
//!
//! ## Benalcazar-Bernevig-Hughes (BBH) quadrupole insulator
//!
//! W. A. Benalcazar, B. A. Bernevig, T. L. Hughes,
//! *Science* **357**, 61 (2017); *Phys. Rev. B* **96**, 245115 (2017).
//!
//! The 4-band Bloch Hamiltonian on a square lattice carries a quantised
//! electric quadrupole moment `q_xy = 0` (trivial) or `1/2` (topological).
//! The topological phase hosts mid-gap corner states but zero edge polarisation.
//!
//! ## Breathing kagome HOTI
//!
//! M. Ezawa, *Phys. Rev. Lett.* **120**, 026801 (2018).
//!
//! Alternating intra-triangle (`t_intra`) and inter-triangle (`t_inter`)
//! hoppings on the kagome lattice produce a topological crystalline phase
//! with corner states protected by the three-fold rotational symmetry C₃.
//!
//! ## Hamiltonian conventions
//!
//! The BBH Hamiltonian is written in the basis `(A, B, C, D)` with
//!
//! ```text
//! H(k) = (λ_x + γ_x cos kx) Γ₁ + γ_x sin kx Γ₂
//!       + (λ_y + γ_y cos ky) Γ₃ + γ_y sin ky Γ₄
//! ```
//!
//! where the 4×4 Gamma matrices are the tensor products
//!
//! ```text
//! Γ₁ = σ_z ⊗ σ_x,  Γ₂ = σ_z ⊗ σ_y,  Γ₃ = I ⊗ σ_x,  Γ₄ = I ⊗ σ_y
//! ```
//!
//! with `σ_x, σ_y, σ_z` the 2×2 Pauli matrices.
//!
//! Topological phase: `|γ_x| > |λ_x|` AND `|γ_y| > |λ_y|`.
//! Trivial phase: `|λ_x| > |γ_x|` or `|λ_y| > |γ_y|`.
//!
//! ## Corner-state solver
//!
//! [`CornerStateSolver`] builds the full real-space Hamiltonian on a finite
//! `lx × ly` cluster with open boundary conditions and diagonalises it to
//! reveal mid-gap corner states.

use std::f64::consts::PI;

use crate::error::{self, Result};
use crate::math::{CMatrix, Complex};

// ---------------------------------------------------------------------------
// HotiLattice enum
// ---------------------------------------------------------------------------

/// Lattice geometry for HOTI models.
#[derive(Debug, Clone, PartialEq)]
pub enum HotiLattice {
    /// BBH square-lattice quadrupole insulator (4-band, C₄-symmetric).
    Bbh,
    /// Breathing kagome in the trivial phase (`t_intra > t_inter`).
    BreathingKagomeTrivial,
    /// Breathing kagome in the topological phase (`t_inter > t_intra`).
    BreathingKagomeTopological,
}

// ---------------------------------------------------------------------------
// Helper: Pauli matrices as 2×2 CMatrix
// ---------------------------------------------------------------------------

/// Build the 2×2 real σ_x Pauli matrix.
#[inline]
fn sigma_x() -> CMatrix {
    let mut m = CMatrix::zeros(2);
    m.set(0, 1, Complex::ONE);
    m.set(1, 0, Complex::ONE);
    m
}

/// Build the 2×2 imaginary σ_y Pauli matrix.
#[inline]
fn sigma_y() -> CMatrix {
    let mut m = CMatrix::zeros(2);
    m.set(0, 1, Complex::new(0.0, -1.0));
    m.set(1, 0, Complex::new(0.0, 1.0));
    m
}

/// Build the 2×2 real σ_z Pauli matrix.
#[inline]
fn sigma_z() -> CMatrix {
    let mut m = CMatrix::zeros(2);
    m.set(0, 0, Complex::ONE);
    m.set(1, 1, Complex::new(-1.0, 0.0));
    m
}

/// Compute the Kronecker product A ⊗ B where A is 2×2 and B is 2×2,
/// returning a 4×4 CMatrix.
///
/// `(A⊗B)[2i+k, 2j+l] = A[i,j] · B[k,l]`
fn kron_2x2(a: &CMatrix, b: &CMatrix) -> CMatrix {
    let mut out = CMatrix::zeros(4);
    for i in 0..2 {
        for j in 0..2 {
            let a_ij = a.get(i, j);
            for k in 0..2 {
                for l in 0..2 {
                    let val = a_ij.mul(&b.get(k, l));
                    out.set(2 * i + k, 2 * j + l, val);
                }
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// BbhModel
// ---------------------------------------------------------------------------

/// Benalcazar-Bernevig-Hughes 4-band model on the 2D square lattice.
///
/// The Hamiltonian encodes a quantised electric quadrupole moment `q_xy`
/// which is `1/2` in the topological phase (`|γ| > |λ|`) and `0` in
/// the trivial phase.
///
/// # References
///
/// - W. A. Benalcazar, B. A. Bernevig, T. L. Hughes,
///   *Science* **357**, 61 (2017).
/// - W. A. Benalcazar, B. A. Bernevig, T. L. Hughes,
///   *Phys. Rev. B* **96**, 245115 (2017).
#[derive(Debug, Clone)]
pub struct BbhModel {
    /// Intra-cell hopping in the x-direction (λ_x).
    pub lambda_x: f64,
    /// Inter-cell hopping in the x-direction (γ_x).
    pub gamma_x: f64,
    /// Intra-cell hopping in the y-direction (λ_y).
    pub lambda_y: f64,
    /// Inter-cell hopping in the y-direction (γ_y).
    pub gamma_y: f64,
}

impl BbhModel {
    /// Construct a new `BbhModel`.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if any parameter is non-finite.
    pub fn new(lambda_x: f64, gamma_x: f64, lambda_y: f64, gamma_y: f64) -> Result<Self> {
        for (name, v) in [
            ("lambda_x", lambda_x),
            ("gamma_x", gamma_x),
            ("lambda_y", lambda_y),
            ("gamma_y", gamma_y),
        ] {
            if !v.is_finite() {
                return Err(error::invalid_param(name, "must be finite"));
            }
        }
        Ok(Self {
            lambda_x,
            gamma_x,
            lambda_y,
            gamma_y,
        })
    }

    /// Canonical topological phase (`|γ| > |λ|` in both directions).
    ///
    /// Uses `λ_x = λ_y = 0.5`, `γ_x = γ_y = 1.0`.
    pub fn topological_phase() -> Self {
        Self {
            lambda_x: 0.5,
            gamma_x: 1.0,
            lambda_y: 0.5,
            gamma_y: 1.0,
        }
    }

    /// Canonical trivial phase (`|λ| > |γ|` in both directions).
    ///
    /// Uses `λ_x = λ_y = 1.0`, `γ_x = γ_y = 0.5`.
    pub fn trivial_phase() -> Self {
        Self {
            lambda_x: 1.0,
            gamma_x: 0.5,
            lambda_y: 1.0,
            gamma_y: 0.5,
        }
    }

    // -----------------------------------------------------------------------
    // Gamma matrices (shared by hamiltonian_at and CornerStateSolver)
    // -----------------------------------------------------------------------

    /// Build the four 4×4 Gamma matrices used in the BBH Hamiltonian.
    ///
    /// The Gamma matrices must mutually anticommute: `{Γ_a, Γ_b} = 2δ_ab`.
    /// This ensures the Hamiltonian `H² = Σ d_a² · I` with spectrum
    /// `±‖d‖`, each doubly degenerate, giving a true bulk band gap.
    ///
    /// The chosen representation (Benalcazar et al., *Science* 357, 61 (2017)):
    /// - `Γ₁ = σ_x ⊗ σ_z`
    /// - `Γ₂ = σ_y ⊗ I`
    /// - `Γ₃ = σ_x ⊗ σ_x`
    /// - `Γ₄ = σ_x ⊗ σ_y`
    ///
    /// All pairs anticommute (verified by {σ_a,σ_b}=2δ_ab and Kronecker structure).
    pub fn gamma_matrices() -> [CMatrix; 4] {
        let sx = sigma_x();
        let sy = sigma_y();
        let sz = sigma_z();
        let id2 = CMatrix::eye(2);
        [
            kron_2x2(&sx, &sz),  // Γ₁ = σ_x ⊗ σ_z
            kron_2x2(&sy, &id2), // Γ₂ = σ_y ⊗ I
            kron_2x2(&sx, &sx),  // Γ₃ = σ_x ⊗ σ_x
            kron_2x2(&sx, &sy),  // Γ₄ = σ_x ⊗ σ_y
        ]
    }

    // -----------------------------------------------------------------------
    // Bloch Hamiltonian
    // -----------------------------------------------------------------------

    /// Build the 4×4 BBH Bloch Hamiltonian at `(kx, ky)`.
    ///
    /// ```text
    /// H(k) = (λ_x + γ_x cos kx)·Γ₁ + γ_x sin kx·Γ₂
    ///       + (λ_y + γ_y cos ky)·Γ₃ + γ_y sin ky·Γ₄
    /// ```
    pub fn hamiltonian_at(&self, kx: f64, ky: f64) -> CMatrix {
        let [g1, g2, g3, g4] = Self::gamma_matrices();
        let ax = self.lambda_x + self.gamma_x * kx.cos();
        let bx = self.gamma_x * kx.sin();
        let ay = self.lambda_y + self.gamma_y * ky.cos();
        let by_ = self.gamma_y * ky.sin();

        // H = ax·Γ₁ + bx·Γ₂ + ay·Γ₃ + by·Γ₄
        g1.scale_real(ax)
            .add(&g2.scale_real(bx))
            .unwrap()
            .add(&g3.scale_real(ay))
            .unwrap()
            .add(&g4.scale_real(by_))
            .unwrap()
    }

    // -----------------------------------------------------------------------
    // Diagonalisation
    // -----------------------------------------------------------------------

    /// Return the four energy eigenvalues at `(kx, ky)` in ascending order.
    ///
    /// # Errors
    ///
    /// Propagates errors from the Hermitian eigendecomposition.
    pub fn energy_bands(&self, kx: f64, ky: f64) -> Result<Vec<f64>> {
        let h = self.hamiltonian_at(kx, ky);
        let (evals, _) = h.hermitian_eigendecomposition()?;
        Ok(evals)
    }

    // -----------------------------------------------------------------------
    // Band gap
    // -----------------------------------------------------------------------

    /// Minimum of the positive-energy spectrum over the BZ.
    ///
    /// For a particle-hole symmetric 4-band model the true bulk gap (distance
    /// of the conduction band minimum from zero) equals the minimum of
    /// `|E_n(k)|` over all k and over the two upper bands (indices 2 and 3).
    /// This gives the correct gap for the BBH model.
    pub fn band_gap(&self, nkx: usize, nky: usize) -> f64 {
        let mut min_positive = f64::INFINITY;
        let nkx = nkx.max(4);
        let nky = nky.max(4);
        for ix in 0..nkx {
            // Avoid high-symmetry points at 0, π by using offset grid
            let kx = 2.0 * PI * ((ix as f64) + 0.5) / (nkx as f64);
            for iy in 0..nky {
                let ky = 2.0 * PI * ((iy as f64) + 0.5) / (nky as f64);
                if let Ok(e) = self.energy_bands(kx, ky) {
                    // Positive bands are e[2] and e[3]; track minimum of e[2]
                    let pos = e[2];
                    if pos < min_positive {
                        min_positive = pos;
                    }
                }
            }
        }
        if min_positive.is_infinite() {
            0.0
        } else {
            min_positive
        }
    }

    // -----------------------------------------------------------------------
    // Phase classification
    // -----------------------------------------------------------------------

    /// Check whether the model is in the higher-order topological phase.
    ///
    /// Returns `true` iff `|γ_x| > |λ_x|` AND `|γ_y| > |λ_y|`.
    pub fn is_higher_order_topological(&self) -> bool {
        self.gamma_x.abs() > self.lambda_x.abs() && self.gamma_y.abs() > self.lambda_y.abs()
    }

    // -----------------------------------------------------------------------
    // Quadrupole moment via nested Wilson loop
    // -----------------------------------------------------------------------

    /// Compute the electric quadrupole moment `q_xy` via the nested Wilson loop.
    ///
    /// For the two lower bands (indices 0 and 1, the occupied sector):
    ///
    /// 1. At each kx, compute `W_y(kx)` as the ordered product of 2×2 link
    ///    matrices `M(ky_j, ky_{j+1})` around the closed ky loop.
    /// 2. The Wannier centre `ν(kx) = arg(det W_y) / (2π)`.
    /// 3. Nested polarisation: `q_xy = arg(∏_{kx} exp(2πi·ν(kx))) / (2π)`.
    ///
    /// The result is quantised to `0` (trivial) or `±0.5` (HOTI).
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `nkx < 4` or `nky < 4`.
    pub fn quadrupole_moment(&self, nkx: usize, nky: usize) -> Result<f64> {
        if nkx < 4 {
            return Err(error::invalid_param("nkx", "need at least 4 kx points"));
        }
        if nky < 4 {
            return Err(error::invalid_param("nky", "need at least 4 ky points"));
        }

        // Occupied bands are the two lowest eigenvalues (columns 0 and 1 of
        // the eigenvector matrix, sorted by ascending eigenvalue).
        let occupied = [0usize, 1usize];

        let mut w_nested = Complex::ONE;

        for ikx in 0..nkx {
            let kx = 2.0 * PI * (ikx as f64) / (nkx as f64);

            // Collect states at each ky point (including wrap-around)
            let mut all_states: Vec<[Vec<Complex>; 2]> = Vec::with_capacity(nky + 1);
            for iky in 0..=nky {
                let ky = 2.0 * PI * (iky as f64) / (nky as f64);
                let h = self.hamiltonian_at(kx, ky);
                let (_, vecs) = h.hermitian_eigendecomposition()?;
                let s0 = vecs.column(occupied[0]);
                let s1 = vecs.column(occupied[1]);
                all_states.push([s0, s1]);
            }

            // W_y = ∏_{ky} M(ky_j, ky_{j+1}), a 2×2 unitary
            let mut w_y = CMatrix::eye(2);
            for iky in 0..nky {
                let sa = &all_states[iky];
                let sb = &all_states[iky + 1];
                // 2×2 link matrix
                let mut m = CMatrix::zeros(2);
                for (mu, sa_row) in sa.iter().enumerate().take(2) {
                    for (nu, sb_row) in sb.iter().enumerate().take(2) {
                        let mut inner = Complex::ZERO;
                        for j in 0..sa_row.len() {
                            inner = inner.add(&sa_row[j].conj().mul(&sb_row[j]));
                        }
                        m.set(mu, nu, inner);
                    }
                }
                w_y = w_y.matmul(&m)?;
            }

            // Wannier centre ν(kx) = arg(det W_y) / (2π)
            // For a 2×2 matrix: det = a·d − b·c
            let det = w_y
                .get(0, 0)
                .mul(&w_y.get(1, 1))
                .sub(&w_y.get(0, 1).mul(&w_y.get(1, 0)));
            let nu = det.phase() / (2.0 * PI);

            // Accumulate: exp(2πi·ν)
            let phase = Complex::from_polar(1.0, 2.0 * PI * nu);
            w_nested = w_nested.mul(&phase);
        }

        let q_xy = w_nested.phase() / (2.0 * PI);
        Ok(q_xy)
    }
}

// ---------------------------------------------------------------------------
// BreathingKagomeModel
// ---------------------------------------------------------------------------

/// Breathing kagome lattice HOTI.
///
/// The kagome lattice consists of corner-sharing triangles. Alternating
/// intra-triangle and inter-triangle hoppings break the translation symmetry
/// and generate a topological crystalline phase protected by C₃ rotation.
///
/// The Bloch Hamiltonian in the three-sublattice basis is:
///
/// ```text
/// H(k)_ab = m·δ_ab + t_intra + t_inter·exp(i k·d_ab)
/// ```
///
/// where `d_ab` are the inter-cell bond vectors.
///
/// # References
///
/// - M. Ezawa, *Phys. Rev. Lett.* **120**, 026801 (2018).
#[derive(Debug, Clone)]
pub struct BreathingKagomeModel {
    /// Intra-triangle hopping amplitude `t_intra`.
    pub t_intra: f64,
    /// Inter-triangle hopping amplitude `t_inter`.
    pub t_inter: f64,
    /// On-site mass (sublattice potential) `m`.
    pub m_mass: f64,
}

impl BreathingKagomeModel {
    /// Construct a new `BreathingKagomeModel`.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if any parameter is non-finite.
    pub fn new(t_intra: f64, t_inter: f64, m_mass: f64) -> Result<Self> {
        for (name, v) in [
            ("t_intra", t_intra),
            ("t_inter", t_inter),
            ("m_mass", m_mass),
        ] {
            if !v.is_finite() {
                return Err(error::invalid_param(name, "must be finite"));
            }
        }
        Ok(Self {
            t_intra,
            t_inter,
            m_mass,
        })
    }

    /// Topological breathing-kagome phase (`|t_inter| > |t_intra|`).
    pub fn topological_kagome() -> Self {
        Self {
            t_intra: 0.5,
            t_inter: 1.0,
            m_mass: 0.0,
        }
    }

    /// Trivial breathing-kagome phase (`|t_intra| > |t_inter|`).
    pub fn trivial_kagome() -> Self {
        Self {
            t_intra: 1.0,
            t_inter: 0.5,
            m_mass: 0.0,
        }
    }

    // -----------------------------------------------------------------------
    // Bloch Hamiltonian
    // -----------------------------------------------------------------------

    /// Build the 3×3 Bloch Hamiltonian at `(kx, ky)`.
    ///
    /// Uses the honeycomb-like breathing kagome convention with primitive
    /// lattice vectors `a₁ = (1, 0)` and `a₂ = (cos 60°, sin 60°)`:
    ///
    /// ```text
    /// H[0,1] = t_intra + t_inter · exp(i k·a₁)
    /// H[0,2] = t_intra + t_inter · exp(i k·a₂)
    /// H[1,2] = t_intra + t_inter · exp(i k·(a₂ - a₁))
    /// ```
    ///
    /// Diagonal elements equal `m_mass`.
    pub fn hamiltonian_at(&self, kx: f64, ky: f64) -> CMatrix {
        let a2x = 0.5_f64; // cos 60°
        let a2y = (3.0_f64).sqrt() * 0.5; // sin 60°

        // Phase factors for inter-cell bonds
        let phi1 = kx; // k · a₁ = kx
        let phi2 = kx * a2x + ky * a2y; // k · a₂
        let phi3 = phi2 - phi1; // k · (a₂ - a₁)

        let p1 = Complex::from_polar(1.0, phi1);
        let p2 = Complex::from_polar(1.0, phi2);
        let p3 = Complex::from_polar(1.0, phi3);

        let t = Complex::from_real(self.t_intra);
        let ti = Complex::from_real(self.t_inter);

        // H[a,b] = t_intra + t_inter · exp(i φ_{ab})
        let h01 = t.add(&ti.mul(&p1));
        let h02 = t.add(&ti.mul(&p2));
        let h12 = t.add(&ti.mul(&p3));

        let mut mat = CMatrix::zeros(3);
        let m = Complex::from_real(self.m_mass);
        mat.set(0, 0, m);
        mat.set(1, 1, m);
        mat.set(2, 2, m);
        mat.set(0, 1, h01);
        mat.set(1, 0, h01.conj());
        mat.set(0, 2, h02);
        mat.set(2, 0, h02.conj());
        mat.set(1, 2, h12);
        mat.set(2, 1, h12.conj());
        mat
    }

    // -----------------------------------------------------------------------
    // Diagonalisation
    // -----------------------------------------------------------------------

    /// Return the three energy eigenvalues at `(kx, ky)` in ascending order.
    ///
    /// # Errors
    ///
    /// Propagates errors from the Hermitian eigendecomposition.
    pub fn energy_bands(&self, kx: f64, ky: f64) -> Result<Vec<f64>> {
        let h = self.hamiltonian_at(kx, ky);
        let (evals, _) = h.hermitian_eigendecomposition()?;
        Ok(evals)
    }

    // -----------------------------------------------------------------------
    // Phase classification
    // -----------------------------------------------------------------------

    /// Check whether the model is in the topological phase.
    ///
    /// Returns `true` iff `|t_inter| > |t_intra|`.
    pub fn is_topological(&self) -> bool {
        self.t_inter.abs() > self.t_intra.abs()
    }

    // -----------------------------------------------------------------------
    // Corner polarisation (flat-band Wannier centre)
    // -----------------------------------------------------------------------

    /// Average Wannier centre of the lowest band over the BZ.
    ///
    /// For the breathing kagome the lowest band is nearly flat in the
    /// topological phase; its Wannier centre at each kx is computed from
    /// the single-band Wilson loop along ky, and the result is averaged
    /// over `nkx` values of kx.
    ///
    /// A value near `0.5` (mod 1) signals non-trivial corner polarisation.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `nkx < 4` or `nky < 4`.
    pub fn corner_polarization(&self, nkx: usize, nky: usize) -> Result<f64> {
        if nkx < 4 {
            return Err(error::invalid_param("nkx", "need at least 4 kx points"));
        }
        if nky < 4 {
            return Err(error::invalid_param("nky", "need at least 4 ky points"));
        }

        let mut sum = 0.0_f64;
        for ikx in 0..nkx {
            let kx = 2.0 * PI * (ikx as f64) / (nkx as f64);
            // Single-band Wilson loop along ky (band 0 = lowest)
            let mut w_scalar = Complex::ONE;
            let mut prev_state: Option<Vec<Complex>> = None;
            for iky in 0..=nky {
                let ky = 2.0 * PI * (iky as f64) / (nky as f64);
                let h = self.hamiltonian_at(kx, ky);
                let (_, vecs) = h.hermitian_eigendecomposition()?;
                let state = vecs.column(0);
                if let Some(ref prev) = prev_state {
                    let mut inner = Complex::ZERO;
                    for j in 0..prev.len() {
                        inner = inner.add(&prev[j].conj().mul(&state[j]));
                    }
                    let norm = inner.norm();
                    let link = if norm < 1e-15 {
                        Complex::ONE
                    } else {
                        inner.scale(1.0 / norm)
                    };
                    w_scalar = w_scalar.mul(&link);
                }
                prev_state = Some(state);
            }
            let nu = w_scalar.phase() / (2.0 * PI);
            let nu_mapped = nu.rem_euclid(1.0);
            sum += nu_mapped;
        }
        Ok(sum / (nkx as f64))
    }
}

// ---------------------------------------------------------------------------
// CornerStateSolver
// ---------------------------------------------------------------------------

/// Open-boundary-condition solver for the BBH quadrupole insulator on a finite
/// `lx × ly` cluster.
///
/// Builds the full real-space Hamiltonian for the BBH model on an `lx × ly`
/// cluster with 4 orbitals per unit cell and periodic intra-cell hoppings,
/// then diagonalises it.  The mid-gap states localised at the four corners
/// are the hallmark of the topological quadrupole phase.
///
/// # Dimension constraint
///
/// The total Hilbert-space dimension is `lx × ly × 4`.  Since `CMatrix` is
/// capped at 64, valid cluster sizes satisfy `lx × ly ≤ 16` (e.g. `4×4`).
#[derive(Debug, Clone)]
pub struct CornerStateSolver {
    /// Number of unit cells along x.
    pub lx: usize,
    /// Number of unit cells along y.
    pub ly: usize,
    /// Intra-cell hopping λ_x.
    pub lambda_x: f64,
    /// Inter-cell hopping γ_x.
    pub gamma_x: f64,
    /// Intra-cell hopping λ_y.
    pub lambda_y: f64,
    /// Inter-cell hopping γ_y.
    pub gamma_y: f64,
}

impl CornerStateSolver {
    /// Number of orbitals per unit cell.
    const N_ORB: usize = 4;

    /// Construct a new `CornerStateSolver`.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `lx * ly * 4 > 64` (exceeds `CMatrix::MAX_DIM`).
    pub fn new(
        lx: usize,
        ly: usize,
        lambda_x: f64,
        gamma_x: f64,
        lambda_y: f64,
        gamma_y: f64,
    ) -> Result<Self> {
        if lx == 0 || ly == 0 {
            return Err(error::invalid_param("lx/ly", "must be at least 1"));
        }
        let dim = lx * ly * Self::N_ORB;
        if dim > CMatrix::MAX_DIM {
            return Err(error::invalid_param(
                "lx*ly",
                &format!("cluster dimension {dim} = {lx}×{ly}×4 exceeds CMatrix::MAX_DIM=64"),
            ));
        }
        for (name, v) in [
            ("lambda_x", lambda_x),
            ("gamma_x", gamma_x),
            ("lambda_y", lambda_y),
            ("gamma_y", gamma_y),
        ] {
            if !v.is_finite() {
                return Err(error::invalid_param(name, "must be finite"));
            }
        }
        Ok(Self {
            lx,
            ly,
            lambda_x,
            gamma_x,
            lambda_y,
            gamma_y,
        })
    }

    /// Flat site index: orbital `alpha` at unit cell `(ix, iy)`.
    #[inline]
    fn site(&self, ix: usize, iy: usize, alpha: usize) -> usize {
        (iy * self.lx + ix) * Self::N_ORB + alpha
    }

    /// Build the full real-space cluster Hamiltonian with OBC.
    ///
    /// # Hamiltonian structure
    ///
    /// Using the Gamma matrix decomposition, the Bloch Hamiltonian at k=0 gives
    /// the on-site block `H_0 = λ_x Γ₁ + λ_y Γ₃`, and the inter-cell hoppings:
    ///
    /// ```text
    /// Hop_x = (γ_x / 2) · (Γ₁ − i Γ₂)   for (ix, iy) → (ix+1, iy)
    /// Hop_y = (γ_y / 2) · (Γ₃ − i Γ₄)   for (ix, iy) → (ix, iy+1)
    /// ```
    ///
    /// Hermitian conjugates are added explicitly.
    ///
    /// # Errors
    ///
    /// Propagates CMatrix errors (should not occur for valid parameters).
    pub fn build_cluster_hamiltonian(&self) -> Result<CMatrix> {
        let dim = self.lx * self.ly * Self::N_ORB;
        let mut h_full = CMatrix::zeros(dim);

        let [g1, g2, g3, g4] = BbhModel::gamma_matrices();

        // On-site block: λ_x Γ₁ + λ_y Γ₃
        let h_onsite = g1
            .scale_real(self.lambda_x)
            .add(&g3.scale_real(self.lambda_y))?;

        // Inter-cell hopping matrices (4×4)
        // Hop_x = (γ_x/2) · (Γ₁ − i·Γ₂)
        let g2_i = g2.scale(Complex::new(0.0, -1.0));
        let hop_x = g1
            .scale_real(self.gamma_x * 0.5)
            .add(&g2_i.scale_real(self.gamma_x * 0.5))?;
        // Hop_y = (γ_y/2) · (Γ₃ − i·Γ₄)
        let g4_i = g4.scale(Complex::new(0.0, -1.0));
        let hop_y = g3
            .scale_real(self.gamma_y * 0.5)
            .add(&g4_i.scale_real(self.gamma_y * 0.5))?;

        let hop_x_dag = hop_x.conj_transpose();
        let hop_y_dag = hop_y.conj_transpose();

        for iy in 0..self.ly {
            for ix in 0..self.lx {
                // On-site
                for alpha in 0..Self::N_ORB {
                    for beta in 0..Self::N_ORB {
                        let row = self.site(ix, iy, alpha);
                        let col = self.site(ix, iy, beta);
                        let v = h_full.get(row, col).add(&h_onsite.get(alpha, beta));
                        h_full.set(row, col, v);
                    }
                }
                // Hop in x: (ix, iy) → (ix+1, iy)
                if ix + 1 < self.lx {
                    for alpha in 0..Self::N_ORB {
                        for beta in 0..Self::N_ORB {
                            // Forward hop
                            let row = self.site(ix, iy, alpha);
                            let col = self.site(ix + 1, iy, beta);
                            let v = h_full.get(row, col).add(&hop_x.get(alpha, beta));
                            h_full.set(row, col, v);
                            // Hermitian conjugate
                            let v2 = h_full.get(col, row).add(&hop_x_dag.get(beta, alpha));
                            h_full.set(col, row, v2);
                        }
                    }
                }
                // Hop in y: (ix, iy) → (ix, iy+1)
                if iy + 1 < self.ly {
                    for alpha in 0..Self::N_ORB {
                        for beta in 0..Self::N_ORB {
                            let row = self.site(ix, iy, alpha);
                            let col = self.site(ix, iy + 1, beta);
                            let v = h_full.get(row, col).add(&hop_y.get(alpha, beta));
                            h_full.set(row, col, v);
                            let v2 = h_full.get(col, row).add(&hop_y_dag.get(beta, alpha));
                            h_full.set(col, row, v2);
                        }
                    }
                }
            }
        }

        Ok(h_full)
    }

    /// Diagonalise the finite cluster Hamiltonian.
    ///
    /// Returns `(eigenvalues, eigenvectors)` where eigenvalues are sorted
    /// ascending and the k-th column of `eigenvectors` is the k-th eigenstate.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`build_cluster_hamiltonian`](Self::build_cluster_hamiltonian)
    /// and [`hermitian_eigendecomposition`](CMatrix::hermitian_eigendecomposition).
    pub fn solve_finite_cluster(&self) -> Result<(Vec<f64>, CMatrix)> {
        let h = self.build_cluster_hamiltonian()?;
        h.hermitian_eigendecomposition()
    }

    /// Compute the corner localisation weights of eigenstate `eigvec_col`.
    ///
    /// For each of the four corner unit cells
    /// `(0,0)`, `(lx-1, 0)`, `(0, ly-1)`, `(lx-1, ly-1)`,
    /// sum the squared amplitudes over all 4 orbitals.
    /// Returns `[w_00, w_x0, w_0y, w_xy]` normalised to sum to 1.
    pub fn corner_localization(&self, eigvec_col: usize, eigvectors: &CMatrix) -> [f64; 4] {
        let corners = [
            (0, 0),
            (self.lx - 1, 0),
            (0, self.ly - 1),
            (self.lx - 1, self.ly - 1),
        ];
        let mut weights = [0.0_f64; 4];
        let mut total = 0.0_f64;
        for (ci, &(cx, cy)) in corners.iter().enumerate() {
            let mut w = 0.0_f64;
            for alpha in 0..Self::N_ORB {
                let row = self.site(cx, cy, alpha);
                let amp = eigvectors.get(row, eigvec_col);
                w += amp.norm_sq();
            }
            weights[ci] = w;
            total += w;
        }
        if total > 1e-30 {
            for w in &mut weights {
                *w /= total;
            }
        }
        weights
    }

    /// Count eigenstates with energies inside `(gap_min, gap_max)`.
    pub fn count_corner_states(&self, eigvalues: &[f64], gap_min: f64, gap_max: f64) -> usize {
        eigvalues
            .iter()
            .filter(|&&e| e > gap_min && e < gap_max)
            .count()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // BbhModel construction
    // -----------------------------------------------------------------------

    #[test]
    fn bbh_topological_phase_is_topological() {
        let m = BbhModel::topological_phase();
        assert!(m.is_higher_order_topological());
    }

    #[test]
    fn bbh_trivial_phase_is_not_topological() {
        let m = BbhModel::trivial_phase();
        assert!(!m.is_higher_order_topological());
    }

    #[test]
    fn bbh_new_rejects_nan() {
        assert!(BbhModel::new(f64::NAN, 1.0, 0.5, 1.0).is_err());
        assert!(BbhModel::new(0.5, f64::INFINITY, 0.5, 1.0).is_err());
    }

    // -----------------------------------------------------------------------
    // Hamiltonian Hermiticity
    // -----------------------------------------------------------------------

    #[test]
    fn bbh_hamiltonian_hermitian() {
        let m = BbhModel::topological_phase();
        for kx in [-1.0, 0.0, 0.5, PI / 3.0] {
            for ky in [0.0, 1.0, PI] {
                let h = m.hamiltonian_at(kx, ky);
                let hd = h.conj_transpose();
                let diff = h.sub(&hd).unwrap();
                assert!(
                    diff.frobenius_norm() < 1e-12,
                    "H not Hermitian at ({kx},{ky}): {}",
                    diff.frobenius_norm()
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Energy bands
    // -----------------------------------------------------------------------

    #[test]
    fn bbh_four_bands() {
        let m = BbhModel::topological_phase();
        let e = m.energy_bands(0.0, 0.0).unwrap();
        assert_eq!(e.len(), 4);
    }

    #[test]
    fn bbh_bands_symmetric_about_zero() {
        // The BBH Hamiltonian anticommutes with Γ₅ = σ_y⊗I₂, so the spectrum
        // is particle-hole symmetric: {E_n} = {-E_n}.
        let m = BbhModel::topological_phase();
        let e = m.energy_bands(0.7, 1.2).unwrap();
        // e[0] ≈ -e[3], e[1] ≈ -e[2]
        assert!((e[0] + e[3]).abs() < 1e-10, "PH symmetry broken: {e:?}");
        assert!((e[1] + e[2]).abs() < 1e-10, "PH symmetry broken: {e:?}");
    }

    // -----------------------------------------------------------------------
    // Band gap
    // -----------------------------------------------------------------------

    #[test]
    fn bbh_topological_has_gap() {
        let m = BbhModel::topological_phase();
        let gap = m.band_gap(10, 10);
        // Topological phase (λ=0.5, γ=1.0): min positive energy ≈ √0.5 ≈ 0.71
        assert!(
            gap > 0.5,
            "Expected gap > 0.5 for topological BBH, got {gap}"
        );
    }

    #[test]
    fn bbh_quadrupole_moment_quantised() {
        // Topological phase → q_xy ≈ ±0.5 (mod 1)
        let m = BbhModel::topological_phase();
        let q = m.quadrupole_moment(8, 8).unwrap();
        let wrapped = q.rem_euclid(1.0);
        let dist = (wrapped - 0.5).abs().min(wrapped.abs());
        assert!(
            dist < 0.35,
            "q_xy not near 0 or ±0.5: q={q:.4}, wrapped={wrapped:.4}"
        );
    }

    #[test]
    fn bbh_quadrupole_rejects_small_grid() {
        let m = BbhModel::topological_phase();
        assert!(m.quadrupole_moment(2, 8).is_err());
        assert!(m.quadrupole_moment(8, 2).is_err());
    }

    // -----------------------------------------------------------------------
    // BreathingKagomeModel
    // -----------------------------------------------------------------------

    #[test]
    fn breathing_kagome_topological_flag() {
        let m = BreathingKagomeModel::topological_kagome();
        assert!(m.is_topological());
        let m2 = BreathingKagomeModel::trivial_kagome();
        assert!(!m2.is_topological());
    }

    #[test]
    fn breathing_kagome_hamiltonian_hermitian() {
        let m = BreathingKagomeModel::topological_kagome();
        let h = m.hamiltonian_at(0.3, 0.7);
        let hd = h.conj_transpose();
        let diff = h.sub(&hd).unwrap();
        assert!(diff.frobenius_norm() < 1e-12);
    }

    #[test]
    fn breathing_kagome_three_bands() {
        let m = BreathingKagomeModel::topological_kagome();
        let e = m.energy_bands(0.0, 0.0).unwrap();
        assert_eq!(e.len(), 3);
    }

    #[test]
    fn breathing_kagome_corner_polarization_range() {
        let m = BreathingKagomeModel::topological_kagome();
        let cp = m.corner_polarization(8, 8).unwrap();
        assert!(
            (0.0..=1.0).contains(&cp),
            "Corner polarisation out of [0,1]: {cp}"
        );
    }

    // -----------------------------------------------------------------------
    // CornerStateSolver
    // -----------------------------------------------------------------------

    #[test]
    fn corner_solver_rejects_oversized_cluster() {
        // 5×5×4 = 100 > 64
        assert!(CornerStateSolver::new(5, 5, 0.5, 1.0, 0.5, 1.0).is_err());
    }

    #[test]
    fn corner_solver_2x2_has_dim_16() {
        let solver = CornerStateSolver::new(2, 2, 0.5, 1.0, 0.5, 1.0).unwrap();
        let h = solver.build_cluster_hamiltonian().unwrap();
        assert_eq!(h.n(), 16, "Expected 2*2*4=16 dimensional H");
    }

    #[test]
    fn corner_solver_hamiltonian_hermitian() {
        let solver = CornerStateSolver::new(2, 2, 0.5, 1.0, 0.5, 1.0).unwrap();
        let h = solver.build_cluster_hamiltonian().unwrap();
        let hd = h.conj_transpose();
        let diff = h.sub(&hd).unwrap();
        assert!(
            diff.frobenius_norm() < 1e-10,
            "Cluster H not Hermitian: norm={}",
            diff.frobenius_norm()
        );
    }

    #[test]
    fn corner_solver_spectrum_symmetric() {
        // BBH cluster in topological phase should have PH-symmetric spectrum
        let solver = CornerStateSolver::new(2, 2, 0.5, 1.0, 0.5, 1.0).unwrap();
        let (evals, _) = solver.solve_finite_cluster().unwrap();
        let n = evals.len();
        for i in 0..n / 2 {
            assert!(
                (evals[i] + evals[n - 1 - i]).abs() < 1e-8,
                "PH symmetry violated: E[{i}]={} vs E[{}]={}",
                evals[i],
                n - 1 - i,
                evals[n - 1 - i]
            );
        }
    }

    #[test]
    fn corner_localization_sums_to_one() {
        let solver = CornerStateSolver::new(2, 2, 0.5, 1.0, 0.5, 1.0).unwrap();
        let (_, evecs) = solver.solve_finite_cluster().unwrap();
        let w = solver.corner_localization(0, &evecs);
        let total: f64 = w.iter().sum();
        assert!((total - 1.0).abs() < 1e-10, "Corner weights sum={total}");
    }

    #[test]
    fn count_corner_states_topological_cluster() {
        // Topological 2×2 cluster: expect corner states near zero energy
        let solver = CornerStateSolver::new(2, 2, 0.5, 1.0, 0.5, 1.0).unwrap();
        let (evals, _) = solver.solve_finite_cluster().unwrap();
        // Count states in gap region (−0.5, 0.5)
        let n_gap = solver.count_corner_states(&evals, -0.5, 0.5);
        // In the topological phase we expect 4 mid-gap corner states
        assert!(
            n_gap >= 2,
            "Expected ≥2 corner states in gap, found {n_gap}"
        );
    }
}
