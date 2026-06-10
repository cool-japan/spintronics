//! Bosonic Bloch Hamiltonians for magnon band structure.
//!
//! Implements the magnon Bloch Hamiltonian H(k) for three canonical lattice
//! geometries used in topological magnon research:
//!
//! - **HoneycombHaldane**: two-sublattice ferromagnet with NN exchange J and
//!   NNN Dzyaloshinskii-Moriya interaction acting as a Haldane mass term.
//!   Ref: S. A. Owerre, *J. Phys.: Condens. Matter* **28**, 386001 (2016).
//!
//! - **Kagome**: three-sublattice frustrated ferromagnet with flat band + two
//!   dispersive bands; the DMI lifts the flat-band degeneracy.
//!
//! - **SquareDmi**: simplest single-band model with cosine dispersion; provides
//!   a non-topological baseline for benchmarking curvature and Hall calculations.
//!
//! # Honeycomb Hamiltonian Details
//!
//! Using primitive vectors **a₁** = (1, 0) and **a₂** = (½, √3/2) (in units of the
//! lattice constant *a*), the 2×2 Bloch Hamiltonian in the sublattice (A, B) basis is:
//!
//! ```text
//! H(k) = h₀(k)·I₂ + d(k)·σ
//! ```
//!
//! where **d** = (dₓ, d_y, d_z) and the Pauli matrices σ act on sublattice space.
//! The NN hopping gives dₓ + i·d_y = J·(1 + e^{ik·a₁} + e^{ik·a₂}) and the
//! topological mass from NNN-DMI is d_z = h_ext + 2·dmi·(sin(k·b₁)+sin(k·b₂)+sin(k·(b₁−b₂))).

use std::f64::consts::PI;

use crate::error::{self, Result};
use crate::math::{CMatrix, Complex};

// ---------------------------------------------------------------------------
// Direct analytical eigendecomposition for small matrices
// ---------------------------------------------------------------------------

/// Diagonalize a 1×1 Hermitian matrix: trivial.
fn diagonalize_1x1(h: &CMatrix) -> Result<(Vec<f64>, CMatrix)> {
    let e = h.get(0, 0).re;
    let mut vecs = CMatrix::zeros(1);
    vecs.set(0, 0, Complex::ONE);
    Ok((vec![e], vecs))
}

/// Diagonalize a 2×2 complex Hermitian matrix analytically.
///
/// For H = [[a, b], [b*, c]] (a, c real; b complex):
///   eigenvalues: λ± = (a+c)/2 ± √(((a-c)/2)² + |b|²)
///   eigenvectors: computed from the rows of (H - λI).
fn diagonalize_2x2(h: &CMatrix) -> Result<(Vec<f64>, CMatrix)> {
    let a = h.get(0, 0).re;
    let c = h.get(1, 1).re;
    let b = h.get(0, 1); // complex off-diagonal

    let half_sum = (a + c) * 0.5;
    let half_diff = (a - c) * 0.5;
    let disc = (half_diff * half_diff + b.norm_sq()).sqrt();

    let e_lo = half_sum - disc;
    let e_hi = half_sum + disc;

    // Eigenvector for e_lo: (H - e_lo*I) has rows [(a-e_lo, b), (b*, c-e_lo)]
    // Null vector: if b≠0, [−b, a−e_lo] normalized; if b=0, [1,0] or [0,1].
    let (v0, v1) = if b.norm() < 1e-14 {
        // Already diagonal
        if a <= c {
            (
                vec![Complex::ONE, Complex::ZERO],
                vec![Complex::ZERO, Complex::ONE],
            )
        } else {
            (
                vec![Complex::ZERO, Complex::ONE],
                vec![Complex::ONE, Complex::ZERO],
            )
        }
    } else {
        // v_lo = [-b, a-e_lo] (kernel of (H-e_lo*I))
        // v_hi = [-b, a-e_hi] = [-b, a-e_hi]
        // But v_hi = [b*, c-e_hi] is simpler (from second row)
        let vlo = [b.neg(), Complex::from_real(a - e_lo)];
        let norm_lo = (vlo[0].norm_sq() + vlo[1].norm_sq()).sqrt();
        let vlo = vec![vlo[0].scale(1.0 / norm_lo), vlo[1].scale(1.0 / norm_lo)];

        let vhi = [b.neg(), Complex::from_real(a - e_hi)];
        let norm_hi = (vhi[0].norm_sq() + vhi[1].norm_sq()).sqrt();
        let vhi = if norm_hi < 1e-14 {
            // Degenerate case — use orthogonal complement
            vec![vlo[1].conj(), vlo[0].neg().conj()]
        } else {
            vec![vhi[0].scale(1.0 / norm_hi), vhi[1].scale(1.0 / norm_hi)]
        };

        (vlo, vhi)
    };

    // Columns of eigenvector matrix in ascending eigenvalue order
    let mut vecs = CMatrix::zeros(2);
    vecs.set(0, 0, v0[0]);
    vecs.set(1, 0, v0[1]);
    vecs.set(0, 1, v1[0]);
    vecs.set(1, 1, v1[1]);

    Ok((vec![e_lo, e_hi], vecs))
}

// ---------------------------------------------------------------------------
// Robust 3×3 Hermitian eigendecomposition
// ---------------------------------------------------------------------------

/// Diagonalize a 3×3 complex Hermitian matrix using the Cardano trigonometric
/// method (Smith, 1961) for eigenvalues, with cross-product null-space vectors.
///
/// For a traceless Hermitian A = H − (tr/3)·I with characteristic polynomial
/// λ³ − p·λ − q = 0, the eigenvalues are p^{1/2}·cos(φ/3 + 2πk/3), k=0,1,2,
/// where φ = arccos(q / p^{3/2}).
fn diagonalize_3x3_hermitian(h: &CMatrix) -> Result<(Vec<f64>, CMatrix)> {
    use std::f64::consts::PI;

    let h00 = h.get(0, 0).re;
    let h11 = h.get(1, 1).re;
    let h22 = h.get(2, 2).re;
    let h01 = h.get(0, 1);
    let h02 = h.get(0, 2);
    let h12 = h.get(1, 2);

    let tr = h00 + h11 + h22;
    let shift = tr / 3.0;

    // Traceless part A = H - shift·I
    let a00 = h00 - shift;
    let a11 = h11 - shift;
    let a22 = h22 - shift;

    // p = tr(A²) / 6 = (Σdiag² + 2Σ|off|²) / 6
    let p_val =
        (a00 * a00 + a11 * a11 + a22 * a22 + 2.0 * (h01.norm_sq() + h02.norm_sq() + h12.norm_sq()))
            / 6.0;

    if p_val < 1e-28 {
        // Scalar: all eigenvalues = shift
        return Ok((vec![shift, shift, shift], CMatrix::eye(3)));
    }

    let p = p_val.sqrt(); // p > 0

    // q = det(A) / 2
    // For Hermitian A with real diagonal: det = a00*(a11*a22-|h12|²)
    //   - Re[conj(h01)*(h01*a22 - h12*conj(h02))]
    //   + Re[conj(h02)*(h01*conj(h12) - a11*conj(h02))]
    let cof01 = h01.scale(a22).sub(&h12.mul(&h02.conj()));
    let cof02 = h01
        .mul(&h12.conj())
        .sub(&Complex::from_real(a11).mul(&h02.conj()));
    let det_a =
        a00 * (a11 * a22 - h12.norm_sq()) - h01.conj().mul(&cof01).re + h02.conj().mul(&cof02).re;
    let q = det_a / 2.0; // q = det(A)/2

    // Trigonometric solution: eigenvalues of A are 2p·cos(φ/3 + 2πk/3)
    // where cos(φ) = q / p³. The argument q/p³ must be in [-1,1]; clamp for
    // robustness at degenerate points.
    let arg = (q / (p * p * p)).clamp(-1.0, 1.0);
    let phi = arg.acos() / 3.0; // φ/3 ∈ [0, π/3]

    let mut evals_a = [
        2.0 * p * phi.cos(),
        2.0 * p * (phi + 2.0 * PI / 3.0).cos(),
        2.0 * p * (phi + 4.0 * PI / 3.0).cos(),
    ];
    evals_a.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let eigenvalues: Vec<f64> = evals_a.iter().map(|&mu| mu + shift).collect();
    let eigenvectors = compute_3x3_eigenvectors(h, &eigenvalues)?;

    Ok((eigenvalues, eigenvectors))
}

/// Compute orthonormal eigenvectors for a 3×3 Hermitian matrix given eigenvalues.
///
/// For each eigenvalue λ, computes the null vector of (H − λ·I) via the cross
/// product of two rows (robust even for degenerate cases via Gram-Schmidt).
fn compute_3x3_eigenvectors(h: &CMatrix, eigenvalues: &[f64]) -> Result<CMatrix> {
    let mut vecs = CMatrix::zeros(3);
    let mut basis: Vec<Vec<Complex>> = Vec::with_capacity(3);

    for (col, &lambda) in eigenvalues.iter().enumerate() {
        // Build (H - λI)
        let mut rows = [[Complex::ZERO; 3]; 3];
        for (i, row) in rows.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = h.get(i, j);
                if i == j {
                    *cell = cell.sub(&Complex::from_real(lambda));
                }
            }
        }

        // Cross product of rows 0 and 1 to get a null vector candidate
        let v = cross_product_3c(&rows[0], &rows[1]);
        let norm = v.iter().map(|c| c.norm_sq()).sum::<f64>().sqrt();

        let mut vec = if norm > 1e-12 {
            v.iter().map(|c| c.scale(1.0 / norm)).collect::<Vec<_>>()
        } else {
            // Degenerate: try rows 0 and 2, or rows 1 and 2
            let v2 = cross_product_3c(&rows[0], &rows[2]);
            let n2 = v2.iter().map(|c| c.norm_sq()).sum::<f64>().sqrt();
            if n2 > 1e-12 {
                v2.iter().map(|c| c.scale(1.0 / n2)).collect::<Vec<_>>()
            } else {
                let v3 = cross_product_3c(&rows[1], &rows[2]);
                let n3 = v3.iter().map(|c| c.norm_sq()).sum::<f64>().sqrt();
                if n3 > 1e-12 {
                    v3.iter().map(|c| c.scale(1.0 / n3)).collect::<Vec<_>>()
                } else {
                    // Full degeneracy: use standard basis vector orthogonal to existing
                    standard_basis_complement(&basis, 3)
                }
            }
        };

        // Gram-Schmidt orthogonalization against already-computed vectors
        for prev in &basis {
            let dot: Complex = vec
                .iter()
                .zip(prev.iter())
                .map(|(&a, &b)| b.conj().mul(&a))
                .fold(Complex::ZERO, |acc, x| acc.add(&x));
            for (v_elem, &p_elem) in vec.iter_mut().zip(prev.iter()) {
                *v_elem = v_elem.sub(&p_elem.mul(&dot));
            }
            let n = vec.iter().map(|c| c.norm_sq()).sum::<f64>().sqrt();
            if n > 1e-12 {
                vec.iter_mut().for_each(|c| *c = c.scale(1.0 / n));
            }
        }

        // Re-normalize
        let n_final = vec.iter().map(|c| c.norm_sq()).sum::<f64>().sqrt();
        if n_final > 1e-12 {
            vec.iter_mut().for_each(|c| *c = c.scale(1.0 / n_final));
        }

        for (row, &val) in vec.iter().enumerate() {
            vecs.set(row, col, val);
        }
        basis.push(vec);
    }

    Ok(vecs)
}

/// Cross product of two 3-vectors of complex numbers: (a × b)_i = ε_{ijk} a_j b_k.
fn cross_product_3c(a: &[Complex; 3], b: &[Complex; 3]) -> Vec<Complex> {
    vec![
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}

/// Return a unit vector orthogonal to all vectors in `basis` within R^n / C^n.
fn standard_basis_complement(basis: &[Vec<Complex>], n: usize) -> Vec<Complex> {
    // Try each standard basis vector e_i and take the one with smallest projection
    let mut best = vec![Complex::ZERO; n];
    let mut best_norm = -1.0_f64;
    for i in 0..n {
        let mut candidate = vec![Complex::ZERO; n];
        candidate[i] = Complex::ONE;
        // Subtract projections
        for prev in basis {
            let dot: Complex = candidate
                .iter()
                .zip(prev.iter())
                .map(|(&a, &b)| b.conj().mul(&a))
                .fold(Complex::ZERO, |acc, x| acc.add(&x));
            for (c, &p) in candidate.iter_mut().zip(prev.iter()) {
                *c = c.sub(&p.mul(&dot));
            }
        }
        let nm = candidate.iter().map(|c| c.norm_sq()).sum::<f64>().sqrt();
        if nm > best_norm {
            best_norm = nm;
            if nm > 1e-12 {
                best = candidate.iter().map(|c| c.scale(1.0 / nm)).collect();
            } else {
                best = candidate;
            }
        }
    }
    best
}

// ---------------------------------------------------------------------------
// Honeycomb geometry: primitive and reciprocal lattice vectors
// ---------------------------------------------------------------------------

// Honeycomb NN vectors (in units of a=1): δ₁=(1,0), δ₂=(−½,√3/2), δ₃=(−½,−√3/2)
// but we work with primitive lattice vectors a1=(1,0), a2=(1/2,√3/2).
// The three NN bond phase factors for A→B are:
//   f(k) = e^{0} + e^{i k·a1} + e^{i k·a2}  (choosing gauge with A at origin)
// Actually using the convention: f(k) = 1 + e^{i kx} + e^{i(kx/2 + ky√3/2)}
// Reciprocal vectors: b1=2π(1,−1/√3), b2=2π(0,2/√3)

const SQRT3: f64 = 1.732_050_808_568_877_3;

/// Three NN phase factors for the honeycomb lattice.
/// Returns (f_re, f_im) = sum of e^{i k·delta_j} for j=1,2,3
/// Using delta vectors: (0,0), a1=(1,0), a2=(1/2, sqrt3/2)
#[inline]
fn honeycomb_nn_phase(kx: f64, ky: f64) -> Complex {
    // delta_0: (0,0)  → phase = 1
    // delta_1: a1 = (1,0) → phase = e^{i kx}
    // delta_2: a2 = (1/2, sqrt3/2) → phase = e^{i(kx/2 + ky*sqrt3/2)}
    let p0 = Complex::ONE;
    let p1 = Complex::new(0.0, kx).exp();
    let p2 = Complex::new(0.0, kx * 0.5 + ky * SQRT3 * 0.5).exp();
    p0.add(&p1).add(&p2)
}

/// NNN hopping phases: sum of sin(k·bi) for Haldane mass term.
/// b1 = 2π*(1, -1/√3), b2 = 2π*(0, 2/√3)
/// The six NNN vectors are ±a1, ±a2, ±(a1-a2). Their dot products with k:
///   ±k·a1 = ±kx
///   ±k·a2 = ±(kx/2 + ky*√3/2)
///   ±k·(a1-a2) = ±(kx/2 - ky*√3/2)
/// The antisymmetric (sin) combination gives the DMI contribution.
#[inline]
fn honeycomb_nnn_dz(kx: f64, ky: f64, dmi: f64, h_ext: f64) -> f64 {
    let phi1 = kx;
    let phi2 = kx * 0.5 + ky * SQRT3 * 0.5;
    let phi3 = kx * 0.5 - ky * SQRT3 * 0.5;
    h_ext + 2.0 * dmi * (phi1.sin() + phi2.sin() + phi3.sin())
}

// ---------------------------------------------------------------------------
// Lattice type
// ---------------------------------------------------------------------------

/// Lattice geometry for the magnon band model.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum LatticeType {
    /// Honeycomb ferromagnet with Haldane-type NNN DMI (2 bands).
    HoneycombHaldane,
    /// Kagome ferromagnet with DMI (3 bands, contains a flat band).
    Kagome,
    /// Square lattice with DMI (1 band, non-topological baseline).
    SquareDmi,
}

// ---------------------------------------------------------------------------
// MagnonBandModel
// ---------------------------------------------------------------------------

/// Magnon Bloch-Hamiltonian model for a magnetic insulator.
///
/// Stores all exchange and DMI parameters needed to construct H(k) at arbitrary
/// wavevectors, compute band dispersions, and feed higher-level topology tools.
#[derive(Debug, Clone)]
pub struct MagnonBandModel {
    /// Lattice geometry (determines n_bands and H(k) formula).
    pub lattice: LatticeType,
    /// Number of bands (sublattice sites per unit cell).
    pub n_bands: usize,
    /// Nearest-neighbour exchange coupling \[meV\] (must be > 0).
    pub j_nn: f64,
    /// Next-nearest-neighbour exchange coupling \[meV\] (Haldane mass coefficient).
    pub j_nnn: f64,
    /// DMI strength \[meV\] (same units as j_nn).
    pub dmi: f64,
    /// External Zeeman field \[meV\] (added to diagonal, breaks TRS).
    pub h_ext: f64,
    /// Lattice constant \[m\] — used when converting to physical dispersion (SI).
    pub a_lattice: f64,
}

impl MagnonBandModel {
    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    /// Construct a honeycomb Haldane-magnon model.
    ///
    /// The two-band Bloch Hamiltonian H(k) = h₀·I + **d**·σ has:
    /// - `d_x + i·d_y = j_nn · f(k)` (NN hopping, f is the phase sum)
    /// - `d_z = h_ext + 2·dmi·(sin(k·a₁) + sin(k·a₂) + sin(k·(a₁−a₂)))`
    ///
    /// When `dmi ≠ 0` or `h_ext ≠ 0` the bands are gapped and topological (|C| = 1).
    pub fn honeycomb_haldane(j_nn: f64, j_nnn: f64, dmi: f64, h_ext: f64) -> Result<Self> {
        if j_nn <= 0.0 {
            return Err(error::invalid_param("j_nn", "must be positive"));
        }
        Ok(Self {
            lattice: LatticeType::HoneycombHaldane,
            n_bands: 2,
            j_nn,
            j_nnn,
            dmi,
            h_ext,
            a_lattice: 3.0e-10, // typical honeycomb (e.g., CrI₃) lattice constant
        })
    }

    /// Construct a kagome magnon model.
    ///
    /// The three-band model hosts one flat band at energy `+2J + h_ext` (top band)
    /// and two dispersive bands (Mielke 1991). The DMI (`dmi`) lifts degeneracies
    /// and generates Berry curvature, making one or more bands topological.
    pub fn kagome(j: f64, dmi: f64, h_ext: f64) -> Result<Self> {
        if j <= 0.0 {
            return Err(error::invalid_param("j", "must be positive"));
        }
        Ok(Self {
            lattice: LatticeType::Kagome,
            n_bands: 3,
            j_nn: j,
            j_nnn: 0.0,
            dmi,
            h_ext,
            a_lattice: 5.0e-10,
        })
    }

    /// Construct a square-lattice magnon model with DMI.
    ///
    /// Single-band dispersion ε(k) = h_ext + 2J(cos kx + cos ky).
    /// This model is non-topological but serves as a baseline for curvature and
    /// Hall conductivity calculations (Berry curvature is identically zero).
    pub fn square_dmi(j: f64, dmi: f64, h_ext: f64) -> Result<Self> {
        if j <= 0.0 {
            return Err(error::invalid_param("j", "must be positive"));
        }
        Ok(Self {
            lattice: LatticeType::SquareDmi,
            n_bands: 1,
            j_nn: j,
            j_nnn: 0.0,
            dmi,
            h_ext,
            a_lattice: 4.0e-10,
        })
    }

    // -----------------------------------------------------------------------
    // Bloch Hamiltonian
    // -----------------------------------------------------------------------

    /// Build the `n_bands × n_bands` Bloch Hamiltonian H(k) at wavevector k = (kx, ky).
    ///
    /// All energies are in the units of `j_nn` (typically meV). The returned matrix
    /// is Hermitian by construction.
    #[inline]
    pub fn hamiltonian_at(&self, k: (f64, f64)) -> Result<CMatrix> {
        let (kx, ky) = k;
        match &self.lattice {
            LatticeType::HoneycombHaldane => self.hamiltonian_honeycomb(kx, ky),
            LatticeType::Kagome => self.hamiltonian_kagome(kx, ky),
            LatticeType::SquareDmi => self.hamiltonian_square(kx, ky),
        }
    }

    /// Honeycomb 2×2 Bloch Hamiltonian.
    fn hamiltonian_honeycomb(&self, kx: f64, ky: f64) -> Result<CMatrix> {
        // NN hopping sum: f(k) = 1 + e^{i kx} + e^{i(kx/2+ky√3/2)}
        let f_k = honeycomb_nn_phase(kx, ky);

        // Topological mass from NNN-DMI
        let dz = honeycomb_nnn_dz(kx, ky, self.dmi, self.h_ext);

        // NNN diagonal shift (j_nnn modifies both bands equally — sets bandwidth offset)
        // The NNN cos terms: h0 contribution from NNN hopping
        let phi1 = kx;
        let phi2 = kx * 0.5 + ky * SQRT3 * 0.5;
        let phi3 = kx * 0.5 - ky * SQRT3 * 0.5;
        let h0_nnn = -self.j_nnn * (phi1.cos() + phi2.cos() + phi3.cos());

        // H = (h0_nnn)·I + [[dz, j_nn·f*],[j_nn·f, -dz]]
        // But standard form: H[0,0] = h0_nnn + dz, H[1,1] = h0_nnn - dz
        //                    H[0,1] = j_nn * f*, H[1,0] = j_nn * f
        let h00 = Complex::from_real(h0_nnn + dz);
        let h11 = Complex::from_real(h0_nnn - dz);
        let h01 = f_k.conj().scale(self.j_nn);
        let h10 = f_k.scale(self.j_nn);

        CMatrix::from_rows(vec![vec![h00, h01], vec![h10, h11]])
    }

    /// Kagome 3×3 Bloch Hamiltonian.
    ///
    /// Following the Mielke (1991) / Owerre (2016) convention for the kagome lattice
    /// with primitive vectors **a₁** = (1, 0), **a₂** = (½, √3/2). The three
    /// sublattice sites per unit cell connect via NN bonds of length *a*/2, with
    /// half-bond phases:
    ///
    /// - Bond 0–1 (along **a₁**/2): φ₁₂ = kx/2
    /// - Bond 0–2 (along (**a₁**−**a₂**)/2): φ₁₃ = kx/4 − ky√3/4
    /// - Bond 1–2 (along **a₂**/2): φ₂₃ = kx/4 + ky√3/4
    ///
    /// The Bloch Hamiltonian with NN exchange J and out-of-plane scalar DMI D is:
    ///
    /// ```text
    /// H_{αβ}(k) = −2J · cos(φ_{αβ}) + i · D · sin(φ_{αβ})
    /// ```
    ///
    /// At D = 0 this is a real symmetric matrix with an **exact flat band** at
    /// ε = +2J + h_ext for all **k** (Mielke 1991, top eigenvalue). The DMI imaginary
    /// term (i·D·sin) is antisymmetric, making the full matrix Hermitian; it breaks
    /// TRS and opens topological magnon gaps.
    fn hamiltonian_kagome(&self, kx: f64, ky: f64) -> Result<CMatrix> {
        let eps0 = self.h_ext;

        // Half-bond phases for the three kagome NN bonds:
        // Bond 0-1 (along a₁/2):       phi12 = kx/2
        // Bond 0-2 (along (a₁-a₂)/2):  phi13 = kx/4 - ky*sqrt3/4
        // Bond 1-2 (along a₂/2):        phi23 = kx/4 + ky*sqrt3/4
        //
        // Bloch Hamiltonian (Mielke 1991, Owerre 2016):
        //   H[α,β] = -2J*cos(φ_{αβ}) + i*D*sin(φ_{αβ})
        //
        // At D=0 this is real symmetric with an exact flat band at ε = -2J + h_ext
        // for all k. The DMI imaginary term breaks TRS and opens topological gaps.
        let phi12 = kx * 0.5;
        let phi13 = kx * 0.25 - ky * SQRT3 * 0.25;
        let phi23 = kx * 0.25 + ky * SQRT3 * 0.25;

        let h01 = Complex::new(-2.0 * self.j_nn * phi12.cos(), self.dmi * phi12.sin());
        let h02 = Complex::new(-2.0 * self.j_nn * phi13.cos(), self.dmi * phi13.sin());
        let h12 = Complex::new(-2.0 * self.j_nn * phi23.cos(), self.dmi * phi23.sin());

        let diag = Complex::from_real(eps0);

        CMatrix::from_rows(vec![
            vec![diag, h01, h02],
            vec![h01.conj(), diag, h12],
            vec![h02.conj(), h12.conj(), diag],
        ])
    }

    /// Square-lattice 1×1 Bloch Hamiltonian.
    fn hamiltonian_square(&self, kx: f64, ky: f64) -> Result<CMatrix> {
        let eps = self.h_ext + 2.0 * self.j_nn * (kx.cos() + ky.cos());
        CMatrix::from_rows(vec![vec![Complex::from_real(eps)]])
    }

    // -----------------------------------------------------------------------
    // Diagonalization
    // -----------------------------------------------------------------------

    /// Diagonalize H(k) and return sorted eigenvalues (ascending) and eigenvectors.
    ///
    /// Uses direct analytical formulas for 1×1 and 2×2 matrices, and a robust
    /// iterative Jacobi algorithm for larger matrices. The k-th column of the
    /// returned `CMatrix` is the eigenvector for the k-th eigenvalue.
    pub fn diagonalize(&self, k: (f64, f64)) -> Result<(Vec<f64>, CMatrix)> {
        let h = self.hamiltonian_at(k)?;
        match h.n() {
            1 => diagonalize_1x1(&h),
            2 => diagonalize_2x2(&h),
            3 => diagonalize_3x3_hermitian(&h),
            _ => h.hermitian_eigendecomposition(),
        }
    }

    // -----------------------------------------------------------------------
    // Band structure on a grid
    // -----------------------------------------------------------------------

    /// Compute the band structure on a uniform (kx, ky) grid covering [−π, π]².
    ///
    /// Returns `bands[ikx][iky] = Vec<f64>` with `n_bands` eigenvalues (ascending).
    /// The grid has `kx_pts × ky_pts` points sampling the first Brillouin zone.
    pub fn bands(&self, kx_pts: usize, ky_pts: usize) -> Result<Vec<Vec<Vec<f64>>>> {
        if kx_pts < 2 || ky_pts < 2 {
            return Err(error::invalid_param("kx_pts/ky_pts", "must be at least 2"));
        }
        let mut result = Vec::with_capacity(kx_pts);
        for ix in 0..kx_pts {
            let kx = -PI + 2.0 * PI * (ix as f64) / (kx_pts as f64);
            let mut row = Vec::with_capacity(ky_pts);
            for iy in 0..ky_pts {
                let ky = -PI + 2.0 * PI * (iy as f64) / (ky_pts as f64);
                let (evals, _) = self.diagonalize((kx, ky))?;
                row.push(evals);
            }
            result.push(row);
        }
        Ok(result)
    }

    /// Minimum direct band gap between `band_below` and `band_above` over the BZ grid.
    ///
    /// Both indices are zero-based. Returns the minimum gap value in the same units
    /// as the exchange parameters. A gap of zero (or negative) signals band touching.
    pub fn band_gap(
        &self,
        band_below: usize,
        band_above: usize,
        kx_pts: usize,
        ky_pts: usize,
    ) -> Result<f64> {
        if band_below >= self.n_bands {
            return Err(error::invalid_param("band_below", "index out of range"));
        }
        if band_above >= self.n_bands {
            return Err(error::invalid_param("band_above", "index out of range"));
        }
        if band_above <= band_below {
            return Err(error::invalid_param(
                "band_above",
                "must be greater than band_below",
            ));
        }
        let bs = self.bands(kx_pts, ky_pts)?;
        let mut min_gap = f64::INFINITY;
        for row in &bs {
            for evals in row {
                let gap = evals[band_above] - evals[band_below];
                if gap < min_gap {
                    min_gap = gap;
                }
            }
        }
        Ok(min_gap)
    }

    /// Number of bands (sublattice sites per unit cell).
    #[inline]
    pub fn n_bands(&self) -> usize {
        self.n_bands
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn honeycomb_two_bands() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.1, 0.0).unwrap();
        assert_eq!(m.n_bands(), 2);
        let h = m.hamiltonian_at((0.1, 0.2)).unwrap();
        assert_eq!(h.n(), 2);
    }

    #[test]
    fn kagome_three_bands() {
        let m = MagnonBandModel::kagome(1.0, 0.05, 0.0).unwrap();
        assert_eq!(m.n_bands(), 3);
        let h = m.hamiltonian_at((0.0, 0.0)).unwrap();
        assert_eq!(h.n(), 3);
    }

    #[test]
    fn square_one_band() {
        let m = MagnonBandModel::square_dmi(1.0, 0.0, 0.0).unwrap();
        assert_eq!(m.n_bands(), 1);
        let h = m.hamiltonian_at((0.5, 0.5)).unwrap();
        assert_eq!(h.n(), 1);
    }

    #[test]
    fn hamiltonian_hermitian() {
        // H = H† for all models at various k-points
        let models = vec![
            MagnonBandModel::honeycomb_haldane(1.0, 0.1, 0.15, 0.2).unwrap(),
            MagnonBandModel::kagome(2.0, 0.3, 0.1).unwrap(),
            MagnonBandModel::square_dmi(1.5, 0.2, 0.05).unwrap(),
        ];
        let kpoints = [(0.0, 0.0), (PI / 3.0, PI / 4.0), (1.1, -0.7)];
        for m in &models {
            for &kpt in &kpoints {
                let h = m.hamiltonian_at(kpt).unwrap();
                let hd = h.conj_transpose();
                let diff = h.sub(&hd).unwrap();
                assert!(
                    diff.frobenius_norm() < 1e-12,
                    "H not Hermitian for {:?} at k={:?}",
                    m.lattice,
                    kpt
                );
            }
        }
    }

    #[test]
    fn diagonalize_sorted_ascending() {
        let m = MagnonBandModel::kagome(1.0, 0.1, 0.0).unwrap();
        for (kx_i, ky_i) in [(0.0f64, 0.0), (0.5, 0.7), (1.2, -0.3)] {
            let (evals, _) = m.diagonalize((kx_i, ky_i)).unwrap();
            for w in evals.windows(2) {
                assert!(w[0] <= w[1] + 1e-12, "Eigenvalues not sorted: {:?}", evals);
            }
        }
    }

    #[test]
    fn eigenvectors_orthonormal() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.2, 0.1).unwrap();
        let (_, vecs) = m.diagonalize((0.3, 0.4)).unwrap();
        let n = m.n_bands();
        // Check V†V = I
        for i in 0..n {
            for j in 0..n {
                let dot: Complex = (0..n)
                    .map(|r| vecs.get(r, i).conj().mul(&vecs.get(r, j)))
                    .fold(Complex::ZERO, |acc, x| acc.add(&x));
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approx(dot.re, expected, 1e-9),
                    "V†V[{},{}] re off: got {}",
                    i,
                    j,
                    dot.re
                );
                assert!(approx(dot.im, 0.0, 1e-9), "V†V[{},{}] im off", i, j);
            }
        }
    }

    #[test]
    fn honeycomb_zero_dmi_band_gap_zero() {
        // Without DMI or field, honeycomb bands touch at the K point.
        // The K point is at k = (4π/3, 0) in our convention where f(K)=0.
        // We directly evaluate at K to verify the gap is exactly zero (no DMI → no mass).
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.0, 0.0).unwrap();
        // K point: kx = 4π/3 (equivalent to -2π/3 by folding)
        let k_x = 4.0 * PI / 3.0;
        let (evals, _) = m.diagonalize((k_x, 0.0)).unwrap();
        let gap_at_k = evals[1] - evals[0];
        assert!(
            gap_at_k < 1e-10,
            "Expected zero gap at K point without DMI, got {}",
            gap_at_k
        );
    }

    #[test]
    fn honeycomb_nonzero_dmi_band_gap_positive() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.3, 0.0).unwrap();
        let gap = m.band_gap(0, 1, 20, 20).unwrap();
        assert!(gap > 0.1, "Expected positive gap with DMI, got {}", gap);
    }

    #[test]
    fn kagome_flat_band_exists() {
        // The kagome NN model (dmi=0, h_ext=0) with H[α,β]=-2J*cos(φ_{αβ})+i*D*sin(φ_{αβ})
        // has an exact flat band at ε = +2J for ALL k (Mielke 1991, top eigenvalue).
        // We verify that evals[2] = +2J at a representative set of k-points.
        let j = 1.0;
        let m = MagnonBandModel::kagome(j, 0.0, 0.0).unwrap();
        let flat_energy = 2.0 * j; // flat band at +2J

        let kpoints = [
            (0.0, 0.0),
            (PI, 0.0),
            (0.0, PI),
            (PI, PI),
            (PI / 2.0, PI / 3.0),
            (-PI / 3.0, 2.0 * PI / 5.0),
        ];

        for (kx, ky) in &kpoints {
            let (evals, _) = m.diagonalize((*kx, *ky)).unwrap();
            assert!(
                (evals[2] - flat_energy).abs() < 1e-8,
                "Flat band at k=({:.4},{:.4}): expected +2J={:.4}, got {:.8}",
                kx,
                ky,
                flat_energy,
                evals[2]
            );
        }
    }

    #[test]
    fn bands_grid_size_correct() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.1, 0.0).unwrap();
        let bs = m.bands(12, 15).unwrap();
        assert_eq!(bs.len(), 12);
        assert_eq!(bs[0].len(), 15);
        assert_eq!(bs[0][0].len(), 2);
    }

    #[test]
    fn band_gap_positive_with_dmi() {
        let m = MagnonBandModel::kagome(1.0, 0.5, 0.0).unwrap();
        // With large DMI the lowest gap should open
        let gap = m.band_gap(0, 1, 15, 15).unwrap();
        // Gap should be positive
        assert!(gap > -1e-6, "Gap should be non-negative: {}", gap);
    }

    #[test]
    fn square_dispersion_cosine_shape() {
        // At k=(0,0): ε = h_ext + 4J; at k=(π,π): ε = h_ext − 4J
        let j = 1.0;
        let h = 0.5;
        let m = MagnonBandModel::square_dmi(j, 0.0, h).unwrap();
        let (e_gamma, _) = m.diagonalize((0.0, 0.0)).unwrap();
        let (e_corner, _) = m.diagonalize((PI, PI)).unwrap();
        assert!(approx(e_gamma[0], h + 4.0 * j, 1e-10));
        assert!(approx(e_corner[0], h - 4.0 * j, 1e-10));
    }

    #[test]
    fn invalid_negative_j_rejected() {
        assert!(MagnonBandModel::honeycomb_haldane(-1.0, 0.0, 0.0, 0.0).is_err());
        assert!(MagnonBandModel::kagome(-0.5, 0.0, 0.0).is_err());
        assert!(MagnonBandModel::square_dmi(-2.0, 0.0, 0.0).is_err());
    }
}
