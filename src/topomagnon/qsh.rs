//! Kane-Mele model for the Quantum Spin Hall (QSH) effect.
//!
//! Implements the original Kane-Mele Hamiltonian for a 2D topological insulator
//! on the honeycomb lattice, as introduced in:
//!
//! - C. L. Kane and E. J. Mele, "Quantum Spin Hall Effect in Graphene",
//!   *Phys. Rev. Lett.* **95**, 226801 (2005).
//! - C. L. Kane and E. J. Mele, "Z₂ Topological Order and the Quantum Spin Hall Effect",
//!   *Phys. Rev. Lett.* **95**, 146802 (2005).
//!
//! # Hamiltonian
//!
//! The full four-band Bloch Hamiltonian in the (A↑, B↑, A↓, B↓) basis is:
//!
//! ```text
//! H(k) = -t  Σ_{<ij>} c†_i c_j                           (NN hopping)
//!       + iλ_SO Σ_{<<ij>>} ν_ij c†_i σ_z c_j             (intrinsic SOC, NNN)
//!       + iλ_R  Σ_{<ij>} c†_i (σ × d̂_ij)_z c_j          (Rashba SOC, NN)
//!       + λ_v Σ_i ξ_i c†_i c_i                            (staggered potential)
//! ```
//!
//! ## 4×4 Block Structure
//!
//! ```text
//! H(k) = | h_uu  h_ud |
//!         | h_du  h_dd |
//! ```
//!
//! where:
//! - `h_uu` (spin-up block, 2×2): NN hopping + intrinsic SOC + staggered potential
//! - `h_dd` (spin-down block, 2×2): time-reversal partner of h_uu
//! - `h_ud`, `h_du`: Rashba spin-orbit coupling connecting up↔down spin
//!
//! # Topological Phase
//!
//! With λ_R = 0 (no Rashba), the model decouples into spin-up and spin-down
//! Haldane models. The QSH phase (Z₂ = 1) occurs when λ_v < 3√3 λ_SO.
//! The topological invariant is the Z₂ index, computed from the spin-up
//! Chern number (C_up) via Z₂ = C_up mod 2 (for λ_R = 0).
//!
//! # Module Contents
//!
//! | Type | Description |
//! |---|---|
//! | [`KaneMeleModel`] | Full Kane-Mele Bloch Hamiltonian and topological analysis |

use std::f64::consts::PI;

use crate::error::{self, Result};
use crate::math::{CMatrix, Complex};

// ---------------------------------------------------------------------------
// Honeycomb geometry constants and helpers
// ---------------------------------------------------------------------------

const SQRT3: f64 = 1.732_050_808_568_877_3;

/// Three NN bond phase factors f(k) for the A→B hoppings on the honeycomb.
///
/// Using primitive vectors a₁ = (1, 0), a₂ = (½, √3/2) in units of lattice
/// constant a=1. The three A→B bond vectors are:
///   δ₀ = (0, 0) relative to A on the same sublattice origin  → but actually
///   we use the A-B offset: gauge choice gives f(k) = 1 + e^{ik·a₁} + e^{ik·a₂}.
///
/// This is the standard Kane-Mele convention where A is at origin and B neighbours
/// are connected via 0, a₁, a₂.
#[inline]
fn nn_phase_factor(kx: f64, ky: f64) -> Complex {
    // f(k) = e^{i·0} + e^{i k·a₁} + e^{i k·a₂}
    //      = 1 + e^{i·kx} + e^{i(kx/2 + ky·√3/2)}
    let p0 = Complex::ONE;
    let p1 = Complex::new(0.0, kx).exp();
    let p2 = Complex::new(0.0, kx * 0.5 + ky * SQRT3 * 0.5).exp();
    p0.add(&p1).add(&p2)
}

/// NNN hopping structure factor g(k) for intrinsic SOC.
///
/// The six NNN vectors are ±a₁, ±a₂, ±(a₁ − a₂). With the clockwise/counterclockwise
/// convention for ν_ij (+1 clockwise, -1 counterclockwise), the antisymmetric
/// combination for A-sublattice (positive ν) gives:
///
///   g(k) = 2[sin(k·a₁) - sin(k·a₂) + sin(k·(a₂ − a₁))]
///
/// This equals the standard Kane-Mele NNN structure factor. For B-sublattice
/// the sign of ν_ij reverses, giving -g(k).
#[inline]
fn nnn_soc_factor(kx: f64, ky: f64) -> f64 {
    // k·a₁ = kx,  k·a₂ = kx/2 + ky·√3/2,  k·(a₁-a₂) = kx/2 - ky·√3/2
    // Kane-Mele NNN structure factor (Eq. 4 of PRL 95, 226801):
    //   g(k) = 2[sin(k·a₁) - sin(k·a₂) + sin(k·(a₂-a₁))]
    let phi1 = kx;
    let phi2 = kx * 0.5 + ky * SQRT3 * 0.5;
    2.0 * (phi1.sin() - phi2.sin() + (phi2 - phi1).sin())
}

/// Rashba coupling vector components from NN bond directions.
///
/// For each NN bond δ = (δx, δy)/|δ|, the Rashba term contributes
/// e^{ik·δ}·(δ_y − i·δ_x)/|δ| to the A→B off-diagonal Rashba matrix element.
///
/// The three NN bond vectors (A→B direction, normalized):
///   δ₀ = (0, 1/√3) → normalized (0, 1)
///   δ₁ = (½, -1/(2√3)) = a₁ - δ₀ direction → normalized (√3/2, -1/2)
///   δ₂ = (-½, -1/(2√3)) = a₂ - δ₀ direction → normalized (-√3/2, -1/2)
///
/// Actually the NN vectors (in the standard honeycomb basis with A at origin and
/// primitive vectors a₁, a₂) are the same as the bond phase factors. We parameterize
/// by the normalized NN bond unit vectors d̂_δ = δ/|δ|.
///
/// Returns the complex scalar R_AB = iλ_R · Σ_δ e^{ik·δ}·(d̂_δ_y - i·d̂_δ_x)
/// which enters the off-diagonal spin block h_ud\[0,1\] = -R_AB, h_ud[1,0] = R_AB*.
#[inline]
fn rashba_ab_coupling(kx: f64, ky: f64, lambda_r: f64) -> Complex {
    if lambda_r.abs() < 1e-15 {
        return Complex::ZERO;
    }
    // NN bond vectors in Cartesian (A at origin, with a lattice const = 1):
    // The B neighbours are at: τ₀ = (0, 1/√3), τ₁ = (½, -1/(2√3)), τ₂ = (-½, -1/(2√3))
    // But the phase factors use the lattice vector convention: bond 0 → (0,0) [intra],
    // bond 1 → a₁ = (1, 0), bond 2 → a₂ = (½, √3/2).
    //
    // The physical bond vectors (A→B) in Cartesian, with A at (0,0) and B at:
    //   B_0: offset (0, 1/√3) (straight up)
    //   B_1: B_0 - a₁ = (-1, 1/√3) ... no, let's use the standard convention:
    //
    // With a₁=(1,0), a₂=(½,√3/2), place A at (0,0) and B at τ = (0, 1/√3).
    // Then B neighbours of A are at: τ (δ₀), τ+a₁ = (1, 1/√3) (δ₁), τ+a₂ = (½, 1/√3+√3/2) ... wrong.
    //
    // Use the simplest consistent convention matching the NN phase factor f(k):
    // f(k) = 1 + e^{ik·a₁} + e^{ik·a₂}
    // Bond 0: lattice shift = (0,0)  → physical direction d̂₀
    // Bond 1: lattice shift = a₁    → physical direction d̂₁
    // Bond 2: lattice shift = a₂    → physical direction d̂₂
    //
    // In the standard honeycomb with A and B in the unit cell, the bond vectors
    // pointing from A to B (in units of a=1) are:
    //   d₀ = (0, 1/√3)                 → |d₀|=1/√3, d̂₀ = (0, 1)
    //   d₁ = a₁ + d₀ - a₁ = wait...
    //
    // More carefully: if A is at (0,0) and the three B neighbours are reached by:
    //   d₀ = (0, 1/√3)                          → phase e^{0}  = 1
    //   d₁ = (1/2, -1/(2√3)) = a₁ - (d₂+d₀)/2  → actually use Bravais convention
    //
    // Standard Kane-Mele: B atoms at positions τ₁=(0,a/√3), τ₂=(-a/2,-a/(2√3)), τ₃=(a/2,-a/(2√3))
    // relative to A. With a=1:
    //   δ₁ = (0,  1/√3),     |δ|=1/√3, d̂ = (0,     1)
    //   δ₂ = (-½, -1/(2√3)), |δ|=1/√3, d̂ = (-√3/2, -½)
    //   δ₃ = ( ½, -1/(2√3)), |δ|=1/√3, d̂ = ( √3/2, -½)
    //
    // Phase factors: e^{ik·δ_j} where k is in reciprocal lattice units.
    // In our convention: phase for bond j is e^{i k·δ_j / a} where a=a_lattice=1.
    //
    // IMPORTANT: The nn_phase_factor function uses f(k)=1+e^{ikx}+e^{i(kx/2+ky√3/2)}.
    // This corresponds to bond vectors: (0,0), a₁=(1,0), a₂=(½,√3/2).
    // These are LATTICE vectors, not the physical bond vectors (A→B).
    // The physical bond vectors depend on where B is within the unit cell.
    //
    // For the Rashba term in the Kane-Mele paper, what matters is the unit vector
    // d̂ along each A-B bond. Using the standard B sublattice position (0,1/√3)
    // relative to A (within unit cell), the three physical A→B bond vectors are:
    //   δ_phys₁ = (0, 1/√3)               [bond within same unit cell]
    //   δ_phys₂ = a₁ + (0,-1/√3) = (1, -1/√3)  ... this doesn't match standard result
    //
    // To avoid inconsistency, use the result derived in the literature: the Rashba
    // coupling R(k) = iλ_R · (Δ_x - i·Δ_y) where:
    //   Δ_x = Σ_δ e^{ik·δ} d̂_x = sum of x-component of bond unit vectors × phase
    //   Δ_y = Σ_δ e^{ik·δ} d̂_y = sum of y-component of bond unit vectors × phase
    //
    // The three normalized bond unit vectors (standard honeycomb, equidistant):
    //   d̂₁ = (0, 1)          with lattice shift (0,0)   → phase = 1
    //   d̂₂ = (-√3/2, -½)    with lattice shift a₁       → phase = e^{ikx}
    //   d̂₃ = ( √3/2, -½)    with lattice shift a₂       → phase = e^{i(kx/2+ky√3/2)}
    //
    // Rashba matrix element (A↑→B↓): h_ud[0,1]
    //   = iλ_R · Σ_j e^{ik·Δ_j} (ẑ × d̂_j)·σ_{↑↓}
    //   = iλ_R · Σ_j e^{ik·Δ_j} (d̂_j_x σ_y - d̂_j_y σ_x)_{↑↓}
    //   For σ_x: (↑↓) element = 1; for σ_y: (↑↓) element = -i
    //   (ẑ × d̂)·σ for ↑→↓: (d_x·σ_y - d_y·σ_x)_{↑↓} = d_x·(-i) - d_y·1 = -i·d_x - d_y
    //
    // So: h_{A↑,B↓} = iλ_R · Σ_j e^{ik·Δ_j} (-i·d̂_j_x - d̂_j_y)
    //               = λ_R · Σ_j e^{ik·Δ_j} (d̂_j_x - i·d̂_j_y) · (-1)
    //               Actually: i·(-i·d_x - d_y) = d_x - i·d_y
    //
    // h_{A↑,B↓} = iλ_R · Σ_j e^{ik·Δ_j} (-i·d_x - d_y)
    //           = λ_R · Σ_j e^{ik·Δ_j} (d_x - i·d_y) · ... let me redo carefully:
    //
    // With i·(-i·d_x - d_y) = i·(-i)·d_x + i·(-1)·d_y = d_x - i·d_y ✓
    //
    // R_AB ≡ h_{A↑,B↓} = λ_R · Σ_j e^{ik·Δ_j} (d̂_j_x - i·d̂_j_y)

    let p1 = Complex::ONE; // bond 1: lattice shift (0,0), d̂₁=(0,1)
    let p2 = Complex::new(0.0, kx).exp(); // bond 2: shift a₁, d̂₂=(-√3/2,-½)
    let p3 = Complex::new(0.0, kx * 0.5 + ky * SQRT3 * 0.5).exp(); // bond 3: shift a₂, d̂₃=(√3/2,-½)

    // (d̂_j_x - i·d̂_j_y) for each bond:
    // Bond 1: d̂=(0,1)        → (0 - i·1) = -i
    // Bond 2: d̂=(-√3/2,-½)  → (-√3/2 - i·(-½)) = -√3/2 + i/2
    // Bond 3: d̂=(√3/2,-½)   → (√3/2 - i·(-½)) = √3/2 + i/2
    let chi1 = Complex::new(0.0, -1.0); // -i
    let chi2 = Complex::new(-SQRT3 * 0.5, 0.5);
    let chi3 = Complex::new(SQRT3 * 0.5, 0.5);

    let sum = p1.mul(&chi1).add(&p2.mul(&chi2)).add(&p3.mul(&chi3));
    sum.scale(lambda_r)
}

// ---------------------------------------------------------------------------
// KaneMeleModel
// ---------------------------------------------------------------------------

/// Kane-Mele model for the Quantum Spin Hall effect on the honeycomb lattice.
///
/// The four-band Bloch Hamiltonian H(k) is a 4×4 matrix in the
/// (A↑, B↑, A↓, B↓) basis. The four physics parameters are:
///
/// - `t`: nearest-neighbour hopping (eV). Gives graphene-like linear Dirac cones.
/// - `lambda_so`: intrinsic spin-orbit coupling (eV), NNN imaginary hopping.
///   Opens a gap Δ ≈ 6√3 λ_SO at the Dirac points.
/// - `lambda_r`: Rashba SOC (eV), NN spin-flipping hopping due to substrate
///   or electric field. Breaks the conservation of S_z.
/// - `lambda_v`: staggered sublattice potential (eV), breaks inversion symmetry.
///   When λ_v > 3√3 λ_SO the gap closes and the system becomes trivial.
///
/// # Topological Phase Diagram (λ_R = 0)
///
/// ```text
/// |λ_v| < 3√3 λ_SO  →  QSH phase (Z₂ = 1, helical edge states)
/// |λ_v| = 3√3 λ_SO  →  Critical point (gapless at K, K')
/// |λ_v| > 3√3 λ_SO  →  Trivial insulator (Z₂ = 0)
/// ```
///
/// The critical value is 3√3 ≈ 5.196; so for λ_SO = 0.05, the boundary is at
/// λ_v ≈ 0.26 eV.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct KaneMeleModel {
    /// Nearest-neighbour hopping amplitude \[eV\]. Must be > 0.
    pub t: f64,
    /// Intrinsic (Kane-Mele) spin-orbit coupling \[eV\]. NNN imaginary hopping.
    pub lambda_so: f64,
    /// Rashba spin-orbit coupling \[eV\]. NN spin-flip hopping.
    pub lambda_r: f64,
    /// Staggered sublattice potential \[eV\]. Breaks inversion symmetry.
    pub lambda_v: f64,
    /// Honeycomb lattice constant \[m\]. Default: 2.46 × 10⁻¹⁰ m (graphene).
    pub a_lattice: f64,
}

impl KaneMeleModel {
    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    /// Build a Kane-Mele model with explicit parameters.
    ///
    /// # Errors
    ///
    /// - `InvalidParameter` if `t ≤ 0` or any parameter is non-finite.
    pub fn new(t: f64, lambda_so: f64, lambda_r: f64, lambda_v: f64) -> Result<Self> {
        if t <= 0.0 {
            return Err(error::invalid_param("t", "NN hopping must be positive"));
        }
        if !t.is_finite() {
            return Err(error::invalid_param("t", "must be finite"));
        }
        if !lambda_so.is_finite() {
            return Err(error::invalid_param("lambda_so", "must be finite"));
        }
        if !lambda_r.is_finite() {
            return Err(error::invalid_param("lambda_r", "must be finite"));
        }
        if !lambda_v.is_finite() {
            return Err(error::invalid_param("lambda_v", "must be finite"));
        }
        Ok(Self {
            t,
            lambda_so,
            lambda_r,
            lambda_v,
            a_lattice: 2.46e-10,
        })
    }

    /// Graphene with intrinsic SOC added.
    ///
    /// Uses graphene parameters: t = 2.8 eV, a = 2.46 Å, λ_R = 0, λ_v = 0.
    /// With λ_SO = 0 this is pristine graphene (gapless Dirac points).
    /// A non-zero λ_SO opens the Kane-Mele topological gap Δ ≈ 6√3 λ_SO.
    pub fn graphene_with_soc(lambda_so: f64) -> Self {
        Self {
            t: 2.8,
            lambda_so,
            lambda_r: 0.0,
            lambda_v: 0.0,
            a_lattice: 2.46e-10,
        }
    }

    /// Pre-built topological phase (QSH insulator, Z₂ = 1).
    ///
    /// Parameters: t = 1.0 eV, λ_SO = 0.1 eV, λ_R = 0.0, λ_v = 0.0.
    /// Well inside the topological phase: λ_v = 0 < 3√3·0.1 ≈ 0.52 eV.
    pub fn topological_phase() -> Self {
        Self {
            t: 1.0,
            lambda_so: 0.1,
            lambda_r: 0.0,
            lambda_v: 0.0,
            a_lattice: 2.46e-10,
        }
    }

    /// Pre-built trivial insulating phase (Z₂ = 0).
    ///
    /// Parameters: t = 1.0 eV, λ_SO = 0.05 eV, λ_R = 0.0, λ_v = 0.5 eV.
    /// Critical value: 3√3·0.05 ≈ 0.26 eV; λ_v = 0.5 > 0.26, so trivial.
    pub fn trivial_phase() -> Self {
        Self {
            t: 1.0,
            lambda_so: 0.05,
            lambda_r: 0.0,
            lambda_v: 0.5,
            a_lattice: 2.46e-10,
        }
    }

    // -----------------------------------------------------------------------
    // Bloch Hamiltonian
    // -----------------------------------------------------------------------

    /// Construct the 4×4 Bloch Hamiltonian H(kx, ky) in the (A↑, B↑, A↓, B↓) basis.
    ///
    /// The matrix is Hermitian by construction and has the block structure:
    ///
    /// ```text
    /// H = | h_uu  h_ud |    (A↑,B↑ | A↓,B↓ blocks)
    ///     | h_du  h_dd |
    /// ```
    ///
    /// ## Spin-up block h_uu (indices 0,1 for A↑,B↑):
    ///
    /// ```text
    /// h_uu = [[  λ_v + λ_SO·g(k),  -t·f*(k)          ],
    ///         [ -t·f(k),           -λ_v - λ_SO·g(k)  ]]
    /// ```
    ///
    /// where f(k) = NN structure factor, g(k) = NNN SOC factor.
    ///
    /// ## Spin-down block h_dd (indices 2,3 for A↓,B↓):
    ///
    /// Time-reversal partner: h_dd(k) = -h_uu*(-k). Since g(-k) = -g(k) and
    /// f(-k) = f*(k):
    ///
    /// ```text
    /// h_dd = [[ -λ_v + λ_SO·g(k),  -t·f*(k)         ],
    ///         [ -t·f(k),            λ_v - λ_SO·g(k)  ]]
    /// ```
    ///
    /// ## Rashba blocks h_ud, h_du (spin-flip, off-diagonal):
    ///
    /// ```text
    /// h_ud = [[  0,      R_AB  ],
    ///         [ -R*_AB,  0     ]]    with R_AB = λ_R · Σ_j e^{ik·δ_j}(d̂_j_x - i·d̂_j_y)
    /// h_du = h_ud†
    /// ```
    ///
    /// # Arguments
    ///
    /// - `kx`, `ky`: wavevector components in units of (rad / a_lattice), i.e., the
    ///   dimensionless wavevector (multiply by a_lattice to get SI units).
    pub fn hamiltonian_at(&self, kx: f64, ky: f64) -> CMatrix {
        // NN structure factor: f(k) = 1 + e^{ikx} + e^{i(kx/2 + ky√3/2)}
        let f_k = nn_phase_factor(kx, ky);

        // NNN SOC structure factor: g(k) = 2[sin(kx) - sin(kx/2+ky√3/2) + sin(-kx/2+ky√3/2)]
        let g_k = nnn_soc_factor(kx, ky);

        // Off-diagonal Rashba coupling
        let r_ab = rashba_ab_coupling(kx, ky, self.lambda_r);

        // ---------- spin-up block (rows/cols 0,1 = A↑,B↑) ----------
        // h_uu[0,0] = λ_v + λ_SO·g(k)     (A↑ on-site: +λ_v from staggered, +λ_SO·g from NNN)
        // h_uu[1,1] = -λ_v - λ_SO·g(k)    (B↑ on-site: -λ_v sublattice flip, -λ_SO·g TRS)
        // h_uu[0,1] = -t·f*(k)             (A↑→B↑ hopping: -t times complex conjugate of f)
        // h_uu[1,0] = -t·f(k)              (B↑→A↑ hopping: Hermitian conjugate)
        let huu_00 = Complex::from_real(self.lambda_v + self.lambda_so * g_k);
        let huu_11 = Complex::from_real(-self.lambda_v - self.lambda_so * g_k);
        let huu_01 = f_k.conj().scale(-self.t);
        let huu_10 = f_k.scale(-self.t);

        // ---------- spin-down block (rows/cols 2,3 = A↓,B↓) ----------
        // Time-reversal: Θ H(k) Θ† = H(-k), with Θ = iσ_y ⊗ I_2 · K (complex conjugation).
        // For the down-block: h_dd(k) corresponds to flipping spin sector.
        // In the Kane-Mele model with spin conservation (λ_R=0), the down block is:
        //   h_dd(k) = conj of h_uu(-k)  (time-reversal partner)
        //
        // Since g(-k) = -g(k) and f(-k) = f*(k):
        //   h_dd[0,0] = -λ_v + λ_SO·g(k)   [A↓ on-site: -λ_v sublattice, same NNN sign due to σ_z flip]
        //   h_dd[1,1] =  λ_v - λ_SO·g(k)   [B↓ on-site]
        //   h_dd[0,1] = -t·f*(k)            [A↓→B↓ same as up, since hopping is spin-diagonal]
        //   h_dd[1,0] = -t·f(k)
        //
        // This is the correct form: h_dd = h_uu but with σ_z term sign-flipped (λ_SO·g → -λ_SO·g)
        // because σ_z = +1 for up, -1 for down. The staggered potential ξ_i also flips sign
        // if it's spin-independent... actually λ_v is spin-independent (it's a lattice potential).
        // So h_dd[0,0] = λ_v - λ_SO·g ... no, let's be careful:
        //
        // The staggered potential ξ_A = +1, ξ_B = -1, spin-independent:
        //   A↑: +λ_v;  B↑: -λ_v;  A↓: +λ_v;  B↓: -λ_v
        //
        // The SOC: iλ_SO ν_ij σ_z — for spin-up σ_z = +1, for spin-down σ_z = -1.
        //   NNN A↑: +λ_SO·g(k);  NNN B↑: -λ_SO·g(k)  (sublattice sign from ν_ij)
        //   NNN A↓: -λ_SO·g(k);  NNN B↓: +λ_SO·g(k)
        //
        // So:
        //   h_dd[0,0] = λ_v - λ_SO·g(k)      (A↓: +λ_v from potential, -λ_SO·g from spin-down SOC)
        //   h_dd[1,1] = -λ_v + λ_SO·g(k)     (B↓: -λ_v from potential, +λ_SO·g from spin-down SOC)
        //   h_dd[0,1] = -t·f*(k)              (same NN hopping, spin-diagonal)
        //   h_dd[1,0] = -t·f(k)
        let hdd_00 = Complex::from_real(self.lambda_v - self.lambda_so * g_k);
        let hdd_11 = Complex::from_real(-self.lambda_v + self.lambda_so * g_k);
        let hdd_01 = f_k.conj().scale(-self.t);
        let hdd_10 = f_k.scale(-self.t);

        // ---------- off-diagonal Rashba blocks h_ud, h_du ----------
        //
        // The Rashba term couples opposite-spin NN sites. In the (A↑,B↑,A↓,B↓)
        // basis the two 2×2 off-diagonal spin blocks satisfy h_du = h_ud†:
        //
        //   h_ud = [[  0,       R_AB  ],     (rows: A↑,B↑; cols: A↓,B↓)
        //           [ -R_AB*,    0    ]]
        //
        //   h_du = h_ud† = [[  0,     -R_AB  ],   (rows: A↓,B↓; cols: A↑,B↑)
        //                   [  R_AB*,   0    ]]
        //
        // Verification that h_du = h_ud†:
        //   h_du[A↓,A↑] = 0             = conj(h_ud[A↑,A↓]) = conj(0)           ✓
        //   h_du[A↓,B↑] = -R_AB         = conj(h_ud[B↑,A↓]) = conj(-R_AB*) = -R_AB ✓
        //   h_du[B↓,A↑] = R_AB*         = conj(h_ud[A↑,B↓]) = conj(R_AB) = R_AB* ✓
        //   h_du[B↓,B↑] = 0             = conj(h_ud[B↑,B↓]) = conj(0)           ✓
        //
        // Matrix positions (global, basis A↑=0,B↑=1,A↓=2,B↓=3):
        //   h_ud block (upper-right):
        //     (0,2) = h_ud[A↑,A↓] = 0
        //     (0,3) = h_ud[A↑,B↓] = R_AB
        //     (1,2) = h_ud[B↑,A↓] = -R_AB*
        //     (1,3) = h_ud[B↑,B↓] = 0
        //   h_du block (lower-left):
        //     (2,0) = h_du[A↓,A↑] = 0
        //     (2,1) = h_du[A↓,B↑] = -R_AB
        //     (3,0) = h_du[B↓,A↑] = R_AB*
        //     (3,1) = h_du[B↓,B↑] = 0
        let r_conj = r_ab.conj(); // R_AB*

        let mut h = CMatrix::zeros(4);

        // Spin-up block
        h.set(0, 0, huu_00);
        h.set(0, 1, huu_01);
        h.set(1, 0, huu_10);
        h.set(1, 1, huu_11);

        // Spin-down block
        h.set(2, 2, hdd_00);
        h.set(2, 3, hdd_01);
        h.set(3, 2, hdd_10);
        h.set(3, 3, hdd_11);

        // Rashba off-diagonal blocks
        // h_ud entries (upper-right):
        h.set(0, 3, r_ab); // h_ud[A↑,B↓] = R_AB
        h.set(1, 2, r_conj.neg()); // h_ud[B↑,A↓] = -R_AB*
                                   // h_du entries (lower-left) = conj-transpose of h_ud:
        h.set(2, 1, r_ab.neg()); // h_du[A↓,B↑] = -R_AB = conj(h_ud[B↑,A↓]) = conj(-R*) = -R
        h.set(3, 0, r_conj); // h_du[B↓,A↑] = R_AB* = conj(h_ud[A↑,B↓]) = conj(R)

        h
    }

    // -----------------------------------------------------------------------
    // Energy bands
    // -----------------------------------------------------------------------

    /// Diagonalize the 4×4 Bloch Hamiltonian at (kx, ky).
    ///
    /// Returns the four eigenvalues in ascending order. The two lowest eigenvalues
    /// correspond to the occupied bands (at half-filling).
    ///
    /// # Errors
    ///
    /// Returns `NumericalError` if the Jacobi eigendecomposition fails to converge.
    pub fn energy_bands(&self, kx: f64, ky: f64) -> Result<Vec<f64>> {
        let h = self.hamiltonian_at(kx, ky);
        let (evals, _) = h.hermitian_eigendecomposition()?;
        Ok(evals)
    }

    // -----------------------------------------------------------------------
    // Band gap
    // -----------------------------------------------------------------------

    /// Minimum direct gap between occupied and empty bands over the BZ.
    ///
    /// Evaluates the 4×4 spectrum on a 30×30 uniform k-grid covering [−π, π]²
    /// and returns the minimum of (E₂ − E₁), i.e., the gap between the 2nd and
    /// 3rd band (zero-indexed: band indices 1 and 2).
    ///
    /// # Errors
    ///
    /// Propagates diagonalization errors.
    pub fn band_gap(&self) -> Result<f64> {
        let n_grid = 30;
        let mut min_gap = f64::INFINITY;
        for ix in 0..n_grid {
            let kx = -PI + 2.0 * PI * (ix as f64) / (n_grid as f64);
            for iy in 0..n_grid {
                let ky = -PI + 2.0 * PI * (iy as f64) / (n_grid as f64);
                let evals = self.energy_bands(kx, ky)?;
                let gap = evals[2] - evals[1];
                if gap < min_gap {
                    min_gap = gap;
                }
            }
        }
        Ok(min_gap)
    }

    // -----------------------------------------------------------------------
    // Z₂ topological invariant
    // -----------------------------------------------------------------------

    /// Compute the Z₂ invariant using the time-reversal invariant TRIM-point method.
    ///
    /// # Method (λ_R = 0 branch — most common case)
    ///
    /// When the Rashba coupling is absent, the model block-diagonalises into
    /// independent spin-up and spin-down 2×2 Haldane models. The Z₂ invariant
    /// satisfies Z₂ = |C_up| mod 2, where C_up is the Chern number of the
    /// spin-up block.
    ///
    /// The spin-up Chern number is computed via the Fukui-Hatsugai discrete method
    /// using a 20×20 k-grid. This is exact (up to grid resolution) and avoids
    /// the need for the Pfaffian.
    ///
    /// # Method (λ_R ≠ 0 branch — general case)
    ///
    /// When Rashba is non-zero, S_z is not conserved. The Z₂ invariant is
    /// computed via the time-reversal polarization at the four TRIM points of the
    /// rectangular BZ: Γ = (0,0), M₁ = (π,0), M₂ = (0,π), M₃ = (π,π).
    ///
    /// At each TRIM Λ_i, we form the sewing matrix:
    ///   m_{mn}(Λ_i) = ⟨u_m(Λ_i) | Θ | u_n(Λ_i)⟩,  m,n ∈ {occupied bands}
    ///
    /// where Θ = iσ_y ⊗ I₂ (acting on spin ⊗ sublattice) is the time-reversal
    /// operator restricted to the Bloch factor. For a 2×2 sewing matrix of an
    /// antisymmetric form [[0, a], [-a, 0]], the Pfaffian is a, and
    ///   δ_i = Pf(m(Λ_i)) / √|det(m)|.
    ///
    /// The Z₂ invariant is: Z₂ = 1 − Π_i sign(Re δ_i) (topological ↔ odd product).
    ///
    /// # Returns
    ///
    /// - `1` if the system is in the QSH (topological) phase (Z₂ = 1)
    /// - `0` if trivial (Z₂ = 0)
    ///
    /// # Errors
    ///
    /// Returns `NumericalError` if the Chern number grid fails, or if a TRIM-point
    /// eigenstate calculation fails.
    pub fn z2_invariant(&self) -> Result<i32> {
        if self.lambda_r.abs() < 1e-12 {
            // Block-diagonal case: Z₂ = |C_up| mod 2
            self.z2_from_spin_chern()
        } else {
            // General case: TRIM-point Pfaffian method
            self.z2_from_trim_pfaffian()
        }
    }

    /// Z₂ via spin Chern number (λ_R = 0 case).
    ///
    /// Computes the Chern number of the spin-up 2×2 block h_uu(k) using the
    /// Fukui-Hatsugai discrete method. Returns Z₂ = |C_up| mod 2.
    ///
    /// # BZ Sampling
    ///
    /// The honeycomb lattice has primitive reciprocal vectors:
    ///   b₁ = 2π·(1, -1/√3),  b₂ = 2π·(0, 2/√3)
    ///
    /// The correct Brillouin zone is the parallelogram spanned by b₁, b₂.
    /// A k-point in this BZ is parameterised as k = s·b₁ + u·b₂ with s,u ∈ \[0,1\].
    /// Using a square [-π,π]² grid is incorrect for the honeycomb because it does
    /// NOT tile the BZ exactly once; it gives a non-integer Chern number.
    fn z2_from_spin_chern(&self) -> Result<i32> {
        // Primitive reciprocal lattice vectors of the honeycomb:
        // b₁ = 2π·(1, -1/√3),  b₂ = 2π·(0, 2/√3)
        // Parameterise: k(s,u) = s·b₁ + u·b₂, with s,u ∈ [0,1].
        let b1x = 2.0 * PI;
        let b1y = -2.0 * PI / SQRT3;
        let b2x = 0.0_f64;
        let b2y = 4.0 * PI / SQRT3;

        let n = 30_usize; // finer grid for honeycomb BZ
        let mut states: Vec<Vec<[Complex; 2]>> = Vec::with_capacity(n + 1);

        for ix in 0..=n {
            let s = (ix as f64) / (n as f64);
            let mut row = Vec::with_capacity(n + 1);
            for iy in 0..=n {
                let u = (iy as f64) / (n as f64);
                let kx = s * b1x + u * b2x;
                let ky = s * b1y + u * b2y;
                let v = self.spin_up_lower_eigenvec(kx, ky)?;
                row.push(v);
            }
            states.push(row);
        }

        let mut flux_sum = 0.0_f64;
        for ix in 0..n {
            let ix1 = ix + 1;
            for iy in 0..n {
                let iy1 = iy + 1;
                let u_x = link_2(states[ix][iy], states[ix1][iy]);
                let u_y = link_2(states[ix][iy], states[ix][iy1]);
                let u_x_py = link_2(states[ix][iy1], states[ix1][iy1]);
                let u_y_px = link_2(states[ix1][iy], states[ix1][iy1]);
                let plaq = u_x.mul(&u_y_px).mul(&u_x_py.conj()).mul(&u_y.conj());
                flux_sum += plaq.phase();
            }
        }

        let chern_real = flux_sum / (2.0 * PI);
        let chern_int = chern_real.round() as i32;
        if (chern_real - chern_int as f64).abs() > 0.15 {
            return Err(error::numerical_error(&format!(
                "Spin Chern number not quantized: {chern_real:.4} (try parameter adjustment)"
            )));
        }

        Ok(chern_int.abs() % 2)
    }

    /// Lower-band eigenvector of the spin-up 2×2 block at (kx, ky).
    ///
    /// Returns the normalized eigenvector [a, b] for the lower eigenvalue of:
    ///   h_uu = [[α, β*], [β, -α]]  where α = λ_v + λ_SO·g(k), β = -t·f(k).
    fn spin_up_lower_eigenvec(&self, kx: f64, ky: f64) -> Result<[Complex; 2]> {
        let g_k = nnn_soc_factor(kx, ky);
        let f_k = nn_phase_factor(kx, ky);

        let alpha = self.lambda_v + self.lambda_so * g_k; // diagonal element
        let beta = f_k.scale(-self.t); // lower-diagonal element h[1,0]
        let beta_conj = f_k.conj().scale(-self.t); // upper-diagonal element h[0,1]

        // Eigenvalues: ±√(α² + |β|²)
        let disc = (alpha * alpha + beta.norm_sq()).sqrt();
        let e_lo = -disc;

        // Eigenvector for e_lo from (h_uu - e_lo·I)|v⟩ = 0:
        //   (α - e_lo)·v₀ + β*·v₁ = 0
        //   β·v₀ + (-α - e_lo)·v₁ = 0
        // If β ≠ 0: v = [-β*, α - e_lo] (null vector of first row)
        let v = if beta.norm_sq() > 1e-28 {
            let raw = [beta_conj.neg(), Complex::from_real(alpha - e_lo)];
            let norm = (raw[0].norm_sq() + raw[1].norm_sq()).sqrt();
            if norm < 1e-15 {
                return Err(error::numerical_error(
                    "degenerate eigenvector in spin-up block",
                ));
            }
            [raw[0].scale(1.0 / norm), raw[1].scale(1.0 / norm)]
        } else {
            // Diagonal case: α ≤ -α → α < 0 for lower band
            if alpha <= 0.0 {
                [Complex::ONE, Complex::ZERO]
            } else {
                [Complex::ZERO, Complex::ONE]
            }
        };

        Ok(v)
    }

    /// Z₂ via TRIM-point Pfaffian method (general λ_R ≠ 0 case).
    ///
    /// Implements the Fu-Kane formula for Z₂ using time-reversal polarization.
    ///
    /// At each TRIM Λ_i, the 2×2 sewing matrix is:
    ///   m_{mn} = ⟨u_m(Λ_i) | iσ_y ⊗ I₂ | u_n(Λ_i)⟩
    ///
    /// For a 2-band occupied subspace (n_occ = 2), the matrix m is antisymmetric
    /// and 2×2. Its Pfaffian is Pf[[0,a],[-a,0]] = a, computed analytically.
    /// The Z₂ parity is δ_i = sign(Pf(m_i) / √|det m_i|) = ±1.
    ///
    /// Z₂ = 1 if Π_i δ_i = -1 (odd number of -1's → topological).
    fn z2_from_trim_pfaffian(&self) -> Result<i32> {
        // TRIM points for the rectangular BZ [−π,π]²: (kx,ky) ∈ {0,π}²
        let trim_points = [(0.0, 0.0), (PI, 0.0), (0.0, PI), (PI, PI)];

        // Time-reversal operator matrix Θ = iσ_y ⊗ I₂ in (A↑,B↑,A↓,B↓) basis:
        //   iσ_y = [[0, 1],[-1,0]] in spin space ⊗ I₂ in sublattice space.
        //
        // In (A↑,B↑,A↓,B↓) ordering:
        //   (A↑,A↓) = spin block for A: row A↑ = 0, row A↓ = 2
        //   (B↑,B↓) = spin block for B: row B↑ = 1, row B↓ = 3
        //
        //   iσ_y ⊗ I₂:
        //   σ_y acts on spin indices {↑,↓}, I₂ acts on sublattice {A,B}.
        //   In our basis (A↑=0, B↑=1, A↓=2, B↓=3):
        //   Row 0 (A↑): (iσ_y)_{↑,↑}·I=0, (iσ_y)_{↑,↑}=0 for spin↑→spin↑,
        //                (iσ_y)_{↑,↓}=1 for spin↑→spin↓ → connects (A↑) to (A↓,B↓)
        //
        //   Actually in the Kronecker product basis ordering, if basis = {A↑,B↑,A↓,B↓}
        //   = {sublattice} ⊗ {spin}, and iσ_y acts on spin part:
        //   In block form: [[0, iσ_y_↑↓·I₂], [iσ_y_↓↑·I₂, 0]]
        //                = [[0·I₂, (+1)·I₂], [(-1)·I₂, 0·I₂]]
        //
        //   But our basis ordering is (A↑,B↑,A↓,B↓), which is (spinup orbitals, spindown orbitals).
        //   In this ordering, iσ_y ⊗ I₂ = [[0, +I₂], [-I₂, 0]] in 2x2 block form:
        //   Matrix:
        //     (A↑→A↓): 0, (A↑→B↓): 0  → row 0, (col 2=+1, col 3=0)
        //     (B↑→A↓): 0, (B↑→B↓): 0  → row 1, (col 2=0, col 3=+1)
        //     (A↓→A↑): -1, (A↓→B↑): 0 → row 2, (col 0=-1, col 1=0)
        //     (B↓→A↑): 0, (B↓→B↑): -1 → row 3, (col 0=0, col 1=-1)
        //
        //   So: theta_mat[i,j]:
        //     (0,2)=+1, (0,3)=0,  (1,2)=0,  (1,3)=+1  (upper-right block = +I₂)
        //     (2,0)=-1, (2,1)=0,  (3,0)=0,  (3,1)=-1  (lower-left block = -I₂)
        //     All other elements = 0.
        //
        // This is the matrix for (iσ_y) ⊗ I₂ in (A↑,B↑,A↓,B↓) = (spin_up block, spin_down block).

        let mut delta_product = 1.0_f64;

        for &(kx, ky) in &trim_points {
            let h = self.hamiltonian_at(kx, ky);
            let (evals, vecs) = h.hermitian_eigendecomposition()?;

            // Two occupied bands (indices 0 and 1, lowest eigenvalues)
            // Extract occupied eigenvectors as 4-component complex vectors
            let u0: Vec<Complex> = (0..4).map(|r| vecs.get(r, 0)).collect();
            let u1: Vec<Complex> = (0..4).map(|r| vecs.get(r, 1)).collect();

            // Verify we're in the gapped regime (gap between bands 1 and 2 > threshold)
            let gap_01 = evals[2] - evals[1];
            if gap_01 < 1e-8 {
                // Near band touching: return error (gapless point)
                return Err(error::numerical_error(&format!(
                    "Band gap too small at TRIM ({kx:.3},{ky:.3}): gap={gap_01:.2e}. \
                     System may be at a topological phase transition."
                )));
            }

            // Apply theta = iσ_y ⊗ I₂ to u0 and u1:
            // theta · u: new[0] = u[2], new[1] = u[3], new[2] = -u[0], new[3] = -u[1]
            let theta_u0 = apply_time_reversal(&u0);
            let theta_u1 = apply_time_reversal(&u1);

            // Sewing matrix m_{mn} = ⟨u_m | theta | u_n⟩ = conj(u_m) · (theta · u_n)
            // m[0,0] = ⟨u0|theta|u0⟩
            // m[0,1] = ⟨u0|theta|u1⟩
            // m[1,0] = ⟨u1|theta|u0⟩
            // m[1,1] = ⟨u1|theta|u1⟩
            //
            // For a time-reversal symmetric Hamiltonian at TRIM, m is antisymmetric:
            // m[0,0] = 0, m[1,1] = 0, m[0,1] = -m[1,0].
            let m01 = inner_product(&u0, &theta_u1);
            let m10 = inner_product(&u1, &theta_u0);

            // Pfaffian of [[0, m01],[-m01, 0]] = m01
            let pf = m01;

            // det(m) = m01 · m10 - 0 = m01·(-m01) = -m01² ... wait:
            // det[[0,a],[-a,0]] = 0*0 - a*(-a) = a²
            // So det(m) = m01 · (-m10_sign) ... let's compute:
            // det = m00*m11 - m01*m10 = 0*0 - m01*m10
            let det_m = m01.mul(&m10).neg(); // = -m01*m10

            // For antisymmetric 2×2: det = pf² = m01² (if m10 = -m01)
            // In our case m10 = -m01 by antisymmetry, so det = m01*m01 = m01²
            // sqrt(|det|) = |m01|
            let sqrt_det = det_m.re.abs().sqrt();
            if sqrt_det < 1e-12 {
                // Singular sewing matrix: system may be gapless
                return Err(error::numerical_error(
                    "Sewing matrix singular at TRIM point; cannot compute Z₂",
                ));
            }

            // δ_i = Pf(m) / √|det(m)| = m01 / |m01| = unit complex on unit circle
            let delta_re = pf.re / sqrt_det;
            delta_product *= delta_re.signum();
        }

        // Z₂ = 1 (topological) if Π_i δ_i = -1 (odd number of -1 factors)
        // Return 1 for topological, 0 for trivial.
        Ok(if delta_product < 0.0 { 1 } else { 0 })
    }

    // -----------------------------------------------------------------------
    // Topological phase test
    // -----------------------------------------------------------------------

    /// Returns `true` if the system is in the QSH (topological) phase.
    ///
    /// Calls `z2_invariant()` and returns `true` iff Z₂ = 1.
    /// If the calculation fails (e.g., at a phase boundary), returns `false`.
    pub fn is_topological(&self) -> bool {
        self.z2_invariant().is_ok_and(|z2| z2 == 1)
    }

    // -----------------------------------------------------------------------
    // Edge spectrum
    // -----------------------------------------------------------------------

    /// Compute the edge spectrum of a zigzag strip with open y-boundary conditions.
    ///
    /// Builds a strip Hamiltonian for `n_cells` unit cells in the y-direction.
    /// Each unit cell contains 4 orbitals (A↑, B↑, A↓, B↓), giving a
    /// `4·n_cells × 4·n_cells` strip matrix. The strip is periodic in x (kx is good)
    /// and has open boundary conditions in y.
    ///
    /// # Strip Hamiltonian Construction
    ///
    /// The strip Hamiltonian at fixed kx is block-tridiagonal:
    /// ```text
    /// H_strip(kx) = diag(H_intra) + superdiag(H_inter) + subdiag(H_inter†)
    /// ```
    /// where:
    /// - `H_intra(kx) = H(kx, ky=0)`: the 4×4 in-cell Hamiltonian (kx-direction hoppings)
    /// - `H_inter`: inter-cell coupling in y, extracted from `(H(kx, δky) − H_intra) · 0.5`
    ///
    /// # Arguments
    ///
    /// - `n_cells`: strip width (number of unit cells in y). Max 16 (→ 64×64 matrix).
    /// - `kx_min`, `kx_max`: range of kx to sample.
    /// - `n_kx`: number of kx points.
    ///
    /// # Returns
    ///
    /// Vec of `(kx, energies)` with `4·n_cells` energies sorted ascending.
    ///
    /// # Errors
    ///
    /// - `InvalidParameter` if `n_cells < 2`, `n_cells > 16`, or `n_kx < 2`.
    /// - `InvalidParameter` if `n_cells * 4 > CMatrix::MAX_DIM` (64).
    pub fn edge_spectrum(
        &self,
        n_cells: usize,
        kx_min: f64,
        kx_max: f64,
        n_kx: usize,
    ) -> Result<Vec<(f64, Vec<f64>)>> {
        if n_cells < 2 {
            return Err(error::invalid_param(
                "n_cells",
                "strip width must be at least 2",
            ));
        }
        let dim = n_cells * 4;
        if dim > CMatrix::MAX_DIM {
            return Err(error::invalid_param(
                "n_cells",
                "n_cells * 4 exceeds CMatrix::MAX_DIM (64); use n_cells ≤ 16",
            ));
        }
        if n_kx < 2 {
            return Err(error::invalid_param("n_kx", "need at least 2 kx points"));
        }

        // Small ky offset for extracting inter-cell hopping numerically
        let delta_ky = PI / (n_cells.max(4) as f64);

        let mut result = Vec::with_capacity(n_kx);
        for i in 0..n_kx {
            let kx = if n_kx == 1 {
                kx_min
            } else {
                kx_min + (kx_max - kx_min) * (i as f64) / ((n_kx - 1) as f64)
            };

            let h_strip = self.build_strip_hamiltonian(kx, n_cells, delta_ky)?;
            let (evals, _) = h_strip.hermitian_eigendecomposition()?;
            result.push((kx, evals));
        }

        Ok(result)
    }

    /// Build the block-tridiagonal strip Hamiltonian at fixed kx.
    fn build_strip_hamiltonian(&self, kx: f64, n_cells: usize, delta_ky: f64) -> Result<CMatrix> {
        let nb = 4; // orbitals per unit cell: A↑,B↑,A↓,B↓
        let dim = n_cells * nb;

        // Intra-cell block: H(kx, ky=0)
        let h_intra = self.hamiltonian_at(kx, 0.0);

        // Inter-cell coupling extracted from finite difference:
        // H_inter ≈ (H(kx, δky) − H(kx, 0)) · 0.5
        // This captures the leading ky-derivative of H(k), representing
        // the inter-cell hopping amplitude in the y-direction.
        let h_at_dky = self.hamiltonian_at(kx, delta_ky);
        let mut h_strip = CMatrix::zeros(dim);

        // Place intra-cell blocks on diagonal
        for iy in 0..n_cells {
            for i in 0..nb {
                for j in 0..nb {
                    let row = iy * nb + i;
                    let col = iy * nb + j;
                    let cur = h_strip.get(row, col);
                    h_strip.set(row, col, cur.add(&h_intra.get(i, j)));
                }
            }
        }

        // Place inter-cell blocks (open boundary: no coupling from cell n_cells-1 to 0)
        for iy in 0..(n_cells - 1) {
            for i in 0..nb {
                for j in 0..nb {
                    // Inter-cell hopping amplitude
                    let t_ij = {
                        let dh = h_at_dky.get(i, j).sub(&h_intra.get(i, j));
                        dh.scale(0.5)
                    };
                    // Super-diagonal: cell iy → cell iy+1
                    let r_up = iy * nb + i;
                    let c_up = (iy + 1) * nb + j;
                    let cur_up = h_strip.get(r_up, c_up);
                    h_strip.set(r_up, c_up, cur_up.add(&t_ij));
                    // Sub-diagonal: cell iy+1 → cell iy (Hermitian conjugate)
                    let r_dn = (iy + 1) * nb + i;
                    let c_dn = iy * nb + j;
                    let cur_dn = h_strip.get(r_dn, c_dn);
                    h_strip.set(r_dn, c_dn, cur_dn.add(&t_ij.conj()));
                }
            }
        }

        Ok(h_strip)
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Apply the time-reversal operator Θ = iσ_y ⊗ I₂ to a 4-component state.
///
/// In the (A↑, B↑, A↓, B↓) basis, the matrix representation of iσ_y ⊗ I₂ is:
/// ```text
/// [[  0,  0, +1,  0 ],
///  [  0,  0,  0, +1 ],
///  [ -1,  0,  0,  0 ],
///  [  0, -1,  0,  0 ]]
/// ```
/// This maps (v₀, v₁, v₂, v₃) → (+v₂, +v₃, −v₀, −v₁).
fn apply_time_reversal(u: &[Complex]) -> Vec<Complex> {
    debug_assert_eq!(u.len(), 4);
    vec![u[2], u[3], u[0].neg(), u[1].neg()]
}

/// Complex inner product ⟨a|b⟩ = Σ_i a_i* · b_i.
fn inner_product(a: &[Complex], b: &[Complex]) -> Complex {
    debug_assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(ai, bi)| ai.conj().mul(bi))
        .fold(Complex::ZERO, |acc, x| acc.add(&x))
}

/// Link variable for 2-component eigenvector array (Fukui-Hatsugai step).
///
/// Returns U = ⟨u_a|u_b⟩ / |⟨u_a|u_b⟩|, a pure phase.
#[inline]
fn link_2(a: [Complex; 2], b: [Complex; 2]) -> Complex {
    let inner = a[0].conj().mul(&b[0]).add(&a[1].conj().mul(&b[1]));
    let norm = inner.norm();
    if norm < 1e-15 {
        Complex::ONE
    } else {
        inner.scale(1.0 / norm)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Relative or absolute tolerance helpers
    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // -----------------------------------------------------------------------
    // Test 1: topological_phase gap > 0
    // -----------------------------------------------------------------------
    #[test]
    fn topological_phase_has_positive_gap() {
        let model = KaneMeleModel::topological_phase();
        let gap = model.band_gap().unwrap();
        assert!(
            gap > 0.0,
            "Topological phase should have positive band gap, got {gap}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 2: trivial_phase also gapped (different topology)
    // -----------------------------------------------------------------------
    #[test]
    fn trivial_phase_has_positive_gap() {
        let model = KaneMeleModel::trivial_phase();
        let gap = model.band_gap().unwrap();
        assert!(
            gap > 0.0,
            "Trivial phase should also have positive band gap, got {gap}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 3: topological and trivial phases differ in Z₂
    // -----------------------------------------------------------------------
    #[test]
    fn topological_and_trivial_differ_in_z2() {
        let topo = KaneMeleModel::topological_phase();
        let trivial = KaneMeleModel::trivial_phase();
        let z2_topo = topo.z2_invariant().unwrap();
        let z2_trivial = trivial.z2_invariant().unwrap();
        assert_ne!(
            z2_topo, z2_trivial,
            "Topological (Z₂={z2_topo}) and trivial (Z₂={z2_trivial}) should differ"
        );
    }

    // -----------------------------------------------------------------------
    // Test 4: H(k) is Hermitian at 3 different k-points
    // -----------------------------------------------------------------------
    #[test]
    fn hamiltonian_is_hermitian() {
        // Test both no-Rashba and with-Rashba models
        let models = [
            KaneMeleModel::new(1.0, 0.1, 0.0, 0.0).unwrap(),
            KaneMeleModel::new(1.0, 0.1, 0.05, 0.1).unwrap(),
            KaneMeleModel::new(2.8, 0.05, 0.02, 0.0).unwrap(),
        ];
        let kpoints = [(0.0, 0.0), (0.5, 0.8), (-1.2, 1.7)];

        for model in &models {
            for &(kx, ky) in &kpoints {
                let h = model.hamiltonian_at(kx, ky);
                let hd = h.conj_transpose();
                let diff = h.sub(&hd).unwrap();
                assert!(
                    diff.frobenius_norm() < 1e-12,
                    "H not Hermitian at k=({kx},{ky}): ||H-H†||={:.2e}",
                    diff.frobenius_norm()
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Test 5: energy_bands returns 4 sorted eigenvalues
    // -----------------------------------------------------------------------
    #[test]
    fn energy_bands_returns_four_sorted_values() {
        let model = KaneMeleModel::topological_phase();
        let evals = model.energy_bands(0.3, 0.7).unwrap();
        assert_eq!(evals.len(), 4, "Must have exactly 4 bands");
        for w in evals.windows(2) {
            assert!(w[0] <= w[1] + 1e-10, "Eigenvalues not sorted: {evals:?}");
        }
    }

    // -----------------------------------------------------------------------
    // Test 6: At Γ (k=0,0), eigenvalues have known spin-degeneracy structure
    // -----------------------------------------------------------------------
    #[test]
    fn gamma_point_kramers_degeneracy() {
        // At k=0 with TRS, eigenvalues should come in Kramers pairs (doubly degenerate)
        // because Θ² = -1 for spin-½. At the Γ point: H(0) and Θ H(0) Θ⁻¹ = H(0)
        // with Θ² = -1 → Kramers degeneracy.
        let model = KaneMeleModel::topological_phase();
        let evals = model.energy_bands(0.0, 0.0).unwrap();
        assert_eq!(evals.len(), 4);
        // Kramers pairs: evals[0]=evals[1] and evals[2]=evals[3]
        // (for λ_R=0: exact due to spin degeneracy; for λ_R≠0: degeneracy at TRIM by TRS)
        assert!(
            approx(evals[0], evals[1], 1e-10),
            "Γ: Kramers pair 1 not degenerate: {} vs {}",
            evals[0],
            evals[1]
        );
        assert!(
            approx(evals[2], evals[3], 1e-10),
            "Γ: Kramers pair 2 not degenerate: {} vs {}",
            evals[2],
            evals[3]
        );
    }

    // -----------------------------------------------------------------------
    // Test 7: is_topological() returns true for topological, false for trivial
    // -----------------------------------------------------------------------
    #[test]
    fn is_topological_classifies_correctly() {
        let topo = KaneMeleModel::topological_phase();
        let trivial = KaneMeleModel::trivial_phase();
        assert!(
            topo.is_topological(),
            "topological_phase() should be classified as topological"
        );
        assert!(
            !trivial.is_topological(),
            "trivial_phase() should be classified as trivial"
        );
    }

    // -----------------------------------------------------------------------
    // Test 8: band_gap decreases as λ_v increases toward critical value
    // -----------------------------------------------------------------------
    #[test]
    fn gap_decreases_as_lambda_v_increases() {
        // For t=1, λ_SO=0.1, λ_R=0: gap closes at λ_v = 3√3·0.1 ≈ 0.5196
        // Test that gap at λ_v=0.0 > gap at λ_v=0.3 > gap at λ_v=0.5 (all below critical)
        let lso = 0.1;
        let lv_values = [0.0, 0.2, 0.4];
        let mut gaps = Vec::new();
        for &lv in &lv_values {
            let m = KaneMeleModel::new(1.0, lso, 0.0, lv).unwrap();
            gaps.push(m.band_gap().unwrap());
        }
        // Gap should be strictly decreasing as λ_v increases (for λ_R=0)
        assert!(
            gaps[0] > gaps[1],
            "Gap should decrease with λ_v: gap(0.0)={} > gap(0.2)={} failed",
            gaps[0],
            gaps[1]
        );
        assert!(
            gaps[1] > gaps[2],
            "Gap should decrease with λ_v: gap(0.2)={} > gap(0.4)={} failed",
            gaps[1],
            gaps[2]
        );
    }

    // -----------------------------------------------------------------------
    // Test 9: graphene_with_soc(0.0) is gapless (Dirac point)
    // -----------------------------------------------------------------------
    #[test]
    fn graphene_no_soc_is_gapless() {
        let model = KaneMeleModel::graphene_with_soc(0.0);
        // Graphene without SOC has Dirac cones — minimum gap should be ~0.
        // We evaluate at the K point: K = (4π/3, 0) in our coordinates.
        let k_x = 4.0 * PI / 3.0;
        let evals = model.energy_bands(k_x, 0.0).unwrap();
        // Bands 1 and 2 should touch at K (or very close): |E₂ - E₁| → 0
        let gap_at_k = evals[2] - evals[1];
        assert!(
            gap_at_k < 1e-6,
            "Graphene without SOC: gap at K should be ~0, got {gap_at_k}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 10: graphene_with_soc(0.1) opens a gap
    // -----------------------------------------------------------------------
    #[test]
    fn graphene_with_soc_opens_gap() {
        let model = KaneMeleModel::graphene_with_soc(0.1);
        let gap = model.band_gap().unwrap();
        assert!(
            gap > 0.0,
            "SOC should open a band gap in graphene, got gap={gap}"
        );
        // Gap ~ 6√3·λ_SO ≈ 6√3·0.1 ≈ 1.04 eV (for the honeycomb KM model)
        // Check it's at least a reasonable fraction of this
        assert!(
            gap > 0.1,
            "Gap with λ_SO=0.1 should be substantial (>0.1 eV), got {gap}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 11: edge_spectrum returns correct number of kx points
    // -----------------------------------------------------------------------
    #[test]
    fn edge_spectrum_correct_kx_count() {
        let model = KaneMeleModel::topological_phase();
        let n_kx = 15;
        let n_cells = 6;
        let result = model.edge_spectrum(n_cells, -PI, PI, n_kx).unwrap();
        assert_eq!(
            result.len(),
            n_kx,
            "Expected {n_kx} kx points, got {}",
            result.len()
        );
        // Each kx point should have 4*n_cells bands
        let expected_bands = 4 * n_cells;
        for (kx, bands) in &result {
            assert_eq!(
                bands.len(),
                expected_bands,
                "At kx={kx:.3}: expected {expected_bands} bands, got {}",
                bands.len()
            );
        }
    }

    // -----------------------------------------------------------------------
    // Test 12: H(k+G) ≈ H(k) (lattice periodicity with reciprocal vector G)
    // -----------------------------------------------------------------------
    #[test]
    fn hamiltonian_is_bz_periodic() {
        // Primitive reciprocal lattice vectors for the honeycomb with
        // a₁=(1,0), a₂=(½,√3/2) satisfy aᵢ·bⱼ = 2π·δᵢⱼ:
        //   b₁ = 2π·(1, -1/√3)
        //   b₂ = 2π·(0,  2/√3)
        //
        // Under k → k + b₁ the NN and NNN structure factors are invariant:
        //   f(k+b₁): e^{i(kx+2π)} = e^{ikx}, e^{i((kx+2π)/2+(ky-2π/√3)√3/2)}
        //            = e^{ikx/2+iπ}·e^{iky√3/2-iπ} = e^{ikx/2+iky√3/2} ✓
        //   g(k+b₁): φ₁→φ₁+2π, φ₂→φ₂, (φ₂-φ₁)→(φ₂-φ₁-2π) → sines unchanged ✓
        //
        // So H(k+b₁) = H(k) and H(k+b₂) = H(k) exactly.
        let model = KaneMeleModel::new(1.0, 0.1, 0.03, 0.05).unwrap();
        let kx_test = 0.4;
        let ky_test = 0.3;

        // Test with b₁ = (2π, -2π/√3)
        let b1x = 2.0 * PI;
        let b1y = -2.0 * PI / SQRT3;
        let h_k = model.hamiltonian_at(kx_test, ky_test);
        let h_kpb1 = model.hamiltonian_at(kx_test + b1x, ky_test + b1y);
        let diff1 = h_k.sub(&h_kpb1).unwrap();
        assert!(
            diff1.frobenius_norm() < 1e-10,
            "b₁ periodicity violated: ||H(k) - H(k+b₁)||={:.2e}",
            diff1.frobenius_norm()
        );

        // Also test with b₂ = (0, 4π/√3)
        let b2x = 0.0;
        let b2y = 4.0 * PI / SQRT3;
        let h_kpb2 = model.hamiltonian_at(kx_test + b2x, ky_test + b2y);
        let diff2 = h_k.sub(&h_kpb2).unwrap();
        assert!(
            diff2.frobenius_norm() < 1e-10,
            "b₂ periodicity violated: ||H(k) - H(k+b₂)||={:.2e}",
            diff2.frobenius_norm()
        );
    }

    // -----------------------------------------------------------------------
    // Test 13: Time-reversal symmetry H(-k) = Θ H(k) Θ†
    // -----------------------------------------------------------------------
    #[test]
    fn time_reversal_symmetry() {
        // The Kane-Mele model satisfies time-reversal symmetry:
        //   Θ H(k) Θ† = H(-k)
        // where Θ = iσ_y ⊗ I₂ acts as complex conjugation in the Bloch factor.
        //
        // In matrix form: Θ H(k) Θ† = Θ_mat · H(k)* · Θ_mat†
        // (since Θ = U·K with U = iσ_y ⊗ I₂ and K = complex conjugation,
        //  Θ A Θ† = U A* U†)
        //
        // Test: Θ_mat · H(kx,ky)* · Θ_mat† = H(-kx,-ky)
        let model = KaneMeleModel::new(1.0, 0.08, 0.0, 0.06).unwrap();
        let test_points = [(0.5, 0.3), (-0.7, 0.9), (1.1, -0.4)];

        // Θ_mat in (A↑,B↑,A↓,B↓) basis: [[0,0,+1,0],[0,0,0,+1],[-1,0,0,0],[0,-1,0,0]]
        // Θ_mat† = Θ_mat⁻¹ = -Θ_mat (since Θ_mat² = -I → Θ_mat† = -Θ_mat for unitary Θ_mat)
        // Let's verify: Θ_mat · H*(k) · Θ_mat† = H(-k)

        for &(kx, ky) in &test_points {
            let h_k = model.hamiltonian_at(kx, ky);
            let h_mk = model.hamiltonian_at(-kx, -ky);

            // Build Θ_mat · H*(k) · Θ_mat†:
            // H*(k) is the elementwise complex conjugate (NOT Hermitian conjugate)
            // Step 1: H_conj = H(k)*
            let mut h_conj = CMatrix::zeros(4);
            for i in 0..4 {
                for j in 0..4 {
                    h_conj.set(i, j, h_k.get(i, j).conj());
                }
            }

            // Step 2: Θ_mat · H_conj · Θ_mat†
            // Θ_mat maps: row 0 → row 2 (× +1), row 1 → row 3 (× +1), row 2 → row 0 (× -1), row 3 → row 1 (× -1)
            // Θ_mat† maps: col 0 → col 2 (× +1), col 1 → col 3 (× +1), col 2 → col 0 (× -1), col 3 → col 1 (× -1)
            // More explicitly: (Θ H* Θ†)[i,j] = Σ_{ab} Θ[i,a] H*[a,b] Θ†[b,j]
            //                                  = Σ_{ab} Θ[i,a] H*[a,b] conj(Θ[j,b])
            // Since Θ_mat is purely real: Θ†[b,j] = Θ[j,b]^T = Θ[b,j] (for Θ_mat real it's the transpose).
            // Wait: Θ_mat is real so Θ_mat† = Θ_mat^T.
            // Θ_mat: rows are standard basis actions:
            //   Θ_mat[0,2]=+1, Θ_mat[1,3]=+1, Θ_mat[2,0]=-1, Θ_mat[3,1]=-1, rest 0
            // Θ_mat^T: Θ_mat^T[2,0]=+1, Θ_mat^T[3,1]=+1, Θ_mat^T[0,2]=-1, Θ_mat^T[1,3]=-1
            //
            // (Θ H* Θ^T)[i,j] = Σ_{ab} Θ[i,a] · H*[a,b] · Θ^T[b,j]
            //                  = Σ_{ab} Θ[i,a] · H*[a,b] · Θ[j,b]
            //
            // Nonzero entries of Θ[i,a]: (i=0,a=2,+1),(i=1,a=3,+1),(i=2,a=0,-1),(i=3,a=1,-1)
            // Nonzero entries of Θ[j,b]: (j=0,b=2,+1),(j=1,b=3,+1),(j=2,b=0,-1),(j=3,b=1,-1)
            //
            // (Θ H* Θ^T)[i,j] = θ_i · H*[a(i), b(j)] · θ_j
            // where θ_0=+1,a(0)=2; θ_1=+1,a(1)=3; θ_2=-1,a(2)=0; θ_3=-1,a(3)=1

            let theta_sign: [f64; 4] = [1.0, 1.0, -1.0, -1.0];
            let theta_map: [usize; 4] = [2, 3, 0, 1];

            let mut thr = CMatrix::zeros(4);
            for i in 0..4 {
                for j in 0..4 {
                    let ai = theta_map[i];
                    let bj = theta_map[j];
                    let val = h_conj.get(ai, bj).scale(theta_sign[i] * theta_sign[j]);
                    thr.set(i, j, val);
                }
            }

            let diff = thr.sub(&h_mk).unwrap();
            assert!(
                diff.frobenius_norm() < 1e-10,
                "TRS violated at k=({kx},{ky}): ||Θ H(k)* Θ^T - H(-k)||={:.2e}",
                diff.frobenius_norm()
            );
        }
    }
}
