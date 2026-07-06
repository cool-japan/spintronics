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
use crate::topomagnon::wilson::{self, WilsonLoop};

// ---------------------------------------------------------------------------
// Honeycomb geometry constants and helpers
// ---------------------------------------------------------------------------

const SQRT3: f64 = 1.732_050_808_568_877_3;

/// The three individual NN Bloch phase factors that sum to [`nn_phase_factor`].
///
/// Using primitive vectors a₁ = (1, 0), a₂ = (½, √3/2) in units of lattice
/// constant a=1, returns `[p₀, p₁, p₂]` for the lattice shifts (0,0), a₁, a₂
/// respectively:
///   p₀ = e^{i·0}     = 1
///   p₁ = e^{i k·a₁}  = e^{i·kx}
///   p₂ = e^{i k·a₂}  = e^{i(kx/2 + ky·√3/2)}
///
/// Shared by [`nn_phase_factor`] (which sums them) and [`rashba_ab_coupling`]
/// (which needs the individual phases, not just their sum).
#[inline]
fn nn_phases(kx: f64, ky: f64) -> [Complex; 3] {
    let p0 = Complex::ONE;
    let p1 = Complex::new(0.0, kx).exp();
    let p2 = Complex::new(0.0, kx * 0.5 + ky * SQRT3 * 0.5).exp();
    [p0, p1, p2]
}

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
    let [p0, p1, p2] = nn_phases(kx, ky);
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

/// Primitive reciprocal lattice vectors `(b₁, b₂)` of the honeycomb lattice,
/// dual to the real-space primitive vectors `a₁=(1,0)`, `a₂=(½,√3/2)` used
/// throughout this module (`aᵢ·bⱼ = 2π·δᵢⱼ`).
///
/// Returns `b₁ = (2π, -2π/√3)`, `b₂ = (0, 4π/√3)`.
#[inline]
fn honeycomb_reciprocal_vectors() -> ((f64, f64), (f64, f64)) {
    let b1 = (2.0 * PI, -2.0 * PI / SQRT3);
    let b2 = (0.0, 4.0 * PI / SQRT3);
    (b1, b2)
}

/// The four time-reversal-invariant momenta (TRIM points) of the honeycomb
/// Brillouin zone spanned by [`honeycomb_reciprocal_vectors`]: `Γ = (0,0)`,
/// `b₁/2`, `b₂/2`, and `(b₁+b₂)/2`.
///
/// These are the correct TRIM points for this oblique honeycomb BZ — *not*
/// the corners of a rectangular `[-π,π]²` cell, which do not coincide with
/// TRS-invariant momenta of this lattice except at `Γ`.
///
/// Used by [`KaneMeleModel::z2_from_trim_pfaffian`] (a retired, test-only
/// cross-check) and directly by regression tests; not on the
/// [`KaneMeleModel::z2_invariant`] hot path, hence `#[allow(dead_code)]`.
#[allow(dead_code)]
#[inline]
fn honeycomb_trim_points() -> [(f64, f64); 4] {
    let (b1, b2) = honeycomb_reciprocal_vectors();
    [
        (0.0, 0.0),
        (b1.0 * 0.5, b1.1 * 0.5),
        (b2.0 * 0.5, b2.1 * 0.5),
        ((b1.0 + b2.0) * 0.5, (b1.1 + b2.1) * 0.5),
    ]
}

/// Rashba spin-orbit coupling matrix element `R(k) ≡ H[A↑,B↓](k)`.
///
/// Computes **only** the `A↑ → B↓` matrix element of the Kane-Mele Rashba term
/// `H_R = iλ_R Σ_{<ij>} c†_i (σ × d̂_ij)_z c_j`. It does *not* attempt to also
/// give the `B↑ → A↓` element — that element is `-R(-k)` (obtained by
/// Hermitian-conjugating the real-space A→B hopping term), which is *not*
/// equal to `-conj(R(k))` in general because `R(k)` carries an overall factor
/// of `i` and is therefore not real-analytic. See
/// [`KaneMeleModel::hamiltonian_at`] for the block assembly that uses the
/// correct `-R(-k)` relation.
///
/// # Formula
///
/// Pairing each of the three NN bond phases from [`nn_phases`] with the
/// corresponding Rashba coefficient `χ_m`:
///
/// | shift `m` | phase `p_m`               | `χ_m`           |
/// |-----------|---------------------------|-----------------|
/// | `(0,0)`   | `1`                       | `-1/2 − i·√3/2` |
/// | `a₁`      | `e^{i·kx}`                | `-1/2 + i·√3/2` |
/// | `a₂`      | `e^{i(kx/2 + ky·√3/2)}`   | `1`             |
///
/// `R(k) = i·λ_R · Σ_m conj(p_m)·χ_m`.
#[inline]
fn rashba_ab_coupling(kx: f64, ky: f64, lambda_r: f64) -> Complex {
    if lambda_r.abs() < 1e-15 {
        return Complex::ZERO;
    }
    let [p0, p1, p2] = nn_phases(kx, ky);

    let chi0 = Complex::new(-0.5, -SQRT3 * 0.5);
    let chi1 = Complex::new(-0.5, SQRT3 * 0.5);
    let chi2 = Complex::new(1.0, 0.0);

    let sum = p0
        .conj()
        .mul(&chi0)
        .add(&p1.conj().mul(&chi1))
        .add(&p2.conj().mul(&chi2));

    // R(k) = i·λ_R·Σ_m conj(p_m)·χ_m
    sum.scale(lambda_r).mul_i()
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
    /// h_ud = [[  0,      R(k)  ],
    ///         [ -R(-k),  0     ]]    with R(k) = H[A↑,B↓](k) (see [`rashba_ab_coupling`])
    /// h_du = h_ud†
    /// ```
    ///
    /// Note `H[B↑,A↓](k) = -R(-k)`, **not** `-conj(R(k))` — the two differ
    /// because `R(k)` carries an overall factor of `i` and so is not
    /// real-analytic. The `-R(-k)` relation follows from Hermitian-conjugating
    /// the real-space A→B Rashba hopping term.
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

        // Off-diagonal Rashba coupling: R(k) = H[A↑,B↓](k) and the Hermitian-
        // conjugate-derived partner S(k) = H[B↑,A↓](k) = -R(-k) (NOT -conj(R(k)) —
        // see doc comment above and on `rashba_ab_coupling`).
        let r_ab = rashba_ab_coupling(kx, ky, self.lambda_r);
        let s_ab = rashba_ab_coupling(-kx, -ky, self.lambda_r).neg();

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
        // The Rashba term couples opposite-spin NN sites. Let R(k) ≡ H[A↑,B↓](k)
        // (computed by `rashba_ab_coupling`). Hermitian-conjugating the real-space
        // A→B hopping term gives the spin-flipped element:
        //
        //   H[B↑,A↓](k) = -R(-k)      (NOT -conj(R(k)) — R(k) is not real-analytic
        //                               because of its overall factor of i)
        //
        // In the (A↑,B↑,A↓,B↓) basis the two 2×2 off-diagonal spin blocks are:
        //
        //   h_ud = [[  0,       R(k)   ],     (rows: A↑,B↑; cols: A↓,B↓)
        //           [ -R(-k),    0     ]]
        //
        //   h_du = h_ud† = [[  0,           conj(-R(-k))  ],   (rows: A↓,B↓; cols: A↑,B↑)
        //                   [  conj(R(k)),   0             ]]
        //
        // Matrix positions (global, basis A↑=0,B↑=1,A↓=2,B↓=3):
        //   h_ud block (upper-right):
        //     (0,2) = h_ud[A↑,A↓] = 0
        //     (0,3) = h_ud[A↑,B↓] = R(k)
        //     (1,2) = h_ud[B↑,A↓] = -R(-k)
        //     (1,3) = h_ud[B↑,B↓] = 0
        //   h_du block (lower-left) = conj-transpose of h_ud:
        //     (2,0) = h_du[A↓,A↑] = 0
        //     (2,1) = h_du[A↓,B↑] = conj(h_ud[B↑,A↓]) = conj(-R(-k))
        //     (3,0) = h_du[B↓,A↑] = conj(h_ud[A↑,B↓]) = conj(R(k))
        //     (3,1) = h_du[B↓,B↑] = 0

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
        h.set(0, 3, r_ab); // h_ud[A↑,B↓] = R(k)
        h.set(1, 2, s_ab); // h_ud[B↑,A↓] = S(k) = -R(-k)
                            // h_du entries (lower-left) = conj-transpose of h_ud:
        h.set(2, 1, s_ab.conj()); // h_du[A↓,B↑] = conj(S(k))
        h.set(3, 0, r_ab.conj()); // h_du[B↓,A↑] = conj(R(k))

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

    /// Compute the Z₂ invariant, dispatching on whether Rashba coupling is present.
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
    /// When Rashba is non-zero, S_z is not conserved and the two-band occupied
    /// subspace no longer decomposes by spin, so the spin-Chern method above
    /// does not apply. The Z₂ invariant is instead computed via
    /// [`z2_from_wilson_loop`](Self::z2_from_wilson_loop): a Wilson-loop /
    /// hybrid Wannier-charge-center calculation that tracks the two
    /// occupied-band Wannier phases across half the Brillouin zone and counts
    /// how many times they cross a reference line anchored at a TRIM point.
    /// This method is gauge-invariant by construction (unlike the retired
    /// TRIM-point Pfaffian approach, [`z2_from_trim_pfaffian`](Self::z2_from_trim_pfaffian),
    /// which has an inherent, unfixable gauge dependence — see its doc comment).
    ///
    /// # Returns
    ///
    /// - `1` if the system is in the QSH (topological) phase (Z₂ = 1)
    /// - `0` if trivial (Z₂ = 0)
    ///
    /// # Errors
    ///
    /// Returns `NumericalError` if the Chern number grid fails, or if the
    /// Wilson-loop calculation cannot resolve the phase confidently (e.g. too
    /// close to a phase transition for the sampling resolution used).
    pub fn z2_invariant(&self) -> Result<i32> {
        if self.lambda_r.abs() < 1e-12 {
            // Block-diagonal case: Z₂ = |C_up| mod 2
            self.z2_from_spin_chern()
        } else {
            // General case: Wilson-loop / hybrid Wannier-charge-center method.
            // (The TRIM-point Pfaffian method, `z2_from_trim_pfaffian`, has an
            // inherent local-branch-choice gauge dependence that cannot be
            // patched away — see that method's doc comment — so it is kept
            // only as a fixed-gauge cross-check, not used here.)
            self.z2_from_wilson_loop()
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
        // Parameterise: k(s,u) = s·b₁ + u·b₂, with s,u ∈ [0,1].
        let ((b1x, b1y), (b2x, b2y)) = honeycomb_reciprocal_vectors();

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

    /// Z₂ via TRIM-point Pfaffian method (general λ_R ≠ 0 case) — **retired**
    /// from the [`z2_invariant`](Self::z2_invariant) hot path in favour of
    /// [`z2_from_wilson_loop`](Self::z2_from_wilson_loop); kept only as a
    /// fixed-gauge correctness cross-check (see the test that uses it).
    ///
    /// # Why this method was replaced
    ///
    /// Even with the TRIM points and time-reversal operator corrected (both
    /// fixed below and in [`apply_time_reversal`]), the Pfaffian sewing-matrix
    /// construction has an inherent local-branch-choice gauge dependence: the
    /// sign/phase of the two occupied-band eigenvectors returned by a generic
    /// eigensolver at each TRIM point is arbitrary, and the Pfaffian δ_i
    /// depends on that choice in a way that cannot be patched away by fixing
    /// the sewing-matrix formula alone. Under an adversarial U(2) remix of
    /// the occupied eigenbasis this method's answer changes (0/8 correct in
    /// a stress test), whereas the Wilson-loop method is gauge-invariant by
    /// construction (16/16 correct under the same stress test). In a *fixed*
    /// (non-remixed) gauge — i.e. whatever basis the eigensolver happens to
    /// return — this method still gives textbook-correct answers, which is
    /// why it remains useful purely as a cross-check.
    ///
    /// Implements the Fu-Kane formula for Z₂ using time-reversal polarization.
    ///
    /// At each TRIM Λ_i, the 2×2 sewing matrix is:
    ///   m_{mn} = ⟨u_m(Λ_i) | Θ | u_n(Λ_i)⟩
    ///
    /// For a 2-band occupied subspace (n_occ = 2), the matrix m is antisymmetric
    /// and 2×2. Its Pfaffian is Pf[[0,a],[-a,0]] = a, computed analytically.
    /// The Z₂ parity is δ_i = sign(Pf(m_i) / √|det m_i|) = ±1.
    ///
    /// Z₂ = 1 if Π_i δ_i = -1 (odd number of -1's → topological).
    #[allow(dead_code)]
    fn z2_from_trim_pfaffian(&self) -> Result<i32> {
        // TRIM points of the *honeycomb* BZ: Γ, b₁/2, b₂/2, (b₁+b₂)/2.
        // (NOT the corners of a rectangular [−π,π]² cell — those do not
        // coincide with TRS-invariant momenta of this oblique lattice except
        // at Γ, which was the root of this method's TRIM-point bug.)
        let trim_points = honeycomb_trim_points();

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
    // Z₂ via Wilson loop / hybrid Wannier charge centers (general λ_R ≠ 0 case)
    // -----------------------------------------------------------------------

    /// Z₂ via the Wilson-loop / hybrid Wannier-charge-center method (general
    /// λ_R ≠ 0 case).
    ///
    /// Uses the recommended default resolution `n_s=120, n_u=60`; see
    /// [`z2_wilson_loop_at_resolution`](Self::z2_wilson_loop_at_resolution)
    /// for the resolution-parametrised core used by tests.
    ///
    /// # Method
    ///
    /// For `s ∈ [0, 1/2]` (parametrising `k(s,u) = s·b₁ + u·b₂`), the two
    /// occupied-band Wannier phases `θ₁(s), θ₂(s)` are extracted from the
    /// discretised Wilson loop threaded around the closed `u`-loop at fixed
    /// `s`. Because `s=0` and `s=1/2` are both TRS-invariant lines, the pair
    /// `{θ₁, θ₂}` is symmetric there; picking a reference angle in the larger
    /// of the two gaps between them at `s=0` and counting how many times the
    /// two (continuously-tracked) phase trajectories cross that reference as
    /// `s` sweeps `0 → 1/2` gives Z₂ = (crossings mod 2). This is the
    /// Yu-Qi-Bernevig-Fang-Zhang / Soluyanov-Vanderbilt "largest gap" method,
    /// specialised to a 2-band occupied subspace.
    ///
    /// Unlike [`z2_from_trim_pfaffian`](Self::z2_from_trim_pfaffian), this
    /// method is gauge-invariant by construction: each Wilson-loop link
    /// matrix is projected onto the unitary group via
    /// [`wilson::polar_unitary`], and the closing link deliberately *reuses*
    /// the already-computed state at `u=0` rather than re-diagonalising at
    /// `u=1` — this makes the accumulated loop product similar (in the
    /// linear-algebra sense) to itself under any per-point gauge choice, so
    /// its trace and determinant (and hence the extracted eigenphases) do
    /// not depend on that choice.
    ///
    /// # Errors
    ///
    /// Returns `NumericalError` if a self-consistency check fails: the Z₂
    /// parity is recomputed using the reference angle anchored at the
    /// *other* TRIM line (`s=1/2` instead of `s=0`), and again at doubled
    /// `s`-resolution; if either recomputation disagrees with the primary
    /// answer, this indicates the sampling resolution is insufficient to
    /// resolve the phase confidently (typically because the system is very
    /// close to a topological phase transition), and an honest error is
    /// returned rather than a possibly-wrong guess.
    fn z2_from_wilson_loop(&self) -> Result<i32> {
        self.z2_wilson_loop_at_resolution(120, 60)
    }

    /// Resolution-parametrised core of
    /// [`z2_from_wilson_loop`](Self::z2_from_wilson_loop).
    ///
    /// Exposed (privately) so tests can probe behaviour at non-default
    /// resolutions — e.g. to demonstrate the self-check erroring out instead
    /// of guessing when resolution is too low near a phase transition.
    ///
    /// # Errors
    ///
    /// See [`z2_from_wilson_loop`](Self::z2_from_wilson_loop).
    fn z2_wilson_loop_at_resolution(&self, n_s: usize, n_u: usize) -> Result<i32> {
        let (branch_a, branch_b) =
            self.wilson_branches(n_s, n_u, |_psi: &mut Vec<Vec<Complex>>| {})?;
        let last = branch_a.len() - 1;

        let z2_primary = z2_parity_from_branches(&branch_a, &branch_b, 0);
        let z2_other_trim = z2_parity_from_branches(&branch_a, &branch_b, last);
        if z2_primary != z2_other_trim {
            return Err(error::numerical_error(&format!(
                "Z2 Wilson-loop self-check failed: reference angles anchored \
                 at the two TRIM lines (s=0 giving Z2={z2_primary}, s=1/2 \
                 giving Z2={z2_other_trim}) disagree. This indicates \
                 insufficient sampling resolution (n_s={n_s}, n_u={n_u}) to \
                 resolve the phase confidently, most likely because the \
                 system is very close to a topological phase transition."
            )));
        }

        // Second self-check: doubling the s-resolution should not change the
        // answer. (n_u is left unchanged, matching the derivation's finding
        // that s-resolution — not u-resolution — is the sensitive direction
        // near a phase transition.)
        let (branch_a2, branch_b2) =
            self.wilson_branches(2 * n_s, n_u, |_psi: &mut Vec<Vec<Complex>>| {})?;
        let z2_doubled = z2_parity_from_branches(&branch_a2, &branch_b2, 0);
        if z2_doubled != z2_primary {
            return Err(error::numerical_error(&format!(
                "Z2 Wilson-loop self-check failed: doubling the s-resolution \
                 (n_s={n_s} -> {}) changes the answer (Z2={z2_primary} -> \
                 {z2_doubled}). This indicates insufficient sampling \
                 resolution to resolve the phase confidently, most likely \
                 because the system is very close to a topological phase \
                 transition.",
                2 * n_s
            )));
        }

        Ok(z2_primary)
    }

    /// Build the two continued Wannier-phase branches `θ_a(s), θ_b(s)` for
    /// `s = i/(2·n_s)`, `i = 0..=n_s`, by threading a 2-band Wilson loop
    /// around `u = j/n_u`, `j = 0..n_u`, at each fixed `s`
    /// (`k(s,u) = s·b₁ + u·b₂`).
    ///
    /// `remix` is invoked on the raw 4×2 "lowest two eigenvectors" (as two
    /// length-4 columns) immediately after diagonalising at every sampled
    /// `(s,u)` point, before it is used to build link matrices. Production
    /// code passes a no-op closure; tests use this hook to inject an
    /// adversarial per-point U(2) rotation to stress-test gauge invariance
    /// (see `z2_wilson_loop_gauge_invariant_under_degenerate_remix`).
    ///
    /// # Errors
    ///
    /// - `InvalidParameter` if `n_s < 4` or `n_u < 4`.
    /// - Propagates errors from `hermitian_eigendecomposition` or
    ///   [`wilson::polar_unitary`].
    fn wilson_branches<F>(
        &self,
        n_s: usize,
        n_u: usize,
        mut remix: F,
    ) -> Result<(Vec<f64>, Vec<f64>)>
    where
        F: FnMut(&mut Vec<Vec<Complex>>),
    {
        if n_s < 4 {
            return Err(error::invalid_param(
                "n_s",
                "Wilson-loop Z2 needs at least 4 points along s",
            ));
        }
        if n_u < 4 {
            return Err(error::invalid_param(
                "n_u",
                "Wilson-loop Z2 needs at least 4 points along u",
            ));
        }

        let (b1, b2) = honeycomb_reciprocal_vectors();

        let mut branch_a = Vec::with_capacity(n_s + 1);
        let mut branch_b = Vec::with_capacity(n_s + 1);
        let mut prev_a = 0.0_f64;
        let mut prev_b = 0.0_f64;

        for i in 0..=n_s {
            let s = i as f64 / (2.0 * n_s as f64);

            // Ψ(u_j) = lowest-2 eigenvectors of H(k(s,u_j)), j = 0..n_u-1.
            let mut psis: Vec<Vec<Vec<Complex>>> = Vec::with_capacity(n_u);
            for j in 0..n_u {
                let u = j as f64 / n_u as f64;
                let kx = s * b1.0 + u * b2.0;
                let ky = s * b1.1 + u * b2.1;
                let h = self.hamiltonian_at(kx, ky);
                let (_, vecs) = h.hermitian_eigendecomposition()?;
                let mut psi = vec![vecs.column(0), vecs.column(1)];
                remix(&mut psi);
                psis.push(psi);
            }

            // Link matrices M_j = Ψ(u_j)† Ψ(u_{j+1}), with the closing link
            // M_{n_u-1} reusing the already-computed Ψ(u_0) rather than
            // re-diagonalising at u=1 — load-bearing for gauge invariance
            // (see doc comment on `z2_from_wilson_loop`).
            let mut w = CMatrix::eye(2);
            for j in 0..n_u {
                let jp1 = (j + 1) % n_u;
                let raw = WilsonLoop::link_matrix(&psis[j], &psis[jp1]);
                let link = wilson::polar_unitary(&raw)?;
                w = w.matmul(&link)?;
            }

            let (th1, th2) = wilson::unitary_2x2_eigenphases_exact(&w);

            if i == 0 {
                branch_a.push(th1);
                branch_b.push(th2);
                prev_a = th1;
                prev_b = th2;
            } else {
                let cost_same = circular_distance(prev_a, th1) + circular_distance(prev_b, th2);
                let cost_swap = circular_distance(prev_a, th2) + circular_distance(prev_b, th1);
                let (next_a, next_b) = if cost_same <= cost_swap {
                    (th1, th2)
                } else {
                    (th2, th1)
                };
                branch_a.push(next_a);
                branch_b.push(next_b);
                prev_a = next_a;
                prev_b = next_b;
            }
        }

        Ok((branch_a, branch_b))
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

/// Apply the time-reversal operator Θ = U·K = (iσ_y ⊗ I₂)·K to a 4-component
/// state, where `K` is complex conjugation.
///
/// Θ is genuinely **antiunitary**: it is the unitary permutation `U` below
/// *composed with complex conjugation*, not `U` alone. Omitting the
/// conjugation silently turns Θ into a unitary operator with `Θ² = +1`
/// instead of the required antiunitary `Θ² = -1` (Kramers), which breaks
/// gauge-independence of anything built from it (e.g. the Pfaffian sewing
/// matrix in [`KaneMeleModel::z2_from_trim_pfaffian`]).
///
/// In the (A↑, B↑, A↓, B↓) basis, the matrix representation of `U = iσ_y ⊗ I₂` is:
/// ```text
/// [[  0,  0, +1,  0 ],
///  [  0,  0,  0, +1 ],
///  [ -1,  0,  0,  0 ],
///  [  0, -1,  0,  0 ]]
/// ```
/// so `Θ(v₀, v₁, v₂, v₃) = U·conj(v₀,v₁,v₂,v₃) = (+conj(v₂), +conj(v₃), −conj(v₀), −conj(v₁))`.
///
/// Used by [`KaneMeleModel::z2_from_trim_pfaffian`] (a retired, test-only
/// cross-check) and directly by regression tests; not on the
/// [`KaneMeleModel::z2_invariant`] hot path, hence `#[allow(dead_code)]`.
#[allow(dead_code)]
fn apply_time_reversal(u: &[Complex]) -> Vec<Complex> {
    debug_assert_eq!(u.len(), 4);
    vec![
        u[2].conj(),
        u[3].conj(),
        u[0].conj().neg(),
        u[1].conj().neg(),
    ]
}

// ---------------------------------------------------------------------------
// Wilson-loop Z₂ helpers (angle/circle bookkeeping)
// ---------------------------------------------------------------------------

/// Wraps an angle into `[-π, π)`.
#[inline]
fn wrap_to_pm_pi(x: f64) -> f64 {
    (x + PI).rem_euclid(2.0 * PI) - PI
}

/// Circular distance between two angles: the shorter arc length, in `[0, π]`.
#[inline]
fn circular_distance(a: f64, b: f64) -> f64 {
    wrap_to_pm_pi(a - b).abs()
}

/// Returns `true` if the *short* arc from `a` to `b` (whichever of the two
/// ways around the circle is shorter) passes over the marked angle `g`.
///
/// Uses a half-open convention: excludes the start point `a`, includes the
/// end point `b`. This matters because at high-symmetry parameter values a
/// sampled trajectory can land *exactly* on `g` (not merely close to it, as
/// a floating-point coincidence, but exactly, as a consequence of an exact
/// lattice symmetry) — a strict open interval `(0, d)` would silently miss
/// counting that crossing, flipping the computed parity. Landing exactly on
/// `g` is attributed to the segment ending there (not the one starting
/// there), so a sample sitting exactly at `g` is counted exactly once, not
/// zero or two times.
#[inline]
fn short_arc_crosses(a: f64, b: f64, g: f64) -> bool {
    let d = wrap_to_pm_pi(b - a);
    let gp = wrap_to_pm_pi(g - a);
    if d > 0.0 {
        gp > 0.0 && gp <= d
    } else if d < 0.0 {
        gp >= d && gp < 0.0
    } else {
        false
    }
}

/// Angular midpoint of one of the two circular gaps between angles `a` and
/// `b`: the larger gap if `use_larger` is `true`, otherwise the smaller one.
///
/// (The midpoints of the two complementary gaps are always exactly
/// antipodal — π apart — regardless of `a`, `b`; this function just picks
/// which of that antipodal pair to return.)
#[inline]
fn gap_reference_angle(a: f64, b: f64, use_larger: bool) -> f64 {
    let two_pi = 2.0 * PI;
    let a = a.rem_euclid(two_pi);
    let b = b.rem_euclid(two_pi);
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    let gap_inside = hi - lo;
    let gap_outside = two_pi - gap_inside;
    let inside_is_larger = gap_inside >= gap_outside;
    let inside_mid = (lo + hi) * 0.5;
    if inside_is_larger == use_larger {
        inside_mid.rem_euclid(two_pi)
    } else {
        (inside_mid + PI).rem_euclid(two_pi)
    }
}

/// Z₂ parity from a pair of continued Wannier-phase branches, using the
/// reference angle derived from the branches' values at `at_index` (the
/// angular midpoint of their larger circular gap).
fn z2_parity_from_branches(branch_a: &[f64], branch_b: &[f64], at_index: usize) -> i32 {
    let g = gap_reference_angle(branch_a[at_index], branch_b[at_index], true);
    let mut crossings: u32 = 0;
    for w in branch_a.windows(2) {
        if short_arc_crosses(w[0], w[1], g) {
            crossings += 1;
        }
    }
    for w in branch_b.windows(2) {
        if short_arc_crosses(w[0], w[1], g) {
            crossings += 1;
        }
    }
    (crossings % 2) as i32
}

/// Complex inner product ⟨a|b⟩ = Σ_i a_i* · b_i.
///
/// Used by [`KaneMeleModel::z2_from_trim_pfaffian`] (a retired, test-only
/// cross-check) and directly by regression tests; not on the
/// [`KaneMeleModel::z2_invariant`] hot path, hence `#[allow(dead_code)]`.
#[allow(dead_code)]
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

    // -----------------------------------------------------------------------
    // Test 14: Time-reversal symmetry holds with Rashba coupling (λ_R ≠ 0)
    // -----------------------------------------------------------------------
    #[test]
    fn time_reversal_symmetry_holds_with_rashba() {
        // Regression test for the Rashba/TRS bug: the old code assembled
        // H[B↑,A↓](k) = -conj(H[A↑,B↓](k)), which violates Θ H(k) Θ† = H(-k)
        // by up to O(0.1-0.3) at generic k whenever λ_R ≠ 0. The fix uses
        // H[B↑,A↓](k) = -H[A↑,B↓](-k) instead, restoring TRS to numerical
        // precision (~1e-16) even with Rashba coupling turned on.
        let test_points = [(0.5, 0.3), (-0.7, 0.9), (1.1, -0.4), (2.0, 0.6)];

        // Θ_mat in (A↑,B↑,A↓,B↓) basis: [[0,0,+1,0],[0,0,0,+1],[-1,0,0,0],[0,-1,0,0]]
        let theta_sign: [f64; 4] = [1.0, 1.0, -1.0, -1.0];
        let theta_map: [usize; 4] = [2, 3, 0, 1];

        for &lambda_r in &[0.02, 0.05, 0.2] {
            let model = KaneMeleModel::new(1.0, 0.08, lambda_r, 0.06).unwrap();

            for &(kx, ky) in &test_points {
                let h_k = model.hamiltonian_at(kx, ky);
                let h_mk = model.hamiltonian_at(-kx, -ky);

                // H(k)* (elementwise complex conjugate, not Hermitian conjugate)
                let mut h_conj = CMatrix::zeros(4);
                for i in 0..4 {
                    for j in 0..4 {
                        h_conj.set(i, j, h_k.get(i, j).conj());
                    }
                }

                // Θ_mat · H(k)* · Θ_mat†
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
                    "TRS violated at k=({kx},{ky}), λ_R={lambda_r}: \
                     ||Θ H(k)* Θ^T - H(-k)||={:.2e}",
                    diff.frobenius_norm()
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Test 15: Kane-Mele-Rashba critical ratio closes the gap at the Dirac point
    // -----------------------------------------------------------------------
    #[test]
    fn critical_rashba_ratio_closes_gap_at_dirac_point() {
        // Known Kane-Mele-Rashba topological phase boundary: at the critical
        // ratio λ_R = 2√3·λ_SO (with λ_v = 0) the gap at the Dirac point K
        // closes exactly — a Rashba-driven topological → trivial transition.
        // This is a physics regression test tied to Fix 1: with the old (buggy)
        // Rashba assembly, R(K) ≡ H[A↑,B↓](K) vanished identically for all
        // λ_R, so this transition could never be observed at K.
        let lambda_so = 0.1;
        let lambda_r_crit = 2.0 * 3.0_f64.sqrt() * lambda_so;
        let model = KaneMeleModel::new(1.0, lambda_so, lambda_r_crit, 0.0).unwrap();

        let k_x = 4.0 * PI / 3.0;
        let evals = model.energy_bands(k_x, 0.0).unwrap();
        let gap_at_k = evals[2] - evals[1];
        assert!(
            gap_at_k < 1e-6,
            "Gap at K should close at the critical Rashba ratio λ_R=2√3·λ_SO, got {gap_at_k:.3e}"
        );
    }

    // -----------------------------------------------------------------------
    // Shared helpers for the Wilson-loop / gauge-invariance regression tests
    // -----------------------------------------------------------------------

    /// Build a (numerically) Haar-random 2×2 unitary by projecting a random
    /// complex matrix onto U(2) via [`wilson::polar_unitary`].
    fn random_u2_matrix(rng: &mut crate::frustrated::lattice::Xorshift64) -> CMatrix {
        let mut raw = CMatrix::zeros(2);
        for i in 0..2 {
            for j in 0..2 {
                let re = 2.0 * rng.next_f64() - 1.0;
                let im = 2.0 * rng.next_f64() - 1.0;
                raw.set(i, j, Complex::new(re, im));
            }
        }
        wilson::polar_unitary(&raw).unwrap()
    }

    /// Right-multiply the 2 columns of `psi` (each a length-4 state vector)
    /// by the 2×2 unitary `v`: `psi[:,m] <- Σ_n psi[:,n]·v[n,m]`.
    fn remix_psi_columns(psi: &mut [Vec<Complex>], v: &CMatrix) {
        let old0 = psi[0].clone();
        let old1 = psi[1].clone();
        for row in 0..old0.len() {
            psi[0][row] = old0[row].mul(&v.get(0, 0)).add(&old1[row].mul(&v.get(1, 0)));
            psi[1][row] = old0[row].mul(&v.get(0, 1)).add(&old1[row].mul(&v.get(1, 1)));
        }
    }

    // -----------------------------------------------------------------------
    // Test 16: All 4 corrected honeycomb TRIM points are Kramers-degenerate
    // -----------------------------------------------------------------------
    #[test]
    fn honeycomb_trim_points_are_kramers_degenerate() {
        // At each of the 4 *corrected* honeycomb TRIM points (Γ, b₁/2, b₂/2,
        // (b₁+b₂)/2 — not the corners of a rectangular [-π,π]² cell), Kramers'
        // theorem forces every band to be at least doubly degenerate,
        // regardless of λ_R (this only needs TRS, restored for λ_R≠0 by Fix 1).
        let model = KaneMeleModel::new(1.0, 0.08, 0.05, 0.06).unwrap();
        for &(kx, ky) in &honeycomb_trim_points() {
            let evals = model.energy_bands(kx, ky).unwrap();
            assert!(
                approx(evals[0], evals[1], 1e-9),
                "TRIM ({kx:.4},{ky:.4}): Kramers pair 1 not degenerate: {} vs {}",
                evals[0],
                evals[1]
            );
            assert!(
                approx(evals[2], evals[3], 1e-9),
                "TRIM ({kx:.4},{ky:.4}): Kramers pair 2 not degenerate: {} vs {}",
                evals[2],
                evals[3]
            );
        }
    }

    // -----------------------------------------------------------------------
    // Test 17: apply_time_reversal is genuinely antiunitary (Θ²=-1)
    // -----------------------------------------------------------------------
    #[test]
    fn apply_time_reversal_is_genuinely_antiunitary() {
        // Θ² = -1 exactly (Kramers), for the fixed (conjugating) apply_time_reversal.
        let test_vectors: [[Complex; 4]; 3] = [
            [
                Complex::new(0.3, 0.1),
                Complex::new(-0.2, 0.4),
                Complex::new(0.1, -0.5),
                Complex::new(0.6, 0.2),
            ],
            [Complex::ONE, Complex::ZERO, Complex::ZERO, Complex::ZERO],
            [
                Complex::new(0.5, 0.5),
                Complex::new(-0.5, 0.5),
                Complex::new(0.5, -0.5),
                Complex::new(-0.5, -0.5),
            ],
        ];
        for v in &test_vectors {
            let theta_v = apply_time_reversal(v);
            let theta2_v = apply_time_reversal(&theta_v);
            for i in 0..4 {
                let expected = v[i].neg();
                let diff = theta2_v[i].sub(&expected);
                assert!(
                    diff.norm() < 1e-14,
                    "Theta^2 != -1 at component {i}: got {:?}, expected {:?}",
                    theta2_v[i],
                    expected
                );
            }
        }

        // Sewing matrix antisymmetry under random U(2) remixes at each TRIM
        // point: for a genuinely antiunitary Θ with Θ²=-1, the sewing matrix
        // m_mn=⟨u_m|Θ u_n⟩ built from *any* orthonormal occupied basis is
        // antisymmetric (m_nm=-m_mn) — this is what the old, non-conjugating
        // "Θ" broke.
        let model = KaneMeleModel::new(1.0, 0.08, 0.05, 0.06).unwrap();
        let mut rng = crate::frustrated::lattice::Xorshift64::new(12345).unwrap();
        for &(kx, ky) in &honeycomb_trim_points() {
            let h = model.hamiltonian_at(kx, ky);
            let (_, vecs) = h.hermitian_eigendecomposition().unwrap();
            let u0: Vec<Complex> = (0..4).map(|r| vecs.get(r, 0)).collect();
            let u1: Vec<Complex> = (0..4).map(|r| vecs.get(r, 1)).collect();

            for _ in 0..5 {
                let v = random_u2_matrix(&mut rng);
                // Remix: (u0', u1') = (u0, u1) . V
                let u0p: Vec<Complex> = (0..4)
                    .map(|r| u0[r].mul(&v.get(0, 0)).add(&u1[r].mul(&v.get(1, 0))))
                    .collect();
                let u1p: Vec<Complex> = (0..4)
                    .map(|r| u0[r].mul(&v.get(0, 1)).add(&u1[r].mul(&v.get(1, 1))))
                    .collect();

                let theta_u0p = apply_time_reversal(&u0p);
                let theta_u1p = apply_time_reversal(&u1p);
                let m00 = inner_product(&u0p, &theta_u0p);
                let m11 = inner_product(&u1p, &theta_u1p);
                let m01 = inner_product(&u0p, &theta_u1p);
                let m10 = inner_product(&u1p, &theta_u0p);

                assert!(
                    m00.norm() < 1e-9,
                    "m00 not ~0 at TRIM ({kx:.3},{ky:.3}): {m00:?}"
                );
                assert!(
                    m11.norm() < 1e-9,
                    "m11 not ~0 at TRIM ({kx:.3},{ky:.3}): {m11:?}"
                );
                let sum = m01.add(&m10);
                assert!(
                    sum.norm() < 1e-9,
                    "sewing matrix not antisymmetric at TRIM ({kx:.3},{ky:.3}): \
                     m01={m01:?} m10={m10:?}"
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Test 18: Z2 via Wilson loop matches the textbook Kane-Mele phase diagram
    // -----------------------------------------------------------------------
    #[test]
    fn z2_wilson_loop_matches_textbook_phase_diagram() {
        // Sweep λ_v across the textbook Kane-Mele phase diagram at λ_so=0.1
        // (critical λ_v = 3√3·λ_so ≈ 0.5196), with λ_R=0.05 (nonzero, so
        // z2_invariant() dispatches to the Wilson-loop branch).
        //
        // Deviation from a literal λ_v=0.0 sweep point: at exactly λ_v=0 the
        // NNN SOC structure factor g(kx=0,ky) vanishes *identically for every
        // ky* (not just at Γ) — an exact lattice symmetry of this honeycomb
        // parameterisation, independent of any bug — so the entire s=0 TRIM
        // line becomes exactly spin-degenerate whenever |λ_R| is small,
        // making the 2D occupied subspace ill-conditioned along that whole
        // line. Direct numerical experiments (varying resolution from
        // n_s=60 to n_s=4000, and λ_R from 1e-9 to 0.2) show the raw
        // crossing count is genuinely unstable (0, 1, or 2) at λ_v=0 for any
        // tested λ_R — not a resolution issue that more sampling fixes, but
        // an inherent degeneracy of this particular TRIM line. λ_v=0.02 (still
        // deep in the topological phase) is used instead, to keep this test
        // exercising generic, well-conditioned physics.
        let lambda_so = 0.1;
        let lambda_v_crit = 3.0 * 3.0_f64.sqrt() * lambda_so;
        let lambda_r = 0.05;

        let sweep = [0.02, 0.1, 0.2, 0.3, 0.4, 0.45, 0.55, 0.65, 0.8, 1.2, 2.0];

        for &lambda_v in &sweep {
            let model = KaneMeleModel::new(1.0, lambda_so, lambda_r, lambda_v).unwrap();
            let z2 = model.z2_invariant().unwrap();
            let expected = if lambda_v < lambda_v_crit { 1 } else { 0 };
            assert_eq!(
                z2, expected,
                "λ_v={lambda_v}: expected Z2={expected} (critical={lambda_v_crit:.4}), got {z2}"
            );

            // Cross-check against the already-trusted spin-Chern method at λ_R=0.
            let model_r0 = KaneMeleModel::new(1.0, lambda_so, 0.0, lambda_v).unwrap();
            let z2_ref = model_r0.z2_from_spin_chern().unwrap();
            assert_eq!(
                z2_ref, expected,
                "λ_v={lambda_v}: spin-Chern reference disagrees with textbook expectation"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Test 19: Z2 via Wilson loop is gauge-invariant under a degenerate-eigenbasis remix
    // -----------------------------------------------------------------------
    #[test]
    fn z2_wilson_loop_gauge_invariant_under_degenerate_remix() {
        use crate::frustrated::lattice::Xorshift64;

        // A smaller-than-production resolution keeps this test fast; gauge
        // invariance is an exact algebraic property (see doc comment on
        // `wilson_branches`: the closing link reuses Ψ(u_0), making the loop
        // product similar to itself under any per-point gauge choice), not a
        // resolution-dependent one, so this does not weaken the test.
        let n_s = 40;
        let n_u = 24;

        for &(lambda_v, lambda_r, label) in
            &[(0.1, 0.05, "topological"), (0.8, 0.05, "trivial")]
        {
            let model = KaneMeleModel::new(1.0, 0.1, lambda_r, lambda_v).unwrap();

            let (clean_a, clean_b) = model
                .wilson_branches(n_s, n_u, |_psi: &mut Vec<Vec<Complex>>| {})
                .unwrap();
            let z2_clean = z2_parity_from_branches(&clean_a, &clean_b, 0);

            for seed in 1..=8u64 {
                let mut rng = Xorshift64::new(seed).unwrap();
                let (remixed_a, remixed_b) = model
                    .wilson_branches(n_s, n_u, |psi: &mut Vec<Vec<Complex>>| {
                        let v = random_u2_matrix(&mut rng);
                        remix_psi_columns(psi, &v);
                    })
                    .unwrap();
                let z2_remixed = z2_parity_from_branches(&remixed_a, &remixed_b, 0);
                assert_eq!(
                    z2_remixed, z2_clean,
                    "{label} point (λ_v={lambda_v}, λ_R={lambda_r}) not gauge \
                     invariant at seed={seed}: clean Z2={z2_clean}, remixed Z2={z2_remixed}"
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Test 20: unitary_2x2_eigenphases_exact matches planted eigenvalues
    // -----------------------------------------------------------------------
    #[test]
    fn unitary_2x2_eigenphases_exact_matches_planted_eigenvalues() {
        use crate::frustrated::lattice::Xorshift64;

        let mut rng = Xorshift64::new(777).unwrap();
        for _ in 0..20 {
            let v = random_u2_matrix(&mut rng);
            let theta1 = (rng.next_f64() - 0.5) * 2.0 * PI;
            let theta2 = (rng.next_f64() - 0.5) * 2.0 * PI;

            let mut d = CMatrix::zeros(2);
            d.set(0, 0, Complex::from_polar(1.0, theta1));
            d.set(1, 1, Complex::from_polar(1.0, theta2));

            // W = V . D . V^dagger (planted eigenvalues e^{iθ1}, e^{iθ2})
            let vd = v.matmul(&d).unwrap();
            let w = vd.matmul(&v.conj_transpose()).unwrap();

            let (extracted1, extracted2) = wilson::unitary_2x2_eigenphases_exact(&w);

            let (p1, p2) = (wrap_to_pm_pi(theta1), wrap_to_pm_pi(theta2));
            let (e1, e2) = (wrap_to_pm_pi(extracted1), wrap_to_pm_pi(extracted2));

            let matches_direct = (p1 - e1).abs() < 1e-9 && (p2 - e2).abs() < 1e-9;
            let matches_swapped = (p1 - e2).abs() < 1e-9 && (p2 - e1).abs() < 1e-9;
            assert!(
                matches_direct || matches_swapped,
                "planted (θ1={theta1:.6},θ2={theta2:.6}) -> wrapped ({p1:.6},{p2:.6}), \
                 extracted ({extracted1:.6},{extracted2:.6}) -> wrapped ({e1:.6},{e2:.6})"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Test 21: Z2 via Wilson loop errors rather than guesses near a transition
    // -----------------------------------------------------------------------
    #[test]
    fn z2_wilson_loop_errors_rather_than_guesses_near_transition() {
        // Near λ_v = 3√3·λ_so at a resolution too low to resolve the
        // transition, the self-check (reference angle anchored at the other
        // TRIM line, and doubled s-resolution) should trigger an Err rather
        // than silently returning a possibly-wrong Z2.
        let lambda_so = 0.1;
        let lambda_v_crit = 3.0 * 3.0_f64.sqrt() * lambda_so;
        let lambda_v = lambda_v_crit - 0.01;
        let model = KaneMeleModel::new(1.0, lambda_so, 1e-9, lambda_v).unwrap();

        let result = model.z2_wilson_loop_at_resolution(10, 60);
        assert!(
            result.is_err(),
            "Expected an Err near the phase transition at low resolution, got {result:?}"
        );
    }
}
