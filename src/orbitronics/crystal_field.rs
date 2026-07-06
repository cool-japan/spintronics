//! Crystal-field operators for 3d transition-metal ions.
//!
//! This module provides the linear-algebra building blocks of ligand-field theory
//! for a single d-electron shell: the real cubic-harmonic d-orbital basis, the
//! orbital angular momentum operators `L_x, L_y, L_z` in that basis, the
//! electrostatic crystal-field Hamiltonian for common coordination environments,
//! and the single-particle spin-orbit Hamiltonian `H_CF + lambda*L.S`.
//!
//! ## Fidelity boundary
//!
//! **This module implements single-configuration ligand-field theory only.** A
//! single d^n configuration is split by an electrostatic crystal field into
//! one-electron orbital energies (t2g/eg for octahedral, e/t2 for tetrahedral, or
//! further split by an axial/tetragonal distortion), and spin-orbit coupling is
//! treated perturbatively as a single-particle `lambda * L.S` term acting on that
//! same one-electron orbital basis (tensored with spin). Concretely, this is:
//!
//! - Atomic Hund's-rule free-ion terms (S, L, J) — see [`super::d_orbital_moment`].
//! - A single crystal-field-split configuration (no configuration interaction
//!   between different free-ion terms of the same d^n, e.g. no ⁴T1g(⁴F)-⁴T1g(⁴P)
//!   mixing).
//! - Perturbative / single-particle spin-orbit coupling (`H_SOC = lambda*L.S` on
//!   the one-electron 10-dimensional orbital⊗spin space).
//!
//! This is **not** full many-electron multiplet / configuration-interaction
//! theory: there are no Tanabe-Sugano diagrams, no Racah `A`/`B`/`C`
//! interelectron-repulsion parameters, and no state mixing between different
//! free-ion terms of the same configuration. This level of theory correctly
//! reproduces ground-term symmetries, ground-state spin (high-spin/low-spin),
//! orbital quenching, and leading-order spin-orbit effects for the *ground*
//! configuration — normally sufficient for magnetism (effective moments,
//! qualitative g-shifts) — but it is not intended for quantitative
//! optical/spectroscopic term-energy prediction (d-d transition energies,
//! Racah-parameter fits), which requires the full multiplet treatment.
//!
//! ## Construction method (important)
//!
//! The real-basis operators `L_x, L_y, L_z` are **not hand-typed**. They are
//! derived by building the exact ladder-operator matrices in the complex
//! `|l=2,m⟩` basis (`m = -2..2`) and rotating them into the real cubic-harmonic
//! basis via a fixed unitary change of basis `U`:
//!
//! ```text
//! L_i(real) = U† L_i(complex) U
//! ```
//!
//! This makes the construction self-verifying: Hermiticity, `L² = 6·I`, and
//! `[L_x, L_y] = i L_z` are properties of the *algebra*, checked in tests against
//! the derived matrices, not against separately hand-copied numbers. The only
//! matrix that is written out explicitly is `U` itself (the standard, textbook
//! change of basis between real cubic harmonics and complex spherical
//! harmonics for `l=2`) — this is unavoidable input data, not a derived result.
//!
//! ## References
//!
//! - C. J. Ballhausen, *Introduction to Ligand Field Theory* (McGraw-Hill, 1962).
//! - J. S. Griffith, *The Theory of Transition-Metal Ions* (Cambridge, 1961).
//! - S. Sugano, Y. Tanabe, H. Kamimura, *Multiplets of Transition-Metal Ions in
//!   Crystals* (Academic Press, 1970).
//! - B. N. Figgis, M. A. Hitchman, *Ligand Field Theory and Its Applications*
//!   (Wiley-VCH, 2000).
//! - A. Abragam, B. Bleaney, *Electron Paramagnetic Resonance of Transition
//!   Ions* (Oxford, 1970).

use std::f64::consts::FRAC_1_SQRT_2;

use crate::error::Result;
use crate::math::{CMatrix, Complex};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// =============================================================================
// Real cubic-harmonic d-orbital basis
// =============================================================================

/// The five real cubic-harmonic d-orbitals, in a fixed enumeration order used
/// throughout this module for matrix rows/columns.
///
/// The first three (`Dxy`, `Dyz`, `Dzx`) form the triply-degenerate cubic set
/// (t2g in O_h, t2 in T_d); the last two (`Dx2y2`, `Dz2`) form the doubly
/// degenerate cubic set (eg in O_h, e in T_d). See [`OrbitalSymmetry`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum DOrbital {
    /// The `d_xy` orbital (lobes between the x and y axes).
    Dxy,
    /// The `d_yz` orbital (lobes between the y and z axes).
    Dyz,
    /// The `d_zx` orbital (lobes between the z and x axes), also written `d_xz`.
    Dzx,
    /// The `d_x²-y²` orbital (lobes along the x and y axes).
    Dx2y2,
    /// The `d_z²` orbital (lobes along the z axis, torus in the xy-plane).
    Dz2,
}

impl DOrbital {
    /// Number of d-orbitals (always 5, for angular momentum `l = 2`).
    pub const COUNT: usize = 5;

    /// All five d-orbitals, in the fixed matrix row/column order used by every
    /// operator/Hamiltonian builder in this module.
    #[inline]
    pub fn all() -> [DOrbital; 5] {
        [
            DOrbital::Dxy,
            DOrbital::Dyz,
            DOrbital::Dzx,
            DOrbital::Dx2y2,
            DOrbital::Dz2,
        ]
    }

    /// The fixed matrix row/column index (0..5) of this orbital.
    #[inline]
    pub fn index(self) -> usize {
        match self {
            DOrbital::Dxy => 0,
            DOrbital::Dyz => 1,
            DOrbital::Dzx => 2,
            DOrbital::Dx2y2 => 3,
            DOrbital::Dz2 => 4,
        }
    }

    /// A short human-readable label (e.g. `"xy"`, `"z2"`).
    #[inline]
    pub fn label(self) -> &'static str {
        match self {
            DOrbital::Dxy => "xy",
            DOrbital::Dyz => "yz",
            DOrbital::Dzx => "zx",
            DOrbital::Dx2y2 => "x2-y2",
            DOrbital::Dz2 => "z2",
        }
    }

    /// The cubic-symmetry set (doubly- or triply-degenerate) this orbital
    /// belongs to, independent of environment or energy ordering.
    #[inline]
    pub fn symmetry(self) -> OrbitalSymmetry {
        match self {
            DOrbital::Dxy | DOrbital::Dyz | DOrbital::Dzx => OrbitalSymmetry::T2,
            DOrbital::Dx2y2 | DOrbital::Dz2 => OrbitalSymmetry::E,
        }
    }
}

/// Cubic-group orbital symmetry label for a set of d-orbitals, independent of
/// which environment (O_h, T_d, ...) or energy ordering is in play.
///
/// The g/u subscript (`T2g`/`Eg` for O_h vs `T2`/`E` for T_d, which lacks an
/// inversion center) is a property of the point group, not of this label.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum OrbitalSymmetry {
    /// The doubly-degenerate cubic set: `d_z²`, `d_x²-y²` (Eg in O_h, E in T_d).
    E,
    /// The triply-degenerate cubic set: `d_xy`, `d_yz`, `d_zx` (T2g in O_h, T2 in T_d).
    T2,
}

/// Coordination environment determining the crystal-field splitting pattern of
/// the five d-orbitals.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum CrystalFieldEnvironment {
    /// Octahedral (six ligands on the +/-x,y,z axes): t2g (xy,yz,zx) lies
    /// `-0.4*10Dq` below, eg (z²,x²-y²) lies `+0.6*10Dq` above the barycenter.
    Octahedral,
    /// Tetrahedral (four ligands at alternate cube corners): the splitting is
    /// inverted relative to octahedral — e (z²,x²-y²) lies `-0.6*10Dq` below,
    /// t2 (xy,yz,zx) lies `+0.4*10Dq` above the barycenter.
    Tetrahedral,
    /// A tetragonal (axial, D4h) distortion superimposed on an octahedral
    /// parent, in Ballhausen's `Ds`/`Dt` parametrization. Splits t2g into
    /// b2g (xy) + eg' (xz,yz) and eg into a1g (z²) + b1g (x²-y²).
    TetragonalDistorted {
        /// Ballhausen's second-order axial crystal-field parameter (energy units
        /// consistent with `ten_dq`). Positive for a tetragonal elongation.
        ds: f64,
        /// Ballhausen's fourth-order axial crystal-field parameter (energy units
        /// consistent with `ten_dq`). Positive for a tetragonal elongation.
        dt: f64,
    },
}

// =============================================================================
// Complex |l=2, m⟩ basis: exact ladder-operator construction
// =============================================================================

/// The five `m` quantum numbers of the `|l=2,m⟩` basis, in ascending
/// row/column order (index `k` <-> `m = k - 2`).
const M_VALUES: [i32; 5] = [-2, -1, 0, 1, 2];

/// `L_z` in the complex `|l=2,m⟩` basis: exactly diagonal with eigenvalues
/// `m = -2, -1, 0, 1, 2` (units of ħ).
fn l_z_complex() -> CMatrix {
    CMatrix::from_diagonal(&[-2.0, -1.0, 0.0, 1.0, 2.0])
}

/// `L_+` in the complex `|l=2,m⟩` basis via the exact ladder-operator formula
/// `L_+|l,m⟩ = sqrt(l(l+1) - m(m+1)) |l,m+1⟩` (here `l=2`, so `l(l+1)=6`).
fn l_plus_complex() -> Result<CMatrix> {
    let mut rows = vec![vec![Complex::ZERO; 5]; 5];
    for (col, &m) in M_VALUES.iter().enumerate().take(4) {
        let m = f64::from(m);
        let coeff = (6.0 - m * (m + 1.0)).sqrt();
        rows[col + 1][col] = Complex::from_real(coeff);
    }
    CMatrix::from_rows(rows)
}

/// `L_- = (L_+)^†` in the complex basis (exact ladder-operator conjugate;
/// `L_+` has purely real entries, so this is also its plain transpose).
fn l_minus_complex() -> Result<CMatrix> {
    Ok(l_plus_complex()?.conj_transpose())
}

/// `L_x = (L_+ + L_-)/2` in the complex `|l=2,m⟩` basis.
fn l_x_complex() -> Result<CMatrix> {
    let lp = l_plus_complex()?;
    let lm = l_minus_complex()?;
    Ok(lp.add(&lm)?.scale_real(0.5))
}

/// `L_y = (L_+ - L_-)/(2i)` in the complex `|l=2,m⟩` basis.
fn l_y_complex() -> Result<CMatrix> {
    let lp = l_plus_complex()?;
    let lm = l_minus_complex()?;
    let diff = lp.sub(&lm)?;
    // 1/(2i) = -i/2
    Ok(diff.scale(Complex::new(0.0, -0.5)))
}

/// Unitary change of basis `U` from the real cubic-harmonic basis (columns, in
/// [`DOrbital::all`] order) to the complex `|l=2,m⟩` basis (rows, `m=-2..2`
/// ascending): `U[m, k] = ⟨2,m | d_k⟩`.
///
/// Standard combinations (e.g. Griffith 1961, Ch. 2):
///
/// ```text
/// |xy⟩    = (i/√2)(|2,-2⟩ - |2,2⟩)      |x²-y²⟩ = (1/√2)(|2,-2⟩ + |2,2⟩)
/// |xz⟩    = (1/√2)(|2,-1⟩ - |2,1⟩)      |yz⟩    = (i/√2)(|2,-1⟩ + |2,1⟩)
/// |z²⟩    = |2,0⟩
/// ```
///
/// Each real orbital is a unit-norm combination confined to a single `|m|`
/// subspace, and distinct orbitals live in disjoint or orthogonal subspaces,
/// so `U` is unitary by construction (verified numerically in tests).
fn u_transform() -> Result<CMatrix> {
    let s = Complex::from_real(FRAC_1_SQRT_2);
    let is = Complex::new(0.0, FRAC_1_SQRT_2);
    let z = Complex::ZERO;
    // Rows: m = -2, -1, 0, 1, 2. Columns: Dxy, Dyz, Dzx, Dx2y2, Dz2.
    let rows = vec![
        vec![is, z, z, s, z],           // m = -2
        vec![z, is, s, z, z],           // m = -1
        vec![z, z, z, z, Complex::ONE], // m =  0
        vec![z, is, s.neg(), z, z],     // m =  1
        vec![is.neg(), z, z, s, z],     // m =  2
    ];
    CMatrix::from_rows(rows)
}

// =============================================================================
// Real-basis orbital angular momentum operators (public, self-verifying)
// =============================================================================

/// The orbital angular momentum operators `(L_x, L_y, L_z)` in the real
/// cubic-harmonic d-orbital basis ([`DOrbital::all`] order), obtained by
/// rotating the exact complex `|l=2,m⟩`-basis ladder-operator matrices:
/// `L_i(real) = U† L_i(complex) U`.
///
/// This is the single source of truth for every `L` operator used elsewhere in
/// this module (the SOC Hamiltonian, orbital quenching diagnostics): the real
/// matrices are never hand-typed, only derived by this basis rotation.
pub fn orbital_angular_momentum_operators() -> Result<(CMatrix, CMatrix, CMatrix)> {
    let u = u_transform()?;
    let u_dag = u.conj_transpose();
    let lx = u_dag.matmul(&l_x_complex()?)?.matmul(&u)?;
    let ly = u_dag.matmul(&l_y_complex()?)?.matmul(&u)?;
    let lz = u_dag.matmul(&l_z_complex())?.matmul(&u)?;
    Ok((lx, ly, lz))
}

/// The orbital angular momentum squared `L² = L_x² + L_y² + L_z²` in the real
/// cubic-harmonic basis. For `l=2` this must equal `l(l+1)·I = 6·I` regardless
/// of basis (a basis-independent Casimir invariant); verified in tests.
pub fn orbital_l_squared() -> Result<CMatrix> {
    let (lx, ly, lz) = orbital_angular_momentum_operators()?;
    let lx2 = lx.matmul(&lx)?;
    let ly2 = ly.matmul(&ly)?;
    let lz2 = lz.matmul(&lz)?;
    lx2.add(&ly2)?.add(&lz2)
}

// =============================================================================
// Crystal-field Hamiltonian
// =============================================================================

/// The one-electron electrostatic crystal-field Hamiltonian, diagonal in the
/// real cubic-harmonic basis ([`DOrbital::all`] order), for the given
/// `environment` and cubic splitting parameter `ten_dq` (the total O_h/T_d gap,
/// "10Dq", in whatever consistent energy unit the caller chooses).
///
/// The trace is always zero (barycenter rule): the crystal field redistributes,
/// but does not shift, the mean one-electron d-orbital energy.
pub fn crystal_field_hamiltonian(
    environment: CrystalFieldEnvironment,
    ten_dq: f64,
) -> Result<CMatrix> {
    let mut energies = [0.0_f64; 5];
    match environment {
        CrystalFieldEnvironment::Octahedral => {
            for orbital in DOrbital::all() {
                energies[orbital.index()] = match orbital.symmetry() {
                    OrbitalSymmetry::T2 => -0.4 * ten_dq,
                    OrbitalSymmetry::E => 0.6 * ten_dq,
                };
            }
        },
        CrystalFieldEnvironment::Tetrahedral => {
            for orbital in DOrbital::all() {
                energies[orbital.index()] = match orbital.symmetry() {
                    OrbitalSymmetry::T2 => 0.4 * ten_dq,
                    OrbitalSymmetry::E => -0.6 * ten_dq,
                };
            }
        },
        CrystalFieldEnvironment::TetragonalDistorted { ds, dt } => {
            energies[DOrbital::Dz2.index()] = 0.6 * ten_dq - 2.0 * ds - 6.0 * dt;
            energies[DOrbital::Dx2y2.index()] = 0.6 * ten_dq + 2.0 * ds - dt;
            energies[DOrbital::Dxy.index()] = -0.4 * ten_dq + 2.0 * ds - dt;
            energies[DOrbital::Dyz.index()] = -0.4 * ten_dq - ds + 4.0 * dt;
            energies[DOrbital::Dzx.index()] = -0.4 * ten_dq - ds + 4.0 * dt;
        },
    }
    Ok(CMatrix::from_diagonal(&energies))
}

// =============================================================================
// Spin-1/2 Pauli matrices and Kronecker product (private linear-algebra glue)
// =============================================================================

pub(crate) fn pauli_x() -> Result<CMatrix> {
    CMatrix::from_rows(vec![
        vec![Complex::ZERO, Complex::ONE],
        vec![Complex::ONE, Complex::ZERO],
    ])
}

pub(crate) fn pauli_y() -> Result<CMatrix> {
    CMatrix::from_rows(vec![
        vec![Complex::ZERO, Complex::new(0.0, -1.0)],
        vec![Complex::new(0.0, 1.0), Complex::ZERO],
    ])
}

pub(crate) fn pauli_z() -> Result<CMatrix> {
    CMatrix::from_rows(vec![
        vec![Complex::ONE, Complex::ZERO],
        vec![Complex::ZERO, Complex::new(-1.0, 0.0)],
    ])
}

/// Kronecker (tensor) product `A ⊗ B`.
pub(crate) fn kron(a: &CMatrix, b: &CMatrix) -> Result<CMatrix> {
    let na = a.n();
    let nb = b.n();
    let n = na * nb;
    let mut rows = vec![vec![Complex::ZERO; n]; n];
    for i in 0..na {
        for j in 0..na {
            let aij = a.get(i, j);
            if aij.re == 0.0 && aij.im == 0.0 {
                continue;
            }
            for p in 0..nb {
                for q in 0..nb {
                    rows[i * nb + p][j * nb + q] = aij.mul(&b.get(p, q));
                }
            }
        }
    }
    CMatrix::from_rows(rows)
}

// =============================================================================
// Single-particle spin-orbit Hamiltonian
// =============================================================================

/// The single-particle crystal-field-plus-spin-orbit Hamiltonian on the
/// 10-dimensional orbital(5)⊗spin(2) space:
///
/// ```text
/// H = H_CF ⊗ I_2  +  lambda * L.S.
///   = H_CF ⊗ I_2  +  lambda * (L_x⊗sigma_x + L_y⊗sigma_y + L_z⊗sigma_z) / 2
/// ```
///
/// where `S_i = sigma_i/2` are the spin-1/2 operators. `soc_lambda` is the
/// single-electron spin-orbit coupling constant (same energy unit as `ten_dq`).
///
/// This is a genuinely single-particle (one electron or, by particle-hole
/// symmetry, one hole) Hamiltonian — see the module-level fidelity boundary.
pub fn soc_hamiltonian(
    environment: CrystalFieldEnvironment,
    ten_dq: f64,
    soc_lambda: f64,
) -> Result<CMatrix> {
    let h_cf = crystal_field_hamiltonian(environment, ten_dq)?;
    let identity2 = CMatrix::eye(2);
    let h_cf_ext = kron(&h_cf, &identity2)?;

    let (lx, ly, lz) = orbital_angular_momentum_operators()?;
    let sx = pauli_x()?;
    let sy = pauli_y()?;
    let sz = pauli_z()?;

    let lx_sx = kron(&lx, &sx)?;
    let ly_sy = kron(&ly, &sy)?;
    let lz_sz = kron(&lz, &sz)?;

    // S_i = sigma_i / 2, so the overall factor is soc_lambda / 2.
    let l_dot_s = lx_sx.add(&ly_sy)?.add(&lz_sz)?.scale_real(soc_lambda * 0.5);

    h_cf_ext.add(&l_dot_s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hermitian_defect(m: &CMatrix) -> f64 {
        let diff = m.sub(&m.conj_transpose()).expect("same dimension");
        diff.frobenius_norm()
    }

    #[test]
    fn test_dorbital_index_roundtrip() {
        for (expected_idx, orbital) in DOrbital::all().into_iter().enumerate() {
            assert_eq!(orbital.index(), expected_idx);
        }
    }

    #[test]
    fn test_dorbital_symmetry_grouping() {
        assert_eq!(DOrbital::Dxy.symmetry(), OrbitalSymmetry::T2);
        assert_eq!(DOrbital::Dyz.symmetry(), OrbitalSymmetry::T2);
        assert_eq!(DOrbital::Dzx.symmetry(), OrbitalSymmetry::T2);
        assert_eq!(DOrbital::Dx2y2.symmetry(), OrbitalSymmetry::E);
        assert_eq!(DOrbital::Dz2.symmetry(), OrbitalSymmetry::E);
    }

    #[test]
    fn test_u_transform_is_unitary() {
        let u = u_transform().expect("u_transform builds");
        let u_dag = u.conj_transpose();
        let should_be_identity = u_dag.matmul(&u).expect("5x5 matmul");
        let identity = CMatrix::eye(5);
        let defect = should_be_identity
            .sub(&identity)
            .expect("same dim")
            .frobenius_norm();
        assert!(defect < 1e-12, "U is not unitary: defect = {}", defect);
    }

    #[test]
    fn test_operators_are_hermitian() {
        let (lx, ly, lz) = orbital_angular_momentum_operators().expect("operators build");
        assert!(hermitian_defect(&lx) < 1e-12, "L_x not Hermitian");
        assert!(hermitian_defect(&ly) < 1e-12, "L_y not Hermitian");
        assert!(hermitian_defect(&lz) < 1e-12, "L_z not Hermitian");
    }

    #[test]
    fn test_commutator_lx_ly_is_i_lz() {
        let (lx, ly, lz) = orbital_angular_momentum_operators().expect("operators build");
        let lxly = lx.matmul(&ly).expect("matmul");
        let lylx = ly.matmul(&lx).expect("matmul");
        let commutator = lxly.sub(&lylx).expect("same dim");
        let i_lz = lz.scale(Complex::I);
        let defect = commutator.sub(&i_lz).expect("same dim").frobenius_norm();
        assert!(defect < 1e-10, "[L_x,L_y] != i L_z: defect = {}", defect);
    }

    #[test]
    fn test_commutator_ly_lz_is_i_lx() {
        let (lx, ly, lz) = orbital_angular_momentum_operators().expect("operators build");
        let lylz = ly.matmul(&lz).expect("matmul");
        let lzly = lz.matmul(&ly).expect("matmul");
        let commutator = lylz.sub(&lzly).expect("same dim");
        let i_lx = lx.scale(Complex::I);
        let defect = commutator.sub(&i_lx).expect("same dim").frobenius_norm();
        assert!(defect < 1e-10, "[L_y,L_z] != i L_x: defect = {}", defect);
    }

    #[test]
    fn test_commutator_lz_lx_is_i_ly() {
        let (lx, ly, lz) = orbital_angular_momentum_operators().expect("operators build");
        let lzlx = lz.matmul(&lx).expect("matmul");
        let lxlz = lx.matmul(&lz).expect("matmul");
        let commutator = lzlx.sub(&lxlz).expect("same dim");
        let i_ly = ly.scale(Complex::I);
        let defect = commutator.sub(&i_ly).expect("same dim").frobenius_norm();
        assert!(defect < 1e-10, "[L_z,L_x] != i L_y: defect = {}", defect);
    }

    #[test]
    fn test_l_squared_is_six_identity() {
        let l2 = orbital_l_squared().expect("l_squared builds");
        let six_identity = CMatrix::eye(5).scale_real(6.0);
        let defect = l2.sub(&six_identity).expect("same dim").frobenius_norm();
        assert!(defect < 1e-9, "L^2 != 6I: defect = {}", defect);
    }

    #[test]
    fn test_l_z_eigenvalue_spectrum() {
        let (_, _, lz) = orbital_angular_momentum_operators().expect("operators build");
        let (mut eigenvalues, _) = lz.hermitian_eigendecomposition().expect("Hermitian eig");
        eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let expected = [-2.0, -1.0, 0.0, 1.0, 2.0];
        for (got, want) in eigenvalues.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-9, "eigenvalue {} != {}", got, want);
        }
    }

    #[test]
    fn test_real_orbitals_have_zero_diagonal_angular_momentum() {
        // General theorem: any real-valued spatial orbital has <psi|L_i|psi> = 0
        // for all i = x,y,z (this is the origin of orbital quenching).
        let (lx, ly, lz) = orbital_angular_momentum_operators().expect("operators build");
        for orbital in DOrbital::all() {
            let idx = orbital.index();
            assert!(
                lx.get(idx, idx).norm() < 1e-10,
                "Lx diagonal nonzero for {:?}",
                orbital
            );
            assert!(
                ly.get(idx, idx).norm() < 1e-10,
                "Ly diagonal nonzero for {:?}",
                orbital
            );
            assert!(
                lz.get(idx, idx).norm() < 1e-10,
                "Lz diagonal nonzero for {:?}",
                orbital
            );
        }
    }

    #[test]
    fn test_l_z_hand_derived_matrix_elements() {
        // Hand-derived fixtures (independent of the U-rotation code path):
        // <xy|Lz|x2-y2> = 2i, <yz|Lz|xz> = i (Griffith 1961-consistent values).
        let (_, _, lz) = orbital_angular_momentum_operators().expect("operators build");
        let xy = DOrbital::Dxy.index();
        let x2y2 = DOrbital::Dx2y2.index();
        let yz = DOrbital::Dyz.index();
        let zx = DOrbital::Dzx.index();

        let element_xy_x2y2 = lz.get(xy, x2y2);
        assert!((element_xy_x2y2.re).abs() < 1e-9);
        assert!(
            (element_xy_x2y2.im - 2.0).abs() < 1e-9,
            "got {:?}",
            element_xy_x2y2
        );

        let element_yz_zx = lz.get(yz, zx);
        assert!((element_yz_zx.re).abs() < 1e-9);
        assert!(
            (element_yz_zx.im - 1.0).abs() < 1e-9,
            "got {:?}",
            element_yz_zx
        );
    }

    #[test]
    fn test_octahedral_splitting_pattern() {
        let ten_dq = 2.0;
        let h =
            crystal_field_hamiltonian(CrystalFieldEnvironment::Octahedral, ten_dq).expect("builds");
        for orbital in DOrbital::all() {
            let e = h.get(orbital.index(), orbital.index()).re;
            let expected = match orbital.symmetry() {
                OrbitalSymmetry::T2 => -0.4 * ten_dq,
                OrbitalSymmetry::E => 0.6 * ten_dq,
            };
            assert!(
                (e - expected).abs() < 1e-12,
                "Oh energy mismatch for {:?}",
                orbital
            );
        }
        assert!(
            h.trace().re.abs() < 1e-12,
            "Oh crystal field must be traceless"
        );
    }

    #[test]
    fn test_tetrahedral_splitting_is_inverted() {
        let ten_dq = 2.0;
        let h_oh =
            crystal_field_hamiltonian(CrystalFieldEnvironment::Octahedral, ten_dq).expect("builds");
        let h_td = crystal_field_hamiltonian(CrystalFieldEnvironment::Tetrahedral, ten_dq)
            .expect("builds");
        for orbital in DOrbital::all() {
            let idx = orbital.index();
            let e_oh = h_oh.get(idx, idx).re;
            let e_td = h_td.get(idx, idx).re;
            assert!(
                (e_td + e_oh).abs() < 1e-12,
                "Td energy should be -1 * Oh energy for {:?}",
                orbital
            );
        }
        assert!(
            h_td.trace().re.abs() < 1e-12,
            "Td crystal field must be traceless"
        );
    }

    #[test]
    fn test_tetragonal_reduces_to_octahedral_at_zero_distortion() {
        let ten_dq = 2.0;
        let h_oh =
            crystal_field_hamiltonian(CrystalFieldEnvironment::Octahedral, ten_dq).expect("builds");
        let h_tet = crystal_field_hamiltonian(
            CrystalFieldEnvironment::TetragonalDistorted { ds: 0.0, dt: 0.0 },
            ten_dq,
        )
        .expect("builds");
        let defect = h_oh.sub(&h_tet).expect("same dim").frobenius_norm();
        assert!(
            defect < 1e-12,
            "tetragonal(ds=dt=0) should equal octahedral"
        );
    }

    #[test]
    fn test_tetragonal_is_traceless() {
        let h = crystal_field_hamiltonian(
            CrystalFieldEnvironment::TetragonalDistorted { ds: 0.3, dt: 0.05 },
            2.0,
        )
        .expect("builds");
        assert!(
            h.trace().re.abs() < 1e-10,
            "tetragonal crystal field must be traceless"
        );
    }

    #[test]
    fn test_tetragonal_splits_xy_from_xz_yz_degeneracy() {
        // A generic tetragonal distortion lifts the t2g degeneracy except for
        // the xz/yz pair, which remains exactly degenerate by C4v symmetry.
        let h = crystal_field_hamiltonian(
            CrystalFieldEnvironment::TetragonalDistorted { ds: 0.3, dt: 0.05 },
            2.0,
        )
        .expect("builds");
        let e_xy = h.get(DOrbital::Dxy.index(), DOrbital::Dxy.index()).re;
        let e_yz = h.get(DOrbital::Dyz.index(), DOrbital::Dyz.index()).re;
        let e_zx = h.get(DOrbital::Dzx.index(), DOrbital::Dzx.index()).re;
        assert!((e_yz - e_zx).abs() < 1e-12, "xz/yz must stay degenerate");
        assert!(
            (e_xy - e_yz).abs() > 1e-6,
            "xy should split away from xz/yz"
        );
    }

    #[test]
    fn test_soc_hamiltonian_is_hermitian_and_10_dimensional() {
        let h = soc_hamiltonian(CrystalFieldEnvironment::Octahedral, 2.0, 0.05).expect("builds");
        assert_eq!(h.n(), 10);
        assert!(hermitian_defect(&h) < 1e-10);
    }

    #[test]
    fn test_soc_hamiltonian_zero_lambda_is_pure_crystal_field() {
        let h = soc_hamiltonian(CrystalFieldEnvironment::Octahedral, 2.0, 0.0).expect("builds");
        let h_cf =
            crystal_field_hamiltonian(CrystalFieldEnvironment::Octahedral, 2.0).expect("builds");
        let h_cf_ext = kron(&h_cf, &CMatrix::eye(2)).expect("kron");
        let defect = h.sub(&h_cf_ext).expect("same dim").frobenius_norm();
        assert!(
            defect < 1e-10,
            "lambda=0 SOC Hamiltonian should reduce to H_CF ⊗ I2"
        );
    }

    #[test]
    fn test_kron_dimensions_and_identity() {
        let a = CMatrix::eye(2);
        let b = CMatrix::eye(3);
        let prod = kron(&a, &b).expect("kron builds");
        assert_eq!(prod.n(), 6);
        let identity6 = CMatrix::eye(6);
        let defect = prod.sub(&identity6).expect("same dim").frobenius_norm();
        assert!(defect < 1e-12, "I2 ⊗ I3 should equal I6");
    }
}
