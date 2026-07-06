//! Local d-orbital magnetic moments: Hund's-rule free-ion terms, crystal-field
//! high-spin/low-spin states, orbital quenching, and effective moments.
//!
//! This module builds on the operators in [`super::crystal_field`] to model the
//! local magnetic moment of a 3d transition-metal ion in a ligand field:
//! atomic Hund's-rule free-ion terms ([`free_ion_term`]), the crystal-field-
//! and-spin-orbit-aware [`CrystalFieldModel`] (high-spin/low-spin, ground-term
//! symmetry, orbital quenching, effective moments), and ~12 preset ions with
//! textbook-matched expected values.
//!
//! ## Fidelity boundary
//!
//! Same single-configuration ligand-field scope as [`super::crystal_field`]:
//! atomic Hund's-rule free-ion terms, a single crystal-field-split
//! configuration (Aufbau + Hund's rules over the one-electron levels, no
//! configuration interaction between different free-ion terms of the same
//! d^n), and perturbative single-particle spin-orbit coupling. **Not** full
//! many-electron multiplet theory: no Tanabe-Sugano diagrams, no Racah `A`,
//! `B`, `C` parameters. `pairing_energy` is a single phenomenological
//! mean-pairing-energy scale per ion (as in introductory ligand-field
//! treatments), not derived from Racah parameters.
//!
//! ## Two distinct notions of spin
//!
//! - [`free_ion_term`] / [`CrystalFieldModel::free_ion_term`]: the *free-ion*
//!   (atomic, Hund's-rule) ground term (S, L, J) — a function of electron
//!   count `n` only, independent of environment. This is Hund's rule 1
//!   (maximum multiplicity) applied with no crystal field at all.
//! - [`CrystalFieldModel::ground_state_spin`]: the *actual* ligand-field
//!   ground-state spin, which equals the free-ion S for high-spin ions but is
//!   *reduced* for low-spin ions (crystal-field-driven pairing beyond Hund's
//!   rule 1, e.g. low-spin d6 Fe(II)/Co(III) are diamagnetic, S=0, while the
//!   free-ion d6 term has S=2).
//!
//! `effective_moment_spin_only` uses the *actual* (environment-aware) spin,
//! matching real 3d-ion magnetochemistry; `effective_moment_lande` uses the
//! free-ion (S, L, J) via the Landé formula, the appropriate limit when
//! spin-orbit coupling dominates the crystal field (as for 4f rare earths, not
//! typically 3d ions — included here for contrast/comparison, and because it
//! cleanly exposes the J=0 free-ion singularity for d4/d6 high-spin ions).

use crate::error::{Error, Result};
use crate::math::{CMatrix, Complex};
use crate::orbitronics::crystal_field::{
    crystal_field_hamiltonian, kron, orbital_angular_momentum_operators, pauli_x, pauli_y, pauli_z,
    soc_hamiltonian, CrystalFieldEnvironment, DOrbital,
};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// `m_l` values in descending order (largest first), used to construct the
/// maximum-`M_L` (Hund's rule 2) state of the maximum-spin manifold.
const M_L_DESCENDING: [i32; 5] = [2, 1, 0, -1, -2];

// =============================================================================
// Free-ion (atomic) Hund's-rule terms
// =============================================================================

/// A free-ion (atomic) Russell-Saunders term: total spin `S`, total orbital
/// angular momentum `L` (always integral for a d-shell, stored as `f64` for
/// uniformity), and total angular momentum `J`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct FreeIonTerm {
    /// Total spin angular momentum quantum number.
    pub s: f64,
    /// Total orbital angular momentum quantum number (integral for a d-shell).
    pub l: f64,
    /// Total angular momentum quantum number (Hund's rule 3).
    pub j: f64,
}

/// The free-ion (atomic) Hund's-rule ground term for a d^n configuration
/// (`n` = 0..=10 electrons), independent of any crystal-field environment.
///
/// Applies Hund's three rules directly:
///
/// 1. **Maximize S**: with `n <= 5` electrons all can occupy distinct `m_l`
///    slots with parallel spin; for `n > 5`, electron-hole symmetry maps the
///    problem onto `10-n` holes.
/// 2. **Maximize L**: among maximum-S arrangements, the maximally "stretched"
///    state (largest `sum(m_l)`) is obtained by greedily occupying the
///    largest-`|m_l|` slots first; this maximally-stretched state is the
///    unique highest-weight state of the `(S_max, L_max)` term.
/// 3. **J = |L-S|** for a less-than-half-filled shell (`n <= 5`), **J = L+S**
///    for a more-than-half-filled shell (`n > 5`) (both formulas agree
///    trivially at `n=5`, where `L=0`).
///
/// # Errors
///
/// Returns [`Error::InvalidParameter`] if `n > 10` (a d-shell holds at most 10
/// electrons).
pub fn free_ion_term(n: u8) -> Result<FreeIonTerm> {
    if n > 10 {
        return Err(Error::InvalidParameter {
            param: "n".to_string(),
            reason: "a d-shell holds at most 10 electrons".to_string(),
        });
    }
    let n_eff = n.min(10 - n);
    let s = f64::from(n_eff) / 2.0;

    let l_sum: i32 = M_L_DESCENDING.iter().take(n_eff as usize).sum();
    let l = f64::from(l_sum);

    let j = if n <= 5 { (l - s).abs() } else { l + s };

    Ok(FreeIonTerm { s, l, j })
}

/// Ground-configuration orbital symmetry label under a crystal field: whether
/// the many-electron ground configuration (built by Hund's-rule/Aufbau filling
/// of the crystal-field-split one-electron levels) carries surviving
/// first-order orbital angular momentum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum GroundTermSymmetry {
    /// Orbitally non-degenerate ground configuration: `L` is quenched to
    /// first order (e.g. `t2g³` or `t2g⁶eg²`).
    A,
    /// Doubly-degenerate ground configuration: `L` is nonetheless quenched to
    /// first order — the cubic `E` representation does not contain the
    /// `T1`-like `L` operator, so `L` has no matrix elements within a pure `E`
    /// term (e.g. `t2g⁶eg³`, Cu²⁺).
    E,
    /// Triply-degenerate ground configuration: `L` survives to first order
    /// (unquenched), e.g. `t2g¹` (Ti³⁺) or `t2g⁵eg²` (high-spin Fe²⁺).
    T,
}

impl GroundTermSymmetry {
    /// Whether this ground-term symmetry has first-order-quenched orbital
    /// angular momentum (`A` and `E`), as opposed to unquenched (`T`).
    #[inline]
    pub fn is_quenched(self) -> bool {
        matches!(self, GroundTermSymmetry::A | GroundTermSymmetry::E)
    }
}

// =============================================================================
// Aufbau-with-pairing filling engine (private)
// =============================================================================

/// Occupation of one crystal-field level (the `members` orbitals, all at a
/// common energy) with `electrons` electrons, filled via Hund's rule 1 within
/// the level (spread out over all member orbitals before pairing).
#[derive(Debug, Clone)]
struct LevelOccupation {
    members: Vec<DOrbital>,
    electrons: usize,
}

impl LevelOccupation {
    /// Degeneracy (number of member orbitals) of this level.
    fn degeneracy(&self) -> usize {
        self.members.len()
    }

    /// Number of distinct ways to arrange `electrons` electrons among the
    /// member orbitals under Hund's rule 1: `C(degeneracy, k)` with
    /// `k = min(electrons, 2*degeneracy - electrons)`. Equals 1 (non-degenerate)
    /// exactly when the level is empty, exactly half-filled, or full.
    fn combinatorial_degeneracy(&self) -> u64 {
        let degeneracy = self.degeneracy();
        let k = self.electrons.min(2 * degeneracy - self.electrons);
        binomial_coefficient(degeneracy as u64, k as u64)
    }

    /// Number of unpaired (parallel-spin) electrons in this level under
    /// Hund's rule 1.
    fn unpaired_electrons(&self) -> usize {
        let degeneracy = self.degeneracy();
        if self.electrons <= degeneracy {
            self.electrons
        } else {
            2 * degeneracy - self.electrons
        }
    }
}

/// Binomial coefficient `C(n, k)` for small non-negative integers.
fn binomial_coefficient(n: u64, k: u64) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result: u64 = 1;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

/// Fills `n` electrons into `levels` (ascending energy, `(energy, member
/// orbitals)` groups) via the standard Aufbau-with-pairing algorithm: each
/// electron goes into the cheapest available slot, where occupying a new
/// (empty) orbital in a level costs that level's energy, and pairing an
/// already-singly-occupied orbital costs the level's energy plus
/// `pairing_energy`. Because the marginal cost of the `k`-th electron in a
/// level is a non-decreasing step function of `k`, this greedy strategy is
/// optimal (a standard exchange-argument result) and reproduces the textbook
/// high-spin/low-spin crossover (10Dq vs pairing energy) for any number of
/// crystal-field levels.
fn fill_levels(
    levels: &[(f64, Vec<DOrbital>)],
    n: u8,
    pairing_energy: f64,
) -> Vec<LevelOccupation> {
    let mut occupations: Vec<LevelOccupation> = levels
        .iter()
        .map(|(_, members)| LevelOccupation {
            members: members.clone(),
            electrons: 0,
        })
        .collect();

    for _ in 0..n {
        let mut best_idx: Option<usize> = None;
        let mut best_cost = f64::INFINITY;
        for (idx, level) in occupations.iter().enumerate() {
            let degeneracy = level.degeneracy();
            if level.electrons >= 2 * degeneracy {
                continue;
            }
            let level_energy = levels[idx].0;
            let cost = if level.electrons < degeneracy {
                level_energy
            } else {
                level_energy + pairing_energy
            };
            if cost < best_cost {
                best_cost = cost;
                best_idx = Some(idx);
            }
        }
        if let Some(idx) = best_idx {
            occupations[idx].electrons += 1;
        }
    }
    occupations
}

/// Extracts the principal sub-block of a 5x5 orbital operator restricted to
/// `members` (in the given order), yielding a `members.len()`-square matrix.
fn orbital_sub_block(op: &CMatrix, members: &[DOrbital]) -> Result<CMatrix> {
    let dim = members.len();
    let mut rows = vec![vec![Complex::ZERO; dim]; dim];
    for (a, &oa) in members.iter().enumerate() {
        for (b, &ob) in members.iter().enumerate() {
            rows[a][b] = op.get(oa.index(), ob.index());
        }
    }
    CMatrix::from_rows(rows)
}

// =============================================================================
// CrystalFieldModel
// =============================================================================

/// A d^n transition-metal ion in a specific crystal-field environment, with
/// spin-orbit coupling, at the single-configuration ligand-field level of
/// theory (see the module-level fidelity boundary).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CrystalFieldModel {
    /// Number of d electrons (1..=9).
    pub n: u8,
    /// Coordination environment (determines the crystal-field splitting).
    pub environment: CrystalFieldEnvironment,
    /// Cubic crystal-field splitting parameter "10Dq" (energy units consistent
    /// with `pairing_energy` and `soc_lambda`, e.g. eV).
    pub ten_dq: f64,
    /// Phenomenological mean spin-pairing energy (same energy units as
    /// `ten_dq`); compared against `ten_dq` to determine high-spin/low-spin.
    pub pairing_energy: f64,
    /// Single-electron spin-orbit coupling constant `lambda` (same energy
    /// units as `ten_dq`), entering `H_SOC = lambda * L.S`.
    pub soc_lambda: f64,
}

impl CrystalFieldModel {
    /// Construct a new crystal-field model, validating physical parameters.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] if `n` is not in `1..=9`, if
    /// `ten_dq` is negative, or if `pairing_energy` is not positive.
    pub fn new(
        n: u8,
        environment: CrystalFieldEnvironment,
        ten_dq: f64,
        pairing_energy: f64,
        soc_lambda: f64,
    ) -> Result<Self> {
        if !(1..=9).contains(&n) {
            return Err(Error::InvalidParameter {
                param: "n".to_string(),
                reason: "d-electron count must be between 1 and 9 for a partially-filled d-shell"
                    .to_string(),
            });
        }
        if ten_dq < 0.0 {
            return Err(Error::InvalidParameter {
                param: "ten_dq".to_string(),
                reason: "crystal-field splitting magnitude must be non-negative".to_string(),
            });
        }
        if pairing_energy <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "pairing_energy".to_string(),
                reason: "mean pairing energy must be positive".to_string(),
            });
        }
        Ok(Self {
            n,
            environment,
            ten_dq,
            pairing_energy,
            soc_lambda,
        })
    }

    /// The free-ion (atomic) Hund's-rule ground term for this ion's `n`
    /// electron count. Environment-independent; see [`free_ion_term`].
    pub fn free_ion_term(&self) -> Result<FreeIonTerm> {
        free_ion_term(self.n)
    }

    /// Distinct one-electron crystal-field energy levels, sorted ascending,
    /// as `(energy, member orbitals)` groups, with numerically-degenerate
    /// orbitals merged (tolerance-based clustering). For
    /// [`CrystalFieldEnvironment::Octahedral`] / [`CrystalFieldEnvironment::Tetrahedral`]
    /// this yields exactly two levels (3-fold and 2-fold); for
    /// [`CrystalFieldEnvironment::TetragonalDistorted`] with generic `ds, dt`
    /// it yields three singlets plus one doublet (`xz`/`yz`, protected by the
    /// residual axial symmetry).
    fn orbital_energy_levels(&self) -> Result<Vec<(f64, Vec<DOrbital>)>> {
        let h_cf = crystal_field_hamiltonian(self.environment, self.ten_dq)?;
        let mut orbital_energies: Vec<(DOrbital, f64)> = DOrbital::all()
            .into_iter()
            .map(|orbital| (orbital, h_cf.get(orbital.index(), orbital.index()).re))
            .collect();
        orbital_energies.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let tolerance = 1e-9 * (self.ten_dq.abs() + 1.0);
        let mut levels: Vec<(f64, Vec<DOrbital>)> = Vec::new();
        for (orbital, e) in orbital_energies {
            if let Some(last) = levels.last_mut() {
                if (e - last.0).abs() < tolerance {
                    last.1.push(orbital);
                    continue;
                }
            }
            levels.push((e, vec![orbital]));
        }
        Ok(levels)
    }

    /// Ligand-field ground-state occupations of the crystal-field levels
    /// (Aufbau + Hund's rule 1 within each level, with the standard 10Dq-vs-
    /// pairing-energy competition between levels).
    fn ground_configuration(&self) -> Result<Vec<LevelOccupation>> {
        let levels = self.orbital_energy_levels()?;
        Ok(fill_levels(&levels, self.n, self.pairing_energy))
    }

    /// The *actual* ligand-field ground-state total spin `S`, accounting for
    /// high-spin/low-spin pairing (as opposed to the free-ion Hund's-rule
    /// maximum-multiplicity spin — see the module-level distinction).
    pub fn ground_state_spin(&self) -> Result<f64> {
        let occupations = self.ground_configuration()?;
        let unpaired: usize = occupations
            .iter()
            .map(LevelOccupation::unpaired_electrons)
            .sum();
        Ok(unpaired as f64 / 2.0)
    }

    /// Whether the ligand-field ground state is high-spin, i.e. whether its
    /// actual spin equals the free-ion (maximum-multiplicity) spin. For `n`
    /// with no high-spin/low-spin ambiguity (1,2,3,8,9), this is always `true`.
    pub fn is_high_spin(&self) -> Result<bool> {
        let s_actual = self.ground_state_spin()?;
        let s_free = free_ion_term(self.n)?.s;
        Ok((s_actual - s_free).abs() < 1e-9)
    }

    /// Combinatorial orbital degeneracy of the ground configuration: the
    /// maximum, over all crystal-field levels, of that level's
    /// `combinatorial_degeneracy()`. At most one level is ever partially
    /// filled away from a "special" occupation (empty/half/full) for the
    /// Aufbau sequence produced by [`fill_levels`], so taking the maximum
    /// (rather than a product over levels) correctly captures the ground
    /// configuration's total orbital degeneracy for every environment
    /// supported by this module.
    fn ground_configuration_degeneracy(&self) -> Result<u64> {
        let occupations = self.ground_configuration()?;
        Ok(occupations
            .iter()
            .map(LevelOccupation::combinatorial_degeneracy)
            .max()
            .unwrap_or(1))
    }

    /// Ground-term symmetry label (A/E/T) — see [`GroundTermSymmetry`].
    pub fn ground_term_symmetry(&self) -> Result<GroundTermSymmetry> {
        let degeneracy = self.ground_configuration_degeneracy()?;
        Ok(match degeneracy {
            1 => GroundTermSymmetry::A,
            2 => GroundTermSymmetry::E,
            _ => GroundTermSymmetry::T,
        })
    }

    /// Whether the ground term's orbital angular momentum is quenched to
    /// first order (`A` or `E` ground-term symmetry).
    pub fn is_orbitally_quenched(&self) -> Result<bool> {
        Ok(self.ground_term_symmetry()?.is_quenched())
    }

    /// The free-ion Landé g-factor, `g_J = 1 + [J(J+1)+S(S+1)-L(L+1)] / [2J(J+1)]`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] if the free-ion term has `J = 0`
    /// (e.g. the high-spin d4 ⁵D₀ term of Cr²⁺/Mn³⁺): the leading-order
    /// magnetic moment vanishes identically for a `J=0` level, so the g-factor
    /// is not physically meaningful at this order.
    pub fn lande_g_factor(&self) -> Result<f64> {
        let term = self.free_ion_term()?;
        if term.j < 1e-9 {
            return Err(Error::InvalidParameter {
                param: "j".to_string(),
                reason: format!(
                    "Lande g-factor is undefined for a J=0 free-ion ground term (d^{}: S={:.1}, \
                     L={:.0}, J=0 - the leading-order magnetic moment vanishes identically, e.g. \
                     Cr2+/Mn3+ high-spin 5D0)",
                    self.n, term.s, term.l
                ),
            });
        }
        let FreeIonTerm { s, l, j } = term;
        Ok(1.0 + (j * (j + 1.0) + s * (s + 1.0) - l * (l + 1.0)) / (2.0 * j * (j + 1.0)))
    }

    /// Spin-only effective moment `mu_eff = 2*sqrt(S(S+1))` (units of `mu_B`),
    /// using the *actual* (high-spin/low-spin-aware) ligand-field ground-state
    /// spin. This is the quantity that matches observed 3d-ion magnetic
    /// moments well, because the crystal field quenches orbital angular
    /// momentum long before spin-orbit coupling could enforce the free-ion
    /// Russell-Saunders `J` coupling (`10Dq >> lambda` for 3d ions).
    pub fn effective_moment_spin_only(&self) -> Result<f64> {
        let s = self.ground_state_spin()?;
        Ok(2.0 * (s * (s + 1.0)).sqrt())
    }

    /// Free-ion Landé effective moment `mu_eff = g_J*sqrt(J(J+1))` (units of
    /// `mu_B`). Appropriate when spin-orbit coupling dominates the crystal
    /// field (as for 4f rare earths); included for contrast with
    /// [`Self::effective_moment_spin_only`], which is normally the physically
    /// relevant quantity for 3d ions.
    ///
    /// # Errors
    ///
    /// Propagates the [`Self::lande_g_factor`] error for `J=0` free-ion terms.
    pub fn effective_moment_lande(&self) -> Result<f64> {
        let g = self.lande_g_factor()?;
        let j = self.free_ion_term()?.j;
        Ok(g * (j * (j + 1.0)).sqrt())
    }

    /// The single-particle crystal-field-plus-spin-orbit Hamiltonian (10x10,
    /// orbital⊗spin); see [`soc_hamiltonian`]. This is `n`-independent (a
    /// generic one-electron-in-the-full-orbital-space building block); for the
    /// `n`-aware demonstration of *this ion's* orbital quenching/unquenching,
    /// see [`Self::ground_manifold_orbital_moment_rms`].
    pub fn single_particle_hamiltonian(&self) -> Result<CMatrix> {
        soc_hamiltonian(self.environment, self.ten_dq, self.soc_lambda)
    }

    /// The orbitally-degenerate "frontier" level of the ground configuration
    /// (the crystal-field level with `combinatorial_degeneracy() > 1`), if one
    /// exists. Returns `None` for a purely `A`-type ground configuration
    /// (every level empty, full, or exactly half-filled) — see
    /// [`Self::ground_term_symmetry`].
    ///
    /// This is the physically relevant reduction for demonstrating orbital
    /// quenching/unquenching: whether a lone electron (or hole) can carry
    /// orbital angular momentum is a property of the many-electron ground
    /// *configuration*, not of a generic lone test-electron in the full
    /// five-orbital space (which would sit in the lower crystal-field level
    /// regardless of `n`, and hence look "unquenched" even for genuinely
    /// quenched ions like Cr³⁺'s exactly-half-filled, Pauli-symmetric t2g³).
    fn frontier_level(&self) -> Result<Option<LevelOccupation>> {
        let occupations = self.ground_configuration()?;
        Ok(occupations
            .into_iter()
            .find(|o| o.combinatorial_degeneracy() > 1))
    }

    /// Reduced single-quasiparticle spin-orbit Hamiltonian and `L_z` operator,
    /// confined to the frontier level's member orbitals only (dimension
    /// `2 * degeneracy`), plus a scale-invariant tolerance for
    /// eigenvalue-degeneracy detection. Returns `None` if there is no frontier
    /// level (see [`Self::frontier_level`]).
    fn frontier_reduced_hamiltonian(&self) -> Result<Option<(CMatrix, CMatrix)>> {
        let Some(level) = self.frontier_level()? else {
            return Ok(None);
        };
        let (lx, ly, lz) = orbital_angular_momentum_operators()?;
        let lx_sub = orbital_sub_block(&lx, &level.members)?;
        let ly_sub = orbital_sub_block(&ly, &level.members)?;
        let lz_sub = orbital_sub_block(&lz, &level.members)?;

        let sx = pauli_x()?;
        let sy = pauli_y()?;
        let sz = pauli_z()?;
        let lx_sx = kron(&lx_sub, &sx)?;
        let ly_sy = kron(&ly_sub, &sy)?;
        let lz_sz = kron(&lz_sub, &sz)?;
        // S_i = sigma_i / 2; the crystal-field energy is a constant additive
        // shift within a single level and is omitted (irrelevant to the
        // eigenvectors / orbital-moment analysis below).
        let h_reduced = lx_sx
            .add(&ly_sy)?
            .add(&lz_sz)?
            .scale_real(self.soc_lambda * 0.5);
        let lz_ext = kron(&lz_sub, &CMatrix::eye(2))?;
        Ok(Some((h_reduced, lz_ext)))
    }

    /// `<L_z>` in one representative ground eigenstate of the reduced
    /// spin-orbit Hamiltonian confined to the ground configuration's frontier
    /// (orbitally-degenerate) level (see `Self::frontier_level`), or exactly
    /// `0.0` if there is no such level (pure `A`-type ground configuration —
    /// no single-particle orbital angular momentum is available at this level
    /// of theory).
    ///
    /// A lone quasiparticle in a degenerate level is a half-integer
    /// (spin-1/2) system, so by Kramers' theorem every eigenvalue is at least
    /// doubly degenerate: within a degenerate eigenspace, this specific
    /// expectation value depends on an arbitrary (but deterministic) choice of
    /// orthonormal basis for that eigenspace. For a basis-independent
    /// diagnostic, use [`Self::ground_manifold_orbital_moment_rms`].
    pub fn expectation_lz_ground_state(&self) -> Result<f64> {
        let Some((h_reduced, lz_ext)) = self.frontier_reduced_hamiltonian()? else {
            return Ok(0.0);
        };
        let (_, eigenvectors) = h_reduced.hermitian_eigendecomposition()?;
        let dim = h_reduced.n();
        let ground: Vec<Complex> = (0..dim).map(|i| eigenvectors.get(i, 0)).collect();
        let mut expectation = Complex::ZERO;
        for (i, gi) in ground.iter().enumerate() {
            let mut row_sum = Complex::ZERO;
            for (j, gj) in ground.iter().enumerate() {
                row_sum = row_sum.add(&lz_ext.get(i, j).mul(gj));
            }
            expectation = expectation.add(&gi.conj().mul(&row_sum));
        }
        Ok(expectation.re)
    }

    /// Basis-independent RMS orbital moment of the (possibly Kramers- or
    /// otherwise degenerate) ground eigenspace of the reduced frontier-level
    /// spin-orbit Hamiltonian (see `Self::frontier_level`):
    /// `sqrt( (1/g) * sum_k <v_k|L_z^2|v_k> )` over the `g` eigenvectors
    /// degenerate with the lowest eigenvalue. Exactly `0.0` if there is no
    /// frontier level (pure `A`-type ground configuration).
    ///
    /// Because this is (up to normalization) `Tr(P L_z^2 P)` for the
    /// ground-eigenspace projector `P`, it is invariant under any unitary
    /// rotation of the degenerate eigenvectors among themselves, unlike
    /// [`Self::expectation_lz_ground_state`]. It is exactly zero — within this
    /// single-configuration model, restricted to the frontier level — for
    /// *both* quenched ground-term symmetries: `A` (no frontier level at all)
    /// and `E` (the cubic `E` representation contains no `L`-operator matrix
    /// elements whatsoever, a rigorous group-theory fact, not merely a small
    /// admixture). It is of order 1 for an unquenched (`T`) ground term — the
    /// robust demonstration that spin-orbit coupling lifts orbital quenching
    /// specifically for `T`-term ions. (Real ions can show a *small*,
    /// second-order-in-`lambda/ten_dq` moment from cross-level `t2g`-`eg`
    /// admixture, but that effect is a form of configuration mixing beyond
    /// this module's single-configuration fidelity boundary, so it is not
    /// modeled here.)
    pub fn ground_manifold_orbital_moment_rms(&self) -> Result<f64> {
        let Some((h_reduced, lz_ext)) = self.frontier_reduced_hamiltonian()? else {
            return Ok(0.0);
        };
        let (energies, eigenvectors) = h_reduced.hermitian_eigendecomposition()?;
        let dim = h_reduced.n();
        let e0 = energies[0];
        let tolerance = 1e-7 * (e0.abs() + 1.0);
        let ground_indices: Vec<usize> = energies
            .iter()
            .enumerate()
            .filter(|&(_, &e)| (e - e0).abs() < tolerance)
            .map(|(idx, _)| idx)
            .collect();

        let mut sum_lz_squared = 0.0_f64;
        for &k in &ground_indices {
            let v: Vec<Complex> = (0..dim).map(|i| eigenvectors.get(i, k)).collect();
            let mut lzv = vec![Complex::ZERO; dim];
            for (i, slot) in lzv.iter_mut().enumerate() {
                let mut acc = Complex::ZERO;
                for (j, vj) in v.iter().enumerate() {
                    acc = acc.add(&lz_ext.get(i, j).mul(vj));
                }
                *slot = acc;
            }
            sum_lz_squared += lzv.iter().map(Complex::norm_sq).sum::<f64>();
        }
        Ok((sum_lz_squared / ground_indices.len() as f64).sqrt())
    }
}

// =============================================================================
// Preset ions
// =============================================================================
//
// Representative (order-of-magnitude, eV) crystal-field parameters informed by
// standard ligand-field-theory references (e.g. Figgis & Hitchman 2000,
// Abragam & Bleaney 1970); not intended for precision spectroscopic fitting
// (see the module-level fidelity boundary). Each preset's 10Dq/pairing-energy
// pair is chosen to reproduce the well-known high-spin/low-spin classification
// and ground-term symmetry of the corresponding real ion.

impl CrystalFieldModel {
    /// Ti³⁺ (3d¹), e.g. `[Ti(H2O)6]³⁺`. No high-spin/low-spin ambiguity.
    /// Free-ion term ²D; ligand-field ground term ²T2g (unquenched).
    pub fn ti3_plus() -> Result<Self> {
        Self::new(1, CrystalFieldEnvironment::Octahedral, 2.5, 3.0, 0.019)
    }

    /// V³⁺ (3d²), e.g. `[V(H2O)6]³⁺`. Free-ion term ³F; ligand-field ground
    /// term ³T1g (unquenched).
    pub fn v3_plus() -> Result<Self> {
        Self::new(2, CrystalFieldEnvironment::Octahedral, 2.3, 3.0, 0.026)
    }

    /// Cr³⁺ (3d³), e.g. `[Cr(H2O)6]³⁺` (chrome alum, ruby dopant). Exactly
    /// half-filled t2g regardless of field strength. Free-ion term ⁴F;
    /// ligand-field ground term ⁴A2g (quenched) — famous for showing
    /// almost exactly the spin-only moment.
    pub fn cr3_plus() -> Result<Self> {
        Self::new(3, CrystalFieldEnvironment::Octahedral, 2.2, 3.0, 0.034)
    }

    /// Mn³⁺ high-spin (3d⁴), e.g. `[Mn(H2O)6]³⁺` (Jahn-Teller active). Free-ion
    /// term ⁵D with `J=0` (Landé g/mu undefined); ligand-field ground term
    /// ⁵Eg (quenched).
    pub fn mn3_plus_high_spin() -> Result<Self> {
        Self::new(4, CrystalFieldEnvironment::Octahedral, 1.8, 3.3, 0.044)
    }

    /// Fe³⁺ high-spin (3d⁵), e.g. `[Fe(H2O)6]³⁺`. Free-ion and ligand-field
    /// term ⁶A1g (quenched); textbook spin-only moment 5.92 `mu_B`.
    pub fn fe3_plus() -> Result<Self> {
        Self::new(5, CrystalFieldEnvironment::Octahedral, 1.7, 3.7, 0.057)
    }

    /// Mn²⁺ high-spin (3d⁵), e.g. `[Mn(H2O)6]²⁺`. Same d⁵ physics as
    /// [`Self::fe3_plus`] (⁶A1g, quenched), with a weaker crystal field.
    pub fn mn2_plus() -> Result<Self> {
        Self::new(5, CrystalFieldEnvironment::Octahedral, 1.0, 3.2, 0.042)
    }

    /// Fe²⁺ high-spin (3d⁶), e.g. `[Fe(H2O)6]²⁺`. Ligand-field ground term
    /// ⁵T2g (unquenched — one of the largest orbital contributions among 3d
    /// ions).
    pub fn fe2_plus_high_spin() -> Result<Self> {
        Self::new(6, CrystalFieldEnvironment::Octahedral, 1.3, 2.6, 0.051)
    }

    /// Fe²⁺ low-spin (3d⁶), e.g. `[Fe(CN)6]⁴⁻` (ferrocyanide, strong field).
    /// Ligand-field ground term ¹A1g, `S=0` (diamagnetic).
    pub fn fe2_plus_low_spin() -> Result<Self> {
        Self::new(6, CrystalFieldEnvironment::Octahedral, 3.8, 2.6, 0.051)
    }

    /// Co²⁺ high-spin (3d⁷), e.g. `[Co(H2O)6]²⁺`. Ligand-field ground term
    /// ⁴T1g (unquenched, large orbital contribution).
    pub fn co2_plus() -> Result<Self> {
        Self::new(7, CrystalFieldEnvironment::Octahedral, 1.1, 2.8, 0.064)
    }

    /// Ni²⁺ (3d⁸), e.g. `[Ni(H2O)6]²⁺`. No high-spin/low-spin ambiguity in
    /// O_h. Ligand-field ground term ³A2g (quenched); textbook spin-only
    /// moment 2.83 `mu_B`.
    pub fn ni2_plus() -> Result<Self> {
        Self::new(8, CrystalFieldEnvironment::Octahedral, 1.05, 2.5, 0.078)
    }

    /// Cu²⁺ (3d⁹), e.g. `[Cu(H2O)6]²⁺` (modeled here as regular octahedral;
    /// real Cu²⁺ complexes are Jahn-Teller distorted — see
    /// [`CrystalFieldEnvironment::TetragonalDistorted`] for that extension).
    /// Ligand-field ground term ²Eg (quenched, Jahn-Teller active); textbook
    /// spin-only moment 1.73 `mu_B`.
    pub fn cu2_plus() -> Result<Self> {
        Self::new(9, CrystalFieldEnvironment::Octahedral, 1.56, 2.5, 0.103)
    }

    /// Co³⁺ low-spin (3d⁶), e.g. `[Co(NH3)6]³⁺` (one of the largest 10Dq
    /// values among 3d ions — virtually all Co³⁺ complexes are low-spin).
    /// Ligand-field ground term ¹A1g, `S=0` (diamagnetic).
    pub fn co3_plus_low_spin() -> Result<Self> {
        Self::new(6, CrystalFieldEnvironment::Octahedral, 3.0, 2.6, 0.069)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_free_ion_term_hund_rule_table() {
        let expected: [(u8, f64, f64); 9] = [
            (1, 0.5, 2.0),
            (2, 1.0, 3.0),
            (3, 1.5, 3.0),
            (4, 2.0, 2.0),
            (5, 2.5, 0.0),
            (6, 2.0, 2.0),
            (7, 1.5, 3.0),
            (8, 1.0, 3.0),
            (9, 0.5, 2.0),
        ];
        for (n, s, l) in expected {
            let term = free_ion_term(n).expect("valid n");
            assert!(
                (term.s - s).abs() < 1e-12,
                "d{}: S expected {}, got {}",
                n,
                s,
                term.s
            );
            assert!(
                (term.l - l).abs() < 1e-12,
                "d{}: L expected {}, got {}",
                n,
                l,
                term.l
            );
        }
    }

    #[test]
    fn test_electron_hole_symmetry_of_s_and_l() {
        for n in 1..=4u8 {
            let term_n = free_ion_term(n).expect("valid");
            let term_hole = free_ion_term(10 - n).expect("valid");
            assert!(
                (term_n.l - term_hole.l).abs() < 1e-12,
                "L({}) = {} != L({}) = {}",
                n,
                term_n.l,
                10 - n,
                term_hole.l
            );
            assert!(
                (term_n.s - term_hole.s).abs() < 1e-12,
                "S({}) != S({})",
                n,
                10 - n
            );
        }
    }

    #[test]
    fn test_hund_third_rule_j_values() {
        // Ti3+ (d1): less than half filled, J = |L-S| = |2-0.5| = 1.5.
        let ti = free_ion_term(1).expect("valid");
        assert!((ti.j - 1.5).abs() < 1e-12);
        // Cu2+ (d9): more than half filled, J = L+S = 2+0.5 = 2.5.
        let cu = free_ion_term(9).expect("valid");
        assert!((cu.j - 2.5).abs() < 1e-12);
        // Fe3+/Mn2+ (d5): L=0, J=S=2.5 regardless of branch.
        let d5 = free_ion_term(5).expect("valid");
        assert!((d5.j - 2.5).abs() < 1e-12);
        // Mn3+/Cr2+ high-spin (d4): L=S=2, less than half filled, J=|2-2|=0.
        let d4 = free_ion_term(4).expect("valid");
        assert!(d4.j.abs() < 1e-12);
    }

    #[test]
    fn test_free_ion_term_rejects_out_of_range_n() {
        assert!(free_ion_term(11).is_err());
    }

    #[test]
    fn test_model_validation_rejects_bad_n() {
        assert!(
            CrystalFieldModel::new(0, CrystalFieldEnvironment::Octahedral, 2.0, 2.5, 0.05).is_err()
        );
        assert!(
            CrystalFieldModel::new(10, CrystalFieldEnvironment::Octahedral, 2.0, 2.5, 0.05)
                .is_err()
        );
    }

    #[test]
    fn test_model_validation_rejects_negative_ten_dq() {
        assert!(
            CrystalFieldModel::new(5, CrystalFieldEnvironment::Octahedral, -1.0, 2.5, 0.05)
                .is_err()
        );
    }

    #[test]
    fn test_model_validation_rejects_non_positive_pairing_energy() {
        assert!(
            CrystalFieldModel::new(5, CrystalFieldEnvironment::Octahedral, 2.0, 0.0, 0.05).is_err()
        );
        assert!(
            CrystalFieldModel::new(5, CrystalFieldEnvironment::Octahedral, 2.0, -1.0, 0.05)
                .is_err()
        );
    }

    #[test]
    fn test_large_ten_dq_forces_low_spin_d4_to_d7() {
        for n in 4..=7u8 {
            let model =
                CrystalFieldModel::new(n, CrystalFieldEnvironment::Octahedral, 1.0e6, 2.5, 0.05)
                    .expect("valid model");
            assert!(
                !model.is_high_spin().expect("computes"),
                "d{} should be forced low-spin",
                n
            );
        }
    }

    #[test]
    fn test_lande_g_factor_ti3_plus() {
        let ti = CrystalFieldModel::ti3_plus().expect("preset builds");
        let g = ti.lande_g_factor().expect("J != 0");
        assert!(
            (g - 0.8).abs() < 1e-9,
            "Ti3+ Lande g should be 0.8, got {}",
            g
        );
    }

    #[test]
    fn test_lande_g_factor_is_two_when_l_is_zero() {
        let fe3 = CrystalFieldModel::fe3_plus().expect("preset builds");
        let g = fe3.lande_g_factor().expect("J != 0");
        assert!(
            (g - 2.0).abs() < 1e-9,
            "L=0 term should have g=2, got {}",
            g
        );
        // When L=0, the Lande and spin-only moments must coincide exactly.
        let mu_lande = fe3.effective_moment_lande().expect("J != 0");
        let mu_so = fe3.effective_moment_spin_only().expect("always defined");
        assert!((mu_lande - mu_so).abs() < 1e-9);
    }

    #[test]
    fn test_lande_g_factor_undefined_for_j_zero_term() {
        let mn3 = CrystalFieldModel::mn3_plus_high_spin().expect("preset builds");
        assert!(mn3.free_ion_term().expect("computes").j.abs() < 1e-9);
        assert!(
            mn3.lande_g_factor().is_err(),
            "g should be undefined for a J=0 term"
        );
        assert!(
            mn3.effective_moment_lande().is_err(),
            "mu_lande should be undefined for J=0"
        );
        let mu_so = mn3
            .effective_moment_spin_only()
            .expect("spin-only always defined");
        assert!(
            (mu_so - 4.899).abs() < 1e-2,
            "Mn3+ HS spin-only moment should be ~4.90, got {}",
            mu_so
        );
    }

    #[test]
    fn test_spin_only_moments_match_textbook_values() {
        let fe3 = CrystalFieldModel::fe3_plus()
            .expect("preset")
            .effective_moment_spin_only()
            .expect("ok");
        let ni2 = CrystalFieldModel::ni2_plus()
            .expect("preset")
            .effective_moment_spin_only()
            .expect("ok");
        let cr3 = CrystalFieldModel::cr3_plus()
            .expect("preset")
            .effective_moment_spin_only()
            .expect("ok");
        let cu2 = CrystalFieldModel::cu2_plus()
            .expect("preset")
            .effective_moment_spin_only()
            .expect("ok");
        assert!(
            (fe3 - 5.92).abs() < 0.01,
            "Fe3+ mu_so expected 5.92, got {}",
            fe3
        );
        assert!(
            (ni2 - 2.83).abs() < 0.01,
            "Ni2+ mu_so expected 2.83, got {}",
            ni2
        );
        assert!(
            (cr3 - 3.87).abs() < 0.01,
            "Cr3+ mu_so expected 3.87, got {}",
            cr3
        );
        assert!(
            (cu2 - 1.73).abs() < 0.01,
            "Cu2+ mu_so expected 1.73, got {}",
            cu2
        );
    }

    #[test]
    fn test_is_orbitally_quenched_matches_symmetry() {
        let ti = CrystalFieldModel::ti3_plus().expect("preset");
        assert!(
            !ti.is_orbitally_quenched().expect("computes"),
            "Ti3+ (T2g) should be unquenched"
        );
        let cr = CrystalFieldModel::cr3_plus().expect("preset");
        assert!(
            cr.is_orbitally_quenched().expect("computes"),
            "Cr3+ (A2g) should be quenched"
        );
        let cu = CrystalFieldModel::cu2_plus().expect("preset");
        assert!(
            cu.is_orbitally_quenched().expect("computes"),
            "Cu2+ (Eg) should be quenched"
        );
    }

    #[test]
    fn test_all_presets_match_expected_physics() {
        struct Expected {
            name: &'static str,
            n: u8,
            s: f64,
            high_spin: bool,
            symmetry: GroundTermSymmetry,
        }
        let cases = [
            Expected {
                name: "Ti3+",
                n: 1,
                s: 0.5,
                high_spin: true,
                symmetry: GroundTermSymmetry::T,
            },
            Expected {
                name: "V3+",
                n: 2,
                s: 1.0,
                high_spin: true,
                symmetry: GroundTermSymmetry::T,
            },
            Expected {
                name: "Cr3+",
                n: 3,
                s: 1.5,
                high_spin: true,
                symmetry: GroundTermSymmetry::A,
            },
            Expected {
                name: "Mn3+ HS",
                n: 4,
                s: 2.0,
                high_spin: true,
                symmetry: GroundTermSymmetry::E,
            },
            Expected {
                name: "Fe3+",
                n: 5,
                s: 2.5,
                high_spin: true,
                symmetry: GroundTermSymmetry::A,
            },
            Expected {
                name: "Mn2+",
                n: 5,
                s: 2.5,
                high_spin: true,
                symmetry: GroundTermSymmetry::A,
            },
            Expected {
                name: "Fe2+ HS",
                n: 6,
                s: 2.0,
                high_spin: true,
                symmetry: GroundTermSymmetry::T,
            },
            Expected {
                name: "Fe2+ LS",
                n: 6,
                s: 0.0,
                high_spin: false,
                symmetry: GroundTermSymmetry::A,
            },
            Expected {
                name: "Co2+",
                n: 7,
                s: 1.5,
                high_spin: true,
                symmetry: GroundTermSymmetry::T,
            },
            Expected {
                name: "Ni2+",
                n: 8,
                s: 1.0,
                high_spin: true,
                symmetry: GroundTermSymmetry::A,
            },
            Expected {
                name: "Cu2+",
                n: 9,
                s: 0.5,
                high_spin: true,
                symmetry: GroundTermSymmetry::E,
            },
            Expected {
                name: "Co3+ LS",
                n: 6,
                s: 0.0,
                high_spin: false,
                symmetry: GroundTermSymmetry::A,
            },
        ];

        let models: Vec<CrystalFieldModel> = vec![
            CrystalFieldModel::ti3_plus().expect("preset"),
            CrystalFieldModel::v3_plus().expect("preset"),
            CrystalFieldModel::cr3_plus().expect("preset"),
            CrystalFieldModel::mn3_plus_high_spin().expect("preset"),
            CrystalFieldModel::fe3_plus().expect("preset"),
            CrystalFieldModel::mn2_plus().expect("preset"),
            CrystalFieldModel::fe2_plus_high_spin().expect("preset"),
            CrystalFieldModel::fe2_plus_low_spin().expect("preset"),
            CrystalFieldModel::co2_plus().expect("preset"),
            CrystalFieldModel::ni2_plus().expect("preset"),
            CrystalFieldModel::cu2_plus().expect("preset"),
            CrystalFieldModel::co3_plus_low_spin().expect("preset"),
        ];

        for (case, model) in cases.iter().zip(models.iter()) {
            assert_eq!(model.n, case.n, "{}: wrong electron count", case.name);
            let s = model.ground_state_spin().expect("computes");
            assert!(
                (s - case.s).abs() < 1e-9,
                "{}: S expected {}, got {}",
                case.name,
                case.s,
                s
            );
            let hs = model.is_high_spin().expect("computes");
            assert_eq!(
                hs, case.high_spin,
                "{}: high-spin classification mismatch",
                case.name
            );
            let sym = model.ground_term_symmetry().expect("computes");
            assert_eq!(
                sym, case.symmetry,
                "{}: ground term symmetry mismatch",
                case.name
            );
        }
    }

    #[test]
    fn test_expectation_lz_ground_state_is_finite_and_bounded() {
        // Ti3+'s frontier level is the 3-fold t2g set, which behaves as an
        // effective l=1 (its Lz sub-block has eigenvalues exactly {-1,0,1}),
        // so <Lz> is bounded by 1, not the full l=2 spectrum's 2.
        let ti = CrystalFieldModel::ti3_plus().expect("preset");
        let lz = ti.expectation_lz_ground_state().expect("computes");
        assert!(lz.is_finite());
        assert!(
            lz.abs() <= 1.001,
            "<Lz> should be bounded by the effective t2g l=1 spectrum, got {}",
            lz
        );
    }

    #[test]
    fn test_soc_lifts_orbital_moment_for_unquenched_t_term() {
        // Ti3+ (d1): ground term T2g (unquenched) -> substantial ground-manifold
        // RMS <Lz>, of order 1, even for small soc_lambda.
        let ti = CrystalFieldModel::ti3_plus().expect("preset");
        assert_eq!(
            ti.ground_term_symmetry().expect("computes"),
            GroundTermSymmetry::T
        );
        let rms_t = ti.ground_manifold_orbital_moment_rms().expect("computes");
        assert!(rms_t.is_finite());
        assert!(
            rms_t > 0.3,
            "T-term ground level should carry substantial orbital moment, got {}",
            rms_t
        );

        // Cr3+ (d3): ground term A2g (quenched, exactly-half-filled t2g^3) ->
        // no frontier level exists at all, so the RMS <Lz> is exactly zero in
        // this single-configuration model.
        let cr = CrystalFieldModel::cr3_plus().expect("preset");
        assert_eq!(
            cr.ground_term_symmetry().expect("computes"),
            GroundTermSymmetry::A
        );
        let rms_a = cr.ground_manifold_orbital_moment_rms().expect("computes");
        assert_eq!(
            rms_a, 0.0,
            "A-term ground level has no frontier level: RMS <Lz> must be exactly 0"
        );
        assert_eq!(cr.expectation_lz_ground_state().expect("computes"), 0.0);
    }

    #[test]
    fn test_e_term_orbital_moment_is_exactly_zero() {
        // Cu2+ (d9): ground term Eg (quenched). The cubic E representation
        // has no L-operator matrix elements at all (not merely "small"), so
        // both diagnostics must be exactly zero within the frontier (eg) level.
        let cu = CrystalFieldModel::cu2_plus().expect("preset");
        assert_eq!(
            cu.ground_term_symmetry().expect("computes"),
            GroundTermSymmetry::E
        );
        assert_eq!(
            cu.ground_manifold_orbital_moment_rms().expect("computes"),
            0.0
        );
        assert_eq!(cu.expectation_lz_ground_state().expect("computes"), 0.0);

        // Mn3+ high-spin (d4): also ground term Eg (quenched, eg^1 hole-free
        // single-electron frontier level), same exact-zero expectation.
        let mn3 = CrystalFieldModel::mn3_plus_high_spin().expect("preset");
        assert_eq!(
            mn3.ground_term_symmetry().expect("computes"),
            GroundTermSymmetry::E
        );
        assert_eq!(
            mn3.ground_manifold_orbital_moment_rms().expect("computes"),
            0.0
        );
    }

    #[test]
    fn test_t_term_presets_all_show_unquenched_orbital_moment() {
        // Every T-term preset (Ti3+, V3+, Fe2+ HS, Co2+) should show a
        // substantial (order-1) ground-manifold RMS orbital moment.
        for model in [
            CrystalFieldModel::ti3_plus().expect("preset"),
            CrystalFieldModel::v3_plus().expect("preset"),
            CrystalFieldModel::fe2_plus_high_spin().expect("preset"),
            CrystalFieldModel::co2_plus().expect("preset"),
        ] {
            assert_eq!(
                model.ground_term_symmetry().expect("computes"),
                GroundTermSymmetry::T
            );
            let rms = model
                .ground_manifold_orbital_moment_rms()
                .expect("computes");
            assert!(
                rms > 0.1,
                "T-term preset (d{}) should show unquenched orbital moment, got {}",
                model.n,
                rms
            );
        }
    }

    #[test]
    fn test_single_particle_hamiltonian_dimension() {
        let ti = CrystalFieldModel::ti3_plus().expect("preset");
        let h = ti.single_particle_hamiltonian().expect("builds");
        assert_eq!(h.n(), 10);
    }
}
