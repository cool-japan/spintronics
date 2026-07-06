//! Valence-bond quantum states: overlaps and Heisenberg matrix elements.
//!
//! A valence-bond (VB) state `|VB_α⟩` for a dimer covering `α` is the tensor product of
//! singlets over its dimers:
//!
//! ```text
//! |VB_α⟩ = ⊗_{(i,j)∈α, i<j} |s_ij⟩,   |s_ij⟩ = (1/√2)(|↑_i↓_j⟩ - |↓_i↑_j⟩)
//! ```
//!
//! The **canonical sign convention** fixes, for every dimer, the *smaller* global site
//! index as the first ("up-first") slot in the singlet formula. This is well-defined for
//! *any* graph (bipartite or not) since it depends only on integer site-index ordering,
//! never on a two-coloring — which matters because non-bipartite clusters (e.g. a
//! triangular plaquette) have no consistent two-coloring to hang a sign convention on.
//!
//! # Overlaps: the Sutherland loop-counting formula
//!
//! For two coverings `α`, `β` of the same `N` sites, superimpose their dimers: every
//! non-monomer site has exactly one `α`-edge and one `β`-edge, so the union graph
//! decomposes into disjoint simple cycles ("transition-graph loops"), each visiting an
//! even number `2p` of sites (a shared dimer, present in both coverings, forms a
//! degenerate 2-cycle with `p = 1`). Sites that are monomers in *both* coverings (see
//! [`super::dimer::DimerCovering::monomers`]) sit outside every loop and contribute
//! nothing (their spin factor is common to both kets and cancels to 1).
//!
//! For a single loop with cyclic vertex sequence `v_1, …, v_{2p}` (alternating `α`- and
//! `β`-edges), direct summation over the two compatible Néel colorings of the loop gives
//!
//! ```text
//! overlap of loop = (-1)^p · S_loop · 2^{1-p},   S_loop = Π_{k=1}^{2p} sign(v_{k+1} - v_k)
//! ```
//!
//! (indices mod `2p`; `S_loop` is independent of traversal direction and starting
//! vertex). Multiplying over all loops of the transition graph gives the **Sutherland
//! overlap formula** used here:
//!
//! ```text
//! ⟨VB_α|VB_β⟩ = 2^{n_loops - N_paired/2} · Π_loops [(-1)^{p_ℓ} · S_ℓ]
//! ```
//!
//! where `N_paired = N - (number of shared monomers)`. This reduces to the textbook
//! bipartite Marshall-sign result on bipartite lattices but — because it never invokes a
//! two-coloring — is equally valid on non-bipartite (frustrated) clusters, which is
//! exactly what is needed for triangular/kagome RVB physics. The overall sign per basis
//! state is a "gauge" choice (rescaling any `|VB_α⟩ → ε_α|VB_α⟩` with `ε_α = ±1` leaves
//! all physical eigenvalues of the generalized eigenproblem in [`super::solver`]
//! invariant); what matters is that the *same* convention is used everywhere, which this
//! module guarantees by construction.
//!
//! # Heisenberg matrix elements
//!
//! Using the spin-1/2 swap identity `S_i·S_j = (1/2)P_ij - (1/4)I`:
//!
//! - If `i, j` are already dimer partners in `β`: `S_i·S_j|VB_β⟩ = -(3/4)|VB_β⟩` exactly
//!   (singlet eigenvalue).
//! - Otherwise, with `k = partner(i)`, `l = partner(j)` (all four sites distinct):
//!   `S_i·S_j|VB_β⟩ = -(1/4)|VB_β⟩ + (1/2)·χ·|VB_{β'}⟩`, where `β'` is `β` with dimers
//!   `(i,k)`, `(j,l)` replaced by the reconnected pair `(i,l)`, `(j,k)`, and
//!   `χ = sign(k−i)·sign(l−j)·sign(l−i)·sign(k−j)` is the sign relating the swap
//!   operator's action (computed in a manifestly convention-free "role" basis) back to
//!   the canonical site-index-ordered singlets. This is an **operator identity** (holds
//!   term-by-term, not just as an inner product), so it assembles the full Hamiltonian
//!   matrix directly: `H_{αβ} = Σ_bonds [diagonal · S_{αβ} + reconnection · S_{α,β'}]`.
//!
//! Note `β'` need not itself be one of the enumerated (nearest-neighbor-restricted) basis
//! coverings — the overlap formula above is fully general and applies to *any* two
//! perfect matchings of the same sites, so `S_{α,β'}` is always well-defined even when
//! `β'` uses a "long-range" (non-bond) dimer. This makes the restricted-basis Hamiltonian
//! matrix element *exact* (a genuine Rayleigh-Ritz projection onto the chosen subspace,
//! not an additional approximation on top of the basis truncation).
//!
//! # References
//!
//! - P.W. Anderson, Mater. Res. Bull. 8, 153 (1973).
//! - B. Sutherland, Phys. Rev. B 37, 3786 (1988).
//! - S. Liang, B. Doucot, P.W. Anderson, Phys. Rev. Lett. 61, 365 (1988).

use super::dimer::DimerCovering;
use crate::error::{self, Result};

/// Sites beyond which brute-force dense-amplitude construction (2^N floats) is refused.
const MAX_DENSE_CONSTRUCTION_SITES: usize = 20;

/// A valence-bond (dimer-covering) quantum state, i.e. a [`DimerCovering`] together with
/// the total site count needed to interpret it as a state in the full spin Hilbert space.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ValenceBondState {
    /// The underlying dimer (or monomer-dimer) covering.
    pub covering: DimerCovering,
    /// Total number of sites (paired + monomer).
    pub num_sites: usize,
}

impl ValenceBondState {
    /// Construct a valence-bond state, validating that the covering matches `num_sites`.
    pub fn new(covering: DimerCovering, num_sites: usize) -> Result<Self> {
        if covering.num_sites() != num_sites {
            return Err(error::invalid_param(
                "covering",
                "covering site count does not match num_sites",
            ));
        }
        Ok(Self {
            covering,
            num_sites,
        })
    }

    /// Overlap `⟨self|other⟩` via the Sutherland loop-counting formula.
    ///
    /// # Errors
    ///
    /// Returns an error if the two states have different site counts or different
    /// (shared) monomer sets.
    pub fn overlap(&self, other: &Self) -> Result<f64> {
        if self.num_sites != other.num_sites {
            return Err(error::invalid_param(
                "other",
                "valence-bond states must share the same num_sites to compute an overlap",
            ));
        }
        overlap(&self.covering, &other.covering, self.num_sites)
    }

    /// Dense amplitude vector (length `2^num_sites`) in the computational S^z basis
    /// (bit `i` of the basis index is 0 for spin-up, 1 for spin-down at site `i`).
    ///
    /// Intended for small clusters (validation/cross-checking), not production-scale use.
    ///
    /// # Errors
    ///
    /// Returns an error if the covering has monomers (a dense single-vector
    /// representation is only well-defined for a perfect matching) or if `num_sites`
    /// exceeds `MAX_DENSE_CONSTRUCTION_SITES`.
    pub fn to_dense_amplitudes(&self) -> Result<Vec<f64>> {
        dense_amplitudes(&self.covering, self.num_sites)
    }
}

/// The action of `coupling_j · S_i·S_j` on a valence-bond covering `β`, decomposed into
/// its (at most two) valence-bond-basis components.
#[derive(Debug, Clone)]
pub(crate) struct BondContribution {
    /// Coefficient of `|VB_β⟩` itself: `-3/4·J` if `(i,j)` are partners in `β`, else
    /// `-1/4·J`.
    pub diagonal: f64,
    /// Present when `(i,j)` are *not* already partners in `β`: the reconnected covering
    /// `β'` and the coefficient of `|VB_{β'}⟩`, equal to `(J/2)·χ`.
    pub reconnection: Option<(DimerCovering, f64)>,
}

/// Compute `⟨VB_a|VB_b⟩` via the Sutherland loop-counting formula (see module docs).
///
/// # Errors
///
/// Returns an error if `a` and `b` have different monomer sets (no consistent shared
/// "spectator" convention would exist), or if `num_sites` does not match either covering.
pub(crate) fn overlap(a: &DimerCovering, b: &DimerCovering, num_sites: usize) -> Result<f64> {
    if a.num_sites() != num_sites || b.num_sites() != num_sites {
        return Err(error::invalid_param(
            "num_sites",
            "covering site count does not match num_sites",
        ));
    }
    if a.monomers != b.monomers {
        return Err(error::invalid_param(
            "b",
            "overlap requires both coverings to share the same (possibly empty) monomer set",
        ));
    }
    let n_paired = num_sites - a.monomers.len();
    if n_paired % 2 != 0 {
        return Err(error::invalid_param(
            "num_sites",
            "the number of paired (non-monomer) sites must be even",
        ));
    }

    let mut visited = vec![false; num_sites];
    for &m in &a.monomers {
        visited[m] = true;
    }

    let mut n_loops: i32 = 0;
    let mut total_sign: f64 = 1.0;

    for start in 0..num_sites {
        if visited[start] {
            continue;
        }
        let mut current = start;
        let mut use_a = true;
        let mut loop_sites: i32 = 0;
        let mut edge_sign_product: f64 = 1.0;
        loop {
            visited[current] = true;
            loop_sites += 1;
            let next = if use_a {
                a.partner_of(current)
            } else {
                b.partner_of(current)
            }
            .ok_or_else(|| {
                error::invalid_param(
                    "a/b",
                    "encountered an unpaired site while tracing a transition-graph loop",
                )
            })?;
            edge_sign_product *= signum_diff(next, current);
            use_a = !use_a;
            current = next;
            if current == start {
                break;
            }
        }
        if loop_sites % 2 != 0 {
            return Err(error::numerical_error(
                "transition-graph loop has odd length; coverings are inconsistent",
            ));
        }
        let p = loop_sites / 2;
        n_loops += 1;
        let loop_sign = if p % 2 == 0 { 1.0 } else { -1.0 };
        total_sign *= loop_sign * edge_sign_product;
    }

    let magnitude = 2.0_f64.powi(n_loops - (n_paired as i32) / 2);
    Ok(total_sign * magnitude)
}

/// `sign(next - current)` as ±1.0 for distinct usize site indices.
fn signum_diff(next: usize, current: usize) -> f64 {
    if next > current {
        1.0
    } else {
        -1.0
    }
}

/// `sign(x)` as ±1.0 for a nonzero `i64`.
fn signum_i64(x: i64) -> f64 {
    if x > 0 {
        1.0
    } else if x < 0 {
        -1.0
    } else {
        0.0
    }
}

/// Given covering `β` and a bond `(i, j)` where `i, j` are *not* already dimer partners
/// in `β`, build the reconnected covering `β'` (dimers `(i,k)`, `(j,l)` replaced by
/// `(i,l)`, `(j,k)`, where `k = partner(i)`, `l = partner(j)`), and the sign `χ` such
/// that `P_ij|VB_β⟩ = χ|VB_{β'}⟩` in the canonical (site-index-ordered) convention.
///
/// Returns `Ok(None)` if `i` and `j` are already partners in `β` (the caller should use
/// the pure diagonal singlet-eigenvalue rule in that case).
pub(crate) fn reconnect(
    covering: &DimerCovering,
    i: usize,
    j: usize,
) -> Result<Option<(DimerCovering, f64)>> {
    let k = covering.partner_of(i).ok_or_else(|| {
        error::invalid_param(
            "i",
            "site is out of range or a monomer (not part of any dimer)",
        )
    })?;
    let l = covering.partner_of(j).ok_or_else(|| {
        error::invalid_param(
            "j",
            "site is out of range or a monomer (not part of any dimer)",
        )
    })?;
    if k == j {
        // (i,j) are already partners (implies l == i too).
        return Ok(None);
    }
    let chi = signum_i64(k as i64 - i as i64)
        * signum_i64(l as i64 - j as i64)
        * signum_i64(l as i64 - i as i64)
        * signum_i64(k as i64 - j as i64);

    let num_sites = covering.num_sites();
    let old_ik = (i.min(k), i.max(k));
    let old_jl = (j.min(l), j.max(l));
    let mut new_dimers: Vec<(usize, usize)> = covering
        .dimers
        .iter()
        .copied()
        .filter(|&d| d != old_ik && d != old_jl)
        .collect();
    new_dimers.push((i.min(l), i.max(l)));
    new_dimers.push((j.min(k), j.max(k)));

    let new_covering = DimerCovering::new(num_sites, new_dimers)?;
    Ok(Some((new_covering, chi)))
}

/// Decompose the action of `coupling_j · S_i·S_j` on covering `β` into VB-basis pieces.
///
/// # Errors
///
/// Propagates errors from [`reconnect`] (out-of-range or monomer sites).
pub(crate) fn bond_contribution(
    covering: &DimerCovering,
    i: usize,
    j: usize,
    coupling_j: f64,
) -> Result<BondContribution> {
    match reconnect(covering, i, j)? {
        None => Ok(BondContribution {
            diagonal: -0.75 * coupling_j,
            reconnection: None,
        }),
        Some((beta_prime, chi)) => Ok(BondContribution {
            diagonal: -0.25 * coupling_j,
            reconnection: Some((beta_prime, 0.5 * chi * coupling_j)),
        }),
    }
}

/// Dense amplitude vector (length `2^num_sites`) for a perfect-matching covering, in the
/// computational S^z basis (bit `i` = 0 for up, 1 for down at site `i`).
pub(crate) fn dense_amplitudes(covering: &DimerCovering, num_sites: usize) -> Result<Vec<f64>> {
    if !covering.monomers.is_empty() {
        return Err(error::invalid_param(
            "covering",
            "dense amplitude construction requires a perfect matching (no monomers)",
        ));
    }
    if num_sites > MAX_DENSE_CONSTRUCTION_SITES {
        return Err(error::invalid_param(
            "num_sites",
            "dense amplitude construction is limited to a small number of sites (2^num_sites amplitudes)",
        ));
    }
    let dim = 1usize << num_sites;
    let inv_sqrt2 = std::f64::consts::FRAC_1_SQRT_2;
    let mut amplitudes = vec![0.0_f64; dim];
    for (b, amp) in amplitudes.iter_mut().enumerate() {
        let mut a = 1.0_f64;
        for &(i, j) in &covering.dimers {
            let bit_i = (b >> i) & 1;
            let bit_j = (b >> j) & 1;
            if bit_i == bit_j {
                a = 0.0;
                break;
            } else if bit_i == 0 {
                a *= inv_sqrt2;
            } else {
                a *= -inv_sqrt2;
            }
        }
        *amp = a;
    }
    Ok(amplitudes)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring4_coverings() -> (DimerCovering, DimerCovering) {
        let a = DimerCovering::new(4, vec![(0, 1), (2, 3)]).expect("valid covering");
        let b = DimerCovering::new(4, vec![(1, 2), (3, 0)]).expect("valid covering");
        (a, b)
    }

    #[test]
    fn test_self_overlap_is_one() {
        let (a, _b) = ring4_coverings();
        let s = overlap(&a, &a, 4).expect("overlap should succeed");
        assert!(
            (s - 1.0).abs() < 1e-12,
            "self-overlap should be exactly 1, got {}",
            s
        );
    }

    #[test]
    fn test_ring4_overlap_magnitude_half() {
        let (a, b) = ring4_coverings();
        let s = overlap(&a, &b, 4).expect("overlap should succeed");
        assert!(
            (s.abs() - 0.5).abs() < 1e-12,
            "expected |S_AB| = 0.5, got {}",
            s
        );
    }

    #[test]
    fn test_overlap_matches_direct_dense_inner_product_ring4() {
        let (a, b) = overlap_direct_setup();
        let loop_overlap = overlap(&a.covering, &b.covering, 4).expect("loop overlap");
        let amp_a = a.to_dense_amplitudes().expect("dense a");
        let amp_b = b.to_dense_amplitudes().expect("dense b");
        let direct: f64 = amp_a.iter().zip(amp_b.iter()).map(|(x, y)| x * y).sum();
        assert!(
            (loop_overlap - direct).abs() < 1e-12,
            "loop-counting overlap {} should match direct dense inner product {}",
            loop_overlap,
            direct
        );
    }

    fn overlap_direct_setup() -> (ValenceBondState, ValenceBondState) {
        let (a, b) = ring4_coverings();
        (
            ValenceBondState::new(a, 4).expect("valid VB state"),
            ValenceBondState::new(b, 4).expect("valid VB state"),
        )
    }

    /// Non-bipartite "rhombus": two triangles (0,1,2) and (0,2,3) sharing edge (0,2).
    /// The extra diagonal bond (0,2) closes an odd (3-site) cycle, so this graph has no
    /// two-coloring — it directly exercises the non-bipartite reconnection sign.
    fn rhombus_bonds() -> Vec<(usize, usize)> {
        vec![(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]
    }

    #[test]
    fn test_rhombus_is_non_bipartite() {
        let coloring = super::super::dimer::bipartite_coloring(4, &rhombus_bonds())
            .expect("bipartite check should not error");
        assert!(
            coloring.is_none(),
            "rhombus graph (contains a triangle) must be non-bipartite"
        );
    }

    #[test]
    fn test_non_bipartite_bond_contribution_matches_direct_computation() {
        // Only two nearest-neighbor perfect matchings use the rhombus's edges:
        // A = {(0,1),(2,3)}, B = {(1,2),(3,0)}. The diagonal bond (0,2) is a valid
        // Hamiltonian exchange bond but participates in neither matching, so evaluating
        // it always triggers the reconnection branch.
        let a = DimerCovering::new(4, vec![(0, 1), (2, 3)]).expect("valid covering");
        let b = DimerCovering::new(4, vec![(1, 2), (3, 0)]).expect("valid covering");

        // Analytic prediction via bond_contribution/reconnect:
        let bc = bond_contribution(&b, 0, 2, 1.0).expect("bond contribution");
        let s_ab = overlap(&a, &b, 4).expect("overlap a,b");
        let (beta_prime, recon_coeff) = bc
            .reconnection
            .expect("bond (0,2) must trigger reconnection (0 and 2 are not partners in B)");
        let s_a_bp = overlap(&a, &beta_prime, 4).expect("overlap a,beta_prime");
        let predicted = bc.diagonal * s_ab + recon_coeff * s_a_bp;

        // Direct computation: apply S_0.S_2 to the dense amplitude vector of B via
        // explicit bit manipulation (S=1/2 Pauli algebra, independent of the VB
        // machinery), then inner-product with A's dense amplitudes.
        let amp_a = dense_amplitudes(&a, 4).expect("dense a");
        let amp_b = dense_amplitudes(&b, 4).expect("dense b");
        let mut applied = [0.0_f64; 16];
        for (bstate, &amp) in amp_b.iter().enumerate() {
            let bit0 = bstate & 1;
            let bit2 = (bstate >> 2) & 1;
            if bit0 == bit2 {
                applied[bstate] += 0.25 * amp;
            } else {
                applied[bstate] -= 0.25 * amp;
                let flipped = bstate ^ 1 ^ (1 << 2);
                applied[flipped] += 0.5 * amp;
            }
        }
        let direct: f64 = amp_a.iter().zip(applied.iter()).map(|(x, y)| x * y).sum();

        assert!(
            (predicted - direct).abs() < 1e-12,
            "VB-basis bond_contribution prediction {} must match direct S_0.S_2 matrix element {}",
            predicted,
            direct
        );
        // Pin the exact numeric value derived analytically (see module docs derivation).
        assert!(
            (direct - (-0.375)).abs() < 1e-12,
            "expected exact value -0.375, got {}",
            direct
        );
    }

    #[test]
    fn test_reconnect_errors_on_monomer_site() {
        let cov = DimerCovering::new(4, vec![(0, 1)]).expect("partial covering");
        let result = reconnect(&cov, 2, 3);
        assert!(result.is_err());
    }

    #[test]
    fn test_reconnect_none_when_already_partners() {
        let (a, _b) = ring4_coverings();
        let result = reconnect(&a, 0, 1).expect("should not error");
        assert!(result.is_none());
    }

    #[test]
    fn test_overlap_rejects_mismatched_monomer_sets() {
        let full = DimerCovering::new(4, vec![(0, 1), (2, 3)]).expect("full covering");
        let partial = DimerCovering::new(4, vec![(0, 1)]).expect("partial covering");
        let result = overlap(&full, &partial, 4);
        assert!(result.is_err());
    }

    #[test]
    fn test_dense_amplitudes_rejects_monomers() {
        let partial = DimerCovering::new(4, vec![(0, 1)]).expect("partial covering");
        let result = dense_amplitudes(&partial, 4);
        assert!(result.is_err());
    }
}
