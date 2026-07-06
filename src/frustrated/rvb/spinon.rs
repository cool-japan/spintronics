//! Spinon-pair energetics and deconfinement diagnostics.
//!
//! A pair of well-separated **spinons** is represented here by a monomer-dimer
//! covering that leaves exactly two designated sites `r1`, `r2` unpaired, with the
//! remaining sites in a valence-bond dimer covering. The two monomer spins are taken in
//! the fully polarized reference state `|↑_{r1}↑_{r2}⟩` (one specific, exactly
//! computable component of what would be a triplet if the two sites were combined),
//! tensored with the singlet background.
//!
//! This choice makes every bond touching a monomer site **exactly** (not just
//! approximately) tractable: for a bond `(r1, k)` with `k` part of a singlet dimer
//! `(k, l)`, `⟨↑_{r1}| S_{r1}^z |↑_{r1}⟩ = 1/2` while `⟨s_{kl}| S_k^a |s_{kl}⟩ = 0` for
//! every spin component `a` (a singlet has zero on-site spin expectation), and the
//! raising/lowering cross-terms vanish because `S_{r1}^+|↑⟩ = 0`. So
//! `⟨ψ|S_{r1}·S_k|ψ⟩ = 0` exactly for this trial-state family — bonds from a monomer to
//! a dimer-covered site simply drop out. A *direct* bond between the two monomers
//! themselves, `(r1, r2)`, is not eliminated this way: both spins are pinned up, so
//! `⟨S_{r1}·S_{r2}⟩ = (1/2)(1/2) = 1/4` exactly (raising/lowering annihilate the fully
//! polarized pair). Both facts are exact consequences of the trial-state symmetry, not
//! approximations.
//!
//! [`spinon_pair_energy`] therefore reduces to: (i) the variational valence-bond ground
//! energy of the `(N-2)`-site sub-lattice (bonds not touching either monomer, reusing
//! `super::solver::generalized_ground_state`), plus (ii) `+J/4` if `(r1,r2)` is itself
//! a bond. [`deconfinement_diagnostic`] sweeps a list of site pairs and reports the
//! energy cost `Δ(r) = E_spinon_pair(r) − E_ground` relative to the true (no-spinon) VB
//! ground state as a function of graph-theoretic separation — the standard diagnostic
//! for spinon confinement (Δ growing without bound) versus deconfinement (Δ saturating).

use super::dimer::DimerCovering;
use super::solver::{generalized_ground_state, RvbSolver};
use crate::error::{self, Error, Result};

/// One point on the spinon-pair energy-cost curve.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpinonSeparationPoint {
    /// First spinon site.
    pub site_a: usize,
    /// Second spinon site.
    pub site_b: usize,
    /// Graph (bond-count) distance between the two sites.
    pub graph_distance: usize,
    /// Energy cost `E_spinon_pair − E_ground` of separating the spinon pair to this
    /// distance, relative to the (no-spinon) valence-bond ground state.
    pub delta_energy: f64,
}

/// Variational energy of a spinon pair fixed at sites `r1`, `r2` (see module docs for
/// the exact trial-state construction and why it is tractable).
///
/// # Errors
///
/// Returns an error if `r1 == r2`, either is out of range, no monomer-dimer matching of
/// the remaining sites exists, or the monomer-dimer basis exceeds
/// [`super::MAX_VB_BASIS`] (named with the actual measured count).
pub fn spinon_pair_energy(solver: &RvbSolver, r1: usize, r2: usize) -> Result<f64> {
    if r1 >= solver.num_sites || r2 >= solver.num_sites {
        return Err(error::invalid_param(
            "r1/r2",
            "spinon site index out of range",
        ));
    }
    if r1 == r2 {
        return Err(error::invalid_param(
            "r1/r2",
            "the two spinon sites must be distinct",
        ));
    }

    let coverings =
        DimerCovering::enumerate_monomer_dimer_matchings(solver.num_sites, &solver.bonds, r1, r2)?;
    if coverings.is_empty() {
        return Err(error::invalid_param(
            "r1/r2",
            "no monomer-dimer matching exists leaving exactly these two sites unpaired",
        ));
    }
    if coverings.len() > super::MAX_VB_BASIS {
        return Err(Error::InvalidParameter {
            param: "r1/r2".to_string(),
            reason: format!(
                "monomer-dimer covering enumeration produced {} basis states, exceeding \
                 MAX_VB_BASIS={} (bounded by CMatrix::MAX_DIM=64)",
                coverings.len(),
                super::MAX_VB_BASIS
            ),
        });
    }

    // Bonds touching either monomer site are exactly zero in expectation for this trial
    // state (see module docs) and are excluded from the sub-lattice Hamiltonian; a
    // direct bond between the two monomers is the one exception, contributing a fixed
    // classical +J/4.
    let filtered_bonds: Vec<(usize, usize)> = solver
        .bonds
        .iter()
        .copied()
        .filter(|&(i, j)| i != r1 && i != r2 && j != r1 && j != r2)
        .collect();
    let direct_bond_correction = if solver
        .bonds
        .iter()
        .any(|&(i, j)| (i == r1 && j == r2) || (i == r2 && j == r1))
    {
        0.25 * solver.coupling_j
    } else {
        0.0
    };

    let ground = generalized_ground_state(
        solver.num_sites,
        &coverings,
        &filtered_bonds,
        solver.coupling_j,
    )?;
    Ok(ground.energy + direct_bond_correction)
}

/// Sweep a list of spinon-separation pairs and report the `Δ(r)` energy-cost curve
/// relative to the true (no-spinon) valence-bond ground state.
///
/// # Errors
///
/// Propagates errors from [`spinon_pair_energy`] and [`RvbSolver::ground_state`], plus
/// an error if a requested pair is disconnected in the bond graph (undefined
/// separation).
pub fn deconfinement_diagnostic(
    solver: &RvbSolver,
    pairs: &[(usize, usize)],
) -> Result<Vec<SpinonSeparationPoint>> {
    let ground = solver.ground_state()?;
    let mut points = Vec::with_capacity(pairs.len());
    for &(r1, r2) in pairs {
        let e_pair = spinon_pair_energy(solver, r1, r2)?;
        let distances = super::dimer::bfs_distances(solver.num_sites, &solver.bonds, r1)?;
        let graph_distance = distances.get(r2).copied().flatten().ok_or_else(|| {
            error::invalid_param(
                "pairs",
                "the two spinon sites must be connected in the bond graph",
            )
        })?;
        points.push(SpinonSeparationPoint {
            site_a: r1,
            site_b: r2,
            graph_distance,
            delta_energy: e_pair - ground.energy,
        });
    }
    Ok(points)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring4_solver() -> RvbSolver {
        RvbSolver::from_bonds(4, vec![(0, 1), (1, 2), (2, 3), (3, 0)], 1.0)
            .expect("valid ring4 solver")
    }

    #[test]
    fn test_ring4_adjacent_spinon_pair_energy_exact() {
        let solver = ring4_solver();
        // Hand-derived: sub-lattice has a single monomer-dimer covering {(2,3)} with
        // monomers {0,1}; diagonal energy -3/4*J (bond (2,3) is a shared dimer); plus
        // the direct-bond correction +J/4 since (0,1) is itself a bond. Total: -0.5.
        let e = spinon_pair_energy(&solver, 0, 1).expect("adjacent spinon pair should solve");
        assert!(
            (e - (-0.5)).abs() < 1e-9,
            "expected exactly -0.5, got {}",
            e
        );

        // By the ring's 4-fold symmetry, every adjacent pair gives the same energy.
        let e12 = spinon_pair_energy(&solver, 1, 2).expect("adjacent pair (1,2)");
        let e23 = spinon_pair_energy(&solver, 2, 3).expect("adjacent pair (2,3)");
        let e30 = spinon_pair_energy(&solver, 3, 0).expect("adjacent pair (3,0)");
        for e_other in [e12, e23, e30] {
            assert!(
                (e_other - e).abs() < 1e-9,
                "ring symmetry should give identical energies"
            );
        }
    }

    #[test]
    fn test_ring4_opposite_corners_have_no_matching() {
        let solver = ring4_solver();
        // Sites 1 and 3 are not directly bonded on the ring, so leaving 0 and 2 as
        // monomers leaves no valid nearest-neighbor matching for {1,3}.
        let result = spinon_pair_energy(&solver, 0, 2);
        assert!(
            result.is_err(),
            "opposite corners should have no valid monomer-dimer matching"
        );
    }

    #[test]
    fn test_spinon_same_site_errors() {
        let solver = ring4_solver();
        assert!(spinon_pair_energy(&solver, 1, 1).is_err());
    }

    #[test]
    fn test_spinon_out_of_range_errors() {
        let solver = ring4_solver();
        assert!(spinon_pair_energy(&solver, 0, 99).is_err());
    }

    #[test]
    fn test_deconfinement_diagnostic_ring4_adjacent_pairs() {
        let solver = ring4_solver();
        let pairs = [(0, 1), (1, 2), (2, 3), (3, 0)];
        let points = deconfinement_diagnostic(&solver, &pairs).expect("diagnostic should solve");
        assert_eq!(points.len(), 4);
        for p in &points {
            assert_eq!(p.graph_distance, 1);
            assert!(p.delta_energy.is_finite());
            // Ground energy is exactly -2.0; adjacent spinon-pair energy is exactly
            // -0.5, so delta must be exactly 1.5.
            assert!(
                (p.delta_energy - 1.5).abs() < 1e-8,
                "expected delta_energy = 1.5, got {}",
                p.delta_energy
            );
        }
    }

    #[test]
    fn test_deconfinement_diagnostic_ladder() {
        let bonds = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (4, 5),
            (5, 6),
            (6, 7),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        let solver = RvbSolver::from_bonds(8, bonds, 1.0).expect("valid ladder solver");
        // Rail-adjacent pair (0,1) (distance 1) and the diagonal corner pair (0,3)
        // (distance 3, shortest path 0-1-2-3) both admit a valid monomer-dimer matching
        // of the remaining 6 sites on the ladder (unlike e.g. (0,2), which strands site 4
        // with no available partner once bonds touching 0 and 2 are excluded).
        let points =
            deconfinement_diagnostic(&solver, &[(0, 1), (0, 3)]).expect("diagnostic should solve");
        assert_eq!(points.len(), 2);
        assert_eq!(points[0].graph_distance, 1);
        assert_eq!(points[1].graph_distance, 3);
        for p in &points {
            assert!(p.delta_energy.is_finite());
        }
    }
}
