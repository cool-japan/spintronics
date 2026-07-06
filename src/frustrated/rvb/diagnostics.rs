//! Spin-liquid diagnostics: real-space correlations, structure factor, and a
//! long-range-order screening report.
//!
//! These operate on a precomputed [`GroundState`] (from [`RvbSolver::ground_state`]) so
//! that sweeping many site pairs (e.g. for [`structure_factor`]) does not repeat the
//! `O(D^3)` Löwdin canonical orthogonalization for every pair.

use super::dimer::bfs_distances;
use super::solver::{build_hamiltonian_matrix, GroundState, RvbSolver};
use crate::error::{self, Result};
use crate::vector3::Vector3;

/// Ground-state expectation value `⟨S_i·S_j⟩` for an arbitrary pair of sites (not
/// necessarily a bonded pair), evaluated in the valence-bond basis via the same
/// Heisenberg matrix-element machinery used to build the full Hamiltonian, restricted
/// to the single "bond" `(i, j)` with unit coupling.
///
/// # Errors
///
/// Returns an error if `i` or `j` is out of range, or if `ground`'s coefficient vector
/// does not match the solver's basis dimension.
pub fn spin_correlation(
    solver: &RvbSolver,
    ground: &GroundState,
    i: usize,
    j: usize,
) -> Result<f64> {
    if i >= solver.num_sites || j >= solver.num_sites {
        return Err(error::invalid_param("i/j", "site index out of range"));
    }
    if ground.coefficients.len() != solver.dim() {
        return Err(error::invalid_param(
            "ground",
            "ground-state coefficient length does not match the solver's basis dimension",
        ));
    }
    if i == j {
        // S_i . S_i = S(S+1) = 3/4 for spin-1/2, independent of the state.
        return Ok(0.75);
    }

    let single_bond = [(i, j)];
    let m = build_hamiltonian_matrix(&solver.coverings, &single_bond, 1.0, solver.num_sites)?;
    let s = solver.overlap_matrix()?;
    let d = solver.dim();
    let c = &ground.coefficients;
    let mut numerator = 0.0_f64;
    let mut denominator = 0.0_f64;
    for a in 0..d {
        for b in 0..d {
            numerator += c[a] * m.get(a, b).re * c[b];
            denominator += c[a] * s.get(a, b).re * c[b];
        }
    }
    if denominator.abs() < 1e-14 {
        return Err(error::numerical_error(
            "ground state has (numerically) zero norm: c^T S c vanished",
        ));
    }
    Ok(numerator / denominator)
}

/// Static structure factor `S(q) = (1/N) Σ_{i,j} cos(q·(r_i − r_j)) ⟨S_i·S_j⟩` at wave
/// vector `q`, given explicit real-space site positions (the abstract, geometry-free
/// [`RvbSolver`] does not itself carry positions — callers using
/// [`RvbSolver::from_lattice`] can supply `lattice.positions` directly).
///
/// # Errors
///
/// Returns an error if `positions.len() != solver.num_sites`, or propagates errors from
/// [`spin_correlation`].
pub fn structure_factor(
    solver: &RvbSolver,
    ground: &GroundState,
    positions: &[Vector3<f64>],
    q: Vector3<f64>,
) -> Result<f64> {
    if positions.len() != solver.num_sites {
        return Err(error::invalid_param(
            "positions",
            "positions length must equal solver.num_sites",
        ));
    }
    let n = solver.num_sites;
    let mut total = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            let corr = spin_correlation(solver, ground, i, j)?;
            let dr = positions[i] - positions[j];
            let phase = q.dot(&dr);
            total += corr * phase.cos();
        }
    }
    Ok(total / n as f64)
}

/// Screening report for spin-liquid behavior: absence of long-range magnetic order,
/// diagnosed via the ground-state spin correlation between the two graph-theoretically
/// most distant sites.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpinLiquidReport {
    /// Valence-bond ground-state energy.
    pub ground_energy: f64,
    /// First site of the most graph-distant pair found.
    pub far_site_a: usize,
    /// Second site of the most graph-distant pair found.
    pub far_site_b: usize,
    /// Graph (bond-count) distance between `far_site_a` and `far_site_b`.
    pub far_graph_distance: usize,
    /// `|⟨S_far_a · S_far_b⟩|` for the ground state.
    pub long_range_correlation: f64,
    /// `true` if `long_range_correlation` is below the supplied threshold — a candidate
    /// signature of a spin liquid (no long-range order) rather than conventional Néel
    /// order (which would show an O(1), non-decaying long-range correlation).
    pub is_candidate: bool,
}

/// Screen an [`RvbSolver`]'s ground state for spin-liquid-like behavior (see
/// [`SpinLiquidReport`]).
///
/// # Errors
///
/// Returns an error if `correlation_threshold <= 0`, or propagates errors from
/// [`RvbSolver::ground_state`] / [`spin_correlation`].
pub fn is_spin_liquid(solver: &RvbSolver, correlation_threshold: f64) -> Result<SpinLiquidReport> {
    if correlation_threshold <= 0.0 {
        return Err(error::invalid_param(
            "correlation_threshold",
            "correlation threshold must be positive",
        ));
    }
    let ground = solver.ground_state()?;
    let n = solver.num_sites;

    let mut far_pair = (0usize, 0usize);
    let mut max_dist = 0usize;
    for source in 0..n {
        let distances = bfs_distances(n, &solver.bonds, source)?;
        for (target, dist_opt) in distances.into_iter().enumerate() {
            if let Some(dist) = dist_opt {
                if dist > max_dist {
                    max_dist = dist;
                    far_pair = (source, target);
                }
            }
        }
    }

    let long_range_correlation = if max_dist == 0 {
        0.0
    } else {
        spin_correlation(solver, &ground, far_pair.0, far_pair.1)?.abs()
    };
    let is_candidate = long_range_correlation < correlation_threshold;

    Ok(SpinLiquidReport {
        ground_energy: ground.energy,
        far_site_a: far_pair.0,
        far_site_b: far_pair.1,
        far_graph_distance: max_dist,
        long_range_correlation,
        is_candidate,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring4_solver() -> RvbSolver {
        RvbSolver::from_bonds(4, vec![(0, 1), (1, 2), (2, 3), (3, 0)], 1.0)
            .expect("valid ring4 solver")
    }

    #[test]
    fn test_self_correlation_is_three_quarters() {
        let solver = ring4_solver();
        let ground = solver.ground_state().expect("ground state");
        let c = spin_correlation(&solver, &ground, 2, 2).expect("self correlation");
        assert!((c - 0.75).abs() < 1e-12);
    }

    #[test]
    fn test_ring4_nearest_neighbor_correlation_exact() {
        let solver = ring4_solver();
        let ground = solver.ground_state().expect("ground state");
        // Ground energy = J * sum_over_4_bonds <S_i.S_j> = 4 * corr_nn (by symmetry) = -2.0.
        let corr = spin_correlation(&solver, &ground, 0, 1).expect("nn correlation");
        assert!(
            (corr - (-0.5)).abs() < 1e-9,
            "expected <S_0.S_1> = -0.5, got {}",
            corr
        );
    }

    #[test]
    fn test_ring4_diagonal_correlation_matches_singlet_constraint() {
        let solver = ring4_solver();
        let ground = solver.ground_state().expect("ground state");
        // S_tot^2 = 0 for the exact ring4 singlet ground state implies (see module/
        // solver docs derivation) <S_0.S_2> = 1/4 exactly.
        let corr = spin_correlation(&solver, &ground, 0, 2).expect("diagonal correlation");
        assert!(
            (corr - 0.25).abs() < 1e-9,
            "expected <S_0.S_2> = 0.25, got {}",
            corr
        );
    }

    #[test]
    fn test_structure_factor_at_q_zero_matches_zero_total_spin() {
        let solver = ring4_solver();
        let ground = solver.ground_state().expect("ground state");
        // At q=0, cos(q.dr)=1 regardless of positions, so S(0) = <S_tot^2>/N = 0 for
        // the exact singlet ground state.
        let positions = vec![Vector3::zero(); 4];
        let s0 = structure_factor(&solver, &ground, &positions, Vector3::zero())
            .expect("structure factor should solve");
        assert!(
            s0.abs() < 1e-9,
            "expected S(q=0) = 0 for a total-spin-0 ground state, got {}",
            s0
        );
    }

    #[test]
    fn test_structure_factor_rejects_wrong_position_count() {
        let solver = ring4_solver();
        let ground = solver.ground_state().expect("ground state");
        let positions = vec![Vector3::zero(); 3];
        assert!(structure_factor(&solver, &ground, &positions, Vector3::zero()).is_err());
    }

    #[test]
    fn test_is_spin_liquid_report_ring4() {
        let solver = ring4_solver();
        let report = is_spin_liquid(&solver, 0.5).expect("report should solve");
        assert_eq!(
            report.far_graph_distance, 2,
            "ring4 diameter is 2 (opposite corners)"
        );
        assert!((report.long_range_correlation - 0.25).abs() < 1e-9);
        assert!((report.ground_energy - (-2.0)).abs() < 1e-9);
        // 0.25 < 0.5 threshold => flagged as a candidate.
        assert!(report.is_candidate);
    }

    #[test]
    fn test_is_spin_liquid_rejects_non_positive_threshold() {
        let solver = ring4_solver();
        assert!(is_spin_liquid(&solver, 0.0).is_err());
        assert!(is_spin_liquid(&solver, -0.1).is_err());
    }
}
