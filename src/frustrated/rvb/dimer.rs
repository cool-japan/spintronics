//! Dimer covering combinatorics for the resonating-valence-bond (RVB) basis.
//!
//! A *dimer covering* (equivalently a *perfect matching*) of a set of `N` sites is a
//! partition of the sites into disjoint pairs ("dimers"), each representing a singlet
//! valence bond. A *monomer-dimer covering* is the same but with exactly two designated
//! sites left unpaired ("monomers"), used to represent a pair of deconfined spinons
//! (see [`super::spinon`]).
//!
//! This module enumerates dimer coverings restricted to a caller-supplied bond list
//! (nearest-neighbor valence bonds only — the standard tractable truncation used in
//! short-range RVB theory, e.g. Liang, Doucot & Anderson, Phys. Rev. Lett. 61, 365 (1988)).
//! Enumeration uses simple backtracking, which is exact (no approximation) and fast
//! enough for the cluster sizes this module targets (covering counts bounded by
//! [`super::MAX_VB_BASIS`]).

use crate::error::{self, Result};

/// A single dimer covering (perfect matching, or monomer-dimer matching) of a finite
/// set of sites.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DimerCovering {
    /// Dimers `(i, j)` with `i < j`, sorted ascending by `i` then `j`.
    pub dimers: Vec<(usize, usize)>,
    /// `partner[i] = Some(j)` if site `i` is paired with site `j` in this covering,
    /// `None` if site `i` is an unpaired monomer.
    pub partner: Vec<Option<usize>>,
    /// Sites left unpaired (monomers), sorted ascending. Empty for a perfect matching.
    pub monomers: Vec<usize>,
}

impl DimerCovering {
    /// Build a covering from an explicit dimer list, validating that it forms a
    /// consistent (monomer-)matching of `num_sites` sites.
    ///
    /// # Errors
    ///
    /// Returns [`crate::error::Error::InvalidParameter`] if any dimer references an
    /// out-of-range site, connects a site to itself, or if any site belongs to more
    /// than one dimer.
    pub fn new(num_sites: usize, dimers: Vec<(usize, usize)>) -> Result<Self> {
        let mut partner: Vec<Option<usize>> = vec![None; num_sites];
        let mut canon_dimers: Vec<(usize, usize)> = Vec::with_capacity(dimers.len());
        for (i, j) in dimers {
            if i >= num_sites || j >= num_sites {
                return Err(error::invalid_param(
                    "dimers",
                    "dimer references a site index out of range",
                ));
            }
            if i == j {
                return Err(error::invalid_param(
                    "dimers",
                    "a dimer cannot connect a site to itself",
                ));
            }
            if partner[i].is_some() || partner[j].is_some() {
                return Err(error::invalid_param(
                    "dimers",
                    "a site cannot belong to more than one dimer in a covering",
                ));
            }
            partner[i] = Some(j);
            partner[j] = Some(i);
            canon_dimers.push((i.min(j), i.max(j)));
        }
        canon_dimers.sort_unstable();
        let monomers: Vec<usize> = (0..num_sites).filter(|&i| partner[i].is_none()).collect();
        Ok(Self {
            dimers: canon_dimers,
            partner,
            monomers,
        })
    }

    /// Number of sites this covering is defined over (paired + monomer).
    pub fn num_sites(&self) -> usize {
        self.partner.len()
    }

    /// Returns `true` if sites `i` and `j` are dimer partners in this covering.
    pub fn is_paired(&self, i: usize, j: usize) -> bool {
        self.partner.get(i).copied().flatten() == Some(j)
    }

    /// The dimer partner of site `i`, or `None` if `i` is a monomer (or out of range).
    pub fn partner_of(&self, i: usize) -> Option<usize> {
        self.partner.get(i).copied().flatten()
    }

    /// Enumerate all perfect matchings of `0..num_sites` using only the given bonds
    /// as candidate dimers.
    ///
    /// # Errors
    ///
    /// Returns an error if `num_sites` is odd (no perfect matching can exist), or if
    /// `bonds` references an out-of-range or self-connecting site.
    pub fn enumerate_perfect_matchings(
        num_sites: usize,
        bonds: &[(usize, usize)],
    ) -> Result<Vec<DimerCovering>> {
        if num_sites == 0 || num_sites % 2 != 0 {
            return Err(error::invalid_param(
                "num_sites",
                "a perfect dimer matching requires a positive, even number of sites",
            ));
        }
        let adjacency = build_adjacency(num_sites, bonds)?;
        let mut used = vec![false; num_sites];
        let mut current = Vec::with_capacity(num_sites / 2);
        let mut results = Vec::new();
        backtrack_matching(&adjacency, &mut used, &mut current, &mut results, &[]);
        results
            .into_iter()
            .map(|dimers| DimerCovering::new(num_sites, dimers))
            .collect()
    }

    /// Enumerate all monomer-dimer matchings of `0..num_sites` that leave exactly the
    /// two designated sites `monomer_a` and `monomer_b` unpaired, matching every other
    /// site using only the given bonds.
    ///
    /// # Errors
    ///
    /// Returns an error if `monomer_a == monomer_b`, either is out of range, if the
    /// remaining site count is odd, or if `bonds` is malformed.
    pub fn enumerate_monomer_dimer_matchings(
        num_sites: usize,
        bonds: &[(usize, usize)],
        monomer_a: usize,
        monomer_b: usize,
    ) -> Result<Vec<DimerCovering>> {
        if monomer_a >= num_sites || monomer_b >= num_sites {
            return Err(error::invalid_param(
                "monomer_a/monomer_b",
                "monomer site index out of range",
            ));
        }
        if monomer_a == monomer_b {
            return Err(error::invalid_param(
                "monomer_a/monomer_b",
                "the two monomer sites must be distinct",
            ));
        }
        if (num_sites - 2) % 2 != 0 {
            return Err(error::invalid_param(
                "num_sites",
                "removing the two monomer sites must leave an even number of sites to match",
            ));
        }
        let adjacency = build_adjacency(num_sites, bonds)?;
        let excluded = [monomer_a, monomer_b];
        let mut used = vec![false; num_sites];
        used[monomer_a] = true;
        used[monomer_b] = true;
        let mut current = Vec::with_capacity((num_sites - 2) / 2);
        let mut results = Vec::new();
        backtrack_matching(&adjacency, &mut used, &mut current, &mut results, &excluded);
        results
            .into_iter()
            .map(|dimers| DimerCovering::new(num_sites, dimers))
            .collect()
    }
}

/// Backtracking enumeration of perfect matchings restricted to `adjacency`.
///
/// `excluded` sites are treated as pre-matched (skipped entirely) so the same routine
/// serves both perfect and monomer-dimer enumeration; `used` must already mark them.
fn backtrack_matching(
    adjacency: &[Vec<usize>],
    used: &mut [bool],
    current: &mut Vec<(usize, usize)>,
    results: &mut Vec<Vec<(usize, usize)>>,
    excluded: &[usize],
) {
    let num_sites = used.len();
    let next = (0..num_sites).find(|&i| !used[i]);
    match next {
        None => results.push(current.clone()),
        Some(site) => {
            // Collect candidate partners first to avoid holding a borrow of `adjacency`
            // across the mutable recursion below.
            let candidates: Vec<usize> = adjacency[site]
                .iter()
                .copied()
                .filter(|&nbr| !used[nbr] && !excluded.contains(&nbr))
                .collect();
            for nbr in candidates {
                used[site] = true;
                used[nbr] = true;
                current.push((site.min(nbr), site.max(nbr)));
                backtrack_matching(adjacency, used, current, results, excluded);
                current.pop();
                used[site] = false;
                used[nbr] = false;
            }
        },
    }
}

/// Build an adjacency list from a bond list, validating indices.
fn build_adjacency(num_sites: usize, bonds: &[(usize, usize)]) -> Result<Vec<Vec<usize>>> {
    let mut adjacency = vec![Vec::new(); num_sites];
    for &(i, j) in bonds {
        if i >= num_sites || j >= num_sites {
            return Err(error::invalid_param(
                "bonds",
                "bond references a site index out of range",
            ));
        }
        if i == j {
            return Err(error::invalid_param(
                "bonds",
                "a bond cannot connect a site to itself",
            ));
        }
        adjacency[i].push(j);
        adjacency[j].push(i);
    }
    Ok(adjacency)
}

/// Attempt a 2-coloring (bipartition) of the graph defined by `bonds` via breadth-first
/// search.
///
/// Returns `Some(coloring)` with `coloring[i] in {0, 1}` if the graph is bipartite (each
/// bond connects different colors), or `None` if an odd cycle is found (non-bipartite).
/// Disconnected components are colored independently (each starting from color 0).
pub(crate) fn bipartite_coloring(
    num_sites: usize,
    bonds: &[(usize, usize)],
) -> Result<Option<Vec<i8>>> {
    let adjacency = build_adjacency(num_sites, bonds)?;
    let mut color: Vec<i8> = vec![-1; num_sites];
    for start in 0..num_sites {
        if color[start] != -1 {
            continue;
        }
        color[start] = 0;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        while let Some(u) = queue.pop_front() {
            for &v in &adjacency[u] {
                if color[v] == -1 {
                    color[v] = 1 - color[u];
                    queue.push_back(v);
                } else if color[v] == color[u] {
                    return Ok(None);
                }
            }
        }
    }
    Ok(Some(color))
}

/// Breadth-first-search graph distances from `source` over the graph defined by `bonds`.
///
/// `distances[i]` is `Some(d)` for the shortest path length in bonds from `source` to
/// site `i`, or `None` if unreachable.
pub(crate) fn bfs_distances(
    num_sites: usize,
    bonds: &[(usize, usize)],
    source: usize,
) -> Result<Vec<Option<usize>>> {
    if source >= num_sites {
        return Err(error::invalid_param(
            "source",
            "source site index out of range",
        ));
    }
    let adjacency = build_adjacency(num_sites, bonds)?;
    let mut distances: Vec<Option<usize>> = vec![None; num_sites];
    distances[source] = Some(0);
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(source);
    while let Some(u) = queue.pop_front() {
        let du = distances[u].unwrap_or(0);
        for &v in &adjacency[u] {
            if distances[v].is_none() {
                distances[v] = Some(du + 1);
                queue.push_back(v);
            }
        }
    }
    Ok(distances)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 4-site ring: 0-1-2-3-0.
    fn ring4_bonds() -> Vec<(usize, usize)> {
        vec![(0, 1), (1, 2), (2, 3), (3, 0)]
    }

    /// 2x4 ladder: rungs (0,4),(1,5),(2,6),(3,7); rails 0-1-2-3 and 4-5-6-7.
    fn ladder_2x4_bonds() -> Vec<(usize, usize)> {
        vec![
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
        ]
    }

    /// 4x4 open (non-periodic) square lattice, sites numbered row-major.
    fn square_4x4_bonds() -> Vec<(usize, usize)> {
        let idx = |x: usize, y: usize| -> usize { y * 4 + x };
        let mut bonds = Vec::new();
        for y in 0..4 {
            for x in 0..4 {
                if x + 1 < 4 {
                    bonds.push((idx(x, y), idx(x + 1, y)));
                }
                if y + 1 < 4 {
                    bonds.push((idx(x, y), idx(x, y + 1)));
                }
            }
        }
        bonds
    }

    #[test]
    fn test_dimer_covering_new_validates_double_use() {
        let result = DimerCovering::new(4, vec![(0, 1), (1, 2)]);
        assert!(result.is_err());
    }

    #[test]
    fn test_dimer_covering_new_self_loop_errors() {
        let result = DimerCovering::new(4, vec![(0, 0)]);
        assert!(result.is_err());
    }

    #[test]
    fn test_dimer_covering_monomers() {
        let cov = DimerCovering::new(4, vec![(0, 1)]).expect("valid partial covering");
        assert_eq!(cov.monomers, vec![2, 3]);
        assert!(cov.is_paired(0, 1));
        assert!(cov.is_paired(1, 0));
        assert!(!cov.is_paired(0, 2));
        assert_eq!(cov.partner_of(2), None);
    }

    #[test]
    fn test_ring4_covering_count_is_two() {
        let coverings = DimerCovering::enumerate_perfect_matchings(4, &ring4_bonds())
            .expect("ring4 enumeration should succeed");
        assert_eq!(
            coverings.len(),
            2,
            "4-ring should have exactly 2 nearest-neighbor dimer coverings"
        );
        for cov in &coverings {
            assert!(cov.monomers.is_empty());
        }
    }

    #[test]
    fn test_ladder_2x4_covering_count_is_five() {
        let coverings = DimerCovering::enumerate_perfect_matchings(8, &ladder_2x4_bonds())
            .expect("ladder enumeration should succeed");
        assert_eq!(
            coverings.len(),
            5,
            "2x4 ladder should have exactly 5 nearest-neighbor dimer coverings"
        );
    }

    #[test]
    fn test_square_4x4_covering_count_is_36() {
        let coverings = DimerCovering::enumerate_perfect_matchings(16, &square_4x4_bonds())
            .expect("4x4 open square enumeration should succeed");
        assert_eq!(
            coverings.len(),
            36,
            "4x4 open square should have exactly 36 nearest-neighbor dimer coverings"
        );
    }

    #[test]
    fn test_odd_num_sites_errors() {
        let result = DimerCovering::enumerate_perfect_matchings(3, &[(0, 1), (1, 2)]);
        assert!(result.is_err());
    }

    #[test]
    fn test_monomer_dimer_ring4() {
        // Leaving sites 0 and 2 (opposite corners) as monomers: the remaining sites
        // 1,3 are not directly bonded on the ring, so no matching should exist.
        let coverings = DimerCovering::enumerate_monomer_dimer_matchings(4, &ring4_bonds(), 0, 2)
            .expect("enumeration should succeed even if empty");
        assert!(coverings.is_empty());

        // Leaving sites 0 and 1 (adjacent) as monomers: the remaining sites 2,3 ARE
        // bonded, so exactly one matching should exist: {(2,3)}.
        let coverings = DimerCovering::enumerate_monomer_dimer_matchings(4, &ring4_bonds(), 0, 1)
            .expect("enumeration should succeed");
        assert_eq!(coverings.len(), 1);
        assert_eq!(coverings[0].monomers, vec![0, 1]);
        assert_eq!(coverings[0].dimers, vec![(2, 3)]);
    }

    #[test]
    fn test_monomer_dimer_invalid_same_site() {
        let result = DimerCovering::enumerate_monomer_dimer_matchings(4, &ring4_bonds(), 1, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_bipartite_coloring_ring4() {
        let coloring = bipartite_coloring(4, &ring4_bonds())
            .expect("bipartite check should not error")
            .expect("4-ring is bipartite");
        assert_ne!(coloring[0], coloring[1]);
        assert_ne!(coloring[1], coloring[2]);
        assert_eq!(coloring[0], coloring[2]);
    }

    #[test]
    fn test_bipartite_coloring_triangle_is_none() {
        let triangle_bonds = vec![(0, 1), (1, 2), (2, 0)];
        let coloring = bipartite_coloring(3, &triangle_bonds).expect("check should not error");
        assert!(
            coloring.is_none(),
            "a triangle (odd cycle) is not bipartite"
        );
    }

    #[test]
    fn test_bfs_distances_ring4() {
        let distances = bfs_distances(4, &ring4_bonds(), 0).expect("bfs should succeed");
        assert_eq!(distances[0], Some(0));
        assert_eq!(distances[1], Some(1));
        assert_eq!(distances[2], Some(2));
        assert_eq!(distances[3], Some(1));
    }

    #[test]
    fn test_bfs_distances_disconnected() {
        // Sites 0-1 connected, site 2 isolated.
        let distances = bfs_distances(3, &[(0, 1)], 0).expect("bfs should succeed");
        assert_eq!(distances[0], Some(0));
        assert_eq!(distances[1], Some(1));
        assert_eq!(distances[2], None);
    }

    #[test]
    fn test_build_adjacency_out_of_range_errors() {
        let result = DimerCovering::enumerate_perfect_matchings(4, &[(0, 5)]);
        assert!(result.is_err());
    }
}
