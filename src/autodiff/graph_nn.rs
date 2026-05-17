//! Graph neural network potentials for spin systems on arbitrary lattices.
//!
//! Builds on the Cartesian-tensor O(3)-equivariant layers from
//! [`crate::autodiff::equivariant`] by adding **message-passing**: each spin's
//! representation is updated by aggregating equivariant messages from its
//! neighbours.  This generalises the fixed-pair `NeuralExchange` (v0.7.0) to
//! arbitrary lattice topologies while preserving SO(3) symmetry by construction.
//!
//! ## Architecture (per layer)
//!
//! For each edge `(i, j)` in the graph:
//!
//! ```text
//! message_ij = f_msg([s_i, s_j, |r_ij|], [v_i, v_j, r̂_ij])
//! ```
//!
//! where `f_msg` is an [`EquivariantLinear`] map producing scalar and vector
//! messages.  The scalar inputs concatenate the two node-scalar features with
//! the edge length `|r_ij|`; the vector inputs concatenate the two node-vector
//! features with the unit edge vector `r̂_ij`.
//!
//! Aggregate per node by summing messages from every incoming edge:
//!
//! ```text
//! m_i_s = Σ_{j ∈ N(i)} message_ij_s
//! m_i_v = Σ_{j ∈ N(i)} message_ij_v        (rotation-equivariant)
//! ```
//!
//! Then update each node with a second equivariant linear map:
//!
//! ```text
//! (s_i', v_i') = f_update([s_i, m_i_s], [v_i, m_i_v])
//! ```
//!
//! For the energy readout, sum the final-layer scalar representations:
//!
//! ```text
//! E_total = Σ_i s_i^(L)
//! ```
//!
//! ## Equivariance
//!
//! Because every node-level operation is built from [`EquivariantLinear`] and
//! summation, each individual message and the aggregated value rotate covariantly
//! with the input spins.  Scalars are computed only from rotation-invariant
//! quantities (vector magnitudes, dot products through `w_sv`), so the total
//! energy is exactly rotation-invariant up to floating-point round-off.
//!
//! ## References
//!
//! - J. Gilmer *et al.*, "Neural Message Passing for Quantum Chemistry",
//!   *ICML* (2017), arXiv:1704.01212.
//! - K. T. Schütt, O. T. Unke & M. Gastegger, "Equivariant Message Passing for
//!   the Prediction of Tensorial Properties and Molecular Spectra",
//!   *ICML* (2021), arXiv:2102.03150.

use crate::autodiff::equivariant::{EquivariantConfig, EquivariantLinear};
use crate::error::{dimension_mismatch, invalid_param, Result};
use crate::vector3::Vector3;

/// Per-node feature pair: one row of scalar features per node, plus one row of
/// vector features per node.  Used as the input and output type of
/// [`GraphMessagePassingLayer::forward`].
pub type NodeFeatures = (Vec<Vec<f64>>, Vec<Vec<Vector3<f64>>>);

// ─── LatticeGraph ────────────────────────────────────────────────────────────

/// Graph topology: an adjacency list of directed edges with their displacement
/// vectors.
///
/// Each entry is `(i, j, r_ij)` meaning a directed edge from `i` to `j` with
/// displacement vector `r_ij = r_j − r_i`.  For *undirected* graphs (the
/// usual physical case), store **both** directions: `(i, j, r_ij)` and
/// `(j, i, −r_ij)`.  The graph helpers below
/// ([`LatticeGraph::chain_1d`], [`LatticeGraph::ring_1d`],
/// [`LatticeGraph::square_lattice_2d`]) do this automatically.
#[derive(Debug, Clone)]
pub struct LatticeGraph {
    /// Number of nodes (sites) in the graph.
    pub n_nodes: usize,
    /// Directed edges `(i, j, r_ij)`.
    pub edges: Vec<(usize, usize, Vector3<f64>)>,
}

impl LatticeGraph {
    /// Build an empty graph with `n_nodes` sites and no edges.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when `n_nodes == 0`.
    pub fn new(n_nodes: usize) -> Result<Self> {
        if n_nodes == 0 {
            return Err(invalid_param("n_nodes", "must be ≥ 1"));
        }
        Ok(Self {
            n_nodes,
            edges: Vec::new(),
        })
    }

    /// Append a directed edge `(i → j)` with displacement vector `r_ij`.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when either node
    /// index is out of range.
    pub fn add_edge(&mut self, i: usize, j: usize, r_ij: Vector3<f64>) -> Result<()> {
        if i >= self.n_nodes {
            return Err(invalid_param(
                "i",
                &format!("node index {i} out of bounds for {} nodes", self.n_nodes),
            ));
        }
        if j >= self.n_nodes {
            return Err(invalid_param(
                "j",
                &format!("node index {j} out of bounds for {} nodes", self.n_nodes),
            ));
        }
        self.edges.push((i, j, r_ij));
        Ok(())
    }

    /// Return every `(j, r_ij)` neighbour of `node` (one entry per outgoing
    /// edge).  Order matches the order in which edges were inserted.
    pub fn neighbors(&self, node: usize) -> Vec<(usize, Vector3<f64>)> {
        let mut out = Vec::new();
        for &(i, j, r) in &self.edges {
            if i == node {
                out.push((j, r));
            }
        }
        out
    }

    /// Total number of (directed) edges in the graph.
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// Build a 1D open chain of `n_sites` nodes along the `x`-axis with NN bonds.
    ///
    /// The resulting graph has `2 (n_sites − 1)` directed edges (each NN bond
    /// stored in both directions).
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when
    /// `n_sites < 2` or `bond_length <= 0`.
    pub fn chain_1d(n_sites: usize, bond_length: f64) -> Result<Self> {
        if n_sites < 2 {
            return Err(invalid_param("n_sites", "chain must have ≥ 2 sites"));
        }
        if !(bond_length > 0.0 && bond_length.is_finite()) {
            return Err(invalid_param(
                "bond_length",
                "must be a positive finite real number",
            ));
        }
        let mut g = Self::new(n_sites)?;
        let dx = Vector3::new(bond_length, 0.0, 0.0);
        let nx = Vector3::new(-bond_length, 0.0, 0.0);
        for i in 0..(n_sites - 1) {
            g.add_edge(i, i + 1, dx)?;
            g.add_edge(i + 1, i, nx)?;
        }
        Ok(g)
    }

    /// Build a 1D periodic ring of `n_sites` nodes with NN bonds (each bond
    /// stored in both directions, including the closing one between site
    /// `n_sites − 1` and site `0`).
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when
    /// `n_sites < 2` or `bond_length <= 0`.
    pub fn ring_1d(n_sites: usize, bond_length: f64) -> Result<Self> {
        if n_sites < 2 {
            return Err(invalid_param("n_sites", "ring must have ≥ 2 sites"));
        }
        if !(bond_length > 0.0 && bond_length.is_finite()) {
            return Err(invalid_param(
                "bond_length",
                "must be a positive finite real number",
            ));
        }
        let mut g = Self::new(n_sites)?;
        let dx = Vector3::new(bond_length, 0.0, 0.0);
        let nx = Vector3::new(-bond_length, 0.0, 0.0);
        for i in 0..n_sites {
            let next = (i + 1) % n_sites;
            g.add_edge(i, next, dx)?;
            g.add_edge(next, i, nx)?;
        }
        Ok(g)
    }

    /// Build a 2D square lattice with `nx × ny` sites and NN bonds on the
    /// `x` and `y` axes.  Each bond is stored in both directions; the lattice
    /// is *not* periodic.
    ///
    /// Nodes are flattened with row-major indexing `node = ix · ny + iy`.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when any side is
    /// `< 1` or `a <= 0`.
    pub fn square_lattice_2d(nx: usize, ny: usize, a: f64) -> Result<Self> {
        if nx == 0 || ny == 0 {
            return Err(invalid_param("nx, ny", "must both be ≥ 1"));
        }
        if !(a > 0.0 && a.is_finite()) {
            return Err(invalid_param(
                "a",
                "lattice constant must be positive and finite",
            ));
        }
        let n_sites = nx * ny;
        let mut g = Self::new(n_sites)?;
        let dx = Vector3::new(a, 0.0, 0.0);
        let nxv = Vector3::new(-a, 0.0, 0.0);
        let dy = Vector3::new(0.0, a, 0.0);
        let nyv = Vector3::new(0.0, -a, 0.0);
        for ix in 0..nx {
            for iy in 0..ny {
                let here = ix * ny + iy;
                if ix + 1 < nx {
                    let right = (ix + 1) * ny + iy;
                    g.add_edge(here, right, dx)?;
                    g.add_edge(right, here, nxv)?;
                }
                if iy + 1 < ny {
                    let up = ix * ny + (iy + 1);
                    g.add_edge(here, up, dy)?;
                    g.add_edge(up, here, nyv)?;
                }
            }
        }
        Ok(g)
    }
}

// ─── GraphMessagePassingLayer ────────────────────────────────────────────────

/// A single message-passing layer that updates node-level scalar/vector
/// features by aggregating equivariant messages from incoming edges.
///
/// The internal layout is:
///   - `message_fn` maps the per-edge concatenated input
///     `([s_i, s_j, |r_ij|], [v_i, v_j, r̂_ij])`
///     (`2·n_scalar + 1` scalars and `2·n_vector + 1` vectors)
///     to a message
///     `([m_s], [m_v])`
///     of the same shape as the per-node features (`n_scalar`, `n_vector`).
///   - `update_fn` maps the post-aggregation per-node input
///     `([s_i, m_i_s], [v_i, m_i_v])`
///     (`2·n_scalar` scalars and `2·n_vector` vectors)
///     back to the per-node feature shape (`n_scalar`, `n_vector`).
///
/// Both blocks use [`EquivariantLinear`] from
/// [`crate::autodiff::equivariant`], which guarantees per-edge equivariance.
/// Summation over neighbours preserves equivariance because each term in the
/// sum rotates the same way.
pub struct GraphMessagePassingLayer {
    /// Per-edge equivariant message function.
    pub message_fn: EquivariantLinear,
    /// Per-node equivariant update function (consumes `[s, m_s]` and `[v, m_v]`).
    pub update_fn: EquivariantLinear,
    /// Node-feature scalar dimension.
    pub n_scalar: usize,
    /// Node-feature vector dimension.
    pub n_vector: usize,
}

impl GraphMessagePassingLayer {
    /// Build a new message-passing layer for node features of shape
    /// `(n_scalar scalars, n_vector vectors)`.
    ///
    /// The required equivariant configs are derived automatically from
    /// `n_scalar`/`n_vector`.  Sub-seeds for the two internal
    /// [`EquivariantLinear`] blocks are spaced by the golden-ratio constant.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when both
    /// `n_scalar` and `n_vector` are zero.
    pub fn new(n_scalar: usize, n_vector: usize, rng_seed: u64) -> Result<Self> {
        if n_scalar == 0 && n_vector == 0 {
            return Err(invalid_param(
                "n_scalar, n_vector",
                "at least one node-feature dimension must be non-zero",
            ));
        }
        // Message-input dims: 2·n_scalar + 1 scalars, 2·n_vector + 1 vectors.
        let msg_cfg =
            EquivariantConfig::new(2 * n_scalar + 1, 2 * n_vector + 1, n_scalar, n_vector);
        let upd_cfg = EquivariantConfig::new(2 * n_scalar, 2 * n_vector, n_scalar, n_vector);
        // The update layer also requires ≥ 1 input on each side; this is
        // satisfied because we just doubled (then added one for messages),
        // but the constructor still validates ≥ 1 output dim per non-zero
        // node side, so we re-check via EquivariantLinear::new.
        let golden = 0x9E37_79B9_7F4A_7C15_u64;
        let message_fn = EquivariantLinear::new(msg_cfg, rng_seed)?;
        let update_fn = EquivariantLinear::new(upd_cfg, rng_seed.wrapping_add(golden))?;
        Ok(Self {
            message_fn,
            update_fn,
            n_scalar,
            n_vector,
        })
    }

    /// Forward pass over the whole graph.
    ///
    /// Inputs:
    ///   - `graph`: topology + per-edge displacement vectors.
    ///   - `node_scalars`: `[n_nodes][n_scalar]`.
    ///   - `node_vectors`: `[n_nodes][n_vector]` (each `Vector3<f64>`).
    ///
    /// Returns updated `(node_scalars', node_vectors')` of the same shape.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the per-node
    /// shapes do not match the layer's configured dimensions, or when an
    /// edge references an out-of-bounds node.
    pub fn forward(
        &self,
        graph: &LatticeGraph,
        node_scalars: &[Vec<f64>],
        node_vectors: &[Vec<Vector3<f64>>],
    ) -> Result<NodeFeatures> {
        if node_scalars.len() != graph.n_nodes {
            return Err(dimension_mismatch(
                &format!("{} node scalar slots", graph.n_nodes),
                &format!("{} node scalar slots", node_scalars.len()),
            ));
        }
        if node_vectors.len() != graph.n_nodes {
            return Err(dimension_mismatch(
                &format!("{} node vector slots", graph.n_nodes),
                &format!("{} node vector slots", node_vectors.len()),
            ));
        }
        for (idx, row) in node_scalars.iter().enumerate() {
            if row.len() != self.n_scalar {
                return Err(dimension_mismatch(
                    &format!("{} scalar dims (node {idx})", self.n_scalar),
                    &format!("{} scalar dims (node {idx})", row.len()),
                ));
            }
        }
        for (idx, row) in node_vectors.iter().enumerate() {
            if row.len() != self.n_vector {
                return Err(dimension_mismatch(
                    &format!("{} vector dims (node {idx})", self.n_vector),
                    &format!("{} vector dims (node {idx})", row.len()),
                ));
            }
        }

        // Aggregator buffers initialised to zero.
        let mut agg_s: Vec<Vec<f64>> = (0..graph.n_nodes)
            .map(|_| vec![0.0_f64; self.n_scalar])
            .collect();
        let mut agg_v: Vec<Vec<Vector3<f64>>> = (0..graph.n_nodes)
            .map(|_| vec![Vector3::<f64>::zero(); self.n_vector])
            .collect();

        // Per-edge: build the message input by concatenation, then add the
        // produced message to the destination node's aggregator.
        let mut msg_scalars: Vec<f64> = Vec::with_capacity(2 * self.n_scalar + 1);
        let mut msg_vectors: Vec<Vector3<f64>> = Vec::with_capacity(2 * self.n_vector + 1);
        for &(i, j, r_ij) in &graph.edges {
            if i >= graph.n_nodes || j >= graph.n_nodes {
                return Err(dimension_mismatch(
                    &format!("edge nodes < {}", graph.n_nodes),
                    &format!("edge ({i}, {j})"),
                ));
            }
            let r_mag = r_ij.magnitude();
            // Unit edge vector — guard against degenerate zero-length edges.
            let r_hat = if r_mag > 0.0 {
                Vector3::new(r_ij.x / r_mag, r_ij.y / r_mag, r_ij.z / r_mag)
            } else {
                Vector3::<f64>::zero()
            };

            // Assemble message inputs: [s_i, s_j, |r_ij|] then [v_i, v_j, r̂_ij].
            msg_scalars.clear();
            msg_scalars.extend_from_slice(&node_scalars[i]);
            msg_scalars.extend_from_slice(&node_scalars[j]);
            msg_scalars.push(r_mag);

            msg_vectors.clear();
            msg_vectors.extend_from_slice(&node_vectors[i]);
            msg_vectors.extend_from_slice(&node_vectors[j]);
            msg_vectors.push(r_hat);

            let (m_s, m_v) = self.message_fn.forward(&msg_scalars, &msg_vectors)?;

            // Accumulate into destination node `i` (the receiver in the
            // convention "j is the source contributing a message to i").
            // Many message-passing papers use either direction; here we pick
            // the more common one used in SchNet-style aggregators: the edge
            // (i, j) carries information *from* j *into* i.
            let dst = &mut agg_s[i];
            for (a, b) in dst.iter_mut().zip(m_s.iter()) {
                *a += *b;
            }
            let dst_v = &mut agg_v[i];
            for (a, b) in dst_v.iter_mut().zip(m_v.iter()) {
                *a = *a + *b;
            }
        }

        // Per-node update: f_update([s_i, m_i_s], [v_i, m_i_v]).
        let mut new_s: Vec<Vec<f64>> = Vec::with_capacity(graph.n_nodes);
        let mut new_v: Vec<Vec<Vector3<f64>>> = Vec::with_capacity(graph.n_nodes);
        let mut upd_scalars: Vec<f64> = Vec::with_capacity(2 * self.n_scalar);
        let mut upd_vectors: Vec<Vector3<f64>> = Vec::with_capacity(2 * self.n_vector);
        for node in 0..graph.n_nodes {
            upd_scalars.clear();
            upd_scalars.extend_from_slice(&node_scalars[node]);
            upd_scalars.extend_from_slice(&agg_s[node]);
            upd_vectors.clear();
            upd_vectors.extend_from_slice(&node_vectors[node]);
            upd_vectors.extend_from_slice(&agg_v[node]);
            let (s_out, v_out) = self.update_fn.forward(&upd_scalars, &upd_vectors)?;
            new_s.push(s_out);
            new_v.push(v_out);
        }

        Ok((new_s, new_v))
    }

    /// Total number of free parameters across both internal blocks.
    pub fn n_params(&self) -> usize {
        self.message_fn.n_params() + self.update_fn.n_params()
    }

    /// Flatten every parameter in a stable order: `message_fn → update_fn`.
    pub fn params_flat(&self) -> Vec<f64> {
        let mut v = self.message_fn.params_flat();
        v.extend(self.update_fn.params_flat());
        v
    }

    /// Inverse of [`Self::params_flat`].
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the slice
    /// length does not equal [`Self::n_params`].
    pub fn set_params(&mut self, flat: &[f64]) -> Result<()> {
        let expected = self.n_params();
        if flat.len() != expected {
            return Err(dimension_mismatch(
                &format!("{expected} params"),
                &format!("{} params", flat.len()),
            ));
        }
        let nm = self.message_fn.n_params();
        self.message_fn.set_params(&flat[..nm])?;
        self.update_fn.set_params(&flat[nm..])?;
        Ok(())
    }
}

// ─── GraphMlp ────────────────────────────────────────────────────────────────

/// Stack of message-passing layers with a sum-over-sites scalar energy readout.
///
/// All layers share the *same* per-node feature shape `(n_scalar, n_vector)`,
/// so the network can be repeated as deep as desired without dimension book-
/// keeping.  A `tanh` non-linearity is applied to scalars between layers
/// (rotation-invariant); vectors flow through unmodified — their non-linear
/// expressivity is contributed by the `w_sv` block of each
/// [`EquivariantLinear`].
pub struct GraphMlp {
    /// Layers in evaluation order.
    pub layers: Vec<GraphMessagePassingLayer>,
    /// Per-node scalar feature dimension shared by every layer.
    pub n_scalar: usize,
    /// Per-node vector feature dimension shared by every layer.
    pub n_vector: usize,
}

impl GraphMlp {
    /// Build a stacked graph MLP with `n_layers` layers, all sharing
    /// `(n_scalar, n_vector)` node-feature dimensions.
    ///
    /// At least one scalar dimension is required because the energy readout
    /// reads `s_i[0]` per site.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] when
    /// `n_layers == 0`, `n_scalar == 0`, or `n_vector == 0`.
    pub fn new(n_layers: usize, n_scalar: usize, n_vector: usize, rng_seed: u64) -> Result<Self> {
        if n_layers == 0 {
            return Err(invalid_param("n_layers", "must be ≥ 1"));
        }
        if n_scalar == 0 {
            return Err(invalid_param(
                "n_scalar",
                "graph MLP requires ≥ 1 scalar feature dim for the energy readout",
            ));
        }
        if n_vector == 0 {
            return Err(invalid_param(
                "n_vector",
                "graph MLP requires ≥ 1 vector feature dim to embed the spins",
            ));
        }
        let golden = 0x9E37_79B9_7F4A_7C15_u64;
        let mut layers = Vec::with_capacity(n_layers);
        for k in 0..n_layers {
            let sub_seed = rng_seed.wrapping_add((k as u64).wrapping_mul(golden));
            layers.push(GraphMessagePassingLayer::new(n_scalar, n_vector, sub_seed)?);
        }
        Ok(Self {
            layers,
            n_scalar,
            n_vector,
        })
    }

    /// Compute the total energy `E = Σ_i s_i^{(L)}[0]` for a spin configuration.
    ///
    /// The initial per-node features are constructed so that the first vector
    /// channel stores the input spin (and is the *only* non-zero channel),
    /// while every scalar channel is initialised to zero.  This makes the
    /// network strictly equivariant from the input layer through to the
    /// readout — any rotation of the spins rotates the vector channels in
    /// the same way and leaves the scalar energy invariant.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the number of
    /// spins does not match `graph.n_nodes`.
    pub fn energy(&self, graph: &LatticeGraph, spins: &[Vector3<f64>]) -> Result<f64> {
        if spins.len() != graph.n_nodes {
            return Err(dimension_mismatch(
                &format!("{} spins", graph.n_nodes),
                &format!("{} spins", spins.len()),
            ));
        }
        // Initialise node features.
        let mut s: Vec<Vec<f64>> = (0..graph.n_nodes)
            .map(|_| vec![0.0_f64; self.n_scalar])
            .collect();
        let mut v: Vec<Vec<Vector3<f64>>> = (0..graph.n_nodes)
            .map(|node| {
                let mut row = vec![Vector3::<f64>::zero(); self.n_vector];
                row[0] = spins[node];
                row
            })
            .collect();

        let last = self.layers.len() - 1;
        for (k, layer) in self.layers.iter().enumerate() {
            let (out_s, out_v) = layer.forward(graph, &s, &v)?;
            s = out_s;
            v = out_v;
            if k != last {
                // Element-wise tanh on every scalar (rotation-invariant).
                for row in s.iter_mut() {
                    for sj in row.iter_mut() {
                        *sj = sj.tanh();
                    }
                }
            }
        }

        // Sum the first scalar feature over every site.
        let mut total = 0.0_f64;
        for row in &s {
            total += row[0];
        }
        Ok(total)
    }

    /// Total number of free parameters across every layer.
    pub fn n_params(&self) -> usize {
        self.layers.iter().map(|l| l.n_params()).sum()
    }

    /// Concatenate every layer's flat parameter vector in evaluation order.
    pub fn params_flat(&self) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.n_params());
        for layer in &self.layers {
            v.extend(layer.params_flat());
        }
        v
    }

    /// Inverse of [`Self::params_flat`].
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] when the slice
    /// length does not equal [`Self::n_params`].
    pub fn set_params(&mut self, flat: &[f64]) -> Result<()> {
        let expected = self.n_params();
        if flat.len() != expected {
            return Err(dimension_mismatch(
                &format!("{expected} params"),
                &format!("{} params", flat.len()),
            ));
        }
        let mut cursor = 0_usize;
        for layer in &mut self.layers {
            let n = layer.n_params();
            layer.set_params(&flat[cursor..cursor + n])?;
            cursor += n;
        }
        Ok(())
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::autodiff::equivariant::{random_so3, rotate_vector};

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // 1. Construct LatticeGraph + add_edge bounds validation.
    #[test]
    fn test_lattice_graph_construct() {
        assert!(LatticeGraph::new(0).is_err());
        let mut g = LatticeGraph::new(4).unwrap();
        assert_eq!(g.n_nodes, 4);
        assert_eq!(g.n_edges(), 0);
        g.add_edge(0, 1, Vector3::unit_x()).unwrap();
        assert_eq!(g.n_edges(), 1);
        // Out-of-bounds rejected.
        assert!(g.add_edge(0, 4, Vector3::unit_x()).is_err());
        assert!(g.add_edge(4, 0, Vector3::unit_x()).is_err());
    }

    // 2. chain_1d has 2·(N−1) directed edges and correct displacements.
    #[test]
    fn test_chain_1d_topology() {
        let g = LatticeGraph::chain_1d(5, 1.0).unwrap();
        assert_eq!(g.n_nodes, 5);
        assert_eq!(g.n_edges(), 2 * (5 - 1));
        let nbrs = g.neighbors(2);
        assert_eq!(nbrs.len(), 2);
        let target_pairs: Vec<usize> = nbrs.iter().map(|(j, _)| *j).collect();
        assert!(target_pairs.contains(&1));
        assert!(target_pairs.contains(&3));
        // Bad parameter rejected.
        assert!(LatticeGraph::chain_1d(1, 1.0).is_err());
        assert!(LatticeGraph::chain_1d(3, -1.0).is_err());
    }

    // 3. ring_1d has 2N directed edges (periodic closing bond).
    #[test]
    fn test_ring_1d_topology() {
        let g = LatticeGraph::ring_1d(6, 2.0).unwrap();
        assert_eq!(g.n_edges(), 2 * 6);
        // Site 0's neighbours are 1 and N − 1 = 5.
        let nbrs0: Vec<usize> = g.neighbors(0).iter().map(|(j, _)| *j).collect();
        assert!(nbrs0.contains(&1));
        assert!(nbrs0.contains(&5));
    }

    // 4. square_lattice_2d has the expected number of edges and valid neighbours.
    #[test]
    fn test_square_lattice_2d_topology() {
        let nx = 3;
        let ny = 4;
        let g = LatticeGraph::square_lattice_2d(nx, ny, 1.0).unwrap();
        assert_eq!(g.n_nodes, nx * ny);
        // Open boundaries: horizontal bonds = (nx − 1)·ny, vertical bonds = nx·(ny − 1).
        let horizontal = (nx - 1) * ny;
        let vertical = nx * (ny - 1);
        let directed = 2 * (horizontal + vertical);
        assert_eq!(g.n_edges(), directed);
        // Bad parameter rejected.
        assert!(LatticeGraph::square_lattice_2d(0, 4, 1.0).is_err());
        assert!(LatticeGraph::square_lattice_2d(3, 4, 0.0).is_err());
    }

    // 5. neighbors() returns correctly-indexed (j, r_ij) pairs.
    #[test]
    fn test_neighbors_returns_correct_list() {
        let mut g = LatticeGraph::new(3).unwrap();
        let r01 = Vector3::new(1.0, 0.0, 0.0);
        let r02 = Vector3::new(0.0, 1.0, 0.0);
        g.add_edge(0, 1, r01).unwrap();
        g.add_edge(0, 2, r02).unwrap();
        g.add_edge(1, 2, Vector3::new(-1.0, 1.0, 0.0)).unwrap();
        let n0 = g.neighbors(0);
        assert_eq!(n0.len(), 2);
        assert_eq!(n0[0].0, 1);
        assert_eq!(n0[1].0, 2);
        assert!(approx(n0[0].1.x, 1.0, 1e-15));
        assert!(approx(n0[1].1.y, 1.0, 1e-15));
        let n2 = g.neighbors(2);
        assert!(n2.is_empty());
    }

    // 6. Rotation invariance of total energy: rotating every spin and every
    //    edge vector by the same R leaves the energy unchanged.  This is the
    //    core correctness check for the architecture.
    #[test]
    fn test_rotation_invariance_of_energy() {
        let mut g = LatticeGraph::new(4).unwrap();
        let edges = [
            (0_usize, 1_usize, Vector3::new(1.0, 0.0, 0.0)),
            (1, 0, Vector3::new(-1.0, 0.0, 0.0)),
            (1, 2, Vector3::new(0.5, 0.8, 0.0)),
            (2, 1, Vector3::new(-0.5, -0.8, 0.0)),
            (2, 3, Vector3::new(0.0, 1.0, 0.3)),
            (3, 2, Vector3::new(0.0, -1.0, -0.3)),
            (3, 0, Vector3::new(-0.4, -0.7, -0.2)),
            (0, 3, Vector3::new(0.4, 0.7, 0.2)),
        ];
        for &(i, j, r) in &edges {
            g.add_edge(i, j, r).unwrap();
        }
        let net = GraphMlp::new(2, 3, 4, 4242).unwrap();
        let spins = vec![
            Vector3::new(0.6, -0.4, 0.5),
            Vector3::new(-0.3, 0.8, 0.1),
            Vector3::new(0.2, 0.2, -0.9),
            Vector3::new(0.7, 0.5, 0.2),
        ];

        let r_mat = random_so3(123);
        let rotated_spins: Vec<Vector3<f64>> =
            spins.iter().map(|s| rotate_vector(&r_mat, *s)).collect();
        let mut rotated_graph = LatticeGraph::new(g.n_nodes).unwrap();
        for &(i, j, r_ij) in &g.edges {
            rotated_graph
                .add_edge(i, j, rotate_vector(&r_mat, r_ij))
                .unwrap();
        }

        let e0 = net.energy(&g, &spins).unwrap();
        let er = net.energy(&rotated_graph, &rotated_spins).unwrap();
        assert!(
            approx(e0, er, 1e-10),
            "energy not invariant under rotation: {} vs {}",
            e0,
            er,
        );
    }

    // 7. n_params() agrees with the block-wise accounting.
    #[test]
    fn test_n_params_consistent() {
        let layer = GraphMessagePassingLayer::new(2, 3, 17).unwrap();
        let expected = layer.message_fn.n_params() + layer.update_fn.n_params();
        assert_eq!(layer.n_params(), expected);
        assert_eq!(layer.params_flat().len(), expected);
        let mlp = GraphMlp::new(3, 2, 3, 17).unwrap();
        let total: usize = mlp.layers.iter().map(|l| l.n_params()).sum();
        assert_eq!(mlp.n_params(), total);
        assert_eq!(mlp.params_flat().len(), total);
    }

    // 8. params_flat / set_params bit-for-bit round-trip.
    #[test]
    fn test_params_roundtrip() {
        let mut mlp = GraphMlp::new(2, 2, 2, 7).unwrap();
        let original = mlp.params_flat();
        let mut perturbed = original.clone();
        for (i, v) in perturbed.iter_mut().enumerate() {
            *v = 0.001 + 0.1 * (i as f64);
        }
        mlp.set_params(&perturbed).unwrap();
        let rt = mlp.params_flat();
        for (a, b) in rt.iter().zip(perturbed.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        // Mismatched-length input is rejected.
        assert!(mlp.set_params(&perturbed[..perturbed.len() - 1]).is_err());
    }

    // 9. Energy is finite for several spin configurations on a 1D ring.
    #[test]
    fn test_energy_finite_for_various_configs() {
        let g = LatticeGraph::ring_1d(6, 1.0).unwrap();
        let net = GraphMlp::new(2, 2, 2, 2024).unwrap();
        let configs: Vec<Vec<Vector3<f64>>> = vec![
            // FM along z.
            (0..6).map(|_| Vector3::unit_z()).collect(),
            // AFM along z.
            (0..6)
                .map(|i| {
                    if i % 2 == 0 {
                        Vector3::unit_z()
                    } else {
                        Vector3::unit_z() * -1.0
                    }
                })
                .collect(),
            // 120° spiral in xy plane.
            (0..6)
                .map(|i| {
                    let theta = 2.0 * std::f64::consts::PI * (i as f64) / 3.0;
                    Vector3::new(theta.cos(), theta.sin(), 0.0)
                })
                .collect(),
        ];
        for cfg in &configs {
            let e = net.energy(&g, cfg).unwrap();
            assert!(
                e.is_finite(),
                "non-finite energy for cfg with len {}",
                cfg.len()
            );
        }
    }

    // 10. Different spin configurations should yield different energies (a
    //     necessary condition for a useful gradient signal).
    #[test]
    fn test_distinct_configs_yield_distinct_energies() {
        let g = LatticeGraph::ring_1d(4, 1.0).unwrap();
        let net = GraphMlp::new(2, 3, 3, 999).unwrap();
        let fm = vec![Vector3::unit_z(); 4];
        let afm = vec![
            Vector3::unit_z(),
            Vector3::unit_z() * -1.0,
            Vector3::unit_z(),
            Vector3::unit_z() * -1.0,
        ];
        let e_fm = net.energy(&g, &fm).unwrap();
        let e_afm = net.energy(&g, &afm).unwrap();
        assert!(
            (e_fm - e_afm).abs() > 1e-8,
            "FM ({e_fm}) and AFM ({e_afm}) should differ",
        );
    }

    // 11. Multi-layer message passing composes: deeper network still passes
    //     the rotation-invariance test.
    #[test]
    fn test_multi_layer_invariance() {
        let g = LatticeGraph::square_lattice_2d(2, 3, 1.0).unwrap();
        let net = GraphMlp::new(3, 2, 3, 314).unwrap();
        let spins: Vec<Vector3<f64>> = (0..g.n_nodes)
            .map(|k| {
                let phase = (k as f64) * 0.7;
                Vector3::new(phase.cos(), phase.sin(), 0.3 * phase.cos())
            })
            .collect();
        let r_mat = random_so3(2718);
        let rotated_spins: Vec<Vector3<f64>> =
            spins.iter().map(|s| rotate_vector(&r_mat, *s)).collect();
        let mut rotated_graph = LatticeGraph::new(g.n_nodes).unwrap();
        for &(i, j, r_ij) in &g.edges {
            rotated_graph
                .add_edge(i, j, rotate_vector(&r_mat, r_ij))
                .unwrap();
        }
        let e0 = net.energy(&g, &spins).unwrap();
        let er = net.energy(&rotated_graph, &rotated_spins).unwrap();
        assert!(approx(e0, er, 1e-10));
    }

    // 12. Seed reproducibility: identical seeds produce bit-identical params.
    #[test]
    fn test_seed_reproducibility() {
        let a = GraphMlp::new(2, 2, 3, 555).unwrap();
        let b = GraphMlp::new(2, 2, 3, 555).unwrap();
        let pa = a.params_flat();
        let pb = b.params_flat();
        assert_eq!(pa.len(), pb.len());
        for (x, y) in pa.iter().zip(pb.iter()) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
        // Independently: a layer with a different seed must differ.
        let c = GraphMlp::new(2, 2, 3, 556).unwrap();
        let pc = c.params_flat();
        let any_diff = pa
            .iter()
            .zip(pc.iter())
            .any(|(x, y)| x.to_bits() != y.to_bits());
        assert!(any_diff, "different seeds should yield distinct params");
    }
}
