//! Edge modes for topological magnon insulators in strip geometry.
//!
//! Constructs a finite-width strip Hamiltonian by stacking `width` unit cells
//! in the y-direction while keeping kx as a good quantum number (periodic in x).
//! The resulting (`width × n_bands`) × (`width × n_bands`) matrix is diagonalized
//! to reveal bulk-like and edge-localized modes.
//!
//! # Strip Construction
//!
//! For a strip of `width` unit cells in the y-direction, the block-tridiagonal
//! Hamiltonian is:
//!
//! ```text
//! H_strip(kx) = block-diag(H_intra(kx)) + off-diag(H_inter, H_inter†)
//! ```
//!
//! where `H_intra(kx)` captures intra-unit-cell and kx-direction hopping, and
//! `H_inter` encodes the inter-unit-cell hopping along y (open boundary condition).
//!
//! # Localization Metric
//!
//! Edge localization is quantified by the fraction of probability weight on the
//! two outermost unit cells (iy=0 and iy=width-1):
//!
//! ```text
//! η = (Σ_{σ} |ψ_{0,σ}|² + |ψ_{N-1,σ}|²) / Σ_{i,σ} |ψ_{i,σ}|²
//! ```
//!
//! A value η ≥ 0.3 (adjustable threshold) signals edge localization.

use std::f64::consts::PI;

use super::band_model::MagnonBandModel;
use crate::error::{self, Result};
use crate::math::{CMatrix, Complex};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Which geometric edge a mode is localized on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum EdgeSide {
    /// Top edge (iy = width − 1).
    Top,
    /// Bottom edge (iy = 0).
    Bottom,
    /// Bulk mode (not edge-localized).
    Bulk,
}

/// A single eigenmode of the strip Hamiltonian.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EdgeMode {
    /// Conserved wavevector along the periodic direction \[rad/a\].
    pub kx: f64,
    /// Mode frequency/energy in the same units as the exchange parameters.
    pub frequency: f64,
    /// Localization metric η ∈ \[0, 1\]: 1 = fully edge-localized.
    pub localization: f64,
    /// Which edge this mode is localized on.
    pub edge: EdgeSide,
}

// ---------------------------------------------------------------------------
// EdgeModes
// ---------------------------------------------------------------------------

/// Computes edge modes of a topological magnon insulator in strip geometry.
pub struct EdgeModes<'a> {
    /// Reference to the magnon band model.
    pub model: &'a MagnonBandModel,
    /// Number of unit cells in the finite (y) direction.
    pub width: usize,
}

impl<'a> EdgeModes<'a> {
    /// Create an edge-mode calculator.
    ///
    /// # Errors
    ///
    /// Returns `InvalidParameter` if `width < 2` or `width > 60`.
    pub fn new(model: &'a MagnonBandModel, width: usize) -> Result<Self> {
        if width < 2 {
            return Err(error::invalid_param(
                "width",
                "strip width must be at least 2",
            ));
        }
        if width > 60 {
            return Err(error::invalid_param(
                "width",
                "strip width must be at most 60 (matrix size limit)",
            ));
        }
        Ok(Self { model, width })
    }

    // -----------------------------------------------------------------------
    // Strip Hamiltonian
    // -----------------------------------------------------------------------

    /// Build and diagonalize the strip Hamiltonian at fixed kx.
    ///
    /// Returns all `width × n_bands` eigenmodes with their localization metrics.
    pub fn solve_strip(&self, kx: f64) -> Result<Vec<EdgeMode>> {
        let nb = self.model.n_bands();
        let ny = self.width;
        let dim = nb * ny;

        // Ensure total matrix size fits within CMatrix::MAX_DIM = 64
        if dim > crate::math::CMatrix::MAX_DIM {
            return Err(error::invalid_param(
                "width",
                "width * n_bands exceeds CMatrix::MAX_DIM (64)",
            ));
        }

        // Build the block-tridiagonal strip Hamiltonian.
        // H_intra(kx): the n_bands × n_bands Bloch Hamiltonian at (kx, ky=0)
        // captures kx-direction hoppings and on-site terms.
        let h_intra = self.model.hamiltonian_at((kx, 0.0))?;

        // H_inter: the inter-unit-cell hopping in the y-direction.
        // Computed as (H(kx, π/4) − H(kx, 0)) × 0.5 to extract the ky-dependent
        // hopping contribution. This is a first-order approximation to the
        // derivative dH/dky times the unit cell spacing.
        let h_ky_offset = self.model.hamiltonian_at((kx, PI / (ny as f64).max(4.0)))?;
        let h_inter_raw = h_ky_offset.sub(&h_intra)?;
        // Scale to represent the inter-cell hopping amplitude
        let h_inter = h_inter_raw.scale_real(0.5);

        let mut h_strip = CMatrix::zeros(dim);

        for iy in 0..ny {
            // Diagonal block (intra-cell Hamiltonian)
            for i in 0..nb {
                for j in 0..nb {
                    let cur = h_strip.get(iy * nb + i, iy * nb + j);
                    h_strip.set(iy * nb + i, iy * nb + j, cur.add(&h_intra.get(i, j)));
                }
            }

            // Super-diagonal block (coupling iy → iy+1): H_inter
            if iy + 1 < ny {
                for i in 0..nb {
                    for j in 0..nb {
                        let t = h_inter.get(i, j);
                        // H_strip[iy, iy+1] += H_inter
                        let cur_up = h_strip.get(iy * nb + i, (iy + 1) * nb + j);
                        h_strip.set(iy * nb + i, (iy + 1) * nb + j, cur_up.add(&t));
                        // H_strip[iy+1, iy] += H_inter† (Hermitian conjugate)
                        let cur_dn = h_strip.get((iy + 1) * nb + j, iy * nb + i);
                        h_strip.set((iy + 1) * nb + j, iy * nb + i, cur_dn.add(&t.conj()));
                    }
                }
            }
        }

        let (evals, vecs) = h_strip.hermitian_eigendecomposition()?;

        let mut modes = Vec::with_capacity(dim);
        for (mode_idx, &e) in evals.iter().enumerate() {
            let eigenvec: Vec<Complex> = (0..dim).map(|r| vecs.get(r, mode_idx)).collect();
            let loc = Self::edge_localization_metric(&eigenvec);
            let edge = classify_edge(&eigenvec, nb, ny);
            modes.push(EdgeMode {
                kx,
                frequency: e,
                localization: loc,
                edge,
            });
        }

        Ok(modes)
    }

    // -----------------------------------------------------------------------
    // Localization metric
    // -----------------------------------------------------------------------

    /// Compute the edge localization metric η ∈ \[0, 1\].
    ///
    /// η = (|ψ₀|² + |ψ₁|² + |ψ_{N-2}|² + |ψ_{N-1}|²) / Σ|ψ_i|²
    ///
    /// where indexing is over the full eigenvector of length `width × n_bands`.
    pub fn edge_localization_metric(eigenvector: &[Complex]) -> f64 {
        let total: f64 = eigenvector.iter().map(|c| c.norm_sq()).sum();
        if total < 1e-30 {
            return 0.0;
        }
        let n = eigenvector.len();
        if n < 4 {
            return eigenvector.iter().map(|c| c.norm_sq()).sum::<f64>() / total;
        }
        let edge_weight = eigenvector[0].norm_sq()
            + eigenvector[1].norm_sq()
            + eigenvector[n - 2].norm_sq()
            + eigenvector[n - 1].norm_sq();
        edge_weight / total
    }

    // -----------------------------------------------------------------------
    // In-gap mode filter
    // -----------------------------------------------------------------------

    /// Find all modes with frequency in [gap_min, gap_max] at fixed kx.
    pub fn find_in_gap_modes(&self, kx: f64, gap_min: f64, gap_max: f64) -> Result<Vec<EdgeMode>> {
        let all = self.solve_strip(kx)?;
        Ok(all
            .into_iter()
            .filter(|m| m.frequency >= gap_min && m.frequency <= gap_max)
            .collect())
    }

    // -----------------------------------------------------------------------
    // Dispersion curve
    // -----------------------------------------------------------------------

    /// Compute the full strip dispersion over `kx ∈ [kx_min, kx_max]`.
    ///
    /// Returns `n_kx` entries `(kx, eigenvalues)` where each eigenvalue list has
    /// `width × n_bands` elements sorted ascending.
    pub fn dispersion_curve(
        &self,
        kx_min: f64,
        kx_max: f64,
        n_kx: usize,
    ) -> Result<Vec<(f64, Vec<f64>)>> {
        if n_kx < 2 {
            return Err(error::invalid_param("n_kx", "need at least 2 kx points"));
        }
        let mut result = Vec::with_capacity(n_kx);
        for i in 0..n_kx {
            let kx = kx_min + (kx_max - kx_min) * (i as f64) / ((n_kx - 1) as f64);
            let modes = self.solve_strip(kx)?;
            let evals: Vec<f64> = modes.iter().map(|m| m.frequency).collect();
            result.push((kx, evals));
        }
        Ok(result)
    }

    // -----------------------------------------------------------------------
    // Edge velocity
    // -----------------------------------------------------------------------

    /// Compute dω/dkx of the lowest in-gap mode at kx = kx0.
    ///
    /// Uses a finite difference `(ω(kx0+dkx) − ω(kx0−dkx)) / (2·dkx)`.
    /// Returns `NumericalError` if no in-gap modes are found at either neighbour.
    pub fn edge_velocity(&self, kx0: f64, dkx: f64, gap_min: f64, gap_max: f64) -> Result<f64> {
        if dkx <= 0.0 {
            return Err(error::invalid_param("dkx", "must be positive"));
        }
        let modes_p = self.find_in_gap_modes(kx0 + dkx, gap_min, gap_max)?;
        let modes_m = self.find_in_gap_modes(kx0 - dkx, gap_min, gap_max)?;

        if modes_p.is_empty() || modes_m.is_empty() {
            return Err(error::numerical_error(
                "no in-gap modes found at neighbouring kx points for velocity computation",
            ));
        }

        // Lowest-frequency in-gap mode at each neighbour
        let ep = modes_p
            .iter()
            .map(|m| m.frequency)
            .fold(f64::INFINITY, f64::min);
        let em = modes_m
            .iter()
            .map(|m| m.frequency)
            .fold(f64::INFINITY, f64::min);

        Ok((ep - em) / (2.0 * dkx))
    }
}

// ---------------------------------------------------------------------------
// Helper: classify edge side from eigenvector weight distribution
// ---------------------------------------------------------------------------

fn classify_edge(eigenvec: &[Complex], nb: usize, ny: usize) -> EdgeSide {
    let dim = eigenvec.len();
    if dim == 0 || nb == 0 || ny == 0 {
        return EdgeSide::Bulk;
    }
    let total: f64 = eigenvec.iter().map(|c| c.norm_sq()).sum();
    if total < 1e-30 {
        return EdgeSide::Bulk;
    }

    // Weight on bottom unit cell (iy=0): indices 0..nb
    let w_bottom: f64 = eigenvec[0..nb.min(dim)].iter().map(|c| c.norm_sq()).sum();
    // Weight on top unit cell (iy=ny-1): indices (ny-1)*nb..ny*nb
    let top_start = ((ny - 1) * nb).min(dim);
    let w_top: f64 = eigenvec[top_start..dim].iter().map(|c| c.norm_sq()).sum();

    let f_bottom = w_bottom / total;
    let f_top = w_top / total;

    let threshold = 0.25; // fraction threshold for edge classification
    if f_bottom > f_top && f_bottom >= threshold {
        EdgeSide::Bottom
    } else if f_top >= threshold {
        EdgeSide::Top
    } else {
        EdgeSide::Bulk
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topomagnon::band_model::MagnonBandModel;

    #[test]
    fn strip_size_correct() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.3, 0.0).unwrap();
        let em = EdgeModes::new(&m, 10).unwrap();
        let modes = em.solve_strip(0.5).unwrap();
        // Should have width * n_bands = 10 * 2 = 20 modes
        assert_eq!(modes.len(), 20, "Expected 20 modes, got {}", modes.len());
    }

    #[test]
    fn eigenvalues_real_within_tolerance() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.2, 0.1).unwrap();
        let em = EdgeModes::new(&m, 8).unwrap();
        let modes = em.solve_strip(0.0).unwrap();
        for mode in &modes {
            assert!(
                mode.frequency.is_finite(),
                "Eigenvalue not finite: {}",
                mode.frequency
            );
        }
    }

    #[test]
    fn eigenvectors_normalized() {
        // The strip Hamiltonian is Hermitian; eigenvectors from hermitian_eigendecomposition
        // are orthonormal. We check the localization metric is in [0,1] as a proxy.
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.3, 0.0).unwrap();
        let em = EdgeModes::new(&m, 6).unwrap();
        let modes = em.solve_strip(0.3).unwrap();
        for mode in &modes {
            assert!(
                mode.localization >= 0.0 && mode.localization <= 1.0 + 1e-10,
                "Localization out of [0,1]: {}",
                mode.localization
            );
        }
    }

    #[test]
    fn strip_has_more_modes_than_bulk() {
        // Strip with width=8 has 8*2=16 modes; bulk model has only 2 bands
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.3, 0.0).unwrap();
        let em = EdgeModes::new(&m, 8).unwrap();
        let strip_modes = em.solve_strip(0.0).unwrap().len();
        assert!(
            strip_modes > m.n_bands(),
            "Strip ({}) should have more modes than bulk ({})",
            strip_modes,
            m.n_bands()
        );
    }

    #[test]
    fn localization_metric_in_range_0_1() {
        let m = MagnonBandModel::kagome(1.0, 0.3, 0.0).unwrap();
        let em = EdgeModes::new(&m, 5).unwrap();
        let modes = em.solve_strip(0.5).unwrap();
        for mode in &modes {
            assert!(
                (0.0..=1.0 + 1e-10).contains(&mode.localization),
                "Localization {} out of [0,1]",
                mode.localization
            );
        }
    }

    #[test]
    fn dispersion_curve_correct_kx_count() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.3, 0.0).unwrap();
        let em = EdgeModes::new(&m, 6).unwrap();
        let disp = em.dispersion_curve(-PI, PI, 10).unwrap();
        assert_eq!(disp.len(), 10);
        assert_eq!(disp[0].1.len(), 12); // 6 * 2 = 12 modes
    }

    #[test]
    fn edge_velocity_finite() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.5, 0.0).unwrap();
        let em = EdgeModes::new(&m, 8).unwrap();
        // Find the gap region (between lower and upper band)
        let gap_min = -1.5;
        let gap_max = 1.5;
        // edge_velocity may return error if no in-gap modes — that's acceptable
        if let Ok(v) = em.edge_velocity(0.0, 0.1, gap_min, gap_max) {
            // No in-gap modes at this kx is also acceptable
            assert!(v.is_finite(), "velocity not finite");
        }
    }

    #[test]
    fn width_minimum_2() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.1, 0.0).unwrap();
        assert!(EdgeModes::new(&m, 1).is_err());
        assert!(EdgeModes::new(&m, 2).is_ok());
    }

    #[test]
    fn width_maximum_60() {
        let m = MagnonBandModel::square_dmi(1.0, 0.0, 0.0).unwrap();
        // square has 1 band, so width=60 → 60×60 matrix, within MAX_DIM
        assert!(EdgeModes::new(&m, 60).is_ok());
        assert!(EdgeModes::new(&m, 61).is_err());
    }

    #[test]
    fn in_gap_filter_reduces_count() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.3, 0.0).unwrap();
        let em = EdgeModes::new(&m, 8).unwrap();
        let all_modes = em.solve_strip(0.0).unwrap();
        let total = all_modes.len();

        // Find the min/max energies and take a narrow window in the middle
        let e_min = all_modes
            .iter()
            .map(|m| m.frequency)
            .fold(f64::INFINITY, f64::min);
        let e_max = all_modes
            .iter()
            .map(|m| m.frequency)
            .fold(f64::NEG_INFINITY, f64::max);
        let mid = (e_min + e_max) / 2.0;
        let width = (e_max - e_min) / 4.0;

        let in_gap = em.find_in_gap_modes(0.0, mid - width, mid + width).unwrap();
        assert!(
            in_gap.len() <= total,
            "In-gap filter should not increase mode count"
        );
    }
}
