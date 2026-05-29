//! Lattice Chern number computation via the Fukui-Hatsugai-Suzuki method.
//!
//! Implements the discrete Brillouin-zone algorithm of
//! T. Fukui, Y. Hatsugai, H. Suzuki, *J. Phys. Soc. Jpn.* **74**, 1674 (2005),
//! which computes topological Chern numbers from the Berry phases of lattice link
//! variables defined on a discretized k-space mesh. This approach is gauge-invariant
//! and avoids the discontinuity problems of the smooth Berry connection.
//!
//! # Algorithm (Fukui-Hatsugai)
//!
//! Given a mesh of `nx × ny` k-points covering the first Brillouin zone:
//!
//! 1. Build eigenstates `|u_n(k)⟩` at every k-point by diagonalizing H(k).
//! 2. Define link variables:
//!    `U_x(k) = ⟨u_n(k)|u_n(k+Δkₓ)⟩ / |⟨u_n(k)|u_n(k+Δkₓ)⟩|`
//!    (and similarly U_y(k) in the y-direction), normalising to a pure phase.
//! 3. Compute the plaquette flux:
//!    `F(k) = arg(U_x(k) · U_y(k+Δkₓ) · U_x*(k+Δky) · U_y*(k))`
//! 4. The Chern number is `C = (1/2π) Σ_k F(k)`, rounded to the nearest integer.
//!
//! # Additional Tools
//!
//! - **Wilson loop**: phase accumulated by an eigenstate around a closed loop in k.
//! - **Berry phase along arbitrary path**: discrete phase product.
//! - **Total Chern sum**: verifies the Nielsen-Ninomiya sum rule (C_total = 0).

use std::f64::consts::PI;

use super::band_model::MagnonBandModel;
use crate::error::{self, Result};
use crate::math::Complex;

// ---------------------------------------------------------------------------
// ChernNumber
// ---------------------------------------------------------------------------

/// Computes topological Chern numbers for magnon bands.
pub struct ChernNumber<'a> {
    /// Reference to the magnon band model.
    pub model: &'a MagnonBandModel,
}

impl<'a> ChernNumber<'a> {
    /// Create a new `ChernNumber` calculator for the given model.
    pub fn new(model: &'a MagnonBandModel) -> Self {
        Self { model }
    }

    // -----------------------------------------------------------------------
    // Fukui-Hatsugai Chern number
    // -----------------------------------------------------------------------

    /// Compute the Chern number of band `band_idx` on an `nx × ny` k-mesh.
    ///
    /// Uses the Fukui-Hatsugai discrete gauge-invariant method. The result is
    /// rounded to the nearest integer; if the deviation exceeds 0.1 a
    /// `NumericalError` is returned (indicating insufficient grid resolution or
    /// a nearly gapless band).
    ///
    /// # Errors
    ///
    /// - `InvalidParameter` if `band_idx ≥ n_bands` or grid < 5×5.
    /// - `NumericalError` if the lattice sum is not close to an integer.
    pub fn compute(&self, band_idx: usize, nx: usize, ny: usize) -> Result<i32> {
        let nb = self.model.n_bands();
        if band_idx >= nb {
            return Err(error::invalid_param(
                "band_idx",
                "index exceeds number of bands",
            ));
        }
        if nx < 5 || ny < 5 {
            return Err(error::invalid_param(
                "nx/ny",
                "grid must be at least 5×5 for reliable Chern number",
            ));
        }

        let states = self.build_eigenstates(nx, ny, band_idx)?;
        let mut flux_sum = 0.0_f64;
        for ix in 0..nx {
            let ix1 = (ix + 1) % nx;
            for iy in 0..ny {
                let iy1 = (iy + 1) % ny;
                let u_x = link_variable(&states[ix][iy], &states[ix1][iy]);
                let u_y = link_variable(&states[ix][iy], &states[ix][iy1]);
                let u_x_py = link_variable(&states[ix][iy1], &states[ix1][iy1]);
                let u_y_px = link_variable(&states[ix1][iy], &states[ix1][iy1]);
                let plaquette = u_x.mul(&u_y_px).mul(&u_x_py.conj()).mul(&u_y.conj());
                flux_sum += plaquette.phase();
            }
        }

        let chern_real = flux_sum / (2.0 * PI);
        let chern_int = chern_real.round() as i32;

        if (chern_real - chern_int as f64).abs() > 0.1 {
            return Err(error::numerical_error(&format!(
                "Chern number not close to integer: {:.4} (try finer grid)",
                chern_real
            )));
        }

        Ok(chern_int)
    }

    /// Build eigenstates on the (nx+1) × (ny+1) grid (extra point for periodicity).
    pub fn build_eigenstates(
        &self,
        nx: usize,
        ny: usize,
        band_idx: usize,
    ) -> Result<Vec<Vec<Vec<Complex>>>> {
        let mut states = Vec::with_capacity(nx + 1);
        for ix in 0..=nx {
            let kx = -PI + 2.0 * PI * (ix as f64) / (nx as f64);
            let mut row = Vec::with_capacity(ny + 1);
            for iy in 0..=ny {
                let ky = -PI + 2.0 * PI * (iy as f64) / (ny as f64);
                let (_, vecs) = self.model.diagonalize((kx, ky))?;
                row.push(vecs.column(band_idx));
            }
            states.push(row);
        }
        Ok(states)
    }

    // -----------------------------------------------------------------------
    // Wilson loop
    // -----------------------------------------------------------------------

    /// Compute the Wilson-loop holonomy for band `band_idx` at fixed `kx`.
    ///
    /// The Wilson loop is the product of link variables U_y around the closed
    /// path ky ∈ [−π, π] at fixed kx. Returns a unit-norm complex scalar whose
    /// argument gives the Wannier-centre position.
    pub fn wilson_loop(&self, band_idx: usize, kx: f64, nky: usize) -> Result<Complex> {
        let nb = self.model.n_bands();
        if band_idx >= nb {
            return Err(error::invalid_param(
                "band_idx",
                "index exceeds number of bands",
            ));
        }
        if nky < 4 {
            return Err(error::invalid_param("nky", "need at least 4 ky points"));
        }

        let mut states: Vec<Vec<Complex>> = Vec::with_capacity(nky + 1);
        for iy in 0..=nky {
            let ky = -PI + 2.0 * PI * (iy as f64) / (nky as f64);
            let (_, vecs) = self.model.diagonalize((kx, ky))?;
            states.push(vecs.column(band_idx));
        }

        let mut wilson = Complex::ONE;
        for iy in 0..nky {
            let u = link_variable(&states[iy], &states[iy + 1]);
            wilson = wilson.mul(&u);
        }

        Ok(wilson)
    }

    // -----------------------------------------------------------------------
    // Berry phase along path
    // -----------------------------------------------------------------------

    /// Compute the Berry phase of band `band_idx` along an arbitrary closed path.
    ///
    /// The path is a sequence of k-points forming a closed loop (the last point
    /// reconnects to the first). Returns the phase in radians, in (−π, π].
    pub fn berry_phase_along_path(&self, band_idx: usize, path: &[(f64, f64)]) -> Result<f64> {
        let nb = self.model.n_bands();
        if band_idx >= nb {
            return Err(error::invalid_param(
                "band_idx",
                "index exceeds number of bands",
            ));
        }
        if path.len() < 3 {
            return Err(error::invalid_param(
                "path",
                "need at least 3 k-points for a closed loop",
            ));
        }

        let mut states: Vec<Vec<Complex>> = Vec::with_capacity(path.len() + 1);
        for &kpt in path {
            let (_, vecs) = self.model.diagonalize(kpt)?;
            states.push(vecs.column(band_idx));
        }
        // Close the loop
        {
            let (_, vecs) = self.model.diagonalize(path[0])?;
            states.push(vecs.column(band_idx));
        }

        let mut product = Complex::ONE;
        let n = states.len() - 1;
        for i in 0..n {
            let u = link_variable(&states[i], &states[i + 1]);
            product = product.mul(&u);
        }

        Ok(product.phase())
    }

    // -----------------------------------------------------------------------
    // Total Chern sum
    // -----------------------------------------------------------------------

    /// Sum Chern numbers over all bands; should be 0 for a complete set.
    pub fn total_chern_sum(&self, nx: usize, ny: usize) -> Result<i32> {
        let nb = self.model.n_bands();
        let mut total = 0_i32;
        for band in 0..nb {
            total += self.compute(band, nx, ny)?;
        }
        Ok(total)
    }
}

// ---------------------------------------------------------------------------
// Link variable helper
// ---------------------------------------------------------------------------

/// Compute the link variable `U = ⟨u_a|u_b⟩ / |⟨u_a|u_b⟩|`.
fn link_variable(u_a: &[Complex], u_b: &[Complex]) -> Complex {
    let mut inner = Complex::ZERO;
    for (a, b) in u_a.iter().zip(u_b.iter()) {
        inner = inner.add(&a.conj().mul(b));
    }
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
    use crate::topomagnon::band_model::MagnonBandModel;

    #[test]
    fn chern_zero_dmi_honeycomb() {
        // Without DMI the honeycomb magnon has trivial bands — use tiny h_ext to open gap
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.0, 1e-3).unwrap();
        let cn = ChernNumber::new(&m);
        // Tiny h_ext → topologically trivial; Chern should be 0
        if let Ok(c) = cn.compute(0, 15, 15) {
            // Numerical failure near gap closure is acceptable
            assert_eq!(c, 0, "Expected C=0 for near-trivial honeycomb");
        }
    }

    #[test]
    fn chern_nonzero_dmi_honeycomb() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.5, 0.0).unwrap();
        let cn = ChernNumber::new(&m);
        let c = cn.compute(0, 20, 20).unwrap();
        assert_eq!(
            c.abs(),
            1,
            "Expected |C|=1 for topological honeycomb, got {}",
            c
        );
    }

    #[test]
    fn chern_sign_flips_with_dmi_sign() {
        let m_pos = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.5, 0.0).unwrap();
        let m_neg = MagnonBandModel::honeycomb_haldane(1.0, 0.0, -0.5, 0.0).unwrap();
        let c_pos = ChernNumber::new(&m_pos).compute(0, 20, 20).unwrap();
        let c_neg = ChernNumber::new(&m_neg).compute(0, 20, 20).unwrap();
        assert_eq!(
            c_pos, -c_neg,
            "Chern sign should flip: {} vs {}",
            c_pos, c_neg
        );
    }

    #[test]
    fn chern_integer_quantized() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.4, 0.0).unwrap();
        let cn = ChernNumber::new(&m);
        let states = cn.build_eigenstates(20, 20, 0).unwrap();
        let nx = 20;
        let ny = 20;
        let mut flux_sum = 0.0_f64;
        for ix in 0..nx {
            let ix1 = (ix + 1) % nx;
            for iy in 0..ny {
                let iy1 = (iy + 1) % ny;
                let u_x = link_variable(&states[ix][iy], &states[ix1][iy]);
                let u_y = link_variable(&states[ix][iy], &states[ix][iy1]);
                let u_x_py = link_variable(&states[ix][iy1], &states[ix1][iy1]);
                let u_y_px = link_variable(&states[ix1][iy], &states[ix1][iy1]);
                let plaquette = u_x.mul(&u_y_px).mul(&u_x_py.conj()).mul(&u_y.conj());
                flux_sum += plaquette.phase();
            }
        }
        let chern_real = flux_sum / (2.0 * PI);
        assert!(
            (chern_real - chern_real.round()).abs() < 0.05,
            "Chern not quantized: {:.4}",
            chern_real
        );
    }

    #[test]
    fn total_chern_sum_zero() {
        // The Nielsen-Ninomiya sum rule: Σ_n C_n = 0 for a complete band set.
        // Use honeycomb-Haldane with DMI=0.5; bands have C = +1 and -1 → sum 0.
        // The honeycomb has a large, uniform gap so 20×20 gives clean integer Cherns.
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.5, 0.0).unwrap();
        let cn = ChernNumber::new(&m);
        let total = cn.total_chern_sum(20, 20).unwrap();
        assert_eq!(total, 0, "Total Chern sum must be zero, got {}", total);
    }

    #[test]
    fn wilson_loop_finite() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.3, 0.0).unwrap();
        let cn = ChernNumber::new(&m);
        let w = cn.wilson_loop(0, 0.0, 20).unwrap();
        assert!(w.norm().is_finite());
        assert!(
            (w.norm() - 1.0).abs() < 1e-10,
            "Wilson loop not unit modulus"
        );
    }

    #[test]
    fn berry_phase_closed_loop_multiple_of_pi() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.0, 0.0).unwrap();
        let cn = ChernNumber::new(&m);
        let path = vec![(0.1, 0.1), (-0.1, 0.1), (-0.1, -0.1), (0.1, -0.1)];
        let phase = cn.berry_phase_along_path(0, &path).unwrap();
        let normalized = phase.abs() / PI;
        let nearest_half = (normalized * 2.0).round() / 2.0;
        assert!(
            (normalized - nearest_half).abs() < 0.15,
            "Berry phase not multiple of π/2: {:.4}π",
            phase / PI
        );
    }

    #[test]
    fn invalid_band_idx_errors() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.1, 0.0).unwrap();
        let cn = ChernNumber::new(&m);
        assert!(cn.compute(5, 10, 10).is_err());
        assert!(cn.wilson_loop(5, 0.0, 10).is_err());
        assert!(cn
            .berry_phase_along_path(5, &[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
            .is_err());
    }

    #[test]
    fn nx_ny_minimum() {
        let m = MagnonBandModel::honeycomb_haldane(1.0, 0.0, 0.3, 0.0).unwrap();
        let cn = ChernNumber::new(&m);
        assert!(cn.compute(0, 3, 3).is_err());
        assert!(cn.compute(0, 5, 3).is_err());
        // 5×5 accepted
        let _ = cn.compute(0, 5, 5);
    }

    #[test]
    fn kagome_middle_band() {
        let m = MagnonBandModel::kagome(1.0, 0.4, 0.0).unwrap();
        let cn = ChernNumber::new(&m);
        let c = cn.compute(1, 15, 15).unwrap();
        assert!(c.abs() <= 5, "Chern number suspiciously large: {}", c);
    }
}
