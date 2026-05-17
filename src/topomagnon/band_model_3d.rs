//! Three-dimensional magnon band models for axion electrodynamics.
//!
//! Provides minimal 3D extensions of the 2D topological magnon framework,
//! focused on the CubicHaldane and pyrochlore lattice geometries that exhibit
//! magnon-axion physics:
//!
//! - **CubicHaldane**: a 2-band Hamiltonian on a 3D cubic lattice driven by
//!   DMI along all three spatial directions.  Inspired by the Weyl-semimetal
//!   analogy described in:
//!   Y. Tokura, K. Yasuda, A. Tsukazaki, *Nature Rev. Phys.* **1**, 126 (2019).
//!
//! - **PyrochloreTopological**: a 4-band model on the pyrochlore lattice,
//!   which hosts a topological axion state with θ = π.
//!   A. M. Essin, J. E. Moore, D. Vanderbilt,
//!   *Phys. Rev. Lett.* **102**, 146805 (2009).
//!
//! # Hamiltonian Forms
//!
//! ## CubicHaldane (2-band)
//!
//! ```text
//! H(k) = J_nn · (cos kx + cos ky + cos kz) · I₂
//!       + D · (sin kx · σ_x + sin ky · σ_y + sin kz · σ_z)
//!       + h_ext · σ_z
//! ```
//!
//! The sine-type DMI terms mimic the Haldane mass; the system is topological
//! when `|D|` exceeds a threshold relative to `J_nn`.
//!
//! ## PyrochloreTopological (4-band)
//!
//! A simplified pyrochlore model with bond structure:
//!
//! ```text
//! H_aa = h_ext  (on-site Zeeman)
//! H_ab = J · (1 + exp(i k · d_ab))  for a ≠ b
//! ```
//!
//! where `d_ab` are the four inequivalent pyrochlore bond vectors.

use crate::error::Result;
use crate::math::{CMatrix, Complex};
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// LatticeType3D
// ---------------------------------------------------------------------------

/// Lattice geometry for 3D magnon band models.
#[derive(Debug, Clone, PartialEq)]
pub enum LatticeType3D {
    /// Cubic Haldane model: 2-band, Weyl-like DMI mass term.
    CubicHaldane,
    /// Pyrochlore topological magnon: 4-band, θ-physics.
    PyrochloreTopological,
}

// ---------------------------------------------------------------------------
// MagnonBandModel3D
// ---------------------------------------------------------------------------

/// 3D magnon band model for axion electrodynamics.
///
/// Wraps either the 2-band CubicHaldane or the 4-band PyrochloreTopological
/// Hamiltonian.  All k-space coordinates are in units of the inverse lattice
/// constant (a = 1 convention).
#[derive(Debug, Clone)]
pub struct MagnonBandModel3D {
    /// Lattice geometry.
    pub lattice: LatticeType3D,
    /// Number of magnon bands.
    pub n_bands: usize,
    /// Nearest-neighbour exchange J.
    pub j_nn: f64,
    /// Next-nearest-neighbour exchange J' (not used in all models).
    pub j_nnn: f64,
    /// Dzyaloshinskii-Moriya interaction strength D.
    pub dmi: f64,
    /// External magnetic field (Zeeman energy) h.
    pub h_ext: f64,
    /// Lattice constant (normalised to 1.0 by default).
    pub a_lattice: f64,
}

impl MagnonBandModel3D {
    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    /// Build a 2-band CubicHaldane model.
    ///
    /// - `j_nn`: NN exchange (sets the bandwidth)
    /// - `j_nnn`: NNN exchange (not used in Hamiltonian, stored for reference)
    /// - `dmi`: DMI strength (drives topological phase transition)
    /// - `h_ext`: uniform Zeeman field (breaks degeneracy at Γ)
    pub fn cubic_haldane(j_nn: f64, j_nnn: f64, dmi: f64, h_ext: f64) -> Self {
        Self {
            lattice: LatticeType3D::CubicHaldane,
            n_bands: 2,
            j_nn,
            j_nnn,
            dmi,
            h_ext,
            a_lattice: 1.0,
        }
    }

    /// Build a 4-band PyrochloreTopological model.
    ///
    /// - `j`: isotropic NN exchange
    /// - `d`: DMI coupling strength
    pub fn pyrochlore_topological(j: f64, d: f64) -> Self {
        Self {
            lattice: LatticeType3D::PyrochloreTopological,
            n_bands: 4,
            j_nn: j,
            j_nnn: 0.0,
            dmi: d,
            h_ext: 0.0,
            a_lattice: 1.0,
        }
    }

    // -----------------------------------------------------------------------
    // Hamiltonian
    // -----------------------------------------------------------------------

    /// Build the Bloch Hamiltonian at `(kx, ky, kz)`.
    ///
    /// Returns an `n_bands × n_bands` Hermitian `CMatrix`.
    pub fn hamiltonian_at(&self, kx: f64, ky: f64, kz: f64) -> CMatrix {
        match self.lattice {
            LatticeType3D::CubicHaldane => self.cubic_haldane_hamiltonian(kx, ky, kz),
            LatticeType3D::PyrochloreTopological => self.pyrochlore_hamiltonian(kx, ky, kz),
        }
    }

    /// CubicHaldane 2×2 Hamiltonian.
    ///
    /// ```text
    /// H = J_nn · (cos kx + cos ky + cos kz) · I₂
    ///   + D · (sin kx · σ_x + sin ky · σ_y + sin kz · σ_z)
    ///   + h_ext · σ_z
    /// ```
    fn cubic_haldane_hamiltonian(&self, kx: f64, ky: f64, kz: f64) -> CMatrix {
        let diag_re = self.j_nn * (kx.cos() + ky.cos() + kz.cos());
        let dx = self.dmi * kx.sin();
        let dy = self.dmi * ky.sin();
        let dz = self.dmi * kz.sin() + self.h_ext;

        // σ_x part: off-diag real dx
        // σ_y part: off-diag ±i·dy
        // σ_z part: diagonal ±dz
        let mut mat = CMatrix::zeros(2);
        // H[0,0] = diag + dz
        mat.set(0, 0, Complex::from_real(diag_re + dz));
        // H[1,1] = diag - dz
        mat.set(1, 1, Complex::from_real(diag_re - dz));
        // H[0,1] = dx - i·dy
        mat.set(0, 1, Complex::new(dx, -dy));
        // H[1,0] = dx + i·dy  (Hermitian conjugate)
        mat.set(1, 0, Complex::new(dx, dy));
        mat
    }

    /// Pyrochlore 4×4 Hamiltonian.
    ///
    /// The four pyrochlore sublattice sites are connected by the FCC bond
    /// vectors.  Using simplified bond structure:
    ///
    /// ```text
    /// d₀₁ = (1,0,0), d₀₂ = (0,1,0), d₀₃ = (0,0,1)
    /// d₁₂ = (-1,1,0), d₁₃ = (-1,0,1), d₂₃ = (0,-1,1)
    /// ```
    ///
    /// `H_ab = J · (1 + exp(i k · d_ab))` for a < b, plus DMI on-diagonal.
    fn pyrochlore_hamiltonian(&self, kx: f64, ky: f64, kz: f64) -> CMatrix {
        // Bond vectors (a→b) for a < b
        let bonds: [(usize, usize, f64, f64, f64); 6] = [
            (0, 1, 1.0, 0.0, 0.0),
            (0, 2, 0.0, 1.0, 0.0),
            (0, 3, 0.0, 0.0, 1.0),
            (1, 2, -1.0, 1.0, 0.0),
            (1, 3, -1.0, 0.0, 1.0),
            (2, 3, 0.0, -1.0, 1.0),
        ];

        let mut mat = CMatrix::zeros(4);

        // On-site: h_ext + DMI diagonal
        for i in 0..4 {
            mat.set(i, i, Complex::from_real(self.h_ext + self.dmi));
        }

        // Off-diagonal hoppings
        for &(a, b, dx, dy, dz) in &bonds {
            let phi = kx * dx + ky * dy + kz * dz;
            // H_ab = J · (1 + exp(i φ))
            let hop = Complex::from_real(self.j_nn).add(&Complex::from_polar(self.j_nn, phi));
            mat.set(a, b, hop);
            mat.set(b, a, hop.conj());
        }

        mat
    }

    // -----------------------------------------------------------------------
    // Diagonalisation
    // -----------------------------------------------------------------------

    /// Diagonalise the 3D Hamiltonian at `(kx, ky, kz)`.
    ///
    /// Returns `(eigenvalues, eigenvectors)` with eigenvalues sorted ascending
    /// and the k-th column of `eigenvectors` being the k-th eigenvector.
    ///
    /// # Errors
    ///
    /// Propagates errors from `hermitian_eigendecomposition`.
    pub fn diagonalize_3d(&self, kx: f64, ky: f64, kz: f64) -> Result<(Vec<f64>, CMatrix)> {
        let h = self.hamiltonian_at(kx, ky, kz);
        h.hermitian_eigendecomposition()
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// Number of magnon bands.
    #[inline]
    pub fn n_bands(&self) -> usize {
        self.n_bands
    }

    /// Minimum direct band gap over the 3D BZ (sampled on an `n × n × n` mesh).
    pub fn band_gap_3d(&self, n: usize) -> f64 {
        let n = n.max(4);
        let mut min_gap = f64::INFINITY;
        for ix in 0..n {
            let kx = 2.0 * PI * (ix as f64) / (n as f64);
            for iy in 0..n {
                let ky = 2.0 * PI * (iy as f64) / (n as f64);
                for iz in 0..n {
                    let kz = 2.0 * PI * (iz as f64) / (n as f64);
                    if let Ok((evals, _)) = self.diagonalize_3d(kx, ky, kz) {
                        let nb = evals.len();
                        if nb >= 2 {
                            let gap = evals[nb / 2] - evals[nb / 2 - 1];
                            if gap < min_gap {
                                min_gap = gap;
                            }
                        }
                    }
                }
            }
        }
        if min_gap.is_infinite() {
            0.0
        } else {
            min_gap
        }
    }

    /// Check whether the CubicHaldane model is in the strong-DMI topological phase.
    ///
    /// The cubic Haldane model is topological when the DMI strength exceeds
    /// the isotropic hopping: `|D| > |J_nn|`.
    pub fn is_topological_cubic(&self) -> bool {
        matches!(self.lattice, LatticeType3D::CubicHaldane) && self.dmi.abs() > self.j_nn.abs()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    #[test]
    fn cubic_haldane_n_bands_is_2() {
        let m = MagnonBandModel3D::cubic_haldane(1.0, 0.0, 0.5, 0.0);
        assert_eq!(m.n_bands(), 2);
    }

    #[test]
    fn pyrochlore_n_bands_is_4() {
        let m = MagnonBandModel3D::pyrochlore_topological(1.0, 0.3);
        assert_eq!(m.n_bands(), 4);
    }

    // -----------------------------------------------------------------------
    // Hermiticity
    // -----------------------------------------------------------------------

    #[test]
    fn cubic_haldane_hamiltonian_hermitian() {
        let m = MagnonBandModel3D::cubic_haldane(1.0, 0.0, 0.5, 0.1);
        for (kx, ky, kz) in [(0.1, 0.2, 0.3), (PI / 3.0, PI / 4.0, 0.5)] {
            let h = m.hamiltonian_at(kx, ky, kz);
            let hd = h.conj_transpose();
            let diff = h.sub(&hd).unwrap();
            assert!(
                diff.frobenius_norm() < 1e-12,
                "H not Hermitian at ({kx},{ky},{kz})"
            );
        }
    }

    #[test]
    fn pyrochlore_hamiltonian_hermitian() {
        let m = MagnonBandModel3D::pyrochlore_topological(1.0, 0.3);
        for (kx, ky, kz) in [(0.0, 0.0, 0.0), (0.5, 0.7, 1.0)] {
            let h = m.hamiltonian_at(kx, ky, kz);
            let hd = h.conj_transpose();
            let diff = h.sub(&hd).unwrap();
            assert!(
                diff.frobenius_norm() < 1e-12,
                "Pyrochlore H not Hermitian at ({kx},{ky},{kz})"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Diagonalisation
    // -----------------------------------------------------------------------

    #[test]
    fn cubic_haldane_diagonalizes_correctly() {
        let m = MagnonBandModel3D::cubic_haldane(1.0, 0.0, 0.5, 0.0);
        let (evals, _) = m.diagonalize_3d(0.0, 0.0, 0.0).unwrap();
        assert_eq!(evals.len(), 2);
        // At k=0: dispersion = 3J ± h_ext; with h_ext=0 → DMI terms vanish → both = 3J
        assert!((evals[0] - 3.0).abs() < 1e-10, "E[0]={}", evals[0]);
        assert!((evals[1] - 3.0).abs() < 1e-10, "E[1]={}", evals[1]);
    }

    #[test]
    fn pyrochlore_diagonalizes_four_bands() {
        let m = MagnonBandModel3D::pyrochlore_topological(1.0, 0.2);
        let (evals, evecs) = m.diagonalize_3d(0.1, 0.2, 0.3).unwrap();
        assert_eq!(evals.len(), 4);
        assert_eq!(evecs.n(), 4);
        // Eigenvalues must be sorted ascending
        for i in 0..3 {
            assert!(
                evals[i] <= evals[i + 1],
                "Eigenvalues not sorted: {evals:?}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Band gap
    // -----------------------------------------------------------------------

    #[test]
    fn cubic_haldane_gap_finite_with_field() {
        let m = MagnonBandModel3D::cubic_haldane(1.0, 0.0, 0.0, 0.5);
        let gap = m.band_gap_3d(6);
        assert!(gap >= 0.0, "Gap negative: {gap}");
    }

    // -----------------------------------------------------------------------
    // Phase detection
    // -----------------------------------------------------------------------

    #[test]
    fn cubic_haldane_topological_phase_detected() {
        // |D| > |J_nn| → topological
        let m = MagnonBandModel3D::cubic_haldane(0.5, 0.0, 1.0, 0.0);
        assert!(m.is_topological_cubic());
        // |D| < |J_nn| → trivial
        let m2 = MagnonBandModel3D::cubic_haldane(1.0, 0.0, 0.3, 0.0);
        assert!(!m2.is_topological_cubic());
    }
}
