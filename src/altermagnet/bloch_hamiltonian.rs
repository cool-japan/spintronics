//! Generic Bloch-Hamiltonian abstraction for the altermagnet electronic band model.
//!
//! [`BlochHamiltonian`] is a small, self-contained trait exposing a momentum-space
//! Hamiltonian matrix `H(kx, ky)` plus a Hermitian-diagonalization convenience
//! method. It is new, non-invasive code: `src/topomagnon/*` (`BerryCurvature`,
//! `MagnonBandModel`, and the Chern-number/Wilson-loop/axion/edge-mode code built
//! on top of them) already has its own internal Hamiltonian-construction and
//! diagonalization pattern and is deliberately left untouched and independent of
//! this trait, so this module's blast radius stays confined to
//! `src/altermagnet/*`.
//!
//! [`AltermagnetSpinHamiltonian`] implements this trait for a single spin sector
//! of [`AltermagnetBandModel`], and
//! [`KuboBerry`](crate::altermagnet::kubo_berry::KuboBerry) is a finite-difference
//! Kubo-formula Berry-curvature calculator built purely against this trait, used
//! as a numerical cross-check of
//! [`AltermagnetBandModel::berry_curvature`]'s closed-form result.

use crate::altermagnet::band_model::{AltermagnetBandModel, Spin};
use crate::error::Result;
use crate::math::CMatrix;

/// A momentum-space (Bloch) Hamiltonian usable by generic diagonalization and
/// Berry-curvature utilities.
///
/// Implementors expose the Hamiltonian matrix `H(kx, ky)` at an arbitrary
/// momentum point and the number of bands (the matrix dimension).
/// [`diagonalize_at`](Self::diagonalize_at) is a provided convenience built
/// directly on top of [`CMatrix::hermitian_eigendecomposition`] -- it is always
/// correct for any Hermitian `H(k)` and does not need to be overridden, though
/// implementors with a faster analytic diagonalization (e.g. a closed-form 2x2
/// solution) are free to do so.
pub trait BlochHamiltonian {
    /// The Hamiltonian matrix `H(kx, ky)` at momentum `(kx, ky)`.
    ///
    /// Implementations must return a Hermitian matrix of dimension
    /// [`n_bands`](Self::n_bands).
    ///
    /// # Errors
    ///
    /// Returns an error if the matrix cannot be constructed at this momentum.
    fn hamiltonian_at(&self, kx: f64, ky: f64) -> Result<CMatrix>;

    /// Matrix dimension (number of bands) of [`hamiltonian_at`](Self::hamiltonian_at).
    fn n_bands(&self) -> usize;

    /// Diagonalize [`hamiltonian_at`](Self::hamiltonian_at) at `(kx, ky)`.
    ///
    /// Returns `(eigenvalues, eigenvectors)` sorted ascending, following the
    /// convention of [`CMatrix::hermitian_eigendecomposition`]: the k-th column
    /// of `eigenvectors` is the eigenvector for the k-th eigenvalue.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`hamiltonian_at`](Self::hamiltonian_at) and from
    /// [`CMatrix::hermitian_eigendecomposition`] (e.g. non-convergence).
    fn diagonalize_at(&self, kx: f64, ky: f64) -> Result<(Vec<f64>, CMatrix)> {
        self.hamiltonian_at(kx, ky)?.hermitian_eigendecomposition()
    }
}

/// A single per-spin `2x2` sector of an [`AltermagnetBandModel`], viewed as a
/// [`BlochHamiltonian`].
///
/// Wraps a `&AltermagnetBandModel` plus a fixed [`Spin`] so the generic
/// [`BlochHamiltonian`] machinery (in particular
/// [`KuboBerry`](crate::altermagnet::kubo_berry::KuboBerry)) can be used without
/// threading a spin argument through the trait itself: `H_sigma(k)` is a
/// genuinely different `2x2` matrix for each spin sector (see the
/// [band model docs](crate::altermagnet::band_model) for the full Hamiltonian).
///
/// # Example
///
/// ```rust
/// use spintronics::altermagnet::{AltermagnetBandModel, AltermagnetSpinHamiltonian, Spin};
/// use spintronics::altermagnet::BlochHamiltonian;
///
/// let model = AltermagnetBandModel::mnte();
/// let up = AltermagnetSpinHamiltonian::new(&model, Spin::Up);
/// assert_eq!(up.n_bands(), 2);
/// let h = up.hamiltonian_at(0.4, 0.15).expect("valid momentum");
/// assert_eq!(h.n(), 2);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct AltermagnetSpinHamiltonian<'a> {
    /// The underlying band model.
    pub model: &'a AltermagnetBandModel,
    /// Which spin sector this Hamiltonian represents.
    pub spin: Spin,
}

impl<'a> AltermagnetSpinHamiltonian<'a> {
    /// Create a per-spin Bloch Hamiltonian view of `model`.
    pub fn new(model: &'a AltermagnetBandModel, spin: Spin) -> Self {
        Self { model, spin }
    }
}

impl BlochHamiltonian for AltermagnetSpinHamiltonian<'_> {
    fn hamiltonian_at(&self, kx: f64, ky: f64) -> Result<CMatrix> {
        self.model.hamiltonian_matrix(kx, ky, self.spin)
    }

    fn n_bands(&self) -> usize {
        2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::altermagnet::materials::AltermagneticSymmetry;

    fn soc_model() -> AltermagnetBandModel {
        AltermagnetBandModel::new(
            AltermagneticSymmetry::DWave,
            1.0,
            0.3,
            0.8,
            0.6,
            0.4,
            0.5,
            0.0,
            1.0,
        )
        .expect("parameters are valid")
    }

    #[test]
    fn n_bands_is_two() {
        let model = soc_model();
        let up = AltermagnetSpinHamiltonian::new(&model, Spin::Up);
        let down = AltermagnetSpinHamiltonian::new(&model, Spin::Down);
        assert_eq!(up.n_bands(), 2);
        assert_eq!(down.n_bands(), 2);
    }

    #[test]
    fn hamiltonian_at_matches_model_hamiltonian_matrix() {
        let model = soc_model();
        let up = AltermagnetSpinHamiltonian::new(&model, Spin::Up);
        let (kx, ky) = (0.4, 0.15);
        let via_trait = up.hamiltonian_at(kx, ky).expect("valid");
        let via_model = model
            .hamiltonian_matrix(kx, ky, Spin::Up)
            .expect("valid");
        for i in 0..2 {
            for j in 0..2 {
                let a = via_trait.get(i, j);
                let b = via_model.get(i, j);
                assert!((a.re - b.re).abs() < 1e-15 && (a.im - b.im).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn diagonalize_at_matches_band_energy() {
        let model = soc_model();
        for spin in [Spin::Up, Spin::Down] {
            let h = AltermagnetSpinHamiltonian::new(&model, spin);
            let (kx, ky) = (0.4, 0.15);
            let (evals, _) = h.diagonalize_at(kx, ky).expect("valid diagonalization");
            assert_eq!(evals.len(), 2);
            let expect_lower =
                model.band_energy(kx, ky, spin, crate::altermagnet::band_model::Band::Lower);
            let expect_upper =
                model.band_energy(kx, ky, spin, crate::altermagnet::band_model::Band::Upper);
            assert!((evals[0] - expect_lower).abs() < 1e-9);
            assert!((evals[1] - expect_upper).abs() < 1e-9);
        }
    }

    // `col` indexes `evals` (via `evals[col]`) and `evecs`/`hv` (via `.get(row,
    // col)`) simultaneously -- a genuinely index-driven loop, not a false
    // positive for `needless_range_loop` (same pattern already suppressed in
    // `src/math/matrix.rs`'s eigendecomposition tests).
    #[allow(clippy::needless_range_loop)]
    #[test]
    fn diagonalize_at_default_impl_is_hermitian_consistent() {
        // The default `diagonalize_at` just calls `hermitian_eigendecomposition`
        // on `hamiltonian_at`; verify H v = lambda v holds for both eigenpairs.
        let model = soc_model();
        let h = AltermagnetSpinHamiltonian::new(&model, Spin::Down);
        let (kx, ky) = (-0.6, 0.9);
        let hmat = h.hamiltonian_at(kx, ky).expect("valid");
        let (evals, evecs) = h.diagonalize_at(kx, ky).expect("valid");
        let hv = hmat.matmul(&evecs).expect("2x2 matmul");
        for col in 0..2 {
            for row in 0..2 {
                let lhs = hv.get(row, col);
                let rhs = evecs.get(row, col).scale(evals[col]);
                assert!((lhs.re - rhs.re).abs() < 1e-9 && (lhs.im - rhs.im).abs() < 1e-9);
            }
        }
    }
}
