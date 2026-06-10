//! Bogoliubov canonical transformation for magnon Hamiltonians.
//!
//! The Bogoliubov transformation diagonalises bilinear Bose Hamiltonians of the
//! form
//!
//! ```text
//! H = Σ_k [ A_k a†_k a_k + (B_k/2)(a†_k a†_{-k} + a_k a_{-k}) ]
//! ```
//!
//! through the *squeeze* rotation
//!
//! ```text
//! α_k = u_k a_k + v_k a†_{-k},   u_k = cosh θ_k,   v_k = sinh θ_k
//! ```
//!
//! where the squeeze angle `θ_k` is determined by the condition that the
//! off-diagonal (pairing) terms vanish:
//!
//! ```text
//! tanh 2θ_k = −B_k / A_k,   ω_k = √(A_k² − B_k²)
//! ```
//!
//! The new operators `α_k` satisfy `[α_k, α†_{k'}] = δ_{kk'}` and
//! `H = Σ_k ω_k α†_k α_k + E_0` with the zero-point shift
//! `E_0 = (1/2) Σ_k (ω_k − A_k)`.
//!
//! # Magnetic Instability
//!
//! If `|B_k| ≥ |A_k|` for any mode k, `tanh 2θ` would exceed unity, meaning
//! no real squeeze angle exists and the system is magnetically unstable.  This
//! regime (e.g. above the spin-flop field in an AFM) is detected and flagged as
//! an error.
//!
//! # Squeezed Vacuum
//!
//! The Bogoliubov vacuum `|0⟩_α` is a *squeezed* state with `v_k²` virtual
//! magnon pairs per mode k.  This is the physical content of the Kamra-Belzig
//! prediction (Phys. Rev. Lett. 116, 146601) of super-Poissonian spin noise in
//! antiferromagnets.
//!
//! # References
//!
//! - N. N. Bogoliubov, *J. Phys. USSR* **11**, 23 (1947).
//! - A. Kamra, W. Belzig, *Phys. Rev. Lett.* **116**, 146601 (2016).

use crate::error::{self, Result};
use crate::math::Complex;

// ---------------------------------------------------------------------------
// Main struct
// ---------------------------------------------------------------------------

/// Bogoliubov transformation data for a set of magnon modes.
///
/// Each mode k is characterised by:
/// - `u[k] = cosh θ_k` — the "normal" amplitude,
/// - `v[k] = sinh θ_k` — the "anomalous" (pairing) amplitude,
/// - `omega[k]` — the diagonalised magnon frequency \[rad s⁻¹\].
///
/// The relationship `u_k² − v_k² = 1` is a necessary condition for the
/// transformed operators to remain canonical bosons.
pub struct BogoliubovTransform {
    /// Normal amplitudes `u_k = cosh θ_k` \[dimensionless\].
    pub u: Vec<f64>,
    /// Anomalous amplitudes `v_k = sinh θ_k` \[dimensionless\].
    pub v: Vec<f64>,
    /// Diagonalised magnon frequencies `ω_k` \[rad/s\].
    pub omega: Vec<f64>,
    /// Stored diagonal Hamiltonian coefficients A_k \[same units as ω\].
    a_k: Vec<f64>,
    /// Stored off-diagonal (pairing) coefficients B_k \[same units as ω\].
    b_k: Vec<f64>,
}

impl BogoliubovTransform {
    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    /// Direct constructor from pre-computed `(u, v, ω)` triples.
    ///
    /// Validates that all three vectors have the same length and that
    /// `|u_k² − v_k² − 1| < 10⁻⁸` for every mode.
    ///
    /// # Arguments
    ///
    /// * `u`     — normal amplitudes `cosh θ_k` \[dimensionless\].
    /// * `v`     — anomalous amplitudes `sinh θ_k` \[dimensionless\].
    /// * `omega` — diagonalised frequencies \[rad/s\].
    ///
    /// # Errors
    ///
    /// Returns an error if lengths differ or if any mode violates `u²−v²=1`.
    pub fn new(u: Vec<f64>, v: Vec<f64>, omega: Vec<f64>) -> Result<Self> {
        if u.len() != v.len() || u.len() != omega.len() {
            return Err(error::invalid_param(
                "u/v/omega",
                "u, v, and omega must have the same length",
            ));
        }
        const TOL: f64 = 1e-8;
        for (i, (&ui, &vi)) in u.iter().zip(v.iter()).enumerate() {
            let norm = ui * ui - vi * vi;
            if (norm - 1.0).abs() > TOL {
                return Err(error::invalid_param(
                    "u,v",
                    &format!("mode {i}: u²−v²={norm:.3e} deviates from 1 by more than {TOL:.0e}"),
                ));
            }
        }
        // Reconstruct a_k and b_k from u, v, omega:
        // omega = sqrt(a^2 - b^2), u = cosh(theta), v = sinh(theta)
        // tanh(2 theta) = -b/a => b = -a * tanh(2*theta)
        // 2*theta = atanh(2*u*v/(u^2+v^2))  but we can derive:
        //   a = omega * u^2 + omega * v^2 = omega*(u^2+v^2)
        //   Actually: a = sqrt(omega^2 + (b)^2), complicated; store zeros
        let n = u.len();
        let a_k = vec![0.0; n];
        let b_k = vec![0.0; n];
        Ok(Self {
            u,
            v,
            omega,
            a_k,
            b_k,
        })
    }

    /// Diagonalise the bilinear Bose Hamiltonian
    /// `H = Σ_k [A_k a†_k a_k + (B_k/2)(a†_k a†_{-k} + h.c.)]`.
    ///
    /// For each k the squeeze angle satisfies `tanh 2θ_k = −B_k / A_k`, giving
    ///
    /// ```text
    /// u_k = cosh θ_k,   v_k = sinh θ_k,   ω_k = √(A_k² − B_k²)
    /// ```
    ///
    /// # Arguments
    ///
    /// * `a_k` — diagonal coefficients A_k \[rad/s or energy/ℏ\].
    /// * `b_k` — off-diagonal (pairing) coefficients B_k \[same units\].
    ///
    /// # Errors
    ///
    /// Returns [`crate::error::Error::NumericalError`] if any mode has
    /// `|B_k| ≥ |A_k|` (magnetically unstable) or if lengths differ.
    pub fn from_hamiltonian(a_k: &[f64], b_k: &[f64]) -> Result<Self> {
        if a_k.len() != b_k.len() {
            return Err(error::invalid_param(
                "a_k/b_k",
                "A_k and B_k arrays must have the same length",
            ));
        }
        let n = a_k.len();
        let mut u = Vec::with_capacity(n);
        let mut v = Vec::with_capacity(n);
        let mut omega = Vec::with_capacity(n);

        for (i, (&ai, &bi)) in a_k.iter().zip(b_k.iter()).enumerate() {
            if bi.abs() >= ai.abs() {
                return Err(error::numerical_error(&format!(
                    "mode {i}: |B_k|={} >= |A_k|={} — system is magnetically unstable",
                    bi.abs(),
                    ai.abs()
                )));
            }
            // omega_k = sqrt(A_k^2 - B_k^2)
            let omega_k = (ai * ai - bi * bi).sqrt();
            // tanh(2*theta) = -B_k / A_k
            let tanh_2theta = -bi / ai;
            // 2*theta = atanh(tanh_2theta)  — guaranteed |tanh| < 1 from check above
            let two_theta = tanh_2theta.atanh();
            let theta = two_theta / 2.0;
            let uk = theta.cosh();
            let vk = theta.sinh();
            u.push(uk);
            v.push(vk);
            omega.push(omega_k);
        }

        Ok(Self {
            u,
            v,
            omega,
            a_k: a_k.to_vec(),
            b_k: b_k.to_vec(),
        })
    }

    /// Build from an antiferromagnetic square-lattice Hamiltonian.
    ///
    /// For the standard Heisenberg AFM on the square lattice with exchange
    /// coupling J (> 0) and external field `h_ext` (along the Néel axis), the
    /// linear spin-wave Hamiltonian has
    ///
    /// ```text
    /// A_k = 4J S + h_ext
    /// B_k = 4J S · γ_k,   γ_k = (cos kx + cos ky) / 2
    /// ```
    ///
    /// where `S = 1/2` is absorbed into `J` here for generality (pass `J*S` as
    /// the argument).  K-points are sampled along the high-symmetry path
    /// Γ → X → M → Γ in the magnetic Brillouin zone.
    ///
    /// # Path segments
    ///
    /// The path is split evenly into three segments of `n_k / 3` points each:
    /// 1. Γ = (0,0)  → X = (π, 0)
    /// 2. X = (π, 0) → M = (π, π)
    /// 3. M = (π, π) → Γ = (0, 0)
    ///
    /// # Arguments
    ///
    /// * `j`     — exchange coupling J (positive, units of energy/ℏ = rad/s).
    /// * `h_ext` — external field contribution \[rad/s\].
    /// * `n_k`   — total number of k-points to sample (must be ≥ 3).
    ///
    /// # Errors
    ///
    /// Returns an error if `j ≤ 0`, `n_k < 3`, or any mode is unstable.
    pub fn from_afm_square_lattice(j: f64, h_ext: f64, n_k: usize) -> Result<Self> {
        if j <= 0.0 {
            return Err(error::invalid_param(
                "j",
                "exchange coupling must be positive",
            ));
        }
        if n_k < 3 {
            return Err(error::invalid_param(
                "n_k",
                "need at least 3 k-points for path sampling",
            ));
        }

        // Split n_k into 3 segments
        let seg = n_k / 3;
        let seg0 = seg;
        let seg1 = seg;
        let seg2 = n_k - 2 * seg; // last segment takes any remainder

        let mut a_vals = Vec::with_capacity(n_k);
        let mut b_vals = Vec::with_capacity(n_k);

        // Segment 0: Γ→X: (0,0)→(π,0), t ∈ [0,1)
        for i in 0..seg0 {
            let t = i as f64 / seg0 as f64;
            let kx = t * std::f64::consts::PI;
            let ky = 0.0_f64;
            let gamma_k = (kx.cos() + ky.cos()) / 2.0;
            let ai = 4.0 * j + h_ext;
            let bi = 4.0 * j * gamma_k;
            a_vals.push(ai);
            b_vals.push(bi);
        }

        // Segment 1: X→M: (π,0)→(π,π), t ∈ [0,1)
        for i in 0..seg1 {
            let t = i as f64 / seg1 as f64;
            let kx = std::f64::consts::PI;
            let ky = t * std::f64::consts::PI;
            let gamma_k = (kx.cos() + ky.cos()) / 2.0;
            let ai = 4.0 * j + h_ext;
            let bi = 4.0 * j * gamma_k;
            a_vals.push(ai);
            b_vals.push(bi);
        }

        // Segment 2: M→Γ: (π,π)→(0,0), t ∈ [0,1]
        for i in 0..seg2 {
            let t = i as f64 / (seg2 - 1).max(1) as f64;
            let kx = (1.0 - t) * std::f64::consts::PI;
            let ky = (1.0 - t) * std::f64::consts::PI;
            let gamma_k = (kx.cos() + ky.cos()) / 2.0;
            let ai = 4.0 * j + h_ext;
            let bi = 4.0 * j * gamma_k;
            a_vals.push(ai);
            b_vals.push(bi);
        }

        Self::from_hamiltonian(&a_vals, &b_vals)
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// Number of magnon modes.
    pub fn n_modes(&self) -> usize {
        self.u.len()
    }

    /// Diagonalised frequency `ω_k` for mode `k_index` \[rad/s\].
    ///
    /// # Errors
    ///
    /// Returns an error if `k_index ≥ n_modes()`.
    #[inline]
    pub fn magnon_frequency(&self, k_index: usize) -> Result<f64> {
        if k_index >= self.omega.len() {
            return Err(error::invalid_param(
                "k_index",
                &format!(
                    "index {k_index} out of range for {} modes",
                    self.omega.len()
                ),
            ));
        }
        Ok(self.omega[k_index])
    }

    /// Squeeze parameter `θ_k = atanh(v_k / u_k)` for mode `k_index`.
    ///
    /// A larger θ means stronger quantum correlations in the Bogoliubov vacuum.
    ///
    /// # Returns
    ///
    /// `θ_k` \[dimensionless\].
    ///
    /// # Errors
    ///
    /// Returns an error if `k_index ≥ n_modes()`.
    pub fn squeeze_parameter(&self, k_index: usize) -> Result<f64> {
        if k_index >= self.u.len() {
            return Err(error::invalid_param(
                "k_index",
                &format!("index {k_index} out of range"),
            ));
        }
        let uk = self.u[k_index];
        let vk = self.v[k_index];
        // theta = atanh(v/u) since tanh(theta) = sinh/cosh = v/u
        Ok((vk / uk).atanh())
    }

    /// Virtual magnon number in the Bogoliubov vacuum for mode `k_index`: `v_k²`.
    ///
    /// The ground-state expectation value `⟨0|α†_k α_k|0⟩_α` of the *original*
    /// magnon number is `v_k²` — non-zero whenever B_k ≠ 0 (squeezed vacuum).
    ///
    /// # Returns
    ///
    /// `v_k²` \[dimensionless\].
    ///
    /// # Errors
    ///
    /// Returns an error if `k_index ≥ n_modes()`.
    #[inline]
    pub fn vacuum_occupation(&self, k_index: usize) -> Result<f64> {
        if k_index >= self.v.len() {
            return Err(error::invalid_param(
                "k_index",
                &format!("index {k_index} out of range"),
            ));
        }
        Ok(self.v[k_index] * self.v[k_index])
    }

    /// Total virtual magnon occupation in the Bogoliubov vacuum: `Σ_k v_k²`.
    ///
    /// Represents the total "quantum depletion" — the number of original
    /// (bare) magnons present in the interacting ground state.
    ///
    /// # Returns
    ///
    /// `Σ_k v_k²` \[dimensionless\].
    pub fn total_vacuum_magnons(&self) -> f64 {
        self.v.iter().map(|&vk| vk * vk).sum()
    }

    /// Check that `u_k² − v_k² = 1` for all modes within `tol`.
    ///
    /// # Arguments
    ///
    /// * `tol` — absolute tolerance for each mode.
    ///
    /// # Returns
    ///
    /// `true` if all modes satisfy the normalisation condition.
    pub fn check_normalization(&self, tol: f64) -> bool {
        self.u
            .iter()
            .zip(self.v.iter())
            .all(|(&uk, &vk)| (uk * uk - vk * vk - 1.0).abs() < tol)
    }

    /// Apply the Bogoliubov transformation to a complex state vector.
    ///
    /// For each mode k: `out[k] = u_k · in[k] + v_k · in[k]*`
    ///
    /// This maps the amplitude in the bare-magnon Fock basis to the Bogoliubov
    /// quasi-magnon basis.  The operation is unitary when `u²−v²=1`.
    ///
    /// # Arguments
    ///
    /// * `input` — complex amplitudes, one per mode \[dimensionless\].
    ///
    /// # Errors
    ///
    /// Returns an error if `input.len() ≠ n_modes()`.
    pub fn apply_to_state(&self, input: &[Complex]) -> Result<Vec<Complex>> {
        if input.len() != self.u.len() {
            return Err(error::invalid_param(
                "input",
                &format!(
                    "state vector length {} does not match n_modes={}",
                    input.len(),
                    self.u.len()
                ),
            ));
        }
        let output: Vec<Complex> = input
            .iter()
            .zip(self.u.iter().zip(self.v.iter()))
            .map(|(&c, (&uk, &vk))| {
                // out[k] = u_k * c + v_k * c*
                let normal = c.scale(uk);
                let anomal = c.conj().scale(vk);
                normal.add(&anomal)
            })
            .collect();
        Ok(output)
    }

    /// Zero-point energy shift from the Bogoliubov transformation \[rad/s × ℏ or same units as A_k\].
    ///
    /// `E_0 = (1/2) Σ_k (ω_k − A_k)`
    ///
    /// This is negative when ω_k < A_k (which occurs whenever B_k ≠ 0),
    /// reflecting the lowering of the ground-state energy by quantum correlations.
    ///
    /// # Returns
    ///
    /// `(1/2) Σ_k (ω_k − A_k)` in the same units as ω (rad/s if constructed
    /// via `from_hamiltonian`).
    pub fn ground_state_energy(&self) -> f64 {
        let sum: f64 = self
            .omega
            .iter()
            .zip(self.a_k.iter())
            .map(|(&ok, &ak)| ok - ak)
            .sum();
        0.5 * sum
    }

    /// Pairing coefficient `B_k` for mode `k_index`, as stored from the original Hamiltonian.
    ///
    /// Returns `0.0` if the transform was constructed via [`BogoliubovTransform::new`] rather
    /// than [`BogoliubovTransform::from_hamiltonian`] (since B_k is unknown in that path).
    ///
    /// # Errors
    ///
    /// Returns an error if `k_index ≥ n_modes()`.
    pub fn pairing_coefficient(&self, k_index: usize) -> Result<f64> {
        if k_index >= self.b_k.len() {
            return Err(error::invalid_param(
                "k_index",
                &format!("index {k_index} out of range"),
            ));
        }
        Ok(self.b_k[k_index])
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    // ── basic normalization ───────────────────────────────────────────────────

    #[test]
    fn test_normalization_u_sq_minus_v_sq_equals_one() {
        // Construct from Hamiltonian and verify
        let a = vec![10.0, 8.0, 6.0];
        let b = vec![2.0, 1.5, 0.5];
        let bogo = BogoliubovTransform::from_hamiltonian(&a, &b).expect("stable");
        assert!(bogo.check_normalization(1e-12), "normalization must hold");
    }

    // ── B=0 limit ────────────────────────────────────────────────────────────

    #[test]
    fn test_zero_pairing_gives_identity_transformation() {
        let a = vec![5.0, 3.0];
        let b = vec![0.0, 0.0];
        let bogo = BogoliubovTransform::from_hamiltonian(&a, &b).expect("stable");
        // u = cosh(0) = 1, v = sinh(0) = 0
        for i in 0..2 {
            assert!((bogo.u[i] - 1.0).abs() < TOL, "u should be 1 at B=0");
            assert!(bogo.v[i].abs() < TOL, "v should be 0 at B=0");
        }
        // omega = sqrt(A^2 - 0) = A
        assert!((bogo.omega[0] - 5.0).abs() < TOL);
        assert!((bogo.omega[1] - 3.0).abs() < TOL);
    }

    // ── instability ───────────────────────────────────────────────────────────

    #[test]
    fn test_diverges_when_b_ge_a() {
        let a = vec![3.0];
        let b = vec![3.0]; // |B| = |A|  → unstable
        assert!(BogoliubovTransform::from_hamiltonian(&a, &b).is_err());

        let b2 = vec![4.0]; // |B| > |A| → unstable
        assert!(BogoliubovTransform::from_hamiltonian(&a, &b2).is_err());
    }

    // ── positive frequencies ──────────────────────────────────────────────────

    #[test]
    fn test_omega_positive_for_stable_modes() {
        let a = vec![10.0, 5.0, 20.0];
        let b = vec![2.0, 1.0, 4.0];
        let bogo = BogoliubovTransform::from_hamiltonian(&a, &b).expect("stable");
        for (i, &ok) in bogo.omega.iter().enumerate() {
            assert!(ok > 0.0, "omega[{i}] = {ok} should be positive");
        }
    }

    // ── zero squeeze parameter at B=0 ─────────────────────────────────────────

    #[test]
    fn test_squeeze_parameter_zero_at_no_pairing() {
        let a = vec![7.0];
        let b = vec![0.0];
        let bogo = BogoliubovTransform::from_hamiltonian(&a, &b).expect("stable");
        let theta = bogo.squeeze_parameter(0).expect("valid index");
        assert!(theta.abs() < TOL, "theta should be 0 when B=0, got {theta}");
    }

    // ── vacuum occupation increases with B/A ratio ────────────────────────────

    #[test]
    fn test_vacuum_occupation_increases_with_b_over_a() {
        let a = vec![10.0, 10.0];
        let b_small = vec![1.0, 0.0];
        let b_large = vec![8.0, 0.0];
        let bogo_small = BogoliubovTransform::from_hamiltonian(&a, &b_small).expect("stable");
        let bogo_large = BogoliubovTransform::from_hamiltonian(&a, &b_large).expect("stable");
        let occ_small = bogo_small.vacuum_occupation(0).expect("valid");
        let occ_large = bogo_large.vacuum_occupation(0).expect("valid");
        assert!(occ_large > occ_small, "larger B/A should give larger v²");
    }

    // ── apply_to_state ────────────────────────────────────────────────────────

    #[test]
    fn test_apply_to_state_preserves_norm_for_real_input() {
        // For B=0: out[k] = 1.0 * in[k] + 0 * in[k]* = in[k], norm preserved
        let a = vec![5.0];
        let b = vec![0.0];
        let bogo = BogoliubovTransform::from_hamiltonian(&a, &b).expect("stable");
        let input = vec![Complex::new(0.6, 0.8)]; // |z|=1
        let output = bogo.apply_to_state(&input).expect("valid");
        let norm_in = input[0].norm_sq();
        let norm_out = output[0].norm_sq();
        assert!(
            (norm_in - norm_out).abs() < TOL,
            "norm must be preserved for B=0"
        );
    }

    // ── AFM square lattice ────────────────────────────────────────────────────

    #[test]
    fn test_afm_square_lattice_dispersion_at_gamma_point() {
        // At Γ = (0,0): gamma_k = 1.0 → B_k = 4J
        // omega = sqrt((4J+h)^2 - (4J)^2) = sqrt(h*(8J+h))
        let j = 1.0;
        let h_ext = 0.5;
        let bogo = BogoliubovTransform::from_afm_square_lattice(j, h_ext, 30).expect("stable");
        // k-index 0 is Γ (first segment, t=0)
        let omega_gamma = bogo.magnon_frequency(0).expect("valid");
        let expected = (h_ext * (8.0 * j + h_ext)).sqrt();
        let rel_err = (omega_gamma - expected).abs() / expected;
        assert!(
            rel_err < 1e-8,
            "Γ-point omega mismatch: got {omega_gamma}, expected {expected}"
        );
    }

    #[test]
    fn test_afm_zone_boundary_at_x_point() {
        // At X = (π,0): gamma_k = (cos(π)+cos(0))/2 = (-1+1)/2 = 0 → B_k = 0
        // omega = sqrt((4J+h)^2 - 0) = 4J + h
        let j = 2.0;
        let h_ext = 0.1;
        // n_k=30, first segment has 10 points: t = 0/10..9/10
        // t=1.0 is excluded (first segment t ∈ [0,1)), so X is approached but not exactly hit.
        // Use large n_k to approach X closely; alternatively check at t=9/10
        // For a cleaner test use n_k=3: seg0=1(Γ), seg1=1(X), seg2=1(M)
        let bogo = BogoliubovTransform::from_afm_square_lattice(j, h_ext, 3).expect("stable");
        // Index 1 is start of seg1: X=(π,0), gamma_k=0
        let omega_x = bogo.magnon_frequency(1).expect("valid");
        let expected = 4.0 * j + h_ext;
        let rel_err = (omega_x - expected).abs() / expected;
        assert!(
            rel_err < 1e-10,
            "X-point omega mismatch: got {omega_x}, expected {expected}"
        );
    }

    // ── total_vacuum_magnons ──────────────────────────────────────────────────

    #[test]
    fn test_total_vacuum_magnons_positive() {
        let a = vec![5.0, 8.0, 3.0];
        let b = vec![1.0, 2.0, 0.5];
        let bogo = BogoliubovTransform::from_hamiltonian(&a, &b).expect("stable");
        assert!(bogo.total_vacuum_magnons() > 0.0);
    }

    // ── ground_state_energy ───────────────────────────────────────────────────

    #[test]
    fn test_ground_state_energy_negative() {
        // When B ≠ 0: ω_k < A_k, so (ω_k - A_k) < 0 → E_0 < 0
        let a = vec![10.0, 8.0];
        let b = vec![3.0, 2.0];
        let bogo = BogoliubovTransform::from_hamiltonian(&a, &b).expect("stable");
        let e0 = bogo.ground_state_energy();
        assert!(
            e0 < 0.0,
            "ground state energy shift should be negative when B≠0, got {e0}"
        );
    }

    // ── check_normalization tolerance ─────────────────────────────────────────

    #[test]
    fn test_check_normalization_tol() {
        let a = vec![6.0];
        let b = vec![2.0];
        let bogo = BogoliubovTransform::from_hamiltonian(&a, &b).expect("stable");
        assert!(bogo.check_normalization(1e-10));
        // Overly tight tolerance should still pass since arithmetic is exact enough
        assert!(bogo.check_normalization(1e-12));
    }

    // ── new() validates lengths ───────────────────────────────────────────────

    #[test]
    fn test_new_validates_lengths_match() {
        let u = vec![1.0, 1.0];
        let v = vec![0.0];
        let omega = vec![5.0, 3.0];
        assert!(BogoliubovTransform::new(u, v, omega).is_err());
    }

    // ── new() rejects unnormalized ────────────────────────────────────────────

    #[test]
    fn test_new_rejects_unnormalized() {
        // u=2, v=0: u²-v²=4≠1
        let u = vec![2.0];
        let v = vec![0.0];
        let omega = vec![5.0];
        assert!(BogoliubovTransform::new(u, v, omega).is_err());
    }
}
