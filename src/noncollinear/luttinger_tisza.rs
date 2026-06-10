//! Luttinger-Tisza method for classical spin ground states
//!
//! The Luttinger-Tisza (LT) method finds the classical ground state of the
//! Heisenberg model H = Σ_{ij} J_{ij} **S**_i · **S**_j by diagonalizing the
//! Fourier-transformed exchange matrix J(**q**) and finding the ordering wavevector
//! **q*** that minimizes the spin energy.
//!
//! # Physics Background
//!
//! ## Fourier Decomposition
//!
//! For a Bravais lattice with periodic exchange interactions, the Hamiltonian in
//! Fourier space is:
//!
//! H = Σ_**q** J(**q**) |**S**(**q**)|²
//!
//! where
//!
//! J(**q**) = Σ_{**δ**} J_{0**δ**} exp(i **q** · **δ**)
//!
//! and the sum runs over all neighbor shells. For a centrosymmetric lattice
//! (where J_{**δ**} = J_{−**δ**}) this is purely real:
//!
//! J(**q**) = Σ_{**δ**} J_{0**δ**} cos(**q** · **δ**)
//!
//! ## Ground State Search
//!
//! The classical ground state minimizes the energy. For the convention
//! H = Σ J_{ij} **S**_i · **S**_j with J > 0 ferromagnetic:
//! - FM (J > 0): J(**q** = 0) is maximum (most positive), so the energy −J(**q**)S²
//!   is most negative at **q** = 0.
//! - AFM (J < 0, single NN): J(**q** = π/a) has the most negative value.
//! - Frustrated (J₁ > 0, J₂ < 0 with |J₂| > J₁/4): intermediate **q*** ≠ 0, π/a.
//!
//! The optimal ordering wavevector minimizes J(**q**), i.e., finds the most negative
//! (lowest energy) Fourier component.
//!
//! ## J₁-J₂ Chain
//!
//! The frustrated J₁-J₂ Heisenberg chain has first-neighbor coupling J₁ > 0
//! (ferromagnetic) and second-neighbor coupling J₂ (typically competing):
//!
//! J(q) = 2J₁ cos(qa) + 2J₂ cos(2qa)
//!
//! The minimum of J(q) moves from q = π/a (AFM, J₂ = 0) to an incommensurate
//! value when J₂ > J₁/4:
//!
//! q* = (1/a) arccos(−J₁ / (4J₂))
//!
//! This is the Lifshitz point separating the commensurate AFM from the spiral phase.
//!
//! ## Mean-Field Ordering Temperature
//!
//! Within mean-field theory the ordering temperature satisfies
//!
//! k_B T_c = S(S+1) |J(**q***)|  / 3
//!
//! which gives a rough but useful estimate of T_c from the exchange parameters.
//!
//! # Sign Convention
//!
//! `ExchangeInteraction::j_exchange > 0` means ferromagnetic:
//! - The Hamiltonian is H = Σ J S_i · S_j so parallel spins (FM) have energy −J < 0
//! - J(**q** = 0) = Σ_δ J_δ = z J for coordination number z and uniform J
//! - The minimum energy state is at the **q** that minimizes J(**q**), which for FM is q = 0
//!
//! # References
//!
//! - J. M. Luttinger and L. Tisza, "Theory of Dipole Interaction in Crystals",
//!   *Phys. Rev.* **70**, 954–964 (1946)
//! - J. Villain, "La structure des substances magnétiques",
//!   *J. Phys. Chem. Solids* **11**, 303–309 (1959)
//! - R. J. Elliott, "Phenomenological Discussion of Magnetic Ordering in the Heavy
//!   Rare-Earth Metals", *Phys. Rev.* **124**, 346–353 (1961)

use std::f64::consts::TAU;

use crate::constants::KB;
use crate::error::{self, Result};
use crate::vector3::Vector3;

/// A single exchange interaction between site 0 and a neighbor at displacement **δ**
///
/// For a centrosymmetric lattice the pair (**δ**, −**δ**) both appear with the same J,
/// so the LT Fourier sum automatically gives a real J(**q**).
#[derive(Debug, Clone)]
pub struct ExchangeInteraction {
    /// Lattice displacement vector to the neighbor \[m\]
    pub delta: Vector3<f64>,
    /// Exchange coupling constant \[J\]; J > 0 is ferromagnetic, J < 0 is antiferromagnetic.
    pub j_exchange: f64,
}

impl ExchangeInteraction {
    /// Construct an exchange interaction.
    ///
    /// # Arguments
    ///
    /// - `delta` — displacement vector to the neighbor \[m\]
    /// - `j_exchange` — coupling constant \[J\]
    pub fn new(delta: Vector3<f64>, j_exchange: f64) -> Self {
        Self { delta, j_exchange }
    }
}

/// Luttinger-Tisza solver for classical spin-lattice ground states
///
/// Finds the ordering wavevector **q*** of the classical ground state of the
/// Heisenberg model H = Σ_{ij} J_{ij} **S**_i · **S**_j by minimizing the
/// Fourier-space exchange matrix J(**q**) over the first Brillouin zone.
#[derive(Debug, Clone)]
pub struct LuttingerTisza {
    /// List of exchange interactions (neighbor vectors and coupling constants)
    interactions: Vec<ExchangeInteraction>,
    /// Spin quantum number S (e.g. 0.5, 1.0, 2.0)
    spin: f64,
    /// Lattice constant of the underlying Bravais lattice \[m\]
    a_lattice: f64,
}

impl LuttingerTisza {
    /// Construct a LT solver with explicit interaction list.
    ///
    /// # Errors
    ///
    /// Returns an error if `spin` ≤ 0 or `a_lattice` ≤ 0 or the interaction list
    /// is empty.
    pub fn new(interactions: Vec<ExchangeInteraction>, spin: f64, a_lattice: f64) -> Result<Self> {
        if spin <= 0.0 {
            return Err(error::invalid_param(
                "spin",
                "spin quantum number must be positive",
            ));
        }
        if a_lattice <= 0.0 {
            return Err(error::invalid_param(
                "a_lattice",
                "lattice constant must be positive",
            ));
        }
        if interactions.is_empty() {
            return Err(error::invalid_param(
                "interactions",
                "at least one exchange interaction is required",
            ));
        }
        Ok(Self {
            interactions,
            spin,
            a_lattice,
        })
    }

    /// Ferromagnet: single nearest-neighbor interaction on a simple chain.
    ///
    /// Sets up a 1D chain with NN coupling `j_nn` > 0 (ferromagnetic convention).
    /// The ground state is FM with **q*** = 0.
    ///
    /// # Arguments
    ///
    /// - `j_nn` — NN exchange constant \[J\]; should be positive for FM
    /// - `a_lattice` — lattice constant \[m\]
    /// - `spin` — spin quantum number S
    pub fn ferromagnet(j_nn: f64, a_lattice: f64, spin: f64) -> Self {
        let interactions = vec![
            ExchangeInteraction::new(Vector3::new(a_lattice, 0.0, 0.0), j_nn),
            ExchangeInteraction::new(Vector3::new(-a_lattice, 0.0, 0.0), j_nn),
        ];
        Self {
            interactions,
            spin,
            a_lattice,
        }
    }

    /// Antiferromagnet: single nearest-neighbor interaction on a simple chain.
    ///
    /// Sets up a 1D chain with NN coupling `j_nn` < 0 (antiferromagnetic convention).
    /// The ground state is AFM with **q*** = π/a.
    ///
    /// # Arguments
    ///
    /// - `j_nn` — NN exchange constant \[J\]; should be negative for AFM
    /// - `a_lattice` — lattice constant \[m\]
    /// - `spin` — spin quantum number S
    pub fn antiferromagnet(j_nn: f64, a_lattice: f64, spin: f64) -> Self {
        let interactions = vec![
            ExchangeInteraction::new(Vector3::new(a_lattice, 0.0, 0.0), j_nn),
            ExchangeInteraction::new(Vector3::new(-a_lattice, 0.0, 0.0), j_nn),
        ];
        Self {
            interactions,
            spin,
            a_lattice,
        }
    }

    /// J₁-J₂ frustrated chain: first- and second-neighbor exchange.
    ///
    /// Implements the 1D J₁-J₂ Heisenberg model:
    ///
    /// H = J₁ Σ_i **S**_i · **S**_{i+1} + J₂ Σ_i **S**_i · **S**_{i+2}
    ///
    /// The frustration condition is J₂ > J₁ / 4 (for J₁ > 0, J₂ < 0):
    /// when satisfied, the ground state is a spiral with
    ///
    /// q* = (1/a) arccos(−J₁ / (4J₂))
    ///
    /// # Arguments
    ///
    /// - `j1` — NN coupling \[J\]; positive = FM
    /// - `j2` — NNN coupling \[J\]; negative = AFM (competing)
    /// - `a_lattice` — lattice constant \[m\]
    /// - `spin` — spin quantum number S
    pub fn j1j2_chain(j1: f64, j2: f64, a_lattice: f64, spin: f64) -> Self {
        let a = a_lattice;
        let interactions = vec![
            // NN interactions
            ExchangeInteraction::new(Vector3::new(a, 0.0, 0.0), j1),
            ExchangeInteraction::new(Vector3::new(-a, 0.0, 0.0), j1),
            // NNN interactions
            ExchangeInteraction::new(Vector3::new(2.0 * a, 0.0, 0.0), j2),
            ExchangeInteraction::new(Vector3::new(-2.0 * a, 0.0, 0.0), j2),
        ];
        Self {
            interactions,
            spin,
            a_lattice,
        }
    }

    /// Evaluate the Fourier-transformed exchange matrix J(**q**).
    ///
    /// For a centrosymmetric lattice (J_{**δ**} = J_{-**δ**}):
    ///
    /// J(**q**) = Σ_{**δ**} J_{0**δ**} cos(**q** · **δ**)
    ///
    /// This is the (real, scalar) exchange eigenvalue for a single-**q** spiral.
    /// The spin energy per site is E = −J(**q**) S² (minus sign because
    /// energy of parallel spins in FM is −J S²).
    ///
    /// **Note**: The sign convention is that the minimum of J(**q**) corresponds to
    /// the lowest-energy ordering (most negative energy per site = −min(J) × S²).
    ///
    /// # Arguments
    ///
    /// - `q` — wavevector in Cartesian coordinates \[1/m\]
    pub fn exchange_fourier(&self, q: &Vector3<f64>) -> f64 {
        self.interactions
            .iter()
            .map(|inter| inter.j_exchange * q.dot(&inter.delta).cos())
            .sum()
    }

    /// Find the ordering wavevector **q*** that maximizes J(**q**) over the BZ.
    ///
    /// Performs a systematic grid search over the first Brillouin zone
    /// [0, 2π/a) × [0, 2π/a) × [0, 2π/a) at n_q³ wavevectors and returns
    /// the **q** vector at which J(**q**) is maximum (most positive).
    ///
    /// The classical ground state energy per spin is E = −J(**q***) S², so the
    /// ground state minimizes E by **maximizing** J(**q**):
    ///
    /// - FM (J > 0): J(**q** = 0) = z J (maximum) → E = −z J S² < 0 ✓
    /// - AFM (J < 0): J(**q** = π/a) = −2J = 2|J| (maximum) → E = 2J S² < 0 ✓
    /// - Spiral: intermediate **q*** where J has an intermediate maximum
    ///
    /// # Arguments
    ///
    /// - `n_q` — number of grid points along each BZ axis (e.g. 100)
    ///
    /// # Note
    ///
    /// For a 1D chain only the q_x component is meaningful; q_y and q_z are
    /// set to zero via the interaction list (no y or z neighbors).
    pub fn find_ground_state_q(&self, n_q: usize) -> Vector3<f64> {
        let n = n_q.max(2);
        let bz_size = TAU / self.a_lattice;
        let dq = bz_size / n as f64;

        let mut max_j = f64::NEG_INFINITY;
        let mut best_q = Vector3::zero();

        // Determine which axes carry interactions
        let has_y = self.interactions.iter().any(|i| i.delta.y.abs() > 1e-30);
        let has_z = self.interactions.iter().any(|i| i.delta.z.abs() > 1e-30);

        let ny = if has_y { n } else { 1 };
        let nz = if has_z { n } else { 1 };

        for ix in 0..n {
            let qx = ix as f64 * dq;
            for iy in 0..ny {
                let qy = if has_y { iy as f64 * dq } else { 0.0 };
                for iz in 0..nz {
                    let qz = if has_z { iz as f64 * dq } else { 0.0 };
                    let q = Vector3::new(qx, qy, qz);
                    let j_q = self.exchange_fourier(&q);
                    if j_q > max_j {
                        max_j = j_q;
                        best_q = q;
                    }
                }
            }
        }

        best_q
    }

    /// Classical energy per spin for a given ordering wavevector.
    ///
    /// E/spin = −J(**q**) × S²
    ///
    /// where J(**q**) is the Fourier exchange. For FM (J > 0, **q** = 0):
    /// E = −z J S² < 0 (negative, correctly stabilizing for FM order).
    ///
    /// # Arguments
    ///
    /// - `q` — ordering wavevector \[1/m\]
    pub fn ground_state_energy_per_spin(&self, q: &Vector3<f64>) -> f64 {
        let j_q = self.exchange_fourier(q);
        -j_q * self.spin * self.spin
    }

    /// Returns `true` if the ground state is a noncollinear (spiral) ordering.
    ///
    /// Finds the optimal **q*** (maximizing J(**q**)) and checks that it is neither
    /// **q** = 0 (ferromagnet) nor at the BZ boundary (antiferromagnet).
    /// The tolerance `tol` (as a fraction of ξ = |**q***| a / (2π)) selects
    /// how close to 0 or 0.5 qualifies as commensurate.
    ///
    /// For the J₁-J₂ chain the analytic crossover is at J₂/J₁ = −1/4 (competing NNN).
    pub fn is_spiral_ground_state(&self, tol: f64) -> bool {
        let q_star = self.find_ground_state_q(200);
        let q_x = q_star.x;
        let pi_over_a = std::f64::consts::PI / self.a_lattice;

        // Check each active component: for 1D chain only q_x matters
        // ξ_x = q_x / (2π/a) = q_x × a / (2π)
        let xi_x = q_x * self.a_lattice / TAU;

        // FM: xi_x ≈ 0; AFM (1D): xi_x ≈ 0.5; spiral: in between
        let not_fm = xi_x > tol;
        let not_afm = (q_x - pi_over_a).abs() > tol * pi_over_a;

        not_fm && not_afm
    }

    /// Estimate the magnetic ordering temperature via mean-field theory.
    ///
    /// T_c ≈ S(S+1) |J(**q***)|  / (3 k_B)
    ///
    /// where **q*** is the ground-state ordering wavevector. This overestimates T_c
    /// for low-dimensional or frustrated systems but gives a useful order-of-magnitude.
    pub fn ordering_temperature_estimate(&self) -> f64 {
        let q_star = self.find_ground_state_q(100);
        let j_q_star = self.exchange_fourier(&q_star).abs();
        self.spin * (self.spin + 1.0) * j_q_star / (3.0 * KB)
    }

    /// Spiral pitch angle [rad/site] for a given optimal wavevector.
    ///
    /// Returns |**q***| × a_lattice, which is the rotation angle (in radians) of the
    /// spin from one lattice site to the next along **q***.
    ///
    /// - For FM: 0 rad/site (no rotation)
    /// - For AFM (1D): π rad/site (full 180° flip)
    /// - For spiral: intermediate angle
    pub fn spiral_pitch_angle(&self, q_optimal: &Vector3<f64>) -> f64 {
        q_optimal.magnitude() * self.a_lattice
    }

    /// Frustration ratio J₂/J₁ for a J₁-J₂ chain.
    ///
    /// Computes the ratio of the NNN to NN exchange magnitude from the stored
    /// interaction list. For models not built with `j1j2_chain()` this may not
    /// be meaningful; it computes the ratio of the longest-range to shortest-range
    /// interaction magnitudes.
    ///
    /// A value > 0.25 indicates classical frustration (spiral ground state) for the
    /// J₁-J₂ chain with J₁ > 0 (FM) and J₂ < 0 (competing AFM NNN).
    pub fn frustration_ratio(&self) -> f64 {
        // Find largest and smallest |delta| among interactions
        let mut min_dist = f64::INFINITY;
        let mut max_dist = 0.0_f64;

        for inter in &self.interactions {
            let d = inter.delta.magnitude();
            if d > 1e-30 {
                min_dist = min_dist.min(d);
                max_dist = max_dist.max(d);
            }
        }

        if max_dist <= min_dist * (1.0 + 1e-6) {
            // All interactions at same range — ratio undefined, return 0
            return 0.0;
        }

        // J1: average exchange for shortest-range interactions
        // J2: average exchange for longest-range interactions
        let j1_avg: f64 = self
            .interactions
            .iter()
            .filter(|i| (i.delta.magnitude() - min_dist).abs() < 1e-10 * min_dist + 1e-30)
            .map(|i| i.j_exchange)
            .sum::<f64>()
            / self
                .interactions
                .iter()
                .filter(|i| (i.delta.magnitude() - min_dist).abs() < 1e-10 * min_dist + 1e-30)
                .count() as f64;

        let j2_avg: f64 = self
            .interactions
            .iter()
            .filter(|i| (i.delta.magnitude() - max_dist).abs() < 1e-10 * max_dist + 1e-30)
            .map(|i| i.j_exchange)
            .sum::<f64>()
            / self
                .interactions
                .iter()
                .filter(|i| (i.delta.magnitude() - max_dist).abs() < 1e-10 * max_dist + 1e-30)
                .count() as f64;

        if j1_avg.abs() < 1e-30 {
            return 0.0;
        }
        j2_avg / j1_avg
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    const TOL: f64 = 1e-6;

    // Test 1: Ferromagnet (J>0): ground state at q=0.
    // Physics: E/spin = -J(q)S², J(q=0) = 2J is maximum → E = -2JS² is minimum.
    #[test]
    fn test_fm_optimal_q_is_zero() {
        let a = 3.0e-10_f64;
        let lt = LuttingerTisza::ferromagnet(1.0e-21, a, 0.5);
        let q_star = lt.find_ground_state_q(200);
        let q_mag = q_star.magnitude();
        let bz = TAU / a;
        let xi = q_mag / bz; // reduced wavevector: 0 = FM, 0.5 = AFM
        assert!(
            xi < 0.05,
            "FM (J>0) ground state should be at q≈0 (xi < 0.05), got xi = {}",
            xi
        );
    }

    // Test 2: Antiferromagnet (J<0): ground state at q=π/a.
    // Physics: J(q=π/a) = -2J = 2|J| (maximum for J<0) → E = 2JS² < 0 (minimum).
    #[test]
    fn test_afm_optimal_q_is_pi_over_a() {
        let a = 3.0e-10_f64;
        let lt = LuttingerTisza::antiferromagnet(-1.0e-21, a, 0.5);
        let q_star = lt.find_ground_state_q(200);
        let q_x = q_star.x;
        let pi_over_a = PI / a;
        let rel_err = (q_x - pi_over_a).abs() / pi_over_a;
        assert!(
            rel_err < 0.05,
            "AFM (J<0) ground state should be at q≈π/a, got q_x = {} (expected {})",
            q_x,
            pi_over_a
        );
    }

    // Test 3: J1-J2 chain with J2=-0.5J1 (|J2|/J1 = 0.5 > 0.25): spiral ground state.
    // The Lifshitz point is at |J2|/J1 = 1/4; above this q* moves incommensurately.
    #[test]
    fn test_j1j2_spiral_ground_state() {
        let a = 3.0e-10_f64;
        let j1 = 1.0e-21_f64;
        let j2 = -0.5 * j1; // competing NNN: |J2|/J1 = 0.5 > 0.25 → spiral
        let lt = LuttingerTisza::j1j2_chain(j1, j2, a, 0.5);
        assert!(
            lt.is_spiral_ground_state(0.05),
            "J1-J2 chain with |J2|/J1=0.5 must have a spiral ground state"
        );
    }

    // Test 4: Ordering temperature estimate is positive and physically reasonable
    #[test]
    fn test_ordering_temperature_positive() {
        let a = 3.0e-10_f64;
        let j_nn = 1.0e-21_f64; // ~72 K / 6 NN ≈ 12 K NN exchange
        let lt = LuttingerTisza::ferromagnet(j_nn, a, 0.5);
        let t_c = lt.ordering_temperature_estimate();
        assert!(t_c > 0.0, "ordering temperature must be positive");
        // For J=1e-21 J, S=0.5: T_c ≈ S(S+1)|J(q)|/(3kB) = 0.5*1.5*2e-21/(3*1.38e-23) ≈ 36 K
        assert!(
            t_c < 1e6,
            "ordering temperature should be physically reasonable, got {} K",
            t_c
        );
    }

    // Test 5: frustration_ratio = J2/J1 correctly
    #[test]
    fn test_frustration_ratio() {
        let a = 3.0e-10_f64;
        let j1 = 1.0e-21_f64;
        let j2 = -0.3 * j1;
        let lt = LuttingerTisza::j1j2_chain(j1, j2, a, 0.5);
        let ratio = lt.frustration_ratio();
        let expected = j2 / j1;
        assert!(
            (ratio - expected).abs() < TOL,
            "frustration ratio = {}, expected {}",
            ratio,
            expected
        );
    }

    // Test 6: exchange_fourier at q=0 equals sum of all J (= z×J for uniform NN)
    #[test]
    fn test_exchange_fourier_at_q_zero() {
        let a = 3.0e-10_f64;
        let j_nn = 2.0e-21_f64;
        let lt = LuttingerTisza::ferromagnet(j_nn, a, 0.5);
        let q0 = Vector3::zero();
        let j_at_0 = lt.exchange_fourier(&q0);
        // For 1D chain with 2 NN: J(0) = 2 × j_nn
        let expected = 2.0 * j_nn;
        assert!(
            (j_at_0 - expected).abs() < TOL * j_nn,
            "J(0) = {}, expected {}",
            j_at_0,
            expected
        );
    }

    // Test 7: ground_state_energy_per_spin is negative for FM (J>0, q=0)
    #[test]
    fn test_fm_energy_negative() {
        let a = 3.0e-10_f64;
        let j_nn = 1.0e-21_f64;
        let lt = LuttingerTisza::ferromagnet(j_nn, a, 0.5);
        let q0 = Vector3::zero();
        let energy = lt.ground_state_energy_per_spin(&q0);
        assert!(
            energy < 0.0,
            "FM ground state energy must be negative, got {}",
            energy
        );
    }

    // Test 8: spiral pitch angle for FM is ~0, for AFM is ~π
    #[test]
    fn test_spiral_pitch_angle() {
        let a = 3.0e-10_f64;
        let lt = LuttingerTisza::ferromagnet(1.0e-21, a, 0.5);
        let q_fm = Vector3::zero();
        let pitch_fm = lt.spiral_pitch_angle(&q_fm);
        assert!(
            pitch_fm.abs() < TOL,
            "FM pitch angle should be 0, got {}",
            pitch_fm
        );

        let q_afm = Vector3::new(PI / a, 0.0, 0.0);
        let pitch_afm = lt.spiral_pitch_angle(&q_afm);
        assert!(
            (pitch_afm - PI).abs() < TOL,
            "AFM pitch angle should be π, got {}",
            pitch_afm
        );
    }

    // Test 9: new() validates parameters correctly
    #[test]
    fn test_new_validates_params() {
        let inter = ExchangeInteraction::new(Vector3::new(3.0e-10, 0.0, 0.0), 1.0e-21);
        assert!(LuttingerTisza::new(vec![inter.clone()], -1.0, 3.0e-10).is_err()); // negative spin
        assert!(LuttingerTisza::new(vec![inter.clone()], 0.5, -3.0e-10).is_err()); // negative a
        assert!(LuttingerTisza::new(vec![], 0.5, 3.0e-10).is_err()); // empty interactions
        assert!(LuttingerTisza::new(vec![inter], 0.5, 3.0e-10).is_ok());
    }

    // Test 10: J1-J2 chain exchange_fourier matches analytic formula
    #[test]
    fn test_j1j2_exchange_fourier_analytic() {
        let a = 3.0e-10_f64;
        let j1 = 1.0e-21_f64;
        let j2 = -0.3e-21_f64;
        let lt = LuttingerTisza::j1j2_chain(j1, j2, a, 0.5);

        // Test at specific q: J(q) = 2J1 cos(qa) + 2J2 cos(2qa)
        let qa = 0.4_f64; // some angle in radians
        let q = Vector3::new(qa / a, 0.0, 0.0);
        let j_numerical = lt.exchange_fourier(&q);
        let j_analytic = 2.0 * j1 * qa.cos() + 2.0 * j2 * (2.0 * qa).cos();
        assert!(
            (j_numerical - j_analytic).abs() < TOL * j1.abs(),
            "J(q) numerical = {}, analytic = {}",
            j_numerical,
            j_analytic
        );
    }
}
