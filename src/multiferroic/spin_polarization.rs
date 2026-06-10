//! Spin-current (KNB) and inverse magnetoelectric mechanisms in spiral magnets
//!
//! This module implements the Katsura-Nagaosa-Balatsky (KNB) spin-current model for
//! magnetically-induced ferroelectricity in spiral magnets (Type-II multiferroics),
//! along with the inverse magnetoelectric effect that allows electric-field control
//! of magnetization.
//!
//! ## KNB Spin-Current Mechanism
//!
//! In a cycloidal spin spiral the spins rotate within a plane containing the
//! propagation vector **q**.  The cross product S_i × S_j — the "spin current" —
//! is perpendicular to both spins and, when combined with the bond direction e_ij,
//! generates a net electric polarization:
//!
//! ```text
//! P_pair = γ · (e_ij × (S_i × S_j))
//! ```
//!
//! where γ [C/m²] is the spin-lattice coupling constant.  For a uniform spiral with
//! propagation vector **q** and spin rotation axis **n̂**, the site-averaged result is:
//!
//! ```text
//! ⟨P⟩ = γ · (q̂ × n̂)   (cycloidal: q ⊥ n̂, non-zero)
//!     = 0                (helical:   q ∥ n̂, zero by symmetry)
//! ```
//!
//! ## Inverse Magnetoelectric Effect
//!
//! By time-reversal symmetry the same coupling that allows H → P also allows E → M:
//!
//! ```text
//! M_induced = α · E
//! ```
//!
//! This electric-field control of magnetism is the basis for proposed ME memory
//! and logic devices that switch magnetic order with voltage rather than current.
//!
//! ## Magnon Drag Contribution
//!
//! At finite temperature, thermally excited magnons (spin waves) carry both spin
//! angular momentum and, in spin-orbit-coupled systems, electric polarization.
//! The magnon-drag contribution to the polarization scales with the Bose-Einstein
//! occupation of the magnon gap mode:
//!
//! ```text
//! P_drag ≈ γ_drag · n_BE(ω_gap, T)
//! ```
//!
//! where n_BE(ω, T) = 1 / (exp(ℏω / k_B T) − 1) at zero chemical potential.
//!
//! ## References
//!
//! - H. Katsura, N. Nagaosa & A. V. Balatsky, "Spin Current and Magnetoelectric Effect
//!   in Noncollinear Magnets", Phys. Rev. Lett. **95**, 057205 (2005) — KNB mechanism
//! - M. Mostovoy, "Ferroelectricity in Spiral Magnets",
//!   Phys. Rev. Lett. **96**, 067601 (2006) — continuum spin-current model
//! - S.-W. Cheong & M. Mostovoy, "Multiferroics: a magnetic twist for ferroelectricity",
//!   Nat. Mater. **6**, 13 (2007) — comprehensive review
//! - T. Kimura et al., "Magnetic control of ferroelectric polarization",
//!   Nature **426**, 55 (2003) — TbMnO₃ experiments

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::constants::{HBAR, KB};
use crate::error::{self, Result};
use crate::vector3::Vector3;

// ============================================================================
// KnbMechanism
// ============================================================================

/// Katsura-Nagaosa-Balatsky (KNB) spin-current mechanism for magnetically-induced
/// ferroelectricity in spiral magnets.
///
/// The KNB model (PRL 95, 057205, 2005) provides a microscopic derivation of the
/// electric polarization induced by a non-collinear spin pair in a dielectric solid.
/// The coupling constant `γ` has units [C/m²] and encodes the strength of the
/// relativistic spin-orbit/lattice coupling; it is typically ~ 10⁻⁴ C/m² for
/// Type-II multiferroics such as TbMnO₃.
///
/// The `crystal_axis` is the normal to the spin-rotation plane of the cycloidal
/// spiral (the axis about which consecutive spins rotate); for TbMnO₃ this is
/// the crystallographic c-axis.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct KnbMechanism {
    /// Spin-lattice coupling constant γ_KNB [C/m²].
    ///
    /// Typical values: ~8 × 10⁻⁴ C/m² for TbMnO₃.
    coupling: f64,

    /// Crystal symmetry axis — the normal to the spin rotation plane (i.e., the
    /// axis about which consecutive spins precess in the spiral).
    ///
    /// For TbMnO₃: c-axis = [0, 0, 1] (orthorhombic setting with b-spiral).
    crystal_axis: Vector3<f64>,
}

impl KnbMechanism {
    // -------------------------------------------------------------------------
    // Constructors
    // -------------------------------------------------------------------------

    /// Create a new KNB mechanism with given coupling constant and crystal axis.
    ///
    /// # Arguments
    /// * `coupling` - Spin-lattice coupling γ [C/m²]; must be strictly positive.
    /// * `crystal_axis` - Normal to spin rotation plane; must be non-zero.
    ///
    /// # Errors
    /// Returns [`crate::error::Error::InvalidParameter`] if `coupling ≤ 0` or
    /// if `crystal_axis` has zero magnitude.
    pub fn new(coupling: f64, crystal_axis: Vector3<f64>) -> Result<Self> {
        if coupling <= 0.0 {
            return Err(error::invalid_param(
                "coupling",
                "KNB spin-lattice coupling constant must be positive",
            ));
        }
        if crystal_axis.magnitude() == 0.0 {
            return Err(error::invalid_param(
                "crystal_axis",
                "crystal axis must be a non-zero vector",
            ));
        }
        Ok(Self {
            coupling,
            crystal_axis,
        })
    }

    /// TbMnO₃ KNB parameters.
    ///
    /// TbMnO₃ is the archetypal Type-II multiferroic.  Below T_N = 41 K the
    /// Mn³⁺ spins form a sinusoidal wave; below T_FE = 27 K this locks into a
    /// cycloidal spiral propagating along the b-axis with spins rotating in the
    /// bc-plane.  The induced polarization points along the c-axis.
    ///
    /// Parameters:
    /// - γ ≈ 8 × 10⁻⁴ C/m² (estimated from P_s ≈ 800 μC/m² and S²≈6)
    /// - crystal_axis = ẑ (c-axis in the orthorhombic frame)
    ///
    /// Reference: T. Kimura et al., Nature 426, 55 (2003)
    pub fn tbmno3() -> Self {
        Self {
            coupling: 8.0e-4,                          // [C/m²]
            crystal_axis: Vector3::new(0.0, 0.0, 1.0), // c-axis
        }
    }

    // -------------------------------------------------------------------------
    // Core physics methods
    // -------------------------------------------------------------------------

    /// KNB pair polarization: P = γ × (e_ij × (S_i × S_j)).
    ///
    /// Computes the electric polarization induced by a single spin pair using the
    /// microscopic KNB formula.  The result is a vector [C/m²] whose direction
    /// is given by e_ij × (S_i × S_j) and whose magnitude is scaled by `coupling`.
    ///
    /// Returns the zero vector if:
    /// - S_i and S_j are parallel (S_i × S_j = 0)
    /// - e_ij is parallel to S_i × S_j (cross product is zero)
    /// - `e_ij` has zero length
    ///
    /// # Arguments
    /// * `s_i` - Spin at site i (normalised or in units of ℏ/2).
    /// * `s_j` - Spin at site j.
    /// * `e_ij` - Bond vector from i to j (does not need to be normalised).
    pub fn polarization_from_spiral(
        &self,
        s_i: &Vector3<f64>,
        s_j: &Vector3<f64>,
        e_ij: &Vector3<f64>,
    ) -> Vector3<f64> {
        let mag = e_ij.magnitude();
        if mag == 0.0 {
            return Vector3::zero();
        }
        let e_hat = Vector3::new(e_ij.x / mag, e_ij.y / mag, e_ij.z / mag);
        let spin_current = s_i.cross(s_j);
        let direction = e_hat.cross(&spin_current);
        Vector3::new(
            self.coupling * direction.x,
            self.coupling * direction.y,
            self.coupling * direction.z,
        )
    }

    /// Total KNB polarization from a spin chain (nearest-neighbour sum).
    ///
    /// Sums the pair contributions over all adjacent spin pairs in the chain:
    ///
    /// ```text
    /// P_total = γ Σ_{<ij>} e_ij × (S_i × S_j)
    /// ```
    ///
    /// The chain has `N` spins and `N-1` bonds.  `lattice_vectors[k]` is the
    /// vector from site k to site k+1.
    ///
    /// # Arguments
    /// * `spins` - Ordered spin vectors S_0 … S_{N-1}.
    /// * `lattice_vectors` - Bond vectors e_{01}, e_{12}, … e_{(N-2)(N-1)}.
    ///   Must have length `spins.len() - 1`.
    ///
    /// # Errors
    /// Returns an error if `spins` has fewer than 2 elements or if
    /// `lattice_vectors.len() != spins.len() - 1`.
    pub fn total_polarization_chain(
        &self,
        spins: &[Vector3<f64>],
        lattice_vectors: &[Vector3<f64>],
    ) -> Result<Vector3<f64>> {
        if spins.len() < 2 {
            return Err(error::invalid_param(
                "spins",
                "chain must have at least 2 spins",
            ));
        }
        let expected_bonds = spins.len() - 1;
        if lattice_vectors.len() != expected_bonds {
            return Err(error::invalid_param(
                "lattice_vectors",
                &format!(
                    "expected {} bond vectors for {} spins, got {}",
                    expected_bonds,
                    spins.len(),
                    lattice_vectors.len()
                ),
            ));
        }
        let mut total = Vector3::zero();
        for k in 0..expected_bonds {
            let p_pair =
                self.polarization_from_spiral(&spins[k], &spins[k + 1], &lattice_vectors[k]);
            total = Vector3::new(total.x + p_pair.x, total.y + p_pair.y, total.z + p_pair.z);
        }
        Ok(total)
    }

    /// Analytically compute the spiral-averaged KNB polarization.
    ///
    /// For a uniform spin spiral with propagation vector **q** and spin-rotation
    /// axis **n̂** the site-averaged polarization is:
    ///
    /// ```text
    /// ⟨P⟩ = γ · (q̂ × n̂)
    /// ```
    ///
    /// This vanishes identically for a helical spiral (**q** ∥ **n̂**) and is
    /// maximal for a cycloidal spiral (**q** ⊥ **n̂**).
    ///
    /// # Arguments
    /// * `spiral_q` - Spiral propagation vector **q** (need not be normalised).
    /// * `rotation_axis` - Axis about which consecutive spins rotate (**n̂**).
    ///
    /// # Returns
    /// Averaged polarization vector [C/m²]; zero vector for helical spirals.
    pub fn polarization_from_spin_spiral(
        &self,
        spiral_q: &Vector3<f64>,
        rotation_axis: &Vector3<f64>,
    ) -> Vector3<f64> {
        let q_mag = spiral_q.magnitude();
        if q_mag == 0.0 {
            return Vector3::zero();
        }
        let q_hat = Vector3::new(spiral_q.x / q_mag, spiral_q.y / q_mag, spiral_q.z / q_mag);
        let cross = q_hat.cross(rotation_axis);
        Vector3::new(
            self.coupling * cross.x,
            self.coupling * cross.y,
            self.coupling * cross.z,
        )
    }

    /// Return the KNB coupling constant γ [C/m²].
    pub fn coupling(&self) -> f64 {
        self.coupling
    }

    /// Return the crystal axis (normalised copy).
    pub fn crystal_axis(&self) -> Vector3<f64> {
        let n = self.crystal_axis.magnitude();
        if n > 0.0 {
            Vector3::new(
                self.crystal_axis.x / n,
                self.crystal_axis.y / n,
                self.crystal_axis.z / n,
            )
        } else {
            self.crystal_axis
        }
    }
}

// ============================================================================
// InverseMagnetoelectric
// ============================================================================

/// Inverse magnetoelectric effect — electric-field control of magnetization.
///
/// The inverse ME couples an applied electric field **E** [V/m] to an induced
/// magnetization **M** [A/m] through the isotropic (scalar) coupling constant α [s/m]:
///
/// ```text
/// M_induced = α · E
/// ```
///
/// The corresponding free-energy coupling term is:
///
/// ```text
/// u_ME = -α (E · H)
/// ```
///
/// which is negative (stabilising) when **E**, **H**, and the ME coupling are
/// aligned, and allows an external electric field to preferentially stabilise one
/// magnetic domain over another.
///
/// This isotropic single-constant model is appropriate when the full tensor α is
/// dominated by a single diagonal component, or as a first-order estimate.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct InverseMagnetoelectric {
    /// Applied electric field **E** [V/m].
    pub e_field: Vector3<f64>,

    /// Isotropic ME coupling constant α [s/m].
    ///
    /// Positive α means that **M** is induced parallel to **E**.
    pub alpha: f64,
}

impl InverseMagnetoelectric {
    // -------------------------------------------------------------------------
    // Constructor
    // -------------------------------------------------------------------------

    /// Create a new `InverseMagnetoelectric` instance.
    ///
    /// # Arguments
    /// * `e_field` - Applied electric field [V/m].
    /// * `alpha`   - ME coupling constant [s/m]; must be non-zero.
    ///
    /// # Errors
    /// Returns an error if `alpha == 0`.
    pub fn new(e_field: Vector3<f64>, alpha: f64) -> Result<Self> {
        if alpha == 0.0 {
            return Err(error::invalid_param(
                "alpha",
                "ME coupling constant must be non-zero",
            ));
        }
        Ok(Self { e_field, alpha })
    }

    // -------------------------------------------------------------------------
    // Response methods
    // -------------------------------------------------------------------------

    /// Compute the magnetization induced by the electric field via the inverse ME effect.
    ///
    /// ```text
    /// M = α · E   [A/m]
    /// ```
    pub fn induced_magnetization(&self) -> Vector3<f64> {
        Vector3::new(
            self.alpha * self.e_field.x,
            self.alpha * self.e_field.y,
            self.alpha * self.e_field.z,
        )
    }

    /// Free-energy density of the ME coupling [J/m³].
    ///
    /// The bimagnetic coupling term in the free energy is:
    ///
    /// ```text
    /// u = -α (E · H)
    /// ```
    ///
    /// Since we may not have the full **H** field, this method approximates it
    /// using the provided magnetization **M** (treating M ~ H for a weakly
    /// magnetic system):
    ///
    /// ```text
    /// u ≈ -α |E| |M| cos θ
    /// ```
    ///
    /// where θ is the angle between **E** and **M**.
    ///
    /// # Arguments
    /// * `magnetization` - Current magnetization **M** [A/m].
    pub fn energy_density(&self, magnetization: &Vector3<f64>) -> f64 {
        -self.alpha * self.e_field.dot(magnetization)
    }

    /// Effective magnetic field arising from the ME coupling (electric-field torque).
    ///
    /// The ME coupling contributes an effective field on the magnetization:
    ///
    /// ```text
    /// H_ME = α · E   [A/m]
    /// ```
    ///
    /// This effective field enters the LLG equation and drives the magnetization
    /// towards alignment with the electric field (for positive α).
    ///
    /// # Arguments
    /// * `_current_m` - Current magnetization state (reserved for future nonlinear extensions).
    pub fn electric_field_control_of_magnetism(&self, _current_m: &Vector3<f64>) -> Vector3<f64> {
        self.induced_magnetization()
    }
}

// ============================================================================
// Free Functions
// ============================================================================

/// Spin current for the KNB picture: J_spin = S_1 × S_2.
///
/// In the Katsura-Nagaosa-Balatsky picture the vector product S_1 × S_2 plays
/// the role of a spin current density flowing between the two sites.  This
/// cross product is the key quantity that, combined with the bond direction,
/// generates the electric polarization.
///
/// Returns a zero vector when S_1 and S_2 are collinear (no spin current for
/// ferromagnetic or antiferromagnetic collinear order).
///
/// # Arguments
/// * `s1` - Spin vector at site 1.
/// * `s2` - Spin vector at site 2.
pub fn spin_current_from_spiral(s1: &Vector3<f64>, s2: &Vector3<f64>) -> Vector3<f64> {
    s1.cross(s2)
}

/// Magnon-drag contribution to the electric polarization at finite temperature.
///
/// Thermally excited magnons propagating in a spin-orbit-coupled multiferroic
/// carry an electric polarization that is proportional to the Bose-Einstein
/// occupation of the magnon gap mode at angular frequency `magnon_gap` \[rad/s\]:
///
/// ```text
/// n_BE(ω_gap, T) = 1 / (exp(ℏ ω_gap / k_B T) − 1)
/// ```
///
/// The polarization contribution scales as:
///
/// ```text
/// P_drag ≈ coupling · n_BE(ω_gap, T)   [C/m²]
/// ```
///
/// At T = 0 there are no thermally excited magnons and P_drag → 0.
/// As T increases above the magnon gap energy (k_B T ≫ ℏ ω_gap), the
/// occupation grows and the polarization increases.
///
/// # Arguments
/// * `temperature` - Temperature \[K\]; must be positive for a physically
///   meaningful result.  Returns 0.0 for T = 0.
/// * `magnon_gap` - Gap angular frequency ω \[rad/s\]; must be non-negative.
/// * `coupling` - Magnon-drag coupling constant [C/m²].
///
/// # Returns
/// Estimated magnon-drag polarization [C/m²].  Returns 0.0 for T ≤ 0 or
/// ω ≤ 0 (gap-less case treated as no magnon contribution).
pub fn magnon_drag_contribution(temperature: f64, magnon_gap: f64, coupling: f64) -> f64 {
    if temperature <= 0.0 || magnon_gap <= 0.0 {
        return 0.0;
    }
    let exponent = HBAR * magnon_gap / (KB * temperature);
    // Guard against overflow in exp
    if exponent > 700.0 {
        // Exponentially suppressed: effectively zero occupation
        return 0.0;
    }
    let n_be = 1.0 / (exponent.exp() - 1.0);
    coupling * n_be
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    const TOL: f64 = 1.0e-10;

    // -----------------------------------------------------------------------
    // KNB: cycloidal spiral gives non-zero P perpendicular to q and n̂
    // -----------------------------------------------------------------------

    #[test]
    fn test_knb_cycloidal_spiral_nonzero() {
        let knb = KnbMechanism::tbmno3();
        // Cycloidal: q ∥ x̂, rotation_axis ∥ ẑ  → P ∝ x̂ × ẑ = -ŷ
        let q = Vector3::new(1.0, 0.0, 0.0);
        let n_hat = Vector3::new(0.0, 0.0, 1.0);
        let p = knb.polarization_from_spin_spiral(&q, &n_hat);
        // P should be non-zero and perpendicular to both q and n_hat
        assert!(p.magnitude() > TOL, "Cycloidal spiral must give non-zero P");
        // P ∝ q̂ × n̂ = x̂ × ẑ = -ŷ
        assert!(p.y < 0.0, "P_y should be negative for q=x̂, n̂=ẑ");
        // Perpendicular to q (x-axis): P·q = 0?  No — P·ẑ = 0 but P can have x-comp.
        // Actually q̂ × n̂ = [1,0,0] × [0,0,1] = [0·1-0·0, 0·0-1·1, 1·0-0·0] = [0,-1,0]
        assert!(p.x.abs() < TOL, "P_x should be 0");
        assert!(p.z.abs() < TOL, "P_z should be 0");
    }

    // -----------------------------------------------------------------------
    // KNB: helical spiral gives zero P
    // -----------------------------------------------------------------------

    #[test]
    fn test_knb_helical_spiral_zero() {
        let knb = KnbMechanism::tbmno3();
        // Helical: q ∥ rotation_axis (both along ẑ)
        let q = Vector3::new(0.0, 0.0, 1.0);
        let n_hat = Vector3::new(0.0, 0.0, 1.0);
        let p = knb.polarization_from_spin_spiral(&q, &n_hat);
        // q̂ × n̂ = ẑ × ẑ = 0
        assert!(p.magnitude() < TOL, "Helical spiral must give zero P");
    }

    // -----------------------------------------------------------------------
    // Chain: AFM (alternating) gives zero total polarization
    // -----------------------------------------------------------------------

    #[test]
    fn test_total_polarization_chain_afm_zero() {
        let knb = KnbMechanism::tbmno3();
        // AFM chain: all spins alternate +x/-x, bond along x
        // S_k = (-1)^k * x̂; e_ij = x̂ for all bonds
        // S_i × S_{i+1} = [1,0,0] × [-1,0,0] = 0 → all pair P = 0
        let spins: Vec<Vector3<f64>> = (0..6)
            .map(|k| {
                if k % 2 == 0 {
                    Vector3::new(1.0, 0.0, 0.0)
                } else {
                    Vector3::new(-1.0, 0.0, 0.0)
                }
            })
            .collect();
        let bonds: Vec<Vector3<f64>> = (0..5).map(|_| Vector3::new(1.0, 0.0, 0.0)).collect();
        let p = knb
            .total_polarization_chain(&spins, &bonds)
            .expect("valid chain");
        assert!(p.magnitude() < TOL, "AFM chain should give zero total P");
    }

    // -----------------------------------------------------------------------
    // Inverse ME: M ∥ E for positive alpha
    // -----------------------------------------------------------------------

    #[test]
    fn test_inverse_me_m_parallel_e() {
        let e = Vector3::new(0.0, 0.0, 1.0e6); // E along z
        let alpha = 4.13e-12_f64;
        let ime = InverseMagnetoelectric::new(e, alpha).expect("valid");
        let m = ime.induced_magnetization();
        // M should be along z
        assert!(
            m.z > 0.0,
            "M_z should be positive for positive alpha and E_z > 0"
        );
        assert!(m.x.abs() < TOL);
        assert!(m.y.abs() < TOL);
    }

    // -----------------------------------------------------------------------
    // Energy density: negative for aligned M and E with positive alpha
    // -----------------------------------------------------------------------

    #[test]
    fn test_energy_density_negative_when_aligned() {
        let e = Vector3::new(0.0, 0.0, 1.0e6);
        let alpha = 4.13e-12_f64;
        let ime = InverseMagnetoelectric::new(e, alpha).expect("valid");
        let m = Vector3::new(0.0, 0.0, 1.0e5); // M along z (same direction as E)
        let u = ime.energy_density(&m);
        // u = -alpha * E·M = -alpha * E_z * M_z < 0
        assert!(
            u < 0.0,
            "Energy should be negative for aligned E and M with α > 0; got {u}"
        );
    }

    // -----------------------------------------------------------------------
    // Spin current: S1 × S2 for known inputs
    // -----------------------------------------------------------------------

    #[test]
    fn test_spin_current_known_result() {
        let s1 = Vector3::new(1.0, 0.0, 0.0);
        let s2 = Vector3::new(0.0, 1.0, 0.0);
        // x̂ × ŷ = ẑ
        let j = spin_current_from_spiral(&s1, &s2);
        assert!((j.x).abs() < TOL);
        assert!((j.y).abs() < TOL);
        assert!((j.z - 1.0).abs() < TOL);
    }

    // -----------------------------------------------------------------------
    // Magnon drag: increases with temperature above gap
    // -----------------------------------------------------------------------

    #[test]
    fn test_magnon_drag_increases_with_temperature() {
        // Gap frequency: THz range, ℏω ≈ k_B * 10 K → ω = k_B * 10 / ℏ
        let omega_gap = KB * 10.0 / HBAR; // ω corresponding to 10 K gap
        let coupling = 1.0e-5;
        let p_low = magnon_drag_contribution(5.0, omega_gap, coupling);
        let p_mid = magnon_drag_contribution(50.0, omega_gap, coupling);
        let p_high = magnon_drag_contribution(300.0, omega_gap, coupling);
        assert!(
            p_low < p_mid && p_mid < p_high,
            "Magnon drag should increase with T: p_low={p_low}, p_mid={p_mid}, p_high={p_high}"
        );
    }

    // -----------------------------------------------------------------------
    // KNB preset: TbMnO3 crystal axis is c-axis (z)
    // -----------------------------------------------------------------------

    #[test]
    fn test_tbmno3_crystal_axis() {
        let knb = KnbMechanism::tbmno3();
        let ax = knb.crystal_axis();
        // normalised c-axis should be [0, 0, 1]
        assert!(ax.x.abs() < TOL);
        assert!(ax.y.abs() < TOL);
        assert!((ax.z - 1.0).abs() < TOL);
    }

    // -----------------------------------------------------------------------
    // polarization_from_spin_spiral: helical returns zero (pair-level check)
    // -----------------------------------------------------------------------

    #[test]
    fn test_polarization_from_spin_spiral_helical_zero() {
        let knb = KnbMechanism::new(1.0e-4, Vector3::new(0.0, 0.0, 1.0)).expect("valid KNB");
        // q ∥ rotation_axis → helical → P = 0
        let q = Vector3::new(1.0, 0.0, 0.0);
        let n = Vector3::new(1.0, 0.0, 0.0); // same direction as q
        let p = knb.polarization_from_spin_spiral(&q, &n);
        assert!(
            p.magnitude() < TOL,
            "Helical (q∥n̂) must give |P|=0; got {}",
            p.magnitude()
        );
    }

    // -----------------------------------------------------------------------
    // polarization_from_spin_spiral: cycloidal returns vector ∥ crystal_axis
    // -----------------------------------------------------------------------

    #[test]
    fn test_polarization_from_spin_spiral_cycloidal_direction() {
        // For TbMnO3: q = ŷ (b-axis), rotation_axis = ẑ (c-axis)
        // P ∝ ŷ × ẑ = x̂ (along a-axis)
        // But the crystal_axis is c = ẑ; P is NOT necessarily ∥ crystal_axis.
        // The KNB formula gives the direction; we just check it is non-zero and
        // along x̂ for this geometry.
        let knb = KnbMechanism::tbmno3(); // crystal_axis = ẑ
        let q = Vector3::new(0.0, 1.0, 0.0); // b-axis propagation
        let n = Vector3::new(0.0, 0.0, 1.0); // rotation about c-axis
                                             // ŷ × ẑ = x̂
        let p = knb.polarization_from_spin_spiral(&q, &n);
        assert!(p.magnitude() > TOL, "Cycloidal spiral must give non-zero P");
        assert!(p.x > 0.0, "P should point along +x̂ for q=ŷ, n̂=ẑ");
        assert!(p.y.abs() < TOL);
        assert!(p.z.abs() < TOL);
    }

    // -----------------------------------------------------------------------
    // Constructor: validation errors
    // -----------------------------------------------------------------------

    #[test]
    fn test_knb_new_rejects_zero_coupling() {
        let result = KnbMechanism::new(0.0, Vector3::new(0.0, 0.0, 1.0));
        assert!(result.is_err(), "Zero coupling should be rejected");
    }

    #[test]
    fn test_knb_new_rejects_zero_axis() {
        let result = KnbMechanism::new(1.0e-4, Vector3::zero());
        assert!(result.is_err(), "Zero crystal axis should be rejected");
    }

    #[test]
    fn test_inverse_me_new_rejects_zero_alpha() {
        let result = InverseMagnetoelectric::new(Vector3::new(1.0, 0.0, 0.0), 0.0);
        assert!(result.is_err(), "Zero alpha should be rejected");
    }

    // -----------------------------------------------------------------------
    // Magnon drag: zero temperature gives zero contribution
    // -----------------------------------------------------------------------

    #[test]
    fn test_magnon_drag_zero_temperature() {
        let omega_gap = KB * 10.0 / HBAR;
        let p = magnon_drag_contribution(0.0, omega_gap, 1.0e-5);
        assert!(p == 0.0, "T=0 should give zero magnon drag");
    }

    // -----------------------------------------------------------------------
    // Chain length validation
    // -----------------------------------------------------------------------

    #[test]
    fn test_chain_wrong_bond_count() {
        let knb = KnbMechanism::tbmno3();
        let spins = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
        ];
        // 3 spins need 2 bonds; provide 3 → error
        let bonds = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
        ];
        let result = knb.total_polarization_chain(&spins, &bonds);
        assert!(result.is_err(), "Wrong bond count should return error");
    }

    // Keep a PI reference to avoid unused import warning
    #[allow(dead_code)]
    const _HALF_PI: f64 = PI / 2.0;
}
