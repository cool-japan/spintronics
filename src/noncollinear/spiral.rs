//! Spin-spiral magnetic structures
//!
//! Implements cycloidal, helical, conical, and fan-type spin spirals with full
//! physical characterization including magnetization profiles, neutron structure
//! factors, chirality determination, and multiferroic polarization via the
//! Katsura-Nagaosa-Balatsky (KNB) mechanism.
//!
//! # Physics Background
//!
//! A spin spiral is described by a single-**q** ordering: at lattice position **r**
//! the spin direction is
//!
//! **S**(**r**) = A [ sin θ (cos(**q**·**r**) **e**₁ + sin(**q**·**r**) **e**₂) + cos θ **ê** ]
//!
//! where **ê** is the rotation axis (the axis about which spins precess), **e**₁ and
//! **e**₂ are unit vectors perpendicular to **ê** forming a right-handed frame with **ê**,
//! θ is the cone half-angle (π/2 for a planar spiral), and A is the spin amplitude.
//!
//! ## Chirality
//!
//! The chirality of a spiral is determined by the sense of rotation of **S** as a
//! function of position along **q**. It is positive (counterclockwise = left-handed
//! in the physics convention) when **S**(0) × **S**(δr) · **q̂** > 0.
//!
//! ## KNB Multiferroic Mechanism
//!
//! Katsura, Nagaosa, and Balatsky, *PRL* **95**, 057205 (2005) showed that for a
//! cycloidal spiral with propagation **q** and rotation axis **ê**:
//!
//! **P** ∝ **ê** × (**q** × **ê**) ∝ **q** − (**q**·**ê**)**ê**
//!
//! For a pure cycloidal spiral (**q** ⊥ **ê**) this simplifies to **P** ∝ **q** × **ê**.
//! Helical spirals (**q** ∥ **ê**) give zero polarization by symmetry.
//!
//! ## Landau-Lifshitz Classical Energy
//!
//! The classical energy density for a single-**q** spiral on a simple cubic lattice
//! (coordination z = 6, NN distance a) is:
//!
//! E/site = −2J [cos(q_x a) + cos(q_y a) + cos(q_z a)] S²
//!
//! which is minimized by **q** = 0 (FM, J > 0) or **q** = (π/a, π/a, π/a) (AFM, J < 0).
//!
//! # References
//!
//! - T. Kimura et al., "Magnetic control of ferroelectric polarization",
//!   *Nature* **426**, 55–58 (2003) — TbMnO₃ multiferroic
//! - H. Katsura, N. Nagaosa, A. V. Balatsky, "Spin current and magnetoelectric effect
//!   in noncollinear magnets", *PRL* **95**, 057205 (2005) — KNB mechanism
//! - M. Mostovoy, "Ferroelectricity in Spiral Magnets",
//!   *PRL* **96**, 067601 (2006) — cycloidal magnetism and ferroelectricity
//! - A. Yoshimori, "A New Type of Antiferromagnetic Structure in the Rutile Type Crystal",
//!   *J. Phys. Soc. Jpn.* **14**, 807–821 (1959)

use std::f64::consts::{PI, TAU};

use crate::error::{self, Result};
use crate::vector3::Vector3;

/// Classification of the spin-spiral geometry
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpiralType {
    /// Cycloidal spiral: the rotation plane **contains** the propagation vector **q**.
    /// Rotation axis is perpendicular to both **q** and the stacking direction.
    /// Prototype: TbMnO₃ (Kimura 2003). Breaks inversion → multiferroic.
    Cycloidal,

    /// Helical (or proper-screw) spiral: the rotation plane is **perpendicular** to **q**.
    /// Rotation axis ∥ **q**. Prototype: MnSi, FeGe. Preserves inversion in the plane
    /// perpendicular to **q** → no electric polarization by KNB mechanism.
    Helical,

    /// Conical spiral: spins trace a cone around a net ferromagnetic axis.
    /// Cone half-angle θ ∈ (0, π/2). At θ = 0 the structure becomes ferromagnetic,
    /// at θ = π/2 it reduces to a planar cycloidal or helical spiral.
    Conical,

    /// Fan structure: spins fan symmetrically about an easy axis, similar to a
    /// conical spiral but with the cone axis in the hard plane. Appears in applied
    /// fields on cycloidal/helical magnets as a metamagnetic phase.
    FanStructure,
}

/// Handedness of the spin rotation along the propagation direction
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpiralChirality {
    /// Spins rotate clockwise (right-handed) as viewed along **q**
    Clockwise,
    /// Spins rotate counter-clockwise (left-handed) as viewed along **q**
    CounterClockwise,
    /// Achiral: the spiral has a mirror symmetry (e.g. pure sinusoidal)
    Achiral,
}

/// A single-**q** spin spiral described by propagation vector, rotation axis, and cone angle
///
/// The magnetization profile is
///
/// **M**(**r**) = A [ sin θ (cos(**q**·**r**) **e**₁ + sin(**q**·**r**) **e**₂) + cos θ **ê** ]
///
/// where **ê** is the `rotation_axis`, θ is the `cone_angle_rad`, A is the `amplitude`,
/// and (**e**₁, **e**₂, **ê**) form a right-handed orthonormal frame.
#[derive(Debug, Clone)]
pub struct SpinSpiral {
    /// Propagation (ordering) wavevector **q** \[1/m\] in Cartesian coordinates
    pub q_vector: Vector3<f64>,

    /// Rotation axis **ê** (normalized): spins precess about this axis.
    /// For cycloidal spirals this is perpendicular to **q**; for helical spirals it equals **q̂**.
    pub rotation_axis: Vector3<f64>,

    /// Cone half-angle θ \[rad\] measured from `rotation_axis`.
    /// θ = π/2 → planar spiral; θ = 0 → ferromagnet; θ ∈ (0, π/2) → conical.
    pub cone_angle_rad: f64,

    /// Spin amplitude A (dimensionless spin quantum number S, e.g. 2.0 for Mn³⁺)
    pub amplitude: f64,

    /// Geometric classification of this spiral
    pub spiral_type: SpiralType,
}

impl SpinSpiral {
    /// Construct a general spin spiral with full parameter control.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `cone_angle_rad` is outside [0, π]
    /// - `amplitude` is not positive
    /// - `q_vector` has zero magnitude (use ferromagnet directly for **q** = 0)
    /// - `rotation_axis` has zero magnitude (cannot normalize)
    pub fn new(
        q_vector: Vector3<f64>,
        rotation_axis: Vector3<f64>,
        cone_angle_rad: f64,
        amplitude: f64,
    ) -> Result<Self> {
        if !(0.0..=PI).contains(&cone_angle_rad) {
            return Err(error::invalid_param(
                "cone_angle_rad",
                "cone angle must be in [0, π]",
            ));
        }
        if amplitude <= 0.0 {
            return Err(error::invalid_param(
                "amplitude",
                "amplitude must be positive",
            ));
        }
        if rotation_axis.magnitude_squared() < 1e-30 {
            return Err(error::invalid_param(
                "rotation_axis",
                "rotation axis must have nonzero magnitude",
            ));
        }

        // Determine spiral type from geometry
        let q_hat = if q_vector.magnitude_squared() > 1e-30 {
            q_vector.normalize()
        } else {
            Vector3::zero()
        };
        let axis_hat = rotation_axis.normalize();
        let q_axis_dot = q_hat.dot(&axis_hat).abs();

        let spiral_type = if (cone_angle_rad - PI / 2.0).abs() < 1e-9 {
            // Planar spiral
            if q_axis_dot < 0.1 {
                SpiralType::Cycloidal
            } else if q_axis_dot > 0.9 {
                SpiralType::Helical
            } else {
                SpiralType::Cycloidal
            }
        } else {
            SpiralType::Conical
        };

        Ok(Self {
            q_vector,
            rotation_axis: axis_hat,
            cone_angle_rad,
            amplitude,
            spiral_type,
        })
    }

    /// Construct a cycloidal spin spiral with rotation axis **q** × **normal**.
    ///
    /// In a cycloidal spiral the rotation plane contains **q**, so spins rotate
    /// within the plane spanned by **q** and (say) the stacking direction.
    /// The rotation axis is then perpendicular to both: **ê** = **q̂** × **n̂**.
    ///
    /// # Arguments
    ///
    /// - `q_vector` — propagation wavevector \[1/m\]
    /// - `normal_axis` — stacking or out-of-plane direction (e.g. `Vector3::unit_z()`)
    /// - `amplitude` — spin quantum number S
    pub fn cycloidal(q_vector: Vector3<f64>, normal_axis: Vector3<f64>, amplitude: f64) -> Self {
        let q_hat = q_vector.normalize();
        let n_hat = normal_axis.normalize();
        // Rotation axis = q × n (perpendicular to the spiral plane)
        let rot_axis = q_hat.cross(&n_hat);
        let rot_axis = if rot_axis.magnitude_squared() > 1e-30 {
            rot_axis.normalize()
        } else {
            // If q ∥ n, fall back to a perpendicular axis
            if q_hat.x.abs() < 0.9 {
                Vector3::new(1.0, 0.0, 0.0).cross(&q_hat).normalize()
            } else {
                Vector3::new(0.0, 1.0, 0.0).cross(&q_hat).normalize()
            }
        };

        Self {
            q_vector,
            rotation_axis: rot_axis,
            cone_angle_rad: PI / 2.0,
            amplitude,
            spiral_type: SpiralType::Cycloidal,
        }
    }

    /// Construct a helical (proper-screw) spin spiral.
    ///
    /// In a helical spiral the rotation axis is parallel to **q**, so spins wind
    /// like a helix. This geometry has no net electric polarization via the
    /// KNB mechanism because the spin current **e**_{ij} × (**S**_i × **S**_j) averages
    /// to zero over one period.
    ///
    /// # Arguments
    ///
    /// - `q_vector` — propagation wavevector \[1/m\]
    /// - `amplitude` — spin quantum number S
    pub fn helical(q_vector: Vector3<f64>, amplitude: f64) -> Self {
        let q_hat = if q_vector.magnitude_squared() > 1e-30 {
            q_vector.normalize()
        } else {
            Vector3::unit_z()
        };

        Self {
            q_vector,
            rotation_axis: q_hat,
            cone_angle_rad: PI / 2.0,
            amplitude,
            spiral_type: SpiralType::Helical,
        }
    }

    /// Construct a conical spin spiral.
    ///
    /// The spin traces a cone of half-angle `cone_angle_rad` about `rotation_axis`.
    /// When `cone_angle_rad` = π/2 this reduces to a planar spiral.
    ///
    /// # Errors
    ///
    /// Returns an error if `cone_angle_rad` is outside [0, π] or `amplitude` ≤ 0.
    pub fn conical(
        q_vector: Vector3<f64>,
        rotation_axis: Vector3<f64>,
        cone_angle_rad: f64,
        amplitude: f64,
    ) -> Result<Self> {
        if !(0.0..=PI).contains(&cone_angle_rad) {
            return Err(error::invalid_param(
                "cone_angle_rad",
                "cone angle must be in [0, π]",
            ));
        }
        if amplitude <= 0.0 {
            return Err(error::invalid_param(
                "amplitude",
                "amplitude must be positive",
            ));
        }
        if rotation_axis.magnitude_squared() < 1e-30 {
            return Err(error::invalid_param(
                "rotation_axis",
                "rotation axis must have nonzero magnitude",
            ));
        }

        Ok(Self {
            q_vector,
            rotation_axis: rotation_axis.normalize(),
            cone_angle_rad,
            amplitude,
            spiral_type: SpiralType::Conical,
        })
    }

    /// Preset for TbMnO₃ (terbium manganite) — the archetypal multiferroic spiral magnet.
    ///
    /// TbMnO₃ exhibits a cycloidal spin spiral with ordering wavevector
    /// **q** = (0, q_b, 0) where q_b ≈ 0.28 × (2π/b) in the *Pbnm* orthorhombic
    /// structure. Below T_c ≈ 27 K the cycloidal order induces a ferroelectric
    /// polarization P ∥ **c** via the KNB mechanism.
    ///
    /// Lattice parameters (approximate): a = 5.30 Å, b = 5.84 Å, c = 7.40 Å.
    /// The spiral is in the *bc* plane with the Mn³⁺ spins confined there.
    /// Mn³⁺ has spin S = 2.
    ///
    /// # Reference
    ///
    /// T. Kimura et al., *Nature* **426**, 55–58 (2003).
    pub fn terbium_manganese_oxide() -> Self {
        // b-axis lattice parameter [m]
        let b_lattice = 5.84e-10_f64;
        // Ordering wavevector along b*: q_b = 0.28 × (2π/b)
        let q_b = 0.28 * TAU / b_lattice;
        let q_vector = Vector3::new(0.0, q_b, 0.0);
        // Rotation axis for cycloidal spiral in the bc-plane is ê = x̂ (a-axis)
        // so the rotation plane spanned by (ŷ, ẑ) contains q ∥ ŷ → cycloidal
        let rotation_axis = Vector3::unit_x();
        // Mn³⁺ spin S = 2
        Self {
            q_vector,
            rotation_axis,
            cone_angle_rad: PI / 2.0,
            amplitude: 2.0,
            spiral_type: SpiralType::Cycloidal,
        }
    }

    // -------------------------------------------------------------------------
    // Derived frame vectors
    // -------------------------------------------------------------------------

    /// Build the orthonormal frame (e₁, e₂, ê) for the spiral.
    ///
    /// **ê** = `rotation_axis` (already normalized in constructor).
    /// **e**₁ and **e**₂ are chosen perpendicular to **ê** and perpendicular to each
    /// other, forming a right-handed triad (**e**₁, **e**₂, **ê**).
    ///
    /// Convention: **e**₁ is chosen as the component of **x̂** perpendicular to **ê**
    /// (falling back to **ŷ** when **ê** ∥ **x̂**).
    fn frame_vectors(&self) -> (Vector3<f64>, Vector3<f64>) {
        let e_hat = self.rotation_axis; // already unit vector
                                        // Choose a reference vector not parallel to ê
        let ref_v = if e_hat.x.abs() < 0.9 {
            Vector3::unit_x()
        } else {
            Vector3::unit_y()
        };
        // e1 = (ref - (ref · ê) ê).normalize()
        let proj = e_hat * ref_v.dot(&e_hat);
        let e1 = (ref_v - proj).normalize();
        // e2 = ê × e1
        let e2 = e_hat.cross(&e1);
        (e1, e2)
    }

    // -------------------------------------------------------------------------
    // Physical observables
    // -------------------------------------------------------------------------

    /// Magnetization vector at position **r** (meters).
    ///
    /// **M**(**r**) = A [ sin θ (cos(**q**·**r**) **e**₁ + sin(**q**·**r**) **e**₂) + cos θ **ê** ]
    ///
    /// For θ = π/2 (planar spiral): **M**(**r**) = A [cos(**q**·**r**) **e**₁ + sin(**q**·**r**) **e**₂].
    /// For θ = 0 (ferromagnet): **M**(**r**) = A **ê** (uniform, independent of **r**).
    pub fn magnetization_at(&self, r: &Vector3<f64>) -> Vector3<f64> {
        let phase = self.q_vector.dot(r); // **q** · **r** [rad]
        let (e1, e2) = self.frame_vectors();
        let sin_theta = self.cone_angle_rad.sin();
        let cos_theta = self.cone_angle_rad.cos();

        let planar_x = e1 * (phase.cos() * sin_theta);
        let planar_y = e2 * (phase.sin() * sin_theta);
        let axial = self.rotation_axis * cos_theta;

        (planar_x + planar_y + axial) * self.amplitude
    }

    /// Chirality of the spiral based on the sign of (**S**(0) × **S**(δ)) · **q̂**.
    ///
    /// Evaluates **M** at two points separated by a quarter wavelength along **q**
    /// and computes the triple product to determine the handedness.
    pub fn chirality(&self) -> SpiralChirality {
        let q_mag = self.q_vector.magnitude();
        if q_mag < 1e-30 {
            return SpiralChirality::Achiral;
        }
        let q_hat = self.q_vector.normalize();
        // Evaluate M at r = 0 and r = (π/2q) q̂ (quarter-period along q)
        let r0 = Vector3::zero();
        let r1 = q_hat * (PI / (2.0 * q_mag));
        let m0 = self.magnetization_at(&r0);
        let m1 = self.magnetization_at(&r1);

        // Triple product: (M(0) × M(λ/4)) · q̂
        let cross = m0.cross(&m1);
        let triple = cross.dot(&q_hat);

        if triple > 1e-10 {
            SpiralChirality::CounterClockwise
        } else if triple < -1e-10 {
            SpiralChirality::Clockwise
        } else {
            SpiralChirality::Achiral
        }
    }

    /// Electric polarization direction from the KNB mechanism (for cycloidal spirals).
    ///
    /// For a cycloidal spiral the polarization is
    ///
    /// **P** ∝ **ê** × (**q** × **ê**) = **q** − (**q**·**ê**)**ê**
    ///
    /// which points along the component of **q** perpendicular to the rotation axis **ê**.
    /// For a helical spiral **q** ∥ **ê**, so **P** = 0 and this method returns `None`.
    ///
    /// Returns `None` for helical spirals (zero polarization) or if **q** is zero.
    pub fn electric_polarization_direction(&self) -> Option<Vector3<f64>> {
        let q_mag = self.q_vector.magnitude();
        if q_mag < 1e-30 {
            return None;
        }
        // P ∝ q - (q·ê) ê  (component of q perpendicular to rotation axis)
        let q_dot_e = self.q_vector.dot(&self.rotation_axis);
        let p_vec = self.q_vector - self.rotation_axis * q_dot_e;
        let p_mag = p_vec.magnitude();
        if p_mag < 1e-30 {
            // q ∥ ê → helical, no polarization
            None
        } else {
            Some(p_vec.normalize())
        }
    }

    /// Wavelength of the spiral [m] = 2π / |**q**|.
    pub fn wavelength(&self) -> f64 {
        let q_mag = self.q_vector.magnitude();
        if q_mag > 1e-30 {
            TAU / q_mag
        } else {
            f64::INFINITY
        }
    }

    /// Pitch of the spiral [m] — synonym for `wavelength`.
    pub fn pitch(&self) -> f64 {
        self.wavelength()
    }

    /// Cone angle in degrees (for human-readable output).
    pub fn cone_angle_deg(&self) -> f64 {
        self.cone_angle_rad.to_degrees()
    }

    /// Returns `true` if the spiral is incommensurate with a lattice of constant `a` [m].
    ///
    /// Tests whether |**q**| × a / (2π) is close to a simple rational fraction p/q
    /// with denominator ≤ 20. An incommensurate wavevector is **not** of this form.
    pub fn is_incommensurate(&self, lattice_constant: f64) -> bool {
        let q_mag = self.q_vector.magnitude();
        if q_mag < 1e-30 || lattice_constant <= 0.0 {
            return false;
        }
        // Reduced wavevector: ξ = |q| a / (2π) — should be irrational for incommensurate
        let xi = q_mag * lattice_constant / TAU;
        // Check against simple fractions p/q with denominator up to 20
        for denom in 1_u32..=20 {
            for numer in 0_u32..=denom {
                let frac = numer as f64 / denom as f64;
                if (xi - frac).abs() < 1e-4 {
                    return false; // commensurate
                }
            }
        }
        true // no simple fraction found → incommensurate
    }

    /// Neutron scattering structure factor: returns spectral weights (S(**q**), S(−**q**)).
    ///
    /// For a single-**q** spiral the magnetic neutron structure factor has delta-function
    /// peaks at **k** = ±**q** (in reduced wavevector from the magnetic Bragg peak).
    /// The spectral weight of the peak at +**q** is proportional to the transverse spin
    /// fluctuation amplitude squared.
    ///
    /// For a planar spiral (θ = π/2) both peaks carry equal weight A²/4.
    /// For a conical spiral the weights are unequal due to the ferromagnetic component.
    ///
    /// Returns `(weight_plus_q, weight_minus_q)` in units of A² (the squared spin amplitude).
    pub fn spin_structure_factor(&self, k: &Vector3<f64>) -> (f64, f64) {
        let (e1, e2) = self.frame_vectors();
        let sin_theta = self.cone_angle_rad.sin();
        let a = self.amplitude;

        // In Fourier space the transverse spin fluctuations at ±q contribute
        // S(+q): weight from the rotating part (e1 + i e2) / 2 component
        // S(-q): weight from the (e1 - i e2) / 2 component
        // Each carries (A sin θ / 2)² to the transverse spin correlation

        // Evaluate how close k is to +q and -q
        let diff_plus = (*k - self.q_vector).magnitude();
        let diff_minus = (*k + self.q_vector).magnitude();
        let q_mag = self.q_vector.magnitude().max(1e-30);

        // Lorentzian broadening with width 0.01 |q| for numerical evaluation
        let width = 0.01 * q_mag;
        let lorentz_plus = width * width / (diff_plus * diff_plus + width * width);
        let lorentz_minus = width * width / (diff_minus * diff_minus + width * width);

        // Transverse weight: (A sin θ)² / 4 per satellite
        // For k perpendicular to q (neutron scattering geometry, full transverse component)
        let perp_e1 = e1.dot(k) / k.magnitude().max(1e-30);
        let perp_e2 = e2.dot(k) / k.magnitude().max(1e-30);
        let transverse_fraction = (1.0 - perp_e1 * perp_e1 - perp_e2 * perp_e2).clamp(0.0, 1.0);
        let base_weight = (a * sin_theta).powi(2) * 0.25 * (1.0 - transverse_fraction);

        let s_plus = base_weight * lorentz_plus;
        let s_minus = base_weight * lorentz_minus;

        (s_plus, s_minus)
    }

    /// Classical Heisenberg energy per spin on a simple cubic lattice (Landau-Lifshitz formula).
    ///
    /// For a simple cubic lattice with NN exchange J (`j_nn`) and lattice constant `a`,
    /// the six nearest-neighbor vectors are ±a **x̂**, ±a **ŷ**, ±a **ẑ**. The classical
    /// energy density for a single-**q** spiral is:
    ///
    /// E/spin = −2J [cos(q_x a) + cos(q_y a) + cos(q_z a)] S²
    ///
    /// - Ferromagnet (**q** = 0, J > 0): E = −6J S²
    /// - Antiferromagnet (**q** = (π/a, π/a, π/a), J < 0 i.e. AFM coupling): E = +6|J| S²
    ///
    /// # Arguments
    ///
    /// - `j_nn` — nearest-neighbor exchange constant [J]; positive = ferromagnetic
    /// - `a_lattice` — cubic lattice constant [m]
    pub fn landau_lifshitz_energy(&self, j_nn: f64, a_lattice: f64) -> f64 {
        let qx = self.q_vector.x;
        let qy = self.q_vector.y;
        let qz = self.q_vector.z;
        let a = a_lattice;
        let s = self.amplitude;

        // Sum over ±a x̂, ±a ŷ, ±a ẑ: Σ cos(q · δ)
        let cos_sum = (qx * a).cos()
            + (-qx * a).cos()
            + (qy * a).cos()
            + (-qy * a).cos()
            + (qz * a).cos()
            + (-qz * a).cos();

        // E/spin = -J Σ_δ cos(q · δ) S² = -J × cos_sum × S²
        -j_nn * cos_sum * s * s
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{PI, TAU};

    use super::*;

    const TOL: f64 = 1e-9;
    const TOL_MEDIUM: f64 = 1e-6;

    // Test 1: Cycloidal spiral gives a finite (nonzero) polarization direction
    #[test]
    fn test_cycloidal_gives_polarization() {
        let q = Vector3::new(1.0e7, 0.0, 0.0);
        let normal = Vector3::unit_z();
        let spiral = SpinSpiral::cycloidal(q, normal, 1.0);
        let p_dir = spiral.electric_polarization_direction();
        assert!(
            p_dir.is_some(),
            "cycloidal spiral must have finite polarization"
        );
        let p = p_dir.unwrap();
        // Polarization should be normalized
        assert!((p.magnitude() - 1.0).abs() < TOL_MEDIUM);
    }

    // Test 2: Helical spiral gives zero polarization direction (returns None)
    #[test]
    fn test_helical_gives_no_polarization() {
        let q = Vector3::new(0.0, 0.0, 1.0e7);
        let spiral = SpinSpiral::helical(q, 1.0);
        let p_dir = spiral.electric_polarization_direction();
        assert!(
            p_dir.is_none(),
            "helical spiral (q ∥ rotation_axis) must have zero polarization"
        );
    }

    // Test 3: magnetization_at returns unit norm for amplitude = 1, planar spiral
    #[test]
    fn test_magnetization_unit_norm() {
        let q = Vector3::new(2.0e7, 1.0e7, 0.0);
        let normal = Vector3::unit_z();
        let spiral = SpinSpiral::cycloidal(q, normal, 1.0);

        let positions = [
            Vector3::zero(),
            Vector3::new(1.0e-9, 0.0, 0.0),
            Vector3::new(0.0, 2.0e-9, 0.0),
            Vector3::new(5.0e-9, 3.0e-9, 1.0e-9),
        ];
        for r in &positions {
            let m = spiral.magnetization_at(r);
            let norm = m.magnitude();
            assert!(
                (norm - 1.0).abs() < 1e-10,
                "magnetization norm = {} at r = ({}, {}, {}), expected 1.0",
                norm,
                r.x,
                r.y,
                r.z
            );
        }
    }

    // Test 4: Conical spiral amplitude = 2 gives |M| = 2 everywhere
    #[test]
    fn test_conical_magnetization_norm() {
        let q = Vector3::new(1.0e7, 0.0, 0.0);
        let axis = Vector3::unit_z();
        let amplitude = 2.0;
        let spiral = SpinSpiral::conical(q, axis, PI / 4.0, amplitude).unwrap();

        let r = Vector3::new(3.14e-10, 0.0, 0.0);
        let m = spiral.magnetization_at(&r);
        assert!(
            (m.magnitude() - amplitude).abs() < 1e-9,
            "conical spiral |M| = {}, expected {}",
            m.magnitude(),
            amplitude
        );
    }

    // Test 5: Wavelength matches 2π/|q|
    #[test]
    fn test_wavelength() {
        let q_mag = 1.0e8_f64;
        let q = Vector3::new(q_mag, 0.0, 0.0);
        let spiral = SpinSpiral::helical(q, 1.0);
        let expected = TAU / q_mag;
        assert!(
            (spiral.wavelength() - expected).abs() < TOL,
            "wavelength = {}, expected {}",
            spiral.wavelength(),
            expected
        );
        assert!((spiral.pitch() - expected).abs() < TOL);
    }

    // Test 6: spin_structure_factor peaks at ±q
    #[test]
    fn test_structure_factor_peaks_at_q() {
        let q_mag = 1.0e8_f64;
        let q = Vector3::new(q_mag, 0.0, 0.0);
        let spiral = SpinSpiral::cycloidal(q, Vector3::unit_z(), 1.0);

        // Evaluate at k = q (should be peak) and at k far from q
        let (s_at_q, _) = spiral.spin_structure_factor(&q);
        let k_far = Vector3::new(0.5e8, 0.0, 0.0);
        let (s_far, _) = spiral.spin_structure_factor(&k_far);

        assert!(
            s_at_q > s_far,
            "structure factor should be larger at k=q ({}) than at k far ({})",
            s_at_q,
            s_far
        );
    }

    // Test 7: TbMnO3 preset has |q| ≈ 0.28 × 2π/b
    #[test]
    fn test_tbmno3_q_vector() {
        let spiral = SpinSpiral::terbium_manganese_oxide();
        let b_lattice = 5.84e-10_f64;
        let expected_q = 0.28 * TAU / b_lattice;
        let q_mag = spiral.q_vector.magnitude();
        let rel_err = (q_mag - expected_q).abs() / expected_q;
        assert!(
            rel_err < 1e-10,
            "|q| = {}, expected {} (0.28 × 2π/b)",
            q_mag,
            expected_q
        );
        // Must be cycloidal
        assert_eq!(spiral.spiral_type, SpiralType::Cycloidal);
        // Must produce a polarization
        assert!(spiral.electric_polarization_direction().is_some());
    }

    // Test 8: FM energy (q=0) on SC lattice = -6J S²
    #[test]
    fn test_fm_energy_sc_lattice() {
        let a = 3.0e-10_f64;
        let j_nn = 1.0e-21;
        let s = 1.0;
        // q = 0 → ferromagnet
        let q = Vector3::zero();
        let spiral = SpinSpiral {
            q_vector: q,
            rotation_axis: Vector3::unit_z(),
            cone_angle_rad: 0.0, // ferromagnet (θ=0)
            amplitude: s,
            spiral_type: SpiralType::Conical,
        };
        let energy = spiral.landau_lifshitz_energy(j_nn, a);
        let expected = -6.0 * j_nn * s * s;
        assert!(
            (energy - expected).abs() < TOL_MEDIUM * j_nn.abs(),
            "FM energy = {}, expected {}",
            energy,
            expected
        );
    }

    // Test 9: AFM energy (q = π/a along all axes) on SC lattice = +6J S²
    #[test]
    fn test_afm_energy_sc_lattice() {
        let a = 3.0e-10_f64;
        let j_nn = 1.0e-21;
        let s = 1.0;
        // q = (π/a, π/a, π/a) for SC AFM
        let q_afm = PI / a;
        let spiral = SpinSpiral {
            q_vector: Vector3::new(q_afm, q_afm, q_afm),
            rotation_axis: Vector3::unit_z(),
            cone_angle_rad: PI / 2.0,
            amplitude: s,
            spiral_type: SpiralType::Cycloidal,
        };
        let energy = spiral.landau_lifshitz_energy(j_nn, a);
        let expected = 6.0 * j_nn * s * s;
        assert!(
            (energy - expected).abs() < TOL_MEDIUM * j_nn.abs(),
            "AFM energy = {}, expected {}",
            energy,
            expected
        );
    }

    // Test 10: cone_angle_deg conversion
    #[test]
    fn test_cone_angle_degrees() {
        let spiral = SpinSpiral::helical(Vector3::new(1.0e7, 0.0, 0.0), 1.0);
        assert!((spiral.cone_angle_deg() - 90.0).abs() < TOL_MEDIUM);
    }

    // Test 11: new() rejects invalid cone angle
    #[test]
    fn test_new_invalid_cone_angle() {
        let q = Vector3::new(1.0e7, 0.0, 0.0);
        let axis = Vector3::unit_z();
        assert!(SpinSpiral::new(q, axis, -0.1, 1.0).is_err());
        assert!(SpinSpiral::new(q, axis, PI + 0.1, 1.0).is_err());
        assert!(SpinSpiral::new(q, axis, PI / 3.0, 1.0).is_ok());
    }

    // Test 12: incommensurate detection
    #[test]
    fn test_incommensurate() {
        let a = 3.0e-10_f64;
        // Commensurate: |q| = 2π/a → reduced wavevector = 1
        let q_comm = Vector3::new(TAU / a, 0.0, 0.0);
        let spiral_comm = SpinSpiral::helical(q_comm, 1.0);
        assert!(
            !spiral_comm.is_incommensurate(a),
            "q = 2π/a must be commensurate"
        );

        // Incommensurate: |q| = 0.28 × 2π/a
        let q_inc = Vector3::new(0.28 * TAU / a, 0.0, 0.0);
        let spiral_inc = SpinSpiral::helical(q_inc, 1.0);
        assert!(
            spiral_inc.is_incommensurate(a),
            "q = 0.28 × 2π/a must be incommensurate"
        );
    }
}
