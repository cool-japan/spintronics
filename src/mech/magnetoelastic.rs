//! Magnetoelastic energy and coupling
//!
//! This module implements the coupling between magnetic order and elastic strain
//! in magnetostrictive materials. Key phenomena include:
//!
//! ## Magnetoelastic Energy
//! The magnetoelastic energy density for cubic crystals:
//! E_me = B₁(ε_xx m_x² + ε_yy m_y² + ε_zz m_z²)
//!      + B₂(ε_xy m_x m_y + ε_yz m_y m_z + ε_xz m_x m_z)
//!
//! where B₁, B₂ are magnetoelastic coupling constants, ε_ij is the strain tensor,
//! and m_i are magnetization direction cosines.
//!
//! ## Magnetostrictive Strain (Inverse Effect)
//! Magnetization induces mechanical strain through magnetostriction constants
//! λ_100 and λ_111 (cubic) or λ_s (isotropic polycrystalline).
//!
//! ## Villari Effect (Direct Effect)
//! Applied stress changes the magnetic anisotropy and effective field:
//! H_σ = (3 λ_s σ) / (μ₀ M_s) · m̂_σ
//!
//! ## References
//! - O'Handley, "Modern Magnetic Materials: Principles and Applications"
//! - Chikazumi, "Physics of Ferromagnetism"

use crate::constants::MU_0;
use crate::error::{self, Result};
use crate::vector3::Vector3;

// =============================================================================
// Strain Tensor
// =============================================================================

/// Symmetric 3x3 strain tensor
///
/// The strain tensor ε_ij describes the deformation of a material.
/// Only the symmetric part is physically meaningful for elasticity:
/// ε_ij = (1/2)(∂u_i/∂x_j + ∂u_j/∂x_i)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrainTensor {
    /// The 3x3 symmetric components [[ε_xx, ε_xy, ε_xz], [ε_yx, ε_yy, ε_yz], [ε_zx, ε_zy, ε_zz]]
    pub components: [[f64; 3]; 3],
}

impl StrainTensor {
    /// Create a new strain tensor from components.
    ///
    /// The input is automatically symmetrized: ε_ij = (c_ij + c_ji) / 2
    pub fn new(components: [[f64; 3]; 3]) -> Self {
        let mut sym = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                sym[i][j] = 0.5 * (components[i][j] + components[j][i]);
            }
        }
        Self { components: sym }
    }

    /// Create a zero strain tensor
    pub fn zero() -> Self {
        Self {
            components: [[0.0; 3]; 3],
        }
    }

    /// Create a uniaxial strain along a given axis
    ///
    /// # Arguments
    /// * `strain` - The strain magnitude
    /// * `axis` - 0 for x, 1 for y, 2 for z
    ///
    /// # Errors
    /// Returns an error if axis is not 0, 1, or 2.
    pub fn uniaxial(strain: f64, axis: usize) -> Result<Self> {
        if axis > 2 {
            return Err(error::invalid_param("axis", "must be 0, 1, or 2"));
        }
        let mut components = [[0.0; 3]; 3];
        components[axis][axis] = strain;
        Ok(Self { components })
    }

    /// Create a biaxial in-plane strain (ε_xx = ε_yy = strain, ε_zz from Poisson)
    ///
    /// # Arguments
    /// * `strain` - The in-plane strain magnitude
    /// * `poisson_ratio` - Poisson's ratio ν
    pub fn biaxial_in_plane(strain: f64, poisson_ratio: f64) -> Self {
        let mut components = [[0.0; 3]; 3];
        components[0][0] = strain;
        components[1][1] = strain;
        // Out-of-plane strain from Poisson effect: ε_zz = -2ν/(1-ν) * ε_in-plane
        components[2][2] = -2.0 * poisson_ratio / (1.0 - poisson_ratio) * strain;
        Self { components }
    }

    /// Create a shear strain in the xy-plane
    pub fn shear_xy(strain: f64) -> Self {
        let mut components = [[0.0; 3]; 3];
        components[0][1] = strain;
        components[1][0] = strain;
        Self { components }
    }

    /// Get diagonal components (normal strains)
    pub fn diagonal(&self) -> (f64, f64, f64) {
        (
            self.components[0][0],
            self.components[1][1],
            self.components[2][2],
        )
    }

    /// Get the trace (volumetric strain)
    pub fn trace(&self) -> f64 {
        self.components[0][0] + self.components[1][1] + self.components[2][2]
    }

    /// Check if the tensor is symmetric (should always be true after construction)
    pub fn is_symmetric(&self, tol: f64) -> bool {
        for i in 0..3 {
            for j in (i + 1)..3 {
                if (self.components[i][j] - self.components[j][i]).abs() > tol {
                    return false;
                }
            }
        }
        true
    }

    /// Compute the elastic energy density: (1/2) C_ijkl ε_ij ε_kl
    /// For isotropic material: E = (λ/2)(tr ε)² + μ ε_ij ε_ij
    ///
    /// # Arguments
    /// * `lame_lambda` - First Lamé parameter λ \[Pa\]
    /// * `shear_modulus` - Shear modulus μ \[Pa\]
    pub fn elastic_energy_density_isotropic(&self, lame_lambda: f64, shear_modulus: f64) -> f64 {
        let tr = self.trace();
        let mut eps_sq = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                eps_sq += self.components[i][j] * self.components[i][j];
            }
        }
        0.5 * lame_lambda * tr * tr + shear_modulus * eps_sq
    }
}

impl Default for StrainTensor {
    fn default() -> Self {
        Self::zero()
    }
}

// =============================================================================
// Stress Tensor
// =============================================================================

/// Symmetric 3x3 stress tensor \[Pa\]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StressTensor {
    /// The 3x3 symmetric components \[Pa\]
    pub components: [[f64; 3]; 3],
}

impl StressTensor {
    /// Create a new stress tensor (automatically symmetrized)
    pub fn new(components: [[f64; 3]; 3]) -> Self {
        let mut sym = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                sym[i][j] = 0.5 * (components[i][j] + components[j][i]);
            }
        }
        Self { components: sym }
    }

    /// Create a zero stress tensor
    pub fn zero() -> Self {
        Self {
            components: [[0.0; 3]; 3],
        }
    }

    /// Create a uniaxial stress along a given axis
    pub fn uniaxial(stress: f64, axis: usize) -> Result<Self> {
        if axis > 2 {
            return Err(error::invalid_param("axis", "must be 0, 1, or 2"));
        }
        let mut components = [[0.0; 3]; 3];
        components[axis][axis] = stress;
        Ok(Self { components })
    }

    /// Create a hydrostatic (isotropic) stress
    pub fn hydrostatic(pressure: f64) -> Self {
        let mut components = [[0.0; 3]; 3];
        components[0][0] = -pressure;
        components[1][1] = -pressure;
        components[2][2] = -pressure;
        Self { components }
    }

    /// Get the trace (negative of 3× pressure for hydrostatic)
    pub fn trace(&self) -> f64 {
        self.components[0][0] + self.components[1][1] + self.components[2][2]
    }

    /// Get the von Mises equivalent stress
    pub fn von_mises(&self) -> f64 {
        let s = self.components;
        let s11 = s[0][0];
        let s22 = s[1][1];
        let s33 = s[2][2];
        let s12 = s[0][1];
        let s23 = s[1][2];
        let s13 = s[0][2];

        let term1 =
            (s11 - s22) * (s11 - s22) + (s22 - s33) * (s22 - s33) + (s33 - s11) * (s33 - s11);
        let term2 = 6.0 * (s12 * s12 + s23 * s23 + s13 * s13);

        (0.5 * (term1 + term2)).sqrt()
    }
}

impl Default for StressTensor {
    fn default() -> Self {
        Self::zero()
    }
}

// =============================================================================
// Magnetoelastic Material
// =============================================================================

/// Magnetoelastic material properties
///
/// Encapsulates the magnetoelastic coupling constants and mechanical properties
/// needed for magnetostriction and inverse magnetostrictive effect calculations.
#[derive(Debug, Clone, Copy)]
pub struct MagnetoelasticMaterial {
    /// Material name
    pub name: &'static str,
    /// B₁ magnetoelastic coupling constant \[J/m³\]
    ///
    /// Couples normal strains to magnetization components: B₁(ε_xx m_x² + ...)
    pub b1: f64,
    /// B₂ magnetoelastic coupling constant \[J/m³\]
    ///
    /// Couples shear strains to magnetization products: B₂(ε_xy m_x m_y + ...)
    pub b2: f64,
    /// Magnetostriction constant λ_100 (dimensionless, often quoted in ppm)
    ///
    /// Strain along \[100\] direction when magnetized along \[100\]
    pub lambda_100: f64,
    /// Magnetostriction constant λ_111 (dimensionless, often quoted in ppm)
    ///
    /// Strain along \[111\] direction when magnetized along \[111\]
    pub lambda_111: f64,
    /// Isotropic (polycrystalline) magnetostriction constant λ_s
    ///
    /// For polycrystalline samples: λ_s = (2/5)λ_100 + (3/5)λ_111
    pub lambda_s: f64,
    /// Young's modulus \[Pa\]
    pub youngs_modulus: f64,
    /// Saturation magnetization \[A/m\]
    pub ms: f64,
}

impl MagnetoelasticMaterial {
    /// Terfenol-D (Tb₀.₃Dy₀.₇Fe₂): giant magnetostrictive material
    ///
    /// The workhorse of magnetostrictive actuators with λ_s ~ 1000-2000 ppm.
    /// B₁ ≈ -9.4 MJ/m³, among the largest known coupling constants.
    pub fn terfenol_d() -> Self {
        Self {
            name: "Terfenol-D",
            b1: -9.38e6,
            b2: -6.16e6,
            lambda_100: 90.0e-6,
            lambda_111: 1640.0e-6,
            lambda_s: 1020.0e-6, // (2/5)*90 + (3/5)*1640 ≈ 1020 ppm
            youngs_modulus: 30.0e9,
            ms: 0.8e6,
        }
    }

    /// Galfenol (Fe₁₋ₓGaₓ, x ≈ 0.19): ductile magnetostrictive alloy
    ///
    /// Moderate magnetostriction (λ_s ~ 200-400 ppm) with good mechanical properties.
    pub fn galfenol() -> Self {
        Self {
            name: "Galfenol",
            b1: -4.0e6,
            b2: -3.0e6,
            lambda_100: 280.0e-6,
            lambda_111: 20.0e-6,
            lambda_s: 124.0e-6, // (2/5)*280 + (3/5)*20 = 124 ppm
            youngs_modulus: 60.0e9,
            ms: 1.3e6,
        }
    }

    /// Nickel: classic negative magnetostriction material
    ///
    /// λ_s ≈ -34 ppm (negative: contracts along magnetization direction)
    pub fn nickel() -> Self {
        Self {
            name: "Nickel",
            b1: 9.38e6,
            b2: 10.0e6,
            lambda_100: -46.0e-6,
            lambda_111: -24.0e-6,
            lambda_s: -32.8e-6, // (2/5)*(-46) + (3/5)*(-24) = -32.8 ppm
            youngs_modulus: 200.0e9,
            ms: 0.49e6,
        }
    }

    /// CoFeB: widely used in spintronics, small tunable magnetostriction
    ///
    /// λ_s ≈ 20-30 ppm, important for MRAM and spin-orbit torque devices
    pub fn cofeb() -> Self {
        Self {
            name: "CoFeB",
            b1: -3.5e6,
            b2: -2.0e6,
            lambda_100: 25.0e-6,
            lambda_111: 25.0e-6,
            lambda_s: 25.0e-6,
            youngs_modulus: 160.0e9,
            ms: 1.0e6,
        }
    }

    /// Iron (bcc Fe): moderate anisotropic magnetostriction
    ///
    /// λ_100 ≈ 21 ppm, λ_111 ≈ -21 ppm (nearly cancel → small λ_s)
    pub fn iron() -> Self {
        Self {
            name: "Iron",
            b1: -3.43e6,
            b2: 7.83e6,
            lambda_100: 21.0e-6,
            lambda_111: -21.0e-6,
            lambda_s: -4.2e-6, // (2/5)*21 + (3/5)*(-21) = -4.2 ppm
            youngs_modulus: 211.0e9,
            ms: 1.71e6,
        }
    }

    /// Compute polycrystalline λ_s from single-crystal constants
    ///
    /// λ_s = (2/5) λ_100 + (3/5) λ_111
    pub fn compute_lambda_s(&self) -> f64 {
        (2.0 / 5.0) * self.lambda_100 + (3.0 / 5.0) * self.lambda_111
    }

    /// Compute the Poisson's ratio assuming ν ≈ 0.3 (typical for metals)
    ///
    /// This is an approximation; for precise work, use measured values.
    pub fn typical_poisson_ratio(&self) -> f64 {
        0.3
    }

    /// Compute the shear modulus from Young's modulus (assuming ν = 0.3)
    ///
    /// G = E / (2(1 + ν))
    pub fn approximate_shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.typical_poisson_ratio()))
    }
}

// =============================================================================
// Magnetoelastic Energy
// =============================================================================

/// Compute the cubic magnetoelastic energy density \[J/m³\]
///
/// E_me = B₁(ε_xx m_x² + ε_yy m_y² + ε_zz m_z²)
///      + B₂(ε_xy m_x m_y + ε_yz m_y m_z + ε_xz m_x m_z)
///
/// # Arguments
/// * `material` - Magnetoelastic material properties
/// * `strain` - Strain tensor
/// * `magnetization_dir` - Unit vector of magnetization direction
///
/// # Errors
/// Returns an error if magnetization_dir has zero magnitude.
pub fn magnetoelastic_energy_density_cubic(
    material: &MagnetoelasticMaterial,
    strain: &StrainTensor,
    magnetization_dir: &Vector3<f64>,
) -> Result<f64> {
    let m_mag = magnetization_dir.magnitude();
    if m_mag < 1e-30 {
        return Err(error::invalid_param(
            "magnetization_dir",
            "must have non-zero magnitude",
        ));
    }

    // Normalize
    let m = *magnetization_dir * (1.0 / m_mag);
    let e = &strain.components;

    // B₁ term: normal strains coupled to direction cosines squared
    let b1_term = material.b1 * (e[0][0] * m.x * m.x + e[1][1] * m.y * m.y + e[2][2] * m.z * m.z);

    // B₂ term: shear strains coupled to products of direction cosines
    let b2_term = material.b2 * (e[0][1] * m.x * m.y + e[1][2] * m.y * m.z + e[0][2] * m.x * m.z);

    Ok(b1_term + b2_term)
}

/// Compute the total magnetoelastic energy for a given volume \[J\]
///
/// E = E_me_density × volume
pub fn magnetoelastic_energy(
    material: &MagnetoelasticMaterial,
    strain: &StrainTensor,
    magnetization_dir: &Vector3<f64>,
    volume: f64,
) -> Result<f64> {
    let density = magnetoelastic_energy_density_cubic(material, strain, magnetization_dir)?;
    Ok(density * volume)
}

// =============================================================================
// Magnetostrictive Strain (Inverse Effect)
// =============================================================================

/// Compute isotropic magnetostrictive strain
///
/// For a polycrystalline material with isotropic magnetostriction:
/// ε_me = (3/2) λ_s (cos²θ - 1/3)
///
/// where θ is the angle between the magnetization and the strain measurement direction.
///
/// # Arguments
/// * `lambda_s` - Isotropic magnetostriction constant
/// * `cos_theta` - Cosine of angle between magnetization and measurement direction
pub fn isotropic_magnetostrictive_strain(lambda_s: f64, cos_theta: f64) -> f64 {
    1.5 * lambda_s * (cos_theta * cos_theta - 1.0 / 3.0)
}

/// Compute the full magnetostrictive strain tensor for cubic crystals
///
/// For a cubic material with magnetization direction (m_x, m_y, m_z):
///
/// ε_xx = (3/2) λ_100 (m_x² - 1/3)
/// ε_yy = (3/2) λ_100 (m_y² - 1/3)
/// ε_zz = (3/2) λ_100 (m_z² - 1/3)
/// ε_xy = (3/2) λ_111 m_x m_y
/// ε_yz = (3/2) λ_111 m_y m_z
/// ε_xz = (3/2) λ_111 m_x m_z
///
/// # Arguments
/// * `material` - Magnetoelastic material
/// * `magnetization_dir` - Unit vector of magnetization direction
///
/// # Errors
/// Returns an error if magnetization_dir has zero magnitude.
pub fn magnetostrictive_strain_tensor(
    material: &MagnetoelasticMaterial,
    magnetization_dir: &Vector3<f64>,
) -> Result<StrainTensor> {
    let m_mag = magnetization_dir.magnitude();
    if m_mag < 1e-30 {
        return Err(error::invalid_param(
            "magnetization_dir",
            "must have non-zero magnitude",
        ));
    }

    let m = *magnetization_dir * (1.0 / m_mag);

    let l100 = material.lambda_100;
    let l111 = material.lambda_111;

    let mut c = [[0.0; 3]; 3];

    // Normal strains from λ_100
    c[0][0] = 1.5 * l100 * (m.x * m.x - 1.0 / 3.0);
    c[1][1] = 1.5 * l100 * (m.y * m.y - 1.0 / 3.0);
    c[2][2] = 1.5 * l100 * (m.z * m.z - 1.0 / 3.0);

    // Shear strains from λ_111
    c[0][1] = 1.5 * l111 * m.x * m.y;
    c[1][0] = c[0][1];
    c[0][2] = 1.5 * l111 * m.x * m.z;
    c[2][0] = c[0][2];
    c[1][2] = 1.5 * l111 * m.y * m.z;
    c[2][1] = c[1][2];

    Ok(StrainTensor { components: c })
}

/// Compute the saturation magnetostrictive strain along the magnetization direction
///
/// For isotropic case: ε_max = (3/2) λ_s × (2/3) = λ_s when θ=0
pub fn saturation_strain(lambda_s: f64) -> f64 {
    // When cos²θ = 1 (θ=0): ε = (3/2)λ_s(1 - 1/3) = λ_s
    1.5 * lambda_s * (1.0 - 1.0 / 3.0)
}

// =============================================================================
// Villari Effect (Direct Magnetoelastic Effect)
// =============================================================================

/// Compute the effective field from applied stress (Villari effect) \[A/m\]
///
/// The Villari effect describes how mechanical stress changes the magnetic
/// anisotropy, creating an effective magnetic field:
///
/// H_σ = (3 λ_s σ) / (μ₀ M_s)
///
/// The field direction depends on the sign of λ_s and σ:
/// - Positive λ_s, tensile σ: field along stress direction
/// - Negative λ_s, tensile σ: field perpendicular to stress
///
/// # Arguments
/// * `material` - Magnetoelastic material
/// * `stress` - Applied stress magnitude \[Pa\]
///
/// # Returns
/// The magnitude of the effective Villari field \[A/m\]
pub fn villari_effective_field_magnitude(
    material: &MagnetoelasticMaterial,
    stress: f64,
) -> Result<f64> {
    if material.ms.abs() < 1e-30 {
        return Err(error::invalid_param(
            "ms",
            "saturation magnetization must be non-zero",
        ));
    }
    let h_eff = (3.0 * material.lambda_s * stress) / (MU_0 * material.ms);
    Ok(h_eff)
}

/// Compute the Villari effective field vector \[A/m\]
///
/// H_σ = (3 λ_s / (μ₀ M_s)) × σ_ij m_j (tensor form)
///
/// For uniaxial stress σ along direction n̂, and magnetization along m̂:
/// The effective field has component along the stress direction weighted by
/// the projection of magnetization.
///
/// # Arguments
/// * `material` - Magnetoelastic material
/// * `stress_tensor` - Applied stress tensor \[Pa\]
/// * `magnetization_dir` - Unit vector of magnetization direction
///
/// # Errors
/// Returns an error if magnetization or M_s is zero.
pub fn villari_effective_field_vector(
    material: &MagnetoelasticMaterial,
    stress_tensor: &StressTensor,
    magnetization_dir: &Vector3<f64>,
) -> Result<Vector3<f64>> {
    let m_mag = magnetization_dir.magnitude();
    if m_mag < 1e-30 {
        return Err(error::invalid_param(
            "magnetization_dir",
            "must have non-zero magnitude",
        ));
    }
    if material.ms.abs() < 1e-30 {
        return Err(error::invalid_param(
            "ms",
            "saturation magnetization must be non-zero",
        ));
    }

    let m = *magnetization_dir * (1.0 / m_mag);
    let s = &stress_tensor.components;
    let prefactor = 3.0 * material.lambda_s / (MU_0 * material.ms);

    // H_i = prefactor * Σ_j σ_ij m_j
    let hx = prefactor * (s[0][0] * m.x + s[0][1] * m.y + s[0][2] * m.z);
    let hy = prefactor * (s[1][0] * m.x + s[1][1] * m.y + s[1][2] * m.z);
    let hz = prefactor * (s[2][0] * m.x + s[2][1] * m.y + s[2][2] * m.z);

    Ok(Vector3::new(hx, hy, hz))
}

/// Compute the stress-induced anisotropy energy density \[J/m³\]
///
/// K_σ = -(3/2) λ_s σ
///
/// For tensile stress (σ > 0) and positive λ_s, the easy axis aligns with stress.
/// For compressive stress (σ < 0) and positive λ_s, easy axis perpendicular to stress.
///
/// # Arguments
/// * `lambda_s` - Isotropic magnetostriction constant
/// * `stress` - Uniaxial stress \[Pa\]
pub fn stress_induced_anisotropy(lambda_s: f64, stress: f64) -> f64 {
    -1.5 * lambda_s * stress
}

/// Compute the magnetoelastic coupling energy between two magnetic layers
/// mediated by an elastic spacer.
///
/// For two magnetostrictive layers coupled through strain transfer:
/// E_coupling = -A_me × (m̂₁ · m̂₂) per unit area
///
/// where A_me depends on elastic properties and geometry.
///
/// # Arguments
/// * `material1` - First magnetic layer
/// * `material2` - Second magnetic layer
/// * `thickness1` - Thickness of first layer \[m\]
/// * `thickness2` - Thickness of second layer \[m\]
/// * `spacer_compliance` - Compliance of spacer (strain per unit stress) [1/Pa]
///
/// # Returns
/// The magnetoelastic interlayer coupling constant A_me \[J/m²\]
pub fn interlayer_magnetoelastic_coupling(
    material1: &MagnetoelasticMaterial,
    material2: &MagnetoelasticMaterial,
    thickness1: f64,
    thickness2: f64,
    spacer_compliance: f64,
) -> f64 {
    // The coupling arises from strain generated by one layer being transferred
    // to the other. The effective coupling constant is:
    // A_me ∝ (3/2)² λ_s1 λ_s2 Y_1 Y_2 t_1 t_2 / (spacer compliance)
    let strain_factor1 = 1.5 * material1.lambda_s * material1.youngs_modulus * thickness1;
    let strain_factor2 = 1.5 * material2.lambda_s * material2.youngs_modulus * thickness2;

    // Coupling mediated by elastic transfer through spacer
    strain_factor1 * strain_factor2 * spacer_compliance
}

// =============================================================================
// Strain from Piezoelectric Substrate
// =============================================================================

/// Compute the strain generated by a piezoelectric substrate
///
/// For a piezoelectric material under applied electric field E:
/// ε = d × E (for direct piezoelectric coefficients d_ij)
///
/// Common piezoelectric substrates:
/// - PZT: d₃₁ ≈ -170 pC/N, d₃₃ ≈ 400 pC/N
/// - PMN-PT: d₃₁ ≈ -1000 pC/N, d₃₃ ≈ 2000 pC/N
///
/// # Arguments
/// * `d31` - Piezoelectric coefficient d₃₁ \[m/V\]
/// * `d33` - Piezoelectric coefficient d₃₃ \[m/V\]
/// * `electric_field` - Applied electric field [V/m]
///
/// # Returns
/// StrainTensor with in-plane (d₃₁ mediated) and out-of-plane (d₃₃) strains
pub fn piezoelectric_strain(d31: f64, d33: f64, electric_field: f64) -> StrainTensor {
    let eps_in_plane = d31 * electric_field;
    let eps_out_of_plane = d33 * electric_field;

    let mut components = [[0.0; 3]; 3];
    components[0][0] = eps_in_plane;
    components[1][1] = eps_in_plane;
    components[2][2] = eps_out_of_plane;

    StrainTensor { components }
}

/// Common piezoelectric substrate parameters
#[derive(Debug, Clone, Copy)]
pub struct PiezoelectricSubstrate {
    /// Material name
    pub name: &'static str,
    /// Piezoelectric coefficient d₃₁ \[m/V\]
    pub d31: f64,
    /// Piezoelectric coefficient d₃₃ \[m/V\]
    pub d33: f64,
    /// Maximum safe electric field [V/m]
    pub max_field: f64,
}

impl PiezoelectricSubstrate {
    /// PZT (Lead Zirconate Titanate)
    pub fn pzt() -> Self {
        Self {
            name: "PZT",
            d31: -170.0e-12,
            d33: 400.0e-12,
            max_field: 2.0e6,
        }
    }

    /// PMN-PT (Lead Magnesium Niobate - Lead Titanate)
    ///
    /// Single crystal with very large piezoelectric coefficients
    pub fn pmn_pt() -> Self {
        Self {
            name: "PMN-PT",
            d31: -1000.0e-12,
            d33: 2000.0e-12,
            max_field: 1.0e6,
        }
    }

    /// BaTiO₃ (Barium Titanate)
    pub fn batio3() -> Self {
        Self {
            name: "BaTiO3",
            d31: -78.0e-12,
            d33: 190.0e-12,
            max_field: 3.0e6,
        }
    }

    /// Compute the strain tensor for a given applied voltage and thickness
    ///
    /// E_field = V / t_substrate
    pub fn strain_from_voltage(
        &self,
        voltage: f64,
        substrate_thickness: f64,
    ) -> Result<StrainTensor> {
        if substrate_thickness <= 0.0 {
            return Err(error::invalid_param(
                "substrate_thickness",
                "must be positive",
            ));
        }
        let e_field = voltage / substrate_thickness;
        Ok(piezoelectric_strain(self.d31, self.d33, e_field))
    }

    /// Maximum achievable in-plane strain |d₃₁ × E_max|
    pub fn max_in_plane_strain(&self) -> f64 {
        (self.d31 * self.max_field).abs()
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nickel_magnetostriction_negative() {
        let ni = MagnetoelasticMaterial::nickel();
        // Nickel has negative magnetostriction (contracts along field)
        assert!(
            ni.lambda_s < 0.0,
            "Nickel should have negative magnetostriction, got {}",
            ni.lambda_s
        );
        // Literature value: λ_s ≈ -34 ppm
        assert!(
            (ni.lambda_s - (-34.0e-6)).abs() < 5.0e-6,
            "Nickel λ_s should be close to -34 ppm, got {} ppm",
            ni.lambda_s * 1e6
        );
    }

    #[test]
    fn test_iron_magnetostriction_anisotropic() {
        let fe = MagnetoelasticMaterial::iron();
        // Iron: λ_100 ≈ 21 ppm, λ_111 ≈ -21 ppm
        assert!(
            (fe.lambda_100 - 21.0e-6).abs() < 2.0e-6,
            "Iron λ_100 should be ~21 ppm, got {} ppm",
            fe.lambda_100 * 1e6
        );
        assert!(
            (fe.lambda_111 - (-21.0e-6)).abs() < 2.0e-6,
            "Iron λ_111 should be ~-21 ppm, got {} ppm",
            fe.lambda_111 * 1e6
        );
        // λ_100 and λ_111 nearly cancel
        assert!(
            fe.lambda_s.abs() < 10.0e-6,
            "Iron λ_s should be small, got {} ppm",
            fe.lambda_s * 1e6
        );
    }

    #[test]
    fn test_terfenol_d_largest_magnetostriction() {
        let terfenol = MagnetoelasticMaterial::terfenol_d();
        let galfenol = MagnetoelasticMaterial::galfenol();
        let nickel = MagnetoelasticMaterial::nickel();
        let cofeb = MagnetoelasticMaterial::cofeb();
        let iron = MagnetoelasticMaterial::iron();

        // Terfenol-D should have the largest |λ_s| among our materials
        assert!(
            terfenol.lambda_s.abs() > galfenol.lambda_s.abs(),
            "Terfenol-D |λ_s| should exceed Galfenol"
        );
        assert!(
            terfenol.lambda_s.abs() > nickel.lambda_s.abs(),
            "Terfenol-D |λ_s| should exceed Nickel"
        );
        assert!(
            terfenol.lambda_s.abs() > cofeb.lambda_s.abs(),
            "Terfenol-D |λ_s| should exceed CoFeB"
        );
        assert!(
            terfenol.lambda_s.abs() > iron.lambda_s.abs(),
            "Terfenol-D |λ_s| should exceed Iron"
        );
    }

    #[test]
    fn test_strain_tensor_symmetry() {
        // Create an asymmetric input; constructor should symmetrize
        let asymmetric = [[0.001, 0.002, 0.0], [0.0, 0.003, 0.004], [0.0, 0.0, 0.005]];
        let tensor = StrainTensor::new(asymmetric);
        assert!(
            tensor.is_symmetric(1e-15),
            "StrainTensor should always be symmetric after construction"
        );
        // Check specific symmetrized values
        let expected_xy = 0.5 * (0.002 + 0.0);
        assert!(
            (tensor.components[0][1] - expected_xy).abs() < 1e-15,
            "ε_xy should be symmetrized"
        );
        assert!(
            (tensor.components[1][0] - expected_xy).abs() < 1e-15,
            "ε_yx should equal ε_xy"
        );
    }

    #[test]
    fn test_magnetoelastic_energy_zero_strain() {
        let ni = MagnetoelasticMaterial::nickel();
        let zero_strain = StrainTensor::zero();
        let m_dir = Vector3::new(1.0, 0.0, 0.0);

        let energy = magnetoelastic_energy_density_cubic(&ni, &zero_strain, &m_dir)
            .expect("should compute energy for zero strain");
        assert!(
            energy.abs() < 1e-20,
            "Magnetoelastic energy should be zero for zero strain, got {}",
            energy
        );
    }

    #[test]
    fn test_villari_effective_field_direction() {
        // For positive λ_s and tensile stress along x,
        // the Villari field should be along x when m is along x
        let cofeb = MagnetoelasticMaterial::cofeb();
        assert!(cofeb.lambda_s > 0.0, "CoFeB should have positive λ_s");

        let tensile_x = StressTensor::uniaxial(100.0e6, 0).expect("should create uniaxial stress");
        let m_along_x = Vector3::new(1.0, 0.0, 0.0);

        let h_vec = villari_effective_field_vector(&cofeb, &tensile_x, &m_along_x)
            .expect("should compute Villari field");

        // For positive λ_s, tensile stress, m along x: field should point along +x
        assert!(
            h_vec.x > 0.0,
            "Villari field x-component should be positive for positive λ_s and tensile stress along x, got {}",
            h_vec.x
        );
        // y and z components should be zero (uniaxial along x, m along x)
        assert!(
            h_vec.y.abs() < 1e-10,
            "Villari field y-component should be ~0, got {}",
            h_vec.y
        );
        assert!(
            h_vec.z.abs() < 1e-10,
            "Villari field z-component should be ~0, got {}",
            h_vec.z
        );
    }

    #[test]
    fn test_magnetostrictive_strain_tensor_cubic() {
        let fe = MagnetoelasticMaterial::iron();

        // Magnetization along [100]: should get ε_xx = λ_100, no shear
        let m_100 = Vector3::new(1.0, 0.0, 0.0);
        let strain =
            magnetostrictive_strain_tensor(&fe, &m_100).expect("should compute strain tensor");

        // ε_xx = (3/2) λ_100 (1 - 1/3) = λ_100
        let expected_exx = fe.lambda_100;
        assert!(
            (strain.components[0][0] - expected_exx).abs() < 1e-10,
            "ε_xx for m along [100] should be λ_100 = {}, got {}",
            expected_exx,
            strain.components[0][0]
        );

        // ε_yy = (3/2) λ_100 (0 - 1/3) = -λ_100/2
        let expected_eyy = -fe.lambda_100 / 2.0;
        assert!(
            (strain.components[1][1] - expected_eyy).abs() < 1e-10,
            "ε_yy for m along [100] should be -λ_100/2 = {}, got {}",
            expected_eyy,
            strain.components[1][1]
        );

        // No shear for magnetization along principal axis
        assert!(
            strain.components[0][1].abs() < 1e-20,
            "Shear ε_xy should be zero for m along [100]"
        );
    }

    #[test]
    fn test_piezoelectric_strain_pmn_pt() {
        let pmn_pt = PiezoelectricSubstrate::pmn_pt();

        // Apply 1 MV/m field
        let e_field = 1.0e6;
        let strain = piezoelectric_strain(pmn_pt.d31, pmn_pt.d33, e_field);

        // In-plane strain: d31 * E = -1000e-12 * 1e6 = -1e-3 (compressive)
        let expected_in_plane = pmn_pt.d31 * e_field;
        assert!(
            (strain.components[0][0] - expected_in_plane).abs() < 1e-10,
            "In-plane strain should be d31*E"
        );

        // PMN-PT should produce larger strains than PZT
        let pzt = PiezoelectricSubstrate::pzt();
        assert!(
            pmn_pt.d31.abs() > pzt.d31.abs(),
            "PMN-PT d31 should be larger than PZT"
        );
    }

    #[test]
    fn test_stress_induced_anisotropy_sign() {
        // Positive λ_s, tensile stress → K_σ < 0 (easy axis along stress)
        let k_sigma = stress_induced_anisotropy(25.0e-6, 100.0e6);
        assert!(
            k_sigma < 0.0,
            "K_σ should be negative for positive λ_s and tensile stress, got {}",
            k_sigma
        );

        // Negative λ_s (Ni), tensile stress → K_σ > 0 (easy axis perpendicular)
        let k_sigma_ni = stress_induced_anisotropy(-34.0e-6, 100.0e6);
        assert!(
            k_sigma_ni > 0.0,
            "K_σ should be positive for negative λ_s and tensile stress, got {}",
            k_sigma_ni
        );
    }

    #[test]
    fn test_polycrystalline_lambda_s_formula() {
        let fe = MagnetoelasticMaterial::iron();
        let computed = fe.compute_lambda_s();
        let expected = (2.0 / 5.0) * fe.lambda_100 + (3.0 / 5.0) * fe.lambda_111;
        assert!(
            (computed - expected).abs() < 1e-15,
            "compute_lambda_s should match formula"
        );
    }

    #[test]
    fn test_von_mises_stress_uniaxial() {
        let sigma = 100.0e6;
        let stress = StressTensor::uniaxial(sigma, 0).expect("should create uniaxial stress");
        let vm = stress.von_mises();
        // For uniaxial stress, von Mises = |σ|
        assert!(
            (vm - sigma.abs()).abs() / sigma.abs() < 1e-10,
            "von Mises for uniaxial should equal |σ|, got {} vs {}",
            vm,
            sigma.abs()
        );
    }
}
