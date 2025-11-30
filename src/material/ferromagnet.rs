//! Ferromagnetic material properties
//!
//! This module defines the parameters needed for the Landau-Lifshitz-Gilbert (LLG)
//! equation, which describes magnetization dynamics in ferromagnetic materials.

use crate::vector3::Vector3;

/// Ferromagnetic material with LLG parameters
///
/// Contains all necessary parameters for simulating magnetization dynamics
/// in ferromagnetic materials like Permalloy (Py), YIG (Yttrium Iron Garnet), etc.
///
/// # Example
/// ```
/// use spintronics::material::Ferromagnet;
///
/// // Use pre-defined materials from database
/// let yig = Ferromagnet::yig();
/// assert_eq!(yig.alpha, 0.0001); // Ultra-low damping
///
/// let permalloy = Ferromagnet::permalloy();
/// assert_eq!(permalloy.ms, 8.0e5); // High saturation magnetization
///
/// // Or create custom material
/// let custom = Ferromagnet {
///     alpha: 0.005,
///     ms: 5.0e5,
///     anisotropy_k: 1.0e4,
///     easy_axis: spintronics::Vector3::new(1.0, 1.0, 0.0).normalize(),
///     exchange_a: 1.5e-11,
/// };
/// ```
#[derive(Debug, Clone)]
pub struct Ferromagnet {
    /// Gilbert damping constant (dimensionless)
    /// Typical values: 0.001-0.1 depending on material
    /// YIG: ~0.0001-0.001, Permalloy: ~0.01
    pub alpha: f64,

    /// Saturation magnetization \[A/m\]
    /// YIG: ~1.4e5, Permalloy: ~8e5
    pub ms: f64,

    /// Uniaxial anisotropy constant \[J/m³\]
    /// Determines the energy barrier for magnetization rotation
    pub anisotropy_k: f64,

    /// Easy axis direction (normalized)
    /// Preferred direction of magnetization due to crystal structure
    pub easy_axis: Vector3<f64>,

    /// Exchange stiffness constant \[J/m\]
    /// Typical values: 1e-11 to 1e-12 for most ferromagnets
    pub exchange_a: f64,
}

impl Default for Ferromagnet {
    fn default() -> Self {
        Self {
            alpha: 0.01,
            ms: 1.4e5,
            anisotropy_k: 0.0,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_a: 1.0e-11,
        }
    }
}

impl Ferromagnet {
    /// Create a YIG (Yttrium Iron Garnet) material
    ///
    /// YIG is known for its extremely low damping, making it ideal for
    /// spin pumping experiments as demonstrated by Saitoh et al.
    pub fn yig() -> Self {
        Self {
            alpha: 0.0001,
            ms: 1.4e5,
            anisotropy_k: 0.0,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_a: 3.7e-12,
        }
    }

    /// Create a Permalloy (Ni80Fe20) material
    pub fn permalloy() -> Self {
        Self {
            alpha: 0.01,
            ms: 8.0e5,
            anisotropy_k: 0.0,
            easy_axis: Vector3::new(1.0, 0.0, 0.0),
            exchange_a: 1.3e-11,
        }
    }

    /// Create a CoFeB material
    ///
    /// Co₂₀Fe₆₀B₂₀ commonly used in magnetic tunnel junctions (MTJs)
    /// Reference: S. Ikeda et al., Nat. Mater. 9, 721 (2010)
    pub fn cofeb() -> Self {
        Self {
            alpha: 0.004,
            ms: 1.0e6,
            anisotropy_k: 1.0e5,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_a: 2.0e-11,
        }
    }

    /// Create pure Iron (Fe) material
    ///
    /// Body-centered cubic (BCC) iron at room temperature
    /// Reference: Kittel, Introduction to Solid State Physics
    pub fn iron() -> Self {
        Self {
            alpha: 0.002,
            ms: 1.71e6,          // Saturation magnetization at room temperature
            anisotropy_k: 4.8e4, // Cubic anisotropy
            easy_axis: Vector3::new(1.0, 0.0, 0.0), // <100> easy axis
            exchange_a: 2.1e-11,
        }
    }

    /// Create pure Cobalt (Co) material
    ///
    /// Hexagonal close-packed (HCP) cobalt at room temperature
    /// Strong perpendicular magnetic anisotropy (PMA)
    pub fn cobalt() -> Self {
        Self {
            alpha: 0.005,
            ms: 1.4e6,
            anisotropy_k: 5.0e5,                    // Strong PMA
            easy_axis: Vector3::new(0.0, 0.0, 1.0), // c-axis (perpendicular)
            exchange_a: 3.0e-11,
        }
    }

    /// Create pure Nickel (Ni) material
    ///
    /// Face-centered cubic (FCC) nickel at room temperature
    pub fn nickel() -> Self {
        Self {
            alpha: 0.045, // Higher damping than Fe or Co
            ms: 4.8e5,
            anisotropy_k: -4.5e3, // Negative (cubic anisotropy)
            easy_axis: Vector3::new(1.0, 1.0, 1.0), // <111> easy axis
            exchange_a: 0.9e-11,
        }
    }

    /// Create CoFe alloy (Co₅₀Fe₅₀)
    ///
    /// Commonly used in spin valves and MTJs for high spin polarization
    /// Reference: J. Magn. Magn. Mater. 159, L1 (1996)
    pub fn cofe() -> Self {
        Self {
            alpha: 0.003,
            ms: 1.95e6, // Very high saturation magnetization
            anisotropy_k: 1.0e4,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
            exchange_a: 2.5e-11,
        }
    }

    /// Create Py (Permalloy, Ni₈₁Fe₁₉) - alias for permalloy()
    ///
    /// Classic soft magnetic material with near-zero anisotropy
    pub fn py() -> Self {
        Self::permalloy()
    }

    /// Builder method to set Gilbert damping
    pub fn with_alpha(mut self, alpha: f64) -> Self {
        self.alpha = alpha;
        self
    }

    /// Builder method to set saturation magnetization
    pub fn with_ms(mut self, ms: f64) -> Self {
        self.ms = ms;
        self
    }

    /// Builder method to set anisotropy constant
    pub fn with_anisotropy(mut self, k: f64) -> Self {
        self.anisotropy_k = k;
        self
    }

    /// Builder method to set easy axis
    pub fn with_easy_axis(mut self, axis: Vector3<f64>) -> Self {
        self.easy_axis = axis.normalize();
        self
    }

    /// Builder method to set exchange stiffness
    pub fn with_exchange(mut self, a: f64) -> Self {
        self.exchange_a = a;
        self
    }

    /// Create a custom ferromagnet with all parameters
    pub fn custom(
        alpha: f64,
        ms: f64,
        anisotropy_k: f64,
        easy_axis: Vector3<f64>,
        exchange_a: f64,
    ) -> Self {
        Self {
            alpha,
            ms,
            anisotropy_k,
            easy_axis: easy_axis.normalize(),
            exchange_a,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_ferromagnet() {
        let fm = Ferromagnet::default();
        assert!(fm.alpha > 0.0);
        assert!(fm.ms > 0.0);
    }

    #[test]
    fn test_yig_properties() {
        let yig = Ferromagnet::yig();
        assert!(yig.alpha < 0.001);
        assert!((yig.ms - 1.4e5).abs() < 1.0);
    }

    #[test]
    fn test_permalloy_properties() {
        let py = Ferromagnet::permalloy();
        assert!((py.ms - 8.0e5).abs() < 1.0);
        assert_eq!(py.anisotropy_k, 0.0); // Soft magnetic material
    }

    #[test]
    fn test_cofeb_properties() {
        let cofeb = Ferromagnet::cofeb();
        assert!(cofeb.ms > 9.0e5);
        assert!(cofeb.anisotropy_k > 0.0); // Has PMA
    }

    #[test]
    fn test_iron_properties() {
        let fe = Ferromagnet::iron();
        assert!(fe.ms > 1.5e6); // High saturation magnetization
        assert!(fe.alpha < 0.01); // Low damping
    }

    #[test]
    fn test_cobalt_properties() {
        let co = Ferromagnet::cobalt();
        assert!(co.anisotropy_k > 1.0e5); // Strong PMA
        assert!((co.easy_axis.z - 1.0).abs() < 1e-10); // Perpendicular easy axis
    }

    #[test]
    fn test_nickel_properties() {
        let ni = Ferromagnet::nickel();
        assert!(ni.alpha > 0.04); // Higher damping
        assert!(ni.ms < 5.0e5);
    }

    #[test]
    fn test_cofe_properties() {
        let cofe = Ferromagnet::cofe();
        assert!(cofe.ms > 1.9e6); // Very high Ms
        assert!(cofe.alpha < 0.01); // Low damping
    }

    #[test]
    fn test_builder_pattern() {
        let custom = Ferromagnet::yig().with_alpha(0.01).with_ms(2.0e5);

        assert_eq!(custom.alpha, 0.01);
        assert_eq!(custom.ms, 2.0e5);
    }

    #[test]
    fn test_custom_ferromagnet() {
        let custom = Ferromagnet::custom(0.015, 1.0e6, 1.0e5, Vector3::new(1.0, 1.0, 0.0), 1.5e-11);

        assert_eq!(custom.alpha, 0.015);
        assert_eq!(custom.ms, 1.0e6);
        // Easy axis should be normalized
        let norm =
            (custom.easy_axis.x.powi(2) + custom.easy_axis.y.powi(2) + custom.easy_axis.z.powi(2))
                .sqrt();
        assert!((norm - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_py_alias() {
        let py1 = Ferromagnet::py();
        let py2 = Ferromagnet::permalloy();

        assert_eq!(py1.alpha, py2.alpha);
        assert_eq!(py1.ms, py2.ms);
    }
}
