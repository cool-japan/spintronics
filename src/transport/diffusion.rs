//! Spin diffusion in normal metals
//!
//! Describes the decay of spin current as it propagates through a normal metal
//! due to spin-flip scattering. The spin current decays exponentially with
//! characteristic length scale λ_sf (spin diffusion length).

use crate::vector3::Vector3;

/// Spin diffusion properties of a normal metal
#[derive(Debug, Clone)]
pub struct SpinDiffusion {
    /// Spin diffusion length \[m\]
    /// Typical values: Pt: ~10 nm, Cu: ~1000 nm, Al: ~500 nm
    pub lambda_sf: f64,

    /// Electrical conductivity \[S/m\]
    pub sigma: f64,
}

impl Default for SpinDiffusion {
    fn default() -> Self {
        Self {
            lambda_sf: 10.0e-9, // 10 nm (typical for Pt)
            sigma: 9.43e6,      // Pt conductivity
        }
    }
}

impl SpinDiffusion {
    /// Create spin diffusion properties for Platinum
    pub fn platinum() -> Self {
        Self {
            lambda_sf: 10.0e-9,
            sigma: 9.43e6,
        }
    }

    /// Create spin diffusion properties for Copper
    pub fn copper() -> Self {
        Self {
            lambda_sf: 1000.0e-9,
            sigma: 5.96e7,
        }
    }

    /// Create spin diffusion properties for Aluminum
    pub fn aluminum() -> Self {
        Self {
            lambda_sf: 500.0e-9,
            sigma: 3.77e7,
        }
    }

    /// Calculate spin current after diffusion over distance x
    ///
    /// J_s(x) = J_s(0) * exp(-x / λ_sf)
    ///
    /// # Arguments
    /// * `js_0` - Initial spin current at x=0
    /// * `distance` - Distance traveled \[m\]
    ///
    /// # Returns
    /// Spin current at distance x
    pub fn decay(&self, js_0: Vector3<f64>, distance: f64) -> Vector3<f64> {
        let decay_factor = (-distance / self.lambda_sf).exp();
        js_0 * decay_factor
    }

    /// Calculate effective decay factor for a given thickness
    pub fn transmission_factor(&self, thickness: f64) -> f64 {
        (-thickness / self.lambda_sf).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spin_diffusion_zero_distance() {
        let diffusion = SpinDiffusion::default();
        let js_0 = Vector3::new(1.0, 0.0, 0.0);
        let js = diffusion.decay(js_0, 0.0);

        assert!((js.x - 1.0).abs() < 1e-10);
        assert!(js.y.abs() < 1e-10);
        assert!(js.z.abs() < 1e-10);
    }

    #[test]
    fn test_spin_diffusion_decay() {
        let diffusion = SpinDiffusion::platinum();
        let js_0 = Vector3::new(1.0, 0.0, 0.0);

        // At one spin diffusion length, current should be 1/e
        let js = diffusion.decay(js_0, diffusion.lambda_sf);
        assert!((js.magnitude() - std::f64::consts::E.recip()).abs() < 1e-10);
    }

    #[test]
    fn test_transmission_factor() {
        let diffusion = SpinDiffusion::platinum();
        let factor = diffusion.transmission_factor(diffusion.lambda_sf);
        assert!((factor - std::f64::consts::E.recip()).abs() < 1e-10);
    }
}
