//! Quantized spin wave modes in confined nanostructures
//!
//! When spin waves are confined to finite geometries (stripes, disks, rectangles),
//! the allowed wavevectors become quantized. This module computes the discrete
//! mode spectrum for common nanostructure geometries.
//!
//! # Physical Background
//!
//! In a stripe of width w, standing spin waves form with quantized wavevectors
//! k_n = n*pi/w (for totally pinned boundary conditions) or with modified boundary
//! conditions accounting for the effective pinning parameter.
//!
//! In circular disk geometries, the modes are described by Bessel functions
//! J_m(k_{m,n} * r), where k_{m,n} = alpha_{m,n} / R with alpha_{m,n} being
//! the zeros of the appropriate Bessel function or its derivative.
//!
//! # References
//!
//! - K. Yu. Guslienko et al., "Effective dipolar boundary conditions for dynamic
//!   magnetization in thin magnetic stripes", Phys. Rev. B 66, 132402 (2002)
//! - C. Bayer et al., "Spin-wave excitations in finite rectangular elements of
//!   Ni80Fe20", Phys. Rev. B 72, 064427 (2005)

use crate::constants::{GAMMA, MU_0};
use crate::error::{self, Result};
use crate::material::Ferromagnet;

/// Geometry of a confined magnetic nanostructure
#[derive(Debug, Clone)]
pub enum NanostructureGeometry {
    /// Rectangular stripe with given width \[m\] and length \[m\]
    Stripe {
        /// Width of the stripe \[m\] (confinement direction)
        width: f64,
        /// Length of the stripe \[m\] (propagation direction)
        length: f64,
    },
    /// Circular disk with given radius \[m\]
    Disk {
        /// Radius of the disk \[m\]
        radius: f64,
    },
    /// Rectangular element with width \[m\] and length \[m\]
    Rectangle {
        /// Width \[m\]
        width: f64,
        /// Length \[m\]
        length: f64,
    },
}

/// Quantized spin wave mode calculator for nanostructures
///
/// Computes the discrete mode spectrum of spin waves confined in finite geometries.
///
/// # Example
///
/// ```
/// use spintronics::spinwave::{QuantizedModes, NanostructureGeometry};
/// use spintronics::material::Ferromagnet;
///
/// let py = Ferromagnet::permalloy();
/// let geo = NanostructureGeometry::Stripe { width: 1e-6, length: 10e-6 };
/// let qm = QuantizedModes::new(&py, geo, 50e-9)
///     .expect("valid parameters");
///
/// let freqs = qm.stripe_mode_frequencies(0.05, 5)
///     .expect("valid calculation");
/// assert_eq!(freqs.len(), 5);
/// ```
#[derive(Debug, Clone)]
pub struct QuantizedModes {
    /// Saturation magnetization \[A/m\]
    ms: f64,
    /// Exchange stiffness \[J/m\]
    #[allow(dead_code)]
    exchange_a: f64,
    /// Film thickness \[m\]
    thickness: f64,
    /// Nanostructure geometry
    geometry: NanostructureGeometry,
    /// Exchange length squared D = 2A/(mu_0*Ms) [m^2]
    exchange_length_sq: f64,
}

impl QuantizedModes {
    /// Create a new quantized mode calculator
    ///
    /// # Arguments
    /// * `material` - Ferromagnetic material parameters
    /// * `geometry` - Nanostructure geometry specification
    /// * `thickness` - Film thickness \[m\]
    ///
    /// # Errors
    /// Returns error for invalid parameters
    pub fn new(
        material: &Ferromagnet,
        geometry: NanostructureGeometry,
        thickness: f64,
    ) -> Result<Self> {
        if material.ms <= 0.0 {
            return Err(error::invalid_param(
                "ms",
                "saturation magnetization must be positive",
            ));
        }
        if thickness <= 0.0 {
            return Err(error::invalid_param(
                "thickness",
                "film thickness must be positive",
            ));
        }

        match &geometry {
            NanostructureGeometry::Stripe { width, length } => {
                if *width <= 0.0 || *length <= 0.0 {
                    return Err(error::invalid_param(
                        "geometry",
                        "stripe width and length must be positive",
                    ));
                }
            },
            NanostructureGeometry::Disk { radius } => {
                if *radius <= 0.0 {
                    return Err(error::invalid_param(
                        "geometry",
                        "disk radius must be positive",
                    ));
                }
            },
            NanostructureGeometry::Rectangle { width, length } => {
                if *width <= 0.0 || *length <= 0.0 {
                    return Err(error::invalid_param(
                        "geometry",
                        "rectangle width and length must be positive",
                    ));
                }
            },
        }

        let exchange_length_sq = 2.0 * material.exchange_a / (MU_0 * material.ms);

        Ok(Self {
            ms: material.ms,
            exchange_a: material.exchange_a,
            thickness,
            geometry,
            exchange_length_sq,
        })
    }

    /// Quantized wavevectors for a stripe geometry: k_n = n * pi / w
    ///
    /// Returns the first `n_modes` quantized wavevectors assuming totally
    /// pinned boundary conditions at the stripe edges.
    ///
    /// # Arguments
    /// * `n_modes` - Number of modes to compute (n = 1, 2, ..., n_modes)
    ///
    /// # Returns
    /// Vector of quantized wavevectors \[rad/m\]
    pub fn stripe_wavevectors(&self, n_modes: usize) -> Result<Vec<f64>> {
        let width = match &self.geometry {
            NanostructureGeometry::Stripe { width, .. } => *width,
            NanostructureGeometry::Rectangle { width, .. } => *width,
            _ => {
                return Err(error::invalid_param(
                    "geometry",
                    "stripe wavevectors require Stripe or Rectangle geometry",
                ))
            },
        };

        if n_modes == 0 {
            return Ok(Vec::new());
        }

        let wavevectors: Vec<f64> = (1..=n_modes)
            .map(|n| n as f64 * std::f64::consts::PI / width)
            .collect();

        Ok(wavevectors)
    }

    /// Compute mode frequencies for a stripe geometry
    ///
    /// Uses the exchange spin wave dispersion omega(k_n) = gamma*(H + D*k_n^2)
    /// combined with the quantized wavevectors k_n = n*pi/w.
    ///
    /// For a more complete treatment including dipolar effects, the Kalinikos-Slavin
    /// formula with the quantized wavevectors would be used.
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `n_modes` - Number of modes to compute
    ///
    /// # Returns
    /// Vector of (mode_number, frequency \[rad/s\]) pairs
    pub fn stripe_mode_frequencies(&self, h_ext: f64, n_modes: usize) -> Result<Vec<(usize, f64)>> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }

        let wavevectors = self.stripe_wavevectors(n_modes)?;
        let d = self.exchange_length_sq;
        let omega_m = GAMMA * MU_0 * self.ms;
        let omega_h = GAMMA * h_ext;

        let mut modes = Vec::with_capacity(n_modes);
        let film_d = self.thickness;

        for (i, &k_n) in wavevectors.iter().enumerate() {
            let n = i + 1;

            // Full dipole-exchange: use Kalinikos-Slavin-like formula
            // For k perpendicular to M (DE-like geometry in stripe):
            let kd = k_n * film_d;
            let p = if kd < 1e-12 {
                1.0 - kd / 2.0
            } else {
                (1.0 - (-kd).exp()) / kd
            };

            let exchange_term = omega_m * d * k_n * k_n;

            // Use the full formula with dipolar contributions
            let term1 = omega_h + exchange_term;
            let term2 = omega_h + exchange_term + omega_m * p;
            let omega_sq = term1 * term2;

            if omega_sq >= 0.0 {
                modes.push((n, omega_sq.sqrt()));
            }
        }

        Ok(modes)
    }

    /// Approximate zeros of Bessel function J_m(x) for disk mode quantization
    ///
    /// Returns the first `n_zeros` positive zeros of J_m(x) using McMahon's
    /// asymptotic expansion, which is accurate for all zeros.
    ///
    /// J_m(alpha_{m,n}) = 0 where alpha_{m,n} are the returned values.
    ///
    /// # Arguments
    /// * `m` - Azimuthal mode number (order of Bessel function)
    /// * `n_zeros` - Number of zeros to compute
    ///
    /// # Returns
    /// Vector of approximate zeros of J_m(x)
    pub fn bessel_zeros(m: u32, n_zeros: usize) -> Vec<f64> {
        if n_zeros == 0 {
            return Vec::new();
        }

        // Use McMahon's asymptotic expansion for zeros of J_m(x)
        // alpha_{m,s} ~ beta - (mu-1)/(8*beta) - 4*(mu-1)*(7*mu-31)/(3*(8*beta)^3) - ...
        // where beta = pi*(s + m/2 - 1/4), mu = 4*m^2
        //
        // For small s, use known tabulated values for better accuracy

        // Known exact zeros for low orders (from Abramowitz & Stegun)
        let known_zeros: &[&[f64]] = &[
            // J_0 zeros
            &[
                2.4048, 5.5201, 8.6537, 11.7915, 14.9309, 18.0711, 21.2116, 24.3525, 27.4935,
                30.6346,
            ],
            // J_1 zeros
            &[
                3.8317, 7.0156, 10.1735, 13.3237, 16.4706, 19.6159, 22.7601, 25.9037, 29.0468,
                32.1897,
            ],
            // J_2 zeros
            &[
                5.1356, 8.4172, 11.6198, 14.7960, 17.9598, 21.1170, 24.2701, 27.4206, 30.5692,
                33.7165,
            ],
            // J_3 zeros
            &[
                6.3802, 9.7610, 13.0152, 16.2235, 19.4094, 22.5828, 25.7482, 28.9084, 32.0649,
                35.2187,
            ],
            // J_4 zeros
            &[
                7.5883, 11.0647, 14.3725, 17.6160, 20.8269, 24.0190, 27.1991, 30.3710, 33.5371,
                36.6990,
            ],
        ];

        let mut zeros = Vec::with_capacity(n_zeros);

        for s in 0..n_zeros {
            let zero = if (m as usize) < known_zeros.len() && s < known_zeros[m as usize].len() {
                known_zeros[m as usize][s]
            } else {
                // McMahon's expansion for higher zeros
                let mu = 4.0 * (m as f64) * (m as f64);
                let beta = std::f64::consts::PI * ((s + 1) as f64 + m as f64 / 2.0 - 0.25);
                let b8 = 8.0 * beta;
                let correction1 = (mu - 1.0) / b8;
                let correction2 = 4.0 * (mu - 1.0) * (7.0 * mu - 31.0) / (3.0 * b8 * b8 * b8);
                beta - correction1 - correction2
            };
            zeros.push(zero);
        }

        zeros
    }

    /// Quantized wavevectors for a circular disk geometry
    ///
    /// The modes in a disk are characterized by radial index n and azimuthal
    /// index m, with wavevectors k_{m,n} = alpha_{m,n} / R where alpha_{m,n}
    /// are zeros of J_m (for pinned boundaries) or J_m' (for free boundaries).
    ///
    /// # Arguments
    /// * `m` - Azimuthal mode number
    /// * `n_radial` - Number of radial modes
    ///
    /// # Returns
    /// Vector of (radial_index, wavevector \[rad/m\]) pairs
    pub fn disk_wavevectors(&self, m: u32, n_radial: usize) -> Result<Vec<(usize, f64)>> {
        let radius = match &self.geometry {
            NanostructureGeometry::Disk { radius } => *radius,
            _ => {
                return Err(error::invalid_param(
                    "geometry",
                    "disk wavevectors require Disk geometry",
                ))
            },
        };

        let zeros = Self::bessel_zeros(m, n_radial);
        let wavevectors: Vec<(usize, f64)> = zeros
            .into_iter()
            .enumerate()
            .map(|(i, alpha)| (i + 1, alpha / radius))
            .collect();

        Ok(wavevectors)
    }

    /// Compute mode frequencies for a circular disk
    ///
    /// Uses exchange spin wave dispersion with quantized disk wavevectors.
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `m` - Azimuthal mode number
    /// * `n_radial` - Number of radial modes
    ///
    /// # Returns
    /// Vector of (radial_index, azimuthal_index, frequency \[rad/s\]) tuples
    pub fn disk_mode_frequencies(
        &self,
        h_ext: f64,
        m: u32,
        n_radial: usize,
    ) -> Result<Vec<(usize, u32, f64)>> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }

        let wavevectors = self.disk_wavevectors(m, n_radial)?;
        let d = self.exchange_length_sq;

        let mut modes = Vec::with_capacity(n_radial);
        for (n, k_mn) in wavevectors {
            let omega = GAMMA * (h_ext + MU_0 * d * k_mn * k_mn);
            modes.push((n, m, omega));
        }

        Ok(modes)
    }

    /// Compute the full mode spectrum for a rectangular element
    ///
    /// In a rectangle of width w and length l, modes are doubly quantized:
    /// k_x = n_x * pi / w, k_y = n_y * pi / l
    /// k^2 = k_x^2 + k_y^2
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[T\]
    /// * `n_max_x` - Maximum mode number in x-direction
    /// * `n_max_y` - Maximum mode number in y-direction
    ///
    /// # Returns
    /// Vector of (n_x, n_y, frequency \[rad/s\]) tuples sorted by frequency
    pub fn rectangle_mode_spectrum(
        &self,
        h_ext: f64,
        n_max_x: usize,
        n_max_y: usize,
    ) -> Result<Vec<(usize, usize, f64)>> {
        if h_ext < 0.0 {
            return Err(error::invalid_param(
                "h_ext",
                "external field must be non-negative",
            ));
        }

        let (width, length) = match &self.geometry {
            NanostructureGeometry::Rectangle { width, length } => (*width, *length),
            NanostructureGeometry::Stripe { width, length } => (*width, *length),
            _ => {
                return Err(error::invalid_param(
                    "geometry",
                    "rectangle mode spectrum requires Rectangle or Stripe geometry",
                ))
            },
        };

        let d = self.exchange_length_sq;
        let omega_m = GAMMA * MU_0 * self.ms;
        let omega_h = GAMMA * h_ext;
        let film_d = self.thickness;

        let mut modes = Vec::new();

        for nx in 1..=n_max_x {
            for ny in 1..=n_max_y {
                let kx = nx as f64 * std::f64::consts::PI / width;
                let ky = ny as f64 * std::f64::consts::PI / length;
                let k_sq = kx * kx + ky * ky;
                let k = k_sq.sqrt();

                let kd = k * film_d;
                let p = if kd < 1e-12 {
                    1.0 - kd / 2.0
                } else {
                    (1.0 - (-kd).exp()) / kd
                };

                // Angle phi between k and M (assume M along x)
                let sin2_phi = if k_sq > 0.0 { ky * ky / k_sq } else { 0.0 };

                let exchange_term = omega_m * d * k_sq;
                let f_nn = 1.0 - p + sin2_phi * p;

                let term1 = omega_h + exchange_term;
                let term2 = omega_h + exchange_term + omega_m * f_nn;
                let omega_sq = term1 * term2;

                if omega_sq >= 0.0 {
                    modes.push((nx, ny, omega_sq.sqrt()));
                }
            }
        }

        modes.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
        Ok(modes)
    }

    /// Standing wave condition: quantized wavevectors k_n = n * pi / L
    ///
    /// Simple utility function for computing standing wave wavevectors.
    ///
    /// # Arguments
    /// * `length` - Confinement length \[m\]
    /// * `n` - Mode number (1, 2, 3, ...)
    ///
    /// # Returns
    /// Quantized wavevector \[rad/m\]
    pub fn standing_wave_k(length: f64, n: usize) -> Result<f64> {
        if length <= 0.0 {
            return Err(error::invalid_param(
                "length",
                "confinement length must be positive",
            ));
        }
        if n == 0 {
            return Err(error::invalid_param(
                "n",
                "mode number must be positive (n >= 1)",
            ));
        }
        Ok(n as f64 * std::f64::consts::PI / length)
    }

    /// Get the geometry of this calculator
    pub fn geometry(&self) -> &NanostructureGeometry {
        &self.geometry
    }

    /// Exchange length squared [m^2]
    pub fn exchange_length_squared(&self) -> f64 {
        self.exchange_length_sq
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    fn stripe_calculator() -> QuantizedModes {
        let py = Ferromagnet::permalloy();
        let geo = NanostructureGeometry::Stripe {
            width: 1e-6,
            length: 10e-6,
        };
        QuantizedModes::new(&py, geo, 50e-9).expect("valid parameters")
    }

    fn disk_calculator() -> QuantizedModes {
        let py = Ferromagnet::permalloy();
        let geo = NanostructureGeometry::Disk { radius: 500e-9 };
        QuantizedModes::new(&py, geo, 50e-9).expect("valid parameters")
    }

    #[test]
    fn test_standing_wave_k() {
        let length = 1e-6; // 1 um
        let k1 = QuantizedModes::standing_wave_k(length, 1).expect("valid");
        let k2 = QuantizedModes::standing_wave_k(length, 2).expect("valid");

        assert!((k1 - PI / length).abs() < 1.0, "k1 should be pi/L");
        assert!((k2 - 2.0 * PI / length).abs() < 1.0, "k2 should be 2*pi/L");
        assert!((k2 / k1 - 2.0).abs() < 1e-10, "k2/k1 should be 2");
    }

    #[test]
    fn test_stripe_wavevectors() {
        let calc = stripe_calculator();
        let kvecs = calc.stripe_wavevectors(5).expect("valid");
        assert_eq!(kvecs.len(), 5);

        // k_n should increase linearly with n
        for i in 1..kvecs.len() {
            assert!(kvecs[i] > kvecs[i - 1], "wavevectors should be increasing");
        }

        // k_n should be n*pi/w
        let width = 1e-6;
        for (i, &k) in kvecs.iter().enumerate() {
            let expected = (i + 1) as f64 * PI / width;
            let rel_diff = (k - expected).abs() / expected;
            assert!(rel_diff < 1e-10, "k[{i}] = {k}, expected {expected}");
        }
    }

    #[test]
    fn test_stripe_mode_frequencies_increasing() {
        let calc = stripe_calculator();
        let modes = calc
            .stripe_mode_frequencies(0.05, 5)
            .expect("valid calculation");
        assert_eq!(modes.len(), 5);

        // Frequencies should increase with mode number
        for i in 1..modes.len() {
            assert!(
                modes[i].1 > modes[i - 1].1,
                "frequencies should increase: f[{}]={} <= f[{}]={}",
                i - 1,
                modes[i - 1].1,
                i,
                modes[i].1
            );
        }
    }

    #[test]
    fn test_disk_wavevectors() {
        let calc = disk_calculator();
        let kvecs = calc.disk_wavevectors(0, 3).expect("valid");
        assert_eq!(kvecs.len(), 3);

        // Wavevectors should increase with radial index
        for i in 1..kvecs.len() {
            assert!(
                kvecs[i].1 > kvecs[i - 1].1,
                "disk wavevectors should increase"
            );
        }
    }

    #[test]
    fn test_disk_mode_frequencies() {
        let calc = disk_calculator();
        let modes = calc
            .disk_mode_frequencies(0.05, 0, 3)
            .expect("valid calculation");
        assert_eq!(modes.len(), 3);

        // Frequencies should increase with radial index
        for i in 1..modes.len() {
            assert!(
                modes[i].2 > modes[i - 1].2,
                "disk frequencies should increase with radial index"
            );
        }
    }

    #[test]
    fn test_bessel_zeros_j0() {
        let zeros = QuantizedModes::bessel_zeros(0, 3);
        assert_eq!(zeros.len(), 3);

        // Known J_0 zeros: 2.4048, 5.5201, 8.6537
        assert!((zeros[0] - 2.4048).abs() < 0.001);
        assert!((zeros[1] - 5.5201).abs() < 0.001);
        assert!((zeros[2] - 8.6537).abs() < 0.001);
    }

    #[test]
    fn test_bessel_zeros_j1() {
        let zeros = QuantizedModes::bessel_zeros(1, 2);
        assert_eq!(zeros.len(), 2);

        // Known J_1 zeros: 3.8317, 7.0156
        assert!((zeros[0] - 3.8317).abs() < 0.001);
        assert!((zeros[1] - 7.0156).abs() < 0.001);
    }

    #[test]
    fn test_rectangle_mode_spectrum() {
        let py = Ferromagnet::permalloy();
        let geo = NanostructureGeometry::Rectangle {
            width: 1e-6,
            length: 2e-6,
        };
        let calc = QuantizedModes::new(&py, geo, 50e-9).expect("valid");
        let modes = calc
            .rectangle_mode_spectrum(0.05, 3, 3)
            .expect("valid calculation");

        // Should have 3x3 = 9 modes
        assert_eq!(modes.len(), 9);

        // Should be sorted by frequency
        for i in 1..modes.len() {
            assert!(
                modes[i].2 >= modes[i - 1].2,
                "modes should be sorted by frequency"
            );
        }
    }

    #[test]
    fn test_quantized_mode_spacing() {
        // For large mode numbers, the spacing should scale as ~n (exchange dominated)
        let calc = stripe_calculator();
        let modes = calc
            .stripe_mode_frequencies(0.05, 10)
            .expect("valid calculation");

        // The frequency should increase faster than linearly (quadratic in k for exchange)
        let delta_low = modes[1].1 - modes[0].1;
        let delta_high = modes[9].1 - modes[8].1;

        assert!(
            delta_high > delta_low,
            "frequency spacing should increase with mode number (exchange-dominated): \
             delta_low={delta_low}, delta_high={delta_high}"
        );
    }

    #[test]
    fn test_invalid_geometry() {
        let py = Ferromagnet::permalloy();
        let bad_geo = NanostructureGeometry::Stripe {
            width: -1.0,
            length: 1e-6,
        };
        assert!(QuantizedModes::new(&py, bad_geo, 50e-9).is_err());

        let bad_disk = NanostructureGeometry::Disk { radius: 0.0 };
        assert!(QuantizedModes::new(&py, bad_disk, 50e-9).is_err());
    }

    #[test]
    fn test_standing_wave_k_errors() {
        assert!(QuantizedModes::standing_wave_k(-1.0, 1).is_err());
        assert!(QuantizedModes::standing_wave_k(1.0, 0).is_err());
    }

    #[test]
    fn test_disk_geometry_mismatch() {
        let calc = stripe_calculator();
        assert!(calc.disk_wavevectors(0, 3).is_err());
    }
}
