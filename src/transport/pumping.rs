//! Spin pumping current calculation
//!
//! Spin pumping is a phenomenon where magnetization precession in a ferromagnet
//! generates a pure spin current into an adjacent normal metal.
//!
//! # Mathematical Formulation
//!
//! The spin current density pumped across a ferromagnet/normal-metal interface
//! is given by Tserkovnyak's formula:
//!
//! $$
//! \mathbf{J}_s = \frac{\hbar}{4\pi} g_r \left(\mathbf{m} \times \frac{d\mathbf{m}}{dt}\right)
//! $$
//!
//! where:
//! - $\mathbf{J}_s$ is the spin current density (J/m²)
//! - $\hbar$ is the reduced Planck constant ($1.055 \times 10^{-34}$ J·s)
//! - $g_r$ is the real part of spin-mixing conductance (1/m²)
//! - $\mathbf{m}$ is the normalized magnetization vector
//! - $d\mathbf{m}/dt$ is the time derivative of magnetization (1/s)
//!
//! The direction of the spin current is perpendicular to both $\mathbf{m}$ and $d\mathbf{m}/dt$,
//! and the spin polarization is along $\mathbf{m}$.
//!
//! # Physical Origin
//!
//! When magnetization precesses, it pumps angular momentum into the adjacent
//! layer through exchange interaction at the interface. The efficiency depends
//! on the spin-mixing conductance $g_r$, which characterizes interface transparency.
//!
//! Typical values:
//! - YIG/Pt: $g_r \sim 10^{14}$ m⁻²
//! - Metallic FM/NM: $g_r \sim 10^{15}$ m⁻²
//!
//! # References
//!
//! - E. Saitoh et al., "Conversion of spin current into charge current at room
//!   temperature: Inverse spin-Hall effect", Appl. Phys. Lett. 88, 182509 (2006)
//! - Y. Tserkovnyak et al., "Enhanced Gilbert Damping in Thin Ferromagnetic Films",
//!   Phys. Rev. Lett. 88, 117601 (2002)
//! - Y. Tserkovnyak et al., "Spin pumping and magnetization dynamics in metallic
//!   multilayers", Rev. Mod. Phys. 77, 1375 (2005)

use std::f64::consts::PI;

use crate::constants::HBAR;
use crate::material::SpinInterface;
use crate::vector3::Vector3;

/// Calculate spin pumping current density using Saitoh's formula
///
/// # Arguments
/// * `interface` - Spin interface properties (contains g_r)
/// * `m` - Normalized magnetization vector
/// * `dm_dt` - Time derivative of magnetization [1/s]
///
/// # Returns
/// Spin current density vector \[J/m²\]
///
/// The direction of the spin current is given by m × dm/dt, and the
/// spin polarization is along the m direction.
///
/// # Physical Interpretation
/// When magnetization precesses (m rotates), angular momentum is transferred
/// from the ferromagnet to conduction electrons at the interface, creating
/// a flow of spin angular momentum (spin current) into the normal metal.
///
/// # Example
/// ```
/// use spintronics::transport::pumping::spin_pumping_current;
/// use spintronics::material::SpinInterface;
/// use spintronics::dynamics::llg::calc_dm_dt;
/// use spintronics::constants::GAMMA;
/// use spintronics::Vector3;
///
/// // YIG/Pt interface (spin pumping experiment)
/// let interface = SpinInterface::yig_pt();
///
/// // Magnetization precessing in xy-plane
/// let m = Vector3::new(1.0, 0.0, 0.0);
/// // External field causes precession
/// let h_ext = Vector3::new(0.0, 0.0, 1.0);
/// let dm_dt = calc_dm_dt(m, h_ext, GAMMA, 0.01);
///
/// // Calculate pumped spin current
/// let js = spin_pumping_current(&interface, m, dm_dt);
///
/// // Spin current should be non-zero when magnetization is precessing
/// assert!(js.magnitude() > 0.0);
/// ```
#[inline]
pub fn spin_pumping_current(
    interface: &SpinInterface,
    m: Vector3<f64>,
    dm_dt: Vector3<f64>,
) -> Vector3<f64> {
    // Cross product: m × dm/dt
    // Physical meaning: The spin current direction is determined by the rate and
    // direction of magnetization precession. When m rotates (precesses), angular
    // momentum is transferred perpendicular to both m and its time derivative.
    // This creates a "pumping" of spin angular momentum into the adjacent layer.
    let cross_term = m.cross(&dm_dt);

    // Prefactor: ℏ / 4π [J·s]
    // Physical meaning: Quantum mechanical coupling strength. The factor ℏ comes
    // from quantization of angular momentum, and 4π from the solid angle integral
    // in the scattering theory (Tserkovnyak et al., Rev. Mod. Phys. 2005).
    let prefactor = HBAR / (4.0 * PI);

    // Spin-mixing conductance: g_r [1/m²]
    // Physical meaning: Measures the transparency of the FM/NM interface to spin
    // current transmission. Higher g_r means more efficient spin injection.
    // Typical values: YIG/Pt ~ 10¹⁴ 1/m², metals ~ 10¹⁵ 1/m²

    // Final result: J_s = (ℏ/4π) g_r (m × dm/dt) \[J/m²\]
    // Unit interpretation: \[J/m²\] = spin angular momentum flux density
    // Alternative unit: Can be expressed as \[A/m²\] when normalized by ℏ/2e
    cross_term * (prefactor * interface.g_r)
}

/// Calculate total spin pumping current (integrated over interface area)
///
/// # Arguments
/// * `interface` - Spin interface properties
/// * `m` - Normalized magnetization vector
/// * `dm_dt` - Time derivative of magnetization [1/s]
///
/// # Returns
/// Total spin current [J·m/s] flowing through the interface
pub fn total_spin_current(
    interface: &SpinInterface,
    m: Vector3<f64>,
    dm_dt: Vector3<f64>,
) -> Vector3<f64> {
    spin_pumping_current(interface, m, dm_dt) * interface.area
}

/// Calculate spin pumping efficiency
///
/// Returns the magnitude of spin current normalized by the precession frequency
pub fn pumping_efficiency(interface: &SpinInterface, m: Vector3<f64>, dm_dt: Vector3<f64>) -> f64 {
    let js = spin_pumping_current(interface, m, dm_dt);
    let omega = dm_dt.magnitude(); // Precession frequency

    if omega > 0.0 {
        js.magnitude() / omega
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spin_pumping_zero_dynamics() {
        let interface = SpinInterface::default();
        let m = Vector3::new(1.0, 0.0, 0.0);
        let dm_dt = Vector3::new(0.0, 0.0, 0.0);

        let js = spin_pumping_current(&interface, m, dm_dt);
        assert!(js.magnitude() < 1e-50);
    }

    #[test]
    fn test_spin_pumping_perpendicular() {
        let interface = SpinInterface::default();
        let m = Vector3::new(1.0, 0.0, 0.0);
        let dm_dt = Vector3::new(0.0, 1.0, 0.0);

        let js = spin_pumping_current(&interface, m, dm_dt);

        // Spin current should be perpendicular to both m and dm/dt
        assert!(js.dot(&m).abs() < 1e-40);
        assert!(js.dot(&dm_dt).abs() < 1e-40);
    }

    #[test]
    fn test_spin_pumping_magnitude() {
        let interface = SpinInterface {
            g_r: 1.0e19,
            ..Default::default()
        };
        let m = Vector3::new(1.0, 0.0, 0.0);
        let dm_dt = Vector3::new(0.0, 1.0e11, 0.0); // Fast precession

        let js = spin_pumping_current(&interface, m, dm_dt);
        assert!(js.magnitude() > 0.0);
    }
}
