//! Runtime physical parameter checks.
//!
//! Functions in this module return `Err(...)` for hard violations (negative
//! temperatures, damping >= 1, etc.) and rely on `debug_assert!` for soft
//! sanity bounds that should never trigger in well-formed simulations but
//! catch obvious unit-conversion mistakes during development.

use crate::error::invalid_param;
use crate::vector3::Vector3;

/// Check that temperature is non-negative
///
/// # Panics (debug builds)
/// Panics if temperature is negative
#[inline]
pub fn check_temperature(t: f64) -> crate::error::Result<()> {
    if t < 0.0 {
        return Err(invalid_param("temperature", "must be non-negative"));
    }
    debug_assert!(t < 10000.0, "Temperature unreasonably high: {} K", t);
    Ok(())
}

/// Check that magnetization magnitude is reasonable
///
/// # Panics (debug builds)
/// Panics if magnetization is negative or unreasonably large
#[inline]
pub fn check_magnetization(ms: f64) -> crate::error::Result<()> {
    if ms < 0.0 {
        return Err(invalid_param("magnetization", "must be non-negative"));
    }
    debug_assert!(ms < 2.0e6, "Magnetization unreasonably high: {} A/m", ms);
    Ok(())
}

/// Check that damping parameter is in valid range
///
/// # Panics (debug builds)
/// Panics if damping is negative or >= 1
#[inline]
pub fn check_damping(alpha: f64) -> crate::error::Result<()> {
    if alpha < 0.0 {
        return Err(invalid_param("damping", "must be non-negative"));
    }
    if alpha >= 1.0 {
        return Err(invalid_param("damping", "must be < 1"));
    }
    debug_assert!(alpha < 0.5, "Damping unusually high: {}", alpha);
    Ok(())
}

/// Check that a vector is normalized (within tolerance)
///
/// # Panics (debug builds)
/// Panics if vector magnitude differs significantly from 1.0
#[inline]
pub fn check_normalized(v: Vector3<f64>, tolerance: f64) {
    let mag = v.magnitude();
    debug_assert!(
        (mag - 1.0).abs() < tolerance,
        "Vector not normalized: magnitude = {}",
        mag
    );
}

/// Check that time step is positive and reasonable
///
/// # Panics (debug builds)
/// Panics if time step is non-positive or unreasonably large
#[inline]
pub fn check_time_step(dt: f64) -> crate::error::Result<()> {
    if dt <= 0.0 {
        return Err(invalid_param("time step", "must be positive"));
    }
    debug_assert!(dt < 1e-6, "Time step unreasonably large: {} s", dt);
    Ok(())
}

/// Check that exchange constant is positive
#[inline]
pub fn check_exchange(a_ex: f64) -> crate::error::Result<()> {
    if a_ex <= 0.0 {
        return Err(invalid_param("exchange constant", "must be positive"));
    }
    debug_assert!(
        a_ex < 1e-10,
        "Exchange constant unreasonably large: {} J/m",
        a_ex
    );
    Ok(())
}

/// Check that current density is reasonable
#[inline]
pub fn check_current_density(j: f64) {
    debug_assert!(
        j.abs() < 1e13,
        "Current density unreasonably high: {} A/m²",
        j
    );
}

/// Check that spin Hall angle is in valid range
#[inline]
pub fn check_spin_hall_angle(theta_sh: f64) -> crate::error::Result<()> {
    if theta_sh.abs() > 1.0 {
        return Err(invalid_param(
            "spin Hall angle",
            "must be <= 1 in magnitude",
        ));
    }
    Ok(())
}

/// Check that mesh dimensions are valid
#[inline]
pub fn check_mesh_dimensions(nx: usize, ny: usize, nz: usize) -> crate::error::Result<()> {
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(invalid_param("mesh dimensions", "must be non-zero"));
    }
    debug_assert!(
        nx * ny * nz < 1_000_000_000,
        "Mesh too large: {}x{}x{}",
        nx,
        ny,
        nz
    );
    Ok(())
}

/// Check that a length scale is positive
#[inline]
pub fn check_length(length: f64, name: &str) -> crate::error::Result<()> {
    if length <= 0.0 {
        return Err(invalid_param(name, "must be positive"));
    }
    Ok(())
}

/// Check that spin polarization is in valid range [0, 1]
#[inline]
pub fn check_spin_polarization(p: f64) -> crate::error::Result<()> {
    if !(0.0..=1.0).contains(&p) {
        return Err(invalid_param("spin polarization", "must be in [0, 1]"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_temperature_validation() {
        assert!(check_temperature(300.0).is_ok());
        assert!(check_temperature(0.0).is_ok());
        assert!(check_temperature(-1.0).is_err());
    }

    #[test]
    fn test_magnetization_validation() {
        assert!(check_magnetization(8e5).is_ok());
        assert!(check_magnetization(0.0).is_ok());
        assert!(check_magnetization(-100.0).is_err());
    }

    #[test]
    fn test_damping_validation() {
        assert!(check_damping(0.01).is_ok());
        assert!(check_damping(0.0).is_ok());
        assert!(check_damping(0.1).is_ok());
        assert!(check_damping(-0.1).is_err());
        assert!(check_damping(1.0).is_err());
    }

    #[test]
    fn test_normalized_check() {
        let v = Vector3::new(1.0, 0.0, 0.0);
        check_normalized(v, 1e-10);

        let v2 = Vector3::new(0.6, 0.8, 0.0);
        check_normalized(v2, 1e-10);
    }

    #[test]
    fn test_time_step_validation() {
        assert!(check_time_step(1e-12).is_ok());
        assert!(check_time_step(0.0).is_err());
        assert!(check_time_step(-1e-12).is_err());
    }

    #[test]
    fn test_spin_hall_angle_validation() {
        assert!(check_spin_hall_angle(0.5).is_ok());
        assert!(check_spin_hall_angle(-0.3).is_ok());
        assert!(check_spin_hall_angle(1.0).is_ok());
        assert!(check_spin_hall_angle(1.5).is_err());
    }

    #[test]
    fn test_mesh_dimensions() {
        assert!(check_mesh_dimensions(10, 10, 1).is_ok());
        assert!(check_mesh_dimensions(0, 10, 1).is_err());
        assert!(check_mesh_dimensions(10, 0, 1).is_err());
    }

    #[test]
    fn test_spin_polarization_validation() {
        assert!(check_spin_polarization(0.5).is_ok());
        assert!(check_spin_polarization(0.0).is_ok());
        assert!(check_spin_polarization(1.0).is_ok());
        assert!(check_spin_polarization(-0.1).is_err());
        assert!(check_spin_polarization(1.5).is_err());
    }
}
