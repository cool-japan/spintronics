//! CPU implementation of the [`Device`] trait.
//!
//! Wraps the existing scalar LLG path
//! ([`crate::dynamics::llg::calc_dm_dt`]) so that any `Device`-based code
//! has a baseline that always works. Performance is identical to direct use
//! of [`crate::dynamics::llg::LlgSolver`] — there is no overhead beyond the
//! trait-object dispatch (and that disappears entirely when the concrete
//! [`CpuDevice`] type is used directly).
//!
//! # RK4 implementation
//!
//! For each spin index `i`, we evaluate the standard explicit RK4 update
//! with the **field frozen** (`h_eff[i]`) over the entire step — i.e. we do
//! not refresh the effective field between sub-stages. This mirrors the
//! "constant-field" semantics that the GPU kernel will use in v1.0.0 and
//! matches the per-step contract of
//! [`crate::dynamics::llg::LlgSolver::step_rk4`] when supplied with a
//! constant `h_eff_fn`.
//!
//! The four slopes are:
//!
//! ```text
//! k1 = f(m,              h)
//! k2 = f(m + k1·dt/2,    h)
//! k3 = f(m + k2·dt/2,    h)
//! k4 = f(m + k3·dt,      h)
//! ```
//!
//! and the update is `m' = m + (k1 + 2 k2 + 2 k3 + k4)·dt/6`, followed by
//! renormalisation to unit length.

use super::Device;
use crate::constants::{GAMMA, MU_0};
use crate::dynamics::llg::calc_dm_dt;
use crate::error::{dimension_mismatch, Result};
use crate::vector3::Vector3;

/// CPU [`Device`] backend.
///
/// Carries an optional thread-count hint reserved for the v1.0.0 rayon-based
/// parallel scheduler. The hint is currently unused — the v0.9.0 CPU path is
/// strictly serial because the existing [`crate::dynamics::llg`] kernels are
/// scalar and `rayon` integration is deferred to v1.0.0.
#[derive(Debug, Clone, Copy, Default)]
pub struct CpuDevice {
    /// Hint for the number of worker threads. `0` means "let the runtime
    /// pick a default". Currently unused; reserved for v1.0.0 rayon
    /// scheduling. Public so users can inspect it.
    pub n_threads: usize,
}

impl CpuDevice {
    /// Construct a CPU device with default settings (`n_threads = 0`,
    /// meaning "runtime default").
    pub fn new() -> Self {
        Self { n_threads: 0 }
    }

    /// Construct a CPU device with an explicit thread-count hint.
    ///
    /// The value is stored for future use; v0.9.0 ignores it because the
    /// scalar LLG path is strictly serial.
    pub fn with_threads(n_threads: usize) -> Self {
        Self { n_threads }
    }
}

/// Convert `[f64; 3]` to [`Vector3<f64>`].
#[inline]
fn to_vec3(a: [f64; 3]) -> Vector3<f64> {
    Vector3::new(a[0], a[1], a[2])
}

/// Convert [`Vector3<f64>`] to `[f64; 3]`.
#[inline]
fn from_vec3(v: Vector3<f64>) -> [f64; 3] {
    [v.x, v.y, v.z]
}

/// Renormalise a magnetisation vector to unit length.
///
/// If the magnitude is below `1e-30` we leave the vector untouched (this
/// only happens in pathological cases where the simulation has already lost
/// the spin direction).
#[inline]
fn renormalize(v: Vector3<f64>) -> Vector3<f64> {
    let mag = v.magnitude();
    if mag > 1.0e-30 {
        Vector3::new(v.x / mag, v.y / mag, v.z / mag)
    } else {
        v
    }
}

/// One explicit RK4 step for a single spin under a frozen field.
#[inline]
fn rk4_single(m: Vector3<f64>, h: Vector3<f64>, alpha: f64, dt: f64) -> Vector3<f64> {
    let k1 = calc_dm_dt(m, h, GAMMA, alpha);
    let k2 = calc_dm_dt(m + k1 * (dt * 0.5), h, GAMMA, alpha);
    let k3 = calc_dm_dt(m + k2 * (dt * 0.5), h, GAMMA, alpha);
    let k4 = calc_dm_dt(m + k3 * dt, h, GAMMA, alpha);
    let m_new = m + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
    renormalize(m_new)
}

impl Device for CpuDevice {
    fn name(&self) -> &'static str {
        "cpu"
    }

    fn is_available(&self) -> bool {
        true
    }

    fn step_llg_rk4(
        &self,
        spins: &mut [[f64; 3]],
        h_eff: &[[f64; 3]],
        alpha: f64,
        dt: f64,
    ) -> Result<()> {
        if spins.len() != h_eff.len() {
            return Err(dimension_mismatch(
                &format!("h_eff.len() = spins.len() = {}", spins.len()),
                &format!("h_eff.len() = {}", h_eff.len()),
            ));
        }
        for (spin, h) in spins.iter_mut().zip(h_eff.iter()) {
            let m = to_vec3(*spin);
            let h_vec = to_vec3(*h);
            *spin = from_vec3(rk4_single(m, h_vec, alpha, dt));
        }
        Ok(())
    }

    fn step_llg_rk4_multi(
        &self,
        spins: &mut [[f64; 3]],
        h_eff: &[[f64; 3]],
        alpha: f64,
        dt: f64,
        n_steps: usize,
    ) -> Result<()> {
        if spins.len() != h_eff.len() {
            return Err(dimension_mismatch(
                &format!("h_eff.len() = spins.len() = {}", spins.len()),
                &format!("h_eff.len() = {}", h_eff.len()),
            ));
        }
        for _ in 0..n_steps {
            for (spin, h) in spins.iter_mut().zip(h_eff.iter()) {
                let m = to_vec3(*spin);
                let h_vec = to_vec3(*h);
                *spin = from_vec3(rk4_single(m, h_vec, alpha, dt));
            }
        }
        Ok(())
    }

    fn zeeman_energy(&self, spins: &[[f64; 3]], h: &[[f64; 3]], ms: f64) -> Result<f64> {
        if spins.len() != h.len() {
            return Err(dimension_mismatch(
                &format!("h.len() = spins.len() = {}", spins.len()),
                &format!("h.len() = {}", h.len()),
            ));
        }
        // E = -μ₀ · M_s · Σᵢ (mᵢ · hᵢ)
        let mut acc = 0.0_f64;
        for (m, hi) in spins.iter().zip(h.iter()) {
            acc += m[0] * hi[0] + m[1] * hi[1] + m[2] * hi[2];
        }
        Ok(-MU_0 * ms * acc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpu_device_new_and_with_threads() {
        let a = CpuDevice::new();
        assert_eq!(a.n_threads, 0);
        let b = CpuDevice::with_threads(8);
        assert_eq!(b.n_threads, 8);
        let c = CpuDevice::default();
        assert_eq!(c.n_threads, 0);
    }

    #[test]
    fn test_cpu_device_name_is_cpu() {
        let dev = CpuDevice::new();
        assert_eq!(dev.name(), "cpu");
    }

    #[test]
    fn test_cpu_device_is_available_true() {
        let dev = CpuDevice::new();
        assert!(dev.is_available());
        // max_spins default is usize::MAX
        assert_eq!(dev.max_spins(), usize::MAX);
    }

    #[test]
    fn test_step_llg_rk4_size_mismatch_returns_error() {
        let dev = CpuDevice::new();
        let mut spins = vec![[1.0, 0.0, 0.0]; 3];
        let h_eff = vec![[0.0, 0.0, 1.0]; 2];
        let result = dev.step_llg_rk4(&mut spins, &h_eff, 0.01, 1.0e-13);
        assert!(result.is_err());
    }

    #[test]
    fn test_step_llg_rk4_multi_size_mismatch_returns_error() {
        let dev = CpuDevice::new();
        let mut spins = vec![[1.0, 0.0, 0.0]; 3];
        let h_eff = vec![[0.0, 0.0, 1.0]; 2];
        let result = dev.step_llg_rk4_multi(&mut spins, &h_eff, 0.01, 1.0e-13, 5);
        assert!(result.is_err());
    }

    #[test]
    fn test_single_spin_larmor_preserves_norm() {
        // A single spin under a constant field should precess (Larmor) and
        // remain unit-length after RK4 + renormalisation.
        let dev = CpuDevice::new();
        let mut spins = vec![[1.0_f64, 0.0, 0.0]];
        // Sub-Tesla field; gamma ~ 1.76e11 rad/(s·T); choose dt so
        // gamma*B*dt ~ 1e-2 (small precession angle per step).
        let b = 1.0_f64; // 1 T
        let h_eff = vec![[0.0, 0.0, b]];
        let dt = 5.0e-14_f64;
        dev.step_llg_rk4(&mut spins, &h_eff, 0.01, dt).unwrap();
        let mag =
            (spins[0][0] * spins[0][0] + spins[0][1] * spins[0][1] + spins[0][2] * spins[0][2])
                .sqrt();
        assert!((mag - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_multi_step_runs_and_norm_preserved() {
        // 100 RK4 steps with a damped LLG: should converge spin towards the
        // field direction (+z) while preserving |m| = 1 throughout.
        let dev = CpuDevice::new();
        let mut spins = vec![[0.7_f64, 0.0, 0.7]];
        // Normalize input
        let n = (spins[0][0] * spins[0][0] + spins[0][1] * spins[0][1] + spins[0][2] * spins[0][2])
            .sqrt();
        spins[0][0] /= n;
        spins[0][1] /= n;
        spins[0][2] /= n;
        let h_eff = vec![[0.0_f64, 0.0, 1.0]];
        let alpha = 0.5_f64; // strong damping → rapid alignment
        let dt = 1.0e-13_f64;
        dev.step_llg_rk4_multi(&mut spins, &h_eff, alpha, dt, 2000)
            .unwrap();
        let mag =
            (spins[0][0] * spins[0][0] + spins[0][1] * spins[0][1] + spins[0][2] * spins[0][2])
                .sqrt();
        assert!((mag - 1.0).abs() < 1.0e-8);
        // m_z should have grown towards +1 due to damping; very loose check.
        assert!(spins[0][2] > 0.7);
    }

    #[test]
    fn test_zeeman_energy_aligned_is_negative() {
        // m and h aligned → energy = -μ₀ · ms · |h|·|m| (most-negative).
        let dev = CpuDevice::new();
        let spins = vec![[0.0, 0.0, 1.0]];
        let h = vec![[0.0, 0.0, 1.0]];
        let ms = 8.0e5_f64; // ~ permalloy
        let e = dev.zeeman_energy(&spins, &h, ms).unwrap();
        let expected = -MU_0 * ms * 1.0_f64;
        assert!((e - expected).abs() < 1.0e-20);
        assert!(e < 0.0);
    }

    #[test]
    fn test_zeeman_energy_orthogonal_is_zero() {
        let dev = CpuDevice::new();
        let spins = vec![[1.0, 0.0, 0.0]];
        let h = vec![[0.0, 0.0, 1.0]];
        let ms = 8.0e5_f64;
        let e = dev.zeeman_energy(&spins, &h, ms).unwrap();
        assert!(e.abs() < 1.0e-25);
    }

    #[test]
    fn test_zeeman_energy_size_mismatch_returns_error() {
        let dev = CpuDevice::new();
        let spins = vec![[1.0, 0.0, 0.0]; 5];
        let h = vec![[0.0, 0.0, 1.0]; 4];
        let ms = 8.0e5_f64;
        let result = dev.zeeman_energy(&spins, &h, ms);
        assert!(result.is_err());
    }

    #[test]
    fn test_large_batch_runs() {
        // 1000-spin batch — verify the loop scales and norms stay preserved.
        let dev = CpuDevice::new();
        let n = 1000;
        let mut spins = vec![[1.0_f64, 0.0, 0.0]; n];
        let h_eff = vec![[0.0_f64, 0.0, 1.0]; n];
        dev.step_llg_rk4_multi(&mut spins, &h_eff, 0.01, 1.0e-13, 10)
            .unwrap();
        for spin in &spins {
            let mag = (spin[0] * spin[0] + spin[1] * spin[1] + spin[2] * spin[2]).sqrt();
            assert!((mag - 1.0).abs() < 1.0e-9);
        }
    }
}
