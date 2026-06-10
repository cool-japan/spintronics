//! GPU acceleration foundation for spintronics (v0.9.0).
//!
//! This module provides a [`Device`] trait abstraction that lets the same LLG
//! simulation code run on CPU or GPU backends. v0.9.0 ships:
//!   - [`CpuDevice`] — always available, drop-in baseline that reuses the
//!     existing [`crate::dynamics::llg::calc_dm_dt`] scalar path
//!   - [`CudaDevice`] (feature = `cuda`) — stub that constructs successfully
//!     but reports `is_available() == false` and returns
//!     [`crate::error::Error::NumericalError`] on every operation. v1.0.0 will
//!     wire this to a real CUDA backend (planned: `cudarc` + PTX kernels)
//!     without breaking existing call sites.
//!
//! The intent is to establish the type signatures and contracts now so that
//! the v1.0.0 CUDA implementation can land as a drop-in replacement.
//!
//! # Design intent
//!
//! Operations work on flat `&mut [[f64; 3]]` arrays (per-spin x, y, z) for
//! FFI-friendliness and zero-copy compatibility with future GPU memory
//! buffers. Implementations decide whether to copy data to GPU and back per
//! call, batch operations, or maintain persistent device buffers.
//!
//! The [`Device`] trait is intentionally **object-safe** so users can store
//! `Box<dyn Device>` and switch backends at runtime via
//! [`select_best_device`].
//!
//! # Example
//!
//! ```ignore
//! use spintronics::gpu::{CpuDevice, Device};
//!
//! let device = CpuDevice::new();
//! let mut spins = vec![[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0]];
//! let h_eff = vec![[0.0_f64, 0.0, 1.0e5]; 2];
//! device.step_llg_rk4(&mut spins, &h_eff, 0.01, 1.0e-13).unwrap();
//! ```
//!
//! # References
//!
//! - A. Vansteenkiste, J. Leliaert, M. Dvornik, M. Helsen, F. Garcia-Sanchez,
//!   and B. Van Waeyenberge, "The design and verification of MuMax3",
//!   *AIP Advances* **4**, 107133 (2014).
//! - MuMAX³ Pro — GPU-accelerated micromagnetics framework (reference design
//!   for batched LLG kernels).

use crate::error::Result;

pub mod cpu;
#[cfg(feature = "cuda")]
pub mod cuda;

pub use cpu::CpuDevice;
#[cfg(feature = "cuda")]
pub use cuda::CudaDevice;

/// Device abstraction for LLG simulation on CPU or GPU backends.
///
/// All operations work on flat `&mut [[f64; 3]]` arrays so the trait can be
/// satisfied identically by CPU code (which mutates in place) and GPU code
/// (which copies to/from device memory). The trait is **object-safe**: no
/// generic methods, no `Self` in arguments other than `&self`, and no
/// `Sized` requirements — so `Box<dyn Device>` works.
pub trait Device: Send + Sync {
    /// Short device kind name, e.g. `"cpu"`, `"cuda"`.
    ///
    /// This is the discriminator used by user code that wants to dispatch on
    /// backend identity (e.g. for logging or feature gating). For richer
    /// descriptions, downcast to the concrete type.
    fn name(&self) -> &'static str;

    /// Maximum number of simultaneously-evolvable spins this device can
    /// accept in a single batch.
    ///
    /// CPU devices return [`usize::MAX`] (limited only by host RAM). GPU
    /// devices may report a lower bound based on device-memory capacity.
    fn max_spins(&self) -> usize {
        usize::MAX
    }

    /// Whether the device backend is currently functional.
    ///
    /// CPU is always `true`. CUDA is `true` only if the runtime is available
    /// and a device handle was successfully obtained. In the v0.9.0
    /// skeleton, [`CudaDevice::is_available`] always returns `false`.
    fn is_available(&self) -> bool;

    /// Evolve `spins` by one RK4 step under the given per-spin effective
    /// field. Each spin is renormalised to unit length after the step.
    ///
    /// # Arguments
    /// * `spins` — per-spin magnetisation direction, mutated in place
    /// * `h_eff` — per-spin effective field \[T\], frozen during the step
    /// * `alpha` — Gilbert damping (dimensionless)
    /// * `dt` — time step \[s\]
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] if `spins.len() !=
    /// h_eff.len()`. GPU backends may additionally return
    /// [`crate::error::Error::NumericalError`] for runtime / kernel-launch
    /// failures.
    fn step_llg_rk4(
        &self,
        spins: &mut [[f64; 3]],
        h_eff: &[[f64; 3]],
        alpha: f64,
        dt: f64,
    ) -> Result<()>;

    /// Evolve `spins` by `n_steps` RK4 steps with the same per-spin `h_eff`
    /// (frozen field — no field update between steps).
    ///
    /// Convenience wrapper for the common "relax under fixed field" pattern;
    /// GPU implementations can amortise device transfers.
    fn step_llg_rk4_multi(
        &self,
        spins: &mut [[f64; 3]],
        h_eff: &[[f64; 3]],
        alpha: f64,
        dt: f64,
        n_steps: usize,
    ) -> Result<()>;

    /// Compute the total Zeeman energy of the system:
    ///
    /// `E = -μ₀ · M_s · Σᵢ (mᵢ · hᵢ)`
    ///
    /// # Arguments
    /// * `spins` — per-spin (unit) magnetisation direction
    /// * `h` — per-spin effective field \[T\]
    /// * `ms` — saturation magnetisation \[A/m\]
    ///
    /// # Errors
    /// Returns [`crate::error::Error::DimensionMismatch`] if `spins.len() !=
    /// h.len()`.
    fn zeeman_energy(&self, spins: &[[f64; 3]], h: &[[f64; 3]], ms: f64) -> Result<f64>;
}

/// Enumerate all devices that the current build can construct.
///
/// The CPU device is always present. If the crate is built with
/// `--features cuda` and [`CudaDevice::new`] succeeds, the resulting handle
/// is appended. Note that in v0.9.0 a successfully-constructed
/// [`CudaDevice`] still reports `is_available() == false`.
pub fn available_devices() -> Vec<Box<dyn Device>> {
    #[cfg(feature = "cuda")]
    {
        let mut devices: Vec<Box<dyn Device>> = vec![Box::new(CpuDevice::new())];
        if let Ok(cuda) = CudaDevice::new() {
            devices.push(Box::new(cuda));
        }
        devices
    }
    #[cfg(not(feature = "cuda"))]
    {
        let devices: Vec<Box<dyn Device>> = vec![Box::new(CpuDevice::new())];
        devices
    }
}

/// Select the best available device.
///
/// Prefers CUDA when the runtime reports `is_available() == true`, otherwise
/// falls back to CPU. In v0.9.0 this always returns the CPU device because
/// the CUDA skeleton is non-functional by design.
pub fn select_best_device() -> Box<dyn Device> {
    #[cfg(feature = "cuda")]
    {
        if let Ok(cuda) = CudaDevice::new() {
            if cuda.is_available() {
                return Box::new(cuda);
            }
        }
    }
    Box::new(CpuDevice::new())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_device_trait_object_construction() {
        // The Device trait must be object-safe: Box<dyn Device> must compile
        // and resolve at runtime. Failure of this test indicates an
        // accidental introduction of generic methods on the trait.
        let cpu: Box<dyn Device> = Box::new(CpuDevice::new());
        assert_eq!(cpu.name(), "cpu");
    }

    #[test]
    fn test_cpu_device_new_works() {
        let cpu = CpuDevice::new();
        assert_eq!(cpu.n_threads, 0);
    }

    #[test]
    fn test_cpu_device_name() {
        let cpu = CpuDevice::new();
        assert_eq!(cpu.name(), "cpu");
    }

    #[test]
    fn test_cpu_device_is_available_true() {
        let cpu = CpuDevice::new();
        assert!(cpu.is_available());
    }

    #[test]
    fn test_available_devices_contains_at_least_one() {
        let devs = available_devices();
        assert!(!devs.is_empty());
        // The first device is always the CPU device.
        assert_eq!(devs[0].name(), "cpu");
    }

    #[test]
    fn test_select_best_device_returns_valid_device() {
        let dev = select_best_device();
        // The selected device must be functional.
        assert!(dev.is_available());
        // And one of the known backend names.
        let name = dev.name();
        assert!(name == "cpu" || name == "cuda");
    }
}
