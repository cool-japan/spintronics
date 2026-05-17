//! CUDA backend skeleton (v0.9.0).
//!
//! This module declares the type signatures and trait implementation for a
//! future CUDA-accelerated LLG path. In v0.9.0 all operations return
//! [`crate::error::Error::NumericalError`] with a `"CUDA backend not yet
//! implemented"` description — this is the public contract that v1.0.0 must
//! preserve when filling in the real kernels.
//!
//! # Why a stub?
//!
//! The point of shipping the skeleton early is to:
//!   - lock in the public API surface ([`CudaDevice::new`],
//!     [`CudaDevice::with_device_id`], [`CudaDevice::device_name`]) so
//!     downstream code can compile against it today
//!   - let the [`super::available_devices`] / [`super::select_best_device`]
//!     plumbing be tested with both backends present
//!   - make it impossible to forget that `is_available()` may legitimately
//!     return `false`
//!
//! # Design intent for v1.0.0
//!
//!   - Use the [`cudarc`](https://crates.io/crates/cudarc) crate for safe
//!     CUDA bindings (chosen for its pure-Rust feel and small surface)
//!   - Persist device memory for `spins` and `h_eff` between calls (avoid
//!     PCIe round-trips per step)
//!   - One PTX kernel per primitive: `cross_product`, `rk4_step`,
//!     `zeeman_energy_reduce`
//!   - Batch all spins in a single launch (block size 256, one thread per
//!     spin) — matches the mumax³ design
//!   - Async copy-back via stream synchronisation; final
//!     `cuda_stream_synchronize` only on borrow back to host slice
//!
//! # Reference
//!
//! - A. Vansteenkiste et al., "The design and verification of MuMax3",
//!   *AIP Adv.* **4**, 107133 (2014). The LLG kernel in section 2 is the
//!   model for the v1.0.0 PTX.

use super::Device;
use crate::error::{numerical_error, Result};

/// Standard error message returned by every operation while the CUDA backend
/// is a skeleton. Centralised so v1.0.0 only has to remove one constant.
const CUDA_NOT_IMPLEMENTED: &str =
    "CUDA backend: not implemented in v0.9.0 (skeleton only — see gpu::cuda)";

/// CUDA [`Device`] handle (v0.9.0 skeleton).
///
/// In v0.9.0 the struct successfully constructs but is permanently
/// non-functional: every trait method returns
/// [`crate::error::Error::NumericalError`], and [`Self::is_available`]
/// returns `false`. The construction itself does not fail so user code can
/// pattern-match on `Result<CudaDevice, _>` and then dispatch on
/// `is_available()` — the same pattern v1.0.0 will use when distinguishing
/// "no CUDA toolkit installed" from "no compatible device found".
#[derive(Debug, Clone, Copy)]
pub struct CudaDevice {
    /// CUDA device ordinal (defaults to `0` — primary GPU).
    pub device_id: usize,
    /// Whether the device handle was successfully bound to a real CUDA
    /// runtime. **Always `false` in v0.9.0** by design.
    available: bool,
}

impl CudaDevice {
    /// Construct a CUDA device handle on the default ordinal (`0`).
    ///
    /// In v0.9.0 this never fails and always returns a non-available stub.
    /// v1.0.0 will probe for a real CUDA runtime and return
    /// [`crate::error::Error::NumericalError`] if none is available.
    pub fn new() -> Result<Self> {
        Ok(Self {
            device_id: 0,
            available: false,
        })
    }

    /// Construct a CUDA device handle on a specific device ordinal.
    ///
    /// In v0.9.0 this never fails and always returns a non-available stub
    /// regardless of `device_id`.
    pub fn with_device_id(device_id: usize) -> Result<Self> {
        Ok(Self {
            device_id,
            available: false,
        })
    }

    /// Reported device name (free-form).
    ///
    /// In the v0.9.0 skeleton this is the literal string
    /// `"CUDA stub (v0.9.0)"`. v1.0.0 will return the actual device name
    /// reported by `cuDeviceGetName`.
    pub fn device_name(&self) -> &'static str {
        "CUDA stub (v0.9.0)"
    }
}

impl Device for CudaDevice {
    fn name(&self) -> &'static str {
        "cuda"
    }

    fn is_available(&self) -> bool {
        self.available
    }

    fn step_llg_rk4(
        &self,
        _spins: &mut [[f64; 3]],
        _h_eff: &[[f64; 3]],
        _alpha: f64,
        _dt: f64,
    ) -> Result<()> {
        Err(numerical_error(CUDA_NOT_IMPLEMENTED))
    }

    fn step_llg_rk4_multi(
        &self,
        _spins: &mut [[f64; 3]],
        _h_eff: &[[f64; 3]],
        _alpha: f64,
        _dt: f64,
        _n_steps: usize,
    ) -> Result<()> {
        Err(numerical_error(CUDA_NOT_IMPLEMENTED))
    }

    fn zeeman_energy(&self, _spins: &[[f64; 3]], _h: &[[f64; 3]], _ms: f64) -> Result<f64> {
        Err(numerical_error(CUDA_NOT_IMPLEMENTED))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Error;

    #[test]
    fn test_cuda_device_construct_returns_ok_unavailable() {
        let dev = CudaDevice::new().expect("CUDA stub construction must succeed");
        assert_eq!(dev.device_id, 0);
        assert!(!dev.available);
        assert!(!dev.is_available());
    }

    #[test]
    fn test_cuda_device_with_device_id_succeeds() {
        let dev = CudaDevice::with_device_id(0).expect("with_device_id(0) must succeed");
        assert_eq!(dev.device_id, 0);
        assert!(!dev.is_available());
        // Also verify a non-zero ordinal works.
        let dev2 = CudaDevice::with_device_id(3).expect("with_device_id(3) must succeed");
        assert_eq!(dev2.device_id, 3);
        assert!(!dev2.is_available());
        // And the device_name surface.
        assert_eq!(dev2.device_name(), "CUDA stub (v0.9.0)");
    }

    #[test]
    fn test_cuda_device_name_is_cuda() {
        let dev = CudaDevice::new().unwrap();
        assert_eq!(dev.name(), "cuda");
    }

    #[test]
    fn test_cuda_is_available_false_in_v09() {
        let dev = CudaDevice::new().unwrap();
        assert!(!dev.is_available());
    }

    #[test]
    fn test_step_llg_rk4_returns_numerical_error() {
        let dev = CudaDevice::new().unwrap();
        let mut spins = vec![[1.0_f64, 0.0, 0.0]; 4];
        let h_eff = vec![[0.0_f64, 0.0, 1.0]; 4];
        let result = dev.step_llg_rk4(&mut spins, &h_eff, 0.01, 1.0e-13);
        assert!(result.is_err());
        match result.unwrap_err() {
            Error::NumericalError { description } => {
                assert!(description.contains("CUDA"));
            },
            other => panic!("expected NumericalError, got {:?}", other),
        }
        // Sanity-check the other two methods also error out.
        assert!(dev
            .step_llg_rk4_multi(&mut spins, &h_eff, 0.01, 1.0e-13, 3)
            .is_err());
        assert!(dev.zeeman_energy(&spins, &h_eff, 8.0e5).is_err());
    }
}
