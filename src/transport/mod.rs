//! Spin transport phenomena

pub mod ac_pumping;
pub mod diffusion;
pub mod pumping;

pub use ac_pumping::{AcSpinPumping, SpinBattery};
pub use diffusion::SpinDiffusion;
pub use pumping::spin_pumping_current;
