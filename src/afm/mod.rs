//! Antiferromagnetic Spintronics - THz Dynamics
//!
//! Antiferromagnets (AFMs) are materials with two (or more) magnetic sublattices
//! that are antiparallel to each other. Unlike ferromagnets, AFMs exhibit:
//!
//! - **THz-order dynamics**: Resonance frequencies in the 100 GHz - THz range
//! - **No stray fields**: Ideal for high-density spintronic devices
//! - **Ultrafast switching**: Potential for sub-ps magnetization reversal
//!
//! ## Physical Background
//!
//! The strong inter-sublattice exchange coupling (typically 100-1000 Tesla equivalent)
//! leads to resonance frequencies orders of magnitude higher than ferromagnets:
//!
//! ```text
//! f_AFM ≈ γ √(H_exchange × H_anisotropy) / (2π)
//!       ≈ 100 GHz - 1 THz
//!
//! Compare with ferromagnets:
//! f_FM ≈ γ H_anisotropy / (2π)
//!      ≈ 1-10 GHz
//! ```
//!
//! ## Applications
//!
//! - Ultrafast magnetic memory
//! - THz spintronic oscillators
//! - Antiferromagnetic tunnel junctions (AFM-TMR)
//! - Spin-orbit torque (SOT) switching of AFMs
//!
//! ## References
//!
//! Based on research from:
//! - T. Jungwirth et al., "Antiferromagnetic spintronics", Nature Nanotechnology (2016)
//! - P. Němec et al., "Antiferromagnetic opto-spintronics", Nature Physics (2018)

pub mod antiferromagnet;

pub use antiferromagnet::Antiferromagnet;
