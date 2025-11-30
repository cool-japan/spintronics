//! Stochastic Magnetization Dynamics - Finite Temperature Effects
//!
//! At finite temperatures, magnetic systems experience thermal fluctuations
//! that can drive important phenomena like thermal activation, noise-assisted
//! switching, and stochastic resonance.
//!
//! ## Fluctuation-Dissipation Theorem
//!
//! The thermal noise field H_th is related to dissipation (damping α) by:
//!
//! ```text
//! ⟨H_th,i(t) H_th,j(t')⟩ = (2 α k_B T) / (γ M_s V Δt) δ_ij δ(t-t')
//! ```
//!
//! This ensures the system reaches thermal equilibrium at temperature T.
//!
//! ## Applications
//!
//! ### 1. Thermally Activated Switching
//! - Energy barrier: ΔE = K V (anisotropy × volume)
//! - Switching rate: Γ = f₀ exp(-ΔE / k_B T)
//! - Crucial for magnetic memory retention time
//!
//! ### 2. Thermal Spin Torque Oscillators
//! - Linewidth broadening due to thermal noise
//! - Phase diffusion: Δφ ~ √(k_B T / Power)
//!
//! ### 3. Spin Seebeck Effect
//! - Temperature gradient drives magnon current
//! - Thermal fluctuations essential for transport
//!
//! ## References
//!
//! - W. F. Brown, "Thermal Fluctuations of a Single-Domain Particle",
//!   Physical Review 130, 1677 (1963)
//! - J. L. García-Palacios and F. J. Lázaro,
//!   "Langevin-dynamics study of the dynamical properties of small magnetic particles",
//!   Physical Review B 58, 14937 (1998)

pub mod thermal;

pub use thermal::{StochasticLLG, ThermalField};
