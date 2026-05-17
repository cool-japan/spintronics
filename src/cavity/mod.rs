//! Cavity Magnonics - Hybrid Quantum Systems
//!
//! Cavity magnonics studies the coupling between magnons (spin wave quanta)
//! and photons (electromagnetic wave quanta) in microwave cavities. When the
//! coupling strength exceeds the dissipation rates, the system enters the
//! **strong coupling regime**, forming hybrid magnon-polariton states.
//!
//! ## Physical Principle
//!
//! A ferromagnetic sphere (typically YIG) placed in a microwave cavity
//! experiences:
//!
//! 1. **Ferromagnetic Resonance (FMR)**: Uniform precession mode (Kittel mode)
//!    - Frequency: ω_FMR = γ H_bias
//!    - Linewidth: Δω_m ≈ α ω_FMR (ultra-narrow for YIG: ~1 MHz)
//!
//! 2. **Cavity Photon Mode**: Standing electromagnetic wave
//!    - Frequency: ω_cavity ≈ 2-20 GHz (X-band, K-band)
//!    - Linewidth: κ ≈ 1-10 MHz
//!
//! 3. **Magnetic Dipole Coupling**: g ≈ 10-100 MHz
//!    - Interaction: oscillating cavity magnetic field ↔ precessing magnetization
//!
//! ## Strong Coupling Regime
//!
//! When cooperativity C = 4g²/(κγ_m) >> 1:
//!
//! - **Level anticrossing**: Avoided crossing in transmission spectrum
//! - **Vacuum Rabi splitting**: Δω = 2g
//! - **Coherent energy exchange**: Photon ↔ Magnon conversion
//!
//! ## Applications
//!
//! ### 1. Quantum Information Processing
//! - Magnons as quantum memory (long coherence time)
//! - Photons for information transfer
//! - Hybrid quantum transducers
//!
//! ### 2. Magnon Amplification
//! - Stimulated emission of magnons
//! - Magnon laser (analogous to photon laser)
//!
//! ### 3. Non-reciprocal Devices
//! - Isolators and circulators
//! - One-way information transfer
//!
//! ## Key Experiments
//!
//! - H. Huebl et al., "High Cooperativity in Coupled Microwave Resonator
//!   Ferrimagnetic Insulator Hybrids", Phys. Rev. Lett. 111, 127003 (2013)
//! - Y. Tabuchi et al., "Coherent coupling between a ferromagnetic magnon
//!   and a superconducting qubit", Science 349, 405 (2015)
//!
//! ## Extended Models
//!
//! Beyond the two-mode Jaynes-Cummings picture this module also contains:
//!
//! - **Tavis-Cummings / Dicke model** ([`TavisCummings`]): N spin-1/2 emitters
//!   coupled collectively to one cavity mode.  The collective coupling g_N = g·√N
//!   gives access to superradiance and large cooperativity with macroscopic spin
//!   ensembles such as YIG spheres.
//!
//! - **Magnon-polariton diagonalisation** ([`MagnonPolariton`], [`MultiModePolariton`]):
//!   analytic 2×2 Hopfield transformation and full (N+1)×(N+1) matrix diagonalisation
//!   for multi-mode photon-magnon hybridisation, Hopfield coefficients, photon/magnon
//!   fractions, and vacuum Rabi splitting.
//!
//! - **Optomagnonics** ([`BrillouinScattering`], [`OptomagnonicCoupling`],
//!   [`MicrowaveToOptical`], [`MagnonicFrequencyComb`]): Brillouin light-scattering
//!   rates and cross-sections, enhanced optomagnonic coupling via intra-cavity field,
//!   three-mode microwave-to-optical quantum transducer efficiency (Hashemi-Mahmoodian
//!   formula), and magnonic frequency-comb spectral properties.
//!
//! ## References
//!
//! Based on research by:
//! - Prof. Y. Nakamura (RIKEN/UTokyo) - Quantum magnonics
//! - Prof. G. E. W. Bauer (Tohoku Univ) - Magnon spintronics
//! - Prof. Y. Otani (RIKEN) - Spin-wave devices

pub mod hybrid;
pub mod optomagnonic;
pub mod polariton;
pub mod tavis_cummings;

pub use hybrid::HybridSystem;
pub use optomagnonic::{
    BrillouinScattering, MagnonicFrequencyComb, MicrowaveToOptical, OptomagnonicCoupling,
};
pub use polariton::{Branch, MagnonPolariton, MultiModePolariton};
pub use tavis_cummings::TavisCummings;
