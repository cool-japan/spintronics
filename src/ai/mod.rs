//! Physical Reservoir Computing with Spintronics
//!
//! Implements neuromorphic computing using the non-linear dynamics of spin waves.
//!
//! ## Physical Reservoir Computing
//!
//! Traditional reservoir computing uses recurrent neural networks where:
//! 1. **Input Layer**: Maps inputs to reservoir
//! 2. **Reservoir Layer**: Fixed non-linear dynamics (not trained)
//! 3. **Readout Layer**: Linear combination (trained)
//!
//! In physical reservoir computing, we replace the RNN with actual physical dynamics:
//!
//! ```text
//! Input Signal → Magnetic Field → Magnon Dynamics → Spin State → Readout
//!     (encoding)                  (non-linear)      (state)     (trained weights)
//! ```
//!
//! ## Magnon-Based Implementation
//!
//! We use [`SpinChain`](crate::magnon::SpinChain) as the computing substrate:
//!
//! - **Non-linearity**: LLG equation with Gilbert damping
//! - **Memory**: Spin wave propagation and interference
//! - **High-dimensionality**: Multiple spins with 3 components each
//!
//! ## Example
//!
//! ```rust
//! use spintronics::ai::MagnonReservoir;
//! use spintronics::vector3::Vector3;
//!
//! // Create reservoir with 20 spins
//! let mut reservoir = MagnonReservoir::uniform_distribution(
//!     20,         // Total spins
//!     3,          // Input nodes
//!     5,          // Readout nodes
//!     1.0e-20,    // Exchange coupling
//!     0.01,       // Damping
//! );
//!
//! // Training data
//! let inputs = vec![
//!     vec![0.1, 0.2],
//!     vec![0.3, 0.1],
//! ];
//! let targets = vec![0.5, 0.8];
//!
//! // Train
//! let h_ext = Vector3::new(0.0, 0.0, 1.0e5);
//! reservoir.train(&inputs, &targets, h_ext, 10, 1.0e-13).unwrap();
//!
//! // Inference
//! let output = reservoir.process(&vec![0.2, 0.15], h_ext, 10, 1.0e-13);
//! ```
//!
//! ## Applications
//!
//! - Time series prediction
//! - Pattern recognition
//! - Signal classification
//! - Non-linear function approximation
//!
//! ## References
//!
//! Based on research by Prof. Eiji Saitoh's group on spin-wave-based neuromorphic computing.

pub mod reservoir;

pub use reservoir::MagnonReservoir;
