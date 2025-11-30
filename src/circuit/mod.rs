//! Spin circuit elements and networks
//!
//! This module implements circuit-level modeling of spin transport:
//! - Spin resistors and conductances
//! - Spin accumulation and relaxation
//! - Multi-terminal spin circuits
//! - Compact models for circuit simulation

pub mod accumulation;
pub mod network;
pub mod resistor;

pub use accumulation::{RelaxationTime, SpinAccumulation};
pub use network::{SpinCircuit, Terminal};
pub use resistor::{SpinConductance, SpinResistor};
