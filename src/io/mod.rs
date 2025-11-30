//! Input/Output Module
//!
//! This module provides file I/O capabilities for interoperability with
//! other simulation tools and data analysis frameworks.

pub mod ovf;

pub use ovf::{OvfData, OvfFormat, OvfReader, OvfWriter};
