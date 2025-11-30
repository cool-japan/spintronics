//! Error Types for Spintronics Library
//!
//! This module defines error types for the spintronics library to provide
//! better error handling and reporting than panics.

use std::fmt;

/// Result type alias for spintronics operations
pub type Result<T> = std::result::Result<T, Error>;

/// Main error type for the spintronics library
#[derive(Debug, Clone)]
pub enum Error {
    /// Invalid physical parameters
    InvalidParameter {
        /// Parameter name
        param: String,
        /// Reason for invalidity
        reason: String,
    },

    /// Dimension mismatch
    DimensionMismatch {
        /// Expected dimension
        expected: String,
        /// Actual dimension
        actual: String,
    },

    /// Numerical error (convergence, overflow, etc.)
    NumericalError {
        /// Description of the numerical issue
        description: String,
    },

    /// I/O error
    IoError {
        /// Description of the I/O error
        description: String,
    },

    /// Configuration error
    ConfigurationError {
        /// Description of the configuration issue
        description: String,
    },

    /// Feature not implemented
    NotImplemented {
        /// Feature description
        feature: String,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidParameter { param, reason } => {
                write!(f, "Invalid parameter '{}': {}", param, reason)
            }
            Error::DimensionMismatch { expected, actual } => {
                write!(
                    f,
                    "Dimension mismatch: expected {}, got {}",
                    expected, actual
                )
            }
            Error::NumericalError { description } => {
                write!(f, "Numerical error: {}", description)
            }
            Error::IoError { description } => {
                write!(f, "I/O error: {}", description)
            }
            Error::ConfigurationError { description } => {
                write!(f, "Configuration error: {}", description)
            }
            Error::NotImplemented { feature } => {
                write!(f, "Feature not implemented: {}", feature)
            }
        }
    }
}

impl std::error::Error for Error {}

impl From<std::io::Error> for Error {
    fn from(err: std::io::Error) -> Self {
        Error::IoError {
            description: err.to_string(),
        }
    }
}

impl From<std::num::ParseFloatError> for Error {
    fn from(err: std::num::ParseFloatError) -> Self {
        Error::NumericalError {
            description: format!("Parse error: {}", err),
        }
    }
}

impl From<std::num::ParseIntError> for Error {
    fn from(err: std::num::ParseIntError) -> Self {
        Error::NumericalError {
            description: format!("Parse error: {}", err),
        }
    }
}

/// Helper function to create invalid parameter error
pub fn invalid_param(param: &str, reason: &str) -> Error {
    Error::InvalidParameter {
        param: param.to_string(),
        reason: reason.to_string(),
    }
}

/// Helper function to create dimension mismatch error
pub fn dimension_mismatch(expected: &str, actual: &str) -> Error {
    Error::DimensionMismatch {
        expected: expected.to_string(),
        actual: actual.to_string(),
    }
}

/// Helper function to create numerical error
pub fn numerical_error(description: &str) -> Error {
    Error::NumericalError {
        description: description.to_string(),
    }
}

/// Helper function to create not implemented error
pub fn not_implemented(feature: &str) -> Error {
    Error::NotImplemented {
        feature: feature.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = invalid_param("temperature", "must be positive");
        assert_eq!(
            err.to_string(),
            "Invalid parameter 'temperature': must be positive"
        );

        let err = dimension_mismatch("3x3", "2x2");
        assert_eq!(err.to_string(), "Dimension mismatch: expected 3x3, got 2x2");
    }

    #[test]
    fn test_io_error_conversion() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file not found");
        let err: Error = io_err.into();
        assert!(matches!(err, Error::IoError { .. }));
    }
}
