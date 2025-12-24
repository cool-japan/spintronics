# Contributing to Spintronics

Thank you for your interest in contributing to the Spintronics library! This document provides guidelines and information for contributors.

## Table of Contents

- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Code Guidelines](#code-guidelines)
- [Physics Guidelines](#physics-guidelines)
- [Testing](#testing)
- [Pull Request Process](#pull-request-process)
- [Code of Conduct](#code-of-conduct)

## Getting Started

### Prerequisites

- Rust 1.70.0 or later
- Git
- Basic knowledge of condensed matter physics (for physics-related contributions)

### Building the Project

```bash
git clone https://github.com/cool-japan/spintronics.git
cd spintronics
cargo build
cargo test
```

## Development Setup

### Running Tests

```bash
# Run all tests
cargo test

# Run tests with all features
cargo test --all-features

# Run specific test module
cargo test material::ferromagnet

# Run tests with output
cargo test -- --nocapture
```

### Linting and Formatting

We enforce code quality with clippy and rustfmt:

```bash
# Check formatting
cargo fmt --all -- --check

# Run clippy (with warnings as errors)
cargo clippy --all-features -- -D warnings

# Auto-format code
cargo fmt --all
```

### Building Documentation

```bash
cargo doc --no-deps --all-features --open
```

## Code Guidelines

### General Principles

1. **Physics First**: Always validate against physical intuition and experimental data
2. **Type Safety**: Use Rust's type system to prevent unphysical states
3. **Performance**: Profile before optimizing; correctness > speed
4. **Simplicity**: Prefer simple solutions over clever ones

### Rust Style

- Follow standard Rust naming conventions:
  - `snake_case` for functions and variables
  - `CamelCase` for types and traits
  - `SCREAMING_SNAKE_CASE` for constants
- Use meaningful variable names that reflect physics (e.g., `magnetization`, `spin_current`)
- Prefer explicit types in public APIs
- Document all public items with `///` doc comments

### API Design

- Implement `Default` for types with sensible defaults
- Implement `Display` for types that users might print
- Use builder pattern for complex configuration
- Prefer `&self` methods over consuming methods where possible
- Return `Result<T, E>` instead of panicking

### Example

```rust
/// Ferromagnetic material properties
///
/// # Example
/// ```
/// use spintronics::material::Ferromagnet;
///
/// let yig = Ferromagnet::yig();
/// println!("Damping: {}", yig.damping());
/// ```
#[derive(Debug, Clone)]
pub struct Ferromagnet {
    /// Gilbert damping constant (dimensionless)
    pub alpha: f64,
    // ... other fields
}

impl Default for Ferromagnet {
    fn default() -> Self {
        Self::yig()
    }
}
```

## Physics Guidelines

### Documentation Standards

Every physics function should include:

1. **LaTeX equations** in doc comments where applicable
2. **Physical units** clearly specified in brackets (e.g., `[A/m]`, `[J/m³]`)
3. **References** to relevant papers
4. **Physical interpretation** explaining the meaning

Example:
```rust
/// Calculate the effective field from anisotropy
///
/// The uniaxial anisotropy field is given by:
///
/// $$
/// \mathbf{H}_{\text{ani}} = \frac{2K}{\mu_0 M_s} (\mathbf{m} \cdot \hat{n}) \hat{n}
/// $$
///
/// # Arguments
/// * `magnetization` - Normalized magnetization direction
/// * `k` - Anisotropy constant [J/m³]
/// * `ms` - Saturation magnetization [A/m]
/// * `easy_axis` - Easy axis direction (normalized)
///
/// # Returns
/// Effective anisotropy field [A/m]
///
/// # References
/// - Kittel, "Introduction to Solid State Physics", Ch. 15
pub fn anisotropy_field(/* ... */) -> Vector3<f64> {
    // ...
}
```

### Physical Constants

- Use constants from `crate::constants` module
- Include NIST reference values with appropriate precision
- Document the physical meaning of each constant

### Validation

- Compare results with published experimental data
- Include validation tests that check against known results
- Document any approximations or limitations

## Testing

### Test Categories

1. **Unit Tests**: Test individual functions in isolation
2. **Integration Tests**: Test module interactions
3. **Doc Tests**: Ensure examples in documentation work
4. **Physics Validation**: Compare with experimental/theoretical results

### Writing Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_physical_quantity() {
        let material = Ferromagnet::yig();
        let result = material.some_calculation();

        // Use appropriate tolerances for floating point
        assert!((result - expected).abs() < 1e-10);
    }

    #[test]
    fn test_against_experiment() {
        // Reference: Saitoh et al., APL 88, 182509 (2006)
        let expected_voltage = 1.0e-6; // 1 μV
        let calculated = spin_pumping_voltage(/* ... */);

        // Within 10% of experimental value
        assert!((calculated - expected_voltage).abs() / expected_voltage < 0.1);
    }
}
```

## Pull Request Process

### Before Submitting

1. Ensure all tests pass: `cargo test --all-features`
2. Run clippy with no warnings: `cargo clippy --all-features -- -D warnings`
3. Format code: `cargo fmt --all`
4. Update documentation if needed
5. Add tests for new functionality

### PR Description

Include in your PR description:

1. **What**: Brief description of changes
2. **Why**: Motivation and context
3. **How**: Technical approach (if non-obvious)
4. **Testing**: How you tested the changes
5. **References**: Any relevant papers or resources

### Review Process

1. PRs require at least one approving review
2. All CI checks must pass
3. Address reviewer feedback constructively
4. Squash commits before merging (if requested)

## Code of Conduct

### Our Standards

- Be respectful and inclusive
- Focus on constructive feedback
- Accept criticism gracefully
- Prioritize the project's best interests

### Reporting Issues

Report unacceptable behavior by opening an issue or contacting the maintainers directly.

## Questions?

- Open a GitHub issue for bugs or feature requests
- Start a GitHub discussion for questions
- Check existing issues before creating new ones

Thank you for contributing to Spintronics!
