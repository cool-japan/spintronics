//! Complex number type for spintronics calculations.
//!
//! A lightweight complex number implementation (`Complex`) designed for use in
//! spin-wave dispersion, Green's function matrices, Bogoliubov transformations,
//! and topological band-structure computations throughout the spintronics library.
//!
//! # Design
//!
//! - No external dependencies — pure Rust arithmetic on `f64`.
//! - All operations use explicit method calls (`a.mul(&b)`, `a.add(&b)`) to avoid
//!   confusion with real arithmetic in physics code.
//! - Constants `ZERO`, `ONE`, `I` are provided as associated constants.
//!
//! # Example
//!
//! ```rust
//! use spintronics::math::Complex;
//!
//! let a = Complex::new(1.0, 2.0);
//! let b = Complex::from_polar(1.0, std::f64::consts::FRAC_PI_4);
//! let product = a.mul(&b);
//! assert!((product.norm() - a.norm() * b.norm()).abs() < 1e-12);
//! ```

/// A complex number `re + i·im` over `f64`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    /// Real part
    pub re: f64,
    /// Imaginary part
    pub im: f64,
}

impl Complex {
    /// The zero complex number `0 + 0i`.
    pub const ZERO: Self = Self { re: 0.0, im: 0.0 };

    /// The unit complex number `1 + 0i`.
    pub const ONE: Self = Self { re: 1.0, im: 0.0 };

    /// The imaginary unit `0 + 1i`.
    pub const I: Self = Self { re: 0.0, im: 1.0 };

    /// Create a new complex number `re + i·im`.
    #[inline]
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// Create from real value `r + 0i`.
    #[inline]
    pub fn from_real(r: f64) -> Self {
        Self { re: r, im: 0.0 }
    }

    /// Create from polar form: `r·e^(iθ) = r·(cos θ + i·sin θ)`.
    #[inline]
    pub fn from_polar(r: f64, theta: f64) -> Self {
        Self {
            re: r * theta.cos(),
            im: r * theta.sin(),
        }
    }

    /// Squared modulus |z|² = re² + im².
    #[inline]
    pub fn norm_sq(&self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// Modulus |z| = √(re² + im²).
    #[inline]
    pub fn norm(&self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// Phase angle arg(z) = atan2(im, re) ∈ (-π, π].
    #[inline]
    pub fn phase(&self) -> f64 {
        self.im.atan2(self.re)
    }

    /// Complex conjugate z* = re − i·im.
    #[inline]
    pub fn conj(&self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    /// Unary negation −z.
    #[inline]
    pub fn neg(&self) -> Self {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }

    /// Addition z + w.
    #[inline]
    pub fn add(&self, other: &Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }

    /// Subtraction z − w.
    #[inline]
    pub fn sub(&self, other: &Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }

    /// Complex multiplication z·w = (ac − bd) + i(ad + bc).
    #[inline]
    pub fn mul(&self, other: &Self) -> Self {
        Self {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }

    /// Scalar multiplication z·s.
    #[inline]
    pub fn scale(&self, s: f64) -> Self {
        Self {
            re: self.re * s,
            im: self.im * s,
        }
    }

    /// Multiplication by the imaginary unit: i·z.
    #[inline]
    pub fn mul_i(&self) -> Self {
        Self {
            re: -self.im,
            im: self.re,
        }
    }

    /// Complex division z/w = z·w* / |w|².
    ///
    /// Returns `NaN + NaN·i` if `w = 0`.
    #[inline]
    pub fn div(&self, other: &Self) -> Self {
        let denom = other.norm_sq();
        let num = self.mul(&other.conj());
        Self {
            re: num.re / denom,
            im: num.im / denom,
        }
    }

    /// Complex exponential e^z = e^re · (cos(im) + i·sin(im)).
    ///
    /// Used for phase factors e^(iφ) = `Complex::from_real(0.0).exp()` for purely imaginary
    /// arguments; for full complex exponents use this method.
    #[inline]
    pub fn exp(&self) -> Self {
        let r = self.re.exp();
        Self {
            re: r * self.im.cos(),
            im: r * self.im.sin(),
        }
    }

    /// Integer power z^n.
    ///
    /// Uses repeated squaring; handles negative exponents via inversion.
    pub fn pow_n(&self, n: i32) -> Self {
        if n == 0 {
            return Self::ONE;
        }
        if n < 0 {
            return Self::ONE.div(&self.pow_n(-n));
        }
        let mut result = Self::ONE;
        let mut base = *self;
        let mut exp = n as u32;
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            exp >>= 1;
        }
        result
    }

    /// Returns `true` if either component is NaN or infinite.
    #[inline]
    pub fn is_finite(&self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
}

impl std::fmt::Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.im >= 0.0 {
            write!(f, "{:.6}+{:.6}i", self.re, self.im)
        } else {
            write!(f, "{:.6}{:.6}i", self.re, self.im)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    use super::*;

    #[test]
    fn test_constants() {
        assert_eq!(Complex::ZERO.re, 0.0);
        assert_eq!(Complex::ZERO.im, 0.0);
        assert_eq!(Complex::ONE.re, 1.0);
        assert_eq!(Complex::ONE.im, 0.0);
        assert_eq!(Complex::I.re, 0.0);
        assert_eq!(Complex::I.im, 1.0);
    }

    #[test]
    fn test_from_real() {
        let z = Complex::from_real(3.5);
        assert!((z.re - 3.5).abs() < 1e-15);
        assert!((z.im).abs() < 1e-15);
    }

    #[test]
    fn test_from_polar() {
        let z = Complex::from_polar(1.0, FRAC_PI_2);
        assert!((z.re).abs() < 1e-14);
        assert!((z.im - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_norm_sq_and_norm() {
        let z = Complex::new(3.0, 4.0);
        assert!((z.norm_sq() - 25.0).abs() < 1e-14);
        assert!((z.norm() - 5.0).abs() < 1e-14);
    }

    #[test]
    fn test_phase() {
        let z = Complex::new(1.0, 1.0);
        assert!((z.phase() - FRAC_PI_4).abs() < 1e-14);
    }

    #[test]
    fn test_conj() {
        let z = Complex::new(1.0, 2.0);
        let zc = z.conj();
        assert!((zc.re - 1.0).abs() < 1e-15);
        assert!((zc.im + 2.0).abs() < 1e-15);
    }

    #[test]
    fn test_neg() {
        let z = Complex::new(1.0, 2.0);
        let zn = z.neg();
        assert!((zn.re + 1.0).abs() < 1e-15);
        assert!((zn.im + 2.0).abs() < 1e-15);
    }

    #[test]
    fn test_mul_i_squared_is_minus_one() {
        let i = Complex::I;
        let i2 = i.mul(&i);
        assert!((i2.re + 1.0).abs() < 1e-14);
        assert!((i2.im).abs() < 1e-14);
    }

    #[test]
    fn test_div_self_is_one() {
        let z = Complex::new(3.0, 4.0);
        let r = z.div(&z);
        assert!((r.re - 1.0).abs() < 1e-14);
        assert!((r.im).abs() < 1e-14);
    }

    #[test]
    fn test_exp_purely_imaginary() {
        // e^(iπ) = -1
        let z = Complex::new(0.0, PI);
        let ez = z.exp();
        assert!((ez.re + 1.0).abs() < 1e-14);
        assert!((ez.im).abs() < 1e-14);
    }

    #[test]
    fn test_pow_n_zero() {
        let z = Complex::new(3.0, 4.0);
        let r = z.pow_n(0);
        assert!((r.re - 1.0).abs() < 1e-14);
        assert!((r.im).abs() < 1e-14);
    }

    #[test]
    fn test_pow_n_negative() {
        let z = Complex::new(2.0, 0.0);
        let r = z.pow_n(-2);
        assert!((r.re - 0.25).abs() < 1e-14);
        assert!((r.im).abs() < 1e-14);
    }

    #[test]
    fn test_is_finite() {
        assert!(Complex::new(1.0, 2.0).is_finite());
        assert!(!Complex::new(f64::NAN, 0.0).is_finite());
        assert!(!Complex::new(0.0, f64::INFINITY).is_finite());
    }
}
