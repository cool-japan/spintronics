//! SIMD-friendly batch processing for Vector3 arrays
//!
//! This module provides a Structure of Arrays (SoA) layout for Vector3 data,
//! enabling auto-vectorization by the compiler on stable Rust without nightly
//! intrinsics. By storing x, y, z components in separate contiguous arrays,
//! the compiler can apply SIMD instructions (SSE, AVX, etc.) to process
//! multiple vectors simultaneously.
//!
//! # Design
//!
//! Traditional Array of Structures (AoS) layout `[{x,y,z}, {x,y,z}, ...]`
//! leads to strided memory accesses that defeat SIMD. The SoA layout used
//! here `{[x,x,...], [y,y,...], [z,z,...]}` enables contiguous loads for
//! each component, which is ideal for vectorized operations.
//!
//! The `#[repr(align(32))]` attribute ensures 32-byte alignment, matching
//! AVX register width for optimal performance.

use crate::error::{self, Result};
use crate::vector3::Vector3;

/// SIMD-friendly batch processing for Vector3 arrays using Structure of Arrays layout.
///
/// By separating x, y, z components into individual contiguous arrays, the compiler
/// can auto-vectorize element-wise operations using SIMD instructions on stable Rust.
///
/// # Example
///
/// ```
/// use spintronics::simd::SimdBatch;
/// use spintronics::Vector3;
///
/// let spins = vec![
///     Vector3::new(1.0, 0.0, 0.0),
///     Vector3::new(0.0, 1.0, 0.0),
/// ];
/// let batch = SimdBatch::from_vector3_slice(&spins);
/// assert_eq!(batch.len(), 2);
/// ```
#[repr(align(32))]
pub struct SimdBatch {
    /// X components
    pub x: Vec<f64>,
    /// Y components
    pub y: Vec<f64>,
    /// Z components
    pub z: Vec<f64>,
}

impl SimdBatch {
    /// Create a new zeroed batch of the given length.
    ///
    /// # Arguments
    /// * `len` - Number of vectors in the batch
    #[inline]
    pub fn new(len: usize) -> Self {
        Self {
            x: vec![0.0; len],
            y: vec![0.0; len],
            z: vec![0.0; len],
        }
    }

    /// Return the number of vectors in the batch.
    #[inline]
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Return whether the batch is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Convert a slice of `Vector3<f64>` into a `SimdBatch` (AoS -> SoA).
    ///
    /// # Arguments
    /// * `spins` - Slice of 3D vectors to convert
    #[inline]
    pub fn from_vector3_slice(spins: &[Vector3<f64>]) -> Self {
        let len = spins.len();
        let mut batch = Self::new(len);
        for (i, v) in spins.iter().enumerate() {
            batch.x[i] = v.x;
            batch.y[i] = v.y;
            batch.z[i] = v.z;
        }
        batch
    }

    /// Convert the batch back to a `Vec<Vector3<f64>>` (SoA -> AoS).
    #[inline]
    pub fn to_vector3_vec(&self) -> Vec<Vector3<f64>> {
        let len = self.len();
        let mut result = Vec::with_capacity(len);
        for i in 0..len {
            result.push(Vector3::new(self.x[i], self.y[i], self.z[i]));
        }
        result
    }
}

/// Check that two batches have the same length, returning a `DimensionMismatch` error
/// if they differ.
#[inline]
fn check_length(a: &SimdBatch, b: &SimdBatch) -> Result<()> {
    if a.len() != b.len() {
        return Err(error::dimension_mismatch(
            &format!("{}", a.len()),
            &format!("{}", b.len()),
        ));
    }
    Ok(())
}

/// Element-wise addition of two batches.
///
/// Returns a new `SimdBatch` where each component is the sum of the
/// corresponding components from `a` and `b`.
///
/// # Errors
/// Returns `DimensionMismatch` if the batches have different lengths.
#[inline]
pub fn batch_add(a: &SimdBatch, b: &SimdBatch) -> Result<SimdBatch> {
    check_length(a, b)?;
    let len = a.len();
    let mut result = SimdBatch::new(len);
    for i in 0..len {
        result.x[i] = a.x[i] + b.x[i];
        result.y[i] = a.y[i] + b.y[i];
        result.z[i] = a.z[i] + b.z[i];
    }
    Ok(result)
}

/// Element-wise subtraction of two batches (a - b).
///
/// Returns a new `SimdBatch` where each component is the difference of the
/// corresponding components from `a` and `b`.
///
/// # Errors
/// Returns `DimensionMismatch` if the batches have different lengths.
#[inline]
pub fn batch_sub(a: &SimdBatch, b: &SimdBatch) -> Result<SimdBatch> {
    check_length(a, b)?;
    let len = a.len();
    let mut result = SimdBatch::new(len);
    for i in 0..len {
        result.x[i] = a.x[i] - b.x[i];
        result.y[i] = a.y[i] - b.y[i];
        result.z[i] = a.z[i] - b.z[i];
    }
    Ok(result)
}

/// Batch cross product of two vector batches.
///
/// Computes `a\[i\] x b\[i\]` for each index `i`, using the standard cross
/// product formula decomposed into component-wise operations for
/// auto-vectorization.
///
/// # Errors
/// Returns `DimensionMismatch` if the batches have different lengths.
#[inline]
pub fn batch_cross(a: &SimdBatch, b: &SimdBatch) -> Result<SimdBatch> {
    check_length(a, b)?;
    let len = a.len();
    let mut result = SimdBatch::new(len);
    for i in 0..len {
        result.x[i] = a.y[i] * b.z[i] - a.z[i] * b.y[i];
        result.y[i] = a.z[i] * b.x[i] - a.x[i] * b.z[i];
        result.z[i] = a.x[i] * b.y[i] - a.y[i] * b.x[i];
    }
    Ok(result)
}

/// Batch dot product of two vector batches.
///
/// Returns a `Vec<f64>` where each element is the dot product of the
/// corresponding vectors from `a` and `b`.
///
/// # Errors
/// Returns `DimensionMismatch` if the batches have different lengths.
#[inline]
pub fn batch_dot(a: &SimdBatch, b: &SimdBatch) -> Result<Vec<f64>> {
    check_length(a, b)?;
    let len = a.len();
    let mut result = vec![0.0; len];
    for (i, val) in result.iter_mut().enumerate().take(len) {
        *val = a.x[i] * b.x[i] + a.y[i] * b.y[i] + a.z[i] * b.z[i];
    }
    Ok(result)
}

/// Scale all vectors in a batch by a scalar factor.
///
/// Returns a new `SimdBatch` where every component is multiplied by `s`.
#[inline]
pub fn batch_scale(a: &SimdBatch, s: f64) -> SimdBatch {
    let len = a.len();
    let mut result = SimdBatch::new(len);
    for i in 0..len {
        result.x[i] = a.x[i] * s;
        result.y[i] = a.y[i] * s;
        result.z[i] = a.z[i] * s;
    }
    result
}

/// Normalize each vector in the batch in-place so that its magnitude is 1.
///
/// Vectors with zero magnitude are left as zero vectors (no division by zero).
#[inline]
pub fn batch_normalize(a: &mut SimdBatch) {
    let len = a.len();
    for i in 0..len {
        let mag_sq = a.x[i] * a.x[i] + a.y[i] * a.y[i] + a.z[i] * a.z[i];
        if mag_sq > 0.0 {
            let inv_mag = 1.0 / mag_sq.sqrt();
            a.x[i] *= inv_mag;
            a.y[i] *= inv_mag;
            a.z[i] *= inv_mag;
        }
    }
}

/// SIMD-friendly batch computation of the LLG equation dm/dt.
///
/// Computes the Landau-Lifshitz-Gilbert equation for an entire array of
/// magnetization vectors simultaneously, enabling auto-vectorization:
///
/// ```text
/// dm/dt = (-gamma / (1 + alpha^2)) * [m x H + alpha * m x (m x H)]
/// ```
///
/// This is mathematically equivalent to `calc_dm_dt` from
/// `crate::dynamics::llg`, but operates on batches for better throughput.
///
/// # Arguments
/// * `m` - Batch of normalized magnetization vectors
/// * `h` - Batch of effective magnetic field vectors \[T\]
/// * `gamma` - Gyromagnetic ratio [rad/(s*T)]
/// * `alpha` - Gilbert damping constant (dimensionless)
///
/// # Errors
/// Returns `DimensionMismatch` if `m` and `h` have different lengths.
///
/// # Example
///
/// ```
/// use spintronics::simd::{SimdBatch, batch_calc_dm_dt};
/// use spintronics::constants::GAMMA;
///
/// let mut m = SimdBatch::new(2);
/// m.x[0] = 1.0; m.x[1] = 0.0;
/// m.y[0] = 0.0; m.y[1] = 1.0;
/// m.z[0] = 0.0; m.z[1] = 0.0;
///
/// let mut h = SimdBatch::new(2);
/// h.z[0] = 1.0; h.z[1] = 1.0;
///
/// let dm_dt = batch_calc_dm_dt(&m, &h, GAMMA, 0.01).expect("dimension match");
/// ```
pub fn batch_calc_dm_dt(m: &SimdBatch, h: &SimdBatch, gamma: f64, alpha: f64) -> Result<SimdBatch> {
    check_length(m, h)?;
    let len = m.len();
    let prefactor = -gamma / (1.0 + alpha * alpha);
    let mut result = SimdBatch::new(len);

    for i in 0..len {
        // Precession: m x H
        let px = m.y[i] * h.z[i] - m.z[i] * h.y[i];
        let py = m.z[i] * h.x[i] - m.x[i] * h.z[i];
        let pz = m.x[i] * h.y[i] - m.y[i] * h.x[i];

        // Double cross product: m x (m x H)
        let dx = m.y[i] * pz - m.z[i] * py;
        let dy = m.z[i] * px - m.x[i] * pz;
        let dz = m.x[i] * py - m.y[i] * px;

        // dm/dt = prefactor * [precession + alpha * damping]
        result.x[i] = prefactor * (px + alpha * dx);
        result.y[i] = prefactor * (py + alpha * dy);
        result.z[i] = prefactor * (pz + alpha * dz);
    }

    Ok(result)
}

/// Element-wise fused multiply-add: result\[i\] = a\[i\] + s * b\[i\].
///
/// Each component of the result is the corresponding component of `a` plus
/// `s` times the corresponding component of `b`. This fused operation avoids
/// an intermediate allocation compared to `batch_add(&a, &batch_scale(&b, s))`.
///
/// # Arguments
/// * `a` - Base batch
/// * `b` - Batch to scale and add
/// * `s` - Scalar multiplier applied to `b`
///
/// # Errors
/// Returns `DimensionMismatch` if the batches have different lengths.
///
/// # Example
///
/// ```
/// use spintronics::simd::{SimdBatch, batch_add_scaled};
///
/// let mut a = SimdBatch::new(2);
/// a.x[0] = 1.0; a.y[0] = 0.0; a.z[0] = 0.0;
/// a.x[1] = 0.0; a.y[1] = 1.0; a.z[1] = 0.0;
///
/// let mut b = SimdBatch::new(2);
/// b.x[0] = 1.0; b.x[1] = 1.0;
///
/// // result[i] = a[i] + 0.5 * b[i]
/// let result = batch_add_scaled(&a, &b, 0.5).expect("same length");
/// assert!((result.x[0] - 1.5).abs() < 1e-12);
/// ```
#[inline]
pub fn batch_add_scaled(a: &SimdBatch, b: &SimdBatch, s: f64) -> Result<SimdBatch> {
    check_length(a, b)?;
    let len = a.len();
    let mut result = SimdBatch::new(len);
    for i in 0..len {
        result.x[i] = a.x[i] + s * b.x[i];
        result.y[i] = a.y[i] + s * b.y[i];
        result.z[i] = a.z[i] + s * b.z[i];
    }
    Ok(result)
}

/// RK4 step for N independent spins in parallel using SoA SIMD layout.
///
/// Each spin evolves under its own `h_eff[i]` with no inter-spin coupling.
/// The standard 4th-order Runge-Kutta scheme is applied to the LLG equation:
///
/// ```text
/// k1 = dt * dm/dt(m,          h_eff)
/// k2 = dt * dm/dt(m + 0.5*k1, h_eff)
/// k3 = dt * dm/dt(m + 0.5*k2, h_eff)
/// k4 = dt * dm/dt(m + k3,     h_eff)
/// m_new = m + (k1 + 2*k2 + 2*k3 + k4) / 6
/// ```
///
/// The result is normalized so that `|m_i| = 1` for every spin.
///
/// # Arguments
/// * `m`     - Batch of normalized magnetization vectors (input)
/// * `h_eff` - Batch of effective magnetic field vectors \[T\]
/// * `alpha` - Gilbert damping constant (dimensionless)
/// * `gamma` - Gyromagnetic ratio \[rad/(s·T)\]
/// * `dt`    - Time step \[s\]
///
/// # Errors
/// Returns `DimensionMismatch` if `m` and `h_eff` have different lengths.
///
/// # Example
///
/// ```
/// use spintronics::simd::{SimdBatch, batch_evolve_rk4};
/// use spintronics::constants::GAMMA;
///
/// let mut m = SimdBatch::new(1);
/// m.x[0] = 1.0; // spin along x
///
/// let mut h = SimdBatch::new(1);
/// h.z[0] = 1.0; // field along z
///
/// let m_new = batch_evolve_rk4(&m, &h, 0.01, GAMMA, 1.0e-12)
///     .expect("same length");
/// let mag = (m_new.x[0].powi(2) + m_new.y[0].powi(2) + m_new.z[0].powi(2)).sqrt();
/// assert!((mag - 1.0).abs() < 1e-10);
/// ```
pub fn batch_evolve_rk4(
    m: &SimdBatch,
    h_eff: &SimdBatch,
    alpha: f64,
    gamma: f64,
    dt: f64,
) -> Result<SimdBatch> {
    check_length(m, h_eff)?;

    // k1 = dm/dt(m, h_eff)
    let k1_raw = batch_calc_dm_dt(m, h_eff, gamma, alpha)?;
    let k1 = batch_scale(&k1_raw, dt);

    // k2 = dm/dt(m + 0.5*k1, h_eff)
    let m2 = batch_add_scaled(m, &k1, 0.5)?;
    let k2_raw = batch_calc_dm_dt(&m2, h_eff, gamma, alpha)?;
    let k2 = batch_scale(&k2_raw, dt);

    // k3 = dm/dt(m + 0.5*k2, h_eff)
    let m3 = batch_add_scaled(m, &k2, 0.5)?;
    let k3_raw = batch_calc_dm_dt(&m3, h_eff, gamma, alpha)?;
    let k3 = batch_scale(&k3_raw, dt);

    // k4 = dm/dt(m + k3, h_eff)
    let m4 = batch_add(m, &k3)?;
    let k4_raw = batch_calc_dm_dt(&m4, h_eff, gamma, alpha)?;
    let k4 = batch_scale(&k4_raw, dt);

    // m_new = m + (k1 + 2*k2 + 2*k3 + k4) / 6
    let len = m.len();
    let mut m_new = SimdBatch::new(len);
    for i in 0..len {
        m_new.x[i] = m.x[i] + (k1.x[i] + 2.0 * k2.x[i] + 2.0 * k3.x[i] + k4.x[i]) / 6.0;
        m_new.y[i] = m.y[i] + (k1.y[i] + 2.0 * k2.y[i] + 2.0 * k3.y[i] + k4.y[i]) / 6.0;
        m_new.z[i] = m.z[i] + (k1.z[i] + 2.0 * k2.z[i] + 2.0 * k3.z[i] + k4.z[i]) / 6.0;
    }

    batch_normalize(&mut m_new);
    Ok(m_new)
}

/// Multi-step evolution: run `num_steps` RK4 iterations on a batch.
///
/// `h_eff` is held constant throughout (static field approximation).
/// This avoids repeated allocation of intermediate batches and is preferred
/// over calling `batch_evolve_rk4` in a loop when the field does not change.
///
/// # Arguments
/// * `m`         - Initial batch of normalized magnetization vectors (consumed)
/// * `h_eff`     - Batch of effective magnetic field vectors \[T\] (constant)
/// * `alpha`     - Gilbert damping constant (dimensionless)
/// * `gamma`     - Gyromagnetic ratio \[rad/(s·T)\]
/// * `dt`        - Time step \[s\]
/// * `num_steps` - Number of RK4 iterations to perform
///
/// # Errors
/// Returns `DimensionMismatch` if `m` and `h_eff` have different lengths.
///
/// # Example
///
/// ```
/// use spintronics::simd::{SimdBatch, batch_evolve_multi_step};
/// use spintronics::constants::GAMMA;
///
/// let mut m = SimdBatch::new(4);
/// for i in 0..4 {
///     m.x[i] = 1.0; // all spins along x
/// }
/// let mut h = SimdBatch::new(4);
/// for i in 0..4 {
///     h.z[i] = 1.0; // field along z for all
/// }
///
/// let m_final = batch_evolve_multi_step(m, &h, 0.01, GAMMA, 1.0e-13, 10)
///     .expect("same length");
/// assert_eq!(m_final.len(), 4);
/// ```
pub fn batch_evolve_multi_step(
    mut m: SimdBatch,
    h_eff: &SimdBatch,
    alpha: f64,
    gamma: f64,
    dt: f64,
    num_steps: usize,
) -> Result<SimdBatch> {
    check_length(&m, h_eff)?;
    for _ in 0..num_steps {
        m = batch_evolve_rk4(&m, h_eff, alpha, gamma, dt)?;
    }
    Ok(m)
}

/// Convert a slice of `Vector3<f64>` into a `SimdBatch` (convenience function).
///
/// This is equivalent to `SimdBatch::from_vector3_slice`.
#[inline]
pub fn from_vector3_slice(spins: &[Vector3<f64>]) -> SimdBatch {
    SimdBatch::from_vector3_slice(spins)
}

/// Convert a `SimdBatch` back to a `Vec<Vector3<f64>>` (convenience function).
///
/// This is equivalent to `SimdBatch::to_vector3_vec`.
#[inline]
pub fn to_vector3_vec(batch: &SimdBatch) -> Vec<Vector3<f64>> {
    batch.to_vector3_vec()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GAMMA;
    use crate::dynamics::llg::calc_dm_dt;

    const TOL: f64 = 1e-12;

    #[test]
    fn test_batch_add_matches_vector3() {
        let a_vecs = vec![
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(-1.0, 0.5, 4.0),
            Vector3::new(0.0, 0.0, 0.0),
        ];
        let b_vecs = vec![
            Vector3::new(4.0, -1.0, 0.5),
            Vector3::new(2.0, 3.0, -1.0),
            Vector3::new(1.0, 1.0, 1.0),
        ];

        let a_batch = SimdBatch::from_vector3_slice(&a_vecs);
        let b_batch = SimdBatch::from_vector3_slice(&b_vecs);
        let result = batch_add(&a_batch, &b_batch).expect("same length batches");
        let result_vecs = result.to_vector3_vec();

        for i in 0..a_vecs.len() {
            let expected = a_vecs[i] + b_vecs[i];
            assert!(
                (result_vecs[i].x - expected.x).abs() < TOL,
                "x mismatch at index {}",
                i
            );
            assert!(
                (result_vecs[i].y - expected.y).abs() < TOL,
                "y mismatch at index {}",
                i
            );
            assert!(
                (result_vecs[i].z - expected.z).abs() < TOL,
                "z mismatch at index {}",
                i
            );
        }
    }

    #[test]
    fn test_batch_cross_matches_vector3() {
        let a_vecs = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 2.0, 3.0),
        ];
        let b_vecs = vec![
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(4.0, 5.0, 6.0),
        ];

        let a_batch = SimdBatch::from_vector3_slice(&a_vecs);
        let b_batch = SimdBatch::from_vector3_slice(&b_vecs);
        let result = batch_cross(&a_batch, &b_batch).expect("same length batches");
        let result_vecs = result.to_vector3_vec();

        for i in 0..a_vecs.len() {
            let expected = a_vecs[i].cross(&b_vecs[i]);
            assert!(
                (result_vecs[i].x - expected.x).abs() < TOL,
                "x mismatch at index {}",
                i
            );
            assert!(
                (result_vecs[i].y - expected.y).abs() < TOL,
                "y mismatch at index {}",
                i
            );
            assert!(
                (result_vecs[i].z - expected.z).abs() < TOL,
                "z mismatch at index {}",
                i
            );
        }
    }

    #[test]
    fn test_batch_calc_dm_dt_matches_scalar() {
        let alpha = 0.01;

        let m_vecs = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(
                1.0 / 3.0_f64.sqrt(),
                1.0 / 3.0_f64.sqrt(),
                1.0 / 3.0_f64.sqrt(),
            ),
        ];
        let h_vecs = vec![
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];

        let m_batch = SimdBatch::from_vector3_slice(&m_vecs);
        let h_batch = SimdBatch::from_vector3_slice(&h_vecs);
        let result =
            batch_calc_dm_dt(&m_batch, &h_batch, GAMMA, alpha).expect("same length batches");
        let result_vecs = result.to_vector3_vec();

        for i in 0..m_vecs.len() {
            let expected = calc_dm_dt(m_vecs[i], h_vecs[i], GAMMA, alpha);
            assert!(
                (result_vecs[i].x - expected.x).abs() < TOL,
                "x mismatch at index {}: got {} expected {}",
                i,
                result_vecs[i].x,
                expected.x
            );
            assert!(
                (result_vecs[i].y - expected.y).abs() < TOL,
                "y mismatch at index {}: got {} expected {}",
                i,
                result_vecs[i].y,
                expected.y
            );
            assert!(
                (result_vecs[i].z - expected.z).abs() < TOL,
                "z mismatch at index {}: got {} expected {}",
                i,
                result_vecs[i].z,
                expected.z
            );
        }
    }

    #[test]
    fn test_roundtrip_vector3_simd_vector3() {
        let original = vec![
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(-0.5, 0.7, -1.3),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1e10, -1e-10, 42.0),
        ];

        let batch = SimdBatch::from_vector3_slice(&original);
        let recovered = batch.to_vector3_vec();

        assert_eq!(original.len(), recovered.len());
        for i in 0..original.len() {
            assert!(
                (recovered[i].x - original[i].x).abs() < TOL,
                "x mismatch at index {}",
                i
            );
            assert!(
                (recovered[i].y - original[i].y).abs() < TOL,
                "y mismatch at index {}",
                i
            );
            assert!(
                (recovered[i].z - original[i].z).abs() < TOL,
                "z mismatch at index {}",
                i
            );
        }
    }

    #[test]
    fn test_dimension_mismatch_error() {
        let a = SimdBatch::new(3);
        let b = SimdBatch::new(5);
        let err = batch_add(&a, &b);
        assert!(err.is_err());
        let err = batch_cross(&a, &b);
        assert!(err.is_err());
        let err = batch_dot(&a, &b);
        assert!(err.is_err());
    }

    // -----------------------------------------------------------------------
    // New tests: batch_add_scaled, batch_evolve_rk4, batch_evolve_multi_step
    // -----------------------------------------------------------------------

    #[test]
    fn test_batch_add_scaled_correctness() {
        // result[i] = a[i] + 3.0 * b[i], verified component by component
        let a_vecs = vec![Vector3::new(1.0, 2.0, 3.0), Vector3::new(-1.0, 0.5, 4.0)];
        let b_vecs = vec![Vector3::new(2.0, -1.0, 0.5), Vector3::new(3.0, 1.0, -2.0)];
        let s = 3.0_f64;

        let a = SimdBatch::from_vector3_slice(&a_vecs);
        let b = SimdBatch::from_vector3_slice(&b_vecs);
        let result = batch_add_scaled(&a, &b, s).expect("same length");

        for i in 0..a_vecs.len() {
            let ex = a_vecs[i].x + s * b_vecs[i].x;
            let ey = a_vecs[i].y + s * b_vecs[i].y;
            let ez = a_vecs[i].z + s * b_vecs[i].z;
            assert!((result.x[i] - ex).abs() < TOL, "x[{}] mismatch", i);
            assert!((result.y[i] - ey).abs() < TOL, "y[{}] mismatch", i);
            assert!((result.z[i] - ez).abs() < TOL, "z[{}] mismatch", i);
        }
    }

    #[test]
    fn test_batch_add_scaled_zero_scalar() {
        // s=0 → result must equal a exactly
        let a_vecs = vec![Vector3::new(5.0, -3.0, 1.5), Vector3::new(0.0, 7.0, -2.0)];
        let b_vecs = vec![
            Vector3::new(100.0, 200.0, 300.0),
            Vector3::new(-50.0, -60.0, -70.0),
        ];

        let a = SimdBatch::from_vector3_slice(&a_vecs);
        let b = SimdBatch::from_vector3_slice(&b_vecs);
        let result = batch_add_scaled(&a, &b, 0.0).expect("same length");

        for (i, av) in a_vecs.iter().enumerate() {
            assert!((result.x[i] - av.x).abs() < TOL, "x[{}] mismatch", i);
            assert!((result.y[i] - av.y).abs() < TOL, "y[{}] mismatch", i);
            assert!((result.z[i] - av.z).abs() < TOL, "z[{}] mismatch", i);
        }
    }

    #[test]
    fn test_batch_evolve_rk4_single_spin_matches_scalar() {
        // A 1-spin batch should produce the same result as the scalar LlgSolver.
        use crate::dynamics::llg::LlgSolver;
        let alpha = 0.01_f64;
        let dt = 1.0e-13_f64;

        let m0 = Vector3::new(1.0, 0.0, 0.0);
        let h_field = Vector3::new(0.0, 0.0, 1.0);

        // Scalar reference
        let solver = LlgSolver::new(alpha, dt);
        let h_fn = |_m: Vector3<f64>| h_field;
        let m_scalar = solver.step_rk4(m0, h_fn);

        // Batch
        let m_batch = SimdBatch::from_vector3_slice(&[m0]);
        let h_batch = SimdBatch::from_vector3_slice(&[h_field]);
        let result = batch_evolve_rk4(&m_batch, &h_batch, alpha, GAMMA, dt).expect("same length");

        // Tolerance is loose because LlgSolver normalizes intermediate stages
        // whereas batch_evolve_rk4 only normalizes at the end.
        let tol = 1e-6_f64;
        assert!(
            (result.x[0] - m_scalar.x).abs() < tol,
            "x mismatch: batch={} scalar={}",
            result.x[0],
            m_scalar.x
        );
        assert!(
            (result.y[0] - m_scalar.y).abs() < tol,
            "y mismatch: batch={} scalar={}",
            result.y[0],
            m_scalar.y
        );
        assert!(
            (result.z[0] - m_scalar.z).abs() < tol,
            "z mismatch: batch={} scalar={}",
            result.z[0],
            m_scalar.z
        );
    }

    #[test]
    fn test_batch_evolve_rk4_n_spins_all_correct() {
        // N=4 spins, each with a distinct initial state, should all evolve
        // and produce unit-length magnetizations.
        let alpha = 0.05_f64;
        let dt = 1.0e-13_f64;

        let m_vecs = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(
                1.0 / 3.0_f64.sqrt(),
                1.0 / 3.0_f64.sqrt(),
                1.0 / 3.0_f64.sqrt(),
            ),
        ];
        let h_vecs = vec![
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];

        let m_batch = SimdBatch::from_vector3_slice(&m_vecs);
        let h_batch = SimdBatch::from_vector3_slice(&h_vecs);
        let result = batch_evolve_rk4(&m_batch, &h_batch, alpha, GAMMA, dt).expect("same length");

        assert_eq!(result.len(), m_vecs.len());
        for i in 0..result.len() {
            let mag = (result.x[i].powi(2) + result.y[i].powi(2) + result.z[i].powi(2)).sqrt();
            assert!(
                (mag - 1.0).abs() < 1e-10,
                "spin {} magnitude={} (expected 1.0)",
                i,
                mag
            );
        }
    }

    #[test]
    fn test_batch_evolve_rk4_normalizes() {
        // |m_i| should be ≈ 1 after one RK4 step, regardless of initial magnitude.
        let alpha = 0.1_f64;
        let dt = 1.0e-12_f64;

        let m_batch = SimdBatch::from_vector3_slice(&[
            Vector3::new(0.6, 0.8, 0.0),    // already unit length
            Vector3::new(0.3, 0.0, 0.9539), // approx unit length
        ]);
        let h_batch = SimdBatch::from_vector3_slice(&[
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.0, 0.0, 1.0),
        ]);

        let result = batch_evolve_rk4(&m_batch, &h_batch, alpha, GAMMA, dt).expect("same length");

        for i in 0..result.len() {
            let mag = (result.x[i].powi(2) + result.y[i].powi(2) + result.z[i].powi(2)).sqrt();
            assert!(
                (mag - 1.0).abs() < 1e-10,
                "spin {} magnitude={} after normalization",
                i,
                mag
            );
        }
    }

    #[test]
    fn test_batch_evolve_rk4_zero_field_no_precession() {
        // H = 0 → dm/dt = 0, so m should remain unchanged.
        let alpha = 0.01_f64;
        let dt = 1.0e-12_f64;

        let m_vecs = vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)];
        let h_vecs = vec![Vector3::zero(), Vector3::zero()];

        let m_batch = SimdBatch::from_vector3_slice(&m_vecs);
        let h_batch = SimdBatch::from_vector3_slice(&h_vecs);
        let result = batch_evolve_rk4(&m_batch, &h_batch, alpha, GAMMA, dt).expect("same length");

        for (i, mv) in m_vecs.iter().enumerate() {
            assert!(
                (result.x[i] - mv.x).abs() < TOL,
                "x[{}] changed under zero field",
                i
            );
            assert!(
                (result.y[i] - mv.y).abs() < TOL,
                "y[{}] changed under zero field",
                i
            );
            assert!(
                (result.z[i] - mv.z).abs() < TOL,
                "z[{}] changed under zero field",
                i
            );
        }
    }

    #[test]
    fn test_batch_evolve_rk4_zero_damping_energy_conserved() {
        // alpha = 0 → no dissipation, Zeeman energy E = -m·H should be conserved.
        let alpha = 0.0_f64;
        let dt = 1.0e-13_f64;

        let m0 = Vector3::new(1.0, 0.0, 0.0); // perpendicular to H
        let h0 = Vector3::new(0.0, 0.0, 1.0);

        let m_batch = SimdBatch::from_vector3_slice(&[m0]);
        let h_batch = SimdBatch::from_vector3_slice(&[h0]);
        let result = batch_evolve_rk4(&m_batch, &h_batch, alpha, GAMMA, dt).expect("same length");

        let e0 = -(m0.x * h0.x + m0.y * h0.y + m0.z * h0.z);
        let e1 = -(result.x[0] * h0.x + result.y[0] * h0.y + result.z[0] * h0.z);

        assert!(
            (e1 - e0).abs() < 1e-8,
            "energy drift without damping: Δe={}",
            (e1 - e0).abs()
        );
    }

    #[test]
    fn test_batch_dimension_mismatch_error() {
        // m and h_eff with different sizes must return Err.
        let m = SimdBatch::new(3);
        let h = SimdBatch::new(5);

        let r1 = batch_add_scaled(&m, &h, 1.0);
        assert!(r1.is_err(), "batch_add_scaled should fail on size mismatch");

        let r2 = batch_evolve_rk4(&m, &h, 0.01, GAMMA, 1.0e-13);
        assert!(r2.is_err(), "batch_evolve_rk4 should fail on size mismatch");

        let r3 = batch_evolve_multi_step(m, &h, 0.01, GAMMA, 1.0e-13, 5);
        assert!(
            r3.is_err(),
            "batch_evolve_multi_step should fail on size mismatch"
        );
    }

    #[test]
    fn test_multi_step_conserves_energy_zero_damping() {
        // alpha = 0, N=50 steps: energy should remain conserved to within RK4 truncation error.
        let alpha = 0.0_f64;
        let dt = 1.0e-13_f64;
        let n_steps = 50_usize;

        let m0 = Vector3::new(1.0, 0.0, 0.0);
        let h0 = Vector3::new(0.0, 0.0, 1.0);

        let m_batch = SimdBatch::from_vector3_slice(&[m0]);
        let h_batch = SimdBatch::from_vector3_slice(&[h0]);

        let m_final = batch_evolve_multi_step(m_batch, &h_batch, alpha, GAMMA, dt, n_steps)
            .expect("same length");

        let e0 = -(m0.x * h0.x + m0.y * h0.y + m0.z * h0.z);
        let e_final = -(m_final.x[0] * h0.x + m_final.y[0] * h0.y + m_final.z[0] * h0.z);

        assert!(
            (e_final - e0).abs() < 1e-7,
            "energy drift over {} steps: Δe={}",
            n_steps,
            (e_final - e0).abs()
        );
    }

    #[test]
    fn test_multi_step_reaches_equilibrium() {
        // alpha > 0, many steps: m should align with H (equilibrium).
        let alpha = 0.5_f64; // strong damping
        let dt = 1.0e-12_f64;
        let n_steps = 2000_usize;

        // Start perpendicular to H
        let m0 = Vector3::new(1.0, 0.0, 0.0);
        let h0 = Vector3::new(0.0, 0.0, 1.0); // field along z

        let m_batch = SimdBatch::from_vector3_slice(&[m0]);
        let h_batch = SimdBatch::from_vector3_slice(&[h0]);

        let m_final = batch_evolve_multi_step(m_batch, &h_batch, alpha, GAMMA, dt, n_steps)
            .expect("same length");

        // At equilibrium, m should be near (0, 0, 1)
        let m_z = m_final.z[0];
        assert!(
            m_z > 0.9,
            "spin did not align with field after {} steps: m_z={}",
            n_steps,
            m_z
        );
    }

    #[test]
    fn test_large_batch_n1024_no_panic() {
        // N=1024 spins, 10 steps should complete without panic or error.
        let n = 1024_usize;
        let alpha = 0.01_f64;
        let dt = 1.0e-13_f64;

        let mut m_batch = SimdBatch::new(n);
        let mut h_batch = SimdBatch::new(n);
        for i in 0..n {
            let angle = (i as f64) * std::f64::consts::TAU / n as f64;
            m_batch.x[i] = angle.cos();
            m_batch.y[i] = angle.sin();
            m_batch.z[i] = 0.0;
            h_batch.z[i] = 1.0; // all fields along z
        }

        let result = batch_evolve_multi_step(m_batch, &h_batch, alpha, GAMMA, dt, 10);
        assert!(
            result.is_ok(),
            "N=1024 multi-step failed: {:?}",
            result.err()
        );
        let m_final = result.expect("N=1024 result");
        assert_eq!(m_final.len(), n);
    }

    #[test]
    fn test_batch_vs_dp45_tolerance() {
        // Compare batch_evolve_rk4 single step with a scalar RK4 reference step
        // (using calc_dm_dt directly) to verify the batch result is within 1e-8.
        // We use calc_dm_dt as the reference since DP45 operates on a different API.
        use crate::dynamics::llg::calc_dm_dt;

        let alpha = 0.01_f64;
        let dt = 1.0e-13_f64;

        let m0 = Vector3::new(1.0 / 2.0_f64.sqrt(), 1.0 / 2.0_f64.sqrt(), 0.0);
        let h0 = Vector3::new(0.0, 0.0, 1.0);

        // Scalar RK4 reference (no intermediate normalization, same as batch)
        let k1 = calc_dm_dt(m0, h0, GAMMA, alpha);
        let m2 = m0 + k1 * (dt * 0.5);
        let k2 = calc_dm_dt(m2, h0, GAMMA, alpha);
        let m3 = m0 + k2 * (dt * 0.5);
        let k3 = calc_dm_dt(m3, h0, GAMMA, alpha);
        let m4 = m0 + k3 * dt;
        let k4 = calc_dm_dt(m4, h0, GAMMA, alpha);
        let dm = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
        let m_ref_raw = m0 + dm;
        let m_ref_mag = m_ref_raw.magnitude();
        let m_ref = if m_ref_mag > 0.0 {
            m_ref_raw * (1.0 / m_ref_mag)
        } else {
            m_ref_raw
        };

        // Batch result
        let m_batch = SimdBatch::from_vector3_slice(&[m0]);
        let h_batch = SimdBatch::from_vector3_slice(&[h0]);
        let result = batch_evolve_rk4(&m_batch, &h_batch, alpha, GAMMA, dt).expect("same length");

        let tol = 1e-10_f64;
        assert!(
            (result.x[0] - m_ref.x).abs() < tol,
            "x mismatch vs scalar RK4: batch={} ref={}",
            result.x[0],
            m_ref.x
        );
        assert!(
            (result.y[0] - m_ref.y).abs() < tol,
            "y mismatch vs scalar RK4: batch={} ref={}",
            result.y[0],
            m_ref.y
        );
        assert!(
            (result.z[0] - m_ref.z).abs() < tol,
            "z mismatch vs scalar RK4: batch={} ref={}",
            result.z[0],
            m_ref.z
        );
    }

    #[test]
    fn test_batch_normalize() {
        let mut batch = SimdBatch::new(3);
        batch.x[0] = 3.0;
        batch.y[0] = 4.0;
        batch.z[0] = 0.0;

        batch.x[1] = 0.0;
        batch.y[1] = 0.0;
        batch.z[1] = 5.0;

        // Index 2 stays zero -- should not cause division by zero
        batch_normalize(&mut batch);

        let mag0 =
            (batch.x[0] * batch.x[0] + batch.y[0] * batch.y[0] + batch.z[0] * batch.z[0]).sqrt();
        assert!(
            (mag0 - 1.0).abs() < TOL,
            "magnitude at 0 should be 1, got {}",
            mag0
        );

        let mag1 =
            (batch.x[1] * batch.x[1] + batch.y[1] * batch.y[1] + batch.z[1] * batch.z[1]).sqrt();
        assert!(
            (mag1 - 1.0).abs() < TOL,
            "magnitude at 1 should be 1, got {}",
            mag1
        );

        // Zero vector stays zero
        assert!((batch.x[2]).abs() < TOL);
        assert!((batch.y[2]).abs() < TOL);
        assert!((batch.z[2]).abs() < TOL);
    }
}
