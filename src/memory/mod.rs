//! Memory Pool Allocator for High-Performance Simulations
//!
//! This module provides memory pooling utilities to reduce allocation overhead
//! in hot paths during spin dynamics simulations. Frequent allocations can
//! significantly impact performance in time-stepping loops.
//!
//! ## Features
//!
//! - **Object pooling**: Reuse `Vec<f64>` and `Vec<Vector3>` allocations
//! - **Thread-local pools**: No synchronization overhead in single-threaded code
//! - **Automatic growth**: Pools grow as needed, shrink when idle
//! - **Zero-copy operations**: Borrow from pool without cloning
//!
//! ## Example Usage
//!
//! ```rust
//! use spintronics::memory::{VectorPool, SpinArrayPool};
//! use spintronics::Vector3;
//!
//! // Create pools
//! let mut f64_pool = VectorPool::<f64>::new();
//! let mut spin_pool = SpinArrayPool::new();
//!
//! // Get a vector from the pool (or allocate new one)
//! let mut work_buffer = f64_pool.get(1000);
//!
//! // Use the buffer
//! for i in 0..1000 {
//!     work_buffer[i] = i as f64;
//! }
//!
//! // Return to pool for reuse
//! f64_pool.put(work_buffer);
//!
//! // Get a spin array
//! let mut spins = spin_pool.get(1000);
//! for spin in spins.iter_mut() {
//!     *spin = Vector3::new(1.0, 0.0, 0.0);
//! }
//! spin_pool.put(spins);
//! ```

use std::cell::RefCell;

use crate::vector3::Vector3;

/// Memory pool for `Vec<T>` allocations
///
/// Maintains a cache of previously allocated vectors to reduce
/// allocation overhead in hot loops.
pub struct VectorPool<T> {
    pool: Vec<Vec<T>>,
    max_cached: usize,
}

impl<T: Default + Clone> VectorPool<T> {
    /// Create a new vector pool
    pub fn new() -> Self {
        Self {
            pool: Vec::new(),
            max_cached: 16,
        }
    }

    /// Create a pool with custom maximum cache size
    pub fn with_capacity(max_cached: usize) -> Self {
        Self {
            pool: Vec::with_capacity(max_cached),
            max_cached,
        }
    }

    /// Get a vector from the pool or allocate a new one
    ///
    /// # Arguments
    /// * `size` - Required minimum size of the vector
    ///
    /// # Returns
    /// A vector with at least `size` elements, all set to default value
    pub fn get(&mut self, size: usize) -> Vec<T> {
        // Try to find a suitable vector in the pool
        if let Some(idx) = self.pool.iter().position(|v| v.capacity() >= size) {
            let mut vec = self.pool.swap_remove(idx);
            vec.clear();
            vec.resize(size, T::default());
            vec
        } else {
            // Allocate a new vector
            vec![T::default(); size]
        }
    }

    /// Return a vector to the pool for reuse
    ///
    /// The vector's contents will be discarded but its capacity preserved.
    pub fn put(&mut self, vec: Vec<T>) {
        if self.pool.len() < self.max_cached && vec.capacity() > 0 {
            self.pool.push(vec);
        }
        // Otherwise, let it drop
    }

    /// Clear the pool, releasing all cached memory
    pub fn clear(&mut self) {
        self.pool.clear();
    }

    /// Number of vectors currently in the pool
    pub fn cached_count(&self) -> usize {
        self.pool.len()
    }

    /// Total capacity of all cached vectors
    pub fn total_capacity(&self) -> usize {
        self.pool.iter().map(|v| v.capacity()).sum()
    }
}

impl<T: Default + Clone> Default for VectorPool<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Specialized pool for spin arrays (`Vec<Vector3<f64>>`)
///
/// Optimized for the common case of magnetization state vectors.
pub struct SpinArrayPool {
    pool: Vec<Vec<Vector3<f64>>>,
    max_cached: usize,
}

impl SpinArrayPool {
    /// Create a new spin array pool
    pub fn new() -> Self {
        Self {
            pool: Vec::new(),
            max_cached: 8,
        }
    }

    /// Get a spin array from the pool or allocate a new one
    pub fn get(&mut self, size: usize) -> Vec<Vector3<f64>> {
        if let Some(idx) = self.pool.iter().position(|v| v.capacity() >= size) {
            let mut vec = self.pool.swap_remove(idx);
            vec.clear();
            vec.resize(size, Vector3::new(0.0, 0.0, 0.0));
            vec
        } else {
            vec![Vector3::new(0.0, 0.0, 0.0); size]
        }
    }

    /// Return a spin array to the pool
    pub fn put(&mut self, vec: Vec<Vector3<f64>>) {
        if self.pool.len() < self.max_cached && vec.capacity() > 0 {
            self.pool.push(vec);
        }
    }

    /// Clear the pool
    pub fn clear(&mut self) {
        self.pool.clear();
    }
}

impl Default for SpinArrayPool {
    fn default() -> Self {
        Self::new()
    }
}

// Thread-local vector pool for convenience
// Uses thread-local storage to provide a global pool without explicit management.
thread_local! {
    static F64_POOL: RefCell<VectorPool<f64>> = RefCell::new(VectorPool::new());
    static SPIN_POOL: RefCell<SpinArrayPool> = RefCell::new(SpinArrayPool::new());
}

/// Get a f64 vector from the thread-local pool
pub fn get_f64_vec(size: usize) -> Vec<f64> {
    F64_POOL.with(|pool| pool.borrow_mut().get(size))
}

/// Return a f64 vector to the thread-local pool
pub fn put_f64_vec(vec: Vec<f64>) {
    F64_POOL.with(|pool| pool.borrow_mut().put(vec));
}

/// Get a spin array from the thread-local pool
pub fn get_spin_array(size: usize) -> Vec<Vector3<f64>> {
    SPIN_POOL.with(|pool| pool.borrow_mut().get(size))
}

/// Return a spin array to the thread-local pool
pub fn put_spin_array(vec: Vec<Vector3<f64>>) {
    SPIN_POOL.with(|pool| pool.borrow_mut().put(vec));
}

/// Workspace buffer for RK4 integration
///
/// Preallocates all temporary buffers needed for RK4 time stepping.
pub struct Rk4Workspace {
    /// k1 coefficient
    pub k1: Vec<Vector3<f64>>,
    /// k2 coefficient
    pub k2: Vec<Vector3<f64>>,
    /// k3 coefficient
    pub k3: Vec<Vector3<f64>>,
    /// k4 coefficient
    pub k4: Vec<Vector3<f64>>,
    /// Temporary state
    pub temp: Vec<Vector3<f64>>,
}

impl Rk4Workspace {
    /// Create a new RK4 workspace for a given system size
    pub fn new(size: usize) -> Self {
        let zero = Vector3::new(0.0, 0.0, 0.0);
        Self {
            k1: vec![zero; size],
            k2: vec![zero; size],
            k3: vec![zero; size],
            k4: vec![zero; size],
            temp: vec![zero; size],
        }
    }

    /// Resize the workspace for a different system size
    pub fn resize(&mut self, size: usize) {
        let zero = Vector3::new(0.0, 0.0, 0.0);
        self.k1.resize(size, zero);
        self.k2.resize(size, zero);
        self.k3.resize(size, zero);
        self.k4.resize(size, zero);
        self.temp.resize(size, zero);
    }

    /// Clear all buffers (set to zero)
    pub fn clear(&mut self) {
        let zero = Vector3::new(0.0, 0.0, 0.0);
        for k in self.k1.iter_mut() {
            *k = zero;
        }
        for k in self.k2.iter_mut() {
            *k = zero;
        }
        for k in self.k3.iter_mut() {
            *k = zero;
        }
        for k in self.k4.iter_mut() {
            *k = zero;
        }
        for t in self.temp.iter_mut() {
            *t = zero;
        }
    }
}

/// Workspace buffer for Heun integration (stochastic LLG)
pub struct HeunWorkspace {
    /// Predictor step
    pub predictor: Vec<Vector3<f64>>,
    /// Corrector derivative
    pub corrector: Vec<Vector3<f64>>,
    /// Noise term
    pub noise: Vec<Vector3<f64>>,
}

impl HeunWorkspace {
    /// Create a new Heun workspace
    pub fn new(size: usize) -> Self {
        let zero = Vector3::new(0.0, 0.0, 0.0);
        Self {
            predictor: vec![zero; size],
            corrector: vec![zero; size],
            noise: vec![zero; size],
        }
    }

    /// Resize the workspace
    pub fn resize(&mut self, size: usize) {
        let zero = Vector3::new(0.0, 0.0, 0.0);
        self.predictor.resize(size, zero);
        self.corrector.resize(size, zero);
        self.noise.resize(size, zero);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vector_pool_basic() {
        let mut pool: VectorPool<f64> = VectorPool::new();

        // Get a vector
        let v1 = pool.get(100);
        assert_eq!(v1.len(), 100);
        assert!(v1.iter().all(|&x| x == 0.0));

        // Return it
        pool.put(v1);
        assert_eq!(pool.cached_count(), 1);

        // Get it back
        let v2 = pool.get(50);
        assert_eq!(v2.len(), 50);
        assert_eq!(pool.cached_count(), 0);
    }

    #[test]
    fn test_spin_array_pool() {
        let mut pool = SpinArrayPool::new();

        let spins = pool.get(100);
        assert_eq!(spins.len(), 100);

        pool.put(spins);

        let spins2 = pool.get(50);
        assert_eq!(spins2.len(), 50);
    }

    #[test]
    fn test_thread_local_pool() {
        let v = get_f64_vec(100);
        assert_eq!(v.len(), 100);
        put_f64_vec(v);

        let spins = get_spin_array(100);
        assert_eq!(spins.len(), 100);
        put_spin_array(spins);
    }

    #[test]
    fn test_rk4_workspace() {
        let ws = Rk4Workspace::new(100);
        assert_eq!(ws.k1.len(), 100);
        assert_eq!(ws.k2.len(), 100);
        assert_eq!(ws.k3.len(), 100);
        assert_eq!(ws.k4.len(), 100);
        assert_eq!(ws.temp.len(), 100);
    }

    #[test]
    fn test_heun_workspace() {
        let ws = HeunWorkspace::new(100);
        assert_eq!(ws.predictor.len(), 100);
        assert_eq!(ws.corrector.len(), 100);
        assert_eq!(ws.noise.len(), 100);
    }
}
