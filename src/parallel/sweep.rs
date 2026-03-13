//! Parallel parameter sweeps for spintronics simulations
//!
//! This module provides utilities for running embarrassingly parallel
//! parameter sweeps using rayon, where independent simulations are
//! executed concurrently across different parameter values.
//!
//! ## Sweep types
//!
//! - **Generic sweep** ([`parallel_sweep`]): map any function over a parameter slice
//! - **Field sweep** ([`field_sweep`]): vary external field, compute equilibrium magnetization
//! - **Temperature sweep** ([`temperature_sweep`]): vary temperature, compute thermal average
//! - **Custom sweep** ([`ParameterSweep`]): builder pattern with progress tracking
//!
//! ## Example
//!
//! ```rust,ignore
//! use spintronics::parallel::sweep::parallel_sweep;
//!
//! let fields: Vec<f64> = (0..100).map(|i| i as f64 * 0.01).collect();
//! let results = parallel_sweep(&fields, |&h| h * h);  // trivial example
//! assert_eq!(results.len(), 100);
//! ```

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use rayon::prelude::*;

use crate::constants::GAMMA;
use crate::error::{Error, Result};
use crate::vector3::Vector3;

// ---------------------------------------------------------------------------
// Generic parallel sweep
// ---------------------------------------------------------------------------

/// Execute a function over each parameter in parallel, collecting results.
///
/// This is the simplest entry point for embarrassingly parallel workloads.
/// The closure `f` is invoked once per element in `params`, and results are
/// returned in the same order as the input parameters.
///
/// # Type parameters
/// * `P` - Parameter type (must be `Send + Sync`)
/// * `R` - Result type (must be `Send`)
/// * `F` - Closure mapping `&P -> R`
///
/// # Example
/// ```rust,ignore
/// let params = vec![1.0_f64, 2.0, 3.0];
/// let squares = parallel_sweep(&params, |&x| x * x);
/// assert_eq!(squares, vec![1.0, 4.0, 9.0]);
/// ```
pub fn parallel_sweep<P, R, F>(params: &[P], f: F) -> Vec<R>
where
    P: Send + Sync,
    R: Send,
    F: Fn(&P) -> R + Send + Sync,
{
    params.par_iter().map(f).collect()
}

/// Execute a function over each parameter in parallel with progress tracking.
///
/// Returns `(results, completed_count)` where `completed_count` is an atomic
/// counter that can be read from another thread to monitor progress.
///
/// # Arguments
/// * `params` - Slice of parameters
/// * `f` - Function to execute per parameter
///
/// # Returns
/// Tuple of (results vector, Arc to progress counter)
pub fn parallel_sweep_with_progress<P, R, F>(params: &[P], f: F) -> (Vec<R>, Arc<AtomicUsize>)
where
    P: Send + Sync,
    R: Send,
    F: Fn(&P) -> R + Send + Sync,
{
    let progress = Arc::new(AtomicUsize::new(0));
    let progress_ref = Arc::clone(&progress);

    let results: Vec<R> = params
        .par_iter()
        .map(|p| {
            let result = f(p);
            progress_ref.fetch_add(1, Ordering::Relaxed);
            result
        })
        .collect();

    (results, progress)
}

// ---------------------------------------------------------------------------
// ParameterSweep builder
// ---------------------------------------------------------------------------

/// A reusable parameter sweep container that stores both parameters and results.
///
/// # Type parameters
/// * `P` - Parameter type
/// * `R` - Result type
#[derive(Debug, Clone)]
pub struct ParameterSweep<P, R> {
    /// Input parameters for the sweep
    pub parameters: Vec<P>,
    /// Results from the sweep (None until executed)
    pub results: Vec<Option<R>>,
}

impl<P, R> ParameterSweep<P, R>
where
    P: Send + Sync + Clone,
    R: Send + Clone,
{
    /// Create a new parameter sweep with the given parameters.
    pub fn new(parameters: Vec<P>) -> Self {
        let n = parameters.len();
        Self {
            parameters,
            results: vec![None; n],
        }
    }

    /// Execute the sweep in parallel, storing results.
    ///
    /// # Arguments
    /// * `f` - Function mapping each parameter to a result
    pub fn execute<F>(&mut self, f: F)
    where
        F: Fn(&P) -> R + Send + Sync,
    {
        let computed: Vec<R> = self.parameters.par_iter().map(f).collect();
        self.results = computed.into_iter().map(Some).collect();
    }

    /// Execute the sweep with an atomic progress counter.
    ///
    /// # Arguments
    /// * `f` - Function mapping each parameter to a result
    ///
    /// # Returns
    /// Arc to the progress counter (reads give number of completed items)
    pub fn execute_with_progress<F>(&mut self, f: F) -> Arc<AtomicUsize>
    where
        F: Fn(&P) -> R + Send + Sync,
    {
        let progress = Arc::new(AtomicUsize::new(0));
        let progress_ref = Arc::clone(&progress);

        let computed: Vec<R> = self
            .parameters
            .par_iter()
            .map(|p| {
                let result = f(p);
                progress_ref.fetch_add(1, Ordering::Relaxed);
                result
            })
            .collect();

        self.results = computed.into_iter().map(Some).collect();
        progress
    }

    /// Return the number of completed results.
    pub fn completed_count(&self) -> usize {
        self.results.iter().filter(|r| r.is_some()).count()
    }

    /// Return true if all results have been computed.
    pub fn is_complete(&self) -> bool {
        self.completed_count() == self.parameters.len()
    }

    /// Collect completed results, returning `None` entries for incomplete ones.
    pub fn results_ref(&self) -> &[Option<R>] {
        &self.results
    }
}

// ---------------------------------------------------------------------------
// Field sweep
// ---------------------------------------------------------------------------

/// Result of a single field-sweep simulation point.
#[derive(Debug, Clone)]
pub struct FieldSweepResult {
    /// Applied external field magnitude \[T\]
    pub field: f64,
    /// Equilibrium magnetization vector (normalized)
    pub magnetization: Vector3<f64>,
    /// Scalar projection of magnetization along field direction
    pub m_parallel: f64,
}

/// Sweep the external magnetic field and compute equilibrium magnetization.
///
/// For each field value, the magnetization is relaxed from `m_init` under
/// the LLG equation until the torque falls below `tol` or `max_steps` is
/// reached.
///
/// # Arguments
/// * `fields` - Slice of field magnitudes \[T\]
/// * `field_direction` - Unit vector for field direction
/// * `m_init` - Initial magnetization (normalized)
/// * `alpha` - Gilbert damping
/// * `dt` - LLG time step \[s\]
/// * `max_steps` - Maximum integration steps per field value
/// * `tol` - Convergence tolerance on |dm/dt|
///
/// # Returns
/// Vector of [`FieldSweepResult`] in the same order as `fields`.
pub fn field_sweep(
    fields: &[f64],
    field_direction: Vector3<f64>,
    m_init: Vector3<f64>,
    alpha: f64,
    dt: f64,
    max_steps: usize,
    tol: f64,
) -> Vec<FieldSweepResult> {
    let dir = field_direction.normalize();

    fields
        .par_iter()
        .map(|&h_mag| {
            let h_ext = dir * h_mag;
            let mut m = m_init.normalize();

            for _ in 0..max_steps {
                let dm_dt = llg_torque_simple(m, h_ext, alpha);
                let torque_mag = dm_dt.magnitude();
                if torque_mag < tol {
                    break;
                }
                m = (m + dm_dt * dt).normalize();
            }

            let m_par = m.dot(&dir);

            FieldSweepResult {
                field: h_mag,
                magnetization: m,
                m_parallel: m_par,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Temperature sweep
// ---------------------------------------------------------------------------

/// Result of a single temperature-sweep simulation point.
#[derive(Debug, Clone)]
pub struct TemperatureSweepResult {
    /// Temperature \[K\]
    pub temperature: f64,
    /// Average magnetization magnitude (0..1)
    pub avg_magnetization: f64,
    /// Standard deviation of magnetization fluctuations
    pub magnetization_std: f64,
}

/// Sweep temperature and compute thermally-averaged magnetization.
///
/// At each temperature, a simplified mean-field Langevin model is used:
///   M(T) = L(mu * H_eff / (k_B * T))
/// where L(x) = coth(x) - 1/x is the Langevin function.
///
/// For a more accurate simulation with stochastic LLG, use the
/// `stochastic` module combined with [`parallel_sweep`].
///
/// # Arguments
/// * `temperatures` - Slice of temperatures \[K\]
/// * `h_eff` - Effective field magnitude \[T\]
/// * `mu` - Magnetic moment per spin \[J/T\]
/// * `curie_temp` - Curie temperature \[K\] (for mean-field correction)
///
/// # Returns
/// Vector of [`TemperatureSweepResult`] in the same order as `temperatures`.
///
/// # Errors
/// Returns error if `curie_temp` is not positive.
pub fn temperature_sweep(
    temperatures: &[f64],
    h_eff: f64,
    mu: f64,
    curie_temp: f64,
) -> Result<Vec<TemperatureSweepResult>> {
    if curie_temp <= 0.0 {
        return Err(Error::InvalidParameter {
            param: "curie_temp".to_string(),
            reason: "Curie temperature must be positive".to_string(),
        });
    }

    let kb = 1.380649e-23; // Boltzmann constant [J/K]

    let results: Vec<TemperatureSweepResult> = temperatures
        .par_iter()
        .map(|&temp| {
            if temp <= 0.0 {
                // At T=0, magnetization is saturated
                return TemperatureSweepResult {
                    temperature: temp,
                    avg_magnetization: 1.0,
                    magnetization_std: 0.0,
                };
            }

            // Mean-field self-consistent Langevin equation:
            // m = L( (mu * h_eff + 3 * T_c * m * kb) / (kb * T) )
            // Solve iteratively
            let mut m = 0.5_f64; // initial guess
            for _ in 0..200 {
                let h_total = mu * h_eff + 3.0 * curie_temp * m * kb;
                let x = h_total / (kb * temp);
                let m_new = langevin(x);
                let dm = (m_new - m).abs();
                m = m_new;
                if dm < 1e-12 {
                    break;
                }
            }

            // Fluctuation estimate from susceptibility
            // chi ~ d<m>/dT, sigma ~ sqrt(kB*T*chi / N) ≈ simplified estimate
            let sigma = if temp < curie_temp {
                ((kb * temp) / (mu * h_eff.max(1e-30))).sqrt() * (1.0 - m).max(0.0)
            } else {
                ((kb * temp) / (mu * h_eff.max(1e-30))).sqrt().min(1.0)
            };

            TemperatureSweepResult {
                temperature: temp,
                avg_magnetization: m.clamp(0.0, 1.0),
                magnetization_std: sigma.clamp(0.0, 1.0),
            }
        })
        .collect();

    Ok(results)
}

// ---------------------------------------------------------------------------
// Multi-parameter sweep
// ---------------------------------------------------------------------------

/// A point in a 2D parameter sweep grid.
#[derive(Debug, Clone)]
pub struct SweepPoint2D<R> {
    /// First parameter value
    pub param1: f64,
    /// Second parameter value
    pub param2: f64,
    /// Result at this point
    pub result: R,
}

/// Execute a 2D parameter sweep (grid of param1 x param2).
///
/// All grid points are evaluated in parallel.
///
/// # Arguments
/// * `params1` - First parameter axis
/// * `params2` - Second parameter axis
/// * `f` - Function `(p1, p2) -> R`
///
/// # Returns
/// Vec of [`SweepPoint2D`] covering the full grid in row-major order.
pub fn parallel_sweep_2d<R, F>(params1: &[f64], params2: &[f64], f: F) -> Vec<SweepPoint2D<R>>
where
    R: Send,
    F: Fn(f64, f64) -> R + Send + Sync,
{
    // Build flat grid of (p1, p2) pairs
    let grid: Vec<(f64, f64)> = params1
        .iter()
        .flat_map(|&p1| params2.iter().map(move |&p2| (p1, p2)))
        .collect();

    grid.par_iter()
        .map(|&(p1, p2)| SweepPoint2D {
            param1: p1,
            param2: p2,
            result: f(p1, p2),
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Simplified LLG torque for single-spin relaxation (no exchange).
fn llg_torque_simple(m: Vector3<f64>, h_ext: Vector3<f64>, alpha: f64) -> Vector3<f64> {
    let m_cross_h = m.cross(&h_ext);
    let m_cross_m_cross_h = m.cross(&m_cross_h);
    let prefactor = -GAMMA / (1.0 + alpha * alpha);
    (m_cross_h + m_cross_m_cross_h * alpha) * prefactor
}

/// Langevin function L(x) = coth(x) - 1/x.
///
/// Handles small x via Taylor expansion to avoid numerical blow-up.
fn langevin(x: f64) -> f64 {
    if x.abs() < 1e-4 {
        // Taylor expansion: L(x) ≈ x/3 - x^3/45 + ...
        x / 3.0 - x * x * x / 45.0
    } else if x.abs() > 20.0 {
        // For large |x|, coth(x) = 1 + 2*exp(-2x)/(1-exp(-2x)) ≈ 1 + 2*exp(-2x)
        // L(x) = coth(x) - 1/x ≈ sign(x) + 2*sign(x)*exp(-2|x|) - 1/x
        let abs_x = x.abs();
        let correction = 2.0 * (-2.0 * abs_x).exp();
        x.signum() * (1.0 + correction) - 1.0 / x
    } else {
        let coth = x.cosh() / x.sinh();
        coth - 1.0 / x
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parallel_sweep_trivial() {
        let params: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let results = parallel_sweep(&params, |&x| x * x);

        assert_eq!(results.len(), 50);
        for (i, &r) in results.iter().enumerate() {
            let expected = (i as f64) * (i as f64);
            assert!(
                (r - expected).abs() < 1e-15,
                "index {}: got {}, expected {}",
                i,
                r,
                expected,
            );
        }
    }

    #[test]
    fn test_parallel_sweep_deterministic() {
        let params: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();

        let r1 = parallel_sweep(&params, |&x| x.sin());
        let r2 = parallel_sweep(&params, |&x| x.sin());

        assert_eq!(r1.len(), r2.len());
        for (a, b) in r1.iter().zip(r2.iter()) {
            assert!(
                (a - b).abs() < 1e-15,
                "determinism violated: {} != {}",
                a,
                b,
            );
        }
    }

    #[test]
    fn test_parameter_sweep_builder() {
        let params: Vec<i32> = (0..20).collect();
        let mut sweep = ParameterSweep::<i32, i64>::new(params);

        assert!(!sweep.is_complete());
        assert_eq!(sweep.completed_count(), 0);

        sweep.execute(|&x| (x as i64) * (x as i64));

        assert!(sweep.is_complete());
        assert_eq!(sweep.completed_count(), 20);

        let results = sweep.results_ref();
        assert_eq!(results[5], Some(25), "5^2 should be 25");
    }

    #[test]
    fn test_field_sweep_saturation() {
        // At very large fields the magnetization should align with the field
        // Use dt small enough for stability: omega*dt < 1 => dt < 1/(GAMMA*H)
        // For H=10T: dt < 1/(1.76e11*10) ~ 5.7e-13, so use 1e-13
        let fields: Vec<f64> = vec![0.01, 0.1, 1.0, 10.0];
        let dir = Vector3::new(0.0, 0.0, 1.0);
        let m_init = Vector3::new(1.0, 0.0, 0.0); // perpendicular to field
        let alpha = 0.5; // heavy damping for fast convergence
        let dt = 1e-13;
        let max_steps = 1_000_000;
        let tol = 1e-10;

        let results = field_sweep(&fields, dir, m_init, alpha, dt, max_steps, tol);

        assert_eq!(results.len(), 4);

        // At largest field, magnetization should be nearly aligned with z
        let last = &results[3];
        assert!(
            last.m_parallel > 0.9,
            "at H=10T, m_parallel={} should be > 0.9",
            last.m_parallel,
        );
    }

    #[test]
    fn test_temperature_sweep_ordering() {
        let temps: Vec<f64> = vec![10.0, 100.0, 300.0, 500.0, 800.0, 1000.0];
        let h_eff = 1.0; // 1 T
        let mu = 9.274e-24; // Bohr magneton
        let tc = 600.0; // Curie temperature

        let results =
            temperature_sweep(&temps, h_eff, mu, tc).expect("temperature sweep should succeed");

        assert_eq!(results.len(), 6);

        // Magnetization should generally decrease with temperature
        // (not strictly monotonic due to mean-field approx, but at high T should be small)
        let m_low = results[0].avg_magnetization;
        let m_high = results[5].avg_magnetization;
        assert!(
            m_low >= m_high,
            "m(10K)={} should be >= m(1000K)={}",
            m_low,
            m_high,
        );
    }

    #[test]
    fn test_temperature_sweep_error_on_bad_curie() {
        let temps = vec![100.0];
        let result = temperature_sweep(&temps, 1.0, 1e-23, 0.0);
        assert!(result.is_err());

        let result = temperature_sweep(&temps, 1.0, 1e-23, -100.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_sweep_with_progress() {
        let params: Vec<u32> = (0..30).collect();
        let (results, progress) = parallel_sweep_with_progress(&params, |&x| x * 2);

        assert_eq!(results.len(), 30);
        assert_eq!(progress.load(Ordering::Relaxed), 30);

        for (i, &r) in results.iter().enumerate() {
            assert_eq!(r, (i as u32) * 2);
        }
    }

    #[test]
    fn test_2d_sweep() {
        let p1: Vec<f64> = vec![1.0, 2.0, 3.0];
        let p2: Vec<f64> = vec![10.0, 20.0];

        let results = parallel_sweep_2d(&p1, &p2, |a, b| a + b);

        assert_eq!(results.len(), 6); // 3 x 2

        // Verify all expected (p1, p2) pairs are present
        let mut found = [false; 6];
        for r in &results {
            let expected = r.param1 + r.param2;
            assert!(
                (r.result - expected).abs() < 1e-15,
                "({}, {}): got {}, expected {}",
                r.param1,
                r.param2,
                r.result,
                expected,
            );
            // Mark as found
            let idx = p1
                .iter()
                .position(|&x| (x - r.param1).abs() < 1e-15)
                .expect("param1 should be in p1");
            let idy = p2
                .iter()
                .position(|&x| (x - r.param2).abs() < 1e-15)
                .expect("param2 should be in p2");
            found[idx * p2.len() + idy] = true;
        }
        assert!(found.iter().all(|&f| f), "not all grid points covered");
    }

    #[test]
    fn test_langevin_function() {
        // L(0) = 0
        assert!((langevin(0.0)).abs() < 1e-10);

        // L(x) -> 1 for large x; L(x) ≈ 1 - 1/x, so L(100) ≈ 0.99
        assert!((langevin(100.0) - 1.0).abs() < 0.011);

        // L(x) ≈ x/3 for small x
        let x = 1e-5;
        assert!((langevin(x) - x / 3.0).abs() < 1e-15);
    }
}
