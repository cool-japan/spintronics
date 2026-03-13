//! Magnon Bose-Einstein Condensation
//!
//! This module implements the physics of magnon Bose-Einstein condensation (BEC),
//! where magnons -- bosonic quasiparticles of spin waves -- accumulate at the
//! lowest energy state under parametric pumping conditions.
//!
//! # Physical Background
//!
//! Magnons obey Bose-Einstein statistics. When the magnon density exceeds a
//! critical value (or equivalently the temperature falls below a critical
//! temperature), the chemical potential reaches the minimum magnon energy and
//! a macroscopic number of magnons occupy the ground state, forming a
//! Bose-Einstein condensate.
//!
//! The condensate is described by a macroscopic wavefunction ψ whose dynamics
//! follow the Gross-Pitaevskii equation (GPE), analogous to atomic BEC but
//! with magnon-specific interactions and dissipation.
//!
//! # Key Features
//!
//! - Bose-Einstein distribution function for magnons
//! - Critical density and BEC transition temperature
//! - Parametric pumping and Suhl instability threshold
//! - Gross-Pitaevskii equation solver (split-step / explicit Euler)
//! - Magnon supercurrent from condensate phase gradient
//!
//! # References
//!
//! - Demokritov et al., Nature 443, 430 (2006) -- first observation of magnon BEC
//! - Serga et al., Nature Communications 5, 3452 (2014) -- Bose-Einstein condensation
//!   of quasi-equilibrium magnons

use std::f64::consts::PI;

use crate::constants::{HBAR, KB};
use crate::error::{self, Result};

// ============================================================================
// Magnon Bose-Einstein Distribution
// ============================================================================

/// Calculate the Bose-Einstein distribution for magnons.
///
/// Returns the average occupation number at energy ℏω for temperature T
/// and chemical potential μ:
///
/// n(ω) = 1 / (exp((ℏω - μ) / (k_B T)) - 1)
///
/// # Arguments
/// * `omega` - Angular frequency of the magnon mode \[rad/s\]
/// * `temperature` - Temperature \[K\]
/// * `chemical_potential` - Chemical potential \[J\]
///
/// # Errors
/// Returns error if temperature is non-positive or if the exponent leads
/// to divergence (μ ≥ ℏω at finite T).
pub fn magnon_distribution(omega: f64, temperature: f64, chemical_potential: f64) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(error::invalid_param(
            "temperature",
            "must be positive for Bose-Einstein distribution",
        ));
    }
    if omega < 0.0 {
        return Err(error::invalid_param("omega", "must be non-negative"));
    }

    let energy = HBAR * omega;
    let exponent = (energy - chemical_potential) / (KB * temperature);

    if exponent <= 0.0 {
        return Err(error::numerical_error(
            "chemical potential exceeds or equals magnon energy; \
             Bose-Einstein distribution diverges",
        ));
    }

    // Guard against overflow in exp
    if exponent > 700.0 {
        // Classical limit: distribution → exp(-(ℏω - μ)/(k_B T))
        return Ok((-exponent).exp());
    }

    let n = 1.0 / (exponent.exp() - 1.0);
    Ok(n)
}

/// Calculate the classical (Boltzmann) limit of the magnon distribution.
///
/// At high temperatures (ℏω - μ) >> k_B T, the Bose-Einstein distribution
/// reduces to the Maxwell-Boltzmann form:
///
/// n(ω) ≈ exp(-(ℏω - μ) / (k_B T))
///
/// # Arguments
/// * `omega` - Angular frequency \[rad/s\]
/// * `temperature` - Temperature \[K\]
/// * `chemical_potential` - Chemical potential \[J\]
pub fn classical_distribution(
    omega: f64,
    temperature: f64,
    chemical_potential: f64,
) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(error::invalid_param("temperature", "must be positive"));
    }

    let energy = HBAR * omega;
    let exponent = (energy - chemical_potential) / (KB * temperature);
    Ok((-exponent).exp())
}

// ============================================================================
// BEC Critical Conditions
// ============================================================================

/// Calculate the critical magnon density for BEC in 3D.
///
/// The critical density is given by:
///
/// n_c = ζ(3/2) · (m* k_B T / (2π ℏ²))^(3/2)
///
/// where ζ(3/2) ≈ 2.612 is the Riemann zeta function.
///
/// # Arguments
/// * `effective_mass` - Effective magnon mass \[kg\]
/// * `temperature` - Temperature \[K\]
///
/// # Errors
/// Returns error if effective mass or temperature is non-positive.
pub fn critical_density(effective_mass: f64, temperature: f64) -> Result<f64> {
    if effective_mass <= 0.0 {
        return Err(error::invalid_param("effective_mass", "must be positive"));
    }
    if temperature <= 0.0 {
        return Err(error::invalid_param(
            "temperature",
            "must be positive for BEC critical density",
        ));
    }

    // Riemann zeta function ζ(3/2)
    let zeta_3_2: f64 = 2.612_375_348_685_488;

    let thermal_factor = effective_mass * KB * temperature / (2.0 * PI * HBAR * HBAR);
    let n_c = zeta_3_2 * thermal_factor.powf(1.5);

    Ok(n_c)
}

/// Calculate the BEC transition temperature for a given magnon density.
///
/// Inverting the critical density formula:
///
/// T_BEC = (2π ℏ² / (m* k_B)) · (n / ζ(3/2))^(2/3)
///
/// # Arguments
/// * `magnon_density` - Magnon number density \[1/m³\]
/// * `effective_mass` - Effective magnon mass \[kg\]
///
/// # Errors
/// Returns error if magnon density or effective mass is non-positive.
pub fn bec_temperature(magnon_density: f64, effective_mass: f64) -> Result<f64> {
    if magnon_density <= 0.0 {
        return Err(error::invalid_param("magnon_density", "must be positive"));
    }
    if effective_mass <= 0.0 {
        return Err(error::invalid_param("effective_mass", "must be positive"));
    }

    let zeta_3_2: f64 = 2.612_375_348_685_488;

    let prefactor = 2.0 * PI * HBAR * HBAR / (effective_mass * KB);
    let density_factor = (magnon_density / zeta_3_2).powf(2.0 / 3.0);

    Ok(prefactor * density_factor)
}

// ============================================================================
// Magnon Condensate
// ============================================================================

/// Magnon Bose-Einstein condensate state.
///
/// Describes the thermodynamic and quantum state of a magnon condensate,
/// including the condensate fraction and interaction parameters needed
/// for Gross-Pitaevskii dynamics.
#[derive(Debug, Clone)]
pub struct MagnonCondensate {
    /// Temperature \[K\]
    pub temperature: f64,

    /// Total magnon number density \[1/m³\]
    pub magnon_density: f64,

    /// Chemical potential \[J\]
    pub chemical_potential: f64,

    /// Effective magnon mass from dispersion curvature \[kg\]
    ///
    /// m* = ℏ² / (d²ε/dk²) evaluated at the band minimum
    pub effective_mass: f64,

    /// Magnon-magnon interaction strength \[J·m³\]
    ///
    /// The s-wave scattering parameter for magnon-magnon interactions
    pub interaction_strength: f64,

    /// Fraction of magnons in the condensate (0 to 1)
    pub condensate_fraction: f64,
}

impl MagnonCondensate {
    /// Create a new magnon condensate state.
    ///
    /// The condensate fraction is computed automatically from the
    /// temperature and density using the 3D ideal Bose gas formula:
    ///
    /// f_c = 1 - (T / T_BEC)^(3/2)   for T < T_BEC
    /// f_c = 0                         for T ≥ T_BEC
    ///
    /// # Arguments
    /// * `temperature` - Temperature \[K\]
    /// * `magnon_density` - Total magnon density \[1/m³\]
    /// * `effective_mass` - Effective magnon mass \[kg\]
    /// * `interaction_strength` - Magnon-magnon interaction \[J·m³\]
    ///
    /// # Errors
    /// Returns error for non-physical parameters.
    pub fn new(
        temperature: f64,
        magnon_density: f64,
        effective_mass: f64,
        interaction_strength: f64,
    ) -> Result<Self> {
        if temperature < 0.0 {
            return Err(error::invalid_param("temperature", "must be non-negative"));
        }
        if magnon_density <= 0.0 {
            return Err(error::invalid_param("magnon_density", "must be positive"));
        }
        if effective_mass <= 0.0 {
            return Err(error::invalid_param("effective_mass", "must be positive"));
        }

        // Calculate BEC temperature
        let t_bec = bec_temperature(magnon_density, effective_mass)?;

        // Condensate fraction
        let condensate_fraction = if temperature < f64::EPSILON {
            // T = 0: all magnons in condensate
            1.0
        } else if temperature < t_bec {
            1.0 - (temperature / t_bec).powf(1.5)
        } else {
            0.0
        };

        // Chemical potential: approaches ℏω_min at BEC
        // Below T_BEC, μ = 0 (in the frame where ω_min = 0)
        // Above T_BEC, μ is negative (below the band minimum)
        let chemical_potential = if temperature < t_bec {
            0.0 // pinned at the band minimum
        } else {
            // Approximate: μ ≈ k_B T ln(1 - (T_BEC/T)^(3/2))
            // This is negative above T_BEC
            let ratio = (t_bec / temperature).powf(1.5);
            if ratio >= 1.0 {
                0.0
            } else {
                KB * temperature * (1.0 - ratio).ln()
            }
        };

        Ok(Self {
            temperature,
            magnon_density,
            chemical_potential,
            effective_mass,
            interaction_strength,
            condensate_fraction,
        })
    }

    /// Get the BEC transition temperature for this system.
    pub fn transition_temperature(&self) -> Result<f64> {
        bec_temperature(self.magnon_density, self.effective_mass)
    }

    /// Get the condensate density n_0 = f_c · n_total.
    pub fn condensate_density(&self) -> f64 {
        self.condensate_fraction * self.magnon_density
    }

    /// Get the thermal (non-condensed) magnon density.
    pub fn thermal_density(&self) -> f64 {
        (1.0 - self.condensate_fraction) * self.magnon_density
    }

    /// Calculate the healing length of the condensate.
    ///
    /// ξ = ℏ / √(2 m* g n_0)
    ///
    /// The healing length sets the scale over which the condensate
    /// wavefunction recovers from a perturbation.
    ///
    /// # Errors
    /// Returns error if condensate density is zero or interaction is non-positive.
    pub fn healing_length(&self) -> Result<f64> {
        let n0 = self.condensate_density();
        if n0 <= 0.0 {
            return Err(error::numerical_error(
                "no condensate present; healing length is undefined",
            ));
        }
        if self.interaction_strength <= 0.0 {
            return Err(error::invalid_param(
                "interaction_strength",
                "must be positive for healing length calculation",
            ));
        }

        let xi = HBAR / (2.0 * self.effective_mass * self.interaction_strength * n0).sqrt();
        Ok(xi)
    }

    /// Calculate the speed of sound in the condensate (Bogoliubov sound).
    ///
    /// c_s = √(g n_0 / m*)
    ///
    /// # Errors
    /// Returns error if condensate density is zero.
    pub fn sound_speed(&self) -> Result<f64> {
        let n0 = self.condensate_density();
        if n0 <= 0.0 {
            return Err(error::numerical_error(
                "no condensate present; sound speed is undefined",
            ));
        }
        if self.interaction_strength <= 0.0 {
            return Err(error::invalid_param(
                "interaction_strength",
                "must be positive for sound speed calculation",
            ));
        }

        let cs = (self.interaction_strength * n0 / self.effective_mass).sqrt();
        Ok(cs)
    }
}

// ============================================================================
// Parametric Pumping
// ============================================================================

/// Parameters for parametric magnon pumping.
///
/// In parametric pumping, an RF field at frequency ω_p creates magnon
/// pairs at ω_p/2 via three-magnon or parametric processes. When the
/// pumping exceeds the Suhl instability threshold, an exponential growth
/// of magnon population occurs.
#[derive(Debug, Clone)]
pub struct ParametricPumping {
    /// Pumping frequency \[rad/s\]
    pub pump_frequency: f64,

    /// Gilbert damping parameter (dimensionless)
    pub damping: f64,

    /// Parametric coupling strength |V_k| \[J\]
    ///
    /// Depends on the magnon wavevector and material parameters
    pub coupling_strength: f64,
}

impl ParametricPumping {
    /// Create new parametric pumping parameters.
    ///
    /// # Arguments
    /// * `pump_frequency` - Pumping frequency \[rad/s\]
    /// * `damping` - Gilbert damping constant
    /// * `coupling_strength` - Parametric coupling |V_k| \[J\]
    ///
    /// # Errors
    /// Returns error for non-physical parameters.
    pub fn new(pump_frequency: f64, damping: f64, coupling_strength: f64) -> Result<Self> {
        if pump_frequency <= 0.0 {
            return Err(error::invalid_param("pump_frequency", "must be positive"));
        }
        if damping <= 0.0 {
            return Err(error::invalid_param("damping", "must be positive"));
        }
        if coupling_strength <= 0.0 {
            return Err(error::invalid_param(
                "coupling_strength",
                "must be positive",
            ));
        }

        Ok(Self {
            pump_frequency,
            damping,
            coupling_strength,
        })
    }

    /// Frequency of the parametrically created magnons.
    ///
    /// Magnon pairs are created at half the pump frequency: ω_k = ω_p / 2
    pub fn magnon_frequency(&self) -> f64 {
        self.pump_frequency / 2.0
    }

    /// Calculate the Suhl instability threshold field.
    ///
    /// h_th = 2 α ω_k / |V_k|
    ///
    /// Above this threshold, parametric magnon amplification occurs.
    /// The threshold is expressed as an RF field amplitude [A/m equivalent].
    pub fn suhl_threshold(&self) -> f64 {
        let omega_k = self.magnon_frequency();
        2.0 * self.damping * omega_k / self.coupling_strength
    }

    /// Calculate the magnon growth rate above threshold.
    ///
    /// For h > h_th, the magnon amplitude grows as exp(Γt) where:
    ///
    /// Γ = |V_k| h / 2 - α ω_k
    ///
    /// # Arguments
    /// * `pump_field` - Applied pump field amplitude [dimensionless coupling·field]
    ///
    /// # Returns
    /// Growth rate \[1/s\]. Positive means amplification.
    pub fn growth_rate(&self, pump_field: f64) -> f64 {
        let omega_k = self.magnon_frequency();
        self.coupling_strength * pump_field / 2.0 - self.damping * omega_k
    }
}

// ============================================================================
// Gross-Pitaevskii Equation Solver
// ============================================================================

/// A complex number representation for the condensate wavefunction.
///
/// ψ = re + i·im
#[derive(Debug, Clone, Copy)]
pub struct Complex {
    /// Real part
    pub re: f64,
    /// Imaginary part
    pub im: f64,
}

impl Complex {
    /// Create a new complex number.
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// Create from polar form: r·e^(iθ)
    pub fn from_polar(r: f64, theta: f64) -> Self {
        Self {
            re: r * theta.cos(),
            im: r * theta.sin(),
        }
    }

    /// Squared modulus |ψ|²
    pub fn norm_sq(&self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// Modulus |ψ|
    pub fn norm(&self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// Phase angle arg(ψ)
    pub fn phase(&self) -> f64 {
        self.im.atan2(self.re)
    }

    /// Complex conjugate
    pub fn conj(&self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    /// Multiplication
    pub fn mul(&self, other: &Self) -> Self {
        Self {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }

    /// Scalar multiplication
    pub fn scale(&self, s: f64) -> Self {
        Self {
            re: self.re * s,
            im: self.im * s,
        }
    }

    /// Addition
    pub fn add(&self, other: &Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }

    /// Subtraction
    pub fn sub(&self, other: &Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }

    /// Multiply by i
    pub fn mul_i(&self) -> Self {
        Self {
            re: -self.im,
            im: self.re,
        }
    }
}

/// 1D Gross-Pitaevskii equation solver for magnon condensate dynamics.
///
/// Solves:
///   iℏ ∂ψ/∂t = (-ℏ²/(2m*) ∂²ψ/∂x² + V(x) + g|ψ|² - μ) ψ
///
/// using explicit Euler time stepping on a uniform spatial grid.
#[derive(Debug, Clone)]
pub struct GrossPitaevskiiSolver {
    /// Condensate wavefunction on the spatial grid
    pub psi: Vec<Complex>,

    /// External potential on the grid \[J\]
    pub potential: Vec<f64>,

    /// Grid spacing \[m\]
    pub dx: f64,

    /// Number of grid points
    pub n_grid: usize,

    /// Effective magnon mass \[kg\]
    pub effective_mass: f64,

    /// Magnon-magnon interaction strength \[J·m\]
    /// (in 1D, the units differ from 3D)
    pub interaction_g: f64,

    /// Chemical potential \[J\]
    pub chemical_potential: f64,

    /// Time step \[s\]
    pub dt: f64,

    /// Current simulation time \[s\]
    pub time: f64,
}

impl GrossPitaevskiiSolver {
    /// Create a new GP solver on a uniform 1D grid.
    ///
    /// # Arguments
    /// * `n_grid` - Number of spatial grid points
    /// * `dx` - Grid spacing \[m\]
    /// * `dt` - Time step \[s\]
    /// * `effective_mass` - Effective magnon mass \[kg\]
    /// * `interaction_g` - Interaction strength \[J·m\] (1D)
    /// * `chemical_potential` - Chemical potential \[J\]
    ///
    /// # Errors
    /// Returns error for invalid parameters.
    pub fn new(
        n_grid: usize,
        dx: f64,
        dt: f64,
        effective_mass: f64,
        interaction_g: f64,
        chemical_potential: f64,
    ) -> Result<Self> {
        if n_grid < 3 {
            return Err(error::invalid_param(
                "n_grid",
                "must be at least 3 for finite differences",
            ));
        }
        if dx <= 0.0 {
            return Err(error::invalid_param("dx", "must be positive"));
        }
        if dt <= 0.0 {
            return Err(error::invalid_param("dt", "must be positive"));
        }
        if effective_mass <= 0.0 {
            return Err(error::invalid_param("effective_mass", "must be positive"));
        }

        Ok(Self {
            psi: vec![Complex::new(0.0, 0.0); n_grid],
            potential: vec![0.0; n_grid],
            dx,
            n_grid,
            effective_mass,
            interaction_g,
            chemical_potential,
            dt,
            time: 0.0,
        })
    }

    /// Initialize the wavefunction with a uniform condensate.
    ///
    /// ψ(x) = √n₀ (real, uniform)
    ///
    /// # Arguments
    /// * `density` - Condensate density \[1/m\] (in 1D)
    pub fn set_uniform_condensate(&mut self, density: f64) {
        let amplitude = density.abs().sqrt();
        for psi_i in &mut self.psi {
            *psi_i = Complex::new(amplitude, 0.0);
        }
    }

    /// Initialize the wavefunction with a Gaussian profile.
    ///
    /// ψ(x) = A · exp(-(x - x₀)² / (2σ²)) · exp(i k₀ x)
    ///
    /// Normalized so that ∫|ψ|² dx = total_particles
    ///
    /// # Arguments
    /// * `center` - Center position index
    /// * `sigma` - Width in grid spacings
    /// * `total_particles` - Total number of particles
    /// * `momentum` - Initial momentum k₀ \[1/m\]
    pub fn set_gaussian_condensate(
        &mut self,
        center: f64,
        sigma: f64,
        total_particles: f64,
        momentum: f64,
    ) {
        // First pass: compute unnormalized Gaussian
        let mut norm_sq_sum = 0.0;
        for i in 0..self.n_grid {
            let x = i as f64 * self.dx;
            let x0 = center * self.dx;
            let gaussian =
                (-((x - x0) * (x - x0)) / (2.0 * sigma * sigma * self.dx * self.dx)).exp();
            let phase = momentum * x;
            self.psi[i] = Complex::from_polar(gaussian, phase);
            norm_sq_sum += gaussian * gaussian * self.dx;
        }

        // Normalize
        if norm_sq_sum > 0.0 {
            let scale = (total_particles / norm_sq_sum).sqrt();
            for psi_i in &mut self.psi {
                *psi_i = psi_i.scale(scale);
            }
        }
    }

    /// Compute the total number of particles: N = ∫|ψ|² dx
    pub fn total_particles(&self) -> f64 {
        self.psi.iter().map(|psi_i| psi_i.norm_sq()).sum::<f64>() * self.dx
    }

    /// Compute the total energy of the condensate.
    ///
    /// E = ∫ [ℏ²/(2m*)|∇ψ|² + V|ψ|² + (g/2)|ψ|⁴] dx
    pub fn total_energy(&self) -> f64 {
        let kinetic_coeff = HBAR * HBAR / (2.0 * self.effective_mass);
        let mut energy = 0.0;

        for i in 0..self.n_grid {
            let psi_i = &self.psi[i];
            let density = psi_i.norm_sq();

            // Kinetic energy via finite difference gradient
            let grad = self.gradient_at(i);
            let grad_sq = grad.norm_sq();
            energy += kinetic_coeff * grad_sq;

            // Potential energy
            energy += self.potential[i] * density;

            // Interaction energy
            energy += 0.5 * self.interaction_g * density * density;
        }

        energy * self.dx
    }

    /// Compute the spatial gradient of ψ at grid point i using central differences.
    fn gradient_at(&self, i: usize) -> Complex {
        // Neumann-like boundary conditions
        let psi_prev = if i == 0 {
            &self.psi[0]
        } else {
            &self.psi[i - 1]
        };
        let psi_next = if i == self.n_grid - 1 {
            &self.psi[self.n_grid - 1]
        } else {
            &self.psi[i + 1]
        };

        psi_next.sub(psi_prev).scale(0.5 / self.dx)
    }

    /// Compute the Laplacian of ψ at grid point i.
    fn laplacian_at(&self, i: usize) -> Complex {
        let psi_i = &self.psi[i];

        // Neumann boundary conditions
        let psi_prev = if i == 0 { psi_i } else { &self.psi[i - 1] };
        let psi_next = if i == self.n_grid - 1 {
            psi_i
        } else {
            &self.psi[i + 1]
        };

        // (ψ_{i+1} - 2ψ_i + ψ_{i-1}) / dx²
        psi_next
            .add(psi_prev)
            .sub(&psi_i.scale(2.0))
            .scale(1.0 / (self.dx * self.dx))
    }

    /// Compute the right-hand side of the GPE: dψ/dt = -i/ℏ · H ψ
    ///
    /// H ψ = (-ℏ²/(2m*) ∇²ψ + V ψ + g|ψ|² ψ - μ ψ)
    fn compute_rhs(&self) -> Vec<Complex> {
        let kinetic_coeff = HBAR * HBAR / (2.0 * self.effective_mass);
        let inv_hbar = 1.0 / HBAR;
        let mut rhs = Vec::with_capacity(self.n_grid);

        for i in 0..self.n_grid {
            let psi_i = &self.psi[i];
            let density = psi_i.norm_sq();

            // Kinetic: -ℏ²/(2m*) ∇²ψ
            let kinetic = self.laplacian_at(i).scale(-kinetic_coeff);

            // Potential + interaction - chemical potential
            let local_potential =
                self.potential[i] + self.interaction_g * density - self.chemical_potential;
            let potential_term = psi_i.scale(local_potential);

            // H ψ = kinetic + potential_term
            let h_psi = kinetic.add(&potential_term);

            // dψ/dt = -i/ℏ · H ψ = (1/ℏ) · (H ψ rotated by -π/2)
            // -i · (a + ib) = b - ia
            let dpsi_dt = Complex::new(h_psi.im, -h_psi.re).scale(inv_hbar);

            rhs.push(dpsi_dt);
        }

        rhs
    }

    /// Advance the wavefunction by one time step using explicit Euler.
    ///
    /// ψ(t + dt) = ψ(t) + dt · dψ/dt
    ///
    /// Note: explicit Euler is only conditionally stable. The time step
    /// must satisfy dt < 2m* dx² / ℏ for stability of the kinetic term.
    pub fn step_euler(&mut self) {
        let rhs = self.compute_rhs();

        for (i, rhs_val) in rhs.iter().enumerate().take(self.n_grid) {
            self.psi[i] = self.psi[i].add(&rhs_val.scale(self.dt));
        }

        self.time += self.dt;
    }

    /// Advance the wavefunction by one time step using Heun's method (RK2).
    ///
    /// Provides better stability and accuracy than explicit Euler.
    pub fn step_heun(&mut self) {
        // k1 = f(t, ψ)
        let k1 = self.compute_rhs();

        // Save original ψ
        let psi_orig: Vec<Complex> = self.psi.clone();

        // Predictor: ψ_pred = ψ + dt·k1
        for (i, k1_val) in k1.iter().enumerate().take(self.n_grid) {
            self.psi[i] = self.psi[i].add(&k1_val.scale(self.dt));
        }

        // k2 = f(t + dt, ψ_pred)
        self.time += self.dt;
        let k2 = self.compute_rhs();

        // Corrector: ψ_new = ψ_orig + (dt/2)(k1 + k2)
        for i in 0..self.n_grid {
            let avg = k1[i].add(&k2[i]).scale(0.5 * self.dt);
            self.psi[i] = psi_orig[i].add(&avg);
        }
    }

    /// Get the density profile |ψ(x)|².
    pub fn density_profile(&self) -> Vec<f64> {
        self.psi.iter().map(|psi_i| psi_i.norm_sq()).collect()
    }

    /// Get the phase profile arg(ψ(x)).
    pub fn phase_profile(&self) -> Vec<f64> {
        self.psi.iter().map(|psi_i| psi_i.phase()).collect()
    }
}

// ============================================================================
// Magnon Supercurrent
// ============================================================================

/// Calculate the magnon supercurrent density at a given grid point.
///
/// J_m = (ℏ / (2 m*)) (ψ* ∇ψ - ψ ∇ψ*)
///     = (ℏ / m*) |ψ|² ∇φ
///
/// where φ is the condensate phase.
///
/// # Arguments
/// * `solver` - The GP solver containing the wavefunction
/// * `index` - Grid point index
///
/// # Returns
/// Magnon current density [1/(m·s)] at the given point.
///
/// # Errors
/// Returns error if index is out of bounds.
pub fn magnon_supercurrent(solver: &GrossPitaevskiiSolver, index: usize) -> Result<f64> {
    if index >= solver.n_grid {
        return Err(error::invalid_param("index", "out of bounds for the grid"));
    }

    let psi_i = &solver.psi[index];
    let grad = solver.gradient_at(index);

    // J = (ℏ / m*) Im(ψ* ∇ψ)
    let psi_conj = psi_i.conj();
    let product = psi_conj.mul(&grad);
    let current = HBAR / solver.effective_mass * product.im;

    Ok(current)
}

/// Calculate the magnon supercurrent profile across the entire grid.
///
/// Returns the supercurrent at each grid point.
pub fn supercurrent_profile(solver: &GrossPitaevskiiSolver) -> Vec<f64> {
    (0..solver.n_grid)
        .map(|i| magnon_supercurrent(solver, i).unwrap_or(0.0))
        .collect()
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_EFFECTIVE_MASS: f64 = 1e-28; // ~100 electron masses
    const TEST_DENSITY: f64 = 1e24; // magnons/m³
    const TEST_INTERACTION: f64 = 1e-50; // J·m³

    #[test]
    fn test_bec_temperature_is_positive_and_physical() {
        let t_bec = bec_temperature(TEST_DENSITY, TEST_EFFECTIVE_MASS)
            .expect("bec_temperature should succeed");
        assert!(t_bec > 0.0, "BEC temperature must be positive");
        // For these parameters, T_BEC should be order of room temperature or below
        assert!(
            t_bec < 1e6,
            "BEC temperature should be physically reasonable, got {}",
            t_bec
        );
    }

    #[test]
    fn test_chemical_potential_below_band_minimum() {
        // Above T_BEC, chemical potential should be negative (below band minimum at 0)
        let condensate = MagnonCondensate::new(
            500.0, // well above any reasonable T_BEC for these params
            TEST_DENSITY,
            TEST_EFFECTIVE_MASS,
            TEST_INTERACTION,
        )
        .expect("MagnonCondensate::new should succeed");

        // If above T_BEC, μ should be ≤ 0 (band minimum)
        if condensate.condensate_fraction == 0.0 {
            assert!(
                condensate.chemical_potential <= 0.0,
                "Chemical potential must be ≤ ℏω_min (0 in our frame)"
            );
        }

        // Below T_BEC, μ is pinned at 0
        let condensate_cold = MagnonCondensate::new(
            0.001, // very cold
            TEST_DENSITY,
            TEST_EFFECTIVE_MASS,
            TEST_INTERACTION,
        )
        .expect("MagnonCondensate::new should succeed for cold system");
        assert!(
            condensate_cold.chemical_potential.abs() < 1e-30,
            "Chemical potential should be pinned at band minimum below T_BEC"
        );
    }

    #[test]
    fn test_condensate_fraction_approaches_one_at_zero_temperature() {
        let condensate = MagnonCondensate::new(
            0.0, // T = 0
            TEST_DENSITY,
            TEST_EFFECTIVE_MASS,
            TEST_INTERACTION,
        )
        .expect("MagnonCondensate::new should succeed at T=0");

        assert!(
            (condensate.condensate_fraction - 1.0).abs() < 1e-12,
            "Condensate fraction must approach 1 as T -> 0, got {}",
            condensate.condensate_fraction
        );
    }

    #[test]
    fn test_condensate_fraction_zero_above_t_bec() {
        let t_bec = bec_temperature(TEST_DENSITY, TEST_EFFECTIVE_MASS)
            .expect("bec_temperature should succeed");

        let condensate = MagnonCondensate::new(
            t_bec * 2.0,
            TEST_DENSITY,
            TEST_EFFECTIVE_MASS,
            TEST_INTERACTION,
        )
        .expect("MagnonCondensate::new should succeed");

        assert!(
            condensate.condensate_fraction.abs() < 1e-12,
            "Condensate fraction must be zero above T_BEC"
        );
    }

    #[test]
    fn test_magnon_distribution_classical_limit() {
        // At high T, Bose-Einstein → classical (Boltzmann)
        let omega = 1e12; // 1 THz magnon
        let temperature = 1000.0; // high T
        let mu = -1e-20; // well below ℏω

        let _n_be = magnon_distribution(omega, temperature, mu)
            .expect("magnon_distribution should succeed");
        let _n_cl = classical_distribution(omega, temperature, mu)
            .expect("classical_distribution should succeed");

        // At high T they should converge (relative error < 10%)
        // For very high T, ℏω/(k_B T) is small, so BE ≈ k_BT/(ℏω - μ) >> 1
        // and classical = exp(-(ℏω - μ)/(k_BT)) which is close to 1
        // The classical limit holds when (ℏω - μ) >> k_BT, so use larger μ gap
        let omega_large = 1e14;
        let mu_far = -1e-18;
        let n_be2 = magnon_distribution(omega_large, temperature, mu_far)
            .expect("magnon_distribution should succeed");
        let n_cl2 = classical_distribution(omega_large, temperature, mu_far)
            .expect("classical_distribution should succeed");

        if n_cl2 > 1e-100 {
            let relative_diff = ((n_be2 - n_cl2) / n_cl2).abs();
            assert!(
                relative_diff < 0.1,
                "BE should approach classical at high T with large energy gap, diff = {}",
                relative_diff
            );
        }
    }

    #[test]
    fn test_suhl_threshold_is_positive() {
        let pumping = ParametricPumping::new(
            2.0 * PI * 10.0e9, // 10 GHz pump
            0.001,             // damping
            1e-24,             // coupling strength
        )
        .expect("ParametricPumping::new should succeed");

        let threshold = pumping.suhl_threshold();
        assert!(
            threshold > 0.0,
            "Suhl instability threshold must be positive"
        );
    }

    #[test]
    fn test_parametric_pumping_frequency() {
        let omega_p = 2.0 * PI * 10.0e9;
        let pumping = ParametricPumping::new(omega_p, 0.001, 1e-24)
            .expect("ParametricPumping::new should succeed");

        let omega_k = pumping.magnon_frequency();
        assert!(
            (omega_k - omega_p / 2.0).abs() < 1e-10,
            "Parametric magnons should be at half the pump frequency"
        );
    }

    #[test]
    fn test_growth_rate_above_threshold() {
        let pumping = ParametricPumping::new(2.0 * PI * 10.0e9, 0.001, 1e-24)
            .expect("ParametricPumping::new should succeed");

        let threshold = pumping.suhl_threshold();

        // Below threshold: negative growth rate
        let rate_below = pumping.growth_rate(threshold * 0.5);
        assert!(
            rate_below < 0.0,
            "Growth rate should be negative below threshold"
        );

        // Above threshold: positive growth rate
        let rate_above = pumping.growth_rate(threshold * 2.0);
        assert!(
            rate_above > 0.0,
            "Growth rate should be positive above threshold"
        );
    }

    #[test]
    fn test_magnon_supercurrent_from_phase_gradient() {
        // Create a condensate with a linear phase gradient -> supercurrent
        let n_grid = 64;
        let dx = 1e-6; // 1 μm
        let dt = 1e-15;
        let m_star = 1e-28;
        let g = 1e-38;
        let mu = 0.0;

        let mut solver = GrossPitaevskiiSolver::new(n_grid, dx, dt, m_star, g, mu)
            .expect("GP solver creation should succeed");

        // Set ψ = √n₀ · exp(i k x) with constant density and linear phase
        // Choose k0 small enough that k0*dx << 1 for accurate finite differences:
        // error ~ (k0*dx)^2/6, so k0*dx = 0.01 gives ~0.002% error
        let n0: f64 = 1e18; // density per meter (1D)
        let k0: f64 = 1e4; // wavevector (k0*dx = 0.01)
        let amplitude = n0.sqrt();

        for i in 0..n_grid {
            let x = i as f64 * dx;
            solver.psi[i] = Complex::from_polar(amplitude, k0 * x);
        }

        // Expected supercurrent: J = (ℏ k₀ / m*) n₀
        let expected_current = HBAR * k0 / m_star * n0;

        // Check current at interior points (away from boundaries)
        let j_mid = magnon_supercurrent(&solver, n_grid / 2)
            .expect("supercurrent calculation should succeed");

        let relative_error = ((j_mid - expected_current) / expected_current).abs();
        assert!(
            relative_error < 0.01,
            "Supercurrent should match ℏk/m* × n₀, relative error = {}",
            relative_error
        );
    }

    #[test]
    fn test_gross_pitaevskii_conserves_particle_number() {
        let n_grid = 128;
        let dx = 1e-7;
        // Use a very small dt for stability of explicit Euler
        let dt = 1e-18;
        let m_star = 1e-28;
        let g = 1e-40;
        let mu = 0.0;

        let mut solver = GrossPitaevskiiSolver::new(n_grid, dx, dt, m_star, g, mu)
            .expect("GP solver creation should succeed");

        // Initialize with a Gaussian wavepacket
        solver.set_gaussian_condensate(n_grid as f64 / 2.0, 10.0, 1000.0, 0.0);

        let n_initial = solver.total_particles();
        assert!(
            n_initial > 0.0,
            "Initial particle number should be positive"
        );

        // Evolve for several steps using Heun (better conservation)
        for _ in 0..50 {
            solver.step_heun();
        }

        let n_final = solver.total_particles();

        // Particle number should be approximately conserved
        // (explicit methods have some drift, but should be small for small dt)
        let relative_change = ((n_final - n_initial) / n_initial).abs();
        assert!(
            relative_change < 0.01,
            "Particle number should be approximately conserved, relative change = {}",
            relative_change
        );
    }

    #[test]
    fn test_critical_density_increases_with_temperature() {
        let m_star = 1e-28;

        let nc_100 =
            critical_density(m_star, 100.0).expect("critical_density should succeed at 100K");
        let nc_200 =
            critical_density(m_star, 200.0).expect("critical_density should succeed at 200K");

        assert!(
            nc_200 > nc_100,
            "Critical density should increase with temperature"
        );
    }

    #[test]
    fn test_healing_length_and_sound_speed() {
        let condensate =
            MagnonCondensate::new(0.001, TEST_DENSITY, TEST_EFFECTIVE_MASS, TEST_INTERACTION)
                .expect("MagnonCondensate::new should succeed");

        let xi = condensate
            .healing_length()
            .expect("healing_length should succeed");
        assert!(xi > 0.0, "Healing length must be positive");

        let cs = condensate
            .sound_speed()
            .expect("sound_speed should succeed");
        assert!(cs > 0.0, "Sound speed must be positive");
    }

    #[test]
    fn test_invalid_parameters_return_errors() {
        // Negative temperature
        assert!(
            MagnonCondensate::new(-1.0, TEST_DENSITY, TEST_EFFECTIVE_MASS, TEST_INTERACTION)
                .is_err()
        );
        // Zero density
        assert!(MagnonCondensate::new(1.0, 0.0, TEST_EFFECTIVE_MASS, TEST_INTERACTION).is_err());
        // Negative mass
        assert!(MagnonCondensate::new(1.0, TEST_DENSITY, -1.0, TEST_INTERACTION).is_err());
        // Invalid distribution: negative temperature
        assert!(magnon_distribution(1e12, -1.0, 0.0).is_err());
        // Invalid pumping: zero frequency
        assert!(ParametricPumping::new(0.0, 0.001, 1e-24).is_err());
    }
}
