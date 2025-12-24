//! WebAssembly bindings for browser-based spintronics simulations
//!
//! This module provides JavaScript-friendly interfaces to the spintronics library,
//! enabling real-time magnetization dynamics, spin current calculations, and
//! interactive physics simulations directly in the web browser.
//!
//! # Features
//!
//! - **Magnetization Dynamics**: Real-time LLG solver for single spins and spin chains
//! - **Spin Current Effects**: Spin Hall effect, spin pumping, spin Seebeck effect
//! - **Material Properties**: Access to comprehensive material database
//! - **Visualization Support**: Export data in JSON format for plotting
//!
//! # Usage
//!
//! Build the WASM package:
//! ```bash
//! wasm-pack build --features wasm --target web
//! ```
//!
//! Use in JavaScript:
//! ```javascript
//! import init, { SpinSimulator } from './pkg/spintronics.js';
//!
//! async function run() {
//!     await init();
//!     const sim = SpinSimulator.new();
//!     sim.step(0.1);  // Advance by 0.1 ns
//!     const mx = sim.get_mx();
//!     console.log(`Magnetization x: ${mx}`);
//! }
//! ```

use wasm_bindgen::prelude::*;

use crate::constants::GAMMA;
use crate::llg::LLGSolver;
use crate::material::Ferromagnet;
use crate::vector3::Vector3;

/// Initialize panic hook for better error messages in browser console
#[wasm_bindgen(start)]
pub fn main() {
    #[cfg(feature = "wasm")]
    console_error_panic_hook::set_once();
}

/// Log a message to the browser console
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

/// JavaScript-accessible Vector3 wrapper
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

#[wasm_bindgen]
impl Vec3 {
    /// Create a new 3D vector
    #[wasm_bindgen(constructor)]
    pub fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3 { x, y, z }
    }

    /// Get x component
    #[wasm_bindgen(getter)]
    pub fn x(&self) -> f64 {
        self.x
    }

    /// Get y component
    #[wasm_bindgen(getter)]
    pub fn y(&self) -> f64 {
        self.y
    }

    /// Get z component
    #[wasm_bindgen(getter)]
    pub fn z(&self) -> f64 {
        self.z
    }

    /// Calculate magnitude
    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Normalize the vector
    pub fn normalized(&self) -> Vec3 {
        let mag = self.magnitude();
        if mag > 1e-10 {
            Vec3 {
                x: self.x / mag,
                y: self.y / mag,
                z: self.z / mag,
            }
        } else {
            *self
        }
    }
}

impl From<Vector3<f64>> for Vec3 {
    fn from(v: Vector3<f64>) -> Self {
        Vec3 {
            x: v.x,
            y: v.y,
            z: v.z,
        }
    }
}

impl From<Vec3> for Vector3<f64> {
    fn from(v: Vec3) -> Self {
        Vector3::new(v.x, v.y, v.z)
    }
}

/// Single-spin magnetization dynamics simulator
///
/// Simulates the time evolution of a magnetic moment using the
/// Landau-Lifshitz-Gilbert (LLG) equation.
#[wasm_bindgen]
pub struct SpinSimulator {
    solver: LLGSolver,
    magnetization: Vector3<f64>,
    time: f64,
}

#[wasm_bindgen]
impl SpinSimulator {
    /// Create a new spin simulator with Permalloy material
    ///
    /// Initial magnetization is along +z direction
    #[wasm_bindgen(constructor)]
    pub fn new() -> SpinSimulator {
        Self::default()
    }

    /// Create simulator with custom material parameters
    ///
    /// # Arguments
    /// * `ms` - Saturation magnetization (A/m)
    /// * `alpha` - Gilbert damping constant
    /// * `gamma` - Gyromagnetic ratio (m/(A·s))
    pub fn with_material(ms: f64, alpha: f64, _gamma: f64) -> SpinSimulator {
        let material = Ferromagnet {
            ms,
            alpha,
            exchange_a: 1.3e-11, // Default Permalloy value
            anisotropy_k: 0.0,
            easy_axis: Vector3::new(0.0, 0.0, 1.0),
        };
        let magnetization = Vector3::new(0.0, 0.0, ms);
        let solver = LLGSolver::new(material);

        SpinSimulator {
            solver,
            magnetization,
            time: 0.0,
        }
    }

    /// Set external magnetic field (A/m)
    pub fn set_field(&mut self, hx: f64, hy: f64, hz: f64) {
        self.solver.h_ext = Vector3::new(hx, hy, hz);
    }

    /// Set magnetization direction (will be normalized)
    pub fn set_magnetization(&mut self, mx: f64, my: f64, mz: f64) {
        let m = Vector3::new(mx, my, mz);
        let m_normalized = m.normalize();
        self.magnetization = m_normalized * self.solver.material.ms;
    }

    /// Advance simulation by time step dt (in nanoseconds)
    ///
    /// Uses RK4 integration for accurate dynamics
    pub fn step(&mut self, dt_ns: f64) {
        let dt = dt_ns * 1e-9; // Convert ns to seconds
        self.magnetization = self.solver.step_rk4(self.magnetization, dt);
        self.time += dt;
    }

    /// Run multiple steps
    pub fn run_steps(&mut self, n_steps: usize, dt_ns: f64) {
        for _ in 0..n_steps {
            self.step(dt_ns);
        }
    }

    /// Get x component of magnetization (normalized)
    pub fn get_mx(&self) -> f64 {
        self.magnetization.x / self.solver.material.ms
    }

    /// Get y component of magnetization (normalized)
    pub fn get_my(&self) -> f64 {
        self.magnetization.y / self.solver.material.ms
    }

    /// Get z component of magnetization (normalized)
    pub fn get_mz(&self) -> f64 {
        self.magnetization.z / self.solver.material.ms
    }

    /// Get current simulation time (seconds)
    pub fn get_time(&self) -> f64 {
        self.time
    }

    /// Get magnetization as Vec3
    pub fn get_magnetization(&self) -> Vec3 {
        Vec3::from(self.magnetization)
    }

    /// Get normalized magnetization direction
    pub fn get_direction(&self) -> Vec3 {
        let m_norm = self.magnetization.normalize();
        Vec3::from(m_norm)
    }

    /// Calculate effective field (A/m)
    pub fn get_effective_field(&self) -> Vec3 {
        let h_eff = self.solver.effective_field(self.magnetization);
        Vec3::from(h_eff)
    }

    /// Reset simulation to initial state
    pub fn reset(&mut self) {
        self.magnetization = Vector3::new(0.0, 0.0, self.solver.material.ms);
        self.time = 0.0;
    }
}

impl Default for SpinSimulator {
    fn default() -> Self {
        let material = Ferromagnet::permalloy();
        let ms = material.ms;
        let solver = LLGSolver::new(material);
        let magnetization = Vector3::new(0.0, 0.0, ms);

        SpinSimulator {
            solver,
            magnetization,
            time: 0.0,
        }
    }
}

/// Spin Hall effect calculator
///
/// Calculates spin current generation from charge current via the spin Hall effect
#[wasm_bindgen]
pub struct SpinHallCalculator {
    theta_sh: f64,
}

#[wasm_bindgen]
impl SpinHallCalculator {
    /// Create calculator for Platinum
    #[wasm_bindgen(constructor)]
    pub fn new() -> SpinHallCalculator {
        Self::default()
    }

    /// Create calculator with custom spin Hall angle
    pub fn with_angle(theta_sh: f64) -> SpinHallCalculator {
        SpinHallCalculator { theta_sh }
    }

    /// Calculate spin current density (A/m²) from charge current density
    ///
    /// # Arguments
    /// * `jc` - Charge current density (A/m²)
    /// * `theta` - Angle of spin polarization (radians)
    /// * `phi` - Azimuthal angle (radians)
    ///
    /// Returns: Spin current magnitude
    pub fn spin_current(&self, jc: f64, _theta: f64, _phi: f64) -> f64 {
        // js = θ_SH * jc (simplified)
        self.theta_sh * jc
    }

    /// Get spin Hall angle
    pub fn get_angle(&self) -> f64 {
        self.theta_sh
    }
}

impl Default for SpinHallCalculator {
    fn default() -> Self {
        SpinHallCalculator {
            theta_sh: 0.07, // Platinum spin Hall angle
        }
    }
}

/// Spin chain simulator for magnon propagation
///
/// Simulates coupled spins in a 1D chain, useful for studying
/// magnon (spin wave) propagation
#[wasm_bindgen]
pub struct SpinChain {
    spins: Vec<Vector3<f64>>,
    material: Ferromagnet,
    exchange_j: f64,
    time: f64,
}

#[wasm_bindgen]
impl SpinChain {
    /// Create a new spin chain
    ///
    /// # Arguments
    /// * `n_spins` - Number of spins in the chain
    /// * `exchange_j` - Exchange coupling (J)
    #[wasm_bindgen(constructor)]
    pub fn new(n_spins: usize, exchange_j: f64) -> SpinChain {
        let material = Ferromagnet::permalloy();
        let spins = vec![Vector3::new(0.0, 0.0, material.ms); n_spins];

        SpinChain {
            spins,
            material,
            exchange_j,
            time: 0.0,
        }
    }

    /// Get number of spins
    pub fn length(&self) -> usize {
        self.spins.len()
    }

    /// Set spin at index i
    pub fn set_spin(&mut self, i: usize, mx: f64, my: f64, mz: f64) {
        if i < self.spins.len() {
            let m = Vector3::new(mx, my, mz);
            let m_normalized = m.normalize();
            self.spins[i] = m_normalized * self.material.ms;
        }
    }

    /// Get spin x-component at index i (normalized)
    pub fn get_mx(&self, i: usize) -> f64 {
        if i < self.spins.len() {
            self.spins[i].x / self.material.ms
        } else {
            0.0
        }
    }

    /// Get spin y-component at index i (normalized)
    pub fn get_my(&self, i: usize) -> f64 {
        if i < self.spins.len() {
            self.spins[i].y / self.material.ms
        } else {
            0.0
        }
    }

    /// Get spin z-component at index i (normalized)
    pub fn get_mz(&self, i: usize) -> f64 {
        if i < self.spins.len() {
            self.spins[i].z / self.material.ms
        } else {
            0.0
        }
    }

    /// Excite a spin wave at position
    ///
    /// Creates a localized perturbation for magnon propagation
    pub fn excite_wave(&mut self, center: usize, amplitude: f64, width: f64) {
        let n = self.spins.len();
        for i in 0..n {
            let x = i as f64;
            let x0 = center as f64;
            let dx = x - x0;
            let perturbation = amplitude * (-dx * dx / (2.0 * width * width)).exp();

            let theta = perturbation;
            let m = Vector3::new(theta.sin(), 0.0, theta.cos());
            self.spins[i] = m.normalize() * self.material.ms;
        }
    }

    /// Advance spin chain dynamics by one time step
    ///
    /// Uses simple Heun method for coupled LLG equations
    pub fn step(&mut self, dt_ns: f64) {
        let dt = dt_ns * 1e-9;
        let alpha = self.material.alpha;
        let n = self.spins.len();

        let mut dmdt = vec![Vector3::new(0.0, 0.0, 0.0); n];

        // Calculate derivatives
        for (i, dm) in dmdt.iter_mut().enumerate().take(n) {
            let m = self.spins[i] * (1.0 / self.material.ms); // Normalized

            // Exchange field from neighbors
            let mut h_ex = Vector3::new(0.0, 0.0, 0.0);
            if i > 0 {
                h_ex = h_ex + self.spins[i - 1];
            }
            if i < n - 1 {
                h_ex = h_ex + self.spins[i + 1];
            }
            h_ex = h_ex * (-2.0 * self.exchange_j / self.material.ms);

            // LLG equation
            let m_cross_h = m.cross(&h_ex);
            let m_cross_m_cross_h = m.cross(&m_cross_h);

            *dm = (m_cross_h * (-GAMMA) + m_cross_m_cross_h * (-GAMMA * alpha)) * self.material.ms;
        }

        // Update spins
        for (i, dm) in dmdt.iter().enumerate().take(n) {
            self.spins[i] = self.spins[i] + *dm * dt;
            // Renormalize
            let mag = self.spins[i].magnitude();
            if mag > 1e-10 {
                self.spins[i] = self.spins[i] * (self.material.ms / mag);
            }
        }

        self.time += dt;
    }

    /// Get current time (seconds)
    pub fn get_time(&self) -> f64 {
        self.time
    }

    /// Export all magnetizations as JSON string
    pub fn export_json(&self) -> String {
        let mut data = String::from("{\"spins\":[");
        for (i, spin) in self.spins.iter().enumerate() {
            if i > 0 {
                data.push(',');
            }
            let mx = spin.x / self.material.ms;
            let my = spin.y / self.material.ms;
            let mz = spin.z / self.material.ms;
            data.push_str(&format!("{{\"mx\":{},\"my\":{},\"mz\":{}}}", mx, my, mz));
        }
        data.push_str(&format!("],\"time\":{}}}", self.time));
        data
    }
}

#[cfg(test)]
mod tests {
    use wasm_bindgen_test::*;

    use super::*;

    #[wasm_bindgen_test]
    fn test_vec3_creation() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(v.x(), 1.0);
        assert_eq!(v.y(), 2.0);
        assert_eq!(v.z(), 3.0);
    }

    #[wasm_bindgen_test]
    fn test_vec3_magnitude() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        assert!((v.magnitude() - 5.0).abs() < 1e-10);
    }

    #[wasm_bindgen_test]
    fn test_spin_simulator_creation() {
        let sim = SpinSimulator::new();
        assert!((sim.get_mz() - 1.0).abs() < 1e-10);
    }

    #[wasm_bindgen_test]
    fn test_spin_simulator_step() {
        let mut sim = SpinSimulator::new();
        sim.set_field(1000.0, 0.0, 0.0); // Field in x
        sim.step(0.1); // 0.1 ns
                       // Magnetization should precess around x-axis
        assert!(sim.get_time() > 0.0);
    }

    #[wasm_bindgen_test]
    fn test_spin_hall_calculator() {
        let calc = SpinHallCalculator::new();
        let jc = 1e10; // 10^10 A/m²
        let js = calc.spin_current(jc, 0.0, 0.0);
        assert!(js > 0.0);
    }

    #[wasm_bindgen_test]
    fn test_spin_chain() {
        let chain = SpinChain::new(10, 1e-21);
        assert_eq!(chain.length(), 10);
        assert!((chain.get_mz(0) - 1.0).abs() < 1e-10);
    }

    #[wasm_bindgen_test]
    fn test_spin_chain_excitation() {
        let mut chain = SpinChain::new(20, 1e-21);
        chain.excite_wave(10, 0.5, 2.0);
        // Center should have tilted magnetization
        assert!(chain.get_mx(10).abs() > 0.1);
    }
}
