//! Cavity Magnonics - Hybrid Magnon-Photon Systems
//!
//! Strong coupling between magnons (spin waves) and microwave photons in
//! cavities leads to formation of hybrid quasiparticles: magnon-polaritons.
//!
//! This enables coherent information transfer between spin and photonic systems.

use crate::constants::GAMMA;
use crate::vector3::Vector3;

/// Coupled magnon-photon system (Jaynes-Cummings-like model)
///
/// The system consists of:
/// - **Magnon mode**: Ferromagnetic resonance (FMR) at Kittel frequency
/// - **Photon mode**: Microwave cavity resonance
/// - **Coupling**: Magnetic dipole interaction
///
/// Hamiltonian (simplified):
/// H = ℏω_m a†a + ℏω_c b†b + ℏg(a†b + ab†)
///
/// where a (magnon), b (photon), g (coupling strength)
#[derive(Debug, Clone)]
pub struct HybridSystem {
    /// Magnetization vector (normalized macro-spin)
    pub magnetization: Vector3<f64>,

    /// Cavity electric field amplitude (normalized)
    /// Represents photon field
    pub cavity_amplitude: f64,

    /// Cavity phase
    pub cavity_phase: f64,

    /// Magnon-photon coupling strength \[rad/s\]
    /// Typical: 1-100 MHz for YIG sphere in cavity
    pub coupling_g: f64,

    /// Cavity resonance frequency \[rad/s\]
    pub omega_cavity: f64,

    /// Magnon (FMR) frequency \[rad/s\]
    pub omega_magnon: f64,

    /// Cavity photon linewidth \[rad/s\]
    pub kappa: f64,

    /// Magnon linewidth (related to Gilbert damping) \[rad/s\]
    pub gamma_m: f64,

    /// Static bias field \[A/m\]
    pub h_bias: Vector3<f64>,

    /// Gyromagnetic ratio
    pub gamma: f64,

    /// Gilbert damping
    pub alpha: f64,
}

impl HybridSystem {
    /// Create a new hybrid magnon-photon system
    ///
    /// # Arguments
    /// * `h_bias` - Static bias field \[A/m\]
    /// * `coupling_g` - Coupling strength \[rad/s\]
    /// * `omega_cavity` - Cavity frequency \[rad/s\]
    /// * `kappa` - Cavity linewidth \[rad/s\]
    /// * `alpha` - Gilbert damping
    pub fn new(
        h_bias: Vector3<f64>,
        coupling_g: f64,
        omega_cavity: f64,
        kappa: f64,
        alpha: f64,
    ) -> Self {
        // Calculate FMR frequency from bias field (H in A/m → B in T)
        let mu0 = 1.256637e-6; // H/m
        let h_bias_magnitude = h_bias.magnitude();
        let omega_magnon = GAMMA * mu0 * h_bias_magnitude;

        Self {
            magnetization: Vector3::new(0.0, 0.0, 1.0), // Initial: aligned with bias
            cavity_amplitude: 0.0,
            cavity_phase: 0.0,
            coupling_g,
            omega_cavity,
            omega_magnon,
            kappa,
            gamma_m: alpha * omega_magnon, // Magnon linewidth
            h_bias,
            gamma: GAMMA,
            alpha,
        }
    }

    /// Create YIG-sphere in microwave cavity (typical experiment)
    ///
    /// Parameters based on canonical cavity magnonics experiments
    pub fn yig_cavity() -> Self {
        let h_bias = Vector3::new(0.0, 0.0, 8.0e4); // 100 mT ≈ 1000 Oe
        let omega_cavity = 2.0 * std::f64::consts::PI * 10.0e9; // 10 GHz
        let coupling_g = 2.0 * std::f64::consts::PI * 100.0e6; // 100 MHz (stronger coupling)
        let kappa = 2.0 * std::f64::consts::PI * 1.0e6; // 1 MHz
        let alpha = 0.0001; // Ultra-low damping for YIG

        Self::new(h_bias, coupling_g, omega_cavity, kappa, alpha)
    }

    /// Evolve coupled system by one time step
    ///
    /// Uses rotating frame approximation and coupled mode theory
    ///
    /// # Arguments
    /// * `h_drive` - External microwave drive field \[A/m\]
    /// * `dt` - Time step \[s\]
    pub fn evolve(&mut self, h_drive: Vector3<f64>, dt: f64) {
        // 1. Photon mode dynamics (harmonic oscillator with drive and coupling)
        //    da/dt = -i ω_c a - κa/2 - i g m_+ + drive
        //
        // Simplified: treat cavity amplitude and phase separately
        let cavity_drive = h_drive.x; // Drive amplitude

        // Magnon contribution to cavity (back-action)
        let m_plus = self.magnetization.x + self.magnetization.y; // Lowering operator component
        let magnon_to_cavity = self.coupling_g * m_plus;

        // Cavity amplitude evolution (envelope)
        let da_dt =
            -self.kappa * 0.5 * self.cavity_amplitude + magnon_to_cavity + cavity_drive * 1.0e9; // Scale drive

        self.cavity_amplitude += da_dt * dt;

        // Cavity phase evolution
        let dphi_dt = self.omega_cavity;
        self.cavity_phase += dphi_dt * dt;

        // 2. Magnon mode dynamics (LLG with cavity feedback)
        //    Cavity provides oscillating transverse field
        let h_cavity_feedback = Vector3::new(
            self.cavity_amplitude * (self.cavity_phase).cos() * self.coupling_g / self.gamma,
            self.cavity_amplitude * (self.cavity_phase).sin() * self.coupling_g / self.gamma,
            0.0,
        );

        let h_total = self.h_bias + h_cavity_feedback + h_drive;

        // Standard LLG
        let m_cross_h = self.magnetization.cross(&h_total);
        let damping = self.magnetization.cross(&m_cross_h) * self.alpha;

        let dm_dt = (m_cross_h + damping) * (-self.gamma / (1.0 + self.alpha * self.alpha));

        self.magnetization = (self.magnetization + dm_dt * dt).normalize();
    }

    /// Calculate detuning between cavity and magnon modes \[rad/s\]
    pub fn detuning(&self) -> f64 {
        self.omega_cavity - self.omega_magnon
    }

    /// Calculate cooperativity parameter
    ///
    /// C = 4g² / (κ γ_m)
    ///
    /// C >> 1: strong coupling regime (can observe level splitting)
    /// C ~ 1: intermediate coupling
    /// C << 1: weak coupling (perturbative regime)
    pub fn cooperativity(&self) -> f64 {
        4.0 * self.coupling_g.powi(2) / (self.kappa * self.gamma_m)
    }

    /// Calculate normal mode splitting (vacuum Rabi splitting) \[rad/s\]
    ///
    /// At zero detuning: Ω_± = ±g
    ///
    /// Returns the frequency splitting between hybridized modes
    pub fn rabi_splitting(&self) -> f64 {
        2.0 * self.coupling_g
    }

    /// Check if system is in strong coupling regime
    pub fn is_strong_coupling(&self) -> bool {
        self.cooperativity() > 1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hybrid_creation() {
        let hybrid = HybridSystem::yig_cavity();

        assert!(hybrid.coupling_g > 0.0);
        assert!(hybrid.omega_cavity > 0.0);
        assert!(hybrid.omega_magnon > 0.0);
    }

    #[test]
    fn test_yig_cavity_strong_coupling() {
        let hybrid = HybridSystem::yig_cavity();

        // YIG in cavity should be in strong coupling regime
        assert!(hybrid.is_strong_coupling());
        assert!(hybrid.cooperativity() > 10.0); // Typical C ~ 100-1000
    }

    #[test]
    fn test_detuning_calculation() {
        let mut hybrid = HybridSystem::yig_cavity();

        // Adjust magnon frequency by changing bias field
        hybrid.h_bias = hybrid.h_bias * 1.1; // +10% field
        hybrid.omega_magnon = GAMMA * hybrid.h_bias.magnitude();

        let detuning = hybrid.detuning();
        assert!(detuning.abs() > 0.0);
    }

    #[test]
    fn test_rabi_splitting() {
        let hybrid = HybridSystem::yig_cavity();
        let splitting = hybrid.rabi_splitting();

        // Splitting should be 2g
        assert!((splitting - 2.0 * hybrid.coupling_g).abs() < 1e-6);
    }

    #[test]
    fn test_evolution() {
        let mut hybrid = HybridSystem::yig_cavity();

        let h_drive = Vector3::new(100.0, 0.0, 0.0); // Small drive

        let m_initial = hybrid.magnetization;

        // Evolve
        for _ in 0..100 {
            hybrid.evolve(h_drive, 1.0e-12); // 1 ps timestep
        }

        // Magnetization should have changed
        assert!((hybrid.magnetization - m_initial).magnitude() > 1e-6);

        // But should remain normalized
        assert!((hybrid.magnetization.magnitude() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_cavity_amplitude_growth() {
        let mut hybrid = HybridSystem::yig_cavity();

        // Apply strong drive
        let h_drive = Vector3::new(1.0e4, 0.0, 0.0);

        for _ in 0..1000 {
            hybrid.evolve(h_drive, 1.0e-12);
        }

        // Cavity amplitude should have grown from zero
        assert!(hybrid.cavity_amplitude.abs() > 1e-6);
    }

    #[test]
    fn test_magnon_frequency() {
        let h_bias = Vector3::new(0.0, 0.0, 8.0e4); // ~100 mT
        let mu0 = 1.256637e-6;
        let omega_expected = GAMMA * mu0 * 8.0e4;

        let hybrid = HybridSystem::new(h_bias, 1.0e8, 1.0e10, 1.0e6, 0.0001);

        assert!((hybrid.omega_magnon - omega_expected).abs() < 1e6);
    }
}
