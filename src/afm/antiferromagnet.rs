//! Antiferromagnetic Dynamics - THz Spintronics
//!
//! Unlike ferromagnets, antiferromagnets (AFMs) have two sublattices (A and B)
//! with antiparallel alignment. This leads to THz-order dynamics due to the
//! strong exchange coupling between sublattices.
//!
//! Key features:
//! - THz frequency magnetization dynamics
//! - Neel vector as order parameter
//! - No stray fields (ideal for high-density devices)
//! - Ultrafast switching potential

use crate::constants::GAMMA;
use crate::vector3::Vector3;

/// Two-sublattice antiferromagnet
///
/// Models coupled LLG equations for sublattices A and B:
/// dm_A/dt = -γ m_A × (H_ext + H_anis - J_ex m_B) + damping
/// dm_B/dt = -γ m_B × (H_ext + H_anis - J_ex m_A) + damping
#[derive(Debug, Clone)]
pub struct Antiferromagnet {
    /// Sublattice A magnetization (normalized)
    pub m_a: Vector3<f64>,

    /// Sublattice B magnetization (normalized)
    pub m_b: Vector3<f64>,

    /// Inter-sublattice exchange field strength \[A/m\]
    /// Typically 100-1000 T equivalent for strong AFMs like NiO, MnF2
    pub h_exchange: f64,

    /// Gyromagnetic ratio \[rad/(s·T)\]
    pub gamma: f64,

    /// Gilbert damping
    pub alpha: f64,

    /// Anisotropy field \[A/m\]
    pub h_anisotropy: f64,
}

impl Antiferromagnet {
    /// Create a new antiferromagnet with default parameters
    ///
    /// Initial state: m_A = +z, m_B = -z (antiparallel)
    pub fn new(h_exchange: f64, alpha: f64, h_anisotropy: f64) -> Self {
        Self {
            m_a: Vector3::new(0.0, 0.0, 1.0),
            m_b: Vector3::new(0.0, 0.0, -1.0),
            h_exchange,
            gamma: GAMMA,
            alpha,
            h_anisotropy,
        }
    }

    /// Create NiO-like antiferromagnet
    ///
    /// NiO is a canonical AFM with strong exchange coupling
    pub fn nio() -> Self {
        Self::new(
            8.0e7, // ~100 T equivalent exchange field
            0.005, // Low damping
            1.0e6, // Anisotropy field
        )
    }

    /// Create MnF2-like antiferromagnet
    pub fn mnf2() -> Self {
        Self::new(
            5.0e7, // Strong exchange
            0.01, 5.0e5,
        )
    }

    /// Calculate Neel vector: n = (m_A - m_B) / 2
    ///
    /// The Neel vector is the order parameter for antiferromagnets
    pub fn neel_vector(&self) -> Vector3<f64> {
        (self.m_a + self.m_b * (-1.0)) * 0.5
    }

    /// Calculate total magnetization: m = (m_A + m_B) / 2
    ///
    /// Ideally zero for perfect AFM, but can be non-zero during dynamics
    pub fn total_magnetization(&self) -> Vector3<f64> {
        (self.m_a + self.m_b) * 0.5
    }

    /// Evolve AFM dynamics by one timestep using Heun's method
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[A/m\]
    /// * `dt` - Time step \[s\]
    pub fn evolve(&mut self, h_ext: Vector3<f64>, dt: f64) {
        // First stage: calculate initial derivatives
        let (dma_dt1, dmb_dt1) = self.calc_derivatives(h_ext);

        // Predictor step
        let m_a_pred = (self.m_a + dma_dt1 * dt).normalize();
        let m_b_pred = (self.m_b + dmb_dt1 * dt).normalize();

        // Save original
        let m_a_orig = self.m_a;
        let m_b_orig = self.m_b;

        // Temporarily update for second stage
        self.m_a = m_a_pred;
        self.m_b = m_b_pred;

        // Second stage: calculate derivatives at predicted point
        let (dma_dt2, dmb_dt2) = self.calc_derivatives(h_ext);

        // Corrector step: average of two derivatives
        self.m_a = (m_a_orig + (dma_dt1 + dma_dt2) * (0.5 * dt)).normalize();
        self.m_b = (m_b_orig + (dmb_dt1 + dmb_dt2) * (0.5 * dt)).normalize();
    }

    /// Calculate time derivatives of both sublattices
    fn calc_derivatives(&self, h_ext: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
        // Effective fields including inter-sublattice exchange
        // H_eff_A = H_ext + H_anis - J_ex * m_B
        // H_eff_B = H_ext + H_anis - J_ex * m_A

        let h_anis_a = Vector3::new(0.0, 0.0, self.h_anisotropy * self.m_a.z);
        let h_anis_b = Vector3::new(0.0, 0.0, self.h_anisotropy * self.m_b.z);

        let h_eff_a = h_ext + h_anis_a + self.m_b * (-self.h_exchange);
        let h_eff_b = h_ext + h_anis_b + self.m_a * (-self.h_exchange);

        // LLG equation for each sublattice
        let dma_dt = self.calc_llg(self.m_a, h_eff_a);
        let dmb_dt = self.calc_llg(self.m_b, h_eff_b);

        (dma_dt, dmb_dt)
    }

    /// Calculate LLG torque for a single magnetization
    fn calc_llg(&self, m: Vector3<f64>, h: Vector3<f64>) -> Vector3<f64> {
        let precession = m.cross(&h);
        let damping = m.cross(&precession) * self.alpha;

        (precession + damping) * (-self.gamma / (1.0 + self.alpha * self.alpha))
    }

    /// Calculate AFM resonance frequency \[Hz\]
    ///
    /// For uniaxial AFM: f_res ≈ γ μ₀ sqrt(H_E * H_A) / (2π)
    /// where H_E is exchange field and H_A is anisotropy field (both in A/m)
    pub fn resonance_frequency(&self) -> f64 {
        let mu0 = 1.256637e-6; // H/m
        let omega = self.gamma * mu0 * (self.h_exchange * self.h_anisotropy).sqrt();
        omega / (2.0 * std::f64::consts::PI)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_afm_creation() {
        let afm = Antiferromagnet::nio();
        assert!(afm.h_exchange > 0.0);
        assert!((afm.m_a.z - 1.0).abs() < 1e-10);
        assert!((afm.m_b.z + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_neel_vector() {
        let afm = Antiferromagnet::nio();
        let neel = afm.neel_vector();

        // Neel vector should be along +z for initial state
        assert!((neel.z - 1.0).abs() < 1e-10);
        assert!((neel.x).abs() < 1e-10);
        assert!((neel.y).abs() < 1e-10);
    }

    #[test]
    fn test_total_magnetization_zero() {
        let afm = Antiferromagnet::nio();
        let m_total = afm.total_magnetization();

        // Perfect AFM should have zero total magnetization
        assert!(m_total.magnitude() < 1e-10);
    }

    #[test]
    fn test_resonance_frequency_thz() {
        let afm = Antiferromagnet::nio();
        let f_res = afm.resonance_frequency();

        // AFM resonance should be in THz range (> 100 GHz)
        assert!(f_res > 1.0e11); // > 100 GHz
        assert!(f_res < 1.0e14); // < 100 THz (reasonable upper bound for strong AFMs)
    }

    #[test]
    fn test_afm_evolution() {
        let mut afm = Antiferromagnet::nio();
        let h_ext = Vector3::new(1.0e6, 0.0, 0.0); // 1 T-equivalent field

        let neel_before = afm.neel_vector();

        // Evolve for a few steps
        for _ in 0..10 {
            afm.evolve(h_ext, 1.0e-14); // 10 fs timestep for THz dynamics
        }

        let neel_after = afm.neel_vector();

        // Neel vector should have changed due to external field
        assert!((neel_before - neel_after).magnitude() > 1e-6);

        // But magnetizations should remain normalized
        assert!((afm.m_a.magnitude() - 1.0).abs() < 1e-6);
        assert!((afm.m_b.magnitude() - 1.0).abs() < 1e-6);
    }
}
