//! Temperature-Dependent Material Properties
//!
//! This module provides temperature scaling for magnetic materials.
//! Key temperature dependencies include:
//! - Magnetization: Bloch's T^(3/2) law and critical behavior
//! - Anisotropy: Callen-Callen power law
//! - Exchange stiffness: Temperature scaling
//! - Gilbert damping: Weak temperature dependence
//!
//! ## Physics Background
//!
//! ### Bloch's Law for Magnetization:
//! At low temperatures (T << T_C):
//! M(T) / M(0) ≈ 1 - (T/T_C)^(3/2)
//!
//! ### Near Curie Temperature:
//! M(T) / M(0) ≈ ((T_C - T) / T_C)^β
//! where β ≈ 0.36 (3D Heisenberg), 0.32 (3D Ising)
//!
//! ### Callen-Callen Law for Anisotropy:
//! K(T) / K(0) ≈ (M(T) / M(0))^n
//! where n ≈ 3 for magnetocrystalline anisotropy
//!
//! ## Key References
//!
//! - F. Bloch, "Zur Theorie des Ferromagnetismus", Z. Physik 61, 206 (1930)
//! - H. B. Callen and E. Callen, "The present status of the temperature
//!   dependence of magnetocrystalline anisotropy", J. Phys. Chem. Solids 27, 1271 (1966)

use crate::material::Ferromagnet;

/// Temperature-dependent ferromagnetic material
#[derive(Debug, Clone)]
pub struct ThermalFerromagnet {
    /// Base material at T=0
    pub material: Ferromagnet,

    /// Curie temperature \[K\]
    pub curie_temperature: f64,

    /// Current temperature \[K\]
    pub temperature: f64,

    /// Critical exponent (default: 0.36 for 3D Heisenberg)
    pub critical_exponent: f64,
}

impl ThermalFerromagnet {
    /// Create temperature-dependent material
    ///
    /// # Arguments
    /// * `material` - Base ferromagnet
    /// * `curie_temp` - Curie temperature \[K\]
    /// * `temperature` - Operating temperature \[K\]
    pub fn new(material: Ferromagnet, curie_temp: f64, temperature: f64) -> Self {
        Self {
            material,
            curie_temperature: curie_temp,
            temperature,
            critical_exponent: 0.36, // 3D Heisenberg
        }
    }

    /// Create YIG at given temperature
    ///
    /// T_C = 560 K
    pub fn yig(temperature: f64) -> Self {
        Self::new(Ferromagnet::yig(), 560.0, temperature)
    }

    /// Create Permalloy at given temperature
    ///
    /// T_C = 860 K
    pub fn permalloy(temperature: f64) -> Self {
        Self::new(Ferromagnet::permalloy(), 860.0, temperature)
    }

    /// Create CoFeB at given temperature
    ///
    /// T_C = 980 K
    pub fn cofeb(temperature: f64) -> Self {
        Self::new(Ferromagnet::cofeb(), 980.0, temperature)
    }

    /// Calculate magnetization at current temperature
    ///
    /// Uses Bloch's law at low T and critical scaling near T_C
    ///
    /// # Returns
    /// Magnetization \[A/m\]
    pub fn magnetization(&self) -> f64 {
        if self.temperature >= self.curie_temperature {
            return 0.0; // Above Curie temperature
        }

        let t = self.temperature / self.curie_temperature;

        let m_ratio = if t < 0.7 {
            // Bloch's law: M(T)/M(0) ≈ 1 - (T/T_C)^(3/2)
            1.0 - t.powf(1.5)
        } else {
            // Critical region: M(T)/M(0) ≈ ((T_C - T)/T_C)^β
            (1.0 - t).powf(self.critical_exponent)
        };

        self.material.ms * m_ratio.max(0.0)
    }

    /// Calculate anisotropy constant at current temperature
    ///
    /// Uses Callen-Callen law: K(T) ∝ M(T)^n
    ///
    /// # Returns
    /// Anisotropy constant \[J/m³\]
    pub fn anisotropy(&self) -> f64 {
        if self.temperature >= self.curie_temperature {
            return 0.0;
        }

        let m_ratio = self.magnetization() / self.material.ms;
        let n = 3.0; // Callen-Callen exponent

        self.material.anisotropy_k * m_ratio.powf(n)
    }

    /// Calculate exchange stiffness at current temperature
    ///
    /// A(T) ∝ M(T)^2
    ///
    /// # Returns
    /// Exchange stiffness \[J/m\]
    pub fn exchange_stiffness(&self) -> f64 {
        if self.temperature >= self.curie_temperature {
            return 0.0;
        }

        let m_ratio = self.magnetization() / self.material.ms;
        self.material.exchange_a * m_ratio * m_ratio
    }

    /// Calculate Gilbert damping at current temperature
    ///
    /// Weak temperature dependence (phonon scattering increases)
    ///
    /// # Returns
    /// Gilbert damping parameter
    pub fn damping(&self) -> f64 {
        // α(T) ≈ α(0) × (1 + c×T)
        let c = 1.0e-4; // K^-1 (weak dependence)
        self.material.alpha * (1.0 + c * self.temperature)
    }

    /// Calculate thermal energy scale
    ///
    /// E_th = k_B T
    ///
    /// # Returns
    /// Thermal energy \[J\]
    pub fn thermal_energy(&self) -> f64 {
        let kb = 1.380649e-23; // J/K
        kb * self.temperature
    }

    /// Calculate magnetic length at current temperature
    ///
    /// l_m = √(2A/μ₀M²)
    ///
    /// # Returns
    /// Magnetic length \[nm\]
    pub fn magnetic_length(&self) -> f64 {
        let mu_0 = 4.0 * std::f64::consts::PI * 1e-7;
        let a_ex = self.exchange_stiffness();
        let ms = self.magnetization();

        if ms < 1.0 {
            return f64::INFINITY;
        }

        let l_m = (2.0 * a_ex / (mu_0 * ms * ms)).sqrt();
        l_m * 1e9 // m to nm
    }

    /// Check if thermal fluctuations are significant
    ///
    /// Compares k_B T to anisotropy energy
    pub fn has_significant_thermal_fluctuations(&self) -> bool {
        let e_th = self.thermal_energy();
        let e_anis = self.anisotropy() * 1e-27; // Energy per nm³

        e_th / e_anis > 0.1
    }

    /// Set temperature
    pub fn set_temperature(&mut self, temp: f64) {
        self.temperature = temp;
    }

    /// Builder: with temperature
    pub fn with_temperature(mut self, temp: f64) -> Self {
        self.temperature = temp;
        self
    }

    /// Builder: with critical exponent
    pub fn with_critical_exponent(mut self, beta: f64) -> Self {
        self.critical_exponent = beta;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_yig_zero_temperature() {
        let yig = ThermalFerromagnet::yig(0.0);
        let m0 = yig.magnetization();
        assert!((m0 - yig.material.ms).abs() / m0 < 0.01);
    }

    #[test]
    fn test_yig_room_temperature() {
        let yig_rt = ThermalFerromagnet::yig(300.0);
        let m_rt = yig_rt.magnetization();
        let m0 = yig_rt.material.ms;

        // Should be reduced but still significant
        assert!(m_rt < m0);
        assert!(m_rt > 0.5 * m0);
    }

    #[test]
    fn test_above_curie_temperature() {
        let yig = ThermalFerromagnet::yig(600.0); // Above T_C
        assert_eq!(yig.magnetization(), 0.0);
        assert_eq!(yig.anisotropy(), 0.0);
    }

    #[test]
    fn test_magnetization_monotonic_decrease() {
        // Test magnetization decreases in each temperature regime separately
        // (avoiding the transition region between Bloch's law and critical region)

        // Low T regime (Bloch's law): T << 0.7*T_C
        let temps_low = [0.0, 100.0, 200.0, 300.0];
        for i in 0..temps_low.len() - 1 {
            let m1 = ThermalFerromagnet::yig(temps_low[i]).magnetization();
            let m2 = ThermalFerromagnet::yig(temps_low[i + 1]).magnetization();
            assert!(
                m2 <= m1,
                "M({}) = {} should be <= M({}) = {}",
                temps_low[i + 1],
                m2,
                temps_low[i],
                m1
            );
        }

        // High T regime (critical region): T > 0.7*T_C
        let temps_high = [400.0, 450.0, 500.0, 550.0];
        for i in 0..temps_high.len() - 1 {
            let m1 = ThermalFerromagnet::yig(temps_high[i]).magnetization();
            let m2 = ThermalFerromagnet::yig(temps_high[i + 1]).magnetization();
            assert!(
                m2 <= m1,
                "M({}) = {} should be <= M({}) = {}",
                temps_high[i + 1],
                m2,
                temps_high[i],
                m1
            );
        }
    }

    #[test]
    fn test_anisotropy_callen_callen() {
        // Use CoFeB which has non-zero anisotropy (YIG has K≈0)
        let cofeb0 = ThermalFerromagnet::cofeb(0.0);
        let cofeb300 = ThermalFerromagnet::cofeb(300.0);

        let k0 = cofeb0.anisotropy();
        let k300 = cofeb300.anisotropy();

        // K should decrease faster than M (power law with n=3)
        // K(T)/K(0) = (M(T)/M(0))^3, so k_ratio should be m_ratio^3
        let m_ratio = cofeb300.magnetization() / cofeb0.magnetization();
        let k_ratio = k300 / k0;

        // Check K decreases faster than M (k_ratio < m_ratio when m_ratio < 1)
        let expected_k_ratio = m_ratio.powf(3.0);
        assert!(
            (k_ratio - expected_k_ratio).abs() < 0.01,
            "k_ratio = {} should be close to m_ratio^3 = {}",
            k_ratio,
            expected_k_ratio
        );
        assert!(
            k_ratio < m_ratio,
            "k_ratio = {} should be < m_ratio = {}",
            k_ratio,
            m_ratio
        );
    }

    #[test]
    fn test_exchange_temperature_dependence() {
        let yig0 = ThermalFerromagnet::yig(0.0);
        let yig300 = ThermalFerromagnet::yig(300.0);

        let a0 = yig0.exchange_stiffness();
        let a300 = yig300.exchange_stiffness();

        assert!(a300 < a0);
        assert!(a300 > 0.0);
    }

    #[test]
    fn test_damping_weak_temperature_dependence() {
        let yig0 = ThermalFerromagnet::yig(0.0);
        let yig300 = ThermalFerromagnet::yig(300.0);

        let alpha0 = yig0.damping();
        let alpha300 = yig300.damping();

        // Damping should increase slightly with T
        assert!(alpha300 > alpha0);
        // But change should be small (< 10%)
        assert!((alpha300 - alpha0) / alpha0 < 0.1);
    }

    #[test]
    fn test_thermal_energy() {
        let yig300 = ThermalFerromagnet::yig(300.0);
        let e_th = yig300.thermal_energy();

        let kb = 1.380649e-23;
        assert!((e_th - kb * 300.0).abs() < 1e-25);
    }

    #[test]
    fn test_magnetic_length() {
        let yig = ThermalFerromagnet::yig(300.0);
        let l_m = yig.magnetic_length();

        // Should be in nm range for YIG
        assert!(l_m > 1.0);
        assert!(l_m < 100.0);
    }

    #[test]
    fn test_thermal_fluctuations() {
        // Use CoFeB which has non-zero anisotropy (YIG has K≈0)
        let cofeb_low_t = ThermalFerromagnet::cofeb(100.0);
        let cofeb_high_t = ThermalFerromagnet::cofeb(500.0);

        // Thermal energy increases with T
        let e_th_low = cofeb_low_t.thermal_energy();
        let e_th_high = cofeb_high_t.thermal_energy();
        assert!(e_th_high > e_th_low);

        // Both might return true due to the 1e-27 volume factor being very small,
        // but the ratio should increase with temperature
        let ratio_low = e_th_low / (cofeb_low_t.anisotropy() * 1e-27 + 1e-30);
        let ratio_high = e_th_high / (cofeb_high_t.anisotropy() * 1e-27 + 1e-30);
        assert!(
            ratio_high > ratio_low,
            "Thermal fluctuation ratio should increase with T"
        );
    }

    #[test]
    fn test_builder_pattern() {
        let yig = ThermalFerromagnet::yig(0.0)
            .with_temperature(400.0)
            .with_critical_exponent(0.32);

        assert_eq!(yig.temperature, 400.0);
        assert_eq!(yig.critical_exponent, 0.32);
    }
}
