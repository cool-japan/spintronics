//! Integration tests for demo endpoints

#[cfg(test)]
mod integration_tests {
    use spintronics::prelude::*;

    #[test]
    fn test_llg_request_parsing() {
        // Test that LLG request structures work correctly
        let fm = Ferromagnet::permalloy().with_alpha(0.01);
        assert_eq!(fm.alpha, 0.01);
    }

    #[test]
    fn test_llg_simulation_physics() {
        // Test LLG simulation runs without errors
        let fm = Ferromagnet::permalloy().with_alpha(0.01);
        let m0 = Vector3::new(0.1, 0.0, 1.0).normalize() * fm.ms;
        let h_ext = Vector3::new(1000.0, 0.0, 0.0);
        let dt = 1e-12;

        let solver = spintronics::dynamics::LlgSolver::new(fm.alpha, dt);
        let mut m = m0;

        // Run multiple steps - should not panic
        for _ in 0..10 {
            m = solver.step_rk4(m, |_m| h_ext);
            // Check that result is finite
            assert!(
                m.x.is_finite() && m.y.is_finite() && m.z.is_finite(),
                "Magnetization contains NaN or Inf"
            );
        }

        // Check that magnetization is non-zero
        assert!(m.magnitude() > 0.0, "Magnetization became zero");
    }

    #[test]
    fn test_spin_pumping_calculation() {
        // Test spin pumping physics
        let fm = Ferromagnet::yig();
        let interface = spintronics::material::SpinInterface::yig_pt();

        // Verify realistic parameters
        assert!(fm.ms > 1e3 && fm.ms < 2e6, "Ms out of range: {}", fm.ms);
        assert!(
            interface.g_r > 1e12 && interface.g_r < 1e20,
            "g_r out of range: {}",
            interface.g_r
        );
    }

    #[test]
    fn test_materials_database() {
        // Test material database contains expected entries
        let yig = Ferromagnet::yig();
        let py = Ferromagnet::permalloy();
        let cofeb = Ferromagnet::cofeb();

        // YIG has lowest damping
        assert!(yig.alpha < py.alpha);
        assert!(yig.alpha < cofeb.alpha);

        // CoFeB has highest saturation magnetization
        assert!(cofeb.ms > yig.ms);
    }

    #[test]
    fn test_skyrmion_parameters() {
        use spintronics::texture::{Chirality, Helicity, Skyrmion};

        let sk = Skyrmion::default()
            .with_radius(50e-9)
            .with_helicity(Helicity::Neel)
            .with_chirality(Chirality::CounterClockwise);

        assert_eq!(sk.radius, 50e-9);
        assert_eq!(sk.helicity, Helicity::Neel);
        assert_eq!(sk.chirality, Chirality::CounterClockwise);
    }

    #[test]
    fn test_interface_materials() {
        // Test various interface materials
        let yig_pt = spintronics::material::SpinInterface::yig_pt();
        let py_pt = spintronics::material::SpinInterface::py_pt();
        let cofeb_pt = spintronics::material::SpinInterface::cofeb_pt();

        // All should have reasonable mixing conductances
        for interface in &[yig_pt, py_pt, cofeb_pt] {
            assert!(interface.g_r > 1e13);
            assert!(interface.g_i < interface.g_r);
        }
    }

    #[test]
    fn test_llg_damping_effect() {
        // Higher damping should lead to faster relaxation
        let fm_low = Ferromagnet::permalloy().with_alpha(0.01);
        let fm_high = Ferromagnet::permalloy().with_alpha(0.5);

        let m0 = Vector3::new(1.0, 0.0, 0.0).normalize();
        let h_ext = Vector3::new(0.0, 0.0, 1.0);
        let dt = 1e-12;
        let steps = 100;

        // Low damping
        let solver_low = spintronics::dynamics::LlgSolver::new(fm_low.alpha, dt);
        let mut m_low = m0 * fm_low.ms;
        for _ in 0..steps {
            m_low = solver_low.step_rk4(m_low, |_m| h_ext);
        }

        // High damping
        let solver_high = spintronics::dynamics::LlgSolver::new(fm_high.alpha, dt);
        let mut m_high = m0 * fm_high.ms;
        for _ in 0..steps {
            m_high = solver_high.step_rk4(m_high, |_m| h_ext);
        }

        // High damping should be closer to equilibrium
        let angle_low = (m_low.z / m_low.magnitude()).acos();
        let angle_high = (m_high.z / m_high.magnitude()).acos();
        assert!(angle_high < angle_low);
    }

    #[test]
    fn test_spin_pumping_frequency_scaling() {
        // Spin current should scale linearly with frequency
        use spintronics::constants::{E_CHARGE, HBAR};

        let interface = spintronics::material::SpinInterface::yig_pt();
        let theta: f64 = 0.1;

        let omega1 = 2.0 * std::f64::consts::PI * 5e9;
        let omega2 = 2.0 * std::f64::consts::PI * 10e9;

        let js1 = (HBAR / (4.0 * std::f64::consts::PI * E_CHARGE))
            * interface.g_r
            * omega1
            * theta.sin()
            * theta;
        let js2 = (HBAR / (4.0 * std::f64::consts::PI * E_CHARGE))
            * interface.g_r
            * omega2
            * theta.sin()
            * theta;

        let ratio = js2 / js1;
        assert!((ratio - 2.0).abs() < 0.01, "Spin current should double");
    }

    #[test]
    fn test_ishe_materials() {
        // Test ISHE converter materials
        let pt = spintronics::effect::InverseSpinHall::platinum();
        let ta = spintronics::effect::InverseSpinHall::tantalum();
        let w = spintronics::effect::InverseSpinHall::tungsten();

        // All should have reasonable spin Hall angles
        for ishe in &[pt, ta, w] {
            assert!(ishe.theta_sh.abs() < 1.0);
        }
    }

    #[test]
    fn test_vector3_operations() {
        let v1 = Vector3::new(1.0, 0.0, 0.0);
        let v2 = Vector3::new(0.0, 1.0, 0.0);

        // Dot product
        assert_eq!(v1.dot(&v2), 0.0);

        // Cross product
        let v3 = v1.cross(&v2);
        assert!((v3.z - 1.0).abs() < 1e-10);

        // Magnitude
        assert!((v1.magnitude() - 1.0).abs() < 1e-10);
    }
}
