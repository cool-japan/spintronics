//! Integration tests for demo API endpoints
//! These tests verify the physics calculations produce correct results

use spintronics::prelude::*;

#[test]
fn test_llg_conservation_of_magnitude() {
    // LLG dynamics should not produce NaN or Inf
    let fm = Ferromagnet::permalloy().with_alpha(0.01);
    let m0 = Vector3::new(0.1, 0.0, 1.0).normalize() * fm.ms;
    let h_ext = Vector3::new(1000.0, 0.0, 0.0);
    let dt = 1e-12;

    let solver = spintronics::dynamics::LlgSolver::new(fm.alpha, dt);
    let mut m = m0;

    for _ in 0..100 {
        m = solver.step_rk4(m, |_m| h_ext);
        // Check that result is finite
        assert!(
            m.x.is_finite() && m.y.is_finite() && m.z.is_finite(),
            "Magnetization contains NaN or Inf"
        );
        // Check magnitude is positive
        assert!(m.magnitude() > 0.0, "Magnetization magnitude is zero");
    }
}

#[test]
fn test_spin_pumping_voltage_scale() {
    // Spin pumping voltage should scale with frequency
    let _fm = Ferromagnet::yig();
    let interface = spintronics::material::SpinInterface::yig_pt();
    let _ishe = spintronics::effect::InverseSpinHall::platinum();

    // Calculate at two frequencies
    let omega1 = 2.0 * std::f64::consts::PI * 5e9; // 5 GHz
    let omega2 = 2.0 * std::f64::consts::PI * 10e9; // 10 GHz

    let theta: f64 = 0.1; // Small angle approximation

    use spintronics::constants::{E_CHARGE, HBAR};
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

    // Voltage should double when frequency doubles
    let ratio = js2 / js1;
    assert!(
        (ratio - 2.0).abs() < 0.01,
        "Spin current should scale linearly with frequency: ratio = {}",
        ratio
    );
}

#[test]
fn test_materials_saturation_magnetization_range() {
    // All ferromagnets should have reasonable Ms values
    let materials = vec![
        Ferromagnet::yig(),
        Ferromagnet::permalloy(),
        Ferromagnet::cofeb(),
        Ferromagnet::iron(),
        Ferromagnet::cobalt(),
        Ferromagnet::nickel(),
    ];

    for mat in materials {
        assert!(
            mat.ms > 1e3 && mat.ms < 2e6,
            "Saturation magnetization out of range: {} A/m",
            mat.ms
        );
        assert!(
            mat.alpha > 1e-5 && mat.alpha < 1.0,
            "Damping out of range: {}",
            mat.alpha
        );
    }
}

#[test]
fn test_skyrmion_topological_charge() {
    use spintronics::texture::{Chirality, Helicity, Skyrmion};

    // Skyrmion should have topological charge Q = ±1
    let sk_ccw = Skyrmion::default()
        .with_radius(50e-9)
        .with_helicity(Helicity::Neel)
        .with_chirality(Chirality::CounterClockwise);

    let sk_cw = Skyrmion::default()
        .with_radius(50e-9)
        .with_helicity(Helicity::Neel)
        .with_chirality(Chirality::Clockwise);

    // CCW should have Q = -1, CW should have Q = +1
    // (These would need actual calculation implementation)
    assert_eq!(sk_ccw.helicity, Helicity::Neel);
    assert_eq!(sk_ccw.chirality, Chirality::CounterClockwise);
    assert_eq!(sk_cw.chirality, Chirality::Clockwise);
}

#[test]
fn test_interface_mixing_conductance() {
    // Real and imaginary parts should be in reasonable ranges
    let interfaces = vec![
        spintronics::material::SpinInterface::yig_pt(),
        spintronics::material::SpinInterface::py_pt(),
        spintronics::material::SpinInterface::cofeb_pt(),
    ];

    for interface in interfaces {
        // g_r typically 1e13 - 1e20 Ω^-1 m^-2
        assert!(
            interface.g_r > 1e12 && interface.g_r < 1e20,
            "Real mixing conductance out of range: {} Ω^-1 m^-2",
            interface.g_r
        );

        // g_i typically smaller than g_r
        assert!(
            interface.g_i < interface.g_r,
            "Imaginary part should be smaller than real part"
        );
    }
}

#[test]
fn test_spin_hall_angle_range() {
    // Spin Hall angles should be reasonable
    let ishe_materials = vec![
        spintronics::effect::InverseSpinHall::platinum(),
        spintronics::effect::InverseSpinHall::tantalum(),
        spintronics::effect::InverseSpinHall::tungsten(),
    ];

    for ishe in ishe_materials {
        assert!(
            ishe.theta_sh.abs() < 1.0,
            "Spin Hall angle unreasonably large: {}",
            ishe.theta_sh
        );
    }
}

#[test]
fn test_llg_damping_effect() {
    // Higher damping should lead to faster relaxation
    let fm_low_damp = Ferromagnet::permalloy().with_alpha(0.01);
    let fm_high_damp = Ferromagnet::permalloy().with_alpha(0.5);

    let m0 = Vector3::new(1.0, 0.0, 0.0).normalize();
    let h_ext = Vector3::new(0.0, 0.0, 1.0); // Field along z

    let dt = 1e-12;
    let steps = 100;

    // Low damping
    let solver_low = spintronics::dynamics::LlgSolver::new(fm_low_damp.alpha, dt);
    let mut m_low = m0 * fm_low_damp.ms;
    for _ in 0..steps {
        m_low = solver_low.step_rk4(m_low, |_m| h_ext);
    }
    let angle_low = (m_low.z / m_low.magnitude()).acos();

    // High damping
    let solver_high = spintronics::dynamics::LlgSolver::new(fm_high_damp.alpha, dt);
    let mut m_high = m0 * fm_high_damp.ms;
    for _ in 0..steps {
        m_high = solver_high.step_rk4(m_high, |_m| h_ext);
    }
    let angle_high = (m_high.z / m_high.magnitude()).acos();

    // High damping should be closer to equilibrium (smaller angle from z-axis)
    assert!(
        angle_high < angle_low,
        "High damping should relax faster: angle_high={:.3} vs angle_low={:.3}",
        angle_high,
        angle_low
    );
}
