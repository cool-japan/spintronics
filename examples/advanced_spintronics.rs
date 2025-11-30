//! Advanced Spintronics Showcase
//!
//! Demonstrates three cutting-edge areas of spintronics research:
//! 1. Antiferromagnetic THz Spintronics
//! 2. Stochastic Thermal Dynamics
//! 3. Cavity Magnonics - Hybrid Quantum Systems
//!
//! Based on Prof. Eiji Saitoh's research group

use spintronics::afm::Antiferromagnet;
use spintronics::cavity::HybridSystem;
use spintronics::stochastic::{StochasticLLG, ThermalField};
use spintronics::vector3::Vector3;

fn main() {
    println!("=== Advanced Spintronics Extensions ===\n");
    println!("Showcasing three frontier research areas:\n");

    // ========================================================================
    // Part 1: Antiferromagnetic THz Spintronics
    // ========================================================================
    println!("--- Part 1: Antiferromagnetic THz Dynamics ---\n");

    let mut afm = Antiferromagnet::nio();

    println!("Material: NiO (canonical antiferromagnet)");
    println!("  Exchange field: {:.2e} A/m", afm.h_exchange);
    println!("  Anisotropy field: {:.2e} A/m", afm.h_anisotropy);

    let f_res = afm.resonance_frequency();
    println!("  Resonance frequency: {:.2} GHz", f_res / 1e9);
    println!();

    // Initial state
    let neel_initial = afm.neel_vector();
    println!("Initial Neel vector (order parameter):");
    println!(
        "  n = [{:.3}, {:.3}, {:.3}]",
        neel_initial.x, neel_initial.y, neel_initial.z
    );

    // Apply THz pulse (transverse field)
    let h_thz_pulse = Vector3::new(1.0e6, 0.0, 0.0); // 1 T-equivalent
    println!();
    println!("Applying THz electromagnetic pulse...");
    println!("  H_pulse = {:.2e} A/m (x-direction)", h_thz_pulse.x);

    // Ultrafast evolution with femtosecond timesteps
    let dt_fs = 10.0e-15; // 10 femtoseconds
    let n_steps = 100;

    println!("  Evolving with Δt = {} fs", dt_fs * 1e15);

    for _ in 0..n_steps {
        afm.evolve(h_thz_pulse, dt_fs);
    }

    let neel_final = afm.neel_vector();
    println!();
    println!("After {} fs:", n_steps as f64 * dt_fs * 1e15);
    println!(
        "  n = [{:.3}, {:.3}, {:.3}]",
        neel_final.x, neel_final.y, neel_final.z
    );
    println!(
        "  Rotation angle: {:.2}°",
        (neel_initial.dot(&neel_final).acos() * 180.0 / std::f64::consts::PI)
    );
    println!();

    // ========================================================================
    // Part 2: Stochastic Thermal Dynamics
    // ========================================================================
    println!("--- Part 2: Finite Temperature Effects ---\n");

    // Compare zero temperature vs 300K dynamics
    let temp_cold = 0.0; // Zero temperature
    let temp_room = 300.0; // Room temperature

    let volume = 1.0e-24; // 1 nm³
    let ms = 1.0e6; // 1 MA/m
    let alpha = 0.01;

    let mut thermal_cold = ThermalField::new(temp_cold, volume, ms, alpha);
    let mut thermal_room = ThermalField::new(temp_room, volume, ms, alpha);

    println!("Thermal field generation:");
    println!("  Temperature: {} K vs {} K", temp_cold, temp_room);
    println!("  Cell volume: {:.1} nm³", volume * 1e27);
    println!();

    let dt = 1.0e-12;
    let h_cold = thermal_cold.generate(dt);
    let h_room = thermal_room.generate(dt);

    println!("Thermal field magnitude:");
    println!("  0 K: {:.3e} A/m (should be zero)", h_cold.magnitude());
    println!(
        "  300 K: {:.3e} A/m (thermal fluctuations)",
        h_room.magnitude()
    );
    println!();

    // Stochastic LLG simulation
    println!("Stochastic LLG at 300 K:");

    let mut sllg = StochasticLLG::new(Vector3::new(1.0, 0.0, 0.0), 300.0, volume, ms, alpha);

    let h_ext = Vector3::new(0.0, 0.0, 1.0e5); // 100 kA/m bias

    let m_initial = sllg.magnetization;

    // Evolve with thermal noise
    for _ in 0..1000 {
        sllg.evolve(h_ext, 1.0e-13);
    }

    let m_final = sllg.magnetization;

    println!(
        "  Initial: m = [{:.3}, {:.3}, {:.3}]",
        m_initial.x, m_initial.y, m_initial.z
    );
    println!(
        "  After 100 ps: m = [{:.3}, {:.3}, {:.3}]",
        m_final.x, m_final.y, m_final.z
    );
    println!(
        "  Thermal deviation: {:.3e}",
        (m_initial - m_final).magnitude()
    );
    println!();

    // Thermally activated switching potential
    println!("Thermally activated switching:");
    let k_b: f64 = 1.38e-23;
    let barrier_energy = 40.0 * k_b * 300.0; // 40 kB T barrier
    let f0 = thermal_room.attempt_frequency();
    let switching_rate = f0 * (-barrier_energy / (k_b * 300.0)).exp();

    println!(
        "  Energy barrier: {:.1} k_B T",
        barrier_energy / (k_b * 300.0)
    );
    println!("  Attempt frequency: {:.2} GHz", f0 / 1e9);
    println!("  Switching rate: {:.3e} Hz", switching_rate);
    println!(
        "  Retention time: {:.2} years",
        1.0 / (switching_rate * 365.25 * 24.0 * 3600.0)
    );
    println!();

    // ========================================================================
    // Part 3: Cavity Magnonics
    // ========================================================================
    println!("--- Part 3: Hybrid Magnon-Photon System ---\n");

    let mut hybrid = HybridSystem::yig_cavity();

    println!("System: YIG sphere in microwave cavity");
    println!(
        "  Cavity frequency: {:.2} GHz",
        hybrid.omega_cavity / (2.0 * std::f64::consts::PI) / 1e9
    );
    println!(
        "  Magnon (FMR) frequency: {:.2} GHz",
        hybrid.omega_magnon / (2.0 * std::f64::consts::PI) / 1e9
    );
    println!(
        "  Coupling strength: {:.1} MHz",
        hybrid.coupling_g / (2.0 * std::f64::consts::PI) / 1e6
    );
    println!();

    // Check coupling regime
    let cooperativity = hybrid.cooperativity();
    let is_strong = hybrid.is_strong_coupling();

    println!("Coupling regime:");
    println!("  Cooperativity C = {:.2}", cooperativity);
    println!(
        "  Status: {}",
        if is_strong {
            "STRONG COUPLING ✓"
        } else {
            "Weak coupling"
        }
    );

    if is_strong {
        let splitting = hybrid.rabi_splitting() / (2.0 * std::f64::consts::PI);
        println!("  Vacuum Rabi splitting: {:.1} MHz", splitting / 1e6);
    }
    println!();

    // Detuning scan
    println!("Detuning scan:");
    let detuning = hybrid.detuning();
    println!(
        "  Δω = ω_cavity - ω_magnon = {:.2} MHz",
        detuning / (2.0 * std::f64::consts::PI) / 1e6
    );
    println!();

    // Drive the system
    println!("Applying microwave drive at cavity frequency...");
    let h_drive = Vector3::new(1000.0, 0.0, 0.0); // Weak drive

    let cavity_amp_initial = hybrid.cavity_amplitude;

    for _ in 0..10000 {
        hybrid.evolve(h_drive, 1.0e-12); // 1 ps timestep
    }

    let cavity_amp_final = hybrid.cavity_amplitude;
    let m_final_cavity = hybrid.magnetization;

    println!("  After 10 ns:");
    println!(
        "  Cavity amplitude: {:.3e} → {:.3e}",
        cavity_amp_initial, cavity_amp_final
    );
    println!(
        "  Magnetization: [{:.3}, {:.3}, {:.3}]",
        m_final_cavity.x, m_final_cavity.y, m_final_cavity.z
    );
    println!();

    // ========================================================================
    // Summary
    // ========================================================================
    println!("=== Summary ===\n");
    println!("Three advanced frontiers demonstrated:");
    println!();
    println!("1. **AFM THz Spintronics**");
    println!(
        "   - Resonance: {:.1} GHz (orders faster than ferromagnets)",
        f_res / 1e9
    );
    println!("   - Ultrafast: femtosecond dynamics");
    println!("   - Application: Ultrafast magnetic memory, THz sources");
    println!();
    println!("2. **Stochastic Thermal Dynamics**");
    println!("   - Thermal noise: crucial for switching, retention");
    println!("   - Fluctuation-Dissipation Theorem enforced");
    println!("   - Application: Magnetic memory stability, spin caloritronics");
    println!();
    println!("3. **Cavity Magnonics**");
    println!("   - Strong coupling: C = {:.2}", cooperativity);
    println!("   - Hybridization: magnon-polariton quasiparticles");
    println!("   - Application: Quantum transducers, magnon lasers");
    println!();
    println!("All three represent active research areas in Prof. Saitoh's group!");
}
