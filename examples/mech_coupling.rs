//! Coupled Spin-Mechanical Dynamics
//!
//! Demonstrates angular momentum exchange between spin and mechanical degrees of freedom.
//! Shows both Barnett effect (rotation → magnetization) and
//! Einstein-de Haas effect (magnetization → rotation).

use spintronics::mech::{BarnettMagnetization, Cantilever, EinsteinDeHaas, SpinMechanicalCoupling};
use spintronics::vector3::Vector3;

fn main() {
    println!("=== Spin-Mechanical Coupling ===\n");

    // Part 1: Barnett Effect
    println!("--- Part 1: Barnett Effect (Rotation → Magnetization) ---\n");

    let barnett = BarnettMagnetization::permalloy();

    println!("Material: Permalloy (Ni₈₀Fe₂₀)");
    println!("  Saturation magnetization: {:.2e} A/m", barnett.ms);
    println!();

    // Moderate rotation
    let omega_moderate = Vector3::new(0.0, 0.0, 1000.0); // 1000 rad/s
    let m_barnett = barnett.magnetization_from_rotation(omega_moderate);

    println!("Rotation: Ω = 1000 rad/s (z-axis)");
    println!(
        "  Frequency: {:.1} Hz",
        1000.0 / (2.0 * std::f64::consts::PI)
    );
    println!();

    println!("Induced magnetization:");
    println!("  M_B = {:.3e} A/m", m_barnett.magnitude());
    println!(
        "  Fractional: M_B/M_s = {:.3e}",
        barnett.fractional_magnetization(omega_moderate)
    );
    println!();

    // Estimate rotation needed for saturation
    let omega_sat = barnett.saturation_rotation_rate();
    println!("Rotation rate for saturation:");
    println!("  Ω_sat = {:.3e} rad/s", omega_sat);
    println!("  (unrealistically high for macroscopic samples)");
    println!();

    // Part 2: Einstein-de Haas Effect
    println!("--- Part 2: Einstein-de Haas Effect (Magnetization → Rotation) ---\n");

    let moment_of_inertia = 1.0e-15; // Small sample: 1 fg·m²
    let volume = 1.0e-18; // 1 μm³
    let density = 8900.0; // Ni density

    let edh = EinsteinDeHaas::new(moment_of_inertia, volume, density);

    println!("Sample: Permalloy nanowire");
    println!("  Volume: 1 μm³");
    println!("  Moment of inertia: {:.2e} kg·m²", moment_of_inertia);
    println!();

    // Magnetization flip
    let ms_permalloy = 8.0e5;
    let rotation_angle = edh.rotation_angle_magnetization_flip(ms_permalloy);

    println!("Magnetization flip: M → -M");
    println!("  ΔM = 2 M_s = {:.2e} A/m", 2.0 * ms_permalloy);
    println!(
        "  Rotation angle: {:.3e} rad = {:.3e} degrees",
        rotation_angle,
        rotation_angle * 180.0 / std::f64::consts::PI
    );
    println!();

    // Instantaneous angular velocity
    let delta_m = Vector3::new(0.0, 0.0, 2.0 * ms_permalloy);
    let omega_edh = edh.angular_velocity_from_magnetization_change(delta_m);

    println!("Instantaneous angular velocity:");
    println!("  Ω = {:.3e} rad/s", omega_edh.magnitude());
    println!();

    // Part 3: Coupled Dynamics Simulation
    println!("--- Part 3: Coupled Dynamics Simulation ---\n");

    let mut coupling = SpinMechanicalCoupling::new(
        Vector3::new(1.0, 0.0, 0.0), // Initial magnetization along x
        moment_of_inertia,
        0.01, // Gilbert damping
        ms_permalloy,
        volume,
    );

    println!("Initial state:");
    println!(
        "  m = [{:.2}, {:.2}, {:.2}]",
        coupling.magnetization.x, coupling.magnetization.y, coupling.magnetization.z
    );
    println!(
        "  Ω = [{:.2e}, {:.2e}, {:.2e}] rad/s",
        coupling.omega.x, coupling.omega.y, coupling.omega.z
    );
    println!();

    // Apply external field along z
    let h_ext = Vector3::new(0.0, 0.0, 1.0e5); // 100 kA/m
    let dt = 1.0e-13; // 0.1 ps
    let n_steps = 1000;

    println!("Applying external field: H_ext = 100 kA/m (z-direction)");
    println!("Time step: {} ps", dt * 1e12);
    println!("Total time: {} ns", dt * n_steps as f64 * 1e9);
    println!();

    println!("Evolving coupled system...");

    for i in 0..n_steps {
        coupling.evolve(h_ext, dt);

        if i % 200 == 0 {
            println!(
                "  t = {:.2} ns: m_z = {:.4}, Ω_z = {:.3e} rad/s",
                i as f64 * dt * 1e9,
                coupling.magnetization.z,
                coupling.omega.z
            );
        }
    }

    println!();
    println!("Final state:");
    println!(
        "  m = [{:.4}, {:.4}, {:.4}]",
        coupling.magnetization.x, coupling.magnetization.y, coupling.magnetization.z
    );
    println!(
        "  |m| = {:.6} (should be 1.0)",
        coupling.magnetization.magnitude()
    );
    println!(
        "  Ω = [{:.3e}, {:.3e}, {:.3e}] rad/s",
        coupling.omega.x, coupling.omega.y, coupling.omega.z
    );
    println!();

    // Part 4: Cantilever Detection
    println!("--- Part 4: Nanomechanical Detection ---\n");

    let mut cantilever = Cantilever::afm_cantilever();

    println!("AFM Cantilever:");
    println!("  Resonance frequency: {:.1} kHz", cantilever.f0 / 1e3);
    println!("  Q factor: {:.0}", cantilever.q_factor);
    println!("  Mass: {:.1} pg", cantilever.mass * 1e15);
    println!();

    // Apply oscillating force from magnetic torque
    let force_amplitude = 1.0e-12; // 1 pN
    let omega_drive = 2.0 * std::f64::consts::PI * cantilever.f0;
    let dt_cant = 1.0e-8; // 10 ns
    let n_cycles = 100;

    println!(
        "Driving force: {:.1} pN at resonance",
        force_amplitude * 1e12
    );
    println!("Simulating {} oscillation cycles...", n_cycles);

    for i in 0..(n_cycles * 20) {
        let t = i as f64 * dt_cant;
        let force = force_amplitude * (omega_drive * t).sin();
        cantilever.evolve(force, dt_cant);
    }

    println!();
    println!("Final cantilever state:");
    println!("  Displacement: {:.3} nm", cantilever.displacement * 1e9);
    println!("  Velocity: {:.3e} m/s", cantilever.velocity);
    println!("  Total energy: {:.3e} J", cantilever.total_energy());
    println!();

    println!("=== Summary ===");
    println!("1. Barnett effect: Rotation induces magnetization (very small for realistic rates)");
    println!(
        "2. Einstein-de Haas: Magnetization change causes rotation (angular momentum conservation)"
    );
    println!("3. Coupled dynamics: Magnetization and rotation exchange angular momentum");
    println!("4. Cantilever: Sensitive detector for mechanical motion induced by spin dynamics");
}
