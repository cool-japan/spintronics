//! Topological Insulator Surface State Transport
//!
//! This example demonstrates:
//! - Surface state band structure (Dirac cone)
//! - Spin-momentum locking
//! - Topological protection of transport
//! - Proximity-induced magnetism
//! - Quantum anomalous Hall effect
//!
//! ## Physics:
//! Topological insulators have insulating bulk but conducting surface states
//! with locked spin and momentum. When magnetized, they exhibit dissipationless
//! edge currents (quantum anomalous Hall effect).
//!
//! ## References:
//! - M. Z. Hasan & C. L. Kane, "Colloquium: Topological insulators",
//!   Rev. Mod. Phys. 82, 3045 (2010)
//! - C.-Z. Chang et al., "Experimental observation of the quantum anomalous
//!   Hall effect in a magnetic topological insulator", Science 340, 167 (2013)
//!
//! Run with: cargo run --example topological_insulator

use std::f64::consts::PI;

use spintronics::prelude::*;

fn main() {
    println!("=== Topological Insulator Surface States ===\n");

    // === 1. Material Selection ===
    println!("=== Material System ===");

    let bi2se3 = TopologicalInsulator::bi2se3();

    // Calculate effective parameters
    let fermi_energy_ev = 0.1; // eV, typical for n-doped Bi2Se3
    let spin_hall_angle = 0.5; // Estimated from spin Hall conductivity

    println!("Material: Bi₂Se₃");
    println!(
        "  Topological class: {}",
        match bi2se3.ti_class {
            TopologicalClass::ThreeDimensional => "3D Strong TI",
            TopologicalClass::TwoDimensional => "2D QSHI",
        }
    );
    println!("  Bulk bandgap: {:.2} eV", bi2se3.bulk_gap);
    println!(
        "  Surface Dirac velocity: {:.2e} m/s",
        bi2se3.fermi_velocity
    );
    println!("  Fermi energy: {:.3} eV (n-doped)", fermi_energy_ev);
    println!("  Spin Hall angle: {:.2}\n", spin_hall_angle);

    // === 2. Surface State Dispersion ===
    println!("=== Dirac Cone Dispersion ===");

    let hbar = HBAR;
    let e = E_CHARGE;

    println!("E(k) = ℏv_F |k| (linear Dirac dispersion)\n");

    println!("Wavevector (nm⁻¹) | Energy (meV) | Group velocity (km/s)");
    println!("--------------------------------------------------------------");

    for k_inv_nm in [0.0, 0.01, 0.02, 0.05, 0.1, 0.2] {
        let k = k_inv_nm * 1e9; // Convert to m⁻¹
        let energy = hbar * bi2se3.fermi_velocity * k / e * 1e3; // meV
        let v_group = bi2se3.fermi_velocity * 1e-3; // km/s

        println!(
            "  {:.3}              | {:.2}         | {:.0}",
            k_inv_nm, energy, v_group
        );
    }

    // === 3. Spin-Momentum Locking ===
    println!("\n=== Spin-Momentum Locking ===");

    println!("Surface states: spin ⊥ momentum");
    println!("  ψ_k ∝ (|↑⟩ + iexp(iθ_k)|↓⟩)");
    println!("  where θ_k = arctan(k_y / k_x)\n");

    // Example: electron moving in +x direction
    let k_x: f64 = 0.1e9; // m⁻¹
    let k_y: f64 = 0.0;
    let theta_k = k_y.atan2(k_x);

    println!(
        "Electron momentum: k = ({:.2}, {:.2}) nm⁻¹",
        k_x * 1e-9,
        k_y * 1e-9
    );
    println!(
        "  → Spin direction: θ = {:.1}° (perpendicular)",
        theta_k * 180.0 / PI + 90.0
    );
    println!("  → Spin polarization: {:.0}% (maximal)\n", 100.0);

    // === 4. Transport Properties ===
    println!("=== Quantum Transport ===");

    // Conductivity from surface states
    let n2d = fermi_wavevector(&bi2se3).powi(2) / (4.0 * PI); // Surface carrier density
    let mu = 2000.0; // cm²/(V·s) typical mobility

    println!("  2D carrier density: {:.2e} cm⁻²", n2d * 1e-4);
    println!("  Mobility: {:.0} cm²/(V·s)", mu);

    let sigma_2d = e * n2d * mu * 1e-4; // S/m (per surface)
    let sheet_resistance = 1.0 / sigma_2d * 1e-3; // kΩ/square

    println!("  Sheet conductance: {:.2e} S", sigma_2d);
    println!("  Sheet resistance: {:.2} kΩ/□\n", sheet_resistance);

    // === 5. Topological Protection ===
    println!("=== Topological Protection ===");

    println!("Backscattering suppression:");
    println!("  • Non-magnetic impurities: FORBIDDEN by time-reversal");
    println!("  • Elastic scattering: Cannot flip spin → cannot backscatter");
    println!("  • Result: High mobility even with disorder\n");

    let tau = mu * 1e-4 * 9.1e-31 / e; // Scattering time
    let mfp = bi2se3.fermi_velocity * tau; // Mean free path

    println!("  Scattering time: {:.2} ps", tau * 1e12);
    println!("  Mean free path: {:.0} nm (long!)\n", mfp * 1e9);

    // === 6. Magnetic Proximity Effect ===
    println!("=== Magnetic Proximity Effect ===");

    println!("Depositing ferromagnet (EuS) on Bi₂Se₃:");
    println!("  • Breaks time-reversal symmetry");
    println!("  • Opens gap at Dirac point");
    println!("  • Enables quantum anomalous Hall effect\n");

    let m_exchange = 0.05; // eV exchange gap
    let gap_magnetic = 2.0 * m_exchange;

    println!("  Exchange coupling: {:.0} meV", m_exchange * 1e3);
    println!("  Magnetic gap: {:.0} meV", gap_magnetic * 1e3);
    println!(
        "  Gap/k_B T (300K): {:.1} (stable)\n",
        gap_magnetic / (KB * 300.0 / e)
    );

    // === 7. Quantum Anomalous Hall Effect ===
    println!("=== Quantum Anomalous Hall Effect ===");

    println!("When magnetized, TI exhibits:");
    println!("  • Chiral edge currents (one direction only)");
    println!("  • Quantized Hall conductance: σ_xy = e²/h");
    println!("  • Zero longitudinal resistance");
    println!("  • No external magnetic field needed!\n");

    let g0 = e * e / HBAR; // Conductance quantum
    let sigma_xy_qahe = g0 * 1.0; // Chern number = 1

    println!("  Quantized Hall conductance: {:.5e} S", sigma_xy_qahe);
    println!("  (= e²/h = {:.5e} S)", g0);
    println!("  Hall resistance: {:.2} kΩ\n", 1.0 / sigma_xy_qahe * 1e-3);

    // Edge state velocity
    let v_edge = gap_magnetic * e / hbar; // From E = ħv k

    println!("  Edge state velocity: {:.2e} m/s", v_edge);
    println!("  Edge current (per channel): {:.2} µA\n", e * v_edge * 1e6);

    // === 8. Spin-Charge Conversion ===
    println!("=== Edelstein Effect ===");

    let j_charge = 1e10; // A/m²
    let lambda_edelstein = bi2se3.edelstein_length * 1e-9; // Convert from nm to m

    println!("  Charge current: {:.0} MA/cm²", j_charge * 1e-10);
    println!("  Edelstein length: {:.2} nm", bi2se3.edelstein_length);

    let spin_density = lambda_edelstein * j_charge;

    println!("  Induced spin density: {:.2e} ħ/m²", spin_density / hbar);
    println!(
        "  Conversion efficiency: {:.1}%\n",
        bi2se3.edelstein_length / 10.0
    );

    // === 9. Comparison with Other Systems ===
    println!("=== Comparison with Conventional Systems ===");

    println!("\nBi₂Se₃ Surface States:");
    println!("  • Spin-momentum locking: 100%");
    println!("  • Spin Hall angle: {:.2}", spin_hall_angle);
    println!("  • Topologically protected");

    println!("\nConventional Rashba (Au surface):");
    println!("  • Spin-momentum locking: ~50%");
    println!("  • Spin Hall angle: ~0.1");
    println!("  • Not protected");

    println!("\nHeavy metals (Pt):");
    println!("  • No momentum locking");
    println!("  • Spin Hall angle: ~0.1");
    println!("  • Bulk effect\n");

    // === 10. Device Applications ===
    println!("=== Device Applications ===");

    println!("\n1. Dissipationless Interconnects:");
    println!("  • QAHE edge currents: zero resistance");
    println!("  • No Joule heating");
    println!("  • Lower power consumption");

    println!("\n2. Spintronic Memory:");
    println!("  • Efficient SOT switching");
    println!("  • θ_SH = {:.2} (high)", spin_hall_angle);
    println!("  • Lower critical current");

    let j_c_ti = 1e10; // A/m² for TI-based SOT
    let j_c_pt = 1e11; // A/m² for Pt-based SOT

    println!("  Critical current (TI): {:.0} MA/cm²", j_c_ti * 1e-10);
    println!("  Critical current (Pt): {:.0} MA/cm²", j_c_pt * 1e-10);
    println!("  → 10× reduction!");

    println!("\n3. Quantum Computing:");
    println!("  • Majorana fermions at TI/SC interface");
    println!("  • Topological qubits");
    println!("  • Protected quantum information");

    // === 11. Experimental Signatures ===
    println!("\n=== Experimental Signatures ===");

    println!("\nARPES (angle-resolved photoemission):");
    println!("  • Direct measurement of Dirac cone");
    println!("  • Confirms linear dispersion");
    println!("  • Shows spin texture");

    println!("\nQuantum oscillations:");
    println!("  • Shubnikov-de Haas at high B");
    println!("  • 2D Fermi surface");
    println!("  • Berry phase = π (topological)");

    println!("\nHall measurements:");
    println!("  • σ_xy = e²/h at low T");
    println!("  • Quantization plateau");
    println!("  • Hysteresis with magnetization");

    // === Summary ===
    println!("\n=== Summary ===");

    println!("\nKey Physics:");
    println!("  ✓ Bulk insulator, surface conductor");
    println!("  ✓ Dirac dispersion: E = ℏv_F |k|");
    println!("  ✓ Spin-momentum locking (100%)");
    println!("  ✓ Topological protection of transport");
    println!("  ✓ QAHE when magnetized");

    println!("\nDevice Advantages:");
    println!("  ✓ Ultra-efficient spin-charge conversion");
    println!("  ✓ Low-power operation");
    println!("  ✓ High spin Hall angle ({:.1})", spin_hall_angle);
    println!("  ✓ Dissipationless edge transport");

    println!("\nChallenges:");
    println!("  • Bulk conductivity from doping");
    println!("  • QAHE requires low temperature (<1K currently)");
    println!("  • Interface quality critical");

    println!("\nFuture Directions:");
    println!("  → Room-temperature QAHE");
    println!("  → Integration with magnetic materials");
    println!("  → Topological quantum computation");
    println!("  → Magneto-electric effects");
}

// Helper functions

fn fermi_wavevector(ti: &TopologicalInsulator) -> f64 {
    // k_F = E_F / (ℏv_F)
    let hbar = HBAR;
    let e = E_CHARGE;
    let fermi_energy_ev = 0.1; // eV, typical for n-doped Bi2Se3

    fermi_energy_ev * e / (hbar * ti.fermi_velocity)
}
