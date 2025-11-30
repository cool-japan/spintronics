//! 2D Magnetic Materials for Spintronics
//!
//! This example demonstrates:
//! - Magnetism in 2D materials (CrI₃, Fe₃GeTe₂)
//! - Layer-dependent magnetic properties
//! - Van der Waals heterostructures
//! - Gate-tunable magnetism
//! - Spin-valley coupling
//!
//! ## Physics:
//! 2D magnetic materials enable atomically thin spintronic devices with
//! unprecedented tunability via electric fields, strain, and stacking.
//! They combine magnetism with strong spin-orbit coupling and valley physics.
//!
//! ## References:
//! - B. Huang et al., "Layer-dependent ferromagnetism in a van der Waals
//!   crystal down to the monolayer limit", Nature 546, 270 (2017)
//! - C. Gong et al., "Discovery of intrinsic ferromagnetism in two-dimensional
//!   van der Waals crystals", Nature 546, 265 (2017)
//!
//! Run with: cargo run --example 2d_material_spintronics

use spintronics::prelude::*;
use std::f64::consts::PI;

fn main() {
    println!("=== 2D Magnetic Materials for Spintronics ===\n");

    // === 1. Material Comparison ===
    println!("=== 2D Magnetic Materials ===");

    let cri3 = Magnetic2D::cri3(1); // Monolayer
    let fgt = Magnetic2D::fe3gete2(5); // Few-layer

    println!("CrI₃ (Chromium triiodide):");
    println!(
        "  Magnetic ordering: {}",
        match cri3.ordering {
            MagneticOrdering::Ferromagnetic => "Ferromagnetic",
            _ => "Other",
        }
    );
    println!(
        "  Curie temperature: {:.0} K (monolayer)",
        cri3.critical_temperature
    );
    println!("  Magnetization: {:.2} µB per atom", cri3.magnetic_moment);
    let cri3_out_of_plane = cri3.easy_axis.z.abs() > 0.5;
    println!(
        "  Easy axis: {}",
        if cri3_out_of_plane {
            "Out-of-plane"
        } else {
            "In-plane"
        }
    );
    println!("  Electronic structure: Insulator");

    println!("\nFe₃GeTe₂:");
    println!(
        "  Magnetic ordering: {}",
        match fgt.ordering {
            MagneticOrdering::Ferromagnetic => "Ferromagnetic",
            _ => "Other",
        }
    );
    println!(
        "  Curie temperature: {:.0} K (bulk)",
        fgt.critical_temperature
    );
    println!("  Magnetization: {:.2} µB per atom", fgt.magnetic_moment);
    let fgt_out_of_plane = fgt.easy_axis.z.abs() > 0.5;
    println!(
        "  Easy axis: {}",
        if fgt_out_of_plane {
            "Out-of-plane"
        } else {
            "In-plane"
        }
    );
    println!("  Electronic structure: Metallic");

    println!("\n→ CrI₃: Insulator, good for tunneling");
    println!("→ Fe₃GeTe₂: Metallic, good for spin injection\n");

    // === 2. Layer Dependence ===
    println!("=== Layer-Dependent Magnetism ===");

    println!("CrI₃ Curie temperature vs. layer number:");
    println!("  Monolayer:  {:.0} K", cri3.critical_temperature);
    println!(
        "  Bilayer:    ~{:.0} K (AFM interlayer)",
        cri3.critical_temperature * 1.2
    );
    println!(
        "  Trilayer:   ~{:.0} K (FM, highest)",
        cri3.critical_temperature * 1.5
    );
    println!("  Bulk:       {:.0} K\n", 61.0);

    println!("Interlayer coupling:");
    println!("  • Bilayer: Antiferromagnetic");
    println!("  • Can be switched by small B field (~0.65T)");
    println!("  • Gate voltage modulation");

    // === 3. Atomically Thin MTJ ===
    println!("\n=== Atomic-Scale Magnetic Tunnel Junction ===");

    println!("Structure: CrI₃ (3 layers) / hBN (2 nm) / CrI₃ (3 layers)");

    let d_barrier = 2e-9; // 2 nm hBN
    let _area = (100e-9) * (100e-9); // 100×100 nm²

    println!("  Barrier: hBN {:.0} nm", d_barrier * 1e9);
    println!("  Junction area: {:.0} × {:.0} nm²\n", 100.0, 100.0);

    // TMR calculation
    let p_cri3 = 0.9; // High spin polarization
    let tmr_max = 2.0 * p_cri3 * p_cri3 / (1.0 - p_cri3 * p_cri3);

    println!("Tunneling Magnetoresistance:");
    println!("  Spin polarization: {:.0}%", p_cri3 * 100.0);
    println!("  TMR (theoretical): {:.0}%", tmr_max * 100.0);

    // Resistance calculation
    let r_parallel = 10e3; // 10 kΩ typical
    let r_antiparallel = r_parallel * (1.0 + tmr_max);

    println!("  R_parallel: {:.0} kΩ", r_parallel * 1e-3);
    println!("  R_antiparallel: {:.0} kΩ", r_antiparallel * 1e-3);
    println!("  MR ratio: {:.0}%\n", tmr_max * 100.0);

    // === 4. Gate Tunability ===
    println!("=== Electric Field Control ===");

    println!("Dual-gate device:");
    println!("  Top gate: Control carrier density");
    println!("  Bottom gate: Tune magnetic anisotropy\n");

    let v_gate_array = [-2.0, -1.0, 0.0, 1.0, 2.0];
    println!("Gate voltage (V) | ΔT_C (K) | TMR change (%)");
    println!("------------------------------------------------");

    for v_gate in v_gate_array.iter() {
        let delta_tc = v_gate * 5.0; // ~5K per volt
        let tmr_change = v_gate * 10.0; // ~10% per volt

        println!(
            "  {:+.1}             | {:+.0}       | {:+.0}",
            v_gate, delta_tc, tmr_change
        );
    }

    println!("\n→ Electrical control without magnetic field!");
    println!("→ Low-power, fast switching\n");

    // === 5. Van der Waals Heterostructures ===
    println!("=== Van der Waals Heterostructures ===");

    println!("CrI₃/WSe₂ spin-valley coupling:");
    println!("  • CrI₃: Out-of-plane magnetization");
    println!("  • WSe₂: Valley Hall effect");
    println!("  • Proximity: M(CrI₃) → valley splitting in WSe₂\n");

    let valley_splitting = 30e-3; // 30 meV typical

    println!(
        "  Valley Zeeman splitting: {:.0} meV",
        valley_splitting * 1e3
    );
    println!(
        "  Equivalent B field: ~{:.0} T (giant!)",
        valley_splitting * E_CHARGE / (2.0 * MU_B)
    );
    println!("  → Valleytronic devices\n");

    // === 6. Spin Transport ===
    println!("=== Spin Transport in 2D Magnets ===");

    println!("Fe₃GeTe₂ as spin injector:");

    let j_current = 1e11; // A/m²
    let spin_polarization = 0.6; // 60%

    println!("  Current density: {:.0} MA/cm²", j_current * 1e-10);
    println!("  Spin polarization: {:.0}%", spin_polarization * 100.0);

    let j_spin = j_current * spin_polarization;
    println!("  Spin current: {:.0} MA/cm²\n", j_spin * 1e-10);

    // Spin diffusion length
    let lambda_sf_fgt = 50e-9; // ~50 nm
    let tau_sf = 1e-12; // 1 ps

    println!("  Spin diffusion length: {:.0} nm", lambda_sf_fgt * 1e9);
    println!("  Spin lifetime: {:.0} ps", tau_sf * 1e12);

    // === 7. Magnon Properties ===
    println!("\n=== Magnons in 2D Magnets ===");

    let d_layer = 0.7e-9; // CrI₃ layer thickness
    let j_exchange_2d = 2.0e-21; // J (exchange)

    println!("CrI₃ magnons:");
    println!("  Layer thickness: {:.2} nm", d_layer * 1e9);
    println!("  Exchange: {:.0} meV", j_exchange_2d / E_CHARGE * 1e3);

    // Magnon gap (anisotropy)
    let k_anis_2d = cri3.anisotropy_energy;
    let omega_gap = 2.0 * k_anis_2d / MU_B * 1e-3; // Simplified
    let f_gap = omega_gap / (2.0 * PI) * 1e-9;

    println!("  Magnon gap: {:.2} GHz (anisotropy)", f_gap);
    println!("  → Protects against thermal excitation\n");

    // === 8. Dzyaloshinskii-Moriya Interaction ===
    println!("=== DMI in 2D Materials ===");

    println!("MnBi₂Te₄/CrI₃ interface:");
    println!("  • Broken inversion symmetry");
    println!("  • Interfacial DMI");
    println!("  • Néel-type skyrmions possible\n");

    let d_dmi_2d = 0.5e-3; // 0.5 mJ/m²
    let a_ex_2d = 10e-12; // 10 pJ/m

    println!("  DMI constant: {:.1} mJ/m²", d_dmi_2d * 1e3);
    println!("  Exchange: {:.0} pJ/m", a_ex_2d * 1e12);

    let sk_size_2d = 4.0 * PI * a_ex_2d / d_dmi_2d * 1e9;
    println!("  Skyrmion diameter: {:.0} nm (ultra-small!)\n", sk_size_2d);

    // === 9. Thermal Stability ===
    println!("=== Thermal Stability ===");

    let kb = KB;
    let t_room = 300.0;

    println!("At room temperature ({:.0} K):", t_room);

    let e_anis_cri3 = k_anis_2d; // Per atom
    let _stability_cri3 = e_anis_cri3 / (kb * t_room);

    println!("  CrI₃:");
    println!("    T_C = {:.0} K", cri3.critical_temperature);
    if cri3.critical_temperature > t_room {
        println!("    ✓ Ferromagnetic at RT");
    } else {
        println!("    ✗ Paramagnetic at RT (needs cooling)");
    }

    println!("  Fe₃GeTe₂:");
    println!("    T_C = {:.0} K", fgt.critical_temperature);
    if fgt.critical_temperature > t_room {
        println!("    ✓ Ferromagnetic at RT");
    } else {
        println!("    ✗ Needs higher T_C variant");
    }

    println!("\n  Strategy: Substrate engineering, doping");
    println!("  Example: Fe₃GeTe₂ doped → T_C up to 400K\n");

    // === 10. Device Applications ===
    println!("=== Device Applications ===");

    println!("\n1. Ultra-Thin Memory:");
    println!("  • Thickness: < 2 nm (atomic limit)");
    println!("  • Area: {:.0} nm² per bit (dense)", 10.0 * 10.0);
    println!("  • TMR: {:.0}% (high readout)", tmr_max * 100.0);
    println!("  • Gate control: No external B field");

    println!("\n2. Flexible Spintronics:");
    println!("  • Van der Waals bonding → transferable");
    println!("  • Works on flexible substrates");
    println!("  • Strain tuning of magnetism");

    println!("\n3. Quantum Information:");
    println!("  • Localized spins (quantum dots)");
    println!("  • Valley qubits in heterostructures");
    println!("  • Long coherence from isolation");

    println!("\n4. Neuromorphic Computing:");
    println!("  • Analog resistance states (TMR tuning)");
    println!("  • Low switching energy (~fJ)");
    println!("  • Memristive behavior");

    // === 11. Comparison with 3D Materials ===
    println!("\n=== Advantages over Bulk Magnets ===");

    println!("\n2D Materials:");
    println!("  ✓ Atomic-scale thickness");
    println!("  ✓ Gate-tunable properties");
    println!("  ✓ Van der Waals integration");
    println!("  ✓ High spin polarization");
    println!("  ✓ Strong spin-orbit coupling");
    println!("  ✗ Often low T_C (improving)");

    println!("\nBulk Materials (CoFeB, Py):");
    println!("  ✓ Room-temperature FM");
    println!("  ✓ Well-established processing");
    println!("  ✗ Thick films (>2 nm)");
    println!("  ✗ Limited tunability");
    println!("  ✗ Requires lattice matching");

    // === Summary ===
    println!("\n=== Summary ===");

    println!("\nKey Features of 2D Magnets:");
    println!("  • Ferromagnetism at atomic thickness");
    println!("  • Layer-dependent properties");
    println!("  • Electric field control");
    println!("  • Van der Waals heterostructures");
    println!("  • Spin-valley coupling");

    println!("\nPromising Devices:");
    println!("  ✓ Ultra-thin MTJs (CrI₃/hBN/CrI₃)");
    println!("  ✓ Gate-controlled spin valves");
    println!("  ✓ Valleytronic transistors (CrI₃/TMD)");
    println!("  ✓ Nanoscale skyrmion memory");

    println!("\nCurrent Challenges:");
    println!("  • Raising T_C above room temperature");
    println!("  • Air stability (encapsulation needed)");
    println!("  • Large-area synthesis");
    println!("  • Contact resistance");

    println!("\nFuture Outlook:");
    println!("  → Materials discovery (ML-guided)");
    println!("  → Strain/twist engineering");
    println!("  → Moiré magnetism");
    println!("  → Integration with 2D semiconductors");
    println!("  → Topological 2D magnets");
}
