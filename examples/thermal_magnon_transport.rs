//! Thermal Magnon Transport and Spin Caloritronics
//!
//! This example demonstrates:
//! - Magnon-mediated thermal transport
//! - Spin Seebeck effect
//! - Magnon drag of electrons
//! - Thermal spin-orbit torques
//!
//! ## Physics:
//! Thermal gradients in magnetic materials excite magnon currents that
//! can transport spin angular momentum and heat. This enables thermoelectric
//! devices without charge current and energy-harvesting applications.
//!
//! ## References:
//! - K. Uchida et al., "Observation of the spin Seebeck effect", Nature 455, 778 (2008)
//! - G. E. W. Bauer et al., "Spin caloritronics", Nature Mater. 11, 391 (2012)
//!
//! Run with: cargo run --example thermal_magnon_transport

use spintronics::prelude::*;
use std::f64::consts::PI;

fn main() {
    println!("=== Thermal Magnon Transport Simulation ===\n");

    // === 1. Material System ===
    println!("=== Material System ===");

    // Using Permalloy (low damping, good magnon conductor)
    let yig = Ferromagnet::permalloy();
    let pt = InverseSpinHall::platinum();
    let _interface = SpinInterface::py_pt();

    println!("Magnetic material: Ni₈₁Fe₁₉ (Permalloy)");
    println!("  M_s = {:.0} kA/m", yig.ms * 1e-3);
    println!("  α = {:.5} (low damping)", yig.alpha);
    println!("  Electrical resistivity: ~16 µΩ·cm");

    println!("\nSpin detector: Pt");
    println!("  Spin Hall angle: {:.2}", pt.theta_sh);
    println!("  Resistivity: {:.1} µΩ·cm\n", pt.rho * 1e8);

    // === 2. Device Geometry ===
    println!("=== Device Structure ===");

    let yig_thickness = 1e-6; // 1 µm YIG film
    let pt_thickness = 5e-9; // 5 nm Pt layer
    let length = 5e-3; // 5 mm
    let width = 2e-3; // 2 mm

    println!("  YIG thickness: {:.1} µm", yig_thickness * 1e6);
    println!("  Pt thickness: {:.0} nm", pt_thickness * 1e9);
    println!("  Sample: {:.0} × {:.0} mm²\n", length * 1e3, width * 1e3);

    // === 3. Temperature Gradient ===
    println!("=== Thermal Conditions ===");

    let t_hot = 320.0; // K
    let t_cold = 300.0; // K
    let delta_t = t_hot - t_cold;
    let grad_t = delta_t / length; // K/m

    println!("  Hot side: {:.0} K", t_hot);
    println!("  Cold side: {:.0} K", t_cold);
    println!("  ΔT = {:.0} K", delta_t);
    println!("  ∇T = {:.2} K/mm\n", grad_t * 1e-3);

    // === 4. Magnon Thermal Conductivity ===
    println!("=== Magnon Thermal Transport ===");

    // YIG magnon thermal conductivity
    let kappa_m = magnon_thermal_conductivity(&yig, 300.0);
    let q_magnon = kappa_m * grad_t; // Heat flux

    println!("  Magnon thermal conductivity: {:.2} W/(m·K)", kappa_m);
    println!("  Magnon heat flux: {:.2e} W/m²", q_magnon);

    // Total heat flow
    let area = width * yig_thickness;
    let q_total = q_magnon * area;

    println!("  Total heat current: {:.2} mW\n", q_total * 1e3);

    // === 5. Spin Seebeck Effect ===
    println!("=== Spin Seebeck Effect ===");

    let sse = SpinSeebeck::yig_pt();

    // Spin current from SSE
    let grad_t_vec = Vector3::new(grad_t, 0.0, 0.0); // Temperature gradient along x
    let js_sse_vec = sse.spin_current(grad_t_vec);
    let js_sse = js_sse_vec.magnitude();

    println!("  Spin Seebeck coefficient: {:.2e} J/(K·m²)", sse.l_s);
    println!("  Spin current density: {:.2e} A/m²", js_sse);

    // ISHE voltage
    let _interface_area = length * width;
    let sigma_direction = Vector3::new(1.0, 0.0, 0.0); // Spin polarization
    let e_ishe = pt.convert(sigma_direction, js_sse_vec);

    let v_sse = e_ishe.magnitude() * width;

    println!("  ISHE electric field: {:.2e} V/m", e_ishe.magnitude());
    println!("  SSE voltage: {:.2} µV\n", v_sse * 1e6);

    // === 6. Magnon Drag ===
    println!("=== Magnon Drag Effect ===");

    // Magnons drag electrons via s-d exchange
    let j_electron = 1e10; // A/m² test current
    let magnon_number = thermal_magnon_density(&yig, 300.0);

    println!("  Thermal magnon density: {:.2e} m⁻³", magnon_number);

    // Drag resistance (simplified)
    let r_drag = magnon_drag_resistivity(&yig, 300.0, grad_t);
    println!("  Magnon drag resistivity: {:.2e} Ω·m", r_drag);

    // Additional voltage from drag
    let v_drag = r_drag * j_electron * length;
    println!(
        "  Drag voltage (at j = {:.0} MA/cm²): {:.2} µV\n",
        j_electron * 1e-10,
        v_drag * 1e6
    );

    // === 7. Thermal Spin-Orbit Torque ===
    println!("=== Thermal Spin-Orbit Torque ===");

    // Temperature gradient induces SOT even without charge current
    let sot_thermal_coeff = 1e-12; // Simplified coefficient
    let torque_thermal = sot_thermal_coeff * grad_t;

    println!("  Thermal SOT coefficient: {:.2e} T·m/K", sot_thermal_coeff);
    println!("  Effective field from ∇T: {:.2e} T", torque_thermal);
    println!("  Can switch magnetization without current!\n");

    // === 8. Spin Nernst Effect ===
    println!("=== Spin Nernst Effect ===");

    let sne = SpinNernst::tungsten();

    // Transverse spin current from temperature gradient
    let grad_t_vec = Vector3::new(grad_t, 0.0, 0.0);
    let spin_dir = Vector3::new(0.0, 1.0, 0.0);
    let js_nernst = sne.spin_current(grad_t_vec, spin_dir);

    println!("  Spin Nernst angle: {:.4}", sne.spin_nernst_angle);
    println!(
        "  Transverse spin current: {:.2e} A/m²",
        js_nernst.magnitude()
    );
    println!("  Direction: perpendicular to ∇T and σ\n");

    // === 9. Energy Conversion Efficiency ===
    println!("=== Thermoelectric Performance ===");

    // Seebeck coefficient (voltage/temperature)
    let s_eff = v_sse / delta_t;
    println!("  Effective Seebeck: {:.2} µV/K", s_eff * 1e6);

    // Figure of merit ZT (simplified)
    let kappa_total = kappa_m + 1.0; // Include phonon contribution
    let sigma_eff = 1.0 / (pt.rho * pt_thickness);
    let zt = (s_eff * s_eff * sigma_eff * 300.0) / kappa_total;

    println!("  ZT figure of merit: {:.4} (low, YIG is insulator)", zt);
    println!("  Better for spin-based devices than power generation\n");

    // === 10. Comparison with Conventional Thermoelectrics ===
    println!("=== Comparison with Conventional TE ===");

    println!("\nConventional (Bi₂Te₃):");
    println!("  Seebeck: ~200 µV/K");
    println!("  ZT ~ 1.0 at 300K");
    println!("  Requires charge current");

    println!("\nSpin Seebeck (YIG/Pt):");
    println!("  Effective Seebeck: {:.1} µV/K", s_eff * 1e6);
    println!("  ZT ~ {:.4}", zt);
    println!("  Pure spin current (insulator)");
    println!("  Advantages:");
    println!("    ✓ No Joule heating");
    println!("    ✓ Works in insulators");
    println!("    ✓ Large-area devices");

    // === 11. Applications ===
    println!("\n=== Applications ===");

    println!("\n1. Waste Heat Recovery:");
    println!("  • Convert waste heat → electricity");
    println!("  • No charge current losses");
    println!("  • Scalable to large areas");

    println!("\n2. Thermal Sensors:");
    let sensitivity = v_sse / delta_t * 1e6; // µV/K
    println!("  • Sensitivity: {:.2} µV/K", sensitivity);
    println!("  • Non-invasive (no current)");
    println!("  • Fast response (magnon dynamics)");

    println!("\n3. Thermal Switching:");
    let switching_field = torque_thermal;
    println!("  • Temperature-controlled magnetization");
    println!("  • Effective field: {:.2e} T", switching_field);
    println!("  • Low-power magnetic memory");

    // === 12. Material Optimization ===
    println!("\n=== Material Engineering ===");

    println!("\nEnhancing spin Seebeck:");
    println!("  • Lower damping → longer magnon lifetime");
    println!("  • Higher M_s → stronger coupling");
    println!("  • Thinner films → larger ∇T");
    println!("  • Better interface → efficient spin injection");

    println!("\nBest materials:");
    println!("  ✓ YIG: α = {:.5} (champion)", yig.alpha);
    println!("  • Permalloy: α = 0.01 (good)");
    println!("  • CoFeB: α = 0.004 (very good)");

    // === Summary ===
    println!("\n=== Summary ===");

    println!("\nKey Results:");
    println!(
        "  Temperature gradient: {:.0} K over {:.0} mm",
        delta_t,
        length * 1e3
    );
    println!("  Spin current (SSE): {:.2e} A/m²", js_sse);
    println!("  ISHE voltage: {:.2} µV", v_sse * 1e6);
    println!("  Heat flux: {:.2e} W/m²", q_magnon);

    println!("\nPhysical Insights:");
    println!("  • Magnons transport heat and spin independently");
    println!("  • Insulators (YIG) avoid charge current losses");
    println!("  • Thermal gradients can switch magnetization");
    println!("  • Promising for energy harvesting and sensing");

    println!("\nFuture Directions:");
    println!("  → Higher ZT materials (doping, nanostructuring)");
    println!("  → Magnonic circuits for heat management");
    println!("  → Thermal magnon logic");
    println!("  → Spin caloritronic quantum devices");
}

// Helper functions

fn magnon_thermal_conductivity(mat: &Ferromagnet, temperature: f64) -> f64 {
    // Simplified magnon thermal conductivity
    // κ_m ~ C_m v_m λ_m
    // where C_m = heat capacity, v_m = magnon velocity, λ_m = mean free path

    let kb = KB;
    let _hbar = HBAR;
    let gamma = GAMMA;

    // Typical magnon velocity ~1 km/s
    let v_magnon = 1000.0;

    // Mean free path from damping
    let lambda_m = v_magnon / (mat.alpha * gamma * mat.ms);

    // Magnon heat capacity (simplified)
    let n_magnons = thermal_magnon_density(mat, temperature);
    let c_magnon = n_magnons * kb;

    // Thermal conductivity
    c_magnon * v_magnon * lambda_m / 3.0
}

fn thermal_magnon_density(mat: &Ferromagnet, temperature: f64) -> f64 {
    // Thermal magnon density from Bose-Einstein statistics
    let kb = KB;
    let hbar = HBAR;
    let gamma = GAMMA;

    // Minimum magnon energy (h * f_FMR)
    let mu0 = 4.0 * PI * 1e-7;
    let h_internal = mat.anisotropy_k / mat.ms;
    let omega_min = gamma * mu0 * h_internal;

    let x = hbar * omega_min / (kb * temperature);

    // For x << 1 (high T), n ~ kT/(ħω)
    if x < 0.1 {
        let n0 = kb * temperature / (hbar * omega_min);
        n0 * 1e20 // Scale to reasonable density
    } else {
        1e18 / (x.exp() - 1.0)
    }
}

fn magnon_drag_resistivity(_mat: &Ferromagnet, _temperature: f64, grad_t: f64) -> f64 {
    // Magnon drag contribution to resistivity
    // ρ_drag ~ (∇T / T)² × ρ₀

    let rho0 = 1e-7; // Base resistivity (Ω·m)
    let t_avg = 300.0;

    rho0 * (grad_t / t_avg).powi(2)
}
