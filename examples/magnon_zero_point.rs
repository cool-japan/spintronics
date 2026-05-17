//! Magnon Zero-Point Fluctuations and Quantum Ground State
//!
//! **Difficulty**: ⭐⭐ Intermediate
//! **Category**: Quantum Magnonics
//!
//! **Physics**: Holstein-Primakoff quantization, zero-point motion,
//! Casimir magnon contribution, quantized spin waves in YIG stripe.
//!
//! This example demonstrates:
//!
//! 1. **Holstein-Primakoff transformation** for YIG (S = 10, ferromagnet):
//!    - S^z expectation value `⟨S^z⟩ = S − n` for n = 0, 1, 2, 3
//!    - Leading-order S^+/S^− coefficients = √(2S)
//!    - Quadratic 1/S correction at n = 2, 3
//!    - Validity fraction `n / (2S)` check
//!
//! 2. **Quantized spin-wave modes** in a YIG stripe (200 nm × 20 nm cross-section,
//!    1 µm long) computed via the Kalinikos-Slavin dipole-exchange formula with
//!    standing-wave quantisation.
//!
//! 3. **Zero-point fluctuations** of those modes: vacuum amplitude `u_n`, total
//!    ground-state energy, and Casimir magnon free energy F(T) at selected
//!    temperatures spanning 0 – 300 K.
//!
//! 4. **Bogoliubov transformation** for each mode, treating a 10% quantum-coupling
//!    fraction as the pairing term B_k = 0.1 A_k.  Reports squeeze parameter θ,
//!    vacuum occupation `v_k²`, and zero-point energy shift per mode.
//!
//! **References**:
//! - T. Holstein & H. Primakoff, Phys. Rev. **58**, 1098 (1940).
//! - C. Kittel, *Quantum Theory of Solids* (Wiley, 1963).
//! - A. V. Chumak et al., Nat. Phys. **11**, 453 (2015).

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Magnon Zero-Point Fluctuations and Quantum Ground State ===\n");

    // -----------------------------------------------------------------------
    // 1. Holstein-Primakoff transformation for YIG (S = 10)
    // -----------------------------------------------------------------------
    println!("=== 1. Holstein-Primakoff Transformation (YIG, S = 10) ===\n");

    // YIG has an effective spin S = 10 per formula unit due to five Fe3+ ions
    // contributing S_eff = 5/2 each, arranged into a net ferrimagnetic S = 10.
    // Holstein-Primakoff::new() requires a half-integer or integer spin.
    let spin_yig = 10.0_f64;
    let hp = HolsteinPrimakoff::new(spin_yig).expect("YIG spin S=10 is a valid integer spin value");

    println!("  Spin S = {}", hp.spin);
    println!("  Expansion order: {:?}", hp.expansion_order);
    println!("  Maximum boson occupation 2S = {}", hp.maximum_n_boson());
    println!(
        "  Validity warning threshold (10%% of S) = {:.1}",
        hp.validity_warning_threshold()
    );
    println!();

    println!(
        "  S^+ coefficient √(2S) = {:.6}  (linear HP)",
        hp.transform_splus_coeff()
    );
    println!(
        "  S^- coefficient √(2S) = {:.6}  (same by Hermitian conjugation)",
        hp.transform_sminus_coeff()
    );
    println!();

    println!(
        "  {:>6}  {:>12}  {:>12}  {:>14}  {:>14}",
        "n", "<S^z>", "n/(2S)", "Quad. corr.", "Valid?"
    );
    println!("  {}", "-".repeat(66));
    for n in [0u64, 1, 2, 3] {
        let n_f = n as f64;
        let sz = hp.transform_sz(n_f);
        let frac = hp.validity_fraction(n_f);
        let quad = hp.quadratic_correction(n_f);
        let valid = if frac < 0.1 { "YES" } else { "WARN" };
        println!(
            "  {:>6}  {:>12.4}  {:>12.6}  {:>14.8}  {:>14}",
            n, sz, frac, quad, valid
        );
    }
    println!();

    // Commutator check at low occupation (n=0, 0)
    let comm_ok = hp.commutator_check(0.0, 0.0, 1e-10);
    println!(
        "  Commutator [S^+, S^-] = 2S^z satisfied at n=0: {}",
        if comm_ok { "PASS" } else { "FAIL" }
    );

    // -----------------------------------------------------------------------
    // 2. Quantized spin-wave modes in a YIG stripe
    // -----------------------------------------------------------------------
    println!("\n=== 2. Quantized Spin-Wave Modes (YIG Stripe, 200 nm × 20 nm) ===\n");

    // YIG stripe geometry: 200 nm wide, 1 µm long, 20 nm thick
    // The thickness is passed to QuantizedModes::new(); Stripe takes width & length.
    let stripe_width = 200e-9_f64; // 200 nm
    let stripe_length = 1e-6_f64; // 1 µm
    let stripe_thickness = 20e-9_f64; // 20 nm

    let yig = Ferromagnet::yig();
    println!(
        "  YIG parameters: Ms = {:.3e} A/m, A_ex = {:.3e} J/m, alpha = {:.4}",
        yig.ms, yig.exchange_a, yig.alpha
    );

    let geometry = NanostructureGeometry::Stripe {
        width: stripe_width,
        length: stripe_length,
    };
    let qm = QuantizedModes::new(&yig, geometry, stripe_thickness)
        .expect("YIG stripe geometry is valid");

    // External field: 100 kA/m = mu_0 * H = 4pi * 1e-7 * 1e5 ≈ 0.12566 T
    // The API expects h_ext in Tesla (omega_H = GAMMA * h_ext)
    let h_ext_ka_per_m = 100e3_f64; // 100 kA/m
    let h_ext_tesla = MU_0 * h_ext_ka_per_m; // ≈ 0.1257 T
    let n_modes = 8_usize;

    println!(
        "  External field: {:.1} kA/m = {:.5} T",
        h_ext_ka_per_m * 1e-3,
        h_ext_tesla
    );
    println!(
        "  Stripe: width = {:.0} nm, length = {:.0} nm, thickness = {:.0} nm",
        stripe_width * 1e9,
        stripe_length * 1e9,
        stripe_thickness * 1e9
    );
    println!("  Number of quantized modes: {}", n_modes);
    println!();

    let mode_freqs: Vec<(usize, f64)> = qm
        .stripe_mode_frequencies(h_ext_tesla, n_modes)
        .expect("stripe mode frequencies computed");

    println!(
        "  {:>6}  {:>18}  {:>16}  {:>12}",
        "Mode n", "omega_n (rad/s)", "f_n (GHz)", "k_n (1/m)"
    );
    println!("  {}", "-".repeat(60));
    let wavevectors = qm
        .stripe_wavevectors(n_modes)
        .expect("stripe wavevectors computed");

    for (i, &(n, omega)) in mode_freqs.iter().enumerate() {
        let f_ghz = omega / (2.0 * std::f64::consts::PI * 1e9);
        let k_n = wavevectors.get(i).copied().unwrap_or(0.0);
        println!(
            "  {:>6}  {:>18.6e}  {:>16.6}  {:>12.4e}",
            n, omega, f_ghz, k_n
        );
    }

    // -----------------------------------------------------------------------
    // 3. Zero-point fluctuations
    // -----------------------------------------------------------------------
    println!("\n=== 3. Zero-Point Fluctuations (Vacuum Amplitudes & Casimir Energy) ===\n");

    // YIG magnon effective mass and mode volume
    // m_eff ≈ 1.5e-28 kg  (from exchange stiffness / group velocity)
    // volume = width × thickness × length (200 nm × 20 nm × 1 µm)
    let m_eff = 1.5e-28_f64;
    let volume = stripe_width * stripe_thickness * stripe_length;
    println!("  Magnon effective mass m_eff = {:.2e} kg", m_eff);
    println!("  Mode volume V = {:.4e} m^3", volume);
    println!("  Spin S = {} (YIG effective spin)", spin_yig);
    println!();

    let zpf = ZeroPointFluctuations::new(mode_freqs.clone(), m_eff, spin_yig, volume)
        .expect("ZeroPointFluctuations constructed from YIG stripe modes");

    // Zero-point amplitude per mode
    println!(
        "  {:>6}  {:>18}  {:>18}  {:>18}",
        "Mode n", "omega_n (rad/s)", "u_n (pm)", "u_n^2 (m^2)"
    );
    println!("  {}", "-".repeat(68));
    for (i, &(n, omega)) in mode_freqs.iter().enumerate() {
        let u_n = zpf
            .zero_point_amplitude(i)
            .expect("valid mode index for zero-point amplitude");
        println!(
            "  {:>6}  {:>18.6e}  {:>18.6e}  {:>18.6e}",
            n,
            omega,
            u_n * 1e12, // picometres
            u_n * u_n
        );
    }
    println!();

    // Ground-state (vacuum) energy
    let e_ground = zpf.ground_state_energy();
    println!(
        "  Ground-state energy E_0 = (ℏ/2) Σ_n ω_n = {:.6e} J",
        e_ground
    );
    println!(
        "  Ground-state energy E_0 = {:.4e} eV",
        e_ground / 1.602_176_634e-19
    );
    println!();

    // Vacuum fluctuation spin density
    let sx_vac = zpf
        .vacuum_fluctuation_sx()
        .expect("vacuum fluctuation Sx density computed");
    println!(
        "  Vacuum fluctuation density Σ u_n^2 / V = {:.4e} m^-1",
        sx_vac
    );
    println!();

    // Casimir magnon free energy at selected temperatures
    let temperatures = [0.0_f64, 1.0, 4.2, 10.0, 77.0, 300.0];
    println!(
        "  {:>10}  {:>20}  {:>20}  {:>12}",
        "T (K)", "F_Casimir (J)", "F_Casimir (eV)", "F / E_0"
    );
    println!("  {}", "-".repeat(70));
    for &temp in &temperatures {
        let f_casimir = zpf
            .casimir_free_energy(temp)
            .expect("Casimir free energy computed");
        let f_ev = f_casimir / 1.602_176_634e-19;
        let ratio = f_casimir / e_ground;
        println!(
            "  {:>10.1}  {:>20.8e}  {:>20.8e}  {:>12.6}",
            temp, f_casimir, f_ev, ratio
        );
    }

    // Casimir pressure (confinement force) at T = 0
    let dl = 1e-12_f64; // 1 pm perturbation
    let casimir_force = zpf
        .casimir_pressure_change(stripe_length, dl)
        .expect("Casimir pressure change computed");
    println!();
    println!(
        "  Casimir confinement force dF/dL ≈ {:.4e} N  (T=0, dL=1 pm)",
        casimir_force
    );

    // -----------------------------------------------------------------------
    // 4. Bogoliubov transformation per mode
    // -----------------------------------------------------------------------
    println!("\n=== 4. Bogoliubov Transformation (10%% Quantum Coupling per Mode) ===\n");

    println!("  Each mode treated as a single-mode system: A_k = ω_n, B_k = 0.1 × A_k");
    println!("  ω_k (Bogoliubov) = √(A_k² − B_k²) = A_k √(1 − 0.01) ≈ 0.99499 A_k");
    println!();

    println!(
        "  {:>6}  {:>16}  {:>14}  {:>16}  {:>16}  {:>16}",
        "Mode n",
        "omega_bogo (rad/s)",
        "theta (squeeze)",
        "vac. occ. v_k^2",
        "E_0 shift (J)",
        "norm ok?"
    );
    println!("  {}", "-".repeat(92));

    // Collect per-mode Bogoliubov data for summary statistics
    let mut total_vac_magnons = 0.0_f64;
    let mut total_e0_shift = 0.0_f64;

    for &(n, omega_n) in &mode_freqs {
        // Single-mode Hamiltonian: a_k = [omega_n], b_k = [0.1 * omega_n]
        let a_k = vec![omega_n];
        let b_k = vec![0.1 * omega_n];

        let bogo = BogoliubovTransform::from_hamiltonian(&a_k, &b_k)
            .expect("Bogoliubov transform: B_k = 0.1 A_k is stable (|B| < |A|)");

        let omega_bogo = bogo.magnon_frequency(0).expect("mode 0 frequency");
        let theta = bogo
            .squeeze_parameter(0)
            .expect("squeeze parameter for mode 0");
        let vac_occ = bogo
            .vacuum_occupation(0)
            .expect("vacuum occupation for mode 0");
        let e0_shift = bogo.ground_state_energy();
        let norm_ok = bogo.check_normalization(1e-10);

        total_vac_magnons += vac_occ;
        total_e0_shift += e0_shift;

        println!(
            "  {:>6}  {:>16.6e}  {:>14.8}  {:>16.10}  {:>16.8e}  {:>16}",
            n,
            omega_bogo,
            theta,
            vac_occ,
            // Convert the dimensionless energy shift (in rad/s) to Joules via ℏ
            e0_shift * HBAR,
            if norm_ok { "PASS" } else { "FAIL" }
        );
    }

    println!();
    println!(
        "  Total vacuum magnon depletion across all modes: {:.6e}",
        total_vac_magnons
    );
    println!(
        "  Total Bogoliubov ground-state energy shift:     {:.6e} J",
        total_e0_shift * HBAR
    );
    println!("  (negative means quantum correlations lower the ground state)");

    // -----------------------------------------------------------------------
    // Summary
    // -----------------------------------------------------------------------
    println!("\n=== Summary ===\n");
    println!("Holstein-Primakoff (YIG, S = 10):");
    println!(
        "  S^+/S^- coefficient = {:.4}  (spin-wave amplitude ∝ √(2S))",
        hp.transform_splus_coeff()
    );
    println!(
        "  Validity fraction at n=1: {:.4} << 1  (linear approximation excellent)",
        hp.validity_fraction(1.0)
    );

    println!("\nQuantized spin-wave modes (n = 1..{}):", n_modes);
    if let Some(&(_, omega_1)) = mode_freqs.first() {
        let f1_ghz = omega_1 / (2.0 * std::f64::consts::PI * 1e9);
        println!("  Lowest mode (n=1): f_1 = {:.4} GHz", f1_ghz);
    }
    if let Some(&(_, omega_n)) = mode_freqs.last() {
        let fn_ghz = omega_n / (2.0 * std::f64::consts::PI * 1e9);
        println!(
            "  Highest mode (n={}): f_{} = {:.4} GHz",
            n_modes, n_modes, fn_ghz
        );
    }

    println!("\nZero-point fluctuations:");
    println!(
        "  Ground-state vacuum energy E_0 = {:.4e} J  ({} modes)",
        e_ground,
        zpf.n_modes()
    );
    println!(
        "  Casimir free energy at 300 K vs 0 K ratio: {:.4}",
        zpf.casimir_free_energy(300.0).expect("300 K free energy")
            / zpf.casimir_free_energy(0.0).expect("0 K free energy")
    );

    println!("\nBogoliubov squeezed vacuum (B_k = 0.1 A_k):");
    println!(
        "  Total virtual magnon depletion (all modes): {:.4e}",
        total_vac_magnons
    );
    println!(
        "  Total ground-state energy lowering:         {:.4e} J",
        total_e0_shift * HBAR
    );

    println!("\nPhysical interpretation:");
    println!("  - Zero-point fluctuations set the quantum noise floor for magnon detection.");
    println!("  - The Casimir-like magnon energy grows with temperature as magnon pairs are");
    println!("    thermally activated, eventually reaching the classical limit F ≈ Nk_BT.");
    println!("  - Bogoliubov squeezing mixes creation and annihilation operators, creating");
    println!("    virtual magnon pairs in the vacuum and enhancing quantum spin noise beyond");
    println!("    the Poisson (shot-noise) limit (Kamra-Belzig effect, PRL 2016).");

    Ok(())
}
