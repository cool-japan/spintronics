//! Axion Electrodynamics and Magnon-Photon Coupling
//!
//! **Difficulty**: ⭐⭐⭐⭐
//! **Category**: Topological Magnon Physics
//! **Physics**: Axion insulators, topological magnetoelectric effect,
//!               magnon-photon conversion, parity violation in hybrid systems
//!
//! ## Background
//!
//! The axion field θ characterises topological insulators via the action term
//! `L_axion = (θ e²/4π²ℏ) E·B`.  For a topological phase θ = π, giving a
//! quantised magnetoelectric polarisability α_TME = e²/(2h).
//!
//! This example demonstrates:
//!
//! 1. Computing the axion angle θ for a CubicHaldane 3D magnon lattice.
//! 2. Repeating the analysis for the PyrochloreTopological 4-band model.
//! 3. Computing the emergent polarisation P and magnetisation M induced by
//!    externally applied E and B fields (axion response).
//! 4. Modelling a hybrid axion-magnon-photon system (YIG cavity) and sweeping
//!    θ from 0 to π to illustrate how coupling efficiency and Faraday angle
//!    vary with the axion angle.
//!
//! ## References
//!
//! - Wilczek, *Phys. Rev. Lett.* **58**, 1799 (1987)
//! - Essin, Moore & Vanderbilt, *Phys. Rev. Lett.* **102**, 146805 (2009)
//! - Tokura, Yasuda & Tsukazaki, *Rev. Mod. Phys.* **91**, 015005 (2019)
//! - Li et al., *Phys. Rev. Lett.* **124**, 167402 (2020)

use spintronics::prelude::*;
use spintronics::topomagnon::{AxionElectrodynamics, AxionMagnonPhoton, MagnonBandModel3D};
use std::f64::consts::PI;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    // ─────────────────────────────────────────────────────────────────────────
    // Section 1: CubicHaldane 3D lattice
    // ─────────────────────────────────────────────────────────────────────────
    println!("======================================================");
    println!(" Section 1: Cubic Haldane 3D lattice (axion angle)");
    println!("======================================================");

    // Build a 2-band CubicHaldane model.
    // Parameters are in normalised units (J_nn = 1 sets the energy scale).
    //   j_nn  = 1.0  (NN exchange, sets bandwidth)
    //   j_nnn = 0.2  (NNN exchange, stored for reference)
    //   dmi   = 0.5  (DMI strength — drives topological phase transition)
    //   h_ext = 0.1  (uniform Zeeman field breaking Γ-point degeneracy)
    let cubic_model = MagnonBandModel3D::cubic_haldane(1.0, 0.2, 0.5, 0.1);

    println!("Model: CubicHaldane");
    println!("  n_bands = {}", cubic_model.n_bands());
    println!("  j_nn    = {:.2}", cubic_model.j_nn);
    println!("  j_nnn   = {:.2}", cubic_model.j_nnn);
    println!("  dmi     = {:.2}", cubic_model.dmi);
    println!("  h_ext   = {:.2}", cubic_model.h_ext);
    println!(
        "  is_topological_cubic = {} (|D|={:.2} vs |J_nn|={:.2})",
        cubic_model.is_topological_cubic(),
        cubic_model.dmi.abs(),
        cubic_model.j_nn.abs()
    );

    // Band gap sampled on a coarse 8-point mesh.
    let gap_cubic = cubic_model.band_gap_3d(8);
    println!("  band_gap (8-pt mesh) = {:.4}", gap_cubic);

    // Construct the axion electrodynamics calculator.
    // Occupied band = 0 (the lower of the two bands).
    // k-mesh: n_kx = n_ky = n_kz = 8 (minimum 4 required; 8 gives reasonable accuracy).
    let axion_cubic = AxionElectrodynamics::new(&cubic_model, vec![0], 8, 8, 8)?;

    let theta_cubic = axion_cubic.axion_angle();
    let alpha_tme_cubic = axion_cubic.topological_magnetoelectric_polarizability();

    println!("\n  -- Axion analysis --");
    println!(
        "  axion angle θ = {:.4} rad  = {:.4}π",
        theta_cubic,
        theta_cubic / PI
    );
    println!(
        "  topological phase: {}",
        if (theta_cubic - PI).abs() < 1e-6 {
            "θ = π  (axion insulator)"
        } else {
            "θ = 0  (trivial)"
        }
    );
    println!(
        "  α_TME = {:.4e} S  (e²/(2h) = {:.4e} S if θ=π)",
        alpha_tme_cubic,
        CONDUCTANCE_QUANTUM * 0.25
    );

    // ─────────────────────────────────────────────────────────────────────────
    // Section 2: Pyrochlore topological model
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n======================================================");
    println!(" Section 2: Pyrochlore topological model (4-band)");
    println!("======================================================");

    // The pyrochlore model hosts a topological axion state with θ = π at strong DMI.
    // j = 1.0 (NN exchange), d = 0.3 (DMI coupling strength).
    let pyro_model = MagnonBandModel3D::pyrochlore_topological(1.0, 0.3);

    println!("Model: PyrochloreTopological");
    println!("  n_bands = {}", pyro_model.n_bands());
    println!("  j_nn    = {:.2}", pyro_model.j_nn);
    println!("  dmi     = {:.2}", pyro_model.dmi);

    // Diagonalise at Γ-point to show the 4-band spectrum.
    let (evals_gamma, _) = pyro_model.diagonalize_3d(0.0, 0.0, 0.0)?;
    println!(
        "  Eigenvalues at Γ: {:?}",
        evals_gamma
            .iter()
            .map(|e| format!("{:.4}", e))
            .collect::<Vec<_>>()
    );

    // For the pyrochlore model (4 bands), we declare the lower half (0, 1) as occupied.
    let axion_pyro = AxionElectrodynamics::new(&pyro_model, vec![0, 1], 8, 8, 8)?;

    let theta_pyro = axion_pyro.axion_angle();
    let alpha_tme_pyro = axion_pyro.topological_magnetoelectric_polarizability();

    println!("\n  -- Axion analysis --");
    println!(
        "  axion angle θ = {:.4} rad  = {:.4}π",
        theta_pyro,
        theta_pyro / PI
    );
    println!(
        "  topological phase: {}",
        if (theta_pyro - PI).abs() < 1e-6 {
            "θ = π  (axion insulator)"
        } else {
            "θ = 0  (trivial)"
        }
    );
    println!("  α_TME = {:.4e} S", alpha_tme_pyro);

    // ─────────────────────────────────────────────────────────────────────────
    // Section 3: Axion electromagnetic response
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n======================================================");
    println!(" Section 3: Axion electromagnetic response");
    println!("======================================================");

    // Apply the axion response using the pyrochlore model (which has θ = π,
    // giving a non-trivial magnetoelectric response).
    // E = (0, 0, 1e5) V/m and B = (0, 0, 1) T along z.
    let e_field = Vector3::new(0.0_f64, 0.0, 1.0e5); // V/m
    let b_field = Vector3::new(0.0_f64, 0.0, 1.0); // T

    println!("Applied fields:");
    println!(
        "  E = ({:.2e}, {:.2e}, {:.2e}) V/m",
        e_field.x, e_field.y, e_field.z
    );
    println!(
        "  B = ({:.2e}, {:.2e}, {:.2e}) T",
        b_field.x, b_field.y, b_field.z
    );
    println!("\n-- Using pyrochlore model (θ = π): --");

    let (polarization, magnetization) = axion_pyro.axion_response(e_field, b_field);

    // Convert to more readable units.
    let p_uc_m2 = polarization.z * 1.0e6; // C/m² → μC/m²
    let m_a_m = magnetization.z; // already in A/m

    println!(
        "Axion response (θ={:.4}π, α_TME={:.4e} S):",
        theta_pyro / PI,
        alpha_tme_pyro
    );
    println!(
        "  Emergent P = ({:.4e}, {:.4e}, {:.4e}) C/m²",
        polarization.x, polarization.y, polarization.z
    );
    println!(
        "    P_z = {:.6e} μC/m²  (= α_TME·B_z = {:.4e}·{:.1})",
        p_uc_m2, alpha_tme_pyro, b_field.z
    );
    println!(
        "  Emergent M = ({:.4e}, {:.4e}, {:.4e}) A/m",
        magnetization.x, magnetization.y, magnetization.z
    );
    println!(
        "    M_z = {:.4} A/m  (= -α_TME·E_z = -{:.4e}·{:.2e})",
        m_a_m, alpha_tme_pyro, e_field.z
    );

    println!("\n-- Using trivial cubic model (θ = 0): --");
    let (p_triv, m_triv) = axion_cubic.axion_response(e_field, b_field);
    println!("  P_z = {:.2e} C/m²  (expected: 0, trivial)", p_triv.z);
    println!("  M_z = {:.2e} A/m  (expected: 0, trivial)", m_triv.z);

    println!("\nConclusion:");
    println!("  Topological (θ=π): P_z = α_TME·B ≠ 0  and  M_z = -α_TME·E ≠ 0.");
    println!("  Trivial     (θ=0): P_z = M_z = 0  — no magnetoelectric response.");

    // ─────────────────────────────────────────────────────────────────────────
    // Section 4: Axion magnon-photon coupling — YIG topological cavity
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n======================================================");
    println!(" Section 4: Axion magnon-photon coupling (YIG cavity)");
    println!("======================================================");

    // Use the built-in YIG topological cavity preset.
    // This models YIG-FMR at ~10 GHz in a microwave cavity with the axion
    // angle tuned to θ = π (topological phase).
    let yig_sys = AxionMagnonPhoton::yig_topological_cavity();

    let g_eff = yig_sys.effective_coupling();
    let detuning = yig_sys.detuning();
    let eta = yig_sys.magnon_photon_conversion_efficiency();
    let faraday_rad = yig_sys.parity_violation_angle();
    let cooperativity = yig_sys.cooperativity();

    println!("YIG topological cavity (preset):");
    println!(
        "  cavity_freq  ω_c = {:.4e} rad/s  ({:.2} GHz)",
        yig_sys.cavity_freq,
        yig_sys.cavity_freq / (2.0 * PI * 1e9)
    );
    println!(
        "  magnon_freq  ω_m = {:.4e} rad/s  ({:.2} GHz)",
        yig_sys.magnon_freq,
        yig_sys.magnon_freq / (2.0 * PI * 1e9)
    );
    println!("  axion_coupling_g = {:.4e}", yig_sys.axion_coupling_g);
    println!(
        "  theta_axion      = {:.4} rad = {:.4}π",
        yig_sys.theta_axion,
        yig_sys.theta_axion / PI
    );
    println!(
        "  kappa_c          = {:.4e} rad/s  ({:.2} MHz)",
        yig_sys.kappa_c,
        yig_sys.kappa_c / (2.0 * PI * 1e6)
    );
    println!(
        "  gamma_m          = {:.4e} rad/s  ({:.2} MHz)",
        yig_sys.gamma_m,
        yig_sys.gamma_m / (2.0 * PI * 1e6)
    );
    println!("\nDerived quantities:");
    println!("  effective coupling g_eff       = {:.4e} rad/s", g_eff);
    println!("  cavity-magnon detuning Δ       = {:.4e} rad/s", detuning);
    println!("  conversion efficiency η        = {:.6e}", eta);
    println!(
        "  Faraday angle (parity-viol.)   = {:.6e} rad  = {:.4} mrad",
        faraday_rad,
        faraday_rad * 1000.0
    );
    println!("  cooperativity C = 4g²/(κγ)    = {:.4e}", cooperativity);

    // ── θ sweep: table of (θ/π, η, Faraday angle mrad)
    println!("\nSweep θ_axion from 0 to π:");
    println!(
        "{:>8}  {:>14}  {:>18}",
        "θ/π", "η (efficiency)", "Faraday (mrad)"
    );
    println!("{}", "-".repeat(46));

    let n_sweep = 11;
    for i in 0..n_sweep {
        let theta = PI * (i as f64) / ((n_sweep - 1) as f64);

        // Build a new AxionMagnonPhoton with the same cavity params but varying θ.
        let sys = AxionMagnonPhoton::new(
            yig_sys.axion_coupling_g,
            yig_sys.cavity_freq,
            yig_sys.magnon_freq,
            yig_sys.kappa_c,
            yig_sys.gamma_m,
            theta,
        )?;
        let eta_i = sys.magnon_photon_conversion_efficiency();
        let f_i = sys.parity_violation_angle() * 1000.0; // mrad

        println!("{:>8.4}  {:>14.6e}  {:>18.6}", theta / PI, eta_i, f_i);
    }

    println!("\nPhysical interpretation:");
    println!("  - At θ=0 (trivial), g_eff=0, η→0, Faraday angle=0.");
    println!("  - At θ=π (topological), g_eff=g_aγγ, η is maximised on resonance.");
    println!("  - The Faraday angle θ_F ≈ θ·α_FS/2 reveals the axion-photon parity violation.");
    println!("  - Strong coupling C>>1 needed for quantum-coherent transduction.");

    Ok(())
}
