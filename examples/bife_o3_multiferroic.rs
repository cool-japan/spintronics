//! BiFeO3 Magnetoelectric Coupling and Electric-Field Control of Magnetism
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Multiferroics / Magnetoelectric Effects
//! **Physics**: Linear magnetoelectric effect; electric polarisation switching;
//! DM-driven polarisation; exchange striction; toroidal moments; Cr₂O₃ comparison.
//!
//! ## Background
//!
//! BiFeO₃ (BFO) is the archetypal room-temperature multiferroic with simultaneous:
//! - Ferroelectric order: P_s ≈ 90 μC/cm² along (111) (T_FE = 1100 K)
//! - G-type antiferromagnetic order (T_N = 643 K) with cycloidal modulation
//!
//! It belongs to Type-I (independent order parameters) — unlike TbMnO₃ (Type-II)
//! where magnetism drives ferroelectricity.  The weak cross-coupling enters through
//! the linear ME tensor α_ij and DM-driven polarisation from the residual canting.
//!
//! ## References
//!
//! - J. Wang et al., "Epitaxial BiFeO3 Multiferroic Thin Film Heterostructures",
//!   *Science* **299**, 1719 (2003).
//! - D. Lebeugle et al., "Electric-field-induced spin flop in BiFeO3 single crystals",
//!   *Phys. Rev. Lett.* **100**, 227602 (2008).
//! - G. Catalan & J. F. Scott, "Physics and Applications of Bismuth Ferrite",
//!   *Adv. Mater.* **21**, 2463 (2009).
//! - I. E. Dzyaloshinskii, "On the magneto-electrical effect in antiferromagnets",
//!   *JETP* **10**, 628 (1960).

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  BiFeO3 Multiferroic: Magnetoelectric Coupling & Field Control");
    println!("=================================================================\n");

    // -------------------------------------------------------------------------
    // 1. Material database: BiFeO3, TbMnO3, Cr2O3
    // -------------------------------------------------------------------------
    let bfo = MagnetoelectricTensor::bife_o3();
    let tbmno3 = MagnetoelectricTensor::tb_mn_o3();
    let cr2o3 = MagnetoelectricTensor::cr2_o3();

    println!("--- Material Overview ---");
    println!(
        "  {:>10}  {:>12}  {:>12}  {:>8}  {:>8}  {:>12}",
        "Material", "Type", "|P_s| (C/m²)", "T_FE (K)", "T_N (K)", "|M_s| (A/m)"
    );
    for (name, mat) in &[("BiFeO3", &bfo), ("TbMnO3", &tbmno3), ("Cr2O3", &cr2o3)] {
        println!(
            "  {:>10}  {:>12}  {:>12.4}  {:>8.0}  {:>8.0}  {:>12.4e}",
            name,
            format!("{:?}", mat.material_type()),
            mat.polarization_s().magnitude(),
            mat.t_ferroelectric(),
            mat.t_magnetic(),
            mat.magnetization_s().magnitude(),
        );
    }
    println!();

    // -------------------------------------------------------------------------
    // 2. Linear magnetoelectric effect in Cr2O3 (canonical ME)
    // -------------------------------------------------------------------------
    println!("--- Linear Magnetoelectric Effect (Cr2O3) ---");
    println!("  Dzyaloshinskii bound:  |α|² ≤ ε₀·μ₀·χ_e·χ_m");
    let bound = cr2o3.bound_on_me_coupling();
    let susc = cr2o3.magnetoelectric_susceptibility();
    println!("  |α|_bound = {:.4e} s/m  (= C/(A·m))", bound);
    println!("  |α|       = {:.4e} s/m", susc);
    println!(
        "  α/α_bound = {:.4}  (must be < 1 by thermodynamics)",
        susc / bound
    );

    // Apply a magnetic field and compute the induced polarization
    let h_z = 1.0e6_f64; // 1 MA/m
    let h_vec = Vector3::new(0.0, 0.0, h_z);
    let p_induced = cr2o3.electric_polarization_from_field(&h_vec);
    let alpha_33 = cr2o3.alpha_tensor()[2][2]; // diagonal element
    println!(
        "\n  Apply H_z = {:.1e} A/m → P_z = {:.4e} C/m²",
        h_z, p_induced.z
    );
    println!("  (α_33·H_z = {:.4e} C/m² ✓)", alpha_33 * h_z);

    // Inverse ME: apply E and compute induced magnetisation
    let e_z = 1.0e6_f64; // 1 MV/m
    let e_vec = Vector3::new(0.0, 0.0, e_z);
    let m_induced = cr2o3.magnetization_from_efield(&e_vec);
    println!(
        "  Apply E_z = {:.1e} V/m → M_z = {:.4e} A/m",
        e_z, m_induced.z
    );
    println!();

    // -------------------------------------------------------------------------
    // 3. BiFeO3 — above/below transition temperatures
    // -------------------------------------------------------------------------
    println!("--- BiFeO3 Phase Diagram ---");
    let temps = [77.0, 300.0, 643.0, 700.0, 1100.0, 1200.0];
    println!(
        "  {:>8}  {:>10}  {:>14}",
        "T (K)", "AFM order?", "Ferroelectric?"
    );
    for t in &temps {
        let afm = !bfo.is_above_magnetic_transition(*t);
        let fe = !bfo.is_above_ferroelectric_transition(*t);
        println!("  {:>8.0}  {:>10}  {:>14}", t, afm, fe);
    }
    println!();

    // -------------------------------------------------------------------------
    // 4. DM-mechanism electric polarization (from spin canting in BFO)
    // -------------------------------------------------------------------------
    println!("--- Dzyaloshinskii-Moriya (KNB) Polarization in BFO ---");
    // Fe sublattice: slightly canted G-type AFM
    // Two adjacent spins with small canting angle θ ≈ 1°
    let theta_rad = 1.0_f64.to_radians();
    let s1 = Vector3::new(theta_rad.sin(), 0.0, theta_rad.cos()); // up + small x
    let s2 = Vector3::new(-theta_rad.sin(), 0.0, theta_rad.cos()); // down + small x
    let _e12 = Vector3::new(1.0, 0.0, 0.0).normalize(); // Fe-Fe bond along a

    // DM vector along b (y) for BFO symmetry
    let d_vec = Vector3::new(0.0, 1.0, 0.0);
    let knb_gamma = 1.0e-3; // representative coupling (C/m²)

    // dzyaloshinskii_moriya_polarization(s_i, s_j, bond_vector): no gamma arg — returns direction
    // We scale by knb_gamma to get SI polarization
    let p_dm_dir = dzyaloshinskii_moriya_polarization(&s1, &s2, &d_vec);
    let p_dm = Vector3::new(
        p_dm_dir.x * knb_gamma,
        p_dm_dir.y * knb_gamma,
        p_dm_dir.z * knb_gamma,
    );
    println!("  Spin canting angle θ = 1°");
    println!("  DM vector D ∥ b-axis");
    println!(
        "  P_DM = ({:>10.4e}, {:>10.4e}, {:>10.4e}) C/m²",
        p_dm.x, p_dm.y, p_dm.z
    );

    // Exchange striction: sensitivity = spin-phonon coupling [C/m²]
    let sensitivity = 5.0e-4_f64; // representative coupling (C/m²)
    let p_es_scalar = exchange_striction_polarization(&s1, &s2, sensitivity);
    println!("  P_exchange-striction = {:.4e} C/m²", p_es_scalar.abs());
    println!();

    // -------------------------------------------------------------------------
    // 5. Toroidal moment — characteristic of BFO's vertex sharing FeO6 octahedra
    // -------------------------------------------------------------------------
    println!("--- Toroidal Moment (Anapole Order) ---");
    // BFO is an anapole magnet: the magnetic order parameter transforms as a toroid
    // T = (1/2V) Σ r_i × m_i
    // Simplified 4-spin ring geometry (planar square)
    let spin_positions = vec![
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, -1.0, 0.0),
    ];
    let spin_moments = vec![
        Vector3::new(0.0, 0.5, 0.0),
        Vector3::new(-0.5, 0.0, 0.0),
        Vector3::new(0.0, -0.5, 0.0),
        Vector3::new(0.5, 0.0, 0.0),
    ];
    let toroid = toroidal_moment(&spin_positions, &spin_moments)?;
    println!("  4-spin vortex (ferroaxial arrangement):");
    println!(
        "  T = ({:>8.4}, {:>8.4}, {:>8.4}) A/m · Å",
        toroid.x, toroid.y, toroid.z
    );
    println!(
        "  |T| = {:.4e}  (non-zero → anapole order)",
        toroid.magnitude()
    );
    println!();

    // -------------------------------------------------------------------------
    // 6. Switching energy under combined E and H fields
    // -------------------------------------------------------------------------
    println!("--- Electric-Field Switching Energy ---");
    println!(
        "  {:>12}  {:>12}  {:>12}  {:>16}",
        "E (MV/m)", "H (kA/m)", "W_switch (eV)", "Description"
    );
    let cases: &[(f64, f64, &str)] = &[
        (1.0, 0.0, "Pure E-field"),
        (0.0, 800.0, "Pure H-field (1 T)"),
        (1.0, 800.0, "Combined E+H"),
        (5.0, 0.0, "High-field E"),
    ];
    for (e_mv, h_ka, label) in cases {
        let w = bfo.switching_energy(*e_mv * 1.0e6, *h_ka * 1.0e3);
        let w_ev = w / 1.602e-19;
        println!(
            "  {:>12.1}  {:>12.0}  {:>12.4e}  {:>16}",
            e_mv, h_ka, w_ev, label
        );
    }
    println!();

    // -------------------------------------------------------------------------
    // 7. Inverse magnetoelectric: E-field-induced magnetisation control
    // -------------------------------------------------------------------------
    println!("--- Inverse ME: Electric Field → Magnetisation Change (Cr2O3) ---");
    println!("  (Underlying principle for electric-field control of magnetic memories)");
    let alpha_me = cr2o3.magnetoelectric_susceptibility();
    let e_fields_mv = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0];
    println!("  {:>12}  {:>16}", "E (MV/m)", "ΔM_z (mA/m)");
    for &e_mv in &e_fields_mv {
        let ime = InverseMagnetoelectric::new(Vector3::new(0.0, 0.0, e_mv * 1.0e6), alpha_me)?;
        let dm = ime.induced_magnetization();
        println!("  {:>12.1}  {:>16.4}", e_mv, dm.z * 1.0e3);
    }

    println!("\n=================================================================");
    println!("  Summary: BiFeO3 Type-I multiferroic — P from polar distortion,");
    println!("  cross-coupling via linear ME + DM + exchange striction terms.");
    println!("=================================================================\n");

    Ok(())
}
