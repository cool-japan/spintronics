//! TbMnO3 Spin Spiral and Electric Polarization via the KNB Mechanism
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Non-Collinear Magnetism / Multiferroics
//! **Physics**: Cycloidal spin spirals in frustrated magnets; Katsura-Nagaosa-Balatsky
//! mechanism for spiral-driven electric polarization; Luttinger-Tisza ground-state search.
//!
//! ## Background
//!
//! TbMnO₃ is a canonical Type-II multiferroic: its electric polarization arises entirely
//! from a cycloidal spin spiral stabilised by competing nearest-neighbour (J₁ < 0) and
//! next-nearest-neighbour (J₂ > 0) exchange on the orthorhombic Mn lattice.
//!
//! The Luttinger-Tisza method locates the ordering wavevector q* by minimising J(q) over
//! the Brillouin zone.  The KNB formula then gives the polarisation direction:
//!
//! ```text
//! P ∝ γ × (e_ij × (S_i × S_j))
//! ```
//!
//! where γ is the crystal-axis coupling constant and e_ij is the Mn–Mn bond vector.
//!
//! ## References
//!
//! - T. Kimura et al., "Magnetic control of ferroelectric polarization",
//!   *Nature* **426**, 55 (2003).
//! - H. Katsura, N. Nagaosa & A. V. Balatsky, "Spin Current and Magnetoelectric Effect
//!   in Noncollinear Magnets", *Phys. Rev. Lett.* **95**, 057205 (2005).
//! - J. M. Luttinger & L. Tisza, "Theory of Dipole Interaction in Crystals",
//!   *Phys. Rev.* **70**, 954 (1946).
//! - M. Mostovoy, "Ferroelectricity in Spiral Magnets",
//!   *Phys. Rev. Lett.* **96**, 067601 (2006).

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  TbMnO3: Spin Spiral Ground State & KNB Electric Polarization");
    println!("=================================================================\n");

    // -------------------------------------------------------------------------
    // 1. Luttinger-Tisza ground-state search for TbMnO3 Mn lattice
    // -------------------------------------------------------------------------
    // Orthorhombic Mn–Mn bonds: a = 5.29 Å, J1 = -2.0 meV (AFM NN),
    // J2 = +0.5 meV (FM NNN) — ratio J2/|J1| = 0.25 drives the spiral.
    let a_mn = 5.29e-10_f64; // m
    let j1 = -2.0e-3 * 1.602e-19; // -2 meV in J
    let j2 = 0.5e-3 * 1.602e-19; //  0.5 meV in J

    // Build J1-J2 frustrated chain model (extends naturally to 3D)
    let lt = LuttingerTisza::j1j2_chain(j1, j2, a_mn, 2.0);

    let q_opt = lt.find_ground_state_q(200);
    let is_spiral = lt.is_spiral_ground_state(1.0e-6);
    let frustration = lt.frustration_ratio();
    let e_per_spin = lt.ground_state_energy_per_spin(&q_opt);
    let pitch_deg = lt.spiral_pitch_angle(&q_opt).to_degrees();
    let t_ord = lt.ordering_temperature_estimate();

    println!("--- Luttinger-Tisza Ground State (J1-J2 Mn chain) ---");
    println!("  J1 = {:>+.2} meV  (AFM NN)", j1 / 1.602e-22);
    println!("  J2 = {:>+.2} meV  (FM NNN)", j2 / 1.602e-22);
    println!(
        "  Frustration ratio J2/|J1| = {:.3} (> 0.25 → spiral phase)",
        frustration
    );
    println!(
        "  Ground-state q* = ({:.4}, {:.4}, {:.4}) × π/a",
        q_opt.x * a_mn / std::f64::consts::PI,
        q_opt.y * a_mn / std::f64::consts::PI,
        q_opt.z * a_mn / std::f64::consts::PI,
    );
    println!("  Energy per spin  = {:.4} meV", e_per_spin / 1.602e-22);
    println!("  Spiral pitch     = {:.1}°", pitch_deg);
    println!("  Is spiral GS?    = {}", is_spiral);
    println!("  T_ord estimate   = {:.1} K\n", t_ord);

    // -------------------------------------------------------------------------
    // 2. Build the TbMnO3 spin spiral using the experimental parameters
    // -------------------------------------------------------------------------
    // Experimental: q ≈ 0.28 × (2π/b) b̂, cycloidal in (a,c) plane
    let spiral = SpinSpiral::terbium_manganese_oxide();
    let chirality = spiral.chirality();
    let wavelength_nm = spiral.wavelength() * 1.0e9;
    let pitch_ang = spiral.cone_angle_deg();

    println!("--- TbMnO3 SpinSpiral Object ---");
    println!("  Spiral type = {:?}", spiral.spiral_type);
    println!("  Chirality   = {:?}", chirality);
    println!("  Wavelength  = {:.2} nm", wavelength_nm);
    println!("  Cone angle  = {:.1}°", pitch_ang);
    println!("  Incommen.?  = {}", spiral.is_incommensurate(a_mn));

    // TbMnO3: spiral propagates along b-axis (ŷ). Sample positions along b.
    let b_lattice = 5.84e-10_f64;
    println!("\n  Magnetization profile along b (y = n·b):");
    println!("  {:>8}  {:>10}  {:>10}  {:>10}", "y (Å)", "Mx", "My", "Mz");
    for n in 0..8 {
        let r = Vector3::new(0.0, n as f64 * b_lattice, 0.0);
        let m = spiral.magnetization_at(&r);
        println!(
            "  {:>8.2}  {:>10.4}  {:>10.4}  {:>10.4}",
            r.y * 1e10,
            m.x,
            m.y,
            m.z
        );
    }

    // -------------------------------------------------------------------------
    // 3. KNB electric polarization from the cycloidal spiral
    // -------------------------------------------------------------------------
    let knb = KnbMechanism::tbmno3();
    let coupling = knb.coupling();

    // TbMnO3 spiral propagates along b-axis (y), so sample nearest Mn neighbors along b.
    // Bond vector e12 along b-axis (y): the Mn-Mn NN bond in the bc-plane.
    let b_mn = 5.84e-10_f64; // b lattice parameter
    let e12 = Vector3::new(0.0, 1.0, 0.0).normalize();
    // Neighbouring spin pair along b: Δr = b (one lattice spacing)
    let r0 = Vector3::new(0.0, 0.0, 0.0);
    let r1 = Vector3::new(0.0, b_mn, 0.0);
    let s1 = spiral.magnetization_at(&r0);
    let s2 = spiral.magnetization_at(&r1);

    let p_knb = knb.polarization_from_spiral(&s1, &s2, &e12);
    let p_mag = p_knb.magnitude();

    println!("\n--- KNB Mechanism Electric Polarization ---");
    println!("  Coupling γ        = {:.4e} C/m²", coupling);
    println!(
        "  S1 (at x=0)       = ({:>8.4}, {:>8.4}, {:>8.4})",
        s1.x, s1.y, s1.z
    );
    println!(
        "  S2 (at x=a)       = ({:>8.4}, {:>8.4}, {:>8.4})",
        s2.x, s2.y, s2.z
    );
    println!(
        "  P direction       = ({:>8.4}, {:>8.4}, {:>8.4})",
        p_knb.x, p_knb.y, p_knb.z
    );
    println!("  |P|               = {:.4e} C/m²", p_mag);
    println!("  |P| (μC/cm²)      = {:.4}", p_mag * 1.0e4);

    // Verify P lies along c-axis (z in our coordinate system)
    let p_polar_frac = if p_mag > 1.0e-30 {
        p_knb.z.abs() / p_mag
    } else {
        0.0
    };
    println!("  P along c-axis?   = {:.1}%", p_polar_frac * 100.0);

    // -------------------------------------------------------------------------
    // 4. Spiral-type comparison: cycloidal vs helical
    // -------------------------------------------------------------------------
    let b_lattice_local = 5.84e-10_f64;
    let q_cycl = Vector3::new(
        0.0,
        2.0 * std::f64::consts::PI * 0.28 / b_lattice_local,
        0.0,
    );
    let spiral_cycl = SpinSpiral::cycloidal(q_cycl, Vector3::new(1.0, 0.0, 0.0), 1.0);
    let spiral_heli = SpinSpiral::helical(q_cycl, 1.0);

    let p_dir_cycl = spiral_cycl.electric_polarization_direction();
    let p_dir_heli = spiral_heli.electric_polarization_direction();

    println!("\n--- Polarization Direction by Spiral Type ---");
    println!(
        "  Cycloidal spiral → P = {}",
        if let Some(p) = p_dir_cycl {
            format!("({:.3}, {:.3}, {:.3})", p.x, p.y, p.z)
        } else {
            "None (forbidden by symmetry)".to_string()
        }
    );
    println!(
        "  Helical spiral   → P = {}",
        if let Some(p) = p_dir_heli {
            format!("({:.3}, {:.3}, {:.3})", p.x, p.y, p.z)
        } else {
            "None (forbidden by symmetry)".to_string()
        }
    );
    println!("  (Helical spirals have P = 0: the DM vector averages to zero)");

    // -------------------------------------------------------------------------
    // 5. Spin structure factor (neutron scattering signal)
    // -------------------------------------------------------------------------
    // Use spiral's own q_vector (along b) to probe the magnetic Bragg peaks
    let q_spiral = spiral.q_vector;
    println!("\n--- Spin Structure Factor S(k) at High-Symmetry k-points ---");
    println!(
        "  {:>20}  {:>12}  {:>12}",
        "k-point", "S+q (wt.)", "S-q (wt.)"
    );
    let q_spiral_neg = q_spiral * (-1.0);
    let q_spiral_dbl = q_spiral * 2.0;
    let k_pts: &[(&str, Vector3<f64>)] = &[
        ("q_TbMnO3", q_spiral),
        ("-q_TbMnO3", q_spiral_neg),
        ("2q_TbMnO3", q_spiral_dbl),
        ("0 (FM peak)", Vector3::new(0.0, 0.0, 0.0)),
    ];
    for (label, k) in k_pts {
        let (w_plus, w_minus) = spiral.spin_structure_factor(k);
        println!("  {:>20}  {:>12.4e}  {:>12.4e}", label, w_plus, w_minus);
    }

    // -------------------------------------------------------------------------
    // 6. Landau-Lifshitz exchange energy
    // -------------------------------------------------------------------------
    let j_nn = -1.5e-3 * 1.602e-19;
    let e_ll = spiral.landau_lifshitz_energy(j_nn, b_lattice);
    println!(
        "\n--- Exchange Energy (Landau-Lifshitz) ---\n  E/spin = {:.4e} meV",
        e_ll / 1.602e-22
    );

    println!("\n=================================================================");
    println!("  Summary: TbMnO3 cycloidal spiral → finite P along c-axis");
    println!("  This demonstrates the KNB mechanism for Type-II multiferroics.");
    println!("=================================================================\n");

    Ok(())
}
