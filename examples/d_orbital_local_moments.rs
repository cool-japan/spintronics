//! Local d-Orbital Magnetic Moments: Crystal-Field Theory for 3d Ions
//!
//! **Difficulty**: ⭐⭐⭐ Advanced
//! **Category**: Orbitronics / Local-Moment Magnetism
//! **Physics**: Crystal-field splitting, Hund's rules, high-spin/low-spin,
//! orbital quenching, spin-orbit-driven unquenching, effective moments
//!
//! This example demonstrates single-configuration ligand-field theory for 3d
//! transition-metal ions: the free-ion Hund's-rule ground term, the
//! octahedral/tetrahedral/tetragonal crystal-field splitting of the five real
//! d-orbitals, high-spin vs low-spin competition (10Dq vs pairing energy),
//! orbital angular-momentum quenching (and its lifting by spin-orbit
//! coupling for orbitally-degenerate T-term ions), and effective magnetic
//! moments compared against textbook values.
//!
//! Fidelity boundary: single-configuration ligand-field theory + atomic
//! Hund's-rule terms + perturbative single-particle spin-orbit coupling —
//! **not** full many-electron multiplet/Tanabe-Sugano theory (see the
//! `crystal_field`/`d_orbital_moment` module docs for details).
//!
//! References:
//! - C. J. Ballhausen, "Introduction to Ligand Field Theory" (McGraw-Hill, 1962)
//! - J. S. Griffith, "The Theory of Transition-Metal Ions" (Cambridge, 1961)

use spintronics::orbitronics::crystal_field::{CrystalFieldEnvironment, DOrbital};
use spintronics::orbitronics::d_orbital_moment::{free_ion_term, GroundTermSymmetry};
use spintronics::orbitronics::CrystalFieldModel;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Local d-Orbital Magnetic Moments: Ligand-Field Theory ===\n");

    // 1. Free-ion Hund's-rule ground terms, d1-d9.
    println!("=== Free-Ion Hund's-Rule Ground Terms (d1-d9) ===");
    println!("{:<6} {:>6} {:>6} {:>6}", "d^n", "S", "L", "J");
    println!("{}", "-".repeat(28));
    for n in 1..=9u8 {
        let term = free_ion_term(n)?;
        println!("d{:<5} {:>6.1} {:>6.0} {:>6.1}", n, term.s, term.l, term.j);
    }
    println!("\n  Electron-hole symmetry: L(n) = L(10-n), e.g. d2 and d8 both have L=3.");
    println!("  Note d4 (Cr2+/Mn3+ high-spin) has J=|L-S|=0: the free-ion moment");
    println!("  formally vanishes at this order (5D0 term).\n");

    // 2. Crystal-field splitting patterns.
    println!("=== Crystal-Field Splitting (10Dq = 2.0 eV) ===");
    let ten_dq = 2.0;
    for (label, env) in [
        ("Octahedral", CrystalFieldEnvironment::Octahedral),
        ("Tetrahedral", CrystalFieldEnvironment::Tetrahedral),
        (
            "Tetragonal (Ds=0.3, Dt=0.05)",
            CrystalFieldEnvironment::TetragonalDistorted { ds: 0.3, dt: 0.05 },
        ),
    ] {
        let h = spintronics::orbitronics::crystal_field_hamiltonian(env, ten_dq)?;
        print!("  {:<30}", label);
        for orbital in DOrbital::all() {
            print!(
                " {}={:+.2}",
                orbital.label(),
                h.get(orbital.index(), orbital.index()).re
            );
        }
        println!();
    }
    println!();

    // 3. Preset ion table: S, ground-term symmetry, quenching, moments.
    println!("=== Preset Ion Table ===");
    println!(
        "  {:<14} {:>4} {:>7} {:>6} {:>10} {:>9} {:>9}",
        "Ion", "d^n", "S", "Term", "Quenched?", "mu_so", "mu_Lande"
    );
    println!("  {}", "-".repeat(70));
    let presets: Vec<(&str, CrystalFieldModel)> = vec![
        ("Ti3+", CrystalFieldModel::ti3_plus()?),
        ("V3+", CrystalFieldModel::v3_plus()?),
        ("Cr3+", CrystalFieldModel::cr3_plus()?),
        ("Mn3+ (HS)", CrystalFieldModel::mn3_plus_high_spin()?),
        ("Fe3+", CrystalFieldModel::fe3_plus()?),
        ("Mn2+", CrystalFieldModel::mn2_plus()?),
        ("Fe2+ (HS)", CrystalFieldModel::fe2_plus_high_spin()?),
        ("Fe2+ (LS)", CrystalFieldModel::fe2_plus_low_spin()?),
        ("Co2+", CrystalFieldModel::co2_plus()?),
        ("Ni2+", CrystalFieldModel::ni2_plus()?),
        ("Cu2+", CrystalFieldModel::cu2_plus()?),
        ("Co3+ (LS)", CrystalFieldModel::co3_plus_low_spin()?),
    ];
    for (name, model) in &presets {
        let s = model.ground_state_spin()?;
        let symmetry = model.ground_term_symmetry()?;
        let symmetry_label = match symmetry {
            GroundTermSymmetry::A => "A",
            GroundTermSymmetry::E => "E",
            GroundTermSymmetry::T => "T",
        };
        let quenched = model.is_orbitally_quenched()?;
        let mu_so = model.effective_moment_spin_only()?;
        let mu_lande = model
            .effective_moment_lande()
            .map(|v| format!("{:.2}", v))
            .unwrap_or_else(|_| "undef.".to_string());
        println!(
            "  {:<14} d{:<3} {:>7.1} {:>6} {:>10} {:>9.2} {:>9}",
            name, model.n, s, symmetry_label, quenched, mu_so, mu_lande
        );
    }

    println!("\n  Key: Cr3+ (A2g) and Ni2+ (A2g) match the spin-only prediction almost");
    println!("  exactly in real complexes; Ti3+/Co2+/Fe2+(HS) (T-term) show sizeable");
    println!("  *positive* deviations in practice from unquenched orbital contributions.");
    println!("  Fe2+/Co3+ low-spin are diamagnetic (S=0) despite sharing d6 free-ion");
    println!("  parentage with high-spin Fe2+ (S=2) -- pure crystal-field pairing effect.\n");

    // 4. Orbital quenching lifted by spin-orbit coupling.
    println!("=== Orbital Quenching vs Spin-Orbit-Driven Unquenching ===");
    println!("  (RMS <Lz> of the ground-configuration's orbitally-degenerate \"frontier\"");
    println!("   level under H_SOC = lambda*L.S; exactly 0 when no such level exists.)\n");
    let ti = CrystalFieldModel::ti3_plus()?;
    let cr = CrystalFieldModel::cr3_plus()?;
    let cu = CrystalFieldModel::cu2_plus()?;
    println!(
        "  Ti3+ (T2g, unquenched): RMS <Lz> = {:.4}",
        ti.ground_manifold_orbital_moment_rms()?
    );
    println!(
        "  Cr3+ (A2g, quenched):   RMS <Lz> = {:.4}  (no frontier level: exactly 0)",
        cr.ground_manifold_orbital_moment_rms()?
    );
    println!(
        "  Cu2+ (Eg, quenched):    RMS <Lz> = {:.4}  (E has zero L matrix elements: exactly 0)\n",
        cu.ground_manifold_orbital_moment_rms()?
    );
    println!("  SOC lifts orbital moment only when the crystal-field ground level is");
    println!("  orbitally degenerate (T); A and E ground terms stay exactly quenched at");
    println!("  this single-configuration level of theory.\n");

    // 5. High-spin -> low-spin crossover as 10Dq grows for a fixed d-count.
    println!("=== High-Spin/Low-Spin Crossover (d6, pairing energy = 2.6 eV) ===");
    println!("  {:>10} {:>12} {:>8}", "10Dq (eV)", "High-spin?", "S");
    for ten_dq in [1.0, 1.8, 2.2, 2.6, 3.0, 3.8, 6.0] {
        let model =
            CrystalFieldModel::new(6, CrystalFieldEnvironment::Octahedral, ten_dq, 2.6, 0.05)?;
        println!(
            "  {:>10.1} {:>12} {:>8.1}",
            ten_dq,
            model.is_high_spin()?,
            model.ground_state_spin()?
        );
    }
    println!(
        "\n  Crossover occurs near 10Dq ~ pairing energy, exactly as textbook theory predicts."
    );

    Ok(())
}
