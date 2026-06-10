//! Path-integral Monte Carlo (PIMC) for a finite-temperature Heisenberg
//! ferromagnetic chain in an applied field.
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Stochastic Methods / Quantum Magnetism
//! **Physics**: Trotterized worldline updates on a classical Heisenberg
//!              chain. The applied field aligns the spins along z; we
//!              measure `⟨M_z⟩`, the per-site energy and the
//!              susceptibility χ as a function of inverse temperature β.
//!
//! ## What this demo shows
//!
//! 1. Build a 1D ferromagnetic chain of `n_spins` sites with periodic
//!    boundary conditions (`Ring1D`).
//! 2. Discretize the imaginary time direction into `n_slices = 4` Trotter
//!    slices.
//! 3. Run PIMC at three temperatures (high T, intermediate T, low T) and
//!    show that the magnetisation grows from ~0 to a polarised state as
//!    β increases.
//!
//! ## References
//! - Suzuki, Prog. Theor. Phys. 56, 1454 (1976)
//! - Sandvik & Kurkijärvi, PRB 43, 5950 (1991)
//! - Bonner & Fisher, PR 135, A640 (1964) — exact Heisenberg chain

#[cfg(feature = "scirs2")]
fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    use spintronics::stochastic::{PimcConfig, PimcLattice, PimcSimulation};

    println!("=============================================================");
    println!("  PIMC: 1D Heisenberg Ring (J = 1, h_z = 0.5) vs β");
    println!("=============================================================");

    let n_spins = 8;
    let j = 1.0_f64;
    let h_z = 0.5_f64;

    println!("\n   β        ⟨M_z⟩        E/N          χ        accept\n");
    for &beta in &[0.05_f64, 0.5, 2.0, 5.0] {
        let lattice = PimcLattice::Ring1D { n_spins, j, h_z };
        let config = PimcConfig {
            n_slices: 4,
            n_steps: 1500,
            n_thermalize: 500,
            beta,
            seed: 2024,
            angular_step: 0.5,
        };
        let mut sim = PimcSimulation::new(lattice, config)?;
        sim.thermalize()?;
        let res = sim.run()?;
        let acc = res.n_accepted as f64 / res.n_proposed.max(1) as f64;
        println!(
            " {beta:5.2}   {:+.4}      {:+.4}     {:.3}     {:.3}",
            res.average_magnetization_z, res.average_energy_per_site, res.susceptibility, acc
        );
    }
    Ok(())
}

#[cfg(not(feature = "scirs2"))]
fn main() {
    println!("Build with `--features scirs2` to run this PIMC example.");
}
