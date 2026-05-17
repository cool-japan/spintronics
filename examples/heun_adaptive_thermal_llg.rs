//! Adaptive Heun integration of the stochastic Landau-Lifshitz-Gilbert
//! equation with frozen-noise step-size control.
//!
//! **Difficulty**: ⭐⭐⭐
//! **Category**: Stochastic Numerical Methods / Thermal LLG
//! **Physics**: Embedded Euler-Heun pair with PI-style adaptive dt;
//!              thermal noise follows the fluctuation-dissipation theorem.
//!
//! ## What this demo shows
//!
//! 1. Build a single-spin SLLG problem with a Tesla-scale external field
//!    along z and Gilbert damping `α = 0.05`.
//! 2. Run the adaptive Heun stepper at room temperature (T = 300 K) for a
//!    few hundred steps.
//! 3. Track the suggested step size, the embedded error estimate and the
//!    z-component of the magnetization to illustrate the controller
//!    behaviour.
//!
//! ## References
//! - Klauder & Petersen, J. Stat. Phys. 39, 53 (1985)
//! - García-Palacios & Lázaro, PRB 58, 14937 (1998)
//! - Hairer & Wanner, "Solving ODEs II" Ch. IV.8 (1996) — PI step control

#[cfg(feature = "scirs2")]
fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    use spintronics::stochastic::{HeunAdaptive, ThermalField};
    use spintronics::vector3::Vector3;

    println!("=============================================================");
    println!("  Adaptive Heun on Stochastic LLG (single spin, T = 300 K)");
    println!("=============================================================");

    let alpha = 0.05_f64;
    let temperature = 300.0_f64;
    let volume = 1.0e-24_f64; // ~1 nm³ cell
    let ms = 1.0e6_f64; // A/m
    let dt0 = 1.0e-13_f64; // initial guess: 0.1 ps
    let h_eff = Vector3::new(0.0, 0.0, 1.0); // T-scale external field

    let mut integ = HeunAdaptive::new(dt0, alpha)?
        .with_tolerances(1.0e-6, 1.0e-4)
        .with_dt_bounds(1.0e-18, 1.0e-10);
    let mut thermal = ThermalField::new(temperature, volume, ms, alpha);

    let mut m = Vector3::new(1.0, 0.0, 0.0);
    let mut t = 0.0_f64;
    let n_steps = 200_usize;

    println!("\n step     t (ps)     dt (fs)    |m|       m_z        err\n");
    for step in 0..n_steps {
        let (m_new, dt_used, err) = integ.step(m, |_| h_eff, &mut thermal)?;
        m = m_new;
        t += dt_used;
        if step % 20 == 0 {
            println!(
                " {step:4}   {:8.3}   {:8.3}   {:.4}   {:+.4}   {:.2e}",
                t * 1.0e12,
                dt_used * 1.0e15,
                m.magnitude(),
                m.z,
                err
            );
        }
    }
    println!(
        "\nFinal: t = {:.3} ps, dt_current = {:.3} fs, |m| = {:.6}",
        t * 1.0e12,
        integ.current_dt * 1.0e15,
        m.magnitude()
    );
    Ok(())
}

#[cfg(not(feature = "scirs2"))]
fn main() {
    println!("Build with `--features scirs2` to run this stochastic example.");
}
