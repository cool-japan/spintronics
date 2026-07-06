//! Streaming Simulation API: Constant-Memory Long-Run Integration
//!
//! **Difficulty**: ⭐⭐ Intermediate
//! **Category**: Simulation Framework
//! **Physics**: LLG dynamics, Zeeman + anisotropy energy, magnetization relaxation
//!
//! `Simulation::run()` collects the *entire* trajectory into
//! `Vec<Vector3<f64>>` / `Vec<f64>` before returning, which is convenient for
//! short runs but becomes memory-prohibitive for very long integrations (or
//! for sweeping many long integrations back to back). `Simulation::run_streaming()`
//! performs the identical integration loop -- same solver dispatch, same
//! finite-value checks, same renormalization -- but instead of collecting a
//! trajectory it invokes a caller-supplied callback once per completed step,
//! `callback(step_index, &magnetization, energy)`, and never retains a
//! growing collection internally.
//!
//! We demonstrate:
//! 1. A long run (500,000 steps -- 500x the typical `num_steps` default of
//!    1000) tracked with a handful of running statistics instead of a
//!    multi-megabyte trajectory buffer.
//! 2. Early termination: returning `Err` from the callback aborts the
//!    integration immediately, which is useful for "stop once converged"
//!    style logic where the required step count isn't known in advance.
//! 3. A back-of-the-envelope comparison of the memory a collected `run()`
//!    trajectory would need at this scale versus the streaming path's
//!    constant footprint.
//!
//! Reference: Gilbert, IEEE Trans. Magn. 40, 3443 (2004)

use spintronics::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== Streaming Simulation API ===\n");

    // ------------------------------------------------------------------
    // 1. Long run tracked with O(1) running statistics (no growing Vec).
    // ------------------------------------------------------------------
    println!("=== Part 1: Long run, constant-memory statistics ===");

    let num_steps = 500_000usize;
    let mut sim = SimulationBuilder::new()
        .material(Ferromagnet::yig())
        .external_field(Vector3::new(0.0, 0.0, 0.1))
        .solver_rk4()
        .initial_magnetization(Vector3::new(1.0, 0.0, 0.0))
        .time_step(1.0e-13)
        .num_steps(num_steps)
        .build()?;

    println!("  Material: YIG, solver: RK4, num_steps = {}", num_steps);
    println!("  (the SimulationBuilder example above uses num_steps = 1000)\n");

    // Running statistics kept in a handful of stack scalars -- this is the
    // *entire* memory footprint of the streaming consumer below, regardless
    // of whether num_steps is 1000 or 1,000,000,000. Contrast with `run()`,
    // whose `SimulationResult` would hold `num_steps + 1` full Vector3/f64
    // entries.
    let mut steps_seen = 0usize;
    let mut energy_min = f64::INFINITY;
    let mut energy_max = f64::NEG_INFINITY;
    let mut energy_sum = 0.0f64;
    let progress_stride = (num_steps / 5).max(1);

    sim.run_streaming(|step_index, m, energy| {
        steps_seen += 1;
        energy_min = energy_min.min(energy);
        energy_max = energy_max.max(energy);
        energy_sum += energy;

        if step_index % progress_stride == 0 {
            println!(
                "  step {:>7}/{}: m = ({:+.6}, {:+.6}, {:+.6}), |m| = {:.8}, E = {:.6e} J/m^3",
                step_index,
                num_steps,
                m.x,
                m.y,
                m.z,
                m.magnitude(),
                energy
            );
        }
        Ok(())
    })?;

    let energy_mean = energy_sum / steps_seen as f64;
    println!(
        "\n  Completed {} reported states (initial state + {} integration steps)",
        steps_seen, num_steps
    );
    println!(
        "  Energy density over the run: min = {:.6e}, max = {:.6e}, mean = {:.6e} J/m^3",
        energy_min, energy_max, energy_mean
    );

    // ------------------------------------------------------------------
    // 2. Early termination via callback error propagation.
    // ------------------------------------------------------------------
    println!("\n=== Part 2: Early stop once magnetization aligns with the field ===");

    let step_budget = 1_000_000usize;
    let mut sim2 = SimulationBuilder::new()
        .material(Ferromagnet::permalloy())
        .external_field(Vector3::new(0.0, 0.0, 0.2))
        .solver_rk4()
        .damping(0.1) // higher damping so relaxation completes within the step budget
        .initial_magnetization(Vector3::new(0.2, 0.0, 0.98))
        .time_step(1.0e-13)
        .num_steps(step_budget)
        .build()?;

    let alignment_threshold = 0.999_9;
    let mut converged_at: Option<usize> = None;
    let outcome = sim2.run_streaming(|step_index, m, _energy| {
        if m.z.abs() >= alignment_threshold {
            converged_at = Some(step_index);
            return Err(Error::ConfigurationError {
                description: format!(
                    "requested stop: |m_z| = {:.8} reached the alignment threshold",
                    m.z.abs()
                ),
            });
        }
        Ok(())
    });

    match (converged_at, outcome) {
        (Some(step), Err(reason)) => {
            println!(
                "  Converged after {} of a {}-step budget -- stopped early: {}",
                step, step_budget, reason
            );
        },
        (None, Ok(())) => {
            println!(
                "  Did not reach |m_z| >= {} within the {}-step budget",
                alignment_threshold, step_budget
            );
        },
        (None, Err(reason)) => {
            // A genuine failure inside the integrator, not our intentional stop.
            return Err(Box::new(reason));
        },
        (Some(_), Ok(())) => unreachable!("callback only returns Err once converged_at is set"),
    }

    // ------------------------------------------------------------------
    // 3. Memory footprint comparison.
    // ------------------------------------------------------------------
    println!("\n=== Part 3: Why streaming matters at this scale ===");
    let bytes_per_state = std::mem::size_of::<Vector3<f64>>() + std::mem::size_of::<f64>();
    let collected_bytes = bytes_per_state * (num_steps + 1);
    let streaming_bytes = std::mem::size_of::<usize>() + 3 * std::mem::size_of::<f64>();
    println!(
        "  A collected run() trajectory at num_steps = {} would retain ~{:.2} MB \
        (trajectory + energies Vecs).",
        num_steps,
        collected_bytes as f64 / (1024.0 * 1024.0)
    );
    println!(
        "  run_streaming() processed the same {} steps using {} bytes of running \
        state -- independent of num_steps.",
        num_steps, streaming_bytes
    );

    println!("\n=== Summary ===");
    println!("  - run_streaming(callback) mirrors run()'s integration loop step-for-step");
    println!("  - callback(step_index, &magnetization, energy) is invoked once per step");
    println!("  - no Vec<Vector3<f64>> / Vec<f64> is ever accumulated internally");
    println!("  - returning Err from the callback aborts the integration immediately");

    Ok(())
}
