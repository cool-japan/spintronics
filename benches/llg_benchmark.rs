use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use spintronics::dynamics::{anisotropy_energy, exchange_energy};
use spintronics::prelude::*;
use spintronics::simd::{batch_evolve_rk4, SimdBatch};

fn bench_calc_dm_dt(c: &mut Criterion) {
    let m = Vector3::new(1.0, 0.0, 0.0);
    let h_eff = Vector3::new(0.0, 0.0, 1.0);
    let alpha = 0.01;

    c.bench_function("calc_dm_dt_single", |b| {
        b.iter(|| black_box(calc_dm_dt(black_box(m), black_box(h_eff), GAMMA, alpha)))
    });
}

fn bench_llg_rk4_step(c: &mut Criterion) {
    let solver = LlgSolver::new(0.01, 1.0e-13);
    let m0 = Vector3::new(0.1, 0.0, 1.0).normalize();
    let h_field = Vector3::new(0.0, 0.0, 1.0);

    c.bench_function("llg_rk4_single_spin", |b| {
        b.iter(|| {
            let h_eff_fn = |_m: Vector3<f64>| h_field;
            black_box(solver.step_rk4(black_box(m0), h_eff_fn))
        })
    });
}

fn bench_llg_evolution_100_spins(c: &mut Criterion) {
    let solver = LlgSolver::new(0.01, 1.0e-13);
    let h_field = Vector3::new(0.0, 0.0, 1.0);

    // Create 100 spins with slight tilts
    let spins: Vec<Vector3<f64>> = (0..100)
        .map(|i| {
            let angle = (i as f64) * 0.01;
            Vector3::new(angle.sin(), 0.0, angle.cos()).normalize()
        })
        .collect();

    c.bench_function("llg_evolution_100spins_10steps", |b| {
        b.iter(|| {
            let mut current_spins = spins.clone();
            for _ in 0..10 {
                for spin in &mut current_spins {
                    let h_eff_fn = |_m: Vector3<f64>| h_field;
                    *spin = solver.step_rk4(*spin, h_eff_fn);
                }
            }
            black_box(&current_spins);
        })
    });
}

fn bench_exchange_energy(c: &mut Criterion) {
    let a_ex = 1.3e-11; // Permalloy exchange stiffness
    let angle = 0.01; // Small rotation angle
    let volume = 1e-24; // 10nm cube
    let length_scale = 1e-8; // 10nm

    c.bench_function("exchange_energy", |b| {
        b.iter(|| {
            black_box(exchange_energy(
                black_box(a_ex),
                black_box(angle),
                black_box(volume),
                black_box(length_scale),
            ))
        })
    });
}

fn bench_anisotropy_energy(c: &mut Criterion) {
    let m = Vector3::new(0.6, 0.0, 0.8);
    let easy_axis = Vector3::new(0.0, 0.0, 1.0);
    let k = 5.0e5; // Uniaxial anisotropy constant

    c.bench_function("anisotropy_energy", |b| {
        b.iter(|| {
            black_box(anisotropy_energy(
                black_box(m),
                black_box(easy_axis),
                black_box(k),
            ))
        })
    });
}

/// Benchmark group "rk4_N1024": scalar loop vs SIMD batch for N=1024 independent spins.
///
/// - `scalar_rk4_N1024`: evolves 1024 spins one at a time via `LlgSolver::step_rk4`.
/// - `simd_rk4_N1024`: evolves 1024 spins simultaneously via `batch_evolve_rk4`.
fn bench_rk4_n1024(c: &mut Criterion) {
    const N: usize = 1024;
    let alpha = 0.01_f64;
    let dt = 1.0e-13_f64;
    let h_field = Vector3::new(0.0, 0.0, 1.0);

    // Pre-build spin arrays
    let spins: Vec<Vector3<f64>> = (0..N)
        .map(|i| {
            let angle = (i as f64) * std::f64::consts::TAU / N as f64;
            Vector3::new(angle.cos(), angle.sin(), 0.0)
        })
        .collect();

    let mut group = c.benchmark_group("rk4_N1024");

    // Scalar: evolve each spin independently using LlgSolver::step_rk4
    group.bench_with_input(
        BenchmarkId::new("scalar_rk4_N1024", N),
        &spins,
        |b, spins| {
            let solver = LlgSolver::new(alpha, dt);
            b.iter(|| {
                let h_fn = |_m: Vector3<f64>| h_field;
                let result: Vec<Vector3<f64>> = spins
                    .iter()
                    .map(|&m| solver.step_rk4(black_box(m), h_fn))
                    .collect();
                black_box(result)
            });
        },
    );

    // SIMD batch: evolve all 1024 spins in one call
    let m_batch = SimdBatch::from_vector3_slice(&spins);
    let mut h_batch = SimdBatch::new(N);
    for i in 0..N {
        h_batch.z[i] = 1.0;
    }

    group.bench_with_input(
        BenchmarkId::new("simd_rk4_N1024", N),
        &(m_batch, h_batch),
        |b, (m_batch, h_batch)| {
            b.iter(|| {
                // Reconstruct from the stored batch each iteration to avoid
                // mutating the benchmark input (batch_evolve_rk4 takes &SimdBatch).
                black_box(batch_evolve_rk4(
                    black_box(m_batch),
                    black_box(h_batch),
                    alpha,
                    GAMMA,
                    dt,
                ))
            });
        },
    );

    group.finish();
}

criterion_group!(
    benches,
    bench_calc_dm_dt,
    bench_llg_rk4_step,
    bench_llg_evolution_100_spins,
    bench_exchange_energy,
    bench_anisotropy_energy,
    bench_rk4_n1024
);
criterion_main!(benches);
