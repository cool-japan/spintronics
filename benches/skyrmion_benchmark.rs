use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use spintronics::builder::SimulationBuilder;

fn bench_skyrmion_dynamics_build(c: &mut Criterion) {
    c.bench_function("skyrmion_dynamics_build", |b| {
        b.iter(|| {
            black_box(
                SimulationBuilder::skyrmion_dynamics()
                    .build()
                    .expect("skyrmion_dynamics preset should build"),
            )
        })
    });
}

fn bench_skyrmion_dynamics_run_50steps(c: &mut Criterion) {
    c.bench_function("skyrmion_dynamics_run_50steps", |b| {
        b.iter(|| {
            let mut sim = SimulationBuilder::skyrmion_dynamics()
                .num_steps(50)
                .build()
                .expect("skyrmion_dynamics preset should build");
            let result = sim
                .run()
                .expect("skyrmion_dynamics preset should run to completion");
            black_box(result.final_magnetization);
        })
    });
}

criterion_group!(
    benches,
    bench_skyrmion_dynamics_build,
    bench_skyrmion_dynamics_run_50steps
);
criterion_main!(benches);
