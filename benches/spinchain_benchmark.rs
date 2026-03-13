use criterion::{black_box, criterion_group, criterion_main, Criterion};
use spintronics::magnon::chain::{ChainParameters, SpinChain};
use spintronics::Vector3;

fn bench_spinchain_creation_100(c: &mut Criterion) {
    let params = ChainParameters::permalloy();

    c.bench_function("spinchain_creation_100", |b| {
        b.iter(|| black_box(SpinChain::new(100, params.clone())))
    });
}

fn bench_spinchain_evolution_100_10steps(c: &mut Criterion) {
    let params = ChainParameters::permalloy();
    let h_ext = Vector3::new(0.0, 0.0, 0.1); // Small external field
    let dt = params.max_stable_dt();

    c.bench_function("spinchain_evolve_100spins_10steps", |b| {
        b.iter(|| {
            let mut chain = SpinChain::new_with_noise(100, params.clone(), 0.01);
            for _ in 0..10 {
                chain.evolve_heun(h_ext, dt);
            }
            black_box(&chain.spins);
        })
    });
}

fn bench_spinchain_evolution_1000_1step(c: &mut Criterion) {
    let params = ChainParameters::permalloy();
    let h_ext = Vector3::new(0.0, 0.0, 0.1);
    let dt = params.max_stable_dt();

    c.bench_function("spinchain_evolve_1000spins_1step", |b| {
        b.iter(|| {
            let mut chain = SpinChain::new_with_noise(1000, params.clone(), 0.01);
            chain.evolve_heun(h_ext, dt);
            black_box(&chain.spins);
        })
    });
}

criterion_group!(
    benches,
    bench_spinchain_creation_100,
    bench_spinchain_evolution_100_10steps,
    bench_spinchain_evolution_1000_1step
);
criterion_main!(benches);
