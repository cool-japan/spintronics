// This benchmark exercises `spintronics::magnon` (spin-chain evolution),
// which is excluded from wasm32 builds (see
// `#[cfg(not(target_arch = "wasm32"))]` on `pub mod magnon;` in
// `src/lib.rs`). `harness = false` benches must provide their own `main`, so
// every item is individually gated (a whole-file `#![cfg(...)]` would also
// erase the wasm32 stub `main` below, since it applies to the enclosing
// crate root as a unit) and a no-op `main` is supplied for wasm32.
#[cfg(not(target_arch = "wasm32"))]
use criterion::{criterion_group, criterion_main, Criterion};
#[cfg(not(target_arch = "wasm32"))]
use std::hint::black_box;

#[cfg(not(target_arch = "wasm32"))]
use spintronics::magnon::chain::{ChainParameters, SpinChain};
#[cfg(not(target_arch = "wasm32"))]
use spintronics::Vector3;

#[cfg(not(target_arch = "wasm32"))]
fn bench_spinchain_creation_100(c: &mut Criterion) {
    let params = ChainParameters::permalloy();

    c.bench_function("spinchain_creation_100", |b| {
        b.iter(|| black_box(SpinChain::new(100, params.clone())))
    });
}

#[cfg(not(target_arch = "wasm32"))]
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

#[cfg(not(target_arch = "wasm32"))]
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

#[cfg(not(target_arch = "wasm32"))]
criterion_group!(
    benches,
    bench_spinchain_creation_100,
    bench_spinchain_evolution_100_10steps,
    bench_spinchain_evolution_1000_1step
);
#[cfg(not(target_arch = "wasm32"))]
criterion_main!(benches);

#[cfg(target_arch = "wasm32")]
fn main() {}
