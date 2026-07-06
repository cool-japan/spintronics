use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use spintronics::material::Ferromagnet;

fn bench_ferromagnet_yig(c: &mut Criterion) {
    c.bench_function("ferromagnet_yig_creation", |b| {
        b.iter(|| black_box(Ferromagnet::yig()))
    });
}

fn bench_ferromagnet_permalloy(c: &mut Criterion) {
    c.bench_function("ferromagnet_permalloy_creation", |b| {
        b.iter(|| black_box(Ferromagnet::permalloy()))
    });
}

fn bench_ferromagnet_cofeb(c: &mut Criterion) {
    c.bench_function("ferromagnet_cofeb_creation", |b| {
        b.iter(|| black_box(Ferromagnet::cofeb()))
    });
}

fn bench_ferromagnet_iron(c: &mut Criterion) {
    c.bench_function("ferromagnet_iron_creation", |b| {
        b.iter(|| black_box(Ferromagnet::iron()))
    });
}

fn bench_ferromagnet_cobalt(c: &mut Criterion) {
    c.bench_function("ferromagnet_cobalt_creation", |b| {
        b.iter(|| black_box(Ferromagnet::cobalt()))
    });
}

fn bench_ferromagnet_nickel(c: &mut Criterion) {
    c.bench_function("ferromagnet_nickel_creation", |b| {
        b.iter(|| black_box(Ferromagnet::nickel()))
    });
}

fn bench_ferromagnet_cofe(c: &mut Criterion) {
    c.bench_function("ferromagnet_cofe_creation", |b| {
        b.iter(|| black_box(Ferromagnet::cofe()))
    });
}

criterion_group!(
    benches,
    bench_ferromagnet_yig,
    bench_ferromagnet_permalloy,
    bench_ferromagnet_cofeb,
    bench_ferromagnet_iron,
    bench_ferromagnet_cobalt,
    bench_ferromagnet_nickel,
    bench_ferromagnet_cofe
);
criterion_main!(benches);
