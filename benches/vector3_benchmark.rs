use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use spintronics::Vector3;

fn bench_vector3_creation(c: &mut Criterion) {
    c.bench_function("vector3_creation", |b| {
        b.iter(|| black_box(Vector3::new(1.0, 2.0, 3.0)))
    });
}

fn bench_vector3_addition(c: &mut Criterion) {
    let v1 = Vector3::new(1.0, 2.0, 3.0);
    let v2 = Vector3::new(4.0, 5.0, 6.0);
    c.bench_function("vector3_add_1000", |b| {
        b.iter(|| {
            let mut result = Vector3::zero();
            for _ in 0..1000 {
                result = result + black_box(v1) + black_box(v2);
            }
            black_box(result)
        })
    });
}

fn bench_vector3_cross(c: &mut Criterion) {
    let v1 = Vector3::new(1.0, 2.0, 3.0);
    let v2 = Vector3::new(4.0, 5.0, 6.0);
    c.bench_function("vector3_cross_product", |b| {
        b.iter(|| black_box(v1.cross(&v2)))
    });
}

fn bench_vector3_dot(c: &mut Criterion) {
    let v1 = Vector3::new(1.0, 2.0, 3.0);
    let v2 = Vector3::new(4.0, 5.0, 6.0);
    c.bench_function("vector3_dot_product", |b| b.iter(|| black_box(v1.dot(&v2))));
}

fn bench_vector3_normalize(c: &mut Criterion) {
    let v = Vector3::new(1.0, 2.0, 3.0);
    c.bench_function("vector3_normalize", |b| b.iter(|| black_box(v.normalize())));
}

fn bench_vector3_magnitude(c: &mut Criterion) {
    let v = Vector3::new(1.0, 2.0, 3.0);
    c.bench_function("vector3_magnitude", |b| b.iter(|| black_box(v.magnitude())));
}

criterion_group!(
    benches,
    bench_vector3_creation,
    bench_vector3_addition,
    bench_vector3_cross,
    bench_vector3_dot,
    bench_vector3_normalize,
    bench_vector3_magnitude
);
criterion_main!(benches);
