use criterion::{black_box, criterion_group, criterion_main, Criterion};

use if97::*;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("pt2h_reg1", |b| b.iter(|| pt(black_box(3.0), black_box(300.0 - 273.15), black_box(OH))));
    c.bench_function("pt2s_reg1", |b| b.iter(|| pt(black_box(3.0), black_box(300.0 - 273.15), black_box(OS))));

    c.bench_function("pt2h_reg2", |b| b.iter(|| pt(black_box(0.0035), black_box(300.0 - 273.15), black_box(OH))));
    c.bench_function("pt2s_reg2", |b| b.iter(|| pt(black_box(0.0035), black_box(300.0 - 273.15), black_box(OS))));

    c.bench_function("tv2h_reg3", |b| b.iter(|| tv(black_box(650.0 - 273.15), black_box(1.0 / 500.0), black_box(OH))));
    c.bench_function("tv2s_reg3", |b| b.iter(|| tv(black_box(650.0 - 273.15), black_box(1.0 / 500.0), black_box(OS))));

    c.bench_function("pT2h_reg5", |b| b.iter(|| pt(black_box(0.5), black_box(1500.0 - 273.15), black_box(OH))));
    c.bench_function("pT2s_reg5", |b| b.iter(|| pt(black_box(0.5), black_box(1500.0 - 273.15), black_box(OS))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
