use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use seuif97::*;
use std::hint::black_box;
use std::time::Duration;

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("IF97");
  
    group.warm_up_time(Duration::from_secs(5));
    group.measurement_time(Duration::from_secs(10));

    group.throughput(Throughput::Elements(1));
    
    //group.bench_function("pt2h_reg1", |b| b.iter(|| pt(black_box(3.0), black_box(300.0 - 273.15), OH)));
    //group.bench_function("pt2s_reg1", |b| b.iter(|| pt(black_box(3.0), black_box(300.0 - 273.15), OS)));
    group.bench_function("pt2v_reg1", |b| b.iter(|| pt(black_box(3.0), black_box(300.0 - 273.15), OV)));
    //group.bench_function("pt2w_reg1", |b| b.iter(|| pt(black_box(3.0), black_box(300.0 - 273.15), OW)));
    
    //group.bench_function("pt2h_reg2", |b| b.iter(|| pt(black_box(0.0035), black_box(300.0 - 273.15), OH)));
    //group.bench_function("pt2s_reg2", |b| b.iter(|| pt(black_box(0.0035), black_box(300.0 - 273.15), OS)));

    //group.bench_function("tv2h_reg3", |b| b.iter(|| tv(black_box(650.0 - 273.15), black_box(1.0 / 500.0), OH)));
    //group.bench_function("tv2s_reg3", |b| b.iter(|| tv(black_box(650.0 - 273.15), black_box(1.0 / 500.0), OS)));
    //group.bench_function("pt2h_reg3", |b| b.iter(|| pt(black_box(50.0), black_box(630.0-273.15), OH)));
    //group.bench_function("pt2s_reg3", |b| b.iter(|| pt(black_box(50.0), black_box(630.0-273.15), OS)));

    //group.bench_function("pt2h_reg5", |b| b.iter(|| pt(black_box(0.5), black_box(1500.0 - 273.15), OH)));
    //group.bench_function("pt2s_reg5", |b| b.iter(|| pt(black_box(0.5), black_box(1500.0-273.15), OS)));
    
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
