use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use seuif97::*;
use std::hint::black_box;
use std::time::Duration;

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("IF97_Region_Speed");
  
    group.warm_up_time(Duration::from_secs(5));
    group.measurement_time(Duration::from_secs(10));

    group.throughput(Throughput::Elements(1));
    
    //group.bench_function("ph2t_reg1", |b| b.iter(|| ph(black_box(3.0), black_box( 0.115331273e+3), OT)));
    //group.bench_function("ph2t_reg2", |b| b.iter(|| ph(black_box(0.001), black_box( 3000.0), OT)));
    //group.bench_function("ph2t_reg3", |b| b.iter(|| ph(black_box(20.0), black_box( 1700.0), OT)));
    //group.bench_function("ph2t_reg5", |b| b.iter(|| ph(black_box(0.5), black_box( 4527.49310), OT)));
    
    //group.bench_function("ps2t_reg1", |b| b.iter(|| ps(black_box(3.0), black_box( 0.392294792), OT)));
   // group.bench_function("ps2t_reg2", |b| b.iter(|| ps(black_box(0.1), black_box( 7.5), OT)));
   // group.bench_function("ps2t_reg3", |b| b.iter(|| ps(black_box(20.0), black_box( 3.8), OT)));
   // group.bench_function("ps2t_reg5", |b| b.iter(|| ps(black_box(0.5), black_box( 9.65408875), OT)));
   
    //group.bench_function("hs2t_reg1", |b| b.iter(|| hs(black_box(0.115331273e+3), black_box( 0.392294792), OT)));
    //group.bench_function("hs2t_reg2", |b| b.iter(|| hs(black_box(3000.01), black_box( 7.5), OT)));
    //group.bench_function("hs2t_reg3", |b| b.iter(|| hs(black_box(1700.0), black_box( 3.8), OT)));
    //group.bench_function("hs2t_reg5", |b| b.iter(|| hs(black_box(4527.49), black_box( 9.65408875), OT)));
   
   // group.bench_function("pv2t_reg1", |b| b.iter(|| pv(black_box(3.0), black_box( 0.100215168e-2), OT)));
   // group.bench_function("pv2t_reg2", |b| b.iter(|| pv(black_box(0.0035), black_box( 0.394913866E+02), OT)));
   // group.bench_function("pv2t_reg3", |b| b.iter(|| pv(black_box(20.0), black_box( 1.761696406e-3), OT)));
   // group.bench_function("pv2t_reg5", |b| b.iter(|| pv(black_box(0.5), black_box( 0.521976855e+4), OT)));
   
    group.bench_function("tv2p_reg1", |b| b.iter(|| tv(black_box(300.0-273.15), black_box( 0.100215168e-2), OP)));
    group.bench_function("tv2p_reg2", |b| b.iter(|| tv(black_box(300.0-273.15), black_box( 0.394913866E+02), OP)));
    group.bench_function("tv2p_reg3", |b| b.iter(|| tv(black_box(650.0-273.15), black_box( 0.0002), OP)));
    group.bench_function("tv2p_reg5", |b| b.iter(|| tv(black_box(1500.0-273.15), black_box( 0.521976855e+4), OP)));
  
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
