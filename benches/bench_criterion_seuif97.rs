#![allow(warnings)]
/// IAPWS-IF97 Comprehensive Performance Benchmark
///
/// This benchmark uses test data from tests/common/mod.rs to ensure
/// consistency between testing and performance evaluation.
///
/// # Usage
/// ```bash
/// # Run all benchmarks
/// cargo bench
///
/// # Run specific region
/// cargo bench -- "IF97_Region1"
///
/// # Run specific function type
/// cargo bench -- "pt2h"
/// ```
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use seuif97::*;
use std::hint::black_box;
use std::time::Duration;

// Include test data directly from tests/common/mod.rs
include!("../tests/common/mod.rs");

/// Benchmark Region 1 using test data from tests/common/mod.rs
fn benchmark_region1(c: &mut Criterion) {
    let mut group = c.benchmark_group("IF97_Region1");
    group.warm_up_time(Duration::from_secs(3));
    group.measurement_time(Duration::from_secs(5));
    group.throughput(Throughput::Elements(1));
    
    // Use first test case from r1_pT_data
    let data = &r1_pT_data[0];
    let t_celsius = data.T - 273.15;
    
    // Forward functions: p, T → properties
    group.bench_with_input(
        BenchmarkId::new("pt2h", format!("p={},T={}", data.p, t_celsius)),
        &(data.p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OH))),
    );

    group.bench_with_input(
        BenchmarkId::new("pt2s", format!("p={},T={}", data.p, t_celsius)),
        &(data.p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OS))),
    );
    group.bench_with_input(
        BenchmarkId::new("pt2v", format!("p={},T={}", data.p, t_celsius)),
        &(data.p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OV))),
    );
    group.bench_with_input(
        BenchmarkId::new("pt2u", format!("p={},T={}", data.p, t_celsius)),
        &(data.p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OU))),
    );
    group.bench_with_input(
        BenchmarkId::new("pt2w", format!("p={},T={}", data.p, t_celsius)),
        &(data.p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OW))),
    );
    
    // Backward functions
    group.bench_function("pv2T", |b| {
        b.iter(|| pv(black_box(data.p), black_box(data.v), black_box(OT)))
    });
    group.bench_function("Tv2p", |b| {
        b.iter(|| tv(black_box(t_celsius), black_box(data.v), black_box(OP)))
    });
    group.bench_function("Th2p", |b| {
        b.iter(|| th(black_box(t_celsius), black_box(data.h), black_box(OP)))
    });
    group.bench_function("Ts2p", |b| {
        b.iter(|| ts(black_box(t_celsius), black_box(data.s), black_box(OP)))
    });
    group.bench_function("ph2t", |b| {
        b.iter(|| ph(black_box(data.p), black_box(data.h), black_box(OT)))
    });
    group.bench_function("ps2t", |b| {
        b.iter(|| ps(black_box(data.p), black_box(data.s), black_box(OT)))
    });
    
    group.bench_with_input(
        BenchmarkId::new("hs2p", format!("h={},s={}", data.h, data.s)),
        &(data.h, data.s),
        |b, (h, s)| b.iter(|| hs(black_box(*h), black_box(*s), black_box(OP))),
    );
    group.finish();
}

/// Benchmark Region 2 using test data from tests/common/mod.rs
fn benchmark_region2(c: &mut Criterion) {
    let mut group = c.benchmark_group("IF97_Region2");
    group.warm_up_time(Duration::from_secs(3));
    group.measurement_time(Duration::from_secs(5));
    group.throughput(Throughput::Elements(1));
    
    let data = &r2_pT_data[0];
    let t_celsius = data.T - 273.15;
    
    // Forward functions
    group.bench_with_input(
        BenchmarkId::new("pt2h", format!("p={},T={}", data.p, t_celsius)),
        &(data.p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OH))),
    );
    group.bench_with_input(
        BenchmarkId::new("pt2s", format!("p={},T={}", data.p, t_celsius)),
        &(data.p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OS))),
    );
    
    // Backward functions
    group.bench_function("pv2T", |b| {
        b.iter(|| pv(black_box(data.p), black_box(data.v), black_box(OT)))
    });
    group.bench_function("ph2t", |b| {
        b.iter(|| ph(black_box(r2_ph_2a[0][0]), black_box(r2_ph_2a[0][1]), black_box(OT)))
    });
    group.bench_function("ps2t", |b| {
        b.iter(|| ps(black_box(r2_ps_2a[0][0]), black_box(r2_ps_2a[0][1]), black_box(OT)))
    });
    group.bench_function("hs2t", |b| {
        b.iter(|| hs(black_box(r2_hs_reg2a[0][0]), black_box(r2_hs_reg2a[0][1]), black_box(OT)))
    });
    
    group.finish();
}

/// Benchmark Region 3 using test data from tests/common/mod.rs
fn benchmark_region3(c: &mut Criterion) {
    let mut group = c.benchmark_group("IF97_Region3");
    group.warm_up_time(Duration::from_secs(3));
    group.measurement_time(Duration::from_secs(5));
    group.throughput(Throughput::Elements(1));
    
    let data = &r3_Td[0];
    let t = data[0];
    let d = data[1];
    let p = data[2];
    let v = 1.0 / d;
    let h = data[3];
    let s = data[5];
    let t_celsius = t - 273.15;
    
    // Forward functions
    group.bench_with_input(
        BenchmarkId::new("pt2h", format!("p={},T={}", p, t_celsius)),
        &(p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OH))),
    );
    group.bench_with_input(
        BenchmarkId::new("pt2s", format!("p={},T={}", p, t_celsius)),
        &(p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OS))),
    );
    
    // Backward functions
    group.bench_function("pv2T", |b| {
        b.iter(|| pv(black_box(p), black_box(v), black_box(OT)))
    });
    group.bench_function("tv2h", |b| {
        b.iter(|| tv(black_box(t_celsius), black_box(v), black_box(OH)))
    });
    group.bench_function("tv2s", |b| {
        b.iter(|| tv(black_box(t_celsius), black_box(v), black_box(OS)))
    });
    group.bench_function("Th2d", |b| {
        b.iter(|| th(black_box(t_celsius), black_box(h), black_box(OD)))
    });
    group.bench_function("Ts2d", |b| {
        b.iter(|| ts(black_box(t_celsius), black_box(s), black_box(OD)))
    });
    group.bench_function("ph2t", |b| {
        b.iter(|| ph(black_box(r3_phTv_3a[0][0]), black_box(r3_phTv_3a[0][1]), black_box(OT)))
    });
    group.bench_function("ps2t", |b| {
        b.iter(|| ps(black_box(r3_psTv_3a[0][0]), black_box(r3_psTv_3a[0][1]), black_box(OT)))
    });
    group.bench_function("hs2t", |b| {
        b.iter(|| hs(black_box(r3_hsp_3a[0][0]), black_box(r3_hsp_3a[0][1]), black_box(OT)))
    });
    
    group.finish();
}

/// Benchmark Region 5 using test data from tests/common/mod.rs
fn benchmark_region5(c: &mut Criterion) {
    let mut group = c.benchmark_group("IF97_Region5");
    group.warm_up_time(Duration::from_secs(3));
    group.measurement_time(Duration::from_secs(5));
    group.throughput(Throughput::Elements(1));
    
    let data = &r5_pT_data[0];
    let t = data[0];
    let p = data[1];
    let v = data[2];
    let h = data[3];
    let s = data[5];
    let t_celsius = t - 273.15;
    
    // Forward functions
    group.bench_with_input(
        BenchmarkId::new("pt2h", format!("p={},T={}", p, t_celsius)),
        &(p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OH))),
    );
    group.bench_with_input(
        BenchmarkId::new("pt2s", format!("p={},T={}", p, t_celsius)),
        &(p, t_celsius),
        |b, (p, t)| b.iter(|| pt(black_box(*p), black_box(*t), black_box(OS))),
    );
    
    // Backward functions
    group.bench_function("pv2T", |b| {
        b.iter(|| pv(black_box(p), black_box(v), black_box(OT)))
    });
    group.bench_function("ph2t", |b| {
        b.iter(|| ph(black_box(p), black_box(h), black_box(OT)))
    });
    group.bench_function("ps2t", |b| {
        b.iter(|| ps(black_box(p), black_box(s), black_box(OT)))
    });
    group.bench_function("hs2t", |b| {
        b.iter(|| hs(black_box(h), black_box(s), black_box(OT)))
    });
    
    group.finish();
}

/// Main benchmark entry point
fn criterion_benchmark(c: &mut Criterion) {
    benchmark_region1(c);
   // benchmark_region2(c);
   // benchmark_region3(c);
    //benchmark_region5(c);
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);