#![allow(warnings)]
/// Benchmark example for SEUIF97 performance testing.
///
/// This example benchmarks the performance of SEUIF97 library across different IAPWS-IF97 regions.
/// It compares the execution time between:
/// - `NoRegion`: Property calculation without region specification (auto-detect)
/// - `WithRegion`: Property calculation with explicit region specification
/// - `Overhead`: The additional time cost introduced by region judgment logic
///
/// Features:
/// - Tests all thermodynamic properties: p, t, V, H, S, U, CP, W
/// - Covers Region 1, 2, 3, 4, and 5
/// - Uses standard test data from tests/common/mod.rs
/// - Known input parameters are marked with `*` in the output table
/// - Smart number formatting: adjusts decimal places based on value magnitude
///
/// Run:
/// ```bash
/// cargo run -r --example bench_seuif97
/// ```
use seuif97::*;
use std::time::Instant;

// Include standard test data
include!("../tests/common/mod.rs");

const prop_map_pt: [(&str, i32); 8] = [("p", OP), ("t", OT),("V", OV), ("H", OH), ("S", OS), ("U", OU), ("CP", OCP), ("W", OW)];

fn format_value(value: f64) -> String {
    let abs_value = value.abs();
    if abs_value >= 10.0 {
        format!("{:.3}", value)
    } else if abs_value >= 1.0 {
        format!("{:.4}", value)
    } else {
        format!("{:.6}", value)
    }
}

/// Benchmark result for a single property
struct BenchmarkRow {
    name: String,
    value: String,
    time_no_region: String,
    time_with_region: String,
    time_diff: String,
}

/// Collect benchmark timing data for all properties
fn collect_benchmark_data<F>(
    v1: f64, v2: f64, func: F, prop_map: &[(&str, i32)], reg: i32, known_props: &[&str],
) -> Vec<BenchmarkRow>
where
    F: Fn(f64, f64, o_id_region_args) -> f64,
{
    const BENCH_ITERATIONS: u128 = 1_000_00;
    let mut rows = Vec::new();
    for e in prop_map {
        let value: f64 = func(v1, v2, (e.1, reg).into());
        let now = Instant::now();
        for _ in 0..BENCH_ITERATIONS {
            std::hint::black_box(func(std::hint::black_box(v1), std::hint::black_box(v2), std::hint::black_box(e.1.into())));
        }
        let elapsed_no_region = now.elapsed();
        let now = Instant::now();
        for _ in 0..BENCH_ITERATIONS {
            std::hint::black_box(func(std::hint::black_box(v1), std::hint::black_box(v2), std::hint::black_box((e.1, reg).into())));
        }
        let elapsed_with_region = now.elapsed();
        let ns_no = elapsed_no_region.as_nanos() as f64 / BENCH_ITERATIONS as f64;
        let ns_with = elapsed_with_region.as_nanos() as f64 / BENCH_ITERATIONS as f64;
        let ns_diff = ns_no - ns_with;
        let is_known = known_props.contains(&e.0);
        let value_str = if is_known {
            format!("{}*", format_value(value))
        } else {
            format_value(value)
        };
        rows.push(BenchmarkRow {
            name: e.0.to_string(),
            value: value_str,
            time_no_region: if is_known { "-".to_string() } else { format!("{:.1}", ns_no) },
            time_with_region: if is_known { "-".to_string() } else { format!("{:.1}", ns_with) },
            time_diff: if is_known { "-".to_string() } else { format!("{:.1}", ns_diff) },
        });
    }
    rows
}

/// Print benchmark results as a formatted table
fn print_benchmark_table(rows: &[BenchmarkRow]) {
    let mut widths: Vec<usize> = (0..rows.len())
        .map(|i| {
            rows[i].name.len()
                .max(rows[i].value.len())
                .max(rows[i].time_no_region.len())
                .max(rows[i].time_with_region.len())
                .max(rows[i].time_diff.len())
        })
        .collect();
    let label_width = "Symbol".len().max("Value".len()).max("NoRegion".len()).max("WithRegion".len()).max("Overhead".len());
    widths.insert(0, label_width);
    let print_separator = |widths: &[usize]| {
        for (i, &w) in widths.iter().enumerate() {
            print!("{}", "-".repeat(w));
            if i < widths.len() - 1 { print!("  "); }
        }
        println!();
    };
    let print_row = |label: &str, items: &[String], center: bool| {
        if center {
            let padding = widths[0] - label.len();
            let left = padding / 2;
            print!("{}{}{}", " ".repeat(left), label, " ".repeat(padding - left));
        } else {
            print!("{:>width$}", label, width = widths[0]);
        }
        for (i, item) in items.iter().enumerate() {
            print!("  ");
            if center {
                let padding = widths[i + 1] - item.len();
                let left = padding / 2;
                print!("{}{}{}", " ".repeat(left), item, " ".repeat(padding - left));
            } else {
                print!("{:>width$}", item, width = widths[i + 1]);
            }
        }
        println!();
    };
    print_row("Symbol", &rows.iter().map(|r| r.name.clone()).collect::<Vec<_>>(), true);
    print_separator(&widths);
    print_row("Value", &rows.iter().map(|r| r.value.clone()).collect::<Vec<_>>(), false);
    print_row("NoRegion", &rows.iter().map(|r| r.time_no_region.clone()).collect::<Vec<_>>(), false);
    print_row("WithRegion", &rows.iter().map(|r| r.time_with_region.clone()).collect::<Vec<_>>(), false);
    print_row("Overhead", &rows.iter().map(|r| r.time_diff.clone()).collect::<Vec<_>>(), false);
    println!();
}

/// Run benchmark and print results
fn benchmark_experiment<F>(v1: f64, v2: f64, func: F, prop_map: &[(&str, i32)], reg: i32, known_props: &[&str])
where
    F: Fn(f64, f64, o_id_region_args) -> f64,
{
    let rows = collect_benchmark_data(v1, v2, func, prop_map, reg, known_props);
    print_benchmark_table(&rows);
}

fn benchmark_region_pt(p:f64, t:f64,region:i32) {
    println!("Region {}: p={:.3} t={:.3} ", region, p, t);
    benchmark_experiment(p, t, pt, &prop_map_pt, region, &["p", "t"]);
    let h = pt(p, t, (OH,region));
    println!("Region {}: p={:.3} h={:.3} ", region, p, h);
    benchmark_experiment(p, h, ph, &prop_map_pt, region, &["p", "H"]);
    let s = pt(p, t, (OS,region));
    println!("Region {}: p={:.3} s={:.3} ", region, p, s);
    benchmark_experiment(p, s, ps, &prop_map_pt, region, &["p", "S"]);
    println!("Region {}: h={:.3} s={:.3} ", region, h, s);
    benchmark_experiment(h, s, hs, &prop_map_pt, region, &["H", "S"]);    
}

fn benchmark_region1() {
    let data = &r1_pT_data[0];
    let p = data.p;
    let t = data.T - 273.15;
    benchmark_region_pt(p, t, 1)
}


fn benchmark_region2() {
    let data = &r2_pT_data[0];
    let p = data.p;
    let t = data.T - 273.15;
    benchmark_region_pt(p, t, 2)
}

fn benchmark_region3() {
    let data = &r3_Td[0];
    let t_kelvin = data[0];
    let d = data[1];
    let t = t_kelvin - 273.15;
    println!("Region 3 t={:.3} d={:.3} ", t, d);
    benchmark_experiment(t, 1.0 / d, tv, &prop_map_pt, 3, &["t", "V"]);
    let p: f64 = tv(t, 1.0 / d, (OP, 3));
    let t: f64 = tv(t, 1.0 / d, (OT, 3));
    benchmark_region_pt(p, t, 3)
}

fn benchmark_region4() {
    let data = &r4_hsT[0];
    let h = data[0];
    let s = data[1];
    println!("Region 4 h={:.3} s={:.3} ", h, s);
    benchmark_experiment(h, s, hs, &prop_map_pt, 4, &["H", "S"]);
}

fn benchmark_region5() {
    let data = &r5_pT_data[0];
    let p = data[1];
    let t = data[0] - 273.15;
    benchmark_region_pt(p, t, 5)
}


fn main() {
    benchmark_region1();
    benchmark_region2();
    benchmark_region3();
    benchmark_region4();
    benchmark_region5();
}