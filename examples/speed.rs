#![allow(warnings)]
/// The example to test the  speed
use if97::*;
use std::time::Instant;

const prop_map_pt: [(&str, i32); 6] = [("V", OV), ("H", OH), ("S", OS), ("U", OU), ("CP", OCP), ("W", OW)];

fn speed_experiment(v1: f64, v2: f64, func: fn(f64, f64, (i32, i32)) -> f64, prop_map: &[(&str, i32)], reg: i32) {
    println!(" v1={:.6} v2={:.6} ", v1, v2);
    for e in prop_map {
        let value: f64 = func(v1, v2, (e.1, reg));
        let now = Instant::now();
        for _ in 0..100000u128 {
            std::hint::black_box(func(v1, v2, (e.1, reg)));
        }
        let elapsed_time = now.elapsed();
        println!("\t {} = {} \n\t\t Running  {:?}", e.0, value, elapsed_time);
    }
}

fn speed_region1() {
    println!("Region 1 ");
    let p: f64 = 3.0;
    let t: f64 = 300.0 - 273.15;
    speed_experiment(p, t, pt, &prop_map_pt, 1);
}

fn speed_region2() {
    println!("Region 2 ");
    let p: f64 = 0.0035;
    let t: f64 = 300.0 - 273.15;
    speed_experiment(p, t, pt, &prop_map_pt, 2);
}

fn speed_region3() {
    let prop_map_Td: [(&str, i32); 6] = [("P", OP), ("H", OH), ("S", OS), ("U", OU), ("CP", OCP), ("W", OW)];
    println!("Region 3 ");
    let T: f64 = 650.0;
    let d: f64 = 500.0;
    speed_experiment(T - 273.15, 1.0 / d, tv, &prop_map_Td, 3);
}

fn speed_region5() {
    println!("Region  5 ");
    let p: f64 = 0.5;
    let t: f64 = 1500.0 - 273.15;
    speed_experiment(p, t, pt, &prop_map_pt, 5);
}

fn main() {
    speed_region1();
    speed_region2();
    speed_region3();
    speed_region5();
}
