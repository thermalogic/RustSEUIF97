#![allow(warnings)]
use std::collections::HashMap;
use std::time::Instant;
pub mod algo;
pub mod common;
pub mod r1;
pub mod r2;
pub mod r3;
pub mod r4;
pub mod r5;

use common::constant::*;
use r1::*;
use r2::*;
use r3::*;
use r4::*;
use r5::*;

pub struct TestData {
    pub p: f64,
    pub T: f64,
    pub v: f64,
    pub h: f64,
    pub s: f64,
    pub u: f64,
    pub cp: f64,
    pub w: f64,
}

const prop_map: [(&str, i32);8] =[
    ("P", OP),
    ("T", OT),
    ("V", OV),
    ("H", OH),
    ("S", OS),
    ("U", OU),
    ("CP", OCP),
    ("W", OW),
];

fn speed_region1() {
    // Table 5，Page 9: p,T,v,h,u,s,cp,w
    const p: [f64; 3] = [3.0, 80.0, 3.0];
    const T: [f64; 3] = [300.0, 300.0, 500.0];

    const v: [f64; 3] = [0.100215168e-2, 0.971180894e-3, 0.120241800e-2];
    const h: [f64; 3] = [0.115331273e+3, 0.184142828e+3, 0.975542239e+3];
    const s: [f64; 3] = [0.392294792, 0.368563852, 0.258041912e1];
    const u: [f64; 3] = [0.112324818e3, 0.106448356e3, 0.971934985e3];
    const cp: [f64; 3] = [0.417301218e1, 0.401008987e1, 0.465580682e1];
    const w: [f64; 3] = [0.150773921e4, 0.163469054e4, 0.124071337e4];

    println!("Region 1 ");

    let mut p1: f64 = p[0];
    let mut t1: f64 = T[0] - 273.15;
    println!(" p={} t={} ", p1, t1);

    for i in 0..8 {
        let value: f64 = pt_reg1(p1, t1,prop_map[i].1);
        println!("\t {} ={}",prop_map[i].0, value);
        let now = Instant::now();
        for _ in 0..100000u128 {
            std::hint::black_box(pt_reg1(p1, t1,prop_map[i].1));
        }
        let elapsed_time = now.elapsed();
        println!("\t\t {} Running  {:?}", prop_map[i].0, elapsed_time);
    }
}

fn speed_region2() {
    const data: [TestData; 3] = [
        TestData {
            p: 0.0035,
            T: 300.0,
            v: 0.394913866E+02,
            h: 0.254991145E+04,
            u: 0.241169160E+04,
            s: 0.852238967E+01,
            cp: 0.191300162E+01,
            w: 0.427920172E+03,
        },
        TestData {
            p: 0.0035,
            T: 700.0,
            v: 0.923015898E+02,
            h: 0.333568375E+04,
            u: 0.301262819E+04,
            s: 0.101749996E+02,
            cp: 0.208141274E+01,
            w: 0.644289068E+03,
        },
        TestData {
            p: 30.0,
            T: 700.0,
            v: 0.542946619E-02,
            h: 0.263149474E+04,
            u: 0.246861076E+04,
            s: 0.517540298E+01,
            cp: 0.103505092E+02,
            w: 0.480386523E+03,
        },
    ];

    println!("Region 2 ");
    let mut p: f64 = data[0].p;
    let mut t: f64 = data[0].T - 273.15;
    println!("\t p={} t={} ", p, t);

    for i in 0..8 {
        let value: f64 = pt_reg2(p, t,prop_map[i].1);
        println!("\t {} ={}", prop_map[i].0, value);
        let now = Instant::now();
        for _ in 0..100000u128 {
            std::hint::black_box(pt_reg2(p, t, prop_map[i].1));
        }
        let elapsed_time = now.elapsed();
        println!("\t\t {} Running  {:?}", prop_map[i].0, elapsed_time);
    }
}

fn speed_region3() {
    // Table 33. Thermodynamic property values calculated from Eq. (28) for selected values of T and  a
    // T,d,p,h,u,s,cp,w
    const tab33: [[f64; 8]; 3] = [
        [
            650.0,
            500.0,
            0.255837018E2,
            0.186343019E4,
            0.181226279E4,
            0.405427273E1,
            0.138935717E2,
            0.502005554E3,
        ],
        [
            650.0,
            200.0,
            0.222930643E2,
            0.237512401E4,
            0.226365868E4,
            0.485438792E1,
            0.446579342E2,
            0.383444594E3,
        ],
        [
            750.0,
            500.0,
            0.783095639E2,
            0.225868845E4,
            0.210206932E4,
            0.446971906E1,
            0.634165359E1,
            0.760696041E3,
        ],
    ];
    println!("Region 3 ");
    let T: f64 = tab33[0][0];
    let d: f64 = tab33[0][1];
    println!("\t T={} d={} ", T, d);

    for i in 0..8 {
        let value: f64 = Td_reg3(T, d, prop_map[i].1);
        println!("\t {}={}", prop_map[i].0, value);
        let now = Instant::now();
        for _ in 0..100000u128 {
            std::hint::black_box(Td_reg3(T, d, prop_map[i].1));
        }
        let elapsed_time = now.elapsed();
        println!("\t\t {} Running  {:?}", prop_map[i].0, elapsed_time);
    }
}

fn speed_region5() {
    // Table 42, page 40 T,p  v,h,u,s,cp,w
    const data: [[f64; 8]; 3] = [
        [
            1500.0,
            0.5,
            1.38455090,
            0.521976855e+4,
            4527.49310,
            9.65408875,
            2.61609445,
            917.068690,
        ],
        [
            1500.0,
            30.0,
            0.0230761299,
            5167.23514,
            4474.95124,
            7.72970133,
            2.72724317,
            928.548002,
        ],
        [
            2000.0,
            30.0,
            0.0311385219,
            6571.22604,
            5637.07038,
            8.53640523,
            2.88569882,
            1067.36948,
        ],
    ];

    println!("Region  5 ");

    let mut p1: f64 = data[0][1];
    let mut t1: f64 = data[0][0] - 273.15;
    println!(" p={} t={} ", p1, t1);
    for i in 0..8 {
        let value: f64 = pt_reg5(p1, t1, prop_map[i].1);
        println!("\t {} ={}", prop_map[i].0, value);
        let now = Instant::now();
        for _ in 0..100000u128 {
            std::hint::black_box(pt_reg5(p1, t1, prop_map[i].1));
        }
        let elapsed_time = now.elapsed();
        println!("\t\t {} Running  {:?}", prop_map[i].0, elapsed_time);
    }
}

fn main() {
    speed_region1();
    speed_region2();
    speed_region3();
    speed_region5();
}
