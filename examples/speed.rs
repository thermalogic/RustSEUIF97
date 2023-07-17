#![allow(warnings)]
/// Example: speeed
use if97::*;
use std::time::Instant;

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

fn main() {
    const prop_map: [(&str, i32); 6] = [
        ("V", OV),
        ("H", OH),
        ("S", OS),
        ("U", OU),
        ("CP", OCP),
        ("W", OW),
    ];
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

    for k in 0..3 {
        let mut p: f64 = data[k].p;
        let mut t: f64 = data[k].T - 273.15;
        println!(" p={:.6} t={:.6} ", p, t);

        for i in 0..6 {
            let value: f64 = pt(p, t, prop_map[i].1);
            println!("\t {} = {}", prop_map[i].0, value);
            let now = Instant::now();
            for _ in 0..100000u128 {
                std::hint::black_box(pt(p, t, prop_map[i].1));
            }
            let elapsed_time = now.elapsed();
            println!("\t\t {} Running  {:?}", prop_map[i].0, elapsed_time);
        }
    }
}
