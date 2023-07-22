#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use if97::*;

#[test]
fn test_hx_region4() {
    let mut h: f64 = 0.0;
    let mut hw: f64 = 0.0;
    let mut hs: f64 = 0.0;
    let x: f64 = 0.35;
    for i in 0..3 {
        let t: f64 = r4_sat_Tp[i][0]-273.15;
        h = if97::tx(t,x, OH);
        assert_approx_eq!(t, if97::hx(h, x, OT), 1.0e-2);
    }
}

#[test]
fn test_sx_region4() {
    let mut s: f64 = 0.0;
    let mut sw: f64 = 0.0;
    let mut ss: f64 = 0.0;
    let x: f64 = 0.35;
    for i in 0..3 {
        let t: f64 = r4_sat_Tp[i][0]-273.15;
        s = if97::tx(t,x, OS);
        assert_approx_eq!(t, if97::sx(s, x, OT), 1.0e-2);
    }
}
