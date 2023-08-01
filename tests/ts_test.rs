#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use if97::*;

#[test]
fn test_ts_region1() {
    for i in 0..3 {
        let T: f64 = r1_pT_data[i].T;
        let p: f64 = r1_pT_data[i].p;
        let s: f64 = pt(p, T - 273.15, OS);
        assert_approx_eq!(r1_pT_data[i].p, ts(T - 273.15, s, OP), 1.0e-4f64);
    }
}

#[test]
fn test_ts_region2() {
    for i in 0..3 {
        let T: f64 = r2_pT_data[i].T;
        let p: f64 = r2_pT_data[i].p;
        let s: f64 = pt(p, T - 273.15, OS);
        assert_approx_eq!(r2_pT_data[i].p, ts(T - 273.15, s, OP), 1.0e-4f64);
    }
}

#[test]
fn test_ts_region3() {
    for i in 0..3 {
        let t: f64 = r3_Td[i][0] - 273.15;
        let d: f64 = r3_Td[i][1];
        let s = tv(t, 1.0 / d, OS);
        assert_approx_eq!(d, if97::ts(t, s, OD), 1.0e-4f64);
    }
}

#[test]
fn test_ts_region4() {
    let mut s: f64 = 0.0;
    let x: f64 = 0.35;
    for i in 0..3 {
        let T: f64 = r4_sat_Tp[i][0];
        s = tx(T - 273.15, x, OS);
        assert_approx_eq!(x, ts(T - 273.15, s, OX));
    }
}

#[test]
fn test_ts_region5() {
    for i in 0..2 {
        let T: f64 = r5_pT_data[i][0];
        let p: f64 = r5_pT_data[i][1];
        let s: f64 = pt(p, T - 273.15, OS);
        assert_approx_eq!(p, ts(T - 273.15, s, OP), 1.0e-4f64);
    }
}
