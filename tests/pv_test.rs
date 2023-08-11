#![allow(warnings)]
///  cargo test --test pv_test
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use if97::*;

#[test]
fn test_pv_region1() {
    for i in 0..3 {
        let t: f64 = r1_pT_data[i].T - 273.15;
        let p: f64 = r1_pT_data[i].p;
        let v: f64 = if97::pt(p, t, OV);
        assert_approx_eq!(r1_pT_data[i].T - 273.15, pv(p, v, OT), 1.0e-4f64);
    }
}

#[test]
fn test_pv_region2() {
    for i in 0..3 {
        let t: f64 = r2_pT_data[i].T - 273.15;
        let p: f64 = r2_pT_data[i].p;
        let v: f64 = pt(p, t, OV);
        assert_approx_eq!(r2_pT_data[i].T - 273.15, pv(p, v, OT), 1.0e-4f64);
    }
}

#[test]
fn test_pv_region3() {
    for i in 0..3 {
        let T: f64 = r3_Td[i][0];
        let d: f64 = r3_Td[i][1];
        let p = tv(T - 273.15, 1.0 / d, OP);
        assert_approx_eq!(r3_Td[i][0] - 273.15, pv(p, 1.0 / d, OT), 1.0e-3f64)
    }
}

#[test]
fn test_pv_region4() {
    let mut v: f64 = 0.0;
    let mut vw: f64 = 0.0;
    let mut vs: f64 = 0.0;
    let x: f64 = 0.35;
    for i in 0..3 {
        let p: f64 = r4_sat_Tp[i][1];
        v = px(p, x, OV);
        assert_approx_eq!(x, pv(p, v, OX));
    }
}

#[test]
fn test_pv_region5() {
    for i in 0..2 {
        let t: f64 = r5_pT_data[i][0] - 273.15;
        let p: f64 = r5_pT_data[i][1];
        let v: f64 = pt(p, t, OV);
        assert_approx_eq!(r5_pT_data[i][0] - 273.15, pv(p, v, OT), 1.0e-4f64);
    }
}
