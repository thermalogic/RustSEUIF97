#![allow(warnings)]
///  cargo test --test tv_test
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use seuif97::*;

#[test]
fn test_tv_region1() {
    for i in 0..3 {
        let T: f64 = r1_pT_data[i].T;
        let p: f64 = r1_pT_data[i].p;
        let v: f64 =  pt(p, T - 273.15, OV);
        assert_approx_eq!(r1_pT_data[i].p, tv(T - 273.15, v, OP), 1.0e-4f64);
    }
}

#[test]
fn test_tv_region2() {
    for i in 0..3 {
        let T: f64 = r2_pT_data[i].T;
        let p: f64 = r2_pT_data[i].p;
        let v: f64 =  pt(p, T - 273.15, OV);
        assert_approx_eq!(r2_pT_data[i].p,  tv(T - 273.15, v, OP), 1.0e-4f64);
    }
}

#[test]
fn test_tv_region3() {
    for i in 0..3 {
        let t: f64 = r3_Td[i][0] - 273.15;
        let v: f64 = 1.0 / r3_Td[i][1];
        assert_approx_eq!(r3_Td[i][2], tv(t, v, OP), 1.0e-5f64);
        assert_approx_eq!(r3_Td[i][3], tv(t, v, OH), 1.0e-5f64);
        assert_approx_eq!(r3_Td[i][4], tv(t, v, OU), 1.0e-5f64);
        assert_approx_eq!(r3_Td[i][5], tv(t, v, OS), 1.0e-5f64);
        assert_approx_eq!(r3_Td[i][6], tv(t, v, OCP), 1.0e-5f64);
        assert_approx_eq!(r3_Td[i][7], tv(t, v, OW), 1.0e-5f64);
    }
}

#[test]
fn test_tv_region4() {
    let mut v: f64 = 0.0;
    let x: f64 = 0.35;
    for i in 0..3 {
        let t: f64 = r4_sat_Tp[i][0] - 273.15;
        v = tx(t, x, OV);
        assert_approx_eq!(x, tv(t, v, OX));
    }
}

#[test]
fn test_tv_region5() {
    for i in 0..2 {
        let t: f64 = r5_pT_data[i][0] - 273.15;
        let p: f64 = r5_pT_data[i][1];
        let v: f64 = pt(p, t, OV);
        assert_approx_eq!(r5_pT_data[i][1], tv(t, v, OP), 1.0e-4f64);
    }
}
