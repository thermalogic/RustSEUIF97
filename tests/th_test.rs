#![allow(warnings)]
///  cargo test --test th_test
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use seuif97::*;

#[test]
fn test_th_region1() {
    for i in 0..3 {
        let T: f64 = r1_pT_data[i].T;
        let p: f64 = r1_pT_data[i].p;
        let h: f64 = pt(p, T - 273.15, OH);
        assert_approx_eq!(r1_pT_data[i].p, th(T - 273.15, h, OP), 1.0e-2f64);
    }
}

#[test]
fn test_th_region2() {
    for i in 0..3 {
        let T: f64 = r2_pT_data[i].T;
        let p: f64 = r2_pT_data[i].p;
        let h: f64 = pt(p, T - 273.15, OH);
        assert_approx_eq!(r2_pT_data[i].p, th(T - 273.15, h, OP), 1.0e-4f64);
    }
}

#[test]
fn test_th_region3() {
    for i in 0..3 {
        let t: f64 = r3_Td[i][0] - 273.15;
        let d: f64 = r3_Td[i][1];
        let h = tv(t, 1.0 / d, OH);
        assert_approx_eq!(d, th(t, h, OD), 1.0e-4f64);
    }
}

#[test]
fn test_th_region4() {
    let mut h: f64 = 0.0;
    let mut hw: f64 = 0.0;
    let mut hs: f64 = 0.0;
    let x: f64 = 0.35;
    for i in 0..3 {
        let T: f64 = r4_sat_Tp[i][0];
        h = tx(T - 273.15, x, OH);
        assert_approx_eq!(x, th(T - 273.15, h, OX));
    }
}

#[test]
fn test_th_region5() {
    for i in 0..2 {
        let T: f64 = r5_pT_data[i][0];
        let p: f64 = r5_pT_data[i][1];
        let h: f64 = pt(p, T - 273.15, OH);
        assert_approx_eq!(p, th(T - 273.15, h, OP), 1.0e-4f64);
    }
}
