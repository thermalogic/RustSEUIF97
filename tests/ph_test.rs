#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use if97::*;

#[test]
fn test_region1_ph() {
    for i in 0..3 {
        assert_approx_eq!(r1_phT[i][2] - 273.15, ph(r1_phT[i][0], r1_phT[i][1], OT), 1.0e-1f64);
    }
}

#[test]
fn test_region2_ph() {
    let mut p: f64 = 0.0;
    let mut h: f64 = 0.0;
    for i in 0..3 {
        p = r2_ph_2a[i][0];
        h = r2_ph_2a[i][1];
        assert_approx_eq!(r2_ph_2a[i][2], ph(p, h, OT) + 273.15, 1.0e-1f64);
        p = r2_ph_2b[i][0];
        h = r2_ph_2b[i][1];
        assert_approx_eq!(r2_ph_2b[i][2], ph(p, h, OT) + 273.15, 1.0e-1f64);
        p = r2_ph_2c[i][0];
        h = r2_ph_2c[i][1];
        assert_approx_eq!(r2_ph_2c[i][2], ph(p, h, OT) + 273.15, 1.0e-1f64);
    }
}

#[test]
fn test_region2_ph_iter() {
    let mut h: f64 = 0.0;
    for i in 0..3 {
        let p: f64 = r2_pT_data[i].p;
        let t: f64 = r2_pT_data[i].T - 273.15;
        h = pt(p, t, OH);
        assert_approx_eq!(t, ph(p, h, OT), 1.0e-2f64);
    }
}

#[test]
fn test_region3_ph() {
    for i in 0..3 {
        let mut p: f64 = r3_phTv_3a[i][0];
        let mut h: f64 = r3_phTv_3a[i][1];
        assert_approx_eq!(r3_phTv_3a[i][2] - 273.15, ph(p, h, OT));
        assert_approx_eq!(r3_phTv_3a[i][3], ph(p, h, OV));
        p = r3_phTv_3b[i][0];
        h = r3_phTv_3b[i][1];
        assert_approx_eq!(r3_phTv_3b[i][2] - 273.15, ph(p, h, OT));
        assert_approx_eq!(r3_phTv_3b[i][3], ph(p, h, OV));
    }
}

#[test]
fn test_region5_ph() {
    for i in 0..2 {
        let p: f64 = r5_pT_data[i][1];
        let h: f64 = r5_pT_data[i][3];
        assert_approx_eq!(r5_pT_data[i][0] - 273.15, ph(p, h, OT));
        assert_approx_eq!(r5_pT_data[i][5], ph(p, h, OS));
    }
}

#[test]
fn test_region4_ph() {
    let mut h: f64 = 0.0;
    let mut h1: f64 = 0.0;
    let mut h2: f64 = 0.0;
    let mut x: f64 = 0.0;
    for i in 0..3 {
        let p: f64 = r4_sat_Tp[i][1];
        x = 0.35;
        h = px(p, x, OH);
        assert_approx_eq!(x, ph(p, h, OX));
    }

    for i in 0..3 {
        let p: f64 = r4_sat_pT[i][0];
        x = 0.55;
        h = px(p, x, OH);
        assert_approx_eq!(x, ph(p, h, OX));
    }
}
