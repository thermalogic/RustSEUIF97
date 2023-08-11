#![allow(warnings)]
///  cargo test --test pt_test
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use seuif97::*;

#[test]
fn test_regino1_pt() {
    for i in 0..3 {
        let t: f64 = r1_pT_data[i].T - 273.15;
        let p: f64 = r1_pT_data[i].p;
        assert_approx_eq!(r1_pT_data[i].v, pt(p, t, OV), 1.0e-9f64);
        assert_approx_eq!(r1_pT_data[i].h, pt(p, t, OH), 1.0e-6f64);
        assert_approx_eq!(r1_pT_data[i].s, pt(p, t, OS), 1.0e-9f64);
        assert_approx_eq!(r1_pT_data[i].u, pt(p, t, OU), 1.0e-6f64);
        assert_approx_eq!(r1_pT_data[i].cp, pt(p, t, OCP), 1.0e-8f64);
        assert_approx_eq!(r1_pT_data[i].w, pt(p, t, OW), 1.0e-5f64);
    }
}

#[test]
fn test_region2_pt() {
    for i in 0..2 {
        let p: f64 = r2_pT_data[i].p;
        let t: f64 = r2_pT_data[i].T - 273.15;
        assert_approx_eq!(r2_pT_data[i].v, pt(p, t, OV), 1.0e-6f64);
        assert_approx_eq!(r2_pT_data[i].h, pt(p, t, OH), 1.0e-5f64);
        assert_approx_eq!(r2_pT_data[i].s, pt(p, t, OS), 1.0e-6f64);
        assert_approx_eq!(r2_pT_data[i].u, pt(p, t, OU), 1.0e-5f64);
        assert_approx_eq!(r2_pT_data[i].cp, pt(p, t, OCP), 1.0e-6f64);
        assert_approx_eq!(r2_pT_data[i].w, pt(p, t, OW), 1.0e-6f64);
    }
}

#[test]
fn test_region3_pt() {
    for i in 0..52 {
        let p: f64 = r3_v_pT[i][1];
        let t: f64 = r3_v_pT[i][2] - 273.15;
        assert_approx_eq!(r3_v_pT[i][0], pt(p, t, OV));
    }
}

#[test]
fn test_region5_pt() {
    for i in 0..2 {
        let t: f64 = r5_pT_data[i][0] - 273.15;
        let p: f64 = r5_pT_data[i][1];
        assert_approx_eq!(r5_pT_data[i][2], pt(p, t, OV), 1.0e-6f64);
        assert_approx_eq!(r5_pT_data[i][3], pt(p, t, OH), 1.0e-4f64);
        assert_approx_eq!(r5_pT_data[i][4], pt(p, t, OU), 1.0e-4f64);
        assert_approx_eq!(r5_pT_data[i][5], pt(p, t, OS), 1.0e-5f64);
        assert_approx_eq!(r5_pT_data[i][6], pt(p, t, OCP), 1.0e-6f64);
        assert_approx_eq!(r5_pT_data[i][7], pt(p, t, OW), 1.0e-6f64);
    }
}
