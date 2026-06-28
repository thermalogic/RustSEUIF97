#![allow(warnings)]
///  cargo test --test hs_test
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use seuif97::*;

#[test]
fn test_region1_hs() {
    for i in 0..=2 {
        assert_approx_eq!(r1_hspT[i][3] - 273.15, hs(r1_hspT[i][0], r1_hspT[i][1], OT), 1.0e-1f64);
        assert_approx_eq!(r1_hspT[i][2], hs(r1_hspT[i][0], r1_hspT[i][1], OP), 1.0e-1f64);
    }
}

#[test]
fn test_region2_hs() {
    let mut s: f64 = 0.0;
    let mut h: f64 = 0.0;
    for i in 0..3 {
        let p: f64 = r2_pT_data[i].p;
        let t: f64 = r2_pT_data[i].T - 273.15;
        h = pt(p, t, OH);
        s = pt(p, t, OS);
        assert_approx_eq!(p, hs(h, s, OP), 1.0e-2f64);
        assert_approx_eq!(t, hs(h, s, OT), 1.0e-2f64);
    }
}

#[test]
fn test_region3_hs() {
    let mut h: f64 = 0.0;
    let mut s: f64 = 0.0;
    for i in 0..3 {
        h = r3_hsp_3a[i][0];
        s = r3_hsp_3a[i][1];
        assert_approx_eq!(r3_hsp_3a[i][2], hs(h, s, OP));
        h = r3_hsp_3b[i][0];
        s = r3_hsp_3b[i][1];
        assert_approx_eq!(r3_hsp_3b[i][2], hs(h, s, OP));
    }
}

#[test]
fn test_region5_hs() {
    for i in 0..2 {
        let h: f64 = r5_pT_data[i][3];
        let s: f64 = r5_pT_data[i][5];
        assert_approx_eq!(r5_pT_data[i][0] - 273.15, hs(h, s, OT), 1.0e-4f64);
        assert_approx_eq!(r5_pT_data[i][1], hs(h, s, OP), 1.0e-4f64);
        assert_approx_eq!(r5_pT_data[i][2], hs(h, s, OV), 1.0e-4f64);
        assert_approx_eq!(r5_pT_data[i][4], hs(h, s, OU), 1.0e-4f64);
        assert_approx_eq!(r5_pT_data[i][6], hs(h, s, OCP), 1.0e-4f64);
        assert_approx_eq!(r5_pT_data[i][7], hs(h, s, OW), 1.0e-4f64);
    }
}

#[test]
fn test_region4_hs() {
    let mut h: f64 = 0.0;
    let mut s: f64 = 0.0;
    for i in 0..2 {
        h = r4_hsT[i][0];
        s = r4_hsT[i][1];
        assert_approx_eq!(r4_hsT[i][2] - 273.15,hs(h, s, OT),1.0e-7f64);
    }
}

#[test]
fn test_region4_x_hs() 
{
    let mut T: f64 = 0.0;
    let mut p: f64 = 0.0;
    let mut h: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut x: f64 = 0.0;
    for i in 0..2 {
        p = r4_sat_Tp[i][1];
        T = r4_sat_Tp[i][0];
        x = 0.6;
        s = px(p, x, OS);
        h = px(p, x, OH);
        assert_approx_eq!(T - 273.15, hs(h, s, OT), 1.0e-1f64);
        assert_approx_eq!(p, hs(h, s, OP), 1.0e-1f64);
        assert_approx_eq!(x, hs(h, s, OX), 1.0e-1f64);
    }
    
//r4_sat_pT: [[f64; 2]; 3] = [[0.1, 372.755919], [1.0, 453.035632], [10.0, 584.149488]];
   for i in 0..0 {
        p = r4_sat_pT[i][0];
        T = r4_sat_pT[i][1];
        x = 0.02;
        s = px(p, x, OS);
        h = px(p, x, OH);
        assert_approx_eq!(T - 273.15, hs(h, s, OT), 1.0e-2f64);
        assert_approx_eq!(p, hs(h, s, OP), 1.0e-3f64);
        assert_approx_eq!(x, hs(h, s, OX), 1.0e-3f64);
    }

    assert_approx_eq!(0.03653974055475902, hs(1800.0, 5.3, OP), 1.0e-4f64);
    assert_approx_eq!(0.270298590472939, hs(1500.0, 4.0, OP), 1.0e-4f64);
}
