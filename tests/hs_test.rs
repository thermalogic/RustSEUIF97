#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use if97::common::*;
use if97::hs;
use if97::pt;
use if97::px;

#[test]
fn test_region1_hs() {
    for i in 0..=2 {
        assert_approx_eq!(
            r1_hspT[i][3] - 273.15,
            hs(r1_hspT[i][0], r1_hspT[i][1], OT),
            1.0e-6f64
        ); //real diff: `0.13
        assert_approx_eq!(
            r1_hspT[i][2],
            hs(r1_hspT[i][0], r1_hspT[i][1], OP),
            1.0e-6f64
        ); //  real diff: `0.181
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
        assert_approx_eq!(p, hs(h, s, OP), 1.0e-6f64);
        assert_approx_eq!(t, hs(h, s, OT), 1.0e-6f64);
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
    const hsT: [[f64; 3]; 3] = [
        [1800.0, 5.3, 346.8475498],
        [2400.0, 6.0, 425.1373305],
        [2500.0, 5.5, 522.5579013],
    ];

    let mut h: f64 = 0.0;
    let mut h1: f64 = 0.0;
    let mut h2: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut s1: f64 = 0.0;
    let mut s2: f64 = 0.0;
    let mut x: f64 = 0.0;

    let mut p: f64 = 0.0;
    let mut T: f64 = 0.0;
    for i in 0..2 {
        h = hsT[i][0];
        s = hsT[i][1];
        assert_approx_eq!(hsT[i][2] - 273.15, hs(h, s, OT));
    }

    for i in 0..2 {
        p = r4_sat_Tp[i][1];
        T = r4_sat_Tp[i][0];
        s1 = px(p, 0.0, OS);
        s2 = px(p, 1.0, OS);
        h1 = px(p, 0.0, OH);
        h2 = px(p, 1.0, OH);
        x = 0.60;
        s = s1 + x * (s2 - s1);
        h = h1 + x * (h2 - h1);
        assert_approx_eq!(T - 273.15, hs(h, s, OT), 1.0e-3f64);
        assert_approx_eq!(p, hs(h, s, OP), 1.0e-3f64);
        assert_approx_eq!(x, hs(h, s, OX), 1.0e-3f64);
    }

    for i in 0..2 {
        p = r4_sat_pT[i][0];
        T = r4_sat_pT[i][1];
        s1 = px(p, 0.0, OS);
        s2 = px(p, 1.0, OS);
        h1 = px(p, 0.0, OH);
        h2 = px(p, 1.0, OH);
        x = 0.05;
        s = s1 + x * (s2 - s1);
        h = h1 + x * (h2 - h1);
        assert_approx_eq!(T - 273.15, hs(h, s, OT), 1.0e-2f64);
        assert_approx_eq!(p, hs(h, s, OP), 1.0e-3f64);
        assert_approx_eq!(x, hs(h, s, OX), 1.0e-3f64);
    }
}
