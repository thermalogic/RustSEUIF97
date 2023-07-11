#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use if97::common::*;
use if97::ps;
use if97::pt;
use if97::region4_pTx::*;

#[test]
fn test_region1_ps() {
    for i in 0..=2 {
        assert_approx_eq!(
            r1_psT[i][2] - 273.15,
            ps(r1_psT[i][0], r1_psT[i][1], OT),
            1.0e-6f64
        );
    }
}

#[test]
fn test_region2_ps() {
    let mut p: f64 = 0.0;
    let mut s: f64 = 0.0;
    for i in 0..3 {
        p = r2_ps_2a[i][0];
        s = r2_ps_2a[i][1];
        assert_approx_eq!(r2_ps_2a[i][2], ps(p, s, OT) + 273.15, 1.0e-1f64);
        p = r2_ps_2b[i][0];
        s = r2_ps_2b[i][1];
        assert_approx_eq!(r2_ps_2b[i][2], ps(p, s, OT) + 273.15, 1.0e-1f64);
        p = r2_ps_2c[i][0];
        s = r2_ps_2c[i][1];
        assert_approx_eq!(r2_ps_2c[i][2], ps(p, s, OT) + 273.15, 1.0e-1f64);
    }
}

#[test]
fn test_region2_ps_iter() {
    let mut s: f64 = 0.0;
    for i in 0..3 {
        let p: f64 = r2_pT_data[i].p;
        let t: f64 = r2_pT_data[i].T - 273.15;
        s = pt(p, t, OS);
        assert_approx_eq!(t, ps(p, s, OT), 1.0e-6f64);
    }
}

#[test]
fn test_region3_ps() {
    for i in 0..3 {
        let mut p: f64 = r3_psTv_3a[i][0];
        let mut s: f64 = r3_psTv_3a[i][1];
        assert_approx_eq!(r3_psTv_3a[i][2] - 273.15, ps(p, s, OT));
        assert_approx_eq!(r3_psTv_3a[i][3], ps(p, s, OV));
        p = r3_psTv_3b[i][0];
        s = r3_psTv_3b[i][1];
        assert_approx_eq!(r3_psTv_3b[i][2] - 273.15, ps(p, s, OT));
        assert_approx_eq!(r3_psTv_3b[i][3], ps(p, s, OV));
    }
}

#[test]
fn test_region5_ps() {
    for i in 0..2 {
        let p: f64 = r5_pT_data[i][1];
        let s: f64 = r5_pT_data[i][5];
        assert_approx_eq!(r5_pT_data[i][0] - 273.15, ps(p, s, OT), 1.0e-4f64);
        assert_approx_eq!(r5_pT_data[i][3], ps(p, s, OH), 1.0e-4f64);
    }
}

#[test]
fn test_region4_ps() {
    let mut s: f64 = 0.0;
    let mut s1: f64 = 0.0;
    let mut s2: f64 = 0.0;
    let mut x: f64 = 0.0;
    for i in 0..3 {
        let p: f64 = r4_sat_Tp[i][1];
        s1 = p2sat_water(p, OS);
        s2 = p2sat_steam(p, OS);
        x = 0.35;
        s = s1 + x * (s2 - s1);
        assert_approx_eq!(x, ps(p, s, OX));
    }

    for i in 0..3 {
        let p: f64 = r4_sat_pT[i][0];
        s1 = p2sat_water(p, OS);
        s2 = p2sat_steam(p, OS);
        x = 0.55;
        s = s1 + x * (s2 - s1);
        assert_approx_eq!(x, ps(p, s, OX));
    }
}
