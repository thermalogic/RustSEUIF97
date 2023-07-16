#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use if97::common::*;
use if97::hs_reg2;
use if97::ph_reg2;
use if97::ps_reg2;
use if97::pt_reg2;

#[test]
fn test_pt() {
    for i in 0..2 {
        let p: f64 = r2_pT_data[i].p;
        let t: f64 = r2_pT_data[i].T - 273.15;
        assert_approx_eq!(r2_pT_data[i].v, pt_reg2(p, t, OV), 1.0e-6f64);
        assert_approx_eq!(r2_pT_data[i].h, pt_reg2(p, t, OH), 1.0e-5f64);
        assert_approx_eq!(r2_pT_data[i].s, pt_reg2(p, t, OS), 1.0e-6f64);
        assert_approx_eq!(r2_pT_data[i].u, pt_reg2(p, t, OU), 1.0e-5f64);
        assert_approx_eq!(r2_pT_data[i].cp, pt_reg2(p, t, OCP), 1.0e-6f64);
        assert_approx_eq!(r2_pT_data[i].w, pt_reg2(p, t, OW), 1.0e-6f64);
    }
}

#[test]
fn test_ph() {
    let mut p: f64 = 0.0;
    let mut h: f64 = 0.0;
    for i in 0..3 {
        p = r2_ph_2a[i][0];
        h = r2_ph_2a[i][1];
        assert_approx_eq!(r2_ph_2a[i][2], ph_reg2(p, h, OT) + 273.15, 1.0e-1f64);
        p = r2_ph_2b[i][0];
        h = r2_ph_2b[i][1];
        assert_approx_eq!(r2_ph_2b[i][2], ph_reg2(p, h, OT) + 273.15, 1.0e-1f64);
        p = r2_ph_2c[i][0];
        h = r2_ph_2c[i][1];
        assert_approx_eq!(r2_ph_2c[i][2], ph_reg2(p, h, OT) + 273.15, 1.0e-1f64);
    }
}

#[test]
fn test_ph_iter() {
    let mut h: f64 = 0.0;
    for i in 0..3 {
        let p: f64 = r2_pT_data[i].p;
        let t: f64 = r2_pT_data[i].T - 273.15;
        h = pt_reg2(p, t, OH);
        assert_approx_eq!(t, ph_reg2(p, h, OT), 1.0e-6f64);
    }
}

#[test]
fn test_ps() {
    let mut p: f64 = 0.0;
    let mut s: f64 = 0.0;
    for i in 0..3 {
        p = r2_ps_2a[i][0];
        s = r2_ps_2a[i][1];
        assert_approx_eq!(r2_ps_2a[i][2], ps_reg2(p, s, OT) + 273.15, 1.0e-2f64);
        p = r2_ps_2b[i][0];
        s = r2_ps_2b[i][1];
        assert_approx_eq!(r2_ps_2b[i][2], ps_reg2(p, s, OT) + 273.15, 1.0e-2f64);
        p = r2_ps_2c[i][0];
        s = r2_ps_2c[i][1];
        assert_approx_eq!(r2_ps_2c[i][2], ps_reg2(p, s, OT) + 273.15, 1.0e-2f64);
    }
}

#[test]
fn test_ps_iter() {
    let mut s: f64 = 0.0;
    for i in 0..3 {
        let p: f64 = r2_pT_data[i].p;
        let t: f64 = r2_pT_data[i].T - 273.15;
        s = pt_reg2(p, t, OS);
        assert_approx_eq!(t, ps_reg2(p, s, OT), 1.0e-6f64);
    }
}

#[test]
fn test_hs() {
    let mut h: f64 = 0.0;
    let mut s: f64 = 0.0;
    for i in 0..3 {
        h = r2_hs_reg2a[i][0];
        s = r2_hs_reg2a[i][1];
        assert_approx_eq!(r2_hs_reg2a[i][2], hs_reg2(h, s, OP), 1.0e-2f64);
        h = r2_hs_reg2b[i][0];
        s = r2_hs_reg2b[i][1];
        assert_approx_eq!(r2_hs_reg2b[i][2], hs_reg2(h, s, OP), 1.0e-2f64);
        h = r2_hs_reg2c[i][0];
        s = r2_hs_reg2c[i][1];
        assert_approx_eq!(r2_hs_reg2c[i][2], hs_reg2(h, s, OP), 1.0e-2f64);
    }
}

#[test]
fn test_hs_iter() {
    let mut s: f64 = 0.0;
    let mut h: f64 = 0.0;
    for i in 0..3 {
        let p: f64 = r2_pT_data[i].p;
        let t: f64 = r2_pT_data[i].T - 273.15;
        h = pt_reg2(p, t, OH);
        s = pt_reg2(p, t, OS);
        assert_approx_eq!(p, hs_reg2(h, s, OP), 1.0e-6f64);
        assert_approx_eq!(t, hs_reg2(h, s, OT), 1.0e-6f64);
    }
}
