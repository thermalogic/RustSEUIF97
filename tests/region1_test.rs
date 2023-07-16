#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use if97::common::*;
use if97::hs_reg1;
use if97::ph_reg1;
use if97::ps_reg1;
use if97::pt_reg1;

#[test]
fn test_pT() {
    for i in 0..3 {
        let t: f64 = r1_pT_data[i].T - 273.15;
        let p: f64 = r1_pT_data[i].p;
        assert_approx_eq!(r1_pT_data[i].v, if97::pt_reg1(p, t, OV), 1.0e-9f64);
        assert_approx_eq!(r1_pT_data[i].h, if97::pt_reg1(p, t, OH), 1.0e-6f64);
        assert_approx_eq!(r1_pT_data[i].s, if97::pt_reg1(p, t, OS), 1.0e-9f64);
        assert_approx_eq!(r1_pT_data[i].u, if97::pt_reg1(p, t, OU), 1.0e-6f64);
        assert_approx_eq!(r1_pT_data[i].cp, if97::pt_reg1(p, t, OCP), 1.0e-8f64);
        assert_approx_eq!(r1_pT_data[i].w, if97::pt_reg1(p, t, OW), 1.0e-5f64);
    }
}

#[test]
fn test_ph() {
    for i in 0..3 {
        assert_approx_eq!(
            r1_phT[i][2] - 273.15,
            ph_reg1(r1_phT[i][0], r1_phT[i][1], OT),
            1.0e-6f64
        ); // real diff: `0.0178
    }
}

#[test]
fn test_ps() {
    for i in 0..=2 {
        assert_approx_eq!(
            r1_psT[i][2] - 273.15,
            ps_reg1(r1_psT[i][0], r1_psT[i][1], OT),
            1.0e-6f64
        );
    }
}

#[test]
pub fn test_hs() {
    for i in 0..=2 {
        assert_approx_eq!(
            r1_hspT[i][3] - 273.15,
            hs_reg1(r1_hspT[i][0], r1_hspT[i][1], OT),
            1.0e-6f64
        ); //real diff: `0.13
        assert_approx_eq!(
            r1_hspT[i][2],
            hs_reg1(r1_hspT[i][0], r1_hspT[i][1], OP),
            1.0e-6f64
        ); //  real diff: `0.181
    }
}

#[test]
fn test_ph_iter() {
    let mut h1: f64 = 0.0;
    let mut v1: f64 = 0.0;
    for i in 0..2 {
        h1 = if97::pt(r1_pT_data[i].p, r1_pT_data[i].T - 273.15, OH);
        v1 = if97::pt(r1_pT_data[i].p, r1_pT_data[i].T - 273.15, OV);
        assert_approx_eq!(v1, ph_reg1(r1_pT_data[i].p, h1, OV), 1.0e-6f64);
    }

    for i in 0..=2 {
        h1 = if97::pt(r1_pT_data[i].p, r1_pT_data[i].T - 273.15, OH);
        assert_approx_eq!(
            r1_pT_data[i].T - 273.15,
            ph_reg1(r1_pT_data[i].p, h1, OT),
            1.0e-6f64
        );
    }
}

#[test]
fn test_ps_iter() {
    for i in 0..=2 {
        assert_approx_eq!(
            r1_pT_data[i].T - 273.15,
            ps_reg1(r1_pT_data[i].p, r1_pT_data[i].s, OT),
            1.0e-6f64
        );
    }
}

#[test]
fn test_hs_iter() {
    for i in 0..=2 {
        assert_approx_eq!(
            r1_pT_data[i].T - 273.15,
            hs_reg1(r1_pT_data[i].h, r1_pT_data[i].s, OT)
        );
        assert_approx_eq!(
            r1_pT_data[i].p,
            hs_reg1(r1_pT_data[i].h, r1_pT_data[i].s, OP)
        );
    }
}
