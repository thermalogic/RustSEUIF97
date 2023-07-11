#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;
mod common;
use common::*;
use if97::common::*;
use if97::hs2p_reg3;
use if97::pT2v_reg3;
use if97::ph2T_reg3;
use if97::ph2v_reg3;
use if97::ps2T_reg3;
use if97::ps2v_reg3;
use if97::Td_reg3;

#[test]
fn test_Td() {
    for i in 0..3 {
        let T: f64 = r3_Td[i][0];
        let d: f64 = r3_Td[i][1];
        assert_approx_eq!(r3_Td[i][2], Td_reg3(T, d, OP), 1.0e-5f64);
        assert_approx_eq!(r3_Td[i][3], Td_reg3(T, d, OH), 1.0e-5f64);
        assert_approx_eq!(r3_Td[i][4], Td_reg3(T, d, OU), 1.0e-5f64);
        assert_approx_eq!(r3_Td[i][5], Td_reg3(T, d, OS), 1.0e-5f64);
        assert_approx_eq!(r3_Td[i][6], Td_reg3(T, d, OCP), 1.0e-5f64);
        assert_approx_eq!(r3_Td[i][7], Td_reg3(T, d, OW), 1.0e-5f64);
    }
}

#[test]
fn test_pT() {
    for i in 0..52 {
        let p: f64 = r3_v_pT[i][1];
        let T: f64 = r3_v_pT[i][2];
        assert_approx_eq!(r3_v_pT[i][0], pT2v_reg3(p, T));
    }
}

#[test]
fn test_ph() {
    for i in 0..3 {
        let mut p: f64 = r3_phTv_3a[i][0];
        let mut h: f64 = r3_phTv_3a[i][1];
        assert_approx_eq!(r3_phTv_3a[i][2], ph2T_reg3(p, h));
        assert_approx_eq!(r3_phTv_3a[i][3], ph2v_reg3(p, h));
        p = r3_phTv_3b[i][0];
        h = r3_phTv_3b[i][1];
        assert_approx_eq!(r3_phTv_3b[i][2], ph2T_reg3(p, h));
        assert_approx_eq!(r3_phTv_3b[i][3], ph2v_reg3(p, h));
    }
}

#[test]
fn test_ps() {
    for i in 0..3 {
        let mut p: f64 = r3_psTv_3a[i][0];
        let mut s: f64 = r3_psTv_3a[i][1];
        assert_approx_eq!(r3_psTv_3a[i][2], ps2T_reg3(p, s));
        assert_approx_eq!(r3_psTv_3a[i][3], ps2v_reg3(p, s));
        p = r3_psTv_3b[i][0];
        s = r3_psTv_3b[i][1];
        assert_approx_eq!(r3_psTv_3b[i][2], ps2T_reg3(p, s));
        assert_approx_eq!(r3_psTv_3b[i][3], ps2v_reg3(p, s));
    }
}

#[test]
fn test_hs() {
    let mut h: f64 = 0.0;
    let mut s: f64 = 0.0;
    for i in 0..3 {
        h = r3_hsp_3a[i][0];
        s = r3_hsp_3a[i][1];
        assert_approx_eq!(r3_hsp_3a[i][2], hs2p_reg3(h, s));
        h = r3_hsp_3b[i][0];
        s = r3_hsp_3b[i][1];
        assert_approx_eq!(r3_hsp_3b[i][2], hs2p_reg3(h, s));
    }
}
