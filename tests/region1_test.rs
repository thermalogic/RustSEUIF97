
// importing common module.
#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;

use if97::pt_reg1;
use if97::ph_reg1;
use if97::ps_reg1;
use if97::hs_reg1;
use if97::common::*;

// Table 5，Page 9: p,T,v,h,u,s,cp,w
const p: [f64; 3] = [3.0, 80.0, 3.0];
const T: [f64; 3] = [300.0, 300.0, 500.0];

const v: [f64; 3] = [0.100215168e-2, 0.971180894e-3, 0.120241800e-2];
const h: [f64; 3] = [0.115331273e+3, 0.184142828e+3, 0.975542239e+3];
const s: [f64; 3] = [0.392294792, 0.368563852, 0.258041912e1];
const u: [f64; 3] = [0.112324818e3, 0.106448356e3, 0.971934985e3];
const cp: [f64; 3] = [0.417301218e1, 0.401008987e1, 0.465580682e1];
//const cv: [f64; 3] = [0.412120160e1, 0.391736606e1, 0.322139223e1];
const w: [f64; 3] = [0.150773921e4, 0.163469054e4, 0.124071337e4];

#[test]
fn test_pT() {
    // using common code.
 
    for i in 0..3 {
        let t:f64 = T[i]- 273.15;
        assert_approx_eq!(v[i], if97::pt_reg1(p[i],t, OV), 1.0e-9f64);
        assert_approx_eq!(h[i], if97::pt_reg1(p[i],t, OH), 1.0e-6f64);
        assert_approx_eq!(s[i], if97::pt_reg1(p[i],t, OS), 1.0e-9f64);
        assert_approx_eq!(u[i], if97::pt_reg1(p[i],t, OU), 1.0e-6f64);
        assert_approx_eq!(cp[i], if97::pt_reg1(p[i], t, OCP), 1.0e-8f64);
        //assert_approx_eq!(cv[i], if97::pt(p[i],t, OCV), 1.0e-8f64);
        assert_approx_eq!(w[i], if97::pt_reg1(p[i], t, OW), 1.0e-5f64);
    }
}

#[test]
fn test_ph() {
    const p: [f64; 3] = [3.0, 80.0, 3.0];
    const h: [f64; 3] = [0.115331273e+3, 0.184142828e+3, 0.975542239e+3];
    const T: [f64; 3] = [300.0, 300.0, 500.0];
    for i in 0..3 {
        assert_approx_eq!(T[i] - 273.15, ph_reg1(p[i], h[i], OT), 1.0e-6f64);
        // real diff: `0.0178
    }
}



#[test]
fn test_ps() {
    const p: [f64; 3] = [3.0, 80.0, 3.0];
    const s: [f64; 3] = [0.392294792, 0.368563852, 0.258041912e1];
    const T: [f64; 3] = [300.0, 300.0, 500.0];
    for i in 0..=2 {
        assert_approx_eq!(T[i] - 273.15, ps_reg1(p[i], s[i], OT), 1.0e-6f64);
    }
}

#[test]
fn test_hs() {
    const h: [f64; 3] = [0.115331273e+3, 0.184142828e+3, 0.975542239e+3];
    const s: [f64; 3] = [0.392294792, 0.368563852, 0.258041912e1];
    const p: [f64; 3] = [3.0, 80.0, 3.0];
    const T: [f64; 3] = [300.0, 300.0, 500.0];
    for i in 0..=2 {
        assert_approx_eq!(T[i] - 273.15, hs_reg1(h[i], s[i], OT), 1.0e-6f64); //real diff: `0.13
        assert_approx_eq!(p[i], hs_reg1(h[i], s[i], OP), 1.0e-6f64); //  real diff: `0.181
    }
}

#[test]
fn test_ph_iter() {
    let mut h1: f64 = 0.0;
    let mut v1: f64 = 0.0;
    for i in 0..2 {
        h1 = if97::pt(p[i], T[i] - 273.15, OH);
        v1 = if97::pt(p[i], T[i] - 273.15, OV);
        assert_approx_eq!(v1, ph_reg1(p[i], h1, OV), 1.0e-6f64);
    }

    for i in 0..=2 {
        h1 = if97::pt(p[i], T[i] - 273.15, OH);
        assert_approx_eq!(T[i] - 273.15, ph_reg1(p[i], h1, OT), 1.0e-6f64);
    }
}

#[test]
fn test_ps_iter() {
    for i in 0..=2 {
        assert_approx_eq!(T[i] - 273.15, ps_reg1(p[i], s[i], OT), 1.0e-6f64);
    }
}

#[test]
fn test_hs_iter() {
    for i in 0..=2 {
        assert_approx_eq!(T[i] - 273.15, hs_reg1(h[i], s[i], OT));
        assert_approx_eq!(p[i], hs_reg1(h[i], s[i], OP));
    }
}
