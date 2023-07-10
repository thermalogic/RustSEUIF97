//cargo test --test region4_test
#![allow(warnings)]
// importing common module.
mod common;
use assert_approx_eq::assert_approx_eq;

use if97::common::*;
use if97::p_sat;
use if97::ph;
use if97::ps;
use if97::hs;
use if97:: hs_reg4;
use if97::region4_T_hs::*;
use if97::region4_pTx::*;
use if97::t_sat;

// saturation T,p
const tab35: [[f64; 2]; 3] = [
    [300.0, 0.00353658941],
    [500.0, 2.63889776],
    [600.0, 12.3443146],
];

// saturation p,T
const tab36: [[f64; 2]; 3] = [[0.1, 372.755919], [1.0, 453.035632], [10.0, 584.149488]];

#[test]
fn test_saturation_pt() {
    //  saturation p(t),t(p)
    for i in 0..3 {
        let t: f64 = tab35[i][0] - 273.15;
        assert_approx_eq!(tab35[i][1], p_sat(t));
    }

    for i in 0..3 {
        let p: f64 = tab36[i][0];
        assert_approx_eq!(tab36[i][1], t_sat(p) + 273.15);
    }
}

#[test]
fn test_region4_ph() {
    let mut h: f64 = 0.0;
    let mut h1: f64 = 0.0;
    let mut h2: f64 = 0.0;
    let mut x: f64 = 0.0;
    for i in 0..3 {
        let p: f64 = tab35[i][1];
        h1 = p2sat_water(p, OH);
        h2 = p2sat_steam(p, OH);
        x = 0.35;
        h = h1 + x * (h2 - h1);
        assert_approx_eq!(x, ph(p, h, OX));
    }

    for i in 0..3 {
        let p: f64 = tab36[i][0];
        h1 = p2sat_water(p, OH);
        h2 = p2sat_steam(p, OH);
        x = 0.55;
        h = h1 + x * (h2 - h1);
        assert_approx_eq!(x, ph(p, h, OX));
    }
}

#[test]
fn test_region4_ps() {
    let mut s: f64 = 0.0;
    let mut s1: f64 = 0.0;
    let mut s2: f64 = 0.0;
    let mut x: f64 = 0.0;
    for i in 0..3 {
        let p: f64 = tab35[i][1];
        s1 = p2sat_water(p, OS);
        s2 = p2sat_steam(p, OS);
        x = 0.35;
        s = s1 + x * (s2 - s1);
        assert_approx_eq!(x, ps(p, s, OX));
    }

    for i in 0..3 {
        let p: f64 = tab36[i][0];
        s1 = p2sat_water(p, OS);
        s2 = p2sat_steam(p, OS);
        x = 0.55;
        s = s1 + x * (s2 - s1);
        assert_approx_eq!(x, ps(p, s, OX));
    }
}

#[test]
fn test_region4_hs()
{
   const hsT:[[f64;3];3] =[
    [1800.0, 5.3, 346.8475498],
    [2400.0, 6.0, 425.1373305],
    [2500.0, 5.5, 522.5579013]];

let mut h:f64=0.0;
let mut h1:f64=0.0;
let mut h2:f64=0.0;
let mut s:f64=0.0;
let mut s1:f64=0.0;
let mut s2:f64=0.0;
let mut x:f64=0.0;

let mut p:f64=0.0;
let mut T:f64=0.0;
for i in 0..2
{
  h = hsT[i][0];
  s = hsT[i][1];
  assert_approx_eq!(hsT[i][2],hs2T_reg4(h, s));
}
}

