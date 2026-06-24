//! IAPWS-IF97 Region5: The extended input pair
//!  *  (p,v)
//!  *  (t,v),(t,h),(t,s)
use crate::algo::root::*;
use crate::common::boundaries::*;
use crate::common::constant::*;
use crate::r5::region5_pT::*;

/// Region 5: (p,v)-> T using the bisection method
pub fn pv2T_reg5(p: f64, v: f64) -> f64 {
    let T1: f64 = T_MIN5;
    let T2: f64 = T_MAX5;
    let func = |T: f64| -> f64 {
        (v - pT2v_reg5(p, T)) / v
    };
    bisection(T1, T2, func,  20000, 1.0e-8, 1.0e-6)       
}

/// Region 5: (T,v)-> p using the bisection method
pub fn Tv2p_reg5(T: f64, v: f64) -> f64 {
    let p1: f64 = P_MIN5;
    let p2: f64 = P_MAX5;
    let func = |p: f64| -> f64 {
        (v - pT2v_reg5(p, T)) / v
    };
    bisection(p1, p2, func,  20000, 1.0e-8, 1.0e-6)      
}

/// Region 5: (T,h)-> p using the secant method
pub fn Th2p_reg5(T: f64, h: f64) -> f64 {
    let mut p1: f64 = P_MIN5;
    let mut h1: f64 = pT2h_reg5(p1, T);
    let mut p2: f64 = P_MAX5; //
    let mut h2: f64 = pT2h_reg5(p2, T);
    let f1: f64 = h - h1;
    let f2: f64 = h - h2;
    p1 = p2 - (p2 - p1) * ((h - h2) / (h1 - h2)).abs();
    rtsec1(pT2h_reg5, T, h, p1, p2, f1, f2, ESP, I_MAX)
}

/// Region 5: (T,s)-> s using the secant method
pub fn Ts2p_reg5(T: f64, s: f64) -> f64 {
    let mut p1: f64 = P_MIN5;
    let mut s1: f64 = pT2s_reg5(p1, T);
    let mut f1: f64 = s - s1;
    let p2: f64 = P_MAX5;
    let s2: f64 = pT2s_reg5(p2, T);
    let f2: f64 = s - s2;
    rtsec1(pT2s_reg5, T, s, p1, p2, f1, f2, ESP, I_MAX)
}
