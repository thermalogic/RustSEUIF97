//! Region 5 - Backward: (p,h)  (p,s) (h,s)
//!
//! Backward of Region 5 using the secant method
//!

use crate::algo::root::I_MAX;
use crate::algo::*;
use crate::common::constant::*;
use crate::r5::region5_pT::*;

pub fn ph2T_reg5(p: f64, h: f64) -> f64 {
    let T1: f64 = T_MIN5;
    let T2: f64 = T_MAX5;
    let func = |T: f64| -> f64 {
        h - pT2h_reg5(p, T)
    };
    bisection(T1, T2, func,  20000, 1.0e-8, 1.0e-6) 
}

pub fn ps2T_reg5(p: f64, s: f64) -> f64 {
    let T1: f64 = T_MIN5;
    let T2: f64 = T_MAX5;
    let func = |T: f64| -> f64 {
        s - pT2s_reg5(p, T)
    };
    bisection(T1, T2, func,  20000, 1.0e-8, 1.0e-6) 
}

pub fn hs2p_reg5(h: f64, s: f64) -> f64 {
    let p1: f64 = P_MIN5;
    let p2: f64 = P_MAX5;
    let func = |p: f64| -> f64 {
        let T: f64 = ph2T_reg5(p, h);
        s-pT2s_reg5(p, T)
    };
    bisection(p1, p2, func,  20000, 1.0e-8, 1.0e-6) 
}
