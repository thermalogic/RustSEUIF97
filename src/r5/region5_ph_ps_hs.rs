//! Region 5 - Backward: (p,h)  (p,s) (h,s)
//!
//！ Backward of Region 5 using the secant method
//！

use crate::algo::root::I_MAX;
use crate::algo::*;
use crate::common::constant::*;
use crate::r5::region5_pT::*;

pub fn ph2T_reg5(p: f64, h: f64) -> f64 {
    let mut T: f64 = -1000.0;
    let mut T1: f64 = 0.5 * (2273.15 + 1073.15);
    let f1: f64 = h - pT2h_reg5(p, T1);
    if f1.abs() > ESP {
        let mut T2: f64 = -1000.0;
        if f1 > 0.0 {
            T2 = (1.0 + f1 / h) * T1;
        } else {
            T2 = (1.0 - f1 / h) * T1;
        }

        let f2: f64 = h - pT2h_reg5(p, T2);

        T = rtsec2(pT2h_reg5, p, h, T1, T2, f1, f2, ESP, I_MAX);
    } else {
        T = T1;
    }

    if T < T_MIN5 {
        T = T_MIN5;
    } else if T > T_MAX5 {
        T = T_MAX5;
    }

    return T;
}

pub fn ps2T_reg5(p: f64, s: f64) -> f64 {
    let mut T: f64 = -1000.0;
    let T1: f64 = 0.5 * (2273.15 + 1073.15); // Get initial value
    let f1: f64 = s - pT2s_reg5(p, T1);

    if f1.abs() > ESP {
        let mut T2: f64 = -1000.0;
        if f1 > 0.0 {
            T2 = (1.0 + f1 / s) * T1;
        } else {
            T2 = (1.0 - f1 / s) * T1;
        }

        let f2: f64 = s - pT2s_reg5(p, T2);
        T = rtsec2(pT2s_reg5, p, s, T1, T2, f1, f2, ESP, I_MAX);
    } else {
        T = T1;
    }

    if T < T_MIN5 {
        T = T_MIN5;
    } else if T > T_MAX5 {
        T = T_MAX5;
    }

    T
}

/// helper for hs2preg5
fn ph2s_reg5(p: f64, h: f64) -> f64 {
    let T: f64 = ph2T_reg5(p, h);
    return pT2s_reg5(p, T);
}

pub fn hs2p_reg5(h: f64, s: f64) -> f64 {
    let mut p: f64 = -1000.0;

    let p1: f64 = 0.5 * (P_MIN5 + P_MAX5);
    let f1: f64 = s - ph2s_reg5(p1, h);
    if f1.abs() > ESP {
        let mut p2: f64 = 0.0;
        if f1 > 0.0 {
            p2 = (1.0 + f1 / s) * p1;
        } else {
            p2 = (1.0 - f1 / s) * p1;
        }

        let f2: f64 = s - ph2s_reg5(p2, h);
        p = rtsec1(ph2s_reg5, h, s, p1, p2, f1, f2, ESP, I_MAX);
    } else {
        p = p1;
    }

    if p < P_MIN5 {
        p = P_MIN5;
    } else if p > P_MAX5 {
        p = P_MAX5;
    }
    p
}
