//! IAPWS-IF97 Region5: The extended input pair
//!  *  (p,v)
//!  *  (t,v),(t,h),(t,s)
use crate::algo::root::*;
use crate::common::boundaries::*;
use crate::common::constant::*;
use crate::r5::region5_pT::*;

/// Region 5: (p,v)-> T using the Brant's method
pub fn pv2T_reg5(p: f64, v: f64) -> f64 {
    let T1: f64 = T_MIN5;
    let T2: f64 = T_MAX5;
   // let func = |T: f64| -> f64 {
   //     (v - pT2v_reg5(p, T)) / v
    //};
    //bisection(T1, T2, func,  20000, 1.0e-8, 1.0e-6)   
    zbrent(pT2v_reg5, p, v, T1, T2, 1, ESP, I_MAX)  
}

/// Region 5: (T,v)-> p using the bisection method
pub fn Tv2p_reg5(T: f64, v: f64) -> f64 {
    let p1: f64 = P_MIN5;
    let p2: f64 = P_MAX5;
    let func = |p: f64| -> f64 {
        (v - pT2v_reg5(p, T)) / v
    };
    bisection(p1, p2, func,  20000, 1.0e-8, 1.0e-6)      
    //zbrent(pT2v_reg5, T, v, p1, p2, 2, ESP, I_MAX)  
}

/// Region 5: (T,h)-> p using the secant method
pub fn Th2p_reg5(T: f64, h: f64) -> f64 {
    let p1: f64 = P_MIN5;
    let p2: f64 = P_MAX5; //
    rtsec(pT2h_reg5, T, h, p1, p2, 2, ESP, I_MAX)
}

/// Region 5: (T,s)-> s using the secant method
pub fn Ts2p_reg5(T: f64, s: f64) -> f64 {
    let p1: f64 = P_MIN5;
    let p2: f64 = P_MAX5;
    rtsec(pT2s_reg5, T, s, p1, p2, 2, ESP, I_MAX)
}
