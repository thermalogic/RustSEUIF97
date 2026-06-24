//!IAPWS-IF97 Region1: The extended input pair
//!    (p,v), (t,h),(t,s),(t,v))
use crate::algo::root::*;
use crate::common::constant::*;
use crate::r1::region1_pT::*;
use crate::r4::region4_sat_pT::*;

/// Region 1  (p,v)->T using the bisection method
/// * p: pressure  MPa
//  * v: specific volume m^3/kg
/// * T: temperature  K
pub fn pv2T_reg1(p: f64, v: f64) -> f64 {
    let mut T1: f64 = T_MIN1;
    let mut T2: f64 = T_MAX1;
    if (p >= 16.5291643) && (p <= 100.0) {
        T2 = T_MAX1;
    } else {
        T2 = T_saturation(p);
    };
    let func = |T: f64| -> f64 {
        (v - pT2v_reg1(p, T)) / v
    };
    bisection(T1, T2, func,  20000, 1.0e-8, 1.0e-6)
}

/// Region 1  (T,v)->p using the secant method
/// * T: temperature  K
/// * v: specific volume m^3/kg
/// * p: pressure  MPa
pub fn Tv2p_reg1(T: f64, v: f64) -> f64 {
    let p1: f64 = 0.3 * (p_saturation(T) + P_MAX1);
    let p2: f64 = 1.05 * p1;
    let f1: f64 = v - pT2v_reg1(p1, T);
    let f: f64 = v - pT2v_reg1(p2, T);
    rtsec1(pT2v_reg1, T, v, p1, p2, f1, f, ESP, I_MAX)
}

/// Region 1  (T,h)->p using the bisection method
///  *  T: temperature  K
///  *  h: specific enthalpy kJ/kg
///  *  p: pressure  MPa
pub fn Th2p_reg1(T: f64, h: f64) -> f64 {
    let mut p1: f64 = p_saturation(T);
    let mut p2: f64 = P_MAX1;
    let func = |p: f64| -> f64 {
        h - pT2h_reg1(p, T)
    };
   bisection(p1, p2, func,  20000, 1.0e-8, 1.0e-6) 
}

/// Region 1  (T,s)->p using the bisection method
///  * T: temperature  K
///  * s: specific entropy  kJ/(kg K)
///  * p: pressure  MPa
pub fn Ts2p_reg1(T: f64, s: f64) -> f64 {
    let mut p1: f64 = p_saturation(T);
    let mut p2: f64 = P_MAX1;
    let func = |p: f64| -> f64 {
        s - pT2s_reg1(p, T)
    };
    bisection(p1, p2, func,  20000, 1.0e-8, 1.0e-6) 
}