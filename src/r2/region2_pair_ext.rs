//!   IAPWS-IF97 Region2: The extended input pair
//!       (p,v)
//!       (t,h),(t,s),(t,v)
use crate::algo::root::*;
use crate::common::boundaries::*;
use crate::common::constant::*;
use crate::r2::region2_pT::*;
use crate::r4::region4_sat_pT::*;

/// the helper for the extended input pair
fn p2Tmin_reg2(p: f64) -> f64 {
    if p > 0.0 && p < 0.000611213 {
        T_MIN2
    } else if p <= p_saturation(623.15) {
        T_saturation(p)
    } else {
        B23_p2T(p)
    }
}

fn T2pmax_reg2(T: f64) -> f64 {
    if T >= 273.15 && T <= 623.15 {
        p_saturation(T)
    } else if T <= 863.15 {
        B23_T2p(T)
    } else if T <= 1073.15 {
        P_MAX2
    } else {
        unreachable!()  
    }
 }
 
/// Region 2  (p,v)->T using the bisection method
///      p: pressure  MPa
//       v: specific volume m^3/kg
///      T: temperature  K
pub fn pv2T_reg2(p: f64, v: f64) -> f64 {
    let Tmin2: f64 = p2Tmin_reg2(p);
    let T1: f64 = Tmin2;
    let T2:f64 = T_MAX2;
    let func = |T: f64| -> f64 {
        (v - pT2v_reg2(p, T)) / v
    };
    bisection(T1, T2, func,  20000, 1.0e-8, 1.0e-6) 
}

/// Region 2  (T,v)->p using the bisection method
///      T: temperature  K
///      p: pressure  MPa
//       v: specific volume m^3/kg
pub fn Tv2p_reg2(T: f64, v: f64) -> f64 {
    let mut p1: f64 = P_MIN2;
    let mut p2: f64 = P_MAX2;
    let func = |p: f64| -> f64 {
        (v - pT2v_reg2(p, T)) / v
    };
    bisection(p1, p2, func,  20000, 1.0e-8, 1.0e-6) 
}

/// Region 2  (T,h)->p using the bisection method
///      T: temperature  K
///      h: specific enthalpy kJ/kg
///      p: pressure  MPa
pub fn Th2p_reg2(T: f64, h: f64) -> f64 {
    let mut p1: f64 = P_MIN2;
    let mut p2: f64 = P_MAX2;
    let func = |p: f64| -> f64 {
        h - pT2h_reg2(p, T)
    };
    bisection(p1, p2, func,  20000, 1.0e-8, 1.0e-6) 
}

/// Region 2  (T,s)->p using the bisection method
///  * T: temperature  K
///  * s: specific entropy  kJ/(kg K)
///  * p: pressure  MPa
pub fn Ts2p_reg2(T: f64, s: f64) -> f64 {
    let mut p1: f64 = P_MIN2;
    let mut p2: f64 =T2pmax_reg2(T);
    let func = |p: f64| -> f64 {
        s - pT2s_reg2(p, T)
    };
    bisection(p1, p2, func,  20000, 1.0e-8, 1.0e-6) 
}
