//!IAPWS-IF97 Region1: The extended input pair
//!    (p,v), (t,h),(t,s),(t,v))
use crate::algo::root::*;
use crate::common::constant::*;
use crate::r1::region1_pT::*;
use crate::r4::region4_sat_pT::*;

/// Region 1  (p,v)->T using the secant method and refine adjust
/// * p: pressure  MPa
//  * v: specific volume m^3/kg
/// * T: temperature  K
pub fn pv2T_reg1(p: f64, v: f64) -> f64 {
    const V_ESPION: f64 = 1.0e-8;
    let mut T1: f64 = T_MIN1;
    let mut v1: f64 = pT2v_reg1(p, T1);
    let mut f1: f64 = v - v1;

    let mut T2: f64 = if p >= 16.5291643 && p <= 100.0 {
        T_MAX1
     } else {
        T_saturation(p)
    };
    
    let mut v2: f64 = pT2v_reg1(p, T2);
    let mut f: f64 = v - v2;
    if (v2 - v1).abs() >ESP {
        T2 = T1 + (T2 - T1) * (v - v1) / (v2 - v1);
    }
    v2 = pT2v_reg1(p, T2);
    f = v - v2;
    if f > 0.0 {
        T2 = T_MAX1;
        v2 = pT2v_reg1(p, T2);
        f = v - v2;
    }
    let mut T: f64 = rtsec2(pT2v_reg1, p, v, T1, T2, f1, f, 1.0e-3, 100);
    let mut r_error:f64 = (v - pT2v_reg1(p, T)) / v;//relative error
    if T >= T_MIN1 && T <= T_MAX1 && r_error.abs() < V_ESPION {
        return T;
    };
    T = T.clamp(T_MIN1, T_MAX1);
    // Region 1 : 
    //  the difference of volume is the very small when the difference T is large
    //  so, we need to adjust the T
    const MAX_STEPS: i32 = 1000000;
    const STEP_UP: f64 = 0.01;
    const STEP_DOWN: f64 = 0.001;
    // T^ -> V^ 
    //  r_error> 0.0, v is bigger than the real value, T need -
    //  r_error< 0.0 , v is smqller than the real value, T need +   
    let mut current_steps:i32=0;
    let step = if r_error > 0.0 { STEP_UP } else { -STEP_DOWN };
    let direction = if r_error > 0.0 { 1.0 } else { -1.0 };
    while r_error.abs() > V_ESPION && current_steps < MAX_STEPS {
          T += step;
          if T < T_MIN1 || T > T_MAX1 {
             return if direction > 0.0 { T_MAX1  } else { T_MIN1 };
          }
         f = v - pT2v_reg1(p, T);
         r_error = f / v;
         current_steps += 1;
   }
   INVALID_VALUE as f64
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
    let p: f64 = rtsec1(pT2v_reg1, T, v, p1, p2, f1, f, ESP, I_MAX);
    return p;
}

/// Region 1  (T,h)->p using the secant method
///  *  T: temperature  K
///  *  h: specific enthalpy kJ/kg
///  *  p: pressure  MPa
pub fn Th2p_reg1(T: f64, h: f64) -> f64 {
    let pmin1: f64 = p_saturation(T);
    let mut p1: f64 = pmin1;
    let mut p2: f64 = P_MAX1; // p1 + stepa
    let mut h1 = pT2h_reg1(p1, T);
    if (h - h1).abs() < ESP {
        return p1;
    };
    let mut h2: f64 = pT2h_reg1(p2, T);
    if (h - h2).abs() < ESP {
        return p2;
    }
    let f1: f64 = h - pT2h_reg1(p1, T);
    let f: f64 = h - pT2h_reg1(p2, T);
    let mut p: f64 = rtsec1(pT2h_reg1, T, h, p1, p2, f1, f, ESP, I_MAX);
    p.clamp(pmin1, P_MAX1)
}

/// Region 1  (T,s)->p using the secant method
///  * T: temperature  K
///  * s: specific entropy  kJ/(kg K)
///  * p: pressure  MPa
pub fn Ts2p_reg1(T: f64, s: f64) -> f64 {
    let pmin1: f64 = p_saturation(T);
    let mut p1: f64 = pmin1; //
    let mut s1: f64 = pT2s_reg1(p1, T);
    let mut f1: f64 = s - s1;
    let p2: f64 = P_MAX1;
    let s2: f64 = pT2s_reg1(p2, T);
    let f2: f64 = s - s2;
    p1 = p2 - (p2 - p1) * (s - s2) / (s1 - s2);
    if p1 < pmin1 {
        p1 = pmin1;
    }
    s1 = pT2s_reg1(p1, T);
    f1 = s - s1;
    if f1.abs() < ESP {
        return p1;
    }
    let mut p: f64 = rtsec1(pT2s_reg1, T, s, p1, p2, f1, f2, ESP, I_MAX);
    p.clamp(pmin1, P_MAX1)
}