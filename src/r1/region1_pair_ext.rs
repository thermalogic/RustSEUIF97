//！IAPWS-IF97 Region1: The extended input pair
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
    let mut T1: f64 = T_MIN1;
    let mut v1: f64 = pT2v_reg1(p, T1);
    let mut f1: f64 = v - pT2v_reg1(p, T1);

    let mut T2: f64 = T_MAX1;
    if (p >= 16.5291643) && (p <= 100.0) {
        T2 = T_MAX1;
    } else {
        T2 = T_saturation(p);
    };
    let mut v2: f64 = pT2v_reg1(p, T2);
    let mut f: f64 = v - v2;
    if (v2 - v1) != 0.0 {
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

    let mut r_error = 0.0; //relative error
    r_error = (v - pT2v_reg1(p, T)) / v;
    if T >= T_MIN1 && T <= T_MAX1 && r_error.abs() < 1.0e-8 {
        return T;
    };
    if T < T_MIN1 {
        T = T_MIN1;
    };
    if T > T_MAX1 {
        T = T_MAX1;
    }

    // Region 1 : the diff volume is the very small when the diff T is large
    // so, we need to adjust the T
    let V_ESPION: f64 = 1.0e-8;
    let MAX_STEPS: i32 = 1000000;
    f = v - pT2v_reg1(p, T);
    r_error = f / v;
    let mut success: bool = false;
    let mut steps: i32 = 0;
    if r_error > 0.0
    // T^ -> V^, the  v of T is samller than the real value,T is  small
    {
        while !success {
            T += 0.01; //the T is smaller than the real value, so  +T
            if T > T_MAX1 {
                return T;
            };
            f = v - pT2v_reg1(p, T);
            r_error = f / v;
            if r_error.abs() < V_ESPION {
                return T;
            } else {
                steps += 1;
                if steps > MAX_STEPS {
                    success = true;
                }
            }
        }
    }

    if r_error < 0.0
    // T^ -> V^, the  v of T is biggerr than the real value,T is  bigger
    {
        success = false;
        steps = 0;
        while !success {
            T -= 0.001; //the T is biggerr than the real value,so  -T
            if T < T_MIN1 {
                return T_MIN1;
            };
            r_error = (v - pT2v_reg1(p, T)) / v;
            if r_error.abs() < V_ESPION {
                return T;
            } else {
                steps += 1;
                if steps > MAX_STEPS {
                    success = true;
                }
            };
        }
    };
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

    if p > P_MAX1 {
        p = P_MAX1;
    }
    if p < pmin1 {
        p = pmin1;
    }
    p
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
    if p > P_MAX1 {
        p = P_MAX1;
    } else if p < pmin1 {
        p = pmin1;
    }
    p
}
