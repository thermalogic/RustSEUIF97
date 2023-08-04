//！IAPWS-IF97 Region3: The extended input pair
//! *  (p,v)
//! *  (T,h),(T,s)

use crate::algo::root::*;
use crate::common::boundaries::*;
use crate::common::constant::*;
use crate::r3::*;

///   (p,v) ->T
pub fn pv2T_reg3(p: f64, v: f64) -> f64 {
    let mut T1: f64 = T_MIN3;
    let mut T2: f64 = B23_p2T(p);
    let d: f64 = 1.0 / v;
    let f1: f64 = p - Td2p_reg3(T1, d);
    let f2: f64 = p - Td2p_reg3(T2, d);
    rtsec1(Td2p_reg3, d, p, T1, T2, f1, f2, ESP, I_MAX)
}

///     (T,h) ->d
pub fn Th2d_reg3(T: f64, h: f64) -> f64 {
    let p1: f64 = B23_T2p(T);
    let d1: f64 = 1.0 / pT2v_reg3(p1, T);
    let p2: f64 = P_MAX3;
    let d2: f64 = 1.0 / pT2v_reg3(p2, T);
    let f1: f64 = h - Td2h_reg3(T, d1);
    let f2: f64 = h - Td2h_reg3(T, d2);
    rtsec2(Td2h_reg3, T, h, d1, d2, f1, f2, ESP, I_MAX)
}

// Region 3  (T,s)->d using the secant method
///  * T: temperature  K
///  * s: specific entropy  kJ/(kg K)
///  * d: density    kg/m^3
pub fn Ts2d_reg3(T: f64, s: f64) -> f64 {
    let d1: f64 = 100.0;
    let d2: f64 = 1.1 * d1;
    let ft: f64 = s - Td2s_reg3(T, d1);
    let f: f64 = s - Td2s_reg3(T, d2);
    rtsec2(Td2s_reg3, T, s, d1, d2, ft, f, ESP, I_MAX)
}
