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
    let mut Tmin = T_MIN2;
    if p > 0.0 && p < 0.000611213 {
        Tmin = T_MIN2;
    } else {
        if p >= 0.000611213 && p <= p_saturation(623.15) {
            Tmin = T_saturation(p);
        } else {
            Tmin = B23_p2T(p);
        }
    }
    Tmin
}

fn T2pmax_reg2(T: f64) -> f64 {
    let mut pmax: f64 = P_MAX2;
    if T >= 273.15 && T <= 623.15 {
        pmax = p_saturation(T);
    } else {
        if T > 623.15 && T <= 863.15 {
            pmax = B23_T2p(T);
        } else {
            if T > 863.15 && T <= 1073.15 {
                pmax = 100.0;
            };
        };
    };
    pmax
}

/// Region 2  (p,v)->T using the secant method and refine adjust
///      p: pressure  MPa
//       v: specific volume m^3/kg
///      T: temperature  K
pub fn pv2T_reg2(p: f64, v: f64) -> f64 {
    let Tmin2: f64 = p2Tmin_reg2(p);
    let mut T1: f64 = Tmin2;
    let v1: f64 = pT2v_reg2(p, T1);
    //let mut T2:f64 = T_MAX2;
    let T2: f64 = 1.1 * T1; //fast ,because the value of volume in region2 ig larger
    let v2 = pT2v_reg2(p, T2);
    T1 = T1 + (T2 - T1) * (v - v1) / (v2 - v1).abs();
    let f1:f64 = v - pT2v_reg2(p, T1);
    if f1.abs() < ESP {
        return T1;
    }
    let f = v - v2;
    let mut T: f64 = rtsec2(pT2v_reg2, p, v, T1, T2, f1, f, ESP, I_MAX);
    let mut v0: f64 = pT2v_reg2(p, T);
    if (v0 - v).abs()< ESP {
        return T;
    }
    if T < Tmin2 {
        T = Tmin2;
    }
    if T > T_MAX2 {
        T = T_MAX2;
    }
    // zoom the solution
    let mut steps: i32 = 0;
    let MAX_STEPS: i32 = 1000;
    while steps < MAX_STEPS {
        T += if v0 > v { -0.001 } else { 0.001  };
        T = T.clamp(Tmin2, T_MAX2);
        v0 = pT2v_reg2(p, T);
        steps += 1;
        if (v0 - v).abs() < ESP {
            break;
        }
    }
    T
}

/// Region 2  (T,v)->p using the secant method
///      T: temperature  K
///      p: pressure  MPa
//       v: specific volume m^3/kg
pub fn Tv2p_reg2(T: f64, v: f64) -> f64 {
    let mut p1: f64 = P_MIN2;

    let mut v1: f64 = pT2v_reg2(p1, T);
    if (v - v1).abs() < ESP {
        return p1;
    }

    let stepa: f64 = 1.0;
    let stepm: f64 = 5.0;
    let mut p: f64 = 0.0;
    let mut p2: f64 = p1 * stepm;
    let mut v2: f64 = pT2v_reg2(p2, T);
    if (v - v2).abs() < ESP {
        return p2;
    }

    let mut bounded: bool = false;
    if (v > v2) && (v < v1) {
        bounded = true;
    }

    let pmax2: f64 = T2pmax_reg2(T);
    while !bounded {
        p1 = p2;
        v1 = v2;
        if p1 > 1.0 {
            p2 = p1 + stepa;
        }
        if p1 < 1.0 {
            p2 = p1 * stepm;
        }
        if p2 >= pmax2 {
            p2 = pmax2;
            bounded = true;
        }
        v2 = pT2v_reg2(p2, T);
        if (v - v2).abs() < ESP {
            return p2;
        }
        if (v > v2) && (v < v1) {
            bounded = true;
        }
    }

    let f2: f64 = v - v2;
    let pmid: f64 = p2 - (p2 - p1) * (v - v2) / (v1 - v2);
    if v < pT2v_reg2(pmid, T) {
        p1 = pmid;
    }
    if p1 < P_MIN2 {
        p1 = P_MIN2;
    }
    let f1: f64 = v - pT2v_reg2(p1, T);
    if f1.abs() < ESP {
        return p1;
    }

    let mut p: f64 = rtsec1(pT2v_reg2, T, v, p1, p2, f1, f2, ESP, I_MAX);
    if p < P_MIN2 {
        p = P_MIN2;
    } else if p > pmax2 {
        p = pmax2;
    }
    p
}

/// Region 2  (T,h)->p using the secant method
///      T: temperature  K
///      h: specific enthalpy kJ/kg
///      p: pressure  MPa
pub fn Th2p_reg2(T: f64, h: f64) -> f64 {
    const stepa: f64 = 1.0;
    const stepm: f64 = 5.0;
    let pmax2: f64 = T2pmax_reg2(T);
    let mut p1: f64 = P_MIN2;
    let mut h1: f64 = pT2h_reg2(p1, T);
    if (h - h1).abs() < ESP {
        return p1;
    }
    let mut p2: f64 = p1 * stepm;
    let mut h2: f64 = pT2h_reg2(p2, T);
    if (h - h2).abs() < ESP {
        return p2;
    }

    let mut bounded: bool = false;
    if (h > h2) && (h < h1) {
        bounded = true;
    }

    while !bounded {
        p1 = p2;
        h1 = h2;
        if p1 > 1.0 {
            p2 = p1 + stepa;
        }
        if p1 < 1.0 {
            p2 = p1 * stepm;
        }
        if p2 >= pmax2 {
            p2 = pmax2;
            bounded = true;
        }
        h2 = pT2h_reg2(p2, T);
        if (h - h2).abs() < ESP {
            return p2;
        }
        if (h > h2) && (h < h1) {
            bounded = true;
        }
    }
    let mut pmid: f64 = p2 - (p2 - p1) * (h - h2) / (h1 - h2);
    if pmid < P_MIN2 {
        pmid = P_MIN2;
    };
    let h_mid: f64 = pT2h_reg2(pmid, T);
    if (h_mid - h).abs() < ESP {
        return pmid;
    }
    if h <= h_mid {
        p1 = pmid;
    }
    h1 = pT2h_reg2(p1, T);
    let f1: f64 = h - h1;
    let f2: f64 = h - h2;
    let mut p: f64 = rtsec1(pT2h_reg2, T, h, p1, p2, f1, f2, ESP, I_MAX);
    if p < P_MIN2 {
        p = P_MIN2;
    }
    if p > pmax2 {
        p = pmax2;
    }
    p
}

/// Region 2  (T,s)->p using the secant method
///  * T: temperature  K
///  * s: specific entropy  kJ/(kg K)
///  * p: pressure  MPa
pub fn Ts2p_reg2(T: f64, s: f64) -> f64 {
    let pmax2: f64 = T2pmax_reg2(T);

    let mut p1: f64 = P_MIN2;
    let mut s1: f64 = pT2s_reg2(p1, T);
    let mut f1: f64 = s - s1;
    let p2: f64 = pmax2;
    let s2: f64 = pT2s_reg2(p2, T);
    let f2: f64 = s - s2;
    p1 = p2 - (p2 - p1) * (s - s2) / (s1 - s2);
    if p1 < P_MIN2 {
        p1 = P_MIN2;
    }
    s1 = pT2s_reg2(p1, T);
    f1 = s - s1;
    if f1.abs() < ESP {
        return p1;
    }
    let mut p: f64 = rtsec1(pT2s_reg2, T, s, p1, p2, f1, f2, ESP, I_MAX);
    if p < P_MIN2 {
        p = P_MIN2;
    } else if p > pmax2 {
        p = pmax2;
    }
    p
}
