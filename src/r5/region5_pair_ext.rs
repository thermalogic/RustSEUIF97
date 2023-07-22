//！   IAPWS-IF97 Region5: The extended input pair
//!       (p,v)
//!       (t,v),(t,h),(t,s)
use crate::algo::root::*;
use crate::common::boundaries::*;
use crate::common::constant::*;
use crate::r5::region5_pT::*;

/// Region 5: (p,v)-> T
pub fn pv2T_reg5(p: f64, v: f64) -> f64 {
    let mut T1: f64 = T_MIN5;
    let mut v1: f64 = pT2v_reg5(p, T1);
    let mut T2: f64 = T_MAX5;
    let mut v2: f64 = pT2v_reg5(p, T2);
    let mut f: f64 = v - pT2v_reg5(p, T2);
    if (v2 - v1) != 0.0 {
        T1 = T1 + (T2 - T1) * (v - v1) / (v2 - v1);
    }
    let mut f1: f64 = v - pT2v_reg5(p, T1);
    let mut T: f64 = rtsec2(pT2v_reg5, p, v, T1, T2, f1, f, ESP, I_MAX);
    if T < T_MIN5 {
        T = T_MIN5;
    } else {
        if T > T_MAX5 {
            T = T_MAX5;
        };
    };
    T
}

/// Region 5: (T,v)-> p
pub fn Tv2p_reg5(T: f64, v: f64) -> f64 {
    let mut p1: f64 = P_MIN5;
    let mut v1: f64 = pT2v_reg5(p1, T);
    if v == v1 {
        return p1;
    }

    let mut p2: f64 = 10.0 * p1;
    if p2 > P_MAX5 {
        p2 = P_MAX5;
    };

    let mut v2: f64 = pT2v_reg5(p2, T);
    if v == v2 {
        return p2;
    }

    let mut vbounded: bool = false;
    if (v > v2) && (v < v1) {
        vbounded = true;
    }

    while !vbounded {
        p1 = p2;
        v1 = v2;
        p2 = 5.0 * p1;
        if p2 >= P_MAX5 {
            p2 = P_MAX5;
            vbounded = true;
        }
        v2 = pT2v_reg5(p2, T);
        if v == v2 {
            return p2;
        }
        if (v > v2) && (v < v1) {
            vbounded = true;
        }
    }

    let mut f2: f64 = v - v2;
    p1 = p2 - (p2 - p1) * (v - v2) / (v1 - v2);
    if p1 < P_MIN5 {
        p1 = P_MIN5;
    }
    let f1: f64 = v - pT2v_reg5(p1, T);
    if f1.abs() < ESP {
        return p1;
    };
    let mut p: f64 = rtsec1(pT2v_reg5, T, v, p1, p2, f1, f2, ESP, I_MAX);
    if p > P_MAX5 {
        p = P_MAX5;
    };
    if p < P_MIN5 {
        p = P_MIN5;
    };
    p
}

/// Region 5: (T,h)-> p
pub fn Th2p_reg5(T: f64, h: f64) -> f64 {
    let mut p1: f64 = P_MIN5;
    let mut h1: f64 = pT2h_reg5(p1, T);
    let mut p2: f64 = P_MAX5; //
    let mut h2: f64 = pT2h_reg5(p2, T);
    let f2: f64 = h - h2;
    p1 = p2 - (p2 - p1) * ((h - h2) / (h1 - h2)).abs();
    let f1: f64 = h - pT2h_reg5(p1, T);
    if f1.abs() < ESP {
        return p1;
    }
    let mut p: f64 = rtsec1(pT2h_reg5, T, h, p1, p2, f1, f2, ESP, I_MAX);
    if p > P_MAX5 {
        p = P_MAX5;
    }
    if p < P_MIN5 {
        p = P_MIN5;
    }
    p
}

/// Region 5: (T,s)-> s
pub fn Ts2p_reg5(T: f64, s: f64) -> f64 {
    let mut p1: f64 = P_MIN5;
    let mut s1: f64 = pT2s_reg5(p1, T);
    let mut f1: f64 = s - s1;
    let p2: f64 = P_MAX5;
    let s2: f64 = pT2s_reg5(p2, T);
    let f2: f64 = s - s2;
    p1 = p2 - (p2 - p1) * (s - s2) / (s1 - s2);
    if p1 < P_MIN5 {
        p1 = P_MIN5;
    }
    s1 = pT2s_reg5(p1, T);
    f1 = s - s1;
    if f1.abs() < ESP {
        return p1;
    }
    let mut p: f64 = rtsec1(pT2s_reg5, T, s, p1, p2, f1, f2, ESP, I_MAX);
    if p > P_MAX5 {
        p = P_MAX5;
    } else if p < P_MIN5 {
        p = P_MIN5;
    }
    p
}
