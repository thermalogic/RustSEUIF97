//! Region 2 - (p,T),(p,h), (p,s),(h,s); (p,v),(t,v),(t,h),(t,s)
//!
//！IAPWS-IF97 and supp release
//！  1. IF97 IAPWS, R7-97(2012) IF97-Rev.pdf: P12-32
//！      * basic: (p,t)->v,u,s,h,cp,cv,w
//！      * backward: (p,h)->T, (p,s)->T
//！  2. IAPWS, SR2-01(2014)
//！       * Supp-PHS12-2014.pdf  (h,s)->p

use crate::common::constant::*;
use crate::common::propertry_id::*;
use crate::r2::*;

/// Region 2 : pT_reg2 - the basic and extended properties
pub fn pT_reg2(p: f64, T: f64, o_id: i32) -> f64 {
    match o_id {
        OV => pT2v_reg2(p, T),
        OD => 1.0 / pT2v_reg2(p, T),
        OH => pT2h_reg2(p, T),
        OS => pT2s_reg2(p, T),
        OU => pT2u_reg2(p, T),
        OCV => pT2cv_reg2(p, T),
        OCP => pT2cp_reg2(p, T),
        OW => pT2w_reg2(p, T),
        _ => pT_ext_reg2(p, T, o_id), // the extended properties
    }
}

pub fn pt_reg2(p: f64, t: f64, o_id: i32) -> f64 {
    let T = t + 273.15;
    pT_reg2(p, T, o_id)
}

pub fn ph_reg2(p: f64, h: f64, o_id: i32) -> f64 {
    let T: f64 = ph2T_reg2(p, h);
    if o_id == OT {
        return T - 273.15;
    };
    pT_reg2(p, T, o_id)
}

pub fn ps_reg2(p: f64, s: f64, o_id: i32) -> f64 {
    let T: f64 = ps2T_reg2(p, s);
    if o_id == OT {
        return T - 273.15;
    };
    pT_reg2(p, T, o_id)
}

pub fn hs_reg2(h: f64, s: f64, o_id: i32) -> f64 {
    let p: f64 = hs2p_reg2(h, s);
    if o_id == OP {
        return p;
    };
    ph_reg2(p, h, o_id)
}

/// IAPWS-IF97 Region 2 : The extended input pair
/// *  (p,v)  (t,v),(t,h),(t,s)

/// Region 2 : (p,v)
pub fn pv_reg2(p: f64, v: f64, o_id: i32) -> f64 {
    let T: f64 = pv2T_reg2(p, v);
    if o_id == OT {
        return T - 273.15;
    };
    pT_reg2(p, T, o_id)
}

///  Region 2 :(t,v)
pub fn tv_reg2(t: f64, v: f64, o_id: i32) -> f64 {
    let T: f64 = t + K;
    let p: f64 = Tv2p_reg2(T, v);
    if o_id == OP {
        return p;
    };
    pT_reg2(p, T, o_id)
}

///  Region 2 :(t,h)
pub fn th_reg2(t: f64, h: f64, o_id: i32) -> f64 {
    let p: f64 = Th2p_reg2(t + 273.15, h);
    if o_id == OP {
        return p;
    };
    pT_reg2(p, t + 273.15, o_id)
}

/// Region 2 : (t, s)
pub fn ts_reg2(t: f64, s: f64, o_id: i32) -> f64 {
    let p: f64 = Ts2p_reg2(t + 273.15, s);
    if o_id == OP {
        return p;
    };
    pT_reg2(p, t + 273.15, o_id)
}
