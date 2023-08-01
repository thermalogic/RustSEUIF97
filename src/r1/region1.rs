//! Region 1 - (p,T),(p,h), (p,s),(h,s); (p,v),(t,v),(t,h),(t,s)

use crate::common::constant::*;
use crate::common::propertry_id::*;

use crate::r1::*;

/// Region 1 : pT_reg1 - basic and extended properties
pub fn pT_reg1(p: f64, T: f64, o_id: i32) -> f64 {
    match o_id {
        OV => pT2v_reg1(p, T),
        OD => 1.0 / pT2v_reg1(p, T),
        OH => pT2h_reg1(p, T),
        OS => pT2s_reg1(p, T),
        OU => pT2u_reg1(p, T),
        OCV => pT2cv_reg1(p, T),
        OCP => pT2cp_reg1(p, T),
        OW => pT2w_reg1(p, T),
        _ => pT_ext_reg1(p, T, o_id), // the extended properties
    }
}

pub fn pt_reg1(p: f64, t: f64, o_id: i32) -> f64 {
    let T: f64 = t + K;
    pT_reg1(p, T, o_id)
}

pub fn ph_reg1(p: f64, h: f64, o_id: i32) -> f64 {
    let T: f64 = ph2T_reg1(p, h);
    if o_id == OT {
        return T - 273.15;
    };
    pT_reg1(p, T, o_id)
}

pub fn ps_reg1(p: f64, s: f64, o_id: i32) -> f64 {
    let T: f64 = ps2T_reg1(p, s);
    if o_id == OT {
        return T - 273.15;
    };
    pT_reg1(p, T, o_id)
}

pub fn hs_reg1(h: f64, s: f64, o_id: i32) -> f64 {
    let p: f64 = hs2p_reg1(h, s);
    if o_id == OP {
        return p;
    }
    ph_reg1(p, h, o_id)
}

// IAPWS-IF97 Region1: The extended input pair
//  (p,v)  (t,h),(t,s),(t,v)

/// Region1: (p,v)
pub fn pv_reg1(p: f64, v: f64, o_id: i32) -> f64 {
    let T: f64 = pv2T_reg1(p, v);
    if o_id == OT {
        return T - 273.15;
    };
    pT_reg1(p, T, o_id)
}

/// Region1: (T,v)
pub fn tv_reg1(t: f64, v: f64, o_id: i32) -> f64 {
    if o_id == OD {
        return 1.0 / v;
    };
    let p: f64 = Tv2p_reg1(t + 273.15, v);
    if o_id == OP {
        return p;
    };
    pT_reg1(p, t + 273.15, o_id)
}

/// Region1: (t,h)
pub fn th_reg1(t: f64, h: f64, o_id: i32) -> f64 {
    let p: f64 = Th2p_reg1(t + 273.15, h);
    if o_id == OP {
        return p;
    };
    pT_reg1(p, t + 273.15, o_id)
}

/// Region1: (t, s)
pub fn ts_reg1(t: f64, s: f64, o_id: i32) -> f64 {
    let p: f64 = Ts2p_reg1(t + 273.15, s);
    if o_id == OP {
        return p;
    };
    pT_reg1(p, t + 273.15, o_id)
}
