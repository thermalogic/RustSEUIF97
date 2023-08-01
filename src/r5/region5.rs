//! Region 5 - (p,T),(p,h), (p,s),(h,s); (p,v),(t,v),(t,h),(t,s)

use crate::common::constant::*;
use crate::common::propertry_id::*;
use crate::r5::*;

/// Region5 : pT_regs - the basic and extended properties
pub fn pT_reg5(p: f64, T: f64, o_id: i32) -> f64 {
    match o_id {
        OV => pT2v_reg5(p, T),
        OD => 1.0 / pT2v_reg5(p, T),
        OH => pT2h_reg5(p, T),
        OS => pT2s_reg5(p, T),
        OU => pT2u_reg5(p, T),
        OCV => pT2cv_reg5(p, T),
        OCP => pT2cp_reg5(p, T),
        OW => pT2w_reg5(p, T),
        _ => pT_ext_reg5(p,T,o_id), //the extended properties
    }
}

pub fn pt_reg5(p: f64, t: f64, o_id: i32) -> f64 {
    let T: f64 = t + K;
    pT_reg5(p, T, o_id)
}

pub fn ph_reg5(p: f64, h: f64, o_id: i32) -> f64 {
    let T: f64 = ph2T_reg5(p, h);
    if o_id == OT {
        return T - 273.15;
    };
    pT_reg5(p, T, o_id)
}

pub fn ps_reg5(p: f64, s: f64, o_id: i32) -> f64 {
    let T: f64 = ps2T_reg5(p, s);
    if o_id == OT {
        return T - 273.15;
    };
    pT_reg5(p, T, o_id)
}

pub fn hs_reg5(h: f64, s: f64, o_id: i32) -> f64 {
    let p: f64 = hs2p_reg5(h, s);
    if o_id == OP {
        return p;
    };
    ph_reg5(p, h, o_id)
}

/// Region 5 : The extended input pair: (p,v), (t,v),(t,h),(t,s)

///  Region 5 :(p,v)
pub fn pv_reg5(p: f64, v: f64, o_id: i32) -> f64 {
    let T: f64 = pv2T_reg5(p, v);
    if o_id == OT {
        return T - 273.15;
    };
    pT_reg5(p, T, o_id)
}

///  Region 5 : (t,v)
pub fn tv_reg5(t: f64, v: f64, o_id: i32) -> f64 {
    if o_id == OD {
        return 1.0 / v;
    };
    let p: f64 = Tv2p_reg5(t + 273.15, v);
    if o_id == OP {
        return p;
    };
    pT_reg5(p, t + 273.15, o_id)
}

///   Region 5 :(t,h)
pub fn th_reg5(t: f64, h: f64, o_id: i32) -> f64 {
    let p: f64 = Th2p_reg5(t + 273.15, h);
    if o_id == OP {
        return p;
    };
    pT_reg5(p, t + 273.15, o_id)
}

///  Region 5 : (t,s)
pub fn ts_reg5(t: f64, s: f64, o_id: i32) -> f64 {
    let p: f64 = Ts2p_reg5(t + 273.15, s);
    if o_id == OP {
        return p;
    };
    pT_reg5(p, t + 273.15, o_id)
}

