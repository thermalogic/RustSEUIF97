//! Region 3 - (p,T),(p,h), (p,s),(h,s); (p,v),(t,v),(t,h),(t,s)
//！* Basic Equation of  IAPWS -IF 97 Region3 IAPWS-IF97, R7-97(2012)
//！    * (T,d)-> p,h,u,s,cp,cv,w
//！* Backward Equation for Region 3.
//!    * (p,h)-> T,v   (p,s)->T,v in IAPWS-IF97-S03rev
//!    * (h,s)-> p  in IAPWS-IF97-S04rev, the region methods of Supp-phs3-2014.pdf in the `boundarie` module
//!    * (p,T)-> d in IAPWS-IF97-S05rev

use crate::common::propertry_id::*;
use crate::common::*;
use crate::r3::*;

/// Region 3 : Td_regs - all  properties
pub fn Td_reg3(T: f64, d: f64, o_id: i32) -> f64 {
    match o_id {
        OST | ODV | OKV | OTC | OSDC | OPR | OTD => pair_transport_reg(T, d, o_id, PAIRS::Td, Td_thermal_reg3),
        _ => Td_thermal_reg3(T, d, o_id), // the  basic and extended properties
    }
}

/// Region 3 : Td_thermal_regs - the basic and extended properties
pub fn Td_thermal_reg3(T: f64, d: f64, o_id: i32) -> f64 {
    match o_id {
        OV => 1.0 / d,
        OP => Td2p_reg3(T, d),
        OH => Td2h_reg3(T, d),
        OS => Td2s_reg3(T, d),
        OU => Td2u_reg3(T, d),
        OCV => Td2cv_reg3(T, d),
        OCP => Td2cp_reg3(T, d),
        OW => Td2w_reg3(T, d),
        _ => Td_ext_reg3(T, d, o_id), // the extended properties
    }
}

pub fn td_reg3(t: f64, d: f64, o_id: i32) -> f64 {
    Td_reg3(t + K, d, o_id)
}

pub fn pT_reg3(p: f64, T: f64, o_id: i32) -> f64 {
    let v: f64 = pT2v_reg3(p, T);
    if o_id == OV {
        return v;
    }
    let d: f64 = 1.0 / v;
    if o_id == OD {
        return d;
    }
    Td_reg3(T, d, o_id)
}

pub fn pt_reg3(p: f64, t: f64, o_id: i32) -> f64 {
    let T: f64 = t + 273.15;
    pT_reg3(p, T, o_id)
}

pub fn ph_reg3(p: f64, h: f64, o_id: i32) -> f64 {
    let v: f64 = ph2v_reg3(p, h);
    if o_id == OV {
        return v;
    };
    let d: f64 = 1.0 / v;
    if o_id == OD {
        return d;
    };
    let T: f64 = ph2T_reg3(p, h);
    if o_id == OT {
        return T - K;
    };
    Td_reg3(T, d, o_id)
}

pub fn ps_reg3(p: f64, s: f64, o_id: i32) -> f64 {
    let v: f64 = ps2v_reg3(p, s);
    if o_id == OV {
        return v;
    };
    let d: f64 = 1.0 / v;
    if o_id == OD {
        return d;
    };
    let T: f64 = ps2T_reg3(p, s);
    if o_id == OT {
        return T - K;
    };
    Td_reg3(T, d, o_id)
}

pub fn hs_reg3(h: f64, s: f64, o_id: i32) -> f64 {
    let p: f64 = hs2p_reg3(h, s);
    if o_id == OP {
        return p;
    }
    let T: f64 = ph2T_reg3(p, h);
    if o_id == OT {
        return T - 273.15;
    }
    let v: f64 = ph2v_reg3(p, h);
    if o_id == OV {
        return v;
    }
    let d = 1.0 / v;
    if o_id == OD {
        return d;
    }
    Td_reg3(T, d, o_id)
}

/// IAPWS-IF97 Region3: The extended input pair
/// * (p,v) (t,v), (t,h),(t,s)

///  Region3:  (p,v)
pub fn pv_reg3(p: f64, v: f64, o_id: i32) -> f64 {
    let d: f64 = 1.0 / v;
    if o_id == OD {
        return d;
    };
    let T: f64 = pv2T_reg3(p, v);
    if o_id == OT {
        return T - K;
    };
    Td_reg3(T, d, o_id)
}

///  Region3:  (t,v)
pub fn tv_reg3(t: f64, v: f64, o_id: i32) -> f64 {
    let d: f64 = 1.0 / v;
    if o_id == OD {
        return d;
    };
    Td_reg3(t + K, d, o_id)
}

///  Region3:  (t,h)
pub fn th_reg3(t: f64, h: f64, o_id: i32) -> f64 {
    let d = Th2d_reg3(t + 273.15, h);
    if o_id == OD {
        return d;
    };
    if o_id == OV {
        return 1.0 / d;
    };
    Td_reg3(t + K, d, o_id)
}

///  Region3:  (t,s)
pub fn ts_reg3(t: f64, s: f64, o_id: i32) -> f64 {
    let d = Ts2d_reg3(t + 273.15, s);
    if o_id == OD {
        return d;
    };
    if o_id == OV {
        return 1.0 / d;
    };
    Td_reg3(t + K, d, o_id)
}
