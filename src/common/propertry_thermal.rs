//! The computing function of basic and extended Thermodynamic  Properties  
//! # Input pairs:
//! *  (p,t) (p,h) (p,s) (h,s)
//! *  (p,v),(t,v),(t,h),(t,s)
//!   
//! # The Basic Thermodynamic Properties
//! *  P: Pressure  MPa
//! *  T: temperature K
//! *  d: Density  kg/m^3
//! *  v: specific volume m^3/kg
//! *  h: specific enthalpy kJ/kg
//! *  u: specific internal energy kJ/kg
//! *  s: specific entropy  kJ/(kg K)
//! *  cp: specific isobaric heat capacity  kJ/(kg K)
//! *  cv: specific isochoric heat capacity kJ/(kg K)
//！*  w:  speed of sound  m/s
//! *  x:  Steam quality
//! *  r:  Region
//! # The extended hermodynamic Properties
//!

use crate::common::*;
use crate::r1::region1::*;
use crate::r2::region2::*;
use crate::r3::region3::*;
use crate::r4::region4::*;
use crate::r5::region5::*;

type REGION_EQ = fn(f64, f64) -> i32;
type PROP_EQ = fn(f64, f64, i32) -> f64;

/// The common method ofcd input pairs except (p,t)
pub fn pair_thermal(
    v1: f64,
    v2: f64,
    o_id: i32,
    fr: REGION_EQ,
    f1: PROP_EQ,
    f2: PROP_EQ,
    f3: PROP_EQ,
    f4: PROP_EQ,
    f5: PROP_EQ,
    reg: i32,
) -> f64 {
    let mut sub_region: i32 = reg;
    if reg == 6 {
        sub_region = fr(v1, v2)
    }
    if o_id == OR {
        return sub_region as f64;
    };
    if o_id == OX {
        match sub_region {
            1 => return 0.0,
            2 | 3 | 5 => return 1.0,
            4 => return f4(v1, v2, OX),
            _ => return INVALID_VALUE as f64,
        }
    }
    match sub_region {
        1 => f1(v1, v2, o_id),
        2 => f2(v1, v2, o_id),
        3 => f3(v1, v2, o_id),
        4 => f4(v1, v2, o_id),
        5 => f5(v1, v2, o_id),
        _ => sub_region as f64,
    }
}

/// pt_thermal(p,t,o_id)：the propertry of o_id (basic and extended thermodynamic)
pub fn pt_thermal(p: f64, t: f64, o_id: i32, reg: i32) -> f64 {
    let mut sub_region: i32 = reg;
    if reg == 6 {
        sub_region =  pT_sub_region(p,t+273.15);
    }
    if o_id == OR {
        return sub_region as f64;
    };
    if o_id == OX {
        match sub_region {
            1 => return 0.0,
            2 | 3 | 5 => return 1.0,
            _ => return INVALID_VALUE as f64,
        }
    }
    match sub_region {
        1 => pt_reg1(p, t, o_id),
        2 => pt_reg2(p, t, o_id),
        3 => pt_reg3(p, t, o_id),
        4 => {
            if (t + 273.15) == TC_WATER && p == PC_WATER {
                return td_reg3(t, DC_WATER, o_id);
            } else {
                return sub_region as f64;
            }
        }
        5 => pt_reg5(p, t, o_id),
        _ => sub_region as f64,
    }
}

/// ph_thermal(p,h,o_id)：the propertry of o_id (basic and extended thermodynamic)
pub fn ph_thermal(p: f64, h: f64, o_id: i32, reg: i32) -> f64 {
    pair_thermal(
        p,
        h,
        o_id,
        ph_sub_region,
        ph_reg1,
        ph_reg2,
        ph_reg3,
        ph_reg4,
        ph_reg5,
        reg,
    )
}

/// ps_thermal(h,s,o_id)：the propertry of o_id (basic and extended thermodynamic)
pub fn ps_thermal(p: f64, s: f64, o_id: i32, reg: i32) -> f64 {
    pair_thermal(
        p,
        s,
        o_id,
        ps_sub_region,
        ps_reg1,
        ps_reg2,
        ps_reg3,
        ps_reg4,
        ps_reg5,
        reg,
    )
}

/// hs_thermal(h,s,o_id)：the propertry of o_id (basic and extended thermodynamic)
pub fn hs_thermal(h: f64, s: f64, o_id: i32, reg: i32) -> f64 {
    pair_thermal(
        h,
        s,
        o_id,
        hs_sub_region,
        hs_reg1,
        hs_reg2,
        hs_reg3,
        hs_reg4,
        hs_reg5,
        reg,
    )
}

//------------------------------------------------------------------
//  Functions of the  extended input pairs (p,v),(t,v),(t,s),(t,h)
//------------------------------------------------------------------

/// The  extended input pair pv_thermal(p,v,o_id)
pub fn pv_thermal(p: f64, v: f64, o_id: i32, reg: i32) -> f64 {
    pair_thermal(
        p,
        v,
        o_id,
        pv_sub_region,
        pv_reg1,
        pv_reg2,
        pv_reg3,
        pv_reg4,
        pv_reg5,
        reg,
    )
}

/// The extended input pair tv_thermal(t,v,o_id)
pub fn tv_thermal(t: f64, v: f64, o_id: i32, reg: i32) -> f64 {
    pair_thermal(
        t,
        v,
        o_id,
        tv_sub_region,
        tv_reg1,
        tv_reg2,
        tv_reg3,
        tv_reg4,
        tv_reg5,
        reg,
    )
}

/// The  extended input pair th_thermal(t,h,o_id)
pub fn th_thermal(t: f64, h: f64, o_id: i32, reg: i32) -> f64 {
    pair_thermal(
        t,
        h,
        o_id,
        th_sub_region,
        th_reg1,
        th_reg2,
        th_reg3,
        th_reg4,
        th_reg5,
        reg,
    )
}

/// The  extended input pair ts_thermal(t,s,o_id)
pub fn ts_thermal(t: f64, s: f64, o_id: i32, reg: i32) -> f64 {
    pair_thermal(
        t,
        s,
        o_id,
        ts_sub_region,
        ts_reg1,
        ts_reg2,
        ts_reg3,
        ts_reg4,
        ts_reg5,
        reg,
    )
}
