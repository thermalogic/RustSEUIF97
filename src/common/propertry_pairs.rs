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
type THERMAL_EQ = fn(f64, f64, i32, i32) -> f64;

pub enum PAIRS {
    pT,
    ph,
    ps,
    pv,
    hs,
    tv,
    th,
    ts,
}

/// The common method of input pairs
pub fn pair_transport(
    v1: f64,
    v2: f64,
    o_id: i32,
    pair: PAIRS,
    fn_thermal: THERMAL_EQ,
    reg: i32,
) -> f64 {
    match o_id {
        OST => match pair {
            PAIRS::tv | PAIRS::th | PAIRS::ts => surface_tension(v1 + 273.15),
            PAIRS::pT => surface_tension(v2),
            _ => surface_tension(fn_thermal(v1, v2, OT, reg) + 273.15),
        },

        ODV => match pair {
            // d,T
            PAIRS::tv => viscosity(1.0 / v2, v1 + 273.15),
            PAIRS::th | PAIRS::ts => viscosity(fn_thermal(v1, v2, OD, reg), v1 + 273.15),
            PAIRS::pT => viscosity(fn_thermal(v1, v2, OD, reg), v2),
            PAIRS::pv => viscosity(1.0 / v2, fn_thermal(v1, v2, OT, reg) + 273.15),
            _ => viscosity(
                fn_thermal(v1, v2, OD, reg),
                fn_thermal(v1, v2, OT, reg) + 273.15,
            ),
        },
        OKV => match pair {// d,T
            PAIRS::tv => v2 * viscosity(1.0 / v2, v1 + 273.15),
            PAIRS::th | PAIRS::ts => {
                let d = fn_thermal(v1, v2, OD, reg);
                viscosity(d, v1 + 273.15) / d
            }
            PAIRS::pT => {
                let d = fn_thermal(v1, v2, OD, reg);
                viscosity(d, v2) / d
            }
            PAIRS::pv => {
                let T = fn_thermal(v1, v2, OT, reg) + 273.15;
                v2 * viscosity(1.0 / v2, T)
            }
            _ => {
                let d = fn_thermal(v1, v2, OD, reg);
                let T = fn_thermal(v1, v2, OT, reg) + 273.15;
                viscosity(d, T) / d
            }
        },

        OTC => match pair { //d,T
            PAIRS::tv => thcond(1.0 / v2, v1 + 273.15),
            PAIRS::th | PAIRS::ts => {
                let d = fn_thermal(v1, v2, OD, reg);
                thcond(d, v1 + 273.15)
            }
            PAIRS::pT => {
                let d = fn_thermal(v1, v2, OD, reg);
                thcond(d, v2)
            }

            PAIRS::pv => {
                let T = fn_thermal(v1, v2, OT, reg) + 273.15;
                thcond(1.0 / v2, T)
            }
            _ => {
                let d = fn_thermal(v1, v2, OD, reg);
                let T = fn_thermal(v1, v2, OT, reg) + 273.15;
                thcond(d, T)
            }
        },
        OSDC => match pair {// d,T
            PAIRS::tv => static_dielectric(1.0 / v2, v1 + 273.15),
            PAIRS::th | PAIRS::ts => {
                let d = fn_thermal(v1, v2, OD, reg);
                static_dielectric(d, v1 + 273.15)
            }
            PAIRS::pT => {
                let d = fn_thermal(v1, v2, OD, reg);
                static_dielectric(d, v2)
            }
            PAIRS::pv => {
                let T = fn_thermal(v1, v2, OT, reg) + 273.15;
                static_dielectric(1.0 / v2, T)
            }
            _ => {
                let d = fn_thermal(v1, v2, OD, reg);
                let T = fn_thermal(v1, v2, OT, reg) + 273.15;
                static_dielectric(d, T)
            }
        },
        OTD => match pair {
            PAIRS::tv => {
                let d = 1.0 / v2;
                let tc: f64 = thcond(d, v1 + 273.15);
                let cp: f64 = fn_thermal(v1, v2, OCP, reg);
                thermal_diffusivity(tc, cp, d)
            }
            PAIRS::th | PAIRS::ts => {
                let d = fn_thermal(v1, v2, OD, reg);
                let tc: f64 = thcond(d, v1 + 273.15);
                let cp: f64 = fn_thermal(v1, v2, OCP, reg);
                thermal_diffusivity(tc, cp, d)
            }
            PAIRS::pT => {
                let d = fn_thermal(v1, v2, OD, reg);
                let tc: f64 = thcond(d, v2);
                let cp: f64 = fn_thermal(v1, v2, OCP, reg);
                thermal_diffusivity(tc, cp, d)
            }
            PAIRS::pv => {
                let d = 1.0 / v2;
                let T = fn_thermal(v1, v2, OT, reg) + 273.15;
                let tc: f64 = thcond(d, T);
                let cp: f64 = fn_thermal(v1, v2, OCP, reg);
                thermal_diffusivity(tc, cp, d)
            }
            _ => {
                let d = fn_thermal(v1, v2, OD, reg);
                let T = fn_thermal(v1, v2, OT, reg) + 273.15;
                let tc: f64 = thcond(d, T);
                let cp: f64 = fn_thermal(v1, v2, OCP, reg);
                thermal_diffusivity(tc, cp, d)
            }
        },
        OPR => match pair {
            PAIRS::tv => {
                let d = 1.0 / v2;
                let T = v1 + 273.15;
                let tc: f64 = thcond(d, T);
                let cp: f64 = fn_thermal(v1, v2, OCP, reg);
                let dv: f64 = viscosity(d, T);
                prandtl_number(dv, cp, tc)
            }
            PAIRS::th | PAIRS::ts => {
                let T = v1 + 273.15;
                let d = fn_thermal(v1, v2, OD, reg);
                let dv: f64 = viscosity(d, T);
                let cp: f64 = fn_thermal(v1, v2, OCP, reg);
                let tc: f64 = thcond(d, T);
                prandtl_number(dv, cp, tc)
            }
            PAIRS::pv => {
                let d = 1.0 / v2;
                let T = fn_thermal(v1, v2, OT, reg) + 273.15;
                let dv: f64 = viscosity(d, T);
                let cp: f64 = fn_thermal(v1, v2, OCP, reg);
                let tc: f64 = thcond(d, T);
                prandtl_number(dv, cp, tc)
            }
            _ => {
                let d = fn_thermal(v1, v2, OD, reg);
                let T = fn_thermal(v1, v2, OT, reg) + 273.15;
                let dv: f64 = viscosity(d, T);
                let cp: f64 = fn_thermal(v1, v2, OCP, reg);
                let tc: f64 = thcond(d, T);
                prandtl_number(dv, cp, tc)
            }
        },
        _ => INVALID_OUTID as f64,
    }
}

/// The common method of input pairs
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

/// pT_thermal(p,T,o_id)：the propertry of o_id (basic and extended thermodynamic)
pub fn pT_thermal(p: f64, T: f64, o_id: i32, reg: i32) -> f64 {
    pair_thermal(
        p,
        T,
        o_id,
        pT_sub_region,
        pT_reg1,
        pT_reg2,
        pT_reg3,
        pT_reg4,
        pT_reg5,
        reg,
    )
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
