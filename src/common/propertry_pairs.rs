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
type THERMAL_REG_EQ = fn(f64, f64, i32) -> f64;

pub enum PAIRS {
    pT,
    Td,
}

/// The common method to get the transport properties of input pairs for each region
/// * (p,T) - 1,2,5; (T,d) -3
pub fn pair_transport_reg(v1: f64, v2: f64, o_id: i32, pair: PAIRS, fn_thermal: THERMAL_REG_EQ) -> f64 {
    match o_id {
        OST => match pair {
            PAIRS::Td => surface_tension(v1),
            PAIRS::pT => surface_tension(v2),
        },

        ODV => match pair {
            // d,T
            PAIRS::Td => viscosity(v2, v1),
            _ => viscosity(fn_thermal(v1, v2, OD), v2),
        },
        OKV => match pair {
            // d,T
            PAIRS::Td => v2 * viscosity(v2, v1),
            _ => {
                let d = fn_thermal(v1, v2, OD);
                viscosity(d, v2) / d
            }
        },

        OTC => match pair {
            //d,T
            PAIRS::Td => thcond(v2, v1),
            _ => {
                let d = fn_thermal(v1, v2, OD);
                thcond(d, v2)
            }
        },
        OSDC => match pair {
            // d,T
            PAIRS::Td => static_dielectric(v2, v1),
            _ => {
                let d = fn_thermal(v1, v2, OD);
                static_dielectric(d, v2)
            }
        },
        OTD => match pair {
            PAIRS::Td => {
                let tc: f64 = thcond(v2, v1);
                let cp: f64 = fn_thermal(v1, v2, OCP);
                thermal_diffusivity(tc, cp, v2)
            }
            _ => {
                let d = fn_thermal(v1, v2, OD);
                let tc: f64 = thcond(d, v2);
                let cp: f64 = fn_thermal(v1, v2, OCP);
                thermal_diffusivity(tc, cp, d)
            }
        },
        OPR => match pair {
            PAIRS::Td => {
                let tc: f64 = thcond(v2, v1);
                let cp: f64 = fn_thermal(v1, v2, OCP);
                let dv: f64 = viscosity(v2, v1);
                prandtl_number(dv, cp, tc)
            }
            _ => {
                let d = 1.0 / v2;
                let dv: f64 = viscosity(d, v2);
                let cp: f64 = fn_thermal(v1, v2, OCP);
                let tc: f64 = thcond(d, v2);
                prandtl_number(dv, cp, tc)
            }
        },
        _ => INVALID_OUTID as f64,
    }
}

/// The common method of input pairs to get all properties
pub fn pair_properties(
    v1: f64, v2: f64, o_id: i32, fr: REGION_EQ, f1: PROP_EQ, f2: PROP_EQ, f3: PROP_EQ, f4: PROP_EQ, f5: PROP_EQ,
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
