#![allow(warnings)]
//use libc;
pub mod algo;
pub mod common;
pub mod r1;
pub mod r2;
pub mod r3;
pub mod r4;
pub mod r5;

pub use common::constant::*;
pub use common::region::*;
pub use common::transport_further::*;
pub use r1::*;
pub use r2::*;
pub use r3::*;
pub use r4::*;
pub use r5::*;

//pub unsafe extern "C" fn pt(p:f64,t:f64) -> f64
pub fn pt_thermal(p: f64, t: f64, o_id: i32) -> f64 {
    let sub_region: i32 = pT_sub_region(p, t + K);
    if o_id == OR {
        return sub_region as f64;
    };
    if  o_id == OX 
    {
       match sub_region{
          1=> {return 0.0},
          2|3|5=> {return 1.0},
          _=> { return INVALID_VALUE as f64}
       }
    }
    match sub_region {
        1 => pt_reg1(p, t, o_id),
        2 => pt_reg2(p, t, o_id),
        3 => pt_reg3(p, t, o_id),
        4 => {
            if ((t + 273.15) == TC_WATER && p == PC_WATER) {
                return td_reg3(t, DC_WATER, o_id);
            } else {
                return sub_region as f64;
            }
        }
        5 => pt_reg5(p, t, o_id),
        _ => sub_region as f64,
    }
}

pub fn pt(p: f64, t: f64, o_id: i32) -> f64 {
    if o_id == OST {
        return surface_tension(t + K);
    }
    match o_id {
        OP | OT | OV | OD | OH | OS | OU | OCP | OCV | OW | OR  | OX => pt_thermal(p, t, o_id),
        ODV | OKV |OTC | OSDC => {
            let d: f64 = pt_thermal(p, t, OD);
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV 
                {
                    value /=d;  // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                }
            };
            value
        }
        OPR | OTD=> {
            let d: f64 = pt_thermal(p, t, OD);
            let cp: f64 = pt_thermal(p, t, OCP);
            let tc: f64 = thcond(d, t + 273.15);
            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }
        _ => INVALID_VALUE as f64,
    }
}

pub fn ph_thermal(p: f64, h: f64, o_id: i32) -> f64 {
    let sub_region: i32 = ph_sub_region(p, h);
    if o_id == OR {
        return sub_region as f64;
    };
    if  o_id == OX 
    {
       match sub_region{
          1=> {return 0.0},
          2|3|5=> {return 1.0},
          4=> { return ph_reg4(p, h, OX)},
          _=> { return INVALID_VALUE as f64}
       }
    }
    match sub_region {
        1 => ph_reg1(p, h, o_id),
        2 => ph_reg2(p, h, o_id),
        3 => ph_reg3(p, h, o_id),
        4 => ph_reg4(p, h, o_id),
        5 => ph_reg5(p, h, o_id),
        _ => sub_region as f64,
    }
}

pub fn ph(p: f64, h: f64, o_id: i32) -> f64 {
    match o_id {
        OP | OT | OV | OD | OH | OS | OU | OCP | OCV | OW | OR | OX => ph_thermal(p, h, o_id),
        OST => {
            let t: f64 = ph_thermal(p, h, OT);
            return surface_tension(t + K);
        }
        ODV | OKV |OTC | OSDC => {
            let d: f64 = ph_thermal(p, h, OD);
            let t: f64 = ph_thermal(p, h, OT);
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV{
                value = viscosity(d, t + 273.15);
                if o_id == OKV 
                {
                    value /=d;  // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                }
            };
            value
        }
        OPR | OTD => {
            let d: f64 = ph_thermal(p, h, OD);
            let t: f64 = ph_thermal(p, h, OT);

            let cp: f64 = ph_thermal(p, h, OCP);
            let tc: f64 = thcond(d, t + 273.15);

            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }
        _ => INVALID_VALUE as f64,
    }
}

pub fn ps_thermal(p: f64, s: f64, o_id: i32) -> f64 {
    let sub_region: i32 = ps_sub_region(p, s);
    if o_id == OR {
        return sub_region as f64;
    };
    if  o_id == OX 
    {
       match sub_region{
          1=> {return 0.0},
          2|3|5=> {return 1.0},
          4=> { return ps_reg4(p, s, OX)},
          _=> { return INVALID_VALUE as f64}
       }
    }
    match sub_region {
        1 => ps_reg1(p, s, o_id),
        2 => ps_reg2(p, s, o_id),
        3 => ps_reg3(p, s, o_id),
        4 => ps_reg4(p, s, o_id),
        5 => ps_reg5(p, s, o_id),
        _ => sub_region as f64,
    }
}

pub fn ps(p: f64, s: f64, o_id: i32) -> f64 {
    match o_id {
        OP | OT | OV | OD | OH | OS | OU | OCP | OCV | OW | OR | OX => ps_thermal(p, s, o_id),
        OST => {
            let t: f64 = ps_thermal(p, s, OT);
            return surface_tension(t + K);
        }
        ODV | OKV | OTC | OSDC => {
            let d: f64 = ps_thermal(p, s, OD);
            let t: f64 = ps_thermal(p, s, OT);
            let mut value: f64 = 0.0;
            if o_id == ODV ||  o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV 
                {
                    value /=d;  // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                };
            };
            value
        }
        OPR | OTD => {
            let d: f64 = ps_thermal(p, s, OD);
            let t: f64 = ps_thermal(p, s, OT);

            let cp: f64 = ps_thermal(p, s, OCP);
            let tc: f64 = thcond(d, t + 273.15);
            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }
        _ => INVALID_VALUE as f64,
    }
}

pub fn hs_thermal(h: f64, s: f64, o_id: i32) -> f64 {
    let sub_region: i32 = hs_sub_region(h, s);
    if o_id == OR {
        return sub_region as f64;
    };
    if  o_id == OX 
    {
       match sub_region{
          1=> {return 0.0},
          2|3|5=> {return 1.0},
          4=> { return hs_reg4(h,s,OX)},
          _=> { return INVALID_VALUE as f64}
       }
    }
    match sub_region {
        1 => hs_reg1(h, s, o_id),
        2 => hs_reg2(h, s, o_id),
        3 => hs_reg3(h, s, o_id),
        4 => hs_reg4(h, s, o_id),
        5 => hs_reg5(h, s, o_id),
        _ => sub_region as f64,
    }
}

pub fn hs(h: f64, s: f64, o_id: i32) -> f64 {
    match o_id {
        OP | OT | OV | OD | OH | OS | OU | OCP | OCV | OW | OR | OX => hs_thermal(h, s, o_id),
        OST => {
            let t: f64 = hs_thermal(h, s, OT);
            return surface_tension(t + K);
        }
        ODV | OKV  | OTC | OSDC => {
            let d: f64 = hs_thermal(h, s, OD);
            let t: f64 = hs_thermal(h, s, OT);
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV 
                {
                    value /=d;  // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                };
            };
            value
        }
        OPR | OTD => {
            let d: f64 = hs_thermal(h, s, OD);
            let t: f64 = hs_thermal(h, s, OT);

            let cp: f64 = hs_thermal(h, s, OCP);
            let tc: f64 = thcond(d, t + 273.15);
            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }

        _ => INVALID_VALUE as f64,
    }
}
