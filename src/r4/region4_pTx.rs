//! Region 4 saturation and wet
//!           p, MPa   T,  K
//!     p2sat_water,   p2sat_steam
//!     T2sat_water,  T2sat_steam

use crate::algo::fast_ipower::sac_pow;
use crate::common::constant::*;
use crate::r1::region1::*;
use crate::r1::region1_pT::*;
use crate::r2::region2::*;
use crate::r2::region2_pT::*;
use crate::r3::region3::*;
use crate::r3::region3_v_pT::*;
use crate::r4::region4_sat_pT::*;

/// saturation water include :region3  pT2v_sat_reg
pub fn p2sat_water(p: f64, o_id: i32) -> f64 {
    if o_id == OP {
        return p;
    }

    let T: f64 = T_saturation(p);
    if o_id == OT {
        return T;
    }

    if (p >= P_MIN && p <= Ps_623) {
        return pT_reg1(p, T, o_id);
    } else {
        // region3
        let mut v: f64 = 0.0;
        if p == PC_WATER {
            v = 1.0 / 322.0;
        } else {
            v = pT2v_sat_reg3(p, T, 0.0);
        }

        if (o_id == OV) {
            return v;
        }
        let d: f64 = 1.0 / v;
        if (o_id == OD) {
            return d;
        }

        return Td_reg3(T, d, o_id);
    }
}

/// saturation  steam
pub fn p2sat_steam(p: f64, o_id: i32) -> f64 {
    if o_id == OP {
        return p;
    }

    let T: f64 = T_saturation(p);
    if o_id == OT {
        return T-273.15;
    }

    if (p >= P_MIN && p <= Ps_623) {
        return pT_reg2(p, T, o_id);
    } else {
        //reg3d =ss
        let mut v: f64 = 0.0;
        if p == PC_WATER {
            let v: f64 = 1.0 / 322.0;
        } else {
            let v: f64 = pT2v_sat_reg3(p, T, 1.0);
        }
        if o_id == OV {
            return v;
        }

        let d: f64 = 1.0 / v;
        if o_id == OD {
            return d;
        }
        return Td_reg3(T, d, o_id);
    }
}

pub fn T2sat_water(T: f64, o_id: i32) -> f64 {
    if o_id == OT {
        return T;
    }

    let mut p: f64 = 0.0;
    if T == TC_WATER {
        p = PC_WATER;
    } else {
        p = p_saturation(T);
    }
    if o_id == OP {
        return p;
    }

    if T <= 623.15 {
        return pT_reg1(p, T, o_id);
    } else {
        // region3
        let mut v: f64 = 0.0;
        if T == TC_WATER {
            if o_id == OS {
                return SC_WATER;
            };
            if o_id == OH {
                return HC_WATER;
            };
            v = 1.0 / 322.0;
        } else {
            v = pT2v_sat_reg3(p, T, 0.0);
        }
        if o_id == OV {
            return v;
        };
        let d: f64 = 1.0 / v;
        if o_id == OD {
            return d;
        }

        return Td_reg3(T, d, o_id);
    }
}

pub fn T2sat_steam(T: f64, o_id: i32) -> f64 {
    if o_id == OT {
        return T;
    }

    let mut p: f64 = 0.0;
    if T == TC_WATER {
        p = PC_WATER;
    } else {
        p = p_saturation(T);
    }
    if o_id == OP {
        return p;
    }

    if T <= 623.15 {
        return pT_reg2(p, T, o_id);
    } else {
        //region3
        let mut v: f64 = 0.0;
        if T == TC_WATER {
            if o_id == OS {
                return SC_WATER;
            }
            if o_id == OH {
                return HC_WATER;
            }

            v = 1.0 / 322.0;
            if o_id == OV {
                return v;
            }
        } else {
            v = pT2v_sat_reg3(p, T, 1.0);
        }
        let d: f64 = 1.0 / v;
        if o_id == OD {
            return d;
        }
        return Td_reg3(T, d, o_id);
    }
}

pub fn px_reg4(p: f64, x: f64, o_id: i32) -> f64 {
    let mut r: f64 = 0.0;
    match o_id {
        OP => return p,
        OX => return x,
        OT => {
            return T_saturation(p) - 273.15;
        }
        _ => {
            if x == 0.0 {
                return p2sat_water(p, o_id);
            } else if x == 1.0 {
                return p2sat_steam(p, o_id);
            } else if (x > 0.0 && x < 1.0) {
                let r1: f64 = p2sat_water(p, o_id);
                let r2: f64 = p2sat_steam(p, o_id);
                r = r1 + x * (r2 - r1);
            }
        }
    }
    return r;
}

pub fn Tx_reg4(T: f64, x: f64, o_id: i32) -> f64 {
    if o_id == OT {
        return T;
    };
    if o_id == OX {
        return x;
    };

    if x == 0.0 {
        return T2sat_water(T, o_id);
    }
    if x == 1.0 {
        return T2sat_steam(T, o_id);
    }

    let mut r: f64 = 0.0;
    if (x > 0.0 && x < 1.0) {
        let p: f64 = p_saturation(T);
        if o_id == OP {
            return p;
        };

        let sw: f64 = T2sat_water(T, o_id);
        let ss: f64 = T2sat_steam(T, o_id);
        r = sw + x * (ss - sw);
    }
    return r;
}
