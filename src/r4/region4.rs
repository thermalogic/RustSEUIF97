//! Region 4
//! * (p,h) (p,s) (h,s) (p,x) (t,x)
//! * (p,v),(t,v),(t,h),(t,s)
//! * (h,x),(s,x)

use crate::common::constant::*;
use crate::common::property_id::*;
use crate::r3::*;
use crate::r4::*;

pub fn p_sat(t: f64) -> f64 {
    p_saturation(t + K)
}

pub fn t_sat(p: f64) -> f64 {
    T_saturation(p) - K
}

/// for TC_WATER && p == PC_WATER
#[inline(always)]
pub fn pT_reg4(p: f64, T: f64, o_id: i32) -> f64 {
    if T == TC_WATER && p == PC_WATER {
        return Td_reg3(TC_WATER, DC_WATER, o_id);
    } else {
        return INVALID_PT as f64 ;// return 4.0 as f64;
    }
}

#[inline(always)]
pub fn ph_reg4(p: f64, h: f64, o_id: i32) -> f64 {
    let h1: f64 = p2sat_water(p, OH);
    let h2: f64 = p2sat_steam(p, OH);
    let x: f64 = (h - h1) / (h2 - h1);
    if o_id == OX {
        return x;
    }
     px_reg4(p, x, o_id)   
}

#[inline(always)]
pub fn ps_reg4(p: f64, s: f64, o_id: i32) -> f64 {
    let s1: f64 = p2sat_water(p, OS);
    let s2: f64 = p2sat_steam(p, OS);
    let x: f64 = (s - s1) / (s2 - s1);
    if o_id == OX {
        return x;
    }
    px_reg4(p, x, o_id)
   
}

#[inline(always)]
pub fn hs_reg4(h: f64, s: f64, o_id: i32) -> f64 {
    // for T<623.15 only
    let T: f64 = hs2T_reg4(h, s);
    if o_id == OT {
        return T - K;
    }
    let p: f64 = p_saturation(T);
    if o_id == OP {
        return p;
    }

    let h1: f64 = p2sat_water(p, OH);
    let h2: f64 = p2sat_steam(p, OH);

    let x: f64 = (h - h1) / (h2 - h1);
    if o_id == OX {
        return x;
    }
    // T is k
   // if o_id == OT {
    //    px_reg4(p, x, o_id) - 273.15
    //} else {
    //    px_reg4(p, x, o_id)
    //}
    px_reg4(p, x, o_id) 
}

/// Region 4 - The extended input pair
///       (p,v) ->x
///       (t,h),(t,s),(t,v)->x
///   x: Steam quality

///  (p,v) ->x
#[inline(always)]
pub fn pv_reg4(p: f64, v: f64, o_id: i32) -> f64 {
    let x: f64 = pv2x_reg4(p, v);
    if o_id == OX {
        return x;
    }
    px_reg4(p, x, o_id)
}

///  (t,v)
#[inline(always)]
pub fn tv_reg4(t: f64, v: f64, o_id: i32) -> f64 {
    let x: f64 = Tv2x_reg4(t + 273.15, v);
    if o_id == OX {
        return x;
    }
    Tx_reg4(t + 273.15, x, o_id)
}

///  (t,h)
#[inline(always)]
pub fn th_reg4(t: f64, h: f64, o_id: i32) -> f64 {
    let x: f64 = Th2x_reg4(t + 273.15, h);
    if o_id == OX {
        return x;
    }
    Tx_reg4(t + 273.15, x, o_id)
}

///  (t,s)
#[inline(always)]
pub fn ts_reg4(t: f64, s: f64, o_id: i32) -> f64 {
    let x: f64 = Ts2x_reg4(t + 273.15, s);
    if o_id == OX {
        return x;
    }
    Tx_reg4(t + 273.15, x, o_id)
}

#[inline(always)]
pub fn hx_reg4(h: f64, x: f64, o_id: i32) -> f64 {
    let mut Tl: f64 = T_MIN4;
    let mut Tr: f64 = T_MAX4;
    let mut T: f64 = 0.5 * (Tl + Tr);
    for _ in 0..200 {
        let pl = p_saturation(Tl);
        let hl1 = p2sat_water(pl, OH);
        let hl2 = p2sat_steam(pl, OH);
        let xl = (h - hl1) / (hl2 - hl1);

        let pr = p_saturation(Tr);
        let hr1 = p2sat_water(pr, OH);
        let hr2 = p2sat_steam(pr, OH);
        let xr = (h - hr1) / (hr2 - hr1);

        let p = p_saturation(T);
        let h1 = p2sat_water(p, OH);
        let h2 = p2sat_steam(p, OH);
        let x_cal = (h - h1) / (h2 - h1);

        if x_cal > x {
            Tl = T;
        } else {
            Tr = T;
        }
        T = 0.5 * (Tl + Tr);
        if (Tr - Tl) < 1.0e-8 {
            break;
        }
    }
    if o_id == OT {
        return T - 273.15;
    }
    Tx_reg4(T, x, o_id)
}

#[inline(always)]
pub fn sx_reg4(s: f64, x: f64, o_id: i32) -> f64 {
    let mut Tl: f64 = T_MIN4;
    let mut Tr: f64 = T_MAX4;
    let mut T: f64 = 0.5 * (Tl + Tr);
    for _ in 0..200 {
        let pl = p_saturation(Tl);
        let sl1 = p2sat_water(pl, OS);
        let sl2 = p2sat_steam(pl, OS);
        let xl = (s - sl1) / (sl2 - sl1);

        let pr = p_saturation(Tr);
        let sr1 = p2sat_water(pr, OS);
        let sr2 = p2sat_steam(pr, OS);
        let xr = (s - sr1) / (sr2 - sr1);

        let p = p_saturation(T);
        let s1 = p2sat_water(p, OS);
        let s2 = p2sat_steam(p, OS);
        let x_cal = (s - s1) / (s2 - s1);

        if x_cal > x {
            Tl = T;
        } else {
            Tr = T;
        }
        T = 0.5 * (Tl + Tr);
        if (Tr - Tl) < 1.0e-8 {
            break;
        }
    }
    if o_id == OT {
        return T - 273.15;
    }
    Tx_reg4(T, x, o_id)
}
