//! Region 4 - The extended input pairs: (p,v)->x, (t,v),(t,h),(t,s) ->x
//! * x: Steam quality

use crate::common::constant::*;
use crate::common::propertry_id::*;

use crate::r4::region4_pTx::*;

/// Region 4: (p,v)-> x
pub fn pv2x_reg4(p: f64, v: f64) -> f64 {
    let v_sat_w: f64 = p2sat_water(p, OV);
    let v_sat_s: f64 = p2sat_steam(p, OV);
    (v - v_sat_w) / (v_sat_s - v_sat_w)
}

/// Region 4: (T,v)-> x
pub fn Tv2x_reg4(T: f64, v: f64) -> f64 {
    let v_sat_w: f64 = T2sat_water(T, OV);
    let v_sat_s: f64 = T2sat_steam(T, OV);
    (v - v_sat_w) / (v_sat_s - v_sat_w)
}

/// Region 4: (T,h)-> x
pub fn Th2x_reg4(T: f64, h: f64) -> f64 {
    let h_sat_w: f64 = T2sat_water(T, OH);
    let h_sat_s: f64 = T2sat_steam(T, OH);
    (h - h_sat_w) / (h_sat_s - h_sat_w)
}

/// Region 4: (T,s)-> x
pub fn Ts2x_reg4(T: f64, s: f64) -> f64 {
    let s_sat_w: f64 = T2sat_water(T, OS);
    let s_sat_s: f64 = T2sat_steam(T, OS);
    (s - s_sat_w) / (s_sat_s - s_sat_w)
}
