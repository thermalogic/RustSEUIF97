//! Region 4 - Saturation  water and steam：
//！ *  p_saturation(T), T_saturation(p)
//!  *  p, MPa   T,K

use crate::algo::*;
use crate::common::constant::*;

//  Initialize coefficients for region 4
const n: [f64; 10] = [
    0.11670521452767E+04,
    -0.72421316703206E+06,
    -0.17073846940092E+02,
    0.12020824702470E+05,
    -0.32325550322333E+07,
    0.14915108613530E+02,
    -0.48232657361591E+04,
    0.40511340542057E+06,
    -0.23855557567849E+00,
    0.65017534844798E+03,
];

/// saturation pressure, MPa
pub fn p_saturation(T: f64) -> f64 {
    let mut ps: f64 = 0.0;
    if T < 273.15 || T > 647.096
    // TC_WATER=647.096
    {
        ps = INVALID_VALUE as f64;
    } else {
        let del: f64 = T + n[8] / (T - n[9]);
        let aco: f64 = del * (del + n[0]) + n[1];
        let bco: f64 = del * (n[2] * del + n[3]) + n[4];
        let cco: f64 = del * (n[5] * del + n[6]) + n[7];
        let a: f64 = 2.0 * cco / (-bco + (bco * bco - 4.0 * aco * cco).sqrt());
        ps = a.powi(4);
    }
    ps
}

/// saturation temperature K
pub fn T_saturation(p: f64) -> f64 {
    let mut Ts: f64 = 0.0;
    if p < 0.000611212677 || p > 22.064 {
        Ts = INVALID_VALUE as f64;
    } else {
        let bet: f64 = p.powf(0.25);
        let eco: f64 = bet * (bet + n[2]) + n[5];
        let fco: f64 = bet * (n[0] * bet + n[3]) + n[6];
        let gco: f64 = bet * (n[1] * bet + n[4]) + n[7];
        let temp: f64 = fco * fco - 4.0 * eco * gco;
        let dco: f64 = 2.0 * gco / (-fco - temp.sqrt());
        let temp: f64 = (n[9] + dco) * (n[9] + dco) - 4.0 * (n[8] + n[9] * dco);
        Ts = 0.5 * (n[9] + dco - temp.sqrt());
    }
    Ts
}
