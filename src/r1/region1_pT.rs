//! Region 1 - Baisic Equation:(p,T)-->v, u,h,s,cp, cv,w
//! *  IAPWS-97(2012) August 2007 Page8 Table 3 : <http://www.iapws.org/relguide/IF97-Rev.html>
//! Thermodynamic  properties(9)
//！ *  T: temperature  K
//！ *  P: pressure  MPa
//！ *  v: specific volume m^3/kg
//！ *  h: specific enthalpy kJ/kg
//！ *  u: specific internal energy kJ/kg
//！ *  s: specific entropy  kJ/(kg K)
//！ *  cp: specific isobaric heat capacity  kJ/(kg K)
//！ *  cv: specific isochoric heat capacity kJ/(kg K)
//！ *   w:  speed of sound  m/s

use crate::algo::*;
use crate::common::constant::*;
use crate::r1::region1_gfe::*;

/// specific volume in region 1
pub fn pT2v_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    0.001 * RGAS_WATER * T * pi * gamma_pi_reg1(pi, tau) / p
}

/// specific internal energy in region 1
pub fn pT2u_reg1(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r1pstar;
    let tau: f64 = r1Tstar / T;
    let (d_pi, d_tau) = polys_i_j_powi_reg1(pi, tau);
    RGAS_WATER * T * (tau * d_tau - pi * d_pi)
}

/// specific entropy in region 1
pub fn pT2s_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    let (v, d_tau) = polys_0_j_powi_reg1(pi, tau);
    RGAS_WATER * (tau * d_tau - v)
}

/// specific enthalpy in region 1
pub fn pT2h_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    RGAS_WATER * T * tau * gamma_tau_reg1(pi, tau)
}

/// specific isobaric heat capacity in region 1
pub fn pT2cp_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    -RGAS_WATER * tau * tau * gamma_tautau_reg1(pi, tau)
}

/// specific isochoric heat capacity in region 1
pub fn pT2cv_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;

    let (d_pi, d_pipi, d_pitau, d_tautau) = polys_i_ii_ij_jj_powi_reg1(pi, tau);
    let a: f64 = -tau * tau * d_tautau;
    let mut b: f64 = d_pi - tau * d_pitau;
    b *= b;
    RGAS_WATER * (a + b / d_pipi)
}

/// speed of sound in region 1:  w in m/s, T in K, p in Mpa
/// w= v*(-(dp/dv)s)^0.5
pub fn pT2w_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi = p / r1pstar;

    let (d_pi, d_pipi, d_pitau, d_tautau) = polys_i_ii_ij_jj_powi_reg1(pi, tau);

    let mut a = d_pi - tau * d_pitau;
    a *= a;
    let b = a / (tau * tau * d_tautau) - d_pipi;
    let temp = 1000.0 * RGAS_WATER * T / b;
    d_pi * temp.sqrt()
}
