//! Region 1 - Baisic Equation:(p,T)-->v, u,h,s,cp, cv,w
//! *  IAPWS-97(2012) August 2007 Page8 Table 3 : <http://www.iapws.org/relguide/IF97-Rev.html>
//! Thermodynamic  properties
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
    0.001 * RGAS_WATER * T * pi * gamma_pi_reg1(tau, pi) / p
}

/// specific internal energy in region 1
pub fn pT2u_reg1(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r1pstar;
    let tau: f64 = r1Tstar / T;
    let (d_pi, d_tau) = sum_i_j_power(7.1 - pi, tau - 1.222, &IJn);
    RGAS_WATER * T * (tau * d_tau - pi * (-d_pi)) ////7.1 - pi1 ,so -d_pi
}

/// specific entropy in region 1
pub fn pT2s_reg1(p: f64, T: f64) -> f64 {
    let tau:f64 = r1Tstar / T;
    let pi :f64= p / r1pstar;
    let (v, d_tau) = sum_0_j_power(7.1 - pi, tau - 1.222, &IJn);
    RGAS_WATER * (tau * d_tau - v)
}

/// specific enthalpy in region 1
pub fn pT2h_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    RGAS_WATER * T * tau * gamma_tau_reg1(tau, pi)
}

/// specific isobaric heat capacity in region 1
pub fn pT2cp_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    -RGAS_WATER * tau * tau * gamma_tautau_reg1(tau, pi)
}

/// specific isochoric heat capacity in region 1
pub fn pT2cv_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    let (d_pi, d_pipi, d_pitau, d_tautau) = sum_i_ii_ij_jj_power(7.1 - pi, tau - 1.222, &IJn);
    let a: f64 = -tau * tau * d_tautau;
    let mut b: f64 = -d_pi - tau * (-d_pitau);// 7.1 - pi1 ,so -d_pi,-d_pitau
    b *= b;
    RGAS_WATER * (a + b / d_pipi)
}

/// speed of sound in region 1:  w in m/s, T in K, p in Mpa
pub fn pT2w_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi = p / r1pstar;
    let (d_pi, d_pipi, d_pitau, d_tautau) = sum_i_ii_ij_jj_power(7.1 - pi, tau - 1.222, &IJn);
    let mut a = -d_pi - tau * (-d_pitau); // 7.1 - pi1,so -d_pi ,-d_pitau
    a *= a;
    let mut b = a / (tau * tau * d_tautau);
    b -= d_pipi;
    let temp = 1000.0 * RGAS_WATER * T / b;
    -d_pi * temp.sqrt() // 7.1 - pi1,so -d_pi
}
