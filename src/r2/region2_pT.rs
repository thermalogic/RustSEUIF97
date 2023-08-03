//! Region 2 - Baisic Equation:(p,T)-->v, u,h,s,cp, cv,w
//!
//! IAPWS-97(2012) August 2007 : <http://www.iapws.org/relguide/IF97-Rev.html>
//! *  T: temperature K
//! *  P: pressure  MPa
//! *  v: specific volume m^3/kg
//! *  h: specific enthalpy kJ/kg
//! *  u: specific internal energy kJ/kg
//! *  s: specific entropy  kJ/(kg K)
//! *  cp: specific isobaric heat capacity  kJ/(kg K)
//! *  cv: specific isochoric heat capacity kJ/(kg K)
//！*   w:  speed of sound  m/s

use crate::algo::*;
use crate::common::constant::*;
use crate::r2::region2_gfe::*;
use crate::r2::*;

pub const r2Tstar: f64 = 540.0;
pub const r2pstar: f64 = 1.0;

pub fn pT2v_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let gammar_pi: f64 = gammar_pi_reg2(pi, tau);
    0.001 * RGAS_WATER * T * (gamma0_pi_reg2(pi) + gammar_pi)
}

pub fn pT2h_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let gammar_tau = gammar_tau_reg2(pi, tau);
    RGAS_WATER * r2Tstar * (gamma0_tau_reg2(tau) + gammar_tau)
}

pub fn pT2s_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    // ----------- fast recursion  get mutil gamma0's items -------------------
    let mut gamma0: f64 = pi.ln();
    let mut gamma0_tau: f64 = 0.0;
    let mut gamma0_item: f64 = 0.0;
    for i in 0..9 {
        gamma0_item = n0[i] * tau.powi(J0[i]);
        gamma0 += gamma0_item;
        gamma0_tau += gamma0_item * J0[i] as f64;
    }
    gamma0_tau /= tau;

    let (gammar, gammar_tau) = polys_0_j_powi_reg2(pi, tau);
    RGAS_WATER * (tau * (gamma0_tau + gammar_tau) - (gamma0 + gammar))
}

pub fn pT2u_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;

    let (gammar_pi, gammar_tau) = polys_i_j_powi_reg2(pi, tau);
    RGAS_WATER * T * (tau * (gamma0_tau_reg2(tau) + gammar_tau) - pi * (gamma0_pi_reg2(pi) + gammar_pi))
}

pub fn pT2cp_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    -RGAS_WATER * tau * tau * (gamma0_tautau_reg2(pi, tau) + gammar_tautau_reg2(pi, tau))
}

pub fn pT2cv_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;

    let gamma0_tautau: f64 = gamma0_tautau_reg2(pi, tau);

    let (gammar_pi, gammar_pipi, gammar_pitau, gammar_tautau) = polys_i_ii_ij_jj_powi_reg2(pi, tau);

    let a: f64 = 1.0 + pi * (gammar_pi - tau * gammar_pitau);
    let b: f64 = a * a / (1.0 - pi * pi * gammar_pipi);
    let c: f64 = -tau * tau * (gamma0_tautau + gammar_tautau);
    RGAS_WATER * (c - b)
}

pub fn pT2w_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;

    let gamma0_tautau: f64 = gamma0_tautau_reg2(pi, tau);

    let (gammar_pi, gammar_pipi, gammar_pitau, gammar_tautau) = polys_i_ii_ij_jj_powi_reg2(pi, tau);

    let mut a: f64 = pi * gammar_pi;
    a *= a;
    let r1: f64 = 1000.0 * RGAS_WATER * T * (1.0 + 2.0 * pi * gammar_pi + a);
    let mut b: f64 = 1.0 + pi * gammar_pi - tau * pi * gammar_pitau;
    b *= b;
    let r2: f64 = 1.0 - pi * pi * gammar_pipi + b / (tau * tau * (gamma0_tautau + gammar_tautau));
    let r: f64 = r1 / r2;
    r.sqrt()
}
