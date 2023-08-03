//! Region 5 - Baisic Equation:(p,T)-->v, u,h,s,cp, cv,w
//!
//!    <http://www.iapws.org/relguide/IF97-Rev.html>, Eq 32-34
//!     P39 (p,T).   T,K   p, MPa        

use crate::algo::*;
use crate::common::constant::*;
use crate::r5::region5_gfe::*;

pub fn pT2v_reg5(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r5Pstar;
    let tau: f64 = r5Tstar / T;
    let a: f64 = pi * (gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau));
    (RGAS_WATER * T / p * a) * 0.001
}

pub fn pT2u_reg5(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r5Pstar;
    let tau: f64 = r5Tstar / T;

    let (sum_gammar_pi, sum_gammar_tau) = polys_i_j_powi(pi, tau, &IJn);
    let a: f64 = tau * (gamma0_tau_reg5(tau) + sum_gammar_tau) - pi * (gamma0_pi_reg5(pi) + sum_gammar_pi);
    RGAS_WATER * T * a
}

pub fn pT2s_reg5(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r5Pstar;
    let tau: f64 = r5Tstar / T;
    let (sum_gammar, sum_gammar_tau) = polys_0_j_powi(pi, tau, &IJn);
    let a: f64 = tau * (gamma0_tau_reg5(tau) + sum_gammar_tau) - (gamma0_reg5(pi, tau) + sum_gammar);
    RGAS_WATER * a
}

pub fn pT2h_reg5(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r5Pstar;
    let tau: f64 = r5Tstar / T;
    let a: f64 = tau * (gamma0_tau_reg5(tau) + gammar_tau_reg5(pi, tau));
    RGAS_WATER * T * a
}

pub fn pT2cp_reg5(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r5Pstar;
    let tau: f64 = r5Tstar / T;
    let a: f64 = -tau * tau * (gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau));
    RGAS_WATER * a
}

pub fn pT2cv_reg5(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r5Pstar;
    let tau: f64 = r5Tstar / T;

    let (gammar_pi, gammar_pipi, gammar_pitau, gammar_tautau) = polys_i_ii_ij_jj_powi(pi, tau, &IJn);

    let a: f64 = -tau * tau * (gamma0_tautau_reg5(tau) + gammar_tautau);
    let b: f64 = 1.0 + pi * gammar_pi - tau * pi * gammar_pitau;
    let c: f64 = 1.0 - pi * pi * gammar_pipi;
    RGAS_WATER * (a - b * b) / c
}

pub fn pT2w_reg5(p: f64, T: f64) -> f64 {
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;

    let (gammar_pi, gammar_pipi, gammar_pitau, gammar_tautau) = polys_i_ii_ij_jj_powi(pi, tau, &IJn);

    let temp: f64 = pi * gammar_pi;
    let a: f64 = 1.0 + temp * (2.0 + temp);
    let b: f64 = 1.0 - pi * pi * gammar_pipi;
    let c: f64 = 1.0 + pi * (gammar_pi - tau * gammar_pitau);
    let d: f64 = tau * tau * (gamma0_tautau_reg5(tau) + gammar_tautau);
    let mut w: f64 = RGAS_WATER * T * a / (b + (c * c) / d) * 1000.0;
    if w < 0.0 {
        w = 0.0;
    }
    w.sqrt()
}
