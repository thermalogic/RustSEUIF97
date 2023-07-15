//! Region 2 - Baisic Equation:(p,T)-->v, u,h,s,cp, cv,w 
//!  
//! IAPWS-97(2012) August 2007 : <http://www.iapws.org/relguide/IF97-Rev.html>
//! 
//! * T: temperature K
//! * P: pressure  Map
//! * v: specific volume m^3/kg
//! * h: specific enthalpy kJ/kg
//! * u: specific internal energy kJ/kg
//! * s: specific entropy  kJ/(kg K)
//! * cp: specific isobaric heat capacity  kJ/(kg K)
//! * cv: specific isochoric heat capacity kJ/(kg K)
//! * w:  speed of sound  m/s

use crate::algo::fast_ipower::sac_pow;
use crate::common::constant::*;
use crate::r2::region2_gfe::*;

const r2Tstar: f64 = 540.0; 
const r2pstar: f64 = 1.0;   

pub fn pT2v_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    return 0.001 * RGAS_WATER * T * pi * (gamma0_pi_reg2(pi) + gammar_pi_reg2(tau, pi)) / p;
}

pub fn pT2h_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    return RGAS_WATER * T * tau * (gamma0_tau_reg2(tau) + gammar_tau_reg2(tau, pi));
}

pub fn pT2s_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    //  return RGAS_WATER * (tau * (gamma0_tau_reg2(tau) + gammar_tau_reg2(tau, pi)) - (gamma0_reg2(tau, pi) + gammar_reg2(tau, pi)));

    // --------------- fast to get mutil gamma's items -------------------
    // gamma0
    let mut gamma0: f64 = pi.ln();
    let mut gamma0_tau: f64 = 0.0;
    let mut gamma0_item: f64 = 0.0;
    for i in 0..9 {
        gamma0_item = n0[i] * sac_pow(tau, J0[i]);
        gamma0 += gamma0_item;
        gamma0_tau += gamma0_item * J0[i] as f64;
    }

    // gammar
    let mut gammar: f64 = 0.0;
    let mut gammar_tau: f64 = 0.0;
    let mut gammar_item: f64 = 0.0;
    let tau1: f64 = tau - 0.5;
    for k in 0..43 {
        gammar_item = IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(tau1, IJn[k].J);
        gammar += gammar_item;
        gammar_tau += gammar_item * IJn[k].J as f64;
    }
    return RGAS_WATER * (tau * (gamma0_tau / tau + gammar_tau / tau1) - (gamma0 + gammar));
}

pub fn pT2u_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    //return RGAS_WATER * T * (tau * (gamma0_tau_reg2(tau) + gammar_tau_reg2(tau, pi)) - pi * (gamma0_pi_reg2(pi) + gammar_pi_reg2(tau, pi)));

    /* ---   fast to get mutil gammar's items  -----------------------------------------------*/
    //  gamma0_pi : 1.0 / pi; only one item of  gamma0: gamma0_tau_reg2(tau), non fast
    //  gammar:  gammar_tau, gammar_pi
    let mut gammar_item: f64 = 0.0;
    let mut gammar_pi: f64 = 0.0;
    let mut gammar_tau: f64 = 0.0;
    let tau1: f64 = tau - 0.5;
    for k in 0..43 {
        gammar_item = IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(tau1, IJn[k].J);
        gammar_pi += IJn[k].I as f64 * gammar_item;
        gammar_tau += IJn[k].J as f64 * gammar_item
    }
    return RGAS_WATER
        * T
        * (tau * (gamma0_tau_reg2(tau) + gammar_tau / tau1)
            - pi * (gamma0_pi_reg2(pi) + gammar_pi / pi));
}

pub fn pT2cp_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    return -RGAS_WATER * tau * tau * (gamma0_tautau_reg2(tau, pi) + gammar_tautau_reg2(tau, pi));
}

pub fn pT2cv_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    //  gamma0_tautau: only one item ,non fast
    let gamma0_tautau: f64 = gamma0_tautau_reg2(tau, pi);

    //let a:f64 = 1.0 + pi *(gammar_pi_reg2(tau, pi) - tau *  gammar_pitau_reg2(tau, pi));
    //let b:f64=a * a / (1.0 - pi * pi * gammar_pipi_reg2(tau, pi));
    //let c:f64= -tau * tau * (gamma0_tautau + gammar_tautau_reg2(tau, pi));
    //return RGAS_WATER *(c-b);

    // --------------------- fast to get mutil gammar's items ----------------------------------
    //  gammar: gammar_tautau, gammar_pi, gammar_pipi,gammar_pitau

    let mut gammar_item: f64 = 0.0;
    let mut gammar_tautau: f64 = 0.0;
    let mut gammar_pi_item: f64 = 0.0;
    let mut gammar_pi: f64 = 0.0;
    let mut gammar_pipi: f64 = 0.0;
    let mut gammar_pitau: f64 = 0.0;

    let tau1: f64 = tau - 0.5;
    for k in 0..43 {
        gammar_item = IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(tau1, IJn[k].J);

        gammar_tautau += (IJn[k].J * (IJn[k].J - 1)) as f64 * gammar_item;

        gammar_pi_item = IJn[k].I as f64 * gammar_item;
        gammar_pi += gammar_pi_item;
        gammar_pipi += (IJn[k].I - 1) as f64 * gammar_pi_item;
        gammar_pitau += IJn[k].J as f64 * gammar_pi_item;
    }

    gammar_pi /= pi;
    gammar_pitau /= pi * tau1;
    gammar_pipi /= pi * pi;
    gammar_tautau /= tau1 * tau1;

    let a: f64 = 1.0 + pi * (gammar_pi - tau * gammar_pitau);
    let b: f64 = a * a / (1.0 - pi * pi * gammar_pipi);
    let c: f64 = -tau * tau * (gamma0_tautau + gammar_tautau);
    return RGAS_WATER * (c - b);
}

pub fn pT2w_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    //  gamma0_tautau: only one item ,non fast
    let gamma0_tautau: f64 = gamma0_tautau_reg2(tau, pi);

    /*
       let mut gammar_pi:f64 =  gammar_pi_reg2(tau, pi);
       let mut a:f64 = pi * gammar_pi;
       a *= a;
       let temp1:f64=1000.0 * RGAS_WATER * T * (1.0 + 2.0 * pi * gammar_pi + a);
       let mut b:f64 = 1.0 + pi * gammar_pi - tau * pi * gammar_pitau_reg2(tau, pi);
       b *= b;
       let temp2:f64=(1.0 - pi * pi * gammar_pipi_reg2(tau, pi)) + b / (tau * tau * (gamma0_tautau + gammar_tautau_reg2(tau, pi)));
       let temp:f64=temp1/temp2;
       return temp.sqrt();
    */

    // ---   fast to get mutil gammar's items  -----------------------------------------------
    //  gammar: gammar_tautau, gammar_pi, gammar_pipi,gammar_pitau
    let mut gammar_item: f64 = 0.0;
    let mut gammar_tautau: f64 = 0.0;
    let mut gammar_pi_item: f64 = 0.0;
    let mut gammar_pi: f64 = 0.0;
    let mut gammar_pipi: f64 = 0.0;
    let mut gammar_pitau: f64 = 0.0;

    let tau1: f64 = tau - 0.5;
    for k in 0..43 {
        gammar_item = IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(tau1, IJn[k].J);

        gammar_tautau += (IJn[k].J * (IJn[k].J - 1)) as f64 * gammar_item;

        gammar_pi_item = IJn[k].I as f64 * gammar_item;
        gammar_pi += gammar_pi_item;
        gammar_pipi += (IJn[k].I - 1) as f64 * gammar_pi_item;
        gammar_pitau += IJn[k].J as f64 * gammar_pi_item;
    }

    gammar_pi /= pi;
    gammar_pitau /= pi * tau1;
    gammar_pipi /= pi * pi;
    gammar_tautau /= tau1 * tau1;
    //println!("gammar_pi {}  {}", gammar_pi_reg2(tau, pi), gammar_pi);
    //println!("gammar_pitau_reg2 {}  {}",gammar_pitau_reg2(tau, pi),gammar_pitau);
    //println!("gammar_pipi_reg2 {}  {}",gammar_pipi_reg2(tau, pi),gammar_pipi);
    //println!("gammar_tautau_reg2 {}  {}",gammar_tautau_reg2(tau, pi),gammar_tautau);
    let mut a: f64 = pi * gammar_pi;
    a *= a;
    let temp1: f64 = 1000.0 * RGAS_WATER * T * (1.0 + 2.0 * pi * gammar_pi + a);
    let mut b: f64 = 1.0 + pi * gammar_pi - tau * pi * gammar_pitau;
    b *= b;
    let temp2: f64 =
        1.0 - pi * pi * gammar_pipi + b / (tau * tau * (gamma0_tautau + gammar_tautau));
    let temp: f64 = temp1 / temp2;
    return temp.sqrt();
}
