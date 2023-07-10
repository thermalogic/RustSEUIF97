//！  Basic Equation of  IAPWS-IF97 Region1
//！       (p,T)-->v, u,h,s,cp, cv,w
//！       T: temperature  K
//！       P: pressure  Map
//！       v: specific volume m^3/kg
//！       h: specific enthalpy kJ/kg
//！       u: specific internal energy kJ/kg
//！       s: specific entropy  kJ/(kg K)
//！       cp: specific isobaric heat capacity  kJ/(kg K)
//！       cv: specific isochoric heat capacity kJ/(kg K)
//！       w:  speed of sound  m/s

use crate::algo::fast_ipower::sac_pow;
use crate::common::constant::*;
use crate::r1::region1_gfe::*;

const r1pstar: f64 = 16.53; // MPa
const r1Tstar: f64 = 1386.0; // K

/// specific volume in region 1
pub fn pT2v_reg1(p: f64, T: f64) -> f64
{
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    return 0.001 * RGAS_WATER * T * pi * gamma_pi_reg1(tau, pi) / p;
}

/// specific internal energy in region 1
pub fn pT2u_reg1(p: f64, T: f64) -> f64
{
   /*
    let pi: f64 = p / r1pstar;
    let tau: f64 = r1Tstar / T;
    return RGAS_WATER * T * (tau * gamma_tau_reg1(tau, pi) - pi * gamma_pi_reg1(tau, pi));
  */
     let pi1: f64 = p / r1pstar;
     let tau1: f64 = r1Tstar / T;
 
     let pi: f64 = 7.1 - pi1;
     let tau: f64 = tau1 - 1.222;
     
     let mut item: f64 = IJn[0].n * sac_pow(pi, IJn[0].I) * sac_pow(tau, IJn[0].J);
     let mut d_pi: f64  = item * IJn[0].I as f64;
     let mut d_tau : f64 = item * IJn[0].J as f64; // n[k] * J , result of gammatau_reg1
     for k in 1..34 {
         item = IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(tau, IJn[k].J); 
         d_pi -= item * IJn[k].I as f64;  // note： - minus because：pi = 7.1 - pi， so have 
         d_tau += item * IJn[k].J as f64;
     }
     return RGAS_WATER *T* (tau1 * d_tau/tau - pi1*d_pi/pi);    
}

/// specific entropy in region 1
pub fn pT2s_reg1(p: f64, T: f64) -> f64
{
    /*
    let tau:f64 = r1Tstar / T;
    let pi :f64= p / r1pstar;
    return RGAS_WATER * (tau * gamma_tau_reg1(tau, pi) - gamma_reg1(tau, pi));
    */
    // -------------- derative  -----------------------------
    let pi1: f64 = p / r1pstar;
    let tau1: f64 = r1Tstar / T;

    let pi: f64 = 7.1 - pi1;
    let tau: f64 = tau1 - 1.222;

    // common item of  gamma_reg1 and  gammatau_reg1
    let mut item: f64 = IJn[0].n * sac_pow(pi, IJn[0].I) * sac_pow(tau, IJn[0].J); // n[k] * sac_pow1(pi, i[k]) * sac_pow2(tau, j[k]);
    let mut v: f64 = item; // result of gamma_reg1
    let mut d_tau = item * IJn[0].J as f64; // n[k] * J , result of gammatau_reg1
    for k in 1..34 {
        item = IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(tau, IJn[k].J); // n[k] * sac_pow1(pi, i[k]) * sac_pow2(tau, j[k]);
        v += item;
        d_tau += item * IJn[k].J as f64;
    }
    return RGAS_WATER * (tau1 * d_tau / tau -v);
}

/// specific enthalpy in region 1
pub fn pT2h_reg1(p: f64, T: f64) -> f64
{
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    return RGAS_WATER * T * tau * gamma_tau_reg1(tau, pi);
}

/// specific isobaric heat capacity in region 1
pub fn pT2cp_reg1(p: f64, T: f64) -> f64
{
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    return -RGAS_WATER * tau * tau * gamma_tautau_reg1(tau, pi);
}

/// specific isochoric heat capacity in region 1
pub fn pT2cv_reg1(p: f64, T: f64) -> f64
{
  /*
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    let a = -tau * tau * gamma_tautau_reg1(tau, pi);
    let mut b = gamma_pi_reg1(tau, pi) - tau * gamma_pitau_reg1(tau, pi);
    b *= b;
    return RGAS_WATER * (a + b / gamma_pipi_reg1(tau, pi));
  */

 // -------------- derative  -----------------------------
  let  pi1:f64 = p / r1pstar;
  let  tau1:f64 = r1Tstar / T;
  let  pi:f64 = 7.1 - pi1;
  let  tau:f64 = tau1 - 1.222;

  let mut item:f64;
  let mut d_piitem:f64;
  let mut d_pi:f64 = 0.0;
  let mut d_pitau:f64 = 0.0;
  let mut d_pipi:f64 = 0.0;
  let mut d_tautau:f64 = 0.0;
  for k in 0..34 
  {
    item = IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(tau, IJn[k].J);
    d_piitem = item * IJn[k].I as f64;
    d_pi -= d_piitem; // note：pi = 7.1 - pi， so -minus
    d_pipi += d_piitem * (IJn[k].I - 1) as f64;
    d_pitau -= d_piitem * IJn[k].J as f64;  // because：pi = 7.1 - pi， so have -minus
    d_tautau += item * IJn[k].J as f64 * (IJn[k].J - 1) as f64;
  }

  let a:f64 = -tau1 * tau1 * d_tautau / tau / tau;
  let mut b:f64 = (d_pi - tau1 * d_pitau / tau) / pi;
  b *= b;
  return  RGAS_WATER * (a + b / (d_pipi / pi / pi));
}

/// speed of sound in region 1:  w in m/s, T in K, p in Mpa
pub fn pT2w_reg1(p: f64, T: f64) -> f64
{
    /*let tau: f64 = r1Tstar / T;
    let pi = p / r1pstar;
    let gamma_pi = gamma_pi_reg1(tau, pi);
    let mut a = gamma_pi - tau * gamma_pitau_reg1(tau, pi);
    a *= a;
    let mut b = a / (tau * tau * gamma_tautau_reg1(tau, pi));
    b = b - gamma_pipi_reg1(tau, pi);
    let temp = 1000.0 * RGAS_WATER * T / b;
    return gamma_pi * temp.sqrt();
    */
 // -------------- derative  -----------------------------
    let pi1:f64 = p / r1pstar;
    let tau1:f64 = r1Tstar / T;
    let pi:f64 = 7.1 - pi1;
    let tau:f64 = tau1 - 1.222;
  
    let mut item:f64;
    let mut d_piitem:f64;
    let mut d_pi:f64 = 0.0;
    let mut d_pitau:f64 = 0.0;
    let mut d_pipi:f64 = 0.0;
    let mut d_tautau:f64 = 0.0;
    for k in 0..34 
    {
      item = IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(tau, IJn[k].J);
      d_piitem = item * IJn[k].I as f64;
      
      d_pi -= d_piitem;// note：pi = 7.1 - pi， so - minus
      d_pipi += d_piitem * (IJn[k].I - 1) as f64;
      d_pitau -= d_piitem * IJn[k].J as f64;  // // note：pi = 7.1 - pi， so have - minus
      d_tautau += item * IJn[k].J  as f64* (IJn[k].J - 1) as f64;
    }
    let mut a = (d_pi - tau1 * d_pitau / tau) / pi;
    a *= a;
    let mut b = a / (tau1 * tau1 * d_tautau / tau / tau);
    b = b - d_pipi / pi / pi;
    let temp=1000.0 *  RGAS_WATER * T / b;
    return d_pi * temp.sqrt() / pi; 
   
}
