//！Region1: The extended Properties:
//!  * pT_ext_reg1(p: f64, T: f64, o_id: i32) -> f64
//! # Properties(17):
//!
//!  *  ks: Isentropic exponent
//!  *  ec: Isobaric cubic expansion coefficient  1/K
//!  *  kt: Isothermal compressibility, [1/MPa]
//! 
//！ *  e: Specific exergy    kJ/kg
//!  *  z: Compressibility factor   -
//!  *  f: Specific Helmholtz free energy kJ/kg
//!  *  g: Specific Gibbs free energy  kJ/kg
//!  *  joule : Joule-Thomson coefficient    K/MPa
//!  *  iJTC: Isothermal Joule-Thomson coefficient kJ/(kg·MPa)
//!  *  pc:  isochoric pressure coefficient  1/K
//!  *  dpdtcv:  Partial derivative (dP/dT)v  MPa/K
//!  *  dvdtcp: Partial derivative (dV/dT)p  m³/(kg·K)
//!  *  dvdpct: Partial derivative (dV/dP)T  m³/(kg·MPa)
//！ *  batap ：Isothermal stress coefficient, kg/m³
//！ *  fi: Fugacity coefficient
//!  *  fu: Fugacity, MPa
//!  * alfap: relative pressure coefficient  1/K

use crate::common::*;
use crate::r1::region1_gfe::*;
use crate::r1::region1_pT::*;

/// Region1: the extended properties
pub fn pT_ext_reg1(p: f64, T: f64, o_id: i32) -> f64 {
    match o_id {
        OF => pT2f_reg1(p, T),
        OG => pT2g_reg1(p, T),
        OZ => pT2z_reg1(p, T),
        OKS => pT2ks_reg1(p, T),
        OKT => pT2kt_reg1(p, T),
        OEC => pT2ec_reg1(p, T),
        OJTC => pT2joule_reg1(p, T),
        OIJTC => pT2iJTC_reg1(p, T),
        OPC => pT2pc_reg1(p, T),
        ODPDT => pT2dpdtcv_reg1(p, T),
        ODVDP => pT2dvdpct_reg1(p, T),
        ODVDT => pT2dvdTcp_reg1(p, T),
        OE => pT2e_reg1(p, T),
        OBETAP => pT2batap_reg1(p, T),
        OFI => pT2fi_reg1(p, T),
        OFU => pT2fu_reg1(p, T),
        OALFAP => pT2alfap_reg1(p, T),
        _ => INVALID_OUTID as f64,
    }
}
/// the specific Gibbs free energy, Page 6, Eq(7)
/// *  g=R*t*gamma
pub fn pT2g_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    RGAS_WATER * T * gamma_reg1(tau, pi)
}

/// the Helmholtz Specific free energy:
/// *  f=u-T*s=R*T*(gamma-gamma_pi)
pub fn pT2f_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    RGAS_WATER * T * (gamma_reg1(tau, pi) - pi * gamma_pi_reg1(tau, pi))
}

/// kt Isothermal compressibility 1/Mpa
// * kt=-(1.0/V)*(dv/dp)T
pub fn pT2kt_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    -gamma_pipi_reg1(tau, pi) / gamma_pi_reg1(tau, pi) / r1pstar
}

/// ec Isobaric volume expansion coefficient  1/K
/// Coefficient of thermal expansion,
/// * α=(1.0/V)*(dv/dT)p
pub fn pT2ec_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    1.0 / T - tau * gamma_pitau_reg1(tau, pi) / gamma_pi_reg1(tau, pi) / T
}

/// Region 1 - (dp/dt)v  MPa/K
pub fn pT2dpdtcv_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    let gamma_pitau: f64 = gamma_pitau_reg1(tau, pi);
    let gamma_pi: f64 = gamma_pi_reg1(tau, pi);
    let gamma_pipi: f64 = gamma_pipi_reg1(tau, pi);
    r1pstar * (gamma_pitau * r1Tstar - gamma_pi * T) / (T * T * gamma_pipi)
}

/// joule ： Joule-Thomson coefficient    K/MPa
/// *  (dt/dp)h
pub fn pT2joule_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    let g1pi: f64 = gamma_pi_reg1(tau, pi);
    let v: f64 = 0.001 * RGAS_WATER * T * g1pi / r1pstar;
    let cp: f64 = -RGAS_WATER * tau * tau * gamma_tautau_reg1(tau, pi);
    let TCex_1: f64 = -tau * gamma_pitau_reg1(tau, pi) / g1pi;
    (v / cp) * TCex_1
}

//  Isothermal throttling coefficient
///  iJTC Isothermal Joule-Thomson coefficient kJ/(kg·MPa)
/// *  (dh/dp)t
pub fn pT2iJTC_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    0.001 * RGAS_WATER * r1Tstar * gamma_pitau_reg1(tau, pi) / r1pstar
}

/// dvdp Partial derivative (dV/dP)T  m3/(kg·MPa)
pub fn pT2dvdpct_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    0.001 * RGAS_WATER * T * gamma_pipi_reg1(tau, pi) / r1pstar / r1pstar
}

// (dv/dT)p m3/(kg.K)
pub fn pT2dvdTcp_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    0.001 * RGAS_WATER * (gamma_pi_reg1(tau, pi) - tau * gamma_pitau_reg1(tau, pi)) / r1pstar
}

/// z: Compressibility factor  
pub fn pT2z_reg1(p: f64, T: f64) -> f64 {
    let v: f64 = pT2v_reg1(p, T);
    1000.0 * p * v / RGAS_WATER / T
}

/// ks:Isentropic exponent OKIE
/// * ks= -(v/p)*/1000*(dp/dv)
pub fn pT2ks_reg1(p: f64, T: f64) -> f64 {
    let v: f64 = pT2v_reg1(p, T);
    let w: f64 = pT2w_reg1(p, T);
    1.0E-6 * w * w / v / p
}

///  isochoric pressure coefficient 1/K
///   * β=(1.0/p)*(dp/dT)v
pub fn pT2pc_reg1(p: f64, T: f64) -> f64 {
    (1.0 / p) * pT2dpdtcv_reg1(p, T)
}

/// e  the specific exergy   kJ/kg  in region 1
pub fn pT2e_reg1(p: f64, T: f64) -> f64 {
    let tau: f64 = r1Tstar / T;
    let pi: f64 = p / r1pstar;
    let gumma: f64 = gamma_reg1(tau, pi);
    let gumma_tau: f64 = gamma_tau_reg1(tau, pi);
    RGAS_WATER * (T * gumma + (T - 273.16) * (tau * gumma_tau - gumma))
}

/// batap Isothermal stress coefficient, kg/m³
///  (-1.0/p)*(dp/dv)T
///  (dp/dv)T=-1.0/(v*kt)
///  (-1.0/p)*(dp/dv)T=(-1.0/p)*(-1.0/(v*kt))=1.0/(p*v*kt)
pub fn pT2batap_reg1(p: f64, T: f64) -> f64 {
    let v: f64 = pT2v_reg1(p, T);
    let kt: f64 = pT2kt_reg1(p, T);
    1.0 / (p * v * kt)
}

///fi Fugacity coefficient  
/// *  fi = exp(g/ R /T)
pub fn pT2fi_reg1(p: f64, T: f64) -> f64 {
    let g: f64 = pT2g_reg1(p, T);
    (g / RGAS_WATER / T).exp()
}

///  fu: Fugacity  MPa
///  * fu =p *fi
pub fn pT2fu_reg1(p: f64, T: f64) -> f64 {
    let fi: f64 = pT2fi_reg1(p, T);
    p * fi
}

/// alfap - Relative pressure coefficient  1/K
///  * alfap=ec/p/kt
pub fn pT2alfap_reg1(p: f64, T: f64) -> f64 {
    let ec: f64 = pT2ec_reg1(p, T);
    let kt: f64 = pT2kt_reg1(p, T);
    ec/p/kt
}