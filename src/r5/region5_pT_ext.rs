//！ Region5: The extended Properties
//! * pT_ext_reg5(p: f64, T: f64, o_id: i32) -> f64

use crate::algo::*;
use crate::common::*;
use crate::r5::region5_gfe::*;
use crate::r5::region5_pT::*;

// Region 5:  the extended properties
pub fn pT_ext_reg5(p: f64, T: f64, o_id: i32) -> f64 {
    match o_id {
        OF => pT2f_reg5(p, T),
        OG => pT2g_reg5(p, T),
        OZ => pT2z_reg5(p, T),
        OKS => pT2ks_reg5(p, T),
        OKT => pT2kt_reg5(p, T),
        OEC => pT2ec_reg5(p, T),
        OJTC => pT2joule_reg5(p, T),
        OIJTC => pT2iJTC_reg5(p, T),
        OPC => pT2pc_reg5(p, T),
        ODPDT => pT2dpdtcv_reg5(p, T),
        ODVDT => pT2dvdtcp_reg5(p, T),
        ODVDP => pT2dvdpct_reg5(p, T),
        OE => pT2e_reg5(p, T),
        OBETAP => pT2batap_reg5(p, T),
        OFI => pT2fi_reg5(p, T),
        OFU => pT2fu_reg5(p, T),
        OALFAP => pT2alfap_reg5(p, T),
        _ => INVALID_OUTID as f64,
    }
}

///  the specific Gibbs free energy, Page 6, Eq(7)
///    g=R*t*gamma
pub fn pT2g_reg5(p: f64, T: f64) -> f64 {
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;
    RGAS_WATER * T * (gamma0_reg5(pi, tau) + gammar_reg5(pi, tau))
}

/// the Helmholtz Specific free energy:
///  f=u-T*s=R*T(gamma-gamma_pi)
pub fn pT2f_reg5(p: f64, T: f64) -> f64 {
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;
    RGAS_WATER
        * T
        * (gamma0_reg5(pi, tau) + gammar_reg5(pi, tau)
            - p * (gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau)))
}

/// kt: Isothermal compressibility, [1/MPa]
/// * kt=-(1.0/V)*(dv/dp)T
pub fn pT2kt_reg5(p: f64, T: f64) -> f64 {
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;
    let gummapipi: f64 = gamma0_pipi_reg5(pi) + gammar_pipi_reg5(pi, tau);
    let gummapi: f64 = gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau);
    -(gummapipi / gummapi)
}

/// dvdt, cp m3/(kg.K)
pub fn pT2dvdtcp_reg5(p: f64, T: f64) -> f64 {
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;
    let gummapi: f64 = gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau);
    let gummapitau: f64 = gamma0_pitau_reg5() + gammar_pitau_reg5(pi, tau);
    0.001 * RGAS_WATER * (gummapi - tau * gummapitau)
}

// coefficient of thermal expansion,
/// ec-Isobaric volume expansion coefficient  1/K
/// α=(1.0/V)*(dv/dT)p
pub fn pT2ec_reg5(p: f64, T: f64) -> f64 {
    let d: f64 = 1.0 / pT2v_reg5(p, T);
    d * pT2dvdtcp_reg5(p, T)
}

/// (dv/dp)T  m3/MPa
pub fn pT2dvdpct_reg5(p: f64, T: f64) -> f64 {
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;
    let gummapipi = gamma0_pipi_reg5(pi) + gammar_pipi_reg5(pi, tau);
    0.001 * RGAS_WATER * T * gummapipi
}

/// joule ： Joule-Thomson coefficient    K/MPa
/// *  (dt/dp)h
pub fn pT2joule_reg5(p: f64, T: f64) -> f64 {
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;
    let gummapi = gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau);
    let gummatautau: f64 = gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau);
    let gummapitau: f64 = gamma0_pitau_reg5() + gammar_pitau_reg5(pi, tau);
    let v = 0.001 * RGAS_WATER * T * gummapi;
    let cp = -RGAS_WATER * tau * tau * gummatautau;
    let TCex_1 = -tau * gummapitau / gummapi;
    1000.0 * (v / cp) * TCex_1
}

//  Isothermal throttling coefficient
///  iJTC Isothermal Joule-Thomson coefficient kJ/(kg·MPa)
/// *  (dh/dp)t
pub fn pT2iJTC_reg5(p: f64, T: f64) -> f64 {
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;
    let gummapitau = gamma0_pitau_reg5() + gammar_pitau_reg5(pi, tau);
    0.001 * RGAS_WATER * r5Tstar * gummapitau
}

/// (dp/dt)v MPa/K
pub fn pT2dpdtcv_reg5(p: f64, T: f64) -> f64 {
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;
    let gummapi = gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau);
    let gummapitau = gamma0_pitau_reg5() + gammar_pitau_reg5(pi, tau);
    let gummapipi = gamma0_pipi_reg5(pi) + gammar_pipi_reg5(pi, tau);
    r5Pstar * (gummapitau * r5Tstar - gummapi * T) / (T * T * gummapipi)
}

///  isochoric pressure coefficient, [1/K]
///   * β=(1.0/p)*(dp/dT)v
pub fn pT2pc_reg5(p: f64, T: f64) -> f64 {
    (1.0 / p) * pT2dpdtcv_reg5(p, T)
}

pub fn pT2e_reg5(p: f64, T: f64) -> f64 {
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;
    let gumma: f64 = gamma0_reg5(pi, tau) + gammar_reg5(pi, tau);
    let gummatau = gamma0_tau_reg5(tau) + gammar_tau_reg5(pi, tau);
    RGAS_WATER * (T * gumma + (T - Tt) * (tau * gummatau - gumma))
}

/// z: Compressibility factor  
pub fn pT2z_reg5(p: f64, T: f64) -> f64 {
    let v: f64 = pT2v_reg5(p, T);
    1000.0 * p * v / RGAS_WATER / T
}

/// ks: Isentropic exponent
/// * ks= -(v/p)*/1000*(dp/dv)
pub fn pT2ks_reg5(p: f64, T: f64) -> f64 {
    let v: f64 = pT2v_reg5(p, T);
    let w: f64 = pT2w_reg5(p, T);
    1.0E-6 * w * w / v / p
}

/// batap Isothermal stress coefficient, kg/m³
///  (-1.0/p)*(dp/dv)T
///  (dp/dv)T=-1.0/(v*kt)
///  (-1.0/p)*(dp/dv)T=(-1.0/p)*(-1.0/(v*kt))=1.0/(p*v*kt)
pub fn pT2batap_reg5(p: f64, T: f64) -> f64 {
    let v: f64 = pT2v_reg5(p, T);
    let kt: f64 = pT2kt_reg5(p, T);
    1.0 / (p * v * kt)
}

///fi Fugacity coefficient  
/// * fi = exp((g-g0)/ R /T)
pub fn pT2fi_reg5(p: f64, T: f64) -> f64 {
    let g: f64 = pT2g_reg5(p, T);
    let tau: f64 = r5Tstar / T;
    let pi: f64 = p / r5Pstar;
    let g0: f64 = T * RGAS_WATER * gamma0_reg5(pi, tau);
    ((g - g0) / RGAS_WATER / T).exp()
}

/// fu: Fugacity  MPa
///   *  fu=p *fi
pub fn pT2fu_reg5(p: f64, T: f64) -> f64 {
    let fi: f64 = pT2fi_reg5(p, T);
    p * fi
}

/// alfap - Relative pressure coefficient  1/K
///  * alfap=ec/p/kt
pub fn pT2alfap_reg5(p: f64, T: f64) -> f64 {
    let ec: f64 = pT2ec_reg5(p, T);
    let kt: f64 = pT2kt_reg5(p, T);
    ec / p / kt
}
