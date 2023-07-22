//！Region2: The extended Properties:
//! pT_ext_reg2(p: f64, T: f64, o_id: i32) -> f64

use crate::algo::*;
use crate::common::*;
use crate::r2::region2_gfe::*;
use crate::r2::region2_pT::*;

/// the extended properties
pub fn pT_ext_reg2(p: f64, T: f64, o_id: i32) -> f64 {
    match o_id {
        OF => pT2f_reg2(p, T),
        OG => pT2g_reg2(p, T),
        OZ => pT2z_reg2(p, T),
        OKS => pT2ks_reg2(p, T),
        OKT => pT2kt_reg2(p, T),
        OEC => pT2ec_reg2(p, T),
        OJTC => pT2joule_reg2(p, T),
        OIJTC => pT2iJTC_reg2(p, T),
        OPC => pT2pc_reg2(p, T),
        ODPDT => pT2dpdtcv_reg2(p, T),
        ODVDP => pT2dvdpct_reg2(p, T),
        ODVDT => pT2dvdtcp_reg2(p, T),
        OE => pT2e_reg2(p, T),
        OBETAP => pT2batap_reg2(p, T),
        OFI => pT2fi_reg2(p, T),
        OFU => pT2fu_reg2(p, T),
        _ => INVALID_OUTID as f64,
    }
}

/// the specific Gibbs free energy
/// *  g=R*t*(gamma0+gammar)
pub fn pT2g_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    RGAS_WATER * T * (gamma0_reg2(pi,tau) + gammar_reg2(pi,tau))
}

/// the Helmholtz Specific free energy:
/// *  f=u-T*s=R*T*(gamma-gamma_pi)
pub fn pT2f_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    RGAS_WATER
        * T
        * (gamma0_reg2(pi,tau) + gammar_reg2(pi,tau)
            - p * (gamma0_pi_reg2(pi) + gammar_pi_reg2(pi,tau)))
}

/// kt Isothermal compressibility 1/MPa
pub fn pT2kt_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let gummapipi: f64 = gamma0_pipi_reg2(pi) + gammar_pipi_reg2(pi,tau);
    let gummapi: f64 = gamma0_pi_reg2(pi) + gammar_pi_reg2(pi,tau);
    -(gummapipi / gummapi)
}

///(dv/dp)t
pub fn pT2dvdpct_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let gummapipi: f64 = gamma0_pipi_reg2(pi) + gammar_pipi_reg2(pi,tau);
    0.001 * RGAS_WATER * T * gummapipi
}

///(dv/dt)p m3/(kg.K)
pub fn pT2dvdtcp_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let gummapi: f64 = gamma0_pi_reg2(pi) + gammar_pi_reg2(pi,tau);
    let gummapitau: f64 = gamma0_pitau_reg2() + gammar_pitau_reg2(pi,tau);
    0.001 * RGAS_WATER * (gummapi - tau * gummapitau)
}

/// region2 ec Isobaric volume expansion coefficient  1/K
pub fn pT2ec_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let d: f64 = 1.0 / pT2v_reg2(p, T);
    d * pT2dvdtcp_reg2(p, T)
}

/// joule ： Joule-Thomson coefficient    K/MPa
/// *  (dt/dp)h
pub fn pT2joule_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let gummapi = gamma0_pi_reg2(pi) + gammar_pi_reg2(pi,tau);
    let gummatautau = gamma0_tautau_reg2(pi,tau) + gammar_tautau_reg2(pi,tau);
    let gummapitau = gamma0_pitau_reg2() + gammar_pitau_reg2(pi,tau);
    let v = RGAS_WATER * T * gummapi;
    let cp = -RGAS_WATER * tau * tau * gummatautau;
    let TCex_1 = -tau * gummapitau / gummapi;
    (v / cp) * TCex_1
}

//  Isothermal throttling coefficient
///  iJTC Isothermal Joule-Thomson coefficient kJ/(kg·MPa)
/// *  (dh/dp)t
pub fn pT2iJTC_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let gummapitau = gamma0_pitau_reg2() + gammar_pitau_reg2(pi,tau);
    0.001 * RGAS_WATER * r2Tstar * gummapitau
}

/// (dp/dt)v MPa/K
///IAPWS, Revised Advisory Note No. 3: Thermodynamic Derivatives from IAPWS
// Formulations,
///* <http://www.iapws.org/relguide/Advise3.pdf>
///
pub fn pT2dpdtcv_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let gummapi = gamma0_pi_reg2(pi) + gammar_pi_reg2(pi,tau);
    let gummapitau = gamma0_pitau_reg2() + gammar_pitau_reg2(pi,tau);
    let gummapipi = gamma0_pipi_reg2(pi) + gammar_pipi_reg2(pi,tau);
    1.0 * (gummapitau * r2Tstar - gummapi * T) / (T * T * gummapipi)
}

///  isochoric pressure coefficient, [1/K]
///   * β=(1.0/p)*(∂p/∂T|v)
pub fn pT2pc_reg2(p: f64, T: f64) -> f64 {
    (1.0 / p) * pT2dpdtcv_reg2(p, T)
}
///  e  Specific exergy    kJ/kg
pub fn pT2e_reg2(p: f64, T: f64) -> f64 {
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let gumma = gamma0_reg2(pi,tau) + gammar_pi_reg2(pi,tau);
    let gummatau = gamma0_tau_reg2(tau) + gammar_tau_reg2(pi,tau);
    RGAS_WATER * (T * gumma + (T - Tt) * (tau * gummatau - gumma))
}

/// z: Compressibility factor  
pub fn pT2z_reg2(p: f64, T: f64) -> f64 {
    let v: f64 = pT2v_reg2(p, T);
    1000.0 * p * v / RGAS_WATER / T
}

/// ks: Isentropic exponent
/// * ks= -(v/p)*/1000*(dp/dv)
pub fn pT2ks_reg2(p: f64, T: f64) -> f64 {
    let v: f64 = pT2v_reg2(p, T);
    let w: f64 = pT2w_reg2(p, T);
    1.0E-6 * w * w / v / p
}

/// batap Isothermal stress coefficient, kg/m³
///  (-1.0/p)*(dp/dv)T
///  (dp/dv)T=-1.0/(v*kt)
///  (-1.0/p)*(dp/dv)T=(-1.0/p)*(-1.0/(v*kt))=1.0/(p*v*kt)
pub fn pT2batap_reg2(p: f64, T: f64) -> f64 {
    let v: f64 = pT2v_reg2(p, T);
    let kt: f64 = pT2kt_reg2(p, T);
    1.0 / (p * v * kt)
}

/// fi Fugacity coefficient  
/// *  fi = exp((g-g0)/ R /T)
pub fn pT2fi_reg2(p: f64, T: f64) -> f64 {
    let g: f64 = pT2g_reg2(p, T);
    let tau: f64 = r2Tstar / T;
    let pi: f64 = p;
    let g0: f64 = T * RGAS_WATER * gamma0_reg2(pi,tau);
    ((g - g0) / RGAS_WATER / T).exp()
}

/// fu: Fugacity  MPa
///  * fu =p *fi
pub fn pT2fu_reg2(p: f64, T: f64) -> f64 {
    let fi: f64 = pT2fi_reg2(p, T);
    p * fi
}
