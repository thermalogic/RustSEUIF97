//！Region3: The extended Properties
//!
//! Td_ext_reg3(T: f64, d: f64, o_id: i32) -> f64

use crate::algo::*;
use crate::common::*;
use crate::r3::region3_hfe::*;
use crate::r3::*;

/// the extended properties
pub fn Td_ext_reg3(T: f64, d: f64, o_id: i32) -> f64 {
    match o_id {
        OF => Td2f_reg3(T, d),
        OG => Td2g_reg3(T, d),
        OZ => Td2z_reg3(T, d),
        OKS => Td2ks_reg3(T, d),
        OKT => Td2kt_reg3(T, d),
        OEC => Td2ec_reg3(T, d),
        OJTC => Td2joule_reg3(T, d),
        OIJTC => Td2iJTC_reg3(T, d),
        OPC => Td2pc_reg3(T, d),
        ODVDP => Td2dvdpct_reg3(T, d),
        ODPDT => Td2dpdtcv_reg3(T, d),
        ODVDT => Td2dvdtcp_reg3(T, d),
        OE => Td2e_reg3(T, d),
        OBETAP => Td2batap_reg3(T, d),
        OFI => Td2fi_reg3(T, d),
        OFU => Td2fu_reg3(T, d),
        OALFAP => Td2alfap_reg3(T, d),
        _ => INVALID_OUTID as f64,
    }
}

///  the specific Gibbs free energy,
///    g=R*T*(phi+delta* phidelta )
pub fn Td2g_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    let phi: f64 = phi_reg3(delta, tau);
    let phidelta = phi_delta_reg3(delta, tau);
    RGAS_WATER * T * (phi + delta * phidelta)
}

/// the Helmholtz Specific free energy:
pub fn Td2f_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    RGAS_WATER * T * phi_reg3(delta, tau)
}

/// (dv/dp)T
pub fn Td2dvdpct_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    let mut ddeltadpi: f64 = 2.0 * phi_delta_reg3(delta, tau) + delta * phi_deltadelta_reg3(delta, tau);
    ddeltadpi = d * RGAS_WATER * T * ddeltadpi;
    -1000.0 * DC_WATER / d / d / ddeltadpi
}

/// kt: Isothermal compressibility, [1/MPa]
pub fn Td2kt_reg3(T: f64, d: f64) -> f64 {
    -d * Td2dvdpct_reg3(T, d)
}

/// (dv/dt)p
pub fn Td2dvdtcp_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    let ddelta: f64 = phi_delta_reg3(delta, tau);
    let d1: f64 = ddelta - tau * phi_deltatau_reg3(delta, tau);
    let d2: f64 = 2.0 * ddelta + delta * phi_deltadelta_reg3(delta, tau);
    d1 / d2 / T / d
}

/// ec-Isobaric volume expansion coefficient  1/K
pub fn Td2ec_reg3(T: f64, d: f64) -> f64 {
    d * Td2dvdtcp_reg3(T, d)
}

/// joule ： Joule-Thomson coefficient    K/MPa
/// *  (dt/dp)h
pub fn Td2joule_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    let v: f64 = 1.0 / d;
    let ddelta = phi_delta_reg3(delta, tau);
    let ddeltatau = phi_deltatau_reg3(delta, tau);
    let ddeltadelta = phi_deltadelta_reg3(delta, tau);

    let mut cp1: f64 = delta * (ddelta - tau * ddeltatau);
    cp1 *= cp1;
    let cp2 = delta * (2.0 * ddelta + delta * ddeltadelta);
    let cp = RGAS_WATER * (-1.0 * tau * tau * phi_tautau_reg3(delta, tau) + cp1 / cp2);

    let d1 = ddelta - tau * ddeltatau; //dpdtcv
    let d2 = 2.0 * ddelta + delta * ddeltadelta; //dvdpct
    let cex = d1 / d2 / T;
    1000.0 * (v / cp) * (T * cex - 1.0)
}

//  Isothermal throttling coefficient
///  iJTC Isothermal Joule-Thomson coefficient kJ/(kg·MPa)
/// *  (dh/dp)t
pub fn Td2iJTC_reg3(T: f64, d: f64) -> f64 {
    -Td2cp_reg3(T, d) * Td2joule_reg3(T, d) / 1000.0
}

///  (dp/dT)v
pub fn Td2dpdtcv_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    0.001 * RGAS_WATER * d * delta * (phi_delta_reg3(delta, tau) - tau * phi_deltatau_reg3(delta, tau))
}

/// alfap: Relative pressure coefficient  
/// * β=(1.0/p)*(dp/dT)
pub fn Td2pc_reg3(T: f64, d: f64) -> f64 {
    let p: f64 = Td2p_reg3(T, d);
    (1.0 / p) * Td2dpdtcv_reg3(T, d)
}

/// e: Specific exergy    kJ/kg  
pub fn Td2e_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    let phi = phi_reg3(delta, tau);
    let phidelta = phi_delta_reg3(delta, tau);
    let phitau = phi_tau_reg3(delta, tau);
    RGAS_WATER * (T * (phi + delta * phidelta) + (T - Tt) * tau * (phitau - phi))
}

/// z: Compressibility factor  
pub fn Td2z_reg3(T: f64, d: f64) -> f64 {
    let v: f64 = 1.0 / d;
    let p: f64 = Td2p_reg3(T, d);
    1000.0 * p * v / RGAS_WATER / T
}

/// ks: Isentropic exponent
/// * ks= -(v/p)*/1000*(dp/dv)
pub fn Td2ks_reg3(T: f64, d: f64) -> f64 {
    let w: f64 = Td2w_reg3(T, d);
    let p: f64 = Td2p_reg3(T, d);
    1.0E-6 * w * w * d / p
}

/// batap Isothermal stress coefficient, kg/m³
///  (-1.0/p)*(dp/dv)T
///  (dp/dv)T=-1.0/(v*kt)
///  (-1.0/p)*(dp/dv)T=(-1.0/p)*(-1.0/(v*kt))=1.0/(p*v*kt)
///  =d/(p*kt)
pub fn Td2batap_reg3(T: f64, d: f64) -> f64 {
    let kt: f64 = Td2kt_reg3(T, d);
    let p: f64 = Td2p_reg3(T, d);
    d / (p * kt)
}

/// fi Fugacity coefficient  
///  *  fi = exp(g/ R /T)
pub fn Td2fi_reg3(T: f64, d: f64) -> f64 {
    let g: f64 = Td2g_reg3(T, d);
    (g / RGAS_WATER / T).exp()
}

///  fu: Fugacity  MPa
///  * fu =p *fi
pub fn Td2fu_reg3(T: f64, d: f64) -> f64 {
    let fi: f64 = Td2fi_reg3(T, d);
    let p: f64 = Td2p_reg3(T, d);
    p * fi
}

/// alfap - Relative pressure coefficient  1/K
///  * alfap=ec/p/kt
pub fn Td2alfap_reg3(T: f64, d: f64) -> f64 {
    let p: f64 = Td2p_reg3(T, d);
    let ec: f64 = Td2ec_reg3(T, d);
    let kt: f64 = Td2kt_reg3(T, d);
    ec / p / kt
}
