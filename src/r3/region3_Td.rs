//! Region 3 - Basic Equation (T,d)-> p,h,u,s,cp,cv,w
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
use crate::r3::region3_hfe::*;

/// pressure in region 3
pub fn Td2p_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    let phi_delta =phi_delta_reg3(delta, tau);
    0.001 * d * RGAS_WATER * T * delta * phi_delta
}

/// speciphic internal energy in region 3 in kJ/kg
pub fn Td2u_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    RGAS_WATER * T * tau * phi_tau_reg3(delta, tau)
}

/// speciphic entropy in region 3 in kJ/(kg K)
pub fn Td2s_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    let (poly_phi, poly_tau) = polys_0_j_powi_reg3(delta, tau);
    RGAS_WATER * (tau * poly_tau - poly_phi)
}

/// speciphic enthalpy in region 3 in kJ/kg
pub fn Td2h_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    let (poly_delta, poly_tau) = polys_i_j_powi_reg3(delta, tau);
    RGAS_WATER * T * (tau * poly_tau + delta * poly_delta)
}

/// speciphic isobaric heat capacity in region 3 in kJ/(kg K)
pub fn Td2cp_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    let (phi_delta, phi_deltadelta, phi_deltatau, phi_tautau) =
        polys_i_ii_ij_jj_powi_reg3(delta, tau);

    let mut a: f64 = delta * (phi_delta - tau * phi_deltatau);
    a *= a;
    let b: f64 = delta * (2.0 * phi_delta + delta * phi_deltadelta);
    RGAS_WATER * (-tau * tau * phi_tautau + a / b)
}

/// speciphic isochoric heat capacity in region 3  in kJ/(kg K)
pub fn Td2cv_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    RGAS_WATER * (-tau * tau * phi_tautau_reg3(delta, tau))
}

/// speed of sound in region 3  in m/s
pub fn Td2w_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    let (phi_delta, phi_deltadelta, phi_deltatau, phi_tautau) =
        polys_i_ii_ij_jj_powi_reg3(delta, tau);

    let mut a: f64 = delta * phi_delta - delta * tau * phi_deltatau;
    a *= a;
    let temp: f64 = 1000.0
        * RGAS_WATER
        * T
        * (2.0 * delta * phi_delta + delta * delta * phi_deltadelta - a / (tau * tau * phi_tautau));
    temp.sqrt()
}
