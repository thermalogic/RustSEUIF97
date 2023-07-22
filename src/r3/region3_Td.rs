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

/// Fast recursion algorithm of phi and its derivatives
///          phi_delta, deltatau, deltadelta, tautau
fn recursion_of_phi_derivatives(delta: f64, tau: f64) -> (f64, f64, f64, f64) {
    let mut phi_delta: f64 = n1 / delta;
    let mut phi_deltadelta = -n1 / delta / delta;

    let (sub_phi_delta, sub_phi_deltadelta, phi_deltatau, phi_tautau) =
        sum_i_ii_ij_jj_power(delta, tau, &IJn);
    phi_delta += sub_phi_delta;
    phi_deltadelta += sub_phi_deltadelta;
    (phi_delta, phi_deltadelta, phi_deltatau, phi_tautau)
}

/// pressure in region 3
pub fn Td2p_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    0.001 * d * RGAS_WATER * T * delta * phi_delta_reg3(tau, delta)
}

/// speciphic internal energy in region 3 in kJ/kg
pub fn Td2u_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    RGAS_WATER * T * tau * phi_tau_reg3(tau, delta)
}

/// speciphic entropy in region 3 in kJ/(kg K)
pub fn Td2s_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
  
    let mut sum_phi: f64 = n1 * delta.ln();
    let (sub_sum_phi, sum_phi_tau) = sum_0_j_power(delta, tau, &IJn);
    sum_phi += sub_sum_phi;
    RGAS_WATER * (tau * sum_phi_tau - sum_phi)
}

/// speciphic enthalpy in region 3 in kJ/kg
pub fn Td2h_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
 
    let mut sum_phi_delta: f64 = n1 / delta;
    let (sub_sum_phi_delta, sum_phi_tau) = sum_i_j_power(delta, tau, &IJn);
    sum_phi_delta += sub_sum_phi_delta;
    RGAS_WATER * T * (tau * sum_phi_tau + delta * sum_phi_delta)
}

/// speciphic isobaric heat capacity in region 3 in kJ/(kg K)
pub fn Td2cp_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;

    let (phi_delta, phi_deltadelta, phi_deltatau, phi_tautau) =
        recursion_of_phi_derivatives(delta, tau);

    let mut a: f64 = delta * (phi_delta - tau * phi_deltatau);
    a *= a;
    let b: f64 = delta * (2.0 * phi_delta + delta * phi_deltadelta);
    RGAS_WATER * (-tau * tau * phi_tautau + a / b)
}

/// speciphic isochoric heat capacity in region 3  in kJ/(kg K)
pub fn Td2cv_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;
    RGAS_WATER * (-tau * tau * phi_tautau_reg3(tau, delta))
}

/// speed of sound in region 3  in m/s
pub fn Td2w_reg3(T: f64, d: f64) -> f64 {
    let tau: f64 = TC_WATER / T;
    let delta: f64 = d / DC_WATER;

    let (phi_delta, phi_deltadelta, phi_deltatau, phi_tautau) =
        recursion_of_phi_derivatives(delta, tau);

    let mut a: f64 = delta * phi_delta - delta * tau * phi_deltatau;
    a *= a;
    let temp: f64 = 1000.0
        * RGAS_WATER
        * T
        * (2.0 * delta * phi_delta + delta * delta * phi_deltadelta - a / (tau * tau * phi_tautau));
    temp.sqrt()
}
