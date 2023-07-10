
//！  Basic Equation of  IAPWS-IF97 Region3
//！  (T,d)-->p,h,u,s,cp,cv,w
//！       T: temperature  K
//！       P: pressure  MPa
//！       v: specific volume m^3/kg
//！       h: specific enthalpy kJ/kg
//！       u: specific internal energy kJ/kg
//！       s: specific entropy  kJ/(kg K)
//！       cp: specific isobaric heat capacity  kJ/(kg K)
//！       cv: specific isochoric heat capacity kJ/(kg K)
//！       w:  speed of sound  m/s

use crate::algo::fast_ipower::sac_pow;
use crate::common::constant::*;
use crate::r3::region3_hfe::*;

/// Fast recursion algorithm of phi and its derivatives 
///          phi_delta, deltatau, deltadelta, tautau
fn recursion_of_phi_derivatives(delta:f64,tau:f64)->(f64,f64,f64,f64) 
{
      let mut item: f64 = IJn[0].n  * sac_pow(delta, IJn[0].I) * sac_pow(tau, IJn[0].J); 
      let mut item_delta:f64= item*IJn[0].I as f64;
      let mut phi_delta: f64 =item_delta;// sum of phi_delta
      let mut phi_deltadelta: f64 =item_delta*(IJn[0].I-1) as f64; // sum of phi_deltadelta
      let mut phi_deltatau: f64 =item_delta*IJn[0].J as f64; // sum of phi_deltatau
      let mut phi_tautau : f64= item * (IJn[0].J *(IJn[0].J-1))as f64; // sum of phi_tau
      
      for k in 1..39 {
          item = IJn[k].n * sac_pow(delta, IJn[k].I) * sac_pow(tau, IJn[k].J); 
          item_delta = item*IJn[k].I as f64;
          phi_delta  += item_delta; 
          phi_deltadelta += item_delta*(IJn[k].I-1) as f64;
          phi_deltatau +=item_delta*IJn[k].J as f64;
          phi_tautau += item * (IJn[k].J *(IJn[k].J-1))as f64
      }
      phi_delta /=delta;
      phi_delta += n1 / delta;

      phi_deltadelta /=(delta*delta);
      phi_deltadelta +=(-n1 / delta / delta);

      phi_deltatau /=(delta*tau);
      phi_tautau /=(tau*tau);
      return (phi_delta,phi_deltadelta,phi_deltatau,phi_tautau);
}

/// pressure in region 3
pub fn Td2p_reg3(T:f64, d:f64)->f64
{
    let tau:f64= TC_WATER / T;
    let delta:f64= d / DC_WATER;
    return 0.001 * d * RGAS_WATER * T * delta * phi_delta_reg3(tau, delta);
}

/// speciphic internal energy in region 3 in kJ/kg
pub fn Td2u_reg3(T:f64, d:f64)->f64
{
    let tau:f64= TC_WATER / T;
    let delta:f64= d / DC_WATER;
    return RGAS_WATER * T * tau * phi_tau_reg3(tau, delta);
}

/// speciphic entropy in region 3 in kJ/(kg K)
pub fn Td2s_reg3(T:f64, d:f64)->f64
{
    let tau:f64= TC_WATER / T;
    let delta:f64= d / DC_WATER;
    //return RGAS_WATER * (tau * phi_tau_reg3(tau, delta) - phi_reg3(tau, delta));

    // common item of  phi and phi_tau
	let mut item: f64 = IJn[0].n  * sac_pow(delta, IJn[0].I) * sac_pow(tau, IJn[0].J); 
	let mut sum_phi: f64 = n1 * delta.ln()+item; // sum of phi
	let mut sum_phi_tau : f64= item * IJn[0].J as f64; // sum of phi_tau
	for k in 1..39 {
		item = IJn[k].n * sac_pow(delta, IJn[k].I) * sac_pow(tau, IJn[k].J); 
		sum_phi += item;
		sum_phi_tau += item * IJn[k].J as f64;
	}
    return RGAS_WATER * (tau * sum_phi_tau/tau - sum_phi);

}

/// speciphic enthalpy in region 3 in kJ/kg
pub fn Td2h_reg3(T:f64, d:f64)->f64
{
    let tau:f64= TC_WATER / T;
    let delta:f64= d / DC_WATER;
    //return RGAS_WATER * T * (tau * phi_tau_reg3(tau, delta) + delta * phi_delta_reg3(tau, delta));

    // common item of  phi_delta and phi_tau
  	let mut item: f64 = IJn[0].n  * sac_pow(delta, IJn[0].I) * sac_pow(tau, IJn[0].J); 
	let mut sum_phi_delta: f64 =item* IJn[0].I as f64; // sum of phi_delta
	let mut sum_phi_tau : f64= item * IJn[0].J as f64; // sum of phi_tau
	for k in 1..39 {
		item = IJn[k].n * sac_pow(delta, IJn[k].I) * sac_pow(tau, IJn[k].J); 
		sum_phi_delta +=item* IJn[k].I as f64;
		sum_phi_tau += item * IJn[k].J as f64;
	}
    sum_phi_delta /=delta;
    sum_phi_delta += n1 / delta;
    return RGAS_WATER * T * (tau * sum_phi_tau/tau + delta * sum_phi_delta);
}

/// speciphic isobaric heat capacity in region 3 in kJ/(kg K)
pub fn Td2cp_reg3(T:f64, d:f64)->f64
{
    let tau:f64= TC_WATER / T;
    let delta:f64= d / DC_WATER;
    
    //let mut a:f64 = delta * (phi_delta_reg3(tau, delta) - tau * phi_deltatau_reg3(tau, delta));
    //a *= a;
    //let b:f64 = delta * (2.0 * phi_delta_reg3(tau, delta) + delta * phi_deltadelta_reg3(tau, delta));
    //return RGAS_WATER * (-tau * tau * phi_tautau_reg3(tau, delta) + a / b);

    let (phi_delta,phi_deltadelta,phi_deltatau, phi_tautau)=recursion_of_phi_derivatives(delta,tau);
     
    let mut a:f64 = delta * (phi_delta - tau * phi_deltatau);
    a *= a;
    let b:f64 = delta * (2.0 * phi_delta + delta * phi_deltadelta);
    return RGAS_WATER * (-tau * tau * phi_tautau + a / b);
}

/// speciphic isochoric heat capacity in region 3  in kJ/(kg K)
pub fn Td2cv_reg3(T:f64, d:f64)->f64
{
    let tau:f64= TC_WATER / T;
    let delta:f64= d / DC_WATER;
    return RGAS_WATER * (-tau * tau * phi_tautau_reg3(tau, delta));
}

/// speed of sound in region 3  in m/s
pub fn Td2w_reg3(T:f64, d:f64)->f64
{
    let tau:f64= TC_WATER / T;
    let delta:f64= d / DC_WATER;

    //let phi_delta:f64=phi_delta_reg3(tau, delta); 
    //let mut a:f64 = delta * phi_delta - delta * tau * phi_deltatau_reg3(tau, delta);
    //a *= a;
    //let temp:f64=1000.0 * RGAS_WATER * T * (2.0 * delta * phi_delta + delta * delta * phi_deltadelta_reg3(tau, delta) - a / (tau * tau * phi_tautau_reg3(tau, delta)));
    //return temp.sqrt();

      
      let (phi_delta,phi_deltadelta,phi_deltatau, phi_tautau)=recursion_of_phi_derivatives(delta,tau);

      let mut a:f64 = delta * phi_delta - delta * tau * phi_deltatau;
      a *= a;
      let temp:f64=1000.0 * RGAS_WATER * T * (2.0 * delta * phi_delta + delta * delta * phi_deltadelta - a / (tau * tau * phi_tautau));
      return temp.sqrt();
}


