//! Region5  Basic Equation
//!    http://www.iapws.org/relguide/IF97-Rev.html, Eq 32-34
//! P39 (p,T)->v,h,s,cp,cv,w,k 
//！    t[K],p[MPa]
//！             k: isentropic exponent
	
use crate::algo::fast_ipower::sac_pow;
use crate::common::constant::*;
use crate::r5::region5_gfe::*;

const r5Pstar:f64 = 1.0;	//[MPa]
const r5Tstar:f64 = 1000.0; //[K]

pub fn pT2v_reg5(p:f64,T:f64)->f64
{ 
	let pi:f64 = p / r5Pstar;
	let tau:f64 = r5Tstar / T;
	let a:f64 = pi * (gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau));
	return  (RGAS_WATER * T / p * a)*0.001;	
}

pub fn pT2u_reg5(p:f64,T:f64)->f64
{ 
	let pi:f64 = p / r5Pstar;
	let tau:f64 = r5Tstar / T;
	let a:f64 = tau * (gamma0_tau_reg5(tau) + gammar_tau_reg5(pi, tau)) - pi * (gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau));
	return RGAS_WATER * T * a;
}

pub fn pT2s_reg5(p:f64,T:f64)->f64
{ 
	let pi:f64 = p / r5Pstar;
	let tau:f64 = r5Tstar / T;
	let a:f64=tau* (gamma0_tau_reg5(tau) + gammar_tau_reg5(pi, tau)) - (gamma0_reg5(pi, tau) + gammar_reg5(pi, tau));
	return RGAS_WATER * a;	
    /*
    // common item of  gammar and  gammar_tau
	let mut item: f64 = nr[0] * sac_pow(pi, Ir[0]) * sac_pow(tau, Jr[0]); 
	let mut sum_gammar: f64 = item; // sum of gammar
	let mut sum_gammar_tau : f64= item * Jr[0] as f64; // sum of gammar_tau
	for k in 1..6 {
		item = nr[k] * sac_pow(pi, Ir[k]) * sac_pow(tau, Jr[k]); 
		sum_gammar += item;
		sum_gammar_tau += item * Jr[k] as f64;
	}
	let a:f64=tau* (gamma0_tau_reg5(tau) + sum_gammar_tau/tau) - (gamma0_reg5(pi, tau) + sum_gammar);
	return RGAS_WATER * a;	*/
}

pub fn pT2h_reg5(p:f64,T:f64)->f64
{ 
	let pi:f64 = p / r5Pstar;
	let tau:f64 = r5Tstar / T;
	let a:f64=tau* (gamma0_tau_reg5(tau) + gammar_tau_reg5(pi, tau));
	return  RGAS_WATER * T * a;	
}

pub fn pT2cp_reg5(p:f64,T:f64)->f64
{ 
	let pi:f64 = p / r5Pstar;
	let tau:f64 = r5Tstar / T;
	let a:f64 = -tau*tau * (gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau));
	return RGAS_WATER * a;
}

pub fn pT2cv_reg5(p:f64,T:f64)->f64
{ 
	let pi:f64 = p / r5Pstar;
	let  tau:f64 = r5Tstar / T;
	let a:f64 = -tau*tau * (gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau));
	let b:f64 = 1.0 + pi * gammar_pi_reg5(pi, tau) - tau * pi * gammar_pitau_reg5(pi, tau);
	let c:f64 = 1.0 - pi*pi * gammar_pipi_reg5(pi, tau);
	return RGAS_WATER * (a - b*b) / c;	
}

pub fn pT2w_reg5(p:f64,T:f64)->f64
{
	let tau:f64 = r5Tstar / T;
	let pi:f64 = p / r5Pstar;
	
	let dgammar_pi:f64 =gammar_pi_reg5(pi, tau);

	let temp:f64=pi*dgammar_pi;
	let a:f64 = 1.0+temp*(2.0+temp);
	
	let b:f64 = 1.0 - pi * pi * gammar_pipi_reg5(pi, tau);
	
	let c:f64 = 1.0 + pi *(dgammar_pi - tau *  gammar_pitau_reg5(pi, tau));
	
	let d:f64 = tau * tau * (gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau));
	
	let mut w:f64 = RGAS_WATER * T * a / (b + (c * c) / d) * 1000.0;
	if (w < 0.0)
	{	w = 0.0;}
	return w.sqrt();
}

/// isentropic exponent
pub fn pT2k_reg5(p:f64,T:f64)->f64
{  	let tau:f64 = r5Tstar / T;
	let pi:f64 = p / r5Pstar;

	let dgammar_pi:f64 = gammar_pi_reg5(pi, tau);
	let dgammar_pi_1:f64 = 1.0+dgammar_pi;

	let a:f64 = dgammar_pi_1*dgammar_pi_1;
	let b:f64 = 1.0 - pi * pi * gammar_pipi_reg5(pi, tau);
	let c:f64 = 1.0 + pi * dgammar_pi - tau * pi * gammar_pitau_reg5(pi, tau);
	let d:f64 = tau * tau * (gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau));
	let e:f64 = pi * (gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau));
	return a / (b + (c * c) / d) / e;
}
