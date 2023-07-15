//！  Region 5 Basic Equation： the dimensionless Gibbs free energy
//！       http://www.iapws.org/relguide/IF97-Rev.html, Eq32-34

use crate::algo::fast_ipower::sac_pow;
use crate::common::constant::*;

/// ideal-gas part  of the  dimensionless Gibbs free energy for region 5 

/*	P37 Table 37: Ideal properties for Region 5
     the coefficients and exponents of ideal-gas part  of the dimensionless Gibbs free energy for region 5, Eq. (33)
*/
const Jo:[i32;6] = [0, 1, -3, -2, -1, 2];

const no:[f64;6] = [-0.13179983674201e+2, 0.68540841634434e+1,-0.24805148933466e-1,
                     0.36901534980333, 	-0.31161318213925e+1, -0.32961626538917];

/// P36 The equation for the ideal-gas part of the dimensionless Gibbs free energy reads eq33
pub fn  gamma0_reg5(pi:f64, tau:f64)->f64
{
	let mut result:f64 =pi.ln();
	for  i in 0..6
	 { result += no[i] * sac_pow(tau, Jo[i]);
	}
	return result;
}

///  region 5 38p
pub fn  gamma0_pi_reg5(pi:f64)->f64
{
	return 1.0 / pi;
}

///  region 5 38p
pub fn  gamma0_pipi_reg5(pi:f64)->f64
{
	return -1.0 / pi / pi;
}

pub fn  gamma0_tau_reg5(tau:f64)->f64
{
	let mut result:f64 = 0.0;
	for  i in 0..6
	{
		result += no[i] * Jo[i] as f64* sac_pow(tau, Jo[i] - 1);
	}
	return result;
}

pub fn  gamma0_tautau_reg5(tau:f64)->f64
{
	let mut result:f64 = 0.0;
	for  i in 0..6
	{
		result += no[i] * (Jo[i] * (Jo[i] - 1)) as f64 * sac_pow(tau, Jo[i] - 2);
	}
	return result;
}

pub fn  gamma0_pitau_reg5()->f64
{
	return 0.0;
}

///----------------residual part r of the dimensionless Gibbs free energy -----------------

// Table 38. coefficients and exponents of the residual part r of the dimensionless Gibbs free energy for region 5, Eq.(34)
pub const Ir:[i32;6] = [1, 1, 1, 2, 2, 3];
pub const Jr:[i32;6] = [1, 2, 3, 3, 9, 7];
pub const nr:[f64;6] = [0.15736404855259e-2, 0.90153761673944e-3, -0.50270077677648e-2,
		    		 0.22440037409485e-5, -0.41163275453471e-5, 0.37919454822955e-7];

pub fn gammar_reg5(pi:f64, tau:f64)->f64
{
	let mut result:f64=0.0;
	for i in 0..6 {
		result += nr[i] * sac_pow(pi, Ir[i]) * sac_pow(tau, Jr[i]);
	}
	return result;
}

// Table 41. The residual part gammar of the dimensionless Gibbs free energy and its derivatives a according to Eq. (34)

/// The residual part gammar of the dimensionless Gibbs free energy
pub fn gammar_pi_reg5(pi:f64, tau:f64)->f64

{
	let mut result:f64=0.0;
	for i in 0..6 {
		result += nr[i] * Ir[i] as f64 * sac_pow(pi, Ir[i] - 1) * sac_pow(tau, Jr[i]);
	}
	return result;
}

pub fn gammar_pipi_reg5(pi:f64, tau:f64)->f64
{
	let mut result:f64=0.0;
	for i in 0..6 {
		result += nr[i] * (Ir[i] * (Ir[i] - 1)) as f64 * sac_pow(pi, Ir[i] - 2) * sac_pow(tau, Jr[i]);
	}
	return result;
}

pub fn gammar_tau_reg5(pi:f64, tau:f64)->f64
{
	let mut result:f64=0.0;
	for i in 0..6 {
		result += nr[i] * sac_pow(pi, Ir[i]) * Jr[i] as f64* sac_pow(tau, Jr[i] - 1);
	}
	return result;
}


/// region5 39p
pub fn gammar_tautau_reg5(pi:f64, tau:f64)->f64
{
	let mut result:f64=0.0;
	for i in 0..6 {
		result += nr[i] * sac_pow(pi, Ir[i]) * (Jr[i] * (Jr[i] - 1)) as f64 * sac_pow(tau, Jr[i] - 2);
	}
	return result;
}

pub fn gammar_pitau_reg5(pi:f64, tau:f64)->f64
{
	let mut result:f64=0.0;
	for i in 0..6 {
		result += nr[i] * Ir[i] as f64 * sac_pow(pi, Ir[i] - 1) * Jr[i] as f64* sac_pow(tau, Jr[i] - 1);
	}
	return result;
}
