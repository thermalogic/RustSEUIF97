//! Region 5 Basic Equation: the dimensionless Gibbs free energy
//!    <http://www.iapws.org/relguide/IF97-Rev.html>, Eq32-34

use crate::algo::*;
use crate::common::constant::*;

pub const r5Pstar: f64 = 1.0; //[MPa]
pub const r5Tstar: f64 = 1000.0; //[K]

/*The ideal-gas part  of the  dimensionless Gibbs free energy for region 5
   P37 Table 37: Ideal properties for Region 5
   the coefficients and exponents of ideal-gas part  of the dimensionless Gibbs free energy for region 5, Eq. (33)
*/
const Jo: [i32; 6] = [0, 1, -3, -2, -1, 2];

const no: [f64; 6] = [
    -0.13179983674201e+2,
    0.68540841634434e+1,
    -0.24805148933466e-1,
    0.36901534980333,
    -0.31161318213925e+1,
    -0.32961626538917,
];

/// P36 Eq33 - The equation for the ideal-gas part of the dimensionless Gibbs free energy
#[inline(always)]
pub fn gamma0_reg5(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = pi.ln();
    let tau_inv: f64 = 1.0 / tau;
	let tau_inv2: f64 = tau_inv * tau_inv;
	// Jo 0, 1, -3, -2, -1, 2
	result += no[0];
    result += no[1]*tau;
	result += no[2]*tau_inv2*tau_inv;
    result += no[3]*tau_inv2;
    result += no[4]*tau_inv;
    result += no[5]*tau*tau;

    //for i in 0..6 {
    //    result += no[i] * tau.powi(Jo[i]);
    //}
    result
}

///  region 5 38p
#[inline(always)]
pub fn gamma0_pi_reg5(pi: f64) -> f64 {
    return 1.0 / pi;
}

///  region 5 38p
#[inline(always)]
pub fn gamma0_pipi_reg5(pi: f64) -> f64 {
    let pi_inv: f64 = 1.0 / pi;
    return pi_inv*pi_inv;
}

#[inline(always)]
pub fn gamma0_tau_reg5(tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    //   Jo   0, 1, -3, -2, -1, 2
    // Jo-1  -1, 0, -4, -3, -2, 1
    let tau_inv:f64 = 1.0 / tau;
	let tau_inv2 = tau_inv * tau_inv;
    result += no[1];
	result += no[2]*(-3.0)*tau_inv2*tau_inv2;
    result += no[3]*(-2.0)*tau_inv2*tau_inv;
    result += no[4]*(-1.0)*tau_inv2;
    result += no[5]*(2.0)*tau;

   // for i in 0..6 {
    //    result += no[i] * Jo[i] as f64 * tau.powi(Jo[i] - 1);
    //}
    result
}

#[inline(always)]
pub fn gamma0_tautau_reg5(tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    //   Jo      0,  1, -3, -2, -1,  2 
    // Jo-1     -1,  0, -4, -3, -2,  1 
    // Jo-2     -2, -1, -5, -4, -3,  0 
    // Jo*(Jo-1) 0,  0, 12,  6,  2,  2 
    
    let tau_inv:f64 = 1.0 / tau;
    let tau_inv3:f64 =  tau_inv *tau_inv * tau_inv;
    let tau_inv4:f64 = tau_inv3 * tau_inv;
    
    // i=0,1: Jo*(Jo-1) = 0 
    result += no[2] * 12.0 * tau_inv4*tau_inv;    // Jo=-3,  12
    result += no[3] * 6.0 * tau_inv4;             // Jo=-2,  6
    result += no[4] * 2.0 * tau_inv3;             // Jo=-1,  2
    result += no[5] * 2.0;                        // Jo=2,  2

   // for i in 0..6 {
    //    result += no[i] * (Jo[i] * (Jo[i] - 1)) as f64 * tau.powi(Jo[i] - 2);
    //}
    result
}

#[inline(always)]
pub fn gamma0_pitau_reg5() -> f64 {
    return 0.0;
}

//----------------residual part r of the dimensionless Gibbs free energy -----------------

// Table 38. coefficients and exponents of the residual part r of the dimensionless Gibbs free energy for region 5, Eq.(34)
pub const IJn: [(i32, i32, f64); 6] = [
    (1, 1, 0.15736404855259e-2),
    (1, 2, 0.90153761673944e-3),
    (1, 3, -0.50270077677648e-2),
    (2, 3, 0.22440037409485e-5),
    (2, 9, -0.41163275453471e-5),
    (3, 7, 0.37919454822955e-7),
];

#[inline(always)]
pub fn gammar_reg5(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = 0.0;
	// I 1 1 1 2 2 3
	// J 1 2 3 3 9 7
	let pi_2:f64 = pi*pi;
	let tau_2:f64=tau*tau;
	let tau_3:f64=tau_2*tau;
	let tau_6:f64=tau_3*tau_3;
	result += IJn[0].2 * pi*tau;
    result += IJn[1].2 * pi*tau_2;
	result += IJn[2].2 * pi*tau_3;
	result += IJn[3].2 * pi_2*tau_3;
	result += IJn[4].2 * pi_2*tau_6*tau_3;
	result += IJn[5].2 * pi_2*pi*tau_6*tau;
    return result
    //poly_powi(pi, tau, &IJn)
}

// Table 41. The residual part gammar of the dimensionless Gibbs free energy and its derivatives a according to Eq. (34)

/// The residual part gammar of the dimensionless Gibbs free energy
#[inline(always)]
pub fn gammar_pi_reg5(pi: f64, tau: f64) -> f64 {
	let mut result: f64 = 0.0;
	//   I  1 1 1 2 2 3
	// I-1  0 0 0 1 1 2
	//   J  1 2 3 3 9 7
	let tau_2:f64=tau*tau;
	let tau_3:f64=tau_2*tau;
	let tau_6:f64=tau_3*tau_3;
	result += IJn[0].2 *tau;
    result += IJn[1].2 * tau_2;
	result += IJn[2].2 * tau_3;
	result += IJn[3].2 * 2.0 * pi*tau_3;
	result += IJn[4].2 * 2.0 *pi*tau_6*tau_3;
	result += IJn[5].2 * 3.0 *pi*pi*tau_6*tau;
    return result
    //poly_i_powi(pi, tau, &IJn)
}

#[inline(always)]
pub fn gammar_pipi_reg5(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = 0.0;
	//   I   1  1  1  2  2  3
	// I-1   0  0  0  1  1  2
	// I-2  -1 -1 -1  0  0  1
	//   J   1  2  3  3  9  7
	let tau_2:f64=tau*tau;
	let tau_3:f64=tau_2*tau;
	let tau_6:f64=tau_3*tau_3;
	result += IJn[3].2 * 2.0 * tau_3;
    result += IJn[4].2 * 2.0 * tau_6*tau_3;
	result += IJn[5].2 *6.0 *pi*tau_6*tau;
    return result
    //poly_ii_powi(pi, tau, &IJn)
}

#[inline(always)]
pub fn gammar_tau_reg5(pi: f64, tau: f64) -> f64 {
	let mut result: f64 = 0.0;
	//   I   1  1  1  2  2  3
	//   J   1  2  3  3  9  7
	// J-1   0  1  2  2  8  6
	let pi_2:f64 = pi * pi;
    let tau_2:f64 = tau * tau;
	let tau_6:f64 = tau_2 * tau_2 * tau_2;
    
    result += IJn[0].2 * pi;                           // J=1
    result += IJn[1].2 * 2.0 * pi * tau;               // J=2
    result += IJn[2].2 * 3.0 * pi * tau_2;             // J=3
    result += IJn[3].2 * 3.0 * pi_2 * tau_2;           // J=3
    result += IJn[4].2 * 9.0 * pi_2 * tau_6* tau_2;    // J=9
    result += IJn[5].2 * 7.0 * pi_2 * pi * tau_6;       // J=7
    return result
   // poly_j_powi(pi, tau, &IJn)
}

/// region5 39p
#[inline(always)]
pub fn gammar_tautau_reg5(pi: f64, tau: f64) -> f64 {
	let mut result: f64 = 0.0;
	//   I   1  1  1  2  2  3
	//   J   1  2  3  3  9  7
	// J-1   0  1  2  2  8  6
	// J-2  -1  0  1  1  7  5

	let pi_2:f64 = pi * pi;
    let tau_2:f64 = tau * tau;
	let tau_5:f64 = tau_2 * tau_2 * tau;
    // i=0: J=1, J*(J-1)=0 
    result += IJn[1].2 * 2.0 * pi;                           // J=2,  2
    result += IJn[2].2 * 6.0 * pi * tau;                      // J=3,  6
    result += IJn[3].2 * 6.0 * pi_2 * tau;                    // J=3,  6
    result += IJn[4].2 * 72.0 * pi_2 * tau_5*tau_2;           // J=9, 72
    result += IJn[5].2 * 42.0 * pi_2 * pi * tau_5; 
	return result
    // poly_jj_powi(pi, tau, &IJn)
}

#[inline(always)]
pub fn gammar_pitau_reg5(pi: f64, tau: f64) -> f64 {

    let mut result: f64 = 0.0;
	//   I   1  1  1  2  2  3
	// I-1   0  0  0  1  1  2
	//   J   1  2  3  3  9  7
	// J-1   0  1  2  2  8  6
	let tau_2:f64 = tau * tau;
    let tau_6:f64 = tau_2 * tau_2 * tau_2;
    
    result += IJn[0].2;                                 // I=1,J=1,  1
    result += IJn[1].2 * 2.0 * tau;                     // I=1,J=2,  2
    result += IJn[2].2 * 3.0 * tau_2;                   // I=1,J=3,  3
    result += IJn[3].2 * 6.0 * pi * tau_2;              // I=2,J=3,  6
    result += IJn[4].2 * 18.0 * pi * tau_6*tau_2;       // I=2,J=9,  18
	result += IJn[5].2 * 21.0 * pi*pi * tau_6;          // I=3,J=27, 21
	return result
   // poly_ij_powi(pi, tau, &IJn)
}

// improve performance 10% ,but code complexity is increased
#[inline(always)]
pub fn gammar_pi_tau_reg5(pi: f64, tau: f64) -> (f64, f64) {
    let mut item:f64=0.0;
    let mut result_pi: f64 = 0.0;
	let mut result_tau: f64 = 0.0;
	// I 1 1 1 2 2 3
	// J 1 2 3 3 9 7
	let pi_2:f64 = pi*pi;
	let tau_2:f64=tau*tau;
	let tau_3:f64=tau_2*tau;
	let tau_6:f64=tau_3*tau_3;
	item = IJn[0].2 * pi*tau;
    result_pi += item;
    result_tau += item;
    item = IJn[1].2 * pi*tau_2;
    result_pi += item;
    result_tau += 2.0*item;
    item = IJn[2].2 * pi*tau_3;
    result_pi += item;
    result_tau += 3.0*item;
    item = IJn[3].2 * pi_2*tau_3;
	result_pi += 2.0*item;
    result_tau += 3.0*item;
    item = IJn[4].2 * pi_2*tau_6*tau_3;
	result_pi += 2.0*item;
    result_tau += 9.0*item;
    item = IJn[5].2 * pi_2*pi*tau_6*tau;
    result_pi += 3.0*item;
    result_tau += 7.0*item;
    (result_pi/pi, result_tau/tau)
}
