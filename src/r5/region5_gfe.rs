//！ Region 5 Basic Equation： the dimensionless Gibbs free energy
//！    <http://www.iapws.org/relguide/IF97-Rev.html>, Eq32-34

use crate::algo::*;
use crate::common::constant::*;

pub const r5Pstar: f64 = 1.0; //[MPa]
pub const r5Tstar: f64 = 1000.0; //[K]

/// ideal-gas part  of the  dimensionless Gibbs free energy for region 5

//	P37 Table 37: Ideal properties for Region 5
//     the coefficients and exponents of ideal-gas part  of the dimensionless Gibbs free energy for region 5, Eq. (33)

const Jo: [i32; 6] = [0, 1, -3, -2, -1, 2];

const no: [f64; 6] = [
    -0.13179983674201e+2,
    0.68540841634434e+1,
    -0.24805148933466e-1,
    0.36901534980333,
    -0.31161318213925e+1,
    -0.32961626538917,
];

/// P36 The equation for the ideal-gas part of the dimensionless Gibbs free energy reads eq33
pub fn gamma0_reg5(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = pi.ln();
    for i in 0..6 {
        result += no[i] * tau.powi(Jo[i]);
    }
    result
}

///  region 5 38p
pub fn gamma0_pi_reg5(pi: f64) -> f64 {
    return 1.0 / pi;
}

///  region 5 38p
pub fn gamma0_pipi_reg5(pi: f64) -> f64 {
    return -1.0 / pi / pi;
}

pub fn gamma0_tau_reg5(tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    for i in 0..6 {
        result += no[i] * Jo[i] as f64 * tau.powi(Jo[i] - 1);
    }
    result
}

pub fn gamma0_tautau_reg5(tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    for i in 0..6 {
        result += no[i] * (Jo[i] * (Jo[i] - 1)) as f64 *tau.powi(Jo[i] - 2);
    }
    result
}

pub fn gamma0_pitau_reg5() -> f64 {
    return 0.0;
}

// Table 38. coefficients and exponents of the residual part r of the dimensionless Gibbs free energy for region 5, Eq.(34)
pub const IJn: [(i32, i32, f64); 6] = [
    (1, 1, 0.15736404855259e-2),
    (1, 2, 0.90153761673944e-3),
    (1, 3, -0.50270077677648e-2),
    (2, 3, 0.22440037409485e-5),
    (2, 9, -0.41163275453471e-5),
    (3, 7, 0.37919454822955e-7),
];

pub fn gammar_reg5(pi: f64, tau: f64) -> f64 {
      poly_powi(pi,tau,&IJn)
}

// Table 41. The residual part gammar of the dimensionless Gibbs free energy and its derivatives a according to Eq. (34)

/// The residual part gammar of the dimensionless Gibbs free energy
pub fn gammar_pi_reg5(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    for i in IJn {
        result += i.2 * i.0 as f64 * pi.powi(i.0 - 1) * tau.powi(i.1);
    }
    result
}

pub fn gammar_pipi_reg5(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = 0.0;
     for i in IJn {
         result += i.2 * (i.0 * (i.0 - 1)) as f64 *pi.powi(i.0 - 2) * tau.powi(i.1);
    }
    result   
}

pub fn gammar_tau_reg5(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    for i in IJn {
       result += i.2 * pi.powi(i.0) * i.1 as f64 * tau.powi(i.1 - 1);
    }
    result
}

/// region5 39p
pub fn gammar_tautau_reg5(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    for i in IJn {
        result += i.2 * pi.powi(i.0) * (i.1 * (i.1 - 1)) as f64 *tau.powi(i.1 - 2);
    }
    result
   
}

pub fn gammar_pitau_reg5(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    for i in IJn {
        result += i.2 * i.0 as f64 * pi.powi(i.0 - 1) * i.1 as f64 * tau.powi(i.1 - 1);
    }
    result
   
}
