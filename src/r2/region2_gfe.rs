//! IAPWS-IF97 Basic Equation for Region 2:
//！  
//！ R7-97(2012) August 2007 : http://www.iapws.org/relguide/IF97-Rev.html
//！
//！  The dimensionless Gibbs free energy gamma and its derivatives
//！       The  basic  equation  Eq.(15), P13
//！            The ideal-gas part ： Eq.(16)
//！            The residual part ： Eq.(17)

use crate::algo::fast_ipower::sac_pow;
use crate::common::constant::*;

//  Initialize coefficients and exponents for region 2
//  Table 10 Page 13：The ideal-gas part gamma'o for region 2, Eq.(16)a
pub const J0: [i32; 9] = [0, 1, -5, -4, -3, -2, -1, 2, 3];
pub const n0: [f64; 9] = [
     -0.96927686500217E+01,
     0.10086655968018E+02,
     -0.5608791128302E-02,
     0.71452738081455E-01,
     -0.40710498223928,
     0.14240819171444E+01,
    -0.43839511319450E+01,
    -0.28408632460772,
     0.21268463753307E-01];

// Table 11 Page 13 The residual part gamma'r for  region 2, Eq. (17)
pub const IJn: [IJnData; 43] = [
    IJnData {
        I: 1,
        J: 0,
        n: -0.17731742473213E-02,
    },
    IJnData {
        I: 1,
        J: 1,
        n: -0.17834862292358E-01,
    },
    IJnData {
        I: 1,
        J: 2,
        n: -0.45996013696365E-01,
    },
    IJnData {
        I: 1,
        J: 3,
        n: -0.57581259083432E-01,
    },
    IJnData {
        I: 1,
        J: 6,
        n: -0.5032527872793E-01,
    },
    IJnData {
        I: 2,
        J: 1,
        n: -0.33032641670203E-04,
    },
    IJnData {
        I: 2,
        J: 2,
        n: -0.18948987516315E-03,
    },
    IJnData {
        I: 2,
        J: 4,
        n: -0.39392777243355E-02,
    },
    IJnData {
        I: 2,
        J: 7,
        n: -0.43797295650573E-01,
    },
    IJnData {
        I: 2,
        J: 36,
        n: -2.6674547914087E-05,
    },
    IJnData {
        I: 3,
        J: 0,
        n: 0.20481737692309E-07,
    },
    IJnData {
        I: 3,
        J: 1,
        n: 0.43870667284435E-06,
    },
    IJnData {
        I: 3,
        J: 3,
        n: -0.3227767723857E-04,
    },
    IJnData {
        I: 3,
        J: 6,
        n: -0.15033924542148E-02,
    },
    IJnData {
        I: 3,
        J: 35,
        n: -0.040668253562649,
    },
    IJnData {
        I: 4,
        J: 1,
        n: -7.8847309559367E-10,
    },
    IJnData {
        I: 4,
        J: 2,
        n: 1.2790717852285E-08,
    },
    IJnData {
        I: 4,
        J: 3,
        n: 4.8225372718507E-07,
    },
    IJnData {
        I: 5,
        J: 7,
        n: 2.2922076337661E-06,
    },
    IJnData {
        I: 6,
        J: 3,
        n: -1.6714766451061E-11,
    },
    IJnData {
        I: 6,
        J: 16,
        n: -2.1171472321355E-03,
    },
    IJnData {
        I: 6,
        J: 35,
        n: -23.895741934104,
    },
    IJnData {
        I: 7,
        J: 0,
        n: -5.905956432427E-18,
    },
    IJnData {
        I: 7,
        J: 11,
        n: -1.2621808899101E-06,
    },
    IJnData {
        I: 7,
        J: 25,
        n: -0.038946842435739,
    },
    IJnData {
        I: 8,
        J: 8,
        n: 1.1256211360459E-11,
    },
    IJnData {
        I: 8,
        J: 36,
        n: -8.2311340897998,
    },
    IJnData {
        I: 9,
        J: 13,
        n: 1.9809712802088E-08,
    },
    IJnData {
        I: 10,
        J: 4,
        n: 1.0406965210174E-19,
    },
    IJnData {
        I: 10,
        J: 10,
        n: -1.0234747095929E-13,
    },
    IJnData {
        I: 10,
        J: 14,
        n: -1.0018179379511E-09,
    },
    IJnData {
        I: 16,
        J: 29,
        n: -8.0882908646985E-11,
    },
    IJnData {
        I: 16,
        J: 50,
        n: 0.10693031879409,
    },
    IJnData {
        I: 18,
        J: 57,
        n: -0.33662250574171,
    },
    IJnData {
        I: 20,
        J: 20,
        n: 8.9185845355421E-25,
    },
    IJnData {
        I: 20,
        J: 35,
        n: 3.0629316876232E-13,
    },
    IJnData {
        I: 20,
        J: 48,
        n: -4.2002467698208E-06,
    },
    IJnData {
        I: 21,
        J: 21,
        n: -5.9056029685639E-26,
    },
    IJnData {
        I: 22,
        J: 53,
        n: 3.7826947613457E-06,
    },
    IJnData {
        I: 23,
        J: 39,
        n: -1.2768608934681E-15,
    },
    //
    IJnData {
        I: 24,
        J: 26,
        n: 7.3087610595061E-29,
    },
    IJnData {
        I: 24,
        J: 40,
        n: 5.5414715350778E-17,
    },
    IJnData {
        I: 24,
        J: 58,
        n: -9.436970724121E-07,
    },
];

/// Eq16 P13 Ideal-gas part of fundamental equation for region 2
pub fn gamma0_reg2(tau: f64, pi: f64) -> f64 {
    let mut result = pi.ln(); //
    for i in 0..9 {
        result += n0[i] * sac_pow(tau, J0[i]);
       // result += n0[i] * tau.powi(J0[i]);
    }
    return result;
}

/// First derivative in pi of ideal-gas part of fundamental equation for region 2
pub fn gamma0_pi_reg2(pi: f64) -> f64
{
    return 1.0 / pi;
}

/// Second derivative in pi of ideal-gas part of fundamental equation for region 2
pub fn gamma0_pipi_reg2(pi: f64) -> f64
{
    return -1.0 / pi / pi;
}

/// First derivative in tau of ideal-gas part of fundamental equation for region 2
pub fn gamma0_tau_reg2(tau: f64) -> f64
{
    let mut result: f64 = 0.0;
    for i in 0..9 {  // 为什么需要：0..9 而不是0..8
       result += n0[i] * J0[i] as f64 * sac_pow(tau, J0[i] - 1);
       //result += n0[i] * J0[i] as f64 * tau.powi(J0[i] - 1);
    }
   return result;
}

/// Second derivative in tau of ideal-gas part of fundamental equation for region 2
pub fn gamma0_tautau_reg2(tau: f64, pi: f64) -> f64
{
    let mut result: f64 = 0.0;
    for i in 0..9 {
        result += n0[i] * (J0[i] * (J0[i] - 1)) as f64 * sac_pow(tau, J0[i] - 2);
        //result += n0[i] * (J0[i] * (J0[i] - 1)) as f64 * tau.powi(J0[i] - 2);
    }
    return result;
}

/// Second derivative in pi and tau of ideal-gas part of fundamental equation for region 2
pub fn gamma0_pitau_reg2() -> f64
{
    return 0.0;
}

//   Eq(17), Page 13   Residual part of fundamental equation for region 2
pub fn gammar_reg2(tau: f64, pi: f64) -> f64 {
    let mut result: f64 = 0.0;
    let tau1: f64 = tau - 0.5;
    for k in 0..43 {
        result += IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(tau1, IJn[k].J);
        //result += IJn[k].n * pi.powi(IJn[k].I) *tau1.powi(IJn[k].J);
    }
    return result;
}

/// First derivative in pi of residual part of fundamental equation for region 2
pub fn gammar_pi_reg2(tau: f64, pi: f64) -> f64
{
    let mut result: f64 = 0.0;
    let tau1: f64 = tau - 0.5;
    for k in 0..43 {
        result += IJn[k].n * IJn[k].I as f64 * sac_pow(pi, IJn[k].I - 1) * sac_pow(tau1, IJn[k].J);
    }
    return result;
}

/// Second derivative in pi of residual part of fundamental equation for region 2
pub fn gammar_pipi_reg2(tau: f64, pi: f64) -> f64
{
    let mut result: f64 = 0.0;
    let tau1: f64 = tau - 0.5;
    for k in 0..43 {
        result += IJn[k].n
            * (IJn[k].I * (IJn[k].I - 1)) as f64
            * sac_pow(pi, IJn[k].I - 2)
            * sac_pow(tau1, IJn[k].J);
    }
    return result;
}

/// First derivative in tau of residual part of fundamental equation for region 2
pub fn gammar_tau_reg2(tau: f64, pi: f64) -> f64
{
    let mut result: f64 = 0.0;
    let tau1: f64 = tau - 0.5;
    for k in 0..43 {
        result +=
           IJn[k].n * sac_pow(pi, IJn[k].I) * IJn[k].J as f64 * sac_pow(tau1, IJn[k].J - 1);
        //result +=
        //    IJn[k].n * pi.powi(IJn[k].I) * IJn[k].J as f64 * tau1.powi(IJn[k].J - 1);
    }
    return result;
}

/// Second derivative in tau of residual part of fundamental equation for region 2
pub fn gammar_tautau_reg2(tau: f64, pi: f64) -> f64
{
    let mut result: f64 = 0.0;
    let tau1: f64 = tau - 0.5;
    for k in 0..43 {
        result += IJn[k].n
            * sac_pow(pi, IJn[k].I)
            * (IJn[k].J * (IJn[k].J - 1)) as f64
            * sac_pow(tau1, IJn[k].J - 2);
    }
    return result;
}

/// Second derivative in pi and tau of residual part of fundamental equation for region 2
pub fn gammar_pitau_reg2(tau: f64, pi: f64) -> f64
{
    let mut result: f64 = 0.0;
    let tau1: f64 = tau - 0.5;
    for k in 0..43 {
        result += IJn[k].n
            * IJn[k].I as f64
            * sac_pow(pi, IJn[k].I - 1)
            * IJn[k].J as f64
            * sac_pow(tau1, IJn[k].J - 1);
    }
    return result;
}
