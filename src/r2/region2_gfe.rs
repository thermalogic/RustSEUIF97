//! IAPWS-IF97 Basic Equation for Region 2:
//！
//！ R7-97(2012) August 2007 : http://www.iapws.org/relguide/IF97-Rev.html
//！
//！ The dimensionless Gibbs free energy gamma and its derivatives
//！ * The  basic  equation  Eq.(15), P13
//！      *  The ideal-gas part： Eq.(16)
//！      *  The residual part： Eq.(17)

use crate::algo::*;
use crate::common::constant::*;
use crate::r2::*;

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
    0.21268463753307E-01,
];

// Table 11 Page 13 The residual part gamma'r for  region 2, Eq. (17)
pub const IJn: [(i32, i32, f64); 43] = [
    (1, 0, -0.17731742473213E-02),
    (1, 1, -0.17834862292358E-01),
    (1, 2, -0.45996013696365E-01),
    (1, 3, -0.57581259083432E-01),
    (1, 6, -0.5032527872793E-01),
    (2, 1, -0.33032641670203E-04),
    (2, 2, -0.18948987516315E-03),
    (2, 4, -0.39392777243355E-02),
    (2, 7, -0.43797295650573E-01),
    (2, 36, -2.6674547914087E-05),
    (3, 0, 0.20481737692309E-07),
    (3, 1, 0.43870667284435E-06),
    (3, 3, -0.3227767723857E-04),
    (3, 6, -0.15033924542148E-02),
    (3, 35, -0.040668253562649),
    (4, 1, -7.8847309559367E-10),
    (4, 2, 1.2790717852285E-08),
    (4, 3, 4.8225372718507E-07),
    (5, 7, 2.2922076337661E-06),
    (6, 3, -1.6714766451061E-11),
    (6, 16, -2.1171472321355E-03),
    (6, 35, -23.895741934104),
    (7, 0, -5.905956432427E-18),
    (7, 11, -1.2621808899101E-06),
    (7, 25, -0.038946842435739),
    (8, 8, 1.1256211360459E-11),
    (8, 36, -8.2311340897998),
    (9, 13, 1.9809712802088E-08),
    (10, 4, 1.0406965210174E-19),
    (10, 10, -1.0234747095929E-13),
    (10, 14, -1.0018179379511E-09),
    (16, 29, -8.0882908646985E-11),
    (16, 50, 0.10693031879409),
    (18, 57, -0.33662250574171),
    (20, 20, 8.9185845355421E-25),
    (20, 35, 3.0629316876232E-13),
    (20, 48, -4.2002467698208E-06),
    (21, 21, -5.9056029685639E-26),
    (22, 53, 3.7826947613457E-06),
    (23, 39, -1.2768608934681E-15),
    //
    (24, 26, 7.3087610595061E-29),
    (24, 40, 5.5414715350778E-17),
    (24, 58, -9.436970724121E-07),
];

/// Eq16 P13 Ideal-gas part of fundamental equation for region 2
pub fn gamma0_reg2(pi: f64, tau: f64) -> f64 {
    let mut result = pi.ln(); //
    for i in 0..9 {
        result += n0[i] * tau.powi(J0[i]);
    }
    result
}

/// First derivative in pi of ideal-gas part of fundamental equation for region 2
pub fn gamma0_pi_reg2(pi: f64) -> f64 {
    1.0 / pi
}

/// Second derivative in pi of ideal-gas part of fundamental equation for region 2
pub fn gamma0_pipi_reg2(pi: f64) -> f64 {
    -1.0 / pi / pi
}

/// First derivative in tau of ideal-gas part of fundamental equation for region 2
pub fn gamma0_tau_reg2(tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    for i in 0..9 {
        result += n0[i] * J0[i] as f64 * tau.powi(J0[i] - 1);
    }
    result
}

/// Second derivative in tau of ideal-gas part of fundamental equation for region 2
pub fn gamma0_tautau_reg2(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    for i in 0..9 {
        result += n0[i] * (J0[i] * (J0[i] - 1)) as f64 * tau.powi(J0[i] - 2);
    }
    result
}

/// Second derivative in pi and tau of ideal-gas part of fundamental equation for region 2
pub fn gamma0_pitau_reg2() -> f64 {
    0.0
}

//   Eq(17), Page 13   Residual part of fundamental equation for region 2
pub fn gammar_reg2(pi: f64, tau: f64) -> f64 {
    let steps: [(usize, usize); 3] = [(0, 19), (19, 38), (38, 43)];
    poly_powi_steps(pi, tau - 0.5, &IJn, &steps)

}

/// First derivative in pi of residual part of fundamental equation for region 2
pub fn gammar_pi_reg2(pi: f64, tau: f64) -> f64 {
    let mut result: f64 = 0.0;
    let tau1: f64 = tau - 0.5;
    let steps: [(usize, usize); 3] = [(0, 16), (16, 32), (32, 43)];
    poly_i_powi_steps(pi, tau - 0.5, &IJn, &steps)
}

/// Second derivative in pi of residual part of fundamental equation for region 2
pub fn gammar_pipi_reg2(pi: f64, tau: f64) -> f64 {
    let steps: [(usize, usize); 3] = [(0, 13), (13, 26), (26, 43)];
    poly_ii_powi_steps(pi, tau - 0.5, &IJn, &steps)
}

/// First derivative in tau of residual part of fundamental equation for region 2
pub fn gammar_tau_reg2(pi: f64, tau: f64) -> f64 {
    let steps: [(usize, usize); 3] = [(0, 13), (13, 26), (26, 43)];
    poly_j_powi_steps(pi, tau - 0.5, &IJn, &steps)
}

/// Second derivative in tau of residual part of fundamental equation for region 2
pub fn gammar_tautau_reg2(pi: f64, tau: f64) -> f64 {
    let steps: [(usize, usize); 3] = [(0, 13), (13, 26), (26, 43)];
    poly_jj_powi_steps(pi, tau - 0.5, &IJn, &steps)
}

/// Second derivative in pi and tau of residual part of fundamental equation for region 2
pub fn gammar_pitau_reg2(pi: f64, tau: f64) -> f64 {
    let steps: [(usize, usize); 3] = [(0, 13), (13, 26), (26, 43)];
    poly_ij_powi_steps(pi, tau - 0.5, &IJn, &steps)
}

// -----------multiple ----------------------------

pub fn polys_0_j_powi_reg2(pi: f64, tau: f64) -> (f64, f64) {
    let steps: [(usize, usize); 3] = [(0, 13), (13, 26), (26, 43)];
    let (gammar, gammar_tau)=polys_0_j_powi_steps(pi, tau - 0.5, &IJn, &steps);
    (gammar, gammar_tau)
}

pub fn polys_i_j_powi_reg2(pi: f64, tau: f64) -> (f64, f64) {
    let steps: [(usize, usize); 3] = [(0, 13), (13, 26), (26, 43)];
    let (gammar, gammar_tau)=polys_i_j_powi_steps(pi, tau - 0.5, &IJn, &steps);
    (gammar, gammar_tau)
}

pub fn polys_i_ii_ij_jj_powi_reg2(pi: f64, tau: f64) -> (f64, f64, f64, f64) {
    let tau1: f64 = tau - 0.5;
    let steps: [(usize, usize); 4] = [(0, 11), (11, 22), (22, 33), (33, 43)];
    let mut item: f64 = 0.0;
    let mut pi_item: f64 = 0.0;

    let mut gammar_pi: f64 = 0.0;
    let mut gammar_pipi: f64 = 0.0;
    let mut gammar_pitau: f64 = 0.0;
    let mut gammar_tautau: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * pi.powi(IJn[k].0) * tau1.powi(IJn[k].1);
            pi_item = IJn[k].0 as f64 * item;
            gammar_pi += pi_item;
            gammar_pipi += (IJn[k].0 - 1) as f64 * pi_item;
            gammar_pitau += IJn[k].1 as f64 * pi_item;
            gammar_tautau += (IJn[k].1 * (IJn[k].1 - 1)) as f64 * item;
        }
    }
    (
        gammar_pi / pi,
        gammar_pipi / pi / pi,
        gammar_pitau / pi / tau1,
        gammar_tautau / tau1 / tau1,
    )
   
}
