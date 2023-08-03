//! Transport and Further Properties(7)
//! *  Dynamic viscosity    Pa.s      dv    
//! *  Kinematic viscosity  m^2/s     kv    
//! *  Thermal conductivity W/(m.K)   tc     
//! *  Thermal diffusivity  m^2/s     td     
//! *  Prandtl number                 pr     
//! *  Static dialectric sonstant     sdc
//！*  Surface tension       N/m       st

use crate::algo::*;
use crate::common::constant::*;

/// Prandtl number=dv*cp/tc
/// * dv: Dynamic viscosity Pa.s
/// * cp: specific isobaric heat capacity
/// * tc: Thermal conductivity   W/(m.K)
pub fn prandtl_number(dv: f64, cp: f64, tc: f64) -> f64 {
    1.0E+3 * dv * cp / tc
}

/// Thermal diffusivity  m²/s
/// * Thermal diffusivity  = 1.0e-3*tc /(cp*d)
///    * cp: specific isobaric heat capacity
///    * tc: Thermal conductivity   W/(m.K)
///    * d  : Density kg/m³
pub fn thermal_diffusivity(tc: f64, cp: f64, d: f64) -> f64 {
    1.0e-3 * tc / (cp * d)
}

/// Viscosity  kg/(m·s)
/// * Parameters
///    * rho :  Density  kg/m³
///    * T :   Temperature K
/// IAPWS, Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary  Water Substance
///         <http://www.iapws.org/relguide/viscosity.html>
pub fn viscosity(rho: f64, T: f64) -> f64 {
    let Tr: f64 = T / TC_WATER;
    let Dr: f64 = rho / DC_WATER;

    let no: [f64; 4] = [1.67752, 2.20462, 0.6366564, -0.241605];
    let mut suma: f64 = 0.0;
    for i in 0..4 {
        suma += no[i] / Tr.powi(i as i32)
    }
    let fi0: f64 = 100.0 * Tr.sqrt() / suma;

    // Page 5 Table2 Coefficients Hij
    const ijH: [(i32, i32, f64); 21] = [
        (0, 0, 5.20094e-1),
        (1, 0, 8.50895e-2),
        (2, 0, -1.08374),
        (3, 0, -2.89555e-1),
        (0, 1, 2.22531e-1),
        (1, 1, 9.99115e-1),
        (2, 1, 1.88797),
        (3, 1, 1.26613),
        (5, 1, 1.20573e-1),
        (0, 2, -2.81378e-1),
        (1, 2, -9.06851e-1),
        (2, 2, -7.72479e-1),
        (3, 2, -4.89837e-1),
        (4, 2, -2.57040e-1),
        (0, 3, 1.61913e-1),
        (1, 3, 2.57399e-1),
        (0, 4, -3.25372e-2),
        (3, 4, 6.98452e-2),
        (4, 5, 8.72102e-3),
        (3, 6, -4.35673e-3),
        (5, 6, -5.93264e-4),
    ];

    suma = poly_powi(1.0 / Tr - 1.0, Dr - 1.0, &ijH);
    let fi1: f64 = (Dr * suma).exp();
    let fi2: f64 = 1.0;
    fi0 * fi1 * fi2 * 1.0e-6
}

/// Thermal conductivity  λ W/mK
///   * rho :  Density kg/m³
///   * T :  Temperature K
/// IAPWS, Release on the IAPWS Formulation 2011 for the Thermal Conductivity  of Ordinary Water Substance
///         <http://www.iapws.org/relguide/ThCond.html>
pub fn thcond(rho: f64, T: f64) -> f64 {
    let d: f64 = rho / DC_WATER;
    let Tr: f64 = T / TC_WATER;
    // Page 6 Table 1.  Coefficients Lk in Eq.(16)
    const no: [f64; 5] = [
        2.443221e-3,
        1.323095e-2,
        6.770357e-3,
        -3.454586e-3,
        4.096266e-4,
    ];
    let mut suma: f64 = 0.0;
    for i in 0..5 {
        suma += no[i] / Tr.powi(i as i32);
    }
    let L0 = Tr.sqrt() / suma;
    // Page 6 Table 2.  Coefficients Lij in Eq(17)
    const nij: [[f64; 6]; 5] = [
        [
            1.60397357,
            -0.646013523,
            0.111443906,
            0.102997357,
            -0.0504123634,
            0.00609859258,
        ],
        [
            2.33771842,
            -2.78843778,
            1.53616167,
            -0.463045512,
            0.0832827019,
            -0.00719201245,
        ],
        [
            2.19650529,
            -4.54580785,
            3.55777244,
            -1.40944978,
            0.275418278,
            -0.0205938816,
        ],
        [
            -1.21051378,
            1.60812989,
            -0.621178141,
            0.0716373224,
            0.0,
            0.0,
        ],
        [
            -2.7203370,
            4.57586331,
            -3.18369245,
            1.1168348,
            -0.19268305,
            0.012913842,
        ],
    ];
    let mut suma: f64 = 0.0;
    for i in 0..5 {
        let mut suma2: f64 = 0.0;
        for j in 0..6 {
            suma2 += nij[i][j] * (d - 1.0).powi(j as i32);
        }
        suma += (1.0 / Tr - 1.0).powi(i as i32) * suma2
    }
    let L1: f64 = (d * suma).exp();
    let L2: f64 = 0.0;
    1e-3 * (L0 * L1 + L2)
}

/// Surface tension  N/m
/// * T :  Temperature K
/// IAPWS, Revised Release on Surface Tension of Ordinary Water Substance June 2014
///   <http://www.iapws.org/relguide/Surf-H2O.html>
pub fn surface_tension(T: f64) -> f64 {
    if 248.15 <= T && T <= TC_WATER {
        let Tr = T / TC_WATER;
        return 1e-3 * (235.8 * (1.0 - Tr).powf(1.256) * (1.0 - 0.625 * (1.0 - Tr)));
    } else {
        return INVALID_VALUE as f64;
    };
}

/// The Static Dielectric Constant of Ordinary Water Substance
///   * rho : Density [kg/m³]
///   * T :  Temperature [K]
/// IAPWS, Release on the Static Dielectric Constant of Ordinary Water
/// Substance for Temperatures from 238 K to 873 K and Pressures up to 1000MPa
//  <http://www.iapws.org/relguide/Dielec.html>
pub fn static_dielectric(rho: f64, T: f64) -> f64 {
    let k: f64 = 1.380658e-23;
    let Na: f64 = 6.0221367e23;
    let alfa: f64 = 1.636e-40;
    let epsilon0: f64 = 8.854187817e-12;
    let mu: f64 = 6.138e-30;
    let M: f64 = 0.018015268;

    let d: f64 = rho / DC_WATER;
    let Tr: f64 = TC_WATER / T;
    const I: [i32; 11] = [1, 1, 1, 2, 3, 3, 4, 5, 6, 7, 10];
    const J: [f64; 11] = [0.25, 1.0, 2.5, 1.5, 1.5, 2.5, 2.0, 2.0, 5.0, 0.5, 10.0];
    const n: [f64; 12] = [
        0.978224486826,
        -0.957771379375,
        0.237511794148,
        0.714692244396,
        -0.298217036956,
        -0.108863472196,
        0.949327488264e-1,
        -0.980469816509e-2,
        0.165167634970e-4,
        0.937359795772e-4,
        -0.12317921872e-9,
        0.196096504426e-2,
    ];

    let mut g: f64 = 1.0 + n[11] * d / (TC_WATER / 228.0 / Tr - 1.0).powf(1.2);
    for i in 0..11 {
        g += n[i] * d.powi(I[i]) * Tr.powf(J[i]);
    }

    let A: f64 = Na * mu * mu * rho * g / M / epsilon0 / k / T;
    let B: f64 = Na * alfa * rho / 3.0 / M / epsilon0;
    let c: f64 = 9.0 + 2.0 * A + 18.0 * B + A * A + 10.0 * A * B + 9.0 * B * B;
    (1.0 + A + 5.0 * B + c.sqrt()) / 4.0 / (1. - B)
}
