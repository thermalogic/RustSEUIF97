//! Constant

/// error code for  input
pub const INVALID_VALUE: i32 = -9999;
pub const INVALID_OUTID: i32 = -1000;
pub const INVALID_P: i32 = -2100;
pub const INVALID_T: i32 = -2101;
pub const INVALID_S: i32 = -2102;
pub const INVALID_H: i32 = -2103;
pub const INVALID_PT: i32 = -2201;
pub const INVALID_HS: i32 = -2202;

///  FLOAT_ERROR for float value with `==`
pub const FLOAT_ERROR: f64 = 1.0e-6;

/// constants  
pub const K: f64 = 273.15;
pub const RGAS_WATER: f64 = 0.461526; //gas constant in KJ/(kg K)

/// critical point
pub const TC_WATER: f64 = 647.096; //critical temperature in K
pub const PC_WATER: f64 = 22.064; //critical p in Mpa
pub const DC_WATER: f64 = 322.0; //critical density in kg/m**3
pub const SC_WATER: f64 = 4.41202148223476; // Critic entropy
pub const HC_WATER: f64 = 2.087546845e+03; // Critic entropy h

/// the triple point of water
pub const Pt: f64 = 611.657e-6; // the triple point
pub const Tt: f64 = 273.16; // the triple point
pub const st_water: f64 = 5.85; // the triple point
pub const ht_water: f64 = 0.611783; // the triple point

/// boundary constants

// T=623.15 region (1,3)
pub const Ps_623: f64 = 16.5291642526045; // P_MIN3 Ps_623 = _PSat_T(623.15)  P Saturation at 623.15 K, boundary region 1-3

pub const P_MIN: f64 = 0.000611212677444; // P_MIN = _PSat_T(273.15)  Mininum pressure
pub const V_MAX: f64 = 1.0E+10;
pub const V_MIN: f64 = 0.00095;
pub const H_MAX: f64 = 7376.99;
pub const H_MIN: f64 = -0.1;
pub const S_MAX: f64 = 18.992; // 1.0E-8, 2273.15
pub const S_MIN: f64 = -0.008583; // 100,273.15

/// Region1  boundary constants
pub const P01: f64 = 16.53;
pub const T01: f64 = 1386.0;
pub const T_MAX1: f64 = 623.15;
pub const T_MIN1: f64 = 273.15;
pub const P_MAX1: f64 = 100.00;
pub const P_MIN1: f64 = 0.000611212677444; // T=273.15

/// Region2  boundary constants
pub const T02: f64 = 540.0;
pub const P02: f64 = 1.0;
pub const T_MAX2: f64 = 1073.15;
pub const T_MIN2: f64 = 273.15;
pub const P_MAX2: f64 = 100.00;
pub const P_MIN2: f64 = 1.0E-8;

/// Region2  boundary constants
pub const D03: f64 = 322.0;
pub const T03: f64 = 647.096;
pub const T_MAX3: f64 = 863.15;
pub const T_MIN3: f64 = 623.15;
pub const P_MIN3: f64 = 16.5291643;
pub const P_MAX3: f64 = 100.0;

/// Region4  boundary constants
pub const T_MAX4: f64 = 647.096;
pub const T_MIN4: f64 = 273.15;
pub const P_MAX4: f64 = 22.064;
pub const P_MIN4: f64 = 0.000611212677444;
pub const H_MAX4: f64 = 2803.2738880203456;
pub const H_MIN4: f64 = 0.38039218566765925;
pub const S_MAX4: f64 = 9.153081306725117;
pub const S_MIN4: f64 = 0.33567511183921145;

/// Region5  boundary constants
pub const T05: f64 = 1000.0;
pub const P05: f64 = 1.0;
pub const T_MAX5: f64 = 2273.15;
pub const T_MIN5: f64 = 1073.15;
pub const P_MAX5: f64 = 50.0;
pub const P_MIN5: f64 = 0.000611212677444; 
