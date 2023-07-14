
/// struct for coefficients and exponents
pub struct IJnData {
   pub I: i32,
   pub J: i32,
   pub n: f64,
}



/// code for error input
pub const INVALID_VALUE:i32=-9999;
pub const INVALID_OUTID:i32=-1000;
pub const INVALID_P:i32=-2100;
pub const INVALID_T:i32=-2101;
pub const INVALID_S:i32=-2102;
pub const INVALID_H:i32=-2103;
pub const INVALID_PT:i32=-2201;
pub const INVALID_HS:i32=-2202;

/// constants  
pub const K:f64=273.15;
pub const RGAS_WATER:f64 = 0.461526; // gas constant in KJ/(kg K)
/// critical point
pub const TC_WATER:f64 = 647.096;          //critical temperature in K
pub const PC_WATER:f64 = 22.064;           //critical p in Mpa
pub const DC_WATER:f64= 322.0;             //critical density in kg/m**3
pub const SC_WATER:f64 = 4.41202148223476; // Critic entropy
pub const HC_WATER:f64= 2.087546845e+03;   // Critic entropy h
/// triple point
pub const Pt:f64 = 611.657e-6;     // the triple point
pub const Tt:f64 = 273.16;         // the triple point
pub const st_water:f64 = 5.85;     // the triple point
pub const ht_water:f64 = 0.611783; // the triple point

// T=623.15 region (1,3)
pub const Ps_623:f64 = 16.5291642526045; // P_MIN3 Ps_623 = _PSat_T(623.15)  P Saturation at 623.15 K, boundary region 1-3
// T=273.15 Tmin
pub const P_MIN:f64 = 0.000611212677444; // P_MIN = _PSat_T(273.15)  Minimum pressure