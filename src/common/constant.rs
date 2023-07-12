
/// struct for coefficients and exponents
pub struct IJnData {
   pub I: i32,
   pub J: i32,
   pub n: f64,
}


/// constants of propertry ID - o_id

// 1 propertry in IF97 Basic Equation
pub const OP:i32=0;   //  p  Pressure   MPa     
pub const OT:i32=1;   //  t  Temperature  °C        
pub const OD:i32=2;   //  d  Density   kg/m^3    
pub const OV:i32=3;   //  v  Specific Volume    m^3/kg   
pub const OH:i32=4;   //  h  Specific enthalpy  kJ/kg 
pub const OS:i32=5;   //  s  Specific entropy   kJ/(kg·K) 
pub const OU:i32=6;   //  u  Specific internal energy  kJ/kg 
pub const OCP:i32=7;  // cp  Specific isobaric heat capacity  kJ/(kg·K) 
pub const OCV:i32=8;  // cv  Specific isochoric heat capacity  kJ/(kg·K) 
pub const OW:i32=9;  //   w  Speed of sound   m/s 

pub const OX:i32=10;   // x  Steam quality   
pub const OR:i32=11;   // r  Region in IF97   

// 2 transport and further propertry 
pub const ODV:i32=12;   // dv Dynamic viscosity  kg/(m·s) or mu(Pa.s)
pub const OKV:i32=13;   // kv Kinematic viscosity  m^2/s 
pub const OTC:i32=14;   // tc Thermal conductivity W/(m.K) 
pub const OTD:i32=15;   // td Thermal diffusivity   um^2/s 
pub const OPR:i32=16;   // pr Prandtl number  
pub const OST:i32=17;   // st Surface tension   mN/m  
pub const OSDC:i32=18;  // sdc Static Dielectric Constant 

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