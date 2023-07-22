//! Propertry ID - o_id
//!
//! 1. The Basic propertry in IAPWS-IF97 Equation - p,t,v,d,h,s,cp,cv,w
//!    * region1_pT.rs, region2_pT.rs, region5_pT.rs
//!    * region3_Td.rs
//! 2. Transport propertry
//!    * common::transport_further.rs
//! 3. The extended propertry
//!    * region1_pT_ext.rs, region2_pT_ext.rs, region5_pT_ext.rs
//!    * region3_Td_ext.rs
//!

/// 0. p - Pressure  MPa   
pub const OP: i32 = 0;
/// 1. t - Temperature  °C
pub const OT: i32 = 1;
/// 2. d - Density   kg/m³
pub const OD: i32 = 2;
/// 3. v - Specific Volume  m³/kg   
pub const OV: i32 = 3;
/// 4. h - Specific enthalpy  kJ/kg  
pub const OH: i32 = 4;
/// 5. s - Specific entropy  kJ/(kg·K) )
pub const OS: i32 = 5;
/// 6. e - Specific exergy  kJ/kg
pub const OE: i32 = 6;
/// 7. u - Specific internal energy  kJ/kg
pub const OU: i32 = 7;
/// 8. cp - Specific isobaric heat capacity  kJ/(kg·K)
pub const OCP: i32 = 8;
/// 9. cv -  Specific isochoric heat capacity  kJ/(kg·K)
pub const OCV: i32 = 9;
/// 10. w -  Speed of sound   m/s
pub const OW: i32 = 10;
/// 11. ks - Isentropic exponent
pub const OKS: i32 = 11;
/// 12. f - Specific Helmholtz free energy   kJ/kg
pub const OF: i32 = 12;
/// 13. g - Specific Gibbs free energy kJ/kg
pub const OG: i32 = 13;
/// 14. z - Compressibility factor
pub const OZ: i32 = 14;
/// 15. x - Steam quality
pub const OX: i32 = 15;
/// 16. r - Region
pub const OR: i32 = 16;
/// 17. ec - Isobaric volume expansion coefficient  1/K
pub const OEC: i32 = 17;
/// 18. kt - Isothermal compressibility   1/MPa  
pub const OKT: i32 = 18;
/// 19. dvdtcp - Partial derivative (∂V/∂T)p  m³/(kg·K)
pub const ODVDT: i32 = 19;
/// 20. dvdpct -  Partial derivative (∂V/∂P)T  m³/(kg·MPa)
pub const ODVDP: i32 = 20;
/// 21. dpdtcv - Partial derivative (∂P/∂T)v   MPa/K
pub const ODPDT: i32 = 21;
/// 22. iJTC -  Isothermal Joule-Thomson coefficient kJ/(kg·MPa)
pub const OIJTC: i32 = 22;
/// 23. joule  - Joule-Thomson coefficient   K/MPa
pub const OJTC: i32 = 23;
/// 24. dv - Dynamic viscosity  Pa.s
pub const ODV: i32 = 24;
/// 25. kv -  Kinematic viscosity  m²/s
pub const OKV: i32 = 25;
/// 26. tc - Thermal conductivity W/(m.K)
pub const OTC: i32 = 26;
/// 27. td - Thermal diffusivity   m²/s
pub const OTD: i32 = 27;
/// 28. pr -  Prandtl number
pub const OPR: i32 = 28;
/// 29. st - Surface tension   N/m
pub const OST: i32 = 29;
/// 30. sdc - Static Dielectric Constant
pub const OSDC: i32 = 30;
/// 31. pc - Isochoric pressure coefficient    1/K
pub const OPC: i32 = 31;
/// 32. betap -  Isothermal stress coefficient, kg/m³
pub const OBETAP: i32 = 32;
/// 33. fi-  Fugacity coefficient 
pub const OFI: i32 = 33;
/// 34. fu - Fugacity MPa
pub const OFU: i32 = 34;
