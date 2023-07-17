//! Propertry ID - o_id

// 1 properties in IAPWS-IF97 Equation

/// 0. p Pressure  MPa      
pub const OP:i32=0;     
/// 1. t - Temperature  °C  
pub const OT:i32=1;     
/// 2. d - Density   kg/m^3     
pub const OD:i32=2;     
/// 3. v - Specific Volume  m^3/kg   
pub const OV:i32=3;  
/// 4. h - Specific enthalpy  kJ/kg    
pub const OH:i32=4;  
/// 5. s - Specific entropy  kJ/(kg·K) 
pub const OS:i32=5;   
/// 6. u - Specific internal energy  kJ/kg 
pub const OU:i32=6;   
/// 7. cp - Specific isobaric heat capacity  kJ/(kg·K)
pub const OCP:i32=7;  
/// 8. cv - Specific isochoric heat capacity  kJ/(kg·K) 
pub const OCV:i32=8; 
/// 9. w - Speed of sound  m/s 
pub const OW:i32=9; 
/// 10. x - Steam quality   
pub const OX:i32=10;  
/// 11. r - Region in IF97   
pub const OR:i32=11;   

// 2 transport and further propertries
/// 12. dv - Dynamic viscosity  kg/(m·s) or Pa.s
pub const ODV:i32=12;   
/// 13. kv - Kinematic viscosity  m^2/s 
pub const OKV:i32=13;  
/// 14. tc - Thermal conductivity W/(m.K) 
pub const OTC:i32=14;   
/// 15. td - Thermal diffusivity  m^2/s 
pub const OTD:i32=15;   
/// 16. pr - Prandtl number 
pub const OPR:i32=16;   
/// 17. st - Surface tension   N/m 
pub const OST:i32=17;    
/// 18. sdc - Static Dielectric constant 
pub const OSDC:i32=18;  