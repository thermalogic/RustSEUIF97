
/*------------------------------------------------------------------

    5 Transport Properties
    
    * ok Dynamic viscosity  Pa.s  dv (mu)   
    * Kinematic viscosity  m^2/s    kv     25 
    * ok Thermal conductivity  W/(m.K)     tc  k  26    
|   * Thermal diffusivity   um^2/s   td        27 
|   * Prandtl number                  pr         28 
|   Further properties-
      * Refractive index -待
      * Relative static dialectric sonstant  -待
      * ok Surface tension     st   mN/m      29 sigma  [N/m]
  
 ---------------------------------------------------------------------------*/
 use crate::algo::fast_ipower::sac_pow;
 use crate::common::constant::*;
 

 /*----------------------------------------------------------------
 
  Transport properties
----------------------------------------------------------------*/
 pub fn viscosity(rho:f64, T:f64)->f64
 {
 /* The Viscosity for IF97
     Parameters
       rho :  Density [kg/m³]
       T :    Temperature [K]
       Returns：  mu : float  Viscosity [Pa·s]
     References
     ----------
     IAPWS, Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary  Water Substance
         http://www.iapws.org/relguide/viscosity.html
*/
     let Tr:f64 = T/TC_WATER;
     let Dr:f64 = rho/DC_WATER;

     let no:[f64;4] = [1.67752, 2.20462, 0.6366564, -0.241605];
     let mut suma:f64 = 0.0;
     for i in 0..4
       {  suma += no[i]/sac_pow(Tr,i as i32)}
     let fi0:f64 = 100.0*Tr.sqrt()/suma;
 
     const I:[i32;21] = [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6];
     const J:[i32;21] = [0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5];
     const nr:[f64;21] = [0.520094, 0.850895e-1, -0.108374e1, -0.289555, 0.222531, 0.999115,
           0.188797e1, 0.126613e1, 0.120573, -0.281378, -0.906851, -0.772479,
           -0.489837, -0.257040, 0.161913, 0.257399, -0.325372e-1, 0.698452e-1,
           0.872102e-2, -0.435673e-2, -0.593264e-3];
     suma = 0.0;
     for i in 0..21
      { 
          suma += nr[i]*sac_pow((Dr-1.0),I[i])*sac_pow((1.0/Tr-1.0),J[i]);
      }
     let fi1:f64 = (Dr*suma).exp();
     let fi2:f64= 1.0;
     return fi0*fi1*fi2*1.0e-6;
    }


pub fn thcond(rho:f64, T:f64)->f64
 /* Equation for the thermal conductivity
        Parameters
          rho :     Density [kg/m³]
          T :    Temperature [K]
       Returns 
           k : Thermal conductivity [W/mK]
    References 
    ----------
    IAPWS, Release on the IAPWS Formulation 2011 for the Thermal Conductivity  of Ordinary Water Substance
         http://www.iapws.org/relguide/ThCond.html
*/
{
    let d:f64 = rho/DC_WATER;
    let Tr:f64 = T/TC_WATER;
    // Page 6 Table 1.  Coefficients Lk in Eq.(16) 
    const no:[f64;5] = [2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4];
    let mut suma:f64 = 0.0;
    for i in 0..5
    {    suma += no[i]/sac_pow(Tr,i as i32);}
    let L0 = Tr.sqrt()/suma;
    // Page 6 Table 2.  Coefficients Lij in Eq. (17) ρ
    const nij:[[f64;6];5] = [
        [1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634,
         0.00609859258],
        [2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019,
         -0.00719201245],
        [2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278,
         -0.0205938816],
        [-1.21051378, 1.60812989, -0.621178141, 0.0716373224, 0.0, 0.0],
        [-2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305,
         0.012913842]];
    let mut suma:f64 = 0.0;
    for i in 0..5
    {   let mut suma2:f64 = 0.0;
        for j in 0..6
        {    suma2 += nij[i][j]*sac_pow(d-1.0,j as i32);}
        suma += sac_pow(1.0/Tr-1.0,i as i32)*suma2
     }
    let L1:f64 =(d*suma).exp();

    let L2:f64 = 0.0;
    return 1e-3*(L0*L1+L2);
 }

 pub fn surface_tension(T:f64)->f64
 {
/* Equation for the surface tension
  Parameters
    ----------
    T :  Temperature [K]
    Returns
    -------
    sigma :   Surface tension [N/m]
   
    References
    ----------
    IAPWS, Revised Release on Surface Tension of Ordinary Water Substance June 2014
           http://www.iapws.org/relguide/Surf-H2O.html
    """
*/
    if 248.15 <= T && T<= TC_WATER
    {    let Tr = T/TC_WATER;
        return 1e-3*(235.8*(1.0-Tr).powf(1.256)*(1.0-0.625*(1.0-Tr)));
    }
    else
       { return INVALID_VALUE as f64};
  }