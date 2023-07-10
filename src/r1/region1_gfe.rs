//！  IAPWS-IF97 Basic Equation for Region 1
//！
//！   The dimensionless Gibbs free energy (p,T) EQUATIONS
//！          gamma and its derivatives
//！   Release : August 2007  Equations for Region 1 P6-9
//！             according to Eq. (7), P8 release if97-rev 2007

use crate::algo::fast_ipower::sac_pow;
use crate::common::constant::*;

//  Initialize coefficients and exponents for region 1
pub const IJn: [IJnData; 34] = [
    IJnData{I:0, J:-2,n: 0.14632971213167E+00},
    IJnData{I:0, J:-1,n: -0.84548187169114E+00},
    IJnData{I:0, J:0, n:-0.37563603672040E+01},
    IJnData{I:0, J:1, n:0.33855169168385E+01},
    IJnData{I:0, J:2, n:-0.95791963387872E+00},
    //
    IJnData{I:0, J:3, n:0.15772038513228E+00},
    IJnData{I:0, J:4, n:-0.16616417199501E-01},
    IJnData{I:0, J:5, n:0.81214629983568E-03},
    IJnData{I:1, J:-9, n:0.28319080123804E-03},
    IJnData{I:1, J:-7, n:-0.60706301565874E-03},
    //
    IJnData{I:1, J:-1, n:-0.18990068218419E-01},
    IJnData{I:1, J:0, n:-0.32529748770505E-01},
    IJnData{I:1, J:1, n:-0.21841717175414E-01},
    IJnData{I:1, J:3, n:-0.52838357969930E-04},
    IJnData{I:2, J:-3, n:-0.47184321073267E-03},
    //
    IJnData{I:2, J:0, n:-0.30001780793026E-03},
    IJnData{I:2, J:1, n:0.47661393906987E-04},
    IJnData{I:2, J:3, n:-0.44141845330846E-05},
    IJnData{I:2, J:17, n:-0.72694996297594E-15},
    IJnData{I:3, J:-4, n:-0.31679644845054E-04},
    //
    IJnData{I:3, J:0, n:-0.28270797985312E-05},
    IJnData{I:3, J:6, n:-0.85205128120103E-09},
    IJnData{I:4, J:-5, n:-0.22425281908000E-05},
    IJnData{I:4, J:-2, n:-0.65171222895601E-06},
    IJnData{I:4, J:10, n:-0.14341729937924E-12},
    //
    IJnData{I:5, J:-8, n:-0.40516996860117E-06},
    IJnData{I:8, J:-11,n: -0.12734301741641E-08},
    IJnData{I:8, J:-6, n:-0.17424871230634E-09},
    IJnData{I:21, J:-29,n: -0.68762131295531E-18},
    IJnData{I:23, J:-31,n: 0.14478307828521E-19},
    //
    IJnData{I:29, J:-38, n:0.26335781662795E-22},
    IJnData{I:30, J:-39, n:-0.11947622640071E-22},
    IJnData{I:31, J:-40, n:0.18228094581404E-23},
    IJnData{I:32, J:-41, n:-0.93537087292458E-25},
];

/// Fundamental equation for region 1
pub fn gamma_reg1(mut tau: f64, mut pi: f64) -> f64
{
    tau = tau - 1.222;
    pi = 7.1 - pi;
    let mut value: f64 = 0.0;
    for k in 0..34 {
        value += IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(tau, IJn[k].J);
        // value += IJn[k].n * pi.powi( IJn[k].I) * tau.powi(IJn[k].J);
    }
    return value;
}

/// First derivative of fundamental equation in pi for region 1
pub fn gamma_pi_reg1(mut tau: f64, mut pi: f64) -> f64
{
    tau = tau - 1.222;
    pi = 7.1 - pi;
    let mut value: f64 = 0.0;
    for k in 0..34 {
        // note: pi = 7.1 - pi， so  - minus
         value -= IJn[k].n * IJn[k].I  as f64 * sac_pow(pi, IJn[k].I - 1) * sac_pow(tau, IJn[k].J);
        // value -= IJn[k].n  as f64 * IJn[k].I  as f64 * pi.powi(IJn[k].I - 1) * tau.powi(IJn[k].J);
    }
    return value;
}

/// Second derivative of fundamental equation in pi for region 1
pub fn gamma_pipi_reg1(mut tau: f64, mut pi: f64) -> f64
{
    tau = tau - 1.222;
    pi = 7.1 - pi;
    let mut value: f64 = 0.0;
    for k in 0..34 {
         value +=
             IJn[k].n  * IJn[k].I as f64  * (IJn[k].I - 1)  as f64 * sac_pow(pi, IJn[k].I - 2)  * sac_pow(tau, IJn[k].J);
       // value +=
        //    IJn[k].n * IJn[k].I as f64  * (IJn[k].I - 1)  as f64 * pi.powi(IJn[k].I - 2)  * tau.powi(IJn[k].J);     
    }
    return value;
}

/// First derivative of fundamental equation in tau for region 1
pub fn gamma_tau_reg1(mut tau: f64, mut pi: f64) -> f64
{
    tau = tau - 1.222;
    pi = 7.1 - pi;
    let mut value: f64 = 0.0;
    for k in 0..34 {
        value += IJn[k].n * sac_pow(pi, IJn[k].I) * IJn[k].J  as f64 * sac_pow(tau, IJn[k].J - 1); // OK!
       // value += IJn[k].n * pi.powi(IJn[k].I) * IJn[k].J as f64 * tau.powi(IJn[k].J - 1);
    }
    return value;
}

/// Second derivative of fundamental equation in tau for region 1
pub fn gamma_tautau_reg1(mut tau: f64, mut pi: f64) -> f64
{
    tau = tau - 1.222;
    pi = 7.1 - pi;
    let mut value = 0.0;
     for k in 0..34 {
          value +=
               IJn[k].n  * sac_pow(pi, IJn[k].I) * IJn[k].J  as f64 * (IJn[k].J - 1)  as f64 * sac_pow(tau, IJn[k].J - 2);
          //value +=
            //  IJn[k].n * pi.powi(IJn[k].I) * IJn[k].J  as f64 * (IJn[k].J - 1)  as f64 * tau.powi(IJn[k].J - 2);     
    }
    return value;
}

/// Second derivative of fundamental equation in pi and tau for region 1
pub fn gamma_pitau_reg1(mut tau: f64, mut pi: f64) -> f64
{
    tau = tau - 1.222;
    pi = 7.1 - pi;
    let mut value: f64 = 0.0;
     for k in 0..34 {
        // note: pi = 7.1 - pi，so - minus
        value -= IJn[k].n  * IJn[k].I  as f64 * sac_pow(pi, IJn[k].I - 1) * IJn[k].J  as f64 * sac_pow(tau, IJn[k].J - 1);
        //value -= IJn[k].n  * IJn[k].I  as f64 * pi.powi(IJn[k].I - 1) * IJn[k].J  as f64 * tau.powi(IJn[k].J - 1);
   
     }
    return value;
}
