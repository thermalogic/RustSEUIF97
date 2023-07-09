/*---------------------------------------------------------------------------
   IAPWS-IF 97 Backware Equation for Region 1:
      August 2007 IF97: IF97-Rev.pdf: P6-9
          (p,h)->T, (p,s)->T
    
  Author: Maohua Cheng
---------------------------------------------------------------------------*/

use crate::common::constant::*;
use crate::algo::fast_ipower::sac_pow;
use crate::algo::root::rtsec2;
use crate::algo::root::ESP;
use crate::algo::root::I_MAX;
use crate::algo::root::IF97_EQ; 
use crate::r1::region1_pT::*;

//-------------------------------------------------------------
// Backward equation T(p,h) for region 1
//--------------------------------------------------------------
pub fn ph2T_reg1(p:f64, h:f64)->f64
{
  // Page 11, Table6 :Initialize coefficients and exponents (P,H)->T for region 1
  const IJn: [IJnData; 20] = [
    IJnData{I:0, J:0,n:-238.72489924521},
      IJnData{I:0, J:1, n:404.21188637945},
      IJnData{I:0, J:2, n:113.49746881718}, 
      IJnData{I:0, J:6, n:-5.8457616048039}, 
      IJnData{I:0, J:22, n:-1.528548241314E-04},

      IJnData{I:0, J:32, n:-1.0866707695377E-06},
      IJnData{I:1, J:0, n:-13.391744872602},
      IJnData{I:1, J:1, n:43.211039183559},
      IJnData{I:1, J:2, n:-54.010067170506},
      IJnData{I:1, J:3, n:30.535892203916},

      IJnData{I:1, J:4, n:-6.5964749423638},
      IJnData{I:1, J:10, n:9.3965400878363E-03},
      IJnData{I:1, J:32, n:1.157364750534E-07},
      IJnData{I:2, J:10, n:-2.5858641282073E-05},
      IJnData{I:2, J:32, n:-4.0644363084799E-09},

      IJnData{I:3, J:10, n:6.6456186191635E-08},
      IJnData{I:3, J:32, n:8.0670734103027E-11},
      IJnData{I:4, J:32, n:-9.3477771213947E-13},
      IJnData{I:5, J:32, n:5.8265442020601E-15},
      IJnData{I:6, J:32, n:-1.5020185953503E-17}];

  let pi:f64 = p / 1.0;
  let eta:f64 = h / 2500.0 + 1.0;
  let mut theta=0.0;
  for k in 0..20
  {
    // theta += IJn[k].n * pi.powi(IJn[k].I) * eta.powi(IJn[k].J);
      theta += IJn[k].n * sac_pow(pi,IJn[k].I) * sac_pow(eta,IJn[k].J);
  }
// return 1.0 *theta;

  // iteration: refine
  let T1:f64 = 1.0 * theta;
  let f1:f64 = h - pT2h_reg1(p, T1);
  let mut T2:f64;
  let mut T:f64;
  if f1.abs() > ESP
  {
    if f1 > 0.0
      {T2 = (1.0 + f1 / h) * T1; }// TODO： 1.05 用 1+f1/h 是不是更快，没有测试
    else
      {T2 = (1.0 - f1 / h) * T1;}

    let f2:f64 = h - pT2h_reg1(p, T2);
    T = rtsec2(pT2h_reg1, p, h, T1, T2, f1, f2, ESP, I_MAX);
  }
  else
  { T = T1;}; 

  return T; 

}

//----------------------------------------------------------------
//  Backward equation T(p,s) for region 1
//----------------------------------------------------------------
pub fn ps2T_reg1(p:f64, s:f64)->f64
// Page 12, Table 8 : Initialize coefficients and exponents (P,S)->T for region 1
{
  const IJn: [IJnData; 20] = [
      IJnData{I:0, J:0, n:0.17478268058307e3},
      IJnData{I:0, J:1, n:0.34806930892873e2},
      IJnData{I:0, J:2, n:0.65292584978455e1},
      IJnData{I:0, J:3, n:0.33039981775489},
      IJnData{I:0, J:11, n:-0.19281382923196e-6},

      IJnData{I:0, J:31, n:-0.24909197244573e-22},
      IJnData{I:1, J:0, n:-0.26107636489332},
      IJnData{I:1, J:1, n:0.22592965981586},
      IJnData{I:1, J:2, n:-0.64256463395226e-1},
      IJnData{I:1, J:3, n:0.78876289270526e-2},

      IJnData{I:1, J:12, n:0.35672110607366e-9},
      IJnData{I:1, J:31, n:0.17332496994895e-23},
      IJnData{I:2, J:0, n:0.56608900654837e-3},
      IJnData{I:2, J:1, n:-0.32635483139717e-3},
      IJnData{I:2, J:2, n:0.44778286690632e-4},

      IJnData{I:2, J:9, n:-0.51322156908507e-9},
      IJnData{I:2, J:31, n:-0.42522657042207e-25},
      IJnData{I:3, J:10, n:0.26400441360689e-12},
      IJnData{I:3, J:32, n:0.78124600459723e-28},
      IJnData{I:4, J:32, n:-0.30732199903668e-30}];

  let pi:f64 = p / 1.0;
  let sigma:f64 = s / 1.0 + 2.0;

  let mut theta = 0.0;
  for k in 0..20
  { 
    // theta += IJn[k].n * pi.powi(IJn[k].I) *sigma.powi(IJn[k].J);
     theta += IJn[k].n * sac_pow(pi,IJn[k].I) *sac_pow(sigma,IJn[k].J);
 
  }
  //return 1.0*theta;

 // iteration: refine
  let T1:f64 = 1.0 * theta;
  let f1:f64 = s - pT2s_reg1(p, T1);
  let T2:f64;
  let mut T:f64;
  if f1.abs() >ESP
  {
    if f1 > 0.0 // pT2s_reg1(p,T1)< s ,the T1< expt T，so， T2=1.05*T1 T（T1,T2)
      {T2 = (1.0 + f1 / s) * T1;}
    else
      {T2 = (1.0 - f1 / s) * T1;}

    let f2:f64 = s - pT2s_reg1(p, T2);
    T = rtsec2(pT2s_reg1, p, s, T1, T2, f1, f2, ESP, I_MAX);
  }
  else
   { T = T1;}
  return T;
}
