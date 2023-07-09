/*---------------------------------------------------------------------------
   IAPWS-IF97 Backware Equation (H,S)->P for Region 1:
     
     IAPWS-IF97-S01: Supp-PHS12-2014.pdf  June 2014

  Author: Cheng Maohua 
--------------------------------------------------------------------- */
use crate::common::constant::*;
use crate::algo::fast_ipower::sac_pow;
use crate::algo::root::rtsec1;
use crate::algo::root::ESP;
use crate::algo::root::I_MAX;
use crate::algo::root::IF97_EQ; 
use crate::r1::region1_pT::*;
use crate::r1::region1_T_phps::*;

// helper for iter (h,s)->p
fn ph2s_reg1(p:f64, h:f64)->f64
{
  let T:f64=ph2T_reg1(p, h);
  return pT2s_reg1(p, T);
}


//----------------------------------------------------------------
//  Backward equation p(h,s) for region 1
//----------------------------------------------------------------
pub fn hs2p_reg1(h:f64, s:f64)->f64
{
  // Page 5, Table 2 :
  // Initialize coefficients and exponents (H,S)->P for region 1
  const IJn: [IJnData; 19] = [
    IJnData{I:0, J:0, n:-0.691997014660582},
      IJnData{I:0, J:1, n:-0.183612548787560e2},
      IJnData{I:0, J:2, n:-0.928332409297335e1},
      IJnData{I:0, J:4, n:0.659639569909906e2},
      IJnData{I:0, J:5, n:-0.162060388912024e2},

      IJnData{I:0, J:6, n:0.450620017338667e3},
      IJnData{I:0, J:8, n:0.854680678224170e3},
      IJnData{I:0, J:14, n:0.607523214001162e4},
      IJnData{I:1, J:0, n:0.326487682621856e2},
      IJnData{I:1, J:1, n:-0.269408844582931e2},

      IJnData{I:1, J:4, n:-0.319947848334300e3},
      IJnData{I:1, J:6, n:-0.928354307043320e3},
      IJnData{I:2, J:0, n:0.303634537455249e2},
      IJnData{I:2, J:1, n:-0.650540422444146e2},
      IJnData{I:2, J:10, n:-0.430991316516130e4},

      IJnData{I:3, J:4, n:-0.747512324096068e3},
      IJnData{I:4, J:1, n:0.730000345529245e3},
      IJnData{I:4, J:4, n:0.114284032569021e4},
      IJnData{I:5, J:0, n:-0.436407041874559e3}];

  let eta:f64 = h / 3400.0 + 0.05;
  let sigma:f64 = s / 7.6 + 0.05;
  let mut pi = 0.0;
  for k in 0..19
  {
      pi += IJn[k].n * sac_pow(eta,IJn[k].I) *sac_pow(sigma,IJn[k].J);
      //pi += IJn[k].n * eta.powi(IJn[k].I) *sigma.powi(IJn[k].J);
  }
  // return 100.0*pi;

  // iteration: refine
  
  let mut p:f64;
  let p1 = 100.0 * pi;
  let f1 = s - ph2s_reg1(p1, h);
  let mut p2:f64;
  if f1.abs() > ESP
  {
    if f1 > 0.0 // pT2sreg1(p,h)< s ,the p1< expt p，so， p2=1.05*p1 p（p1,p2)
      { p2 = (1.0 + f1 / s) * p1;}
    else
      { p2 = (1.0 - f1 / s) * p1;}

    let f2:f64 = s - ph2s_reg1(p2, h);
    p = rtsec1(ph2s_reg1, h, s, p1, p2, f1, f2,ESP, I_MAX);
  }
  else
  {  p = p1;};
  return p; 
}
