//！ Region 2 API
//！  IAPWS-IF97 and supp release
//！      1: IF97 IAPWS, R7-97(2012)
//！            IF97-Rev.pdf: P12-32
//！               1)fundamental: (p,t)->v,u,s,h,cp,cv,w
//！               2)backward: (p,h)->T, (p,s)->T
//！      2: IAPWS, SR2-01(2014)
//！               Supp-PHS12-2014.pdf  (h,s)->p

use crate::common::constant::*;
use crate::common::propertry_id::*;

use crate::r2::region2_T_ph::*;
use crate::r2::region2_T_ps::*;
use crate::r2::region2_pT::*;
use crate::r2::region2_p_hs::*;

pub fn pT_reg2(p:f64,T:f64, o_id:i32)->f64
{
    match o_id
    {
      OV=> pT2v_reg2(p, T),
      OD=> 1.0/pT2v_reg2(p, T),
      OH=> pT2h_reg2(p, T),
      OS=> pT2s_reg2(p, T),
      OU=> pT2u_reg2(p, T),
      OCV=> pT2cv_reg2(p, T),
      OCP=> pT2cp_reg2(p, T),
      OW=> pT2w_reg2(p, T),
      _=> INVALID_OUTID as f64,
    }
    
}

pub fn pt_reg2(p:f64,t:f64,o_id:i32) -> f64 {
    let T=t+273.15;
    pT_reg2(p,T,o_id)
}

pub fn ph_reg2(p:f64,h:f64, o_id:i32)->f64
{
   let T:f64=ph2T_reg2(p,h);
   if o_id==OT 
   {   return T-273.15;   };
   return  pT_reg2(p, T, o_id); 
}

pub fn ps_reg2(p:f64,s:f64, o_id:i32)->f64
{
    let T:f64=ps2T_reg2(p,s);
    if o_id==OT 
    {   return T-273.15;   };
    return  pT_reg2(p, T, o_id);
}

pub fn hs_reg2(h:f64,s:f64, o_id:i32)->f64
{
   let p:f64=hs2p_reg2(h,s);
   if o_id==OP
   {
       return p;
   }
   return ph_reg2(p,h, o_id);  
}

