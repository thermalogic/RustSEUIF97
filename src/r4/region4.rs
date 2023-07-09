
/*---------------------------------------------------------------------------
  Region 4 API: 
    (p,x) (t,x)  (p,h) (p,s) (h,s)
      
 Author: Maohua Cheng
---------------------------------------------------------------------------*/
use crate::common::constant::*;
use crate::r4::region4_sat_pT::*;
use crate::r4::region4_pTx::*;
use crate::r4::region4_T_hs::*;

pub fn p_sat(t:f64)->f64
{
   p_saturation(t+K)
}

pub fn t_sat(p:f64)->f64
{
   T_saturation(p)-K
}

pub fn ph_reg4(p:f64,h:f64, o_id:i32)->f64
{
   let h1:f64 = p2sat_water(p, OH);
   let h2:f64 = p2sat_steam(p, OH);
   let x:f64 = (h - h1) / (h2 - h1);
   if o_id == OX
   {  return x;}
   return px(p, x, o_id);
}

pub fn  ps_reg4(p:f64,s:f64, o_id:i32)->f64
{
   let s1:f64 = p2sat_water(p, OS);
   let s2:f64 = p2sat_steam(p, OS);
   let x:f64 = (s - s1) / (s2 - s1);
   if o_id == OX
   {   return x;}
   return px(p, x, o_id);
}

pub fn hs_reg4(h:f64,s:f64, o_id:i32)->f64
{

   // python版本（h,s)  只能计算T<623.15,
   let T:f64 = hs2T_reg4(h, s);
   if o_id == OT
    {  return T-K;}
   let p:f64 = p_saturation(T);
   if o_id == OP
   {   return p;}

   let h1:f64 = p2sat_water(p, OH);
   let h2:f64 = p2sat_steam(p, OH);

   let x:f64 = (h - h1) / (h2 - h1);
   if o_id == OX
   {   return x;}

   return px(p, x, o_id);
}