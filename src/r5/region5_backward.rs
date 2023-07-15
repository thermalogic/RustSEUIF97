//! Region 5 - Backward: (p,h)  (p,s) (h,s)
//! 
//！ Backward of Region 5 using the secant method
//！      (p,h)  (p,s) (h,s)

use crate::algo::fast_ipower::sac_pow;
use crate::algo::root::rtsec1;
use crate::algo::root::rtsec2;
use crate::algo::root::ESP;
use crate::algo::root::I_MAX;
use crate::algo::root::IF97_EQ; 
use crate::common::constant::*;
use crate::r5::region5_pT::*;

const TMAX5:f64 = 2273.15;
const TMIN5:f64 = 1073.15;

const PMAX5:f64 = 50.0;
const PMIN5:f64 = 0.000611212677444;

pub fn ph2T_reg5(p:f64, h:f64)->f64
{
   // TODO:  插值求迭代初值方法
   // double  hmin =pT2hreg2(p,1073.15);
   // double  hmax = pT2hreg5(p,2273.15);
   //  T1=1073.15+(2273.15-1073.15)*(h-hmin)/(hmax-hmin);
   //  求迭代初值方法2： 以三相点为起点,采用与理想气体的对比关系求迭代初值

   let mut T:f64 = -1000.0;
   let mut T1:f64 = 0.5 * (2273.15 + 1073.15);
   let f1:f64 = h - pT2h_reg5(p, T1);
   if f1.abs() > ESP
   {
      let mut T2:f64 =-1000.0;
      if f1 > 0.0
       {  T2 = (1.0 + f1 / h) * T1;}
      else
       {  T2 = (1.0 - f1 / h) * T1;}

      let f2:f64 = h - pT2h_reg5(p, T2);

      T = rtsec2(pT2h_reg5, p, h, T1, T2, f1, f2,ESP, I_MAX);
   }
   else
    {  T = T1;}

   if T < TMIN5
       { T = TMIN5;}
   else  if T > TMAX5
       { T = TMAX5;}

   return T;
}

pub fn ps2T_reg5(p:f64, s:f64)->f64
{
   //TODO: 插值求迭代初值方法1
   // double  smin =pT2sreg2(p,1073.15);
   // double  smax = pT2sreg5(p,2273.15);
   // T1=1073.15+(2273.15-1073.15)*(s-smin)/(smax-smin);

   let mut T:f64 = -1000.0;
   let T1:f64 = 0.5 * (2273.15 + 1073.15); // Get initial value
   let f1:f64 = s - pT2s_reg5(p, T1);
   
   if f1.abs() > ESP
   {
      let mut T2:f64 =-1000.0;
       if f1 > 0.0
      {   T2 = (1.0 + f1 / s) * T1;}
      else
      {   T2 = (1.0 - f1 / s) * T1;}

      let f2:f64 = s - pT2s_reg5(p, T2);
      T = rtsec2(pT2s_reg5, p, s, T1, T2, f1, f2, ESP, I_MAX);
   }
   else
    {  T = T1;}

   if T < TMIN5
      {T = TMIN5;}
   else if T > TMAX5
      { T = TMAX5;}

   return T;
}

/// helper for hs2preg5
fn ph2s_reg5(p:f64, h:f64)->f64
{
   let T:f64 = ph2T_reg5(p, h);
   return pT2s_reg5(p, T);
 }

pub fn hs2p_reg5(h:f64,s:f64)->f64
{
   // TODO:  迭代初始值，可测试那个更好?
   // 也可以计算smin,smax,2元插值得到更接近的p1
  // 测试表明：更复杂的方法计算迭代初始数值并不更好
   //double  hmin =pT2hreg5(PMIN5,1073.15);
   //double  hmax =pT2hreg5(PMAX5,2273.15);
   //p1=PMIN5+(PMAX5-PMIN5)*(h-hmin)/(hmax-hmin);
   
   let mut p:f64=-1000.0;

   let p1:f64 = 0.5 * (PMIN5 + PMAX5);
   let f1:f64 = s - ph2s_reg5(p1, h);
   if f1.abs() > ESP
   {
      let mut p2:f64=0.0; 
      if f1 > 0.0
       {  p2 = (1.0 + f1 / s) * p1;}
      else
       {  p2 = (1.0 - f1 / s) * p1;}

      let f2:f64 = s - ph2s_reg5(p2, h);
      p = rtsec1(ph2s_reg5, h, s, p1, p2, f1, f2, ESP, I_MAX);
   }
   else
    {  p = p1;}
  
   if p < PMIN5
      { p = PMIN5;}
   else if p > PMAX5
      { p = PMAX5;}
   return p;
}
