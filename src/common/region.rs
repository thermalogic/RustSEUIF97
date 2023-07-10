//!
//! Check the Region
//!    pT_sub_region
//!    ph_sub_region
//!    ps_sub_region
//!    hs_sub_region

use crate::common::boundaries::*;
use crate::common::constant::*;
use crate::r1::region1_pT::*;
use crate::r1::region1_T_phps::*;
use crate::r2::region2_pT::*;
use crate::r2::region2_T_ps::*;
use crate::r2::region2_p_hs::*;
use crate::r3::region3_Td::*;
use crate::r3::region3_Tv_phps::*;
use crate::r4::region4_sat_pT::*;
use crate::r5::region5_pT::*;
use crate::r5::region5_backward::*;

/// T in up-order  to check region，
///    p in MPa
///    T in K
/// returns the region
pub fn pT_sub_region(p: f64, T: f64) -> i32
{
    if (p < P_MIN || p > 100.0) {
        return INVALID_P;
    }
    if (T < 273.15 || T > 2273.15) {
        return INVALID_T;
    }
    if (T > 1073.15 && T <= 2273.15 && p > 50.0) {
        return INVALID_P;
    }

    // ON TOP: Saturaton lines、critical point等特殊点判断在前
    // 以后区域判断中没有这些特殊，减少区域判断的复杂度

    //TODO: Saturaton Pressure Tolerance
    const psatTol: f64 = 1.0e-6;
    if (T >= 273.15 && T < TC_WATER) {
        let ps: f64 = p_saturation(T);
        if ((p - ps).abs() / ps < psatTol) {
            return 4;
        }
    }
    // TODO: critical point 放在4区,如果放在3区？
    if (T == TC_WATER && p == PC_WATER) {
        return 4;
    }

    if (T >= 273.15 && T <= 623.15) {
        if (p >= p_saturation(T) && p <= 100.0) {
            return 1;
        }
        if (p < p_saturation(T) && p > P_MIN) {
            return 2;
        }
    };

    // T（623.15,tc_water)之间的饱和线，critical point情况，已在前面处理，这里无需判断，简化了区域判断
    if (T > 623.15 && T <= 863.15) {
        if (p >= P_MIN && p <= B23_T2p(T)) {
            return 2;
        }
        if (p > B23_T2p(T) && p <= 100.0) {
            return 3;
        }
    };

    if (T > 863.15 && T <= 1073.15 && p >= P_MIN && p <= 100.0) {
        return 2;
    }

    if (1073.15 < T && T <= 2273.15 && P_MIN <= p && p <= 50.0) {
        return 5;
    }
    return INVALID_VALUE;
}

/// Pmin -> Ps_623-> Pc-> 100MP ，3 range to check region
///  in each sub region use(hmin, hmax) to chack region
pub fn ph_sub_region(p:f64,h:f64)->i32
{
    let hmin:f64 = pT2h_reg1(p, 273.15);
    let hmax:f64 = pT2h_reg5(p, 2273.15);

    if (P_MIN <= p && p <= Ps_623) // Ps_623 
    {
        let T_sat:f64=T_saturation(p);
        let h14:f64 = pT2h_reg1(p, T_sat);
        let h24:f64 = pT2h_reg2(p, T_sat);
        let h25:f64 = pT2h_reg2(p, 1073.15);

        if (hmin <= h && h <= h14)
            {return 1;}
        else if (h14 < h && h < h24)
            { return 4; }
        else if (h24 <= h && h <= h25)
            { return 2;}
        else if (h25 < h && h <= hmax)
            { return 5;}
    };

    if (Ps_623 < p && p < PC_WATER)
    {

        let h13:f64 = pT2h_reg1(p, 623.15);
        let h32:f64 = pT2h_reg2(p, B23_p2T(p)); //boundaries
        let h25:f64 = pT2h_reg2(p, 1073.15);

        if (hmin <= h && h <= h13)
         {   return 1;}

        if (h13 < h && h < h32)
        {
            let p34:f64 = h2p_sat_reg3(h);//boundaries
            if p < p34
               { return 4;}
            else
               { return 3;}
        };

        if (h32 <= h && h <= h25)
           { return 2;}
        if (h25 < h && h <= hmax)
           { return 5;}
    };

    if (PC_WATER <= p && p <= 100.0)
    {
        let h13:f64 = pT2h_reg1(p, 623.15);
        let h32:f64 = pT2h_reg2(p, B23_p2T(p));
        let h25:f64 = pT2h_reg2(p, 1073.15);

        if (hmin <= h && h <= h13)
           { return 1;}
        if (h13 < h && h < h32)
           {return 3;}
        if (h32 <= h && h <= h25)
        {   return 2;}
        if ((p <= 50.0) && (h25 <= h && h <= hmax))
        {   return 5;}
    };
    return INVALID_VALUE;
}

/// Pmin -> Ps_623-> Pc-> 100MP ，3 range to check region
///  in each sub region use(smin ,smax) to chack region
pub fn ps_sub_region(p:f64,s:f64)->i32
{
    let smin:f64 = pT2s_reg1(p, 273.15);
    let smax:f64 = pT2s_reg5(p, 2273.15);
    
    // 1. First Range: [P_MIN ,Ps_623]
    if (P_MIN <= p && p <= Ps_623)
    {
        let Tsat:f64=T_saturation(p);
        let s14:f64 = pT2s_reg1(p, Tsat);
        let s24:f64 = pT2s_reg2(p,  Tsat);
        let s25:f64 = pT2s_reg2(p, 1073.15);

        if (smin <= s && s <= s14)
           { return 1;}
        if (s14 < s && s < s24)
            {return 4;}
        if (s24 <= s && s <= s25)
           { return 2;}
        if (s25 < s && s <= smax)
           { return 5;}
    };
   
    // 2. Secode Range: (Ps_623,PC_WATER)
    if (Ps_623 < p && p < PC_WATER)
    {
        let s13:f64 = pT2s_reg1(p, 623.15);
        let s32:f64 = pT2s_reg2(p, B23_p2T(p));
        let s25 :f64= pT2s_reg2(p, 1073.15);
        if (smin <= s && s <= s13)
         {   return 1;}
        if (s13 < s && s < s32)
        {
            let p34:f64 = s2p_sat_reg3(s); // boundaries;
            if (p < p34)
              { return 4;}
            else
             {   return 3;}
        };

        if (s32 <= s && s <= s25)
          {  return 2;}
        if (s25 < s && s <= smax)
           { return 5;}
    };
    // 3. Third: [PC_WATER,100.0]
    if (PC_WATER <= p && p <= 100.0)
    {
        let s13:f64 = pT2s_reg1(p, 623.15);
        let s32 :f64= pT2s_reg2(p, B23_p2T(p));
        let s25 :f64= pT2s_reg2(p, 1073.15);
        if (smin <= s && s <= s13)
           { return 1;}
        if (s13 < s && s < s32)
          {  return 3;}
        if (s32 <= s && s <= s25)
         {   return 2;}

        if (p <= 50.0 && s25 <= s && s <= smax)
          {  return 5;}
    };
    return INVALID_VALUE;
}

///  region 1,2,3,4 (smin ->smax), region 5
pub fn hs_sub_region(h:f64, s:f64)->i32
{
    let mut T:f64=0.0;
    let mut p:f64=0.0;
    let mut v:f64=0.0;
    let mut hs:f64=0.0;
    let s13:f64 = pT2s_reg1(100.0, 623.15);
    let s13s:f64 = pT2s_reg1(Ps_623, 623.15);
    let sTPmax:f64 = pT2s_reg2(100.0, 1073.15);
    let s2ab:f64 = pT2s_reg2(4.0, 1073.15); // TODO： p=4 2ab s2ab 的意义？

    // Left point in h-s plot
    let mut smin:f64 = pT2s_reg1(100.0, 273.15);
    let mut hmin:f64 = pT2h_reg1(P_MIN, 273.15);

    // Right point in h-s plot
    let mut hmax:f64 = pT2h_reg2(P_MIN, 1073.15);
    let mut smax:f64 = pT2s_reg2(P_MIN, 1073.15);

    // Region 4 left and right point
    let h4l:f64 = pT2h_reg1(P_MIN, 273.15);
    let s4l:f64 = pT2s_reg1(P_MIN, 273.15);

    let h4v:f64 = pT2h_reg2(P_MIN, 273.15);
    let s4v:f64 = pT2s_reg2(P_MIN, 273.15);

    // !!!! Check region 5 MUST On TOP !!!
    // if （s4v <= s && s<= smax） (h,s)may be setup to error region2
    if (pT2s_reg5(50.0, 1073.15) < s && s <= pT2s_reg5(P_MIN, 2273.15) &&
        pT2h_reg5(50.0, 1073.15) < h && h <= pT2h_reg5(P_MIN, 2273.15))
    {
        p = hs2p_reg5(h, s);
        T = ph2T_reg5(p, h);
        if (1073.15 < T && T <= 2273.15 && P_MIN <= p && p <= 50.0)
          {  return 5; }
    };
    if (smin <= s && s <= s13)
    {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs= hs_region_h1_s(s);
        T = ps2T_reg1(100.0, s) - 0.0218;
        hmax = pT2h_reg1(100.0, T);
        if (hmin <= h && h < hs)
         {  return 4;}
        if (hs <= h && h <= hmax)
          {  return 1;}
    };

 
    if (s13 < s && s <= s13s)
    {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h1_s(s);
        let h13:f64 = hs_region_h13_s(s);
        v = ps2v_reg3(100.0, s) * (1.0 + 9.6e-5);
        T = ps2T_reg3(100.0, s) - 0.0248;
        hmax = Td2h_reg3(T, 1.0/ v);
        if (hmin <= h && h < hs)
          {   return 4;}
        if (hs <= h && h < h13)
           { return 1;}
        if (h13 <= h && h <= hmax)
           { return 3;}
    };

    if (s13s < s && s <= SC_WATER)
    {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h3a_s(s);
        v = ps2v_reg3(100.0, s) * (1.0 + 9.6e-5);
        T = ps2T_reg3(100.0, s) - 0.0248;
        hmax = Td2h_reg3(T, 1.0 / v);
        if (hmin <= h && h < hs)
          {  return 4;}
        if (hs <= h && h <= hmax)
           { return 3;}
    };

    if (SC_WATER < s && s < 5.049096828)
    {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2c3b_s(s);
        v = ps2v_reg3(100.0, s) * (1.0 + 9.6e-5);
        T = ps2T_reg3(100.0, s) - 0.0248;
        hmax = Td2h_reg3(T, 1.0 / v);
        if (hmin <= h && h < hs)
          {  return 4;}
        if (hs <= h && h <= hmax)
          {  return 3;}
    };

    if (5.049096828 <= s && s < 5.260578707)
    {
        // Specific zone with 2-3 boundary in s shape
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2c3b_s(s);
        let h23max:f64 = pT2h_reg2(100.0, 863.15);
        let h23min:f64 = pT2h_reg2(Ps_623, 623.15);
        T = ps2T_reg2(100.0, s) - 0.019;
        hmax = pT2h_reg2(100.0, T);

        if (hmin <= h && h < hs)
           { return 4;}
        if (hs <= h && h < h23min)
            {return 3;}

        if (h23min <= h && h < h23max)
          {  if (hs2p_reg2c(h, s) <= B23_T2p(hs_region_t_hs(h, s))) //hs2p_reg2c r2::region2_p_hs
             {   return 2; }
            else
             {   return 3;}
          }

        if (h23max <= h && h <= hmax)
           { return 2;}
    };

    if (5.260578707 <= s && s < 5.85)
    {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2c3b_s(s);
        T = ps2T_reg2(100.0, s) - 0.019;
        hmax = pT2h_reg2(100.0, T);
        if (hmin <= h && h < hs)
         {   return 4;}
        if (hs <= h && h <= hmax)
          {  return 2;}
    }

    if (5.85 <= s && s < sTPmax)
    {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2ab_s(s);
        T = ps2T_reg2(100.0, s) - 0.019;
        hmax = pT2h_reg2(100.0, T);
        if (hmin <= h && h < hs)
         {   return 4;}
        if (hs <= h && h <= hmax)
          {  return 2;}
    };
 
    if (sTPmax <= s && s < s2ab)
    {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2ab_s(s);
        p = hs2p_reg2(h, s);
        hmax = pT2h_reg2(p, 1073.15);
        if (hmin <= h && h < hs)
          {  return 4;}
        if (hs <= h && h <= hmax)
           { return 2;}
    };

    if (s2ab <= s && s < s4v)
    {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2ab_s(s);
        p = hs2p_reg2(h, s);
        hmax = pT2h_reg2(p, 1073.15);
        if (hmin <= h && h < hs)
          {  return 4;}
        if (hs <= h && h <= hmax)
           { return 2;}
    }

    if (s4v <= s && s <= smax)
    {
        hmin = pT2h_reg2(P_MIN, 273.15);
        p = hs2p_reg2a(h, s);   //hs2p_reg2a r2::region2_p_hs
        hmax = pT2h_reg2(p, 1073.15);
        if (P_MIN <= p && p <= 100.0 && hmin <= h && h <= hmax)
         {   return 2;}
    }
    return INVALID_VALUE;
}
