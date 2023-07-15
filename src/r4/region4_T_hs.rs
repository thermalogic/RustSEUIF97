
//! Region 4 - Backward Equation Tsat(h,s)
//! 
//！Page25,Page30
//！    http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 9
//！    5.3 Backward Equation Tsat(h,s)
//！          s> 5.210 887 825
//！    Temperature range is T（273.15,623.15）

use crate::algo::fast_ipower::sac_pow;
use crate::common::constant::*;
use crate::common::propertry_id::*;

use crate::r1::region1_pT::*;
use crate::r1::region1::*;
use crate::r2::region2_pT::*;
use crate::r3::region3_v_pT::*;
use crate::r3::region3::*;
use crate::r3::region3_Td::*;
use crate::r4::region4_sat_pT::*;
use crate::r4::region4_pTx::*;
use crate::r4::region4::*;
use crate::common::boundaries::*;


const s4L_623:f64 = 3.778281340;       // T=623.16 ，Sature Liquid s
const s4V_623:f64 = 5.210887825;       // T=623.16 ，Sature steam   s
const h4V_623:f64 = 2.5635920043e+03; // Page 25  Sature steam   h
const s4L_273:f64 = -1.545495919e-04;  // T=273.15 ，Sature Liquid s
const s4V_273:f64= 9.155759395;       // T=273.15 ，Sature steam s

///  Page30  5.3 Backward Functions Tsat(h,s),
pub fn  hs2T_reg43(h:f64, s:f64)->f64
{
    // Table 28 
    const IJn:[(i32,i32,f64);36] = [
        (0,0, 0.179882673606601),
        (0,3, -0.267507455199603),
        (0,12, 0.116276722612600e1),
        (1,0, 0.147545428713616),
        (1,1, -0.512871635973248),

        (1,2, 0.421333567697984),
        (1,5, 0.563749522189870),
        (2,0, 0.429274443819153),
        (2,5, -0.335704552142140e1),
        (2,8, 0.108890916499278e2),

        (3,0, -0.248483390456012),
        (3,2, 0.304153221906390),
        (3,3, -0.494819763939905),
        (3,4, 0.107551674933261e1),
        (4,0, 0.733888415457688e-1),

        (4,1, 0.140170545411085e-1),
        (5,1, -0.106110975998808),
        (5,2, 0.168324361811875e-1),
        (5,4, 0.125028363714877e1),
        (5,16, 0.101316840309509e4),

        (6,6, -0.151791558000712e1),
        (6,8, 0.524277865990866e2),
        (6,22, 0.230495545563912e5),
        (8,1, 0.249459806365456e-1),
        (10,20, 0.210796467412137e+07),

        (10, 36, 0.366836848613065e9),
        (12, 24, -0.144814105365163e9),
        (14, 1, -0.179276373003590e-2),
        (14, 28, 0.489955602100459e10),
        (16, 12, 0.471262212070518e3),

        (16, 32, -0.829294390198652e11),
        (18, 14, -0.171545662263191e4),
        (18, 22, 0.355777682973575e7),
        (18, 36, 0.586062760258436e12),
        (20, 24, -0.129887635078195e8),

        (28, 36, 0.317247449371057e11)];

    if s < s4V_623
        { return INVALID_S as f64;}

    let nu1:f64 = h / 2800.0-0.119;
    let sigma1:f64 = s / 9.2-1.07;
    let mut suma:f64 = 0.0;
    for i in 0..36
     {   suma += IJn[i].2 * sac_pow(nu1, IJn[i].0) * sac_pow(sigma1, IJn[i].1);}
    return 550.0 * suma;
}

pub fn hs2T_reg4(h:f64, s:f64)->f64
{
    if s > s4V_623 && s < s4V_273
    {
        return  hs2T_reg43(h, s);
    };

    //The if97 function hs2Treg43 is only valid for part of region4. Use iteration outsida.
    let mut p_Low_Bound:f64=0.0;
    let mut p_High_Bound:f64=0.0;
    let mut PL:f64=0.0;
    let mut Tsat:f64=0.0;

    if s > s4L_273 && s <= s4L_623
    {
        p_Low_Bound = P_MIN;
        p_High_Bound =Ps_623;

        let mut hL:f64 = -1000.0;
        while (hL - h).abs() > 1.0e-04 && (p_High_Bound - p_Low_Bound).abs() > 1.0e-4
        {
            PL = 0.5*(p_Low_Bound + p_High_Bound); 
            Tsat = T_saturation(PL);
            hL = pT2h_reg1(PL, Tsat);
            if hL > h
               { p_High_Bound = PL;}
            else
               { p_Low_Bound = PL;}
        }
    };

    if s > s4L_623 && s <= SC_WATER
    {
        PL = h2p_sat_reg3(h); // liquid , boundaries
        p_Low_Bound = P_MIN;
        p_High_Bound = PL;
    }
    if s > SC_WATER && s <= s4V_623
    {
        PL = h2p_sat_reg3(h); // steam,boundaries
        p_Low_Bound = P_MIN;
        p_High_Bound = PL;
    }

    let mut sss:f64 = -1000.0;
    let mut p:f64=0.0;
    let mut xs:f64;
    let mut s4v:f64;
    let mut s4L:f64;
    let mut v4v:f64;
    let mut v4L:f64;

    while (s - sss).abs() > 1.0e-6
    {
        p = 0.5 * (p_Low_Bound + p_High_Bound);

        Tsat = T_saturation(p);
        xs = ph_reg4(p, h, OX); // region4_phps 

        if p < Ps_623
        {
            s4v = pT2s_reg2(p, Tsat);
            s4L = pT2s_reg1(p, Tsat);
        }
        else
        {
            v4v = ph_reg3(p, p2sat_steam(p, OH), OV);
            s4v = Td2s_reg3(Tsat, 1.0/ v4v);
            v4L = ph_reg3(p, p2sat_water(p, OH), OV);
            s4L = Td2s_reg3(Tsat, 1.0/v4L);
        };

        sss = xs * s4v + (1.0 - xs) * s4L;

        if sss < s
        {
            p_High_Bound = p;
            p_Low_Bound = (1.0 + (sss - s) / s) * p;
        }
        else
        {
            p_Low_Bound = p;
            p_High_Bound = (1.0 + (sss - s) / s) * p;
        }
    } // end of while  ( abs(s - sss).abd() > 1.0e-6)

    return  T_saturation(p);
}