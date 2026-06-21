//! Check the Region
//! （h,s)

use crate::common::boundaries::*;
use crate::common::constant::*;
use crate::common::property_id::*;

use crate::r1::region1_T_phps::*;
use crate::r1::region1_pT::*;

use crate::r2::region2_T_ps::*;
use crate::r2::region2_pT::*;
use crate::r2::region2_p_hs::*;

use crate::r3::region3::*;
use crate::r3::region3_Td::*;
use crate::r3::region3_Tv_phps::*;
use crate::r3::region3_v_pT::*;

use crate::r4::region4_pTx::*;
use crate::r4::region4_sat_pT::*;

use crate::r5::region5_pT::*;
use crate::r5::region5_ph_ps_hs::*;

///  region 1,2,3,4 (smin ->smax), region 5
pub fn hs_sub_region(h: f64, s: f64) -> i32 {
    let mut T: f64 = 0.0;
    let mut p: f64 = 0.0;
    let mut v: f64 = 0.0;
    let mut hs: f64 = 0.0;
    let s13: f64 = 3.3977829547018907;// pT2s_reg1(100.0, 623.15);
    let s13s: f64 =  3.7782813395443466;//pT2s_reg1(Ps_623, 623.15);
    let sTPmax: f64 = 6.040483671712382;//pT2s_reg2(100.0, 1073.15);
    let s2ab: f64 = 7.85234039987851;//pT2s_reg2(4.0, 1073.15); // TODO： p=4 2ab s2ab

    // Left point in h-s plot
    let mut smin: f64 =-0.00858228709262268;// pT2s_reg1(100.0, 273.15);
    let mut hmin: f64 = -0.0415878259881163;//pT2h_reg1(P_MIN, 273.15);

    // Right point in h-s plot
    let mut hmax: f64 =4160.660928250124;// pT2h_reg2(P_MIN, 1073.15);
    let mut smax: f64 =  11.921055068613795;//pT2s_reg2(P_MIN, 1073.15);

    // Region 4 left and right point
    let h4l: f64 = -0.0415878259881163 ;//pT2h_reg1(P_MIN, 273.15);
    let s4l: f64 = -0.00015454959194230702;//pT2s_reg1(P_MIN, 273.15);

    let h4v: f64 = 2500.892617817172;//pT2h_reg2(P_MIN, 273.15);
    let s4v: f64 = 9.155759395224662;//pT2s_reg2(P_MIN, 273.15);

    // !!!! Check region 5 MUST On TOP !!!
    // if （s4v <= s && s<= smax） (h,s)may be setup to error region2
    if pT2s_reg5(50.0, 1073.15) < s
        && s <= pT2s_reg5(P_MIN, 2273.15)
        && pT2h_reg5(50.0, 1073.15) < h
        && h <= pT2h_reg5(P_MIN, 2273.15)
    {
        p = hs2p_reg5(h, s);
        T = ph2T_reg5(p, h);
        if 1073.15 < T && T <= 2273.15 && P_MIN <= p && p <= 50.0 {
            return 5;
        }
    };
   
    if smin <= s && s <= s13 {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h1_s(s);
        T = ps2T_reg1(100.0, s) - 0.0218;
        hmax = pT2h_reg1(100.0, T);
        if hmin <= h && h < hs {
            return 4;
        }
        if hs <= h && h <= hmax {
            return 1;
        }
    };

    if s13 < s && s <= s13s {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h1_s(s);
        let h13: f64 = hs_region_h13_s(s);
        v = ps2v_reg3(100.0, s) * (1.0 + 9.6e-5);
        T = ps2T_reg3(100.0, s) - 0.0248;
        hmax = Td2h_reg3(T, 1.0 / v);
        if hmin <= h && h < hs {
            return 4;
        }
        if hs <= h && h < h13 {
            return 1;
        }
        if h13 <= h && h <= hmax {
            return 3;
        }
    };

    if s13s < s && s <= SC_WATER {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h3a_s(s);
        v = ps2v_reg3(100.0, s) * (1.0 + 9.6e-5);
        T = ps2T_reg3(100.0, s) - 0.0248;
        hmax = Td2h_reg3(T, 1.0 / v);
      //  println!("hmin={:.6}",hmin);
       // println!("hs={:.6}",hs);
       // println!("hmax={:.6}",hmax);
       // println!("v={:.6}",v);
       // println!("T={:.6}",T);
        if hmin <= h && h < hs {
         //   println!("4");
            return 4;
        }
        if hs <= h && h <= hmax {
            return 3;
        }
    };

    if SC_WATER < s && s < 5.049096828 {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2c3b_s(s);
        v = ps2v_reg3(100.0, s) * (1.0 + 9.6e-5);
        T = ps2T_reg3(100.0, s) - 0.0248;
        hmax = Td2h_reg3(T, 1.0 / v);
        if hmin <= h && h < hs {
            return 4;
        }
        if hs <= h && h <= hmax {
            return 3;
        }
    };

    if 5.049096828 <= s && s < 5.260578707 {
        // Specific zone with 2-3 boundary in s shape
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2c3b_s(s);
        let h23max: f64 = pT2h_reg2(100.0, 863.15);
        let h23min: f64 = pT2h_reg2(Ps_623, 623.15);
        T = ps2T_reg2(100.0, s) - 0.019;
        hmax = pT2h_reg2(100.0, T);

        if hmin <= h && h < hs {
            return 4;
        }
        if hs <= h && h < h23min {
            return 3;
        }

        if h23min <= h && h < h23max {
            if hs2p_reg2c(h, s) <= B23_T2p(hs_region_t_hs(h, s))
            //hs2p_reg2c r2::region2_p_hs
            {
                return 2;
            } else {
                return 3;
            }
        }
        if h23max <= h && h <= hmax {
            return 2;
        }
    };

    if 5.260578707 <= s && s < 5.85 {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2c3b_s(s);
        T = ps2T_reg2(100.0, s) - 0.019;
        hmax = pT2h_reg2(100.0, T);
        if hmin <= h && h < hs {
            return 4;
        }
        if hs <= h && h <= hmax {
            return 2;
        }
    }

    if 5.85 <= s && s < sTPmax {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2ab_s(s);
        T = ps2T_reg2(100.0, s) - 0.019;
        hmax = pT2h_reg2(100.0, T);
        if hmin <= h && h < hs {
            return 4;
        }
        if hs <= h && h <= hmax {
            return 2;
        }
    };

    if sTPmax <= s && s < s2ab {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2ab_s(s);
        p = hs2p_reg2(h, s);
        hmax = pT2h_reg2(p, 1073.15);
        if hmin <= h && h < hs {
            return 4;
        }
        if hs <= h && h <= hmax {
            return 2;
        }
    };

    if s2ab <= s && s < s4v {
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        hs = hs_region_h2ab_s(s);
        p = hs2p_reg2(h, s);
        hmax = pT2h_reg2(p, 1073.15);
        if hmin <= h && h < hs {
            return 4;
        }
        if hs <= h && h <= hmax {
            return 2;
        }
    }

    if s4v <= s && s <= smax {
        hmin = pT2h_reg2(P_MIN, 273.15);
        p = hs2p_reg2a(h, s); //hs2p_reg2a r2::region2_p_hs
        hmax = pT2h_reg2(p, 1073.15);
        if P_MIN <= p && p <= 100.0 && hmin <= h && h <= hmax {
            return 2;
        }
    }
    
    INVALID_HS
}
