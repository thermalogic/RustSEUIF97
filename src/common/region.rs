//! Check the Region
//! * Basic input pairs :  (p,T) (p,h) (p,s) (h,s)  
//! * Extended input pairs:: (p,v) (t,v) (t,h) (t,s)
//! * lazy version: (p,h),(p,s),h,s)
//!    ph - reg 1 +47%,reg2 +5% reg3 +22% reg5 same
//!    ps - reg1 +71%, reg2 +3%, reg3 +25%,reg 5 same
//!    hs - reg +96%, reg 2 20% ,reg 3 69%, reg5 same

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

/// T in the order to fast check region,
///    p in MPa ,  T in K, returns the region
pub fn pT_sub_region(p: f64, T: f64) -> i32 {
    if p < P_MIN || p > P_MAX {
        return INVALID_P;
    }
    if T < T_MIN || T > T_MAX {
        return INVALID_T;
    }
    if T > T_MIN5 && T <= T_MAX5 && p > P_MAX5 {
        return INVALID_PT;
    }
    
  
    // up to fast check region
    if T > T_MAX3 && T <= T_MAX2 && p >= P_MIN && p <= P_MAX2 {
        return 2;
    }

    if T_MIN5 < T && T <= T_MAX5 && P_MIN5 <= p && p <= P_MAX5 {
        return 5;
    }

    if T >= T_MIN1 && T <= T_MIN3 {
        if p > P_MIN3 && p < P_MAX1 {
            return 1;
        }
        let p_s = p_saturation(T);
        if p >= p_s && p <= P_MAX1 {
            return 1;
        }
        if p < p_s && p > P_MIN {
            return 2;
        }
    };

    // T（623.15,T_MAX3)
    if T > T_MIN3 && T <= T_MAX3 {
        if p < P_MIN3 && p > P_MIN {
            return 2;
        }
        let p23 = B23_T2p(T);
        if p >= P_MIN && p <= p23 {
            return 2;
        }
        if p > p23 && p <= P_MAX3 {
            return 3;
        }
    };

    // Bottom of the  Saturation lines、critical point  to fast check the region
    // the special point (P_MIN, T_MIN) to region 1
    if (T -T_MIN).abs() < T_TOL && (p-P_MIN).abs() < P_TOL {
        return 1;
    }
   
    if T >= T_MIN && T < TC_WATER {
        let p_s: f64 = p_saturation(T);
        if (p - p_s).abs() < P_TOL {
            return 4;
        }
    }
   
    // the critical point in region 3
    if (T - TC_WATER).abs() < T_TOL && (p - PC_WATER).abs() < P_TOL {
        return 3;
    }

    INVALID_PT
}

/// Pmin -> Ps_623-> Pc-> 100MP, 3 range to check region
///  in each sub region use(hmin, hmax) to check region
///   - lazy version: -reg1 +47%, reg2 +5%  reg3 +22% reg5 same
pub fn ph_sub_region(p: f64, h: f64) -> i32 {
    if p < P_MIN || p > P_MAX {
        return INVALID_P;
    }

    if h < H_MIN || h > H_MAX {
        return INVALID_H;
    }

    let hmin: f64 = pT2h_reg1(p, 273.15);

    if P_MIN <= p && p <= Ps_623
    // Ps_623
    {
        let T_sat: f64 = T_saturation(p);
        let h_sw: f64 = pT2h_reg1(p, T_sat);
        if hmin <= h && h <= h_sw {
            return 1;
        };

        let h_ss: f64 = pT2h_reg2(p, T_sat);
        if h_sw < h && h < h_ss {
            return 4;
        };
        let h_25: f64 = pT2h_reg2(p, 1073.15);
        if h_ss <= h && h <= h_25 {
            return 2;
        };
        let hmax: f64 = pT2h_reg5(p, 2273.15);
        if h_25 < h && h <= hmax {
            return 5;
        };
    };

    if Ps_623 < p && p < PC_WATER {
        let h_13: f64 = pT2h_reg1(p, 623.15);
        if hmin <= h && h <= h_13 {
            return 1;
        }
        let h_32: f64 = pT2h_reg2(p, B23_p2T(p)); //boundaries
        if h_13 < h && h < h_32 {
            let p_34: f64 = h2p_sat_reg3(h); //boundaries
            if p < p_34 {
                return 4;
            } else {
                return 3;
            }
        };
        let h_25: f64 = pT2h_reg2(p, 1073.15);
        if h_32 <= h && h <= h_25 {
            return 2;
        }
        let hmax: f64 = pT2h_reg5(p, 2273.15);
        if h_25 < h && h <= hmax {
            return 5;
        }
    };

    if PC_WATER <= p && p <= 100.0 {
        let h_13: f64 = pT2h_reg1(p, 623.15);
        if hmin <= h && h <= h_13 {
            return 1;
        }

        let h_32: f64 = pT2h_reg2(p, B23_p2T(p));
        if h_13 < h && h < h_32 {
            return 3;
        }

        let h_25: f64 = pT2h_reg2(p, 1073.15);
        if h_32 <= h && h <= h_25 {
            return 2;
        }

        let hmax: f64 = pT2h_reg5(p, 2273.15);
        if (p <= 50.0) && (h_25 <= h && h <= hmax) {
            return 5;
        }
    };
    INVALID_VALUE
}

/// Pmin -> Ps_623-> Pc-> 100MP, 3 range to check region
///  in each sub region use(smin ,smax) to check region
///  - lazy version: -reg1 +71%, reg2 +3%, reg3 +25%, reg4%, reg5 same
pub fn ps_sub_region(p: f64, s: f64) -> i32 {
    
     if p < P_MIN || p > P_MAX {
        return INVALID_P;
    }

    if s < S_MIN || s > S_MAX {
        return INVALID_S;
    }
    let smin: f64 = pT2s_reg1(p, 273.15);

    // 1. First Range: [P_MIN ,Ps_623]
    if P_MIN <= p && p <= Ps_623 {
        let Tsat: f64 = T_saturation(p);
        let s14: f64 = pT2s_reg1(p, Tsat);
        if smin <= s && s <= s14 {
            return 1;
        }
        let s24: f64 = pT2s_reg2(p, Tsat);
        if s14 < s && s < s24 {
            return 4;
        }
        let s25: f64 = pT2s_reg2(p, 1073.15);
        if s24 <= s && s <= s25 {
            return 2;
        }
        let smax: f64 = pT2s_reg5(p, 2273.15);
        if s25 < s && s <= smax {
            return 5;
        }
    };

    // 2. Second Range: (Ps_623,PC_WATER)
    if Ps_623 < p && p < PC_WATER {
        let s13: f64 = pT2s_reg1(p, 623.15);
        if smin <= s && s <= s13 {
            return 1;
        }
        let s32: f64 = pT2s_reg2(p, B23_p2T(p));
        if s13 < s && s < s32 {
            let p34: f64 = s2p_sat_reg3(s); // boundaries;
            if p < p34 {
                return 4;
            } else {
                return 3;
            }
        };
        let s25: f64 = pT2s_reg2(p, 1073.15);
        if s32 <= s && s <= s25 {
            return 2;
        }
        let smax: f64 = pT2s_reg5(p, 2273.15);
        if s25 < s && s <= smax {
            return 5;
        }
    };

    // 3. Third: [PC_WATER,100.0]
    if PC_WATER <= p && p <= 100.0 {
        let s13: f64 = pT2s_reg1(p, 623.15);
        if smin <= s && s <= s13 {
            return 1;
        }
        let s32: f64 = pT2s_reg2(p, B23_p2T(p));
        if s13 < s && s < s32 {
            return 3;
        }
        let s25: f64 = pT2s_reg2(p, 1073.15);
        if s32 <= s && s <= s25 {
            return 2;
        }
        let smax: f64 = pT2s_reg5(p, 2273.15);
        if p <= 50.0 && s25 <= s && s <= smax {
            return 5;
        }
    };
    INVALID_VALUE
}

///  region 1,2,3,4 (smin ->smax), region 5
///  lazy version: -reg1 +96%, reg 2 same ,reg 3 69%, reg5 same  
macro_rules! define_region3_hmax_boundary_points {
    ($s:expr, $hmax:ident) => {
        let v = ps2v_reg3(100.0, $s) * (1.0 + 9.6e-5);
        let T = ps2T_reg3(100.0, $s) - 0.0248;
        $hmax = Td2h_reg3(T, 1.0 / v);
    };
}

// Start at the low entropy curves and work our way up
pub fn hs_sub_region(h: f64, s: f64) -> i32 {

    if h < H_MIN || h > H_MAX {
        return INVALID_H;
    }

    if s < S_MIN || s > S_MAX {
        return INVALID_S;
    }

    let mut T: f64 = 0.0;
    let mut p: f64 = 0.0;

    const h4l: f64 = 0.0; // pT2h_reg1(P_MIN, 273.15);
    const s4l: f64 = 0.0; // pT2s_reg1(P_MIN, 273.15);
    const h4v: f64 = 2500.89261782; //pT2h_reg2(P_MIN, 273.15);
    const s4v: f64 = 9.15575940; // pT2s_reg2(P_MIN, 273.15);

    // !!!! Check region 5 MUST BE On TOP !!!
    // if (s4v <= s && s<= smax) (h,s)may be setup to error region2
    let s_r5_1 = 6.51708829; //pT2s_reg5(50.0, 1073.15);
    let s_r5_2 = 13.90495608; //pT2s_reg5(P_MIN, 2273.15);
    let h_r5_1 = 3926.05014007; //pT2h_reg5(50.0, 1073.15);
    let h_r5_2 = 7376.98026360; //pT2h_reg5(P_MIN, 2273.15);
    if s_r5_1 < s && s <= s_r5_2 && h_r5_1 < h && h <= h_r5_2 {
        p = hs2p_reg5(h, s);
        T = ph2T_reg5(p, h);
        if 1073.15 < T && T <= 2273.15 && P_MIN <= p && p <= 50.0 {
            return 5;
        }
    };

    let mut v: f64 = 0.0;
    let mut hs: f64 = 0.0;
    // Left point in h-s plot
    let smin: f64 = 0.0; // pT2s_reg1(100.0, 273.15);
    let mut hmin: f64 = 0.0;
    // Right point in h-s plot
    let smax: f64 = 11.92105507; // pT2s_reg2(P_MIN, 1073.15);
    let mut hmax: f64 = 0.0;

    let s13: f64 = 3.39778295; //pT2s_reg1(100.0, 623.15);
    if smin <= s && s <= s13 {
        T = ps2T_reg1(100.0, s) - 0.0218;
        hmax = pT2h_reg1(100.0, T);
        hs = hs_region_h1_s(s);
        if hs <= h && h <= hmax {
            return 1;
        }
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        if hmin <= h && h < hs {
            return 4;
        };
    };

    let s13s: f64 = 3.77828134; //pT2s_reg1(Ps_623, 623.15);
    if s13 < s && s <= s13s {
        hs = hs_region_h1_s(s);
        let h13: f64 = hs_region_h13_s(s);
        if hs <= h && h < h13 {
            return 1;
        }
        define_region3_hmax_boundary_points!(s, hmax);
        if h13 <= h && h <= hmax {
            return 3;
        }
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        if hmin <= h && h < hs {
            return 4;
        }
    };

    if s13s < s && s <= SC_WATER {
        hs = hs_region_h3a_s(s);
        define_region3_hmax_boundary_points!(s, hmax);
        if hs <= h && h <= hmax {
            return 3;
        }
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        if hmin <= h && h < hs {
            return 4;
        }
    };

    if SC_WATER < s && s < 5.049096828 {
        hs = hs_region_h2c3b_s(s);
        define_region3_hmax_boundary_points!(s, hmax);
        if hs <= h && h <= hmax {
            return 3;
        }
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        if hmin <= h && h < hs {
            return 4;
        }
    };

    if 5.049096828 <= s && s < 5.260578707 {
        // Specific zone with 2-3 boundary in s shape
        hs = hs_region_h2c3b_s(s);
        let h23min: f64 = 2563.59200389; //pT2h_reg2(Ps_623, 623.15);
        if hs <= h && h < h23min {
            return 3;
        }
        T = ps2T_reg2(100.0, s) - 0.019;
        hmax = pT2h_reg2(100.0, T);
        let h23max: f64 = 2812.94206060; //pT2h_reg2(100.0, 863.15);
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
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        if hmin <= h && h < hs {
            return 4;
        }
    };

    if 5.260578707 <= s && s < 5.85 {
        hs = hs_region_h2c3b_s(s);
        T = ps2T_reg2(100.0, s) - 0.019;
        hmax = pT2h_reg2(100.0, T);
        if hs <= h && h <= hmax {
            return 2;
        }
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        if hmin <= h && h < hs {
            return 4;
        }
    }

    let sTPmax: f64 = 6.04048367; // pT2s_reg2(100.0, 1073.15);
    if 5.85 <= s && s < sTPmax {
        hs = hs_region_h2ab_s(s);
        T = ps2T_reg2(100.0, s) - 0.019;
        hmax = pT2h_reg2(100.0, T);
        if hs <= h && h <= hmax {
            return 2;
        }
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        if hmin <= h && h < hs {
            return 4;
        }
    };

    let s2ab: f64 = 7.85234040; // pT2s_reg2(4.0, 1073.15); // TODO： p=4 2ab s2ab
    if sTPmax <= s && s < s2ab {
        hs = hs_region_h2ab_s(s);
        p = hs2p_reg2(h, s);
        hmax = pT2h_reg2(p, 1073.15);
        if hs <= h && h <= hmax {
            return 2;
        }
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        if hmin <= h && h < hs {
            return 4;
        }
    };

    if s2ab <= s && s < s4v {
        hs = hs_region_h2ab_s(s);
        p = hs2p_reg2(h, s);
        hmax = pT2h_reg2(p, 1073.15);
        if hs <= h && h <= hmax {
            return 2;
        }
        hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
        if hmin <= h && h < hs {
            return 4;
        }
    }

    if s4v <= s && s <= smax {
        hmin = 2500.89261782; //pT2h_reg2(P_MIN, 273.15);
        p = hs2p_reg2a(h, s); //hs2p_reg2a r2::region2_p_hs
        hmax = pT2h_reg2(p, 1073.15);
        if P_MIN <= p && p <= 100.0 && hmin <= h && h <= hmax {
            return 2;
        }
    }
    INVALID_VALUE
}

///  Region for extended input pairs (p,v),(t,v),(t,s),(t,h)
///     Region (p,v)
pub fn pv_sub_region(p: f64, v: f64) -> i32 {
    if (p < P_MIN) || (p > P_MAX) {
        return INVALID_P;
    }

   if (v < V_MIN) || (v > V_MAX) {
        return INVALID_V;
    }

    let vt273: f64 = pT2v_reg1(p, 273.15);
    let vt1073: f64 = pT2v_reg2(p, 1073.15);
    let mut vt2273: f64 = 0.0;
    if (p >= P_MIN5) && (p <= P_MAX5) {
        vt2273 = pT2v_reg5(p, 2273.15);
    }

    let mut T1: f64 = 0.0;
    let mut vsw: f64 = 0.0;
    let mut vss: f64 = 0.0;
    let mut vt623: f64 = 0.0;
    let mut vB23: f64 = 0.0;

    if (p > 0.000611213) && (p <= P_MIN3) {
        T1 = T_saturation(p);
        vsw = pT2v_reg1(p, T1);
        vss = pT2v_reg2(p, T1);
    } else if (p > P_MIN3) && (p <= P_MAX3) {
        vt623 = pT2v_reg1(p, 623.15);
        T1 = B23_p2T(p);
        vB23 = pT2v_reg2(p, T1); //
        if p <= PC_WATER {
            T1 = T_saturation(p);
            vsw = pT2v_sat_reg3(p, T1, 0.0);
            vss = pT2v_sat_reg3(p, T1, 1.0);
        }
    }
    if ((p >= P_MIN1) && (p <= P_MIN3)) && ((v >= vt273) && (v <= vsw))
        || ((p > P_MIN3) && (p <= P_MAX1)) && ((v >= vt273) && (v <= vt623))
    {
        return 1;
    };
    if ((p >= P_MIN2) && (p <= P_MIN3) && (v > vss) && (v <= vt1073))
        || ((p > P_MIN3) && (p <= P_MAX2) && (v >= vB23) && (v <= vt1073))
    {
        return 2;
    };
    if ((p > P_MIN3) && (p <= PC_WATER) && (((v > vt623) && (v < vsw)) || ((v > vss) && (v < vB23))))
        || ((p > PC_WATER) && (p <= P_MAX3) && (v > vt623) && (v < vB23))
    {
        return 3;
    };
    if (p > P_MIN) && (p <= PC_WATER) && (v >= vsw) && (v <= vss) {
        return 4; // x is obtained in r4::region4_pair_ext::pv2x_reg4
    };
    if (p > P_MIN5) && (p <= P_MAX5) && (v > vt1073) && (v <= vt2273) {
        return 5;
    };
    INVALID_VALUE
}

///  Region (t,v)
pub fn tv_sub_region(t: f64, v: f64) -> i32 {
    let T: f64 = t + 273.15;
    if (T < 273.15) || (T > 2273.15) {
        return INVALID_T;
    }

    if (v < V_MIN) || (v > V_MAX) {
        return INVALID_V;
    }

    let mut vpmax2: f64 = 0.0;
    if (T >= T_MIN2) && (T <= T_MAX2) {
        vpmax2 = pT2v_reg2(P_MIN2, T);
    }

    let mut vpmax5: f64 = 0.0;
    let mut vp50: f64 = 0.0;
    if (T >= T_MIN5) && (T <= T_MAX5) {
        vpmax5 = pT2v_reg5(P_MIN5, T);
        vp50 = pT2v_reg5(50.0, T);
    }

    // p=100MPa
    let mut vp100: f64 = 0.0;
    if (T >= T_MIN1) && (T <= T_MAX1) {
        vp100 = pT2v_reg1(100.0, T);
    } else if (T > T_MIN3) && (T <= T_MAX3) {
        vp100 = pT2v_reg3(100.0, T);
    } else if (T > T_MAX3) && (T <= T_MAX2) {
        vp100 = pT2v_reg2(100.0, T);
    }

    let mut vsw: f64 = 0.0;
    let mut vss: f64 = 0.0;
    let mut vB23: f64 = 0.0;

    if (T >= T_MIN1) && (T <= T_MAX1) {
        let p_s = p_saturation(T);
        vsw = pT2v_reg1(p_s, T);
        vss = pT2v_reg2(p_s, T);
    } else if (T > T_MIN3) && (T <= T_MAX3) {
        let p_b23 = B23_T2p(T);
        vB23 = pT2v_reg2(p_b23, T); //
        if T <= TC_WATER {
            let p_s = p_saturation(T);
            vsw = pT2v_sat_reg3(p_s, T, 0.0);
            vss = pT2v_sat_reg3(p_s, T, 1.0);
        };
    };

    //
    if (T >= T_MIN1) && (T <= T_MAX1) && (v < vsw) && (v > vp100) {
        return 1;
    };
    if ((T >= T_MIN2) && (T <= T_MAX1) && (v > vss) && (v < vpmax2))
        || ((T > T_MAX1) && (T <= T_MAX3) && (v >= vB23) && (v < vpmax2))
    {
        return 2;
    };
    if ((T > T_MIN3) && (T <= TC_WATER) && (((v > vss) && (v < vB23)) || ((v >= vp100) && (v < vsw))))
        || ((T > TC_WATER) && (T <= T_MAX3) && (v >= vp100) && (v < vB23))
    {
        return 3;
    };

    if (T >= T_MIN1) && (T <= TC_WATER) && (v >= vsw) && (v <= vss) {
        return 4; // x is obtained in r4::region4_pair_ext::Tv2x_reg4
    };

    if (T > T_MIN5) && (T <= T_MAX5) && (v >= vp50) && (v <= vpmax5) {
        return 5;
    };
    INVALID_VALUE
}

/// Region (t,h)
pub fn th_sub_region(t: f64, h: f64) -> i32 {
    let T: f64 = t + 273.15;
    if T < 273.15 || T > 2273.15 {
        return INVALID_T;
    }
    if h < H_MIN || h > H_MAX {
        return INVALID_H;
    }

    let mut hpmax2: f64 = 0.0;
    if (T >= T_MIN2) && (T <= T_MAX2) {
        hpmax2 = pT2h_reg2(P_MIN2, T);
    }

    let mut hpm50: f64 = 0.0;
    let mut hpmax5: f64 = 0.0;
    if (T > T_MIN5) && (T <= T_MAX5) {
        hpm50 = pT2h_reg5(P_MAX5, T);
        hpmax5 = pT2h_reg5(P_MIN5, T);
    }
    // p=100MPa
    let mut hp100: f64 = 0.0;
    if (T >= T_MIN1) && (T <= T_MAX1) {
        hp100 = pT2h_reg1(100.0, T);
    } else if (T > T_MIN3) && (T <= T_MAX3) {
        hp100 = pt_reg3(100.0, T - 273.15, OH);
    }

    let mut hsw: f64 = 0.0;
    let mut hss: f64 = 0.0;
    let mut p_s: f64 = 0.0;
    let mut hB23: f64 = 0.0;
    if (T >= T_MIN1) && (T <= T_MAX1) {
        p_s = p_saturation(T);
        hsw = pT2h_reg1(p_s, T);
        hss = pT2h_reg2(p_s, T);
    } else if (T > T_MIN3) && (T <= T_MAX3) {
        let p_b23 = B23_T2p(T);
        hB23 = pT2h_reg2(p_b23, T);
        if T <= TC_WATER {
            p_s = p_saturation(T);
            let mut v: f64 = pT2v_sat_reg3(p_s, T, 0.0);
            hsw = Td2h_reg3(T, 1.0 / v);
            v = pT2v_sat_reg3(p_s, T, 1.0);
            hss = Td2h_reg3(T, 1.0 / v);
        };
    };

    // if T is very small , P^ -> h^, h1 > hsat_water(T)
    if (T >= T_MIN1) && (T < (250.0 + T_MIN1)) && (h > hsw) && (h <= hp100) {
        return 1;
    }
    //  if p>Ps, p^, h--, Hmin,p^,h++, h<hsat_water(T)
    if (T >= (250.0 + 273.15)) && (T <= T_MAX1) && (h < hsw) && (h <= hp100) {
        let mut Hmin: f64 = hsw;
        let mut p11: f64 = p_s + 0.1;
        let mut Hmid: f64 = pT2h_reg1(p11, T);
        if Hmid < Hmin {
            let mut HminFounded: bool = false;

            while !HminFounded {
                p11 += 0.1;
                Hmid = pT2h_reg1(p11, T);
                if Hmid < Hmin {
                    Hmin = Hmid;
                } else {
                    HminFounded = true;
                }
            }
            if h >= Hmin {
                return 1;
            } else {
                return INVALID_VALUE;
            };
        } else {
            return INVALID_VALUE;
        }
    };

    if (T >= T_MIN2) && (T <= T_MAX1) && (h > hss) && (h < hpmax2) {
        return 2;
    };
    if (T > T_MAX1) && (T <= T_MAX3) && (h >= hB23) && (h < hpmax2) {
        return 2;
    };
    if (T > T_MIN3) && (T <= TC_WATER) && (((h > hss) && (h < hB23)) || ((h >= hp100) && (h < hsw)))
        || ((T > TC_WATER) && (T <= T_MAX3) && (h >= hp100) && (h < hB23))
    {
        return 3;
    };

    //  (t,h)
    if (T >= T_MIN1) && (T <= TC_WATER) && (h >= hsw) && (h <= hss) {
        return 4;
    };

    if (T > T_MIN5) && (T <= T_MAX5) && (h >= hpm50) && (h <= hpmax5) {
        return 5;
    };
    INVALID_VALUE
}

///  Region (t,s)
pub fn ts_sub_region(t: f64, s: f64) -> i32 {
    let T: f64 = t + 273.15;
    if (T < 273.15) || (T > 2273.15) {
        return INVALID_T;
    }

    if (s > S_MAX) || (s < S_MIN) {
        return INVALID_S;
    }

    let mut spmax2: f64 = 0.0;
    if (T >= T_MIN2) && (T <= T_MAX2) {
        spmax2 = pT2s_reg2(P_MIN2, T);
    }

    let mut spmax5: f64 = 0.0;
    let mut sp50: f64 = 0.0;
    if (T >= T_MIN5) && (T <= T_MAX5) {
        spmax5 = pT2s_reg5(P_MIN5, T);
        sp50 = pT2s_reg5(50.0, T);
    }
    // p=100MPa
    let mut sp100: f64 = 0.0;
    if (T >= T_MIN1) && (T <= T_MAX1) {
        sp100 = pT2s_reg1(100.0, T);
    } else if (T > T_MIN3) && (T <= T_MAX3) {
        sp100 = pT_reg3(100.0, T, OS);
    } else if (T > T_MAX3) && (T <= T_MAX2) {
        sp100 = pT2s_reg2(100.0, T);
    }

    let mut p_s: f64 = 0.0;
    let mut p_b23: f64 = 0.0;
    let mut ssw: f64 = 0.0;
    let mut sss: f64 = 0.0;
    let mut sB23: f64 = 0.0;
    if (T >= T_MIN1) && (T <= T_MAX1) {
        p_s = p_saturation(T);
        ssw = pT2s_reg1(p_s, T);
        sss = pT2s_reg2(p_s, T);
    } else if (T > T_MIN3) && (T <= T_MAX3) {
        p_b23 = B23_T2p(T);
        sB23 = pT2s_reg2(p_b23, T); //
        if T <= TC_WATER {
            p_s = p_saturation(T);
            let mut v: f64 = pT2v_sat_reg3(p_s, T, 0.0);
            ssw = Td2s_reg3(T, 1.0 / v);
            v = pT2v_sat_reg3(p_s, T, 1.0);
            sss = Td2s_reg3(T, 1.0 / v);
        };
    };
    //
    if (T >= T_MIN1) && (T <= T_MAX1) && (s < ssw) && (s > sp100) {
        return 1;
    };
    if ((T >= T_MIN2) && (T <= T_MAX1) && (s > sss) && (s < spmax2))
        || ((T > T_MAX1) && (T <= T_MAX3) && (s >= sB23) && (s < spmax2))
    {
        return 2;
    };
    if ((T > T_MIN3) && (T <= TC_WATER) && (((s > sss) && (s < sB23)) || ((s >= sp100) && (s < ssw))))
        || ((T > TC_WATER) && (T <= T_MAX3) && (s >= sp100) && (s < sB23))
    {
        return 3;
    };
    if (T >= T_MIN1) && (T <= TC_WATER) && (s >= ssw) && (s <= sss) {
        return 4;
    };
    if (T > T_MIN5) && (T <= T_MAX5) && (s >= sp50) && (s <= spmax5) {
        return 5;
    };
    INVALID_VALUE
}
