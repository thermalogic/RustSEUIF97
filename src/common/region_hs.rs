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

   if s < S_MIN || s > S_MAX {
        return INVALID_S;
    }

    if h < H_MIN || h > H_MAX {
        return INVALID_H;
    }

    let s13: f64 = 3.39778295; //pT2s_reg1(100.0, 623.15);
    let s13s: f64 = 3.77828134; //pT2s_reg1(Ps_623, 623.15);
    let sTPmax: f64 = 6.04048367; // pT2s_reg2(100.0, 1073.15);
    let s2ab: f64 = 7.85234040; // pT2s_reg2(4.0, 1073.15); // TODO： p=4 2ab s2ab
    // Left point in h-s plot
    let smin: f64 = 0.0; // pT2s_reg1(100.0, 273.15);
    let hmin: f64 = pT2h_reg1(P_MIN,273.15)
    // Right point in h-s plot
    hmax = pT2h_reg2(P_MIN,1073.15)
    smax = pT2s_reg2(P_MIN,1073.15)
    // Region 4 left and right point
    const h4l: f64 = 0.0; // pT2h_reg1(P_MIN, 273.15);
    const s4l: f64 = 0.0; // pT2s_reg1(P_MIN, 273.15);
    const h4v: f64 = 2500.89261782; //pT2h_reg2(P_MIN, 273.15);
    const s4v: f64 = 9.15575940; // pT2s_reg2(P_MIN, 273.15);



    INVALID_VALUE
}

