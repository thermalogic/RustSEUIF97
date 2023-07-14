//！  Region 5 API
//！     (p,T) (p,h) (p,s) (h,s)
   
use crate::common::constant::*;

use crate::r5::region5_pT::*;
use crate::r5::region5_backward::*;

pub fn pT_reg5(p: f64, T: f64, o_id: i32) -> f64
{
    match o_id {
        OV => pT2v_reg5(p, T),
        OD => 1.0 / pT2v_reg5(p, T),
        OH => pT2h_reg5(p, T),
        OS => pT2s_reg5(p, T),
        OU => pT2u_reg5(p, T),
        OCV => pT2cv_reg5(p, T),
        OCP => pT2cp_reg5(p, T),
        OW => pT2w_reg5(p, T),
        // OK => pT2k_reg5(p, T), //isentropic exponent
        _ => INVALID_OUTID as f64,
    }
}

pub fn pt_reg5(p:f64,t:f64,o_id:i32) -> f64 {
    let T:f64=t+K;
    pT_reg5(p,T,o_id)
}

pub fn ph_reg5(p: f64, h: f64, o_id: i32) -> f64
{
    let T: f64 = ph2T_reg5(p, h);
    if o_id == OT {
        return T-273.15;
    };
    return pT_reg5(p, T, o_id);
}

pub fn ps_reg5(p: f64, s: f64, o_id: i32) -> f64
{
    let T: f64 = ps2T_reg5(p, s);
    if o_id == OT {
        return T-273.15;
    };
    return pT_reg5(p, T, o_id);
}


pub fn hs_reg5(h: f64, s: f64, o_id: i32) -> f64
{
    let p: f64 = hs2p_reg5(h, s);
    if o_id==OP
    {
        return p;
    }
    return ph_reg5(p, h, o_id);
}
