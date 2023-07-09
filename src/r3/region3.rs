/*---------------------------------------------------------------------------
  Region 3 API:   IAPWS-IF97 and supp release

    Basic Equation of  IAPWS -IF 97 Region3
      *  IAPWS-IF97, R7-97(2012)
          (T,d)-->p,h,u,s,cp,w

    Backward Equation for Region 3:
       * IAPWS-IF97-S03rev
           (p,h) -T,v
            (p,s)->T,v

       *  IAPWS-IF97-S04rev
             (h,s)->p

           Supp-phs3-2014.pdf文档中其他区域判断公式在boundarie中

        * IAPWS-IF97-S05rev
              (p,T)->d

 Author: Maohua Cheng
---------------------------------------------------------------------------*/
use crate::common::constant::*;
use crate::r3::region3_Td::*;
use crate::r3::region3_Tv_phps::*;
use crate::r3::region3_v_pT::*;
use crate::r3::region3_p_hs::*;


pub fn Td_reg3(T: f64, d: f64, o_id: i32) -> f64 {
    match o_id {
        OT => T,
        OD => d,
        OV => 1.0 / d,
        OP => Td2p_reg3(T, d),
        OH => Td2h_reg3(T, d),
        OS => Td2s_reg3(T, d),
        OU => Td2u_reg3(T, d),
        OCV => Td2cv_reg3(T, d),
        OCP => Td2cp_reg3(T, d),
        OW => Td2w_reg3(T, d),
        _ => INVALID_OUTID as f64,
    }
}

// --------------  region3 --------------------------
pub fn td_reg3(t: f64, d: f64, o_id: i32) -> f64 {
    Td_reg3(t + K, d, o_id)
}

pub fn pT_reg3(p: f64, T: f64, o_id: i32) -> f64 {
    if o_id == OP {
        return p;
    }
    if o_id == OT {
        return T;
    }
    let v: f64 = pT2v_reg3(p, T);
    if o_id == OV {
        return v;
    }
    let d: f64 = 1.0 / v;
    if o_id == OD {
        return d;
    }
    return Td_reg3(T, d, o_id);
}

pub fn pt_reg3(p: f64, t: f64, o_id: i32) -> f64 {
    if o_id == OP {
        return p;
    }
    if o_id == OT {
        return t;
    }
    let T: f64 = t + 273.15;
    pT_reg3(p, T, o_id)
}

pub fn ph_reg3(p: f64, h: f64, o_id: i32) -> f64 {
    if o_id == OP {
        return p;
    };
    if o_id == OH {
        return h;
    };
    let v: f64 = ph2v_reg3(p, h);
    if o_id == OV {
        return v;
    };
    let d: f64 = 1.0 / v;
    if o_id == OD {
        return d;
    };
    let T: f64 = ph2T_reg3(p, h);
    if o_id == OT {
        return T - K;
    };
    return Td_reg3(T, d, o_id);
}

pub fn ps_reg3(p:f64,s:f64, o_id:i32)->f64
{
    if o_id==OP
    {   return p;};
    if o_id==OS 
    {   return s;   };

    let v: f64 = ps2v_reg3(p, s);
    if o_id == OV {
        return v;
    };
    let d: f64 = 1.0 / v;
    if o_id == OD {
        return d;
    };
    let T: f64 = ps2T_reg3(p, s);
    if o_id == OT {
        return T - K;
    };
    return  Td_reg3(T, d, o_id);
}

pub fn hs_reg3(h:f64, s:f64,o_id:i32)->f64
{
    if o_id==OH
     {  return h;}
    if o_id==OS
    {   return s;}

    let p:f64= hs2p_reg3(h,s);
    if o_id==OP
    {   return p;}
    let T:f64=ph2T_reg3(p,h);
    if o_id==OT
    {   return T-273.15;}
    let v:f64=ph2v_reg3(p,h);
    if o_id==OV
    {  return v;}
    let d=1.0/v;
    if o_id==OD
    {  return d;}
    return Td_reg3(T,d,o_id);
}
