//!  The C API: stdcall - Win32 API functions
//!
use crate::algo::*;
use crate::common::*;
use crate::r1::*;
use crate::r2::*;
use crate::r3::*;
use crate::r4::*;
use crate::r5::*;

/// double pt(double p,double t,short o_id)- the property of `o_id` (thermodynamic,transport,etc)  
#[no_mangle]
pub unsafe extern "stdcall" fn pt(p: f64, t: f64, o_id: i32) -> f64 {
    let T: f64 = t + 273.15;
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OT => return t,
        _ => pair_properties(p, T, o_id, pT_sub_region, pT_reg1, pT_reg2, pT_reg3, pT_reg4, pT_reg5, reg),
    }
}

/// double ph(double p,double h,short o_id)- the property of `o_id` (thermodynamic,transport,etc)  
#[no_mangle]
pub unsafe extern "stdcall" fn ph(p: f64, h: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OH => return h,
        _ => pair_properties(p, h, o_id, ph_sub_region, ph_reg1, ph_reg2, ph_reg3, ph_reg4, ph_reg5, reg),
    }
}

/// double ps(double p,double s,short o_id)- the property of `o_id` (thermodynamic,transport,etc)  
#[no_mangle]
pub unsafe extern "stdcall" fn ps(p: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OS => return s,
        _ => pair_properties(p, s, o_id, ps_sub_region, ps_reg1, ps_reg2, ps_reg3, ps_reg4, ps_reg5, reg),
    }
}

/// double hs(double h,double s,short o_id)- the property of `o_id` (thermodynamic,transport,etc)  
#[no_mangle]
pub unsafe extern "stdcall" fn hs(h: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OH => return h,
        OS => return s,
        _ => pair_properties(h, s, o_id, hs_sub_region, hs_reg1, hs_reg2, hs_reg3, hs_reg4, hs_reg5, reg),
    }
}

/// double px(double p,double x,short o_id) - the property of `o_id` (thermodynamic)
#[no_mangle]
pub unsafe extern "stdcall" fn px(p: f64, x: f64, o_id: i32) -> f64 {
    if p > P_MAX4 || p < P_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OP => return p,
        OX => return x,
        _ => return px_reg4(p, x, o_id),
    }
}

/// double tx(double t,double x,short o_id) - the property of `o_id` (thermodynamic)
#[no_mangle]
pub unsafe extern "stdcall" fn tx(t: f64, x: f64, o_id: i32) -> f64 {
    match o_id {
        OT => return t,
        OX => return x,
        _ => {
            let T: f64 = t + 273.15;
            if T > T_MAX4 || T < T_MIN4 || x > 1.0 || x < 0.0 {
                return INVALID_VALUE as f64;
            }
            Tx_reg4(T, x, o_id)
        }
    }
}

/// double pv(double p,double v,short o_id)- the property of `o_id` (thermodynamic,transport,etc)  
#[no_mangle]
pub unsafe extern "stdcall" fn pv(p: f64, v: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OV => return v,
        _ => pair_properties(p, v, o_id, pv_sub_region, pv_reg1, pv_reg2, pv_reg3, pv_reg4, pv_reg5, reg),
    }
}

/// double tv(double t,double v,short o_id)- the property of `o_id` (thermodynamic,transport,etc)  
#[no_mangle]
pub unsafe extern "stdcall" fn tv(t: f64, v: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => return t,
        OV => return v,
        _ => pair_properties(t, v, o_id, tv_sub_region, tv_reg1, tv_reg2, tv_reg3, tv_reg4, tv_reg5, reg),
    }
}

/// double th(double t,double h,short o_id)- the property of `o_id` (thermodynamic,transport,etc)  
#[no_mangle]
pub unsafe extern "stdcall" fn th(t: f64, h: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => return t,
        OH => return h,
        _ => pair_properties(t, h, o_id, th_sub_region, th_reg1, th_reg2, th_reg3, th_reg4, th_reg5, reg),
    }
}

/// double ts(double t,double s,short o_id)- the property of `o_id` (thermodynamic,transport,etc)  
#[no_mangle]
pub unsafe extern "stdcall" fn ts(t: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => return t,
        OS => return s,
        _ => pair_properties(t, s, o_id, ts_sub_region, ts_reg1, ts_reg2, ts_reg3, ts_reg4, ts_reg5, reg),
    }
}

/// double hx(double h,double x,short o_id)- the property of `o_id` (thermodynamic)  
#[no_mangle]
pub unsafe extern "stdcall" fn hx(h: f64, x: f64, o_id: i32) -> f64 {
    if h > H_MAX4 || h < H_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OH => return h,
        OX => return x,
        _ => hx_reg4(h, x, o_id),
    }
}

/// double sx(double s,double x,short o_id)- the property of `o_id` (thermodynamic)
#[no_mangle]
pub unsafe extern "stdcall" fn sx(s: f64, x: f64, o_id: i32) -> f64 {
    if s > S_MAX4 || s < S_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OS => return s,
        OX => return x,
        _ => sx_reg4(s, x, o_id),
    }
}

// Convenience Functions (Direct Output)

/// double pt2h(double p,double t) - Calculate specific enthalpy from pressure and temperature
#[no_mangle]
pub unsafe extern "stdcall" fn pt2h(p: f64, t: f64) -> f64 {
    pt(p, t, OH)
}

/// double pt2s(double p,double t) - Calculate specific entropy from pressure and temperature
#[no_mangle]
pub unsafe extern "stdcall" fn pt2s(p: f64, t: f64) -> f64 {
    pt(p, t, OS)
}

/// double pt2v(double p,double t) - Calculate specific volume from pressure and temperature
#[no_mangle]
pub unsafe extern "stdcall" fn pt2v(p: f64, t: f64) -> f64 {
    pt(p, t, OV)
}

/// double pt2x(double p,double t) - Calculate steam quality from pressure and temperature
#[no_mangle]
pub unsafe extern "stdcall" fn pt2x(p: f64, t: f64) -> f64 {
    pt(p, t, OX)
}

/// double ph2t(double p,double h) - Calculate temperature from pressure and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn ph2t(p: f64, h: f64) -> f64 {
    ph(p, h, OT)
}

/// double ph2s(double p,double h) - Calculate specific entropy from pressure and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn ph2s(p: f64, h: f64) -> f64 {
    ph(p, h, OS)
}

/// double ph2v(double p,double h) - Calculate specific volume from pressure and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn ph2v(p: f64, h: f64) -> f64 {
    ph(p, h, OV)
}

/// double ph2x(double p,double h) - Calculate steam quality from pressure and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn ph2x(p: f64, h: f64) -> f64 {
    ph(p, h, OX)
}

/// double ps2t(double p,double s) - Calculate temperature from pressure and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ps2t(p: f64, s: f64) -> f64 {
    ps(p, s, OT)
}

/// double ps2h(double p,double s) - Calculate specific enthalpy from pressure and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ps2h(p: f64, s: f64) -> f64 {
    ps(p, s, OH)
}

/// double ps2v(double p,double s) - Calculate specific volume from pressure and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ps2v(p: f64, s: f64) -> f64 {
    ps(p, s, OV)
}

/// double ps2x(double p,double s) - Calculate steam quality from pressure and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ps2x(p: f64, s: f64) -> f64 {
    ps(p, s, OX)
}

/// double pv2t(double p,double v) - Calculate temperature from pressure and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn pv2t(p: f64, v: f64) -> f64 {
    pv(p, v, OT)
}

/// double pv2h(double p,double v) - Calculate specific enthalpy from pressure and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn pv2h(p: f64, v: f64) -> f64 {
    pv(p, v, OH)
}

/// double pv2s(double p,double v) - Calculate specific entropy from pressure and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn pv2s(p: f64, v: f64) -> f64 {
    pv(p, v, OS)
}

/// double pv2x(double p,double v) - Calculate steam quality from pressure and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn pv2x(p: f64, v: f64) -> f64 {
    pv(p, v, OX)
}

/// double hs2p(double h,double s) - Calculate pressure from enthalpy and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn hs2p(h: f64, s: f64) -> f64 {
    hs(h, s, OP)
}

/// double hs2t(double h,double s) - Calculate temperature from enthalpy and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn hs2t(h: f64, s: f64) -> f64 {
    hs(h, s, OT)
}

/// double hs2v(double h,double s) - Calculate specific volume from enthalpy and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn hs2v(h: f64, s: f64) -> f64 {
    hs(h, s, OV)
}

/// double hs2x(double h,double s) - Calculate steam quality from enthalpy and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn hs2x(h: f64, s: f64) -> f64 {
    hs(h, s, OX)
}

/// double th2p(double t,double h) - Calculate pressure from temperature and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn th2p(t: f64, h: f64) -> f64 {
    th(t, h, OP)
}

/// double th2s(double t,double h) - Calculate specific entropy from temperature and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn th2s(t: f64, h: f64) -> f64 {
    th(t, h, OS)
}

/// double th2v(double t,double h) - Calculate specific volume from temperature and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn th2v(t: f64, h: f64) -> f64 {
    th(t, h, OV)
}

/// double th2x(double t,double h) - Calculate steam quality from temperature and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn th2x(t: f64, h: f64) -> f64 {
    th(t, h, OX)
}

/// double ts2p(double t,double s) - Calculate pressure from temperature and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ts2p(t: f64, s: f64) -> f64 {
    ts(t, s, OP)
}

/// double ts2h(double t,double s) - Calculate specific enthalpy from temperature and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ts2h(t: f64, s: f64) -> f64 {
    ts(t, s, OH)
}

/// double ts2v(double t,double s) - Calculate specific volume from temperature and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ts2v(t: f64, s: f64) -> f64 {
    ts(t, s, OV)
}

/// double ts2x(double t,double s) - Calculate steam quality from temperature and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ts2x(t: f64, s: f64) -> f64 {
    ts(t, s, OX)
}

/// double tv2p(double t,double v) - Calculate pressure from temperature and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn tv2p(t: f64, v: f64) -> f64 {
    tv(t, v, OP)
}

/// double tv2h(double t,double v) - Calculate specific enthalpy from temperature and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn tv2h(t: f64, v: f64) -> f64 {
    tv(t, v, OH)
}

/// double tv2s(double t,double v) - Calculate specific entropy from temperature and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn tv2s(t: f64, v: f64) -> f64 {
    tv(t, v, OS)
}

/// double tv2x(double t,double v) - Calculate steam quality from temperature and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn tv2x(t: f64, v: f64) -> f64 {
    tv(t, v, OX)
}

/// double px2t(double p,double x) - Calculate temperature from pressure and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn px2t(p: f64, x: f64) -> f64 {
    px(p, x, OT)
}

/// double px2h(double p,double x) - Calculate specific enthalpy from pressure and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn px2h(p: f64, x: f64) -> f64 {
    px(p, x, OH)
}

/// double px2s(double p,double x) - Calculate specific entropy from pressure and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn px2s(p: f64, x: f64) -> f64 {
    px(p, x, OS)
}

/// double px2v(double p,double x) - Calculate specific volume from pressure and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn px2v(p: f64, x: f64) -> f64 {
    px(p, x, OV)
}

/// double tx2p(double t,double x) - Calculate pressure from temperature and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn tx2p(t: f64, x: f64) -> f64 {
    tx(t, x, OP)
}

/// double tx2h(double t,double x) - Calculate specific enthalpy from temperature and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn tx2h(t: f64, x: f64) -> f64 {
    tx(t, x, OH)
}

/// double tx2s(double t,double x) - Calculate specific entropy from temperature and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn tx2s(t: f64, x: f64) -> f64 {
    tx(t, x, OS)
}

/// double tx2v(double t,double x) - Calculate specific volume from temperature and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn tx2v(t: f64, x: f64) -> f64 {
    tx(t, x, OV)
}

/// double hx2p(double h,double x) - Calculate pressure from enthalpy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn hx2p(h: f64, x: f64) -> f64 {
    hx(h, x, OP)
}

/// double hx2t(double h,double x) - Calculate temperature from enthalpy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn hx2t(h: f64, x: f64) -> f64 {
    hx(h, x, OT)
}

/// double hx2s(double h,double x) - Calculate specific entropy from enthalpy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn hx2s(h: f64, x: f64) -> f64 {
    hx(h, x, OS)
}

/// double hx2v(double h,double x) - Calculate specific volume from enthalpy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn hx2v(h: f64, x: f64) -> f64 {
    hx(h, x, OV)
}

/// double sx2p(double s,double x) - Calculate pressure from entropy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn sx2p(s: f64, x: f64) -> f64 {
    sx(s, x, OP)
}

/// double sx2t(double s,double x) - Calculate temperature from entropy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn sx2t(s: f64, x: f64) -> f64 {
    sx(s, x, OT)
}

/// double sx2h(double s,double x) - Calculate specific enthalpy from entropy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn sx2h(s: f64, x: f64) -> f64 {
    sx(s, x, OH)
}

/// double sx2v(double s,double x) - Calculate specific volume from entropy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn sx2v(s: f64, x: f64) -> f64 {
    sx(s, x, OV)
}
