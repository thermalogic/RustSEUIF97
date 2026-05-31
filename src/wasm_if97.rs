use wasm_bindgen::prelude::*;

use crate::algo::*;
use crate::common::*;
use crate::r1::*;
use crate::r2::*;
use crate::r3::*;
use crate::r4::*;
use crate::r5::*;

#[wasm_bindgen]
pub fn pt(p: f64, t: f64, o_id: i32) -> f64 {
    let T: f64 = t + 273.15;
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OT => return t,
        _ => pair_properties(p, T, o_id, pT_sub_region, pT_reg1, pT_reg2, pT_reg3, pT_reg4, pT_reg5, reg),
    }
}

#[wasm_bindgen]
pub fn ph(p: f64, h: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OH => return h,
        _ => pair_properties(p, h, o_id, ph_sub_region, ph_reg1, ph_reg2, ph_reg3, ph_reg4, ph_reg5, reg),
    }
}

#[wasm_bindgen]
pub fn ps(p: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OS => return s,
        _ => pair_properties(p, s, o_id, ps_sub_region, ps_reg1, ps_reg2, ps_reg3, ps_reg4, ps_reg5, reg),
    }
}

#[wasm_bindgen]
pub fn hs(h: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OH => return h,
        OS => return s,
        _ => pair_properties(h, s, o_id, hs_sub_region, hs_reg1, hs_reg2, hs_reg3, hs_reg4, hs_reg5, reg),
    }
}

#[wasm_bindgen]
pub fn px(p: f64, x: f64, o_id: i32) -> f64 {
    if p > P_MAX4 || p < P_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OP => return p,
        OX => return x,
        _ => return px_reg4(p, x, o_id),
    }
}

#[wasm_bindgen]
pub fn tx(t: f64, x: f64, o_id: i32) -> f64 {
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

#[wasm_bindgen]
pub fn pv(p: f64, v: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OV => return v,
        _ => pair_properties(p, v, o_id, pv_sub_region, pv_reg1, pv_reg2, pv_reg3, pv_reg4, pv_reg5, reg),
    }
}

#[wasm_bindgen]
pub fn tv(t: f64, v: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => return t,
        OV => return v,
        _ => pair_properties(t, v, o_id, tv_sub_region, tv_reg1, tv_reg2, tv_reg3, tv_reg4, tv_reg5, reg),
    }
}

#[wasm_bindgen]
pub fn th(t: f64, h: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => return t,
        OH => return h,
        _ => pair_properties(t, h, o_id, th_sub_region, th_reg1, th_reg2, th_reg3, th_reg4, th_reg5, reg),
    }
}

#[wasm_bindgen]
pub fn ts(t: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => return t,
        OS => return s,
        _ => pair_properties(t, s, o_id, ts_sub_region, ts_reg1, ts_reg2, ts_reg3, ts_reg4, ts_reg5, reg),
    }
}

#[wasm_bindgen]
pub fn hx(h: f64, x: f64, o_id: i32) -> f64 {
    if h > H_MAX4 || h < H_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OH => return h,
        OX => return x,
        _ => hx_reg4(h, x, o_id),
    }
}

#[wasm_bindgen]
pub fn sx(s: f64, x: f64, o_id: i32) -> f64 {
    if s > S_MAX4 || s < S_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OS => return s,
        OX => return x,
        _ => sx_reg4(s, x, o_id),
    }
}

// the convenience functions (direct output)
// p,t->h,s,v,x
#[wasm_bindgen]
pub fn pt2h(p: f64, t: f64) -> f64 {
    crate::rust_if97::pt(p, t, OH)
}

#[wasm_bindgen]
pub fn pt2s(p: f64, t: f64) -> f64 {
    crate::rust_if97::pt(p, t, OS)
}

#[wasm_bindgen]
pub fn pt2v(p: f64, t: f64) -> f64 {
    crate::rust_if97::pt(p, t, OV)
}

#[wasm_bindgen]
pub fn pt2x(p: f64, t: f64) -> f64 {
    crate::rust_if97::pt(p, t, OX)
}

// p,h->t,s,v,x
#[wasm_bindgen]
pub fn ph2t(p: f64, h: f64) -> f64 {
    crate::rust_if97::ph(p, h, OT)
}

#[wasm_bindgen]
pub fn ph2s(p: f64, h: f64) -> f64 {
    crate::rust_if97::ph(p, h, OS)
}

#[wasm_bindgen]
pub fn ph2v(p: f64, h: f64) -> f64 {
    crate::rust_if97::ph(p, h, OV)
}

#[wasm_bindgen]
pub fn ph2x(p: f64, h: f64) -> f64 {
    crate::rust_if97::ph(p, h, OX)
}
// p,s->t,h,v,x
#[wasm_bindgen]
pub fn ps2t(p: f64, s: f64) -> f64 {
    crate::rust_if97::ps(p, s, OT)
}

#[wasm_bindgen]
pub fn ps2h(p: f64, s: f64) -> f64 {
    crate::rust_if97::ps(p, s, OH)
}

#[wasm_bindgen]
pub fn ps2v(p: f64, s: f64) -> f64 {
    crate::rust_if97::ps(p, s, OV)
}

#[wasm_bindgen]
pub fn ps2x(p: f64, s: f64) -> f64 {
    crate::rust_if97::ps(p, s, OX)
}

// p,v->t,h,s,x
#[wasm_bindgen]
pub fn pv2t(p: f64, v: f64) -> f64 {
    crate::rust_if97::pv(p, v, OT)
}

#[wasm_bindgen]
pub fn pv2h(p: f64, v: f64) -> f64 {
    crate::rust_if97::pv(p, v, OH)
}

#[wasm_bindgen]
pub fn pv2s(p: f64, v: f64) -> f64 {
    crate::rust_if97::pv(p, v, OS)
}

#[wasm_bindgen]
pub fn pv2x(p: f64, v: f64) -> f64 {
    crate::rust_if97::pv(p, v, OX)
}

// h,s->p,t,v,x
#[#[wasm_bindgen]
pub fn hs2p(h: f64, s: f64) -> f64 {
    crate::rust_if97::hs(h, s, OP)
}

#[wasm_bindgen]
pub fn hs2t(h: f64, s: f64) -> f64 {
    crate::rust_if97::hs(h, s, OT)
}

#[wasm_bindgen]
pub fn hs2v(h: f64, s: f64) -> f64 {
    crate::rust_if97::hs(h, s, OV)
}

#[wasm_bindgen]
pub fn hs2x(h: f64, s: f64) -> f64 {
    crate::rust_if97::hs(h, s, OX)
}

// t,h->p,s,v,x
#[wasm_bindgen]
pub fn th2p(t: f64, h: f64) -> f64 {
    crate::rust_if97::th(t, h, OP)
}

#[wasm_bindgen]
pub fn th2s(t: f64, h: f64) -> f64 {
    crate::rust_if97::th(t, h, OS)
}

#[wasm_bindgen]
pub fn th2v(t: f64, h: f64) -> f64 {
    crate::rust_if97::th(t, h, OV)
}

#[wasm_bindgen]
pub fn th2x(t: f64, h: f64) -> f64 {
    crate::rust_if97::th(t, h, OX)
}

// t,s->p,h,v,x
#[wasm_bindgen]
pub fn ts2p(t: f64, s: f64) -> f64 {
    crate::rust_if97::ts(t, s, OP)
}

#[wasm_bindgen]
pub fn ts2h(t: f64, s: f64) -> f64 {
    crate::rust_if97::ts(t, s, OH)
}

#[wasm_bindgen]
pub fn ts2v(t: f64, s: f64) -> f64 {
    crate::rust_if97::ts(t, s, OV)
}

#[wasm_bindgen]
pub fn ts2x(t: f64, s: f64) -> f64 {
    crate::rust_if97::ts(t, s, OX)
}

// t,v->p,h,s,x
#[wasm_bindgen]
pub fn tv2p(t: f64, v: f64) -> f64 {
    crate::rust_if97::tv(t, v, OP)
}

#[wasm_bindgen]
pub fn tv2h(t: f64, v: f64) -> f64 {
    crate::rust_if97::tv(t, v, OH)
}

#[wasm_bindgen]
pub fn tv2s(t: f64, v: f64) -> f64 {
    crate::rust_if97::tv(t, v, OS)
}

#[wasm_bindgen]
pub fn tv2x(t: f64, v: f64) -> f64 {
    crate::rust_if97::tv(t, v, OX)
}

// p,x->t,h,s,x
#[wasm_bindgen]
pub fn px2t(p: f64, x: f64) -> f64 {
    crate::rust_if97::px(p, x, OT)
}

#[wasm_bindgen]
pub fn px2h(p: f64, x: f64) -> f64 {
    crate::rust_if97::px(p, x, OH)
}

#[wasm_bindgen]
pub fn px2s(p: f64, x: f64) -> f64 {
    crate::rust_if97::px(p, x, OS)
}

#[wasm_bindgen]
pub fn px2v(p: f64, x: f64) -> f64 {
    crate::rust_if97::px(p, x, OV)
}

// t,x->p,h,s,v
#[wasm_bindgen]
pub fn tx2p(t: f64, x: f64) -> f64 {
    crate::rust_if97::tx(t, x, OP)
}

#[wasm_bindgen]
pub fn tx2h(t: f64, x: f64) -> f64 {
    crate::rust_if97::tx(t, x, OH)
}

#[wasm_bindgen]
pub fn tx2s(t: f64, x: f64) -> f64 {
    crate::rust_if97::tx(t, x, OS)
}

#[wasm_bindgen]
pub fn tx2v(t: f64, x: f64) -> f64 {
    crate::rust_if97::tx(t, x, OV)
}

// h,x->p,t,s,v
#[wasm_bindgen]
pub fn hx2p(h: f64, x: f64) -> f64 {
    crate::rust_if97::hx(h, x, OP)
}

#[wasm_bindgen]
pub fn hx2t(h: f64, x: f64) -> f64 {
    crate::rust_if97::hx(h, x, OT)
}

#[wasm_bindgen]
pub fn hx2s(h: f64, x: f64) -> f64 {
    crate::rust_if97::hx(h, x, OS)
}

#[wasm_bindgen]
pub fn hx2v(h: f64, x: f64) -> f64 {
    crate::rust_if97::hx(h, x, OV)
}

// s,x->p,t,h,v
#[wasm_bindgen]
pub fn sx2p(s: f64, x: f64) -> f64 {
    crate::rust_if97::sx(s, x, OP)
}

#[wasm_bindgen]
pub fn sx2t(s: f64, x: f64) -> f64 {
    crate::rust_if97::sx(s, x, OT)
}

#[wasm_bindgen]
pub fn sx2h(s: f64, x: f64) -> f64 {
    crate::rust_if97::sx(s, x, OH)
}

#[wasm_bindgen]
pub fn sx2v(s: f64, x: f64) -> f64 {
    crate::rust_if97::sx(s, x, OV)
}
