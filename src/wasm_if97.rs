use wasm_bindgen::prelude::*;

use crate::if97_core::*;

#[wasm_bindgen]
pub fn pt(p: f64, t: f64, o_id: i32) -> f64 {
    core_pt(p, t, o_id)
}

#[wasm_bindgen]
pub fn ph(p: f64, h: f64, o_id: i32) -> f64 {
    core_ph(p, h, o_id)
}

#[wasm_bindgen]
pub fn ps(p: f64, s: f64, o_id: i32) -> f64 {
    core_ps(p, s, o_id)
}

#[wasm_bindgen]
pub fn hs(h: f64, s: f64, o_id: i32) -> f64 {
    core_hs(h, s, o_id)
}

#[wasm_bindgen]
pub fn px(p: f64, x: f64, o_id: i32) -> f64 {
    core_px(p, x, o_id)
}

#[wasm_bindgen]
pub fn tx(t: f64, x: f64, o_id: i32) -> f64 {
    core_tx(t, x, o_id)
}

#[wasm_bindgen]
pub fn pv(p: f64, v: f64, o_id: i32) -> f64 {
    core_pv(p, v, o_id)
}

#[wasm_bindgen]
pub fn tv(t: f64, v: f64, o_id: i32) -> f64 {
    core_tv(t, v, o_id)
}

#[wasm_bindgen]
pub fn th(t: f64, h: f64, o_id: i32) -> f64 {
    core_th(t, h, o_id)
}

#[wasm_bindgen]
pub fn ts(t: f64, s: f64, o_id: i32) -> f64 {
    core_ts(t, s, o_id)
}

#[wasm_bindgen]
pub fn hx(h: f64, x: f64, o_id: i32) -> f64 {
    core_hx(h, x, o_id)
}

#[wasm_bindgen]
pub fn sx(s: f64, x: f64, o_id: i32) -> f64 {
    core_sx(s, x, o_id)
}

#[wasm_bindgen]
pub fn pt2h(p: f64, t: f64) -> f64 {
    core_pt2h(p, t)
}

#[wasm_bindgen]
pub fn pt2s(p: f64, t: f64) -> f64 {
    core_pt2s(p, t)
}

#[wasm_bindgen]
pub fn pt2v(p: f64, t: f64) -> f64 {
    core_pt2v(p, t)
}

#[wasm_bindgen]
pub fn pt2x(p: f64, t: f64) -> f64 {
    core_pt2x(p, t)
}

#[wasm_bindgen]
pub fn ph2t(p: f64, h: f64) -> f64 {
    core_ph2t(p, h)
}

#[wasm_bindgen]
pub fn ph2s(p: f64, h: f64) -> f64 {
    core_ph2s(p, h)
}

#[wasm_bindgen]
pub fn ph2v(p: f64, h: f64) -> f64 {
    core_ph2v(p, h)
}

#[wasm_bindgen]
pub fn ph2x(p: f64, h: f64) -> f64 {
    core_ph2x(p, h)
}

#[wasm_bindgen]
pub fn ps2t(p: f64, s: f64) -> f64 {
    core_ps2t(p, s)
}

#[wasm_bindgen]
pub fn ps2h(p: f64, s: f64) -> f64 {
    core_ps2h(p, s)
}

#[wasm_bindgen]
pub fn ps2v(p: f64, s: f64) -> f64 {
    core_ps2v(p, s)
}

#[wasm_bindgen]
pub fn ps2x(p: f64, s: f64) -> f64 {
    core_ps2x(p, s)
}

#[wasm_bindgen]
pub fn pv2t(p: f64, v: f64) -> f64 {
    core_pv2t(p, v)
}

#[wasm_bindgen]
pub fn pv2h(p: f64, v: f64) -> f64 {
    core_pv2h(p, v)
}

#[wasm_bindgen]
pub fn pv2s(p: f64, v: f64) -> f64 {
    core_pv2s(p, v)
}

#[wasm_bindgen]
pub fn pv2x(p: f64, v: f64) -> f64 {
    core_pv2x(p, v)
}

#[wasm_bindgen]
pub fn hs2p(h: f64, s: f64) -> f64 {
    core_hs2p(h, s)
}

#[wasm_bindgen]
pub fn hs2t(h: f64, s: f64) -> f64 {
    core_hs2t(h, s)
}

#[wasm_bindgen]
pub fn hs2v(h: f64, s: f64) -> f64 {
    core_hs2v(h, s)
}

#[wasm_bindgen]
pub fn hs2x(h: f64, s: f64) -> f64 {
    core_hs2x(h, s)
}

#[wasm_bindgen]
pub fn th2p(t: f64, h: f64) -> f64 {
    core_th2p(t, h)
}

#[wasm_bindgen]
pub fn th2s(t: f64, h: f64) -> f64 {
    core_th2s(t, h)
}

#[wasm_bindgen]
pub fn th2v(t: f64, h: f64) -> f64 {
    core_th2v(t, h)
}

#[wasm_bindgen]
pub fn th2x(t: f64, h: f64) -> f64 {
    core_th2x(t, h)
}

#[wasm_bindgen]
pub fn ts2p(t: f64, s: f64) -> f64 {
    core_ts2p(t, s)
}

#[wasm_bindgen]
pub fn ts2h(t: f64, s: f64) -> f64 {
    core_ts2h(t, s)
}

#[wasm_bindgen]
pub fn ts2v(t: f64, s: f64) -> f64 {
    core_ts2v(t, s)
}

#[wasm_bindgen]
pub fn ts2x(t: f64, s: f64) -> f64 {
    core_ts2x(t, s)
}

#[wasm_bindgen]
pub fn tv2p(t: f64, v: f64) -> f64 {
    core_tv2p(t, v)
}

#[wasm_bindgen]
pub fn tv2h(t: f64, v: f64) -> f64 {
    core_tv2h(t, v)
}

#[wasm_bindgen]
pub fn tv2s(t: f64, v: f64) -> f64 {
    core_tv2s(t, v)
}

#[wasm_bindgen]
pub fn tv2x(t: f64, v: f64) -> f64 {
    core_tv2x(t, v)
}

#[wasm_bindgen]
pub fn px2t(p: f64, x: f64) -> f64 {
    core_px2t(p, x)
}

#[wasm_bindgen]
pub fn px2h(p: f64, x: f64) -> f64 {
    core_px2h(p, x)
}

#[wasm_bindgen]
pub fn px2s(p: f64, x: f64) -> f64 {
    core_px2s(p, x)
}

#[wasm_bindgen]
pub fn px2v(p: f64, x: f64) -> f64 {
    core_px2v(p, x)
}

#[wasm_bindgen]
pub fn tx2p(t: f64, x: f64) -> f64 {
    core_tx2p(t, x)
}

#[wasm_bindgen]
pub fn tx2h(t: f64, x: f64) -> f64 {
    core_tx2h(t, x)
}

#[wasm_bindgen]
pub fn tx2s(t: f64, x: f64) -> f64 {
    core_tx2s(t, x)
}

#[wasm_bindgen]
pub fn tx2v(t: f64, x: f64) -> f64 {
    core_tx2v(t, x)
}

#[wasm_bindgen]
pub fn hx2p(h: f64, x: f64) -> f64 {
    core_hx2p(h, x)
}

#[wasm_bindgen]
pub fn hx2t(h: f64, x: f64) -> f64 {
    core_hx2t(h, x)
}

#[wasm_bindgen]
pub fn hx2s(h: f64, x: f64) -> f64 {
    core_hx2s(h, x)
}

#[wasm_bindgen]
pub fn hx2v(h: f64, x: f64) -> f64 {
    core_hx2v(h, x)
}

#[wasm_bindgen]
pub fn sx2p(s: f64, x: f64) -> f64 {
    core_sx2p(s, x)
}

#[wasm_bindgen]
pub fn sx2t(s: f64, x: f64) -> f64 {
    core_sx2t(s, x)
}

#[wasm_bindgen]
pub fn sx2h(s: f64, x: f64) -> f64 {
    core_sx2h(s, x)
}

#[wasm_bindgen]
pub fn sx2v(s: f64, x: f64) -> f64 {
    core_sx2v(s, x)
}

#[wasm_bindgen]
pub fn ishd(pi: f64, ti: f64, pe: f64) -> f64 {
    core_ishd(pi, ti, pe)
}

#[wasm_bindgen]
pub fn ief(pi: f64, ti: f64, pe: f64, te: f64) -> f64 {
    core_ief(pi, ti, pe, te)
}