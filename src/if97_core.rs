//! Core IAPWS-IF97 implementation with i32-based property IDs
//! Used by all language bindings (Python, WASM, C)

use crate::algo::*;
use crate::common::*;
use crate::r1::*;
use crate::r2::*;
use crate::r3::*;
use crate::r4::*;
use crate::r5::*;

#[inline(always)]
pub fn core_pt(p: f64, t: f64, o_id: i32) -> f64 {
    let T: f64 = t + 273.15;
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => p,
        OT => t,
        _ => pair_properties(p, T, o_id, pT_sub_region, pT_reg1, pT_reg2, pT_reg3, pT_reg4, pT_reg5, reg),
    }
}

#[inline(always)]
pub fn core_ph(p: f64, h: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => p,
        OH => h,
        _ => pair_properties(p, h, o_id, ph_sub_region, ph_reg1, ph_reg2, ph_reg3, ph_reg4, ph_reg5, reg),
    }
}

#[inline(always)]
pub fn core_ps(p: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => p,
        OS => s,
        _ => pair_properties(p, s, o_id, ps_sub_region, ps_reg1, ps_reg2, ps_reg3, ps_reg4, ps_reg5, reg),
    }
}

#[inline(always)]
pub fn core_hs(h: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OH => h,
        OS => s,
        _ => pair_properties(h, s, o_id, hs_sub_region, hs_reg1, hs_reg2, hs_reg3, hs_reg4, hs_reg5, reg),
    }
}

#[inline(always)]
pub fn core_pv(p: f64, v: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => p,
        OV => v,
        _ => pair_properties(p, v, o_id, pv_sub_region, pv_reg1, pv_reg2, pv_reg3, pv_reg4, pv_reg5, reg),
    }
}

#[inline(always)]
pub fn core_tv(t: f64, v: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => t,
        OV => v,
        _ => pair_properties(t, v, o_id, tv_sub_region, tv_reg1, tv_reg2, tv_reg3, tv_reg4, tv_reg5, reg),
    }
}

#[inline(always)]
pub fn core_th(t: f64, h: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => t,
        OH => h,
        _ => pair_properties(t, h, o_id, th_sub_region, th_reg1, th_reg2, th_reg3, th_reg4, th_reg5, reg),
    }
}

#[inline(always)]
pub fn core_ts(t: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => t,
        OS => s,
        _ => pair_properties(t, s, o_id, ts_sub_region, ts_reg1, ts_reg2, ts_reg3, ts_reg4, ts_reg5, reg),
    }
}

#[inline(always)]
pub fn core_px(p: f64, x: f64, o_id: i32) -> f64 {
    if p > P_MAX4 || p < P_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OP => p,
        OX => x,
        _ => px_reg4(p, x, o_id),
    }
}

#[inline(always)]
pub fn core_tx(t: f64, x: f64, o_id: i32) -> f64 {
    match o_id {
        OT => t,
        OX => x,
        _ => {
            let T: f64 = t + 273.15;
            if T > T_MAX4 || T < T_MIN4 || x > 1.0 || x < 0.0 {
                return INVALID_VALUE as f64;
            }
            Tx_reg4(T, x, o_id)
        }
    }
}

#[inline(always)]
pub fn core_hx(h: f64, x: f64, o_id: i32) -> f64 {
    if h > H_MAX4 || h < H_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OH => h,
        OX => x,
        _ => hx_reg4(h, x, o_id),
    }
}

#[inline(always)]
pub fn core_sx(s: f64, x: f64, o_id: i32) -> f64 {
    if s > S_MAX4 || s < S_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OS => s,
        OX => x,
        _ => sx_reg4(s, x, o_id),
    }
}

#[inline(always)]
pub fn core_pt2h(p: f64, t: f64) -> f64 { core_pt(p, t, OH) }
#[inline(always)]
pub fn core_pt2s(p: f64, t: f64) -> f64 { core_pt(p, t, OS) }
#[inline(always)]
pub fn core_pt2v(p: f64, t: f64) -> f64 { core_pt(p, t, OV) }
#[inline(always)]
pub fn core_pt2x(p: f64, t: f64) -> f64 { core_pt(p, t, OX) }

#[inline(always)]
pub fn core_ph2t(p: f64, h: f64) -> f64 { core_ph(p, h, OT) }
#[inline(always)]
pub fn core_ph2s(p: f64, h: f64) -> f64 { core_ph(p, h, OS) }
#[inline(always)]
pub fn core_ph2v(p: f64, h: f64) -> f64 { core_ph(p, h, OV) }
#[inline(always)]
pub fn core_ph2x(p: f64, h: f64) -> f64 { core_ph(p, h, OX) }

#[inline(always)]
pub fn core_ps2t(p: f64, s: f64) -> f64 { core_ps(p, s, OT) }
#[inline(always)]
pub fn core_ps2h(p: f64, s: f64) -> f64 { core_ps(p, s, OH) }
#[inline(always)]
pub fn core_ps2v(p: f64, s: f64) -> f64 { core_ps(p, s, OV) }
#[inline(always)]
pub fn core_ps2x(p: f64, s: f64) -> f64 { core_ps(p, s, OX) }

#[inline(always)]
pub fn core_pv2t(p: f64, v: f64) -> f64 { core_pv(p, v, OT) }
#[inline(always)]
pub fn core_pv2h(p: f64, v: f64) -> f64 { core_pv(p, v, OH) }
#[inline(always)]
pub fn core_pv2s(p: f64, v: f64) -> f64 { core_pv(p, v, OS) }
#[inline(always)]
pub fn core_pv2x(p: f64, v: f64) -> f64 { core_pv(p, v, OX) }

#[inline(always)]
pub fn core_hs2p(h: f64, s: f64) -> f64 { core_hs(h, s, OP) }
#[inline(always)]
pub fn core_hs2t(h: f64, s: f64) -> f64 { core_hs(h, s, OT) }
#[inline(always)]
pub fn core_hs2v(h: f64, s: f64) -> f64 { core_hs(h, s, OV) }
#[inline(always)]
pub fn core_hs2x(h: f64, s: f64) -> f64 { core_hs(h, s, OX) }

#[inline(always)]
pub fn core_th2p(t: f64, h: f64) -> f64 { core_th(t, h, OP) }
#[inline(always)]
pub fn core_th2s(t: f64, h: f64) -> f64 { core_th(t, h, OS) }
#[inline(always)]
pub fn core_th2v(t: f64, h: f64) -> f64 { core_th(t, h, OV) }
#[inline(always)]
pub fn core_th2x(t: f64, h: f64) -> f64 { core_th(t, h, OX) }

#[inline(always)]
pub fn core_ts2p(t: f64, s: f64) -> f64 { core_ts(t, s, OP) }
#[inline(always)]
pub fn core_ts2h(t: f64, s: f64) -> f64 { core_ts(t, s, OH) }
#[inline(always)]
pub fn core_ts2v(t: f64, s: f64) -> f64 { core_ts(t, s, OV) }
#[inline(always)]
pub fn core_ts2x(t: f64, s: f64) -> f64 { core_ts(t, s, OX) }

#[inline(always)]
pub fn core_tv2p(t: f64, v: f64) -> f64 { core_tv(t, v, OP) }
#[inline(always)]
pub fn core_tv2h(t: f64, v: f64) -> f64 { core_tv(t, v, OH) }
#[inline(always)]
pub fn core_tv2s(t: f64, v: f64) -> f64 { core_tv(t, v, OS) }
#[inline(always)]
pub fn core_tv2x(t: f64, v: f64) -> f64 { core_tv(t, v, OX) }

#[inline(always)]
pub fn core_px2t(p: f64, x: f64) -> f64 { core_px(p, x, OT) }
#[inline(always)]
pub fn core_px2h(p: f64, x: f64) -> f64 { core_px(p, x, OH) }
#[inline(always)]
pub fn core_px2s(p: f64, x: f64) -> f64 { core_px(p, x, OS) }
#[inline(always)]
pub fn core_px2v(p: f64, x: f64) -> f64 { core_px(p, x, OV) }

#[inline(always)]
pub fn core_tx2p(t: f64, x: f64) -> f64 { core_tx(t, x, OP) }
#[inline(always)]
pub fn core_tx2h(t: f64, x: f64) -> f64 { core_tx(t, x, OH) }
#[inline(always)]
pub fn core_tx2s(t: f64, x: f64) -> f64 { core_tx(t, x, OS) }
#[inline(always)]
pub fn core_tx2v(t: f64, x: f64) -> f64 { core_tx(t, x, OV) }

#[inline(always)]
pub fn core_hx2p(h: f64, x: f64) -> f64 { core_hx(h, x, OP) }
#[inline(always)]
pub fn core_hx2t(h: f64, x: f64) -> f64 { core_hx(h, x, OT) }
#[inline(always)]
pub fn core_hx2s(h: f64, x: f64) -> f64 { core_hx(h, x, OS) }
#[inline(always)]
pub fn core_hx2v(h: f64, x: f64) -> f64 { core_hx(h, x, OV) }

#[inline(always)]
pub fn core_sx2p(s: f64, x: f64) -> f64 { core_sx(s, x, OP) }
#[inline(always)]
pub fn core_sx2t(s: f64, x: f64) -> f64 { core_sx(s, x, OT) }
#[inline(always)]
pub fn core_sx2h(s: f64, x: f64) -> f64 { core_sx(s, x, OH) }
#[inline(always)]
pub fn core_sx2v(s: f64, x: f64) -> f64 { core_sx(s, x, OV) }

/// Isentropic enthalpy drop: h(pi,ti) - h(pe, s=const)
///   ishd(pi,ti,pe)
///
/// # Examples
///
/// ```
/// use seuif97::*;
///
/// let pi:f64 = 16.0;
/// let ti:f64 = 535.1;
/// let pe:f64 = 5.0;
/// let delta_h = ishd(pi, ti, pe);
/// println!("pi={pi} ti={ti} pe={pe} ishd={delta_h:.3}");
/// ```
/// Returns INVALID_VALUE if input is invalid
pub fn core_ishd(pi: f64, ti: f64, pe: f64) -> f64 {
    if pi <= pe {
        return INVALID_VALUE as f64;
    }
    let hi = core_pt(pi, ti, OH);
    if hi <= 0.0 {
        return INVALID_VALUE as f64;
    }
    let si = core_pt(pi, ti, OS);
    if si < 0.0 {
        return INVALID_VALUE as f64;
    }
    let he_isos = core_ps(pe, si, OH);
    if he_isos < 0.0 {
        return INVALID_VALUE as f64;
    }
    hi - he_isos
}

/// ief(pi,ti,pe,te) - Isentropic efficiency (%) for superheated steam expansion
///
/// # Examples
///
/// ```
///  use seuif97::*;
///
/// let pi:f64 = 16.0;
/// let ti:f64 = 535.1;
/// let pe:f64 = 5.0;
/// let te:f64 = 350.0;
/// let eff = ief(pi, ti, pe, te);
/// println!("pi={pi} ti={ti} pe={pe} te={te} ief={eff:.2}%");
/// ```
/// Returns INVALID_VALUE if input is invalid
pub fn core_ief(pi: f64, ti: f64, pe: f64, te: f64) -> f64 {
    if pi <= pe || ti <= te {
        return INVALID_VALUE as f64;
    }
    let hi = core_pt(pi, ti, OH);
    if hi < 0.0 {
        return INVALID_VALUE as f64;
    }
    let si = core_pt(pi, ti, OS);
    if si < 0.0 {
        return INVALID_VALUE as f64;
    }
    let he_isos = core_ps(pe, si, OH);
    if he_isos < 0.0 {
        return INVALID_VALUE as f64;
    }
    let ishd_val = hi - he_isos;

    let he = core_pt(pe, te, OH);
    if he < 0.0 {
        return INVALID_VALUE as f64;
    }
    let se = core_pt(pe, te, OS);
    if se < 0.0 {
        return INVALID_VALUE as f64;
    }
    if (se - si) <= 0.0 {
        return INVALID_VALUE as f64;
    }

    let ahd = hi - he;
    100.0 * ahd / ishd_val
}