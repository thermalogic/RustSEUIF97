//!  The Python API

use pyo3::prelude::*;

use crate::algo::*;
use crate::common::*;
use crate::r1::*;
use crate::r2::*;
use crate::r3::*;
use crate::r4::*;
use crate::r5::*;

#[pyfunction]
fn pt(p: f64, t: f64, o_id: i32) -> f64 {
    let T: f64 = t + 273.15;
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OT => return t,
        _ => pair_properties(p, T, o_id, pT_sub_region, pT_reg1, pT_reg2, pT_reg3, pT_reg4, pT_reg5, reg),
    }
}

#[pyfunction]
fn ph(p: f64, h: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OH => return h,
        _ => pair_properties(p, h, o_id, ph_sub_region, ph_reg1, ph_reg2, ph_reg3, ph_reg4, ph_reg5, reg),
    }
}

#[pyfunction]
fn ps(p: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OS => return s,
        _ => pair_properties(p, s, o_id, ps_sub_region, ps_reg1, ps_reg2, ps_reg3, ps_reg4, ps_reg5, reg),
    }
}

#[pyfunction]
fn hs(h: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OH => return h,
        OS => return s,
        _ => pair_properties(h, s, o_id, hs_sub_region, hs_reg1, hs_reg2, hs_reg3, hs_reg4, hs_reg5, reg),
    }
}

#[pyfunction]
fn px(p: f64, x: f64, o_id: i32) -> f64 {
    if p > P_MAX4 || p < P_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OP => return p,
        OX => return x,
        _ => return px_reg4(p, x, o_id),
    }
}

#[pyfunction]
fn tx(t: f64, x: f64, o_id: i32) -> f64 {
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

#[pyfunction]
fn pv(p: f64, v: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OP => return p,
        OV => return v,
        _ => pair_properties(p, v, o_id, pv_sub_region, pv_reg1, pv_reg2, pv_reg3, pv_reg4, pv_reg5, reg),
    }
}

#[pyfunction]
fn tv(t: f64, v: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => return t,
        OV => return v,
        _ => pair_properties(t, v, o_id, tv_sub_region, tv_reg1, tv_reg2, tv_reg3, tv_reg4, tv_reg5, reg),
    }
}

#[pyfunction]
fn th(t: f64, h: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => return t,
        OH => return h,
        _ => pair_properties(t, h, o_id, th_sub_region, th_reg1, th_reg2, th_reg3, th_reg4, th_reg5, reg),
    }
}

#[pyfunction]
fn ts(t: f64, s: f64, o_id: i32) -> f64 {
    let reg: i32 = REGION_NONE;
    match o_id {
        OT => return t,
        OS => return s,
        _ => pair_properties(t, s, o_id, ts_sub_region, ts_reg1, ts_reg2, ts_reg3, ts_reg4, ts_reg5, reg),
    }
}

#[pyfunction]
fn hx(h: f64, x: f64, o_id: i32) -> f64 {
    if h > H_MAX4 || h < H_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OH => return h,
        OX => return x,
        _ => hx_reg4(h, x, o_id),
    }
}

#[pyfunction]
fn sx(s: f64, x: f64, o_id: i32) -> f64 {
    if s > S_MAX4 || s < S_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OS => return s,
        OX => return x,
        _ => sx_reg4(s, x, o_id),
    }
}

#[pymodule]
fn seuif97(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pt, m)?)?;
    m.add_function(wrap_pyfunction!(ph, m)?)?;
    m.add_function(wrap_pyfunction!(ps, m)?)?;
    m.add_function(wrap_pyfunction!(pv, m)?)?;

    m.add_function(wrap_pyfunction!(tv, m)?)?;
    m.add_function(wrap_pyfunction!(th, m)?)?;
    m.add_function(wrap_pyfunction!(ts, m)?)?;

    m.add_function(wrap_pyfunction!(hs, m)?)?;

    m.add_function(wrap_pyfunction!(px, m)?)?;
    m.add_function(wrap_pyfunction!(tx, m)?)?;
    m.add_function(wrap_pyfunction!(hx, m)?)?;
    m.add_function(wrap_pyfunction!(sx, m)?)?;
    Ok(())
}
