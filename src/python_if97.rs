//!  The Python API
#[pyo3::pymodule]
mod seuif97 {
    use pyo3::prelude::*;

    use crate::if97_core::*;

    #[pyfunction]
    fn pt(p: f64, t: f64, o_id: i32) -> f64 {
        core_pt(p, t, o_id)
    }

    #[pyfunction]
    fn ph(p: f64, h: f64, o_id: i32) -> f64 {
        core_ph(p, h, o_id)
    }

    #[pyfunction]
    fn ps(p: f64, s: f64, o_id: i32) -> f64 {
        core_ps(p, s, o_id)
    }

    #[pyfunction]
    fn hs(h: f64, s: f64, o_id: i32) -> f64 {
        core_hs(h, s, o_id)
    }

    #[pyfunction]
    fn px(p: f64, x: f64, o_id: i32) -> f64 {
        core_px(p, x, o_id)
    }

    #[pyfunction]
    fn tx(t: f64, x: f64, o_id: i32) -> f64 {
        core_tx(t, x, o_id)
    }

    #[pyfunction]
    fn pv(p: f64, v: f64, o_id: i32) -> f64 {
        core_pv(p, v, o_id)
    }

    #[pyfunction]
    fn tv(t: f64, v: f64, o_id: i32) -> f64 {
        core_tv(t, v, o_id)
    }

    #[pyfunction]
    fn th(t: f64, h: f64, o_id: i32) -> f64 {
        core_th(t, h, o_id)
    }

    #[pyfunction]
    fn ts(t: f64, s: f64, o_id: i32) -> f64 {
        core_ts(t, s, o_id)
    }

    #[pyfunction]
    fn hx(h: f64, x: f64, o_id: i32) -> f64 {
        core_hx(h, x, o_id)
    }

    #[pyfunction]
    fn sx(s: f64, x: f64, o_id: i32) -> f64 {
        core_sx(s, x, o_id)
    }

    #[pyfunction]
    fn pt2h(p: f64, t: f64) -> f64 {
        core_pt2h(p, t)
    }

    #[pyfunction]
    fn pt2s(p: f64, t: f64) -> f64 {
        core_pt2s(p, t)
    }

    #[pyfunction]
    fn pt2v(p: f64, t: f64) -> f64 {
        core_pt2v(p, t)
    }

    #[pyfunction]
    fn pt2x(p: f64, t: f64) -> f64 {
        core_pt2x(p, t)
    }

    #[pyfunction]
    fn ph2t(p: f64, h: f64) -> f64 {
        core_ph2t(p, h)
    }

    #[pyfunction]
    fn ph2s(p: f64, h: f64) -> f64 {
        core_ph2s(p, h)
    }

    #[pyfunction]
    fn ph2v(p: f64, h: f64) -> f64 {
        core_ph2v(p, h)
    }

    #[pyfunction]
    fn ph2x(p: f64, h: f64) -> f64 {
        core_ph2x(p, h)
    }

    #[pyfunction]
    fn ps2t(p: f64, s: f64) -> f64 {
        core_ps2t(p, s)
    }

    #[pyfunction]
    fn ps2h(p: f64, s: f64) -> f64 {
        core_ps2h(p, s)
    }

    #[pyfunction]
    fn ps2v(p: f64, s: f64) -> f64 {
        core_ps2v(p, s)
    }

    #[pyfunction]
    fn ps2x(p: f64, s: f64) -> f64 {
        core_ps2x(p, s)
    }

    #[pyfunction]
    fn pv2t(p: f64, v: f64) -> f64 {
        core_pv2t(p, v)
    }

    #[pyfunction]
    fn pv2h(p: f64, v: f64) -> f64 {
        core_pv2h(p, v)
    }

    #[pyfunction]
    fn pv2s(p: f64, v: f64) -> f64 {
        core_pv2s(p, v)
    }

    #[pyfunction]
    fn pv2x(p: f64, v: f64) -> f64 {
        core_pv2x(p, v)
    }

    #[pyfunction]
    fn hs2p(h: f64, s: f64) -> f64 {
        core_hs2p(h, s)
    }

    #[pyfunction]
    fn hs2t(h: f64, s: f64) -> f64 {
        core_hs2t(h, s)
    }

    #[pyfunction]
    fn hs2v(h: f64, s: f64) -> f64 {
        core_hs2v(h, s)
    }

    #[pyfunction]
    fn hs2x(h: f64, s: f64) -> f64 {
        core_hs2x(h, s)
    }

    #[pyfunction]
    fn th2p(t: f64, h: f64) -> f64 {
        core_th2p(t, h)
    }

    #[pyfunction]
    fn th2s(t: f64, h: f64) -> f64 {
        core_th2s(t, h)
    }

    #[pyfunction]
    fn th2v(t: f64, h: f64) -> f64 {
        core_th2v(t, h)
    }

    #[pyfunction]
    fn th2x(t: f64, h: f64) -> f64 {
        core_th2x(t, h)
    }

    #[pyfunction]
    fn ts2p(t: f64, s: f64) -> f64 {
        core_ts2p(t, s)
    }

    #[pyfunction]
    fn ts2h(t: f64, s: f64) -> f64 {
        core_ts2h(t, s)
    }

    #[pyfunction]
    fn ts2v(t: f64, s: f64) -> f64 {
        core_ts2v(t, s)
    }

    #[pyfunction]
    fn ts2x(t: f64, s: f64) -> f64 {
        core_ts2x(t, s)
    }

    #[pyfunction]
    fn tv2p(t: f64, v: f64) -> f64 {
        core_tv2p(t, v)
    }

    #[pyfunction]
    fn tv2h(t: f64, v: f64) -> f64 {
        core_tv2h(t, v)
    }

    #[pyfunction]
    fn tv2s(t: f64, v: f64) -> f64 {
        core_tv2s(t, v)
    }

    #[pyfunction]
    fn tv2x(t: f64, v: f64) -> f64 {
        core_tv2x(t, v)
    }

    #[pyfunction]
    fn px2t(p: f64, x: f64) -> f64 {
        core_px2t(p, x)
    }

    #[pyfunction]
    fn px2h(p: f64, x: f64) -> f64 {
        core_px2h(p, x)
    }

    #[pyfunction]
    fn px2s(p: f64, x: f64) -> f64 {
        core_px2s(p, x)
    }

    #[pyfunction]
    fn px2v(p: f64, x: f64) -> f64 {
        core_px2v(p, x)
    }

    #[pyfunction]
    fn tx2p(t: f64, x: f64) -> f64 {
        core_tx2p(t, x)
    }

    #[pyfunction]
    fn tx2h(t: f64, x: f64) -> f64 {
        core_tx2h(t, x)
    }

    #[pyfunction]
    fn tx2s(t: f64, x: f64) -> f64 {
        core_tx2s(t, x)
    }

    #[pyfunction]
    fn tx2v(t: f64, x: f64) -> f64 {
        core_tx2v(t, x)
    }

    #[pyfunction]
    fn hx2p(h: f64, x: f64) -> f64 {
        core_hx2p(h, x)
    }

    #[pyfunction]
    fn hx2t(h: f64, x: f64) -> f64 {
        core_hx2t(h, x)
    }

    #[pyfunction]
    fn hx2s(h: f64, x: f64) -> f64 {
        core_hx2s(h, x)
    }

    #[pyfunction]
    fn hx2v(h: f64, x: f64) -> f64 {
        core_hx2v(h, x)
    }

    #[pyfunction]
    fn sx2p(s: f64, x: f64) -> f64 {
        core_sx2p(s, x)
    }

    #[pyfunction]
    fn sx2t(s: f64, x: f64) -> f64 {
        core_sx2t(s, x)
    }

    #[pyfunction]
    fn sx2h(s: f64, x: f64) -> f64 {
        core_sx2h(s, x)
    }

    #[pyfunction]
    fn sx2v(s: f64, x: f64) -> f64 {
        core_sx2v(s, x)
    }

    #[pyfunction]
    fn ishd(pi: f64, ti: f64, pe: f64) -> f64 {
        core_ishd(pi, ti, pe)
    }

    #[pyfunction]
    fn ief(pi: f64, ti: f64, pe: f64, te: f64) -> f64 {
        core_ief(pi, ti, pe, te)
    }
}