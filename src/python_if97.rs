//!  The Python API
#[pyo3::pymodule]
mod seuif97 {
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

    // the convenience functions (direct output)
    //p,t -> h,s,v,x
   #[pyfunction]
    fn pt2h(p: f64, t: f64) -> f64 {
         crate::rust_if97::pt(p,t,OH)
    }

   #[pyfunction]
    fn pt2s(p: f64, t: f64) -> f64 {
        crate::rust_if97::pt(p,t,OS)
    }

   #[pyfunction]
    fn pt2v(p: f64, t: f64) -> f64 {
        crate::rust_if97::pt(p,t,OV)
    }

    #[pyfunction]
    fn pt2x(p: f64, t: f64) -> f64 {
        crate::rust_if97::pt(p,t,OX)
    }
    // p,h->t,s,v,x
     #[pyfunction]
    #[pyfunction]
    fn ph2t(p: f64, h: f64) -> f64 {
        crate::rust_if97::ph(p,h,OT)    
    }
    
    #[pyfunction]
    fn ph2s(p: f64, h: f64) -> f64 {
        crate::rust_if97::ph(p,h,OS)
    }

    #[pyfunction]
    fn ph2v(p: f64, h: f64) -> f64 {
        crate::rust_if97::ph(p,h,OV)
    }

   #[pyfunction]
   fn ph2x(p: f64, h: f64) -> f64 {
        crate::rust_if97::ph(p,h,OX)
    }

    // p,s->t,h,v,x
    #[pyfunction]
    fn ps2t(p: f64, s: f64) -> f64 {
        crate::rust_if97::ps(p,s,OT)
    }
  
    #[pyfunction]
    fn ps2h(p: f64, s: f64) -> f64 {
        crate::rust_if97::ps(p,s,OH)
    }
  
    #[pyfunction]
    fn ps2v(p: f64, s: f64) -> f64 {
        crate::rust_if97::ps(p,s,OV)
    }
  
    #[pyfunction]
    fn ps2x(p: f64, s: f64) -> f64 {
        crate::rust_if97::ps(p,s,OX)
    }

    //p,v  -> t,h,s,x
    #[pyfunction]
    fn pv2t(p: f64, v: f64) -> f64 {
        crate::rust_if97::pv(p,v,OT)
    }
   
    #[pyfunction]
    fn pv2h(p: f64, v: f64) -> f64 {
        crate::rust_if97::pv(p,v,OH)
    }
   
    #[pyfunction]
    fn pv2s(p: f64, v: f64) -> f64 {
        crate::rust_if97::pv(p,v,OS)
    }
  
    #[pyfunction]
    fn pv2x(p: f64, v: f64) -> f64 {
        crate::rust_if97::pv(p,v,OX)
    }
    
    // h,s->p,t,v,x
    #[pyfunction]
    fn hs2p(h: f64, s: f64) -> f64 {
        crate::rust_if97::hs(h,s,OP)
    }
  
    #[pyfunction]
    fn hs2t(h: f64, s: f64) -> f64 {
        crate::rust_if97::hs(h,s,OT)
    }
  
    #[pyfunction]
    fn hs2v(h: f64, s: f64) -> f64 {
        crate::rust_if97::hs(h,s,OV)
    }
  
    #[pyfunction]
    fn hs2x(h: f64, s: f64) -> f64 {
        crate::rust_if97::hs(h,s,OX)
    }
       
    // t,h -> p,s,v,x
    #[pyfunction]
    fn th2p(t: f64, h: f64) -> f64 {
        crate::rust_if97::th(t,h,OP)
    }

    #[pyfunction]
    fn th2s(t: f64, h: f64) -> f64 {
        crate::rust_if97::th(t,h,OS)
    }

   #[pyfunction]
    fn th2v(t: f64, h: f64) -> f64 {
        crate::rust_if97::th(t,h,OV)
    }

    #[pyfunction]
    fn th2x(t: f64, h: f64) -> f64 {
        crate::rust_if97::th(t,h,OX)
    }

    // t,s -> p,h,v,x
    #[pyfunction]
    fn ts2p(t: f64, s: f64) -> f64 {
        crate::rust_if97::ts(t,s,OP)
    }

    #[pyfunction]
    fn ts2h(t: f64, s: f64) -> f64 {
        crate::rust_if97::ts(t,s,OH)
    }

    #[pyfunction]
    fn ts2v(t: f64, s: f64) -> f64 {
        crate::rust_if97::ts(t,s,OV)
    }

    #[pyfunction]
    fn ts2x(t: f64, s: f64) -> f64 {
        crate::rust_if97::ts(t,s,OX)
    }

    // t,v -> p,h,s,x
    #[pyfunction]
    fn tv2p(t: f64, v: f64) -> f64 {
        crate::rust_if97::tv(t,v,OP)
    }
    
    #[pyfunction]
    fn tv2h(t: f64, v: f64) -> f64 {
        crate::rust_if97::tv(t,v,OH)
    }

    #[pyfunction]
    fn tv2s(t: f64, v: f64) -> f64 {
        crate::rust_if97::tv(t,v,OS)
    }

    #[pyfunction]
    fn tv2x(t: f64, v: f64) -> f64 {
        crate::rust_if97::tv(t,v,OX)
    }

    // p,x-> t,h,s,v
    #[pyfunction]
    fn px2t(p: f64, x: f64) -> f64 {
        crate::rust_if97::px(p,x,OT)
    }

    #[pyfunction]
    fn px2h(p: f64, x: f64) -> f64 {
        crate::rust_if97::px(p,x,OH)
    }

    #[pyfunction]
    fn px2s(p: f64, x: f64) -> f64 {
        crate::rust_if97::px(p,x,OS)
    }

    #[pyfunction]
    fn px2v(p: f64, x: f64) -> f64 {
        crate::rust_if97::px(p,x,OV)
    }

    // t,x -> p,h,s,v
     #[pyfunction]
    fn tx2p(t: f64, x: f64) -> f64 {
        crate::rust_if97::tx(t,x,OP)
    }

    #[pyfunction]
    fn tx2h(t: f64, x: f64) -> f64 {
        crate::rust_if97::tx(t,x,OH)
    }
   
    #[pyfunction]
    fn tx2s(t: f64, x: f64) -> f64 {
        crate::rust_if97::tx(t,x,OS)
    }

    #[pyfunction]
    fn tx2v(t: f64, x: f64) -> f64 {
        crate::rust_if97::tx(t,x,OV)
    }

    //h,x->p,t,s,v
    #[pyfunction]
    fn hx2p(h: f64, x: f64) -> f64 {
        crate::rust_if97::hx(h,x,OP)
    }
  
    #[pyfunction]
    fn hx2t(h: f64, x: f64) -> f64 {
        crate::rust_if97::hx(h,x,OT)
    }

    #[pyfunction]
    fn hx2s(h: f64, x: f64) -> f64 {
        crate::rust_if97::hx(h,x,OS)
    }

    #[pyfunction]
    fn hx2v(h: f64, x: f64) -> f64 {
        crate::rust_if97::hx(h,x,OV)
    }

    // s,x->p,t,h,v
    #[pyfunction]
    fn sx2p(s: f64, x: f64) -> f64 {    
       crate::rust_if97::sx(s,x,OP)
    }
    
    #[pyfunction]
    fn sx2t(s: f64, x: f64) -> f64 {    
       crate::rust_if97::sx(s,x,OT)
    }
   
    #[pyfunction]
    fn sx2h(s: f64, x: f64) -> f64 {    
       crate::rust_if97::sx(s,x,OH)
    }
    #[pyfunction]
    fn sx2v(s: f64, x: f64) -> f64 {    
       crate::rust_if97::sx(s,x,OV)
    }

}
    /*
//  py03 : 0.24
#[pymodule]
fn seuif97(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pt, m)?)?;
    m.add_function(wrap_pyfunction!(ph, m)?)?;
    m.add_function(wrap_pyfunction!(ps, m)?)?;
    m.add_function(wrap_pyfunction!(pv, m)?)?;//

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
*/

