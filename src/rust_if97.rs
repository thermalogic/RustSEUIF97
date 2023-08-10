//! The Rust API
//!  
use crate::algo::*;
use crate::common::*;
use crate::r1::*;
use crate::r2::*;
use crate::r3::*;
use crate::r4::*;
use crate::r5::*;

/// the paramters: <br/>
///   `o_id`: the propertry of id;<br/>
///   `region`: the region in the IAPWS-IF97,**optional**
pub struct o_id_region_args {
    o_id: i32,
    region: i32,
}

impl Default for o_id_region_args {
    fn default() -> Self {
        o_id_region_args {
            o_id: 0,
            region: REGION_NONE,
        }
    }
}

impl From<i32> for o_id_region_args {
    fn from(o_id: i32) -> Self {
        Self {
            o_id,
            ..Self::default()
        }
    }
}

impl From<(i32, i32)> for o_id_region_args {
    fn from((o_id, region): (i32, i32)) -> Self {
        Self {
            o_id,
            region,
        }
    }
}

fn o_id_reg_info<R>(o_id_reg: R) -> (i32, i32)
where
    R: Into<o_id_region_args>,
{
    let args = o_id_reg.into();
    let o_id: i32 = args.o_id;
    let reg: i32 = args.region;
    (o_id, reg)
}

/// pt(p,t,o_id) - the propertry of `o_id` (thermodynamic,transport,etc) <br/>
/// pt(p,t,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///
/// # Examples
///
/// ```
/// use if97::*;
///
/// let p:f64 = 3.0;
/// let t:f64= 300.0-273.15;
/// let h=pt(p,t,OH);
/// // set the region
/// let s=pt(p,t,(OS,1));
/// println!("p={p:.6} t={t:.6} h={h:.6} s={s:.6}");    
/// ```
///
pub fn pt<R>(p: f64, t: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    let T: f64 = t + 273.15;
    match o_id {
        OP => return p,
        OT => return t,
        _ => pair_properties(p, T, o_id, pT_sub_region, pT_reg1, pT_reg2, pT_reg3, pT_reg4, pT_reg5, reg),
    }
}

/// ph(p,h,o_id) - the propertry of `o_id` (thermodynamic,transport,etc)<br/>
/// ph(p,h,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///
/// # Examples
///
/// ```
/// use if97::*;
///
/// let p:f64 = 3.0;
/// let h:f64= 0.115331273e+3;
/// let t=ph(p,h,OT);
/// // set the region
/// let s=ph(p,h,(OS,1));
/// println!("p={p:.6} h={h:.6} t={t:.6} s={s:.6}");    
/// ```
///
pub fn ph<R>(p: f64, h: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OP => return p,
        OH => return h,
        _ => pair_properties(p, h, o_id, ph_sub_region, ph_reg1, ph_reg2, ph_reg3, ph_reg4, ph_reg5, reg),
    }
}

/// ps(p,s,o_id) - the propertry of `o_id` (thermodynamic,transport,etc)<br/>
/// ps(p,s,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
/// # Examples
///
///```
/// use if97::*;
///
/// let p:f64= 3.0;
/// let s:f64= 0.392294792;
/// let t=ps(p,s,OT);
/// // set the region
/// let h=ps(p,s,(OH,1));
/// println!("p={p:.6} s={s:.6} t={t:.6} h={h:6}");    
/// ```
pub fn ps<R>(p: f64, s: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OP => return p,
        OS => return s,
        _ => pair_properties(p, s, o_id, ps_sub_region, ps_reg1, ps_reg2, ps_reg3, ps_reg4, ps_reg5, reg),
    }
}

/// hs(h,s,o_id) - the propertry of `o_id` (thermodynamic,transport,etc)<br/>
/// hs(h,s,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///
/// # Examples
///
///```
/// use if97::*;
///
/// let h:f64= 0.115331273e+3;;
/// let s:f64= 0.392294792;
/// let p=hs(h,s,OP);
/// // set the region
/// let t=hs(h,s,(OT,1));
/// println!("h={h:.6} s={s:.6} p={p:.6} t={t:.6}");
/// ```
///   
pub fn hs<R>(h: f64, s: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OH => return h,
        OS => return s,
        _ => pair_properties(h, s, o_id, hs_sub_region, hs_reg1, hs_reg2, hs_reg3, hs_reg4, hs_reg5, reg),
    }
}

///  px(p,x,o_id) - the propertry of `o_id` (thermodynamic)
///
///  # Examples
///
///```
///  use if97::*;
///
///  let p: f64 = 0.1;
///  let x: f64 = 0.0; // x= 0.0 saturation water ,x=1.0 saturation steam
///  let t: f64 = px(p, x, OT);
///  let h: f64 = px(p, x, OH);
///  println!("px: p={p:.6} x={x:.6} t={t:.6} h={h:.6}");
///```

// #[no_mangle]
// pub unsafe extern "C"  fn px(p: f64, x: f64, o_id: i32) -> f64 {
pub fn px(p: f64, x: f64, o_id: i32) -> f64 {
    if p > P_MAX4 || p < P_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OP => return p,
        OX => return x,
        OT => return px_reg4(p, x, o_id) - 273.15,
        _ => return px_reg4(p, x, o_id),
    }
}

///  tx(t,x,o_id) - the propertry of `o_id` (thermodynamic)
///
///  # Examples
///
///```
///  use if97::*;
///
///  let t: f64 =372.755919-273.15;
///  let x: f64 = 0.0; // x= 0.0 saturation water ,x=1.0 saturation steam
///  let p: f64 = tx(t, x, OP);
///  let h: f64 = tx(t, x, OH);
///  println!("tx: p={p:.6} x={x:.6} t={t:.6} h={h:.6}");
/// ```
///
pub fn tx(t: f64, x: f64, o_id: i32) -> f64 {
    if (o_id==OT)
    {
        return t;
    } 
    if (o_id==OX)
    {
        return x;
    } 
    let T: f64 = t + 273.15;
    if T > T_MAX4 || T < T_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    Tx_reg4(T, x, o_id)
}

//  Functions of the  extended input pairs (p,v),(t,v),(t,s),(t,h)

/// pv(p,v,o_id) - the propertry of `o_id`(thermodynamic,transport,etc)<br/>
/// pv(p,v,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///
/// # Examples
///
///```
/// use if97::*;
///
/// let p:f64= 3.0;
/// let v:f64= 0.100215168e-2;
/// let t=pv(p,v,OT);
/// //set the region
/// let h=pv(p,v,(OH,1));
/// println!("p={p:.6} v={v:.6} t={t:.6}");    
/// ```
pub fn pv<R>(p: f64, v: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OP => return p,
        OV => return v,
        _ => pair_properties(p, v, o_id, pv_sub_region, pv_reg1, pv_reg2, pv_reg3, pv_reg4, pv_reg5, reg),
    }
}

/// tv(t,v,o_id) - the propertry of `o_id`(thermodynamic,transport,etc)<br/>
/// tv(t,v,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///
/// # Examples
///
///```
/// use if97::*;
///
/// let t:f64=300.0-273.15;
/// let v:f64= 0.100215168e-2;
/// let p=tv(t,v,OP);
/// //set the regiion
/// let s=tv(t,v,(OS,1));
/// println!("t={p:.6} v={v:.6} p={p:.6} s={s:.6}");    
/// ```
pub fn tv<R>(t: f64, v: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OT => return t,
        OV => return v,
        _ => pair_properties(t, v, o_id, tv_sub_region, tv_reg1, tv_reg2, tv_reg3, tv_reg4, tv_reg5, reg),
    }
}

/// th(t,h,o_id) - the propertry of `o_id`(thermodynamic,transport,etc)<br/>
/// th(t,h,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///  
/// # Examples
///
///```
/// use if97::*;
///
/// let t:f64=300.0-273.15;
/// let h:f64=0.115331273e+3;
/// let p=th(t,h,OP);
/// // set the region
/// let s=th(t,h,(OS,1));
/// println!("t={p:.6} h={h:.6} p={p:.6} s={s:.6}");    
/// ```
pub fn th<R>(t: f64, h: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OT => return t,
        OH => return h,
        _ => pair_properties(t, h, o_id, th_sub_region, th_reg1, th_reg2, th_reg3, th_reg4, th_reg5, reg),
    }
}

/// ts(t,s,o_id) - the propertry of `o_id`(thermodynamic,transport,etc)<br/>
/// ts(t,s,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///   
/// # Examples
///
///```
/// use if97::*;
///
/// let t:f64=300.0-273.15;
/// let s:f64=0.392294792;
/// let p=ts(t,s,OP);
/// // set the region
/// let h=ts(t,s,(OH,1));
/// println!("t={p:.6} s={s:.6} p={p:.6} h={h:.6}");    
/// ```
pub fn ts<R>(t: f64, s: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OT => return t,
        OS => return s,
        _ => pair_properties(t, s, o_id, ts_sub_region, ts_reg1, ts_reg2, ts_reg3, ts_reg4, ts_reg5, reg),
    }
}

/// hx(h,x,o_id) - the propertry of `o_id`(thermodynamic)
///
///  # Examples
///
///```
///  use if97::*;
///
///  let h: f64 =1094.690434;
///  let x: f64 = 0.3;
///  let p: f64 = hx(h, x, OP);
///  let t: f64 = hx(h, x, OT);
///  println!("hx: h={p:.6} x={x:.6} p={p:.6} t={t:.6}");
/// ```
///
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

/// sx(s,x,o_id) - the propertry of `o_id`(thermodynamic)
///
///  # Examples
///
///```
///  use if97::*;
///
///  let s: f64 =3.119434;
///  let x: f64 = 0.3;
///  let p: f64 = sx(s, x, OP);
///  let t: f64 = sx(s, x, OT);
///  println!("sx: s={p:.6} x={x:.6} p={p:.6} t={t:.6}");
/// ```
///
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
