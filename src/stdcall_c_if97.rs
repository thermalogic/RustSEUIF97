//!  The C API: stdcall - Win32 API functions

/// double pt(double p,double t,short o_id)- the property of `o_id` (thermodynamic,transport,etc)
#[no_mangle]
pub unsafe extern "stdcall" fn pt(p: f64, t: f64, o_id: i32) -> f64 {
    crate::if97_core::core_pt(p, t, o_id)
}

/// double ph(double p,double h,short o_id)- the property of `o_id` (thermodynamic,transport,etc)
#[no_mangle]
pub unsafe extern "stdcall" fn ph(p: f64, h: f64, o_id: i32) -> f64 {
    crate::if97_core::core_ph(p, h, o_id)
}

/// double ps(double p,double s,short o_id)- the property of `o_id` (thermodynamic,transport,etc)
#[no_mangle]
pub unsafe extern "stdcall" fn ps(p: f64, s: f64, o_id: i32) -> f64 {
    crate::if97_core::core_ps(p, s, o_id)
}

/// double hs(double h,double s,short o_id)- the property of `o_id` (thermodynamic,transport,etc)
#[no_mangle]
pub unsafe extern "stdcall" fn hs(h: f64, s: f64, o_id: i32) -> f64 {
    crate::if97_core::core_hs(h, s, o_id)
}

/// double px(double p,double x,short o_id) - the property of `o_id` (thermodynamic)
#[no_mangle]
pub unsafe extern "stdcall" fn px(p: f64, x: f64, o_id: i32) -> f64 {
    crate::if97_core::core_px(p, x, o_id)
}

/// double tx(double t,double x,short o_id) - the property of `o_id` (thermodynamic)
#[no_mangle]
pub unsafe extern "stdcall" fn tx(t: f64, x: f64, o_id: i32) -> f64 {
    crate::if97_core::core_tx(t, x, o_id)
}

/// double pv(double p,double v,short o_id)- the property of `o_id` (thermodynamic,transport,etc)
#[no_mangle]
pub unsafe extern "stdcall" fn pv(p: f64, v: f64, o_id: i32) -> f64 {
    crate::if97_core::core_pv(p, v, o_id)
}

/// double tv(double t,double v,short o_id)- the property of `o_id` (thermodynamic,transport,etc)
#[no_mangle]
pub unsafe extern "stdcall" fn tv(t: f64, v: f64, o_id: i32) -> f64 {
    crate::if97_core::core_tv(t, v, o_id)
}

/// double th(double t,double h,short o_id)- the property of `o_id` (thermodynamic,transport,etc)
#[no_mangle]
pub unsafe extern "stdcall" fn th(t: f64, h: f64, o_id: i32) -> f64 {
    crate::if97_core::core_th(t, h, o_id)
}

/// double ts(double t,double s,short o_id)- the property of `o_id` (thermodynamic,transport,etc)
#[no_mangle]
pub unsafe extern "stdcall" fn ts(t: f64, s: f64, o_id: i32) -> f64 {
    crate::if97_core::core_ts(t, s, o_id)
}

/// double hx(double h,double x,short o_id)- the property of `o_id` (thermodynamic)
#[no_mangle]
pub unsafe extern "stdcall" fn hx(h: f64, x: f64, o_id: i32) -> f64 {
    crate::if97_core::core_hx(h, x, o_id)
}

/// double sx(double s,double x,short o_id)- the property of `o_id` (thermodynamic)
#[no_mangle]
pub unsafe extern "stdcall" fn sx(s: f64, x: f64, o_id: i32) -> f64 {
    crate::if97_core::core_sx(s, x, o_id)
}

/// double pt2h(double p,double t) - Calculate specific enthalpy from pressure and temperature
#[no_mangle]
pub unsafe extern "stdcall" fn pt2h(p: f64, t: f64) -> f64 {
    crate::if97_core::core_pt2h(p, t)
}

/// double pt2s(double p,double t) - Calculate specific entropy from pressure and temperature
#[no_mangle]
pub unsafe extern "stdcall" fn pt2s(p: f64, t: f64) -> f64 {
    crate::if97_core::core_pt2s(p, t)
}

/// double pt2v(double p,double t) - Calculate specific volume from pressure and temperature
#[no_mangle]
pub unsafe extern "stdcall" fn pt2v(p: f64, t: f64) -> f64 {
    crate::if97_core::core_pt2v(p, t)
}

/// double pt2x(double p,double t) - Calculate steam quality from pressure and temperature
#[no_mangle]
pub unsafe extern "stdcall" fn pt2x(p: f64, t: f64) -> f64 {
    crate::if97_core::core_pt2x(p, t)
}

/// double ph2t(double p,double h) - Calculate temperature from pressure and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn ph2t(p: f64, h: f64) -> f64 {
    crate::if97_core::core_ph2t(p, h)
}

/// double ph2s(double p,double h) - Calculate specific entropy from pressure and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn ph2s(p: f64, h: f64) -> f64 {
    crate::if97_core::core_ph2s(p, h)
}

/// double ph2v(double p,double h) - Calculate specific volume from pressure and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn ph2v(p: f64, h: f64) -> f64 {
    crate::if97_core::core_ph2v(p, h)
}

/// double ph2x(double p,double h) - Calculate steam quality from pressure and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn ph2x(p: f64, h: f64) -> f64 {
    crate::if97_core::core_ph2x(p, h)
}

/// double ps2t(double p,double s) - Calculate temperature from pressure and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ps2t(p: f64, s: f64) -> f64 {
    crate::if97_core::core_ps2t(p, s)
}

/// double ps2h(double p,double s) - Calculate specific enthalpy from pressure and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ps2h(p: f64, s: f64) -> f64 {
    crate::if97_core::core_ps2h(p, s)
}

/// double ps2v(double p,double s) - Calculate specific volume from pressure and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ps2v(p: f64, s: f64) -> f64 {
    crate::if97_core::core_ps2v(p, s)
}

/// double ps2x(double p,double s) - Calculate steam quality from pressure and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ps2x(p: f64, s: f64) -> f64 {
    crate::if97_core::core_ps2x(p, s)
}

/// double pv2t(double p,double v) - Calculate temperature from pressure and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn pv2t(p: f64, v: f64) -> f64 {
    crate::if97_core::core_pv2t(p, v)
}

/// double pv2h(double p,double v) - Calculate specific enthalpy from pressure and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn pv2h(p: f64, v: f64) -> f64 {
    crate::if97_core::core_pv2h(p, v)
}

/// double pv2s(double p,double v) - Calculate specific entropy from pressure and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn pv2s(p: f64, v: f64) -> f64 {
    crate::if97_core::core_pv2s(p, v)
}

/// double pv2x(double p,double v) - Calculate steam quality from pressure and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn pv2x(p: f64, v: f64) -> f64 {
    crate::if97_core::core_pv2x(p, v)
}

/// double hs2p(double h,double s) - Calculate pressure from enthalpy and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn hs2p(h: f64, s: f64) -> f64 {
    crate::if97_core::core_hs2p(h, s)
}

/// double hs2t(double h,double s) - Calculate temperature from enthalpy and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn hs2t(h: f64, s: f64) -> f64 {
    crate::if97_core::core_hs2t(h, s)
}

/// double hs2v(double h,double s) - Calculate specific volume from enthalpy and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn hs2v(h: f64, s: f64) -> f64 {
    crate::if97_core::core_hs2v(h, s)
}

/// double hs2x(double h,double s) - Calculate steam quality from enthalpy and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn hs2x(h: f64, s: f64) -> f64 {
    crate::if97_core::core_hs2x(h, s)
}

/// double th2p(double t,double h) - Calculate pressure from temperature and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn th2p(t: f64, h: f64) -> f64 {
    crate::if97_core::core_th2p(t, h)
}

/// double th2s(double t,double h) - Calculate specific entropy from temperature and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn th2s(t: f64, h: f64) -> f64 {
    crate::if97_core::core_th2s(t, h)
}

/// double th2v(double t,double h) - Calculate specific volume from temperature and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn th2v(t: f64, h: f64) -> f64 {
    crate::if97_core::core_th2v(t, h)
}

/// double th2x(double t,double h) - Calculate steam quality from temperature and enthalpy
#[no_mangle]
pub unsafe extern "stdcall" fn th2x(t: f64, h: f64) -> f64 {
    crate::if97_core::core_th2x(t, h)
}

/// double ts2p(double t,double s) - Calculate pressure from temperature and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ts2p(t: f64, s: f64) -> f64 {
    crate::if97_core::core_ts2p(t, s)
}

/// double ts2h(double t,double s) - Calculate specific enthalpy from temperature and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ts2h(t: f64, s: f64) -> f64 {
    crate::if97_core::core_ts2h(t, s)
}

/// double ts2v(double t,double s) - Calculate specific volume from temperature and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ts2v(t: f64, s: f64) -> f64 {
    crate::if97_core::core_ts2v(t, s)
}

/// double ts2x(double t,double s) - Calculate steam quality from temperature and entropy
#[no_mangle]
pub unsafe extern "stdcall" fn ts2x(t: f64, s: f64) -> f64 {
    crate::if97_core::core_ts2x(t, s)
}

/// double tv2p(double t,double v) - Calculate pressure from temperature and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn tv2p(t: f64, v: f64) -> f64 {
    crate::if97_core::core_tv2p(t, v)
}

/// double tv2h(double t,double v) - Calculate specific enthalpy from temperature and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn tv2h(t: f64, v: f64) -> f64 {
    crate::if97_core::core_tv2h(t, v)
}

/// double tv2s(double t,double v) - Calculate specific entropy from temperature and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn tv2s(t: f64, v: f64) -> f64 {
    crate::if97_core::core_tv2s(t, v)
}

/// double tv2x(double t,double v) - Calculate steam quality from temperature and specific volume
#[no_mangle]
pub unsafe extern "stdcall" fn tv2x(t: f64, v: f64) -> f64 {
    crate::if97_core::core_tv2x(t, v)
}

/// double px2t(double p,double x) - Calculate temperature from pressure and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn px2t(p: f64, x: f64) -> f64 {
    crate::if97_core::core_px2t(p, x)
}

/// double px2h(double p,double x) - Calculate specific enthalpy from pressure and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn px2h(p: f64, x: f64) -> f64 {
    crate::if97_core::core_px2h(p, x)
}

/// double px2s(double p,double x) - Calculate specific entropy from pressure and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn px2s(p: f64, x: f64) -> f64 {
    crate::if97_core::core_px2s(p, x)
}

/// double px2v(double p,double x) - Calculate specific volume from pressure and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn px2v(p: f64, x: f64) -> f64 {
    crate::if97_core::core_px2v(p, x)
}

/// double tx2p(double t,double x) - Calculate pressure from temperature and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn tx2p(t: f64, x: f64) -> f64 {
    crate::if97_core::core_tx2p(t, x)
}

/// double tx2h(double t,double x) - Calculate specific enthalpy from temperature and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn tx2h(t: f64, x: f64) -> f64 {
    crate::if97_core::core_tx2h(t, x)
}

/// double tx2s(double t,double x) - Calculate specific entropy from temperature and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn tx2s(t: f64, x: f64) -> f64 {
    crate::if97_core::core_tx2s(t, x)
}

/// double tx2v(double t,double x) - Calculate specific volume from temperature and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn tx2v(t: f64, x: f64) -> f64 {
    crate::if97_core::core_tx2v(t, x)
}

/// double hx2p(double h,double x) - Calculate pressure from enthalpy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn hx2p(h: f64, x: f64) -> f64 {
    crate::if97_core::core_hx2p(h, x)
}

/// double hx2t(double h,double x) - Calculate temperature from enthalpy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn hx2t(h: f64, x: f64) -> f64 {
    crate::if97_core::core_hx2t(h, x)
}

/// double hx2s(double h,double x) - Calculate specific entropy from enthalpy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn hx2s(h: f64, x: f64) -> f64 {
    crate::if97_core::core_hx2s(h, x)
}

/// double hx2v(double h,double x) - Calculate specific volume from enthalpy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn hx2v(h: f64, x: f64) -> f64 {
    crate::if97_core::core_hx2v(h, x)
}

/// double sx2p(double s,double x) - Calculate pressure from entropy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn sx2p(s: f64, x: f64) -> f64 {
    crate::if97_core::core_sx2p(s, x)
}

/// double sx2t(double s,double x) - Calculate temperature from entropy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn sx2t(s: f64, x: f64) -> f64 {
    crate::if97_core::core_sx2t(s, x)
}

/// double sx2h(double s,double x) - Calculate specific enthalpy from entropy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn sx2h(s: f64, x: f64) -> f64 {
    crate::if97_core::core_sx2h(s, x)
}

/// double sx2v(double s,double x) - Calculate specific volume from entropy and steam quality
#[no_mangle]
pub unsafe extern "stdcall" fn sx2v(s: f64, x: f64) -> f64 {
    crate::if97_core::core_sx2v(s, x)
}