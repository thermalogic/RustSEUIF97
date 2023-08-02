//! Region 4
//! * (p,h) (p,s) (h,s) (p,x) (t,x)
//! * (p,v),(t,v),(t,h),(t,s)
//! * (h,x),(s,x)

use crate::common::constant::*;
use crate::common::propertry_id::*;
use crate::r3::*;
use crate::r4::*;

pub fn p_sat(t: f64) -> f64 {
    p_saturation(t + K)
}

pub fn t_sat(p: f64) -> f64 {
    T_saturation(p) - K
}

/// for TC_WATER && p == PC_WATER
pub fn pT_reg4(p: f64, T: f64, o_id: i32) -> f64 {
    if T == TC_WATER && p == PC_WATER {
        return Td_reg3(TC_WATER, DC_WATER, o_id);
    } else {
        return 4.0 as f64;
    }
}

pub fn ph_reg4(p: f64, h: f64, o_id: i32) -> f64 {
    let h1: f64 = p2sat_water(p, OH);
    let h2: f64 = p2sat_steam(p, OH);
    let x: f64 = (h - h1) / (h2 - h1);
    if o_id == OX {
        return x;
    }
    px_reg4(p, x, o_id)
}

pub fn ps_reg4(p: f64, s: f64, o_id: i32) -> f64 {
    let s1: f64 = p2sat_water(p, OS);
    let s2: f64 = p2sat_steam(p, OS);
    let x: f64 = (s - s1) / (s2 - s1);
    if o_id == OX {
        return x;
    }
    px_reg4(p, x, o_id)
}

pub fn hs_reg4(h: f64, s: f64, o_id: i32) -> f64 {
    // for T<623.15 only
    let T: f64 = hs2T_reg4(h, s);
    if o_id == OT {
        return T - K;
    }
    let p: f64 = p_saturation(T);
    if o_id == OP {
        return p;
    }

    let h1: f64 = p2sat_water(p, OH);
    let h2: f64 = p2sat_steam(p, OH);

    let x: f64 = (h - h1) / (h2 - h1);
    if o_id == OX {
        return x;
    }

    px_reg4(p, x, o_id)
}

/// Region 4 - The extended input pair
///       (p,v) ->x
///       (t,h),(t,s),(t,v)->x
///   x: Steam quality

///  (p,v) ->x
pub fn pv_reg4(p: f64, v: f64, o_id: i32) -> f64 {
    let x: f64 = pv2x_reg4(p, v);
    if o_id == OX {
        return x;
    }
    px_reg4(p, x, o_id)
}

///  (t,v)
pub fn tv_reg4(t: f64, v: f64, o_id: i32) -> f64 {
    let x: f64 = Tv2x_reg4(t + 273.15, v);
    if o_id == OX {
        return x;
    }
    Tx_reg4(t + 273.15, x, o_id)
}

///  (t,h)
pub fn th_reg4(t: f64, h: f64, o_id: i32) -> f64 {
    let x: f64 = Th2x_reg4(t + 273.15, h);
    if o_id == OX {
        return x;
    }
    Tx_reg4(t + 273.15, x, o_id)
}

///  (t,s)
pub fn ts_reg4(t: f64, s: f64, o_id: i32) -> f64 {
    let x: f64 = Ts2x_reg4(t + 273.15, s);
    if o_id == OX {
        return x;
    }
    Tx_reg4(t + 273.15, x, o_id)
}

/// helper of (h,x),(s,x) -failed!
fn helper_hxsx_reg4(y: f64, x: f64, y_id: i32) -> f64 {
    let mut T: f64 = 0.0;
    let mut T1: f64 = T_MIN4 + 0.01;
    let mut T2: f64 = 0.5 * (T_MIN4 + T_MAX4);
    let sw1: f64 = T2sat_water(T1, y_id);
    let ss1: f64 = T2sat_steam(T1, y_id);
    let sw2: f64 = T2sat_water(T2, y_id);
    let ss2: f64 = T2sat_steam(T2, y_id);
    let mut x1: f64 = (y - sw1) / (ss1 - sw1);
    let mut x2: f64 = (y - sw2) / (ss2 - sw2);
    let mut ft: f64 = x - x1;
    let mut f: f64 = x - x2;
    let mut swap: f64 = 0.0;

    let mut pt: f64 = 0.0;
    if ft.abs() < f.abs() {
        T = T1;
        pt = T2;
        swap = ft;
        ft = f;
        f = swap;
    } else {
        pt = T1;
        T = T2;
    };
    let mut again_t: bool = true;
    let mut n: i32 = 0;
    let mut dt: f64 = 0.0;
    let mut xx: f64 = 0.0;
    let mut sw: f64 = 0.0;
    let mut ss: f64 = 0.0;

    while again_t {
        dt = (pt - T) * f / (f - ft);
        pt = T;
        ft = f;
        T += dt;
        if (T > T_MIN4) && (T < T_MAX4) {
            sw = T2sat_water(T, y_id);
            ss = T2sat_steam(T, y_id);
            xx = (y - sw) / (ss - sw);
            f = x - xx;
        } else {
            if T < T_MIN4 {
                T = T_MIN4;
            } else if T > T_MAX4 {
                T = T_MAX4;
            }
        };
        if f.abs() > 0.001 {
            n += 1;
        } else {
            again_t = false;
        }
        if n > 10000 {
            again_t = false;
        };
    }
    if n > 10000 {
        return INVALID_VALUE as f64;
    }
    T
}

/// function for getting the steam quality,residuals: x(T,y)-x
fn Ty2x_residuals(T: f64, x: f64, y: f64, y_id: i32) -> f64 {
    let sw = T2sat_water(T, y_id);
    let ss = T2sat_steam(T, y_id);
    (y - sw) / (ss - sw) - x
}

/// Bisection for the root : Ty2x_residuals(T,x, y, y_id)=0
fn bisection_reg4(
    x: f64,
    y: f64,
    y_id: i32,
    mut Tl: f64,
    mut Tr: f64,
    tol: f64,
    maxiter: i32,
) -> f64 {
    let mut T: f64 = 0.0;
    let mut fl: f64 = Ty2x_residuals(Tl, x, y, y_id); // residual for left  bound
    let mut fr: f64 = Ty2x_residuals(Tr, x, y, y_id); //resdiual for right bound
    let mut f: f64 = 0.0;
    let mut numIters: i32 = 0;

    for i in 0..maxiter {
        numIters += 1;
        // get midpoint
        T = 0.5 * (Tl + Tr);
        // evaluate resdiual at midpoint
        f = Ty2x_residuals(T, x, y, y_id);
        //  check for convergence
        if f.abs() < tol {
            break;
        };

        // reset the bounds
        if f * fl < 0.0 {
            // move right bound info to mid
            Tr = T;
            fr = f;
        } else {
            // move left bound info to mid
            Tl = T;
            fl = f;
        }
    }
    T
}

/// (h,x,o_id)
pub fn hx_reg4(h: f64, x: f64, o_id: i32) -> f64 {
    //let T: f64 = helper_hxsx_reg4(h, x, OH);
    let mut Tl: f64 = T_MIN4;
    let mut Tr: f64 = T_MAX4;
    let T: f64 = bisection_reg4(x, h, OH, Tl, Tr, 0.00001, 1000);
    if o_id == OT {
        return T - 273.15;
    };
    Tx_reg4(T, x, o_id)
}

/// (s,x,o_id)
pub fn sx_reg4(s: f64, x: f64, o_id: i32) -> f64 {
    //let T: f64 = helper_hxsx_reg4(s, x, OS);
    let mut Tl: f64 = T_MIN4;
    let mut Tr: f64 = T_MAX4;
    let T: f64 = bisection_reg4(x, s, OS, Tl, Tr, 0.00001, 1000);
    if o_id == OT {
        return T - 273.15;
    };
    Tx_reg4(T, x, o_id)
}
