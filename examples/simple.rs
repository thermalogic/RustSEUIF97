#![allow(warnings)]
/// The simple example of the input pairs
/// *  (p,t) (p,h) (p,s)
/// *  (h,s)
///
///   Run:
///    cargo run -r --example simple
///   
use seuif97::*;
fn main() {
    let mut p: f64 = 3.0;
    let mut t: f64 = 300.0 - 273.15;

    let mut h = pt(p, t, OH);
    let mut s = pt(p, t, OS);
    let mut v = pt(p, t, OV);
    // set the region 1
    let u = pt(p, t, (OU, 1));
    let w = pt(p, t, (OW, 1));

    println!("pt: p={p:.6} t={t:.6} h={h:.6} s={s:.6} v={v:.6},u={u:.6} w={w:.6}");

    t = ph(p, h, OT);
    println!("ph: p={p:.6} h={h:.6} t={t:.6}");

    t = ps(p, s, OT);
    println!("ps: p={p:.6} s={s:.6} t={t:.6}");

    p = hs(h, s, OP);
    t = hs(h, s, OT);
    println!("hs: h={h:.6} s={s:.6} p={p:.6} t={t:.6}");

    let p: f64 = 0.1;
    let t_sat: f64 = px(p, 0.0, OT);
    println!("t={t_sat:.6}");

    let p: f64 = 0.5;
    let t_sat: f64 = px(p, 0.5, OT);
    println!("t={t_sat:.6}");
}
