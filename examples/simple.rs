#![allow(warnings)]
/// The simple example
use if97::*;
fn main() {
    let p: f64 = 3.0;
    let t: f64 = 300.0 - 273.15;

    let h = pt(p, t, OH);
    let s = pt(p, t, OS);
    let v = pt(p, t, OV);
    let u = pt(p, t, OU);
    let w = pt(p, t, OW);
    println!("pt: p={p:.6} t={t:.6} h={h:.6} s={s:.6} v={v:.6},u={u:.6} w={w:.6}");

    let t = ph(p, h, OT);
    println!("ph: p={p:.6} h={h:.6} t={t:.6}");

    let t = ps(p, s, OT);
    println!("ps: p={p:.6} s={s:.6} t={t:.6}");

    let p = hs(h, s, OP);
    let t = hs(h, s, OT);
    println!("hs: h={h:.6} s={s:.6} p={p:.6} t={t:.6}");

    let p: f64 = 0.1;
    let x: f64 = 0.0; // x= 0.0 saturation water ,x=1.0 saturation steam
    let t: f64 = px(p, x, OT);
    let h: f64 = px(p, x, OH);
    println!("px: p={p:.6} x={x:.6} t={t:.6} h={h:.6}");

    let t: f64 = 372.755919 - 273.15;
    let x: f64 = 1.0; // x= 0.0 saturation water ,x=1.0 saturation steam
    let p: f64 = tx(t, x, OP);
    let h: f64 = tx(t, x, OH);
    println!("tx: p={p:.6} x={x:.6} t={t:.6} h={h:.6}");
}
