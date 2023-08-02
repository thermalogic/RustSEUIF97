#![allow(warnings)]
/// The simple example for all input pairs
/// *  (p,t) (p,h) (p,s) (p,v) 
/// *  (t,h) (t,s) (t,v) 
/// *  (p,x) (t,x) (h,x) (t,s) 
/// *  (h,s)
///   
use if97::*;
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
    // extended 
    t=pv(p,v,OT);
    println!("pv: p={p:.6} v={v:.6} t={t:.6}");  
    p=tv(t,v,OP);
    println!("tv t={t:.6} v={v:.6} p={p:.6}");
    
    p=th(t,h,OP);
    println!("th: t={p:.6} h={h:.6} p={p:.6}");    
     
    p=ts(t,s,OP);
    println!("ts: t={p:.6} s={s:.6} p={p:.6}"); 

    // px,tx
    let x: f64 = 0.3; // x= 0.0 saturation water ,x=1.0 saturation steam
    p = 0.1;
    t = px(p, x, OT);
    h = px(p, x, OH);
    s = px(p, x, OS);
    println!("px: p={p:.6} x={x:.6} t={t:.6} h={h:.6} s={s:.6}");

    t = 372.755919 - 273.15;
    p = tx(t, x, OP);
    h = tx(t, x, OH);
    s = tx(t, x, OS);
    println!("tx: p={p:.6} x={x:.6} t={t:.6} h={h:.6}  s={s:.6}");

    p = hx(h, x, OP);
    t = hx(h, x, OT);
    println!("hx: h={h:.6} x={x:.6} p={p:.6} t={t:.6}");

    p = sx(s, x, OP);
    t = sx(s, x, OT);
    println!("sx: s={s:.6} x={x:.6} p={p:.6} t={t:.6}");
}
