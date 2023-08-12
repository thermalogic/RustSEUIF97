#![allow(warnings)]
/// Example: input pairs
/// *  (p,t) (p,h) (p,s) ( (p,v)
/// *  (t,v),(t,h),(t,s)
/// *  (h,s)
/// *  (p,x),(t,x),(h,x),(s,x)  
///  
///  Run:
///    cargo run -r --example input_pair
///
use seuif97::*;
fn main() {
    let p: f64 = 3.0;
    let t: f64 = 300.0 - 273.15;
    // pt
    let h = pt(p, t, OH);
    let s = pt(p, t, OS);

    println!("\tpt\n\t   p={:.6} t={:.6} h={:.6} s={:6}", p, t, h, s);
    // ph -> t
    println!("\tph\n\t   p={:.6} t={:.6} h={:.6} s={:6}", p, ph(p, h, OT), h, s);
    // ps -> t
    println!("\tps\n\t   p={:.6} t={:.6} h={:.6} s={:6}", p, ps(p, s, OT), h, s);
    // hs ->t,p
    println!("\ths\n\t   p={:.6} t={:.6} h={:.6} s={:6}", hs(h, s, OP), hs(h, s, OT), h, s);

    //extened (t,s) ->p
    println!("\tpt\n\t   p={:.6} t={:.6}  s={:.6}", p, t, s);
    println!("\tts\n\t   p={:.6} t={:.6}  s={:.6}", ts(t, s, OP), t, s);
    //extened (t,h) ->p
    println!("\tpt\n\t   p={:.6} t={:.6}  h={:.6}", p, t, h);
    println!("\tth\n\t   p={:.6} t={:.6}  h={:.6}", th(t, h, OP), t, h);

    //extened (p,v),(t,v)
    let v = pt(p, t, OV);

    println!("\tpt\n\t   p={:.6} t={:.6} v={:.6}", p, t, v);
    //extened (p,v) ->t
    println!("\tpv\n\t   p={:.6} t={:.6} v={:.6}", p, pv(p, v, OT), v);
    //extened (t,v) ->p
    println!("\ttv\n\t   p={:.6} t={:.6} v={:.6}", tv(t, v, OP), t, v);

    // px x= 0.0 saturation water ,x=1.0 saturation steam
    let p: f64 = 0.1;
    let t_sat: f64 = px(p, 0.0, OT);
    let h_sat_w: f64 = px(p, 0.0, OH);
    let h_sat_s: f64 = px(p, 1.0, OH);
    println!("\tpx\n\t   p_sat={:.6} t_sat={:.6} h_sat_w={:.6} h_sat_s={:6}", p, t_sat, h_sat_w, h_sat_s);

    // tx x= 0.0 saturation water ,x=1.0 saturation steam
    let p_sat: f64 = tx(t_sat, 0.0, OP);
    let h_sat_w: f64 = tx(t_sat, 0.0, OH);
    let h_sat_s: f64 = tx(t_sat, 1.0, OH);
    println!("\ttx\n\t   p_sat={:.6} t_sat={:.6} h_sat_w={:.6} h_sat_s={:6}", p_sat, t_sat, h_sat_w, h_sat_s);

    let x: f64 = 0.6;
    // (h,x)->p,t
    let hxv: f64 = px(p, x, OH);
    println!("\tpx\n\t   p={:.6} t_sat={:.6} x={:.6} h={:.6}", p, t_sat, x, hxv);
    println!("\thx\n\t   p={:.6} t_sat={:.6} x={:.6} h={:.6}", hx(hxv, x, OP), hx(hxv, x, OT), x, hxv);
    // (s,x)->p,t
    let sxv: f64 = px(p, x, OS);
    println!("\tpx\n\t   p_sat={:.6} t_sat={:.6} x={:.6} s={:.6}", p, t_sat, x, sxv);
    println!("\tsx\n\t   p_sat={:.6} t_sat={:.6} x={:.6} s={:.6}", sx(sxv, x, OP), sx(sxv, x, OT), x, sxv);
}
