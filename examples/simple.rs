#![allow(warnings)]
/// The example to test the types of funtion
/// 
///  * pt(p, t, OH); 
///  * pt(p, t, (OS, 1));
///  * pt2v
///
///  Run:
///    cargo run -r --example simple
///   
use seuif97::*;
fn main() {
    let mut p: f64 = 3.0;
    let mut t: f64 = 300.0 - 273.15;
    // Universal Functions (with o_id parameter)
    let mut h = pt(p, t, OH);
    //  Universal Functions (with o_id parameter and region parameter)
    let s = pt(p, t, (OS, 1));
    // Direct Property Functions
    let mut v = pt2v(p, t);
    println!("pt: p={p:.6} t={t:.6} h={h:.6} s={s:.6} v={v:.6}");   
}
