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
    
    let P_MIN: f64 = 0.000611212677444;
    let s_r5_1=pt(50.0, 1073.15-275.15,(OS,5));
    let s_r5_2=pt(P_MIN, 2273.15-273.15,(OS,5));
    let h_r5_1=pt(50.0, 1073.15-273.15,(OH,5)); 
    let h_r5_2= pt(P_MIN, 2273.15-273.15,(OH,5));
    println!(" s_r5_1={s_r5_1:.8} s_r5_2={s_r5_2:.8} h_r5_1={h_r5_1:.8} h_r5_2={h_r5_2:.8}");   

    let s_min= pt(100.0, 273.15-273.15,(OS,1));
    let s_max= pt(P_MIN, 1073.15-273.15,(OS,2));
    println!(" s_min={s_min:.8} s_max={s_max:.8}");   

     let s13: f64 = pt(100.0, 623.15-273.15,(OS,1));
     println!(" s13={s13:.8}");   
     let Ps_623: f64 = 16.5291642526045;
     let s13s: f64 = pt(Ps_623, 623.15-273.15,(OS,1));
     println!(" s13s={s13s:.8}");   
     let h23min: f64 = pt(Ps_623, 623.15-273.15,(OH,2));
     println!(" h23min={h23min:.8}");  
     let h23max: f64 = pt(100.0, 863.15-273.15,(OH,2)); 
     println!(" h23max={h23max:.8}");  
     let sTPmax: f64 = pt(100.0, 1073.15-273.15,(OS,2));
     println!(" sTPmax={sTPmax:.8}");  
     let s2ab: f64 = pt(4.0, 1073.15-273.15,(OS,2)); 
     println!(" s2ab={s2ab:.8}"); 
     let hmin: f64 = pt(P_MIN, 273.15-273.15,(OH,2)); 
     println!(" hmin={hmin:.8}");  

    let h4l: f64 = pt(P_MIN, 273.15-273.15,(OH,1));
    let s4l: f64 = pt(P_MIN, 273.15-273.15,(OS,1));
    let h4v: f64 = pt(P_MIN, 273.15-273.15,(OH,2));
    let s4v: f64 = pt(P_MIN, 273.15-273.15,(OS,2));
    println!(" h4l={h4l:.8} s4l={s4l:.8} h4v={h4v:.8} s4v={s4v:.8}");  
}
