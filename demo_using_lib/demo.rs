//  The Rust example of using the rust library 
//  rustc -O demo.rs --extern seuif97=../target/release/libseuif97.rlib
#![allow(warnings)]
fn main() {
    
    let p:f64 = 3.0;
    let t:f64= 300.0-273.15;
   
    let h=seuif97::pt(p,t,seuif97::OH);
    // the funnction with o_id and region for faster calculation
    let s=seuif97::pt(p,t,(seuif97::OS,1));
    let v=seuif97::pt2v(p,t);
    println!("p={p:.6} t={t:.6} h={t:.6} s={s:.6} v={v:.6}");   
}
