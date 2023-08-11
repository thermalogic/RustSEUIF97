//  The Rust example of using the shared library 
//  rustc -O demo.rs --extern seuif97=../target/release/libseuif97.rlib
#![allow(warnings)]
fn main() {
    
    let p:f64 = 3.0;
    let t:f64= 300.0-273.15;
   
    let h=seuif97::pt(p,t,seuif97::OH);
    let s=seuif97::pt(p,t,seuif97::OS);
    // set the region
    let v=seuif97::pt(p,t,(seuif97::OV,1));
    println!("p={p:.6} t={t:.6} h={t:.6} s={s:.6} v={v:.6}");   
}
