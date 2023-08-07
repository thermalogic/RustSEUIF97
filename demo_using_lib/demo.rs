//  The Rust example of using a Library 
//  rustc -O demo.rs --extern if97=../target/release/libif97.rlib
fn main() {
    
    let p:f64 = 3.0;
    let t:f64= 300.0-273.15;
   
    let h=if97::pt(p,t,if97::OH);
    let s=if97::pt(p,t,if97::OS);
    // set the region
    let v=if97::pt(p,t,(if97::OV,1));
    println!("p={p:.6} t={t:.6} h={t:.6} s={s:.6} v={v:.6}");   
}
