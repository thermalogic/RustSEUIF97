#![allow(warnings)]
/// The simple example
 
fn main() {
    
    let p:f64 = 3.0;
    let t:f64= 300.0-273.15;
    
    let h=if97::pt(p,t,if97::OH);
    let s=if97::pt(p,t,if97::OS);
    let v=if97::pt(p,t,if97::OV);
    let u=if97::pt(p,t,if97::OU);
    let w=if97::pt(p,t,if97::OW);
    println!("p={} t={} h={} s={} v={},u={} w={}",p,t,h,s,v,u,w);      

}
    