#![allow(warnings)]
fn main() {
    
    let p:f64 = 3.0;
    let t:f64= 300.0-273.15;
    // pt
    let h=if97::pt(p,t,if97::OH);
    let s=if97::pt(p,t,if97::OS);
    println!("pt p={:.6} t={:.6} h={:.6} s={:6}",p,t,h,s);    
    // ph
    let t=if97::ph(p,h,if97::OT);
    println!("ph p={:.6} t={:.6} h={:.6} s={:6}",p,t,h,s);    
    // ps
    let t=if97::ps(p,s,if97::OT);
    println!("ps p={:.6} t={:.6} h={:.6} s={:6}",p,t,h,s);    
    // hs
    let t=if97::hs(h,s,if97::OT);
    let p=if97::hs(h,s,if97::OP);
    println!("hs p={:.6} t={:.6} h={:.6} s={:6}",p,t,h,s);    
    // px
    let p:f64=0.1;  
    let t_sat:f64 =if97::px(p,0.0,if97::OT);
    let h_sat_w:f64 =if97::px(p,0.0,if97::OH);
    let h_sat_s:f64 =if97::px(p,1.0,if97::OH);
    println!("px p_sat={:.6} t_sat={:.6} h_sat_w={:.6} h_sat_s={:6}", p,t_sat,h_sat_w,h_sat_s);   
    //tx
    let p_sat:f64 =if97::tx(t_sat,0.0,if97::OP);
    let h_sat_w:f64 =if97::tx(t_sat,0.0,if97::OH);
    let h_sat_s:f64 =if97::tx(t_sat,1.0,if97::OH);
    println!("tx p_sat={:.6} t_sat={:.6} h_sat_w={:.6} h_sat_s={:6}",p_sat,t_sat,h_sat_w,h_sat_s);
   

   

    

}
    