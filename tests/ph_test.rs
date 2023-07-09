


// importing common module.
#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;

use if97::ph;
use if97::pt;
use if97::region4_pTx::*;
use if97::common::*;

#[test]
fn test_region1_ph() {

  const p: [f64; 3] = [3.0, 80.0, 3.0];
    const h: [f64; 3] = [0.115331273e+3, 0.184142828e+3, 0.975542239e+3];
    const T: [f64; 3] = [300.0, 300.0, 500.0];
    for i in 0..3 {
      assert_approx_eq!(T[i] - 273.15, ph(p[i], h[i], OT), 1.0e-6f64);
      // real diff: `0.0178
  }
}
/*
#[test]
fn test_region2_ph() {
   //  table15，Page17: p,T,v,h,u,s,cp,ws
 // Table24 a，Page25: p,h,T
 const tab24_2a:[[f64;3];3] = [[0.001, 3000.0, 0.534433241E+03],
 [3.0, 3000.0, 0.575373370E+03],
 [3.0, 4000.0, 0.101077577E+04]];

// table24 b，Page25: p,h,T
const tab24_2b:[[f64;3];3] = [[5.0, 3500.0, 0.801299102E+03],
   [5.0, 4000.0, 0.101531583E+04],
   [25.0, 3500.0, 0.875279054E+03]];

// table24 c，Page25: p,h,T
const tab24_2c:[[f64;3];3]  = [[40.0, 2700.0, 0.743056411E+03],
   [60.0, 2700.0, 0.791137067E+03],
   [60.0, 3200.0, 0.882756860E+03]];
let mut p:f64=0.0;
let mut h:f64=0.0;
for  i in 0..3
{
p = tab24_2a[i][0];
h = tab24_2a[i][1];
assert_approx_eq!(tab24_2a[i][2], ph(p, h,OT)+273.15, 1.0e-5f64);
p = tab24_2b[i][0];
h = tab24_2b[i][1];
assert_approx_eq!(tab24_2b[i][2], ph(p, h,OT)+273.15, 1.0e-5f64);
p = tab24_2c[i][0];
h = tab24_2c[i][1];
assert_approx_eq!(tab24_2c[i][2], ph(p, h,OT)+273.15, 1.0e-5f64);
}
}
*/
#[test] 
fn test_region2_ph_iter()
{
  pub struct TestData {
    pub p: f64,
    pub T: f64,
    pub v: f64,
    pub h: f64,
    pub s: f64,
    pub u: f64,
    pub cp: f64,
    pub w: f64  
  }
  //  Table15，Page17: p,T,v,h,u,s,cp,w
  const data:[TestData;3] = [
    TestData{p:0.0035,T:300.0, v:0.394913866E+02, h:0.254991145E+04, u:0.241169160E+04, s:0.852238967E+01, cp:0.191300162E+01, w:0.427920172E+03},
    TestData{p:0.0035, T:700.0, v:0.923015898E+02, h:0.333568375E+04, u:0.301262819E+04, s:0.101749996E+02, cp:0.208141274E+01,w:0.644289068E+03},
    TestData{p:30.0, T:700.0, v:0.542946619E-02, h:0.263149474E+04, u:0.246861076E+04, s:0.517540298E+01, cp:0.103505092E+02,w:0.480386523E+03}];
  let mut h: f64 = 0.0;
  for i in 0..3 {
      let p:f64 = data[i].p;
      let t:f64 = data[i].T- 273.15;
      h = pt(p, t, OH);
      assert_approx_eq!(t, ph(p, h,OT), 1.0e-6f64);
  }
}



#[test]
fn test_region3_ph() {
// Supp-Tv(ph,ps)3-2014.pdf
  //      Page 8 Table 5. Selected temperature values calculated from Eqs. (2) and (3) a
  //                  p,h->T,
  //       Page 10Table 8. Selected specific volume values calculated from Eqs. (4) and (5) a
  //                 p,h,->v     p,h,T,v
  const tab5_3a:[[f64;4];3] = [
      [20.0, 1700.0, 6.293083892e2, 1.749903962e-3],
      [50.0, 2000.0, 6.905718338e2, 1.908139035e-3],
      [100.0, 2100.0, 7.336163014e2, 1.676229776e-3]];
  const tab5_3b:[[f64;4];3]  = [
      [20.0, 2500.0, 6.418418053e2, 6.670547043e-3],
      [50.0, 2400.0, 7.351848618e2, 2.801244590e-3],
      [100.0, 2700.0, 8.420460876e2, 2.404234998e-3]];

  for i in 0..3
  {
    let mut p:f64 = tab5_3a[i][0];
    let mut h:f64 = tab5_3a[i][1];
    assert_approx_eq!(tab5_3a[i][2]-273.15, ph(p, h,OT));
    assert_approx_eq!(tab5_3a[i][3], ph(p, h,OV));
    p = tab5_3b[i][0];
    h = tab5_3b[i][1];
    assert_approx_eq!(tab5_3b[i][2]-273.15, ph(p, h,OT));
    assert_approx_eq!(tab5_3b[i][3], ph(p, h,OV));
  }
 
}

#[test] 
fn test_region5_ph()
{
  // Table 42, page 40 T,p  v,h,u,s,cp,w
const data:[[f64;8];3] = [
  [1500.0, 0.5, 1.38455090, 0.521976855e+4, 4527.49310, 9.65408875, 2.61609445, 917.068690],
  [1500.0, 30.0, 0.0230761299, 5167.23514, 4474.95124, 7.72970133, 2.72724317, 928.548002],
  [2000.0, 30.0, 0.0311385219, 6571.22604, 5637.07038, 8.53640523, 2.88569882, 1067.36948]];

  for  i in 0..2
  {
    let p:f64 = data[i][1];
    let h:f64 = data[i][3];
    assert_approx_eq!(data[i][0]-273.15 ,ph(p, h,OT));
    assert_approx_eq!(data[i][5], ph(p, h, OS));
  }
}

#[test] 
fn test_region4_ph() {
 
     // saturation T,p
const tab35:[[f64;2];3] = [[300.0, 0.00353658941],
[500.0, 2.63889776],
[600.0, 12.3443146]];

// saturation p,T
const tab36:[[f64;2];3]=[[0.1, 372.755919],
[1.0, 453.035632],
[10.0, 584.149488]];

  let mut h:f64=0.0;
  let mut h1:f64=0.0;
  let mut h2:f64=0.0;
  let mut x:f64=0.0;
  for i in 0..3 
  {
    let p:f64=tab35[i][1];
    h1=p2sat_water(p,OH);
    h2=p2sat_steam(p,OH);
    x=0.35;
    h=h1+x*(h2-h1);
    assert_approx_eq!(x, ph(p,h,OX));
  }

  for i in 0..3 
  {
    let p:f64=tab36[i][0];
    h1=p2sat_water(p,OH);
    h2=p2sat_steam(p,OH);
    x=0.55;
    h=h1+x*(h2-h1);
    assert_approx_eq!(x, ph(p,h,OX));    
  }
   
}


