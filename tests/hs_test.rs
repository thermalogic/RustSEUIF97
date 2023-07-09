


// importing common module.
#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;

use if97::hs;
use if97::pt;
use if97::px;
use if97::common::*;

#[test]
fn test_region1_hs() {
    const h: [f64; 3] = [0.115331273e+3, 0.184142828e+3, 0.975542239e+3];
    const s: [f64; 3] = [0.392294792, 0.368563852, 0.258041912e1];
    const p: [f64; 3] = [3.0, 80.0, 3.0];
    const T: [f64; 3] = [300.0, 300.0, 500.0];
    for i in 0..=2 {
        assert_approx_eq!(T[i] - 273.15, hs(h[i], s[i], OT), 1.0e-6f64); //real diff: `0.13
        assert_approx_eq!(p[i], hs(h[i], s[i], OP), 1.0e-6f64); //  real diff: `0.181
    }
}

#[test] 
fn test_region2_hs()
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
  let mut s: f64 = 0.0;
  let mut h: f64 = 0.0;
  for i in 0..3 {
      let p:f64 = data[i].p;
      let t:f64 = data[i].T- 273.15;
      h = pt(p, t, OH);
      s = pt(p, t, OS);
      assert_approx_eq!(p, hs(h, s,OP), 1.0e-6f64);
      assert_approx_eq!(t, hs(h, s,OT), 1.0e-6f64);
  }
}

#[test] 
fn test_region3_hs()
{
  //_Backward3_P_hs Table 5 pag 10
  const tab_3a:[[f64;3];3]=[
      [1700.0, 3.8, 25.5570324620],
      [2000.0, 4.2, 45.40873468],
      [2100.0, 4.3, 60.78123340100]];

  const tab_3b:[[f64;3];3]=[
      [2600.0, 5.1, 34.34999263],
      [2400.0, 4.7, 63.6392488750],
      [2700.0, 5.0, 88.3904328]];

  let mut h:f64=0.0;
  let mut s:f64=0.0;
  for i in 0..3
  {
    h = tab_3a[i][0];
    s = tab_3a[i][1];
    assert_approx_eq!(tab_3a[i][2], hs(h, s,OP));
    h = tab_3b[i][0];
    s = tab_3b[i][1];
    assert_approx_eq!(tab_3b[i][2], hs(h, s,OP));
  }
}

#[test] 
fn test_region5_hs()
{
  // Table 42, page 40 T,p  v,h,u,s,cp,w
  const data:[[f64;8];3] = [
          [1500.0, 0.5, 1.38455090, 0.521976855e+4, 4527.49310, 9.65408875, 2.61609445, 917.068690],
          [1500.0, 30.0, 0.0230761299, 5167.23514, 4474.95124, 7.72970133, 2.72724317, 928.548002],
          [2000.0, 30.0, 0.0311385219, 6571.22604, 5637.07038, 8.53640523, 2.88569882, 1067.36948]];

  for  i in 0..2
  {
    let h:f64 = data[i][3];
    let s:f64 = data[i][5];
    assert_approx_eq!(data[i][0]-273.15, hs(h, s, OT), 1.0e-4f64);
    assert_approx_eq!(data[i][1], hs(h, s,OP), 1.0e-4f64);
    assert_approx_eq!(data[i][2], hs(h, s, OV), 1.0e-4f64);
    assert_approx_eq!(data[i][4], hs(h, s, OU), 1.0e-4f64);
    assert_approx_eq!(data[i][6], hs(h, s, OCP), 1.0e-4f64);
    assert_approx_eq!(data[i][7], hs(h, s, OW), 1.0e-4f64);
    
  }
}

#[test]
fn test_region4_hs()
{
   const hsT:[[f64;3];3] =[
    [1800.0, 5.3, 346.8475498],
    [2400.0, 6.0, 425.1373305],
    [2500.0, 5.5, 522.5579013]];

let mut h:f64=0.0;
let mut h1:f64=0.0;
let mut h2:f64=0.0;
let mut s:f64=0.0;
let mut s1:f64=0.0;
let mut s2:f64=0.0;
let mut x:f64=0.0;

let mut p:f64=0.0;
let mut T:f64=0.0;
for i in 0..2
{
  h = hsT[i][0];
  s = hsT[i][1];
  assert_approx_eq!(hsT[i][2]-273.15,hs(h, s,OT));
}

// saturation T,p
const tab35: [[f64; 2]; 3] = [
    [300.0, 0.00353658941],
    [500.0, 2.63889776],
    [600.0, 12.3443146],
];

// saturation p,T
const tab36: [[f64; 2]; 3] = [[0.1, 372.755919], [1.0, 453.035632], [10.0, 584.149488]];

for i in 0..2
{
  p = tab35[i][1];
  T = tab35[i][0];
  s1 = px(p,0.0, OS);
  s2 = px(p,1.0, OS);
  h1 = px(p,0.0, OH);
  h2 = px(p,1.0 ,OH);
  x = 0.60;
  s = s1 + x * (s2 - s1);
  h = h1 + x * (h2 - h1);
  assert_approx_eq!(T-273.15, hs(h, s,OT),1.0e-3f64);
  assert_approx_eq!(p, hs(h, s, OP),1.0e-3f64);
  assert_approx_eq!(x, hs(h, s, OX),1.0e-3f64);
}

for i in 0..2
{
  p = tab36[i][0];
  T = tab36[i][1];
  s1 = px(p,0.0,OS);
  s2 = px(p,1.0,OS);
  h1 = px(p,0.0,OH);
  h2 = px(p,1.0,OH);
  x = 0.05;
  s = s1 + x * (s2 - s1);
  h = h1 + x * (h2 - h1);
  assert_approx_eq!(T-273.15, hs(h, s,OT),1.0e-2f64);
  assert_approx_eq!(p, hs(h, s, OP),1.0e-3f64);
  assert_approx_eq!(x, hs(h, s, OX),1.0e-3f64);
}
}

