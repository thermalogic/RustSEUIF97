


// importing common module.
#![allow(warnings)]
mod common;
use assert_approx_eq::assert_approx_eq;

use if97::ps;
use if97::pt;
use if97::region4_pTx::*;
use if97::common::*;

#[test]
fn test_region1_ps() {
  const p: [f64; 3] = [3.0, 80.0, 3.0];
  const s: [f64; 3] = [0.392294792, 0.368563852, 0.258041912e1];
  const T: [f64; 3] = [300.0, 300.0, 500.0];
  for i in 0..=2 {
      assert_approx_eq!(T[i] - 273.15, ps(p[i], s[i], OT), 1.0e-6f64);
  }
  
}
/*
#[test]
fn test_region2_ps() {
 // table29 a，Page29: p,s,T
 const tab29_2a:[[f64;3];3] = [[0.1, 7.5, 0.399517097E+03],
 [0.1, 8.0, 0.514127081E+03],
 [2.5, 8.0, 0.103984917E+04]];

// table29 b，Page29: p,s,T
const tab29_2b:[[f64;3];3] = [[8.0, 6.0, 0.600484040E+03],
 [8.0, 7.5, 0.106495556E+04],
 [90.0, 6.0, 0.103801126E+04]];

// table29 c，Page29: p,s,T
const tab29_2c:[[f64;3];3] = [[20.0, 5.75, 0.697992849E+03],
 [80.0, 5.25, 0.854011484E+03],
 [80.0, 5.75, 0.949017998E+03]];
let mut p:f64=0.0;
let mut s:f64=0.0;
for i in  0..3
{
p = tab29_2a[i][0];
s = tab29_2a[i][1];
assert_approx_eq!(tab29_2a[i][2], ps(p, s,OT)+273.15, 1.0e-5f64);
p = tab29_2b[i][0];
s = tab29_2b[i][1];
assert_approx_eq!(tab29_2b[i][2], ps(p, s,OT)+273.15, 1.0e-5f64);
p = tab29_2c[i][0];
s = tab29_2c[i][1];
assert_approx_eq!(tab29_2c[i][2], ps(p, s,OT)+273.15,1.0e-5f64);
}
}
*/
#[test] 
fn test_region2_ps_iter()
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
  for i in 0..3 {
      let p:f64 = data[i].p;
      let t:f64 = data[i].T- 273.15;
      s = pt(p, t, OS);
      assert_approx_eq!(t, ps(p, s,OT), 1.0e-6f64);
  }
}

#[test]
fn test_region3_ps() {
// Page13 Table 12. Selected temperature values calculated from Eqs. (6) and (7)a
  // Page15 Table 15. Selected specific volume values calculated from Eqs. (8) and (9) a
  const tab12_3a:[[f64;4];3] =[
      [20.0, 3.8, 6.282959869e2, 1.733791463e-3],
      [50.0, 3.6, 6.297158726e2, 1.469680170e-3],
      [100.0, 4.0, 7.056880237e2, 1.555893131e-3]];

  const tab12_3b:[[f64;4];3] = [
      [20.0, 5.0, 6.401176443e2, 6.262101987e-3],
      [50.0, 4.5, 7.163687517e2, 2.332634294e-3],
      [100.0, 5.0, 8.474332825e2, 2.449610757e-3]];

  for i in 0..3
  {
    let mut p:f64 = tab12_3a[i][0];
    let mut s:f64 = tab12_3a[i][1];
    assert_approx_eq!(tab12_3a[i][2]-273.15, ps(p, s,OT));
    assert_approx_eq!(tab12_3a[i][3], ps(p, s,OV));
    p = tab12_3b[i][0];
    s = tab12_3b[i][1];
    assert_approx_eq!(tab12_3b[i][2]-273.15, ps(p, s,OT));
    assert_approx_eq!(tab12_3b[i][3], ps(p, s,OV));
  }
 
}

#[test] 
fn test_region5_ps()
{
  // Table 42, page 40 T,p  v,h,u,s,cp,w
 // Table 42, page 40 T,p  v,h,u,s,cp,w
const data:[[f64;8];3] = [
  [1500.0, 0.5, 1.38455090, 0.521976855e+4, 4527.49310, 9.65408875, 2.61609445, 917.068690],
  [1500.0, 30.0, 0.0230761299, 5167.23514, 4474.95124, 7.72970133, 2.72724317, 928.548002],
  [2000.0, 30.0, 0.0311385219, 6571.22604, 5637.07038, 8.53640523, 2.88569882, 1067.36948]];
  for  i in 0..2
  {
    let p:f64 = data[i][1];
    let s:f64 = data[i][5];
    assert_approx_eq!(data[i][0]-273.15, ps(p, s,OT), 1.0e-4f64);
    assert_approx_eq!(data[i][3], ps(p, s, OH), 1.0e-4f64);
  }
}

#[test] 
fn test_region4_ps() {
 
     // saturation T,p
const tab35:[[f64;2];3] = [[300.0, 0.00353658941],
[500.0, 2.63889776],
[600.0, 12.3443146]];

// saturation p,T
const tab36:[[f64;2];3]=[[0.1, 372.755919],
[1.0, 453.035632],
[10.0, 584.149488]];

let mut s:f64=0.0;
let mut s1:f64=0.0;
let mut s2:f64=0.0;
let mut x:f64=0.0;
for i in 0..3
{
  let p:f64 = tab35[i][1];
  s1 = p2sat_water(p, OS);
  s2 = p2sat_steam(p, OS);
  x = 0.35;
  s = s1 + x * (s2 - s1);
  assert_approx_eq!(x, ps(p, s, OX));
}

for i in 0..3
{
  let p:f64 = tab36[i][0];
  s1 = p2sat_water(p, OS);
  s2 = p2sat_steam(p, OS);
  x = 0.55;
  s = s1 + x * (s2 - s1);
  assert_approx_eq!(x, ps(p, s, OX));
}
}
   



