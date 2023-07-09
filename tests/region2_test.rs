
#![allow(warnings)]
// importing common module.
use assert_approx_eq::assert_approx_eq;

use if97::pt_reg2;
use if97::ph_reg2;
use if97::ps_reg2;
use if97::hs_reg2;
use if97::common::*;

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

#[test]
fn test_pt()
{

  for i in 0..2
  {
    let p:f64 = data[i].p;
    let t:f64 = data[i].T- 273.15;
    assert_approx_eq!(data[i].v, pt_reg2(p, t, OV), 1.0e-6f64); 
    assert_approx_eq!(data[i].h, pt_reg2(p, t, OH), 1.0e-5f64); 
    assert_approx_eq!(data[i].s, pt_reg2(p, t, OS), 1.0e-6f64);
    assert_approx_eq!(data[i].u, pt_reg2(p, t, OU), 1.0e-5f64);
    assert_approx_eq!(data[i].cp,pt_reg2(p,t, OCP), 1.0e-6f64);
    assert_approx_eq!(data[i].w, pt_reg2(p,t, OW), 1.0e-6f64);
   
  }
}

#[test] 
fn test_ph()
{
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
    assert_approx_eq!(tab24_2a[i][2], ph_reg2(p, h,OT)+273.15, 1.0e-1f64);
    p = tab24_2b[i][0];
    h = tab24_2b[i][1];
    assert_approx_eq!(tab24_2b[i][2], ph_reg2(p, h,OT)+273.15, 1.0e-1f64);
    p = tab24_2c[i][0];
    h = tab24_2c[i][1];
    assert_approx_eq!(tab24_2c[i][2], ph_reg2(p, h,OT)+273.15, 1.0e-1f64);
  }
}

#[test] 
fn test_ph_iter()
{
  let mut h: f64 = 0.0;
  for i in 0..3 {
      let p:f64 = data[i].p;
      let t:f64 = data[i].T- 273.15;
      h = pt_reg2(p, t, OH);
      assert_approx_eq!(t, ph_reg2(p, h,OT), 1.0e-6f64);
  }
}


#[test] 
fn test_ps()
{
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
    assert_approx_eq!(tab29_2a[i][2], ps_reg2(p, s,OT)+273.15, 1.0e-2f64);
    p = tab29_2b[i][0];
    s = tab29_2b[i][1];
    assert_approx_eq!(tab29_2b[i][2], ps_reg2(p, s,OT)+273.15, 1.0e-2f64);
    p = tab29_2c[i][0];
    s = tab29_2c[i][1];
    assert_approx_eq!(tab29_2c[i][2], ps_reg2(p, s,OT)+273.15,1.0e-2f64);
   }
}

#[test] 
fn test_ps_iter()
{
  let mut s: f64 = 0.0;
  for i in 0..3 {
      let p:f64 = data[i].p;
      let t:f64 = data[i].T- 273.15;
      s = pt_reg2(p, t, OS);
      assert_approx_eq!(t, ps_reg2(p, s,OT), 1.0e-6f64);
  }
}


#[test] 
fn test_hs()
{
  // Supp-PHS12-2014 ,table9 a，Page10: h,s,p
  const tab9_hs_reg2a:[[f64;3];3]= [[2800.0, 6.5, 1.371012767],
                                 [2800.0, 9.5, 1.879743844e-03],
                                 [4100.0, 9.5, 1.024788997e-01]];

  // Supp-PHS12-2014 ,table9 b，Page10: h,s,p
  const tab9_hs_reg2b:[[f64;3];3] = [[2800.0, 6.0, 4.793911442],
                                 [3600.0, 6.0, 8.395519209e+01],  // `75.59971677713756?????
                                 [3600.0, 7.0, 7.527161441]];

  // Supp-PHS12-2014 ,table9 c Page10: h,s,p
  const tab9_hs_reg2c:[[f64;3];3]  = [[2800.0, 5.1, 9.439202060E+01],
                                 [2800.0, 5.8, 8.414574124],
                                 [3400.0, 5.8, 8.376903879e+01]];
  let mut h:f64=0.0;
  let mut s:f64=0.0;
  for i in 0..3
  {
    h = tab9_hs_reg2a[i][0];
    s = tab9_hs_reg2a[i][1];
    assert_approx_eq!(tab9_hs_reg2a[i][2], hs_reg2(h, s,OP),1.0e-2f64);
    h = tab9_hs_reg2b[i][0];
    s = tab9_hs_reg2b[i][1];
    assert_approx_eq!(tab9_hs_reg2b[i][2], hs_reg2(h, s,OP),1.0e-2f64);
    h = tab9_hs_reg2c[i][0];
    s = tab9_hs_reg2c[i][1];
    assert_approx_eq!(tab9_hs_reg2c[i][2], hs_reg2(h, s,OP),1.0e-2f64);
  }
}

#[test] 
fn test_hs_iter()
{
  let mut s: f64 = 0.0;
  let mut h: f64 = 0.0;
  for i in 0..3 {
      let p:f64 = data[i].p;
      let t:f64 = data[i].T- 273.15;
      h = pt_reg2(p, t, OH);
      s = pt_reg2(p, t, OS);
      assert_approx_eq!(p, hs_reg2(h, s,OP), 1.0e-6f64);
      assert_approx_eq!(t, hs_reg2(h, s,OT), 1.0e-6f64);
  }
}



