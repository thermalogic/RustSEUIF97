
//cargo test --test region3_test
#![allow(warnings)]
// importing common module.
mod common;
use assert_approx_eq::assert_approx_eq;

use if97::Td_reg3;
use if97::pT2v_reg3;
use if97:: ph2T_reg3;
use if97:: ph2v_reg3;
use if97:: ps2T_reg3;
use if97:: ps2v_reg3;
use if97::hs2p_reg3;
use if97::common::*;


#[test] 
fn test_Td()
{
  // Table 33. Thermodynamic property values calculated from Eq. (28) for selected values of T and  a
  // T,d,p,h,u,s,cp,w
  const tab33:[[f64;8];3]=[[650., 500., 0.255837018E2, 0.186343019E4, 0.181226279E4, 0.405427273E1, 0.138935717E2, 0.502005554E3],
      [650., 200., 0.222930643E2, 0.237512401E4, 0.226365868E4, 0.485438792E1, 0.446579342E2, 0.383444594E3],
      [750., 500., 0.783095639E2, 0.225868845E4, 0.210206932E4, 0.446971906E1, 0.634165359E1, 0.760696041E3]];

  for i in 0..3 
  {
    let T:f64 = tab33[i][0];
    let d:f64 = tab33[i][1];
    assert_approx_eq!(tab33[i][2], Td_reg3(T, d,OP),1.0e-5f64);
    assert_approx_eq!(tab33[i][3], Td_reg3(T, d,OH),1.0e-5f64);
    assert_approx_eq!(tab33[i][4], Td_reg3(T, d,OU),1.0e-5f64);
    assert_approx_eq!(tab33[i][5], Td_reg3(T, d,OS),1.0e-5f64);
    assert_approx_eq!(tab33[i][6], Td_reg3(T, d,OCP),1.0e-5f64);
    assert_approx_eq!(tab33[i][7], Td_reg3(T, d,OW),1.0e-5f64);
  }
}

#[test] 
fn test_pT()
{
  //_Backward3_PT v， p,T
  const tab:[[f64;3];52] = [
      [1.470853100e-3, 50.0, 630.0], // a
      [1.503831359e-3, 80.0, 670.0],
      [2.204728587e-3, 50.0, 710.0], // b
      [1.973692940e-3, 80.0, 750.0],
      [1.761696406e-3, 20.0, 630.0], // c
      [1.819560617e-3, 30.0, 650.0],
      [2.245587720e-3, 26.0, 656.0], // d
      [2.506897702e-3, 30.0, 670.0],
      [2.970225962e-3, 26.0, 661.0], // e
      [3.004627086e-3, 30.0, 675.0],
      [5.019029401e-3, 26.0, 671.0], // f
      [4.656470142e-3, 30.0, 690.0],
      [2.163198378e-3, 23.6, 649.0], // g
      [2.166044161e-3, 24.0, 650.0],
      [2.651081407e-3, 23.6, 652.0], // h
      [2.967802335e-3, 24.0, 654.0],
      [3.273916816e-3, 23.6, 653.0], // i
      [3.550329864e-3, 24.0, 655.0],
      [4.545001142e-3, 23.5, 655.0], // j
      [5.100267704e-3, 24.0, 660.0],
      [6.109525997e-3, 23.0, 660.0], // k
      [6.427325645e-3, 24.0, 670.0],
      [2.117860851e-3, 22.6, 646.0], // l
      [2.062374674e-3, 23.0, 646.0],
      [2.533063780e-3, 22.6, 648.6], // m
      [2.572971781e-3, 22.8, 649.3],
      [2.923432711e-3, 22.6, 649.0], // n
      [2.913311494e-3, 22.8, 649.7],
      [3.131208996e-3, 22.6, 649.1], // o
      [3.221160278e-3, 22.8, 649.9],
      [3.715596186e-3, 22.6, 649.4], // p
      [3.664754790e-3, 22.8, 650.2],
      [1.970999272e-3, 21.1, 640.0], // q
      [2.043919161e-3, 21.8, 643.0],
      [5.251009921e-3, 21.1, 644.0], // r
      [5.256844741e-3, 21.8, 648.0],
      [1.932829079e-3, 19.1, 635.0], // s
      [1.985387227e-3, 20.0, 638.0],
      [8.483262001e-3, 17.0, 626.0], // t
      [6.227528101e-3, 20.0, 640.0],
      [2.268366647e-3, 21.5, 644.6], // u
      [2.296350553e-3, 22.0, 646.1],
      [2.832373260e-3, 22.5, 648.6], // v
      [2.811424405e-3, 22.3, 647.9],
      [3.694032281e-3, 22.15, 647.5], // w
      [3.622226305e-3, 22.3, 648.1],
      [4.528072649e-3, 22.11, 648.0], // x
      [4.556905799e-3, 22.3, 649.0],
      [2.698354719e-3, 22.0, 646.84], // y
      [2.717655648e-3, 22.064, 647.05],
      [3.798732962e-3, 22.0, 646.89], // z
      [3.701940010e-3, 22.064, 647.15]];

  for i in 0..52
  {
    let p:f64 = tab[i][1];
    let T:f64 = tab[i][2];
    assert_approx_eq!(tab[i][0], pT2v_reg3(p, T));
  }
}

#[test] 
fn test_ph()
{
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
    assert_approx_eq!(tab5_3a[i][2], ph2T_reg3(p, h));
    assert_approx_eq!(tab5_3a[i][3], ph2v_reg3(p, h));
    p = tab5_3b[i][0];
    h = tab5_3b[i][1];
    assert_approx_eq!(tab5_3b[i][2], ph2T_reg3(p, h));
    assert_approx_eq!(tab5_3b[i][3], ph2v_reg3(p, h));
  }
}

#[test] 
fn test_ps()
{
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
    assert_approx_eq!(tab12_3a[i][2], ps2T_reg3(p, s));
    assert_approx_eq!(tab12_3a[i][3], ps2v_reg3(p, s));
    p = tab12_3b[i][0];
    s = tab12_3b[i][1];
    assert_approx_eq!(tab12_3b[i][2], ps2T_reg3(p, s));
    assert_approx_eq!(tab12_3b[i][3], ps2v_reg3(p, s));
  }
}

#[test] 
fn test_hs()
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
    assert_approx_eq!(tab_3a[i][2], hs2p_reg3(h, s));
    h = tab_3b[i][0];
    s = tab_3b[i][1];
    assert_approx_eq!(tab_3b[i][2], hs2p_reg3(h, s));
  }
}

