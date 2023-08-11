#![allow(warnings)]
///  cross test
///     cargo test --test cross_test
use assert_approx_eq::assert_approx_eq;
use if97::*;
#[test]
fn test_cross() {
    let p: f64 = 3.0;
    let t: f64 = 300.0 - 273.15;
    // pt,ph,ps,hs
    let h = pt(p, t, if97::OH);
    let s = pt(p, t, if97::OS);
    //
    assert_approx_eq!(t, ph(p, h, if97::OT), 1.0e-1f64);
    assert_approx_eq!(t, ps(p, s, if97::OT), 1.0e-1f64);

    assert_approx_eq!(t, hs(h, s, if97::OT), 1.0e-1f64);
    assert_approx_eq!(p, hs(h, s, if97::OP), 1.0e-1f64);

    // px,tx
    let p: f64 = 0.1;
    let t_sat: f64 = px(p, 0.0, if97::OT);
    assert_approx_eq!(p, if97::tx(t_sat, 0.0, if97::OP), 1.0e-6f64);

    assert_approx_eq!(px(p, 0.0, if97::OH), tx(t_sat, 0.0, if97::OH), 1.0e-6f64);
    assert_approx_eq!(px(p, 0.0, if97::OS), tx(t_sat, 0.0, if97::OS), 1.0e-6f64);
}
