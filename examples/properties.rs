//! The example to get all properties
#![allow(warnings)]
use if97::*;
fn main() {
    let prop_map: [(&str, i32); 19] = [
        // 1 propertry in IF97 Basic Equation
        ("p  Pressure(MPa)", OP),
        ("t  Temperature(°C)", OT),
        ("d  Density(kg/m^3)", OD),
        ("v  Specific Volume(m^3/kg)", OV),
        ("h  Specific enthalpy(kJ/kg)", OH),
        ("s  Specific entropy(kJ/(kg·K))", OS),
        ("u  Specific internal energy(kJ/kg)", OU),
        ("cp Specific isobaric heat capacity(kJ/(kg·K))", OCP),
        ("cv Specific isochoric heat capacity(kJ/(kg·K))", OCV),
        ("w  Speed of sound(m/s)", OW),
        ("x  Steam quality", OX),
        ("r  Region", OR),
        // 2 transport and further propertry
        ("dv Dynamic viscosity(kg/(m·s))", ODV),
        ("kv Kinematic viscosity(m^2/s)", OKV),
        ("tc Thermal conductivity(W/(m.K))", OTC),
        ("td Thermal diffusivity(m^2/s)", OTD),
        ("pr Prandtl number", OPR),
        ("st Surface tension(N/m)", OST),
        ("sdc Static Dielectric Constant", OSDC),
    ];
    let p: f64 = 3.0;
    let t: f64 = 300.0 - 273.15;
    for i in 0..19 {
        let value: f64 = pt(p, t, prop_map[i].1);
        println!("\t {} = {:.6e}", prop_map[i].0, value);
    }
}
