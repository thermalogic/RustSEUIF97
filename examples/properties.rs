//! The example to get all properties
#![allow(warnings)]
fn main() {
    let prop_map: [(&str, i32); 19] = [
        // 1 propertry in IF97 Basic Equation
        ("p  Pressure(MPa)", if97::OP),
        ("t  Temperature(°C)", if97::OT),
        ("d  Density(kg/m^3)", if97::OD),
        ("v  Specific Volume(m^3/kg)", if97::OV),
        ("h  Specific enthalpy(kJ/kg)", if97::OH),
        ("s  Specific entropy(kJ/(kg·K))", if97::OS),
        ("u  Specific internal energy(kJ/kg)", if97::OU),
        ("cp Specific isobaric heat capacity(kJ/(kg·K))", if97::OCP),
        ("cv Specific isochoric heat capacity(kJ/(kg·K))", if97::OCV),
        ("w  Speed of sound(m/s)", if97::OW),
        ("x  Steam quality", if97::OX),
        ("r  Steam quality", if97::OR),
        // 2 transport and further propertry
        ("dv Dynamic viscosity(kg/(m·s))", if97::ODV),
        ("kv Kinematic viscosity(m^2/s)", if97::OKV),
        ("tc Thermal conductivity(W/(m.K))", if97::OTC),
        ("td Thermal diffusivity(m^2/s)", if97::OTD),
        ("pr Prandtl number", if97::OPR),
        ("st Surface tension(N/m)", if97::OST),
        ("sdc Static Dielectric Constant", if97::OSDC),
    ];
    let p: f64 = 3.0;
    let t: f64 = 300.0 - 273.15;
    for i in 0..19 {
        let value: f64 = if97::pt(p, t, prop_map[i].1);
        println!("\t {} = {:.6e}", prop_map[i].0, value);
    }
}
