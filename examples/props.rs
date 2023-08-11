/// The example to test all properties
/// 
///   Run: 
///    cargo run -r --example props
/// 
#![allow(warnings)]
use if97::*;

const prop_map: [(&str, &str, i32); 36] = [
    ("Pressure", "MPa", OP),
    ("Temperature", "°C", OT),
    ("Density", "kg/m³", OD),
    ("Specific Volume", "m³/kg", OV),
    ("Specific enthalpy", "kJ/kg", OH),
    ("Specific entropy", "kJ/(kg·K)", OS),
    ("Specific exergy", "kJ/kg", OE),
    ("Specific internal energy", "kJ/kg", OU),
    ("Specific isobaric heat capacity", "kJ/(kg·K)", OCP),
    ("Specific isochoric heat capacity", "kJ/(kg·K)", OCV),
    ("Speed of sound", "m/s", OW),
    ("Isentropic exponent", "-", OKS),
    ("Specific Helmholtz free energy", "kJ/kg", OF),
    ("Specific Gibbs free energy", "kJ/kg", OG),
    ("Compressibility factor", "-", OZ),
    ("Steam quality", "-", OX),
    ("Region", "-", OR),
    ("Isobaric volume expansion coefficient", "1/K", OEC),
    ("Isothermal compressibility", "1/MPa", OKT),
    ("Partial derivative (∂V/∂T)p", "m³/(kg·K)", ODVDT),
    ("Partial derivative (∂V/∂P)T", "m³/(kg·MPa)", ODVDP),
    ("Partial derivative (∂P/∂T)v", "MPa/K", ODPDT),
    ("Isothermal Joule-Thomson coefficient", "kJ/(kg·MPa)", OIJTC),
    ("Joule-Thomson coefficient", "K/MPa", OJTC),
    ("Dynamic viscosity", "kg/(m·s)", ODV),
    ("Kinematic viscosity", "m²/s", OKV),
    ("Thermal conductivity", "W/(m.K)", OTC),
    ("Thermal diffusivity", "m²/s", OTD),
    ("Prandtl number", "-", OPR),
    ("Surface tension", "N/m", OST),
    ("Static Dielectric Constant", "-", OSDC),
    ("Isochoric pressure coefficient", "-", OPC),
    ("Isothermal stress coefficient", "kg/m³", OBETAP),
    ("Fugacity coefficient", "-", OFI),
    ("Fugacity", "Mpa ", OFU),
    ("Relative pressure coefficient", "1/K", OALFAP),
];

/// v1 :the first input value
/// v2 :the second input value
/// func: the function to get properties
fn get_props(v1: f64, v2: f64, func: fn(f64, f64, i32) -> f64) {
    println!("\t┌{}┐", "─".repeat(78));
    println!("\t| {:2} |  {:^40} | {:^13}| {:^12} |", "No", "Propertry", "Unit", "Value");
    println!("\t|{}|", "─".repeat(78));
    for e in prop_map {
        let value: f64 = func(v1, v2, e.2);
        println!("\t| {:2} |  {:40} | {:^12} |{value:>12.5e}  |", e.2, e.0, e.1);
    }
    println!("\t└{}┘ ", "─".repeat(78));
}

fn main() {
    println!("{:^50}", "Region 1");
    let mut p: f64 = 3.0;
    let mut t: f64 = 300.0 - 273.15;
    get_props(p, t, pt);

    println!("\n {:^50}", "Region 2");
    p = 24.0;
    t = 540.0;
    get_props(p, t, pt);

    println!("\n {:^50}", "Region 3");
    p = 50.0;
    t = 630.0 - 273.15;
    get_props(p, t, pt);

    println!("\n {:^50}", "Region 5");
    p = 30.0;
    t = 1500.0 - 273.15;
    get_props(p, t, pt);
}
