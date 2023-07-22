//! The example to test all properties
#![allow(warnings)]
use if97::*;

const prop_map: [(&str, &str, &str, i32); 35] = [
    ("p", "Pressure", "MPa", OP),
    ("t", "Temperature", "°C", OT),
    ("d", "Density", "kg/m³", OD),
    ("v", "Specific Volume", "m³/kg", OV),
    ("h", "Specific enthalpy", "kJ/kg", OH),
    ("s", "Specific entropy", "kJ/(kg·K)", OS),
    ("e", "Specific exergy", "kJ/kg", OE),
    ("u", "Specific internal energy", "kJ/kg", OU),
    ("cp", "Specific isobaric heat capacity", "kJ/(kg·K)", OCP),
    ("cv", "Specific isochoric heat capacity", "kJ/(kg·K)", OCV),
    ("w", "Speed of sound", "m/s", OW),
    ("ks", "Isentropic exponent", "-", OKS),
    ("f", "Specific Helmholtz free energy", "kJ/kg", OF),
    ("g", "Specific Gibbs free energy", "kJ/kg", OG),
    ("z", "Compressibility factor", "-", OZ),
    ("x", "Steam quality", "-", OX),
    ("r", "Region", "-", OR),
    ("ec", "Isobaric volume expansion coefficient", "1/K", OEC),
    ("kt", "Isothermal compressibility", "1/MPa", OKT),
    ("dvctcp", "Partial derivative (∂V/∂T)p", "m³/(kg·K)", ODVDT),
    (
        "dvdpct",
        "Partial derivative (∂V/∂P)T",
        "m³/(kg·MPa)",
        ODVDP,
    ),
    ("dpdtcv", "Partial derivative (∂P/∂T)v", "MPa/K", ODPDT),
    (
        "iJTC",
        "Isothermal Joule-Thomson coefficient",
        "kJ/(kg·MPa)",
        OIJTC,
    ),
    ("joule", "Joule-Thomson coefficient", "K/MPa", OJTC),
    ("dv", "Dynamic viscosity", "kg/(m·s)", ODV),
    ("kv", "Kinematic viscosity", "m²/s", OKV),
    ("tc", "Thermal conductivity", "W/(m.K)", OTC),
    ("td", "Thermal diffusivity", "m²/s", OTD),
    ("pr", "Prandtl number", "-", OPR),
    ("st", "Surface tension", "N/m", OST),
    ("sdc", "Static Dielectric Constant", "-", OSDC),
    ("pc", "Isochoric pressure coefficient", "-", OPC),
    ("betap", "Isothermal stress coefficient", "kg/m³", OBETAP),
    ("fi", "Fugacity coefficient", "-", OFI),
    ("fu", "Fugacity", "Mpa ", OFU),
];

/// v1 :the first input value
/// v2 :the second input value
/// func: the function to get properties
fn get_props(v1: f64, v2: f64, func: fn(f64, f64, i32) -> f64) {
    println!("\t┌{}┐", "─".repeat(78));
    println!(
        "\t| {:2} |  {:^40} | {:^13}| {:^12} |",
        "No", "Propertry", "Unit", "Value"
    );
    println!("\t|{}|", "─".repeat(78));
    for e in prop_map {
        let value: f64 = func(v1, v2, e.3);
        println!(
            "\t| {:2} |  {:40} | {:^12} |{value:>12.5e}  |",
            e.3, e.1, e.2
        );
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
