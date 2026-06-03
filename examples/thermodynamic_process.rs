#![allow(warnings)]
/// The example to test the thermodynamic process functions
///
///  * ishd(pi, ti, pe) - Isentropic enthalpy drop
///  * ief(pi, ti, pe, te) - Isentropic efficiency (%)
///
///  Run:
///    cargo run -r --example thermodynamic_process
///
use seuif97::*;

fn main() {
    // Inlet conditions
    let pi: f64 = 16.0;  // MPa
    let ti: f64 = 535.1; // °C

    // Outlet pressure
    let pe: f64 = 5.0;   // MPa

    // Isentropic enthalpy drop
    let delta_h = ishd(pi, ti, pe);
    println!("Inlet:  pi={pi:.1} MPa, ti={ti:.1} °C");
    println!("Outlet: pe={pe:.1} MPa");
    println!("Isentropic enthalpy drop: ishd={delta_h:.3} kJ/kg");

    // Isentropic efficiency
    let te: f64 = 350.0; // °C (actual outlet temperature)
    let efficiency = ief(pi, ti, pe, te);
    println!("\nActual outlet temperature: te={te:.1} °C");
    println!("Isentropic efficiency: ief={efficiency:.2}%");
}
