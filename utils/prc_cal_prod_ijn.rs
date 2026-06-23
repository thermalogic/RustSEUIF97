#![allow(warnings)]
///  prc_cal_prod_ijn:
///  Precomputing the product of IJn yields marginal speed gains 
///  while increasing code complexity; this optimization is not adopted.
///     cargo run -r --example prc_cal_prod_ijn
///  
use seuif97::r1::IJn;
use std::fs::File;
use std::io::Write;

fn main() {
    let mut n_I = Vec::new();
    let mut n_I_I1 = Vec::new();
    let mut n_J = Vec::new();
    let mut n_I_J = Vec::new();
    let mut n_J_J1 = Vec::new();

    for e in IJn {
        n_I.push(e.2 * e.0 as f64);
        n_I_I1.push(e.2 * (e.0 * (e.0 - 1)) as f64);
        n_J.push(e.2 * e.1 as f64);
        n_I_J.push(e.2 * e.0 as f64 * e.1 as f64);
        n_J_J1.push(e.2 * (e.1 * (e.1 - 1)) as f64);
    }

    let output = format!(
        "pub const n_I_prod: [f64; 34] = [ {} ];\n\
         pub const n_I_I1_prod: [f64; 34] = [{}; ]\n\
         pub const n_J_prod: [f64; 34] = [ {} ];\n\
         pub const n_I_J_prod: [f64; 34] = [ {} ];\n\
         pub const n_J_J1_prod: [f64; 34] = [ {} ];\n",
        n_I.iter().map(|v| format!("{:.15e}", v)).collect::<Vec<_>>().join(", "),
        n_I_I1.iter().map(|v| format!("{:.15e}", v)).collect::<Vec<_>>().join(", "),
        n_J.iter().map(|v| format!("{:.15e}", v)).collect::<Vec<_>>().join(", "),
        n_I_J.iter().map(|v| format!("{:.15e}", v)).collect::<Vec<_>>().join(", "),
        n_J_J1.iter().map(|v| format!("{:.15e}", v)).collect::<Vec<_>>().join(", "),
    );

    let filename = "prc_cal_prod_ijn_results.txt";
    match File::create(filename) {
        Ok(mut file) => {
            if let Err(e) = file.write_all(output.as_bytes()) {
                eprintln!("Error writing to file: {}", e);
            } else {
                println!("Results saved to: {}", filename);
            }
        }
        Err(e) => {
            eprintln!("Error creating file: {}", e);
        }
    }

    println!("{}", output);
}