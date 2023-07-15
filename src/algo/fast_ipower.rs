//!  The fast integer power using the shortest addition chains,n[0,58]
  
fn possac(x: f64, n: i32) -> f64 {
    let mut x2: f64;
    let mut x3: f64;
    let mut x4: f64;
    let mut x5: f64;
    let mut x6: f64;
    let mut x7: f64;
    let mut x8: f64;
    let mut x9: f64;

    let mut x10: f64;
    let mut x11: f64;
    let mut x12: f64;
    let mut x13: f64;
    let mut x14: f64;
    let mut x15: f64;
    let mut x16: f64;
    let mut x17: f64;
    let mut x18: f64;
    let mut x19: f64;

    let mut x20: f64;
    let mut x21: f64;
    let mut x22: f64;
    let mut x23: f64;
    let mut x24: f64;
    let mut x25: f64;
    let mut x26: f64;
    let mut x27: f64;
    let mut x28: f64;
    let mut x29: f64;

    let mut x32: f64;
    let mut x33: f64;
    let mut x36: f64;
    let mut x37: f64;

    let mut x46: f64;
    let mut x48: f64;
    let mut x49: f64;
    let mut x54: f64;

    match n {
        0 => 1.0,
        1 => x,
        2 => x * x,
        3 => {
            x2 = x * x;
            return x2 * x;
        }
        4 => {
            x2 = x * x;
            return x2 * x2;
        }
        5 => {
            x2 = x * x;
            x3 = x2 * x;
            return x2 * x3;
        }
        6 => {
            x2 = x * x;
            x3 = x2 * x;
            return x3 * x3;
        }
        7 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            return x5 * x2;
        }
        8 => {
            x2 = x * x;
            x4 = x2 * x2;
            return x4 * x4;
        }
        9 => {
            x2 = x * x;
            x4 = x2 * x2;
            return x * x4 * x4;
        }
        10 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            return x5 * x5;
        }
        11 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            return x5 * x5 * x;
        }
        12 => {
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            return x6 * x6;
        }
        13 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            return x8 * x4 * x;
        }
        14 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            x7 = x5 * x2;
            return x7 * x7;
        }
        15 => {
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            return x6 * x6 * x3;
        }
        16 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            return x8 * x8;
        }
        17 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            return x8 * x8 * x;
        }
        18 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x16 = x8 * x8;
            return x16 * x2;
        }
        19 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x16 = x8 * x8;
            return x16 * x2 * x;
        }
        20 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            x10 = x5 * x5;
            return x10 * x10;
        }
        21 => {
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            return x12 * x6 * x3;
        }
        22 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            x10 = x5 * x5;
            return x10 * x10 * x2;
        }
        23 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            x10 = x5 * x5;
            return x10 * x10 * x3;
        }
        24 => {
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            return x12 * x12;
        }
        25 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x9 = x8 * x;
            return x9 * x8 * x8;
        }
        26 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x9 = x8 * x;
            x13 = x9 * x4;
            return x13 * x13;
        }
        27 => {
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            return x12 * x12 * x3;
        }
        28 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            x7 = x5 * x2;
            x14 = x7 * x7;
            return x14 * x14;
        }
        29 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x9 = x8 * x;
            x17 = x8 * x9;
            return x17 * x8 * x4;
        }

        30 => {
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            x15 = x12 * x3;
            return x15 * x15;
        }
        31 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            x7 = x5 * x2;
            x14 = x7 * x7;
            return x14 * x14 * x3;
        }
        32 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x16 = x8 * x8;
            return x16 * x16;
        }
        33 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x16 = x8 * x8;
            return x16 * x16 * x;
        }
        34 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x9 = x8 * x;
            x17 = x8 * x9;
            return x17 * x17;
        }
        35 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x9 = x8 * x;
            x13 = x4 * x9;
            return x13 * x13 * x9;
        }
        36 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x18 = x8 * x8 * x2;
            return x18 * x18;
        }
        37 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x18 = x8 * x8 * x2;
            return x18 * x18 * x;
        }
        38 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x16 = x8 * x8;
            x18 = x16 * x2;
            x19 = x18 * x;
            return x19 * x19;
        }

        39 => {
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            x24 = x12 * x12;
            return x24 * x12 * x3;
        }

        40 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            x10 = x5 * x5;
            x20 = x10 * x10;
            return x20 * x20;
        }
        41 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x3 * x2;
            x10 = x5 * x5;
            x20 = x10 * x10;
            return x20 * x20 * x;
        }
        42 => {
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            x21 = x12 * x3 * x6;
            return x21 * x21;
        }
        43 => {
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x9 = x8 * x;
            x17 = x8 * x9;
            return x17 * x17 * x9;
        }
        44 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x2 * x3;
            x11 = x5 * x5 * x;
            x22 = x11 * x11;
            return x22 * x22;
        }
        45 => {
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            x15 = x12 * x3;
            return x15 * x15 * x15;
        }
        46 => {
            x2 = x * x;
            x3 = x2 * x;
            x5 = x2 * x3;
            x10 = x5 * x5;
            x20 = x10 * x10;
            x23 = x20 * x3;
            return x23 * x23;
        }
        47=>{
            x2 = x * x;
            x3 = x2 * x;
            x5 = x2 * x3;
            x10 = x5 * x5;
            x20 = x10 * x10;
            x23 = x20 * x3;
            x46 = x23 * x23;
            return x46 * x;
        }  
        48=>{
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            x24 = x12 * x12;
            return x24 * x24;
        }
        49=>{
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x16 = x8 * x8;
            x32 = x16 * x16;
            x33 = x32 * x;
            return x33 * x16;
        }
        50=>{
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x9 = x8 * x;
            x17 = x8 * x9;
            x25 = x17 * x8;
            return x25 * x25
        }
        51=>{
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            x24 = x12 * x12;
            x48 = x24 * x24;
            return x48 * x3;
        }
        52=>{
            x2 = x * x;
            x4 = x2 * x2;
            x5 = x4 * x;
            x8 = x4 * x4;
            x13 = x8 * x5;
            x26 = x13 * x13;
            return x26 * x26;
        }
        53=>{	x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x16 = x8 * x8;
            x18 = x16 * x2;
            x36 = x18 * x18;
            x37 = x36 * x;
            return x37 * x16;
        }
        54=>{
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            x24 = x12 * x12;
            x27 = x24 * x3;
            return x27 * x27;
        }
        55=>{
            x2 = x * x;
            x3 = x2 * x;
            x6 = x3 * x3;
            x12 = x6 * x6;
            x24 = x12 * x12;
            x27 = x24 * x3;
            x54 = x27 * x27;
            return x54 * x;
        }
        56=>{
            x2 = x * x;
            x3 = x2 * x;
            x4 = x2 * x2;
            x7 = x4 * x3;
            x14 = x7 * x7;
            x28 = x14 * x14;
            return x28 * x28;
        }
        57=>{	x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x16 = x8 * x8;
            x32 = x16 * x16;
            x33 = x32 * x;
            x49 = x33 * x16;
            return x49 * x8;
        }
        58=>{
            x2 = x * x;
            x4 = x2 * x2;
            x8 = x4 * x4;
            x9 = x8 * x;
            x17 = x8 * x9;
            x25 = x17 * x8;
            x29 = x25 * x4;
            return x29 * x29;
        }
      _ => x.powi(n)
    }
}

/// The fast integer power using the shortest addition chains,n[0,58]
pub fn sac_pow(x: f64, n: i32) -> f64 {
    if n >= 0 {
        return possac(x, n);
    } else {
        return 1.0 / possac(x, -n);
    }
}
