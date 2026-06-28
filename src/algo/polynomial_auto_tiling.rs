//! The functions compute the polynomial and its derivatives values 
//!    using
//!     * shared-power scaling method 
//!     * auto tiling
//! * IJn[(i32,i32,f64)]
//!   * I - IJn[k].0
//!   * J - IJn[k].1
//!   * n - IJn[k].2
//! * vi - the base of I
//! * vj - the base of J
//!    polynomial = n * vi^I * vj^J =  IJn[k].2 * vi^ IJn[k].0 * vj^ IJn[k].1

use crate::algo::*;

/// Auto tiling strategy:
/// - len < 36: 2 segments, equally distributed
/// - len >= 36: 3 segments, equally distributed
fn auto_steps(len: usize) -> Vec<(usize, usize)> {
    if len < 36 {
        let mid = len / 2;
        vec![(0, mid), (mid, len)]
    } else {
        let third = len / 3;
        let third_2=third+third;
        vec![(0, third), (third, third_2), (third_2, len)]
    }
}

/// Dedicated for four-output functions: 4 segments equally distributed
fn auto_steps_4(len: usize) -> Vec<(usize, usize)> {
    let quarter = len / 4;
    let quarter_2=quarter+quarter;
    let quarter_3=quarter_2+quarter;
    vec![
        (0, quarter),
        (quarter, quarter_2),
        (quarter_2, quarter_3),
        (quarter_3, len)
    ]
}

///  the polynomial:  n*vi^i* vj^j
pub fn poly_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps(len);
       poly_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        poly_powi(vi, vj, &IJn)
    }
}

/// the polynomial of the derivative (∂f/∂vi) and vj  
/// * n * (i-1)*vi^(i-1) * vj^j
pub fn poly_i_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps(len);
       poly_i_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        poly_i_powi(vi, vj, &IJn)
    }
}

/// the polynomial of the derivative (∂²f/∂²vi) and vj
/// * n*i*(i-1)*vi^(i-2) * vj^j
pub fn poly_ii_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps(len);
       poly_ii_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        poly_ii_powi(vi, vj, &IJn)
    }
}

/// the polynomial of vi and  the derivative (∂f/∂vj)
///  n* vi^i  *j* vj^(j-1)
pub fn poly_j_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps(len);
       poly_j_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        poly_j_powi(vi, vj, &IJn)
    }
}

/// the polynomial of vi and the derivative (∂²f/∂²vj)
/// *  n* vi^i  *j*(j-1)* vj^(j-2)
pub fn poly_jj_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps(len);
       poly_jj_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        poly_jj_powi(vi, vj, &IJn)
    }
}

/// the polynomial of the derivative (∂²f/∂vi∂vj) 
/// * n*i*vi^(i-1) *j*vj^(j-1)
pub fn poly_ij_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps(len);
       poly_ij_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        poly_ij_powi(vi, vj, &IJn)
    }
}

/// The shared-power scaling method to get the polynomials
///  * the power of vi and vj  
///  * the power of vi and the derivative (∂f/∂vj)
pub fn polys_0_j_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64) {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps(len);
       polys_0_j_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        polys_0_j_powi(vi, vj, &IJn)
    }
}

/// The shared-power scaling method to get the polynomials
///  * the power of the derivative (∂f/∂vi) and vj
///  * the power of vi and the derivative (∂f/∂vj)
pub fn polys_i_j_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64) {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps(len);
       polys_i_j_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        polys_i_j_powi(vi, vj, &IJn)
    }
}


/// The shared-power scaling method to get the polynomials
///  * the power of the derivative (∂f/∂vi) and vj
///  * the power of the derivative (∂²f/∂vi∂vj)
pub fn polys_i_ij_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64) {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps(len);
       polys_i_ij_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        polys_i_ij_powi(vi, vj, &IJn)
    }
}

/// The shared-power scaling method to get the polynomials
///  * the power of the derivative (∂f/∂vi) and  (∂f/∂²vi)
pub fn polys_i_ii_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64) {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps(len);
       polys_i_ii_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        polys_i_ii_powi(vi, vj, &IJn)
    }
}

/// The shared-power scaling method to compute four derivatives simultaneously
///  * the power of the derivative (∂f/∂vi) and vj
///  * the power of the derivative (∂²f/∂²vi) and vj
///  * the power of the derivative (∂²f/∂vi∂vj)
///  * the power of vi and the derivative (∂²f/∂²vj)
pub fn polys_i_ii_ij_jj_powi_auto_tiling(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64, f64, f64) {
    let len:usize = IJn.len();
    if len >= 29  
    { 
       let steps = auto_steps_4(len);
       polys_i_ii_ij_jj_powi_steps(vi, vj, &IJn, &steps)
    }
    else
    {
        polys_i_ii_ij_jj_powi(vi, vj, &IJn)
    }
}