//! The functions compute the polynomial values of the base variable and its derivatives
//!  1. To sum the one product of the powers of base variable and its derivatives
//!  2. To sum the multiple products of the powers of base variable and its derivatives recursively
//！# Variables
//! * IJn[(i32,i32,f64)]
//!     *  e - the element in IJn[k]
//!     *  n = e.2, i=e.0, j=e.1  
//! * vi - the base of i=IJn[k][0]
//! * vj - the base of j=IJn[k][1]
//!     * power = n * vi^i * i^j = e.2 * vi^e.0 * vj^e.1

use crate::algo::sac_pow;


/// using powi() to sum the product of the powers of vi and vj
///  * n*vi^i *  vj^j
pub fn poly_powi(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * vi.powi(e.0) * vj.powi(e.1);
    }
    value
}

/// sum the product of the powers of vi and vj
///  * n*vi^i *  vj^j
pub fn sum_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * sac_pow(vi, e.0) * sac_pow(vj, e.1);
    }
    value
}


/// sum the product of the power of the derivative (∂f/∂vi) and vj  
/// * n * (i-1)*vi^(i-1) * vj^j
pub fn sum_di_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * e.0 as f64 * sac_pow(vi, e.0 - 1) * sac_pow(vj, e.1);
    }
    value
}

/// sum the product of the power of the derivative (∂²f/∂²vi) and vj
/// * n*i*(i-1)*vi^(i-2) * vj^j
pub fn sum_dii_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * (e.0 * (e.0 - 1)) as f64 * sac_pow(vi, e.0 - 2) * sac_pow(vj, e.1);
    }
    value
}

/// sum the product of the power of vi and  the derivative (∂f/∂vj)
///  n* vi^i  *j* vj^(j-1)
pub fn sum_dj_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * sac_pow(vi, e.0) * e.1 as f64 * sac_pow(vj, e.1 - 1);
    }
    value
}

/// sum the product of the power of vi and the derivative (∂²f/∂²vj)
/// *  n* vi^i  *j*(j-1)* vi^(j-2)
pub fn sum_djj_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value = 0.0;
    for e in IJn {
        value += e.2 * sac_pow(vi, e.0) * (e.1 * (e.1 - 1)) as f64 * sac_pow(vj, e.1 - 2);
    }
    value
}

/// sum the product of the power of the derivative (∂f/∂vi) and (∂f/∂vj)
///   * n* i*vi^(i-1) *j** vi^(j-1)
pub fn sum_dij_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * e.0 as f64 * sac_pow(vi, e.0 - 1) * e.1 as f64 * sac_pow(vj, e.1 - 1);
    }
    value
}

/// The recursive method to sum the product of
///  * the power of vi and vj  
///  * the power of vi and the derivative (∂f/∂vj)
pub fn sum_0_j_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut sum_0: f64 = 0.0;
    let mut sum_j: f64 = 0.0;

    for e in IJn {
        item = e.2 * sac_pow(vi, e.0) * sac_pow(vj, e.1);
        sum_0 += item;
        sum_j += e.1 as f64 * item;
    }
    sum_j /= vj;
    (sum_0, sum_j)
}

/// The recursive method to sum the product of
///  * the power of the derivative (∂f/∂vi) and vj
///  * the power of vi and the derivative (∂f/∂vj)
pub fn sum_i_j_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut sum_i: f64 = 0.0;
    let mut sum_j: f64 = 0.0;

    for e in IJn {
        item = e.2 * sac_pow(vi, e.0) * sac_pow(vj, e.1);
        sum_i += e.0 as f64 * item;
        sum_j += e.1 as f64 * item;
    }
    sum_i /= vi;
    sum_j /= vj;
    (sum_i, sum_j)
}

/// recursion the powers of the derivatives of (∂f/∂vi),(∂²f/∂²vi),(∂²f/∂vi∂vj) and (∂²f/∂²vj)
/// The recursive method to sum the product of
///  * the power of the derivative (∂f/∂vi) and vj
///  * the power of the derivative (∂²f/∂²vi) and vj
///  * the power of the derivative (∂f/∂vi) and (∂f/∂vi)
///  * the power of the derivative vi and (∂²f/∂²vj)
pub fn sum_i_ii_ij_jj_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64, f64, f64) {
    let mut item: f64 = 0.0;
    let mut i_item: f64 = 0.0;
    let mut sum_i: f64 = 0.0;
    let mut sum_ii: f64 = 0.0;
    let mut sum_ij: f64 = 0.0;
    let mut sum_jj: f64 = 0.0;

    for e in IJn {
        item = e.2 * sac_pow(vi, e.0) * sac_pow(vj, e.1);

        sum_jj += (e.1 * (e.1 - 1)) as f64 * item;

        i_item = e.0 as f64 * item;
        sum_i += i_item;
        sum_ii += (e.0 - 1) as f64 * i_item;
        sum_ij += e.1 as f64 * i_item;
    }

    sum_i /= vi;
    sum_ii /= vi * vi;
    sum_ij /= vi * vj;
    sum_jj /= vj * vj;
    (sum_i, sum_ii, sum_ij, sum_jj)
}
