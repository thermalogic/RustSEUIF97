//! sum the power of function and it's derivative on IJn[(i32,i32,f64)]
//! * e - the element in IJn[(i32,i32,f64)]
//! * vi - the item powered by i(index is 0) in IJn[(i32,i32,f64) :vi^e.0
//! * vj - the item powered by j(index is 1) in IJn[(i32,i32,f64) :vj^e.1

use crate::algo::sac_pow;

/// sum the powers of function
pub fn sum_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * sac_pow(vi, e.0) * sac_pow(vj, e.1);
    }
    value
}

/// sum the power of the first derivative in vi
pub fn sum_di_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * e.0 as f64 * sac_pow(vi, e.0 - 1) * sac_pow(vj, e.1);
    }
    value
}

/// sum the power of the second derivative in vi
pub fn sum_dii_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * (e.0 * (e.0 - 1)) as f64 * sac_pow(vi, e.0 - 2) * sac_pow(vj, e.1);
    }
    value
}

/// sum the power of the first derivative in vj
pub fn sum_dj_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * sac_pow(vi, e.0) * e.1 as f64 * sac_pow(vj, e.1 - 1);
    }
    value
}

/// sum the power of the second derivative in vj
pub fn sum_djj_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value = 0.0;
    for e in IJn {
        value += e.2 * sac_pow(vi, e.0) * (e.1 * (e.1 - 1)) as f64 * sac_pow(vj, e.1 - 2);
    }
    value
}

/// sum the power of the second derivative in vi and vj
pub fn sum_dij_power(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * e.0 as f64 * sac_pow(vi, e.0 - 1) * e.1 as f64 * sac_pow(vj, e.1 - 1);
    }
    value
}

///  Fast recursion algorithm of 0 and j derivatives
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

///  Fast recursion algorithm of i and j derivatives
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

/// Fast recursion algorithm of i,ii,ij,jj derivatives
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
