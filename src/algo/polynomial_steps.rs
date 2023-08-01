//! The multi-step method to enable the compiler optimizations for using powi() within `for` loop
//!
//! The functions compute the polynomial values of the base variable and its derivatives
//!  1. To the polynomial of base variable and its derivatives
//!  2. To the polynomial of base variable and its derivatives recursively
//！# Variables
//! * IJn[(i32,i32,f64)]
//!   * vi - the base of i=IJn[k][0]
//!   * vj - the base of j=IJn[k][1]
//!   * power = n * vi^i * i^j =  IJn[k].2 * vi^ IJn[k].0 * vj^ IJn[k].1
//!

#[inline(always)]
pub fn poly_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> f64 {
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            value += IJn[k].2 * vi.powi(IJn[k].0) * vj.powi(IJn[k].1);
        }
    }
    value
}

#[inline(always)]
pub fn poly_i_powi_steps(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    steps: &[(usize, usize)],
) -> f64 {
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            value += IJn[k].2 * IJn[k].0 as f64 * vi.powi(IJn[k].0 - 1) * vj.powi(IJn[k].1);
        }
    }
    value
}

#[inline(always)]
pub fn poly_ii_powi_steps(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    steps: &[(usize, usize)],
) -> f64 {
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            value += IJn[k].2
                * (IJn[k].0 * (IJn[k].0 - 1)) as f64
                * vi.powi(IJn[k].0 - 2)
                * vj.powi(IJn[k].1);
        }
    }
    value
}

#[inline(always)]
pub fn poly_ij_powi_steps(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    steps: &[(usize, usize)],
) -> f64 {
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            value += IJn[k].2
                * IJn[k].0 as f64
                * vi.powi(IJn[k].0)
                * IJn[k].1 as f64
                * vj.powi(IJn[k].1);
        }
    }
    value / vi / vj
}

#[inline(always)]
pub fn poly_j_powi_steps(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    steps: &[(usize, usize)],
) -> f64 {
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            value += IJn[k].2 * vi.powi(IJn[k].0) * IJn[k].1 as f64 * vj.powi(IJn[k].1 - 1);
        }
    }
    value
}

/// the polynomial of vi and the derivative (∂²f/∂²vj)
/// *  n* vi^i  *j*(j-1)* vi^(j-2)
#[inline(always)]
pub fn poly_jj_powi_steps(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    steps: &[(usize, usize)],
) -> f64 {
    let mut value = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            value += IJn[k].2
                * vi.powi(IJn[k].0)
                * (IJn[k].1 * (IJn[k].1 - 1)) as f64
                * vj.powi(IJn[k].1 - 2);
        }
    }
    value
}

//---------------------- The recursive method to compute the multiple polynomials ---------------------------

/// The recursive method to get the polynomials
///  * the power of vi and vj  
///  * the power of vi and the derivative (∂f/∂vj)
#[inline(always)]
pub fn polys_0_j_powi_steps(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_0: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * vi.powi(IJn[k].0) * vj.powi(IJn[k].1);
            poly_0 += item;
            poly_j += IJn[k].1 as f64 * item;
        }
    }
    poly_j /= vj;
    (poly_0, poly_j)
}

#[inline(always)]
pub fn polys_i_j_powi_steps(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * vi.powi(IJn[k].0) * vj.powi(IJn[k].1);
            poly_i += IJn[k].0 as f64 * item;
            poly_j += IJn[k].1 as f64 * item;
        }
    }
    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}

#[inline(always)]
pub fn polys_i_ij_powi_steps(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_ij: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * IJn[k].0 as f64 * vi.powi(IJn[k].0 - 1) * vj.powi(IJn[k].1);
            poly_i += item;
            poly_ij += IJn[k].1 as f64 * item;
        }
    }
    poly_ij /= vj;
    (poly_i, poly_ij)
}

#[inline(always)]
pub fn polys_i_ii_powi_steps(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_ii: f64 = 0.0;

    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * IJn[k].0 as f64 * vi.powi(IJn[k].0 - 1) * vj.powi(IJn[k].1);
            poly_i += item;
            poly_ii += (IJn[k].0 - 1) as f64 * item;
        }
    }
    poly_ii /= vi;
    (poly_i, poly_ii)
}

#[inline(always)]
pub fn polys_i_ii_ij_jj_powi_steps(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    steps: &[(usize, usize)],
) -> (f64, f64, f64, f64) {
    let mut item: f64 = 0.0;
    let mut i_item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_ii: f64 = 0.0;
    let mut poly_ij: f64 = 0.0;
    let mut poly_jj: f64 = 0.0;

    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * vi.powi(IJn[k].0) * vj.powi(IJn[k].1);
            i_item = IJn[k].0 as f64 * item;
            poly_i += i_item;
            poly_ii += (IJn[k].0 - 1) as f64 * i_item;
            poly_ij += IJn[k].1 as f64 * i_item;
            poly_jj += (IJn[k].1 * (IJn[k].1 - 1)) as f64 * item;
        }
    }

    poly_i /= vi;
    poly_ii /= (vi * vi);
    poly_ij /= (vi * vj);
    poly_jj /= (vj * vj);
    (poly_i, poly_ii, poly_ij, poly_jj)
}
