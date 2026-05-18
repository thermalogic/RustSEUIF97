//! The multi-step method to enable the compiler optimizations for using powi() within `for` loop
//!
//! The functions compute the polynomial values of the base variable and its derivatives
//!  1. To the polynomial of base variable and its derivatives
//!  2. To the polynomial of base variable and its derivatives recursively
//! # Variables
//! * IJn[(i32,i32,f64)]
//!   * vi - the base of i=IJn[k][0]
//!   * vj - the base of j=IJn[k][1]
//!
//!

/// 使用预计算的幂值计算多项式
/// 参数：
/// - IJn: 系数数组
/// - x_powers: 预计算的 x^i 数组
/// - y_powers: 预计算的 y^j 数组
/// - steps: 计算步骤
#[inline(always)]
pub fn poly_powi_steps_precomputed(
    IJn: &[(i32, i32, f64)], 
    x_powers: &[f64], 
    y_powers: &[f64], 
    steps: &[(usize, usize)]
) -> f64 {
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            value += IJn[k].2 * x_powers[k] * y_powers[k];
        }
    }
    value
}

/// 使用预计算的幂值计算对 x 的一阶导数
#[inline(always)]
pub fn poly_i_powi_steps_precomputed(
    vi: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> f64 {
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            if IJn[k].0 != 0 {
                value += IJn[k].2 * IJn[k].0 as f64 * (x_powers[k] / vi) * y_powers[k];
            }
        }
    }
    value
}

/// 使用预计算的幂值计算对 x 的二阶导数
#[inline(always)]
pub fn poly_ii_powi_steps_precomputed(
    vi: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> f64 {
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            if IJn[k].0 >= 2 {
                value += IJn[k].2 * (IJn[k].0 * (IJn[k].0 - 1)) as f64 * (x_powers[k] / (vi * vi)) * y_powers[k];
            }
        }
    }
    value
}

/// 使用预计算的幂值计算对 x 和 y 的混合导数
#[inline(always)]
pub fn poly_ij_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> f64 {
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            if IJn[k].0 != 0 && IJn[k].1 != 0 {
                value += IJn[k].2 * IJn[k].0 as f64 * IJn[k].1 as f64 * (x_powers[k] / vi) * (y_powers[k] / vj);
            }
        }
    }
    value
}

/// 使用预计算的幂值计算对 y 的一阶导数
#[inline(always)]
pub fn poly_j_powi_steps_precomputed(
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> f64 {
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            if IJn[k].1 != 0 {
                value += IJn[k].2 * x_powers[k] * IJn[k].1 as f64 * (y_powers[k] / vj);
            }
        }
    }
    value
}

/// 使用预计算的幂值计算对 y 的二阶导数
#[inline(always)]
pub fn poly_jj_powi_steps_precomputed(
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> f64 {
    let mut value = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            if IJn[k].1 >= 2 {
                value += IJn[k].2 * x_powers[k] * (IJn[k].1 * (IJn[k].1 - 1)) as f64 * (y_powers[k] / (vj * vj));
            }
        }
    }
    value
}

//---------------------- The recursive method to compute the multiple polynomials ---------------------------

/// The recursive method to get the polynomials
///  * the power of vi and vj  
///  * the power of vi and the derivative (∂f/∂vj)

/// 使用预计算的幂值计算多项式和对 y 的一阶导数
#[inline(always)]
pub fn polys_0_j_powi_steps_precomputed(
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_0: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * x_powers[k] * y_powers[k];
            poly_0 += item;
            if IJn[k].1 != 0 {
                poly_j += IJn[k].1 as f64 * item;
            }
        }
    }
    poly_j /= vj;
    (poly_0, poly_j)
}

/// 使用预计算的幂值计算对 x 和 y 的一阶导数
#[inline(always)]
pub fn polys_i_j_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * x_powers[k] * y_powers[k];
            if IJn[k].0 != 0 {
                poly_i += IJn[k].0 as f64 * item;
            }
            if IJn[k].1 != 0 {
                poly_j += IJn[k].1 as f64 * item;
            }
        }
    }
    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}

/// 使用预计算的幂值计算对 x 的一阶导数和混合导数
#[inline(always)]
pub fn polys_i_ij_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_ij: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            if IJn[k].0 != 0 {
                item = IJn[k].2 * IJn[k].0 as f64 * (x_powers[k] / vi) * y_powers[k];
                poly_i += item;
                if IJn[k].1 != 0 {
                    poly_ij += IJn[k].1 as f64 * item;
                }
            }
        }
    }
    poly_ij /= vj;
    (poly_i, poly_ij)
}

/// 使用预计算的幂值计算对 x 的一阶和二阶导数
#[inline(always)]
pub fn polys_i_ii_powi_steps_precomputed(
    vi: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_ii: f64 = 0.0;

    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            if IJn[k].0 != 0 {
                item = IJn[k].2 * IJn[k].0 as f64 * (x_powers[k] / vi) * y_powers[k];
                poly_i += item;
                if IJn[k].0 >= 2 {
                    poly_ii += (IJn[k].0 - 1) as f64 * item;
                }
            }
        }
    }
    poly_ii /= vi;
    (poly_i, poly_ii)
}

/// 使用预计算的幂值计算所有一阶和二阶导数
#[inline(always)]
pub fn polys_i_ii_ij_jj_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
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
            item = IJn[k].2 * x_powers[k] * y_powers[k];
            if IJn[k].0 != 0 {
                i_item = IJn[k].0 as f64 * item;
                poly_i += i_item;
                if IJn[k].0 >= 2 {
                    poly_ii += (IJn[k].0 - 1) as f64 * i_item;
                }
                if IJn[k].1 != 0 {
                    poly_ij += IJn[k].1 as f64 * i_item;
                }
            }
            if IJn[k].1 >= 2 {
                poly_jj += (IJn[k].1 * (IJn[k].1 - 1)) as f64 * item;
            }
        }
    }

    poly_i /= vi;
    poly_ii /= (vi * vi);
    poly_ij /= (vi * vj);
    poly_jj /= (vj * vj);
    (poly_i, poly_ii, poly_ij, poly_jj)
}