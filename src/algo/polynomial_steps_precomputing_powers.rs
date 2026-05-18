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
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn poly_powi_steps_precomputed(
    IJn: &[(i32, i32, f64)], 
    x_powers: &[f64], 
    y_powers: &[f64], 
    steps: &[(usize, usize)]
) -> f64 {
    let mut value: f64 = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                let x = x_powers.get_unchecked(k);
                let y = y_powers.get_unchecked(k);
                value += coeff.2 * x * y;
            }
        }
    }
    value
}

/// 使用预计算的幂值计算对 x 的一阶导数
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn poly_i_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> f64 {
    let vi_inv = 1.0 / vi;
    let mut value: f64 = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                if coeff.0 != 0 {
                    let x = x_powers.get_unchecked(k);
                    let y = y_powers.get_unchecked(k);
                    value += coeff.2 * coeff.0 as f64 * x * vi_inv * y;
                }
            }
        }
    }
    value
}

/// 使用预计算的幂值计算对 x 的二阶导数
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn poly_ii_powi_steps_precomputed(
    vi: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> f64 {
    let vi_inv = 1.0 / vi;
    let vi_inv2 = vi_inv * vi_inv;
    let mut value: f64 = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                if coeff.0 >= 2 {
                    let x = x_powers.get_unchecked(k);
                    let y = y_powers.get_unchecked(k);
                    value += coeff.2 * (coeff.0 * (coeff.0 - 1)) as f64 * x * vi_inv2 * y;
                }
            }
        }
    }
    value
}

/// 使用预计算的幂值计算对 x 和 y 的混合导数
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn poly_ij_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> f64 {
    let vi_inv = 1.0 / vi;
    let vj_inv = 1.0 / vj;
    let mut value: f64 = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                if coeff.0 != 0 && coeff.1 != 0 {
                    let x = x_powers.get_unchecked(k);
                    let y = y_powers.get_unchecked(k);
                    value += coeff.2 * coeff.0 as f64 * coeff.1 as f64 * x * vi_inv * y * vj_inv;
                }
            }
        }
    }
    value
}

/// 使用预计算的幂值计算对 y 的一阶导数
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn poly_j_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> f64 {
    let vj_inv = 1.0 / vj;
    let mut value: f64 = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                if coeff.1 != 0 {
                    let x = x_powers.get_unchecked(k);
                    let y = y_powers.get_unchecked(k);
                    value += coeff.2 * x * coeff.1 as f64 * y * vj_inv;
                }
            }
        }
    }
    value
}

/// 使用预计算的幂值计算对 y 的二阶导数
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn poly_jj_powi_steps_precomputed(
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> f64 {
    let vj_inv = 1.0 / vj;
    let vj_inv2 = vj_inv * vj_inv;
    let mut value = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                if coeff.1 >= 2 {
                    let x = x_powers.get_unchecked(k);
                    let y = y_powers.get_unchecked(k);
                    value += coeff.2 * x * (coeff.1 * (coeff.1 - 1)) as f64 * y * vj_inv2;
                }
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
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn polys_0_j_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let vj_inv = 1.0 / vj;
    let mut item: f64 = 0.0;
    let mut poly_0: f64 = 0.0;
    let mut poly_j: f64 = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                let x = x_powers.get_unchecked(k);
                let y = y_powers.get_unchecked(k);
                item = coeff.2 * x * y;
                poly_0 += item;
                if coeff.1 != 0 {
                    poly_j += coeff.1 as f64 * item;
                }
            }
        }
    }
    poly_j *= vj_inv;
    (poly_0, poly_j)
}

/// 使用预计算的幂值计算对 x 和 y 的一阶导数
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn polys_i_j_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let vi_inv = 1.0 / vi;
    let vj_inv = 1.0 / vj;
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                let x = x_powers.get_unchecked(k);
                let y = y_powers.get_unchecked(k);
                item = coeff.2 * x * y;
                if coeff.0 != 0 {
                    poly_i += coeff.0 as f64 * item;
                }
                if coeff.1 != 0 {
                    poly_j += coeff.1 as f64 * item;
                }
            }
        }
    }
    poly_i *= vi_inv;
    poly_j *= vj_inv;
    (poly_i, poly_j)
}

/// 使用预计算的幂值计算对 x 的一阶导数和混合导数
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn polys_i_ij_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let vi_inv = 1.0 / vi;
    let vj_inv = 1.0 / vj;
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_ij: f64 = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                if coeff.0 != 0 {
                    let x = x_powers.get_unchecked(k);
                    let y = y_powers.get_unchecked(k);
                    item = coeff.2 * coeff.0 as f64 * x * vi_inv * y;
                    poly_i += item;
                    if coeff.1 != 0 {
                        poly_ij += coeff.1 as f64 * item;
                    }
                }
            }
        }
    }
    poly_ij *= vj_inv;
    (poly_i, poly_ij)
}

/// 使用预计算的幂值计算对 x 的一阶和二阶导数
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn polys_i_ii_powi_steps_precomputed(
    vi: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> (f64, f64) {
    let vi_inv = 1.0 / vi;
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_ii: f64 = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                if coeff.0 != 0 {
                    let x = x_powers.get_unchecked(k);
                    let y = y_powers.get_unchecked(k);
                    item = coeff.2 * coeff.0 as f64 * x * vi_inv * y;
                    poly_i += item;
                    if coeff.0 >= 2 {
                        poly_ii += (coeff.0 - 1) as f64 * item;
                    }
                }
            }
        }
    }
    poly_ii *= vi_inv;
    (poly_i, poly_ii)
}

/// 使用预计算的幂值计算所有一阶和二阶导数
/// 优化：使用 unsafe 消除数组边界检查
#[inline(always)]
pub fn polys_i_ii_ij_jj_powi_steps_precomputed(
    vi: f64,
    vj: f64,
    IJn: &[(i32, i32, f64)],
    x_powers: &[f64],
    y_powers: &[f64],
    steps: &[(usize, usize)],
) -> (f64, f64, f64, f64) {
    let vi_inv = 1.0 / vi;
    let vi_inv2 = vi_inv * vi_inv;
    let vj_inv = 1.0 / vj;
    let vj_inv2 = vj_inv * vj_inv;
    let mut item: f64 = 0.0;
    let mut i_item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_ii: f64 = 0.0;
    let mut poly_ij: f64 = 0.0;
    let mut poly_jj: f64 = 0.0;
    unsafe {
        for m in 0..steps.len() {
            let range = steps.get_unchecked(m);
            for k in range.0..range.1 {
                let coeff = IJn.get_unchecked(k);
                let x = x_powers.get_unchecked(k);
                let y = y_powers.get_unchecked(k);
                item = coeff.2 * x * y;
                if coeff.0 != 0 {
                    i_item = coeff.0 as f64 * item;
                    poly_i += i_item;
                    if coeff.0 >= 2 {
                        poly_ii += (coeff.0 - 1) as f64 * i_item;
                    }
                    if coeff.1 != 0 {
                        poly_ij += coeff.1 as f64 * i_item;
                    }
                }
                if coeff.1 >= 2 {
                    poly_jj += (coeff.1 * (coeff.1 - 1)) as f64 * item;
                }
            }
        }
    }
    poly_i *= vi_inv;
    poly_ii *= vi_inv2;
    poly_ij *= vi_inv * vj_inv;
    poly_jj *= vj_inv2;
    (poly_i, poly_ii, poly_ij, poly_jj)
}