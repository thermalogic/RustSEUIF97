//! The functions compute the polynomial values of the base variable and its derivatives
//!  1. To the polynomial of base variable and its derivatives
//!  2. To the polynomial of base variable and its derivatives recursively
//！# Variables
//! * IJn[(i32,i32,f64)]
//!     *  e - the element in IJn[k]
//!     *  n = e.2, i=e.0, j=e.1  
//! * vi - the base of i=IJn[k][0]
//! * vj - the base of j=IJn[k][1]
//!     * power = n * vi^i * i^j = e.2 * vi^e.0 * vj^e.1

///  powi()
///  * n*vi^i *  vj^j
pub fn poly_powi(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * vi.powi(e.0) * vj.powi(e.1);
    }
    value
}

/// the polynomial of the derivative (∂f/∂vi) and vj  
/// * n * (i-1)*vi^(i-1) * vj^j
pub fn poly_i_powi(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * e.0 as f64 * vi.powi(e.0 - 1) * vj.powi(e.1);
    }
    value
}

/// the polynomial of the derivative (∂²f/∂²vi) and vj
/// * n*i*(i-1)*vi^(i-2) * vj^j
pub fn poly_ii(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * (e.0 * (e.0 - 1)) as f64 * vi.powi(e.0 - 2) * vj.powi(e.1);
    }
    value
}

/// the polynomial of the derivative (∂²f/∂²vi) and vj
/// * n*i*(i-1)*vi^(i-2) * vj^j
pub fn poly_ii_powi(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * (e.0 * (e.0 - 1)) as f64 * vi.powi(e.0) * vj.powi(e.1);
    }
    value
}

/// the polynomial of vi and  the derivative (∂f/∂vj)
///  n* vi^i  *j* vj^(j-1)
pub fn poly_j(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * vi.powi(e.0) * e.1 as f64 * vj.powi(e.1 - 1);
    }
    value
}

/// the polynomial of vi and  the derivative (∂f/∂vj)
///  n* vi^i  *j* vj^(j-1)
pub fn poly_j_powi(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * vi.powi(e.0) * e.1 as f64 * vj.powi(e.1 - 1);
    }
    value
}

/// the polynomial of vi and the derivative (∂²f/∂²vj)
/// *  n* vi^i  *j*(j-1)* vi^(j-2)
pub fn poly_jj(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value = 0.0;
    for e in IJn {
        value += e.2 * vi.powi(e.0) * (e.1 * (e.1 - 1)) as f64 * vj.powi(e.1 - 2);
    }
    value
}

/// the polynomial of vi and the derivative (∂²f/∂²vj)
/// *  n* vi^i  *j*(j-1)* vi^(j-2)
pub fn poly_jj_powi(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value = 0.0;
    for e in IJn {
        value += e.2 * vi.powi(e.0) * (e.1 * (e.1 - 1)) as f64 * vj.powi(e.1 - 2);
    }
    value
}

/// the polynomial of the derivative (∂f/∂vi) and (∂f/∂vj)
///   * n* i*vi^(i-1) *j** vi^(j-1)
pub fn poly_ij(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * e.0 as f64 * vi.powi(e.0 - 1) * e.1 as f64 * vj.powi(e.1 - 1);
    }
    value
}

pub fn poly_ij_powi(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> f64 {
    let mut value: f64 = 0.0;
    for e in IJn {
        value += e.2 * e.0 as f64 * vi.powi(e.0 - 1) * e.1 as f64 * vj.powi(e.1 - 1);
    }
    value
}

/// The recursive method to get the polynomials
///  * the power of vi and vj  
///  * the power of vi and the derivative (∂f/∂vj)
pub fn polys_0_j(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_0: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for e in IJn {
        item = e.2 * vi.powi(e.0) * vj.powi(e.1);
        poly_0 += item;
        poly_j += e.1 as f64 * item;
    }
    poly_j /= vj;
    (poly_0, poly_j)
}

/// The recursive method to get the polynomials
///  * the power of vi and vj  
///  * the power of vi and the derivative (∂f/∂vj)
pub fn polys_0_j_powi(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_0: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for e in IJn {
        item = e.2 * vi.powi(e.0) * vj.powi(e.1);
        poly_0 += item;
        poly_j += e.1 as f64 * item;
    }
    poly_j /= vj;
    (poly_0, poly_j)
}

/// The recursive method to get the polynomials
///  * the power of the derivative (∂f/∂vi) and vj
///  * the power of vi and the derivative (∂f/∂vj)
pub fn polys_i_j(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for e in IJn {
        item = e.2 * vi.powi(e.0) * vj.powi(e.1);
        poly_i += e.0 as f64 * item;
        poly_j += e.1 as f64 * item;
    }
    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}

/// The recursive method to get the polynomials
///  * the power of the derivative (∂f/∂vi) and vj
///  * the power of vi and the derivative (∂f/∂vj)
pub fn polys_i_j_powi(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for e in IJn {
        item = e.2 * vi.powi(e.0) * vj.powi(e.1);
        poly_i += e.0 as f64 * item;
        poly_j += e.1 as f64 * item;
    }
    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}

/// recursion the powers of the derivatives of (∂f/∂vi),(∂²f/∂²vi),(∂²f/∂vi∂vj) and (∂²f/∂²vj)
/// The recursive method to get the polynomials
///  * the power of the derivative (∂f/∂vi) and vj
///  * the power of the derivative (∂²f/∂²vi) and vj
///  * the power of  the derivative (∂f/∂vi) and (∂f/∂vi)
///  * the power of the derivative vi and (∂²f/∂²vj)the polynomiapoly
pub fn polys_i_ii_ij_jj(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64, f64, f64) {
    let mut item: f64 = 0.0;
    let mut i_item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_ii: f64 = 0.0;
    let mut poly_ij: f64 = 0.0;
    let mut poly_jj: f64 = 0.0;

    for e in IJn {
        item = e.2 * vi.powi(e.0) * vj.powi(e.1);

        poly_jj += (e.1 * (e.1 - 1)) as f64 * item;

        i_item = e.0 as f64 * item;
        poly_i += i_item;
        poly_ii += (e.0 - 1) as f64 * i_item;
        poly_ij += e.1 as f64 * i_item;
    }

    poly_i /= vi;
    poly_ii /= vi * vi;
    poly_ij /= vi * vj;
    poly_jj /= vj * vj;
    (poly_i, poly_ii, poly_ij, poly_jj)
}

pub fn polys_i_ii_ij_jj_powi(vi: f64, vj: f64, IJn: &[(i32, i32, f64)]) -> (f64, f64, f64, f64) {
    let mut item: f64 = 0.0;
    let mut i_item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_ii: f64 = 0.0;
    let mut poly_ij: f64 = 0.0;
    let mut poly_jj: f64 = 0.0;

    for e in IJn {
        item = e.2 * vi.powi(e.0) * vj.powi(e.1);

        poly_jj += (e.1 * (e.1 - 1)) as f64 * item;

        i_item = e.0 as f64 * item;
        poly_i += i_item;
        poly_ii += (e.0 - 1) as f64 * i_item;
        poly_ij += e.1 as f64 * i_item;
    }

    poly_i /= vi;
    poly_ii /= vi * vi;
    poly_ij /= vi * vj;
    poly_jj /= vj * vj;
    (poly_i, poly_ii, poly_ij, poly_jj)
}
