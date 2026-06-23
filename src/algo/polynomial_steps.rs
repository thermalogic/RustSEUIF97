//! The functions compute the polynomial and its derivatives values 
//!    using loop tiling and shared-power scaling method 
//! 
//! * IJn[(i32,i32,f64)]
//!   * I - IJn[k].0
//!   * J - IJn[k].1
//!   * n - IJn[k].2
//!       let (I, J, n) = IJn[k];
//! * vi - the base of I
//! * vj - the base of J
//!    polynomial = n * vi^I * vj^J =  IJn[k].2 * vi^ IJn[k].0 * vj^ IJn[k].1
//!

///  the polynomial:  n*vi^i* vj^j
#[inline(always)]
pub fn poly_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> f64 {
    let mut value: f64 = 0.0;
    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            value += n * vi.powi(I) * vj.powi(J);  
        }
    }
    value
}

/// the polynomial of the derivative (∂f/∂vi)   
/// * n * i*vi^(i-1) * vj^j
#[inline(always)]
pub fn poly_i_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> f64 {
    let mut value: f64 = 0.0;
    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            value += n * I as f64* vi.powi(I - 1) * vj.powi(J);
        }
    }
    value
}

/// the polynomial of the derivative (∂²f/∂²vi) 
/// * n*i*(i-1)*vi^(i-2) * vj^j
#[inline(always)]
pub fn poly_ii_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> f64 {
    let mut value: f64 = 0.0;
    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            value += n * (I*(I-1))  as f64 *vi.powi(I - 2) * vj.powi(J);
        }
    }
    value
}

/// the polynomial of the derivative (∂²f/∂vi∂vj) 
/// * n*i*vi^(i-1) *j*vj^(j-1)
#[inline(always)]
pub fn poly_ij_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> f64 {
    let mut value: f64 = 0.0;
    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            value += n * I as f64* vi.powi(I - 1) *J as f64* vj.powi(J - 1);
        }
    }
    value
}

/// the polynomial of derivative (∂f/∂vj)
///  n* vi^i  *j* vj^(j-1)
#[inline(always)]
pub fn poly_j_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> f64 {
    let mut value: f64 = 0.0;
    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            value += n * vi.powi(I) *J as f64*vj.powi(J - 1);
        }
    }
    value
}

/// the polynomial of the derivative (∂²f/∂²vj)
/// *  n* vi^i  *j*(j-1)* vj^(j-2)
#[inline(always)]
pub fn poly_jj_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> f64 {
    let mut value = 0.0;
    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            value += n * vi.powi(I) * (J*(J-1)) as f64*vj.powi(J - 2);
        }
    }
    value
}

//---------------------- The recursive method to compute the multiple polynomials ---------------------------

/// The recursive method to get the polynomials
///  * the power of vi and vj  
///  * the power of vi and the derivative (∂f/∂vj)
#[inline(always)]
pub fn polys_0_j_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> (f64, f64) {
    let mut poly_0: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            let item = n * vi.powi(I) * vj.powi(J);
            poly_0 += item;
            poly_j += J as f64 * item;
        }
    }
    poly_j /= vj;
    (poly_0, poly_j)
}

#[inline(always)]
pub fn polys_i_j_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> (f64, f64) {
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            let item = n * vi.powi(I) * vj.powi(J);
            poly_i += I as f64 * item;
            poly_j += J as f64 * item;
        }
    }
    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}

#[inline(always)]
pub fn polys_i_ij_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> (f64, f64) {
    let mut poly_i: f64 = 0.0;
    let mut poly_ij: f64 = 0.0;
    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            let item = n *I as f64 * vi.powi(I - 1) * vj.powi(J);
            poly_i += item;
            poly_ij += J as f64 * item;
        }
    }
    poly_ij /= vj;
    (poly_i, poly_ij)
}

#[inline(always)]
pub fn polys_i_ii_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> (f64, f64) {
    let mut poly_i: f64 = 0.0;
    let mut poly_ii: f64 = 0.0;
    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            let item = n * I as f64 * vi.powi(I - 1) * vj.powi(J);
            poly_i += item;
            poly_ii += (I - 1) as f64 * item;
        }
    }
    poly_ii /= vi;
    (poly_i, poly_ii)
}

#[inline(always)]
pub fn polys_i_ii_ij_jj_powi_steps(
    vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)],
) -> (f64, f64, f64, f64) {
    let mut poly_i: f64 = 0.0;
    let mut poly_ii: f64 = 0.0;
    let mut poly_ij: f64 = 0.0;
    let mut poly_jj: f64 = 0.0;
    for &(start, end) in steps {
        for k in start..end {
            let (I, J, n) = IJn[k];
            let item = n * vi.powi(I) * vj.powi(J);
            let i_item = I as f64 * item;
            poly_i += i_item;
            poly_ii += (I - 1) as f64 * i_item;
            poly_ij += J as f64 * i_item;
            poly_jj += (J * (J - 1)) as f64 * item;
        }
    }

    poly_i /= vi;
    poly_ii /= (vi * vi);
    poly_ij /= (vi * vj);
    poly_jj /= (vj * vj);
    (poly_i, poly_ii, poly_ij, poly_jj)
}
