# The code snippets of the acceleration methods

[![DOI](https://img.shields.io/badge/DOI-10.20944/preprints202606.0793.v1)](https://doi.org/10.20944/preprints202606.0793.v1)

The code snippets demonstrate the acceleration methods to calculate the specific internal energy $u$ in region 1, illustrating the flow from the optimized kernel to the final physical property calculation:

## The IAPWS-IF97 Equations

The basic equation for this region 1 is a fundamental equation for the specific **Gibbs free energy** $g$. 

$$\frac{g(p,T)}{RT} = \gamma(\pi,\tau) = \sum_{i=1}^{34} n_i (7.1-\pi)^{I_i} (\tau-1.222)^{J_i}$$

where 

$\pi = p/p^{*}$

$\tau = T^{*}/T$

To derive the **specific internal energy**

$$u = g - T \left( \frac{\partial g}{\partial T} \right)_p - p \left( \frac{\partial g}{\partial p} \right)_T$$

$$\frac{u(\pi, \tau)}{RT} = \tau \gamma_{\tau} - \pi \gamma_{\pi}$$

## Code snippets

```rust
// --- Module: algo/polynomial_tile.rs ---
// 1. The optimized kernel: Uses loop tiling and Shared-Power Scaling
#[inline(always)]
pub fn polys_i_j_powi_tile(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], tiles: &[(usize, usize)]) -> (f64, f64) {
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for &(start, end) in tiles {
       for &(I, J, n) in &IJn[start..end] {
            let item = n * vi.powi(I) * vj.powi(J);
            poly_i += I as f64 * item;
            poly_j += J as f64 * item;
        }
    }
    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}


// --- Module: r1/region1_gre.rs ---
// 2. Region 1 Wrapper: Handles specific region 1 equations for Gibbs free energy
// and returns the partial derivatives of Gibbs free energy with respect to pi and tau
pub fn polys_i_j_powi_reg1(pi: f64, tau: f64) -> (f64, f64) {
    // profiling-guided loop tiling: define calculation tiles explicitly to assist compiler optimization
    let tiles: [(usize, usize); 3] = [(0, 16), (16, 26), (26, 34)];
    let (d_pi, d_tau) = polys_i_j_powi_tile(7.1 - pi, tau - 1.222, &IJn, &tiles);
    (-d_pi, d_tau)
}

// --- Module: r1/region_pT.rs ---
// 3. API: Calculates specific internal energy
pub fn pT2u_reg1(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r1pstar;
    let tau: f64 = r1Tstar / T;
    let (d_pi, d_tau) = polys_i_j_powi_reg1(pi, tau);
    RGAS_WATER * T * (tau * d_tau - pi * d_pi)
}
```
