# The code snippets of the acceleration methods

The code snippets demonstrate the acceleration methods to calculate the specific internal energy in region 1, illustrating the flow from the optimized kernel to the final physical property calculation:

## Key Acceleration Methods

* Loop Tiling Method: Unleashes the full power of compiler optimizations, surpassing the performance of the single loop.

* Recurrence Method for Multi-Polynomial Evaluation: By leveraging the relationship between polynomials and their derivatives, only a single polynomial needs to be computed directly. The remaining values are derived via multiplication or division by the base. This approach eliminates redundant calculations and significantly improves performance.

## The IAPWS-IF97 Equations

The basic equation for this region 1 is a fundamental equation for the specific **Gibbs free energy** $g$. 

$$\frac{g(p,T)}{RT} = \gamma(\pi,\tau) = \sum_{i=1}^{34} n_i (7.1-\pi)^{I_i} (\tau-1.222)^{J_i}$$

where $\pi = p/p^{*}\quad p^{*}=16.53$MPa

To derive the **specific internal energy**

$$u = g - T \left( \frac{\partial g}{\partial T} \right)_p - p \left( \frac{\partial g}{\partial p} \right)_T$$

$$\frac{u(\pi, \tau)}{RT} = \tau \gamma_{\tau} - \pi \gamma_{\pi}$$

## Implementation Details

```rust
// --- Module: algo/polynomial_steps.rs ---
// 1. The optimized kernel: Uses loop splitting (steps) and aggressive inlining
#[inline(always)]
pub fn polys_i_j_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    // Loop splitting for better cache locality or SIMD potential
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * vi.powi(IJn[k].0) * vj.powi(IJn[k].1);
            poly_i += IJn[k].0 as f64 * item;
            poly_j += IJn[k].1 as f64 * item;
        }
    }

    // Multi-polynomial evaluation derived via base scaling (multiplication/division)
    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}

// --- Module: r1/region1_gre.rs ---
// 2. Region 1 Wrapper: Handles specific coordinate transformations (pi, tau)
pub fn polys_i_j_powi_reg1(pi: f64, tau: f64) -> (f64, f64) {
    // Define calculation steps explicitly to assist compiler optimization
    let steps: [(usize, usize); 3] = [(0, 16), (16, 26), (26, 34)];
    let (d_pi, d_tau) = polys_i_j_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps);
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
