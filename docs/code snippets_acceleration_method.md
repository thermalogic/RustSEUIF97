# The code snippets of the acceleration methods

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
// --- Module: algo/polynomial_steps.rs ---
// 1. The optimized kernel: Uses loop splitting (steps) and Shared-Power Scaling
#[inline(always)]
pub fn polys_i_j_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for m in 0..steps.len() {
       // Loop splitting for better cache locality or SIMD potential
        for k in steps[m].0..steps[m].1 {
            //Shared-Power Scaling: Compute shared power terms only once
            item = IJn[k].2 * vi.powi(IJn[k].0) * vj.powi(IJn[k].1);
            poly_i += IJn[k].0 as f64 * item;
            poly_j += IJn[k].1 as f64 * item;
        }
    }

    // the base scalingvision)
    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}

// --- Module: r1/region1_gre.rs ---
// 2. Region 1 Wrapper: Handles specific region 1 equations for Gibbs free energy
// and returns the partial derivatives of Gibbs free energy with respect to pi and tau
pub fn polys_i_j_powi_reg1(pi: f64, tau: f64) -> (f64, f64) {
    // profiling-guided loop tiling: define calculation steps explicitly to assist compiler optimization
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

## Performance Comparison: SEUIF97 vs CoolProp IF97

The comparison evaluates SEUIF97 against CoolProp IF97, which employs its own repeated-squaring fast integer power algorithm. 

**Test Environment:** 
- CPU: Intel Core i7-1165G7 @ 2.80GHz
- RAM: 8GB DDR4 
- OS: Windows 11

**Compilation:**
- SEUIF97: Rust 1.96.0, release mode (opt-level = 3, target-cpu = native)
- CoolProp IF97: MSVC 19.50.35719.0 (-O3, -march=native)

Measurements were performed using `clock_t`.

|Case|Region|Input| CoolProp IF97(ns) | SEUIF97(ns) |Speedup|
|---|---|---|---|---|---|
|(p, T) → h| 1 |3.0MPa, 300K |151.1|41.4|3.7x|
|(p, T) → s |1| 3.0MPa, 300K|246.3|46.1|5.3x|
|(p, T) → h |2| 0.0035MPa, 300K |207.9|47.4|4.4x|
|(p, T) → s |2| 0.0035MPa, 300K |377.4|57.0|6.6x|
|(p, T) → h |3| 50MPa, 630K |389.5|89.3|4.4x|
|(p, T) → s |3 |50MPa, 630K |422.1|83.1|5.1x|
|(p, T) → h |5| 0.5MPa, 1500K | 54.2| 9.9| 5.5x|
|(p, T) → s |5| 0.5MPa, 1500K | 74.8 |16.8| 4.5x|

SEUIF97 achieves 3.7–6.6x speedups over CoolProp IF97. These gains stem from the proposed algorithmic optimizations.
