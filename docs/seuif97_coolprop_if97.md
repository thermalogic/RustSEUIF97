# Performance Comparison: SEUIF97 vs CoolProp IF97

The comparison evaluates SEUIF97 against [CoolProp IF97](https://github.com/CoolProp/IF97), which employs its own repeated-squaring fast integer power algorithm. 

SEUIF97 achieves 1.9–6.9x speedups over CoolProp IF97. These gains stem from the proposed algorithmic optimizations.

**Note:** SEUIF97 is accessed via C FFI (foreign function interface). The FFI call overhead (~1.5 ns/call) was measured using a dummy function and subtracted from the reported times to ensure a fair comparison with the natively compiled CoolProp IF97.

**Test Environment:** 
- CPU: Intel Core i7-1165G7 @ 2.80GHz
- RAM: 8GB DDR4 
- OS: Windows 11

**Compilation:**
- SEUIF97: Rust 1.96.0, release mode (opt-level = 3, target-cpu = native)
- CoolProp IF97: MSVC 19.50.35719.0 (-O3, -march=native)

Measurements were performed using `std::chrono::high_resolution_clock`.

|Case|Reg. | Input | CoolProp IF97(ns) | SEUIF97(ns) |Speedup|
|---|:---:|---:|:------:|:---:|:---:|
|(p, T) → h| 1 |3.0 MPa, 300 K |134.7|38.4|3.5x|
|(p, T) → s |1| 3.0 MPa, 300 K|255.8|42.8|6.0x|
|(p, T) → h |2| 0.0035 MPa, 300 K |188.5|49.6|3.8x|
|(p, T) → s |2| 0.0035 MPa, 300 K |374.2|54.1|6.9x|
|(T, v) → h |3| 650 K, 0.002 m³/kg |213.0|29.6|7.2x|
|(T, v) → s |3 |650 K, 0.002 m³/kg |213.1|28.0|7.6x|
|(p, T) → h |5| 0.5 MPa, 1500 K | 30.9| 9.1| 3.4x|
|(p, T) → s |5| 0.5 MPa, 1500 K | 56.7 |15.1| 3.8x|
