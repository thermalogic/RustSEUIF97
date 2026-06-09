# Performance Comparison: SEUIF97 vs CoolProp IF97

The comparison evaluates SEUIF97 against [CoolProp IF97](https://github.com/CoolProp/IF97), which employs its own repeated-squaring fast integer power algorithm. 

SEUIF97 achieves 1.9–6.7x speedups over CoolProp IF97. These gains stem from the proposed algorithmic optimizations.

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
|(p, T) → h| 1 |3.0Mpa, 300K |136.1|42.6|3.2x|
|(p, T) → s |1| 3.0Mpa, 300K|255.5|45.8|5.6x|
|(p, T) → h |2| 0.0035Mpa, 300K |188.0|48.8|3.9x|
|(p, T) → s |2| 0.0035Mpa, 300K |381.9|56.8|6.7x|
|(T, v) → h |3| 630K,0.002m^3/kg |203.6|106.7|1.9x|
|(T, v) → s |3 |630K,0.002m^3/kg |212.8|107.4|2.0x|
|(p, T) → h |5| 0.5Mpa, 1500K | 32.0| 10.9| 2.9x|
|(p, T) → s |5| 0.5Mpa, 1500K | 56.9 |17.4| 3.3x|
