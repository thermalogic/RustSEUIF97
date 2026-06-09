# Performance Comparison: SEUIF97 vs CoolProp IF97SEUIF97 vs CoolProp IF9

The comparison evaluates SEUIF97 against [CoolProp IF97](https://github.com/CoolProp/IF97), which employs its own repeated-squaring fast integer power algorithm. 

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
