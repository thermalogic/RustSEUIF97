# SEUIF97: Rust vs C Version Comparison

The Rust version of SEUIF97 is a major upgrade over the original C implementation, delivering significant improvements in performance, functionality, and ecosystem support.

## Summary Comparison Table

| Feature | C Version | Rust Version |
|---------|-----------|--------------|
| **Calculation Speed** | Baseline | **~2× speedup** |
| **Supported Properties** | 30 properties | **36 properties** (+6 new) |
| **Package Distribution** | PyPI only | **Crates.io, PyPI, npm** |
| **Universal Functions** | ✓ | ✓ |
| **Direct Property Functions** | ✗ | **✓** (new) |

## Key Improvements

### 1. ~2× Speedup in Calculation

The Rust implementation achieves a ~2× speedup over the C version through:
- Native `powi()` function: After compiler optimization, integer power calculations are extremely fast
- Loop tiling optimization that unleashes full compiler auto-vectorization
- Zero-cost abstractions and compile-time optimizations (LTO, fat codegen)

### 2. Extended Property Support: 30 → 36 Properties

The Rust version adds 6 new thermodynamic properties:
- Static Dielectric Constant (ε)
- Isochoric Pressure Coefficient (β)
- Isothermal Stress Coefficient (βp)
- Fugacity Coefficient (fi)
- Fugacity (f*)
- Relative Pressure Coefficient (αp)

### 3. Multi-Ecosystem Package Distribution

| Ecosystem | C Version | Rust Version |
|-----------|-----------|--------------|
| Rust (Crates.io) | ✗ | ✓ |
| Python (PyPI) | ✓ | ✓ |
| JavaScript/TypeScript (npm) | ✗ | ✓ |

### 4. Dual API Design: Universal + Direct Functions

The Rust version introduces two complementary calculation approaches:

| API Type | Description | Example |
|----------|-------------|---------|
| **Universal Functions** | Single function with property ID parameter | `pt(p, t, OH)` |
| **Direct Property Functions** | Dedicated function for each property | `pt2h(p, t)` |

The C version only supports universal functions. The Rust version's direct property functions provide a more convenient way to calculate commonly used properties.

### 5. Thermodynamic Process Calculation

The C version supports thermodynamic process calculations (e.g., isentropic expansion, polytropic processes). This feature is not yet implemented in the Rust version but is planned for future releases.
