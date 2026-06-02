# SEUIF97: Rust vs C Version Comparison

The Rust version of SEUIF97 is a major upgrade over the original C implementation, delivering significant improvements in performance, functionality, and ecosystem support.

## Summary Comparison Table

| Feature | C Version | Rust Version |
|---------|-----------|--------------|
| **Calculation Speed** | Baseline | **3× faster** |
| **Supported Properties** | 30 properties | **36 properties** (+6 new) |
| **Package Distribution** | PyPI only | **Crates.io, PyPI, npm** |
| **Supported OS (Pre-built)** | Windows, Linux | **Windows, Linux, macOS** |
| **Universal Functions** | ✓ | ✓ |
| **Direct Property Functions** | ✗ | **✓** (new) |
| **Thermodynamic Process Calculation** | ✓ | ✗ |

## Key Improvements

### 1. 3× Faster Calculation Speed

The Rust implementation achieves a 3× speedup over the C version through:
- Loop tiling optimization that unleashes full compiler auto-vectorization
- Recurrence method for multi-polynomial evaluation, eliminating redundant calculations
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

### 4. Cross-Platform Shared Library & Package Support

| Platform | C Version | Rust Version |
|----------|-----------|--------------|
| **Windows** | ✓ Pre-built | ✓ Pre-built |
| **Linux** | ✓ Pre-built | ✓ Pre-built |
| **macOS** | ✗ Compile from source | ✓ Pre-built |

The Rust version provides pre-compiled shared libraries for all three major operating systems (Windows, Linux, macOS) via GitHub Actions CI/CD. The C version requires macOS users to compile from source themselves.

### 5. Dual API Design: Universal + Direct Functions

The Rust version introduces two complementary calculation approaches:

| API Type | Description | Example |
|----------|-------------|---------|
| **Universal Functions** | Single function with property ID parameter | `pt(p, t, OH)` |
| **Direct Property Functions** | Dedicated function for each property | `pt2h(p, t)` |

The C version only supports universal functions. The Rust version's direct property functions provide:
- Cleaner, more readable code
- Compile-time type safety
- Better IDE autocomplete support

### 6. Thermodynamic Process Calculation

The C version supports thermodynamic process calculations (e.g., isentropic expansion, polytropic processes). This feature is not yet implemented in the Rust version but is planned for future releases.
