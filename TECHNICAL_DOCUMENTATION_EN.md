# SEUIF97 Technical Documentation

## 1. Project Overview

SEUIF97 is a high-speed IAPWS-IF97 water and steam property calculation library implemented in Rust. It is specifically designed for computation-intensive tasks such as non-stationary process simulation, on-line process monitoring, and optimization.

### 1.1 Core Performance Advantages

- **Performance Improvement**: Achieves **5x to 20x speedup** compared to direct implementations using Rust standard library's `powi()` for Region 1, 2, and 3 basic equations
- **Significantly outperforms** various approximate equations and algorithms for fast water and steam property calculations

### 1.2 Supported Input Parameter Pairs

```txt
(p,t) (p,h) (p,s) (p,v)
(t,h) (t,s) (t,v)
(p,x) (t,x) (h,x) (s,x)
(h,s)
```

### 1.3 Calculatable Properties

Supports **36 thermodynamic, transport, and derived properties**.

## 2. Technical Architecture

### 2.1 Project Structure

```
RustSEUIF97/
├── src/
│   ├── algo/              # Core algorithm modules
│   │   ├── polynomial.rs          # Polynomial calculation
│   │   ├── polynomial_steps.rs    # Step-by-step polynomial (acceleration core)
│   │   └── root.rs                # Root-finding algorithms
│   ├── common/            # Common components
│   │   ├── boundaries.rs          # Region boundary determination
│   │   ├── constant.rs            # Physical constants
│   │   ├── property_id.rs         # Property ID definitions
│   │   ├── property_pairs.rs      # Property pair calculations
│   │   ├── region.rs              # Region determination logic
│   │   └── transport_further.rs   # Transport property calculations
│   ├── r1/                # Region 1 (Liquid water region)
│   ├── r2/                # Region 2 (Superheated vapor region)
│   ├── r3/                # Region 3 (Near critical point region)
│   ├── r4/                # Region 4 (Saturation region)
│   ├── r5/                # Region 5 (High temperature region)
│   ├── cdecl_c_if97.rs    # C interface (cdecl)
│   ├── stdcall_c_if97.rs  # C interface (stdcall)
│   ├── python_if97.rs     # Python bindings
│   ├── rust_if97.rs       # Rust API
│   └── lib.rs             # Library entry point
├── demo_using_lib/        # Multi-language example code
├── dynamic_lib/           # Pre-compiled dynamic libraries
├── benches/              # Performance benchmarks
├── examples/              # Rust examples
└── tests/                 # Test suite
```

### 2.2 Module Responsibilities

| Module | Responsibility | Key Files |
|--------|----------------|-----------|
| **algo** | Core mathematical algorithm implementation | polynomial_steps.rs, root.rs |
| **common** | Common constants, boundary determination, property definitions | constant.rs, boundaries.rs, region.rs |
| **r1-r5** | Property calculations for each thermodynamic region | regionX_pT.rs, regionX_ph_ps_hs.rs, etc. |
| **cdecl_c_if97** | C interface (cdecl calling convention) | cdecl_c_if97.rs |
| **stdcall_c_if97** | C interface (stdcall calling convention) | stdcall_c_if97.rs |
| **python_if97** | Python bindings | python_if97.rs |
| **rust_if97** | Rust high-level API | rust_if97.rs |

### 2.3 Thermodynamic Region Modules

The project implements the 5 thermodynamic regions of the IAPWS-IF97 standard. The structure of each region module is as follows:

#### Region 1 - Liquid Water Region

**Temperature Range**: 273.15K ≤ T ≤ 623.15K
**Pressure Range**: p ≥ saturation pressure (up to 100MPa)

| File | Responsibility |
|------|----------------|
| `region1.rs` | Region 1 main entry, property dispatch |
| `region1_pT.rs` | Basic equation: (p,T) → v,u,s,h,cp,cv,w |
| `region1_T_phps.rs` | Temperature calculation: (T,p/h/s) → p/h/s |
| `region1_p_hs.rs` | Backward equation: (p,h/s) → T |
| `region1_gfe.rs` | Gibbs free energy equation implementation |
| `region1_pT_ext.rs` | Extended property calculations |
| `region1_pair_ext.rs` | Extended input pair implementation |

#### Region 2 - Superheated Vapor Region

**Temperature Range**: 273.15K ≤ T ≤ 1073.15K
**Pressure Range**: p < saturation pressure (up to 100MPa)

| File | Responsibility |
|------|----------------|
| `region2.rs` | Region 2 main entry, property dispatch |
| `region2_pT.rs` | Basic equation: (p,T) → v,u,s,h,cp,cv,w |
| `region2_T_ph.rs` | Backward equation: (p,h) → T |
| `region2_T_ps.rs` | Backward equation: (p,s) → T |
| `region2_p_hs.rs` | Backward equation: (h,s) → p,T |
| `region2_gfe.rs` | Gibbs free energy equation implementation |
| `region2_pT_ext.rs` | Extended property calculations |
| `region2_pair_ext.rs` | Extended input pair implementation |

#### Region 3 - Near Critical Point Region

**Temperature Range**: 623.15K < T ≤ 863.15K
**Pressure Range**: 16.529MPa ≤ p ≤ 100MPa

| File | Responsibility |
|------|----------------|
| `region3.rs` | Region 3 main entry, property dispatch |
| `region3_Td.rs` | Basic equation: (T,d) → p,h,u,s,cp,cv,w |
| `region3_v_pT.rs` | Backward equation: (p,T) → v |
| `region3_v_subregion_pT.rs` | Subregion subdivision calculation |
| `region3_Tv_phps.rs` | Temperature calculation: (T,v/p/h/s) → p/h/s |
| `region3_p_hs.rs` | Backward equation: (p,h/s) → T,v |
| `region3_hfe.rs` | Helmholtz free energy equation implementation |
| `region3_Td_ext.rs` | Extended property calculations |
| `region3_pair_ext.rs` | Extended input pair implementation |

#### Region 4 - Saturation Region

**Temperature Range**: 273.15K ≤ T ≤ 647.096K (Critical temperature)
**Pressure Range**: 0.000611MPa ≤ p ≤ 22.064MPa (Critical pressure)

| File | Responsibility |
|------|----------------|
| `region4.rs` | Region 4 main entry, saturation property calculation |
| `region4_sat_pT.rs` | Saturation line equation: p=f(T), T=f(p) |
| `region4_pTx.rs` | Saturation region properties: (p,T,x) → h,s,v |
| `region4_T_hs.rs` | Saturation region backward equation |
| `region4_pair_ext.rs` | Extended input pair implementation |

#### Region 5 - High Temperature Region

**Temperature Range**: 1073.15K < T ≤ 2273.15K
**Pressure Range**: p ≤ 50MPa

| File | Responsibility |
|------|----------------|
| `region5.rs` | Region 5 main entry, property dispatch |
| `region5_pT.rs` | Basic equation: (p,T) → v,u,s,h,cp,cv,w |
| `region5_ph_ps_hs.rs` | Backward equation implementation |
| `region5_gfe.rs` | Gibbs free energy equation implementation |
| `region5_pT_ext.rs` | Extended property calculations |
| `region5_pair_ext.rs` | Extended input pair implementation |

---

## 3. Core Algorithms

### 3.1 Acceleration Techniques

SEUIF97 employs two core acceleration techniques:

**Design Principles**:
1. **Loop Tiling**: Splits a single loop into multiple steps, improving cache locality and SIMD vectorization potential
2. **Recursive Evaluation**: Directly derives derivative values by dividing by base values vi/vj, computing polynomial values and their derivatives simultaneously, avoiding redundant calculations

#### 3.1.1 Loop Tiling Method

Splits a single loop into multiple small steps, fully leveraging compiler optimization capabilities and surpassing single-loop performance.

```rust
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

    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}
```

#### 3.1.2 Recurrence Method for Multi-Polynomial Evaluation

By utilizing the relationship between polynomials and their derivatives, only a single polynomial needs to be computed directly. The remaining values are derived via multiplication or division by the base, eliminating redundant calculations and significantly improving computational performance.

##### 3.1.2.1 IAPWS-IF97 Basic Equations

The basic equation for Region 1 is based on the **Gibbs free energy**:

$$\frac{g(p,T)}{RT} = \gamma(\pi,\tau) = \sum_{i=1}^{34} n_i (7.1-\pi)^{I_i} (\tau-1.222)^{J_i}$$

where:
- $\pi = p/p^{*}$
- $\tau = T^{*}/T$

Derivation of specific internal energy:

$$u = g - T \left( \frac{\partial g}{\partial T} \right)_p - p \left( \frac{\partial g}{\partial p} \right)_T$$

$$\frac{u(\pi, \tau)}{RT} = \tau \gamma_{\tau} - \pi \gamma_{\pi}$$

##### 3.1.2.2 Implementation Code

```rust
// --- Module: algo/polynomial_steps.rs ---
// 1. Optimized kernel: Uses loop splitting (steps) and aggressive inlining
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

// --- Module: r1/region1_gfe.rs ---
// 2. Region 1 Wrapper: Handles specific coordinate transformations (pi, tau)
pub fn polys_i_j_powi_reg1(pi: f64, tau: f64) -> (f64, f64) {
    // Define calculation steps explicitly to assist compiler optimization
    let steps: [(usize, usize); 3] = [(0, 16), (16, 26), (26, 34)];
    let (d_pi, d_tau) = polys_i_j_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps);
    (-d_pi, d_tau)
}

// --- Module: r1/region1_pT.rs ---
// 3. API: Calculates specific internal energy
pub fn pT2u_reg1(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r1pstar;
    let tau: f64 = r1Tstar / T;
    let (d_pi, d_tau) = polys_i_j_powi_reg1(pi, tau);
    RGAS_WATER * T * (tau * d_tau - pi * d_pi)
}
```

**Design Principles**:
1. **Loop Splitting**: Splits a single loop into multiple steps, improving cache locality and SIMD vectorization potential
2. **Recursive Evaluation**: Computes polynomial values and their derivatives simultaneously in a single traversal, avoiding redundant calculations
3. **Base Scaling**: Directly derives derivative values by dividing by base values vi/vj, rather than recalculating

### 3.2 Region Determination Logic

Region determination is a critical step in the calculation flow, selecting different determination strategies based on input parameter pairs.

**(p, T) Determination Flow**:
1. **Parameter Boundary Validation**: Checks if pressure and temperature are within valid ranges
2. **Saturation Line Detection**: Determines if located in the saturation region (Region 4)
3. **Region Determination**: Determines the specific region based on temperature and pressure ranges

| Region | Range | State |
|--------|-------|-------|
| **Region 1** | 273.15K ≤ T ≤ 623.15K, p ≥ saturation pressure | Liquid water |
| **Region 2** | T ≥ 273.15K, p < saturation pressure | Superheated vapor |
| **Region 3** | 623.15K < T ≤ 863.15K, high pressure region | Near critical point |
| **Region 4** | 273.15K ≤ T ≤ 647.096K | Saturated two-phase region |
| **Region 5** | 1073.15K < T ≤ 2273.15K, p ≤ 50MPa | High temperature region |

```rust
pub fn pT_sub_region(p: f64, T: f64) -> i32 {
    // 1. Boundary check
    if p < P_MIN || p > 100.0 { return INVALID_P; }
    if T < 273.15 || T > 2273.15 { return INVALID_T; }

    // 2. Saturation line check
    if T >= 273.15 && T < TC_WATER {
        let ps: f64 = p_saturation(T);
        if (p - ps).abs() / ps < psatTol { return 4; }
    }

    // 3. Region determination
    if T >= 273.15 && T <= 623.15 {
        if p >= p_saturation(T) && p <= 100.0 { return 1; }
        if p < p_saturation(T) && p > P_MIN { return 2; }
    }
    // ... other region determination logic
}
```

## 4. API Interface Design

### 4.1 Rust API

#### 4.1.1 Function Signature

```rust
pub fn pt<R>(p: f64, t: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
```

**Parameter Description**:
- `p`: Pressure (MPa)
- `t`: Temperature (°C)
- `o_id_reg`: Property ID (optional with region specification)

#### 4.1.2 Input Parameter Pair Functions

| Function | Input Parameters | Description |
|----------|------------------|--------------|
| `pt(p, t, o_id)` | Pressure, Temperature | Most commonly used input method |
| `ph(p, h, o_id)` | Pressure, Enthalpy | Applicable when enthalpy is known |
| `ps(p, s, o_id)` | Pressure, Entropy | Applicable when entropy is known |
| `pv(p, v, o_id)` | Pressure, Specific Volume | Extended input pair |
| `th(t, h, o_id)` | Temperature, Enthalpy | Extended input pair |
| `ts(t, s, o_id)` | Temperature, Entropy | Extended input pair |
| `tv(t, v, o_id)` | Temperature, Specific Volume | Extended input pair |
| `hs(h, s, o_id)` | Enthalpy, Entropy | Common thermodynamic process input |
| `px(p, x, o_id)` | Pressure, Dryness | Wet steam region only |
| `tx(t, x, o_id)` | Temperature, Dryness | Wet steam region only |
| `hx(h, x, o_id)` | Enthalpy, Dryness | Wet steam region only |
| `sx(s, x, o_id)` | Entropy, Dryness | Wet steam region only |

#### 4.1.3 Property ID Constants

| Property | Symbol | o_id | Unit |
|----------|--------|------|------|
| Pressure | p | OP(0) | MPa |
| Temperature | t | OT(1) | °C |
| Density | ρ | OD(2) | kg/m³ |
| Specific Volume | v | OV(3) | m³/kg |
| Specific Enthalpy | h | OH(4) | kJ/kg |
| Specific Entropy | s | OS(5) | kJ/(kg·K) |
| Specific Exergy | e | OE(6) | kJ/kg |
| Specific Internal Energy | u | OU(7) | kJ/kg |
| Isobaric Heat Capacity | cp | OCP(8) | kJ/(kg·K) |
| Isochoric Heat Capacity | cv | OCV(9) | kJ/(kg·K) |
| Speed of Sound | w | OW(10) | m/s |
| Isentropic Exponent | k | OKS(11) | - |
| Helmholtz Free Energy | f | OF(12) | kJ/kg |
| Gibbs Free Energy | g | OG(13) | kJ/kg |
| Compressibility Factor | z | OZ(14) | - |
| Steam Quality | x | OX(15) | - |
| Region | r | OR(16) | - |
| Isobaric Cubic Expansion Coefficient | αv | OEC(17) | 1/K |
| Isothermal Compressibility | kT | OKT(18) | 1/MPa |
| Partial Derivative (∂V/∂T)p | - | ODVDT(19) | m³/(kg·K) |
| Partial Derivative (∂V/∂p)T | - | ODVDP(20) | m³/(kg·MPa) |
| Partial Derivative (∂P/∂T)v | - | ODPDT(21) | MPa/K |
| Isothermal Throttling Coefficient | δt | OIJTC(22) | kJ/(kg·MPa) |
| Joule-Thomson Coefficient | μ | OJTC(23) | K/MPa |
| Dynamic Viscosity | η | ODV(24) | Pa·s |
| Kinematic Viscosity | ν | OKV(25) | m²/s |
| Thermal Conductivity | λ | OTC(26) | W/(m·K) |
| Thermal Diffusivity | a | OTD(27) | m²/s |
| Prandtl Number | Pr | OPR(28) | - |
| Surface Tension | σ | OST(29) | N/m |
| Static Dielectric Constant | ε | OSDC(30) | - |
| Isochoric Pressure Coefficient | β | OPC(31) | 1/K |
| Isothermal Stress Coefficient | βp | OBETAP(32) | kg/m³ |
| Fugacity Coefficient | fi | OFI(33) | - |
| Fugacity | f* | OFU(34) | MPa |
| Relative Pressure Coefficient | αp | OALFAP(35) | 1/K |

#### 4.1.4 Usage Example

```rust
use seuif97::*;

fn main() {
    let p: f64 = 16.0;    // MPa
    let t: f64 = 535.1;   // °C

    // Calculate enthalpy
    let h = pt(p, t, OH);
    println!("h = {:.3f} kJ/kg", h);

    // Calculate entropy with specified region
    let s = pt(p, t, (OS, 1));
    println!("s = {:.3f} kJ/(kg·K)", s);
}
```

### 4.2 C Language Interface

#### 4.2.1 Build Options

```bash
# cdecl calling convention
cargo build -r --features cdecl

# stdcall calling convention (Windows 64-bit)
cargo build -r --features stdcall

# stdcall calling convention (Windows 32-bit)
cargo build -r --target=i686-pc-windows-msvc --features stdcall
```

#### 4.2.2 Function Declarations

```c
double pt(double p, double t, short o_id);
double ph(double p, double h, short o_id);
double ps(double p, double s, short o_id);
double pv(double p, double v, short o_id);
double tv(double t, double v, short o_id);
double th(double t, double h, short o_id);
double ts(double t, double s, short o_id);
double hs(double h, double s, short o_id);
double px(double p, double x, short o_id);
double tx(double t, double x, short o_id);
double hx(double h, double x, short o_id);
double sx(double s, double x, short o_id);
```

#### 4.2.3 Usage Example

```c
#include <stdio.h>

#define OH 4
#define OS 5

extern double pt(double p, double t, short o_id);

int main(void) {
    double p = 16.0;
    double t = 530.0;
    double h = pt(p, t, OH);
    double s = pt(p, t, OS);
    printf("p,t %f,%f h= %f s= %f\n", p, t, h, s);
    return 0;
}
```

### 4.3 Python Interface

#### 4.3.1 Installation

```bash
# Install from PyPI
pip install seuif97

# Build and install locally
python ./setup.py install
```

#### 4.3.2 Usage Example

```python
from seuif97 import *

OH = 4
p = 16.0
t = 535.1

# Calculate using property ID
h = pt(p, t, OH)

# Calculate specific property directly
s = pt2s(p, t)

print(f"p={p}, t={t} h={h:.3f} s={s:.3f}")
```

#### 4.3.3 Python Convenience Functions

| Function | Description |
|----------|-------------|
| `pt2h(p, t)` | Calculate enthalpy from (p,t) |
| `pt2s(p, t)` | Calculate entropy from (p,t) |
| `pt2v(p, t)` | Calculate specific volume from (p,t) |
| `ph2t(p, h)` | Calculate temperature from (p,h) |
| `ps2t(p, s)` | Calculate temperature from (p,s) |
| `hs2p(h, s)` | Calculate pressure from (h,s) |
| `hs2t(h, s)` | Calculate temperature from (h,s) |


## 5. Physical Constants

### 5.1 Key Constants

```rust
pub const K: f64 = 273.15;                    // Celsius temperature conversion constant
pub const RGAS_WATER: f64 = 0.461526;        // Gas constant kJ/(kg·K)

// Critical point parameters
pub const TC_WATER: f64 = 647.096;           // Critical temperature K
pub const PC_WATER: f64 = 22.064;            // Critical pressure MPa
pub const DC_WATER: f64 = 322.0;             // Critical density kg/m³

// Boundary constants
pub const P_MIN: f64 = 0.000611212677444;   // Minimum pressure MPa
pub const P_MAX1: f64 = 100.0;               // Region 1 maximum pressure
pub const T_MAX1: f64 = 623.15;              // Region 1 maximum temperature
```

## 6 Error Handling

### 6.1 Error Code Definitions

| Error Code | Constant | Meaning |
|------------|----------|---------|
| -9999 | `INVALID_VALUE` | Invalid value |
| -1000 | `INVALID_OUTID` | Invalid property ID |
| -2100 | `INVALID_P` | Invalid pressure |
| -2101 | `INVALID_T` | Invalid temperature |
| -2102 | `INVALID_S` | Invalid entropy |
| -2103 | `INVALID_H` | Invalid enthalpy |
| -2201 | `INVALID_PT` | Invalid (p,T) combination |
| -2202 | `INVALID_HS` | Invalid (h,s) combination |

### 6.2 Input Validation Flow

```
Input Parameters → Range Check → Region Determination → Property Calculation → Return Result
                        ↓
                   Out of Range?
                        ↓
                   Return Error Code
```

## 7. Build and Testing

### 7.1 Build Commands

```bash
# Development build
cargo build

# Release build (optimized)
cargo build --release

# Build Python extension
cargo build --release --features python

# Run tests
cargo test

# Run benchmarks
cargo bench
```

### 7.2 Test Suite

The project includes comprehensive test cases:

| Test File | Test Content |
|-----------|--------------|
| `pt_test.rs` | (p,T) input pair tests |
| `ph_test.rs` | (p,h) input pair tests |
| `ps_test.rs` | (p,s) input pair tests |
| `pv_test.rs` | (p,v) input pair tests |
| `th_test.rs` | (t,h) input pair tests |
| `ts_test.rs` | (t,s) input pair tests |
| `tv_test.rs` | (t,v) input pair tests |
| `hs_test.rs` | (h,s) input pair tests |
| `hxsx_test.rs` | Wet steam region tests |
| `cross_test.rs` | Cross-region boundary tests |

### 7.3 Performance Benchmarks

The project uses the [Criterion](https://crates.io/crates/criterion) framework for performance benchmarking. The benchmark file is located at `benches/speed_benchmark.rs`.

#### 7.3.1 Test Content

Performance benchmarks cover typical calculation scenarios for each thermodynamic region:

| Test Function | Region | Input Parameters | Calculated Property |
|--------------|--------|------------------|---------------------|
| `pt2h_reg1` | Region 1 | p=3.0MPa, t=26.85°C | Specific Enthalpy |
| `pt2s_reg1` | Region 1 | p=3.0MPa, t=26.85°C | Specific Entropy |
| `pt2h_reg2` | Region 2 | p=0.0035MPa, t=26.85°C | Specific Enthalpy |
| `pt2s_reg2` | Region 2 | p=0.0035MPa, t=26.85°C | Specific Entropy |
| `tv2h_reg3` | Region 3 | t=376.85°C, v=0.002 m³/kg | Specific Enthalpy |
| `tv2s_reg3` | Region 3 | t=376.85°C, v=0.002 m³/kg | Specific Entropy |
| `pT2h_reg5` | Region 5 | p=0.5MPa, t=1226.85°C | Specific Enthalpy |
| `pT2s_reg5` | Region 5 | p=0.5MPa, t=1226.85°C | Specific Entropy |

#### 7.3.2 Test Code Implementation

```rust
use criterion::{black_box, criterion_group, criterion_main, Criterion};

use seuif97::*;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("pt2h_reg1", |b| {
        b.iter(|| pt(black_box(3.0), black_box(300.0 - 273.15), black_box(OH)))
    });
    c.bench_function("pt2s_reg1", |b| {
        b.iter(|| pt(black_box(3.0), black_box(300.0 - 273.15), black_box(OS)))
    });
    // ... other region tests
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
```

#### 7.3.3 Running Performance Tests

```bash
# Run all benchmarks
cargo bench

# View generated HTML report
# Report is located at target/criterion/index.html
```

#### 7.3.4 Performance Test Configuration

Criterion dependency and benchmark settings are configured in `Cargo.toml`:

```toml
[dev-dependencies]
criterion = { version = "0.5.1", features = ["html_reports"] }

[[bench]]
name = "speed_benchmark"
harness = false
```

**Configuration Notes**:
- `features = ["html_reports"]`: Generate HTML format performance reports
- `harness = false`: Use Criterion's benchmark interface instead of built-in harness

---

## 8. Deployment and Integration

### 8.1 Dynamic Library Deployment

Pre-compiled dynamic libraries are located in the `dynamic_lib/` directory:

| Platform | File | Path |
|----------|------|------|
| Windows 64-bit | `seuif97.dll` | `dynamic_lib/windows_x64/` |
| Windows 32-bit | `seuif97.dll` | `dynamic_lib/windows_x86/` |
| Linux 64-bit | `libseuif97.so` | `dynamic_lib/linux_x64/` |

### 8.2 Multi-language Integration Examples

Supported programming languages:
- ✅ Rust
- ✅ C / C++
- ✅ Python
- ✅ C#
- ✅ Java
- ✅ Fortran
- ✅ Go
- ✅ Excel VBA

## 9. Performance Optimization Recommendations

### 9.1 Usage Recommendations

1. **Batch Calculations**: For large numbers of calculations, use loop tiling techniques to fully utilize cache
2. **Pre-determine Region**: If the calculation region is known, specifying region parameters directly avoids region determination overhead
3. **Avoid Redundant Calculations**: For multiple queries with the same input parameter pairs, consider caching results

### 9.2 Performance Comparison

| Implementation | Performance | Description |
|----------------|-------------|-------------|
| SEUIF97 | Baseline | Optimized high-speed implementation |
| Rust standard library `powi()` | 5-20x slower | Unoptimized direct implementation |

## 10. Maintenance and Contribution

### 10.1 Code Standards

- Use Rust 2021 edition
- Follow `rustfmt` code formatting rules
- Use `clippy` for code checks

### 10.2 Contribution Process

1. Fork the repository
2. Create a feature branch
3. Write code and tests
4. Run tests to ensure they pass
5. Submit a Pull Request

### 10.3 Version Management

Version format: `MAJOR.MINOR.PATCH`

- **MAJOR**: API incompatible changes
- **MINOR**: New features, backward compatible
- **PATCH**: Bug fixes, backward compatible

## 11. References

* https://iapws.org/documents/release/IF97-Rev


**Document Version**: v1.2.2

**Generated Date**: 2024

**Author**: Cheng Maohua <cmh@seu.edu.cn>

**Project URL**: https://github.com/thermalogic/RustSEUIF97