# SEUIF97

![docs.rs](https://img.shields.io/docsrs/seuif97)  [![Build test](https://github.com/thermalogic/RustSEUIF97/actions/workflows/rust.yml/badge.svg)](https://github.com/thermalogic/RustSEUIF97/actions/workflows/rust.yml) ![Crates.io Version](https://img.shields.io/crates/v/seuif97) ![Crates.io Total Downloads](https://img.shields.io/crates/d/seuif97) ![Crates.io Downloads (recent)](https://img.shields.io/crates/dr/seuif97)![PyPI](https://img.shields.io/pypi/v/seuif97) [![Downloads](https://static.pepy.tech/badge/seuif97)](https://pepy.tech/project/seuif97) [![Downloads](https://static.pepy.tech/badge/seuif97/month)](https://pepy.tech/project/seuif97) ![npm version](https://img.shields.io/npm/v/seuif97)![NPM Downloads](https://img.shields.io/npm/dt/seuif97)![NPM Downloads](https://img.shields.io/npm/dm/seuif97)

This is the Rust implementation of the high-speed IAPWS-IF97 package **SEUIF97** with **C, Python and WASM** bindings. It is designed for computation-intensive tasks, such as simulating non-stationary processes, on-line process monitoring, and optimization.

Through the high-speed package, IAPWS-IF97 calculations achieve a **5-20x speedup** compared to direct implementations using the Rust standard library's `powi()` within loops for the basic equations of Regions 1, 2 and 3.

**SEUIF97** also significantly outperforms various approximate equations and algorithms typically used for fast water and steam property calculations.

This package supports **12 distinct input state pairs** for calculating **36 thermodynamic, transport, and derived properties** (see [Properties](#properties)).

## Acceleration Methods

- Loop Tiling Method: Unleashes the full power of compiler optimizations, surpassing the performance of the single loop.
- Recurrence Method for Multi-Polynomial Evaluation: By utilizing the relationship between polynomials and their derivatives, only a single polynomial needs to be computed directly. The remaining values are derived via multiplication or division by the base. This approach eliminates redundant calculations and significantly boosts computational performance.

Please refer to [The acceleration methods](./docs/the_acceleration_methods.md) for more details on the algorithm

## What's New in the Rust Version

The Rust version of SEUIF97 is a major upgrade over [the original C implementation](https://github.com/thermalogic/SEUIF97), delivering significant improvements in performance, functionality, and ecosystem support.

| Feature                               | C Version      | Rust Version               |
| ------------------------------------- | -------------- | -------------------------- |
| **Calculation Speed**                 | Baseline       | **2-3× speedup**           |
| **Supported Properties**              | 30 properties  | **36 properties** (+6 new) |
| **Package Distribution**              | PyPI only      | **Crates.io, PyPI, npm**   |
| **Universal Functions**               | ✓              | ✓                          |
| **Direct Property Functions**         | ✗              | **✓** (new)                |
| **Thermodynamic Process Calculation** | ✓              | ✗ (planned)                |

For detailed comparison and key improvements, see [Rust vs C Version Comparison](./docs/RUST_VS_C_COMPARISON.md).

## Install the crate

```bash
cargo add seuif97
```

## API Reference

The package provides two types of API.

1. Universal Functions (with o\_id and optional region parameter)
   - These functions accept an input property pair plus a property ID([o\_id](#properties)) to calculate the desired output property. For example: `pt(p,t,o_id,<region>)`, where `o_id` specifies the output property, and `region` is optional.
2. Direct Property Functions
   - These functions directly calculate a specific property `(p,t,h,s,v,x)` from the input property pairs without requiring the property ID parameter. For example: `pt2h(p,t)`.

**C, Python and WASM** bindings support all functions of both types in Rust, **except for the** **`optional region`** **parameter**.

### Universal Functions (with o\_id and `<optional>` region parameter)

The following function signature is provided:

```txt
struct o_id_region_args {
   o_id: i32,
   region: i32,
}

fn<R>(f64,f64,R) -> f64
where
    R: Into<o_id_region_args>,
```

- the first,second input parameters(f64) : the input property pairs
- the third and fourth input parameters<R>:
  - the third : the property ID of the calculated property - [o\_id](#properties)
  - the fourth (`optional`) parameter: IAPWS-IF97 region specification
- the return(f64): the calculated property value of o\_id

The following 12 input pairs are implemented:

```txt
pt<R>(p:f64,t:f64,o_id_region:R)->f64
ph<R>(p:f64,h:f64,o_id_region:R)->f64
ps<R>(p:f64,s:f64,o_id_region:R)->f64
pv<R>(p:f64,v:f64,o_id_region:R)->f64

th<R>(t:f64,h:f64,o_id_region:R)->f64
ts<R>(t:f64,s:f64,o_id_region:R)->f64
tv<R>(t:f64,v:f64,o_id_region:R)->f64

hs<R>(h:f64,s:f64,o_id_region:R)->f64

px(p:f64,x:f64,o_id:i32)->f64
tx(p:f64,x:f64,o_id:i32)->f64
hx(h:f64,x:f64,o_id:i32)->f64
sx(s:f64,x:f64,o_id:i32)->f64
```

### Direct Property Functions

The following 12 input pairs are implemented:

```txt
pt2h(p, t)  pt2s(p, t)  pt2v(p, t)  pt2x(p, t)
ph2t(p, h)  ph2s(p, h)  ph2v(p, h)  ph2x(p, h)   
ps2t(p, s)  ps2h(p, s)  ps2v(p, s)  ps2x(p, s)  
pv2t(p, v)  pv2h(p, v)  pv2s(p, v)  pv2x(p, v)  

th2p(t, h)  th2s(t, h)  th2v(t, h)  th2x(t, h)   
ts2p(t, s)  ts2h(t, s)  ts2v(t, s)  ts2x(t, s)  
tv2p(t, v)  tv2h(t, v)  tv2s(t, v)  tv2x(t, v)  

hs2p(h, s)  hs2t(h, s)  hs2v(h, s)  hs2x(h, s)    

px2t(p, x)  px2h(p, x)  px2s(p, x)  px2v(p, x)
tx2p(t, x)  tx2h(t, x)  tx2s(t, x)  tx2v(t, x)

hx2p(h, x)  hx2t(h, x)  hx2s(h, x)  hx2v(h, x)
sx2p(s, x)  sx2t(s, x)  sx2h(s, x)  sx2v(s, x)
```

### Usage

```rust
use seuif97::*;
fn main() {
    
    let p:f64 = 3.0;
    let t:f64= 300.0-273.15;
    // universal functions (with o_id parameter only)
    let h=pt(p,t,OH);
    // universal functions with explicit region for faster calculation
    let s=pt(p,t,(OS,1));
    // direct property functions
    let v=pt2v(p,t);

    println!("p={p:.6} t={t:.6} h={h:.6} s={s:.6} v={v:.6}");   
}
```

## C Shared Library

**Building the dynamic link library**

- cdecl

```bash
cargo build -r --features cdecl
```

- stdcall: Windows API functions(64bit)

```bash
cargo build -r --features stdcall
```

- stdcall: Windows API functions(32bit)

```bash
cargo build -r  --target=i686-pc-windows-msvc --features stdcall
```

Pre-compiled dynamic link libraries for Windows, Linux and macOS are available in [GitHub Releases](https://github.com/thermalogic/RustSEUIF97/releases).

Legacy pre-compiled libraries are also provided in the [./dynamic\_lib/](./dynamic_lib/) directory.

- `seuif97.dll`: [windows\_x64](./dynamic_lib/windows_x64/)  and [windows\_x86](./dynamic_lib/windows_x86/)
- `libseuif97.so`: [linux\_x64](./dynamic_lib/linux_x64/)

The shared library supports all functions of both types in Rust, except for the `optional region` parameter. For example in C:

```c
double pt(double p,double t,short o_id);
double pt2s(double p,double t);
```

Interfaces and examples are provided in the [./demo\_using\_lib/](./demo_using_lib/) directory, supporting a wide range of languages and environments

- C/C++, Python, C#, Java, Excel VBA, Rust, Fortran, Golang, JavaScript/TypeScript

```c
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define OH 4

extern double pt(double p,double t,short o_id);
extern double pt2s(double p,double t);

int main(void)
{
    double p = 16.0;
    double t = 530.0;
    // universal functions (with o_id parameter only)
    double h = pt(p, t, OH);
    // direct property functions
    double s = pt2s(p, t);
    printf("p,t %f,%f h= %f s= %f\n", p, t, h, s);
    return EXIT_SUCCESS;
}
```

**Comprehensive Cross-language Examples**

- [The Rankine Cycle Steady-state Simulator in Python, C++, Rust and Modelica](https://github.com/thermalogic/SimRankine)

## Python binding

**Install from pypi**

- <https://pypi.org/project/seuif97/>

```bash
pip install seuif97
```

```python
from seuif97 import *

OH=4

p=16.0
t=535.1
# universal functions (with o_id parameter only)
h=pt(p,t,OH)
# direct property functions
s=pt2s(p,t)
print(f"p={p}, t={t} h={h:.3f} s={s:.3f}")
```

**The Comprehensive Examples in Python**

- [T-S Diagram](./demo_using_lib/Diagram_T-S.py)
- [H-S Diagram](./demo_using_lib/Diagram_H-S.py)
- [H-S Diagram of Steam Turbine Expansion](./demo_using_lib/Turbine_H-S.py)
- [The Hybrid Steady-state Simulator of Rankine Cycle in Python](https://github.com/thermalogic/PyRankine)

![T-S Diagram](./img/T-S.jpg)

## WASM binding

- WASM - [README\_WASM.md](./README_WASM.md)
- NPM package: [seuif97](https://www.npmjs.com/seuif97)

```javascript
import init, { pt, pt2s } from 'seuif97';

await init();

const p = 16.0;  // MPa
const t = 535.1; // °C
// universal functions (with o_id parameter only)
const h = pt(p, t, 4);     // kJ/kg
// direct property functions
const s = pt2s(p, t);      // kJ/(kg·K)

console.log('Properties at p = 16.0 MPa, t = 535.1 °C:');
console.log(`H: ${h.toFixed(3)} kJ/kg`);
console.log(`S: ${s.toFixed(5)} kJ/(kg·K)`);
```

## Properties

| Property                             |     Unit    |  Symbol  |  o\_id | o\_id(i32) |
| ------------------------------------ | :---------: | :------: | -----: | :--------: |
| Pressure                             |     MPa     |     p    |     OP |      0     |
| Temperature                          |      °C     |     t    |     OT |      1     |
| Density                              |    kg/m³    |     ρ    |     OD |      2     |
| Specific Volume                      |    m³/kg    |     v    |     OV |      3     |
| Specific enthalpy                    |    kJ/kg    |     h    |     OH |      4     |
| Specific entropy                     |  kJ/(kg·K)  |     s    |     OS |      5     |
| Specific exergy                      |    kJ/kg    |     e    |     OE |      6     |
| Specific internal energy             |    kJ/kg    |     u    |     OU |      7     |
| Specific isobaric heat capacity      |  kJ/(kg·K)  |    cp    |    OCP |      8     |
| Specific isochoric heat capacity     |  kJ/(kg·K)  |    cv    |    OCV |      9     |
| Speed of sound                       |     m/s     |     w    |     OW |     10     |
| Isentropic exponent                  |    <br />   |     k    |    OKS |     11     |
| Specific Helmholtz free energy       |    kJ/kg    |     f    |     OF |     12     |
| Specific Gibbs free energy           |    kJ/kg    |     g    |     OG |     13     |
| Compressibility factor               |    <br />   |     z    |     OZ |     14     |
| Steam quality                        |    <br />   |     x    |     OX |     15     |
| Region                               |    <br />   |     r    |     OR |     16     |
| Isobaric cubic expansion coefficient |     1/K     |    ɑv    |    OEC |     17     |
| Isothermal compressibility           |    1/MPa    |    kT    |    OKT |     18     |
| Partial derivative (∂V/∂T)p          |  m³/(kg·K)  | (∂V/∂T)p |  ODVDT |     19     |
| Partial derivative (∂V/∂p)T          | m³/(kg·MPa) | (∂v/∂p)t |  ODVDP |     20     |
| Partial derivative (∂P/∂T)v          |    MPa/K    | (∂p/∂t)v |  ODPDT |     21     |
| Isothermal throttling coefficient    | kJ/(kg·MPa) |    δt    |  OIJTC |     22     |
| Joule-Thomson coefficient            |    K/MPa    |     μ    |   OJTC |     23     |
| Dynamic viscosity                    |     Pa·s    |     η    |    ODV |     24     |
| Kinematic viscosity                  |     m²/s    |     ν    |    OKV |     25     |
| Thermal conductivity                 |   W/(m.K)   |     λ    |    OTC |     26     |
| Thermal diffusivity                  |     m²/s    |     a    |    OTD |     27     |
| Prandtl number                       |    <br />   |    Pr    |    OPR |     28     |
| Surface tension                      |     N/m     |     σ    |    OST |     29     |
| Static Dielectric Constant           |    <br />   |     ε    |   OSDC |     30     |
| Isochoric pressure coefficient       |     1/K     |     β    |    OPC |     31     |
| Isothermal stress coefficient        |    kg/m³    |    βp    | OBETAP |     32     |
| Fugacity coefficient                 |    <br />   |    fi    |    OFI |     33     |
| Fugacity                             |     MPa     |    f\*   |    OFU |     34     |
| Relative pressure coefficient        |     1/K     |    αp    | OAFLAP |     35     |

