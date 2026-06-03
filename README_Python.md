# SEUIF97

 ![PyPI](https://img.shields.io/pypi/v/seuif97) [![Downloads](https://static.pepy.tech/badge/seuif97)](https://pepy.tech/project/seuif97) [![Downloads](https://static.pepy.tech/badge/seuif97/month)](https://pepy.tech/project/seuif97)

**SEUIF97 Version 2** is the Python API of the high-speed IAPWS-IF97 package in Rust. 

SEUIF97 2, built on Rust, is a major upgrade over SEUIF97 1.* (built on C), delivering significant improvements in performance, functionality and ecosystem support.

It is suitable for computation-intensive calculations, such as heat cycle calculations, simulations of non-stationary processes, real-time process monitoring and optimizations.

Through the high-speed package, IAPWS-IF97 calculations achieve a **5-20x speedup** compared to direct implementations using the Rust standard library's `powi()` within loops for the basic equations of Regions 1, 2 and 3.

This package supports **12 distinct input state pairs** for calculating **36 thermodynamic, transport, and derived properties** (see [Properties](#properties)).

## Acceleration Methods

* Loop Tiling Method: Unleashes the full power of compiler optimizations, surpassing the performance of the single loop.

* Recurrence Method for Multi-Polynomial Evaluation: By leveraging the relationship between polynomials and their derivatives, only a single polynomial needs to be computed directly. The remaining values are derived via multiplication or division by the base. This approach eliminates redundant calculations and significantly boosts performance.

## What's New in SEUIF97 2 

| Feature                               | 1.*            | 2.*                        |
| ------------------------------------- | -------------- | -------------------------- |
| **Implementation**                    | C              | **Rust**                   |
| **Calculation Speed**                 | Baseline       | **~2x speedup**            |
| **Supported Properties**              | 30 properties  | **36 properties** (+6 new) |
| **Supported OS**                      | Windows, Linux | **Windows, Linux, macOS**  |
| **Universal Functions**               | ✓              | ✓                          |
| **Direct Property Functions**         | ✗              | **✓** (new)                |
| **Thermodynamic Process Calculation** | ✓              | ✗ (planned)                |

## API Reference

The package provides two types of APIs.

 1.  Universal Functions (with o_id parameter)
     - These functions accept an input property pair plus a property ID([o_id](#properties)) to calculate the desired output property. For example: `pt(p,t,o_id)`, where `o_id` specifies the output property.

 2. Direct Property Functions
    -  These functions directly calculate a specific property `(p,t,h,s,v,x)` without requiring the property ID parameter. For example: `pt2h(p,t)`.

### Universal Functions (with o_id parameter) 

The following 12 input pairs are implemented:

```python
pt(p,t,o_id) ph(p,h,o_id) ps(p,s,o_id) pv(p,v,o_id)

th(t,h,o_id) ts(t,s,o_id) tv(t,v,o_id)

hs(h,s,o_id)

px(p,x,o_id) tx(p,x,o_id) hx(h,x,o_id) sx(s,x,o_id)
```

### Direct Property Functions

The following 12 input pairs are implemented:

```python
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

## Usage 

```python
from seuif97 import *

OH=4

p=16.0
t=535.1
# universal functions (with o_id parameter)
h=pt(p,t,OH)
# direct property functions
s=pt2s(p,t)
print(f"p={p}, t={t} h={h:.3f} s={s:.3f}")
```

## Examples

* [T-S Diagram](https://github.com/thermalogic/RustSEUIF97/blob/seuif97-pypi-multi-platform/demo_using_lib/Diagram_T-S.py)

* [H-S Diagram](https://github.com/thermalogic/RustSEUIF97/blob/seuif97-pypi-multi-platform/demo_using_lib/Diagram_H-S.py)

* [H-S Diagram of Steam Turbine Expansion](https://github.com/thermalogic/RustSEUIF97/blob/seuif97-pypi-multi-platform/demo_using_lib/Turbine_H-S.py)

* [The Hybrid Steady-state Simulator of Rankine Cycle in Python](https://github.com/thermalogic/PyRankine)

![T-S Diagram](https://github.com/thermalogic/RustSEUIF97/raw/seuif97-pypi-multi-platform/img/T-S.jpg)

## Properties

| Property                              |    Unit     | Symbol | o_id  | o_id(i32)|
| ------------------------------------- | :---------: |:------:|------:|:--------:|
| Pressure                              |     MPa     |      p |   OP  |       0  |
| Temperature                           |     °C      |      t |   OT  |       1  |
| Density                               |   kg/m³     |      ρ |   OD  |       2  |
| Specific Volume                       |   m³/kg     |      v |   OV  |       3  |
| Specific enthalpy                     |    kJ/kg    |      h |   OH  |       4  |
| Specific entropy                      |  kJ/(kg·K)  |      s |   OS  |       5  |
| Specific exergy                       |    kJ/kg    |      e |   OE  |       6  |
| Specific internal energy              |    kJ/kg    |      u |   OU  |       7  |
| Specific isobaric heat capacity       |  kJ/(kg·K)  |     cp |  OCP  |       8  |
| Specific isochoric heat capacity      |  kJ/(kg·K)  |     cv |  OCV  |       9  |
| Speed of sound                        |     m/s     |      w |   OW  |       10 |
| Isentropic exponent                   |             |     k  |  OKS  |       11 |
| Specific Helmholtz free energy        |    kJ/kg    |     f  |   OF  |       12 |
| Specific Gibbs free energy            |    kJ/kg    |     g  |   OG  |       13 |
| Compressibility factor                |             |     z  |   OZ  |       14 |
| Steam quality                         |             |     x  |   OX  |       15 |
| Region                                |             |     r  |   OR  |       16 |
| Isobaric cubic expansion coefficient  |     1/K     |   ɑv   |  OEC  |       17 |
| Isothermal compressibility            |    1/MPa    |    kT  |  OKT  |       18 |
| Partial derivative (∂V/∂T)p           |  m³/(kg·K)  |(∂V/∂T)p| ODVDT |       19 |
| Partial derivative (∂V/∂p)T           | m³/(kg·MPa) |(∂v/∂p)T| ODVDP |       20 |
| Partial derivative (∂P/∂T)v           |    MPa/K    |(∂p/∂T)v| ODPDT |       21 |
| Isothermal throttling coefficient     | kJ/(kg·MPa) |   δt   | OIJTC |       22 |
| Joule-Thomson coefficient             |    K/MPa    |    μ   | OJTC  |       23 |
| Dynamic viscosity                     |   Pa·s      |    η   |  ODV  |       24 |
| Kinematic viscosity                   |    m²/s     |    ν   |  OKV  |       25 |
| Thermal conductivity                  |   W/(m.K)   |    λ   |  OTC  |       26 |
| Thermal diffusivity                   |    m²/s     |    a   |  OTD  |       27 |
| Prandtl number                        |             |    Pr  |  OPR  |       28 |
| Surface tension                       |    N/m      |    σ   |  OST  |       29 |
| Static Dielectric Constant            |             |    ε   | OSDC  |       30 |
| Isochoric pressure coefficient        |    1/K      |    β   | OPC   |       31 |
| Isothermal stress coefficient         |   kg/m³     |    βp  | OBETAP|       32 |
| Fugacity coefficient                  |             |    fi  |   OFI |       33 |
| Fugacity                              |     MPa     |     f* |   OFU |       34 |
| Relative pressure coefficient         |     1/K     |    αp  | OAFLAP|       35 |