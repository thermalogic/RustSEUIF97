# SEUIF97

This package is the WebAssembly implementation of the high-speed IAPWS-IF97 package SEUIF97 in Rust. enabling fast and accurate thermodynamic property calculations for water and steam directly in the browser or Node.js. 
 
Through the high-speed package, the results of the IAPWS-IF97 are accurately produced at about 5-20x speed-up compared to using the `powi()` of the Rust standard library in the `for` loop directly when computing the basic equations of Regions 1, 2, and 3.

## Acceleration Methods

* Loop Tiling Method: Unleashes the full power of compiler optimizations, surpassing the performance of the single loop.

* Recurrence Method for Multi-Polynomial Evaluation: By leveraging the relationship between polynomials and their derivatives, only a single polynomial needs to be computed directly. The remaining values are derived via multiplication or division by the base. This approach eliminates redundant calculations and significantly improves performance.

## Related Packages

- [seuif97 (Python)](https://pypi.org/project/seuif97/) - Python version
- [seuif97 (Rust)](https://crates.io/crates/seuif97) - Rust version

## Installation

```bash
npm install seuif97
```

## Input Pairs and Properties

In the package, [36 thermodynamic, transport and further properties](#properties) can be calculated. 

The following 12 input pairs are implemented:

```txt
  (p,t) (p,h) (p,s) (p,v) 
  
  (p,x) (t,x) (h,x) (s,x) 

  (t,h) (t,s) (t,v) 

  (h,s)
```

## Functions

The two types of functions are provided in the package.

 1. the input property pairs and the property ID([o_id](#properties)) to get the value of the specified property

 2. the input property pairs to get one of  `p`,`t`,`h`,`s`,`v` or `x` directly

### The input property pairs and the property ID 

```python 
  ??(in1,in2,o_id)
```

* the first, second input parameters: the input property pairs
* the third input parameters: the property ID of the calculated property - [o_id](#properties)
* the return: the calculated property value of o_id

```javascript
pt(p,t,o_id)  ph(p,h,o_id) ps(p,s,o_id) pv(p,v,o_id)

th(t,h,o_id)  ts(t,s,o_id) v(t,v,o_id)

hs(h,s,o_id)

px(p,x,o_id) tx(p,x,o_id) hx(h,x,o_id)sx(s,x,o_id)
```

```javascript
import init, { pt } from 'seuif97';

await init();

const p = 16.0;  // MPa
const t = 535.1; // °C

// Calculate all properties
const enthalpy = pt(p, t, 4);     // kJ/kg
const entropy = pt(p, t, 5);      // kJ/(kg·K)

console.log('Properties at p = 16.0 MPa, t = 535.1 °C:');
console.log(`Enthalpy: ${enthalpy.toFixed(3)} kJ/kg`);
console.log(`Entropy: ${entropy.toFixed(5)} kJ/(kg·K)`);
```

### The input property pairs 

```python 
  ??2?(in1,in2)
```

* the `?` in `2?` is the one of `p`,`t`,`h`,`s`,`v` or `x`

```javascript
pt2h(p, t)  pt2s(p, t)  pt2v(p, t)  pt2x(p, t)
ph2t(p, h)  ph2s(p, h)  ph2v(p, h)  ph2x(p, h)   
ps2t(p, s)  ps2h(p, s)  ps2v(p, s)  ps2x(p, s)  
pv2t(p, v)  pv2h(p, v)  pv2s(p, v)  pv2x(p, v)  

hs2p(h, s)  hs2t(h, s)  hs2v(h, s)  hs2x(h, s)    

th2p(t, h)  th2s(t, h)  th2v(t, h)  th2x(t, h)   
ts2p(t, s)  ts2h(t, s)  ts2v(t, s)  ts2x(t, s)  
tv2p(t, v)  tv2h(t, v)  tv2s(t, v)  tv2x(t, v)  

px2t(p, x)  px2h(p, x)  px2s(p, x)  px2v(p, x)
tx2p(t, x)  tx2h(t, x)  tx2s(t, x)  tx2v(t, x)

hx2p(h, x)  hx2t(h, x)  hx2s(h, x)  hx2v(h, x)
sx2p(s, x)  sx2t(s, x)  sx2h(s, x)  sx2v(s, x)
```

```javascript
import init, { pt2h, pt2s, pt2v } from 'seuif97';

// Initialize the WASM module
await init();

// Calculate properties
const p = 16.0;  // Pressure in MPa
const t = 535.1; // Temperature in °C

const h = pt2h(p, t);
const s = pt2s(p, t);
const v = pt2v(p, t);

console.log(`h = ${h.toFixed(3)} kJ/kg`);
console.log(`s = ${s.toFixed(5)} kJ/(kg·K)`);
console.log(`v = ${v.toFixed(6)} m³/kg`);
```
![](https://raw.githubusercontent.com/thermalogic/RustSEUIF97/seuif97-pypi-multi-platform/img/turbine_hs.jpg)

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
| Isobaric cubic expansion coefficient   |     1/K     |   ɑv   |  OEC  |       17 |
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

