# SEUIF97

![npm version](https://img.shields.io/npm/v/seuif97)![NPM Downloads](https://img.shields.io/npm/dt/seuif97)![NPM Downloads](https://img.shields.io/npm/dm/seuif97)

The WebAssembly implementation of the high-speed IAPWS-IF97 package SEUIF97 in Rust, enabling fast and accurate thermodynamic property calculations for water and steam in the browser or Node.js. 
 
Through the high-speed package, the results of the IAPWS-IF97 are accurately produced at about **5-20x** speed-up compared to using the `powi()` of the Rust standard library in the `for` loop directly when computing the basic equations of Regions 1, 2 and 3.

This package supports **12 distinct input state pairs** for calculating **36 thermodynamic, transport, and derived properties**(see [Properties](#properties)).

## Acceleration Methods

* Loop Tiling Method: Unleashes the full power of compiler optimizations, surpassing the performance of the single loop.

* Recurrence Method for Multi-Polynomial Evaluation: By leveraging the relationship between polynomials and their derivatives, only a single polynomial needs to be computed directly. The remaining values are derived via multiplication or division by the base. This approach eliminates redundant calculations and significantly improves performance.

## Related Packages

-  [Rust crate: seuif97](https://crates.io/crates/seuif97)  
-  [Python package: seuif97](https://pypi.org/project/seuif97/)

## Installation

```bash
npm install seuif97
```

## API Reference

The two types of API are provided in the package.

 1.  Universal Functions (with o_id parameter)
     - These functions accept an input property pair plus a property ID([o_id](#properties)) to calculate the desired output property. For example: `pt(p,t,o_id)`, where `o_id` is the property ID of the calculated property.

 2. Convenience Functions (Direct Output)
    -  These functions directly calculate a specific property `(p,t,h,s,v,x)` without requiring the property ID parameter. For example: `pt2h（p,t)`

### Universal Functions (with o_id parameter)

The following 12 input pairs are implemented:

```javascript
pt(p,t,o_id)  ph(p,h,o_id) ps(p,s,o_id) pv(p,v,o_id)

th(t,h,o_id)  ts(t,s,o_id) tv(t,v,o_id)

hs(h,s,o_id)

px(p,x,o_id) tx(p,x,o_id) hx(h,x,o_id) sx(s,x,o_id)
```

### Convenience Functions (Direct Output)

The following 12 input pairs are implemented:

```javascript
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

### Example

```javascript
import init, { pt,pt2s } from 'seuif97';

await init();

const p = 16.0;  // Pressure in MPa
const t = 535.1; // Temperature in °C

// universal function (with o_id parameter)
const h = pt(p, t, 4);   
// convenience function (Direct Output)
const s = pt2s(p, t);

console.log('Properties at p = 16.0 MPa, t = 535.1 °C:');
console.log(`h: ${h.toFixed(3)} kJ/kg`);
console.log(`s: ${s.toFixed(5)} kJ/(kg·K)`);
```

**T-s Diagram**

![](https://raw.githubusercontent.com/thermalogic/RustSEUIF97/seuif97-pypi-multi-platform/img/ts_wasm.jpg)

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

