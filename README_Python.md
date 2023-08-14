# SEUIF97

The version of **seuif97 2**  is the Python API of the high-speed IAPWS-IF97 package in Rust. It is suitable for computation-intensive calculations, such as heat cycle calculations, simulations of non-stationary processes, real-time process monitoring and optimizations.   
 
Through the high-speed package, the results of the IAPWS-IF97 are accurately produced at about 5-20x speed-up compared to  using the `powi()` of the Rust standard library in the `for`loop directly when computing the basic equations of Region 1,2,3.

**The Fast Methods**

1. The multi-step method unleashes the full power of the compiler optimizations while using `powi()` with the `for` loop
2. The recursive  method computes the polynomial values of the base variable and its derivatives

In seuif97, [36 thermodynamic, transport and  further properties](#properties) can be calculated. 

The following 12 input pairs are implemented:

```txt
  (p,t) (p,h) (p,s) (p,v) 
  
  (p,x) (t,x) (h,x) (s,x) 

  (t,h) (t,s) (t,v) 

  (h,s)
```

## The functions 

The two types of functions are provided in the package

 1. the input propertry pairs and the property ID([o_id](#properties)) to get the value of the specified property

 2. the input propertry pairs to get the one of  `p`,`t`,`h`,`s`,`v` or `x` directly

### The input propertry pairs and the property ID 

```python 
  ??(in1,in2,o_id)
```

* the first,second input parameters : the input propertry pairs
* the third input parametes: the property ID of the calculated property - [o_id](#properties)
* the return: the calculated property value of o_id

```python
pt(p,t,o_id)
ph(p,h,o_id)
ps(p,s,o_id)
pv(p,v,o_id)

th(t,h,o_id)
ts(t,s,o_id)
tv(t,v,o_id)

hs(h,s,o_id)

px(p,x,o_id)
tx(p,x,o_id)
hx(h,x,o_id)
sx(s,x,o_id)
```

### The input propertry pairs 

```python 
  ??2?(in1,in2)
```

* the `?` in `2?` is the one of `p`,`t`,`h`,`s`,`v` or `x`

```python
pt2h(p, t)  pt2s(p, t)  pt2v(p, t)  pt2x(p, t)
ph2t(p, h)  ph2s(p, h)  ph2v(p, h)  ph2x(p, h)   
ps2t(p, s)  ps2h(p, s)  ps2v(p, s)  ps2x(p, s)  
pv2t(p, v)  pv2h(p, v)  pv2s(p, v)  pv2x(p, v)  

hs2t(h, s)  hs2p(h, s)  hs2v(p, s)  hs2x(h, s)    

th2p(t, h)  th2s(t, h)  th2v(t, h)  th2x(t, h)   
ts2p(t, s)  ts2h(t, s)  th2v(t, s)  ts2x(t, s)  
tv2p(t, v)  tv2h(t, v)  tv2s(t, v)  tv2x(t, v)  

px2t(p, x)  px2h(p, x)  px2s(p, x)  px2v(p, x)
tx2p(t, x)  tx2h(t, x)  tx2s(t, x)  tx2v(t, x)

hx2p(h, x)  hx2t(h, x)  hx2s(h, x)  hx2v(h, x)
sx2p(s, x)  sx2t(s, x)  sx2h(s, x)  sx2v(s, x)
```

## Examples

```python
from seuif97 import *

OH=4

p=16.0
t=535.1
# ??(in1,in2,o_id)
h=pt(p,t,OH)
# ??2?(in1,in2)
s=pt2s(p,t)
print(f"p={p}, t={t} h={h:.3f} s={s:.3f}")
```
    
## Properties

| Propertry                             |    Unit     | Symbol | o_id  | o_id(i32)|
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
| Isobari cubic expansion coefficient   |     1/K     |   ɑv   |  OEC  |       17 |
| Isothermal compressibility            |    1/MPa    |    kT  |  OKT  |       18 |
| Partial derivative (∂V/∂T)p           |  m³/(kg·K)  |(∂V/∂T)p| ODVDT |       19 |
| Partial derivative (∂V/∂p)T           | m³/(kg·MPa) |(∂v/∂p)t| ODVDP |       20 |
| Partial derivative (∂P/∂T)v           |    MPa/K    |(∂p/∂t)v| ODPDT |       21 |
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
| Relative pressure coefficient         |     1/K     |    αp  | OAFLAP|        35|







