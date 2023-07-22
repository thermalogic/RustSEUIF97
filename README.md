# IF97

IF97 is the high-speed package of IAPWS-IF97 in Rust. It is suitable for computation-intensive calculations，such as heat cycle calculations, simulations of non-stationary processes, real-time process monitoring and optimizations.
 
Through the high-speed package, the results of the IAPWS-IF97 are accurately produced at about 5-15x speed-up compared to  the `powi()` of the Rust standard library，2-3x speed-up compared to the C implementation using the same fast algorithms.

* The comparison results of the computing-time are obtained using the [criterion.rs](https://bheisler.github.io/criterion.rs/book/index.html). 

The following input pairs are implemented: 

```
(p,t) (p,h) (p,s) (p,v) 

(t,h) (t,s) (t,v) 

(p,x) (t,x) (h,x) (s,x) 

(h,s)  
```
## Usage

Install the crate

```bash
cargo add if97
```

The type of functions are provided in the if97 package:

```rust
fn(f64,f64,i32) -> f64
``````

* the first,second input parameters: the input propertry pairs
* the third input parameter: the property ID of the calculated property - o_id
* the return: the calculated property value of o_id

```rust
pt(p:f64,t:f64,o_id:i32)->f64
ph(p:f64,h:f64,o_id:i32)->f64
ps(p:f64,s:f64,o_id:i32)->f64
pv(p:f64,v:f64,o_id:i32)->f64

th(t:f64,h:f64,o_id:i32)->f64
ts(t:f64,s:f64,o_id:i32)->f64
tv(t:f64,v:f64,o_id:i32)->f64

px(p:f64,x:f64,o_id:i32)->f64
tx(p:f64,x:f64,o_id:i32)->f64
hx(h:f64,x:f64,o_id:i32)->f64
sx(s:f64,x:f64,o_id:i32)->f64

hs(h:f64,s:f64,o_id:i32)->f64
```
Example

```rust
use if97::*;
fn main() {
    
    let p:f64 = 3.0;
    let t:f64= 300.0-273.15;
   
    let h=pt(p,t,OH);
    let s=pt(p,t,OS);
    let v=pt(p,t,OV);
    println!("p={p:.6} t={t:.6} h={t:.6} s={s:.6} v={v:.6}");   
}
```
    
## Properties

| Propertry                             |    Unit     | Symbol | o_id  | o_id(i32)|
| ------------------------------------- | :---------: |:------:|-----:|:--------:|
| Pressure                              |     MPa     |      p |   OP  |       0  |
| Temperature                           |     °C      |      t |   OT  |       1  |
| Density                               |   kg/m³     |      d |   OD  |       2  |
| Specific Volume                       |   m³/kg     |      v |   OV  |       3  |
| Specific enthalpy                     |    kJ/kg    |      h |   OH  |       4  |
| Specific entropy                      |  kJ/(kg·K)  |      s |   OS  |       5  |
| Specific exergy                       |    kJ/kg    |      e |   OE  |       6  |
| Specific internal energy              |    kJ/kg    |      u |   OU  |       7  |
| Specific isobaric heat capacity       |  kJ/(kg·K)  |     cp |  OCP  |       8  |
| Specific isochoric heat capacity      |  kJ/(kg·K)  |     cv |  OCV  |       9  |
| Speed of sound                        |     m/s     |      w |   OW  |       10 |
| Isentropic exponent                   |             |     ks |  OKS  |       11 |
| Specific Helmholtz free energy        |    kJ/kg    |     f  |   OF  |       12 |
| Specific Gibbs free energy            |    kJ/kg    |     g  |   OG  |       13 |
| Compressibility factor                |             |     z  |   OZ  |       14 |
| Steam quality                         |             |     x  |   OX  |       15 |
| Region                                |             |      r |   OR  |       16 |
| Isobaric volume expansion coefficient |     1/K     |     ec |  OEC  |       17 |
| Isothermal compressibility            |    1/MPa    |     kt |  OKT  |       18 |
| Partial derivative (∂V/∂T)p           |  m³/(kg·K)  | dvdtcp | ODVDT |       19 |
| Partial derivative (∂V/∂p)T           | m³/(kg·MPa) | dvdpct | ODVDP |       20 |
| Partial derivative (∂P/∂T)v           |    MPa/K    | dpdtcv | ODPDT |       21 |
| Isothermal Joule-Thomson coefficient  | kJ/(kg·MPa) |   iJTC | OIJTC |       22 |
| Joule-Thomson coefficient             |    K/MPa    |   joule| OJTC  |       23 |
| Dynamic viscosity                     |   Pa·s      |     dv |  ODV  |       24 |
| Kinematic viscosity                   |    m²/s    |     kv |  OKV  |       25 |
| Thermal conductivity                  |   W/(m.K)   |     tc |  OTC  |       26 |
| Thermal diffusivity                   |    m²/s    |     td |  OTD  |       27 |
| Prandtl number                        |             |     pr |  OPR  |       28 |
| Surface tension                       |    N/m      |     st |  OST  |       29 |
| Static Dielectric Constant            |             |    sdc | OSDC  |       30 |
| Isochoric pressure coefficient        |    1/K      |    pc  | OPC   |       31 |
| Isothermal stress coefficient         |   kg/m³     |   betap| OBETAP|       32 |
| Fugacity coefficient                  |             |      fi|   OFI |       33 |
| Fugacity                              |     Mpa     |      fu|   OFU |       34 |


