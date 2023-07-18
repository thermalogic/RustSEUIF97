# IF97

IF97 is the high-speed package of IAPWS-IF97 in Rust. It is suitable for computation-intensive calculations，such as heat cycle calculations, simulations of non-stationary processes, real-time process monitoring and optimizations.
 
Through the high-speed package, the results of the IAPWS-IF97 are accurately produced at about 5-10x speed-up compared to  the `powi()` of the Rust standard library，2-3x speed-up compared to the C implementation using the same fast algorithms.

The following input pairs are implemented: 

```
(p,t) (p,h) (p,s) (p,x) 

(t,x) 

(h,s)  
```
## Usage

Install the crate

```bash
cargo add if97
```

The type of functions are provided in the if97 package:

```rust
fn(f64,f64,i32) -> f64;
``````

* first,second input parameters: the input propertry pairs(f64)
* third input parameter: o_id - the property ID of the calculated property(i32)
* the return: the calculated property value(f64)

```rust
pt(p:f64,t:f64,o_id::i32)->f64;
ph(p:f64,h:f64,o_id::i32)->f64;
ps(p:f64,s:f64,o_id::i32)->f64;
hs(h:f64,s:f64,o_id::i32)->f64;
px(p:f64,x:f64,o_id::i32)->f64;
tx(p:f64,x:f64,o_id::i32)->f64;
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
    
## Properties in if97

| Propertry                             |    Unit     | Symbol | o_id |  
| ------------------------------------- | :---------: | ------:|-----:|
| Pressure                              |     MPa     |      p |  OP  |
| Temperature                           |     °C      |      t |  OT  |
| Density                               |   kg/m^3    |      d |  OD  |
| Specific Volume                       |   m^3/kg    |      v |  OV  |
| Specific enthalpy                     |    kJ/kg    |      h |  OH  |
| Specific entropy                      |  kJ/(kg·K)  |      s |  OS  |
| Specific internal energy              |    kJ/kg    |      u |  OU  |
| Specific isobaric heat capacity       |  kJ/(kg·K)  |     cp | OCP  |
| Specific isochoric heat capacity      |  kJ/(kg·K)  |     cv | OCV  |
| Speed of sound                        |     m/s     |      w |  OW  |
| Steam quality                         |             |      x |  OX  |
| Region                                |             |      r |  OR  |
| Dynamic viscosity                     |  kg/(m·s)   |     dv |  ODV |
| Kinematic viscosity                   |    m^2/s    |     kv |  OKV |
| Thermal conductivity                  |   W/(m.K)   |     tc |  OTC |
| Thermal diffusivity                   |   m^2/s     |     td |  OTD |
| Prandtl number                        |             |     pr |  OPR |
| Surface tension                       |    N/m      |     st |  OST |
| Static Dielectric Constant            |             |    sdc | OSDC | 

