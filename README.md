# IF97

IF97 is the high-speed package of IAPWS-IF97 in Rust. It is suitable for computation-intensive calculations，such as heat cycle calculations, simulations of non-stationary processes, real-time process monitoring and optimizations.
 
Through the high-speed package, the results of the IAPWS-IF97 are accurately produced at about 5~10x speed-up compared to  the `powi()` of the Rust standard library，5x speed-up compared to [the shared library in C using the same fast algorithms](https://github.com/thermalogic/SEUIF97).  

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
Examples

```rust
use if97::pt;
fn main() {
    
    let p:f64 = 3.0;
    let t:f64= 300.0-273.15;
   
    let h=pt(p,t,if97::OH);
    let s=pt(p,t,if97::OS);
    let v=pt(p,t,if97::OV);
    
    println!("p={} t={} h={} s={} v={}",p,t,h,s,v);    
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
| Thermal diffusivity                   |   um^2/s    |     td |  OTD |
| Prandtl number                        |             |     pr |  OPR |
| Surface tension                       |    mN/m     |     st |  OST |
| Static Dielectric Constant            |             |    sdc | OSDC | 

