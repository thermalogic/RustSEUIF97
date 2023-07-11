# IF97

IF97 is the high-speed package of IAPWS-IF97 in Rust. It is suitable for computation-intensive calculations，such as heat cycle calculations, simulations of non-stationary processes, real-time process monitoring and optimizations.
 
Through the high-speed package, the results of the IAPWS-IF97 are accurately produced at about 5~10x speed-up compared to  the `powi()` of the Rust standard library，5x speed-up compared to  the shared library in C using the same fast algorithms.  

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

| Propertry                             |    Unit     | symbol | o_id | o_id value | 
| ------------------------------------- | :---------: | ------:|-----:|-----------:|
| Pressure                              |     MPa     |      p |  OH  |   0        |
| Temperature                           |     °C      |      t |  OT  |   1        |
| Density                               |   kg/m^3    |      d |  OD  |   2        |
| Specific Volume                       |   m^3/kg    |      v |  OV  |   3        |
| Specific enthalpy                     |    kJ/kg    |      h |  OH  |   4        |
| Specific entropy                      |  kJ/(kg·K)  |      s |  OS  |   5        |
| Specific internal energy              |    kJ/kg    |      u |  OU  |   7        |
| Specific isobaric heat capacity       |  kJ/(kg·K)  |     cp | OCP  |   8        |
| Specific isochoric heat capacity      |  kJ/(kg·K)  |     cv | OCV  |   9        |
| Speed of sound                        |     m/s     |      w |  OW  |   10       |
| Steam quality                         |             |      x |  OX  |   15       |
| Region                                |             |      r |  OR  |   16       |
| Dynamic viscosity                     |  kg/(m·s)   |     dv |  ODV |   24       |
| Kinematic viscosity                   |    m^2/s    |     kv |  OKV |   25       |
| Thermal conductivity                  |   W/(m.K)   |     tc |  OTC |   26       |
| Thermal diffusivity                   |   um^2/s    |     td |  OTD |   27       |
| Prandtl number                        |             |     pr |  OPR |   28       |
| Surface tension                       |    mN/m     |     st |  OST |   29       |
| Static Dielectric Constant            |             |    sdc | OSDC |   30       | 

