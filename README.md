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
* third input parameter: the propertyID of the calculated property(i32)
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

| Properties                            |    Unit     | symbol | propertyID |
| ------------------------------------- | :---------: | -----: | ---------: |
| Pressure                              |     MPa     |      p |          0 |
| Temperature                           |     °C      |      t |          1 |
| Density                               |   kg/m^3    |      d |          2 |
| Specific Volume                       |   m^3/kg    |      v |          3 |
| Specific enthalpy                     |    kJ/kg    |      h |          4 |
| Specific entropy                      |  kJ/(kg·K)  |      s |          5 |
| Specific exergy                       |    kJ/kg    |      e |          6 |
| Specific internal energy              |    kJ/kg    |      u |          7 |
| Specific isobaric heat capacity       |  kJ/(kg·K)  |     cp |          8 |
| Specific isochoric heat capacity      |  kJ/(kg·K)  |     cv |          9 |
| Speed of sound                        |     m/s     |      w |         10 |
| Isentropic exponent                   |             |     ks |         11 |
| Specific Helmholtz free energy        |    kJ/kg    |      f |         12 |
| Specific Gibbs free energy            |    kJ/kg    |      g |         13 |
| Compressibility factor                |             |      z |         14 |
| Steam quality                         |             |      x |         15 |
| Region                                |             |      r |         16 |
| Isobaric volume expansion coefficient |     1/K     |     ec |         17 |
| Isothermal compressibility            |    1/MPa    |     kt |         18 |
| Partial derivative (dV/dT)p           |  m3/(kg·K)  |   dvdt |         19 |
| Partial derivative (dV/dP)T           | m3/(kg·MPa) |   dvdp |         20 |
| Partial derivative (dP/dT)v           |    MPa/K    |   dpdt |         21 |
| Isothermal Joule-Thomson coefficient  | kJ/(kg·MPa) |   iJTC |         22 |
| Joule-Thomson coefficient             |    K/MPa    |    JTC |         23 |
| Dynamic viscosity                     |  kg/(m·s)   |     dv |         24 |
| Kinematic viscosity                   |    m^2/s    |     kv |         25 |
| Thermal conductivity                  |   W/(m.K)   |     tc |         26 |
| Thermal diffusivity                   |   um^2/s    |     td |         27 |
| Prandtl number                        |             |     pr |         28 |
| Surface tension                       |    mN/m     |     st |         29 |
