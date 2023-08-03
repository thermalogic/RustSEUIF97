#![allow(warnings)]
//#![warn(missing_docs)]
/*!
# IF97

IF97 is the high-speed package of IAPWS-IF97 in Rust. It is suitable for computation-intensive calculations，such as heat cycle calculations, simulations of non-stationary processes, real-time process monitoring and optimizations.

Through the high-speed package, the results of the IAPWS-IF97 are accurately produced at about 5-20x speed-up compared to  using the `powi()` of
the Rust standard library in the `for`loop directly when computing the basic equations of Region 1,2,3.

* The comparison results of the computing-time are obtained using the [criterion.rs](https://bheisler.github.io/criterion.rs/book/index.html).

**The Fast Methods**

1. The multi-step method unleashes the full power of the compiler optimizations while using `powi()` with the `for` loop
2. The recursive method computes the polynomial values of the base variable and its derivatives

In IF97, [36 thermodynamic, transport and  further properties](#properties) can be calculated.

The following input pairs are implemented:

```txt
(p,t) (p,h) (p,s) (p,v)

(t,h) (t,s) (t,v)

(p,x) (t,x) (h,x) (s,x)

(h,s)
```
## Release Notes


* [1.1.7](https://crates.io/crates/if97/1.1.7) -  The multi-step method unleashes the full power of the compiler optimizations while using `powi()` with the `for` loop
   * Add the optional parameter of `region` to computer the properties of the specified region quickly

* [1.0.9](https://crates.io/crates/if97/1.0.9) - The shortest addition chain computes integer powers of a number to provide the stable speed performance under the various hardware and software platforms

## Usage

The type of functions are provided in the if97 package:

```rust
struct o_id_region_args {
    o_id: i32,
    region: i32,
}

fn<R>(f64,f64,R) -> f64
where
    R: Into<o_id_region_args>,
``````

* the first,second input parameters(f64) : the input propertry pairs
* the third and fourth input parametes<R>:
   * the third : the property ID of the calculated property - [o_id](#properties)
   * the fourth `option` parameter: the region of IAPWS-IF97
* the return(f64): the calculated property value of o_id

```rust
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
**Example**

```rust
use if97::*;
fn main() {

    let p:f64 = 3.0;
    let t:f64= 300.0-273.15;

    let h=pt(p,t,OH);
    let s=pt(p,t,OS);
    // set the region
    let v=pt(p,t,(OV,1));
    println!("p={p:.6} t={t:.6} h={t:.6} s={s:.6} v={v:.6}");

}
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

*/

mod algo;
mod common;
mod r1;
mod r2;
mod r3;
mod r4;
mod r5;

pub use common::propertry_id::*;
use common::*;
use r1::*;
use r2::*;
use r3::*;
use r4::*;
use r5::*;

/// the paramters: <br/>
///   `o_id`: the propertry of id;<br/>
///   `region`: the region in the IAPWS-IF97,**optional**
pub struct o_id_region_args {
    o_id: i32,
    region: i32,
}

impl Default for o_id_region_args {
    fn default() -> Self {
        o_id_region_args { o_id: 0, region: 6 }
    }
}

impl From<i32> for o_id_region_args {
    fn from(o_id: i32) -> Self {
        Self {
            o_id: o_id,
            ..Self::default()
        }
    }
}

impl From<(i32, i32)> for o_id_region_args {
    fn from((o_id, region): (i32, i32)) -> Self {
        Self {
            o_id: o_id,
            region: region,
        }
    }
}

fn o_id_reg_info<R>(o_id_reg: R) -> (i32, i32)
where
    R: Into<o_id_region_args>,
{
    let args = o_id_reg.into();
    let o_id: i32 = args.o_id;
    let reg: i32 = args.region;
    (o_id, reg)
}

/// pt(p,t,o_id) - the propertry of `o_id` (thermodynamic,transport,etc) <br/>
/// pt(p,t,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///
/// # Examples
///
/// ```
/// use if97::*;
///
/// let p:f64 = 3.0;
/// let t:f64= 300.0-273.15;
/// let h=pt(p,t,OH);
/// // set the region
/// let s=pt(p,t,(OS,1));
/// println!("p={p:.6} t={t:.6} h={h:.6} s={s:.6}");    
/// ```
///
pub fn pt<R>(p: f64, t: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    let T: f64 = t + 273.15;
    match o_id {
        OP => return p,
        OT => return t,
        OST | ODV | OKV | OTC | OSDC | OPR | OTD => {
            pair_transport(p, T, o_id, PAIRS::pT, pT_thermal, reg)
        }
        _ => pT_thermal(p, T, o_id, reg),
    }
}

/// ph(p,h,o_id) - the propertry of `o_id` (thermodynamic,transport,etc)<br/>
/// ph(p,h,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///
/// # Examples
///
/// ```
/// use if97::*;
///
/// let p:f64 = 3.0;
/// let h:f64= 0.115331273e+3;
/// let t=ph(p,h,OT);
/// // set the region
/// let s=ph(p,h,(OS,1));
/// println!("p={p:.6} h={h:.6} t={t:.6} s={s:.6}");    
/// ```
///
pub fn ph<R>(p: f64, h: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OP => return p,
        OH => return h,
        OST | ODV | OKV | OTC | OSDC | OPR | OTD => {
            pair_transport(p, h, o_id, PAIRS::ph, ph_thermal, reg)
        }
        _ => ph_thermal(p, h, o_id, reg),
    }
}

/// ps(p,s,o_id) - the propertry of `o_id` (thermodynamic,transport,etc)<br/>
/// ps(p,s,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
/// # Examples
///
///```
/// use if97::*;
///
/// let p:f64= 3.0;
/// let s:f64= 0.392294792;
/// let t=ps(p,s,OT);
/// // set the region
/// let h=ps(p,s,(OH,1));
/// println!("p={p:.6} s={s:.6} t={t:.6} h={h:6}");    
/// ```
pub fn ps<R>(p: f64, s: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OP => return p,
        OS => return s,
        OST | ODV | OKV | OTC | OSDC | OPR | OTD => {
            pair_transport(p, s, o_id, PAIRS::ps, ps_thermal, reg)
        }
        _ => ps_thermal(p, s, o_id, reg),
    }
}

/// hs(h,s,o_id) - the propertry of `o_id` (thermodynamic,transport,etc)<br/>
/// hs(h,s,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///
/// # Examples
///
///```
/// use if97::*;
///
/// let h:f64= 0.115331273e+3;;
/// let s:f64= 0.392294792;
/// let p=hs(h,s,OP);
/// // set the region
/// let t=hs(h,s,(OT,1));
/// println!("h={h:.6} s={s:.6} p={p:.6} t={t:.6}");
/// ```
///   
pub fn hs<R>(h: f64, s: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OH => return h,
        OS => return s,
        OST | ODV | OKV | OTC | OSDC | OPR | OTD => {
            pair_transport(h, s, o_id, PAIRS::hs, hs_thermal, reg)
        }
        _ => hs_thermal(h, s, o_id, reg),
    }
}

///  px(p,x,o_id) - the propertry of `o_id` (thermodynamic)
///
///  # Examples
///
///```
///  use if97::*;
///
///  let p: f64 = 0.1;
///  let x: f64 = 0.0; // x= 0.0 saturation water ,x=1.0 saturation steam
///  let t: f64 = px(p, x, OT);
///  let h: f64 = px(p, x, OH);
///  println!("px: p={p:.6} x={x:.6} t={t:.6} h={h:.6}");
///```

// #[no_mangle]
// pub unsafe extern "C"  fn px(p: f64, x: f64, o_id: i32) -> f64 {
pub fn px(p: f64, x: f64, o_id: i32) -> f64 {
    if p > P_MAX4 || p < P_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OP => return p,
        OX => return x,
        OT => return px_reg4(p, x, o_id) - 273.15,
        _ => return px_reg4(p, x, o_id),
    }
}

///  tx(t,x,o_id) - the propertry of `o_id` (thermodynamic)
///
///  # Examples
///
///```
///  use if97::*;
///
///  let t: f64 =372.755919-273.15;
///  let x: f64 = 0.0; // x= 0.0 saturation water ,x=1.0 saturation steam
///  let p: f64 = tx(t, x, OP);
///  let h: f64 = tx(t, x, OH);
///  println!("tx: p={p:.6} x={x:.6} t={t:.6} h={h:.6}");
/// ```
///
pub fn tx(t: f64, x: f64, o_id: i32) -> f64 {
    let T: f64 = t + 273.15;
    if T > T_MAX4 || T < T_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    Tx_reg4(T, x, o_id)
}

//  Functions of the  extended input pairs (p,v),(t,v),(t,s),(t,h)

/// pv(p,v,o_id) - the propertry of `o_id`(thermodynamic,transport,etc)<br/>
/// pv(p,v,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///
/// # Examples
///
///```
/// use if97::*;
///
/// let p:f64= 3.0;
/// let v:f64= 0.100215168e-2;
/// let t=pv(p,v,OT);
/// //set the region
/// let h=pv(p,v,(OH,1));
/// println!("p={p:.6} v={v:.6} t={t:.6}");    
/// ```
pub fn pv<R>(p: f64, v: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OP => return p,
        OV => return v,
        OST | ODV | OKV | OTC | OSDC | OPR | OTD => {
            pair_transport(p, v, o_id, PAIRS::pv, pv_thermal, reg)
        }
        _ => pv_thermal(p, v, o_id, reg),
    }
}

/// tv(t,v,o_id) - the propertry of `o_id`(thermodynamic,transport,etc)<br/>
/// tv(t,v,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///
/// # Examples
///
///```
/// use if97::*;
///
/// let t:f64=300.0-273.15;
/// let v:f64= 0.100215168e-2;
/// let p=tv(t,v,OP);
/// //set the regiion
/// let s=tv(t,v,(OS,1));
/// println!("t={p:.6} v={v:.6} p={p:.6} s={s:.6}");    
/// ```
pub fn tv<R>(t: f64, v: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    match o_id {
        OT => return t,
        OV => return v,
        OST | ODV | OKV | OTC | OSDC | OPR | OTD => {
            pair_transport(t, v, o_id, PAIRS::tv, tv_thermal, reg)
        }
        _ => tv_thermal(t, v, o_id, reg),
    }
}

/// th(t,h,o_id) - the propertry of `o_id`(thermodynamic,transport,etc)<br/>
/// th(t,h,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///  
/// # Examples
///
///```
/// use if97::*;
///
/// let t:f64=300.0-273.15;
/// let h:f64=0.115331273e+3;
/// let p=th(t,h,OP);
/// // set the region
/// let s=th(t,h,(OS,1));
/// println!("t={p:.6} h={h:.6} p={p:.6} s={s:.6}");    
/// ```
pub fn th<R>(t: f64, h: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id, reg) = o_id_reg_info(o_id_reg);
    let T: f64 = t + 273.15;
    match o_id {
        OT => return t,
        OH => return h,
        OST | ODV | OKV | OTC | OSDC | OPR | OTD => {
            pair_transport(t, h, o_id, PAIRS::th, th_thermal, reg)
        }
        _ => th_thermal(t, h, o_id, reg),
    }
}

/// ts(t,s,o_id) - the propertry of `o_id`(thermodynamic,transport,etc)<br/>
/// ts(t,s,(o_id,reg)) - the propertry of `o_id` in the region of `reg`(thermodynamic,transport,etc)
///   
/// # Examples
///
///```
/// use if97::*;
///
/// let t:f64=300.0-273.15;
/// let s:f64=0.392294792;
/// let p=ts(t,s,OP);
/// // set the region
/// let h=ts(t,s,(OH,1));
/// println!("t={p:.6} s={s:.6} p={p:.6} h={h:.6}");    
/// ```
pub fn ts<R>(t: f64, s: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
{
    let (o_id,reg)=o_id_reg_info( o_id_reg);
    let T: f64 = t + 273.15;
    match o_id {
        OT => return t,
        OS => return s,
        OST | ODV | OKV | OTC | OSDC | OPR | OTD => {
            pair_transport(t, s, o_id, PAIRS::ts, ts_thermal, reg)
        }
        _ => ts_thermal(t, s, o_id, reg),
    }
}

/// hx(h,x,o_id) - the propertry of `o_id`(thermodynamic)
///
///  # Examples
///
///```
///  use if97::*;
///
///  let h: f64 =1094.690434;
///  let x: f64 = 0.3;
///  let p: f64 = hx(h, x, OP);
///  let t: f64 = hx(h, x, OT);
///  println!("hx: h={p:.6} x={x:.6} p={p:.6} t={t:.6}");
/// ```
///
pub fn hx(h: f64, x: f64, o_id: i32) -> f64 {
    if h > H_MAX4 || h < H_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OH => return h,
        OX => return x,
        _ => hx_reg4(h, x, o_id),
    }
}

/// sx(s,x,o_id) - the propertry of `o_id`(thermodynamic)
///
///  # Examples
///
///```
///  use if97::*;
///
///  let s: f64 =3.119434;
///  let x: f64 = 0.3;
///  let p: f64 = sx(s, x, OP);
///  let t: f64 = sx(s, x, OT);
///  println!("sx: s={p:.6} x={x:.6} p={p:.6} t={t:.6}");
/// ```
///
pub fn sx(s: f64, x: f64, o_id: i32) -> f64 {
    if s > S_MAX4 || s < S_MIN4 || x > 1.0 || x < 0.0 {
        return INVALID_VALUE as f64;
    }
    match o_id {
        OS => return s,
        OX => return x,
        _ => sx_reg4(s, x, o_id),
    }
}
