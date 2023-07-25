#![allow(warnings)]
//#![warn(missing_docs)]
/*!
# IF97

IF97 is the high-speed package of IAPWS-IF97 in Rust. It is suitable for computation-intensive calculations，such as heat cycle calculations, simulations of non-stationary processes, real-time process monitoring and optimizations.
 
Through the high-speed package, the results of the IAPWS-IF97 are accurately produced at about 5-15x speed-up compared to  the `powi()` of the Rust standard library.

* The comparison results of the computing-time are obtained using the [criterion.rs](https://bheisler.github.io/criterion.rs/book/index.html). 

**The Fast Algorithm**

1. The shortest addition chain computes integer powers of a number.([the paper in chinese](https://github.com/thermalogic/SEUIF97/blob/master/doc/水和水蒸汽热力性质IAPWS-IF97公式的通用计算模型.pdf))
2. The recursive algorithm computes the polynomial values of the base variable and it's derivatives

In IF97, [36 thermodynamic, transport and  further properties](#properties) can be calculated. 

The following input pairs are implemented: 

```txt
(p,t) (p,h) (p,s) (p,v) 

(t,h) (t,s) (t,v) 

(p,x) (t,x) (h,x) (s,x) 

(h,s)  
```
## Usage

The type of functions are provided in the if97 package:

```rust
fn(f64,f64,i32) -> f64;
``````

* the first,second input parameters: the input propertry pairs
* the third input parameter: the property ID of the calculated property - [o_id](#properties)
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
**Example**

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


/// pt(p,t,o_id) - the propertry of o_id(thermodynamic,transport,etc)
///
/// # Examples
///
/// ```
/// use if97::*;
///
/// let p:f64 = 3.0;
/// let t:f64= 300.0-273.15;
/// let h=pt(p,t,OH);
/// let s=pt(p,t,OS);
/// let v=pt(p,t,OV);
/// println!("p={p:.6} t={t:.6} h={h:.6} s={s:.6} v={v:.6}");    
/// ```
///
pub fn pt(p: f64, t: f64, o_id: i32) -> f64 {
    match o_id {
        OP => return p,
        OT => return t,
        OST => return surface_tension(t + K),
        ODV | OKV | OTC | OSDC => {
            let d: f64 = pt_thermal(p, t, OD);
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV {
                    value /= d; // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                }
            };
            value
        }
        OPR | OTD => {
            let d: f64 = pt_thermal(p, t, OD);
            let cp: f64 = pt_thermal(p, t, OCP);
            let tc: f64 = thcond(d, t + 273.15);
            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }
        _ => pt_thermal(p, t, o_id),
    }
}

/// ph(p,h,o_id) - the propertry of o_id (thermodynamic,transport,etc)
///
/// # Examples
///
/// ```
/// use if97::*;
///
/// let p:f64 = 3.0;
/// let h:f64= 0.115331273e+3;
/// let t=ph(p,h,OT);
/// println!("p={p:.6} h={h:.6} t={t:.6}");    
/// ```
///
pub fn ph(p: f64, h: f64, o_id: i32) -> f64 {
    match o_id {
        OP => return p,
        OH => return h,
        OST => {
            let t: f64 = ph_thermal(p, h, OT);
            return surface_tension(t + K);
        }
        ODV | OKV | OTC | OSDC => {
            let d: f64 = ph_thermal(p, h, OD);
            let t: f64 = ph_thermal(p, h, OT);
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV {
                    value /= d; // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                }
            };
            value
        }
        OPR | OTD => {
            let d: f64 = ph_thermal(p, h, OD);
            let t: f64 = ph_thermal(p, h, OT);

            let cp: f64 = ph_thermal(p, h, OCP);
            let tc: f64 = thcond(d, t + 273.15);

            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }
        _ => ph_thermal(p, h, o_id),
    }
}

/// ps(p,s,o_id) - the propertry of o_id (thermodynamic,transport,etc)
///
/// # Examples
///
///```
/// use if97::*;
///
/// let p:f64= 3.0;
/// let s:f64= 0.392294792;
/// let t=ps(p,s,OT);
/// println!("p={p:.6} s={s:.6} t={t:.6}");    
/// ```
pub fn ps(p: f64, s: f64, o_id: i32) -> f64 {
    match o_id {
        OP => return p,
        OS => return s,
        OST => {
            let t: f64 = ps_thermal(p, s, OT);
            return surface_tension(t + K);
        }
        ODV | OKV | OTC | OSDC => {
            let d: f64 = ps_thermal(p, s, OD);
            let t: f64 = ps_thermal(p, s, OT);
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV {
                    value /= d; // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                };
            };
            value
        }
        OPR | OTD => {
            let d: f64 = ps_thermal(p, s, OD);
            let t: f64 = ps_thermal(p, s, OT);

            let cp: f64 = ps_thermal(p, s, OCP);
            let tc: f64 = thcond(d, t + 273.15);
            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }
        _ => ps_thermal(p, s, o_id),
    }
}

/// hs(h,s,o_id) - the propertry of o_id (thermodynamic,transport,etc)
///
/// # Examples
///
///```
/// use if97::*;
///
/// let h:f64= 0.115331273e+3;;
/// let s:f64= 0.392294792;
/// let p=hs(h,s,OP);
/// let t=hs(h,s,OT);
/// println!("h={h:.6} s={s:.6} p={p:.6} t={t:.6}");
/// ```
///   
pub fn hs(h: f64, s: f64, o_id: i32) -> f64 {
    match o_id {
        OH => return h,
        OS => return s,
        OST => {
            let t: f64 = hs_thermal(h, s, OT);
            return surface_tension(t + K);
        }
        ODV | OKV | OTC | OSDC => {
            let d: f64 = hs_thermal(h, s, OD);
            let t: f64 = hs_thermal(h, s, OT);
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV {
                    value /= d; // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                };
            };
            value
        }
        OPR | OTD => {
            let d: f64 = hs_thermal(h, s, OD);
            let t: f64 = hs_thermal(h, s, OT);

            let cp: f64 = hs_thermal(h, s, OCP);
            let tc: f64 = thcond(d, t + 273.15);
            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }

        _ => hs_thermal(h, s, o_id),
    }
}

///  px(p,x,o_id) - the propertry of o_id (thermodynamic)
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

///  tx(t,x,o_id) - the propertry of o_id (thermodynamic)
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

/// pv(p,v,o_id) - the propertry of o_id (thermodynamic,transport,etc)
///
/// # Examples
///
///```
/// use if97::*;
///
/// let p:f64= 3.0;
/// let v:f64= 0.100215168e-2;
/// let t=pv(p,v,OT);
/// println!("p={p:.6} v={v:.6} t={t:.6}");    
/// ```
pub fn pv(p: f64, v: f64, o_id: i32) -> f64 {
    match o_id {
        OP => return p,
        OV => return v,
        OST => {
            let t: f64 = pv_thermal(p, v, OT);
            return surface_tension(t + K);
        }
        ODV | OKV | OTC | OSDC => {
            let d: f64 = 1.0 / v;
            let t: f64 = pv_thermal(p, v, OT);
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV {
                    value /= d; // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                };
            };
            value
        }
        OPR | OTD => {
            let d: f64 = 1.0 / v;
            let t: f64 = pv_thermal(p, v, OT);

            let cp: f64 = pv_thermal(p, v, OCP);
            let tc: f64 = thcond(d, t + 273.15);
            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }
        _ => pv_thermal(p, v, o_id),
    }
}

/// tv(t,v,o_id) - the propertry of o_id (thermodynamic,transport,etc)
///
/// # Examples
///
///```
/// use if97::*;
///
/// let t:f64=300.0-273.15;
/// let v:f64= 0.100215168e-2;
/// let p=tv(t,v,OP);
/// println!("t={p:.6} v={v:.6} p={p:.6}");    
/// ```
pub fn tv(t: f64, v: f64, o_id: i32) -> f64 {
    match o_id {
        OT => return t,
        OV => return v,
        OST => return surface_tension(t + K),
        ODV | OKV | OTC | OSDC => {
            let d: f64 = 1.0 / v;
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV {
                    value /= d; // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                };
            };
            value
        }
        OPR | OTD => {
            let d: f64 = 1.0 / v;
            let cp: f64 = tv_thermal(t, v, OCP);
            let tc: f64 = thcond(d, t + 273.15);
            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }
        _ => tv_thermal(t, v, o_id),
    }
}

/// th(t,h,o_id) - the propertry of o_id(thermodynamic,transport,etc)
///  
/// # Examples
///
///```
/// use if97::*;
///
/// let t:f64=300.0-273.15;
/// let h:f64=0.115331273e+3;
/// let p=th(t,h,OP);
/// println!("t={p:.6} h={h:.6} p={p:.6}");    
/// ```
pub fn th(t: f64, h: f64, o_id: i32) -> f64 {
    match o_id {
        OT => return t,
        OH => return h,
        OST => return surface_tension(t + K),
        ODV | OKV | OTC | OSDC => {
            let d: f64 = th_thermal(t, h, OD);
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV {
                    value /= d; // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                };
            };
            value
        }
        OPR | OTD => {
            let d: f64 = th_thermal(t, h, OD);
            let cp: f64 = th_thermal(t, h, OCP);
            let tc: f64 = thcond(d, t + 273.15);
            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }
        _ => th_thermal(t, h, o_id),
    }
}

/// ts(t,s,o_id) - the propertry of o_id (thermodynamic,transport,etc)
///   
/// # Examples
///
///```
/// use if97::*;
///
/// let t:f64=300.0-273.15;
/// let s:f64=0.392294792;
/// let p=ts(t,s,OP);
/// println!("t={p:.6} s={s:.6} p={p:.6}");    
/// ```
pub fn ts(t: f64, s: f64, o_id: i32) -> f64 {
    match o_id {
        OT => return t,
        OS => return s,
        OST => return surface_tension(t + K),
        ODV | OKV | OTC | OSDC => {
            let d: f64 = ts_thermal(t, s, OD);
            let mut value: f64 = 0.0;
            if o_id == ODV || o_id == OKV {
                value = viscosity(d, t + 273.15);
                if o_id == OKV {
                    value /= d; // Kinematic viscosity=Dynamic viscosity/density
                }
            } else {
                if o_id == OTC {
                    value = thcond(d, t + 273.15)
                } else {
                    if o_id == OSDC {
                        value = static_dielectric(d, t + 273.15);
                    };
                };
            };
            value
        }
        OPR | OTD => {
            let d: f64 = ts_thermal(t, s, OD);
            let cp: f64 = ts_thermal(t, s, OCP);
            let tc: f64 = thcond(d, t + 273.15);
            let mut value: f64 = 0.0;
            if o_id == OTD {
                value = thermal_diffusivity(tc, cp, d);
            } else if o_id == OPR {
                let dv: f64 = viscosity(d, t + 273.15);
                value = prandtl_number(dv, cp, tc);
            }
            value
        }
        _ => ts_thermal(t, s, o_id),
    }
}

/// hx(h,x,o_id) - the propertry of o_id (thermodynamic)
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

/// sx(s,x,o_id) - the propertry of o_id (thermodynamic)
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
