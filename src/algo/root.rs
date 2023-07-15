//! The secant method to find the root
//！
//！Numerical Reciples  Ch.9.2: Pages 357

pub type IF97_EQ = fn(f64,f64) -> f64;

pub const ESP:f64=1.0E-08;
pub const I_MAX:i32=100;

///  rtsec1 :fun(x,con_var2)： the result is the first x：x，x={x1,x2}
///  the second con_var2 in constant 
///    fr: the value of fun(x,con_var2)
///    fl: fr-fun(x1,con_var2)
///    f: fr-fun(x2,con_var2)
/// rts：the solution of x
pub fn  rtsec1(fun:IF97_EQ, con_var2:f64, fr:f64,x1:f64, x2:f64, mut fl:f64,mut f:f64, xacc:f64,i_max:i32)->f64
{
   let mut xl:f64;
   let mut rts:f64;
   let mut swap:f64;
   let mut dx:f64=0.0;
  // pick the bound with the smaller function value as the most recent guess
  if fl.abs() < f.abs() { 
      rts=x1;
      xl=x2;
      swap=fl;
      fl=f;
      f=swap;
  } 
  else {
    xl=x1;
    rts=x2;
  };

//scant loop
    let mut i:i32=0;    
    if (f-fl)!=0.0
    {
      dx=(xl-rts)*f/(f-fl); 
      while (dx.abs() > xacc)&&(i<i_max)&&(f!=0.0)&&((f-fl)!=0.0)
      {
        dx=(xl-rts)*f/(f-fl); // increment with respect to latest value
        xl=rts;
        fl=f;
        rts += dx;
        //rts may out-bounded in region X
        if rts<=0.0
        { 
          rts=0.000001;
        };
        //if (rts>100) rts=100;
        f=fr-fun(rts,con_var2);
        i+=1;
       };
   };
   return rts;
}

///  rtsec2 :fun(con_var1，x),the result is the second x，x={x1,x2}
///     the first con_var1 is constant
///     fr: the value of func(con_var1,x)
///     fl: fr-func(con_var1,x1)
///     f: fr-func(con_var1,x2)
///     rts：the solution of x
pub fn  rtsec2(fun:IF97_EQ, con_var1:f64, fr:f64, x1:f64, x2:f64,mut fl:f64,mut f:f64, xacc:f64, i_max:i32)->f64
{
  let mut xl:f64;
  let mut rts:f64;
  let mut swap:f64;
  let mut dx:f64=0.0;
  // pick the bound with the smaller function value as the most recent guess
  if fl.abs() <f.abs() {
        rts=x1;
        xl=x2;
        swap=fl;
        fl=f;
        f=swap;
  } else {
        xl=x1;
        rts=x2;
  }
  //scant loop
  let mut i:i32=0;
  if (f-fl)!=0.0
  {  dx=(xl-rts)*f/(f-fl);
    while (dx.abs() > xacc) && (i<i_max) && (f!=0.0) && (f-fl)!=0.0  // Covergence
      {
          dx=(xl-rts)*f/(f-fl); // increment with respect to latest value
          xl=rts;
          fl=f;
          rts += dx;
          //rts must bounded in region X
          f=fr-fun(con_var1,rts);
          i+=1;
      }
  }
  return rts;
}

