//! The secant method to find the root
//! * Numerical Reciples  Ch.9.2

type IF97_EQ = fn(f64, f64) -> f64;

pub const ESP: f64 = 1.0E-08;
pub const I_MAX: i32 = 100;

/// Finds the root of the equation fun(x) = target using the secant method.
///
/// # Arguments
/// * `fun` - Target function fn(f64, f64) -> f64
/// * `var` - Fixed parameter value
/// * `target` - Target value (solve for fun = target)
/// * `x1` - Left boundary of the search interval
/// * `x2` - Right boundary of the search interval
/// * `var_position` - Position of the fixed parameter: 1=fun(var,x), 2=fun(x,var)
/// * `xacc` - Convergence precision
/// * `i_max` - Maximum number of iterations
///
/// # Returns
/// Approximate root satisfying the precision requirement
pub fn rtsec(
    fun: IF97_EQ, var: f64, target: f64, x1: f64, x2: f64,
    var_position: i32, xacc: f64, i_max: i32,
) -> f64 {
    let mut xl: f64;
    let mut rts: f64;
    let mut swap: f64;
    let mut dx: f64 = 0.0;
    
    // Calculate function values based on variable position
    let (mut fl, mut f) = if var_position == 1 {
        // fun(var, x) - first parameter is fixed
        (target - fun(var, x1), target - fun(var, x2))
    } else {
        // fun(x, var) - second parameter is fixed
        (target - fun(x1, var), target - fun(x2, var))
    };
    
    // pick the bound with the smaller function value as the most recent guess
    if fl.abs() < f.abs() {
        rts = x1;
        xl = x2;
        swap = fl;
        fl = f;
        f = swap;
    } else {
        xl = x1;
        rts = x2;
    };

    // secant loop
    let mut i: i32 = 0;
    if (f - fl) != 0.0 {
        dx = (xl - rts) * f / (f - fl);
        while (dx.abs() > xacc) && (i < i_max) && (f != 0.0) && ((f - fl) != 0.0) {
            dx = (xl - rts) * f / (f - fl); // increment with respect to latest value
            xl = rts;
            fl = f;
            rts += dx;
            // rts <=0.0 out of bounds
            if rts <= 0.0 {
                 rts = 0.000001;
             }
            // Calculate function value based on variable position
            if var_position == 1 {
                f = target - fun(var, rts);
            } else {
                f = target - fun(rts, var);
            }
            i += 1;
        }
    };
    
    rts
}


/// Finds the root of the equation f(x) = 0 using the bisection method.
///
/// # Arguments
/// * `t1` - Left boundary of the search interval
/// * `t2` - Right boundary of the search interval
/// * `f` - Target function (returns f64, we need to find x where f(x) = 0)
/// * `max_iter` - Maximum number of iterations
/// * `tol` - Function value tolerance (|f(x)| < tol)
/// * `x_tol` - Interval length tolerance (|t1 - t2| < x_tol)
///
/// # Returns
/// Approximate root satisfying the precision requirement
pub fn bisection<F>(
    mut t1: f64,
    mut t2: f64,
    f: F,
    max_iter: usize,
    tol: f64,
    x_tol: f64,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let mut r_t1 = f(t1);
    let mut r_t2 = f(t2);
     if r_t1 * r_t2 > 0.0 {
         panic!("Bisection failed: f(t1) and f(t2) must have opposite signs!");
     }

    for _ in 0..max_iter {
        let tm = 0.5 * (t1 + t2);
        let r_tm = f(tm);
       if r_tm.abs() < tol || (t1 - t2).abs() < x_tol {
            return tm;
        }
       if r_t1 * r_tm < 0.0 {
            t2 = tm;
            r_t2 = r_tm;
        } else {
            t1 = tm;
            r_t1 = r_tm;
        }
    }
    0.5 * (t1 + t2)
}