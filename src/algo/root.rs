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

// Brent's method: combines bisection, secant, and inverse quadratic interpolation
/// to find the root of fun(x) = target within bracket [x1, x2].
///
/// # Arguments
/// * `fun` - Target function fn(f64, f64) -> f64
/// * `var` - Fixed parameter value
/// * `target` - Target value (solve for fun = target)
/// * `x1` - Left boundary of the search interval
/// * `x2` - Right boundary of the search interval
/// * `var_position` - Position of the fixed parameter: 1=fun(var,x), 2=fun(x,var)
/// * `tol` - Convergence precision (on x)
/// * `i_max` - Maximum number of iterations
///
/// # Returns
/// Approximate root satisfying the precision requirement
///
/// # Panics
/// Panics if f(x1) and f(x2) do not have opposite signs (no bracket).
pub fn zbrent(
    fun: IF97_EQ,
    var: f64,
    target: f64,
    x1: f64,
    x2: f64,
    var_position: i32,
    tol: f64,
    i_max: i32,
) -> f64 {
    // Helper closure to evaluate f(x) = target - fun(x)
    let mut eval = |x: f64| -> f64 {
        if var_position == 1 {
            target - fun(var, x)
        } else {
            target - fun(x, var)
        }
    };

    let mut a = x1;
    let mut b = x2;
    let mut c = x2;

    let mut fa = eval(a);
    let mut fb = eval(b);
    let mut fc = fb;

    // Verify bracket
    if fa * fb > 0.0 {
        panic!("zbrent: root must be bracketed");
    }

    let mut d = b - a;
    let mut e = d;

    let mut s: f64;
    let mut p: f64;
    let mut q: f64;
    let mut r: f64;

    for _ in 0..i_max {
        // Ensure b is the best estimate so far, and c is the previous value of b
        if fb * fc > 0.0 {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }

        // Reorder so that |fc| > |fb| > |fa| is NOT maintained;
        // we want b to be the best root estimate (smallest |f|)
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        // Convergence criterion on x
        let tol1 = 2.0 * f64::EPSILON * b.abs() + 0.5 * tol;
        let xm = 0.5 * (c - b);

        if xm.abs() <= tol1 || fb == 0.0 {
            return b;
        }

        // Attempt interpolation if we have a good bracket
        if xm.abs() < e.abs() && fa != fb {
            let mut use_secant = false;

            // Attempt inverse quadratic interpolation
            if a == c {
                // Secant method (only two points distinct)
                s = fb / fa;
                p = 2.0 * xm * s;
                q = 1.0 - s;
                use_secant = true;
            } else {
                // Inverse quadratic interpolation
                q = fa / fc;
                r = fb / fc;
                s = fb / fa;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if p > 0.0 {
                q = -q;
            }
            p = p.abs();

            let min1 = 3.0 * xm * q - (tol1 * q).abs();
            let min2 = (e * q).abs();

            // Accept interpolation only if it falls within the bracket
            // and is making sufficient progress
            if 2.0 * p < min1.min(min2) {
                e = d;
                d = p / q;
            } else {
                // Fall back to bisection
                d = xm;
                e = d;
            }
        } else {
            // Bounds decreasing too slowly, use bisection
            d = xm;
            e = d;
        }

        // Update a to be the previous best estimate
        a = b;
        fa = fb;

        // Compute new estimate
        if d.abs() > tol1 {
            b += d;
        } else {
            b += tol1.copysign(xm);
        }

        // Prevent non-positive values if needed
        if b <= 0.0 {
            b = 0.000001;
        }

        fb = eval(b);
    }

    // Max iterations reached
    b
}