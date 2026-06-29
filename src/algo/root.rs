//! Find the root of the equation f(x)
//! Numerical Recipes  Ch.9
//!   * Bisection method
//!   * Secant method
//!   * Brent's method

type IF97_EQ = fn(f64, f64) -> f64;

pub const FIRST_FIXED: i32 = 1;        /* f(fvar, x): first parameter fixed, search second */
pub const SECOND_FIXED: i32 = 2;        /* f(x, fvar): second parameter fixed, search first */
pub const CONVERGENCE_PRECISION: f64 = 1.0E-08; /* rtsec xacc: convergence threshold for search variable */
pub const FN_TOLERANCE: f64 = 1.0E-08;   /* bisection tol: function value tolerance */
pub const INTERVAL_TOLERANCE: f64 = 1.0E-08;   /* bisection x_tol: interval length tolerance */
pub const MAX_ITER: i32 = 20000;     /* maximum number of iterations */
pub const ESP: f64 = 3.0E-08; /* machine floating-point precision */

/// Finds the root of the equation f(x) = 0 using the bisection method.
///
/// # Arguments
/// * `x1` - Left boundary of the search interval
/// * `x2` - Right boundary of the search interval
/// * `f` - Target function (returns f64, we need to find x where f(x) = 0)
/// * `max_iter` - Maximum number of iterations
/// * `tol` - Function value tolerance (|f(x)| < tol)
/// * `x_tol` - Interval length tolerance (|x1 - x2| < x_tol)
///
/// # Returns
/// Approximate root satisfying the precision requirement
pub fn bisection<F>(
    mut x1: f64,
    mut x2: f64,
    f: F,
    max_iter: usize,
    tol: f64,
    x_tol: f64,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let mut r_x1 = f(x1);
    let mut r_x2 = f(x2);
     if r_x1 * r_x2 > 0.0 {
         panic!("Bisection failed: f(x1) and f(x2) must have opposite signs!");
     }

    for _ in 0..max_iter {
        let xm = 0.5 * (x1 + x2);
        let r_xm = f(xm);
       if r_xm.abs() < tol || (x1 - x2).abs() < x_tol {
            return xm;
        }
       if r_x1 * r_xm < 0.0 {
            x2 = xm;
            r_x2 = r_xm;
        } else {
            x1 = xm;
            r_x1 = r_xm;
        }
    }
    0.5 * (x1 + x2)
}

/// Finds the root of the equation fun(f64, f64) = target using the secant method.
/// Numerical Recipes  Ch.9.2
/// # Arguments
/// * `fun` - Target function fun(f64, f64) -> f64
/// * `fvar` - Fixed parameter value
/// * `target` - Target value (solve for fun(f64, f64) = target)
/// * `x1` - Left boundary of the search interval
/// * `x2` - Right boundary of the search interval
/// * `fvar_position` - Position of the fixed parameter: 1=fun(fvar, x), 2=fun(x, fvar)
/// * `i_max` - Maximum number of iterations
/// * `xacc` - Convergence precision
///
/// # Returns
///   Approximate root satisfying the precision requirement
pub fn rtsec(fun: IF97_EQ, fvar: f64, target: f64, x1: f64, x2: f64,
            fvar_position: i32, i_max: i32, xacc: f64) -> f64 {
    let mut xl: f64;
    let mut rts: f64;
    let mut swap: f64;
    let mut dx: f64 = 0.0;
    
    // Helper closure to evaluate f(x) = target - fun(x)
    let mut eval = |x: f64| -> f64 {
        if fvar_position == 1 {
            target - fun(fvar, x)
        } else {
            target - fun(x, fvar)
        }
    };
    
    // Calculate function values
    let (mut fl, mut f) = (eval(x1), eval(x2));
    
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
            // Prevent non-positive root (physical quantities like T/P must be positive)
            if rts <= 0.0 {
                 rts = 0.000001;
             }
            // Calculate function value
            f = eval(rts);
            i += 1;
        }
    };
    rts
}


/// Brent's method: combines bisection, secant, and inverse quadratic interpolation
/// to find the root of fun(f64, f64) = target within bracket [x1, x2].
///
/// # Arguments
/// * `fun` - Target function fn(f64, f64) -> f64
/// * `fvar` - Fixed parameter value
/// * `target` - Target value (solve for fun(f64, f64) = target)
/// * `x1` - Left boundary of the search interval
/// * `x2` - Right boundary of the search interval
/// * `fvar_position` - Position of the fixed parameter: 1=fun(fvar,x), 2=fun(x,fvar)
/// * `i_max` - Maximum number of iterations
/// * `tol` - Convergence precision (on x)
///
/// # Returns
///   Approximate root satisfying the precision requirement
///
/// # Panics
///   Panics if f(x1) and f(x2) do not have opposite signs (no bracket).
pub fn zbrent(
    fun: IF97_EQ,
    fvar: f64,
    target: f64,
    x1: f64,
    x2: f64,
    fvar_position: i32,
    i_max: i32,
    tol: f64,
) -> f64 {
    // Helper closure to evaluate f(x) = target - fun(x)
    let mut eval = |x: f64| -> f64 {
        if fvar_position == 1 {
            target - fun(fvar, x)
        } else {
            target - fun(x, fvar)
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
        let tol1 = 2.0 * ESP * b.abs() + 0.5 * tol;
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