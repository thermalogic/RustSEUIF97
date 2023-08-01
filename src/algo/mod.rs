//! Algorithms
//! * The secant method to find the root
//! * The function to sum the power of function and its derivative on IJn[(i32,i32,f64)

pub mod root;
pub mod polynomial;
pub mod polynomial_steps;

pub use self::root::*;
pub use self::polynomial::*;
pub use self::polynomial_steps::*;

