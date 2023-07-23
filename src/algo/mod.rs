//! Algorithms
//! * The fast integer power using the shortest addition chains
//! * The secant method to find the root
//! * The helper function to sum the power of function and it's derivative on IJn[(i32,i32,f64)
pub mod fast_ipow;
pub mod root;
pub mod sum_product_power;

pub use self::fast_ipow::*;
pub use self::root::*;
pub use self::sum_product_power::*;
