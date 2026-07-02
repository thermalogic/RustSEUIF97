#![allow(warnings)]
// allow snake case for using the thermodynamics notation
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]  
#![allow(clippy::approx_constant)]
#![doc=include_str!("../README.md")]
//#![warn(missing_docs)]

mod algo;
mod common;
mod if97_core;
pub mod r1;
pub mod r2;
mod r3;
mod r4;
pub mod r5;

pub use common::property_id::*;
use common::*;
pub use r1::*;
pub use r2::*;
use r3::*;
use r4::*;
pub use r5::*;

mod rust_if97;
pub use rust_if97::*;

#[cfg(feature = "cdecl")]
pub mod cdecl_c_if97;

#[cfg(feature = "stdcall")]
pub mod stdcall_c_if97;

#[cfg(feature = "python")]
pub mod python_if97;

#[cfg(feature = "wasm")]
pub mod wasm_if97;
