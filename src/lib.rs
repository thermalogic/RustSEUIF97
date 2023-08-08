#![allow(warnings)]
// allow snake case for using the thermodynamics notation
#![allow(non_snake_case)]
#![allow(clippy::approx_constant)]
#![doc=include_str!("../README.md")]
//#![warn(missing_docs)]

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

mod rust_if97;
pub use rust_if97::*;

#[cfg(feature = "cdecl")]
pub mod cdecl_c_if97;

#[cfg(feature = "stdcall")]
pub mod stdcall_c_if97;
