//!     
//! Region 2 - (p,T),(p,h), (p,s),(h,s)
//! 
mod region2_gfe;
pub mod region2_pT;
pub mod region2_p_hs;
mod region2_T_ph;
pub mod region2_T_ps;
pub mod region2;

pub use self::region2_pT::*;
pub use self::region2_p_hs::*;
pub use self::region2_T_ph::*;
pub use self::region2_T_ps::*;
pub use self::region2::*;
