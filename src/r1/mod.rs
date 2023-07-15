//!     
//! Region 1 - (p,T),(p,h), (p,s),(h,s)
//! 
mod region1_gfe;
pub mod region1_pT;
pub mod region1_p_hs;
pub mod region1_T_phps;
pub mod region1;

pub use self::region1_pT::*;
pub use self::region1_p_hs::*;
pub use self::region1_T_phps::*;
pub use self::region1::*;
