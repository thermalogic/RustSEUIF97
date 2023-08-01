//! Region 1 -  (p,T),(p,h),(p,s),(h,s); (p,v),(t,v),(t,h),(t,s)
pub mod region1;
pub mod region1_T_phps;
mod region1_gfe;
pub mod region1_pT;
pub mod region1_pT_ext;
pub mod region1_p_hs;
pub mod region1_pair_ext;

pub use self::region1::*;
pub use self::region1_T_phps::*;
pub use self::region1_pT::*;
pub use self::region1_pT_ext::*;
pub use self::region1_p_hs::*;
pub use self::region1_pair_ext::*;
