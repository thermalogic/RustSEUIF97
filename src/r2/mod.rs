//! Region 2 - (p,T),(p,h), (p,s),(h,s); (p,v),(t,v),(t,h),(t,s)
pub mod region2;
mod region2_T_ph;
pub mod region2_T_ps;
mod region2_gfe;
pub mod region2_pT;
pub mod region2_pT_ext;
pub mod region2_p_hs;
pub mod region2_pair_ext;

pub use self::region2::*;
pub use self::region2_T_ph::*;
pub use self::region2_T_ps::*;
pub use self::region2_pT::*;
pub use self::region2_pT_ext::*;
pub use self::region2_p_hs::*;
pub use self::region2_pair_ext::*;
