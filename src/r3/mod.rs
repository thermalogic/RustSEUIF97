//! Region 3 - (p,T),(p,h), (p,s),(h,s); (p,v),(t,v),(t,h),(t,s)
pub mod region3;
pub mod region3_Td;
pub mod region3_Tv_phps;
mod region3_hfe;
pub mod region3_p_hs;
pub mod region3_pair_ext;
pub mod region3_v_pT;
pub mod region3_v_subregion_pT;
pub mod region3_Td_ext;

pub use self::region3::*;
pub use self::region3_Td::*;
pub use self::region3_Tv_phps::*;
pub use self::region3_p_hs::*;
pub use self::region3_pair_ext::*;
pub use self::region3_v_pT::*;
pub use self::region3_v_subregion_pT::*;
pub use self::region3_Td_ext::*;
