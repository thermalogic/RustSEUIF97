//! Region 5 - (p,T),(p,h), (p,s),(h,s); (p,v),(t,v),(t,h),(t,s)

pub mod region5;
mod region5_gfe;
pub mod region5_pT;
pub mod region5_pT_ext;
pub mod region5_pair_ext;
pub mod region5_ph_ps_hs;

pub use self::region5::*;
pub use self::region5_pT::*;
pub use self::region5_pT_ext::*;
pub use self::region5_pair_ext::*;
pub use self::region5_ph_ps_hs::*;
