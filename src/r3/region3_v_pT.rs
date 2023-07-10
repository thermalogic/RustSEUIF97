//! Backward Equation for Region 3:
//!    IAPWS-IF97-S05rev 
//!     (p,T)->v 
//!        p  pressure in MPa
//!        T  temperature in K
//!         v  density in kg/m3

use crate::common::constant::*;
use crate::algo::fast_ipower::sac_pow;
use crate::r4::region4_sat_pT::*;
use crate::r3::region3_v_subregion_pT::*;

/// boundary is one of 3ab, 3cd, ...
/// p is pressure in MPa
/// returns the temperature at the boundary in K
pub fn T_atRegionBoundary(p:f64, boundary:String)->f64
{
    const I3ab:[i32;6] =[0, 0, 1, 2, -1, -2];
    const I3cd:[i32;5] = [0, 0, 1, 2, 3];
    const I3gh:[i32;6] = [0,0, 1, 2, 3, 4];
    const I3ij:[i32;6] = [0,0, 1, 2, 3, 4];
    const I3jk:[i32;6] = [0,0, 1, 2, 3, 4];
    const I3mn:[i32;5] = [0,0, 1, 2, 3];
    const I3op:[i32;6] = [0,0, 1, 2, -1, -2];
    const I3qu:[i32;5] = [0,0, 1, 2, 3];
    const I3rx:[i32;5] = [0,0, 1, 2, 3];
    const I3uv:[i32;5] = [0,0, 1, 2, 3];
    const I3wx:[i32;6] = [0,0, 1, 2, -1, -2];

    const n3ab:[f64;6] = [0.0, 0.154793642129415e4, -0.187661219490113e3, 0.213144632222113e2,
                     -0.191887498864292e4, 0.918419702359447e3];
    const n3cd:[f64;5] = [0.0, 0.585276966696349e3, 0.278233532206915e1, -0.127283549295878e-1,
                     0.159090746562729e-3];
    const n3gh:[f64;6]=[0.0, -0.249284240900418e5, 0.428143584791546e4, -0.269029173140130e3,
                     0.751608051114157e1, -0.787105249910383e-1];
    const n3ij:[f64;6]=[0.0, 0.584814781649163e3, -0.616179320924617, 0.260763050899562,
                     -0.587071076864459e-2, 0.515308185433082e-4];
    const n3jk:[f64;6]=[0.0, 0.617229772068439e3, -0.770600270141675e1, 0.697072596851896,
                     -0.157391839848015e-1, 0.137897492684194e-3];
    const n3mn:[f64;5]=[0.0, 0.535339483742384e3, 0.761978122720128e1, -0.158365725441648,
                     0.192871054508108e-2];
    const n3op:[f64;6]=[0.0, 0.969461372400213e3, -0.332500170441278e3, 0.642859598466067e2,
                     0.773845935768222e3, -0.152313732937084e4];
    const n3qu:[f64;5]=[0.0, 0.565603648239126e3, 0.529062258221222e1, -0.102020639611016,
                     0.122240301070145e-2];
    const n3rx:[f64;5]=[0.0, 0.584561202520006e3, -0.102961025163669e1, 0.243293362700452,
                     -0.294905044740799e-2];
    const n3uv:[f64;5]=[0.0, 0.528199646263062e3, 0.890579602135307e1, -0.222814134903755,
                     0.286791682263697e-2];
    const n3wx:[f64;6]=[0.0, 0.728052609145380e1, 0.973505869861952e2, 0.147370491183191e2,
                     0.329196213998375e3, 0.873371668682417e3];

    let mut T:f64 = 0.0; // temperature at the boundary in K
    match boundary.as_str() {
            "3ab"=>{ for i in 1..=5  {  T += n3ab[i] * sac_pow(p.ln(), I3ab[i]);}},
            "3op"=>{ for i in 1..=5  {  T += n3op[i] * sac_pow(p.ln(), I3op[i]);}},
            "3ef"=>{ T = 3.727888004 * (p - 22.064) + 647.096;},
            "3cd"=>{ for i  in 1..=4 {  T += n3cd[i] * sac_pow(p, I3cd[i]);}},
            "3gh"=>{  for i in 1..=5 {  T += n3gh[i] * sac_pow(p, I3gh[i]);}},
            "3ij"=>{  for i in 1..=5 {  T += n3ij[i] * sac_pow(p, I3ij[i]);}},
            "3jk"=>{  for i in 1..=5 {  T += n3jk[i] * sac_pow(p, I3jk[i]);}},
            "3mn"=>{  for i in 1..=4 {  T += n3mn[i] * sac_pow(p, I3mn[i]);}},
            "3qu"=>{  for i in 1..=4 {  T += n3qu[i] * sac_pow(p, I3qu[i]);}},
            "3rx"=>{  for i in 1..=4 {  T += n3rx[i] * sac_pow(p, I3rx[i]);}},
            "3uv"=>{  for i in 1..=4 {  T += n3uv[i] * sac_pow(p, I3uv[i]);}},
            "3wx"=>{  for i in 1..=5 {  T += n3wx[i] * sac_pow(p.ln(), I3wx[i]);}},
            _=>(),
    };
    return T
}

/// p is pressure in MPa
/// t is temperature in K
/// x = 0 is liquid, =1 is steam]  x : integer    Vapor quality [-]
/// returns density in kg/m3
pub fn pT2v_sat_reg3(p:f64, T:f64,x:f64)->f64
{
    let mut v:f64=0.0;
    let mut subRegion:char;
    // set the region
    // force the region if density of saturated water is asked for
    if x == 0.0
    {
        if p < 19.0088
        {
            subRegion = 'c';
        }
        else
        {
            if p < 21.0434
            {
                subRegion = 's';
            }
            else
            {
                if (p < 21.9316)
                 { subRegion = 'u';}
                else
                 {   subRegion = 'y';}
            }
        }
    }
    else
    {
        if (p < 20.5)
        {
            subRegion = 't';
        }
        else
        {
            if (p < 21.0434)
               { subRegion = 'r';}
            else if (p < 21.9009)
               { subRegion = 'x';}
            else
               { subRegion = 'z';}
        }
    }

    return pT2v_3subreg(p, T, subRegion);
}

/// sets the subregion in region3
pub fn sub_region3_pT(p:f64,T:f64)->char
{

    let mut subRegion:char;

    if (p > 40.0 && p <= 100.0)
    {
        let tB:f64 = T_atRegionBoundary(p, "3ab".to_string());
        if T > tB
         {   subRegion = 'b';}
        else
        {    subRegion = 'a';}
    }
    else if p > 25.0
    {
        let tBab:f64 = T_atRegionBoundary(p, "3ab".to_string());
        let tBcd:f64 = T_atRegionBoundary(p, "3cd".to_string());
        let tBef:f64 = T_atRegionBoundary(p, "3ef".to_string());

        if T <= tBcd
         {   subRegion = 'c';}
        else if T <= tBab
         {   subRegion = 'd';}
        else if (T <= tBef)
         {   subRegion = 'e';}
        else
          {  subRegion = 'f';}
    }
    else if p > 23.5
    {
        let tBcd:f64 = T_atRegionBoundary(p, "3cd".to_string());
        let tBef:f64 = T_atRegionBoundary(p, "3ef".to_string());
        let tBgh:f64 = T_atRegionBoundary(p, "3gh".to_string());
        let tBij:f64 = T_atRegionBoundary(p, "3ij".to_string());
        let tBjk:f64 = T_atRegionBoundary(p, "3jk".to_string());

        if T <= tBcd
         {  subRegion = 'c';}
        else if T <= tBgh
          {  subRegion = 'g';}
        else if (T <= tBef)
           { subRegion = 'h';}
        else if (T <= tBij)
            { subRegion = 'i';}
        else if (T <= tBjk)
            { subRegion = 'j';}
        else
             { subRegion = 'k';}
    }
    else if p > 23.0
    {
        let tBcd:f64 = T_atRegionBoundary(p, "3cd".to_string());
        let tBef:f64 = T_atRegionBoundary(p, "3ef".to_string());
        let tBgh:f64 = T_atRegionBoundary(p, "3gh".to_string());
        let tBij:f64 = T_atRegionBoundary(p, "3ij".to_string());
        let tBjk:f64 = T_atRegionBoundary(p, "3jk".to_string());

        if T <= tBcd
           { subRegion = 'c';}
        else if T <= tBgh
           { subRegion = 'l';}
        else if T <= tBef
            {subRegion = 'h';}
        else if T <= tBij
           { subRegion = 'i';}
        else if T <= tBjk
           { subRegion = 'j';}
        else
            {subRegion = 'k';}
    }
    else if p > 22.5
    {
        let tBcd:f64 = T_atRegionBoundary(p, "3cd".to_string());
        let tBef:f64 = T_atRegionBoundary(p, "3ef".to_string());
        let tBgh:f64 = T_atRegionBoundary(p, "3gh".to_string());
        let tBij:f64 = T_atRegionBoundary(p, "3ij".to_string());
        let tBjk:f64 = T_atRegionBoundary(p, "3jk".to_string());
        let tBmn:f64 = T_atRegionBoundary(p, "3mn".to_string());
        let tBop:f64 = T_atRegionBoundary(p, "3op".to_string());

        if T <= tBcd
           { subRegion = 'c';}
        else if T <= tBgh
          {  subRegion = 'l';}
        else if T <= tBmn
          {  subRegion = 'm';}
        else if T <= tBef
           { subRegion = 'n';}
        else if T <= tBop
           { subRegion = 'o';}
        else if T <= tBij
           { subRegion = 'p';}
        else if T <= tBjk
            {subRegion = 'j';}
        else
            {subRegion = 'k';}
    }
    else
    {
        let p_sat:f64 = p_saturation(643.15); //需要补充4区，饱和方程
        if p > p_sat
        {
            let tBcd:f64 = T_atRegionBoundary(p, "3cd".to_string());
            let tBqu:f64 = T_atRegionBoundary(p, "3qu".to_string());
            let tBrx:f64 = T_atRegionBoundary(p, "3rx".to_string());
            let tBjk:f64 = T_atRegionBoundary(p, "3jk".to_string());

            if T <= tBcd
             {   subRegion = 'c';}
            else if T <= tBqu
            {
                subRegion = 'q';
                //这里的判断算法有点问题，下面这个 subRegion会被重新设定u
                // 简便处理，这个判断后，立即返回
                return subRegion;
            }

            if (T > tBrx && T <= tBjk)
            {
                subRegion = 'r';
                // 这里的判断算法有点问题，下面这个 subRegion会被重新设定u
                // 简便处理，这个判断后，立即返回
                return subRegion;
            }

            if T > tBjk
               { subRegion = 'k';}

            // tBqu < T <= tBrx
            // small regions right around critical point
            // 3u, 3x, 3y, 3z, 3v and 3w

            let tBuv:f64 = T_atRegionBoundary(p, "3uv".to_string());
            let tBwx:f64 = T_atRegionBoundary(p, "3wx".to_string());
            let tBef:f64 = T_atRegionBoundary(p, "3ef".to_string());

            if p > 22.11
            {
                if (T <= tBuv)
                    {subRegion = 'u';}
                else if (T < tBef)
                  {  subRegion = 'v';}
                else if (T < tBwx)
                   { subRegion = 'w';}
                else
                   { subRegion = 'x';}
            }
            else if p > 22.064
            {
                let tBwx:f64 = T_atRegionBoundary(p, "3wx".to_string());
                let tBuv:f64 = T_atRegionBoundary(p, "3uv".to_string());
                if T <= tBuv
                   { subRegion = 'u';}
                else if T <= tBef
                  {  subRegion = 'y';}
                else if (T <= tBwx)
                  {  subRegion = 'z';}
                else
                   { subRegion = 'x';}
            }
            else
            {
                let T_sat:f64 =T_saturation(p);

                if (T <= T_sat)
                {
                    if p > 21.93161551
                    {
                        if T < tBuv
                           { subRegion = 'u';}
                        else
                           { subRegion = 'y';}
                    }
                    else
                      {  subRegion = 'u';}
                }
                else
                {
                    if (p > 21.90096265)
                    {
                        if T < tBwx
                          {  subRegion = 'z';}
                        else
                          {  subRegion = 'x';}
                    }
                    else
                      {  subRegion = 'x';}
                }
            }
        }
        else if (p > 20.5)
        {
            let tBcd = T_atRegionBoundary(p, "3cd".to_string());
            let tBjk = T_atRegionBoundary(p, "3jk".to_string());
            let T_sat:f64 = T_saturation(p);

            if T <= tBcd
              {  subRegion = 'c';}
            else if T <= T_sat
              {  subRegion = 's';}
            else if T <= tBjk
             {   subRegion = 'r';}
            else
               { subRegion = 'k';}
        }
        else if (p > 19.00881189173929)
        {
            let tBcd = T_atRegionBoundary(p, "3cd".to_string());
            let T_sat:f64 = T_saturation(p);

            if T <= tBcd
                {subRegion = 'c';}
            else if T <= T_sat
             {   subRegion = 's';}
            else
            {    subRegion = 't';}
        }
        else
        {
            let T_sat:f64 = T_saturation(p);
            if T <= T_sat
               { subRegion = 'c';}
            else
               { subRegion = 't';}
        }
    }
    return subRegion;
}

pub  fn pT2v_reg3(p:f64,T:f64)->f64
{
    let sub_region:char = sub_region3_pT(p, T);
    let v:f64=pT2v_3subreg(p, T, sub_region);
    return v;
}
