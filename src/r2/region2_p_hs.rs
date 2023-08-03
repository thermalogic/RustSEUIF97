//! Region 2 - Backward Equation(h,s)->p
//！
//！ <http://www.iapws.org/relguide/Supp-PHS12-2014.pdf>
//!
//！ Backward Equations p(h,s) for Region 2:  hs2p_reg2(h,s)

use crate::algo::*;
use crate::common::constant::*;
use crate::r2::region2_T_ph::*;
use crate::r2::region2_pT::*;

/// for iter (h,s)->p
pub fn ph2s_reg2(p: f64, h: f64) -> f64 {
    let T: f64 = ph2T_reg2(p, h);
    return pT2s_reg2(p, T);
}

pub fn hs2p_reg2a(h: f64, s: f64) -> f64 {
    // Initialize coefficients and exponents (h,s)->P for region 2a
    const IJn: [(i32, i32, f64); 29] = [
        (0, 1, -0.182575361923032E-01),
        (0, 3, -0.125229548799536),
        (0, 6, 0.592290437320145),
        (0, 16, 0.604769706185122E+01),
        (0, 20, 0.238624965444474E+03),
        //
        (0, 22, -0.298639090222922E+03),
        (1, 0, 0.512250813040750E-01),
        (1, 1, -0.437266515606486),
        (1, 2, 0.413336902999504),
        (1, 3, -0.516468254574773E+01),
        //
        (1, 5, -0.557014838445711E+01),
        (1, 6, 0.128555037824478E+02),
        (1, 10, 0.114144108953290E+02),
        (1, 16, -0.119504225652714E+03),
        (1, 20, -0.284777985961560E+04),
        //
        (1, 22, 0.431757846408006E+04),
        (2, 3, 0.112894040802650E+01),
        (2, 16, 0.197409186206319E+04),
        (2, 20, 0.151612444706087E+04),
        (3, 0, 0.141324451421235E-01),
        //
        (3, 2, 0.585501282219601),
        (3, 3, -0.297258075863012E+01),
        (3, 6, 0.594567314847319E+01),
        (3, 16, -0.623656565798905E+04),
        (4, 16, 0.965986235133332E+04),
        //
        (5, 3, 0.681500934948134E+01),
        (5, 16, -0.633207286824489E+04),
        (6, 3, -0.558919224465760E+01),
        (7, 1, 0.400645798472063E-01),
    ];

    let eta: f64 = h / 4200.0 - 0.5;
    let sigma: f64 = s / 12.0 - 1.2;
    let pi: f64 = poly_powi(eta, sigma, &IJn);
    let pi2: f64 = pi * pi;
    4.0 * pi2 * pi2
}

pub fn hs2p_reg2b(h: f64, s: f64) -> f64 {
    // Table 7 Initialize coefficients and exponents (H,S)->P for region 2b
    const IJn: [(i32, i32, f64); 33] = [
        (0, 0, 0.801496989929495E-01),
        (0, 1, -0.543862807146111),
        (0, 2, 0.337455597421283),
        (0, 4, 0.890555451157450E+01),
        (0, 8, 0.313840736431485E+03),
        //
        (1, 0, 0.797367065977789),
        (1, 1, -0.121616973556240E+01),
        (1, 2, 0.872803386937477E+01),
        (1, 3, -0.169769781757602E+02),
        (1, 5, -0.186552827328416E+03),
        //
        (1, 12, 0.951159274344237E+05),
        (2, 1, -0.189168510120494E+02),
        (2, 6, -0.433407037194840E+04),
        (2, 18, 0.543212633012715E+09),
        (3, 0, 0.144793408386013),
        //
        (3, 1, 0.128024559637516E+03),
        (3, 7, -0.672309534071268E+05),
        (3, 12, 0.336972380095287E+08),
        (4, 1, -0.586634196762720E+03),
        (4, 16, -0.221403224769889E+11),
        //
        (5, 1, 0.171606668708389E+04),
        (5, 12, -0.570817595806302E+09),
        (6, 1, -0.312109693178482E+04),
        (6, 8, -0.207841384633010E+07),
        (6, 18, 0.305605946157786E+13),
        //
        (7, 1, 0.322157004314333E+04),
        (7, 16, 0.326810259797295E+12),
        (8, 1, -0.144104158934487E+04),
        (8, 3, 0.410694867802691E+03),
        (8, 14, 0.109077066873024E+12),
        //
        (8, 18, -0.247964654258893E+14),
        (12, 10, 0.188801906865134E+10),
        (14, 16, -0.123651009018773E+15),
    ];
    let eta: f64 = h / 4100.0 - 0.6;
    let sigma: f64 = s / 7.9 - 1.01;
    let pi: f64 = poly_powi(eta, sigma, &IJn);
    let pi2: f64 = pi * pi;
    100.0 * pi2 * pi2
}

pub fn hs2p_reg2c(h: f64, s: f64) -> f64 {
    //   Table 8 Initialize coefficients and exponents (H,S)->P for region 2c
    const IJn: [(i32, i32, f64); 31] = [
        (0, 0, 0.112225607199012E+00),
        (0, 1, -0.339005953606712E+01),
        (0, 2, -0.320503911730094E+02),
        (0, 3, -0.197597305104900E+03),
        (0, 4, -0.407693861553446E+03),
        //
        (0, 8, 0.132943775222331E+05),
        (1, 0, 0.170846839774007E+01),
        (1, 2, 0.373694198142245E+02),
        (1, 5, 0.358144365815434E+04),
        (1, 8, 0.423014446424664E+06),
        //
        (1, 14, -0.751071025760063E+09),
        (2, 2, 0.523446127607898E+02),
        (2, 3, -0.228351290812417E+03),
        (2, 7, -0.960652417056937E+06),
        (2, 10, -0.807059292526074E+08),
        //
        (2, 18, 0.162698017225669E+13),
        (3, 0, 0.772465073604171),
        (3, 5, 0.463929973837746E+05),
        (3, 8, -0.137317885134128E+08),
        (3, 16, 0.170470392630512E+13),
        //
        (3, 18, -0.251104628187308E+14),
        (4, 18, 0.317748830835520E+14),
        (5, 1, 0.538685623675312E+02),
        (5, 4, -0.553089094625169E+05),
        (5, 6, -0.102861522421405E+07),
        //
        (5, 14, 0.204249418756234E+13),
        (6, 8, 0.273918446626977E+09),
        (6, 18, -0.263963146312685E+16),
        (10, 7, -0.107890854108088E+10),
        (12, 7, -0.296492620980124E+11),
        //
        (16, 10, -0.111754907323424E+16),
    ];

    let eta: f64 = h / 3500.0 - 0.7;
    let sigma: f64 = s / 5.9 - 1.1;
    let pi: f64 = poly_powi(eta, sigma, &IJn);
    let pi2: f64 = pi * pi;
    100.0 * pi2 * pi2
}

//------------------------------------------------------------
// Region 2 (h,s)
//------------------------------------------------------------

/// Define the boundary between Region 2a and 2b, h=f(s)
///    Supp-PHS12-2014.pdf, Eq 2
pub fn s2h_reg2_ab(s: f64) -> f64 {
    let n: [f64; 4] = [-0.349898083432139E+04, 0.257560716905876E+04, -0.421073558227969E+03, 0.276349063799944E+02];
    let sigma: f64 = s / 1.0;
    let eta: f64 = n[0] + sigma * (n[1] + sigma * (n[2] + n[3] * sigma));
    eta * 1.0
}

pub fn hs2p_reg2(h: f64, s: f64) -> f64 {
    let h2ab: f64 = s2h_reg2_ab(s);
    let mut p: f64 = 0.0;
    if h > h2ab {
        if s >= 5.85 {
            p = hs2p_reg2b(h, s);
        } else {
            p = hs2p_reg2c(h, s);
        };
    } else {
        p = hs2p_reg2a(h, s);
    };
    p
}
