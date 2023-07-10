//！ Backward Equation for Region 2:
//！   IAPWS-IF97-Rev : (P,H)->T
//！       6.3.1 The Backward Equations T( p, h ) for Subregions 2a, 2b, and 2c.
//!               ph2T_reg2(p,h)

use crate::algo::fast_ipower::sac_pow;
use crate::algo::root::rtsec2;
use crate::algo::root::ESP;
use crate::algo::root::IF97_EQ;
use crate::algo::root::I_MAX;
use crate::common::constant::*;
use crate::r2::region2_pT::*;

pub fn ph2T_reg2a(p: f64, h: f64) -> f64 {
// Table 20 the coefficients and exponents of  the backward equation T ( p,h ) for subregion 2a, Eq.(22)
    const IJn: [IJnData; 34] = [
        IJnData {
            I: 0,
            J: 0,
            n: 0.10898952318288E+04,
        },
        IJnData {
            I: 0,
            J: 1,
            n: 0.84951654495535E+03,
        },
        IJnData {
            I: 0,
            J: 2,
            n: -0.10781748091826E+03,
        },
        IJnData {
            I: 0,
            J: 3,
            n: 0.33153654801263E+02,
        },
        IJnData {
            I: 0,
            J: 7,
            n: -0.74232016790248E+01,
        },
        IJnData {
            I: 0,
            J: 20,
            n: 0.11765048724356E+02,
        },
        IJnData {
            I: 1,
            J: 0,
            n: 0.18445749355790E+01,
        },
        IJnData {
            I: 1,
            J: 1,
            n: -0.41792700549624E+01,
        },
        IJnData {
            I: 1,
            J: 2,
            n: 0.62478196935812E+01,
        },
        IJnData {
            I: 1,
            J: 3,
            n: -0.17344563108114E+02,
        },
        IJnData {
            I: 1,
            J: 7,
            n: -0.20058176862096E+03,
        },
        IJnData {
            I: 1,
            J: 9,
            n: 0.27196065473796E+03,
        },
        IJnData {
            I: 1,
            J: 11,
            n: -0.45511318285818E+03,
        },
        IJnData {
            I: 1,
            J: 18,
            n: 0.30919688604755E+04,
        },
        IJnData {
            I: 1,
            J: 44,
            n: 0.25226640357872E+06,
        },
        IJnData {
            I: 2,
            J: 0,
            n: -0.61707422868339E-02,
        },
        IJnData {
            I: 2,
            J: 2,
            n: -0.31078046629583E+00,
        },
        IJnData {
            I: 2,
            J: 7,
            n: 0.11670873077107E+02,
        },
        IJnData {
            I: 2,
            J: 36,
            n: 0.12812798404046E+09,
        },
        IJnData {
            I: 2,
            J: 38,
            n: -0.98554909623276E+09,
        },
        //
        IJnData {
            I: 2,
            J: 40,
            n: 0.28224546973002E+10,
        },
        IJnData {
            I: 2,
            J: 42,
            n: -0.35948971410703E+10,
        },
        IJnData {
            I: 2,
            J: 44,
            n: 0.17227349913197E+10,
        },
        IJnData {
            I: 3,
            J: 24,
            n: -0.13551334240775E+05,
        },
        IJnData {
            I: 3,
            J: 44,
            n: 0.12848734664650E+08,
        },
        IJnData {
            I: 4,
            J: 12,
            n: 0.13865724283226E+01,
        },
        IJnData {
            I: 4,
            J: 32,
            n: 0.23598832556514E+06,
        },
        IJnData {
            I: 4,
            J: 44,
            n: -0.13105236545054E+08,
        },
        IJnData {
            I: 5,
            J: 32,
            n: 0.73999835474766E+04,
        },
        IJnData {
            I: 5,
            J: 36,
            n: -0.55196697030060E+06,
        },
        IJnData {
            I: 5,
            J: 42,
            n: 0.37154085996233E+07,
        },
        IJnData {
            I: 6,
            J: 34,
            n: 0.19127729239660E+05,
        },
        IJnData {
            I: 6,
            J: 44,
            n: -0.41535164835634E+06,
        },
        IJnData {
            I: 7,
            J: 28,
            n: -0.62459855192507E+02,
        },
    ];

    let pi: f64 = p / 1.0;
    let eta: f64 = h / 2000.0 - 2.1;
    let mut theta: f64 = 0.0;
    for k in 0..34 {
        theta += IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(eta, IJn[k].J);
    }
    return 1.0 * theta;
}

pub fn ph2T_reg2b(p: f64, h: f64) -> f64 {
// Table 21. the coefficients and exponents of the backward equation T ( p,h ) for  subregion 2b Eq.(23)
    const IJn: [IJnData; 38] = [
        IJnData {
            I: 0,
            J: 0,
            n: 1.4895041079516e3,
        },
        IJnData {
            I: 0,
            J: 1,
            n: 7.4307798314034e2,
        },
        IJnData {
            I: 0,
            J: 2,
            n: -9.7708318797837e1,
        },
        IJnData {
            I: 0,
            J: 12,
            n: 2.4742464705674,
        },
        IJnData {
            I: 0,
            J: 18,
            n: -6.3281320016026e-1,
        },
        IJnData {
            I: 0,
            J: 24,
            n: 1.1385952129658,
        },
        IJnData {
            I: 0,
            J: 28,
            n: -4.7811863648625e-1,
        },
        IJnData {
            I: 0,
            J: 40,
            n: 8.5208123431544e-3,
        },
        IJnData {
            I: 1,
            J: 0,
            n: 9.3747147377932e-1,
        },
        IJnData {
            I: 1,
            J: 2,
            n: 3.3593118604916,
        },
        IJnData {
            I: 1,
            J: 6,
            n: 3.3809355601454,
        },
        IJnData {
            I: 1,
            J: 12,
            n: 1.6844539671904e-1,
        },
        IJnData {
            I: 1,
            J: 18,
            n: 7.3875745236695e-1,
        },
        IJnData {
            I: 1,
            J: 24,
            n: -4.7128737436186e-1,
        },
        IJnData {
            I: 1,
            J: 28,
            n: 1.5020273139707e-1,
        },
        IJnData {
            I: 1,
            J: 40,
            n: -2.1764114219750e-3,
        },
        IJnData {
            I: 2,
            J: 2,
            n: -2.1810755324761e-2,
        },
        IJnData {
            I: 2,
            J: 8,
            n: -1.0829784403677e-1,
        },
        IJnData {
            I: 2,
            J: 18,
            n: -4.6333324635812e-2,
        },
        IJnData {
            I: 2,
            J: 40,
            n: 7.1280351959551e-5,
        },
        IJnData {
            I: 3,
            J: 1,
            n: 1.1032831789999e-4,
        },
        IJnData {
            I: 3,
            J: 2,
            n: 1.8955248387902e-4,
        },
        IJnData {
            I: 3,
            J: 12,
            n: 3.0891541160537e-3,
        },
        IJnData {
            I: 3,
            J: 24,
            n: 1.3555504554949e-3,
        },
        IJnData {
            I: 4,
            J: 2,
            n: 2.8640237477456e-7,
        },
        IJnData {
            I: 4,
            J: 12,
            n: -1.0779857357512e-5,
        },
        IJnData {
            I: 4,
            J: 18,
            n: -7.6462712454814e-5,
        },
        IJnData {
            I: 4,
            J: 24,
            n: 1.4052392818316e-5,
        },
        IJnData {
            I: 4,
            J: 28,
            n: -3.1083814331434e-5,
        },
        IJnData {
            I: 4,
            J: 40,
            n: -1.0302738212103e-6,
        },
        IJnData {
            I: 5,
            J: 18,
            n: 2.8217281635040e-7,
        },
        IJnData {
            I: 5,
            J: 24,
            n: 1.2704902271945e-6,
        },
        IJnData {
            I: 5,
            J: 40,
            n: 7.3803353468292e-8,
        },
        IJnData {
            I: 6,
            J: 28,
            n: -1.1030139238909e-8,
        },
        IJnData {
            I: 7,
            J: 2,
            n: -8.1456365207833e-14,
        },
        IJnData {
            I: 7,
            J: 28,
            n: -2.5180545682962e-11,
        },
        IJnData {
            I: 9,
            J: 1,
            n: -1.7565233969407e-18,
        },
        IJnData {
            I: 9,
            J: 40,
            n: 8.6934156344163e-15,
        },
    ];

    let pi: f64 = p / 1.0 - 2.0;
    let eta: f64 = h / 2000.0 - 2.6;
    let mut theta = 0.0;
    for k in 0..38 {
        theta += IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(eta, IJn[k].J);
    }
    return 1.0 * theta;
}

pub fn ph2T_reg2c(p: f64, h: f64) -> f64 {
// Table 22. the coefficients and exponents of the backward  equation T(p,h) for subregion 2c, Eq.(24)
    const IJn: [IJnData; 23] = [
        IJnData {
            I: -7,
            J: 0,
            n: -3236839855524.2,
        },
        IJnData {
            I: -7,
            J: 4,
            n: 7326335090218.1,
        },
        IJnData {
            I: -6,
            J: 0,
            n: 358250899454.47,
        },
        IJnData {
            I: -6,
            J: 2,
            n: -583401318515.9,
        },
        IJnData {
            I: -5,
            J: 0,
            n: -10783068217.47,
        },
        //
        IJnData {
            I: -5,
            J: 2,
            n: 20825544563.171,
        },
        IJnData {
            I: -2,
            J: 0,
            n: 610747.83564516,
        },
        IJnData {
            I: -2,
            J: 1,
            n: 859777.2253558,
        },
        IJnData {
            I: -1,
            J: 0,
            n: -25745.72360417,
        },
        IJnData {
            I: -1,
            J: 2,
            n: 31081.088422714,
        },
        IJnData {
            I: 0,
            J: 0,
            n: 1208.2315865936,
        },
        IJnData {
            I: 0,
            J: 1,
            n: 482.19755109255,
        },
        IJnData {
            I: 1,
            J: 4,
            n: 3.7966001272486,
        },
        IJnData {
            I: 1,
            J: 8,
            n: -10.842984880077,
        },
        IJnData {
            I: 2,
            J: 4,
            n: -0.04536417267666,
        },
        //
        IJnData {
            I: 6,
            J: 0,
            n: 1.4559115658698E-13,
        },
        IJnData {
            I: 6,
            J: 1,
            n: 1.126159740723E-12,
        },
        IJnData {
            I: 6,
            J: 4,
            n: -1.7804982240686E-11,
        },
        IJnData {
            I: 6,
            J: 10,
            n: 1.2324579690832E-07,
        },
        IJnData {
            I: 6,
            J: 12,
            n: -1.1606921130984E-06,
        },
        IJnData {
            I: 6,
            J: 16,
            n: 2.7846367088554E-05,
        },
        IJnData {
            I: 6,
            J: 20,
            n: -5.9270038474176E-04,
        },
        IJnData {
            I: 6,
            J: 22,
            n: 1.2918582991878E-03,
        },
    ];

    let pi: f64 = p / 1.0 + 25.0;
    let eta: f64 = h / 2000.0 - 1.8;
    let mut theta: f64 = 0.0;
    for k in 0..23 {
        theta += IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(eta, IJn[k].J);
    }
    return 1.0 * theta;
}

/// http://www.iapws.org/relguide/IF97-Rev.html, Eq 21
pub fn enthalpy_2bc(p: f64) -> f64
{
    const n: [f64; 3] = [0.12809002730136e-3, 0.26526571908428e4, 0.45257578905948e1];
    return n[1] + ((p - n[2]) / n[0]).sqrt();
}

pub fn ph2T_reg2(p: f64, h: f64) -> f64 {
    let mut T: f64 = 0.0;

    if p > 4.0 {
        if h < enthalpy_2bc(p) {
            T = ph2T_reg2c(p, h);
        } else {
            T = ph2T_reg2b(p, h);
        }
    } else {
        T = ph2T_reg2a(p, h);
    }
    //println!(" non iter to improve h - pT2h_reg2(p, T) {}",h-pT2h_reg2(p, T));
    // return T;

    let T1: f64 = T;
    let f1: f64 = h - pT2h_reg2(p, T1);
    //不满足精度要求才迭代改进
    // println!("h - pT2h_reg2(p, T1) {}",f1);
    if f1.abs() > ESP
    {
       let mut T2: f64 = 0.0;
       if f1 > 0.0 {
            T2 = (1.0 + f1 / h) * T1;
        }
        // TODO： 1.05 用 1+f1/h 是不是更快，没有测试
        else {
            T2 = (1.0 - f1 / h) * T1;
        }
        let f2: f64 = h - pT2h_reg2(p, T2);
        T = rtsec2(pT2h_reg2, p, h, T1, T2, f1, f2, ESP,I_MAX);
        // println!("h - pT2h_reg2(p, T) {}",h-pT2h_reg2(p, T));
    }; 
    return T;
}
