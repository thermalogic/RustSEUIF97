/*--------------------------------------------------------------------------
 Backward Equation for Region 2:
   IAPWS-IF97-Rev : (P,s)->T
   Page 25:
       6.3.2 The Backward Equations T( p, s ) for Subregions 2a, 2b, and 2c.
          ps2T_reg2(p,s)

    ps2T_reg2a, I 不是整数，不能使用SAC方法！
 Author： Maohua Cheng
-------------------------------------------------------------------------*/
use crate::algo::fast_ipower::sac_pow;
use crate::algo::root::rtsec2;
use crate::algo::root::ESP;
use crate::algo::root::IF97_EQ;
use crate::algo::root::I_MAX;
use crate::common::constant::*;
use crate::r2::region2_pT::*;

pub struct fIJnData
// region 2: ps2T_reg2a,  I is float
{
    pub I: f64,
    pub J: i32,
    pub n: f64,
}

//-----  the backward equation T( p,s ) ---------------------
pub fn ps2T_reg2a(p: f64, s: f64) -> f64 {
    // Table 25. Page 26  Numerical values of the coefficients and exponents of
    // the backward equation T( p,s ) for
    // subregion 2a, Eq. (25)
    const IJn: [fIJnData; 46] = [
        fIJnData {
            I: -1.5,
            J: -24,
            n: -0.39235983861984E+06,
        },
        fIJnData {
            I: -1.5,
            J: -23,
            n: 0.51526573827270E+06,
        },
        fIJnData {
            I: -1.5,
            J: -19,
            n: 0.40482443161048E+05,
        },
        fIJnData {
            I: -1.5,
            J: -13,
            n: -0.32193790923902E+03,
        },
        fIJnData {
            I: -1.5,
            J: -11,
            n: 0.96961424218694E+02,
        },
        //6
        fIJnData {
            I: -1.5,
            J: -10,
            n: -0.22867846371773E+02,
        },
        fIJnData {
            I: -1.25,
            J: -19,
            n: -0.44942914124357E+06,
        },
        fIJnData {
            I: -1.25,
            J: -15,
            n: -0.50118336020166E+04,
        },
        fIJnData {
            I: -1.25,
            J: -6,
            n: 0.35684463560015E+00,
        },
        fIJnData {
            I: -1.0,
            J: -26,
            n: 0.44235335848190E+05,
        },
        // 11
        fIJnData {
            I: -1.0,
            J: -21,
            n: -0.13673388811708E+05,
        },
        fIJnData {
            I: -1.0,
            J: -17,
            n: 0.42163260207864E+06,
        },
        fIJnData {
            I: -1.0,
            J: -16,
            n: 0.22516925837475E+05,
        },
        fIJnData {
            I: -1.0,
            J: -9,
            n: 0.47442144865646E+03,
        },
        fIJnData {
            I: -1.0,
            J: -8,
            n: -0.14931130797647E+03,
        },
        // 16
        fIJnData {
            I: -0.75,
            J: -15,
            n: -0.19781126320452E+06,
        },
        fIJnData {
            I: -0.75,
            J: -14,
            n: -0.23554399470760E+05,
        },
        fIJnData {
            I: -0.5,
            J: -26,
            n: -0.19070616302076E+05,
        },
        fIJnData {
            I: -0.5,
            J: -13,
            n: 0.55375669883164E+05,
        },
        fIJnData {
            I: -0.5,
            J: -9,
            n: 0.38293691437363E+04,
        },
        //21
        fIJnData {
            I: -0.5,
            J: -7,
            n: -0.60391860580567E+03,
        },
        fIJnData {
            I: -0.25,
            J: -27,
            n: 0.19363102620331E+04,
        },
        fIJnData {
            I: -0.25,
            J: -25,
            n: 0.42660643698610E+04,
        },
        fIJnData {
            I: -0.25,
            J: -11,
            n: -0.59780638872718E+04,
        },
        fIJnData {
            I: -0.25,
            J: -6,
            n: -0.70401463926862E+03,
        },
        //26
        fIJnData {
            I: 0.25,
            J: 1,
            n: 0.33836784107553E+03,
        },
        fIJnData {
            I: 0.25,
            J: 4,
            n: 0.20862786635187E+02,
        },
        fIJnData {
            I: 0.25,
            J: 8,
            n: 0.33834172656196E-01,
        },
        fIJnData {
            I: 0.25,
            J: 11,
            n: -0.43124428414893E-04,
        },
        fIJnData {
            I: 0.5,
            J: 0,
            n: 0.16653791356412E+03,
        },
        // 31
        fIJnData {
            I: 0.5,
            J: 1,
            n: -0.13986292055898E+03,
        },
        fIJnData {
            I: 0.5,
            J: 5,
            n: -0.78849547999872E+00,
        },
        fIJnData {
            I: 0.5,
            J: 6,
            n: 0.72132411753872E-01,
        },
        fIJnData {
            I: 0.5,
            J: 10,
            n: -0.59754839398283E-02,
        },
        fIJnData {
            I: 0.5,
            J: 14,
            n: -0.12141358953904E-04,
        },
        // 36
        fIJnData {
            I: 0.5,
            J: 16,
            n: 0.23227096733871E-06,
        },
        fIJnData {
            I: 0.75,
            J: 0,
            n: -0.10538463566194E+02,
        },
        fIJnData {
            I: 0.75,
            J: 4,
            n: 0.20718925496502E+01,
        },
        fIJnData {
            I: 0.75,
            J: 9,
            n: -0.72193155260427E-01,
        },
        fIJnData {
            I: 0.75,
            J: 17,
            n: 0.20749887081120E-06,
        },
        // 41
        fIJnData {
            I: 1.0,
            J: 7,
            n: -0.18340657911379E-01,
        },
        fIJnData {
            I: 1.0,
            J: 18,
            n: 0.29036272348696E-06,
        },
        fIJnData {
            I: 1.25,
            J: 3,
            n: 0.21037527893619E+00,
        },
        fIJnData {
            I: 1.25,
            J: 15,
            n: 0.2568123972999E-03,
        },
        fIJnData {
            I: 1.5,
            J: 5,
            n: -0.1279900293381E-01,
        },
        // 46
        fIJnData {
            I: 1.5,
            J: 18,
            n: -0.82198102652018E-05,
        },
    ];

    let pi: f64 = p / 1.0;
    let sigma: f64 = s / 2.0 - 2.0;

    let mut theta: f64 = 0.0;
    println!("pi {} sigma {}", pi, sigma);
    for k in 0..46 {
        theta += IJn[k].n * pi.powf(IJn[k].I) * sac_pow(sigma, IJn[k].J); // IJn[k].I) is float
                                                                          // theta += IJn[k].n * pi.powf(IJn[k].I) * sigma.powi(IJn[k].J); // IJn[k].I) is float
    }
    return 1.0 * theta;
}

pub fn ps2T_reg2b(p: f64, s: f64) -> f64 {
    // Table 26. Page 27
    //   Numerical values of the coefficients and exponents of
    //     the backward equation T( p,s ) for subregion 2b, Eq. (26)
    const IJn: [IJnData; 44] = [
        IJnData {
            I: -6,
            J: 0,
            n: 0.31687665083497e6,
        },
        IJnData {
            I: -6,
            J: 11,
            n: 0.20864175881858e2,
        },
        IJnData {
            I: -5,
            J: 0,
            n: -0.39859399803599e6,
        },
        IJnData {
            I: -5,
            J: 11,
            n: -0.21816058518877e2,
        },
        IJnData {
            I: -4,
            J: 0,
            n: 0.22369785194242e6,
        },
        IJnData {
            I: -4,
            J: 1,
            n: -0.27841703445817e4,
        },
        IJnData {
            I: -4,
            J: 11,
            n: 0.99207436071480e1,
        },
        IJnData {
            I: -3,
            J: 0,
            n: -0.75197512299157e5,
        },
        IJnData {
            I: -3,
            J: 1,
            n: 0.29708605951158e4,
        },
        IJnData {
            I: -3,
            J: 11,
            n: -0.34406878548526e1,
        },
        IJnData {
            I: -3,
            J: 12,
            n: 0.38815564249115,
        },
        IJnData {
            I: -2,
            J: 0,
            n: 0.17511295085750e5,
        },
        IJnData {
            I: -2,
            J: 1,
            n: -0.14237112854449e4,
        },
        IJnData {
            I: -2,
            J: 6,
            n: 0.10943803364167e1,
        },
        IJnData {
            I: -2,
            J: 10,
            n: 0.89971619308495,
        },
        IJnData {
            I: -1,
            J: 0,
            n: -0.33759740098958e4,
        },
        IJnData {
            I: -1,
            J: 1,
            n: 0.47162885818355e3,
        },
        IJnData {
            I: -1,
            J: 5,
            n: -0.19188241993679e1,
        },
        IJnData {
            I: -1,
            J: 8,
            n: 0.41078580492196,
        },
        IJnData {
            I: -1,
            J: 9,
            n: -0.33465378172097,
        },
        IJnData {
            I: 0,
            J: 0,
            n: 0.13870034777505e4,
        },
        IJnData {
            I: 0,
            J: 1,
            n: -0.40663326195838e3,
        },
        IJnData {
            I: 0,
            J: 2,
            n: 0.41727347159610e2,
        },
        IJnData {
            I: 0,
            J: 4,
            n: 0.21932549434532e1,
        },
        IJnData {
            I: 0,
            J: 5,
            n: -0.10320050009077e1,
        },
        IJnData {
            I: 0,
            J: 6,
            n: 0.35882943516703,
        },
        IJnData {
            I: 0,
            J: 9,
            n: 0.52511453726066e-2,
        },
        IJnData {
            I: 1,
            J: 0,
            n: 0.12838916450705e2,
        },
        IJnData {
            I: 1,
            J: 1,
            n: -0.28642437219381e1,
        },
        IJnData {
            I: 1,
            J: 2,
            n: 0.56912683664855,
        },
        IJnData {
            I: 1,
            J: 3,
            n: -0.99962954584931e-1,
        },
        IJnData {
            I: 1,
            J: 7,
            n: -0.32632037778459e-2,
        },
        IJnData {
            I: 1,
            J: 8,
            n: 0.23320922576723e-3,
        },
        IJnData {
            I: 2,
            J: 0,
            n: -0.15334809857450,
        },
        IJnData {
            I: 2,
            J: 1,
            n: 0.29072288239902e-1,
        },
        IJnData {
            I: 2,
            J: 5,
            n: 0.37534702741167e-3,
        },
        IJnData {
            I: 3,
            J: 0,
            n: 0.17296691702411e-2,
        },
        IJnData {
            I: 3,
            J: 1,
            n: -0.38556050844504e-3,
        },
        IJnData {
            I: 3,
            J: 3,
            n: -0.35017712292608e-4,
        },
        IJnData {
            I: 4,
            J: 0,
            n: -0.14566393631492e-4,
        },
        IJnData {
            I: 4,
            J: 1,
            n: 0.56420857267269e-5,
        },
        IJnData {
            I: 5,
            J: 0,
            n: 0.41286150074605e-7,
        },
        IJnData {
            I: 5,
            J: 1,
            n: -0.20684671118824e-7,
        },
        IJnData {
            I: 5,
            J: 2,
            n: 0.16409393674725e-8,
        },
    ];

    let pi: f64 = p / 1.0;
    let sigma: f64 = 10.0 - s / 0.7853;

    let mut theta: f64 = 0.0;
    for k in 0..44 {
        theta += IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(sigma, IJn[k].J);
    }
    return 1.0 * theta;
}

pub fn ps2T_reg2c(p: f64, s: f64) -> f64 {
    // Table 27. Page 28
    // Numerical values of the coefficient s and exponents of
    //  the backward equation T( p,s ) for subregion 2c, Eq. (27)
    const IJn: [IJnData; 30] = [
        IJnData {
            I: -2,
            J: 0,
            n: 0.90968501005365e3,
        },
        IJnData {
            I: -2,
            J: 1,
            n: 0.24045667088420e4,
        },
        IJnData {
            I: -1,
            J: 0,
            n: -0.59162326387130e3,
        },
        IJnData {
            I: 0,
            J: 0,
            n: 0.54145404128074e3,
        },
        IJnData {
            I: 0,
            J: 1,
            n: -0.27098308411192e3,
        },
        IJnData {
            I: 0,
            J: 2,
            n: 0.97976525097926e3,
        },
        IJnData {
            I: 0,
            J: 3,
            n: -0.46966772959435e3,
        },
        IJnData {
            I: 1,
            J: 0,
            n: 0.14399274604723e2,
        },
        IJnData {
            I: 1,
            J: 1,
            n: -0.19104204230429e2,
        },
        IJnData {
            I: 1,
            J: 3,
            n: 0.53299167111971e1,
        },
        IJnData {
            I: 1,
            J: 4,
            n: -0.21252975375934e2,
        },
        IJnData {
            I: 2,
            J: 0,
            n: -0.31147334413760,
        },
        IJnData {
            I: 2,
            J: 1,
            n: 0.60334840894623,
        },
        IJnData {
            I: 2,
            J: 2,
            n: -0.42764839702509e-1,
        },
        IJnData {
            I: 3,
            J: 0,
            n: 0.58185597255259e-2,
        },
        IJnData {
            I: 3,
            J: 1,
            n: -0.14597008284753e-1,
        },
        IJnData {
            I: 3,
            J: 5,
            n: 0.56631175631027e-2,
        },
        IJnData {
            I: 4,
            J: 0,
            n: -0.76155864584577e-4,
        },
        IJnData {
            I: 4,
            J: 1,
            n: 0.22440342919332e-3,
        },
        IJnData {
            I: 4,
            J: 4,
            n: -0.12561095013413e-4,
        },
        IJnData {
            I: 5,
            J: 0,
            n: 0.63323132660934e-6,
        },
        IJnData {
            I: 5,
            J: 1,
            n: -0.20541989675375e-5,
        },
        IJnData {
            I: 5,
            J: 2,
            n: 0.36405370390082e-7,
        },
        IJnData {
            I: 6,
            J: 0,
            n: -0.29759897789215e-8,
        },
        IJnData {
            I: 6,
            J: 1,
            n: 0.10136618529763e-7,
        },
        IJnData {
            I: 7,
            J: 0,
            n: 0.59925719692351e-11,
        },
        IJnData {
            I: 7,
            J: 1,
            n: -0.20677870105164e-10,
        },
        IJnData {
            I: 7,
            J: 3,
            n: -0.20874278181886e-10,
        },
        IJnData {
            I: 7,
            J: 4,
            n: 0.10162166825089e-9,
        },
        IJnData {
            I: 7,
            J: 5,
            n: -0.16429828281347e-9,
        },
    ];

    let pi: f64 = p / 1.0;
    let sigma: f64 = 2.0 - s / 2.9251;

    let mut theta: f64 = 0.0;
    for k in 0..30 {
        theta += IJn[k].n * sac_pow(pi, IJn[k].I) * sac_pow(sigma, IJn[k].J);
    }
    return 1.0 * theta;
}

pub fn ps2T_reg2(p: f64, s: f64) -> f64 {
    let mut T: f64 = 0.0;
    if p > 4.0 {
        if s < 5.85 {
            T = ps2T_reg2c(p, s);
        } else {
            T = ps2T_reg2b(p, s);
        }
    } else {
        T = ps2T_reg2a(p, s);
    }
    //return T;

    let mut T1: f64 = T;
    let f1: f64 = s - pT2s_reg2(p, T1);

    if f1.abs() > ESP {
        let mut T2: f64 = 0.0;
        if f1 > 0.0
        // pT2sreg1(p,T1)< s ,the T1< expt T，so， T2=1.05*T1 T（T1,T2)
        {
            T2 = (1.0 + f1 / s) * T1;
        } else {
            T2 = (1.0 - f1 / s) * T1;
        }

        let f2: f64 = s - pT2s_reg2(p, T2);

        T = rtsec2(pT2s_reg2, p, s, T1, T2, f1, f2, ESP, I_MAX);
    };
    return T;
}
