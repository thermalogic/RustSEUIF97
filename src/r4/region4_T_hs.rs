//! Region 4 - Backward Equation Tsat(h,s)
//! Page25,Page30 Eq 9: <http://www.iapws.org/relguide/Supp-phs3-2014.pdf>
//! *  5.3 Backward Equation Tsat(h,s)
//!   * s> 5.210 887 825
//!   * Temperature range is T(273.15,623.15)

use crate::algo::*;
use crate::common::boundaries::*;
use crate::common::constant::*;
use crate::common::property_id::*;
use crate::r1::*;
use crate::r2::region2_pT::*;
use crate::r3::*;
use crate::r4::*;

const s4V_623: f64 = 5.210887825; // T=623.16 ，Sature steam   s
const h4V_623: f64 = 2.5635920043e+03; // Page 25  Sature steam   h
const s4L_273: f64 = -1.545495919e-04; // T=273.15 ，Sature Liquid s
const s4V_273: f64 = 9.155759395; // T=273.15 ，Sature steam s

///  Page30  5.3 Backward Functions Tsat(h,s),
pub fn hs2T_reg43(h: f64, s: f64) -> f64 {
    // Table 28
    const IJn: [(i32, i32, f64); 36] = [
        (0, 0, 0.179882673606601),
        (0, 3, -0.267507455199603),
        (0, 12, 0.116276722612600e1),
        (1, 0, 0.147545428713616),
        (1, 1, -0.512871635973248),
        (1, 2, 0.421333567697984),
        (1, 5, 0.563749522189870),
        (2, 0, 0.429274443819153),
        (2, 5, -0.335704552142140e1),
        (2, 8, 0.108890916499278e2),
        (3, 0, -0.248483390456012),
        (3, 2, 0.304153221906390),
        (3, 3, -0.494819763939905),
        (3, 4, 0.107551674933261e1),
        (4, 0, 0.73388415457688e-1),
        (4, 1, 0.140170545411085e-1),
        (5, 1, -0.106110975998808),
        (5, 2, 0.168324361811875e-1),
        (5, 4, 0.125028363714877e1),
        (5, 16, 0.101316840309509e4),
        (6, 6, -0.151791158000712e1),
        (6, 8, 0.524277865990866e2),
        (6, 22, 0.230495545563912e5),
        (8, 1, 0.249459806365456e-1),
        (10, 20, 0.210796467412137e+07),
        (10, 36, 0.366836848613065e9),
        (12, 24, -0.1448105565163e9),
        (14, 1, -0.179276373003590e-2),
        (14, 28, 0.489955602100459e10),
        (16, 12, 0.471262212070518e3),
        (16, 32, -0.829294390198652e11),
        (18, 14, -0.17154662263191e4),
        (18, 22, 0.355777682973575e7),
        (18, 36, 0.586062760258436e12),
        (20, 24, -0.129887635078195e8),
        (28, 36, 0.317247449371057e11),
    ];

    let nu: f64 = h / 2800.0 - 0.119;
    let sigma: f64 = s / 9.2 - 1.07;
    let suma: f64 = poly_powi(nu, sigma, &IJn);
    550.0 * suma
}

/// 1D 残差: r(T) = (h-hl)*(sv-sl) - (s-sl)*(hv-hl)
/// 在正确的饱和温度下, 由 h 和 s 分别算出的干度 x 应当一致
fn residual_T(T: f64, h: f64, s: f64) -> f64 {
    let p = p_saturation(T);
    let hl = pT2h_reg1(p, T);
    let hv = pT2h_reg2(p, T);
    let sl = pT2s_reg1(p, T);
    let sv = pT2s_reg2(p, T);
    (h - hl)/(hv - hl) - (s - sl)/ (sv - sl)
}

fn bisect_T(h: f64, s: f64, mut a: f64, mut b: f64, tol: f64) -> f64 {
    let mut ra = residual_T(a, h, s);
    let mut rb = residual_T(b, h, s);

    // 端点符号相同: 向外扩展搜索区间
    let mut n_expand = 0;
    while ra * rb > 0.0 && n_expand < 200 {
        let width = b - a;
        a = (a - width).max(T_MIN4);
        b = (b + width).min(TC_WATER - 1.0);
        ra = residual_T(a, h, s);
        rb = residual_T(b, h, s);
        n_expand += 1;
        if a <= T_MIN4 + 0.5 && b >= TC_WATER - 1.5 {
            break;
        }
    }

    if ra * rb > 0.0 {
        // 全区间扫描找符号变化
        let n_steps = 60;
        let step = (TC_WATER - 1.0 - T_MIN4) / n_steps as f64;
        let mut prev_r = residual_T(T_MIN4, h, s);
        let mut found_lo = T_MIN4;
        let mut found_hi = T_MIN4 + step;
        for i in 1..=n_steps {
            let Ti = T_MIN4 + step * i as f64;
            let ri = residual_T(Ti, h, s);
            if prev_r * ri < 0.0 {
                found_lo = Ti - step;
                found_hi = Ti;
                break;
            }
            prev_r = ri;
        }
        a = found_lo;
        b = found_hi;
        ra = residual_T(a, h, s);
    }

    // 二分
    for _ in 0..200 {
        let mid = 0.5 * (a + b);
        let rm = residual_T(mid, h, s);
        if rm.abs() < tol || (b - a) < 1.0e-10 {
            return mid;
        }
        if ra * rm < 0.0 {
            b = mid;
            rb = rm;
        } else {
            a = mid;
            ra = rm;
        }
    }
    0.5 * (a + b)
}

pub fn hs2T_reg4(h: f64, s: f64) -> f64 {
    let mut T = hs2T_reg43(h, s);
    if s > s4V_623 && s < s4V_273 {
       return T;
    };
    if T < T_MIN4 || T > TC_WATER - 1.0 {
        T = 300.0;
    }
    bisect_T(h, s, T - 0.5, T + 0.5, 1.0e-8)   
}
