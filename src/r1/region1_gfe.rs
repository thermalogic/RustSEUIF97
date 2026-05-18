//！IAPWS-IF97 Basic Equation for Region 1
//! * The basic equation of the dimensionless Gibbs free energy(p,T) and its derivatives
//！   Eq.(7), P6 http://www.iapws.org/relguide/IF97-Rev.html

use crate::algo::*;
use crate::common::constant::*;
use crate::r1::*;

pub const r1pstar: f64 = 16.53; // MPa
pub const r1Tstar: f64 = 1386.0; // K

// 宏：预计算所有幂次项 x^i 和 y^j
// 使用纯乘法计算，避免 powi 函数调用
// 优化：先计算所有独立的幂次，再组装结果数组，避免重复计算
// 输入：x, y 为变量值，coeffs 为系数数组引用
// 输出：返回两个 Vec<f64>，分别包含每个系数对应的 x^i 和 y^j 值
macro_rules! precompute_powers {
    ($x:expr, $y:expr, $coeffs:expr) => {{
        let x_val = $x;
        let y_val = $y;
        
        // 预计算所有独立的 x 幂次
        // IJn 中出现的 i 值：0, 1, 2, 3, 4, 5, 8, 21, 23, 29, 30, 31, 32, -1, -2, -3, -4, -5, -6, -8, -9, -11
        let x0 = 1.0;
        let x1 = x_val;
        let x2 = x_val * x_val;
        let x3 = x2 * x_val;
        let x4 = x2 * x2;
        let x5 = x4 * x_val;
        let x8 = x4 * x4;
        let x16 = x8 * x8;
        let x21 = x16 * x8 * x4 * x2 * x1;  // x^21 = x^16 * x^4 * x^1
        let x23 = x16 * x4 * x2 * x1;        // x^23 = x^16 * x^4 * x^2 * x^1
        let x29 = x16 * x8 * x4 * x1;        // x^29 = x^16 * x^8 * x^4 * x^1
        let x30 = x16 * x8 * x4 * x2;        // x^30 = x^16 * x^8 * x^4 * x^2
        let x31 = x16 * x8 * x4 * x2 * x1;   // x^31 = x^16 * x^8 * x^4 * x^2 * x^1
        let x32 = x16 * x16;                 // x^32 = x^16 * x^16
        
        // 负数幂次
        let x_1 = 1.0 / x_val;
        let x_2 = x_1 * x_1;
        let x_3 = x_2 * x_1;
        let x_4 = x_2 * x_2;
        let x_5 = x_4 * x_1;
        let x_6 = x_3 * x_3;
        let x_8 = 1.0 / x8;
        let x_9 = x_8 * x_1;
        let x_11 = 1.0 / (x8 * x2 * x1);
        
        // 预计算所有独立的 y 幂次
        // IJn 中出现的 j 值：-2, -1, 0, 1, 2, 3, 4, 5, 6, 10, 17, -5, -3, -6, -7, -8, -9, -11, -29, -31, -38, -39, -40, -41
        let y0 = 1.0;
        let y1 = y_val;
        let y2 = y_val * y_val;
        let y3 = y2 * y1;
        let y4 = y2 * y2;
        let y5 = y4 * y1;
        let y6 = y3 * y3;
        let y8 = y4 * y4;
        let y16 = y8 * y8;
        let y10 = y8 * y2;                   // y^10 = y^8 * y^2
        let y17 = y16 * y1;                  // y^17 = y^16 * y^1
        
        // 负数幂次
        let y_1 = 1.0 / y_val;
        let y_2 = y_1 * y_1;
        let y_3 = y_2 * y_1;
        let y_5 = 1.0 / y5;
        let y_6 = 1.0 / y6;
        let y_7 = y_6 * y_1;
        let y_8 = 1.0 / y8;
        let y_9 = y_8 * y_1;
        let y_11 = 1.0 / (y8 * y2 * y1);
        let y_29 = 1.0 / (y16 * y8 * y4 * y1);   // y^-29 = 1/(y^16 * y^8 * y^4 * y^1)
        let y_31 = 1.0 / (y16 * y8 * y4 * y2 * y1); // y^-31 = 1/(y^16 * y^8 * y^4 * y^2 * y^1)
        let y_38 = 1.0 / (y16 * y16 * y4 * y2);   // y^-38 = 1/(y^16 * y^16 * y^4 * y^2)
        let y_39 = 1.0 / (y16 * y16 * y4 * y2 * y1); // y^-39 = 1/(y^16 * y^16 * y^4 * y^2 * y^1)
        let y_40 = 1.0 / (y16 * y16 * y8);        // y^-40 = 1/(y^16 * y^16 * y^8)
        let y_41 = 1.0 / (y16 * y16 * y8 * y1);   // y^-41 = 1/(y^16 * y^16 * y^8 * y^1)
        
        // 组装结果数组
        let mut x_powers = Vec::with_capacity($coeffs.len());
        let mut y_powers = Vec::with_capacity($coeffs.len());
        for &(i, j, _) in $coeffs {
            // 根据 i 值选择对应的幂次
            let x_power = match i {
                0 => x0,
                1 => x1,
                2 => x2,
                3 => x3,
                4 => x4,
                5 => x5,
                8 => x8,
                21 => x21,
                23 => x23,
                29 => x29,
                30 => x30,
                31 => x31,
                32 => x32,
                -1 => x_1,
                -2 => x_2,
                -3 => x_3,
                -4 => x_4,
                -5 => x_5,
                -6 => x_6,
                -8 => x_8,
                -9 => x_9,
                -11 => x_11,
                _ => x_val.powi(i), // 兜底使用 powi
            };
            
            // 根据 j 值选择对应的幂次
            let y_power = match j {
                0 => y0,
                1 => y1,
                2 => y2,
                3 => y3,
                4 => y4,
                5 => y5,
                6 => y6,
                10 => y10,
                17 => y17,
                -1 => y_1,
                -2 => y_2,
                -3 => y_3,
                -5 => y_5,
                -6 => y_6,
                -7 => y_7,
                -8 => y_8,
                -9 => y_9,
                -11 => y_11,
                -29 => y_29,
                -31 => y_31,
                -38 => y_38,
                -39 => y_39,
                -40 => y_40,
                -41 => y_41,
                _ => y_val.powi(j), // 兜底使用 powi
            };
            
            x_powers.push(x_power);
            y_powers.push(y_power);
        }
        
        (x_powers, y_powers)
    }};
}

//  Initialize coefficients and exponents for region 1
pub const IJn: [(i32, i32, f64); 34] = [
    (0, -2, 0.14632971213167E+00),
    (0, -1, -0.84548187169114E+00),
    (0, 0, -0.37563603672040E+01),
    (0, 1, 0.33855169168385E+01),
    (0, 2, -0.95791963387872E+00),
    (0, 3, 0.15772038513228E+00),
    (0, 4, -0.16616417199501E-01),
    (0, 5, 0.81214629983568E-03),
    (1, -9, 0.28319080123804E-03),
    (1, -7, -0.60706301565874E-03),
    (1, -1, -0.18990068218419E-01),
    (1, 0, -0.32529748770505E-01),
    (1, 1, -0.21841717175414E-01),
    (1, 3, -0.52838357969930E-04),
    (2, -3, -0.47184321073267E-03),
    (2, 0, -0.30001780793026E-03),
    (2, 1, 0.47661393906987E-04),
    (2, 3, -0.44141845330846E-05),
    (2, 17, -0.72694996297594E-15),
    (3, -4, -0.31679644845054E-04),
    (3, 0, -0.28270797985312E-05),
    (3, 6, -0.85205128120103E-09),
    (4, -5, -0.22425281908000E-05),
    (4, -2, -0.65171222895601E-06),
    (4, 10, -0.14341729937924E-12),
    (5, -8, -0.40516996860117E-06),
    (8, -11, -0.12734301741641E-08),
    (8, -6, -0.17424871230634E-09),
    (21, -29, -0.68762131295531E-18),
    (23, -31, 0.14478307828521E-19),
    (29, -38, 0.26335781662795E-22),
    (30, -39, -0.11947622640071E-22),
    (31, -40, 0.18228094581404E-23),
    (32, -41, -0.93537087292458E-25),
];

/// Fundamental equation for region 1
pub fn gamma_reg1(pi: f64, tau: f64) -> f64 {
    let (x_powers, y_powers) = precompute_powers!(pi, tau, &IJn);
     // x_powers[i] = x^i (第 i 个系数的 i 值对应的幂次)
    // y_powers[i] = y^j (第 i 个系数的 j 值对应的幂次)
    let steps: [(usize, usize); 2] = [(0, 19), (19, 34)];
    // poly_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps)
    poly_powi_steps_precomputed(&IJn, &x_powers, &y_powers, &steps)
}

/// First derivative of fundamental equation in pi for region 1
pub fn gamma_pi_reg1(pi: f64, tau: f64) -> f64 {
    let steps: [(usize, usize); 2] = [(0, 17), (17, 34)];
    -poly_i_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps)
}

/// Second derivative of fundamental equation in pi for region 1
pub fn gamma_pipi_reg1(pi: f64, tau: f64) -> f64 {
    let steps: [(usize, usize); 2] = [(0, 17), (17, 34)];
    poly_ii_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps)
}

/// First derivative of fundamental equation in tau for region 1
pub fn gamma_tau_reg1(pi: f64, tau: f64) -> f64 {
    let steps: [(usize, usize); 2] = [(0, 17), (17, 34)];
    poly_j_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps)
}

/// Second derivative of fundamental equation in tau for region 1
pub fn gamma_tautau_reg1(pi: f64, tau: f64) -> f64 {
    // let steps: [(usize, usize); 2] = [(0, 17), (17, 34)];
    let steps: [(usize, usize); 3] = [(0, 15), (15, 26), (26, 34)];
    poly_jj_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps)
}

/// Second derivative of fundamental equation in pi and tau for region 1
pub fn gamma_pitau_reg1(pi: f64, tau: f64) -> f64 {
    let steps: [(usize, usize); 2] = [(0, 17), (17, 34)];
    -poly_ij_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps)
}

// ------------------- multiple -------------------------------------
pub fn polys_i_j_powi_reg1(pi: f64, tau: f64) -> (f64, f64) {
    let steps: [(usize, usize); 3] = [(0, 16), (16, 26), (26, 34)];
    let (d_pi, d_tau) = polys_i_j_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps);
    (-d_pi, d_tau)
}

pub fn polys_i_ii_powi_reg1(pi: f64, tau: f64) -> (f64, f64) {
    let steps: [(usize, usize); 3] = [(0, 14), (14, 26), (26, 34)];
    let (poly_pi, poly_pipi) = polys_i_ii_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps);
    (-poly_pi, poly_pipi) // 7.1 - pi1,-> -d_pi
}

pub fn polys_0_j_powi_reg1(pi: f64, tau: f64) -> (f64, f64) {
    let steps: [(usize, usize); 2] = [(0, 16), (16, 34)];
    polys_0_j_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps)
}

pub fn polys_i_ij_powi_reg1(pi: f64, tau: f64) -> (f64, f64) {
    let steps: [(usize, usize); 2] = [(0, 17), (17, 34)];
    let (ploy_pi, ploy_pitau) = polys_i_ij_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps);
    (-ploy_pi, -ploy_pitau) // 7.1 - pi1,so -> -d_pi ,-d_pitau
}

/// Fast recursion algorithm
pub fn polys_i_ii_ij_jj_powi_reg1(pi: f64, tau: f64) -> (f64, f64, f64, f64) {
    let steps: [(usize, usize); 4] = [(0, 11), (11, 20), (20, 28), (28, 34)];
    let (ploy_pi, ploy_pipi, ploy_pitau, ploy_tautau) =
        polys_i_ii_ij_jj_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps);
    (-ploy_pi, ploy_pipi, -ploy_pitau, ploy_tautau)
}