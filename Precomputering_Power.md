# 预先计算幂性能分析

## 1. 概述

本文分析预计算幂次对IF97多项式求值的加速效果：**预计算幂性能不升反降**

预计算幂性能不升反降的真正原因：

1. 数组索引边界检查开销 - 每次访问 x_powers[k] 编译器要检查 k 是否越界
2. 预计算幂本身的开销 - 34项的预计算可能超过实际运行的收益
3. powi 编译器优化 - 标准库 powi 已经是高度优化的

直接使用原 powi 方法 （ polynomial_steps.rs ），它已经是最优的，因为：

1. 循环内无除法
2. powi 编译器高度优化
3. 无数组边界检查开销

## 2. 性能问题根源

### 2.1 除法开销分析

| 操作 | 相对开销 |
|------|----------|
| 乘法 | 1x |
| 除法 | 3-5x |
| powi (高次) | 5-10x |

### 2.2 原预计算方案的问题

原方案在导数计算中每项都进行除法：

```rust
// 原代码 - 循环内除法（慢）
value += coeff * x_powers[k] * (y_powers[k] / vj);  // 每项都要除
poly_j /= vj;  // 循环后还要除
```

## 3. 高性能预计算方案

### 3.1 核心优化：预计算倒数

**关键思想**：在循环外只做一次除法得到倒数，后续全部用乘法

```rust
// 优化后 - 循环外预计算倒数
let vj_inv = 1.0 / vj;  // 1次除法
let vi_inv = 1.0 / vi;

// 循环内全部用乘法（快）
for k in steps[m].0..steps[m].1 {
    value += coeff * x_powers[k] * y_powers[k] * vj_inv;  // 乘法
}
poly_j *= vj_inv;  // 乘法
```

### 3.2 Region 1 宏优化

**precompute_pi_powers 宏**：优化大整数幂次计算

```rust
macro_rules! precompute_pi_powers {
    ($x:ident) => {{
        let x0 = 1.0;
        let x1 = $x;
        let x2 = x1 * x1;
        let x3 = x2 * x1;
        let x4 = x2 * x2;
        let x5 = x4 * x1;
        let x8 = x4 * x4;
        // 优化：先计算 x16，减少乘法次数
        let x16 = x8 * x8;
        let x21 = x16 * x5;   // 2次乘法（原：x8*x8*x5 = 3次）
        let x23 = x21 * x2;   // 1次乘法
        let x29 = x21 * x8;   // 1次乘法
        let x30 = x29 * x1;
        let x31 = x30 * x1;
        let x32 = x31 * x1;
        // ... 数组
    }};
}
```

**precompute_tau_powers 宏**：只用一次除法

```rust
macro_rules! precompute_tau_powers {
    ($y:ident) => {{
        // 正幂次 - 全用乘法
        let y1 = $y;
        let y2 = $y * $y;
        // ... 其他正幂次

        // 关键优化：只做一次除法
        let y_inv = 1.0 / $y;      // 1次除法
        let y_1 = y_inv;            // 0次额外运算
        let y_2 = y_inv * y_inv;    // 1次乘法
        let y_3 = y_2 * y_inv;     // 1次乘法
        // ... 其余负幂次全用乘法
        let y41_inv = y40_inv * y_inv;  // 1次乘法

        // ... 数组
    }};
}
```

### 3.3 多项式函数优化

所有导数计算函数都已优化：

```rust
// poly_j_powi_steps_precomputed
pub fn poly_j_powi_steps_precomputed(...) -> f64 {
    let vj_inv = 1.0 / vj;  // 循环外1次除法
    let mut value: f64 = 0.0;
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            if IJn[k].1 != 0 {
                // 循环内全用乘法
                value += IJn[k].2 * x_powers[k] * IJn[k].1 as f64 * y_powers[k] * vj_inv;
            }
        }
    }
    value
}

// polys_0_j_powi_steps_precomputed
pub fn polys_0_j_powi_steps_precomputed(...) -> (f64, f64) {
    let vj_inv = 1.0 / vj;  // 预计算倒数
    let mut item: f64 = 0.0;
    let mut poly_0: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * x_powers[k] * y_powers[k];
            poly_0 += item;
            if IJn[k].1 != 0 {
                poly_j += IJn[k].1 as f64 * item;
            }
        }
    }
    poly_j *= vj_inv;  // 乘法代替除法
    (poly_0, poly_j)
}
```

## 4. 优化效果对比

| 函数 | 原方法除法次数 | 优化后除法次数 | 改进 |
|------|----------------|----------------|------|
| poly_j_powi_steps | 每项1次 | 1次（循环外） | 大幅减少 |
| polys_0_j_powi_steps | 每项+1次循环后 | 1次（循环外） | 大幅减少 |
| polys_i_j_powi_steps | 每项2次+循环后2次 | 2次（循环外） | 大幅减少 |
| polys_i_ii_ij_jj_powi_steps | 每项4次+循环后4次 | 4次（循环外） | 大幅减少 |

## 5. 进一步优化方向

### 5.1 避免数组边界检查

使用切片时，编译器可能无法证明索引不会越界。可以考虑：

```rust
// 使用裸指针访问（需 unsafe）
// 或确保 steps 是编译期常量
```

### 5.2 SIMD 向量化

当前循环结构可能阻碍自动向量化：

```rust
// 分离计算以利于向量化
let coeff = IJn[k].2;
let power = x_powers[k] * y_powers[k];
item = coeff * power;
```

### 5.3 循环展开

对关键路径进行手动循环展开：

```rust
// 对 steps[0] 展开
for k in 0..11 {
    // ... 计算
}
```

## 6. 结论

**主要优化**：
1. ✅ **消除循环内除法**：预计算 `vj_inv = 1.0 / vj`
2. ✅ **优化大幂次计算**：减少乘法次数
3. ✅ **循环外预计算**：所有导数函数都已优化

**性能提升**：消除每项的除法操作后，预计算方法的性能应显著提升。

**建议**：运行基准测试验证优化效果，对比原 `powi` 方法和预计算方法的实际性能。