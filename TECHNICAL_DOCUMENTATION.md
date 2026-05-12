# SEUIF97 技术文档

## 1. 项目概述

SEUIF97 是基于Rust语言实现的高速IAPWS-IF97水和水蒸汽性质计算库。该库专为计算密集型任务设计，如非稳态过程模拟、在线过程监控和优化等场景。

### 1.1 核心性能优势

- **计算速度提升**：相比直接使用 Rust 标准库 `powi()` 的实现，Region 1、2、3 的基本方程计算速度提升 **5x 至 20x**
- **显著优于**各类近似方程和快速水蒸汽性质计算算法

### 1.2 支持的输入参数对

```txt
(p,t) (p,h) (p,s) (p,v)
(t,h) (t,s) (t,v)
(p,x) (t,x) (h,x) (s,x)
(h,s)
```

### 1.3 计算属性数量

支持 **36 种**热力学、传输和衍生属性的计算。

## 2. 技术架构

### 2.1 项目结构

```
RustSEUIF97/
├── src/
│   ├── algo/              # 核心算法模块
│   │   ├── polynomial.rs          # 多项式计算
│   │   ├── polynomial_steps.rs    # 分步骤多项式（加速核心）
│   │   └── root.rs                # 根求解算法
│   ├── common/            # 通用组件
│   │   ├── boundaries.rs          # 区域边界判定
│   │   ├── constant.rs            # 物理常数
│   │   ├── property_id.rs         # 属性ID定义
│   │   ├── property_pairs.rs      # 属性对计算
│   │   ├── region.rs              # 区域判断逻辑
│   │   └── transport_further.rs   # 传输属性计算
│   ├── r1/                # Region 1 (液态区)
│   ├── r2/                # Region 2 (过热蒸汽区)
│   ├── r3/                # Region 3 (临界点附近区)
│   ├── r4/                # Region 4 (饱和区)
│   ├── r5/                # Region 5 (高温区)
│   ├── cdecl_c_if97.rs    # C接口 (cdecl)
│   ├── stdcall_c_if97.rs  # C接口 (stdcall)
│   ├── python_if97.rs     # Python绑定
│   ├── rust_if97.rs       # Rust API
│   └── lib.rs             # 库入口
├── demo_using_lib/        # 多语言示例代码
├── dynamic_lib/           # 预编译动态库
├── benches/              # 性能基准测试
├── examples/              # Rust示例
└── tests/                 # 测试套件
```

### 2.2 模块职责说明

| 模块 | 职责 | 关键文件 |
|------|------|----------|
| **algo** | 核心数学算法实现 | polynomial_steps.rs, root.rs |
| **common** | 通用常量、边界判定、属性定义 | constant.rs, boundaries.rs, region.rs |
| **r1-r5** | 各热力学区域的属性计算 | regionX_pT.rs, regionX_ph_ps_hs.rs 等 |
| **cdecl_c_if97** | C语言接口 (cdecl调用约定) | cdecl_c_if97.rs |
| **stdcall_c_if97** | C语言接口 (stdcall调用约定) | stdcall_c_if97.rs |
| **python_if97** | Python绑定 | python_if97.rs |
| **rust_if97** | Rust高级API | rust_if97.rs |

### 2.3 热力学区域模块

项目实现了IAPWS-IF97标准的5个热力学区域，各区域模块结构如下：

#### Region 1 - 液态水区

**温度范围**：273.15K ≤ T ≤ 623.15K  
**压力范围**：p ≥ 饱和压力 (最高100MPa)

| 文件 | 职责 |
|------|------|
| `region1.rs` | Region 1主入口，属性分发 |
| `region1_pT.rs` | 基本方程：(p,T) → v,u,s,h,cp,cv,w |
| `region1_T_phps.rs` | 温度计算：(T,p/h/s) → p/h/s |
| `region1_p_hs.rs` | 后向方程：(p,h/s) → T |
| `region1_gfe.rs` | 吉布斯自由能方程实现 |
| `region1_pT_ext.rs` | 扩展属性计算 |
| `region1_pair_ext.rs` | 扩展输入对实现 |

#### Region 2 - 过热蒸汽区

**温度范围**：273.15K ≤ T ≤ 1073.15K  
**压力范围**：p < 饱和压力 (最高100MPa)

| 文件 | 职责 |
|------|------|
| `region2.rs` | Region 2主入口，属性分发 |
| `region2_pT.rs` | 基本方程：(p,T) → v,u,s,h,cp,cv,w |
| `region2_T_ph.rs` | 后向方程：(p,h) → T |
| `region2_T_ps.rs` | 后向方程：(p,s) → T |
| `region2_p_hs.rs` | 后向方程：(h,s) → p,T |
| `region2_gfe.rs` | 吉布斯自由能方程实现 |
| `region2_pT_ext.rs` | 扩展属性计算 |
| `region2_pair_ext.rs` | 扩展输入对实现 |

#### Region 3 - 临界点附近区

**温度范围**：623.15K < T ≤ 863.15K  
**压力范围**：16.529MPa ≤ p ≤ 100MPa

| 文件 | 职责 |
|------|------|
| `region3.rs` | Region 3主入口，属性分发 |
| `region3_Td.rs` | 基本方程：(T,d) → p,h,u,s,cp,cv,w |
| `region3_v_pT.rs` | 后向方程：(p,T) → v |
| `region3_v_subregion_pT.rs` | 子区域细分计算 |
| `region3_Tv_phps.rs` | 温度计算：(T,v/p/h/s) → p/h/s |
| `region3_p_hs.rs` | 后向方程：(p,h/s) → T,v |
| `region3_hfe.rs` | Helmholtz自由能方程实现 |
| `region3_Td_ext.rs` | 扩展属性计算 |
| `region3_pair_ext.rs` | 扩展输入对实现 |

#### Region 4 - 饱和区

**温度范围**：273.15K ≤ T ≤ 647.096K (临界温度)  
**压力范围**：0.000611MPa ≤ p ≤ 22.064MPa (临界压力)

| 文件 | 职责 |
|------|------|
| `region4.rs` | Region 4主入口，饱和属性计算 |
| `region4_sat_pT.rs` | 饱和线方程：p=f(T), T=f(p) |
| `region4_pTx.rs` | 饱和区属性：(p,T,x) → h,s,v |
| `region4_T_hs.rs` | 饱和区后向方程 |
| `region4_pair_ext.rs` | 扩展输入对实现 |

#### Region 5 - 高温区

**温度范围**：1073.15K < T ≤ 2273.15K  
**压力范围**：p ≤ 50MPa

| 文件 | 职责 |
|------|------|
| `region5.rs` | Region 5主入口，属性分发 |
| `region5_pT.rs` | 基本方程：(p,T) → v,u,s,h,cp,cv,w |
| `region5_ph_ps_hs.rs` | 后向方程实现 |
| `region5_gfe.rs` | 吉布斯自由能方程实现 |
| `region5_pT_ext.rs` | 扩展属性计算 |
| `region5_pair_ext.rs` | 扩展输入对实现 |

---

## 3. 核心算法

### 3.1 加速技术

SEUIF97 采用两种核心加速方法：

**设计原理**：
1. **循环拆分**：将单一循环拆分为多个步骤，提升缓存局部性和SIMD向量化潜力
2. **递归求值**：通过除以基值 vi/vj 直接导出导数值，计算多项式值及其导数，避免重复计算

#### 3.1.1 Loop Tiling 方法

将单个循环拆分为多个小步骤，充分释放编译器优化能力，超越单循环性能。

```rust
#[inline(always)]
pub fn polys_i_j_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    // Loop splitting for better cache locality or SIMD potential
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * vi.powi(IJn[k].0) * vj.powi(IJn[k].1);
            poly_i += IJn[k].0 as f64 * item;
            poly_j += IJn[k].1 as f64 * item;
        }
    }

    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}
```

#### 3.1.2 递归多项式求值方法（Recurrence Method for Multi-Polynomial Evaluation）

通过利用多项式与其导数之间的关系，只需直接计算单个多项式，其余值通过乘以或除以基值导出，消除冗余计算，显著提升计算性能。

##### 3.1.2.1 IAPWS-IF97 基本方程

Region 1 的基本方程基于**吉布斯自由能**：

$$\frac{g(p,T)}{RT} = \gamma(\pi,\tau) = \sum_{i=1}^{34} n_i (7.1-\pi)^{I_i} (\tau-1.222)^{J_i}$$

其中：
- $\pi = p/p^{*}$
- $\tau = T^{*}/T$

比内能的推导：

$$u = g - T \left( \frac{\partial g}{\partial T} \right)_p - p \left( \frac{\partial g}{\partial p} \right)_T$$

$$\frac{u(\pi, \tau)}{RT} = \tau \gamma_{\tau} - \pi \gamma_{\pi}$$

##### 3.1.2.2 实现代码

```rust
// --- Module: algo/polynomial_steps.rs ---
// 1. 优化内核：使用循环拆分（steps）和强制内联
#[inline(always)]
pub fn polys_i_j_powi_steps(vi: f64, vj: f64, IJn: &[(i32, i32, f64)], steps: &[(usize, usize)]) -> (f64, f64) {
    let mut item: f64 = 0.0;
    let mut poly_i: f64 = 0.0;
    let mut poly_j: f64 = 0.0;

    // Loop splitting for better cache locality or SIMD potential
    for m in 0..steps.len() {
        for k in steps[m].0..steps[m].1 {
            item = IJn[k].2 * vi.powi(IJn[k].0) * vj.powi(IJn[k].1);
            poly_i += IJn[k].0 as f64 * item;
            poly_j += IJn[k].1 as f64 * item;
        }
    }

    // Multi-polynomial evaluation derived via base scaling (multiplication/division)
    poly_i /= vi;
    poly_j /= vj;
    (poly_i, poly_j)
}

// --- Module: r1/region1_gfe.rs ---
// 2. Region 1 包装器：处理特定坐标变换 (pi, tau)
pub fn polys_i_j_powi_reg1(pi: f64, tau: f64) -> (f64, f64) {
    // 显式定义计算步骤以辅助编译器优化
    let steps: [(usize, usize); 3] = [(0, 16), (16, 26), (26, 34)];
    let (d_pi, d_tau) = polys_i_j_powi_steps(7.1 - pi, tau - 1.222, &IJn, &steps);
    (-d_pi, d_tau)
}

// --- Module: r1/region1_pT.rs ---
// 3. API：计算比内能
pub fn pT2u_reg1(p: f64, T: f64) -> f64 {
    let pi: f64 = p / r1pstar;
    let tau: f64 = r1Tstar / T;
    let (d_pi, d_tau) = polys_i_j_powi_reg1(pi, tau);
    RGAS_WATER * T * (tau * d_tau - pi * d_pi)
}
```

### 3.2 区域判定逻辑

区域判定是计算流程的关键环节，根据输入参数对选择不同的判定策略：

**判定流程**：
1. 参数边界校验
2. 饱和线检测（Region 4）
3. 根据温度压力范围判定具体区域

```rust
pub fn pT_sub_region(p: f64, T: f64) -> i32 {
    // 1. 边界检查
    if p < P_MIN || p > 100.0 { return INVALID_P; }
    if T < 273.15 || T > 2273.15 { return INVALID_T; }
    
    // 2. 饱和线检查
    if T >= 273.15 && T < TC_WATER {
        let ps: f64 = p_saturation(T);
        if (p - ps).abs() / ps < psatTol { return 4; }
    }
    
    // 3. 区域判定
    if T >= 273.15 && T <= 623.15 {
        if p >= p_saturation(T) && p <= 100.0 { return 1; }
        if p < p_saturation(T) && p > P_MIN { return 2; }
    }
    // ... 其他区域判定逻辑
}
```

## 4. API 接口设计

### 4.1 Rust API

#### 4.1.1 函数签名

```rust
pub fn pt<R>(p: f64, t: f64, o_id_reg: R) -> f64
where
    R: Into<o_id_region_args>,
```

**参数说明**：
- `p`: 压力 (MPa)
- `t`: 温度 (°C)
- `o_id_reg`: 属性ID（可选带区域指定）

#### 4.1.2 输入参数对函数

| 函数 | 输入参数 | 说明 |
|------|----------|------|
| `pt(p, t, o_id)` | 压力, 温度 | 最常用的输入方式 |
| `ph(p, h, o_id)` | 压力, 焓 | 适用于焓已知的场景 |
| `ps(p, s, o_id)` | 压力, 熵 | 适用于熵已知的场景 |
| `pv(p, v, o_id)` | 压力, 比容 | 扩展输入对 |
| `th(t, h, o_id)` | 温度, 焓 | 扩展输入对 |
| `ts(t, s, o_id)` | 温度, 熵 | 扩展输入对 |
| `tv(t, v, o_id)` | 温度, 比容 | 扩展输入对 |
| `hs(h, s, o_id)` | 焓, 熵 | 常用热力过程输入 |
| `px(p, x, o_id)` | 压力, 干度 | 湿蒸汽区专用 |
| `tx(t, x, o_id)` | 温度, 干度 | 湿蒸汽区专用 |
| `hx(h, x, o_id)` | 焓, 干度 | 湿蒸汽区专用 |
| `sx(s, x, o_id)` | 熵, 干度 | 湿蒸汽区专用 |

#### 4.1.3 属性ID常量

| 属性 | 符号 | o_id | 单位 |
|------|------|------|------|
| Pressure | p | OP(0) | MPa |
| Temperature | t | OT(1) | °C |
| Density | ρ | OD(2) | kg/m³ |
| Specific Volume | v | OV(3) | m³/kg |
| Specific enthalpy | h | OH(4) | kJ/kg |
| Specific entropy | s | OS(5) | kJ/(kg·K) |
| Specific exergy | e | OE(6) | kJ/kg |
| Specific internal energy | u | OU(7) | kJ/kg |
| Isobaric heat capacity | cp | OCP(8) | kJ/(kg·K) |
| Isochoric heat capacity | cv | OCV(9) | kJ/(kg·K) |
| Speed of sound | w | OW(10) | m/s |
| Isentropic exponent | k | OKS(11) | - |
| Helmholtz free energy | f | OF(12) | kJ/kg |
| Gibbs free energy | g | OG(13) | kJ/kg |
| Compressibility factor | z | OZ(14) | - |
| Steam quality | x | OX(15) | - |
| Region | r | OR(16) | - |
| Isobaric cubic expansion coefficient | αv | OEC(17) | 1/K |
| Isothermal compressibility | kT | OKT(18) | 1/MPa |
| Partial derivative (∂V/∂T)p | - | ODVDT(19) | m³/(kg·K) |
| Partial derivative (∂V/∂p)T | - | ODVDP(20) | m³/(kg·MPa) |
| Partial derivative (∂P/∂T)v | - | ODPDT(21) | MPa/K |
| Isothermal throttling coefficient | δt | OIJTC(22) | kJ/(kg·MPa) |
| Joule-Thomson coefficient | μ | OJTC(23) | K/MPa |
| Dynamic viscosity | η | ODV(24) | Pa·s |
| Kinematic viscosity | ν | OKV(25) | m²/s |
| Thermal conductivity | λ | OTC(26) | W/(m·K) |
| Thermal diffusivity | a | OTD(27) | m²/s |
| Prandtl number | Pr | OPR(28) | - |
| Surface tension | σ | OST(29) | N/m |
| Static Dielectric Constant | ε | OSDC(30) | - |
| Isochoric pressure coefficient | β | OPC(31) | 1/K |
| Isothermal stress coefficient | βp | OBETAP(32) | kg/m³ |
| Fugacity coefficient | fi | OFI(33) | - |
| Fugacity | f* | OFU(34) | MPa |
| Relative pressure coefficient | αp | OALFAP(35) | 1/K |

#### 4.1.4 使用示例

```rust
use seuif97::*;

fn main() {
    let p: f64 = 16.0;    // MPa
    let t: f64 = 535.1;   // °C
    
    // 计算焓值
    let h = pt(p, t, OH);
    println!("h = {:.3f} kJ/kg", h);
    
    // 指定区域计算熵
    let s = pt(p, t, (OS, 1));
    println!("s = {:.3f} kJ/(kg·K)", s);
}
```

### 4.2 C 语言接口

#### 4.2.1 编译选项

```bash
# cdecl 调用约定
cargo build -r --features cdecl

# stdcall 调用约定 (Windows 64位)
cargo build -r --features stdcall

# stdcall 调用约定 (Windows 32位)
cargo build -r --target=i686-pc-windows-msvc --features stdcall
```

#### 4.2.2 函数声明

```c
double pt(double p, double t, short o_id);
double ph(double p, double h, short o_id);
double ps(double p, double s, short o_id);
double pv(double p, double v, short o_id);
double tv(double t, double v, short o_id);
double th(double t, double h, short o_id);
double ts(double t, double s, short o_id);
double hs(double h, double s, short o_id);
double px(double p, double x, short o_id);
double tx(double t, double x, short o_id);
double hx(double h, double x, short o_id);
double sx(double s, double x, short o_id);
```

#### 4.2.3 使用示例

```c
#include <stdio.h>

#define OH 4
#define OS 5

extern double pt(double p, double t, short o_id);

int main(void) {
    double p = 16.0;
    double t = 530.0;
    double h = pt(p, t, OH);
    double s = pt(p, t, OS);
    printf("p,t %f,%f h= %f s= %f\n", p, t, h, s);
    return 0;
}
```

### 4.3 Python 接口

#### 4.3.1 安装方式

```bash
# 从 PyPI 安装
pip install seuif97

# 本地构建安装
python ./setup.py install
```

#### 4.3.2 使用示例

```python
from seuif97 import *

OH = 4
p = 16.0
t = 535.1

# 使用属性ID计算
h = pt(p, t, OH)

# 直接计算特定属性
s = pt2s(p, t)

print(f"p={p}, t={t} h={h:.3f} s={s:.3f}")
```

#### 4.3.3 Python 快捷函数

| 函数 | 说明 |
|------|------|
| `pt2h(p, t)` | 从(p,t)计算焓 |
| `pt2s(p, t)` | 从(p,t)计算熵 |
| `pt2v(p, t)` | 从(p,t)计算比容 |
| `ph2t(p, h)` | 从(p,h)计算温度 |
| `ps2t(p, s)` | 从(p,s)计算温度 |
| `hs2p(h, s)` | 从(h,s)计算压力 |
| `hs2t(h, s)` | 从(h,s)计算温度 |

---

## 5. 热力学区域说明

### 5.1 区域划分

| 区域 | 范围 | 状态 |
|------|------|------|
| **Region 1** | 273.15K ≤ T ≤ 623.15K, p ≥ 饱和压力 | 液态水 |
| **Region 2** | T ≥ 273.15K, p < 饱和压力 | 过热蒸汽 |
| **Region 3** | 623.15K < T ≤ 863.15K, 高压区 | 临界点附近 |
| **Region 4** | 273.15K ≤ T ≤ 647.096K | 饱和两相区 |
| **Region 5** | 1073.15K < T ≤ 2273.15K, p ≤ 50MPa | 高温区 |

### 5.2 区域边界

```
                    T (K)
                     ^
              2273.15|          Region 5
                     |            (p ≤ 50MPa)
                     |         +------------------+
                     |         |                  |
              1073.15|    +----+                  +----+
                     |    |    |                  |    |
                     |    |    |   Region 2       |    |
                     |    |    |                  |    |
               863.15|    |    |    +----------+  |    |
                     |    |    |    | Region 3  |  |    |
               647.10|    |    |    |(critical)|  |    |
                     |    |    |    +----------+  |    |
               623.15|    |    +------------------+    |
                     |    |        ^                  |
                     |    |   Region 4                |
                     |    |  (saturation)             |
               273.15|----+--------+------------------+----> p (MPa)
                     |   Region 1  |
                     |  (liquid)   |  Region 2
                     +-------------+ (vapor)
```

---

## 6. 物理常量定义

### 6.1 关键常量

```rust
pub const K: f64 = 273.15;                    // 摄氏温度转换常数
pub const RGAS_WATER: f64 = 0.461526;         // 气体常数 kJ/(kg·K)

// 临界点参数
pub const TC_WATER: f64 = 647.096;            // 临界温度 K
pub const PC_WATER: f64 = 22.064;             // 临界压力 MPa
pub const DC_WATER: f64 = 322.0;              // 临界密度 kg/m³

// 边界常量
pub const P_MIN: f64 = 0.000611212677444;     // 最小压力 MPa
pub const P_MAX1: f64 = 100.0;                // Region 1 最大压力
pub const T_MAX1: f64 = 623.15;               // Region 1 最大温度
```

---

## 7. 错误处理机制

### 7.1 错误码定义

| 错误码 | 常量 | 含义 |
|--------|------|------|
| -9999 | `INVALID_VALUE` | 无效值 |
| -1000 | `INVALID_OUTID` | 无效属性ID |
| -2100 | `INVALID_P` | 无效压力 |
| -2101 | `INVALID_T` | 无效温度 |
| -2102 | `INVALID_S` | 无效熵 |
| -2103 | `INVALID_H` | 无效焓 |
| -2201 | `INVALID_PT` | 无效(p,T)组合 |
| -2202 | `INVALID_HS` | 无效(h,s)组合 |

### 7.2 输入验证流程

```
输入参数 → 范围检查 → 区域判定 → 属性计算 → 返回结果
              ↓
         超出范围?
              ↓
         返回错误码
```

---

## 8. 构建与测试

### 8.1 构建命令

```bash
# 开发构建
cargo build

# 发布构建 (优化)
cargo build --release

# 构建 Python 扩展
cargo build --release --features python

# 运行测试
cargo test

# 运行基准测试
cargo bench
```

### 8.2 测试套件

项目包含全面的测试用例：

| 测试文件 | 测试内容 |
|----------|----------|
| `pt_test.rs` | (p,T)输入对测试 |
| `ph_test.rs` | (p,h)输入对测试 |
| `ps_test.rs` | (p,s)输入对测试 |
| `pv_test.rs` | (p,v)输入对测试 |
| `th_test.rs` | (t,h)输入对测试 |
| `ts_test.rs` | (t,s)输入对测试 |
| `tv_test.rs` | (t,v)输入对测试 |
| `hs_test.rs` | (h,s)输入对测试 |
| `hxsx_test.rs` | 湿蒸汽区测试 |
| `cross_test.rs` | 跨区域边界测试 |

### 8.3 性能基准测试

项目使用 [Criterion](https://crates.io/crates/criterion) 框架进行性能基准测试，测试文件位于 `benches/speed_benchmark.rs`。

#### 8.3.1 测试内容

性能基准测试覆盖了各热力学区域的典型计算场景：

| 测试函数 | 区域 | 输入参数 | 计算属性 |
|----------|------|----------|----------|
| `pt2h_reg1` | Region 1 | p=3.0MPa, t=26.85°C | 比焓 |
| `pt2s_reg1` | Region 1 | p=3.0MPa, t=26.85°C | 比熵 |
| `pt2h_reg2` | Region 2 | p=0.0035MPa, t=26.85°C | 比焓 |
| `pt2s_reg2` | Region 2 | p=0.0035MPa, t=26.85°C | 比熵 |
| `tv2h_reg3` | Region 3 | t=376.85°C, v=0.002 m³/kg | 比焓 |
| `tv2s_reg3` | Region 3 | t=376.85°C, v=0.002 m³/kg | 比熵 |
| `pT2h_reg5` | Region 5 | p=0.5MPa, t=1226.85°C | 比焓 |
| `pT2s_reg5` | Region 5 | p=0.5MPa, t=1226.85°C | 比熵 |

#### 8.3.2 测试代码实现

```rust
use criterion::{black_box, criterion_group, criterion_main, Criterion};

use seuif97::*;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("pt2h_reg1", |b| {
        b.iter(|| pt(black_box(3.0), black_box(300.0 - 273.15), black_box(OH)))
    });
    c.bench_function("pt2s_reg1", |b| {
        b.iter(|| pt(black_box(3.0), black_box(300.0 - 273.15), black_box(OS)))
    });
    // ... 其他区域测试
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
```

#### 8.3.3 运行性能测试

```bash
# 运行所有基准测试
cargo bench

# 查看生成的HTML报告
# 报告位于 target/criterion/index.html
```

#### 8.3.4 性能测试配置

在 `Cargo.toml` 中配置了 Criterion 依赖和基准测试设置：

```toml
[dev-dependencies]
criterion = { version = "0.5.1", features = ["html_reports"] }

[[bench]]
name = "speed_benchmark"
harness = false
```

**配置说明**：
- `features = ["html_reports"]`：生成HTML格式的性能报告
- `harness = false`：使用Criterion的基准测试接口而非内置harness

---

## 9. 部署与集成

### 9.1 动态库部署

预编译的动态库位于 `dynamic_lib/` 目录：

| 平台 | 文件 | 路径 |
|------|------|------|
| Windows 64位 | `seuif97.dll` | `dynamic_lib/windows_x64/` |
| Windows 32位 | `seuif97.dll` | `dynamic_lib/windows_x86/` |
| Linux 64位 | `libseuif97.so` | `dynamic_lib/linux_x64/` |

### 9.2 多语言集成示例

支持的编程语言：
- ✅ Rust
- ✅ C / C++
- ✅ Python
- ✅ C#
- ✅ Java
- ✅ Fortran
- ✅ Go
- ✅ Excel VBA

---

## 10. 性能优化建议

### 10.1 使用建议

1. **批量计算**：对于大量计算任务，建议使用循环分块技术，充分利用缓存
2. **区域预判定**：如果已知计算区域，直接指定区域参数可避免区域判定开销
3. **避免重复计算**：对于相同输入参数对的多次查询，考虑缓存结果

### 10.2 性能对比

| 实现方式 | 性能 | 说明 |
|----------|------|------|
| SEUIF97 | 基准 | 优化后的高速实现 |
| Rust标准库 `powi()` | 慢 5-20x | 无优化的直接实现 |
| 其他近似算法 | 精度损失 | 速度可能更快但精度不足 |

---

## 11. 维护与贡献

### 11.1 代码规范

- 使用 Rust 2021 edition
- 遵循 `rustfmt` 代码格式化规则
- 使用 `clippy` 进行代码检查

### 11.2 贡献流程

1. Fork 仓库
2. 创建特性分支
3. 编写代码和测试
4. 运行测试确保通过
5. 提交 Pull Request

### 11.3 版本管理

版本格式：`MAJOR.MINOR.PATCH`

- **MAJOR**：API 不兼容变更
- **MINOR**：新增功能，向后兼容
- **PATCH**：Bug 修复，向后兼容

---

## 12. 参考文献

1. IAPWS-IF97: "Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam"
2. IAPWS Supplementary Release: "Supplementary Release on Backward Equations for the Properties of Water and Steam"
3. IAPWS Supp-Tv(ph,ps)-2014: "Supplementary Release for the Region 3 Boundaries"
4. IAPWS Supp-phs3-2014: "Supplementary Release for the (h,s) Region Boundaries"

---

## 附录：单位转换表

| 物理量 | SI单位 | 工程单位 | 转换关系 |
|--------|--------|----------|----------|
| 压力 | Pa | MPa | 1 MPa = 10^6 Pa |
| 温度 | K | °C | T(K) = t(°C) + 273.15 |
| 焓 | J/kg | kJ/kg | 1 kJ/kg = 10^3 J/kg |
| 熵 | J/(kg·K) | kJ/(kg·K) | 1 kJ/(kg·K) = 10^3 J/(kg·K) |
| 比容 | m³/kg | m³/kg | - |
| 密度 | kg/m³ | kg/m³ | - |

---

**文档版本**: v1.2.2  
**生成日期**: 2024年  
**作者**: Cheng Maohua <cmh@seu.edu.cn>  
**项目地址**: https://github.com/thermalogic/RustSEUIF97