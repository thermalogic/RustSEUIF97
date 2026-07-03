"""
B142边界: h = f(p)，5段多项式，端点约束，h归一化
Segments: (P_MIN,0.1), (0.1,5.0), (5.0,Ps_623), (Ps_623,20.0), (20.0,PC_WATER)
Fit form: h/HC_WATER = a3*(x-x0)^3 + a2*(x-x0)^2 + a1*(x-x0) + a0
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.abspath(r'D:\sim_rankine_workbench\iapws-if97\RustSEUIF97\target\release'))
from seuif97 import tx, px, pt

OP = 0; OH = 4
P_MIN = 0.000611212677444
Ps_623 = 16.5291642526045
PC_WATER = 22.064
HC_WATER = 2.087546845e+03

Tcrit_K = 647.096
Tc = Tcrit_K - 273.15
T_sat = np.linspace(0.01, Tc, 800)

p_liq = np.array([tx(t, 0.0, OP) for t in T_sat])
h_liq = np.array([tx(t, 0.0, OH) for t in T_sat])
p_vap = np.array([tx(t, 1.0, OP) for t in T_sat])
h_vap = np.array([tx(t, 1.0, OH) for t in T_sat])

def fit_segment_endpoints(x_data, y_data, x0, y0, x1, y1):
    """用(x-x0)基，左端点y0固定，然后调整a1使右端=y1"""
    dx = x_data - x0
    dx1 = x1 - x0
    A = np.column_stack([dx**3, dx**2, dx])
    b = y_data - y0
    coeffs, *_ = np.linalg.lstsq(A, b, rcond=None)
    a3, a2, a1 = coeffs
    unconstrained_val = a3*dx1**3 + a2*dx1**2 + a1*dx1
    a1_adj = a1 + (y1 - y0 - unconstrained_val) / dx1 if dx1 != 0 else a1
    coeffs_final = np.array([a3, a2, a1_adj, y0])
    h_fit_norm = np.polyval(coeffs_final, dx)
    max_err = np.max(np.abs(y_data - h_fit_norm))
    return coeffs_final, max_err

def h_at_p(p, p_data, h_data):
    if p <= p_data[0]: return h_data[0]
    if p >= p_data[-1]: return h_data[-1]
    idx = np.searchsorted(p_data, p)
    if idx >= len(p_data): idx = len(p_data) - 1
    if idx == 0: idx = 1
    t_ratio = (p - p_data[idx-1]) / (p_data[idx] - p_data[idx-1])
    return h_data[idx-1] + t_ratio * (h_data[idx] - h_data[idx-1])

def fit_line(p_data, h_data, seg_bounds):
    """Fit h(p) with segments. Returns list of (p0, p1, x0, xform, coeffs, err).
    Shared boundary h-values ensure continuity (jump=0 at boundaries)."""
    seg_coeffs = []
    total_max_err = 0.0
    for i in range(len(seg_bounds) - 1):
        p0 = seg_bounds[i]
        p1 = seg_bounds[i+1]
        mask = (p_data >= p0) & (p_data <= p1)
        n_points = int(np.sum(mask))
        if n_points < 4:  # cubic needs at least 4 points, endpoints may have only on a short segment, skip if too few
            continue
        p_seg = p_data[mask]
        hn_seg = h_data[mask] / HC_WATER
        y0 = h_at_p(p0, p_data, h_data) / HC_WATER
        y1 = h_at_p(p1, p_data, h_data) / HC_WATER
        if p1 <= Ps_623:
            x_data = np.log10(p_seg)
            x0 = np.log10(p0)
            x1_val = np.log10(p1)
            xform = "log10"
        else:
            x_data = p_seg / Ps_623
            x0 = p0 / Ps_623
            x1_val = p1 / Ps_623
            xform = "p/Ps_623"
        coeffs, err = fit_segment_endpoints(x_data, hn_seg, x0, y0, x1_val, y1)
        seg_coeffs.append((p0, p1, x0, xform, coeffs, err * HC_WATER))
        total_max_err = max(total_max_err, err * HC_WATER)
    return seg_coeffs, total_max_err

def h_eval(p, seg_coeffs):
    """Evaluate h(p) by finding right segment and evaluating."""
    for (p0, p1, x0, xform, coeffs, _err) in seg_coeffs:
        if p0 <= p <= p1:
            if xform == "log10":
                x = np.log10(p)
            else:
                x = p / Ps_623
            dx = x - x0
            hn = np.polyval(coeffs, dx)
            return hn * HC_WATER
    if p < seg_coeffs[0][0]:
        _, _, x0, xform, coeffs, _ = seg_coeffs[0]
        x = np.log10(p) if xform == "log10" else p / Ps_623
        return np.polyval(coeffs, x - x0) * HC_WATER
    else:
        _, _, x0, xform, coeffs, _ = seg_coeffs[-1]
        x = np.log10(p) if xform == "log10" else p / Ps_623
        return np.polyval(coeffs, x - x0) * HC_WATER

# ========== 主程序 ==========
seg_bounds = [P_MIN, 0.1, 5.0, Ps_623, 20.0, 21.5, 21.8, 22.0, PC_WATER]

coeffs_liq, err_liq = fit_line(p_liq, h_liq, seg_bounds)
coeffs_vap, err_vap = fit_line(p_vap, h_vap, seg_bounds)

# 打印每段
for name, sc, err_total in [("液体", coeffs_liq, err_liq), ("蒸汽", coeffs_vap, err_vap)]:
    print(f"\n=== {name}线  Max err: {err_total:.3f} kJ/kg ===")
    for i, (p0, p1, x0, xform, coeffs, err) in enumerate(sc):
        print(f"  段{i+1}: {p0:.5f}~{p1:.5f} MPa, x={xform}, x0={x0:.6f}, err={err:.3f} kJ/kg")
        print(f"    a3={coeffs[0]:.8e}, a2={coeffs[1]:.8e}, a1={coeffs[2]:.8e}, a0={coeffs[3]:.8e}")

# 验证
for name, p_data, h_data, sc in [("Liquid", p_liq, h_liq, coeffs_liq), ("Vapor", p_vap, h_vap, coeffs_vap)]:
    print(f"\n=== Verification - Saturated {name} ===")
    print(f"{'p(MPa)':>10} {'h_true':>12} {'h_fit':>12} {'error':>12}")
    for p in [P_MIN, 0.001, 0.01, 0.1, 0.5, 1.0, 3.0, 5.0, 8.0, 10.0, 13.0, Ps_623, 18.0, 19.0, 20.0, 21.0, 21.5, 21.8, 22.0, PC_WATER]:
        if p > p_data[-1]: continue
        h_t = h_at_p(p, p_data, h_data)
        h_f = h_eval(p, sc)
        print(f"{p:>10.6f} {h_t:>12.4f} {h_f:>12.4f} {h_t-h_f:>+12.4f}")

# 连续性检查
print("\n=== Continuity ===")
for pb in seg_bounds[1:-1]:
    hl = h_eval(pb - 1e-10, coeffs_liq)
    hr = h_eval(pb + 1e-10, coeffs_liq)
    print(f"p={pb:.6f}: 液体 jump={hr-hl:.6f}")
for pb in seg_bounds[1:-1]:
    hl = h_eval(pb - 1e-10, coeffs_vap)
    hr = h_eval(pb + 1e-10, coeffs_vap)
    print(f"p={pb:.6f}: 蒸汽 jump={hr-hl:.6f}")

# Debug: 每段数据点数
print(f"\n=== Data point counts ===")
for i in range(len(seg_bounds)-1):
    p0, p1 = seg_bounds[i], seg_bounds[i+1]
    n_liq = int(np.sum((p_liq >= p0) & (p_liq <= p1)))
    n_vap = int(np.sum((p_vap >= p0) & (p_vap <= p1)))
    print(f"  {p0:.5f}~{p1:.5f}: 液体 {n_liq}点, 蒸汽 {n_vap}点")