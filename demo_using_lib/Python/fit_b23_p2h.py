"""
Fit B23 boundary curve h = f(p) as 3-segment cubic polynomials.

B23 curve: IAPWS-IF97 Eq.(6) T = n4 + sqrt((p - n5) / n3)
Range: 623.15K (16.5292 MPa) to 863.15K (100 MPa)

Segments:
  1: 16.529 ≤ p ≤ 30 MPa  (covers non-monotonic region)
  2: 30 ≤ p ≤ 60 MPa
  3: 60 ≤ p ≤ 100 MPa
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.abspath(r'D:\sim_rankine_workbench\iapws-if97\RustSEUIF97\target\release'))
from seuif97 import tx, pt
OH = 4

n3 = 0.10192970039326e-2
n4 = 0.57254459862746e3
n5 = 0.13918839778870e2
Ps_623 = 16.5291642526045

# Generate B23 data points
p_b23 = np.linspace(Ps_623 + 0.001, 100.0, 500)
T_b23 = np.array([n4 + ((p - n5) / n3) ** 0.5 for p in p_b23])
h_b23 = np.array([pt(p, T - 273.15, OH) for p, T in zip(p_b23, T_b23)])

# Prepend start point from saturated vapor line
h_start = tx(623.15 - 273.15, 1.0, OH)
p_b23 = np.concatenate(([Ps_623], p_b23))
T_b23 = np.concatenate(([623.15], T_b23))
h_b23 = np.concatenate(([h_start], h_b23))

# ============================================================
# 3-segment cubic polynomial fit
# ============================================================
segments = [
    (16.5291642526045, 30.0),
    (30.0, 60.0),
    (60.0, 100.0),
]

print("B23 curve h = f(p) — 3-segment cubic polynomial fit")
print("=" * 60)

all_coeffs = []
for i, (p_lo, p_hi) in enumerate(segments):
    mask = (p_b23 >= p_lo) & (p_b23 <= p_hi)
    p_seg = p_b23[mask]
    h_seg = h_b23[mask]
    coeffs = np.polyfit(p_seg, h_seg, 3)
    all_coeffs.append(coeffs)

    h_fit = np.polyval(coeffs, p_seg)
    max_err = np.max(np.abs(h_seg - h_fit))
    rel_err = np.max(np.abs((h_seg - h_fit) / h_seg)) * 100

    print(f"\nSegment {i+1}: {p_lo:.4f} ≤ p ≤ {p_hi:.1f} MPa")
    print(f"  h = a3*p^3 + a2*p^2 + a1*p + a0")
    print(f"  a3 = {coeffs[0]:.15e}")
    print(f"  a2 = {coeffs[1]:.15e}")
    print(f"  a1 = {coeffs[2]:.15e}")
    print(f"  a0 = {coeffs[3]:.15e}")
    print(f"  max_abs_err = {max_err:.6f} kJ/kg")
    print(f"  max_rel_err = {rel_err:.6f}%")

# ============================================================
# Verification at key points
# ============================================================
print("\n" + "=" * 60)
print("Verification at key points:")
print(f"{'p (MPa)':>10} {'h_true':>14} {'h_fit':>14} {'error':>12}")

def h_b23_eval(p):
    for i, (p_lo, p_hi) in enumerate(segments):
        if p_lo - 1e-9 <= p <= p_hi + 1e-9:
            return np.polyval(all_coeffs[i], p)
    return None

check_points = [16.5292, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
for p in check_points:
    T = n4 + ((p - n5) / n3) ** 0.5
    # At p=Ps_623, pt() returns liquid-side value; use tx() for vapor-side
    if abs(p - Ps_623) < 0.001:
        h_true = tx(623.15 - 273.15, 1.0, OH)
    else:
        h_true = pt(p, T - 273.15, OH)
    h_fit = h_b23_eval(p)
    err = h_true - h_fit
    print(f"{p:>10.4f} {h_true:>14.6f} {h_fit:>14.6f} {err:>+12.6f}")

# ============================================================
# Boundary continuity check
# ============================================================
print("\n" + "=" * 60)
print("Continuity at segment boundaries:")
for p_boundary in [30.0, 60.0]:
    h_left = h_b23_eval(p_boundary - 1e-10)
    h_right = h_b23_eval(p_boundary + 1e-10)
    print(f"  p={p_boundary} MPa: h_left={h_left:.6f}, h_right={h_right:.6f}, jump={h_right-h_left:.6f}")
