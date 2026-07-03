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
from seuif97 import tx, pt, px
OH = 4

n3 = 0.10192970039326e-2
n4 = 0.57254459862746e3
n5 = 0.13918839778870e2
Ps_623 = 16.5291642526045

# Generate B23 data points from exact boundary
# Use px(p, x=1) for saturated vapor enthalpy at boundary
h_ref = px(Ps_623, 1.0, OH)            # reference h at p = Ps_623

# Start exactly from Ps_623
p_b23 = np.linspace(Ps_623, 100.0, 501)
T_b23 = np.array([n4 + ((p - n5) / n3) ** 0.5 for p in p_b23])
h_b23 = np.empty_like(p_b23)
h_b23[0] = h_ref                       # boundary: saturated vapor
for i in range(1, len(p_b23)):
    h_b23[i] = pt(p_b23[i], T_b23[i] - 273.15, OH)

# ============================================================
# 3-segment cubic polynomial fit
# ============================================================
# ============================================================
# Constrained cubic fit: force exact match at segment start
# ============================================================

def fit_segment_relative(p_data, h_data, p_fix, h_fix, degree=3):
    """Fit using relative pressure pr = p / Ps_623.
    h = a3*pr^3 + a2*pr^2 + a1*pr + a0
    Constraint: force exact match at pr_fix = p_fix / Ps_623
    """
    pr = p_data / Ps_623
    pr_fix = p_fix / Ps_623
    h = h_data
    A = np.column_stack([
        pr**3 - pr_fix**3,
        pr**2 - pr_fix**2,
        pr - pr_fix
    ])
    b = h - h_fix
    coeffs3, *_ = np.linalg.lstsq(A, b, rcond=None)
    a3, a2, a1 = coeffs3
    a0 = h_fix - a3*pr_fix**3 - a2*pr_fix**2 - a1*pr_fix
    coeffs = np.array([a3, a2, a1, a0])
    h_fit = np.polyval(coeffs, pr)
    max_err = np.max(np.abs(h_data - h_fit))
    return coeffs, max_err

def evaluate_segments_relative(segs):
    all_coeffs = []
    total_max_err = 0.0
    for p_lo, p_hi in segs:
        mask = (p_b23 >= p_lo) & (p_b23 <= p_hi)
        p_seg = p_b23[mask]
        h_seg = h_b23[mask]
        p_fix = p_seg[0]
        h_fix = h_seg[0]
        coeffs, max_err = fit_segment_relative(p_seg, h_seg, p_fix, h_fix, degree=3)
        all_coeffs.append(coeffs)
        total_max_err = max(total_max_err, max_err)
    return all_coeffs, total_max_err

# ============================================================
# Test different segment schemes
# ============================================================
PC_WATER = 22.064

def test_segments(name, segs):
    coeffs, err = evaluate_segments_relative(segs)
    print(f"\n{'='*60}")
    print(f"{name}")
    print(f"Overall max_abs_err = {err:.6f} kJ/kg")
    
    def h_eval(p):
        for i, (p_lo, p_hi) in enumerate(segs):
            if p_lo - 1e-9 <= p <= p_hi + 1e-9:
                return np.polyval(coeffs[i], p / Ps_623)
        return None
    
    print(f"{'p (MPa)':>10} {'h_true':>14} {'h_fit':>14} {'error':>12}")
    for p in [Ps_623, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]:
        if abs(p - Ps_623) < 0.001:
            h_true = px(Ps_623, 1.0, OH)
        else:
            T = n4 + ((p - n5) / n3) ** 0.5
            h_true = pt(p, T - 273.15, OH)
        h_fit = h_eval(p)
        print(f"{p:>10.4f} {h_true:>14.6f} {h_fit:>14.6f} {h_true-h_fit:>+12.6f}")
    
    return err, coeffs, segs

# Test 5-segment (original best)
best_err = float('inf')
best_b2 = 40.0
for b2 in np.linspace(30.0, 50.0, 101):
    segs = [(Ps_623, PC_WATER), (PC_WATER, b2), (b2, 60.0), (60.0, 80.0), (80.0, 100.0)]
    coeffs, err = evaluate_segments_relative(segs)
    if err < best_err:
        best_err = err
        best_b2 = b2
        best_coeffs = coeffs
        best_segs = segs

err5, coeffs5, segs5 = test_segments(
    f"5-segment: optimal b2={best_b2:.2f}",
    best_segs
)

all_coeffs, segments = coeffs5, segs5
print("  -> Using 5-segment scheme")

# Final detailed output
print(f"\n{'='*60}")
print("Final coefficients:")
for i, (p_lo, p_hi) in enumerate(segments):
    mask = (p_b23 >= p_lo) & (p_b23 <= p_hi)
    coeffs = all_coeffs[i]
    pr = p_b23[mask] / Ps_623
    h_fit = np.polyval(coeffs, pr)
    max_err = np.max(np.abs(h_b23[mask] - h_fit))
    rel_err = np.max(np.abs((h_b23[mask] - h_fit) / h_b23[mask])) * 100

    print(f"\nSegment {i+1}: {p_lo:.4f} ≤ p ≤ {p_hi:.1f} MPa")
    print(f"  h = a3*pr^3 + a2*pr^2 + a1*pr + a0,  pr = p / {Ps_623}")
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
print("Final verification at key points:")
print(f"{'p (MPa)':>10} {'h_true':>14} {'h_fit':>14} {'error':>12}")

def h_b23_eval(p):
    for i, (p_lo, p_hi) in enumerate(segments):
        if p_lo - 1e-9 <= p <= p_hi + 1e-9:
            return np.polyval(all_coeffs[i], p / Ps_623)
    return None

check_points = [Ps_623, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
for p in check_points:
    if abs(p - Ps_623) < 0.001:
        h_true = px(Ps_623, 1.0, OH)
    else:
       T = n4 + ((p - n5) / n3) ** 0.5
       h_true = pt(p, T - 273.15, OH)
    h_fit = h_b23_eval(p)
    err = h_true - h_fit
    print(f"{p:>10.4f} {h_true:>14.6f} {h_fit:>14.6f} {err:>+12.6f}")

# ============================================================
# Boundary continuity check
# ============================================================
print("\n" + "=" * 60)
print("Continuity at segment boundaries:")
for p_boundary in [p_hi for p_lo, p_hi in segments[:-1]]:
    h_left = h_b23_eval(p_boundary - 1e-10)
    h_right = h_b23_eval(p_boundary + 1e-10)
    print(f"  p={p_boundary} MPa: h_left={h_left:.6f}, h_right={h_right:.6f}, jump={h_right-h_left:.6f}")

# ============================================================
# Output Rust code for boundaries_supp.rs
# ============================================================
print("\n" + "=" * 60)
print("Rust code for boundaries_supp.rs:")
print("=" * 60)
print("// B23 boundary: h = f(p) using 5-segment cubic polynomials")
print("// Relative pressure: pr = p / Ps_623")
print(f"const PS_623: f64 = {Ps_623:.15e};")
print(f"const PC_WATER: f64 = {PC_WATER:.15e};")
print()
print("// Segment boundaries: [start, end] for each segment")
print("const B23_P_BOUNDARIES: [(f64, f64); 5] = [")
for i, (p_lo, p_hi) in enumerate(segments):
    print(f"    ({p_lo:.15e}, {p_hi:.15e}), // Segment {i+1}")
print("];")
print()
print("// Coefficients for h = a3*pr^3 + a2*pr^2 + a1*pr + a0")
print("const B23_COEFFS: [[f64; 4]; 5] = [")
for i, coeffs in enumerate(all_coeffs):
    print(f"    // Segment {i+1}: {segments[i][0]:.4f} ~ {segments[i][1]:.1f} MPa")
    print(f"    [{coeffs[0]:.15e}, {coeffs[1]:.15e}, {coeffs[2]:.15e}, {coeffs[3]:.15e}],")
print("];")
print()
print("/// B23 boundary enthalpy h = f(p) in kJ/kg")
print("/// p: pressure in MPa")
print("pub fn b23_p2h(p: f64) -> f64 {")
print("    let pr = p / PS_623;")
print("    for i in 0..5 {")
print("        if p >= B23_P_BOUNDARIES[i].0 && p <= B23_P_BOUNDARIES[i].1 {")
print("            let c = &B23_COEFFS[i];")
print("            return ((c[0]*pr + c[1])*pr + c[2])*pr + c[3];")
print("        }")
print("    }")
print("    f64::NAN")
print("}")